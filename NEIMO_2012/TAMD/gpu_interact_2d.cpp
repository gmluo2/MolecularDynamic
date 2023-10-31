#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "def.h"

using namespace std;

#include "ranlib.h"

#include "show.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"

#include "ZMatrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "cluster.h"

#include "var.h"

#define INTERACT_DEF 1
#define _INTERACTION 0
#include "Interaction.h"
#include "Interact1.h"
#include "Interact2.h"

#include <fftw3.h>
#include "complex.h"

#include "NEIMO.h"
#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"
#include "MD_POLYMER.h"
#include "MD_SM.h"
#include "md_cell.h"

#include "vmd_cell.h"

//#include "EwaldSumVar.h"
#include "EwaldSum.h"
#include "spme.h"
#include "spme_2d.h"
#include "spme_interact.h"
#include "spme_interact_2d.h"

using namespace _evatom_;
using namespace _EwaldSum_real_;
using namespace _EwaldSum_2d_;

#include "atom_sr.h"
#include "atom_sr_interact.h"


#include "gpu_vector.h"
#include "gpu_spme.h"
#include "gpu_spme_2d.h"
#include "gpu_interact.h"
#include "gpu_interact_2d.h"

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);

extern double max_diff(double *a1, double *a2, int n);

#ifdef __GPU__

namespace _gpu_md_ {
	bool setup_gpu_spme(_gpu_md_::GPU_SPME_2d *spme_gpu, BSplineFunc<2> *bsp, bool bAutoExpandZ, float* xl, float *xd, int *dim_spme, double kappa, int natoms, int mMP) {
		spme_gpu->vspme_host._init(); // important

		float c[3] = {xl[0], xl[1], xl[2]}, w[3] = {xd[0], xd[1], xd[2]};
		int dim[3] = {dim_spme[0], dim_spme[1], dim_spme[2]};
		if (bAutoExpandZ) spme_gpu->vspme_host.zExpandDim(c[2], w[2], dim[2]);

		size_t s1 = sizeof(cuDoubleComplex), s2 = sizeof(fftw_complex);
		if (s1 == s2) {
#ifdef __CPU_FFT__
			if (dim_spme[0] * dim_spme[1] * dim_spme[2] < 32 * 32 * 32) {
				spme_gpu->vspme_host.bUseHostFFT = 1; // use host multi-thread FFTW
				spme_gpu->init_fft_plans(dim[0], dim[1], dim[2]);
				spme_gpu->vspme_host.ftc_host = (cuDoubleComplex*)(spme_gpu->ftc.m_dev);
				spme_gpu->vspme_host.ftd_host = spme_gpu->ftd.m_dev;

				show_log("use multi-thread CPU-FFTW with GPU-SPME", true);
			}
#endif
		}
		else show_infor("WARNNING: cuDoubleComplex is not same as fftw_complex. can not use CPU-FFTW with GPU SPME!", true);

		// do this before setup_vspme_host(...)
		int ndbf = 4096 * 16;
		spme_gpu->dbuf.use_GPU_MapMem(true);
		spme_gpu->dbuf.set_array(ndbf);
		spme_gpu->vspme_host.dbuff_host = spme_gpu->dbuf.m;
		spme_gpu->vspme_host.dbuff_dev = spme_gpu->dbuf.m_dev;
		spme_gpu->vspme_host.ndbuff = spme_gpu->dbuf.n;

		bool status = _gpu_md_::_spme_::setup_vspme_host(&(spme_gpu->vspme_host), &(spme_gpu->vspme_dev), false, spme_gpu->status.m_dev, 
			bsp->indx, bsp->x, bsp->y[0], bsp->y[1], bsp->y[2], bsp->npts, 
			c, w, kappa, dim[0], dim[1], dim[2], bsp->indx, mMP, spme_gpu->r.m, spme_gpu->q.m, spme_gpu->mu.m_dev, 
			spme_gpu->E.m_dev, spme_gpu->F.m_dev, natoms);

		return status;
	}
}

namespace _gpu_md_ {
} // end of namespace _gpu_md_

namespace _atomic_cmm_ {
	extern void cudaEwaldSum2d_K0(EAtCELL *acell);
	extern void cpuEwaldSum2d_K0(EAtCELL *acell);
	extern void gpu_reset_esr_neighbors(EAtCELL *acell);

	/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
	void SPME_EF(EAtCELL &cell, _gpu_md_::GPU_SPME_2d &vspme, int nThreads, _spme_2d_::SPME_VAR& var, InteractRes& res, bool bInitCoordination) {
		res.reset();
		InteractRes tres;

		int i;
		JobIndx job[MAX_THREADS];
		for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);

		cell.bCalE = (var.spme_3d_var.bEfieldOnly ? 1 : 2);
	
		// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
		int nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.atom.n - 1, 400, nThreads);
		if (cell.bCalE == 1 || cell.bCalE == 11) MTOperate1<JobIndx, _gpu_md_::GPU_SPME >((void*)(&_gpu_md_::reset_global_efield), job, nthreads, (_gpu_md_::GPU_SPME*)(&vspme), true);
		else if (cell.bCalE == 2 || cell.bCalE == 12) MTOperate1<JobIndx, _gpu_md_::GPU_SPME >((void*)(&_gpu_md_::reset_global_force), job, nthreads, (_gpu_md_::GPU_SPME*)(&vspme), true);

		// ***************** following should be done in GPU  ***************
		if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
		else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

		// calculate short-range interaction with GPU, and
		// 2d Ewald Sum has K0 term calculated in real space
		cell.k0var = var.K0var;
		cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.spme_3d_var.esv1));
		AsynRun<EAtCELL> asynRun;
		asynRun.set_par(&cell);
		if (cell.bUseGPU) asynRun.start_func((void*)&cudaEwaldSum2d_K0);
		else asynRun.start_func((void*)&cpuEwaldSum2d_K0);
		//asynRun.wait();
		// ***************** end of jobs done in GPU  ***************

#if _CHECK_TIME_STAMP_ == 1
		char msg[256] = "\0";
		TIME t; int dt;
		t.start();
#endif
		// accumulate the electric field from the image part
		bool bQ = true, bMU = (vspme.vspme_host.mMP > 0 ? true : false);
		if (bInitCoordination) {
			// cubic cell in meshed cell could be not enough, or overflow, status_host[0] is the amount overflowed
			_gpu_md_::_spme_::init_normalized_coordinates(&(vspme.vspme_host), vspme.vspme_dev);
			if (vspme.vspme_host.status_host[0] > 0) {
				show_log("WARNING: SPME meshed cell on GPU has atoms overflowed from the cubic, SPME calculation is not correct!", true);
			}
		}

#if _CHECK_TIME_STAMP_ == 1
		dt = t.elapse(); sprintf(msg, "time stamp 1: %d ms", dt); show_log(msg, true); t.start();
#endif
		//bool bstatus1 = _gpu_dstruct_::_gpu_util_::check_gpu();
		_gpu_md_::_spme_::cal_Q_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU); // SPME Q matrix

#if _CHECK_TIME_STAMP_ == 1
		dt = t.elapse(); sprintf(msg, "time stamp 2 Q: %d ms", dt); show_log(msg, true);  t.start();
#endif

		//bool bstatus2 = _gpu_dstruct_::_gpu_util_::check_gpu();
		if (vspme.vspme_host.bUseHostFFT == 0 || vspme.vspme_host.ftd_host == NULL) {
			_gpu_md_::_spme_::vspme_FT_IFT(&(vspme.vspme_host), vspme.vspme_dev); //Fourier transform and inverse Fourier Transform in SPME 
		}
		else {
			_gpu_md_::_spme_::vspme_FT_prior(&(vspme.vspme_host), vspme.vspme_dev); //Fourier transform and inverse Fourier Transform in SPME 
			fftw_execute(vspme.ft_plan); // ftd ==FFT==> ftc
			_gpu_md_::_spme_::vspme_FT_post(&(vspme.vspme_host), vspme.vspme_dev); //Fourier transform and inverse Fourier Transform in SPME 
			fftw_execute(vspme.ift_plan); // ftc ==FFT==> ftd
			_gpu_md_::_spme_::vspme_IFT_post(&(vspme.vspme_host), vspme.vspme_dev); //Fourier transform and inverse Fourier Transform in SPME 
		}

#if _CHECK_TIME_STAMP_ == 1
		dt = t.elapse(); sprintf(msg, "time stamp 3 FFT: %d ms", dt); show_log(msg, true);  t.start();
#endif
		double U = 0, STxx, STyy, STzz;
		//bool bstatus = _gpu_dstruct_::_gpu_util_::check_gpu();
		if (!var.spme_3d_var.bEfieldOnly) { // force only
			U = _gpu_md_::_spme_::cal_U_recp(&(vspme.vspme_host), vspme.vspme_dev) * ::eUnit_estat_kT; // calculate the energy
			res.U += U;

			// calculate the electric field, or force for each atom in reciprocal space, without dipole
			_gpu_md_::_spme_::cal_EF_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU, var.spme_3d_var.bEfieldOnly, !var.spme_3d_var.bEfieldOnly); 
			_gpu_md_::dump_atomic_F(&vspme, nThreads);

			
			// stress tensor TQ
			// this has to be done after U is calculated, because the calculation of U used ftd, and TQ ==> FTQ overwrite ftd
			_gpu_md_::_spme_::cal_TQ_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU);
			if (vspme.vspme_host.bUseHostFFT == 0 || vspme.vspme_host.ftd_host == NULL) {
				_gpu_md_::_spme_::cal_FTQ(&(vspme.vspme_host), vspme.vspme_dev); // used GPU-FFTW for Fourier transform
			}
			else {
				// stress tensor TQ ==> FTQ
				_gpu_md_::_spme_::cal_FTQ_prior(&(vspme.vspme_host), vspme.vspme_dev); // TQ ==> ftd
				fftw_execute(vspme.ft_plan); // ftd ==FFT==> ftc
				_gpu_md_::_spme_::cal_FTQ_post(&(vspme.vspme_host), vspme.vspme_dev); // ftc ==> FTQ
			}
			
			_gpu_md_::_spme_::cal_ST_recp(&(vspme.vspme_host), vspme.vspme_dev, STxx, STyy, STzz);
			res.STxx += STxx * ::fUnit_estat_mass;
			res.STyy += STyy * ::fUnit_estat_mass;
			res.STzz += STzz * ::fUnit_estat_mass;
			res.STtr += (STxx + STyy + STzz) * ::fUnit_estat_mass;
		}
		else { // electric field only
			// calculate the electric field, or force for each atom in reciprocal space, without dipole
			_gpu_md_::_spme_::cal_EF_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU, var.spme_3d_var.bEfieldOnly, !var.spme_3d_var.bEfieldOnly); 
			_gpu_md_::dump_atomic_E(&vspme, nThreads);
		}

#if _CHECK_TIME_STAMP_ == 1
		dt = t.elapse(); sprintf(msg, "time stamp 4 EF: %d ms", dt); show_log(msg, true); t.start();
#endif

		// virial coming from the surface dipole, U(V), from receiprcal part
		if (!var.spme_3d_var.bEfieldOnly && var.spme_3d_var.bVirial) {
			EwaldSum_SurfaceDipole(cell.mu_surf, (EwaldSumRealVars0*)(&(var.spme_3d_var.esv1)), tres);
			res += tres;
		}

		//sprintf(errmsg, "After surface dipole correction, virial changes to %f", res.STtr); show_log(errmsg, true);

		// we have to wait for the GPU result on real part ewald-sum
		asynRun.wait(); 
		// accumulating the electric field or force in CELL, or GPU
		nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
		if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
		else if (cell.bCalE == 2) { // force
			MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
		}


		// then we needs to dump the efield or force of local interaction
		SVECTOR<double, 4> sr_res[MAX_THREADS];
		nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
		if (var.spme_3d_var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _gpu_md_::GPU_SPME >((void*)(&_gpu_md_::collect_efield_sr_job), job, nthreads, &cell, (_gpu_md_::GPU_SPME*)(&vspme), true);
		else {
			MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _gpu_md_::GPU_SPME >((void*)(&_gpu_md_::collect_eforce_sr_job), job, nthreads, &cell, (_gpu_md_::GPU_SPME*)(&vspme), true);
			MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, SVECTOR<double, 4> >((void*)(&collect_Uvirial), job, nthreads, &cell, sr_res, true);
			for (i = 0; i < nthreads; i++) {
				res.U += sr_res[i].v[0];
				res.STxx += sr_res[i].v[1];
				res.STyy += sr_res[i].v[2];
				res.STzz += sr_res[i].v[3];
				res.STtr += sr_res[i].v[1] + sr_res[i].v[2] + sr_res[i].v[3];
			}
		}
	}


extern void dmax(double &x1, double &x2, double &res);
//void dmax(double &x1, double &x2, double &res) {
//	res = (x1 > x2 ? x1 : x2);
//}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
	bool Polarize_SPME_Interact(MMOL_MD_CELL &mdcell, EAtCELL &cell, _gpu_md_::GPU_SPME_2d &vspme, int nThreads, bool bPolarize, _spme_2d_::SPME_VAR &var, InteractRes &res) {
		// we have to calculate the real part of electro-static interaction, and LJ interaction
		// in which case, the electric field of each atom will be initialized/refreshed at the beginning
		// then we can calculate the image part of SPME Ewald Sum
		// and then to accumulate the force from electric field

		bool bEfieldOnly = var.spme_3d_var.bEfieldOnly;
		res.reset();
		InteractRes tres;

#if _CHECK_TIME_STAMP_ == 1
		char msg[256] = "\0";
		TIME t; int dt;
		t.start();
#endif

		_gpu_md_::_spme_::init_normalized_coordinates(&(vspme.vspme_host), vspme.vspme_dev);
		if (vspme.vspme_host.status_host[0] > 0) {
			show_log("WARNING: SPME meshed cell on GPU has atoms overflowed from the cubic, SPME calculation is not correct!", true);
		}

#if _CHECK_TIME_STAMP_ == 1
		dt = t.elapse(); sprintf(msg, "time stamp -- initialize coordination: %d ms", dt); show_log(msg, true); t.start();
#endif

		double dmu_max = 0, vt = 0, muConvg = ::muConvg; // eA
		int iloop = 0, max_loops = ::nmax_loops_polarize;
		bool status = false;
		var.spme_3d_var.bEfieldOnly = true;
		cell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

		gpu_reset_esr_neighbors(&cell);
		while (bPolarize && iloop < max_loops) {
			cal_dipole(mdcell, nThreads, mdcell.mu_surf); // accumulate dipole of each atom, and the whole mdcell
			memcpy(&(cell.mu_surf), &(mdcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
			_atomic_cmm_::update_atomic_charge_dipole(&cell); // update charge and dipole on cell

			memcpy(vspme.mu.m, cell.mu.m, cell.mu.n * sizeof(double)); // update dipole on vspme, we do not need to update charge
			_gpu_md_::_spme_::update_dipole_device_spme(&(vspme.vspme_host));

			SPME_EF(cell, vspme, nThreads, var, tres, false); // calculate the dipole only

			// polarization
			if (mdcell.mm.n > 0) MOperate<MMOLECULE>((void*)(&polarize_MM), mdcell.mm.m, mdcell.mm.n, nThreads);
			if (mdcell.sm.n > 0) MOperate<SMOLECULE>((void*)(&polarize_SM), mdcell.sm.m, mdcell.sm.n, nThreads);
			if (mdcell.pm.n > 0) MOperate<PMOLECULE>((void*)(&polarize_PM), mdcell.pm.m, mdcell.pm.n, nThreads);
			dmu_max = 0; vt = 0;
			if (mdcell.ind_mu_mm.n > 0) MOperate2<MOL_POL<MMOLECULE>, double>((void*)(&MM_maxdiff_backup_polarization), mdcell.ind_mu_mm.m, mdcell.ind_mu_mm.n, nThreads, vt, false, (void*)(&dmax));
			dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
			if (mdcell.ind_mu_sm.n > 0) MOperate2<MOL_POL<SMOLECULE>, double>((void*)(&SM_maxdiff_backup_polarization), mdcell.ind_mu_sm.m, mdcell.ind_mu_sm.n, nThreads, vt, false, (void*)(&dmax));
			dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
			if (mdcell.ind_mu_pm.n > 0) MOperate2<MOL_POL<PMOLECULE>, double>((void*)(&PM_maxdiff_backup_polarization), mdcell.ind_mu_pm.m, mdcell.ind_mu_pm.n, nThreads, vt, false, (void*)(&dmax));
			dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;

			iloop++;
			if (dmu_max <= muConvg && iloop > 0) {status = true; break;}
		}

		//char msg[256] = "\0";
		//if (status) sprintf(msg, "Convergent polarized dipole is obtainted with %d loops", iloop); 
		//else sprintf(msg, "No convergent polarized dipole is obtained with %d loops", iloop);
		//show_log(msg, true);

		cal_dipole(mdcell, nThreads, mdcell.mu_surf); // accumulate dipole of each atom, and the whole mdcell
		memcpy(&(cell.mu_surf), &(mdcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
		_atomic_cmm_::update_atomic_charge_dipole(&cell);

		memcpy(vspme.mu.m, cell.mu.m, cell.mu.n * sizeof(double)); // update dipole on vspme, we do not need to update charge
		_gpu_md_::_spme_::update_dipole_device_spme(&(vspme.vspme_host));

		if (bEfieldOnly) { // calculate f
			SPME_EF(cell, vspme, nThreads, var, res, false); // calculate the dipole and force
		}
		else { // calculate f
			var.spme_3d_var.bEfieldOnly = bEfieldOnly;
			//cell.bCalE = (bEfieldOnly ? 1 : 2);

			SPME_EF(cell, vspme, nThreads, var, res, false); // calculate the dipole and force
		}

		return status;
	}


// polarized SPME, with the interaction between induced dipole excluded
// in some polarized model, the interaction between induced dipole is excluded, e.g. only q-induced dipole is considered.
	void VSPME_F_exclude_dipole(_gpu_md_::GPU_SPME_2d &vspme, int nThreads, _spme_2d_::SPME_VAR &var, InteractRes &res, bool bInitCoordination) {
		// we have to calculate the real part of electro-static interaction, and LJ interaction
		// in which case, the electric field of each atom will be initialized/refreshed at the beginning
		// then we can calculate the image part of SPME Ewald Sum
		// and then to accumulate the force from electric field
		res.reset();
		InteractRes tres;

		// accumulate the electric field from the image part
		if (bInitCoordination) {
			// cubic cell in meshed cell could be not enough, or overflow, status_host[0] is the amount overflowed
			_gpu_md_::_spme_::init_normalized_coordinates(&(vspme.vspme_host), vspme.vspme_dev);
			if (vspme.vspme_host.status_host[0] > 0) {
				show_log("WARNING: SPME meshed cell on GPU has atoms overflowed from the cubic, SPME calculation is not correct!", true);
			}
		}

		bool bQ = false, bMU = true; // do not calculate Q
		//bool bstatus1 = _gpu_dstruct_::_gpu_util_::check_gpu();
		_gpu_md_::_spme_::cal_Q_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU); // SPME Q matrix
		if (!var.spme_3d_var.bEfieldOnly) _gpu_md_::_spme_::cal_TQ_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU); // TQ

		//bool bstatus2 = _gpu_dstruct_::_gpu_util_::check_gpu();
		if (vspme.vspme_host.bUseHostFFT == 0 || vspme.vspme_host.ftd_host == NULL) {
			_gpu_md_::_spme_::vspme_FT_IFT(&(vspme.vspme_host), vspme.vspme_dev); //Fourier transform and inverse Fourier Transform in SPME 
		}
		else {
			_gpu_md_::_spme_::vspme_FT_prior(&(vspme.vspme_host), vspme.vspme_dev); //Fourier transform and inverse Fourier Transform in SPME 
			fftw_execute(vspme.ft_plan); // ftd ==FFT==> ftc
			_gpu_md_::_spme_::vspme_FT_post(&(vspme.vspme_host), vspme.vspme_dev); //Fourier transform and inverse Fourier Transform in SPME 
			fftw_execute(vspme.ift_plan); // ftc ==FFT==> ftd
			_gpu_md_::_spme_::vspme_IFT_post(&(vspme.vspme_host), vspme.vspme_dev); //Fourier transform and inverse Fourier Transform in SPME 
		}

		
		double U = 0, STxx, STyy, STzz;
		//bool bstatus = _gpu_dstruct_::_gpu_util_::check_gpu();
		if (!var.spme_3d_var.bEfieldOnly) { // force only
			U = _gpu_md_::_spme_::cal_U_recp(&(vspme.vspme_host), vspme.vspme_dev) * ::eUnit_estat_kT; // calculate the energy
			res.U += U;

			// calculate the electric field, or force for each atom in reciprocal space, without dipole
			_gpu_md_::_spme_::cal_EF_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU, var.spme_3d_var.bEfieldOnly, !var.spme_3d_var.bEfieldOnly); 
			_gpu_md_::dump_subtract_atomic_F(&vspme, nThreads);

			
			// stress tensor TQ
			// this has to be done after U is calculated, because the calculation of U used ftd, and TQ ==> FTQ overwrite ftd
			_gpu_md_::_spme_::cal_TQ_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU);
			if (vspme.vspme_host.bUseHostFFT == 0 || vspme.vspme_host.ftd_host == NULL) {
				_gpu_md_::_spme_::cal_FTQ(&(vspme.vspme_host), vspme.vspme_dev); // used GPU-FFTW for Fourier transform
			}
			else {
				// stress tensor TQ ==> FTQ
				_gpu_md_::_spme_::cal_FTQ_prior(&(vspme.vspme_host), vspme.vspme_dev); // TQ ==> ftd
				fftw_execute(vspme.ft_plan); // ftd ==FFT==> ftc
				_gpu_md_::_spme_::cal_FTQ_post(&(vspme.vspme_host), vspme.vspme_dev); // ftc ==> FTQ
			}
			
			_gpu_md_::_spme_::cal_ST_recp(&(vspme.vspme_host), vspme.vspme_dev, STxx, STyy, STzz);
			res.STxx += STxx * ::fUnit_estat_mass;
			res.STyy += STyy * ::fUnit_estat_mass;
			res.STzz += STzz * ::fUnit_estat_mass;
			res.STtr += (STxx + STyy + STzz) * ::fUnit_estat_mass;
		}
		else { // electric field only
			// calculate the electric field, or force for each atom in reciprocal space, without dipole
			_gpu_md_::_spme_::cal_EF_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU, var.spme_3d_var.bEfieldOnly, !var.spme_3d_var.bEfieldOnly); 
			_gpu_md_::dump_subtract_atomic_E(&vspme, nThreads);
		}
	}


	bool Polarize_SPME_Interact2(MMOL_MD_CELL &mdcell, EAtCELL &cell, _gpu_md_::GPU_SPME_2d &vspme, int nThreads, bool bPolarize, _spme_2d_::SPME_VAR &var, InteractRes &res) {
		bool status = Polarize_SPME_Interact(mdcell, cell, vspme, nThreads, bPolarize, var, res);
		InteractRes tres;
		if (bPolarize && var.spme_3d_var.esv1.iExIndDipole == 1) {
			VSPME_F_exclude_dipole(vspme, nThreads, var, tres, false);
			res -= tres;
		}
		return status;
	}

} // end of namespace _atomic_cmm_

#endif
