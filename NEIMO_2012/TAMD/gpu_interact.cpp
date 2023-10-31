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
#include "gpu_interact.h"

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);

double max_diff(double *a1, double *a2, int n) {
	double vd, vdmax = 0;
	for (int i = 0; i < n; i++) {
		vd = a2[i] - a1[i]; vd = FABS(vd);
		vdmax = (vd > vdmax ? vd : vdmax);
	}
	return vdmax;
};

#ifdef __GPU__

namespace _gpu_md_ {
	bool setup_gpu_spme(_gpu_md_::GPU_SPME_3d *spme_gpu, BSplineFunc<2> *bsp, float* xl, float *xd, int *dim_spme, double kappa, int natoms, int mMP) {
		spme_gpu->vspme_host._init(); // important

		size_t s1 = sizeof(cuDoubleComplex), s2 = sizeof(fftw_complex);
		if (s1 == s2) {
#ifdef __CPU_FFT__
			if (dim_spme[0] * dim_spme[1] * dim_spme[2] < 32 * 32 * 32) { //< 64 * 64 * 128) {
				spme_gpu->vspme_host.bUseHostFFT = 1; // use host multi-thread FFTW
				spme_gpu->init_fft_plans(dim_spme[0], dim_spme[1], dim_spme[2]);
				spme_gpu->vspme_host.ftc_host = (cuDoubleComplex*)(spme_gpu->ftc.m_dev);
				spme_gpu->vspme_host.ftd_host = spme_gpu->ftd.m_dev;

				show_log("use multi-thread CPU-FFTW with GPU-SPME", true);
			}
			else show_log("use cuFFT for GPU-SPME", true);
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
		bsp->indx, bsp->x, bsp->y[0], bsp->y[1], bsp->y[2], bsp->npts, xl, xd, kappa, 
		dim_spme[0], dim_spme[1], dim_spme[2], bsp->indx, mMP, spme_gpu->r.m, spme_gpu->q.m, spme_gpu->mu.m_dev, 
		spme_gpu->E.m_dev, spme_gpu->F.m_dev, natoms);

		if (!status) show_infor("failure on initialize memory", true);

		return status;
	}
}

namespace _gpu_md_ {
	void init_atom(ARRAY<MATOM*> &matom, ARRAY<AtomGPUSPME> &atom) {
		if (matom.n != atom.n) return;
		int i;
		for (i = 0; i < matom.n; i++) {
			atom.m[i].r = matom.m[i]->rg.v;
			atom.m[i].q = &(matom.m[i]->c);
			atom.m[i].mu = matom.m[i]->mu.v;
			atom.m[i].E = matom.m[i]->E.v;
			atom.m[i].F = matom.m[i]->F.v;
		}
	};

	void update_atom_coordinate(JobIndx *job, GPU_SPME *gspme) {
		int i, i3;
		for (i = job->n1; i <= job->n2; i++) {
			if (i >= gspme->atom.n) break;
			i3 = i * 3;
			memcpy(gspme->r.m + i3, gspme->atom.m[i].r, SIZE_V3);
		}
	};
	void update_atom_dipole(JobIndx *job, GPU_SPME *gspme) {
		int i, i3;
		for (i = job->n1; i <= job->n2; i++) {
			if (i >= gspme->atom.n) break;
			i3 = i * 3;
			memcpy(gspme->mu.m + i3, gspme->atom.m[i].mu, SIZE_V3);
		}
	};

	void update_atomic_coordinate(GPU_SPME *gspme, int nThreads) {
		int i;
		JobIndx job[MAX_THREADS];
		for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
		int nthreads = -1; assignJobIndx1(job, nthreads, 0, gspme->atom.n - 1, 200, nThreads);
		MTOperate1<JobIndx, GPU_SPME>((void*)(&update_atom_coordinate), job, nthreads, gspme, true);
	};

	void update_atomic_dipole(GPU_SPME *gspme, int nThreads) {
		int i;
		JobIndx job[MAX_THREADS];
		for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
		int nthreads = -1; assignJobIndx1(job, nthreads, 0, gspme->atom.n - 1, 200, nThreads);
		MTOperate1<JobIndx, GPU_SPME>((void*)(&update_atom_dipole), job, nthreads, gspme, true);
	};

	void dump_atom_E(JobIndx *job, GPU_SPME *gspme) {
		int i, i3;
		for (i = job->n1; i <= job->n2; i++) {
			if (i >= gspme->atom.n) break;
			i3 = i * 3;
			gspme->atom.m[i].E[0] += gspme->E.m[i3];
			gspme->atom.m[i].E[1] += gspme->E.m[i3 + 1];
			gspme->atom.m[i].E[2] += gspme->E.m[i3 + 2];
		}
	};
	void dump_atom_F(JobIndx *job, GPU_SPME *gspme) {
		int i, i3;
		for (i = job->n1; i <= job->n2; i++) {
			if (i >= gspme->atom.n) break;
			i3 = i * 3;
			gspme->atom.m[i].F[0] += gspme->F.m[i3] * ::fUnit_estat_mass;
			gspme->atom.m[i].F[1] += gspme->F.m[i3 + 1] * ::fUnit_estat_mass;
			gspme->atom.m[i].F[2] += gspme->F.m[i3 + 2] * ::fUnit_estat_mass;
		}
	};

	void dump_atomic_E(GPU_SPME *gspme, int nThreads) {
		int i;
		JobIndx job[MAX_THREADS];
		for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
		int nthreads = -1; assignJobIndx1(job, nthreads, 0, gspme->atom.n - 1, 200, nThreads);
		MTOperate1<JobIndx, GPU_SPME>((void*)(&dump_atom_E), job, nthreads, gspme, true);
	};
	void dump_atomic_F(GPU_SPME *gspme, int nThreads) {
		int i;
		JobIndx job[MAX_THREADS];
		for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
		int nthreads = -1; assignJobIndx1(job, nthreads, 0, gspme->atom.n - 1, 200, nThreads);
		MTOperate1<JobIndx, GPU_SPME>((void*)(&dump_atom_F), job, nthreads, gspme, true);
	};


	void dump_subtract_atom_E(JobIndx *job, GPU_SPME *gspme) {
		int i, i3;
		for (i = job->n1; i <= job->n2; i++) {
			if (i >= gspme->atom.n) break;
			i3 = i * 3;
			gspme->atom.m[i].E[0] -= gspme->E.m[i3];
			gspme->atom.m[i].E[1] -= gspme->E.m[i3 + 1];
			gspme->atom.m[i].E[2] -= gspme->E.m[i3 + 2];
		}
	};
	void dump_subtract_atom_F(JobIndx *job, GPU_SPME *gspme) {
		int i, i3;
		for (i = job->n1; i <= job->n2; i++) {
			if (i >= gspme->atom.n) break;
			i3 = i * 3;
			gspme->atom.m[i].F[0] -= gspme->F.m[i3] * ::fUnit_estat_mass;
			gspme->atom.m[i].F[1] -= gspme->F.m[i3 + 1] * ::fUnit_estat_mass;
			gspme->atom.m[i].F[2] -= gspme->F.m[i3 + 2] * ::fUnit_estat_mass;
		}
	};
	void dump_subtract_atomic_E(GPU_SPME *gspme, int nThreads) {
		int i;
		JobIndx job[MAX_THREADS];
		for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
		int nthreads = -1; assignJobIndx1(job, nthreads, 0, gspme->atom.n - 1, 200, nThreads);
		MTOperate1<JobIndx, GPU_SPME>((void*)(&dump_subtract_atom_E), job, nthreads, gspme, true);
	};
	void dump_subtract_atomic_F(GPU_SPME *gspme, int nThreads) {
		int i;
		JobIndx job[MAX_THREADS];
		for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
		int nthreads = -1; assignJobIndx1(job, nthreads, 0, gspme->atom.n - 1, 200, nThreads);
		MTOperate1<JobIndx, GPU_SPME>((void*)(&dump_subtract_atom_F), job, nthreads, gspme, true);
	};

	void collect_efield_sr_job(JobIndx* job, _atomic_cmm_::EAtCELL *cell, GPU_SPME* vspme) {
		double *E;
		int ia, i3a;
		for (ia = job->n1; ia <= job->n2; ia++) {
			i3a = ia * 3;
			if (ia >= vspme->atom.n) break;
			E = vspme->atom.m[ia].E;
			E[0] += cell->E.m[i3a];
			E[1] += cell->E.m[i3a + 1];
			E[2] += cell->E.m[i3a + 2];
		}
	}

	void collect_eforce_sr_job(JobIndx* job, _atomic_cmm_::EAtCELL *cell, GPU_SPME* vspme) {
		double *F;
		int ia, i3a;
		for (ia = job->n1; ia <= job->n2; ia++) {
			i3a = ia * 3;
			if (ia >= vspme->atom.n) break;
			F = vspme->atom.m[ia].F;
			F[0] += cell->Fe.m[i3a] * ::fUnit_estat_mass;
			F[1] += cell->Fe.m[i3a + 1] * ::fUnit_estat_mass;
			F[2] += cell->Fe.m[i3a + 2] * ::fUnit_estat_mass;
		}
	}

	void reset_global_efield(JobIndx *job, GPU_SPME* vspme) {
		int i;
		double *E;
		for (i = job->n1; i <= job->n2; i++) {
			if (i >= vspme->atom.n) break;
			E = vspme->atom.m[i].E;
			memset(E, 0, SIZE_V3);
		}
	};

	void reset_global_force(JobIndx *job, GPU_SPME* vspme) {
		int i;
		double *F;
		for (i = job->n1; i <= job->n2; i++) {
			if (i >= vspme->atom.n) break;
			F = vspme->atom.m[i].F;
			memset(F, 0, SIZE_V3);
		}
	};

	extern __host__ void reset_esr_neighbors(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev, cudaStream_t &stream);

} // _gpu_md_


namespace _atomic_cmm_ {
	extern void cudaEwaldSumReal_3d(EAtCELL *acell);
	extern void cpuEwaldSumReal_3d(EAtCELL *acell);

	void gpu_reset_esr_neighbors(EAtCELL *acell) {
		cudaStream_t cudaStream;
		cudaStreamCreate(&cudaStream);
		_gpu_md_::reset_esr_neighbors(&(acell->gpu_cell), acell->gpu_cell_dev, cudaStream);
		cudaStreamSynchronize(cudaStream);
		cudaStreamDestroy(cudaStream);
	}

	/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
	void SPME_EF(EAtCELL &cell, _gpu_md_::GPU_SPME_3d &vspme, int nThreads, SPME_VAR& var, InteractRes& res, bool bInitCoordination) {
		res.reset();
		InteractRes tres;

		char msg[256] = "\0";

		int i;
		JobIndx job[MAX_THREADS];
		for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);

		cell.bCalE = (var.bEfieldOnly ? 1 : 2);
	
		// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
		int nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.atom.n - 1, 400, nThreads);
		if (cell.bCalE == 1 || cell.bCalE == 11) MTOperate1<JobIndx, _gpu_md_::GPU_SPME >((void*)(&_gpu_md_::reset_global_efield), job, nthreads, (_gpu_md_::GPU_SPME*)(&vspme), true);
		else if (cell.bCalE == 2 || cell.bCalE == 12) MTOperate1<JobIndx, _gpu_md_::GPU_SPME >((void*)(&_gpu_md_::reset_global_force), job, nthreads, (_gpu_md_::GPU_SPME*)(&vspme), true);

		// ***************** following should be done in GPU  ***************
		if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
		else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

		// calculate short-range interaction with GPU
		cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.esv1));
		AsynRun<EAtCELL> asynRun;
		asynRun.set_par(&cell);
		if (cell.bCalUseGPU) asynRun.start_func((void*)&cudaEwaldSumReal_3d);
		else asynRun.start_func((void*)&cpuEwaldSumReal_3d);

		//asynRun.wait();
		// ***************** end of jobs done in GPU  ***************

#if _CHECK_TIME_STAMP_ == 1
		TIME t; int dt;
		t.start();
#endif
		cufftResult cuffstatus;
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
		if (!var.bEfieldOnly) { // force only
			U = _gpu_md_::_spme_::cal_U_recp(&(vspme.vspme_host), vspme.vspme_dev) * ::eUnit_estat_kT; // calculate the energy
			res.U += U;
			
			_gpu_md_::_spme_::cal_ST_recp(&(vspme.vspme_host), vspme.vspme_dev, STxx, STyy, STzz);
			res.STxx += STxx * ::fUnit_estat_mass;
			res.STyy += STyy * ::fUnit_estat_mass;
			res.STzz += STzz * ::fUnit_estat_mass;
			res.STtr += (STxx + STyy + STzz) * ::fUnit_estat_mass;

			// calculate the electric field, or force for each atom in reciprocal space, without dipole
			_gpu_md_::_spme_::cal_EF_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU, var.bEfieldOnly, !var.bEfieldOnly); 
			_gpu_md_::dump_atomic_F(&vspme, nThreads);
		}
		else { // electric field only
			// calculate the electric field, or force for each atom in reciprocal space, without dipole
			_gpu_md_::_spme_::cal_EF_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU, var.bEfieldOnly, !var.bEfieldOnly);
			_gpu_md_::dump_atomic_E(&vspme, nThreads);
		}


#if _CHECK_TIME_STAMP_ == 1
		dt = t.elapse(); sprintf(msg, "time stamp 4 EF: %d ms", dt); show_log(msg, true); t.start();
#endif

		// virial coming from the surface dipole, U(V), from receiprcal part
		if (!var.bEfieldOnly && var.bVirial) {
			EwaldSum_SurfaceDipole(cell.mu_surf, (EwaldSumRealVars0*)(&(var.esv1)), tres);
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

#if _CHECK_TIME_STAMP_ == 1
		dt = t.elapse(); sprintf(msg, "time stamp 5 wating to finish real-EwaldSum: %d ms", dt); show_log(msg, true); t.start();
#endif

		// then we needs to dump the efield or force of local interaction
		SVECTOR<double, 4> sr_res[MAX_THREADS];
		nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
		if (var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _gpu_md_::GPU_SPME >((void*)(&_gpu_md_::collect_efield_sr_job), job, nthreads, &cell, (_gpu_md_::GPU_SPME*)(&vspme), true);
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
	bool Polarize_SPME_Interact(MMOL_MD_CELL &mdcell, EAtCELL &cell, _gpu_md_::GPU_SPME_3d &vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
		// we have to calculate the real part of electro-static interaction, and LJ interaction
		// in which case, the electric field of each atom will be initialized/refreshed at the beginning
		// then we can calculate the image part of SPME Ewald Sum
		// and then to accumulate the force from electric field

		bool bEfieldOnly = var.bEfieldOnly;
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
		var.bEfieldOnly = true;
		cell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

		//gpu_reset_esr_neighbors(&cell);
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

#if _CHECK_TIME_STAMP_ == 1
		dt = t.elapse(); sprintf(msg, "time stamp -- consistent polarization takes: %d ms", dt); show_log(msg, true); t.start();
#endif

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
			var.bEfieldOnly = bEfieldOnly;
			//cell.bCalE = (bEfieldOnly ? 1 : 2);

			SPME_EF(cell, vspme, nThreads, var, res, false); // calculate the dipole and force
		}

#if _CHECK_TIME_STAMP_ == 1
		dt = t.elapse(); sprintf(msg, "time stamp -- calculating eforce takes: %d ms", dt); show_log(msg, true); t.start();
#endif

		return status;
	}


// polarized SPME, with the interaction between induced dipole excluded
// in some polarized model, the interaction between induced dipole is excluded, e.g. only q-induced dipole is considered.
	void VSPME_F_exclude_dipole(_gpu_md_::GPU_SPME_3d &vspme, int nThreads, SPME_VAR &var, InteractRes &res, bool bInitCoordination) {
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
		if (!var.bEfieldOnly) { // force only
			U = _gpu_md_::_spme_::cal_U_recp(&(vspme.vspme_host), vspme.vspme_dev) * ::eUnit_estat_kT; // calculate the energy
			res.U += U;
			_gpu_md_::_spme_::cal_ST_recp(&(vspme.vspme_host), vspme.vspme_dev, STxx, STyy, STzz);
			res.STxx += STxx * ::fUnit_estat_mass;
			res.STyy += STyy * ::fUnit_estat_mass;
			res.STzz += STzz * ::fUnit_estat_mass;
			res.STtr += (STxx + STyy + STzz) * ::fUnit_estat_mass;

			// calculate the electric field, or force for each atom in reciprocal space, without dipole
			_gpu_md_::_spme_::cal_EF_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU, var.bEfieldOnly, !var.bEfieldOnly); 
			_gpu_md_::dump_subtract_atomic_F(&vspme, nThreads);
		}
		else { // electric field only
			// calculate the electric field, or force for each atom in reciprocal space, without dipole
			_gpu_md_::_spme_::cal_EF_spme(&(vspme.vspme_host), vspme.vspme_dev, bQ, bMU, var.bEfieldOnly, !var.bEfieldOnly); 
			_gpu_md_::dump_subtract_atomic_E(&vspme, nThreads);
		}
	}


	bool Polarize_SPME_Interact2(MMOL_MD_CELL &mdcell, EAtCELL &cell, _gpu_md_::GPU_SPME_3d &vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
		bool status = Polarize_SPME_Interact(mdcell, cell, vspme, nThreads, bPolarize, var, res);
		InteractRes tres;
		if (bPolarize && var.esv1.iExIndDipole == 1) {
			VSPME_F_exclude_dipole(vspme, nThreads, var, tres, false);
			res -= tres;
		}
		return status;
	}

} // end of namespace _atomic_cmm_

#endif //__GPU__
