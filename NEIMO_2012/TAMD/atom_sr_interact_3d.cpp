#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "ranlib.h"

#include "def.h"
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

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
//long t_real = 0, t_recp = 0, t_emp = 0, t_cluster = 0, t_LJ = 0;
#endif


extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
extern void cal_dipole(MMOL_VMD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);

extern PAIRWISE_DB< PAIRWISE<EF_DPROF> > sr_db; // pair-wise short-range interaction database

#ifdef __GPU__
namespace _gpu_md_ {
	extern __host__ void GPU_EwaldSum_Real_3d(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev, cudaStream_t &stream);
}
#endif

namespace _atomic_cmm_ {
#ifdef __GPU__
extern void gpu_reset_esr_neighbors(EAtCELL *acell);
#endif

#ifdef __GPU__
void cudaEwaldSumReal_3d(EAtCELL *acell) { // asynchrotron thread-function
	if (cudaSetDevice(acell->igpuDev) != cudaSuccess) {
		show_infor("GPU is not enabled!", true); return;
	}
	acell->gpu_prior_cal(3);
	cudaStream_t cudaStream;
	cudaStreamCreate(&cudaStream);
	_gpu_md_::GPU_EwaldSum_Real_3d(&(acell->gpu_cell), acell->gpu_cell_dev, cudaStream);
	cudaStreamSynchronize(cudaStream);
	cudaStreamDestroy(cudaStream);
}

void cpuEwaldSumReal_3d(EAtCELL *cell) {
	int nThreads = MAX_THREADS, nthreads = nThreads;

	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);

	// ***************** following should be done in GPU  ***************
	//if (cell->bCalE == 1 || cell->bCalE == 11) cell->reset_efield(); // efield calculation
	//else if (cell->bCalE == 2 || cell->bCalE == 12) cell->reset_eforce(); // force and energy calculation

	assignJobIndx(job, nThreads, 0, cell->Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell->Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, cell, true);
	
	if (cell->Nb > 0) { // bound
		if (cell->Nb > MAX_THREADS * 10) {
			assignJobIndx(job, nThreads, 0, cell->Nb - 1);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nThreads, cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell->Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, cell);}
	}
	// need to do 1-4 neighbors here

	/*
	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell->Na - 1, 400, nThreads);
	if (cell->bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, cell, true);
	else if (cell->bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, cell, true);
	}
	*/
	// ***************** end of jobs done in GPU  ***************
}

#endif

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void LocalF(EAtCELL &cell, _spme_::VSPME<_evatom_::_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	res.reset();
	InteractRes tres;

	cell.bCalE = 2; // force only

	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	
	// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
	int nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
	if (cell.bCalE == 1 || cell.bCalE == 11) MTOperate1<JobIndx, _spme_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::reset_global_efield_q_3d), job, nthreads, &vspme, true);
	else if (cell.bCalE == 2 || cell.bCalE == 12) MTOperate1<JobIndx, _spme_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::reset_global_force_q_3d), job, nthreads, &vspme, true);

	// ***************** following should be done in GPU  ***************

	if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
	else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

#ifdef __GPU__
	// calculate short-range interaction with GPU
	cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.esv1));
	AsynRun<EAtCELL> asynRun;
	asynRun.set_par(&cell);
	if (cell.bCalUseGPU) asynRun.start_func((void*)&cudaEwaldSumReal_3d);
	else asynRun.start_func((void*)&cpuEwaldSumReal_3d);
	//asynRun.wait();
#else
	assignJobIndx(job, nThreads, 0, cell.Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, &cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell.Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, &cell, true);
	
	if (cell.Nb > 0) { // bound
		if (cell.Nb > MAX_THREADS * 10) {
			assignJobIndx(job, nThreads, 0, cell.Nb - 1);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nThreads, &cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell.Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, &cell);}
	}
	// need to do 1-4 neighbors here

	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif
	// ***************** end of jobs done in GPU  ***************

#ifdef __GPU__ // we have to wait for the GPU result on real part ewald-sum
	asynRun.wait();
	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif


	// then we needs to dump the efield or force of local interaction
	SVECTOR<double, 4> sr_res[MAX_THREADS];
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::collect_efield_sr_q_3d), job, nthreads, &cell, &vspme, true);
	else {
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::collect_fresult_sr_q_3d), job, nthreads, &cell, &vspme, true);
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

void LocalF(EAtCELL &cell, _spme_::VSPME<_evatom_::_VMUatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	res.reset();
	InteractRes tres;

	cell.bCalE = 2; // force only

	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);

	// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
	int nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
	if (cell.bCalE == 1 || cell.bCalE == 11) MTOperate1<JobIndx, _spme_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::reset_global_efield_mu_3d), job, nthreads, &vspme, true);
	else if (cell.bCalE == 2 || cell.bCalE == 12) MTOperate1<JobIndx, _spme_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::reset_global_force_mu_3d), job, nthreads, &vspme, true);

	// ***************** following should be done in GPU  ***************

	if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
	else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

#ifdef __GPU__
	// calculate short-range interaction with GPU
	cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.esv1));
	AsynRun<EAtCELL> asynRun;
	asynRun.set_par(&cell);
	if (cell.bCalUseGPU) asynRun.start_func((void*)&cudaEwaldSumReal_3d);
	else asynRun.start_func((void*)&cpuEwaldSumReal_3d);
	//asynRun.wait();
#else
	assignJobIndx(job, nThreads, 0, cell.Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, &cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell.Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, &cell, true);
	
	if (cell.Nb > 0) { // bound
		if (cell.Nb > MAX_THREADS * 10) {
			assignJobIndx(job, nThreads, 0, cell.Nb - 1);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nThreads, &cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell.Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, &cell);}
	}
	// need to do 1-4 neighbors here

	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif
	// ***************** end of jobs done in GPU  ***************

#ifdef __GPU__ // we have to wait for the GPU result on real part ewald-sum
	asynRun.wait();
	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif


	// then we needs to dump the efield or force of local interaction
	SVECTOR<double, 4> sr_res[MAX_THREADS];
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, nThreads);
	if (var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::collect_efield_sr_mu_3d), job, nthreads, &cell, &vspme, true);
	else {
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::collect_fresult_sr_mu_3d), job, nthreads, &cell, &vspme, true);
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


/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(/*MMOL_MD_CELL &mdcell,*/ EAtCELL &cell, _spme_::VSPME<_evatom_::_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	res.reset();
	InteractRes tres;

	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);

	cell.bCalE = (var.bEfieldOnly ? 1 : 2);
	
	// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
	int nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
	if (cell.bCalE == 1 || cell.bCalE == 11) MTOperate1<JobIndx, _spme_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::reset_global_efield_q_3d), job, nthreads, &vspme, true);
	else if (cell.bCalE == 2 || cell.bCalE == 12) MTOperate1<JobIndx, _spme_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::reset_global_force_q_3d), job, nthreads, &vspme, true);

	// ***************** following should be done in GPU  ***************
	if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
	else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

#ifdef __GPU__
	// calculate short-range interaction with GPU
	cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.esv1));
	AsynRun<EAtCELL> asynRun;
	asynRun.set_par(&cell);
	if (cell.bCalUseGPU) asynRun.start_func((void*)&cudaEwaldSumReal_3d);
	else asynRun.start_func((void*)&cpuEwaldSumReal_3d);
	//asynRun.wait();
#else
	assignJobIndx(job, nThreads, 0, cell.Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, &cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell.Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, &cell, true);
	
	if (cell.Nb > 0) { // bound
		if (cell.Nb > MAX_THREADS * 10) {
			assignJobIndx(job, nThreads, 0, cell.Nb - 1);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nThreads, &cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell.Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, &cell);}
	}
	// need to do 1-4 neighbors here

	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif
	// ***************** end of jobs done in GPU  ***************


	// accumulate the electric field from the image part
	//vspme.reset_local_EF(); // reset the variable for E & F
	init_normalized_coordinates(vspme, nThreads);
	cal_Q_spme(vspme, nThreads); // SPME Q matrix
	double U = vspme.cal_U_recp(); //calculate the energy from the image part 
	res.U += U;
	cal_E_recp(vspme, !var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.bEfieldOnly && var.bVirial) {
		vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {res.STxx += vspme.STxx; res.STyy += vspme.STyy; res.STzz += vspme.STzz;}
		res.STtr += vspme.STtr;
	}
	//sprintf(errmsg, "Recip : %f, %f", U, vspme.STtr); show_log(errmsg, true);

	// virial coming from the surface dipole, U(V), from receiprcal part
	if (!var.bEfieldOnly && var.bVirial) {
		EwaldSum_SurfaceDipole(cell.mu_surf, (EwaldSumRealVars0*)(&(var.esv1)), tres);
		res += tres;
	}
	
	//sprintf(errmsg, "After surface dipole correction, virial changes to %f", res.STtr); show_log(errmsg, true);

	bool bDumpF = (cell.bCalE == 2 ? true : false);
	bool bDumpE = (cell.bCalE == 1 ? true : false);
	if (nThreads < 2 || vspme.va.n < N4PARALLEL) vspme.dump_EF(bDumpF, true, bDumpE); // transfer the calculated E & F to the variables in real atom
	else {
		assignJobIndx(job, nThreads, 0, vspme.va.n - 1);
		MTOperate3< JobIndx, _spme_::VSPME<_evatom_::_VQatom>, bool, bool>((void*)(&_spme_::dump_EF_Q_Job), job, nThreads, &vspme, &bDumpE, &bDumpF, true);
	}


#ifdef __GPU__ // we have to wait for the GPU result on real part ewald-sum
	asynRun.wait();
	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif


	// then we needs to dump the efield or force of local interaction
	SVECTOR<double, 4> sr_res[MAX_THREADS];
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::collect_efield_sr_q_3d), job, nthreads, &cell, &vspme, true);
	else {
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::collect_fresult_sr_q_3d), job, nthreads, &cell, &vspme, true);
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

// with dipole
/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(/*MMOL_MD_CELL &mdcell, */EAtCELL &cell, _spme_::VSPME<_evatom_::_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res, bool bInitCoordination) {
	res.reset();
	InteractRes tres;

	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);

	cell.bCalE = (var.bEfieldOnly ? 1 : 2);
	
	// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
	int nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _spme_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::reset_global_efield_mu_3d), job, nthreads, &vspme, true);
	else if (cell.bCalE == 2) MTOperate1<JobIndx, _spme_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::reset_global_force_mu_3d), job, nthreads, &vspme, true);

	// ***************** following should be done in GPU  ***************

	if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
	else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

#ifdef __GPU__
	// calculate short-range interaction with GPU
	cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.esv1));
	AsynRun<EAtCELL> asynRun;
	asynRun.set_par(&cell);
	if (cell.bCalUseGPU) asynRun.start_func((void*)&cudaEwaldSumReal_3d);
	else asynRun.start_func((void*)&cpuEwaldSumReal_3d);
#else
	assignJobIndx(job, nThreads, 0, cell.Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, &cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell.Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, &cell, true);
	
	if (cell.Nb > 0) { // bound
		if (cell.Nb > MAX_THREADS * 10) {
			assignJobIndx(job, nThreads, 0, cell.Nb - 1);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nThreads, &cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell.Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, &cell);}
	}
	// need to do 1-4 neighbors here

	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif
	// ***************** end of jobs done in GPU  ***************



	// accumulate the electric field from the image part
	//vspme.reset_local_EF(); // reset the variable for E & F
	if (bInitCoordination) init_normalized_coordinates(vspme, nThreads); //-- vspme is initilized outside before this funciton is called
	cal_Q_spme(vspme, nThreads); // SPME Q matrix
	double U = vspme.cal_U_recp(); //calculate the energy from the image part 
	res.U += U;
	cal_E_recp(vspme, !var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.bEfieldOnly && var.bVirial) {
		cal_StressTensor(vspme, nThreads);
		if (var.bST_diagonal) {
			res.STxx += vspme.STxx; res.STyy += vspme.STyy; res.STzz += vspme.STzz;
		}
		res.STtr += vspme.STtr;
	}

	// virial coming from the surface dipole, U(V), from receiprical part
	if (!var.bEfieldOnly && var.bVirial) {
		EwaldSum_SurfaceDipole(cell.mu_surf, (EwaldSumRealVars0*)(&(var.esv1)), tres);
		res += tres;
	}


	bool bDumpF = (cell.bCalE == 2 ? true : false);
	bool bDumpE = (cell.bCalE == 1 ? true : false);
	if (nThreads < 2 || vspme.va.n < N4PARALLEL) vspme.dump_EF(bDumpF, true, bDumpE); // transfer the calculated E & F to the variables in real atom
	else {
		assignJobIndx(job, nThreads, 0, vspme.va.n - 1);
		MTOperate3< JobIndx, _spme_::VSPME<_evatom_::_VMUatom>, bool, bool>((void*)(&_spme_::dump_EF_MU_Job), job, nThreads, &vspme, &bDumpE, &bDumpF, true);
	}



#ifdef __GPU__ // we have to wait for the GPU result on real part ewald-sum
	char msg[256] = "\0";
	asynRun.wait();
	/*
	sprintf(msg, "==>  total neighbors %d", cell.gpu_cell.ivar_map[1]); show_log(msg, true);
	sprintf(msg, "==>  neighbor overflowed %d", cell.gpu_cell.ivar_map[0]); show_log(msg, true);
	sprintf(msg, "==>  ignored atoms %d", cell.gpu_cell.ivar_map[2]); show_log(msg, true);
	sprintf(msg, "==>  center atom on device locates at [%d, %d, %d]", cell.gpu_cell.ivar_map[5], cell.gpu_cell.ivar_map[6], cell.gpu_cell.ivar_map[7]); show_log(msg, true);
	sprintf(msg, "==>  center atom on host locates at [%d, %d, %d]", cell.gpu_cell.ir_map[3], cell.gpu_cell.ir_map[4], cell.gpu_cell.ir_map[5]); show_log(msg, true);

	{
		int i, j, k, ia, *atom;
		for (i = 0; i < cell.Nx; i++) {for (j = 0; j < cell.Ny; j++) {for (k = 0; k < cell.Nz; k++) {
			atom = cell.cubic(i, j, k);
			if (atom[0] != 0) {
				sprintf(msg, "[%d, %d, %d]: ", i, j, k); show_log(msg, false);
				sprintf(msg, "%d -- ", atom[0]); show_log(msg, false);

				for (ia = 0; ia < atom[0]; ia++) {sprintf(msg, "%d [%f, %f, %f], ", atom[ia + 1], cell.gpu_cell.r_map[atom[ia+1] * 3], cell.gpu_cell.r_map[atom[ia+1] * 3 + 1], cell.gpu_cell.r_map[atom[ia+1] * 3 + 2]); show_log(msg, false);}
				show_log("", true);
			}
		}}}
	}
	*/

	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif


	// then we needs to dump the efield or force of local interaction
	SVECTOR<double, 4> sr_res[MAX_THREADS];
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::collect_efield_sr_mu_3d), job, nthreads, &cell, &vspme, true);
	else {
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::collect_fresult_sr_mu_3d), job, nthreads, &cell, &vspme, true);
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
bool Polarize_SPME_Interact(MMOL_MD_CELL &mdcell, EAtCELL &cell, _spme_::VSPME<_evatom_::_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field

	bool bEfieldOnly = var.bEfieldOnly;
	res.reset();
	InteractRes tres;

	init_normalized_coordinates(vspme, nThreads); // vspme has to be initialized before SPME_EF for atom with dipole
	double dmu_max = 0, vt = 0, muConvg = ::muConvg; // eA
	int iloop = 0, max_loops = ::nmax_loops_polarize;
	bool status = false;
	var.bEfieldOnly = true;
	cell.bCalE = 1; // electric field only
	cell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

	while (bPolarize && iloop < max_loops) {
		cal_dipole(mdcell, nThreads, mdcell.mu_surf); // accumulate dipole of each atom, and the whole mdcell
		memcpy(&(cell.mu_surf), &(mdcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
		_atomic_cmm_::update_atomic_charge_dipole(&cell);

		SPME_EF(cell, vspme, nThreads, var, tres, false); // calculate the dipole only

		// add the electric field from the atom inside the cluster
		//if (mdcell.sm.n > 0) MOperate<SMOLECULE>((void*)(&AtomicEFieldCorrect_SM), mdcell.sm.m, mdcell.sm.n, nThreads);

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

	if (bEfieldOnly) { // calculate f
		SPME_EF(cell, vspme, nThreads, var, res, false); // calculate the dipole and force
	}
	else { // calculate f
		var.bEfieldOnly = bEfieldOnly;
		cell.bCalE = (bEfieldOnly ? 1 : 2);

		SPME_EF(cell, vspme, nThreads, var, res, false); // calculate the dipole and force
	}

	return status;
}

// polarized SPME, with the interaction between induced dipole excluded
// in some polarized model, the interaction between induced dipole is excluded, e.g. only q-induced dipole is considered.
void VSPME_F_exclude_induced(MMOL_MD_CELL &mdcell, EAtCELL &cmm_cell, _spme_::VSPME<_evatom_::_VMUatom>& vspme_induced, int nThreads, SPME_VAR &var, InteractRes &res, bool bInitABSP_buff) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	// accumulate the electric field from the image part
	//vspme_induced.reset_local_EF(); // reset the variable for E & F
	init_normalized_coordinates(vspme_induced, nThreads); //-- vspme is initilized outside before this funciton is called
	if (!bInitABSP_buff) vspme_induced.bReady_absp = true; // copied from vspme
	cal_Q_spme(vspme_induced, nThreads); // SPME Q matrix
	double U = vspme_induced.cal_U_recp(!var.bEfieldOnly); //calculate the energy from the image part 
	res.U = U;
	cal_E_recp(vspme_induced, true, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (var.bVirial) {
		cal_StressTensor(vspme_induced, nThreads);
		if (var.bST_diagonal) {
			res.STxx += vspme_induced.STxx; res.STyy += vspme_induced.STyy; res.STzz += vspme_induced.STzz;
		}
		res.STtr += vspme_induced.STtr;
	}
/*
	if (nThreads > 1) {
		vk0.thread_var.WaitUntilThreadOver();
		res += tres_k0;
	}
*/
	//vspme.dump_EF(!var.K0var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
	//bool bDumpF = !var.K0var.bEfieldOnly;

	JobIndx job[MAX_THREADS];
	for (int i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	int nthreads;

	bool bDumpF = true; // useless
	bool bDumpE = false;
	//if (nThreads < 2 || vspme_induced.va.n < N4PARALLEL) Slab_dump_F_exclude_induced_mu(vspme_induced, 0, vspme_induced.va.n - 1, bDumpF); // transfer the calculated E & F to the variables in real atom
	//else MRunSys2< VSPME<_evatom_::_VMUatom>, bool >((void*)(&Slab_dump_F_exclude_induced_mu), vspme_induced, 0, vspme_induced.va.n, nThreads, bDumpF);
	if (nThreads < 2 || vspme_induced.va.n < N4PARALLEL) _spme_::dump_F_exclude_induced_mu(vspme_induced, 0, vspme_induced.va.n - 1, bDumpF); // transfer the calculated E & F to the variables in real atom
	else {
		//MRunSys2< _spme_2d_::VSPME<_evatom_::_VMUatom>, bool>((void*)(&_spme_2d_::dump_EF_MU), vspme_induced, 0, vspme_induced.va.n, nThreads, bDumpF);
		nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme_induced.va.n - 1, 400, nThreads);
		MTOperate1< JobIndx, _spme_::VSPME<_evatom_::_VMUatom> >((void*)(&_spme_::dumpJob_F_exclude_induced_mu), job, nthreads, &vspme_induced, true);
	}
}


bool Polarize_SPME_Interact2(MMOL_MD_CELL &mdcell, EAtCELL &cell, _spme_::VSPME<_evatom_::_VMUatom>& vspme, _spme_::VSPME<_evatom_::_VMUatom>& vspme_induced, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
	/*
	vspme_induced.release_plans();
	vspme.init_fft_plans();
	bool status = Polarize_SPME_Interact(mdcell, cmm_cell, vspme, nThreads, bPolarize, var, res);
	InteractRes tres;
	if (bPolarize && ::bExcludeInducedDipole) {
		vspme.release_plans();
		vspme_induced.init_fft_plans();
		VSPME_F_exclude_induced(mdcell, cmm_cell, vspme_induced, nThreads, var, tres);
		res += tres;
	}
	return status;
	*/

	vspme_induced.release_plans();
	vspme.init_fft_plans();
#if _CHECK_TIME_STAMP_
	show_log("Polarize: ", false);
	show_vlog<int>(mt.elapse(), false);
#endif
	bool status = Polarize_SPME_Interact(mdcell, cell, vspme, nThreads, bPolarize, var, res);
	InteractRes tres;
	bool bABSP_Buff = false;
	if (bPolarize && var.esv1.iExIndDipole == 1) {
		vspme.release_plans();
		vspme_induced.init_fft_plans();
		if (vspme.absp.n == vspme_induced.absp.n) {
			cp_absp_buff((__AtomBSP_BUFF*)(&vspme), (__AtomBSP_BUFF*)(&vspme_induced));
			bABSP_Buff = true;
		}
		else bABSP_Buff = false;
		VSPME_F_exclude_induced(mdcell, cell, vspme_induced, nThreads, var, tres, !bABSP_Buff);
		res -= tres;
	}
#if _CHECK_TIME_STAMP_
	show_log("", true);
#endif
	return status;
}


} // end of namespace


