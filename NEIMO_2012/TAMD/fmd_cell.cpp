#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
#include "show.h"
#include "ranlib.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "var.h"

#include "Interaction.h"
#include "Interact1.h"
#include "Interact2.h"
#include "Interact3.h"
#include "NEIMO.h"
#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"
#include "MD_POLYMER.h"
#include "MD_SM.h"
#include "md_cell.h"

#include "fmd_cell.h"

#include "complex.h"
#include "fftw3.h"
#include "EwaldSum.h"
#include "spme.h"
#include "spme_interact.h"
extern BSplineFunc<2> bsp;

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

#if _SYS_ == _WINDOWS_SYS_
int PM_SpeedVerlet_FMD_thread_func(LPVOID *par) {
#elif _SYS_ == _LINUX_SYS_
void *PM_SpeedVerlet_FMD_thread_func(void *par) {
#endif
	PM_FMD_THREAD_PAR *pm_mdpar = (PM_FMD_THREAD_PAR*)par;
	LFVERLET_PM_MD *lfmdpar = NULL;
	PMOLECULE *pm = NULL;
	int nm = 0;

	bool speed_verlet = (pm_mdpar->md_loop == 0 ? false : true);

	strcpy(pm_mdpar->msg_show, "\0"); strcpy(pm_mdpar->msg_log, "\0");

	pm_mdpar->Ek_total = 0;
	for (nm = 0; nm < pm_mdpar->nM; nm++) {
		pm = pm_mdpar->pm + nm;

		lfmdpar = pm_mdpar->lfpmdpar + nm;
		if (pm_mdpar->res[nm] == true) pm_mdpar->res[nm] = PM_SpeedVerlet(pm_mdpar->local_relax, speed_verlet,
			*pm, *lfmdpar, pm_mdpar->dv[nm], pm_mdpar->msg_show, pm_mdpar->msg_log, 5120);
		else pm_mdpar->dv[nm] = 0;
		pm_mdpar->Ek_total += pm->Ek;
	}
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(pm_mdpar->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void PM_SPEEDVERLET_MULTI_THREAD_MD(PM_FMD_THREAD_PAR* pm_mdpar, int nThreads) { // we assume nThread > 1
	int i = 0;
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nThreads; i++) {
		ResetEvent(pm_mdpar[i].hMMThread);
		pm_mdpar[i].thread = AfxBeginThread((AFX_THREADPROC)PM_SpeedVerlet_FMD_thread_func, (LPVOID)(pm_mdpar + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nThreads; i++) {
		WaitForSingleObject(pm_mdpar[i].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nThreads; i++) {
		pthread_create(&(pm_mdpar[i].thread), NULL, &PM_SpeedVerlet_FMD_thread_func, (void *)(pm_mdpar + i));
	}
	for (i = 0; i < nThreads; i++) {
		if (pm_mdpar[i].thread != 0) pthread_join(pm_mdpar[i].thread, NULL); // wait the thread over
	}
#endif
}


#if _SYS_ == _WINDOWS_SYS_
int SM_SpeedVerlet_FMD_thread_func(LPVOID *par) {
#elif _SYS_ == _LINUX_SYS_
void *SM_SpeedVerlet_FMD_thread_func(void *par) {
#endif
	SM_FMD_THREAD_PAR *sm_mdpar = (SM_FMD_THREAD_PAR*)par;
	LFVERLET_SM_MD *lfsmdpar = NULL;
	SMOLECULE *sm = NULL;
	int nm = 0;

	bool speed_verlet = (sm_mdpar->md_loop == 0 ? false : true);

	strcpy(sm_mdpar->msg_show, "\0"); strcpy(sm_mdpar->msg_log, "\0");

	sm_mdpar->Ek_total = 0;
	for (nm = 0; nm < sm_mdpar->nM; nm++) {
		sm = sm_mdpar->sm + nm;
		sm->calInertiaMoment();

		lfsmdpar = sm_mdpar->lfsmdpar + nm;
		if (sm_mdpar->res[nm] == true) sm_mdpar->res[nm] = SM_SpeedVerlet(sm_mdpar->local_relax, speed_verlet,
			*sm, *lfsmdpar, sm_mdpar->dv[nm], sm_mdpar->dw[nm], sm_mdpar->msg_show, sm_mdpar->msg_log, 5120);
		else {sm_mdpar->dv[nm] = 0; sm_mdpar->dw[nm] = 0;}
		sm_mdpar->Ek_total += sm->Ek;
	}
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(sm_mdpar->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void SM_SPEEDVERLET_MULTI_THREAD_MD(SM_FMD_THREAD_PAR* sm_mdpar, int nThreads) { // we assume nThread > 1
	int i = 0;
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nThreads; i++) {
		ResetEvent(sm_mdpar[i].hMMThread);
		sm_mdpar[i].thread = AfxBeginThread((AFX_THREADPROC)SM_SpeedVerlet_FMD_thread_func, (LPVOID)(sm_mdpar + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nThreads; i++) {
		WaitForSingleObject(sm_mdpar[i].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nThreads; i++) {
		pthread_create(&(sm_mdpar[i].thread), NULL, &SM_SpeedVerlet_FMD_thread_func, (void *)(sm_mdpar + i));
	}
	for (i = 0; i < nThreads; i++) {
		if (sm_mdpar[i].thread != 0) pthread_join(sm_mdpar[i].thread, NULL); // wait the thread over
	}
#endif
}

#if _SYS_ == _WINDOWS_SYS_
int MM_SpeedVerlet_FMD_thread_func(LPVOID *par) {
#elif _SYS_ == _LINUX_SYS_
void *MM_SpeedVerlet_FMD_thread_func(void *par) {
#endif
	MM_FMD_THREAD_PAR *mm_mdpar = (MM_FMD_THREAD_PAR*)par;
	LFVERLET_MM_MD *lfmdpar = NULL;
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	int nm = 0;

	bool speed_verlet = (mm_mdpar->md_loop == 0 ? false : true);

	strcpy(mm_mdpar->msg_show, "\0"); strcpy(mm_mdpar->msg_log, "\0");
	mm_mdpar->Ek_total = 0;
	for (nm = 0; nm < mm_mdpar->nM; nm++) {
		mm = mm_mdpar->mm + nm;
		//if (bRandomForce) gen_add_random_force(*mm, 0, mm->nCluster, randf);
		//CalNetForce(*mm, 0, mm->nCluster); // the net force on molecular mass center

		MakeFreeBaseVelocityConsistent(*mm);

		lfmdpar = mm_mdpar->lfmdpar + nm;
		// if mm_mdpar->res[nm] is false, the speed was rescaled in last iteration, so we stop the iteration
		if (mm_mdpar->res[nm] == true) mm_mdpar->res[nm] = MM_SpeedVerlet(mm_mdpar->local_relax, speed_verlet, 
			*mm, *lfmdpar, mm_mdpar->vect1, mm_mdpar->vect2, mm_mdpar->dv[nm], mm_mdpar->dw[nm], mm_mdpar->msg_show, mm_mdpar->msg_log, 5120); // translation and rotation of whole molecule are 1.5kT
		else {mm_mdpar->dv[nm] = 0; mm_mdpar->dw[nm] = 0;}
		mm_mdpar->Ek_total += mm->Ek;
	}
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(mm_mdpar->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void MM_SPEEDVERLET_MULTI_THREAD_MD(MM_FMD_THREAD_PAR* mm_mdpar, int nThreads) { // we assume nThread > 1
	int i = 0;
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nThreads; i++) {
		ResetEvent(mm_mdpar[i].hMMThread);
		mm_mdpar[i].thread = AfxBeginThread((AFX_THREADPROC)MM_SpeedVerlet_FMD_thread_func, (LPVOID)(mm_mdpar + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nThreads; i++) {
		WaitForSingleObject(mm_mdpar[i].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nThreads; i++) {
		pthread_create(&(mm_mdpar[i].thread), NULL, &MM_SpeedVerlet_FMD_thread_func, (void *)(mm_mdpar + i));
	}
	for (i = 0; i < nThreads; i++) {
		if (mm_mdpar[i].thread != 0) pthread_join(mm_mdpar[i].thread, NULL); // wait the thread over
	}
#endif
}

bool SelfConsistent_VelocityAcceleration(MMOL_MD_CELL &mm_mdcell, MM_FMD_THREAD_PAR *mm_mdpar, int nMMThreads, SM_FMD_THREAD_PAR *sm_mdpar, int nSMThreads, PM_FMD_THREAD_PAR *pm_mdpar, int nPMThreads, double vConvg, double wConvg, float dt) {
	int nm, i;
	float Ek_total = 0;
	for (nm = 0; nm < mm_mdcell.mm.n; nm++) Ek_total += mm_mdcell.mm.m[nm].Ek;
	for (nm = 0; nm < mm_mdcell.sm.n; nm++) Ek_total += mm_mdcell.sm.m[nm].Ek;
	for (nm = 0; nm < mm_mdcell.pm.n; nm++) Ek_total += mm_mdcell.pm.m[nm].Ek;

	// important : the res in each thread for each molecule, needs to be initialized to true so that we can try to get consistent velocity and acceleration
	// res could be changed during the iteration, if the velocity and acceleration can not be consistent
	for (i = 0; i < nMMThreads; i++) {for (nm = 0; nm < mm_mdpar[i].nM; nm++) mm_mdpar[i].res[nm] = true;}
	for (i = 0; i < nSMThreads; i++) {for (nm = 0; nm < sm_mdpar[i].nM; nm++) sm_mdpar[i].res[nm] = true;}
	for (i = 0; i < nPMThreads; i++) {for (nm = 0; nm < pm_mdpar[i].nM; nm++) pm_mdpar[i].res[nm] = true;}

	double ksi = 0, vksi = 0, T = 0;
	HOOVER_DYN *hdyn = &(mm_mdcell.hdyn.m[0]);

	bool convergent = false;
	double dv_max = 0, dw_max = 0;
	int max_loops = 10, loop = 0;

	while (!convergent) {
		T = Ek_total / mm_mdcell.Ethermal;
		hdyn->nhc->set_dT(T - 1);
		hdyn->nhc->verlet_NHC();
		ksi = hdyn->nhc->ksi.v[0];

		for (i = 0; i < nMMThreads; i++) mm_mdpar[i].ksi = ksi;
		for (i = 0; i < nSMThreads; i++) sm_mdpar[i].ksi = ksi;
		for (i = 0; i < nPMThreads; i++) pm_mdpar[i].ksi = ksi;

		Ek_total = 0; dv_max = 0;
		if (mm_mdcell.mm.n > 0) {
			MM_SPEEDVERLET_MULTI_THREAD_MD(mm_mdpar, nMMThreads);
			for (i = 0; i < nMMThreads; i++) {
				Ek_total += mm_mdpar[i].Ek_total;
				for (nm = 0; nm < mm_mdpar[i].nM; nm++) {
					if (dv_max < mm_mdpar[i].dv[nm]) dv_max = mm_mdpar[i].dv[nm];
					if (dw_max < mm_mdpar[i].dw[nm]) dw_max = mm_mdpar[i].dw[nm];
				}
			}
		}
		if (mm_mdcell.sm.n > 0) {
			SM_SPEEDVERLET_MULTI_THREAD_MD(sm_mdpar, nSMThreads);
			for (i = 0; i < nSMThreads; i++) {
				Ek_total += sm_mdpar[i].Ek_total;
				for (nm = 0; nm < sm_mdpar[i].nM; nm++) {
					if (dv_max < sm_mdpar[i].dv[nm]) dv_max = sm_mdpar[i].dv[nm];
					if (dw_max < sm_mdpar[i].dw[nm]) dw_max = sm_mdpar[i].dw[nm];
				}
			}
		}
		if (mm_mdcell.pm.n > 0) {
			PM_SPEEDVERLET_MULTI_THREAD_MD(pm_mdpar, nPMThreads);
			for (i = 0; i < nPMThreads; i++) {
				Ek_total += pm_mdpar[i].Ek_total;
				for (nm = 0; nm < pm_mdpar[i].nM; nm++) {
					if (dv_max < pm_mdpar[i].dv[nm]) dv_max = sm_mdpar[i].dv[nm];
				}
			}
		}
		if (dv_max < vConvg && dw_max < wConvg) convergent = true;
		loop++;
		if (loop > max_loops) {convergent = false; break;}
	}
	hdyn->nhc->accept_nextstep();
	return convergent;
}

// calculate effective dipole moment of multi-macromolecules
void MM_cal_dipole(MMOLECULE *mm, int n, _EwaldSum_real_::SURF_DIPOLE &mu_surf) {
	V3zero(mu_surf.qr) V3zero(mu_surf.mu0) V3zero(mu_surf.mu) mu_surf.q = 0; V3zero(mu_surf.mu_ind)
	int imol, nc, iatom;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	float q = 0;
	VECTOR3 mu0, qr;
	for (imol = 0; imol < n; imol++) {
		for (nc = 0; nc < mm[imol].nCluster; nc++) {
			pc = mm[imol].cluster + nc;
			if (pc->eNeutral) continue;
			q = 0; V3zero(mu0) V3zero(qr)
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom;
				q += patom->c;
				qr.v[0] += patom->c * patom->rg.v[0];
				qr.v[1] += patom->c * patom->rg.v[1];
				qr.v[2] += patom->c * patom->rg.v[2];
				if (patom->mMP > 0) {
					V3plusV3(patom->intr_mu, patom->ind_mu, patom->mu) // sum the intrinsic dipole with the induced dipole
					mu0.v[0] += patom->mu.v[0];
					mu0.v[1] += patom->mu.v[1];
					mu0.v[2] += patom->mu.v[2];

					V3plusV3(mu_surf.mu_ind, patom->ind_mu, mu_surf.mu_ind)
				}
			}
			mu_surf.q += q;
			V3plusV3(mu_surf.qr, qr, mu_surf.qr)
			V3plusV3(mu_surf.mu0, mu0, mu_surf.mu0)
		}
	}
	V3plusV3(mu_surf.qr, mu_surf.mu0, mu_surf.mu)
}

// calculate effective dipole moment of multi-macromolecules
void SM_cal_dipole(SMOLECULE *sm, int n, _EwaldSum_real_::SURF_DIPOLE &mu_surf) {
	V3zero(mu_surf.qr) V3zero(mu_surf.mu0) V3zero(mu_surf.mu) mu_surf.q = 0; V3zero(mu_surf.mu_ind)
	int imol, iatom;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	float q = 0;
	VECTOR3 qr, mu0;
	for (imol = 0; imol < n; imol++) {
		pc = sm[imol].c;
		if (pc->eNeutral) continue;
		q = 0; V3zero(qr) V3zero(mu0)
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom;
			q += patom->c;
			qr.v[0] += patom->c * patom->rg.v[0];
			qr.v[1] += patom->c * patom->rg.v[1];
			qr.v[2] += patom->c * patom->rg.v[2];
			if (patom->mMP > 0) {
				V3plusV3(patom->intr_mu, patom->ind_mu, patom->mu) // sum the intrinsic dipole with the induced dipole
				mu0.v[0] += patom->mu.v[0];
				mu0.v[1] += patom->mu.v[1];
				mu0.v[2] += patom->mu.v[2];

				V3plusV3(mu_surf.mu_ind, patom->ind_mu, mu_surf.mu_ind)
			}
		}
		mu_surf.q += q;
		V3plusV3(mu_surf.qr, qr, mu_surf.qr)
		V3plusV3(mu_surf.mu0, mu0, mu_surf.mu0)
	}
	V3plusV3(mu_surf.qr, mu_surf.mu0, mu_surf.mu)
}

// calculate effective dipole moment of multi-macromolecules
void PM_cal_dipole(PMOLECULE *pm, int n, _EwaldSum_real_::SURF_DIPOLE &mu_surf) {
	V3zero(mu_surf.qr) V3zero(mu_surf.mu0) V3zero(mu_surf.mu) mu_surf.q = 0; V3zero(mu_surf.mu_ind)
	int imol;
	float c;
	double *r;
	MATOM *patom = NULL;
	float q = 0;
	VECTOR3 qr, mu0;
	for (imol = 0; imol < n; imol++) {
		patom = pm[imol].c->atom; c = pm[imol].c->atom->c;
		r = pm[imol].r.v;
		q += c;
		qr.v[0] += c * r[0];
		qr.v[1] += c * r[1];
		qr.v[2] += c * r[2];
		if (patom->mMP > 0) {
			patom = pm[imol].c->atom;
			V3plusV3(patom->intr_mu, patom->ind_mu, patom->mu) // sum the intrinsic dipole with the induced dipole
			mu0.v[0] += patom->mu.v[0];
			mu0.v[1] += patom->mu.v[1];
			mu0.v[2] += patom->mu.v[2];

			V3plusV3(mu_surf.mu_ind, patom->ind_mu, mu_surf.mu_ind)
		}
	}
	mu_surf.q = q;
	V32V3(qr, mu_surf.qr)
	V32V3(mu0, mu_surf.mu0)
	V3plusV3(mu_surf.qr, mu_surf.mu0, mu_surf.mu)
}

void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, _EwaldSum_real_::SURF_DIPOLE &mu_surf) {
	V3zero(mu_surf.qr) V3zero(mu_surf.mu0) V3zero(mu_surf.mu) V3zero(mu_surf.mu_ind) mu_surf.q = 0;
	if (iSurfaceBoundary == 1) return;
	_EwaldSum_real_::SURF_DIPOLE mu_sf;
	if (mmcell.mm.m != NULL) {
		MOperate2<MMOLECULE, _EwaldSum_real_::SURF_DIPOLE>((void*)(&MM_cal_dipole), mmcell.mm.m, mmcell.mm.n, nThreads, mu_sf);
		mu_surf.q += mu_sf.q; 
		V3plusV3(mu_surf.qr, mu_sf.qr, mu_surf.qr)
		V3plusV3(mu_surf.mu0, mu_sf.mu0, mu_surf.mu0)
		V3plusV3(mu_surf.mu_ind, mu_sf.mu_ind, mu_surf.mu_ind)
	}
	if (mmcell.sm.m != NULL) {
		mu_sf.reset();
		MOperate2<SMOLECULE, _EwaldSum_real_::SURF_DIPOLE>((void*)(&SM_cal_dipole), mmcell.sm.m, mmcell.sm.n, nThreads, mu_sf);
		mu_surf.q += mu_sf.q; 
		V3plusV3(mu_surf.qr, mu_sf.qr, mu_surf.qr)
		V3plusV3(mu_surf.mu0, mu_sf.mu0, mu_surf.mu0)
		V3plusV3(mu_surf.mu_ind, mu_sf.mu_ind, mu_surf.mu_ind)
	}
	if (mmcell.pm.m != NULL) {
		mu_sf.reset();
		MOperate2<PMOLECULE, _EwaldSum_real_::SURF_DIPOLE>((void*)(&PM_cal_dipole), mmcell.pm.m, mmcell.pm.n, nThreads, mu_sf);
		mu_surf.q += mu_sf.q; 
		V3plusV3(mu_surf.qr, mu_sf.qr, mu_surf.qr)
		V3plusV3(mu_surf.mu0, mu_sf.mu0, mu_surf.mu0)
		V3plusV3(mu_surf.mu_ind, mu_sf.mu_ind, mu_surf.mu_ind)
	}
	V3plusV3(mu_surf.qr, mu_surf.mu0, mu_surf.mu)

	if (::iSurfaceBoundary == 2) { // disable z direction
		mu_surf.qr.v[2] = 0; mu_surf.mu0.v[2] = 0; mu_surf.mu.v[2] = 0; mu_surf.mu_ind.v[2] = 0;
	}
}

void polarize(MMOL_MD_CELL &mmcell, int nThreads) {
	if (mmcell.mm.m != NULL) {
		MOperate<MMOLECULE>((void*)(&polarize_MM), mmcell.mm.m, mmcell.mm.n, nThreads);
	}
	if (mmcell.sm.m != NULL) {
		MOperate<SMOLECULE>((void*)(&polarize_SM), mmcell.sm.m, mmcell.sm.n, nThreads);
	}
	if (mmcell.pm.m != NULL) {
		MOperate<PMOLECULE>((void*)(&polarize_PM), mmcell.pm.m, mmcell.pm.n, nThreads);
	}
}

double cal_charge(MMOL_MD_CELL &mmcell) {
	double q = 0;
	double qm = 0, qc = 0;
	int imol, nc, iatom;
	BASIC_CLUSTER *pc = NULL;
	for (imol = 0; imol < mmcell.mm.n; imol++) {
		qm = 0;
		for (nc = 0; nc < mmcell.mm.m[imol].nCluster; nc++) {
			pc = mmcell.mm.m[imol].cluster + nc;
			qc = 0;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) qc += pc->atom[iatom].c;
			qm += qc;
		}
		q += qm;
	}
	for (imol = 0; imol < mmcell.sm.n; imol++) {
		pc = mmcell.sm.m[imol].c;
		qc = 0;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) qc += pc->atom[iatom].c;
		q += qc;
	}
	for (imol = 0; imol < mmcell.pm.n; imol++) {
		pc = mmcell.pm.m[imol].c;
		q += pc->atom[0].c;
	}
	return q;
}

extern void check_atom_rg(MMOL_MD_CELL &mmcell, int nThreads);

#if _DISABLE_ == 0
void FMD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm) {
	using namespace _spme_;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	InitEwaldSumPars(cmm->xd[0], cmm->xd[1], cmm->xd[2]);

	init_BSplineFunc<2>(::bsp, ::indx_bspline);
	VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	if (bPolarize) {
		vspme_mu.bsp = &(::bsp);
		init_spme(&mmcell, vspme_mu); // init the viral atoms
		vspme_mu.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_mu.xl[0] = cmm->xl[0]; vspme_mu.xl[1] = cmm->xl[1]; vspme_mu.xl[2] = cmm->xl[2];
		vspme_mu.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		vspme_mu.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();

		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
	}
	else { // charge only
		vspme_q.bsp = &(::bsp);
		init_spme(&mmcell, vspme_q); // init the viral atoms
		vspme_q.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_q.xl[0] = cmm->xl[0]; vspme_q.xl[1] = cmm->xl[1]; vspme_q.xl[2] = cmm->xl[2];
		vspme_q.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		vspme_q.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}
#endif

	mmcell.init_big_small_mm();

	mmcell.check_freedom();
	mmcell.Ethermal = ::Ek0Internal * mmcell.NF; 

	VECTOR3 mass_center;
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}

	molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell
	check_atom_rg(mmcell, MAX_THREADS);
	// use SPME for Ewald Sum

	int i = 0, j = 0, nloop = 0, nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0;

	for (nm = 0; nm < mmcell.mm.n; nm++) {
		if (nm >= mmcell.mm.n) break;
		mmcell.mm.m[nm].setup_kin_dyn();
	}
	if (mmcell.sm.m != NULL) {
		for (nm = 0; nm <= mmcell.sm.n; nm++) {
			if (nm >= mmcell.sm.n) break;
			mmcell.sm.m[nm].setup_kin_dyn();
		}
	}
	if (mmcell.pm.m != NULL) {
		for (nm = 0; nm <= mmcell.pm.n; nm++) {
			if (nm >= mmcell.pm.n) break;
			mmcell.pm.m[nm].setup_dyn();
		}
	}

	// since molecule mass center position could be shifted, distribute cluster in CMM again 
	//for (nm = 0; nm < mmcell.mm.n; nm++) CalMolInertiaMoment(mmcell.mm.m[nm]); // here we will calculate mass center
	//for (nm = 0; nm < mmcell.sm.n; nm++) mmcell.sm.m[nm].calInertiaMoment();
	bool reset_cell_chain = true;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		reset_cell_chain = (nm == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.mm.m[nm], *cmm, reset_cell_chain); // set base cluster at the center cell
	}
	reset_cell_chain = (mmcell.mm.n == 0 ? true : false);
	cmm_init_distribute_cluster(mmcell.sm.m, mmcell.sm.n, *cmm, reset_cell_chain); // add simple solid molecule into cmm
	cmm_init_subcell<BASIC_CLUSTER>(*cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
	// here we have to check the periodity of the cluster position
	cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************

	show_all_molecules();

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;

	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;

	bool Overlap = false;
	int nstart = 0, nEachThread = 0, nmol = 0;

	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	LFVERLET_MM_MD *lfmdpar = NULL;
	int nMaxCluster = mmcell.max_ncluster();

	int nMMThreads = 0, nMMEachThread = 0;
	MM_FMD_THREAD_PAR mm_mdpar[MAX_THREADS];

	float md_dt = 0;

	if (mmcell.mm.n > 0) {
		nMMEachThread = mmcell.mm.n / MAX_THREADS;
		if (nMMEachThread > 0) nMMThreads = MAX_THREADS;
		else {nMMThreads = mmcell.mm.n; nMMEachThread = 1;}

		nstart = 0;
		for (i = 0; i < nMMThreads; i++) {
			if (i == nMMThreads - 1) nmol = mmcell.mm.n - nstart;
			else nmol = nMMEachThread;
			mm_mdpar[i].set(mmcell.mm.m + nstart, mmcell.lfmd.m + nstart, ::local_relax, nmol, i);
			mm_mdpar[i].set_buffer_memory(nMaxCluster);
			nstart += nmol;
		}
		if (md_dt == 0) md_dt = mmcell.lfmd.m[0].dt;
	}

	SMOLECULE *sm = NULL;
	LFVERLET_SM_MD *lfsmdpar = NULL;

	int nSMThreads = 0, nSMEachThread = 0;
	SM_FMD_THREAD_PAR sm_mdpar[MAX_THREADS];

	if (mmcell.sm.n > 0) {
		nSMEachThread = mmcell.sm.n / MAX_THREADS;
		if (nSMEachThread > 0) nSMThreads = MAX_THREADS;
		else {nSMThreads = mmcell.sm.n; nSMEachThread = 1;}

		nstart = 0;
		for (i = 0; i < nSMThreads; i++) {
			if (i == nSMThreads - 1) nmol = mmcell.sm.n - nstart;
			else nmol = nSMEachThread;
			sm_mdpar[i].set(mmcell.sm.m + nstart, mmcell.lfsmd.m + nstart, nmol, i);
			nstart += nmol;
		}
		if (md_dt == 0) md_dt = mmcell.lfsmd.m[0].dt;
	}

	PMOLECULE *pm = NULL;
	LFVERLET_PM_MD *lfpmdpar = NULL;

	int nPMThreads = 0, nPMEachThread = 0;
	PM_FMD_THREAD_PAR pm_mdpar[MAX_THREADS];

	if (mmcell.sm.n > 0) {
		nPMEachThread = mmcell.pm.n / MAX_THREADS;
		if (nPMEachThread > 0) nPMThreads = MAX_THREADS;
		else {nPMThreads = mmcell.pm.n; nPMEachThread = 1;}

		nstart = 0;
		for (i = 0; i < nPMThreads; i++) {
			if (i == nPMThreads - 1) nmol = mmcell.pm.n - nstart;
			else nmol = nPMEachThread;
			pm_mdpar[i].set(mmcell.pm.m + nstart, mmcell.lfpmd.m + nstart, nmol, i);
			nstart += nmol;
		}
		if (md_dt == 0) md_dt = mmcell.lfpmd.m[0].dt;
	}

	double vConvg = 0.005 / md_dt, wConvg = 0.005 / md_dt; //0.017 is 1 arc degree;
	if (vConvg < 2e-4) vConvg = 2e-4;
	if (wConvg < 2e-4) wConvg = 2e-4;


	VECTOR<6> Iw0;

	int nConvg = 0;

	//char init_struct_fname[256] = "\0";
	mmcell.Ethermal = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mm = mmcell.mm.m + nm;
		//lfmdpar = mmcell.lfmd.m + nm;
		//sprintf(init_struct_fname, "%s.xyz", lfmdpar->msave.title);
		//show_molecule_struct(*mm, init_struct_fname);
		//sprintf(errmsg, "initial structure of molecule %d is saved in %s", nm, init_struct_fname);
		//show_infor(errmsg);
		//strcpy(errmsg, "\0");

		InitVelocity(*mm);
		mm->calClusterSolvationRadius(); // it uses the cluster's mass center position
		mm->resetDynPars(); // set force & torque to 0
	}

	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		InitVelocity(*sm);
		// set force & torque to 0
		V6zero(sm->c->dyn->fc);
	}

	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m + nm;
		InitVelocity(*pm);
		// set force & torque to 0
		V6zero(pm->c->dyn->fc);
	}

	CHAIN<short> *ch_cluster_id = NULL;
	strcpy(::errmsg, "\n"); show_log(::errmsg, true);
	sprintf(::errmsg, "--------------------------------------------------------"); show_log(::errmsg, true);
	sprintf(::errmsg, "  CLUSTER UID: ita, cd_translation, cd_rotation"); show_log(::errmsg, true);
	sprintf(::errmsg, "  Brownian force: f = -ita * V - cd_trans/rot * |V| * V"); show_log(::errmsg, true);
	sprintf(::errmsg, "  Brownian force unit: mass * Angs/fs^2"); show_log(::errmsg, true);
	sprintf(::errmsg, "--------------------------------------------------------"); show_log(::errmsg, true);
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mm = mmcell.mm.m + nm;
		mm->show_cluster_BrownianForcePars(&ch_cluster_id);
	}
	release_chain<short>(&ch_cluster_id, true);
	sprintf(::errmsg, "------------------------- OVER -------------------------"); show_log(::errmsg, true);
	strcpy(::errmsg, "\n"); show_log(::errmsg, true);


	mmcell.Init_MD(); //copy the initial torsion angle, speed to the mdpar

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	long LOOPS = 0;
	int save_indx = 0;
	MD_SAVE *mdsave = NULL;
	if (mmcell.mm.n > 0) {LOOPS = mmcell.lfmd.m[0].LOOPS; mdsave = mmcell.lfmd.m[0].mmsave.mdsave;}
	else if (mmcell.sm.n > 0) {LOOPS = mmcell.lfsmd.m[0].LOOPS; mdsave = mmcell.lfsmd.m[0].smsave.mdsave;}

	NEIMO_THREAD_VARS *neimo_matrix_cal_var = new NEIMO_THREAD_VARS[mmcell.mm.n];

	bool bCheckClusterInCMM = true;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	cmm->mu_surf.q = cal_charge(mmcell);
#endif

	double U_LJ = 0, U_estat = 0;
	InteractVar LJ_var; 
	LJ_var.bVirial = ::bVirial; LJ_var.bST_diagonal = true;
	SPME_VAR spme_var;
	spme_var.bVirial = ::bVirial; spme_var.bST_diagonal = true; spme_var.bEfieldOnly = false;
	spme_var.esv1.bEwald = true; spme_var.esv1.init_EwaldSum(cmm->xd[0] * cmm->xd[1] * cmm->xd[2], ::rcut_Ewd);
	spme_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	InteractRes LJ_res, spme_res;

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command
		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			// here we have to check the periodity of the molecule position
			molecule_check_periodic_cell(mmcell);
		}

		check_atom_rg(mmcell, MAX_THREADS);

		// calculate the inertia moment of each cluster and the whole molecule
		for (nm = 0; nm < mmcell.bmm.n; nm++) {
			mm = mmcell.bmm.m[nm];
			CalMolInertiaMoment(*mm);
			//neimo_matrix_cal_var[nm].set_pars(mm, 0);
			//Start_NEIMO_MATRIX_CALC(*mm, neimo_matrix_cal_var[nm]);
		}
		MOperate<MMOLECULE*>((void*)(&MultiMMP_NEIMO_MOMENT), mmcell.smm.m, mmcell.smm.n, MAX_THREADS);
		// sm inertia moment will be calculated in SM_LFVERLET_MD_thread_func before Newton law calculation

		/*
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			for (nc = 0; nc < mmcell.mm.m[nm].nCluster; nc++) {
				V6zero(mmcell.mm.m[nm].cluster[nc].dyn->fc)
				mmcell.mm.m[nm].cluster[nc].dyn->overlap = false;
			}
		}
		for (nm = 0; nm < mmcell.sm.n; nm++) {
			sm = mmcell.sm.m + nm;
			V6zero(sm->c->dyn->fc)
		}
		*/
		CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);

	// re-distribute / check clusters in CMM, use vp of each cluster only
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
		else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		show_infor(errmsg);
	#endif

		Ep = 0;
		// check relationship of each cluster, electrostatic and LJ
		if (iCheckCMM == 0) {
		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
			ParallelCheckRelshipInCells(*cmm, 0, cmm->acell.n - 1, MAX_THREADS); // setup relationship for each atom with CMM
		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
			show_infor(errmsg);
		#endif
		}

		// important : calculate LJ interaction first, here cluster force & torque are reset
		LJ_Interact(&mmcell, cmm, MAX_THREADS, LJ_var, LJ_res);
		U_LJ = LJ_res.U;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

		spme_var.bEfieldOnly = false;
		// important : calculate Coulomb secondly, here cluster force & torque are accumulated
		if (::bPolarize) {
			if (nloop > 3) Guess_induced_dipole(mmcell, MAX_THREADS);
			// calculat the total dipole of the whole cell
			cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
			if (!Polarize_SPME_Interact(&mmcell, cmm, vspme_mu, MAX_THREADS, true, spme_var, spme_res)) {
				sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
			}
			Backup_induced_dipole_hist(mmcell, MAX_THREADS);
		}
		else {
			// calculat the total dipole of the whole cell
			cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
			SPME_Interact(&mmcell, cmm, vspme_q, MAX_THREADS, spme_var, spme_res);
		}
		U_estat = spme_res.U;
#endif
		if (mmcell.mm.m != NULL) MOperate<MMOLECULE>(MM_TOTAL_FORCE, mmcell.mm.m, mmcell.mm.n, MAX_THREADS);
		if (mmcell.sm.m != NULL) MOperate<SMOLECULE>(SM_TOTAL_FORCE, mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
		if (mmcell.pm.m != NULL) MOperate<PMOLECULE>(PM_TOTAL_FORCE, mmcell.pm.m, mmcell.pm.n, MAX_THREADS);
		Ep = U_LJ + U_estat;

		/*
		if (scale_force) {
			CellOperate<BASIC_CLUSTER>((void*)(&SCALE_FORCE_CMM), *cmm, 0, cmm->acell.n, MAX_THREADS);
		}
		*/

		for (nm = 0; nm < mmcell.mm.n; nm++) {
			mm = mmcell.mm.m + nm;
			//mm->calTorsion(Etorsion);
			MacroMoleculeOperate2<MMOLECULE, double>((void*)(&calTorsion), *mm, 0, mm->nCluster, MAX_THREADS, Etorsion);
			Ep += Etorsion;
		}

		if (bImplicitSolvent) {
			// neighbors of cluster for implicit solvation 
			if (bCheckClusterInCMM) CellOperate<BASIC_CLUSTER>((void*)(&CMMCheckImpSolvNeighbors_Cell), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			// solvation energy
			CellOperate<BASIC_CLUSTER>((void*)(&CMM_ImplicitSolvationForce), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			for (nm = 0; nm < mmcell.mm.n; nm++) {
				mm = mmcell.mm.m + nm;
				for (nc = 0; nc < mm->nCluster; nc++) Ep += mm->cluster[nc].E_ImplicitSolvation;
			}
		}

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		for (i = 0; i < nMMThreads; i++) mm_mdpar[i].set_loop(nloop);
		for (i = 0; i < nSMThreads; i++) sm_mdpar[i].set_loop(nloop);
		for (i = 0; i < nPMThreads; i++) pm_mdpar[i].set_loop(nloop);

		mmcell.MDSave(nloop, (float)Ep, false);

		// wait until the neimo-matrix calculation is done
		for (nm = 0; nm < mmcell.mm.n; nm++) {
#if _SYS_ == _WINDOWS_SYS_
			WaitForSingleObject(neimo_matrix_cal_var[nm].hMMThread, INFINITE);
#elif _SYS_ == _LINUX_SYS_
			pthread_join(neimo_matrix_cal_var[nm].thread, NULL);
#endif
		}

		//MM_MULTI_THREAD_MD(mm_mdpar, nMMThreads);
		//SM_MULTI_THREAD_MD(sm_mdpar, nSMThreads);
		// local relax was set at beginning
		for (i = 0; i < nMMThreads; i++) mm_mdpar[i].md_loop = nloop;
		for (i = 0; i < nSMThreads; i++) sm_mdpar[i].md_loop = nloop;
		for (i = 0; i < nPMThreads; i++) pm_mdpar[i].md_loop = nloop;
		SelfConsistent_VelocityAcceleration(mmcell, mm_mdpar, nMMThreads, sm_mdpar, nSMThreads, pm_mdpar, nPMThreads, vConvg, wConvg, md_dt);

		for (i = 0; i < nMMThreads; i++) {
			if (strlen(mm_mdpar[i].msg_show) > 1) show_infor(mm_mdpar[i].msg_show, true);
			if (strlen(mm_mdpar[i].msg_log) > 1) mlog.show(mm_mdpar[i].msg_log, true);
		}
		for (i = 0; i < nSMThreads; i++) {
			if (strlen(sm_mdpar[i].msg_show) > 1) show_infor(sm_mdpar[i].msg_show, true);
			if (strlen(sm_mdpar[i].msg_log) > 1) mlog.show(sm_mdpar[i].msg_log, true);
		}

		Ep_t = Ep; Ek_t = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			Ek_t += mmcell.mm.m[nm].Ek;
			Verlet(mmcell.lfmd.m[nm], mmcell.mm.m[nm]); // Verlet algorithm to calculate the ta and tp for next timestep
		}
		for (nm = 0; nm < mmcell.sm.n; nm++) {
			Ek_t += mmcell.sm.m[nm].Ek;
			Verlet(mmcell.lfsmd.m[nm], mmcell.sm.m[nm]); // Verlet algorithm to calculate the ta and tp for next timestep
		}
		for (nm = 0; nm < mmcell.pm.n; nm++) {
			Ek_t += mmcell.pm.m[nm].Ek;
			Verlet(mmcell.lfpmd.m[nm], mmcell.pm.m[nm]); // Verlet algorithm to calculate the ta and tp for next timestep
		}

		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			for (nm = 0; nm < mmcell.mm.n; nm++) {
				mm = mmcell.mm.m + nm;
				lfmdpar = mmcell.lfmd.m + nm;
				check_angle(*mm, *lfmdpar);
			}
		}
		iCheckAngle++; if (iCheckAngle >= nCheckAngle) iCheckAngle = 0;

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mmcell.CalibrateTorsionAngle();

	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "checking procedure takes %d ms", time.glance());
		show_infor(errmsg);
	#endif

		nloop++;
	}
	for (i = 0; i < nMMThreads; i++) mm_mdpar[i].reset();
	for (i = 0; i < nSMThreads; i++) sm_mdpar[i].reset();
	if (neimo_matrix_cal_var != NULL) {delete[] neimo_matrix_cal_var; neimo_matrix_cal_var = NULL;}

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	vspme_q.release_plans();
	vspme_mu.release_plans();
#endif
}

#endif // _DISABLE_
