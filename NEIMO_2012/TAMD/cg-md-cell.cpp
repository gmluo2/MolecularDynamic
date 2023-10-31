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
#include "cg-mm.h"
#include "CMM.h"
using namespace _cmm_3d_;
#include "cg-cmm.h"
#include "var.h"

#include "Interaction.h"
#include "Interact1.h"
#include "Interact2.h"
#include "Interact3.h"
#include "Interact4.h"
#include "MD.h"
#include "cg-md.h"
//#include "MD_POLYMER.h"
#include "cgmd_polymer.h"
#include "md_cell.h"
#include "cg-md-cell.h"

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

extern void show_all_molecules();

namespace _coarse_grain_ {

void CGMM_MD_CELL::init_nhc_mols() { // group molecules for each HOOVER_DYN
	int nMM = 0, nSM = 0, nPM = 0;
	int imol = 0, ihdyn = 0, indx = 0;
	for (ihdyn = 0; ihdyn < hdyn.n; ihdyn++) {
		nMM = 0; nSM = 0; nPM = 0;
		for (imol = 0; imol < lfmd.n; imol++) {
			if (lfmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) nMM++;
		}
		for (imol = 0; imol < lfsmd.n; imol++) {
			if (lfsmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) nSM++;
		}
		for (imol = 0; imol < lfpmd.n; imol++) {
			if (lfpmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) nPM++;
		}
		hdyn.m[ihdyn].mmIndx.SetArray(nMM);
		if (nMM > 0) {
			indx = 0;
			for (imol = 0; imol < lfmd.n; imol++) {
				if (lfmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) {
					hdyn.m[ihdyn].mmIndx.m[indx] = imol; indx++;
				}
			}
		}
		hdyn.m[ihdyn].smIndx.SetArray(nSM);
		if (nSM > 0) {
			indx = 0;
			for (imol = 0; imol < lfsmd.n; imol++) {
				if (lfsmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) {
					hdyn.m[ihdyn].smIndx.m[indx] = imol; indx++;
				}
			}
		}
		hdyn.m[ihdyn].pmIndx.SetArray(nPM);
		if (nPM > 0) {
			indx = 0;
			for (imol = 0; imol < lfpmd.n; imol++) {
				if (lfpmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) {
					hdyn.m[ihdyn].pmIndx.m[indx] = imol; indx++;
				}
			}
		}
	}
}

void CGMM_MD_CELL::check_freedom() {
	int NF = 0;
	int imol = 0, ihdyn = 0, i;
	MDCELL_HDYN *pdyn = NULL;
	this->NF = 0;
	for (ihdyn = 0; ihdyn < hdyn.n; ihdyn++) {
		pdyn = hdyn.m + ihdyn;
		NF = 0;
		for (i = 0; i < pdyn->mmIndx.n; i++) {
			imol = pdyn->mmIndx.m[i];
			NF += mm.m[imol].nCluster * 3;
		}
		// each sm has 3 freedomes, translation
		NF += pdyn->smIndx.n *3;
		// each sm has 3 freedomes, translation only
		NF += pdyn->pmIndx.n * 3;
		//for (i = 0; i < pm.n; i++) NF += 3; // translation only

		pdyn->nhc->set_freedom(NF);
		this->NF += NF;
	}
	this->Ethermal = this->NF * ::Ek0Internal;
}

bool CGMM_MD_CELL::MDSave(long nloop, float Ep) {
	int nsave = 0, n = 0;
	MD_SAVE *mdsave = NULL;
	if (mm.n > 0) mdsave = lfmd.m[0].mmsave.mdsave;
	else if (sm.n > 0) mdsave = lfsmd.m[0].smsave.mdsave;
	else return false;
	nsave = mdsave->nsave;

	char title[20] = "\0";
	bool save_md_title = true;
	bool save_status = false;

	if (!mdsave->mKinSave.one_more()) return false;
	if (!mdsave->mESave.one_more()) return false;
	if (!mdsave->mDynSave.one_more()) return false;

	if (mdsave->mKinSave.bsave()) {
		save_md_title = true;
		for (n = 0; n < sm.n; n++) {
			sprintf(title, "SM %d", n);
			lfsmd.m[n].smsave.save_kinetic_pars(sm.m[n], nloop, false, title, save_md_title);
			lfsmd.m[n].smsave.mdsave->save_energy(sm.m[n].Ek, nloop, false, NULL, save_md_title);
			lfsmd.m[n].save_mol_dyn_pars(sm.m[n], nloop, false, title, save_md_title);
			save_md_title = false;
		}
		for (n = 0; n < mm.n; n++) {
			sprintf(title, "MM %d", n);
			lfmd.m[n].mmsave.save_kinetic_pars(mm.m[n], nloop, false, title, save_md_title);
			lfmd.m[n].mmsave.mdsave->save_energy(mm.m[n].Ek, nloop, false, NULL, save_md_title);
			lfmd.m[n].save_mol_dyn_pars(mm.m[n], nloop, false, title, save_md_title);
			save_md_title = false;
		}
		mdsave->save_energy(Ep, nloop, false, NULL, false);
		mdsave->mKinSave.save_over();
		mdsave->mESave.save_over();
		mdsave->mDynSave.save_over();
	}
	return true;
}

void set_free_cell_mass_center(CGMM_MD_CELL &mmcell, VECTOR3 &O) {
	VECTOR3 mc;
	float M = 0;
	CG_MMOLECULE *m = NULL;
	int nm = 0, i = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		m->calMassCenter();
		for (i = 0; i < 3; i++) {
			mc.v[i] += m->M * m->cm.v[i];
			M += m->M;
		}
	}
	CGSM *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 3; i++) {
			mc.v[i] += *(sm->M) * sm->r->v[i];
			M += *(sm->M);
		}
	}
	for (i = 0; i < 3; i++) mc.v[i] /= M;
	
	VECTOR3 ds;
	for (i = 0; i < 3; i++) ds.v[i] = O.v[i] - mc.v[i];
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		m->shiftMol(ds);
	}
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		sm->shift(ds); 
	}
}

void cluster_check_periodic_cell(CGMM_MD_CELL &mmcell) {
	CG_MMOLECULE *m = NULL;
	CG_CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 3; i++) {
				xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
				x = (float)(pc->r->v[i]);
				PERIOD_CHECK(x, xl, xr, period, pc->dr_cmm.v[i])
				if (pc->dr_cmm.v[i] != 0) pc->bCMMShift = true;
			}
			pc->update_vp();
		}
	}
	/* not necessary, because sm has only one cluster
	CGSM *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(sm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shift(ds);
	}
	*/
}

void molecule_check_periodic_cell(CGMM_MD_CELL &mmcell) {
	CG_MMOLECULE *m = NULL;
	int nm = 0, i = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		m->calMassCenter();
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(m->cm.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) m->shiftMol(ds);
	}
	CGSM *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(sm->r->v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shift(ds);
	}

	cluster_check_periodic_cell(mmcell);
}

#if _SYS_ == _WINDOWS_SYS_
int CGSM_LFVERLET_MD_thread_func(LPVOID *par) {
#elif _SYS_ == _LINUX_SYS_
void* CGSM_LFVERLET_MD_thread_func(void *par) {
#endif
	CGSM_MD_THREAD_PAR *sm_mdpar = (CGSM_MD_THREAD_PAR*)par;
	LFVERLET_CGSM_MD *lfsmdpar = NULL;
	CGSM *sm = NULL;
	int nm = 0;

	strcpy(sm_mdpar->msg_show, "\0"); strcpy(sm_mdpar->msg_log, "\0");

	double Ethermal = 3 * ::Ek0Internal;

	for (nm = 0; nm < sm_mdpar->nM; nm++) {
		sm = sm_mdpar->sm + nm;

		lfsmdpar = sm_mdpar->lfsmdpar + nm;
		CGSMLeapFrog_Verlet(*sm, *lfsmdpar, Ethermal, (sm_mdpar->md_loop == 0 ? true : false), sm_mdpar->msg_show, sm_mdpar->msg_log, 5120);
	}
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(sm_mdpar->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void CGSM_MULTI_THREAD_MD(CGSM_MD_THREAD_PAR* sm_mdpar, int nThreads) { // we assume nThread > 1
	int i = 0;
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nThreads; i++) {
		ResetEvent(sm_mdpar[i].hMMThread);
		sm_mdpar[i].thread = AfxBeginThread((AFX_THREADPROC)CGSM_LFVERLET_MD_thread_func, (LPVOID)(sm_mdpar + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nThreads; i++) {
		WaitForSingleObject(sm_mdpar[i].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nThreads; i++) {
		pthread_create(&(sm_mdpar[i].thread), NULL, &CGSM_LFVERLET_MD_thread_func, (void *)(sm_mdpar + i));
	}
	for (i = 0; i < nThreads; i++) {
		if (sm_mdpar[i].thread != 0) pthread_join(sm_mdpar[i].thread, NULL); // wait the thread over
	}
#endif
}

#if _SYS_ == _WINDOWS_SYS_
int CGMM_LFVERLET_MD_thread_func(LPVOID *par) {
#elif _SYS_ == _LINUX_SYS_
void* CGMM_LFVERLET_MD_thread_func(void *par) {
#endif
	CGMM_MD_THREAD_PAR *mm_mdpar = (CGMM_MD_THREAD_PAR*)par;
	LFVERLET_CGMM_MD *lfmdpar = NULL;
	CG_MMOLECULE *mm = NULL;
	CG_CLUSTER *pc = NULL;
	int nm = 0, nc = 0, j = 0;

	double Ethermal = 3 * ::Ek0Internal; // kT
	
	strcpy(mm_mdpar->msg_show, "\0"); strcpy(mm_mdpar->msg_log, "\0");
	for (nm = 0; nm < mm_mdpar->nM; nm++) {
		lfmdpar = mm_mdpar->lfmdpar + nm;
		mm = mm_mdpar->mm + nm;
		CGLeapFrog_Verlet(*mm, *lfmdpar, mm_mdpar->vbf, Ethermal, (mm_mdpar->md_loop == 0 ? true : false), 
			mm_mdpar->msg_show, mm_mdpar->msg_log, 5120); // translation and rotation of whole molecule are 1.5kT
	}
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(mm_mdpar->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void CGMM_MULTI_THREAD_MD(CGMM_MD_THREAD_PAR* mm_mdpar, int nThreads) { // we assume nThread > 1
	int i = 0;
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nThreads; i++) {
		ResetEvent(mm_mdpar[i].hMMThread);
		mm_mdpar[i].thread = AfxBeginThread((AFX_THREADPROC)CGMM_LFVERLET_MD_thread_func, (LPVOID)(mm_mdpar + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nThreads; i++) {
		WaitForSingleObject(mm_mdpar[i].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nThreads; i++) {
		pthread_create(&(mm_mdpar[i].thread), NULL, &CGMM_LFVERLET_MD_thread_func, (void *)(mm_mdpar + i));
	}
	for (i = 0; i < nThreads; i++) {
		if (mm_mdpar[i].thread != 0) pthread_join(mm_mdpar[i].thread, NULL); // wait the thread over
	}
#endif
}

void InitCGClusterForceInCMM(CMM_CELL3D<CG_CLUSTER> &cell, int n1Cell, int n2Cell) {
	CMM_CELL<CG_CLUSTER>* pcell = NULL;
	CG_CLUSTER *pc = NULL;

	ITERATOR<CG_CLUSTER, CMM_CELL<CG_CLUSTER> > it;

	for (int nCell = n1Cell; nCell <= n2Cell; nCell++) {
		if (nCell >= cell.acell.n) break;
		pcell = cell.acell.m[nCell]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			if (pc->dyn != NULL) {
				V6zero(pc->dyn->fc)
				pc->dyn->overlap = false;
			}
			it.next_iterate();
		}
	}
}

void MD_PROC(CGMM_MD_CELL &mmcell, CMM_CELL3D<CG_CLUSTER> *cmm, int md_mode) {
	if (::bCGdipole) {
		sprintf(errmsg, "external E-field : %8.2f kV/mm", ::Eext); show_msg(errmsg);
	}

	VECTOR3 mass_center;
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0) {
		show_msg(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}
	if (md_mode == MULTI_CGMM_FREE_CELL_MD) {
		set_free_cell_mass_center(mmcell, mass_center);
	}
	else if (md_mode == MULTI_CGMM_PERIDIC_CELL_MD) {
		molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell
		// using Ewald_Sum to calculate interaction
	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		//_EwaldSum_recp_::InitEwaldSumPars(cmm->xd[0]* cmm->xd[1]* cmm->xd[2]);
//#error needs to be completed for EwaldSum, if it is required
	#endif
	}

	int i = 0, j = 0, nloop = 0, nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0;

	for (nm = 0; nm < mmcell.mm.n; nm++) {
		if (nm >= mmcell.mm.n) break;
		mmcell.mm.m[nm].setup_kin_dyn();
		mmcell.mm.m[nm].calBrownianPars();
		// bound-dipole
		if (::bCGdipole) mmcell.mm.m[nm].setup_bounds();
	}
	if (mmcell.sm.m != NULL) {
		for (nm = 0; nm <= mmcell.sm.n; nm++) {
			if (nm >= mmcell.sm.n) break;
			mmcell.sm.m[nm].setup_kin_dyn();
			mmcell.sm.m[nm].calBrownianPars();
		}
	}

	// since molecule mass center position could be shifted, distribute cluster in CMM again 
	bool reset_cell_chain = true;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		reset_cell_chain = (nm == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.mm.m[nm], *cmm, reset_cell_chain); // set base cluster at the center cell
	}
	reset_cell_chain = (mmcell.mm.n == 0 ? true : false);
	cmm_init_distribute_cluster(mmcell.sm.m, mmcell.sm.n, *cmm, reset_cell_chain); // add simple solid molecule into cmm
	cmm_init_subcell<CG_CLUSTER>(*cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
	// here we have to check the periodity of the cluster position
	if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************

	show_all_molecules();

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;

	bool Overlap = false;
	int nstart = 0, nEachThread = 0, nmol = 0;

	CG_MMOLECULE *mm = NULL;
	CG_CLUSTER *pc = NULL;
	LFVERLET_CGMM_MD *lfmdpar = NULL;
	int nMaxCluster = mmcell.max_ncluster();

	int nMMThreads = 0, nMMEachThread = 0;
	CGMM_MD_THREAD_PAR mm_mdpar[MAX_THREADS];

	if (mmcell.mm.n > 0) {
		nMMEachThread = mmcell.mm.n / MAX_THREADS;
		if (nMMEachThread > 0) nMMThreads = MAX_THREADS;
		else {nMMThreads = mmcell.mm.n; nMMEachThread = 1;}

		nstart = 0;
		for (i = 0; i < nMMThreads; i++) {
			if (i == nMMThreads - 1) nmol = mmcell.mm.n - nstart;
			else nmol = nMMEachThread;
			mm_mdpar[i].set(mmcell.mm.m + nstart, mmcell.lfmd.m + nstart, nmol, i);
			mm_mdpar[i].set_buffer_memory(nMaxCluster);
			nstart += nmol;
		}
	}

	CGSM *sm = NULL;
	LFVERLET_CGSM_MD *lfsmdpar = NULL;

	int nSMThreads = 0, nSMEachThread = 0;
	CGSM_MD_THREAD_PAR sm_mdpar[MAX_THREADS];

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
	}

	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mm = mmcell.mm.m + nm;
		InitVelocity(*mm, 3 * ::Ek0Internal);
		for (nc = 0; nc < mm->nCluster; nc++) {
			V3zero(mm->cluster[nc].dyn->fc)
			mm->cluster[nc].dyn->overlap = false;
		}
	}

	if (mmcell.sm.n > 0) {
		Ek =  3 * ::Ek0Internal;
		MOperate1<CGSM, double>((void*)(&InitCGSMVelocity), mmcell.sm.m, mmcell.sm.n, Ek, MAX_THREADS);
	}

	mmcell.Init(); //copy the initial torsion angle, speed to the mdpar

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	long LOOPS = 0;
	int save_indx = 0;
	MD_SAVE *mdsave = NULL;
	if (mmcell.mm.n > 0) {LOOPS = mmcell.lfmd.m[0].LOOPS; mdsave = mmcell.lfmd.m[0].mmsave.mdsave;}
	else if (mmcell.sm.n > 0) {LOOPS = mmcell.lfsmd.m[0].LOOPS; mdsave = mmcell.lfsmd.m[0].smsave.mdsave;}

	bool bCheckClusterInCMM = true;

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command
		if (iCheckCMM == 0) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) {
			switch(md_mode) {
			case MULTI_CGMM_PERIDIC_CELL_MD:
				// here we have to check the periodity of the molecule position
				molecule_check_periodic_cell(mmcell);
				break;
			case MULTI_CGMM_FREE_CELL_MD:
				set_free_cell_mass_center(mmcell, mass_center);
				break;
			}
		}

		for (nm = 0; nm < mmcell.mm.n; nm++) {
			mm = mmcell.mm.m + nm;
			calInertiaMoment(mm);
		}

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

		CellOperate<CG_CLUSTER>((void*)(&InitCGClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);

	// re-distribute / check clusters in CMM, use vp of each cluster only
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		if (nloop == 0) cmm_check_cluster<CG_CLUSTER>(*cmm);
		else if (bCheckClusterInCMM) cmm_check<CG_CLUSTER>(*cmm, MAX_THREADS);
	#if _INTRA_INTER_INTERACTS == 1
		if (nloop == 0 || bCheckClusterInCMM) {
			for (nm = 0; nm < mmcell.mm.n; nm++)
				mmcell.mm.m[nm].check_non_neighbor_sm(::r2cut_intramol);
		}
	#endif
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		show_infor(errmsg);
	#endif

		Ep = 0;
		switch(md_mode) {
		case MULTI_CGMM_PERIDIC_CELL_MD:
		#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		#if _EWALD_SUM == _MultiPole_EWALD_SUM
			calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, true, MAX_THREADS);
			// calculate the multipole of each cell
			calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, false, true, MAX_THREADS);
		#elif _EWALD_SUM == _SPME_EWALD_SUM
			//cal_dipole(mmcell, cmm->mu);
		#endif
		#endif

			// check relationship of each cluster, electrostatic and LJ
			if (iCheckCMM == 0) {
			#if _CHECK_TIME_STAMP_ == 1
				time.start();
			#endif
				// imethod == 0 : LJ interaction for all clusters, ignoring the neighbors
				// imethod == 1 : LJ interaction for the clusters in different neighbors
				// here we use imethod = 0
				CheckLocalRelship(*cmm, 0, cmm->bcell.nx, 0, MAX_THREADS); // setup relationship for each atom with CMM, with the bcell
			#if _CHECK_TIME_STAMP_ == 1
				sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
				show_infor(errmsg);
			#endif
			}
			Ep = 0;
			for (nm = 0; nm < mmcell.mm.n; nm++) {
				mm = mmcell.mm.m + nm;
			#if _DISABLE_ == 0
				//MMClusterInteraction<CG_MMOLECULE, CG_CLUSTER>((void*)(&FreeCellCMMInteraction_SingleCGMM), *cmm, *mm, 0, mm->nCluster, Ep_t, Overlap, false, MAX_THREADS);
				Ep += Ep_t;
			#endif
			}
			if (mmcell.sm.n > 0) {
			#if _DISABLE_ == 0
				MultiMolInteraction<CGSM, CG_CLUSTER>((void*)(&FreeCellCMMInteraction_CGSM), *cmm, mmcell.sm.m, mmcell.sm.n, Ep_t, Overlap, MAX_THREADS);
				Ep += Ep_t;
			#endif
			}

			break;

		case MULTI_CGMM_FREE_CELL_MD:
		#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		#if _EWALD_SUM == _MultiPole_EWALD_SUM
			calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, false, MAX_THREADS);
		#elif _EWALD_SUM == _SPME_EWALD_SUM
			//cal_dipole(mmcell, cmm->mu);
		#endif
		#endif
			// check relationship of each cluster, electrostatic and LJ
			if (iCheckCMM == 0) {
			#if _CHECK_TIME_STAMP_ == 1
				time.start();
			#endif
				CheckLocalRelship(*cmm, 0, cmm->bcell.nx, 0, MAX_THREADS); // setup relationship for each atom with CMM, with the bcell
			#if _CHECK_TIME_STAMP_ == 1
				sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
				show_infor(errmsg);
			#endif
			}

			Ep = 0;
			for (nm = 0; nm < mmcell.mm.n; nm++) {
				mm = mmcell.mm.m + nm;
			#if _DISABLE_ == 0
				MMClusterInteraction<CG_MMOLECULE, CG_CLUSTER>((void*)(&FreeCellCMMInteraction_SingleCGMM), *cmm, *mm, 0, mm->nCluster, Ep_t, Overlap, false, MAX_THREADS);
				Ep += Ep_t;
			#endif
			}
			if (mmcell.sm.n > 0) {
			#if _DISABLE_ == 0
				MultiMolInteraction<CGSM, CG_CLUSTER>((void*)(&FreeCellCMMInteraction_CGSM), *cmm, mmcell.sm.m, mmcell.sm.n, Ep_t, Overlap, MAX_THREADS);
				Ep += Ep_t;
			#endif
			}

			break;
		}
		// intra-molecule non-bound interaction
	#if _INTRA_INTER_INTERACTS == 1
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			mm = mmcell.mm.m + nm;
			MacroMoleculeOperate2<CG_MMOLECULE, double>((void*)(&cgmm_IntraMolInteract), *mm, 0, mm->nCluster, MAX_THREADS, Ep_t);
			Ep += Ep_t;
		}
	#endif

		// bound-interaction
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			mm = mmcell.mm.m + nm;
			MacroMoleculeOperate2<CG_MMOLECULE, double>((void*)(&cgbound_interact), *mm, 0, mm->nCluster, MAX_THREADS, Ep_t);
			Ep += Ep_t;
			if (::bCGdipole && ::Eext != 0) {
				mm->bound_dipole();
				MacroMoleculeOperate2<CG_MMOLECULE, double>((void*)(&calDipoleEffect), *mm, 0, mm->nBound, MAX_THREADS, Ep_t);
				Ep += Ep_t;
			}
		}

		if (bImplicitSolvent) {
			// neighbors of cluster for implicit solvation 
			for (nm = 0; nm < mmcell.mm.n; nm++) {
				CGMMClusterImpSolvNeighborsInCell(cmm, mmcell.mm.m + nm, 0, mmcell.mm.m[nm].nCluster - 1);
				// solvation energy
				MacroMolCellOperate((void*)(&CGMM_ImplicitSolvationForce), *cmm, mmcell.mm.m[nm], 0, mmcell.mm.m[nm].nCluster - 1, MAX_THREADS);
				for (nc = 0; nc < mm->nCluster; nc++) Ep += mm->cluster[nc].E_ImplicitSolvation;
			}
		}

		/*
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			mm = mmcell.mm.m + nm;
			//mm->calTorsion(Etorsion);
			MacroMoleculeOperate2<MMOLECULE, double>((void*)(&calTorsion), *mm, 0, mm->nCluster, MAX_THREADS, Etorsion);
			Ep += Etorsion;
		}
		*/
		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		if (scale_force) {
			for (nm = 0; nm < mmcell.mm.n; nm++) {
				mm = mmcell.mm.m + nm;
				MacroMoleculeOperate<CG_MMOLECULE>((void*)(&scale_force_cgcluster), *mm, 0, mm->nCluster, MAX_THREADS);
			}
			if (mmcell.sm.n > 0) {
				MOperate<CGSM>((void*)(&scale_force_cgsm), mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
			}
		}

		for (i = 0; i < nMMThreads; i++) mm_mdpar[i].set_loop(nloop);
		for (i = 0; i < nSMThreads; i++) sm_mdpar[i].set_loop(nloop);

		mmcell.MDSave(nloop, (float)Ep);

		CGMM_MULTI_THREAD_MD(mm_mdpar, nMMThreads);
		CGSM_MULTI_THREAD_MD(sm_mdpar, nSMThreads);

		for (i = 0; i < nMMThreads; i++) {
			if (strlen(mm_mdpar[i].msg_show) > 1) show_infor(mm_mdpar[i].msg_show, true);
			if (strlen(mm_mdpar[i].msg_log) > 1) mlog.show(mm_mdpar[i].msg_log, true);
		}
		for (i = 0; i < nSMThreads; i++) {
			if (strlen(sm_mdpar[i].msg_show) > 1) show_infor(sm_mdpar[i].msg_show, true);
			if (strlen(sm_mdpar[i].msg_log) > 1) mlog.show(sm_mdpar[i].msg_log, true);
		}

		Ek = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) Ek += mmcell.mm.m[nm].Ek;
		for (nm = 0; nm < mmcell.sm.n; nm++) Ek += mmcell.sm.m[nm].Ek;

		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek, (float)Ep);

		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "checking procedure takes %d ms", time.glance());
			show_infor(errmsg);
		#endif

		nloop++;
	}
	for (i = 0; i < nMMThreads; i++) mm_mdpar[i].reset();
	for (i = 0; i < nSMThreads; i++) sm_mdpar[i].reset();
}

} // end of namespace _coarse_grain_ 


