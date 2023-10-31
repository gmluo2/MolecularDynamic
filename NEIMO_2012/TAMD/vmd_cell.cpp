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
using namespace _cmm_3d_;
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
#include "vmd_cell.h"

#include "complex.h"
#include "fftw3.h"
#include "EwaldSum.h"
#include "spme.h"
#include "spme_interact.h"
extern BSplineFunc<2> bsp;

using namespace _EwaldSum_real_;

void MMOL_VMD_CELL::check_freedom() {
	int NF = 0;
	int imol = 0, i;

	NF = 0; 
	for (i = 0; i < mm.n; i++) {
		NF += (mm.m[imol]->nCluster - 1);
		NF += 6;
	}
	// each sm has 6 freedomes, translation + rotation
	NF += sm.n * 6;
	// each sm has 3 freedomes, translation only
	NF += pm.n * 3;
	//for (i = 0; i < pm.n; i++) NF += 3; // translation only

	hdyn.nhc->set_freedom(NF);
	this->NF = NF;
	this->Ethermal = this->NF * ::Ek0Internal;
}

bool MMOL_VMD_CELL::MDSave(long nloop, float Ep, bool bEachMolDyn) {
	int nsave = 0, n = 0;
	MD_SAVE *mdsave = NULL;
	if (mm.n > 0) mdsave = lfmd.m[0]->mmsave.mdsave;
	else if (sm.n > 0) mdsave = lfsmd.m[0]->smsave.mdsave;
	else if (pm.n > 0) mdsave = lfpmd.m[0]->pmsave.mdsave;
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
		for (n = 0; n < pm.n; n++) {
			sprintf(title, "PM %d", n);
			lfpmd.m[n]->pmsave.save_kinetic_pars(*(pm.m[n]), nloop, false, title, save_md_title);
			lfpmd.m[n]->pmsave.mdsave->save_energy(pm.m[n]->Ek, nloop, false, NULL, save_md_title);
			//if (bEachMolDyn) lfpmd.m[n].save_mol_dyn_pars(pm.m[n], nloop, false, title, save_md_title);
			save_md_title = false;
		}
		for (n = 0; n < sm.n; n++) {
			sprintf(title, "SM %d", n);
			lfsmd.m[n]->smsave.save_kinetic_pars(*(sm.m[n]), nloop, false, title, save_md_title);
			lfsmd.m[n]->smsave.mdsave->save_energy(sm.m[n]->Ek, nloop, false, NULL, save_md_title);
			//if (bEachMolDyn) lfsmd.m[n].save_mol_dyn_pars(sm.m[n], nloop, false, title, save_md_title);
			save_md_title = false;
		}
		for (n = 0; n < mm.n; n++) {
			sprintf(title, "MM %d", n);
			lfmd.m[n]->mmsave.save_kinetic_pars(*(mm.m[n]), nloop, false, title, save_md_title);
			//lfmd.m[n].mmsave.mdsave->save_energy(mm.m[n].Ek, nloop, false, NULL, save_md_title);
			lfmd.m[n]->mmsave.save_energy(*(mm.m[n]), nloop, false, NULL, save_md_title);
			//if (bEachMolDyn) lfmd.m[n].save_mol_dyn_pars(mm.m[n], nloop, false, title, save_md_title);
			save_md_title = false;
		}
		mdsave->save_energy(Ep, nloop, false, NULL, false);

		if (!bEachMolDyn && mdsave->mDynSave.bsave()) {
			*(mdsave->mDynSave.out)<<"[MD] "<<mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
			//*(mdsave->mDynSave.out)<<hdyn.ksi_Hoover<<" "<<hdyn.vksi_Hoover<<endl<<endl;
			if (hdyn.nhc != NULL) {
				*(mdsave->mDynSave.out)<<n<<":  "<<hdyn.nhc->ksi.v[0]<<" "<<hdyn.nhc->vksi.v[0]<<endl;
			}
			*(mdsave->mDynSave.out)<<endl;
		}
		mdsave->mKinSave.save_over();
		mdsave->mESave.save_over();
		mdsave->mDynSave.save_over();
	}
	return true;
}

void MMOL_VMD_CELL::init_mdcell_velocity_memory() {
	int imol;
	mvb.mmv1.SetArray(mm.n); mvb.mmv2.SetArray(mm.n); 
	mvb.mm_dvmax.SetArray(mm.n); mvb.mm_dwmax.SetArray(mm.n); 
	mvb.mm_status.SetArray(mm.n); mvb.mm_sc_status.SetArray(mm.n);

	for (imol = 0; imol < mm.n; imol++) {
		mvb.mmv1.m[imol].cv.SetArray(mm.m[imol]->nCluster);
		mvb.mmv2.m[imol].cv.SetArray(mm.m[imol]->nCluster);
	}

	mvb.sm_dvmax.SetArray(sm.n); mvb.sm_dwmax.SetArray(sm.n);
	mvb.sm_status.SetArray(sm.n); mvb.sm_sc_status.SetArray(sm.n);

	mvb.pm_dvmax.SetArray(pm.n); 
	mvb.pm_status.SetArray(pm.n); mvb.pm_sc_status.SetArray(pm.n);
}

void MMOL_VMD_CELL::init_polarization_buff() {
	int imol, iatom, natom = 0, ic;
	MMOLECULE *pmm = NULL;
	ind_mu_mm.SetArray(mm.n);
	for (imol = 0; imol < mm.n; imol++) {
		pmm = mm.m[imol];
		ind_mu_mm.m[imol].m = pmm;
		natom = 0;
		for (ic = 0; ic < pmm->nCluster; ic++) natom += pmm->cluster[ic].nAtoms;
		ind_mu_mm.m[imol].mu.SetArray(natom);
	}
	ind_mu_sm.SetArray(sm.n);
	for (imol = 0; imol < sm.n; imol++) {
		ind_mu_sm.m[imol].m = sm.m[imol];
		ind_mu_sm.m[imol].mu.SetArray(sm.m[imol]->c->nAtoms);
	}
	ind_mu_pm.SetArray(pm.n);
	for (imol = 0; imol < pm.n; imol++) {
		ind_mu_pm.m[imol].mu.SetArray(pm.m[imol]->c->nAtoms);
	}
}

int nAtoms(MMOL_VMD_CELL &mdcell) {
	int na = 0;
	int im, ic;
	for (im = 0; im < mdcell.mm.n; im++) {
		for (ic = 0; ic < mdcell.mm.m[im]->nCluster; ic++) na += mdcell.mm.m[im]->cluster[ic].nAtoms;
	}
	for (im = 0; im < mdcell.sm.n; im++) na += mdcell.sm.m[im]->c->nAtoms;
	for (im = 0; im < mdcell.pm.n; im++) na += 1;
	return na;
}

void MMOL_VMD_CELL::set_atom_multipole(int mMP) {
	int imol, iatom, ic, imu;
	MMOLECULE *pmm = NULL;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	for (imol = 0; imol < mm.n; imol++) {
		pmm = mm.m[imol];
		for (ic = 0; ic < pmm->nCluster; ic++) {
			pc = pmm->cluster + ic;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom; 
				if (patom->par->alpha != 0) {
					patom->mMP = mMP;
					if (patom->eNeutral) patom->eNeutral = false;
				}
				else if (patom->c != 0) {patom->mMP = mMP; patom->eNeutral = false;}
				else {patom->mMP = 0; patom->eNeutral = true;}
			}
			pc->check_eneutral();
		}
	}
	for (imol = 0; imol < sm.n; imol++) {
		pc = sm.m[imol]->c;
		for (iatom = 0; iatom < sm.m[imol]->c->nAtoms; iatom++) {
			patom = pc->atom + iatom; 
			if (patom->par->alpha != 0) {
				patom->mMP = mMP;
				if (patom->eNeutral) patom->eNeutral = false;
			}
			else if (patom->c != 0) {patom->mMP = mMP; patom->eNeutral = false;}
			else {patom->mMP = 0; patom->eNeutral = true;}
		}
		pc->check_eneutral();
	}
	for (imol = 0; imol < pm.n; imol++) {
		pc = pm.m[imol]->c;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom; 
			if (patom->par->alpha != 0) {
				patom->mMP = mMP;
				if (patom->eNeutral) patom->eNeutral = false;
			}
			else if (patom->c != 0) {patom->mMP = mMP; patom->eNeutral = false;}
			else {patom->mMP = 0; patom->eNeutral = true;}
		}
		pc->check_eneutral();
	}
}

void cluster_check_periodic_cell(MMOL_VMD_CELL &mmcell) {
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m[nm];
		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 3; i++) {
				xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
				x = (float)(pc->cm.v[i]);
				PERIOD_CHECK(x, xl, xr, period, pc->dr_cmm.v[i])
				if (pc->dr_cmm.v[i] != 0) pc->bCMMShift = true;
			}
			pc->update_vp();
		}
	}
	/* not necessary, because sm has only one cluster
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m[nm];
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(sm->r0.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
	}
	*/
}

void molecule_check_periodic_cell(MMOL_VMD_CELL &mmcell) {
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m[nm];
		m->calMassCenter(false);
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(m->mDynKin.cm.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) m->shiftMol(ds);

		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 3; i++) {
				xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
				x = (float)(pc->cm.v[i]);
				PERIOD_CHECK(x, xl, xr, period, pc->dr_cmm.v[i])
				if (pc->dr_cmm.v[i] != 0) pc->bCMMShift = true;
			}
			pc->update_vp();
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m[nm];
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(sm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m[nm];
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(pm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) pm->shiftMol(ds);
	}

	//cluster_check_periodic_cell(mmcell);
}


void VMDCELL_molecule_PeriodicCell_check(MDCELL_ThreadJob<MMOL_VMD_CELL> *job) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0, imol = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		m = mdcell->mm.m[imol];
		m->calMassCenter(false);
		for (i = 0; i < 3; i++) {
			xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
			x = (float)(m->mDynKin.cm.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) m->shiftMol(ds);

		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 3; i++) {
				xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
				x = (float)(pc->cm.v[i]);
				PERIOD_CHECK(x, xl, xr, period, pc->dr_cmm.v[i])
				if (pc->dr_cmm.v[i] != 0) pc->bCMMShift = true;
			}
			pc->update_vp();
		}
	}

	char msg[256] = "\0";
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol]; 
		for (i = 0; i < 3; i++) {
			xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
			x = (float)(sm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
		/*
		if (sm->r.v[0] > mdcell->h[0] || sm->r.v[0] < -mdcell->h[0] ||
			sm->r.v[1] > mdcell->h[1] || sm->r.v[1] < -mdcell->h[1] ||
			sm->r.v[2] > mdcell->h[2] || sm->r.v[2] < -mdcell->h[2]){
				sprintf(msg, "strange: sm %d has coordinate out of cell [%f, %f, %f]", sm->r.v[0], sm->r.v[1], sm->r.v[2]);
				show_log(msg, true);
		}
		*/
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
		for (i = 0; i < 3; i++) {
			xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
			x = (float)(pm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) pm->shiftMol(ds);
	}
}


// calculate effective dipole moment of multi-macromolecules
void MM_cal_dipole(MMOLECULE **mm, int n, SURF_DIPOLE &mu_surf) {
	V3zero(mu_surf.qr) V3zero(mu_surf.mu0) V3zero(mu_surf.mu) mu_surf.q = 0;
	int imol, nc, iatom;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	float q = 0;
	VECTOR3 mu0, qr;
	for (imol = 0; imol < n; imol++) {
		for (nc = 0; nc < mm[imol]->nCluster; nc++) {
			pc = mm[imol]->cluster + nc;
			if (pc->eNeutral) continue;
			q = 0; V3zero(mu0) V3zero(qr)
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom;
				q += patom->c;
				qr.v[0] += patom->c * patom->r.v[0];
				qr.v[1] += patom->c * patom->r.v[1];
				qr.v[2] += patom->c * patom->r.v[2];
				if (patom->mMP > 0) {
					V3plusV3(patom->intr_mu, patom->ind_mu, patom->mu) // sum the intrinsic dipole with the induced dipole
					mu0.v[0] += patom->mu.v[0];
					mu0.v[1] += patom->mu.v[1];
					mu0.v[2] += patom->mu.v[2];
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
void SM_cal_dipole(SMOLECULE **sm, int n, SURF_DIPOLE &mu_surf) {
	V3zero(mu_surf.qr) V3zero(mu_surf.mu0) V3zero(mu_surf.mu) mu_surf.q = 0;
	int imol, iatom;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	float q = 0;
	VECTOR3 qr, mu0;
	for (imol = 0; imol < n; imol++) {
		pc = sm[imol]->c;
		if (pc->eNeutral) continue;
		q = 0; V3zero(qr) V3zero(mu0)
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom;
			q += patom->c;
			qr.v[0] += patom->c * patom->r.v[0];
			qr.v[1] += patom->c * patom->r.v[1];
			qr.v[2] += patom->c * patom->r.v[2];
			if (patom->mMP > 0) {
				V3plusV3(patom->intr_mu, patom->ind_mu, patom->mu) // sum the intrinsic dipole with the induced dipole
				mu0.v[0] += patom->mu.v[0];
				mu0.v[1] += patom->mu.v[1];
				mu0.v[2] += patom->mu.v[2];
			}
		}
		mu_surf.q += q;
		V3plusV3(mu_surf.qr, qr, mu_surf.qr)
		V3plusV3(mu_surf.mu0, mu0, mu_surf.mu0)
	}
	V3plusV3(mu_surf.qr, mu_surf.mu0, mu_surf.mu)
}

// calculate effective dipole moment of multi-macromolecules
void PM_cal_dipole(PMOLECULE **pm, int n, SURF_DIPOLE &mu_surf) {
	V3zero(mu_surf.qr) V3zero(mu_surf.mu0) V3zero(mu_surf.mu) mu_surf.q = 0;
	int imol;
	float c;
	double *r;
	MATOM *patom = NULL;
	float q = 0;
	VECTOR3 qr, mu0;
	for (imol = 0; imol < n; imol++) {
		patom = pm[imol]->c->atom; c = pm[imol]->c->atom->c;
		r = pm[imol]->r.v;
		q += c;
		qr.v[0] += c * r[0];
		qr.v[1] += c * r[1];
		qr.v[2] += c * r[2];
		if (patom->mMP > 0) {
			patom = pm[imol]->c->atom;
			V3plusV3(patom->intr_mu, patom->ind_mu, patom->mu) // sum the intrinsic dipole with the induced dipole
			mu0.v[0] += patom->mu.v[0];
			mu0.v[1] += patom->mu.v[1];
			mu0.v[2] += patom->mu.v[2];
		}
	}
	mu_surf.q = q;
	V32V3(qr, mu_surf.qr)
	V32V3(mu0, mu_surf.mu0)
	V3plusV3(mu_surf.qr, mu_surf.mu0, mu_surf.mu)
}

void cal_dipole(MMOL_VMD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf) {
	V3zero(mu_surf.qr) V3zero(mu_surf.mu0) V3zero(mu_surf.mu) mu_surf.q = 0;
	if (::iSurfaceBoundary == 1) return;
	SURF_DIPOLE mu_sf;
	if (mmcell.mm.m != NULL) {
		MOperate2<MMOLECULE*, SURF_DIPOLE>((void*)(&MM_cal_dipole), mmcell.mm.m, mmcell.mm.n, nThreads, mu_sf);
		mu_surf.q += mu_sf.q; 
		V3plusV3(mu_surf.qr, mu_sf.qr, mu_surf.qr)
		V3plusV3(mu_surf.mu0, mu_sf.mu0, mu_surf.mu0)
	}
	if (mmcell.sm.m != NULL) {
		MOperate2<SMOLECULE*, SURF_DIPOLE>((void*)(&SM_cal_dipole), mmcell.sm.m, mmcell.sm.n, nThreads, mu_sf);
		mu_surf.q += mu_sf.q; 
		V3plusV3(mu_surf.qr, mu_sf.qr, mu_surf.qr)
		V3plusV3(mu_surf.mu0, mu_sf.mu0, mu_surf.mu0)
	}
	if (mmcell.pm.m != NULL) {
		MOperate2<PMOLECULE*, SURF_DIPOLE>((void*)(&PM_cal_dipole), mmcell.pm.m, mmcell.pm.n, nThreads, mu_sf);
		mu_surf.q += mu_sf.q; 
		V3plusV3(mu_surf.qr, mu_sf.qr, mu_surf.qr)
		V3plusV3(mu_surf.mu0, mu_sf.mu0, mu_surf.mu0)
	}
	V3plusV3(mu_surf.qr, mu_surf.mu0, mu_surf.mu)

	if (::iSurfaceBoundary == 2) { // disable z direction
		mu_surf.qr.v[2] = 0; mu_surf.mu0.v[2] = 0; mu_surf.mu.v[2] = 0; mu_surf.mu_ind.v[2] = 0;
	}
}

void polarize(MMOL_VMD_CELL &mmcell, int nThreads) {
	if (mmcell.mm.m != NULL) {
		MOperate<MMOLECULE*>((void*)(&polarize_MM), mmcell.mm.m, mmcell.mm.n, nThreads);
	}
	if (mmcell.sm.m != NULL) {
		MOperate<SMOLECULE*>((void*)(&polarize_SM), mmcell.sm.m, mmcell.sm.n, nThreads);
	}
	if (mmcell.pm.m != NULL) {
		MOperate<PMOLECULE*>((void*)(&polarize_PM), mmcell.pm.m, mmcell.pm.n, nThreads);
	}
}

double cal_charge(MMOL_VMD_CELL &mmcell) {
	double q = 0;
	double qm = 0, qc = 0;
	int imol, nc, iatom;
	BASIC_CLUSTER *pc = NULL;
	for (imol = 0; imol < mmcell.mm.n; imol++) {
		qm = 0;
		for (nc = 0; nc < mmcell.mm.m[imol]->nCluster; nc++) {
			pc = mmcell.mm.m[imol]->cluster + nc;
			qc = 0;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) qc += pc->atom[iatom].c;
			qm += qc;
		}
		q += qm;
	}
	for (imol = 0; imol < mmcell.sm.n; imol++) {
		pc = mmcell.sm.m[imol]->c;
		qc = 0;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) qc += pc->atom[iatom].c;
		q += qc;
	}
	for (imol = 0; imol < mmcell.pm.n; imol++) {
		pc = mmcell.pm.m[imol]->c;
		q += pc->atom[0].c;
	}
	return q;
}
