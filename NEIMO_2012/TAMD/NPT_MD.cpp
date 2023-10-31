#include "project.h"

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
#include "NPT_MD.h"

#include "complex.h"
#include "fftw3.h"
#include "EwaldSum.h"
#include "spme.h"
#include "spme_interact.h"
extern BSplineFunc<2> bsp;

using namespace _spme_;

void InitVelocity(ISOTROPIC_PDYN &pdyn, float Ek) {
	pdyn.v = 2 * Ek / pdyn.W;
	pdyn.v = sqrt(pdyn.v);
	pdyn.Ek = Ek;
}

int check_NF_pressure(MMOL_MD_CELL *mdcell) {
	int NF = 0;
	NF += mdcell->mm.n; // each macromolecule has only one translation
	NF += mdcell->sm.n + mdcell->pm.n;
	return NF * 3;
}

void MDCELL_Ekt(MDCELL_ThreadJob<ISOTROPIC_NPT_CELL> *job, double *Ekt) {
	int nm = 0, ic = 0, imol = 0, iatom = 0;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	MMOL_MD_CELL *mdcell = job->job->mdcell->mdcell;
	// RF is accumulated based on molecule
	double Vk = 0;

	// IMPORTANT: in this function, we did not subtract the local average speed for each molecule

	// As a rigid cluster, only the translation energy at mass center, and the force at mass center
	// contribute to the internal pressure. Rotation at mass center has no contribution to internal pressure

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		Vk += mm->Ekt; // only the translation energy at mass center has contribution to the internal pressure
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; 
		Vk += sm->Ekt;  // only the translation energy at mass center has contribution to the internal pressure
	}
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol; //patom = pm->c->atom;
		Vk += pm->Ekt;  // point molecule has only one atom, no correction is required on virial term
	}
	*Ekt = Vk;
}

void MDCELL_PressureForce(MDCELL_ThreadJob<ISOTROPIC_NPT_CELL> *job, double *v, int nJob) {
	ISOTROPIC_NPT_CELL *npt_cell = job->job->mdcell;
	int i;
	float Ek = 0, RF = 0;
	MTOperate1< MDCELL_ThreadJob<ISOTROPIC_NPT_CELL>, double >((void*)(&MDCELL_Ekt), job, nJob, v);
	Ek = 0;
	for (i = 0; i < nJob; i++) Ek += v[i]; // has unit in kT
	npt_cell->Ekt = Ek;
	RF = npt_cell->RF;
	//char msg[256] = "\0";
	//sprintf(msg, "RF : %f, Ekt : %f", RF * unit_Ek_kT, Ek); show_log(msg, true);

	float Vol = npt_cell->pdyn.V;
	npt_cell->pdyn.Pint = (RF + Ek * 2 / unit_Ek_kT) / (::fUnit_atm_u * Vol * 3); // unit -- atmosphere
	//npt_cell->pdyn.Pint = Ek * 2 / (unit_Ek_kT * ::fUnit_atm_u * Vol * 3); // unit -- atmosphere

	double dP = npt_cell->pdyn.Pint - npt_cell->Pext;
	/*
	double dPmax = 50000; // atmosphere
	if (FABS(dP) > dPmax) dP = sign(dP) * dPmax;
	if (dP > dPmax) dP = dPmax;
	*/

	// pdyn.fc has unit in kT, matching to the unit of barostat mass, kT * fs^2
	npt_cell->pdyn.fc = 3 * Vol * dP * ::fUnit_atm_u * unit_Ek_kT + (Ek + Ek) / npt_cell->NF_pressure;
	//npt_cell->pdyn.fc = 0;
	// the last term above is tricky, in the publication of MTK, sometime is 1/N * SUM(p^2/m), sometime is d/N * SUM(p^2/m)
	// see J. Chem. Phys. 101(5) 4177 (1994) and J. Chem. Phys. 115(4)1678 (2001)
}

void PressureForce(ISOTROPIC_NPT_CELL *npt_cell) {
	int i;
	float Ek = 0, RF = 0;
	Ek = npt_cell->Ekt;
	RF = npt_cell->RF;
	//char msg[256] = "\0";
	//sprintf(msg, "RF : %f, Ekt : %f", RF * unit_Ek_kT, Ek); show_log(msg, true);

	float Vol = npt_cell->pdyn.V;
	npt_cell->pdyn.Pint = (RF + Ek * 2 / unit_Ek_kT) / (::fUnit_atm_u * Vol * 3); // unit -- atmosphere
	//npt_cell->pdyn.Pint = Ek * 2 / (unit_Ek_kT * ::fUnit_atm_u * Vol * 3); // unit -- atmosphere
	
	double dP = npt_cell->pdyn.Pint - npt_cell->Pext;
	/*
	double dPmax = 50000; // atmosphere
	if (FABS(dP) > dPmax) dP = sign(dP) * dPmax;
	if (dP > dPmax) dP = dPmax;
	*/

	// pdyn.fc has unit in kT, matching to the unit of barostat mass, kT * fs^2
	npt_cell->pdyn.fc = 3 * Vol * dP * ::fUnit_atm_u * unit_Ek_kT + (Ek + Ek) / npt_cell->NF_pressure;
	//npt_cell->pdyn.fc = 0;
	// the last term above is tricky, in the publication of MTK, sometime is 1/N * SUM(p^2/m), sometime is d/N * SUM(p^2/m)
	// see J. Chem. Phys. 101(5) 4177 (1994) and J. Chem. Phys. 115(4)1678 (2001)
}

void cal_dyn(ISOTROPIC_PDYN &pdyn, float ksi) {
	pdyn.alpha = pdyn.fc / pdyn.W - ksi * pdyn.v;
	//double amax = sqrt(0.2 / pdyn.W) / 20;
	//if (FABS(pdyn.alpha) > amax) pdyn.alpha = sign(pdyn.alpha) * amax;
}

void P_NHC(ISOTROPIC_PDYN &pdyn, HOOVER_DYN &phdyn) {
	calKineticEnergy(pdyn);
	double T = pdyn.Ek / ::Ek0Internal;
	phdyn.nhc->set_dT(T - 1);
	phdyn.nhc->verlet_NHC();
}

// MD functions for NPT ISOTROPIC_PDYN

bool ISOTROPIC_P_MD_SAVE::save_kinetic_pars(ISOTROPIC_PDYN &pdyn, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	if (mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!mdsave->mKinSave.one_more()) return false;
	}
	if (!mdsave->mKinSave.bsave()) return true; // this loop will be not saved

	if (bAsync) {
		mdsave->init_buf();
		if (save_md_title) mdsave->mKinSave.sbuf<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			mdsave->mKinSave.sbuf<<additional_title<<endl;
		mdsave->mKinSave.sbuf<<pdyn.V<<"    "<<pdyn.Pint<<endl;
		mdsave->mKinSave.sbuf<<pdyn.eps<<" "<<pdyn.v<<" "<<pdyn.alpha<<"  "<<pdyn.Ek<<endl;

		mdsave->mKinSave.sbuf<<endl;

		mdsave->mKinSave.flush_buf();
	}
	else {
		if (save_md_title) *(mdsave->mKinSave.out)<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			*(mdsave->mKinSave.out)<<additional_title<<endl;
		*(mdsave->mKinSave.out)<<pdyn.V<<"    "<<pdyn.Pint<<endl;
		*(mdsave->mKinSave.out)<<pdyn.eps<<" "<<pdyn.v<<" "<<pdyn.alpha<<"  "<<pdyn.Ek<<endl;

		*(mdsave->mKinSave.out)<<endl;
	}

	if (check_status) mdsave->mKinSave.save_over();
	return true;
}

void Speed_Verlet(LFVERLET_ISOTROPIC_P_MD &mdpar, ISOTROPIC_PDYN &pdyn) {
	mdpar.v_ph = mdpar.v_mh + pdyn.alpha * mdpar.dt;
	mdpar.v0 = 0.5 * (mdpar.v_mh + mdpar.v_ph);
	pdyn.v = mdpar.v0;

	//double dt = mdpar.dt * 0.5; // half step
	//mdpar.eps_ph = mdpar.eps0 + (pdyn.v + pdyn.alpha * dt * 0.5) * dt; // half step of eps
	//mdpar.eps_p = mdpar.eps0 + (pdyn.v + pdyn.alpha * mdpar.dt * 0.5) * mdpar.dt; // full step of eps
}

void Speed_Kinetic(LFVERLET_ISOTROPIC_P_MD &mdpar, ISOTROPIC_PDYN &pdyn) {
	double dv = 0, dt = mdpar.dt * 0.5;

	mdpar.v0 = pdyn.v;
	dv = pdyn.alpha * dt;
	mdpar.v_mh = mdpar.v0 - dv;
	mdpar.v_ph = mdpar.v0 + dv;

	//dv = (pdyn.v + pdyn.alpha * dt * 0.5) * dt; // change of eps within half step
	//mdpar.eps0 = pdyn.eps;
	//mdpar.eps_ph = mdpar.eps0 + dv;
	//mdpar.eps_mh = mdpar.eps0 - dv;
	//mdpar.eps_p = mdpar.eps0 + (pdyn.v + pdyn.alpha * mdpar.dt * 0.5) * mdpar.dt; // full step of eps
}

void NPT_Verlet(LFVERLET_ISOTROPIC_P_MD &mdpar, ISOTROPIC_PDYN &pdyn) {
	double d = 0, hdt = mdpar.dt * 0.5;
	mdpar.v_ph = mdpar.v_mh + pdyn.alpha * mdpar.dt;
	mdpar.v0 = (mdpar.v_mh + mdpar.v_ph) * 0.5;

	d = mdpar.v_ph * mdpar.dt;
	pdyn.eps += d;
	pdyn.V = pdyn.V0 * exp(3 * pdyn.eps); // volumn in next step
	
	// guess the value for next step
	// speed of next step
	mdpar.v0 = 1.5 * mdpar.v_ph - 0.5 * mdpar.v_mh; 
	pdyn.v = mdpar.v0; mdpar.v_mh = mdpar.v_ph; 
}

bool LFVERLET_ISOTROPIC_P_MD::save_dyn_pars(ISOTROPIC_PDYN &pdyn, long nloop, bool check_status, char* additional_title, bool save_md_title) {
	if (psave.mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!psave.mdsave->mDynSave.one_more()) return false;
	}
	if (!psave.mdsave->mDynSave.bsave()) return true; // this loop will be not saved

	if (save_md_title) *(psave.mdsave->mDynSave.out)<<"[MD] "<<psave.mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
	if (additional_title != NULL && strlen(additional_title) > 0) 
		*(psave.mdsave->mDynSave.out)<<additional_title<<endl;

	if (hdyn != NULL) {
		*(psave.mdsave->mDynSave.out)<<this->hdyn->nhc->ksi.v[0]<<" "<<this->hdyn->nhc->vksi.v[0]<<endl;
	}
	*(psave.mdsave->mDynSave.out)<<endl;

	if (check_status) psave.mdsave->mDynSave.save_over();
	return true;
}



void NPT_Speed_Verlet(LFVERLET_MM_MD &mdpar, MMOLECULE &mm) {
	// use the velocity of the base cluster
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *pc = NULL;

	double hdt = mdpar.dt * 0.5;

	for (i = 0; i < 6; i++) {
		mdpar.v0.v[i] = mm.V0.v[i];
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + mm.alpha0.v[i] * mdpar.dt;
		mdpar.v0.v[i] = mdpar.v0_mh.v[i] + mm.alpha0.v[i] * hdt;
		mm.V0.v[i] = mdpar.v0.v[i];
		V62V6(mm.V0, mm.base->bkm->V0)
	}

	// check velocity of clusters except base 
	for (i = 0; i < mm.nCluster; i++) {
		// Leapfrog Verlet algorithm
		pc = mm.cluster + i;
		if (mm.free_base && mm.base == pc) {
			mdpar.tp_mh[i] = 0; mdpar.tp_ph[i] = 0; pc->km.tpp = 0; pc->km.tp = 0;
			continue;
		}
		mdpar.tp_ph[i] = mdpar.tp_mh[i] + pc->km.tpp * mdpar.dt;
		mdpar.tp[i] = mdpar.tp_mh[i] + pc->km.tpp * hdt;
		pc->km.tp = mdpar.tp[i];
	}

	calClusterVelocity0(mm);
	calSpatialMomentum0(mm, true);
	mm.CalMassCenterVelocity_MassCenterMomentum();
}

void NPT_Speed_Kinetic(LFVERLET_MM_MD &mdpar, MMOLECULE &mm) {
	// use the velocity of base cluster
	// suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *pc = NULL;

	double hdt = mdpar.dt * 0.5;

	double dv = 0;

	for (i = 0; i < 6; i++) {
		mdpar.v0.v[i] = mm.V0.v[i];
		dv = mm.alpha0.v[i] * hdt;
		mdpar.v0_mh.v[i] = mm.V0.v[i] - dv;
		mdpar.v0_ph.v[i] = mm.V0.v[i] + dv;
	}

	for (i = 0; i < mm.nCluster; i++) {
		// Leapfrog Verlet algorithm
		pc = mm.cluster + i;
		if (pc = mm.base) continue;
		mdpar.tp[i] = pc->km.tp;
		dv = pc->km.tpp * hdt;
		mdpar.tp_mh[i] = mdpar.tp[i] - dv;
		mdpar.tp_ph[i] = mdpar.tp[i] - dv;
	}
}

void NPT_Verlet(LFVERLET_MM_MD &mdpar, MMOLECULE &mm, double e_eps) {
	// use the velocity of base cluster
	// suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *pc = NULL;

	double hdt = mdpar.dt * 0.5;

	double delta = 0;
	VECTOR<6> d;
	VECTOR3 dr, drot;

	// move base cluster and guess velocity for next step
	// because the base velocity is re-configured with the spacial moment at mass center,
	// the base velocity is hard to follow the leaflog algorithem

	for (i = 0; i < 6; i++) {
		// guess base cluster velocity
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + mm.alpha0.v[i] * mdpar.dt;
		mdpar.v0.v[i] = mdpar.v0_mh.v[i] + mm.alpha0.v[i] * hdt;
		d.v[i] = (mdpar.v0.v[i] + mm.alpha0.v[i] * hdt) * mdpar.dt;
		
		// guess base cluster velocity for next step
		mdpar.v0.v[i] = (mdpar.v0_mh.v[i] + 1.5 * mm.alpha0.v[i] * mdpar.dt) * e_eps; // rescale the speed for next volumn
		mm.V0.v[i] = mdpar.v0.v[i];
		V62V6(mm.V0, mm.base->bkm->V0)
		mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i] * e_eps; //for next step, and rescale
	}
	V2wv(drot, dr, d)

	mm.set_base_move(dr, drot);

	for (i = 0; i < mm.nCluster; i++) {
		pc = mm.cluster + i;
		if (mm.free_base && pc == mm.base) continue;
		//delta = torsion_angle(p->parent, (BASIC_CLUSTER*)p); // check torsion angle
		// Leapfrog Verlet algorithm
		delta = mdpar.tp_ph[i] * mdpar.dt;
		pc->rot.dta = delta;
	}
	
	MacroMolMove(&mm);

	// guess the rotation speed for next step
	// rotation speed is NOT scaled
	for (i = 0; i < mm.nCluster; i++) {
		pc = mm.cluster + i;
		if (pc == mm.base) continue;
		mdpar.t[i] += pc->rot.dta;
		pc->km.ta = mdpar.t[i];
		
		//delta = torsion_angle(p->parent, p); // check torsion angle

		// guess new position velocity
		mdpar.tp[i] = 1.5 * mdpar.tp_ph[i] - 0.5 * mdpar.tp_mh[i];
		pc->km.tp = mdpar.tp[i];
		mdpar.tp_mh[i] = mdpar.tp_ph[i];
	}

	// should we rescale the coordinates of macromolecule's mass center, or the the coordinates of base cluster ?
	// it might be more reasonable to rescale the coordinates on the old coordinate, while not the new coordinates
	// So, rescaling the base-cluster's coordinate sounds more reasonable
	
	// rescale the molecular mass center
	/*
	mm.calMassCenter(true);
	dr.v[0] = (e_eps - 1) * mm.mDynKin.cm.v[0];
	dr.v[1] = (e_eps - 1) * mm.mDynKin.cm.v[1];
	dr.v[2] = (e_eps - 1) * mm.mDynKin.cm.v[2];
	mm.shiftMol(dr);
	*/
	
	// rescale the base cluster of macromolecule
	dr.v[0] = (e_eps - 1) * mm.base->cm.v[0];
	dr.v[1] = (e_eps - 1) * mm.base->cm.v[1];
	dr.v[2] = (e_eps - 1) * mm.base->cm.v[2];
	mm.shiftMol(dr);
	
}






void NPT_Speed_Verlet(LFVERLET_SM_MD &mdpar, SMOLECULE &sm) {
	// use the velocity of the base cluster
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	double hdt = mdpar.dt * 0.5;

	double *alpha = sm.c->bkm->alpha.v;
	double *V = sm.c->bkm->V0.v;

	for (i = 0; i < 6; i++) {
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + alpha[i] * mdpar.dt;
		mdpar.v0.v[i] = mdpar.v0_mh.v[i] + alpha[i] * hdt;
		V[i] = mdpar.v0.v[i];
	}
}

void NPT_Speed_Kinetic(LFVERLET_SM_MD &mdpar, SMOLECULE &sm) {
	// use the velocity of base cluster
	// suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	double hdt = mdpar.dt * 0.5;

	double *alpha = sm.c->bkm->alpha.v;
	double *V = sm.c->bkm->V0.v;

	double dv = 0;

	for (i = 0; i < 6; i++) {
		dv = alpha[i] * hdt;
		mdpar.v0_mh.v[i] = V[i] - dv;
		mdpar.v0_ph.v[i] = V[i] + dv;
		mdpar.v0.v[i] = V[i];
	}
}

void NPT_Verlet(LFVERLET_SM_MD &mdpar, SMOLECULE &sm, double e_eps) {
	// use the velocity of base cluster
	// suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *pc = NULL;

	double hdt = mdpar.dt * 0.5;

	double delta = 0;
	VECTOR<6> d;
	VECTOR3 dr, drot;

	double *alpha = sm.c->bkm->alpha.v;
	double *V = sm.c->bkm->V0.v;
	double *R = sm.r.v;

	// move base cluster and guess velocity for next step
	// because the base velocity is re-configured with the spacial moment at mass center,
	// the base velocity is hard to follow the leaflog algorithem

	for (i = 0; i < 6; i++) {
		// guess base cluster velocity
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + alpha[i] * mdpar.dt;
		mdpar.v0.v[i] = mdpar.v0_mh.v[i] + alpha[i] * hdt;
		d.v[i] = (mdpar.v0.v[i] + alpha[i] * hdt) * mdpar.dt;
		
		// guess base cluster velocity for next step
		mdpar.v0.v[i] = (mdpar.v0_mh.v[i] + 1.5 * alpha[i] * mdpar.dt) * e_eps; // rescale the velocity for next step
		V[i] = mdpar.v0.v[i];
		mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i] * e_eps; //for next step, and rescale the speed
	}
	V2wv(drot, dr, d)
	// rescale the cordinate for next step of cell volumn, and calculate the translation again
	dr.v[0] = e_eps * (sm.r.v[0] + dr.v[0]) - sm.r.v[0];
	dr.v[1] = e_eps * (sm.r.v[1] + dr.v[1]) - sm.r.v[1];
	dr.v[2] = e_eps * (sm.r.v[2] + dr.v[2]) - sm.r.v[2];

	sm.shiftMol(dr);
	rotate_SM(sm, drot, true);
}







void NPT_Speed_Verlet(LFVERLET_PM_MD &mdpar, PMOLECULE &pm) {
	// use the velocity of the base cluster
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	double hdt = mdpar.dt * 0.5;

	double *alpha = pm.alpha.v;
	double *V = pm.v.v;

	for (i = 0; i < 3; i++) {
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + alpha[i] * mdpar.dt;
		mdpar.v0.v[i] = mdpar.v0_mh.v[i] + alpha[i] * hdt;
		V[i] = mdpar.v0.v[i];
	}
}

void NPT_Speed_Kinetic(LFVERLET_PM_MD &mdpar, PMOLECULE &pm) {
	// use the velocity of base cluster
	// suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	double hdt = mdpar.dt * 0.5;

	double *alpha = pm.alpha.v;
	double *V = pm.v.v;

	double dv = 0;

	for (i = 0; i < 3; i++) {
		dv = alpha[i] * hdt;
		mdpar.v0_mh.v[i] = V[i] - dv;
		mdpar.v0_ph.v[i] = V[i] + dv;
		mdpar.v0.v[i] = V[i];
	}
}

void NPT_Verlet(LFVERLET_PM_MD &mdpar, PMOLECULE &pm, double e_eps) {
	// use the velocity of base cluster
	// suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	double hdt = mdpar.dt * 0.5;
	double delta = 0;
	VECTOR3 d;
	double *alpha = pm.alpha.v;
	double *V = pm.v.v;
	double *R = pm.r.v;

	// move base cluster and guess velocity for next step
	// because the base velocity is re-configured with the spacial moment at mass center,
	// the base velocity is hard to follow the leaflog algorithem

	for (i = 0; i < 3; i++) {
		// guess base cluster velocity
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + alpha[i] * mdpar.dt;
		mdpar.v0.v[i] = mdpar.v0_mh.v[i] + alpha[i] * hdt;
		d.v[i] = (mdpar.v0.v[i] + alpha[i] * hdt) * mdpar.dt;
		
		// guess base cluster velocity for next step
		mdpar.v0.v[i] = (mdpar.v0_mh.v[i] + 1.5 * alpha[i] * mdpar.dt) * e_eps; // rescale the velocity for next step
		V[i] = mdpar.v0.v[i];
		mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i] * e_eps; //for next step, and rescale it
	}

	// rescale the cordinate for next step, and calculate the translation again
	d.v[0] = e_eps * (pm.r.v[0] + d.v[0]) - pm.r.v[0];
	d.v[1] = e_eps * (pm.r.v[1] + d.v[1]) - pm.r.v[1];
	d.v[2] = e_eps * (pm.r.v[2] + d.v[2]) - pm.r.v[2];

	pm.shiftMol(d);
}

void calTransSpatialMomentum(MMOLECULE &mm, VECTOR3 &P) {
	CLUSTER *pc = NULL;
	VECTOR3 *v, *w, *rc;
	VECTOR<6> *V;
	SVECTOR<double, 3> dv;
	V3zero(P)
	int ic;
	for (ic = 0; ic < mm.nCluster; ic++) {
		pc = mm.cluster + ic;
		v = &(pc->mcDyn->v); w = &(pc->mcDyn->w); rc = &(pc->cm0); V = &(pc->bkm->V0);
		//translate the velocity to the mass center
		V2wv((*w), (*v), (*V))
		V3PV3((*w), (*rc), dv)
		V3plusV3((*v), dv, (*v))

		// mass center
		P.v[0] += pc->M * v->v[0];
		P.v[1] += pc->M * v->v[1];
		P.v[2] += pc->M * v->v[2];
	}
}

// **************************************************************************************************************************
// VERY IMPORTANT: in NPT calculation using Anderson's method, rescaling the translation coordinates with volumn change, 
// the thermokinetic energy contributing to pressure is only V = exp(eps) * V0, while not dr/dt = V + v_eps * r (dr/dt is
// the general speed with the cell volumn change, while NOT the speed within specific volumn). V is same as the speed in NVT.
// **************************************************************************************************************************

void NPT_calKineticEnergy_Cluster(MMOLECULE &mm, bool bTransOnly = false) {
	double Ekt = 0, Ekr = 0, Ektx, Ekty, Ektz;
	int ic = 0, iatom;
	CLUSTER *pc = NULL;
	VECTOR3 *v, *w;
	SVECTOR<double, 3> Iw;
	MATRIX<3> *I = NULL;
	mm.Ek = 0; mm.Ekt = 0; mm.Ektx = 0; mm.Ekty = 0; mm.Ektz = 0;
	double uEk2kT = ::unit_Ek_kT * 0.5;
	for (ic = 0; ic < mm.nCluster; ic++) {
		pc = mm.cluster + ic; v = &(pc->mcDyn->v); w = &(pc->mcDyn->w);
		Ektx = pc->M * v->v[0] * v->v[0] * uEk2kT;
		Ekty = pc->M * v->v[1] * v->v[1] * uEk2kT;
		Ektz = pc->M * v->v[2] * v->v[2] * uEk2kT;
		Ekt = Ektx + Ekty + Ektz;
		if (!bTransOnly) {
			I = &(pc->mcDyn->I);
			MpV3((*I), (*w), Iw)
			Ekr = Iw.v[0] * w->v[0] + Iw.v[1] * w->v[1] + Iw.v[2] * w->v[2];
			mm.Ek += Ekt + Ekr * uEk2kT;
		}
		mm.Ekt += Ekt;
		mm.Ektx += Ektx;
		mm.Ekty += Ekty;
		mm.Ektz += Ektz;
	}
}

void NPT_calKineticEnergy_Mol(MMOLECULE &mm, bool bTransOnly = false) {
	// calculate the translation energy at mass center of macromolecule
	double Ekt = 0, Ekr = 0, Ektx, Ekty, Ektz;
	int ic = 0, iatom;
	CLUSTER *pc = NULL;
	VECTOR3 *v, *w;
	SVECTOR<double, 3> Iw;
	VECTOR3 P;
	MATRIX<3> *I = NULL;
	mm.Ek = 0; mm.Ekt = 0;
	double uEk2kT = ::unit_Ek_kT * 0.5;
	for (ic = 0; ic < mm.nCluster; ic++) {
		pc = mm.cluster + ic; v = &(pc->mcDyn->v); w = &(pc->mcDyn->w);
		P.v[0] += pc->M * v->v[0];
		P.v[1] += pc->M * v->v[1];
		P.v[2] += pc->M * v->v[2];

		Ektx = pc->M * v->v[0] * v->v[0] * uEk2kT;
		Ekty = pc->M * v->v[1] * v->v[1] * uEk2kT;
		Ektz = pc->M * v->v[2] * v->v[2] * uEk2kT;
		Ekt = Ektx + Ekty + Ektz;
		if (!bTransOnly) {
			Ekt = pc->M * (v->v[0] * v->v[0] + v->v[1] * v->v[1] + v->v[2] * v->v[2]);
			I = &(pc->mcDyn->I);
			MpV3((*I), (*w), Iw)
			Ekr = Iw.v[0] * w->v[0] + Iw.v[1] * w->v[1] + Iw.v[2] * w->v[2];
			mm.Ek += (Ekt + Ekr) * uEk2kT;
		}
	}
	mm.Ektx = P.v[0] * P.v[0] / mm.mDynKin.M * uEk2kT;
	mm.Ekty = P.v[1] * P.v[1] / mm.mDynKin.M * uEk2kT;
	mm.Ektz = P.v[2] * P.v[2] / mm.mDynKin.M * uEk2kT;
	mm.Ekt = mm.Ektx + mm.Ekty + mm.Ektz;
}

void NPT_calKineticEnergy(SMOLECULE &sm, bool bTransOnly = false) {
	double Ekr = 0;
	VECTOR<6> *v6 = NULL;
	SVECTOR<double, 3> v, w, dv, Iw;
	MATRIX<3> *I = &(sm.I);
	double uEk2kT = ::unit_Ek_kT * 0.5;

	v6 = &(sm.c->bkm->V0);
	V2wv(w, v, (*v6))
		
	sm.Ektx = sm.c->M * v.v[0] * v.v[0] * uEk2kT;
	sm.Ekty = sm.c->M * v.v[1] * v.v[1] * uEk2kT;
	sm.Ektz = sm.c->M * v.v[2] * v.v[2] * uEk2kT;
	sm.Ekt = sm.Ektx + sm.Ekty + sm.Ektz;
	if (!bTransOnly) {
		MpV3((*I), w, Iw)
		Ekr = Iw.v[0] * w.v[0] + Iw.v[1] * w.v[1] + Iw.v[2] * w.v[2];
		sm.Ek = sm.Ekt + Ekr * uEk2kT;
	}
}

void NPT_calKineticEnergy(PMOLECULE &pm) {
	double uEk2kT = ::unit_Ek_kT * 0.5;
	double* v = pm.v.v;
		
	pm.Ektx = pm.c->M * v[0] * v[0] * uEk2kT;
	pm.Ekty = pm.c->M * v[1] * v[1] * uEk2kT;
	pm.Ektz = pm.c->M * v[2] * v[2] * uEk2kT;
	pm.Ekt = pm.Ektx + pm.Ekty + pm.Ektz;
}

bool NPT_MM_SpeedVerlet_Simplified(bool speed_verlet, MMOLECULE &mm, LFVERLET_MM_MD& mdpar, LFVERLET_ISOTROPIC_P_MD &pmdpar, VECTOR<6> *vect1, VECTOR<6> *vect2, double &vmax_diff, double &wmax_diff) {
	int nc = 0;
	CLUSTER *pc = NULL;

	VECTOR<6> *v_last = vect1, *v_new = vect2;

	// speed is varying to get consisten speed and accleration based on leap-frog velert MD 
	// calculate Coriolis acceleration term and gyroscopic spatial force
	CalMol_ab(mm);

	for (nc = 0; nc < mm.nCluster; nc++) {V62V6(mm.cluster[nc].bkm->V0, v_last[nc])}
	
	NEIMO(mm, true);  // it is always true even if velocity changed, total force is to be calculated in NEIMO(...)

	if (speed_verlet) NPT_Speed_Verlet(mdpar, mm); // get new speed at current position and +1/2 time step, and reset the speed in mm to the new speed
	else NPT_Speed_Kinetic(mdpar, mm);
	
	for (nc = 0; nc < mm.nCluster; nc++) {V62V6(mm.cluster[nc].bkm->V0, v_new[nc])}
	max_diff(v_new, v_last, mm.nCluster, vmax_diff, wmax_diff);

	//calKineticEnergy(mm); // kinetic energy 
	//calFrameKineticEnergy(mm);
	//mm.Ek_internal = mm.Ek - mm.Ek_frame;

	return true;
}

extern double max_diff(VECTOR3 &v1, VECTOR3 &v2);
bool NPT_SM_SpeedVerlet_Simplified(bool speed_verlet, SMOLECULE &m, LFVERLET_SM_MD& mdpar, LFVERLET_ISOTROPIC_P_MD &pmdpar, double& vmax_diff, double& wmax_diff) {
	VECTOR3 v_new, v_old, w_new, w_old;

	V2wv(w_old, v_old, m.c->bkm->V0)

	m.calTotalForce();
	calAccel_f(m);

	if (speed_verlet) NPT_Speed_Verlet(mdpar, m); // get new speed at current position and +1/2 time step
	else NPT_Speed_Kinetic(mdpar, m);

	V2wv(w_new, v_new, m.c->bkm->V0);

	vmax_diff = max_diff(v_new, v_old);
	wmax_diff = max_diff(w_new, w_old);

	//NPT_calKineticEnergy(m); // kinetic energy 

	return true;
}

bool NPT_PM_SpeedVerlet_Simplified(bool speed_verlet, PMOLECULE &m, LFVERLET_PM_MD& mdpar, LFVERLET_ISOTROPIC_P_MD &pmdpar, double& vmax_diff) {
	VECTOR3 v_new, v_old;
	
	V32V3(m.v, v_old)

	m.calTotalForce();
	calAccel_f(m);

	if (speed_verlet) NPT_Speed_Verlet(mdpar, m); // get new speed at current position and +1/2 time step
	else NPT_Speed_Kinetic(mdpar, m);

	V32V3(m.v, v_new)
	vmax_diff = max_diff(v_new, v_old);

	//NPT_calKineticEnergy(m); // kinetic energy 

	return true;
}

// **********************************************************************************************************************
// In NPT, V = dr/dt - veps * r, is the speed within specific cell-volumn, like that in NVT. 
// For rigid cluster, virtual force from a couple with volumn change, is applied to the mass center only of cluster, for translation.
// In effect, a virtual force F = ... is applied to the mass center of cluster
// Rotation of the cluster (at mass center) is released from the coupling with barostat
// and has no contribution to the internal pressure
// **********************************************************************************************************************
void cal_npt_force(MMOLECULE *mm, double b, double ksi, double ksi_mc, VECTOR3 &vmc) { // b is the coefficient for speed
	int method = ::MM_SpeedVerlet_method;

	int ic;
	CLUSTER *pc;
	SVECTOR<double, 3> force, torq, dr;
	SVECTOR<double, 6> f, f1;
	MATRIX<6> *phai;
	float mass;
	double *fh; // Hoover force
	double *v, *v1, *v2;

	if (method == 0) {
		V6zero(mm->mDynKin.f_Hoover)
	}
	NEIMO_CalHooverForce(*mm, ksi, mm->mDynKin.f_Hoover);
	for (ic = 0; ic < mm->nCluster; ic++) {
		pc = mm->cluster + ic;
		mass = pc->M; 

		v = pc->mcDyn->v.v;
		force.v[0] = mass * (b * v[0] - ksi_mc * vmc.v[0]);
		force.v[1] = mass * (b * v[1] - ksi_mc * vmc.v[1]);
		force.v[2] = mass * (b * v[2] - ksi_mc * vmc.v[2]);
		VECT3((*(pc->vp)), pc->cm, dr)
		V3PV3(dr, force, torq);
		
		fh = pc->dyn->f_Hoover.v;
		fh[0] += torq.v[0];
		fh[1] += torq.v[1];
		fh[2] += torq.v[2];
		fh[3] += force.v[0];
		fh[4] += force.v[1];
		fh[5] += force.v[2];

		if (method == 0) {
			phai = &(pc->invNE->phai_cm2me);
			wv2V(torq, force, f)
			MpV6((*phai), f, f1) // f1 = phai(cm=>pc) * f_Hoover

			v1 = mm->mDynKin.f_Hoover.v; v2 = f1.v; V6plus(v1, v2, v1) // mDynKin.f += f1
		}
	}
}

void cal_npt_force(SMOLECULE *sm, double b, double ksi, double ksi_mc, VECTOR3 &vmc) { // b is the coefficient for speed
	SVECTOR<double, 3> force, torque;
	SVECTOR<double, 3> w, v;
	V2wv(w, v, sm->c->bkm->V0)
	MATRIX<3> *I = &(sm->I);
	float mass = sm->c->M;
	double *v1, *v2;
	double *fh = sm->c->dyn->f_Hoover.v;

	memset(fh, 0, SIZE_V6);

	force.v[0] = mass * (b * v.v[0] - ksi * v.v[0] - ksi_mc * vmc.v[0]);
	force.v[1] = mass * (b * v.v[1] - ksi * v.v[1] - ksi_mc * vmc.v[1]);
	force.v[2] = mass * (b * v.v[2] - ksi * v.v[2] - ksi_mc * vmc.v[2]);

	torque.v[0] = -ksi * sm->c->bkm->P0.v[0];
	torque.v[1] = -ksi * sm->c->bkm->P0.v[1];
	torque.v[2] = -ksi * sm->c->bkm->P0.v[2];

	// total virtual force
	fh[0] = torque.v[0];
	fh[1] = torque.v[1];
	fh[2] = torque.v[2];
	fh[3] = force.v[0];
	fh[4] = force.v[1];
	fh[5] = force.v[2];
}

void cal_npt_force(PMOLECULE *pm, double b, double ksi, double ksi_mc, VECTOR3 &vmc) { // b is the coefficient for speed
	double *v = pm->v.v;
	double *v1, *v2;
	double *fh = pm->c->dyn->f_Hoover.v;
	float mass = pm->c->M;

	memset(fh, 0, SIZE_V6);

	// total virtual force
	fh[3] = mass *(v[0] * b - ksi * v[0] - ksi_mc * vmc.v[0]);
	fh[4] = mass *(v[1] * b - ksi * v[1] - ksi_mc * vmc.v[1]);
	fh[5] = mass *(v[2] * b - ksi * v[2] - ksi_mc * vmc.v[2]);
}

void MDCELL_NPT_SpeedVerlet(MDCELL_ThreadJob< ISOTROPIC_NPT_CELL> *job, bool *bSpeedVerlet) {
	int nm = 0, imol = 0, ihdyn = 0;
	char msg[256] = "\0";

	ISOTROPIC_NPT_CELL *npt_cell = job->job->mdcell;
	MMOL_MD_CELL *mdcell = npt_cell->mdcell;

	double dt =  npt_cell->pmdpar.dt;
	double uConvg = 0.005 / dt, wConvg = 0.005 / dt; //0.017 is 1 arc degree;
	if (uConvg < 2e-4) uConvg = 2e-4;
	if (wConvg < 2e-4) wConvg = 2e-4;

	double b, ksi;
	double veps = npt_cell->pdyn.v, alpha_eps = npt_cell->pdyn.alpha;

	double ksi_mc = 0;
	if (mdcell->bHDynMC_T) ksi_mc = mdcell->mcHDyn.nhc->ksi.v[0];
	
	VECTOR3* Vmc = &(mdcell->Vmc);
	VECTOR3 P;

	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mdcell->lfmd.m[imol].fetch_ksi();

		// calculate the NPT force on each cluster
		calTransSpatialMomentum(mdcell->mm.m[imol], P);
		ksi = mdcell->lfmd.m[imol].ksi; // actually the ksi_internal is same as ksi_frame in NPT, we do not support independent ksi_frame
		b = (2 + 1.0 / npt_cell->NF_pressure) * veps;
		cal_npt_force(mdcell->mm.m + imol, -b, ksi, ksi_mc, *Vmc);

		if (::bBrownian) Brownian_force(mdcell->mm.m + imol, 0, mdcell->mm.m[imol].nCluster - 1, mdcell->mm.m[imol].mDynKin.f_Brown);

		mdcell->mvb.mm_status.m[imol] = NPT_MM_SpeedVerlet_Simplified(*bSpeedVerlet, mdcell->mm.m[imol], mdcell->lfmd.m[imol], npt_cell->pmdpar,
			mdcell->mvb.mmv1.m[imol].cv.m, mdcell->mvb.mmv2.m[imol].cv.m, 
			mdcell->mvb.mm_dvmax.m[imol], mdcell->mvb.mm_dwmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.mm_status.m[imol]) {
			sprintf(msg, "Macromol # %d [%s] : fails for speed-verlet", imol, mdcell->mm.m[imol].mol);
			show_log(msg, true); //show_infor(msg_show);
		}
		if (mdcell->mvb.mm_dvmax.m[imol] <= uConvg && mdcell->mvb.mm_dwmax.m[imol] <= wConvg) mdcell->mvb.mm_sc_status.m[imol] = true;
		else mdcell->mvb.mm_sc_status.m[imol] = false;
	}

	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		if (mdcell->lfsmd.m[imol].hdyn->hdyn_status) continue;  // convergent already

		// calculate the NPT force on each cluster
		ksi = mdcell->lfsmd.m[imol].hdyn->nhc->ksi.v[0]; // actually the ksi_internal is same as ksi_frame in NPT, we do not support independent ksi_frame
		b = (2 + 1.0 / npt_cell->NF_pressure) * veps;
		cal_npt_force(mdcell->sm.m + imol, -b, ksi, ksi_mc, *Vmc);

		if (::bBrownian) SM_Brownian_force(mdcell->sm.m + imol);

		mdcell->mvb.sm_status.m[imol] = NPT_SM_SpeedVerlet_Simplified(*bSpeedVerlet, mdcell->sm.m[imol], mdcell->lfsmd.m[imol], npt_cell->pmdpar,
			mdcell->mvb.sm_dvmax.m[imol], mdcell->mvb.sm_dwmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.sm_status.m[imol]) {
			sprintf(msg, "mol # %d [%s] : fails for speed-verlet", imol, mdcell->sm.m[imol].mol);
			show_log(msg, true); //show_infor(msg_show);
		}
		if (mdcell->mvb.sm_dvmax.m[imol] <= uConvg && mdcell->mvb.sm_dwmax.m[imol] < wConvg) mdcell->mvb.sm_sc_status.m[imol] = true;
		else mdcell->mvb.sm_sc_status.m[imol] = false;
	}

	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		if (mdcell->lfpmd.m[imol].hdyn->hdyn_status) continue;  // convergent already

		// calculate the NPT force on each cluster
		ksi = mdcell->lfpmd.m[imol].hdyn->nhc->ksi.v[0]; // actually the ksi_internal is same as ksi_frame in NPT, we do not support independent ksi_frame
		b = (2 + 1.0 / npt_cell->NF_pressure) * veps;
		cal_npt_force(mdcell->pm.m + imol, -b, ksi, ksi_mc, *Vmc);

		if (::bBrownian) PM_Brownian_force(mdcell->pm.m + imol);

		mdcell->mvb.pm_status.m[imol] = NPT_PM_SpeedVerlet_Simplified(*bSpeedVerlet, mdcell->pm.m[imol], mdcell->lfpmd.m[imol], npt_cell->pmdpar,
			mdcell->mvb.pm_dvmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.pm_status.m[imol]) {
			sprintf(msg, "mol # %d [%s] : fails for speed-verlet", imol, mdcell->pm.m[imol].mol);
			show_log(msg, true); //show_infor(msg_show);
		}
		if (mdcell->mvb.pm_dvmax.m[imol] <= uConvg) mdcell->mvb.pm_sc_status.m[imol] = true;
		else mdcell->mvb.pm_sc_status.m[imol] = false;
	}
}

void MDCELL_NPT_Verlet(MDCELL_ThreadJob< ISOTROPIC_NPT_CELL > *job) {
	int nm = 0, imol = 0, ihdyn = 0;

	ISOTROPIC_NPT_CELL *npt_cell = job->job->mdcell;
	MMOL_MD_CELL *mdcell = npt_cell->mdcell;

	double e_eps = job->job->mdcell->pmdpar.e_eps;

	//char msg_show[1024] = "\0", msg_log[1024] = "\0";

	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		NPT_Verlet(mdcell->lfmd.m[imol], mdcell->mm.m[imol], e_eps); // Verlet algorithm to calculate the ta and tp for next timestep
	}

	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		NPT_Verlet(mdcell->lfsmd.m[imol], mdcell->sm.m[imol], e_eps); // Verlet algorithm to calculate the ta and tp for next timestep
	}

	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		NPT_Verlet(mdcell->lfpmd.m[imol], mdcell->pm.m[imol], e_eps); // Verlet algorithm to calculate the ta and tp for next timestep
	}
}

void NPT_MDCELL_CalMolKineticEnergy(MDCELL_ThreadJob<ISOTROPIC_NPT_CELL> *job, double *Ekt) {
	MMOL_MD_CELL *mdcell = job->job->mdcell->mdcell;
	int nm = 0, imol = 0;
	//double veps = job->job->mdcell->pdyn.v;
	double Ek = 0;
	VECTOR3 P;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		calTransSpatialMomentum(*mm, P); // velocity of each cluster at its mass-center
		if (::method_Virial == CLUSTER_VIRIAL) {
			NPT_calKineticEnergy_Cluster(*mm); // translating kinetic energy of all clusters
		}
		else if (::method_Virial == MOLECULE_VIRIAL) {
			NPT_calKineticEnergy_Mol(*mm); // translating kinetic energy of mass center of the molecule
		}
		Ek += mm->Ekt;
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		NPT_calKineticEnergy(*sm);
		Ek += sm->Ekt;
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		NPT_calKineticEnergy(*pm);
		Ek += pm->Ekt;
	}
	*Ekt = Ek;
}

void MDCELL_MMClusterVelocity(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell; 
	int nm = 0, imol = 0;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		calClusterVelocity0(*mm);
	}
}

/*
void NPT_MDCELL_SpatialMomentum(MDCELL_ThreadJob< ISOTROPIC_NPT_CELL> *job, SVECTOR<double, 6> *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	double veps = job->job->mdcell->pdyn.v;
	VECTOR3 P, v, p;
	VECTOR3 Pw, w, pw, pw1;
	VECTOR3 Rcm;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		calTransSpatialMomentum(*mm, P); // velocity of each cluster at its mass-center
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			Rcm.v[0] = pc->cm.v[0] + pc->dr_cmm.v[0];
			Rcm.v[1] = pc->cm.v[1] + pc->dr_cmm.v[1];
			Rcm.v[2] = pc->cm.v[2] + pc->dr_cmm.v[2];
			p.v[0] = pc->M * (pc->mcDyn->v.v[0] - veps * Rcm.v[0]);
			p.v[1] = pc->M * (pc->mcDyn->v.v[1] - veps * Rcm.v[1]);
			p.v[2] = pc->M * (pc->mcDyn->v.v[2] - veps * Rcm.v[2]);

			V3plusV3(p, P, P)

			// for rotation speed, Rcm x v = Rcm x (v - veps * Rcm)
			// inertia momentum from the translation of particle
			V3PV3(Rcm, p, pw)
			// inertia momentum from the rotation of particle
			MpV3(pc->mcDyn->I, pc->mcDyn->w, pw1)
			Pw.v[0] += pw.v[0] + pw1.v[0];
			Pw.v[1] += pw.v[1] + pw1.v[1];
			Pw.v[2] += pw.v[2] + pw1.v[2];
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		V32V3(sm->r, Rcm)
		//V32V3(sm->Rg, Rcm)
		p.v[0] = sm->m * (sm->c->bkm->V0.v[3] - veps * Rcm.v[0]);
		p.v[1] = sm->m * (sm->c->bkm->V0.v[4] - veps * Rcm.v[1]);
		p.v[2] = sm->m * (sm->c->bkm->V0.v[5] - veps * Rcm.v[2]);

		V3plusV3(p, P, P)

		// for rotation speed, Rcm x v = Rcm x (v - veps * Rcm)
		// inertia momentum from the translation of particle
		V3PV3(Rcm, p, pw)
		// inertia momentum from the rotation of particle
		w.v[0] = sm->c->bkm->V0.v[0];
		w.v[1] = sm->c->bkm->V0.v[1];
		w.v[2] = sm->c->bkm->V0.v[2];
		MpV3(sm->I, w, pw1)
		Pw.v[0] += pw.v[0] + pw1.v[0];
		Pw.v[1] += pw.v[1] + pw1.v[1];
		Pw.v[2] += pw.v[2] + pw1.v[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		V32V3(pm->r, Rcm)
		//V32V3(pm->Rg, Rcm)
		p.v[0] = pm->m * (pm->v.v[0] - veps * Rcm.v[0]);
		p.v[1] = pm->m * (pm->v.v[1] - veps * Rcm.v[1]);  
		p.v[2] = pm->m * (pm->v.v[2] - veps * Rcm.v[2]);

		V3plusV3(p, P, P)

		// for rotation speed, Rcm x v = Rcm x (v - veps * Rcm)
		// inertia momentum from the translation of particle
		V3PV3(Rcm, p, pw)
		// PM has not rotation
		Pw.v[0] += pw.v[0];
		Pw.v[1] += pw.v[1];
		Pw.v[2] += pw.v[2];
	}
	mp->v[0] = Pw.v[0]; mp->v[1] = Pw.v[1]; mp->v[2] = Pw.v[2];
	mp->v[3] = P.v[0]; mp->v[4] = P.v[1]; mp->v[5] = P.v[2];
}

void NPT_MDCELL_ZeroSpatialMomentum(MDCELL_ThreadJob< ISOTROPIC_NPT_CELL > *job, SVECTOR<double, 6> *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	double veps = job->job->mdcell->pdyn.v;

	VECTOR3 v, w, Vw, P;
	V2wv(w, v, (*mp))

	if (FABS(v.v[0]) < 1e-8 && FABS(v.v[1]) < 1e-8 && FABS(v.v[2]) < 1e-8) return;
	else show_log("zero the total spatial momentum of the cell", true);

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		mm->V0.v[3] -= v.v[0]; mm->V0.v[4] -= v.v[0]; mm->V0.v[5] -= v.v[0];
		// Update all the requirement change relating to the velocity change
		calClusterVelocity0(*mm);
		calSpatialMomentum0(*mm, true);
		mm->CalMassCenterVelocity_MassCenterMomentum();
		calTransSpatialMomentum(*mm, P); // velocity of each cluster at its mass-center
		NPT_calKineticEnergy(*mm, veps);

		// the rotation of whole cell on the macromolecule is not corrected !!!
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		// compensate the rotation of the whole cell with the translation of the particle
		V3PV3(w, sm->r, Vw)
		sm->c->bkm->V0.v[3] -= v.v[0] + Vw.v[0];
		sm->c->bkm->V0.v[4] -= v.v[1] + Vw.v[1];
		sm->c->bkm->V0.v[5] -= v.v[2] + Vw.v[2];
		NPT_calKineticEnergy(*sm, veps);
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		// compensate the rotation of the whole cell with the translation of the particle
		V3PV3(w, pm->r, Vw)
		pm->v.v[0] -= v.v[0] + Vw.v[0]; 
		pm->v.v[1] -= v.v[1] + Vw.v[1]; 
		pm->v.v[2] -= v.v[2] + Vw.v[2];
		NPT_calKineticEnergy(*pm, veps);
	}
}

void NPT_ZeroSpatialMomentum(MDCELL_ThreadJob< ISOTROPIC_NPT_CELL> *job, SVECTOR<double, RF_NVAR> *mp, int nJob) {
	MTOperate1<MDCELL_ThreadJob< ISOTROPIC_NPT_CELL>, SVECTOR<double, 6> >((void*)(&NPT_MDCELL_SpatialMomentum), job, nJob, mp);
	VECTOR3 P, v, Pw, w;
	float M = job->job->mdcell->mdcell->M;
	MATRIX<3> *invI = &(job->job->mdcell->mdcell->invI);
	int i;
	for (i = 0; i < nJob; i++) {
		Pw.v[0] += mp[i].v[0]; Pw.v[1] += mp[i].v[1]; Pw.v[2] += mp[i].v[2];
		P.v[0] += mp[i].v[3]; P.v[1] += mp[i].v[4]; P.v[2] += mp[i].v[5];
	}
	MpV3((*invI), Pw, w)
	v.v[0] = P.v[0] / M; v.v[1] = P.v[1] / M; v.v[2] = P.v[2] / M;
	
	for (i = 0; i < nJob; i++) {
		mp[i].v[0] = w.v[0]; mp[i].v[1] = w.v[1]; mp[i].v[2] = w.v[2];
		mp[i].v[3] = v.v[0]; mp[i].v[4] = v.v[1]; mp[i].v[5] = v.v[2];
	}
	MTOperate1<MDCELL_ThreadJob< ISOTROPIC_NPT_CELL>, SVECTOR<double, 6> >((void*)(&NPT_MDCELL_ZeroSpatialMomentum), job, nJob, mp);
}
*/
bool NPT_LeapFrogVerlet_SelfConsistent_SpeedAccel(ISOTROPIC_NPT_CELL *npt_cell, MDCELL_ThreadJob<MMOL_MD_CELL> *job_mdcell, MDCELL_ThreadJob<ISOTROPIC_NPT_CELL> *job_nptcell, MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> *hdyn_thread_var, HingeConstraintVar *hinge_constraint_var, int nThreads, bool bSpeedVerlet) {
	MMOL_MD_CELL *mdcell = npt_cell->mdcell;
	bool bDynConvergent = false;
	int iloop = 0, max_loops = 8, i, imol;
	ARRAY<int> *indx;

	double RF0 = npt_cell->RF;
	SVECTOR<double, 3> virial0;
	virial0.v[0] = npt_cell->virial.v[0]; virial0.v[1] = npt_cell->virial.v[1]; virial0.v[2] = npt_cell->virial.v[2];
	SVECTOR<double, 4> virial;

	ARRAY<InteractRes> interactRes; interactRes.SetArray(nThreads);

	float Ek = 0;
	VECTOR3 w, Pw;

	ARRAY<bool> status; status.SetArray(nThreads);
	ARRAY<double> av; av.SetArray(nThreads);
	ARRAY<SVECTOR<double, 3> > sv3; sv3.SetArray(nThreads);
	ARRAY<SVECTOR<double, 4> > sv4; sv4.SetArray(nThreads);

	extern void MDCELL_SpatialMomentum(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 3> *P);

	MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_MMClusterVelocity), job_mdcell, nThreads);

	for (i = 0; i < mdcell->hdyn.n; i++) mdcell->hdyn.m[i].hdyn_status = false;
	while (!bDynConvergent && iloop < max_loops) {
		 // new ksi of barostat, dependent on the kinetic energy of barostat
		P_NHC(npt_cell->pdyn, npt_cell->phdyn);
		// calculate the kinetic energy of each molecule, with the old veps, also the cluster's translation & rotation speed at mass center
		MTOperate1< MDCELL_ThreadJob<ISOTROPIC_NPT_CELL>, double>((void*)(&NPT_MDCELL_CalMolKineticEnergy), job_nptcell, nThreads, av.m); 
		npt_cell->Ekt = 0;
		for (i = 0; i < nThreads; i++) npt_cell->Ekt += av.m[i];
		// run NHC to get a new ksi for molecules
		MDCELL_HDyn_MultiThreadOperate<MMOL_MD_CELL>((void*)(&MDCELL_HDyn), hdyn_thread_var, nThreads);

		// run NHC for the mass-center velocity of the whole cell
		if (mdcell->bHDynMC_T) {
			//MTOperate1< MDCELL_ThreadJob< ISOTROPIC_NPT_CELL>, SVECTOR<double, 6> >((void*)(&NPT_MDCELL_SpatialMomentum), job_nptcell, nThreads, sv6.m);
			MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 3> >((void*)(&MDCELL_SpatialMomentum), job_mdcell, nThreads, sv3.m);
			V3zero(mdcell->Vmc) V3zero(Pw)
			for (i = 0; i < nThreads; i++) {
				Pw.v[0] += sv3.m[i].v[0];
				Pw.v[1] += sv3.m[i].v[1];
				Pw.v[2] += sv3.m[i].v[2];
			}
			mdcell->Vmc.v[0] = Pw.v[0] / mdcell->M; mdcell->Vmc.v[1] = Pw.v[1] / mdcell->M; mdcell->Vmc.v[2] = Pw.v[2] / mdcell->M;

			if (mdcell->bHDynMC_T) {
				V3ABS2(mdcell->Vmc, Ek) Ek *= mdcell->M * unit_Ek_kT;
				mdcell->Ekt_mc = Ek * 0.5;
				mdcell->mcHDyn.nhc->set_dT(Ek); // supposedly the spatial momentum should be 0
				mdcell->mcHDyn.nhc->verlet_NHC();
			}
		}

		// run Speed-Verlet of molecules
		for (i = 0; i < nThreads; i++) status.m[i] = bSpeedVerlet;
		MTOperate1< MDCELL_ThreadJob<ISOTROPIC_NPT_CELL>, bool>((void*)(&MDCELL_NPT_SpeedVerlet), job_nptcell, nThreads, status.m); 

		// then run pressure's speed-verlet
		// important: this verlet will change veps, which is used to for speed-verlet of molecules, so run pressure's speed-verlet after molecules'
		//MDCELL_PressureForce(job_nptcell, av.m, nThreads); // calculate the internal pressure and force 
		// Ekt was calculated by NPT_MDCELL_CalMolKineticEnergy above

		if (::method_Virial == CLUSTER_VIRIAL) {
			// Virial correction is required for the force on hinges of MM, which is velocity/acceleration related
			virial.v[0] = 0; virial.v[1] = 0; virial.v[2] = 0; virial.v[3] = 0; 
			if (mdcell->bHingeConstraint && mdcell->mm.n > 0 && ::method_Virial == CLUSTER_VIRIAL) {
				MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, HingeConstraintVar, InteractRes >((void*)(&MDCELL_ClusterVirialCorrect_HingeConstraint), job_mdcell, nThreads, hinge_constraint_var, interactRes.m);
				for (i = 0; i < nThreads; i++) {
					virial.v[0] += interactRes.m[i].STxx;
					virial.v[1] += interactRes.m[i].STyy;
					virial.v[2] += interactRes.m[i].STzz;
					virial.v[3] += interactRes.m[i].STtr;
				}
			}
			//sprintf(errmsg, "virial from hinge force: %f", virial.v[3] * ::unit_Ek_kT); show_log(errmsg, true);

			npt_cell->RF = RF0 + virial.v[3];
			npt_cell->virial.v[0] = virial0.v[0] + virial.v[0];
			npt_cell->virial.v[1] = virial0.v[1] + virial.v[1];
			npt_cell->virial.v[2] = virial0.v[2] + virial.v[2];
		}
		else {}

		PressureForce(npt_cell); // calculate the internal pressure and force 
		cal_dyn(npt_cell->pdyn, npt_cell->phdyn.nhc->ksi.v[0]); //calculate the acceleration of eps
		if (bSpeedVerlet) Speed_Verlet(npt_cell->pmdpar, npt_cell->pdyn);
		else Speed_Kinetic(npt_cell->pmdpar, npt_cell->pdyn);

		// check whether the speed and accel are convergent
		MDCELL_HDyn_MultiThreadOperate<MMOL_MD_CELL>((void*)(&MDCELL_CheckHDynCovergent), hdyn_thread_var, nThreads);
		bDynConvergent = true;
		for (i = 0; i < mdcell->hdyn.n; i++) {
			if (!mdcell->hdyn.m[i].hdyn_status) {bDynConvergent = false; break;}
		}
		// we assume that eps is always convergent. Is this reasonable ?
		iloop += 1;
	}
	npt_cell->pmdpar.cal_e(); // calculation the rescale-ratio for next step, based on v_ph on the barostat
	return bDynConvergent;
}

void MDCellVolUpdate(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, double V) {
	float x = cmm->xd[0], y = cmm->xd[1], z = cmm->xd[2], v0 = x * y * z;
	float c = powf(V / v0, 1./3);
	x *= c; y *= c; z *= c;
	float hx = x * 0.5, hy = y * 0.5, hz = z * 0.5;
	mdcell->h[0] = hx; mdcell->h[1] = hy; mdcell->h[2] = hz;
	cmm->set_cell_range(cmm->periodic_cell, -hx, hx, -hy, hy, -hz, hz);
	cmm->set_cell_pos();
}

template <class VATOM> void spmeCellSizeUpdate(VSPME<VATOM> &vspme, float dx, float dy, float dz) {
	if (bPolarize) {
		vspme.xl[0] = -dx * 0.5; vspme.xl[1] = -dy * 0.5; vspme.xl[2] = -dz * 0.5;
		vspme.cell_dim(dx, dy, dz); // important, after init();
		vspme.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		vspme.init_C(); 
		//vspme.init_b(); vspme.init_vars();
	}
}

void NPT_MD_PROC(ISOTROPIC_NPT_CELL &npt_cell) {
	MMOL_MD_CELL *mdcell = npt_cell.mdcell;
	CMM_CELL3D<BASIC_CLUSTER> *cmm = npt_cell.cmm;

	if (!cmm->periodic_cell) {
		show_msg("NPT has to work on a periodical cell"); return;
	}

	int i;

	// we need to make sure the mass center of the cell is at zero
	VECTOR3 mass_center; // zero
	for (i = 0; i < 5; i++) {
		set_free_cell_mass_center(*mdcell, mass_center);
		molecule_check_periodic_cell(*mdcell); // check the molecule mass center is inside the cell
		// after the periodical check, the mass center could be changed again. So, do it iteratively
	}

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************
	show_all_molecules();

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

	init_BSplineFunc<2>(::bsp, ::indx_bspline);
	VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	vspme_q.bST = true; vspme_q.bST_diagonal = true;
	VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	vspme_mu.bST = true; vspme_mu.bST_diagonal = true;
	if (bPolarize) {
		mdcell->init_polarization_buff();
		mdcell->init_dipole_hist();
		mdcell->set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		vspme_mu.bsp = &(::bsp);
		init_spme(mdcell, vspme_mu); // init the viral atoms
		vspme_mu.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_mu.xl[0] = cmm->xl[0]; vspme_mu.xl[1] = cmm->xl[1]; vspme_mu.xl[2] = cmm->xl[2];
		vspme_mu.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		vspme_mu.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();
	}
	else { // charge only
		mdcell->set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		vspme_q.bsp = &(::bsp);
		init_spme(mdcell, vspme_q); // init the viral atoms
		vspme_q.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_q.xl[0] = cmm->xl[0]; vspme_q.xl[1] = cmm->xl[1]; vspme_q.xl[2] = cmm->xl[2];
		vspme_q.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		vspme_q.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
#endif

extern void MDCELL_mass_center(MMOL_MD_CELL &mmcell);

	bool bStorageChain = false;
	InitMD(*mdcell, cmm, bStorageChain);

	npt_cell.pdyn.Ek = 0; npt_cell.pdyn.v = 0; 
	InitVelocity(npt_cell.pdyn, ::Ek0Internal); 
	npt_cell.pdyn.alpha = 0; npt_cell.pdyn.fc = 0;
	npt_cell.pdyn.V0 = cmm->xd[0] * cmm->xd[1] * cmm->xd[2]; npt_cell.pdyn.V = npt_cell.pdyn.V0;
	npt_cell.pdyn.eps = 0; npt_cell.pdyn.Pint = 0; 
	//npt_cell.Pext = 1; // external pressure : 1 atmosphere


	long nloop = 0;
	int nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0, U_LJ = 0, U_estat = 0;
	
	InteractVar LJ_var; 
	LJ_var.bVirial = true; LJ_var.bST_diagonal = true;

	SPME_VAR spme_var;
	spme_var.bVirial = true; spme_var.bST_diagonal = true; spme_var.bEfieldOnly = false;
	spme_var.esv1.bEwald = true; spme_var.esv1.init_EwaldSum(cmm->xd[0] * cmm->xd[1] * cmm->xd[2], ::rcut_Ewd);
	spme_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;

	InteractRes LJ_res, spme_res;

	TorsionVar torsion_var[MAX_THREADS];
	HingeConstraintVar hinge_constraint_var[MAX_THREADS];
	InteractRes mires[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		torsion_var[i].bVirial = true; torsion_var[i].bST_diagonal = true;
		hinge_constraint_var[i].bVirial = true; hinge_constraint_var[i].bST_diagonal = true;
	}

	double Uext = 0; // external electric field

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	int nMDThreads = MAX_THREADS;
	MDCELL_Job<MMOL_MD_CELL> mdcellJob[MAX_THREADS];
	AssignJob_MDCELL<MMOL_MD_CELL>(mdcell, mdcellJob, nMDThreads);
	sprintf(errmsg, "Job # : %d; working with molecules:", MAX_THREADS); show_log(errmsg, true);
	for (i = 0; i < MAX_THREADS; i++) {
		sprintf(errmsg, "%d -- MM %d, SM %d, PM %d", i, mdcellJob[i].mmIndx.n, mdcellJob[i].smIndx.n, mdcellJob[i].pmIndx.n);
		show_log(errmsg, true);
		//for (nm = 0; nm < mdcellJob[i].mmIndx.n; nm++) {
		//	sprintf(errmsg, "%d ", mdcellJob[i].mmIndx.m[nm]); show_log(errmsg, false);
		//}
		//show_log("", true);
	}

	MDCELL_Job<ISOTROPIC_NPT_CELL> nptcellJob[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		nptcellJob[i].mmIndx.SetArray(mdcellJob[i].mmIndx.n);
		for (nc = 0; nc < mdcellJob[i].mmIndx.n; nc++) nptcellJob[i].mmIndx.m[nc] = mdcellJob[i].mmIndx.m[nc];
		nptcellJob[i].smIndx.SetArray(mdcellJob[i].smIndx.n);
		for (nc = 0; nc < mdcellJob[i].smIndx.n; nc++) nptcellJob[i].smIndx.m[nc] = mdcellJob[i].smIndx.m[nc];
		nptcellJob[i].pmIndx.SetArray(mdcellJob[i].pmIndx.n);
		for (nc = 0; nc < mdcellJob[i].pmIndx.n; nc++) nptcellJob[i].pmIndx.m[nc] = mdcellJob[i].pmIndx.m[nc];

		nptcellJob[i].mdcell = &npt_cell;
	}

	MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> mdcellHDynThreadVar[MAX_THREADS];
	AssignJob_MDCELL_HDYN<MMOL_MD_CELL>(mdcell, mdcellHDynThreadVar, nMDThreads, false);
	// if Barostat has an isolated NHC, the last HDYN in mdcell is independent to the barostat
	// else, barostat use the last NHC in HDYN, the default one

	// assuming each cluster has minimum dimension 2.5^3, and we expand the array dimension for each subcell 50% bigger than estimation
	float vcluster = 1.5 * 1.5 * 1.5;
	int ncs = int(cmm->fw[0] * cmm->fw[1] * cmm->fw[2] / vcluster * 1.2 + 0.5);
	if (ncs < 2) ncs = 2;
	CMM_array_set_storage<BASIC_CLUSTER>(cmm, ncs);
	CMM_CELL3D<BASIC_CLUSTER> job_cmm[MAX_THREADS];

	INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation

	MDCELL_ThreadJob<MMOL_MD_CELL> job_mdcell[MAX_THREADS];
	MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > job_mdcellCMM[MAX_THREADS];

	bool status_job[MAX_THREADS];
	double djob[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	MATRIX<3> mjob[MAX_THREADS];

	MDCELL_ThreadJob<ISOTROPIC_NPT_CELL> job_nptcell[MAX_THREADS];

	SVECTOR<double, 4> virial;

	for (i = 0; i < MAX_THREADS; i++) {
		if (MAX_THREADS > 1) {
			CMM_array_construct_for_storage<BASIC_CLUSTER>(job_cmm + i, cmm, ncs);
		}

		job_mdcell[i].set_pars(i, mdcellJob + i);
		job_mdcellCMM[i].set_pars(i, mdcellJob + i);
		job_mdcellCMM[i].cmm = cmm;
		job_nptcell[i].set_pars(i, nptcellJob + i);

		// using for explicit interaction calculation
		interact_job[i].var = LJ_var;
		interact_job[i].var.bVirial = true; interact_job[i].var.bST_diagonal = true;
		interact_job[i].evar.esv1.bEwald = false; // explicit direct interaction, non-EwaldSum
		interact_job[i].evar.bEfieldOnly = false; interact_job[i].evar.bVirial = true; interact_job[i].evar.bST_diagonal = true;
		interact_job[i].evar.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
		interact_job[i].evar.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	}

	if (mdcell->bHDynMC_T) {
		MDCELL_mass_center(*mdcell);
		// make sure the spatial momentum of the whole cell is zero
		//NPT_ZeroSpatialMomentum(job_nptcell, sv6job, MAX_THREADS);
		extern void ZeroSpatialMomentum(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob);
		ZeroSpatialMomentum(job_mdcell, sv3job, MAX_THREADS);
	}

	cmm->nCheckCluster = 1;
	int nCheckAngle = 15, iCheckAngle = 0;
	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;
	bool bSpeedVerlet = true;

	mdcell->Init_MD(); //copy the initial torsion angle, speed to the mdpar
	if (bPolarize) mdcell->init_dipole_hist();

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	long LOOPS = 0;
	int save_indx = 0;
	MD_SAVE *mdsave = NULL;
	if (mdcell->mm.n > 0) {LOOPS = mdcell->lfmd.m[0].LOOPS; mdsave = mdcell->lfmd.m[0].mmsave.mdsave;}
	else if (mdcell->sm.n > 0) {LOOPS = mdcell->lfsmd.m[0].LOOPS; mdsave = mdcell->lfsmd.m[0].smsave.mdsave;}
	else if (mdcell->pm.m > 0) {LOOPS = mdcell->lfpmd.m[0].LOOPS; mdsave = mdcell->lfpmd.m[0].pmsave.mdsave;}

	ITER_SAVE PVsave, trace_save;
	PVsave.set_file(mdsave->title, "pv");
	PVsave.set_iter(mdsave->nsave, mdsave->max_save);
	trace_save.set_file(mdsave->title, "trace");
	trace_save.set_iter(mdsave->nsave, mdsave->max_save);

	_rand_::RAND_ARRAY<int> randa_m, randa_c;

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command

		//molecule_check_periodic_cell(*mdcell); // check the molecule mass center is inside the cell
		//set_free_cell_mass_center(*mdcell, mass_center);
		//molecule_check_periodic_cell(*mdcell); // check the molecule mass center is inside the cell

		// we do re-position of mmcell coupling with CMM position check, relationship check
		// here we have to check the periodity of the molecule position
		//molecule_check_periodic_cell(mmcell);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_molecule_PeriodicCell_check), job_mdcell, nMDThreads);

		check_atom_rg(*mdcell, nMDThreads);

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalInertiaMoment), job_mdcell, nMDThreads);
		MDCELL_mass_center(*mdcell);

		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_InitClusterForce), job_mdcell, nMDThreads);

	// re-distribute / check clusters in CMM, use vp of each cluster only
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		if (bStorageChain) {
			if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
			else cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
		}
		else {
			if (MAX_THREADS == 1) {
				job_mdcellCMM[0].cmm = cmm;
				cmm->reset_cell_chain();
				MMCELL_CMM_cluster_allocate(job_mdcellCMM, status_job);
			}
			else {
				for (i = 0; i < nMDThreads; i++) job_mdcellCMM[i].cmm = job_cmm + i; // set the JOB-related CMM is the job-related CMM
				MTOperate1<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, bool>((void*)(&MMCELL_CMM_cluster_allocate), job_mdcellCMM, nMDThreads, status_job);
				for (i = 0; i < nMDThreads; i++) {
					if (!status_job[i]) {
						sprintf(errmsg, "loop %d thread %d: buffer of CMM cell leaks when allocating cluster", nloop, i);
						show_log(errmsg, true);
					}
				}
				if (!combine_AtomAllocatedCMM_array(job_cmm, nMDThreads, cmm, true, true)) {
					sprintf(errmsg, "loop %d : buffer of CMM cell leaks when combining allocated cluster", nloop);
					show_log(errmsg, true);
				}
			}
			for (i = 0; i < MAX_THREADS; i++) job_mdcellCMM[i].cmm = cmm; // reset the JOB-related CMM is the main CMM
		}
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		show_infor(errmsg);
	#endif

		Ep = 0;
		// check relationship of each cluster, electrostatic and LJ
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		if (bStorageChain) {
			ParallelCheckRelshipInCells(*cmm, 0, cmm->acell.n - 1, MAX_THREADS); // setup relationship for each atom with CMM
		}
		else {
			MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >((void*)(&MMCELL_CMM_CheckRelship), job_mdcellCMM, nMDThreads);
		}
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
		show_infor(errmsg);
	#endif
		
		// important : calculate LJ interaction first, here cluster force & torque are reset
		LJ_Interact(mdcell, cmm, MAX_THREADS, LJ_var, LJ_res);
		U_LJ = LJ_res.U;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!::local_interact && !mdcell->eNeutral) {
			spme_var.bEfieldOnly = false;
			// important : calculate Coulomb secondly, here cluster force & torque are accumulated
			if (::bPolarize) {
				if (nloop > 3) Guess_induced_dipole(*mdcell, MAX_THREADS);
				// calculat the total dipole of the whole cell
				cal_dipole(*mdcell, MAX_THREADS, mdcell->mu_surf); memcpy(&(cmm->mu_surf), &(mdcell->mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				if (!Polarize_SPME_Interact(mdcell, cmm, vspme_mu, MAX_THREADS, true, spme_var, spme_res)) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
				}
			}
			else {
				// calculat the total dipole of the whole cell
				cal_dipole(*mdcell, MAX_THREADS, mdcell->mu_surf); memcpy(&(cmm->mu_surf), &(mdcell->mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				SPME_Interact(mdcell, cmm, vspme_q, MAX_THREADS, spme_var, spme_res);
			}

			U_estat = spme_res.U;
		}
#endif
		Ep = U_LJ + U_estat;

		virial.v[0] = LJ_res.STxx + spme_res.STxx;
		virial.v[1] = LJ_res.STyy + spme_res.STyy;
		virial.v[2] = LJ_res.STzz + spme_res.STzz;
		virial.v[3] = LJ_res.STtr + spme_res.STtr;
		
		PVsave.one_more(); trace_save.one_more();

		
		// Explicit calculation of electrostatic interaction, iterating over image cell
		/*
		MTOperate1< MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, INTERACT>((void*)(&MMCELL_ExplicitInteract), job_mdcellCMM, nMDThreads, interact_job);
		Ep = 0; virial.v[0] = 0; virial.v[1] = 0; virial.v[2] = 0; virial.v[3] = 0;
		for (i = 0; i < nMDThreads; i++) {
			Ep += interact_job[i].LJ_res.U + interact_job[i].eres.U; 
			virial.v[0] += interact_job[i].LJ_res.STxx + interact_job[i].eres.STxx;
			virial.v[1] += interact_job[i].LJ_res.STyy + interact_job[i].eres.STyy;
			virial.v[2] += interact_job[i].LJ_res.STzz + interact_job[i].eres.STzz;
			virial.v[3] += interact_job[i].LJ_res.STtr + interact_job[i].eres.STtr;
		}
		*/

		//sprintf(errmsg, "U : %f, Virial : %f", Ep, virial.v[3] * ::unit_Ek_kT); show_log(errmsg, true);

		if (mdcell->mm.n > 0) {
			MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, TorsionVar, InteractRes>((void*)(&MMCELL_MM_Torsion), job_mdcell, nMDThreads, torsion_var, mires);
			for (i = 0; i < nMDThreads; i++) {
				Ep += mires[i].U; Etorsion += mires[i].U;
				if (torsion_var[i].bVirial) {
					virial.v[0] += mires[i].STxx;
					virial.v[1] += mires[i].STyy;
					virial.v[2] += mires[i].STzz;
					virial.v[3] += mires[i].STtr;
				}
			}
			//sprintf(errmsg, "%d: Torsion %f kT", nloop, Etorsion); show_log(errmsg, true);
		}

		if (::bEext) {
			if (mdcell->mm.m != NULL) {
				MOperate3<MMOLECULE, float, double>((void*)(&MM_ExternalEField), mdcell->mm.m, mdcell->mm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mdcell->sm.m != NULL) {
				MOperate3<SMOLECULE, float, double>((void*)(&SM_ExternalEField), mdcell->sm.m, mdcell->sm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mdcell->pm.m != NULL) {
				MOperate3<PMOLECULE, float, double>((void*)(&PM_ExternalEField), mdcell->pm.m, mdcell->pm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
		}

		//if (mmcell.mm.m != NULL) MOperate<MMOLECULE>(MM_TOTAL_FORCE, mmcell.mm.m, mmcell.mm.n, MAX_THREADS);
		//if (mmcell.sm.m != NULL) MOperate<SMOLECULE>(SM_TOTAL_FORCE, mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
		//if (mmcell.pm.m != NULL) MOperate<PMOLECULE>(PM_TOTAL_FORCE, mmcell.pm.m, mmcell.pm.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ClusterForce), job_mdcell, nMDThreads);

		// correction of virial, for rigid molecules
		
		if (mdcell->mm.n > 0 || mdcell->sm.n > 0) {
			for (i = 0; i < nMDThreads; i++) {
				sv4job[i].v[0] = 0; sv4job[i].v[1] = 0; sv4job[i].v[2] = 0; sv4job[i].v[3] = 0;
			}
			switch (::method_Virial) {
			case CLUSTER_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_VirialCorrect_RigidCluster), job_mdcell, nMDThreads, sv4job);
				break;
			case MOLECULE_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_VirialCorrect_Molecule), job_mdcell, nMDThreads, sv4job);
				break;
			}
			for (i = 0; i < nMDThreads; i++) {
				virial.v[0] -= sv4job[i].v[0];
				virial.v[1] -= sv4job[i].v[1];
				virial.v[2] -= sv4job[i].v[2];
				virial.v[3] -= sv4job[i].v[0] + sv4job[i].v[1] + sv4job[i].v[2];
			}
		}
		
		npt_cell.RF = virial.v[3];
		npt_cell.virial.v[0] = virial.v[0];
		npt_cell.virial.v[1] = virial.v[1];
		npt_cell.virial.v[2] = virial.v[2];
		
		// IMPORTANT: the forces from here on needs to calculate the induced torque explicitly
		// randome force and ImplicitSolvent force, and the induced torque, will be calculated explicitly
		// so we do not need to call MDCELL_ClusterForce.

		if (bRandomForce) gen_random_force(*mdcell, randa_m, randa_c);

		if (bImplicitSolvent) {
			// neighbors of cluster for implicit solvation 
			//if (bCheckClusterInCMM) CellOperate<BASIC_CLUSTER>((void*)(&CMMCheckImpSolvNeighbors_Cell), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			// solvation energy
			//CellOperate<BASIC_CLUSTER>((void*)(&CMM_ImplicitSolvationForce), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			//for (nm = 0; nm < mmcell.mm.n; nm++) {
			//	mm = mmcell.mm.m + nm;
			//	for (nc = 0; nc < mm->nCluster; nc++) Ep += mm->cluster[nc].E_ImplicitSolvation;
			//}
			MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >((void*)(&MMCELL_CMM_CheckImpSolvRelship), job_mdcellCMM, nMDThreads);
			MTOperate1<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, double >((void*)(&MMCELL_cluster_impsolv_interact), job_mdcellCMM, nMDThreads, djob);
			for (i = 0; i < nMDThreads; i++) Ep += djob[i];
		}

		//MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ScaleClusterForce), mdvar.job_mdcell, nMDThreads);
		SCALE_FORCE_CMM(*cmm, 0, cmm->acell.n - 1);

		// self-consistent velocity and acceleration, Nose-Hoover Chain
		bSpeedVerlet = (nloop == 0 ? false : true);
		if (!NPT_LeapFrogVerlet_SelfConsistent_SpeedAccel(&npt_cell, job_mdcell, job_nptcell, mdcellHDynThreadVar, hinge_constraint_var, nMDThreads, bSpeedVerlet)) {
			sprintf(errmsg, "failure to get consistent velocity and acceleration at loop %d", nloop);
			show_log(errmsg, true);
		}

		Ep_t = Ep; Ek_t = 0;
		for (nm = 0; nm < mdcell->mm.n; nm++) Ek_t += mdcell->mm.m[nm].Ek;
		for (nm = 0; nm < mdcell->sm.n; nm++) Ek_t += mdcell->sm.m[nm].Ek;
		for (nm = 0; nm < mdcell->pm.n; nm++) Ek_t += mdcell->pm.m[nm].Ek;

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);

		mdcell->MDSave(nloop, (float)Ep_t, false);
		// save npt
		npt_cell.pmdpar.save_dyn_pars(npt_cell.pdyn, nloop, true, NULL, true);
		npt_cell.pmdpar.psave.save_kinetic_pars(npt_cell.pdyn, nloop, true, NULL, true);
		
		{	
			float Ekt = npt_cell.Ekt, RF = npt_cell.RF;
			if (PVsave.bsave()) {

			//if (npt_cell.pmdpar.psave.mdsave->mKinSave.out != NULL) {
			//	*(npt_cell.pmdpar.psave.mdsave->mKinSave.out)<<Ekt * 2<<" + "<<RF * ::unit_Ek_kT<<",  U = "<<Ep_t<<endl<<endl;
			//}

			/*
			MTOperate1<MDCELL_ThreadJob< ISOTROPIC_NPT_CELL>, SVECTOR<double, 6> >(&NPT_MDCELL_SpatialMomentum, job_nptcell, MAX_THREADS, sv3job);
			VECTOR3 P, Pw;
			for (i = 0; i < MAX_THREADS; i++) {
				Pw.v[0] += sv6job[i].v[0];
				Pw.v[1] += sv6job[i].v[1];
				Pw.v[2] += sv6job[i].v[2];
				P.v[0] += sv6job[i].v[3];
				P.v[1] += sv6job[i].v[4];
				P.v[2] += sv6job[i].v[5];
			}
			//if (npt_cell.pmdpar.psave.mdsave->mKinSave.out != NULL) {
			//	*(npt_cell.pmdpar.psave.mdsave->mKinSave.out)<<"Pw = ["<<Pw.v[0]<<", "<<Pw.v[1]<<", "<<Pw.v[2]<<"]"<<endl;
			//	*(npt_cell.pmdpar.psave.mdsave->mKinSave.out)<<"P = ["<<P.v[0]<<", "<<P.v[1]<<", "<<P.v[2]<<"]"<<endl<<endl;
			//}

			//VECTOR3 Rcm;
			//MDCELL_mass_center(*mdcell, Rcm);
			(*PVsave.out)<<nloop<<"  "<<npt_cell.pdyn.Pint<<"  "<<npt_cell.pdyn.V<<"  "<<RF * ::unit_Ek_kT<<"  "<<Rcm.v[0]<<"  "<<Rcm.v[1]<<"  "<<Rcm.v[2]<<endl;

			MpV3(mdcell->I, mdcell->Wmc, Pw)
			Ekt = mdcell->M * (mdcell->Vmc.v[0] * mdcell->Vmc.v[0] + mdcell->Vmc.v[1] * mdcell->Vmc.v[1] + mdcell->Vmc.v[2] * mdcell->Vmc.v[2]);
			Ekt *= ::unit_Ek_kT;
			float Ekr = 0;
			Ekr = mdcell->Wmc.v[0] * Pw.v[0] + mdcell->Wmc.v[1] * Pw.v[1] + mdcell->Wmc.v[2] * Pw.v[2];
			Ekr *= ::unit_Ek_kT;
			*/

			(*PVsave.out)<<nloop<<"  "<<npt_cell.pdyn.Pint<<"  "<<npt_cell.pdyn.V;//<<endl;
			(*PVsave.out)<<"  "<<npt_cell.virial.v[0] * ::unit_Ek_kT;
			(*PVsave.out)<<"  "<<npt_cell.virial.v[1] * ::unit_Ek_kT;
			(*PVsave.out)<<"  "<<npt_cell.virial.v[2] * ::unit_Ek_kT;
			(*PVsave.out)<<"  "<<npt_cell.Ekt<<endl;
			/*
			if (trace_save.bsave()) {
				(*trace_save.out)<<nloop;
				double *v;
				for (i = 15; i < 30; i++) {
					if (i >= mdcell->pm.n) break;
					//RF = mdcell->sm.m[i].Rg.v[0] * mdcell->sm.m[i].c->dyn->fc.v[3];
					//RF += mdcell->sm.m[i].Rg.v[1] * mdcell->sm.m[i].c->dyn->fc.v[4];
					//RF += mdcell->sm.m[i].Rg.v[2] * mdcell->sm.m[i].c->dyn->fc.v[5];
					v = mdcell->pm.m[i].c->dyn->fc.v + 3;
					//v = mdcell->pm.m[i].Rg.v;
					//RF = v[0];
					RF = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; RF = sqrt(RF);
					(*trace_save.out)<<"  "<<RF;
			}
			(*trace_save.out)<<endl;
			*/
			}
		}
		

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob< ISOTROPIC_NPT_CELL > >((void*)(&MDCELL_NPT_Verlet), job_nptcell, nMDThreads);
		NPT_Verlet(npt_cell.pmdpar, npt_cell.pdyn);

		// V of the cell is changed now, we need to reset the cmm and mdcell, and clusters distributed in CMM
		MDCellVolUpdate(mdcell, cmm, npt_cell.pdyn.V);
		for (i = 0; i < MAX_THREADS; i++) CMM_construct<BASIC_CLUSTER>(job_cmm + i, cmm); // Update the dimesion of job_cmm[i], which is used to distribute the cluster in CMM
	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!::local_interact) {
			if (bPolarize) spmeCellSizeUpdate<_VMUatom>(vspme_mu, cmm->xd[0], cmm->xd[1], cmm->xd[2]);
			else spmeCellSizeUpdate<_VQatom>(vspme_q, cmm->xd[0], cmm->xd[1], cmm->xd[2]);
			spme_var.esv1.init_EwaldSum(cmm->xd[0] * cmm->xd[1] * cmm->xd[2], ::rcut_Ewd);
		}
	#endif

		// accept current NHC vars and prepair for next step
		for (i = 0; i < mdcell->hdyn.n; i++) mdcell->hdyn.m[i].nhc->accept_nextstep();
		npt_cell.pmdpar.hdyn->nhc->accept_nextstep();
		mdcell->mcHDyn.nhc->accept_nextstep();

		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			if (mdcell->mm.n > 2) {
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CheckRotAngle), job_mdcell, nMDThreads);
			}
			else {
				for (nm = 0; nm < mdcell->mm.n; nm++) {
					mm = mdcell->mm.m + nm;
					check_angle(*mm, mdcell->lfmd.m[nm]);
				}
			}
		}
		iCheckAngle++; if (iCheckAngle >= nCheckAngle) iCheckAngle = 0;

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mdcell->CalibrateTorsionAngle();

		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "checking procedure takes %d ms", time.glance());
			show_infor(errmsg);
		#endif

		PVsave.save_over(); trace_save.save_over();

		nloop++;
	}
	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) job_cmm[i].reset_basic_cell_matrix();
	}
	PVsave.reset(); trace_save.reset();
}
