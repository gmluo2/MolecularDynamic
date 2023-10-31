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

#include "Interact_gc.h"
#include "mdvar.h"
#include "gcmd.h"

using namespace _EwaldSum_real_;


#include "interface.h"
using namespace _interface_constraint_;
extern void InterfaceConstraints_CMM(CMM_CELL3D<BASIC_CLUSTER> *cmm, int nThreads, double &Ut);


extern void ExplicitEInteract(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, EInteractVar& evar, InteractRes &res);

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

// **********************************************************************************************************************
// calculating the Hoover Force on each cluster
// **********************************************************************************************************************
void cal_Hoover_force(MMOLECULE *mm, double ksi) {
	int ic;
	CLUSTER *pc;
	VECTOR<6> *V;
	MATRIX<6> *M, *phai;
	SVECTOR<double, 6> f;
	double *v1, *v2;
	float mass;
	double *fh;

	int method = ::MM_SpeedVerlet_method;

	if (method == 0) {
		V6zero(mm->mDynKin.f_Hoover)
	}

	for (ic = 0; ic < mm->nCluster; ic++) {
		pc = mm->cluster + ic;
		mass = pc->M; 
		V = &(pc->bkm->V0); M = &(pc->invNE->M);
		MpV6((*M), (*V), f)

		// accumulate Hoover force
		fh = pc->dyn->f_Hoover.v;
		memset(fh, 0, SIZE_V6);

		fh[0] = -ksi * f.v[0];
		fh[1] = -ksi * f.v[1];
		fh[2] = -ksi * f.v[2];
		fh[3] = -ksi * f.v[3];
		fh[4] = -ksi * f.v[4];
		fh[5] = -ksi * f.v[5];
		/*
		v1 = pc->dyn->f.v;
		v2 = pc->dyn->fc.v;
		V6plus(v2, fh, v1)
		*/
		if (method == 0) {
			phai = &(pc->invNE->phai_cm2me);
			MpV6((*phai), pc->dyn->f_Hoover, f) // f = phai(cm=>pc) * f

			v1 = mm->mDynKin.f_Hoover.v; v2 = f.v; V6plus(v1, v2, v1) // mDynKin.f += f
		}
	}
}

void cal_Hoover_force(SMOLECULE *sm, double ksi) {
	VECTOR<6> *V = &(sm->c->bkm->V0);
	SVECTOR<double, 3> w;
	SVECTOR<double, 3> force, torque;
	MATRIX<3> *I = &(sm->I);
	float mass = sm->c->M;
	double *fh = sm->c->dyn->f_Hoover.v;

	w.v[0] = -ksi * V->v[0];
	w.v[1] = -ksi * V->v[1];
	w.v[2] = -ksi * V->v[2];
	MpV3((*I), w, torque)

	force.v[0] = -mass * (ksi * V->v[3]);
	force.v[1] = -mass * (ksi * V->v[4]);
	force.v[2] = -mass * (ksi * V->v[5]);
	
	// total virtual force
	fh[0] = torque.v[0];
	fh[1] = torque.v[1];
	fh[2] = torque.v[2];
	fh[3] = force.v[0];
	fh[4] = force.v[1];
	fh[5] = force.v[2];

	/*
	v1 = sm->c->dyn->f.v;
	v2 = sm->c->dyn->fc.v;
	V6plus(v2, fh, v1) // f = fc + f_Hoover
	*/
}

void cal_Hoover_force(PMOLECULE *pm, double ksi) {
	double *v = pm->v.v;
	double *fh = pm->c->dyn->f_Hoover.v;
	float mass = pm->c->M;

	memset(fh, 0, SIZE_V6);

	// total virtual force
	fh[3] = -mass *(ksi * v[0]);
	fh[4] = -mass *(ksi * v[1]);
	fh[5] = -mass *(ksi * v[2]);

	/*
	v1 = pm->c->dyn->f.v;
	v2 = pm->c->dyn->fc.v;
	v1[3] = v2[3] + fh[3];
	v1[4] = v2[4] + fh[4];
	v1[5] = v2[5] + fh[5];
	*/
}

void VMDCELL_SpeedVerlet(MDCELL_ThreadJob<MMOL_VMD_CELL> *job, bool *status) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	VECTOR3 f, torq;
	int nm = 0, imol = 0;

	double dt =  mdcell->hdyn.nhc->dt;
	double uConvg = 0.005 / dt, wConvg = 0.005 / dt; //0.017 is 1 arc degree;
	if (uConvg < 2e-4) uConvg = 2e-4;
	if (wConvg < 2e-4) wConvg = 2e-4;

	char msg_show[1024] = "\0", msg_log[1024] = "\0", msg[1024] = "\0";

	bool bSpeedVerlet = *status;

	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		if (mdcell->lfmd.m[imol]->iHDyn->hdyn_status) continue; // convergent already

		mdcell->lfmd.m[imol]->fetch_ksi();

		cal_Hoover_force(mm, mdcell->lfmd.m[imol]->ksi);

		if (::bBrownian) Brownian_force(mm, 0, mm->nCluster - 1, mm->mDynKin.f_Brown);

		mdcell->mvb.mm_status.m[imol] = MM_SpeedVerlet_Simplified(bSpeedVerlet, *(mdcell->mm.m[imol]), *(mdcell->lfmd.m[imol]), 
			mdcell->mvb.mmv1.m[imol].cv.m, mdcell->mvb.mmv2.m[imol].cv.m, 
			mdcell->mvb.mm_dvmax.m[imol], mdcell->mvb.mm_dwmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.mm_status.m[imol]) {
			sprintf(msg, "Macromol # %d [%s] : Speed-Verlet fails", imol, mdcell->mm.m[imol]->mol);
			show_log(msg, true); //show_infor(msg_show);
		}
		if (mdcell->mvb.mm_dvmax.m[imol] <= uConvg && mdcell->mvb.mm_dwmax.m[imol] <= wConvg) mdcell->mvb.mm_sc_status.m[imol] = true;
		else mdcell->mvb.mm_sc_status.m[imol] = false;
	}

	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		if (mdcell->lfsmd.m[imol]->hdyn->hdyn_status) continue;  // convergent already

		cal_Hoover_force(sm, mdcell->lfsmd.m[imol]->hdyn->nhc->ksi.v[0]);

		if (::bBrownian) SM_Brownian_force(sm);

		mdcell->mvb.sm_status.m[imol] = SM_SpeedVerlet_Simplified(bSpeedVerlet, *(mdcell->sm.m[imol]), *(mdcell->lfsmd.m[imol]), 
			mdcell->mvb.sm_dvmax.m[imol], mdcell->mvb.sm_dwmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.sm_status.m[imol]) {
			sprintf(msg, "mol # %d [%s] : Speed-Verlet fails", imol, mdcell->sm.m[imol]->mol);
			show_log(msg_show, true); //show_log(msg_log, true);
		}
		if (mdcell->mvb.sm_dvmax.m[imol] <= uConvg && mdcell->mvb.sm_dwmax.m[imol] < wConvg) mdcell->mvb.sm_sc_status.m[imol] = true;
		else mdcell->mvb.sm_sc_status.m[imol] = false;
	}

	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
		if (mdcell->lfpmd.m[imol]->hdyn->hdyn_status) continue;  // convergent already

		cal_Hoover_force(pm, mdcell->lfpmd.m[imol]->hdyn->nhc->ksi.v[0]);

		if (::bBrownian) PM_Brownian_force(pm);

		mdcell->mvb.pm_status.m[imol] = PM_SpeedVerlet_Simplified(bSpeedVerlet, *(mdcell->pm.m[imol]), *(mdcell->lfpmd.m[imol]), 
			mdcell->mvb.pm_dvmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.pm_status.m[imol]) {
			sprintf(msg, "mol # %d [%s] : Speed-Verlet fails", imol, mdcell->pm.m[imol]->mol);
			show_log(msg_show, true); //show_log(msg_log, true);
		}
		if (mdcell->mvb.pm_dvmax.m[imol] <= uConvg) mdcell->mvb.pm_sc_status.m[imol] = true;
		else mdcell->mvb.pm_sc_status.m[imol] = false;
	}
}

void VMDCELL_Verlet(MDCELL_ThreadJob<MMOL_VMD_CELL> *job) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0, ihdyn = 0;

	//char msg_show[1024] = "\0", msg_log[1024] = "\0";

	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		Verlet(*(mdcell->lfmd.m[imol]), *(mdcell->mm.m[imol])); // Verlet algorithm to calculate the ta and tp for next timestep
	}

	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		Verlet(*(mdcell->lfsmd.m[imol]), *(mdcell->sm.m[imol])); // Verlet algorithm to calculate the ta and tp for next timestep
	}

	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		Verlet(*(mdcell->lfpmd.m[imol]), *(mdcell->pm.m[imol])); // Verlet algorithm to calculate the ta and tp for next timestep
	}
}

void VMDCELL_CalMolKineticEnergy(MDCELL_ThreadJob<MMOL_VMD_CELL> *job) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0;
	int method = ::MM_SpeedVerlet_method;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		if (method == 0) MakeFreeBaseVelocityConsistent(*mm);
		else if (method == 1) {calClusterVelocity0(*mm); calSpatialMomentum0(*mm, true);}
		calKineticEnergy(*mm);
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		calKineticEnergy(*sm);
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
		calKineticEnergy(*pm);
	}
}

void VMDCELL_CalMolNetForce(MDCELL_ThreadJob<MMOL_VMD_CELL> *job) {
	int nm = 0, imol = 0;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= job->job->mdcell->mm.n) continue;
		mm = job->job->mdcell->mm.m[imol];
		CalNetForce(*mm, 0, mm->nCluster); // calculate the net force also here !
	}
	/*
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < smIndx.n; nm++) {
		imol = smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		sm->Ek = calKineticEnergy(*sm);
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < pmIndx.n; nm++) {
		imol = pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
		pm->Ek = calKineticEnergy(*pm);
	}
	*/
}


void VMDCELL_HDyn(MMOL_VMD_CELL* mdcell) {
	int i = 0, im;
	double Ek1;
	double T = 0;

	char msg[256] = "\0";
	bool status = true;

	Ek1 = 0;
	for (im = 0; im < mdcell->mm.n; im++) {
		Ek1 += mdcell->mm.m[im]->Ek;
		if (!mdcell->mvb.mm_status.m[im]) status = false;
	}
	for (im = 0; im < mdcell->sm.n; im++) {
		Ek1 += mdcell->sm.m[im]->Ek;
		if (!mdcell->mvb.sm_status.m[im]) status = false;
	}
	for (im = 0; im < mdcell->pm.n; im++) {
		Ek1 += mdcell->pm.m[im]->Ek;
		if (!mdcell->mvb.pm_status.m[im]) status = false;
	}

	T = Ek1 / (mdcell->hdyn.nhc->N * ::Ek0Internal);
	mdcell->hdyn.nhc->set_dT(T - 1);
	mdcell->hdyn.nhc->verlet_NHC();

	//sprintf(msg, "Hdyn : T[%f], ksi[%f], vksi[%f]", T, mdcell->hdyn.nhc->ksi.v[0], mdcell->hdyn.nhc->vksi.v[0]);
	//show_log(msg, true);

	//if (!status) {
	//	sprintf(msg, "NHC #%d is reset", ihdyn); show_log(msg, true);
	//}
}

void VMDCELL_CheckHDynCovergent(MMOL_VMD_CELL* mdcell) {
	int i = 0, im;

	mdcell->hdyn.hdyn_status = true;
	for (im = 0; im < mdcell->mm.n; im++) {
		if (!mdcell->mvb.mm_status.m[im] || mdcell->mvb.mm_sc_status.m[im]) continue; // LFMD speed is rescaled, or the speed/accel is covergent
		else {mdcell->hdyn.hdyn_status = false; return;}
	}
	for (im = 0; im < mdcell->sm.n; im++) {
		if (!mdcell->mvb.sm_status.m[im] || mdcell->mvb.sm_sc_status.m[im]) continue; // LFMD speed is rescaled, or the speed/accel is covergent
		else {mdcell->hdyn.hdyn_status = false; return;}
	}
	for (im = 0; im < mdcell->pm.n; im++) {
		if (!mdcell->mvb.pm_status.m[im] || mdcell->mvb.pm_sc_status.m[im]) continue; // LFMD speed is rescaled, or the speed/accel is covergent
		else {mdcell->hdyn.hdyn_status = false; return;}
	}
}

void VMDCELL_AcceptCurrentDynForNextStep(MMOL_VMD_CELL* mdcell) {
	mdcell->hdyn.nhc->accept_nextstep();
}

bool LeapFrogVerlet_SelfConsistent_SpeedAccel(MMOL_VMD_CELL *mdcell, MDCELL_ThreadJob<MMOL_VMD_CELL> *job_mdcell, bool bSpeedVerlet, int nThreads) {
	//extern void MDCELL_SpatialMomentum(MDCELL_ThreadJob<MMOL_VMD_CELL> *job, SVECTOR<double, 3> *P);
	/*
	ARRAY< SVECTOR<double, 6> > aP; aP.SetArray(nThreads);
	*/
	ARRAY<bool> status; status.SetArray(nThreads);

	// before velocity verlet, backup the real force on the cluster, because the force on cluster (fc) will be adjusted by mass center velocity constraint
	//MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&MDCELL_ClusterRealForceBackup), job_mdcell, nThreads);

	extern void VMDCELL_MMClusterVelocity(MDCELL_ThreadJob<MMOL_VMD_CELL> *job);
	if (mdcell->mm.n > 0) MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_MMClusterVelocity), job_mdcell, nThreads);
	// calculate the kinetic energy of each molecule, 
	// initializing the velocities of each cluster and the mass-center speed based on the Spatial Moment of frame
	MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_CalMolKineticEnergy), job_mdcell, nThreads); 

	bool bDynConvergent = false;
	int iloop = 0, max_loops = 8, i;

	float Ek = 0;

	//VECTOR3 w, Pw;

	mdcell->hdyn.hdyn_status = false;
	while (!bDynConvergent && iloop < max_loops) {
		// run NHC to get a new ksi
		VMDCELL_HDyn(mdcell);

		// run NHC for the mass-center velocity of the whole cell
		/*
		if (mdcell->bHDynMC_T) {
			MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 3> >(&MDCELL_SpatialMomentum, job_mdcell, nThreads, aP.m);
			V3zero(mdcell->Vmc) V3zero(Pw)
			for (i = 0; i < nThreads; i++) {
				Pw.v[0] += aP.m[i].v[0];
				Pw.v[1] += aP.m[i].v[1];
				Pw.v[2] += aP.m[i].v[2];
				mdcell->Vmc.v[0] += aP.m[i].v[3];
				mdcell->Vmc.v[1] += aP.m[i].v[4];
				mdcell->Vmc.v[2] += aP.m[i].v[5];
			}
			MpV3(mdcell->invI, Pw, w)
			V32V3(w, mdcell->Wmc)
			mdcell->Vmc.v[0] /= mdcell->M; mdcell->Vmc.v[1] /= mdcell->M; mdcell->Vmc.v[2] /= mdcell->M;

			if (mdcell->bHDynMC_T) {
				V3ABS2(mdcell->Vmc, Ek) Ek *= mdcell->M * 0.5; Ek *= unit_Ek_kT;
				mdcell->Ekt_mc = Ek * 0.5;
				mdcell->mcHDyn.nhc->set_dT(Ek); // supposedly the spatial momentum should be 0
				mdcell->mcHDyn.nhc->verlet_NHC(1e-6);
			}

		}
		*/

		// Do the Speed-Verlet
		for (i = 0; i < nThreads; i++) status.m[i] = bSpeedVerlet;
		MTOperate1< MDCELL_ThreadJob<MMOL_VMD_CELL>, bool>((void*)(&VMDCELL_SpeedVerlet), job_mdcell, nThreads, status.m); 

		// check whether the speed and accel are convergent
		VMDCELL_CheckHDynCovergent(mdcell);
		bDynConvergent = mdcell->hdyn.hdyn_status;

		iloop += 1;
	}
	return bDynConvergent;
}

void VMDCELL_InitClusterForce(MDCELL_ThreadJob<MMOL_VMD_CELL> *job) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0, nc = 0, iatom, i;
	CLUSTER *pc = NULL;

	MMOLECULE *mm = NULL;
	MATOM *patom = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = mm->cluster + nc;
			V6zero(pc->dyn->fc) V6zero(pc->dyn->f_Hoover)
			if (pc->fphinge != NULL) {V6zero((*(pc->fphinge)))}
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom; 
				V3zero(patom->F) V3zero(patom->E)
				for (i = 0; i < patom->f.n; i++) {V3zero(patom->f.m[i])}
				patom->bShocked = false;
			}
			pc->bShocked = false;
		}
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		V6zero(sm->c->dyn->fc) V6zero(sm->c->dyn->f_Hoover) 
		for (iatom = 0; iatom < sm->c->nAtoms; iatom++) {
			patom = sm->c->atom + iatom; 
			V3zero(patom->F) V3zero(patom->E)
			for (i = 0; i < patom->f.n; i++) {V3zero(patom->f.m[i])}
			patom->bShocked = false;
		}
		sm->c->bShocked = false;
	}
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
		V6zero(pm->c->dyn->fc) V6zero(pm->c->dyn->f_Hoover) 
		for (iatom = 0; iatom < pm->c->nAtoms; iatom++) {
			patom = pm->c->atom + iatom; 
			V3zero(patom->F) V3zero(patom->E)
			for (i = 0; i < patom->f.n; i++) {V3zero(patom->f.m[i])}
			patom->bShocked = false;
		}
		pm->c->bShocked = false;
	}
}

void VMDCELL_ClusterForce(MDCELL_ThreadJob<MMOL_VMD_CELL> *job) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0, nc = 0;
	CLUSTER *pc = NULL;

	MMOLECULE *mm = NULL;
	MATOM *patom = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = mm->cluster + nc;
			V6zero(pc->dyn->fc)
			CLUSTER_TOTAL_FORCE(pc);
		}
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		V6zero(sm->c->dyn->fc)
		CLUSTER_TOTAL_FORCE(sm->c);
	}
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
		//V6zero(pm->c->dyn->fc) 
		//CLUSTER_TOTAL_FORCE(pm->c);
		pm->c->dyn->fc.v[3] = pm->c->atom[0].F.v[0];
		pm->c->dyn->fc.v[4] = pm->c->atom[0].F.v[1];
		pm->c->dyn->fc.v[5] = pm->c->atom[0].F.v[2];
	}
}


void MM_atomic_rg(MMOLECULE **mm, int nMM) {
	int imol, ic, iatom;
	CLUSTER *pc;
	MATOM *patom;
	for (imol = 0; imol < nMM; imol++) {
		for (ic = 0; ic < mm[imol]->nCluster; ic++) {
			pc = mm[imol]->cluster + ic;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = mm[imol]->cluster[ic].atom + iatom;
				V3plusV3(patom->r, pc->dr_cmm, patom->rg)
			}
		}
	}
}

void SM_atomic_rg(SMOLECULE**sm, int nSM) {
	int imol, iatom;
	MATOM *patom;
	for (imol = 0; imol < nSM; imol++) {
		for (iatom = 0; iatom < sm[imol]->c->nAtoms; iatom++) {
			patom = sm[imol]->c->atom + iatom;
			V3plusV3(patom->r, sm[imol]->c->dr_cmm, patom->rg)
		}
	}
}

void PM_atomic_rg(PMOLECULE**pm, int nPM) {
	int imol;
	MATOM *patom;
	for (imol = 0; imol < nPM; imol++) {
		patom = pm[imol]->c->atom;
		V32V3(patom->r, patom->rg)
	}
}

void check_atom_rg(MMOL_VMD_CELL &mmcell, int nThreads) {
	if (mmcell.mm.n > 0) {
		MOperate<MMOLECULE*>((void*)(&MM_atomic_rg), mmcell.mm.m, mmcell.mm.n, (mmcell.mm.n > nThreads ? nThreads : mmcell.mm.n));
	}
	if (mmcell.sm.n > 0) {
		MOperate<SMOLECULE*>((void*)(&SM_atomic_rg), mmcell.sm.m, mmcell.sm.n, (mmcell.sm.n > nThreads ? nThreads : mmcell.sm.n));
	}
	if (mmcell.pm.n > 0) {
		MOperate<PMOLECULE*>((void*)(&PM_atomic_rg), mmcell.pm.m, mmcell.pm.n, (mmcell.pm.n > nThreads ? nThreads : mmcell.pm.n));
	}

	/*
	char msg[256] = "\0";
	int nm, ia;
	MATOM *patom;
	float xl[3] = {-mmcell.h[0], -mmcell.h[1], -mmcell.h[2]};
	float xr[3] = {mmcell.h[0], mmcell.h[1], mmcell.h[2]};
	xl[0] *= 1.5; xl[1] *= 1.5; xl[2] *= 1.5;
	xr[0] *= 1.5; xr[1] *= 1.5; xr[2] *= 1.5;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		for (ia = 0; ia < mmcell.sm.m[nm].c->nAtoms; ia++) {
			patom = mmcell.sm.m[nm].c->atom + ia;
			if (patom->rg.v[0] < xl[0] || patom->rg.v[0] > xr[0] ||
				patom->rg.v[1] < xl[1] || patom->rg.v[1] > xr[1] ||
				patom->rg.v[2] < xl[2] || patom->rg.v[2] > xr[2]) {
					sprintf(msg, "strange: SM %d atom %d has unit coordinates [%f, %f, %f]", nm, ia, patom->rg.v[0], patom->rg.v[1], patom->rg.v[2]);
					show_log(msg, true);
			}
		}
	}
	*/
}

void update_atom_rg(MMOL_VMD_CELL &mmcell) {
	if (mmcell.mm.n > 0) MM_atomic_rg(mmcell.mm.m, mmcell.mm.n);
	if (mmcell.sm.n > 0) SM_atomic_rg(mmcell.sm.m, mmcell.sm.n);
	if (mmcell.pm.n > 0) PM_atomic_rg(mmcell.pm.m, mmcell.pm.n);
}


void V_MultiMM_NEIMO_MOMENT(MMOLECULE **mm, int nmm) {
	for (int imol = 0; imol < nmm; imol++) {
		CalMolInertiaMoment(*(mm[imol]));
		NEIMO_MOMENT_MATRIX(*(mm[imol]));
	}
}

void VMDCELL_CheckRotAngle(MDCELL_ThreadJob<MMOL_VMD_CELL> *job) {
	int nm = 0, imol = 0, ihdyn = 0;

	char msg_show[1024] = "\0", msg_log[1024] = "\0";

	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= job->job->mdcell->mm.n) continue;
		check_angle(*(job->job->mdcell->mm.m[imol]), *(job->job->mdcell->lfmd.m[imol]));
	}
	/*
	for (nm = 0; nm < smIndx.n; nm++) {
		imol = smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		// work to be done
	}

	for (nm = 0; nm < pmIndx.n; nm++) {
		imol = pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		// work to be done
	}
	*/
}

void VMMCELL_CMM_cluster_allocate(MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator, bool *status) {
	int i, imol;
	bool tstatus;
	*status = true;
	CMM_array_reset_storage<BASIC_CLUSTER>(allocator->cmm);
	for (i = 0; i < allocator->job->mmIndx.n; i++) {
		imol = allocator->job->mmIndx.m[i];
		if (imol >= allocator->job->mdcell->mm.n) continue;
		tstatus = MM_AllocateCMM_array<BASIC_CLUSTER, MMOLECULE>(allocator->job->mdcell->mm.m[imol], allocator->cmm); 
		if (!tstatus) *status = false;
	}
	for (i = 0; i < allocator->job->smIndx.n; i++) {
		imol = allocator->job->smIndx.m[i];
		if (imol >= allocator->job->mdcell->sm.n) continue;
		tstatus = SM_AllocateCluster_CMM_array(allocator->job->mdcell->sm.m[imol], 1, allocator->cmm); 
		if (!tstatus) *status = false;
	}
	for (i = 0; i < allocator->job->pmIndx.n; i++) {
		imol = allocator->job->pmIndx.m[i];
		if (imol >= allocator->job->mdcell->pm.n) continue;
		tstatus = PM_AllocateCluster_CMM_array(allocator->job->mdcell->pm.m[imol], 1, allocator->cmm); 
		if (!tstatus) *status = false;
	}
}

void VMMCELL_CMM_CheckRelship(MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator) {
	int i, imol, nc;
	MMOLECULE *mm = NULL;
	BASIC_CLUSTER *pc = NULL;
	for (i = 0; i < allocator->job->mmIndx.n; i++) {
		imol = allocator->job->mmIndx.m[i];
		if (imol >= allocator->job->mdcell->mm.n) continue;
		mm = allocator->job->mdcell->mm.m[imol];
		for (nc = 0; nc < mm->nCluster; nc++) {
			CheckClusterRelship((BASIC_CLUSTER*)(mm->cluster + nc), *(allocator->cmm)); 
		}
	}
	for (i = 0; i < allocator->job->smIndx.n; i++) {
		imol = allocator->job->smIndx.m[i];
		if (imol >= allocator->job->mdcell->sm.n) continue;
		CheckClusterRelship(allocator->job->mdcell->sm.m[imol]->c, *(allocator->cmm)); 
	}
	for (i = 0; i < allocator->job->pmIndx.n; i++) {
		imol = allocator->job->pmIndx.m[i];
		if (imol >= allocator->job->mdcell->pm.n) continue;
		CheckClusterRelship(allocator->job->mdcell->pm.m[imol]->c, *(allocator->cmm)); 
	}
}

void VMMCELL_CMM_CheckImpSolvRelship(MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator) {
	extern float rSolvent, dSolv_max;
	//float r0_solvent = rSolvent;
	CMM_CELL3D<BASIC_CLUSTER> *cmm = allocator->cmm;
	float xd[3] = {cmm->xd[0] / cmm->bcell.nx, cmm->xd[1] / cmm->bcell.ny, cmm->xd[2] / cmm->bcell.nz};
	int dnc[3];
	float rcut = dSolv_max;
	
	if (::local_relax) {dnc[0] = 1; dnc[1] = 1; dnc[2] = 1;}
	else {
		dnc[0] = int(rcut / xd[0] + 0.5f); if (dnc[0] == 0) dnc[0] = 1; if (dnc[0] * xd[0] < rcut) dnc[0] += 1;
		dnc[1] = int(rcut / xd[1] + 0.5f); if (dnc[1] == 0) dnc[1] = 1; if (dnc[1] * xd[1] < rcut) dnc[1] += 1;
		dnc[2] = int(rcut / xd[2] + 0.5f); if (dnc[2] == 0) dnc[2] = 1; if (dnc[2] * xd[2] < rcut) dnc[2] += 1;
		dnc[0] += 1; dnc[1] += 1; dnc[2] += 1;
	}

	int i, imol, nc;
	MMOLECULE *mm = NULL;
	BASIC_CLUSTER *pc = NULL;
	for (i = 0; i < allocator->job->mmIndx.n; i++) {
		imol = allocator->job->mmIndx.m[i];
		if (imol >= allocator->job->mdcell->mm.n) continue;
		mm = allocator->job->mdcell->mm.m[imol];
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = (BASIC_CLUSTER*)(mm->cluster + nc);
			ClusterImpSolvNeighborsInCell<BASIC_CLUSTER>(cmm, pc, dnc); 
		}
	}
	for (i = 0; i < allocator->job->smIndx.n; i++) {
		imol = allocator->job->smIndx.m[i];
		if (imol >= allocator->job->mdcell->sm.n) continue;
		pc = allocator->job->mdcell->sm.m[imol]->c;
		ClusterImpSolvNeighborsInCell<BASIC_CLUSTER>(cmm, pc, dnc); 
	}
	for (i = 0; i < allocator->job->pmIndx.n; i++) {
		imol = allocator->job->pmIndx.m[i];
		if (imol >= allocator->job->mdcell->pm.n) continue;
		pc = allocator->job->mdcell->pm.m[imol]->c;
		ClusterImpSolvNeighborsInCell<BASIC_CLUSTER>(cmm, pc, dnc); 
	}
}

void VMMCELL_cluster_impsolv_interact(MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *mdcell_impsolv, double *U) {
	int i, imol, nc;
	MMOLECULE *mm = NULL;
	BASIC_CLUSTER *pc = NULL;
	CMM_CELL3D<BASIC_CLUSTER> *cmm = mdcell_impsolv->cmm;
	MMOL_VMD_CELL *mdcell = mdcell_impsolv->job->mdcell;
	// calculate SASA of each cluster
	for (i = 0; i < mdcell_impsolv->job->mmIndx.n; i++) {
		imol = mdcell_impsolv->job->mmIndx.m[i];
		if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = (BASIC_CLUSTER*)(mm->cluster + nc);
			Cluster_SASA(pc, cmm);
		}
	}
	for (i = 0; i < mdcell_impsolv->job->smIndx.n; i++) {
		imol = mdcell_impsolv->job->smIndx.m[i];
		if (imol >= mdcell->sm.n) continue;
		pc = mdcell->sm.m[imol]->c;
		Cluster_SASA(pc, cmm);
	}
	for (i = 0; i < mdcell_impsolv->job->pmIndx.n; i++) {
		imol = mdcell_impsolv->job->pmIndx.m[i];
		if (imol >= mdcell->pm.n) continue;
		pc = mdcell->pm.m[imol]->c;
		Cluster_SASA(pc, cmm);
	}

	// calculate the force & energy
	double E = 0;
	for (i = 0; i < mdcell_impsolv->job->mmIndx.n; i++) {
		imol = mdcell_impsolv->job->mmIndx.m[i];
		if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = (BASIC_CLUSTER*)(mm->cluster + nc);
			ClusterImpSolvForce(pc);
			E += pc->E_ImplicitSolvation;
		}
	}
	for (i = 0; i < mdcell_impsolv->job->smIndx.n; i++) {
		imol = mdcell_impsolv->job->smIndx.m[i];
		if (imol >= mdcell->sm.n) continue;
		pc = mdcell->sm.m[imol]->c;
		ClusterImpSolvForce(pc);
		E += pc->E_ImplicitSolvation;
	}
	for (i = 0; i < mdcell_impsolv->job->pmIndx.n; i++) {
		imol = mdcell_impsolv->job->pmIndx.m[i];
		if (imol >= mdcell->pm.n) continue;
		pc = mdcell->pm.m[imol]->c;
		ClusterImpSolvForce(pc);
		E += pc->E_ImplicitSolvation;
	}
	*U = E;
}

//#error the function can work with reference as variables, while not pointer for variables, how can it work with MTOperate2 ?
//void MMCELL_MM_Torsion(MDCELL_ThreadJob<MMOL_MD_CELL> *job, TorsionVar &var, InteractRes &ires) {
void VMMCELL_MM_Torsion(MDCELL_ThreadJob<MMOL_VMD_CELL> *job, TorsionVar *var, InteractRes *ires) {
	int i, imol;
	MMOLECULE *mm = NULL;
	ires->reset();
	InteractRes tres;
	for (i = 0; i < job->job->mmIndx.n; i++) {
		imol = job->job->mmIndx.m[i];
		if (imol >= job->job->mdcell->mm.n) continue;
		mm = job->job->mdcell->mm.m[imol];
		TorsionInteract(mm, 0, mm->nCluster, *var, tres);
		(*ires) += tres;
	}
}

void VMMCELL_ExplicitInteract(MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *job, INTERACT *interact) {
	int i, imol, ic;
	MMOLECULE *mm = NULL;
	CMM_CELL3D<BASIC_CLUSTER> *cmm = job->cmm;
	BASIC_CLUSTER * bpc = NULL;
	InteractRes tres;
	interact->LJ_res.reset(); interact->eres.reset();
	interact->evar.esv1.bEwald = false; // NON-EwaldSum, but direct interaction
	for (i = 0; i < job->job->mmIndx.n; i++) {
		imol = job->job->mmIndx.m[i];
		if (imol >= job->job->mdcell->mm.n) continue;
		mm = job->job->mdcell->mm.m[imol];
		for (ic = 0; ic < mm->nCluster; ic++) {
			bpc = (BASIC_CLUSTER*)(mm->cluster + ic);
			CLUSTER_LJ_INTERACT(*cmm, bpc, interact->var, tres);
			interact->LJ_res += tres;
			if (!bpc->eNeutral) {
				ExplicitEInteract(*cmm, bpc, interact->evar, tres);
				interact->eres += tres;
			}
		}
	}
	for (i = 0; i < job->job->smIndx.n; i++) {
		imol = job->job->smIndx.m[i];
		if (imol >= job->job->mdcell->sm.n) continue;
		bpc = job->job->mdcell->sm.m[imol]->c;
		CLUSTER_LJ_INTERACT(*cmm, bpc, interact->var, tres);
		interact->LJ_res += tres;
		if (!bpc->eNeutral) {
			ExplicitEInteract(*cmm, bpc, interact->evar, tres);
			interact->eres += tres;
		}
	}
	for (i = 0; i < job->job->pmIndx.n; i++) {
		imol = job->job->pmIndx.m[i];
		if (imol >= job->job->mdcell->pm.n) continue;
		bpc = job->job->mdcell->pm.m[imol]->c;
		CLUSTER_LJ_INTERACT(*cmm, bpc, interact->var, tres);
		interact->LJ_res += tres;
		if (!bpc->eNeutral) {
			ExplicitEInteract(*cmm, bpc, interact->evar, tres);
			interact->eres += tres;
		}
	}
}


void VMDCELL_VirialCorrect_RigidCluster(MDCELL_ThreadJob<MMOL_VMD_CELL> *job, SVECTOR<double, 4>* virial) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	int nm = 0, imol = 0, ic, ia;
	double STxx = 0, STyy = 0, STzz = 0;
	SVECTOR<double, 3> dr;
	SVECTOR<double, 3> F;
	double *f, *r;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			if (pc->fatom.m == NULL) {
				for (ia = 0; ia < pc->nAtoms; ia++) {
					patom = pc->atom + ia;
					VECT3(pc->cm0, patom->r0, dr) r = dr.v;
					if (::bEext) { // subtract the external force
						F.v[0] = patom->F.v[0];
						F.v[1] = patom->F.v[1];
						F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
						f = F.v;
					}
					else f = patom->F.v;
					STxx += r[0] * f[0];
					STyy += r[1] * f[1];
					STzz += r[2] * f[2];
				}
			}
			else {
				for (ia = 0; ia < pc->fatom.n; ia++) {
					patom = pc->fatom.m[ia];
					VECT3(pc->cm, patom->r, dr) r = dr.v;
					if (::bEext) { // subtract the external force
						F.v[0] = patom->F.v[0];
						F.v[1] = patom->F.v[1];
						F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
						f = F.v;
					}
					else f = patom->F.v;
					STxx += r[0] * f[0];
					STyy += r[1] * f[1];
					STzz += r[2] * f[2];
				}
			}
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		if (sm->c->fatom.m == NULL) {
			for (ia = 0; ia < sm->c->nAtoms; ia++) {
				patom = sm->c->atom + ia;
				// SM has mass center at r of molecule, so r0 is the displacement from mass center
				r = patom->r0.v;
				if (::bEext) { // subtract the external force
					F.v[0] = patom->F.v[0];
					F.v[1] = patom->F.v[1];
					F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
					f = F.v;
				}
				else f = patom->F.v;
				STxx += r[0] * f[0];
				STyy += r[1] * f[1];
				STzz += r[2] * f[2];
			}
		}
		else {
			for (ia = 0; ia < sm->c->fatom.n; ia++) {
				patom = sm->c->fatom.m[ia];
				// SM has mass center at r of molecule, so r0 is the displacement from mass center
				r = patom->r0.v;
				if (::bEext) { // subtract the external force
					F.v[0] = patom->F.v[0];
					F.v[1] = patom->F.v[1];
					F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
					f = F.v;
				}
				else f = patom->F.v;
				STxx += r[0] * f[0];
				STyy += r[1] * f[1];
				STzz += r[2] * f[2];
			}
		}
	}

	// point molecule has only one atom, no correction is required for virial
	/*
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
	}
	*/
	//STxx *= 0.5; STyy *= 0.5; STzz *= 0.5;
	virial->v[0] = STxx; virial->v[1] = STyy; virial->v[2] = STzz;
	virial->v[3] = STxx + STyy + STzz;
}

void VMDCELL_ClusterVirialCorrect_HingeConstraint(MDCELL_ThreadJob<MMOL_VMD_CELL> *job, HingeConstraintVar &var, InteractRes &ires) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL, *parent;
	HINGE_CONSTRAINT *hc;
	int nm = 0, imol = 0, ic, ia;
	double STxx = 0, STyy = 0, STzz = 0;
	SVECTOR<double, 3> dr;
	double *f, *r;
	double *v1, *v2, *v3;

	ires.reset();
	InteractRes tres;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		NEIMO_HingeForce(*mm);
		if (!HingeConstraintForce(*mm, var, tres)) continue;
		ires += tres;
	}
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			if (pc->hinge_constraint == NULL) continue;
			hc = pc->hinge_constraint;

			v1 = pc->cm.v; v3 = dr.v; r = dr.v;
			// 00
			v2 = hc->patom_00->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc0[0].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];
			
			// 01
			v2 = hc->patom_01->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc0[1].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];

			// 02
			v2 = hc->patom_02->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc0[2].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];

			// parent cluster
			parent = (CLUSTER*)(pc->parent);
			v1 = parent->cm.v; v3 = dr.v; r = dr.v;
			// 10
			v2 = hc->patom_10->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc1[0].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];
			
			// 11
			v2 = hc->patom_11->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc1[1].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];

			// 12
			v2 = hc->patom_12->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc1[2].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];
		}
	}
	ires.STxx -= STxx; 
	ires.STyy -= STyy;
	ires.STzz -= STzz;
	ires.STtr -= STxx + STyy + STzz;
}

void VMDCELL_VirialCorrect_Molecule(MDCELL_ThreadJob<MMOL_VMD_CELL> *job, SVECTOR<double, 4>* virial) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	int nm = 0, imol = 0, ic, ia;
	double STxx = 0, STyy = 0, STzz = 0;
	SVECTOR<double, 3> dr;
	SVECTOR<double, 3> F;
	double *f, *r;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			for (ia = 0; ia < pc->nAtoms; ia++) {
				patom = pc->atom + ia;
				VECT3(mm->mDynKin.cm, patom->r, dr) r = dr.v;
				if (::bEext) { // subtract the external force
					F.v[0] = patom->F.v[0];
					F.v[1] = patom->F.v[1];
					F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
					f = F.v;
				}
				else f = patom->F.v;
				STxx += r[0] * f[0];
				STyy += r[1] * f[1];
				STzz += r[2] * f[2];
			}
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		for (ia = 0; ia < sm->c->nAtoms; ia++) {
			patom = sm->c->atom + ia;
			// SM has mass center at r of molecule, so r0 is the displacement from mass center
			r = patom->r0.v;
			if (::bEext) { // subtract the external force
				F.v[0] = patom->F.v[0];
				F.v[1] = patom->F.v[1];
				F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
				f = F.v;
			}
			else f = patom->F.v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];
		}
	}

	// point molecule has only one atom, no correction is required for virial
	/*
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
	}
	*/
	//STxx *= 0.5; STyy *= 0.5; STzz *= 0.5;
	virial->v[0] = STxx; virial->v[1] = STyy; virial->v[2] = STzz;
	virial->v[3] = STxx + STyy + STzz;
}

void gen_random_force(MMOL_VMD_CELL &mdcell, _rand_::RAND_ARRAY<int> &randa_m, _rand_::RAND_ARRAY<int> &randa_c) {
	int imol, nc, i, i_m, i_c;
	SVECTOR<double, 3> f;
	BASIC_CLUSTER *pc = NULL;
	if (mdcell.mm.n > 0) {
		_rand_::init_random_array(randa_m, mdcell.mm.n);
		for (i_m = 0; i_m < randa_m.nw; i_m++) {
			imol = randa_m.a.m[i_m]; if (imol >= mdcell.mm.n) continue;
			_rand_::init_random_array(randa_c, mdcell.mm.m[imol]->nCluster);
			for (i_c = 0; i_c < mdcell.mm.m[imol]->nCluster; i_c++) {
				nc = randa_c.a.m[i_c]; if (nc >= mdcell.mm.m[imol]->nCluster) continue;
				pc = mdcell.mm.m[imol]->cluster + nc;
				gauss_random_force(randf, fwhm_randf, f.v, 1);
				for (i = 0; i < 3; i++) {
					pc->dyn->f.v[i+3] += f.v[i];
					pc->dyn->fc.v[i+3] += f.v[i];
				}
			}
		}
	}

	if (mdcell.sm.n > 0) {
		_rand_::init_random_array(randa_m, mdcell.sm.n);
		for (i_m = 0; i_m < randa_m.nw; i_m++) {
			imol = randa_m.a.m[i_m]; if (imol >= mdcell.sm.n) continue;
			pc = mdcell.sm.m[imol]->c;
			gauss_random_force(randf, fwhm_randf, f.v, 1);
			for (i = 0; i < 3; i++) {
				pc->dyn->f.v[i+3] += f.v[i];
				pc->dyn->fc.v[i+3] += f.v[i];
			}
		}
	}

	if (mdcell.pm.n > 0) {
		_rand_::init_random_array(randa_m, mdcell.pm.n);
		for (i_m = 0; i_m < randa_m.nw; i_m++) {
			imol = randa_m.a.m[i_m]; if (imol >= mdcell.pm.n) continue;
			pc = mdcell.pm.m[imol]->c;
			gauss_random_force(randf, fwhm_randf, f.v, 1);
			for (i = 0; i < 3; i++) {
				pc->dyn->fc.v[i+3] += f.v[i];
				pc->dyn->f.v[i+3] += f.v[i];
			}
		}
	}
}

extern void calTransSpatialMomentum(MMOLECULE &mm, VECTOR3 &P);
void VMDCELL_SpatialMomentum(MDCELL_ThreadJob<MMOL_VMD_CELL> *job, SVECTOR<double, 3> *mp) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0;
// do we need to check the rotation momentum ?
	VECTOR3 P, p;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		calTransSpatialMomentum(*mm, p); // velocity of each cluster at its mass-center
		V3plusV3(p, P, P)
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol]; M = sm->c->M;
		P.v[0] += sm->c->bkm->P0.v[3];
		P.v[1] += sm->c->bkm->P0.v[4];
		P.v[2] += sm->c->bkm->P0.v[5];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol]; M = pm->c->M;
		P.v[0] += M * pm->v.v[0];
		P.v[1] += M * pm->v.v[1];  
		P.v[2] += M * pm->v.v[2];
	}
	mp->v[0] = P.v[0]; mp->v[1] = P.v[1]; mp->v[2] = P.v[2];
}

void VMDCELL_ZeroSpatialMomentum(MDCELL_ThreadJob< MMOL_VMD_CELL> *job, SVECTOR<double, 3> *vp) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	VECTOR3 *v = &(mdcell->Vmc);

	if (FABS(v->v[0]) < 1e-8 && FABS(v->v[1]) < 1e-8 && FABS(v->v[2]) < 1e-8) return;
	else show_log("zero the total spatial momentum of the cell", true);

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		for (ic = 0; ic < mm->nCluster; ic++) {
			mm->cluster[ic].bkm->V0.v[3] -= vp->v[0];
			mm->cluster[ic].bkm->V0.v[4] -= vp->v[1];
			mm->cluster[ic].bkm->V0.v[5] -= vp->v[2];

			mm->cluster[ic].bkm->P0.v[3] -= mm->cluster[ic].M * vp->v[0];
			mm->cluster[ic].bkm->P0.v[4] -= mm->cluster[ic].M * vp->v[1];
			mm->cluster[ic].bkm->P0.v[5] -= mm->cluster[ic].M * vp->v[2];

			mm->cluster[ic].mcDyn->v.v[0] -= vp->v[0];
			mm->cluster[ic].mcDyn->v.v[1] -= vp->v[1];
			mm->cluster[ic].mcDyn->v.v[2] -= vp->v[2];
		}
		mm->V0.v[3] -= vp->v[0];
		mm->V0.v[4] -= vp->v[1];
		mm->V0.v[5] -= vp->v[2];
		mm->mDynKin.V.v[3] -= vp->v[0];
		mm->mDynKin.V.v[4] -= vp->v[1];
		mm->mDynKin.V.v[5] -= vp->v[2];

		mm->mDynKin.P.v[3] -= mm->mDynKin.M * vp->v[0];
		mm->mDynKin.P.v[4] -= mm->mDynKin.M * vp->v[1];
		mm->mDynKin.P.v[5] -= mm->mDynKin.M * vp->v[2];
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		sm->c->bkm->V0.v[3] -= vp->v[0];
		sm->c->bkm->V0.v[4] -= vp->v[1];
		sm->c->bkm->V0.v[5] -= vp->v[2];
		
		sm->c->bkm->P0.v[3] -= sm->c->M * vp->v[0];
		sm->c->bkm->P0.v[4] -= sm->c->M * vp->v[1];
		sm->c->bkm->P0.v[5] -= sm->c->M * vp->v[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
		pm->v.v[0] -= -vp->v[0];
		pm->v.v[1] -= vp->v[1];
		pm->v.v[2] -= vp->v[2];

		pm->c->bkm->P0.v[3] -= pm->c->M * vp->v[0];
		pm->c->bkm->P0.v[4] -= pm->c->M * vp->v[1];
		pm->c->bkm->P0.v[5] -= pm->c->M * vp->v[2];
	}
}

void ZeroSpatialMomentum(MDCELL_ThreadJob< MMOL_VMD_CELL> *job, SVECTOR<double, 3> *mp, int nJob) {
	MTOperate1<MDCELL_ThreadJob<MMOL_VMD_CELL>,SVECTOR<double, 3> >((void*)(&VMDCELL_SpatialMomentum), job, nJob, mp);
	VECTOR3 P, v;
	float M = job->job->mdcell->M;
	int i;
	for (i = 0; i < nJob; i++) {
		P.v[0] += mp[i].v[1]; P.v[1] += mp[i].v[2]; P.v[2] += mp[i].v[3];
	}
	v.v[0] = P.v[0] / M; v.v[1] = P.v[1] / M; v.v[2] = P.v[2] / M;
	for (i = 0; i < nJob; i++) memcpy(mp[i].v, v.v, SIZE_V3);

	MTOperate1<MDCELL_ThreadJob<MMOL_VMD_CELL>, SVECTOR<double, 3> >((void*)(&VMDCELL_ZeroSpatialMomentum), job, nJob, mp);
}

void VMDCELL_MMClusterVelocity(MDCELL_ThreadJob<MMOL_VMD_CELL> *job) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell; 
	int nm = 0, imol = 0;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		calClusterVelocity0(*mm);
	}
}

void cmm_init_distribute_cluster(SMOLECULE** sm, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain) {
	if (reset_cell_chain) cell.reset_cell_chain();

	int nm = 0, nc = 0;
	int ni = 0, nj = 0, nk = 0;
	VECTOR<3> dr, rc;
	int vi = 0;

	BASIC_CLUSTER* pc = NULL;
	for (nm = 0; nm < n; nm++) {
		pc = sm[nm]->c;
		V32V3((*(pc->vp)), rc)  // use vp position

		CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
		//cell_n(cell, ni, nj, nk, nc)
		//if (nc >= cell.nCells) nc = cell.nCells - 1;
		//cell.cell[nc].attach(pc);
		cell.bcell.m[ni].m[nj].m[nk].attach(pc);
	}
}

void cmm_init_distribute_cluster(PMOLECULE** pm, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain) {
	if (reset_cell_chain) cell.reset_cell_chain();

	int nm = 0, nc = 0;
	int ni = 0, nj = 0, nk = 0;
	VECTOR<3> dr, rc;
	int vi = 0;

	BASIC_CLUSTER* pc = NULL;
	for (nm = 0; nm < n; nm++) {
		pc = pm[nm]->c;
		V32V3((*(pc->vp)), rc)  // use vp position

		CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
		//cell_n(cell, ni, nj, nk, nc)
		//if (nc >= cell.nCells) nc = cell.nCells - 1;
		//cell.cell[nc].attach(pc);
		cell.bcell.m[ni].m[nj].m[nk].attach(pc);
	}
}

void InitMD(MMOL_VMD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, bool bStorageChain) {
	//mmcell.init_big_small_mm();

	::bVirial = true;

	int nc = 0, nm = 0;
	double Ek = 0;

	float vcluster = 1.5 * 1.5 * 1.5, rcut = 0;
	bool bRelshipStorageChain = false;
	rcut = ::rcut_LJ + ::size_cluster * 2;
	int ncluster_LJ = int(rcut * rcut * rcut / vcluster * 1.2 + 0.5);
	rcut = ::rcut_Ewd + ::size_cluster * 2;
	int ncluster_SPME = int(rcut * rcut * rcut / vcluster * 1.2 + 0.5);
	int ncluster_ImpSolv = int(::dSolv_max * ::dSolv_max * ::dSolv_max / vcluster * 1.2 + 0.5);

	mmcell.bHingeConstraint = false;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mmcell.mm.m[nm]->setup_cluster_fatom();// be sure the eneutral is checked
		mmcell.mm.m[nm]->setup_mcDyn();
		mmcell.mm.m[nm]->setup_kin_dyn();
		mmcell.mm.m[nm]->init_hinge_force();
		if (mmcell.mm.m[nm]->setup_hinge_constraint()) {
			mmcell.bHingeConstraint = false; show_log("disable hinge constraint", true);
		}
		mmcell.mm.m[nm]->setup_hinge_dihedral();
		if (!bRelshipStorageChain) {
			for (nc = 0; nc < mmcell.mm.m[nm]->nCluster; nc++) {
				mmcell.mm.m[nm]->cluster[nc].LJRship.use_array(ncluster_LJ);
				mmcell.mm.m[nm]->cluster[nc].spmeRship.use_array(ncluster_SPME);
				if (::bImplicitSolvent) mmcell.mm.m[nm]->cluster[nc].impsolv_rship.use_array(ncluster_ImpSolv);
			}
		}
	}
	if (mmcell.sm.m != NULL) {
		for (nm = 0; nm < mmcell.sm.n; nm++) {
			mmcell.sm.m[nm]->c->setup_cluster_fatom();// be sure the eneutral is checked
			mmcell.sm.m[nm]->setup_kin_dyn();
			if (!bRelshipStorageChain) {
				mmcell.sm.m[nm]->c->LJRship.use_array(ncluster_LJ);
				mmcell.sm.m[nm]->c->spmeRship.use_array(ncluster_SPME);
				if (::bImplicitSolvent) mmcell.sm.m[nm]->c->impsolv_rship.use_array(ncluster_ImpSolv);
			}
		}
	}
	if (mmcell.pm.m != NULL) {
		for (nm = 0; nm < mmcell.pm.n; nm++) {
			mmcell.pm.m[nm]->c->setup_cluster_fatom();// be sure the eneutral is checked
			mmcell.pm.m[nm]->setup_dyn();
			if (!bRelshipStorageChain) {
				mmcell.pm.m[nm]->c->LJRship.use_array(ncluster_LJ);
				mmcell.pm.m[nm]->c->spmeRship.use_array(ncluster_SPME);
				if (::bImplicitSolvent) mmcell.pm.m[nm]->c->impsolv_rship.use_array(ncluster_ImpSolv);
			}
		}
	}

	check_atom_rg(mmcell, MAX_THREADS);
	// since molecule mass center position could be shifted, distribute cluster in CMM again 
	//for (nm = 0; nm < mmcell.mm.n; nm++) CalMolInertiaMoment(mmcell.mm.m[nm]); // here we will calculate mass center
	//for (nm = 0; nm < mmcell.sm.n; nm++) mmcell.sm.m[nm].calInertiaMoment();

	extern void cmm_init_distribute_cluster(MMOLECULE &mm, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain);
	extern void cmm_init_distribute_cluster(SMOLECULE** sm, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain);
	extern void cmm_init_distribute_cluster(PMOLECULE** pm, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain);

	bool bCMM_Storage_Chain = bStorageChain;
	bool reset_cell_chain = true;
	if (bCMM_Storage_Chain) {
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			reset_cell_chain = (nm == 0 ? true : false);
			cmm_init_distribute_cluster(*(mmcell.mm.m[nm]), *cmm, reset_cell_chain); // set base cluster at the center cell
		}
		reset_cell_chain = (mmcell.mm.n == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.sm.m, mmcell.sm.n, *cmm, reset_cell_chain); // add simple solid molecule into cmm
		reset_cell_chain = (mmcell.mm.n + mmcell.sm.n == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.pm.m, mmcell.pm.n, *cmm, reset_cell_chain); // add point molecule into cmm

		//cmm_init_subcell<BASIC_CLUSTER>(*cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
	}

	mmcell.init_mdcell_velocity_memory();

	// ************ SETUP the KinDyn variables *******************
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mm = mmcell.mm.m[nm];
		InitVelocity(*mm);
		mm->calClusterSolvationRadius(); // it uses the cluster's mass center position
		//mm->resetDynPars(); // set force & torque to 0
	}

	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m[nm];
		InitVelocity(*sm);
		sm->c->calSolvationRadius();
	}

	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m[nm];
		InitVelocity(*pm);
		pm->c->calSolvationRadius();
	}

	CHAIN<short> *ch_cluster_id = NULL;
	if (::bBrownian) {
		strcpy(::errmsg, "\n"); show_log(::errmsg, true);
		sprintf(::errmsg, "----------------------  Brownian Force  ---------------------"); show_log(::errmsg, true);
		sprintf(::errmsg, "  general parameters :"); show_log(::errmsg, true);
		sprintf(::errmsg, "  ita / radiu = %f", ::ita_Brownian); show_log(::errmsg, true);
		sprintf(::errmsg, "  cd_trans / d^2 = %f", ::cd_drag_trans); show_log(::errmsg, true);
		sprintf(::errmsg, "  cd_rot / d^2 = %f", ::cd_drag_rot); show_log(::errmsg, true);
		show_log("", true);
		sprintf(::errmsg, "Eistein-type diffusion coefficient D = kT / ita of cluster"); show_log(::errmsg, true);
		sprintf(::errmsg, "  For each cluster, ita = ita * radiu"); show_log(::errmsg, true);
		sprintf(::errmsg, "                    cd_trans = cd_trans * d^2"); show_log(::errmsg, true);
		sprintf(::errmsg, "                    cd_rot = cd_rot * d^2"); show_log(::errmsg, true);
		sprintf(::errmsg, ""); show_log(::errmsg, true);

		sprintf(::errmsg, "  CLUSTER UID: ita, cd_translation, cd_rotation"); show_log(::errmsg, true);
		sprintf(::errmsg, "  Brownian force: f = -ita * V - cd_trans/rot * |V| * V"); show_log(::errmsg, true);
		sprintf(::errmsg, "  Brownian force unit: mass * Angs/fs^2"); show_log(::errmsg, true);
		sprintf(::errmsg, "--------------------------------------------------------"); show_log(::errmsg, true);
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			mm = mmcell.mm.m[nm];
			mm->show_cluster_BrownianForcePars(&ch_cluster_id);
		}
		for (nm = 0; nm < mmcell.sm.n; nm++) {
			sm = mmcell.sm.m[nm];
			sm->c->show_BrownianForcePars(&ch_cluster_id);
		}
		for (nm = 0; nm < mmcell.pm.n; nm++) {
			pm = mmcell.pm.m[nm];
			pm->c->show_BrownianForcePars(&ch_cluster_id);
		}
		release_chain<short>(&ch_cluster_id, true);
		sprintf(::errmsg, "------------------------- OVER -------------------------"); show_log(::errmsg, true);
		show_log("", true);
	}

	/*
	show_log("-------  NHC ---------", true);
	sprintf(errmsg, "  dt -- %f, tao -- %f", nc, mmcell.hdyn.nhc->dt, 1 / sqrt(mmcell.hdyn.nhc->inv_taos));
	show_log(errmsg, true);

	for (nm = 0; nm < mmcell.lfpmd.n; nm++) {
		sprintf(errmsg, "LFPMD # %d : dt = %f", nm, mmcell.lfpmd.m[nm]->dt); show_log(errmsg, true);
	}
	show_log("-------- OVER --------", true);
	show_log("", true);
	*/
}

void VMDCELL_CalMolTransKineticEnergy(MDCELL_ThreadJob<MMOL_VMD_CELL> *job, SVECTOR<double, 4> *Ekt) {
	extern void NPT_calKineticEnergy_Cluster(MMOLECULE& mm, bool bTransOnly = false);
	extern void NPT_calKineticEnergy_Mol(MMOLECULE& mm, bool bTransOnly = false);
	extern void NPT_calKineticEnergy(SMOLECULE& sm, bool bTransOnly = false);
	extern void NPT_calKineticEnergy(PMOLECULE& pm);

	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0;
	double Ek = 0, Ekx = 0, Eky = 0, Ekz = 0;
	VECTOR3 P;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		calTransSpatialMomentum(*mm, P); // update velocity of each cluster at its mass-center
		if (::method_Virial == CLUSTER_VIRIAL) {
			NPT_calKineticEnergy_Cluster(*mm, true); // translation kinetic energy of all clusters
		}
		else if (::method_Virial == MOLECULE_VIRIAL) {
			NPT_calKineticEnergy_Mol(*mm, true); // translation kinetic energy of the mass center of the molecule
		}
		Ek += mm->Ekt;
		Ekx += mm->Ektx;
		Eky += mm->Ekty;
		Ekz += mm->Ektz;
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		NPT_calKineticEnergy(*sm, true);
		Ek += sm->Ekt;
		Ekx += sm->Ektx;
		Eky += sm->Ekty;
		Ekz += sm->Ektz;
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
		NPT_calKineticEnergy(*pm);
		Ek += pm->Ekt;
		Ekx += pm->Ektx;
		Eky += pm->Ekty;
		Ekz += pm->Ektz;
	}
	Ekt->v[0] = Ekx;
	Ekt->v[1] = Eky;
	Ekt->v[2] = Ekz;
	Ekt->v[3] = Ek;
}

void log_dipole(MMOL_VMD_CELL *mdcell) {
	char msg[256] = "\0", buffer[256] = "\0";
	VECTOR3 mu, mut;
	BASIC_CLUSTER *pc;
	int im, ia;
	for (im = 0; im < mdcell->sm.n; im++) {
		pc = mdcell->sm.m[im]->c; V3zero(mu) V3zero(mut)
		sprintf(msg, "sm %d : ", im);
		for (ia = 0; ia < pc->nAtoms; ia++) {
			mu += pc->atom[ia].mu;
			mut.v[0] += pc->atom[ia].c * pc->atom[ia].r0.v[0];
			mut.v[1] += pc->atom[ia].c * pc->atom[ia].r0.v[1];
			mut.v[2] += pc->atom[ia].c * pc->atom[ia].r0.v[2];
			sprintf(buffer, "%f, ", pc->atom[ia].mu.abs());
			strcat(msg, buffer);
		}
		mut += mu;
		sprintf(buffer, "induced dipole %f eA, total dipole %f eA", mu.abs(), mut.abs());
		strcat(msg, buffer);
		show_log(msg, true);
	}
	for (im = 0; im < mdcell->pm.n; im++) {
		pc = mdcell->pm.m[im]->c; V3zero(mu) V3zero(mut)
		sprintf(msg, "pm %d : ", im);
		for (ia = 0; ia < pc->nAtoms; ia++) {
			mu += pc->atom[ia].mu;
			mut.v[0] += pc->atom[ia].c * pc->atom[ia].r0.v[0];
			mut.v[1] += pc->atom[ia].c * pc->atom[ia].r0.v[1];
			mut.v[2] += pc->atom[ia].c * pc->atom[ia].r0.v[2];
			sprintf(buffer, "%f, ", pc->atom[ia].mu.abs());
			strcat(msg, buffer);
		}
		mut += mu;
		sprintf(buffer, "induced dipole %f eA, total dipole %f eA", mu.abs(), mut.abs());
		strcat(msg, buffer);
		show_log(msg, true);
	}
	show_log("", true);
}


void MDCELL_mass_center(MMOL_VMD_CELL &mmcell) {
	VECTOR3 mc;
	float M = 0;
	MMOLECULE *m = NULL;
	int nm = 0, i = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m[nm];
		m->calMassCenter(true);
		for (i = 0; i < 3; i++) mc.v[i] += m->mDynKin.M * m->mDynKin.cm.v[i];
		M += m->mDynKin.M;
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m[nm];
		for (i = 0; i < 3; i++) mc.v[i] += sm->c->M * sm->r.v[i];
		M += sm->c->M;
	}
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m[nm];
		for (i = 0; i < 3; i++) mc.v[i] += pm->c->M * pm->r.v[i];
		M += pm->c->M;
	}
	for (i = 0; i < 3; i++) mc.v[i] /= M;
	V32V3(mc, mmcell.rc)
}

void VMDCELL_CalCellInertiaTensor(MDCELL_ThreadJob<MMOL_VMD_CELL> *job, MATRIX<3> *mI) {
	MMOL_VMD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	MATRIX<3> Iw;
	VECTOR3 Rcm, rc;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m[imol];
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			Rcm.v[0] = pc->cm.v[0] + pc->dr_cmm.v[0];
			Rcm.v[1] = pc->cm.v[1] + pc->dr_cmm.v[1];
			Rcm.v[2] = pc->cm.v[2] + pc->dr_cmm.v[2];

			VECT3(mdcell->rc, Rcm, rc)

			Iw.m[0][0] += pc->M * (rc.v[1] * rc.v[1] + rc.v[2] * rc.v[2]);
			Iw.m[1][1] += pc->M * (rc.v[0] * rc.v[0] + rc.v[2] * rc.v[2]);
			Iw.m[2][2] += pc->M * (rc.v[0] * rc.v[0] + rc.v[1] * rc.v[1]);
			Iw.m[0][1] -= pc->M * rc.v[0] * rc.v[1];
			Iw.m[0][2] -= pc->M * rc.v[0] * rc.v[2];
			Iw.m[1][2] -= pc->M * rc.v[1] * rc.v[2];
		}
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m[imol];
		M = sm->c->M;

		VECT3(mdcell->rc, sm->r, rc)

		Iw.m[0][0] += M * (rc.v[1] * rc.v[1] + rc.v[2] * rc.v[2]);
		Iw.m[1][1] += M * (rc.v[0] * rc.v[0] + rc.v[2] * rc.v[2]);
		Iw.m[2][2] += M * (rc.v[0] * rc.v[0] + rc.v[1] * rc.v[1]);
		Iw.m[0][1] -= M * rc.v[0] * rc.v[1];
		Iw.m[0][2] -= M * rc.v[0] * rc.v[2];
		Iw.m[1][2] -= M * rc.v[1] * rc.v[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m[imol];
		M = pm->c->M;

		VECT3(mdcell->rc, pm->r, rc)

		Iw.m[0][0] += M * (rc.v[1] * rc.v[1] + rc.v[2] * rc.v[2]);
		Iw.m[1][1] += M * (rc.v[0] * rc.v[0] + rc.v[2] * rc.v[2]);
		Iw.m[2][2] += M * (rc.v[0] * rc.v[0] + rc.v[1] * rc.v[1]);
		Iw.m[0][1] -= M * rc.v[0] * rc.v[1];
		Iw.m[0][2] -= M * rc.v[0] * rc.v[2];
		Iw.m[1][2] -= M * rc.v[1] * rc.v[2];
	}
	mI->m[0][0] = Iw.m[0][0];
	mI->m[1][1] = Iw.m[1][1];
	mI->m[2][2] = Iw.m[2][2];
	mI->m[0][1] = Iw.m[0][1];
	mI->m[0][2] = Iw.m[0][2];
	mI->m[1][2] = Iw.m[1][2];
	mI->m[1][0] = Iw.m[0][1];
	mI->m[2][0] = Iw.m[0][2];
	mI->m[2][1] = Iw.m[1][2];
}

void NVT_MD(MMOL_VMD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, bool bPolarize, NVT_VMDVar<MMOL_VMD_CELL, BASIC_CLUSTER> &mdvar, long nloops, int iloop_start) {
	using namespace _spme_;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
extern void cal_dipole(MMOL_VMD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
#endif

	int i = 0;

	long nloop = 0;
	int nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0, U_LJ = 0, U_estat = 0;

	double Uext = 0; // external electric field

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	double djob[MAX_THREADS];
	MATRIX<3> mjob[MAX_THREADS];
	bool status_job[MAX_THREADS];
	
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];

	
	// make sure the spatial momentum of the whole cell is zero

	int nCheckAngle = 50, iCheckAngle = 0;
	//int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;
	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;
	bool bCheckClusterInCMM = true;
	bool Overlap = false;
	int nConvg = 0;
	bool bSpeedVerlet = true;

	// the amount of molecules could be changed, and the job assignment needs to be refreshed
	mdvar.init_job(mmcell, cmm); 
	mdvar.init_SPME_SYS(mmcell, bPolarize); // charges in SPME needs to be updated also

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	long LOOPS = 0;
	int save_indx = 0;
	MD_SAVE *mdsave = mmcell.msave;

	int nMDThreads = MAX_THREADS;

	int ncs_mmcell = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) ncs_mmcell += mmcell.mm.m[nm]->nCluster;
	ncs_mmcell += mmcell.sm.n + mmcell.pm.n;

	while (nloop < nloops) {
		if (COMMAND_MD_STOP) break; // stop with given command

		iCheckCMM = 0;

		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			// here we have to check the periodity of the molecule position
			//molecule_check_periodic_cell(mmcell);
			MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_molecule_PeriodicCell_check), mdvar.job_mdcell, nMDThreads);
		}

		check_atom_rg(mmcell, nMDThreads);

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&MDCELL_CalInertiaMoment), mdvar.job_mdcell, nMDThreads);
		
		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&MDCELL_InitClusterForce), mdvar.job_mdcell, nMDThreads);

	// re-distribute / check clusters in CMM, use vp of each cluster only
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		if (mdvar.bStorageChain) {
			if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
			else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
		}
		else {
			if (nloop == 0 || bCheckClusterInCMM) {
				if (MAX_THREADS == 1) {
					mdvar.job_mdcellCMM[0].cmm = cmm;
					cmm->reset_cell_chain();
					VMMCELL_CMM_cluster_allocate(mdvar.job_mdcellCMM, status_job);
				}
				else {
					for (i = 0; i < MAX_THREADS; i++) mdvar.job_mdcellCMM[i].cmm = mdvar.job_cmm + i; // set the JOB-related CMM is the job-related CMM
					MTOperate1<MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, bool>((void*)(&VMMCELL_CMM_cluster_allocate), mdvar.job_mdcellCMM, nMDThreads, status_job);
					for (i = 0; i < nMDThreads; i++) {
						if (!status_job[i]) {
							sprintf(errmsg, "loop %d thread %d: buffer of CMM cell leaks when allocating cluster", nloop, i);
							show_log(errmsg, true);
						}
					}
					if (!combine_AtomAllocatedCMM_array<BASIC_CLUSTER>(mdvar.job_cmm, nMDThreads, cmm, true, true)) {
						sprintf(errmsg, "loop %d : buffer of CMM cell leaks when combining allocated cluster", nloop);
						show_log(errmsg, true);
					}
				}
				for (i = 0; i < MAX_THREADS; i++) mdvar.job_mdcellCMM[i].cmm = cmm; // reset the JOB-related CMM is the main CMM
			}

			if (bCheckClusterInCMM) {
			{
				nc = 0;
				for (i = 0; i < cmm->acell.n; i++) {
					nc += cmm->acell.m[i]->length();
					//sprintf(errmsg, "%d %d", i, cmm->acell.m[i]->length()); show_log(errmsg, true);
				}
				if (nc != ncs_mmcell) {
					sprintf(errmsg, "%d : cluster in CMM -- %d, real -- %d", nloop, nc, ncs_mmcell); show_log(errmsg, true);
				}
			}
			}
			
		}
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
			if (mdvar.bStorageChain) {
				ParallelCheckRelshipInCells(*cmm, 0, cmm->acell.n - 1, MAX_THREADS); // setup relationship for each atom with CMM
			}
			else {
				MTOperate<MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >((void*)(&VMMCELL_CMM_CheckRelship), mdvar.job_mdcellCMM, nMDThreads);
			}
		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
			show_infor(errmsg);
		#endif
		}

		// important : calculate LJ interaction first, here cluster force & torque are reset
		
		LJ_Interact(&mmcell, cmm, MAX_THREADS, *(mdvar.LJ_var), *(mdvar.LJ_res));
		U_LJ = mdvar.LJ_res->U;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!::local_interact && !mmcell.eNeutral) {
			// calculat the total dipole of the whole cell
			cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));

			mdvar.spme_var->bEfieldOnly = false;
			// important : calculate Coulomb secondly, here cluster force & torque are accumulated
			if (::bPolarize) {
				if (!Polarize_SPME_Interact(&mmcell, cmm, *(mdvar.vspme_mu), MAX_THREADS, true, *(mdvar.spme_var), *(mdvar.spme_res))) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
				}
			}
			else SPME_Interact(&mmcell, cmm, *(mdvar.vspme_q), MAX_THREADS, *(mdvar.spme_var), *(mdvar.spme_res));

			U_estat = mdvar.spme_res->U;
		}
#endif
		
		Ep = U_LJ + U_estat;
		if (::bVirial) {
			mdvar.virial->v[0] = mdvar.LJ_res->STxx + mdvar.spme_res->STxx;
			mdvar.virial->v[1] = mdvar.LJ_res->STyy + mdvar.spme_res->STyy;
			mdvar.virial->v[2] = mdvar.LJ_res->STzz + mdvar.spme_res->STzz;
			mdvar.virial->v[3] = mdvar.LJ_res->STtr + mdvar.spme_res->STtr;

			if (!mdvar.virial_save->one_more()) break;
			mdvar.bSaveVirial = mdvar.virial_save->bsave(); 
		}
		
		//mdvar.virial.v[0] = 0; mdvar.virial.v[1] = 0; mdvar.virial.v[2] = 0; mdvar.virial.v[3] = 0; Ep = 0;

		if (mmcell.mm.n > 0) {
			MTOperate2< MDCELL_ThreadJob<MMOL_VMD_CELL>, TorsionVar, InteractRes>((void*)(&VMMCELL_MM_Torsion), mdvar.job_mdcell, nMDThreads, mdvar.torsion_var, mdvar.mires);
			Etorsion = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ep += mdvar.mires[i].U; Etorsion += mdvar.mires[i].U;
				if (mdvar.torsion_var[i].bVirial) {
					mdvar.virial->v[0] += mdvar.mires[i].STxx;
					mdvar.virial->v[1] += mdvar.mires[i].STyy;
					mdvar.virial->v[2] += mdvar.mires[i].STzz;
					mdvar.virial->v[3] += mdvar.mires[i].STtr;
				}
			}
			//sprintf(errmsg, "%d: Torsion %f kT", nloop, Etorsion); show_log(errmsg, true);
		}

		if (::bEext) {
			if (mmcell.mm.m != NULL) {
				MOperate3<MMOLECULE*, float, double>((void*)(&VMM_ExternalEField), mmcell.mm.m, mmcell.mm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.sm.m != NULL) {
				MOperate3<SMOLECULE*, float, double>((void*)(&VSM_ExternalEField), mmcell.sm.m, mmcell.sm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.pm.m != NULL) {
				MOperate3<PMOLECULE*, float, double>((void*)(&VPM_ExternalEField), mmcell.pm.m, mmcell.pm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
		}


		//if (mmcell.mm.m != NULL) MOperate<MMOLECULE*>(VMM_TOTAL_FORCE, mmcell.mm.m, mmcell.mm.n, MAX_THREADS);
		//if (mmcell.sm.m != NULL) MOperate<SMOLECULE*>(VSM_TOTAL_FORCE, mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
		//if (mmcell.pm.m != NULL) MOperate<PMOLECULE*>(VPM_TOTAL_FORCE, mmcell.pm.m, mmcell.pm.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_ClusterForce), mdvar.job_mdcell, nMDThreads);

		if (::bVirial && mdvar.bSaveVirial && (mmcell.mm.n > 0 || mmcell.sm.n > 0)) {
			for (i = 0; i < nMDThreads; i++) {
				sv4job[i].v[0] = 0; sv4job[i].v[1] = 0; sv4job[i].v[2] = 0; sv4job[i].v[3] = 0;
			}
			switch (::method_Virial) {
			case CLUSTER_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_VMD_CELL>, SVECTOR<double, 4> >((void*)(&VMDCELL_VirialCorrect_RigidCluster), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			case MOLECULE_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_VMD_CELL>, SVECTOR<double, 4> >((void*)(&VMDCELL_VirialCorrect_Molecule), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			}
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial->v[0] -= sv4job[i].v[0];
				mdvar.virial->v[1] -= sv4job[i].v[1];
				mdvar.virial->v[2] -= sv4job[i].v[2];
				mdvar.virial->v[3] -= sv4job[i].v[0] + sv4job[i].v[1] + sv4job[i].v[2];
			}
		}
		
		if (::bVirial) {
			if (mdvar.virial_save->bsave()) {
				(*mdvar.virial_save->out)<<nloop<<"  "<<mdvar.virial->v[0] * ::unit_Ek_kT<<"  "<<mdvar.virial->v[1] * ::unit_Ek_kT<<"  "<<mdvar.virial->v[2] * ::unit_Ek_kT<<"  "<<mdvar.virial->v[3] * ::unit_Ek_kT;
			}
		}

		// IMPORTANT: the forces from here on needs to calculate the induced torque explicitly
		// randome force and ImplicitSolvent force, and the induced torque, will be calculated explicitly
		// so we do not need to call MDCELL_ClusterForce.

		if (bImplicitSolvent) {
			// neighbors of cluster for implicit solvation 
			//if (bCheckClusterInCMM) CellOperate<BASIC_CLUSTER>((void*)(&CMMCheckImpSolvNeighbors_Cell), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			// solvation energy
			//CellOperate<BASIC_CLUSTER>((void*)(&CMM_ImplicitSolvationForce), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			//for (nm = 0; nm < mmcell.mm.n; nm++) {
			//	mm = mmcell.mm.m + nm;
			//	for (nc = 0; nc < mm->nCluster; nc++) Ep += mm->cluster[nc].E_ImplicitSolvation;
			//}
			MTOperate<MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >((void*)(&VMMCELL_CMM_CheckImpSolvRelship), mdvar.job_mdcellCMM, nMDThreads);
			MTOperate1< MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, double>((void*)(&VMMCELL_cluster_impsolv_interact), mdvar.job_mdcellCMM, nMDThreads, djob);
			for (i = 0; i < nMDThreads; i++) Ep += djob[i];
		}

		if (bRandomForce) gen_random_force(mmcell, *(mdvar.rand0), *(mdvar.rand1));

		/*
		if (scale_force) { // scaled in LJ interaction calculation
			CellOperate<BASIC_CLUSTER>((void*)(&SCALE_FORCE_CMM), *cmm, 0, cmm->acell.n, MAX_THREADS);
		}
		*/

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		// self-consistent velocity and acceleration, Nose-Hoover Chain
		bSpeedVerlet = (nloop == 0 ? false : true);
		if (!LeapFrogVerlet_SelfConsistent_SpeedAccel(&mmcell, mdvar.job_mdcell, bSpeedVerlet, nMDThreads)) {
			sprintf(errmsg, "failure to get consistent velocity and acceleration at loop %d", nloop);
			show_log(errmsg, true);
		}

		Ep_t = Ep; Ek_t = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) Ek_t += mmcell.mm.m[nm]->Ek;
		for (nm = 0; nm < mmcell.sm.n; nm++) Ek_t += mmcell.sm.m[nm]->Ek;
		for (nm = 0; nm < mmcell.pm.n; nm++) Ek_t += mmcell.pm.m[nm]->Ek;

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop + iloop_start, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);
		mmcell.MDSave(nloop + iloop_start, (float)Ep_t, false);

		//if (::bPolarize && mdsave->mKinSave.bsave()) log_dipole(&mmcell); // write the dipole of cluster into log file

		if (::bVirial && mdvar.bSaveVirial && ::method_Virial == CLUSTER_VIRIAL) {
			if (mmcell.mm.n > 0 && mmcell.bHingeConstraint) {
				MTOperate2< MDCELL_ThreadJob<MMOL_VMD_CELL>, HingeConstraintVar, InteractRes>((void*)(&MDCELL_ClusterVirialCorrect_HingeConstraint), mdvar.job_mdcell, nMDThreads, mdvar.hinge_constraint_var, mdvar.mires);
				for (i = 0; i < nMDThreads; i++) {
					mdvar.virial->v[0] += mdvar.mires[i].STxx;
					mdvar.virial->v[1] += mdvar.mires[i].STyy;
					mdvar.virial->v[2] += mdvar.mires[i].STzz;
					mdvar.virial->v[3] += mdvar.mires[i].STtr;
				}
			}

			MTOperate1<MDCELL_ThreadJob<MMOL_VMD_CELL>, SVECTOR<double, 4> >((void*)(&VMDCELL_CalMolTransKineticEnergy), mdvar.job_mdcell, nMDThreads, sv4job);
			Ek_t = 0;
			for (i = 0; i < nMDThreads; i++) Ek_t += sv4job[i].v[3];
			*(mdvar.Pint) = (mdvar.virial->v[3] + Ek_t * 2 / unit_Ek_kT) / (::fUnit_atm_u * mmcell.h[0] * mmcell.h[1] * mmcell.h[2] * 8 * 3);
			if (mdvar.bSaveVirial) {
				(*mdvar.virial_save->out)<<"          "<<mdvar.virial->v[3] * ::unit_Ek_kT<<"  "<<Ek_t<<"  "<<mdvar.Pint<<endl;
			}
		}

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_Verlet), mdvar.job_mdcell, nMDThreads);
		// accept current NHC vars and prepair for next step
		mmcell.hdyn.nhc->accept_nextstep();

		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			if (mmcell.mm.n > 2) {
				MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_CheckRotAngle), mdvar.job_mdcell, nMDThreads);
			}
			else {
				for (nm = 0; nm < mmcell.mm.n; nm++) {
					mm = mmcell.mm.m[nm];
					check_angle(*mm, *mmcell.lfmd.m[nm]);
				}
			}
		}
		iCheckAngle++; if (iCheckAngle >= nCheckAngle) iCheckAngle = 0;

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mmcell.CalibrateTorsionAngle();

		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "checking procedure takes %d ms", time.glance());
			show_infor(errmsg);
		#endif

		if (::bVirial) mdvar.virial_save->save_over();

		nloop++;
	}
	/*
	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
	if (::bVirial) mdvar.virial_save.reset();
	*/
}



#ifdef _DISABLE_
void GC_MD(MMOL_VMD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, bool bPolarize, GC_VMDVar<MMOL_VMD_CELL, BASIC_CLUSTER> &mdvar, long nloops, int iloop_start) {
	using namespace _spme_;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
extern void cal_dipole(MMOL_VMD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
#endif

	int i = 0;

	long nloop = 0;
	int nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0, U_LJ = 0, U_estat = 0;

	double Uext = 0; // external electric field

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	double djob[MAX_THREADS];
	//MATRIX<3> mjob[MAX_THREADS];
	bool status_job[MAX_THREADS];
	
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];

	//InteractRes tres;

	
	// make sure the spatial momentum of the whole cell is zero

	int nCheckAngle = 50, iCheckAngle = 0;
	//int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;
	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;
	bool bCheckClusterInCMM = true;
	bool Overlap = false;
	int nConvg = 0;
	bool bSpeedVerlet = true;

	// the amount of molecules could be changed, and the job assignment needs to be refreshed
	mdvar.init_job(mmcell, cmm); 
	mdvar.init_SPME_SYS(mmcell, bPolarize); // charges in SPME needs to be updated also

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	long LOOPS = 0;
	int save_indx = 0;
	MD_SAVE *mdsave = mmcell.msave;

	int nMDThreads = MAX_THREADS;

	int ncs_mmcell = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) ncs_mmcell += mmcell.mm.m[nm]->nCluster;
	ncs_mmcell += mmcell.sm.n + mmcell.pm.n;

	while (nloop < nloops) {
		if (COMMAND_MD_STOP) break; // stop with given command

		iCheckCMM = 0;

		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			// here we have to check the periodity of the molecule position
			//molecule_check_periodic_cell(mmcell);
			MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_molecule_PeriodicCell_check), mdvar.job_mdcell, nMDThreads);
		}

		check_atom_rg(mmcell, nMDThreads);

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&MDCELL_CalInertiaMoment), mdvar.job_mdcell, nMDThreads);
		
		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&MDCELL_InitClusterForce), mdvar.job_mdcell, nMDThreads);

	// re-distribute / check clusters in CMM, use vp of each cluster only
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		if (mdvar.bStorageChain) {
			if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
			else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
		}
		else {
			if (nloop == 0 || bCheckClusterInCMM) {
				if (MAX_THREADS == 1) {
					mdvar.job_mdcellCMM[0].cmm = cmm;
					cmm->reset_cell_chain();
					VMMCELL_CMM_cluster_allocate(mdvar.job_mdcellCMM, status_job);
				}
				else {
					for (i = 0; i < MAX_THREADS; i++) mdvar.job_mdcellCMM[i].cmm = mdvar.job_cmm + i; // set the JOB-related CMM is the job-related CMM
					MTOperate1<MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, bool>((void*)(&VMMCELL_CMM_cluster_allocate), mdvar.job_mdcellCMM, nMDThreads, status_job);
					for (i = 0; i < nMDThreads; i++) {
						if (!status_job[i]) {
							sprintf(errmsg, "loop %d thread %d: buffer of CMM cell leaks when allocating cluster", nloop, i);
							show_log(errmsg, true);
						}
					}
					if (!combine_AtomAllocatedCMM_array<BASIC_CLUSTER>(mdvar.job_cmm, nMDThreads, cmm, true, true)) {
						sprintf(errmsg, "loop %d : buffer of CMM cell leaks when combining allocated cluster", nloop);
						show_log(errmsg, true);
					}
				}
				for (i = 0; i < MAX_THREADS; i++) mdvar.job_mdcellCMM[i].cmm = cmm; // reset the JOB-related CMM is the main CMM
			}

			if (bCheckClusterInCMM) {
			{
				nc = 0;
				for (i = 0; i < cmm->acell.n; i++) {
					nc += cmm->acell.m[i]->length();
					//sprintf(errmsg, "%d %d", i, cmm->acell.m[i]->length()); show_log(errmsg, true);
				}
				if (nc != ncs_mmcell) {
					sprintf(errmsg, "%d : cluster in CMM -- %d, real -- %d", nloop, nc, ncs_mmcell); show_log(errmsg, true);
				}
			}
			}
			
		}
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
			if (mdvar.bStorageChain) {
				ParallelCheckRelshipInCells(*cmm, 0, cmm->acell.n - 1, MAX_THREADS); // setup relationship for each atom with CMM
			}
			else {
				MTOperate<MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >((void*)(&VMMCELL_CMM_CheckRelship), mdvar.job_mdcellCMM, nMDThreads);
			}
		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
			show_infor(errmsg);
		#endif
		}

		// important : calculate LJ interaction first, here cluster force & torque are reset
		
		LJ_Interact(&mmcell, cmm, MAX_THREADS, *(mdvar.LJ_var), *(mdvar.LJ_res));
		GC_LJ_Interact_Correct(&mmcell, cmm, *(mdvar.gc_InteractVar), *(mdvar.gc_LJRes));
		*(mdvar.LJ_res) += *((InteractRes*)(&mdvar.gc_LJRes));
		U_LJ = mdvar.LJ_res->U;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!::local_interact && !mmcell.eNeutral) {
			// calculat the total dipole of the whole cell
			cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));

			mdvar.spme_var->bEfieldOnly = false;
			// important : calculate Coulomb secondly, here cluster force & torque are accumulated
			if (::bPolarize) {
				if (!Polarize_SPME_Interact(&mmcell, cmm, *(mdvar.vspme_mu), MAX_THREADS, true, *(mdvar.spme_var), *(mdvar.spme_res))) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
				}
				GC_Electrostatic_Interact_Correct(&mmcell, cmm, *(mdvar.gc_EVar), *(mdvar.gc_ERes));
				*(mdvar.spme_res) += *((InteractRes*)(&mdvar.gc_ERes));
			}
			else {
				SPME_Interact(&mmcell, cmm, *(mdvar.vspme_q), MAX_THREADS, *(mdvar.spme_var), *(mdvar.spme_res));
				GC_Electrostatic_Interact_Correct(&mmcell, cmm, *(mdvar.gc_EVar), *(mdvar.gc_ERes));
				*(mdvar.spme_res) += *((InteractRes*)(&mdvar.gc_ERes));
			}


			U_estat = mdvar.spme_res->U;
		}
#endif
		
		Ep = U_LJ + U_estat;
		if (::bVirial) {
			mdvar.virial->v[0] = mdvar.LJ_res->STxx + mdvar.spme_res->STxx;
			mdvar.virial->v[1] = mdvar.LJ_res->STyy + mdvar.spme_res->STyy;
			mdvar.virial->v[2] = mdvar.LJ_res->STzz + mdvar.spme_res->STzz;
			mdvar.virial->v[3] = mdvar.LJ_res->STtr + mdvar.spme_res->STtr;

			if (!mdvar.virial_save->one_more()) break;
			mdvar.bSaveVirial = mdvar.virial_save->bsave(); 
		}
		
		//mdvar.virial.v[0] = 0; mdvar.virial.v[1] = 0; mdvar.virial.v[2] = 0; mdvar.virial.v[3] = 0; Ep = 0;

		if (mmcell.mm.n > 0) {
			MTOperate2< MDCELL_ThreadJob<MMOL_VMD_CELL>, TorsionVar, InteractRes>((void*)(&VMMCELL_MM_Torsion), mdvar.job_mdcell, nMDThreads, mdvar.torsion_var, mdvar.mires);
			Etorsion = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ep += mdvar.mires[i].U; Etorsion += mdvar.mires[i].U;
				if (mdvar.torsion_var[i].bVirial) {
					mdvar.virial->v[0] += mdvar.mires[i].STxx;
					mdvar.virial->v[1] += mdvar.mires[i].STyy;
					mdvar.virial->v[2] += mdvar.mires[i].STzz;
					mdvar.virial->v[3] += mdvar.mires[i].STtr;
				}
			}
			//sprintf(errmsg, "%d: Torsion %f kT", nloop, Etorsion); show_log(errmsg, true);
		}

		if (::bEext) {
			if (mmcell.mm.m != NULL) {
				MOperate3<MMOLECULE*, float, double>((void*)(&VMM_ExternalEField), mmcell.mm.m, mmcell.mm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.sm.m != NULL) {
				MOperate3<SMOLECULE*, float, double>((void*)(&VSM_ExternalEField), mmcell.sm.m, mmcell.sm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.pm.m != NULL) {
				MOperate3<PMOLECULE*, float, double>((void*)(&VPM_ExternalEField), mmcell.pm.m, mmcell.pm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
		}


		//if (mmcell.mm.m != NULL) MOperate<MMOLECULE*>(VMM_TOTAL_FORCE, mmcell.mm.m, mmcell.mm.n, MAX_THREADS);
		//if (mmcell.sm.m != NULL) MOperate<SMOLECULE*>(VSM_TOTAL_FORCE, mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
		//if (mmcell.pm.m != NULL) MOperate<PMOLECULE*>(VPM_TOTAL_FORCE, mmcell.pm.m, mmcell.pm.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_ClusterForce), mdvar.job_mdcell, nMDThreads);

		if (::bVirial && mdvar.bSaveVirial && (mmcell.mm.n > 0 || mmcell.sm.n > 0)) {
			for (i = 0; i < nMDThreads; i++) {
				sv4job[i].v[0] = 0; sv4job[i].v[1] = 0; sv4job[i].v[2] = 0; sv4job[i].v[3] = 0;
			}
			switch (::method_Virial) {
			case CLUSTER_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_VMD_CELL>, SVECTOR<double, 4> >((void*)(&VMDCELL_VirialCorrect_RigidCluster), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			case MOLECULE_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_VMD_CELL>, SVECTOR<double, 4> >((void*)(&VMDCELL_VirialCorrect_Molecule), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			}
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial->v[0] -= sv4job[i].v[0];
				mdvar.virial->v[1] -= sv4job[i].v[1];
				mdvar.virial->v[2] -= sv4job[i].v[2];
				mdvar.virial->v[3] -= sv4job[i].v[0] + sv4job[i].v[1] + sv4job[i].v[2];
			}
		}
		
		if (::bVirial) {
			if (mdvar.virial_save->bsave()) {
				(*mdvar.virial_save->out)<<nloop<<"  "<<mdvar.virial->v[0] * ::unit_Ek_kT<<"  "<<mdvar.virial->v[1] * ::unit_Ek_kT<<"  "<<mdvar.virial->v[2] * ::unit_Ek_kT<<"  "<<mdvar.virial->v[3] * ::unit_Ek_kT;
			}
		}

		// IMPORTANT: the forces from here on needs to calculate the induced torque explicitly
		// randome force and ImplicitSolvent force, and the induced torque, will be calculated explicitly
		// so we do not need to call MDCELL_ClusterForce.

		if (bImplicitSolvent) {
			// neighbors of cluster for implicit solvation 
			//if (bCheckClusterInCMM) CellOperate<BASIC_CLUSTER>((void*)(&CMMCheckImpSolvNeighbors_Cell), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			// solvation energy
			//CellOperate<BASIC_CLUSTER>((void*)(&CMM_ImplicitSolvationForce), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			//for (nm = 0; nm < mmcell.mm.n; nm++) {
			//	mm = mmcell.mm.m + nm;
			//	for (nc = 0; nc < mm->nCluster; nc++) Ep += mm->cluster[nc].E_ImplicitSolvation;
			//}
			MTOperate<MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >((void*)(&VMMCELL_CMM_CheckImpSolvRelship), mdvar.job_mdcellCMM, nMDThreads);
			MTOperate1< MDCELL_CMM_ThreadJob<MMOL_VMD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, double>((void*)(&VMMCELL_cluster_impsolv_interact), mdvar.job_mdcellCMM, nMDThreads, djob);
			for (i = 0; i < nMDThreads; i++) Ep += djob[i];
		}

		if (bRandomForce) gen_random_force(mmcell, *(mdvar.rand0), *(mdvar.rand1));

		/*
		if (scale_force) { // scaled in LJ interaction calculation
			CellOperate<BASIC_CLUSTER>((void*)(&SCALE_FORCE_CMM), *cmm, 0, cmm->acell.n, MAX_THREADS);
		}
		*/

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		// self-consistent velocity and acceleration, Nose-Hoover Chain
		bSpeedVerlet = (nloop == 0 ? false : true);
		if (!LeapFrogVerlet_SelfConsistent_SpeedAccel(&mmcell, mdvar.job_mdcell, bSpeedVerlet, nMDThreads)) {
			sprintf(errmsg, "failure to get consistent velocity and acceleration at loop %d", nloop);
			show_log(errmsg, true);
		}

		Ep_t = Ep; Ek_t = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) Ek_t += mmcell.mm.m[nm]->Ek;
		for (nm = 0; nm < mmcell.sm.n; nm++) Ek_t += mmcell.sm.m[nm]->Ek;
		for (nm = 0; nm < mmcell.pm.n; nm++) Ek_t += mmcell.pm.m[nm]->Ek;

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop + iloop_start, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);
		mmcell.MDSave(nloop + iloop_start, (float)Ep_t, false);

		//if (::bPolarize && mdsave->mKinSave.bsave()) log_dipole(&mmcell); // write the dipole of cluster into log file

		if (::bVirial && mdvar.bSaveVirial && ::method_Virial == CLUSTER_VIRIAL) {
			if (mmcell.mm.n > 0 && mmcell.bHingeConstraint) {
				MTOperate2< MDCELL_ThreadJob<MMOL_VMD_CELL>, HingeConstraintVar, InteractRes>((void*)(&MDCELL_ClusterVirialCorrect_HingeConstraint), mdvar.job_mdcell, nMDThreads, mdvar.hinge_constraint_var, mdvar.mires);
				for (i = 0; i < nMDThreads; i++) {
					mdvar.virial->v[0] += mdvar.mires[i].STxx;
					mdvar.virial->v[1] += mdvar.mires[i].STyy;
					mdvar.virial->v[2] += mdvar.mires[i].STzz;
					mdvar.virial->v[3] += mdvar.mires[i].STtr;
				}
			}

			MTOperate1<MDCELL_ThreadJob<MMOL_VMD_CELL>, SVECTOR<double, 4> >((void*)(&VMDCELL_CalMolTransKineticEnergy), mdvar.job_mdcell, nMDThreads, sv4job);
			Ek_t = 0;
			for (i = 0; i < nMDThreads; i++) Ek_t += sv4job[i].v[3];
			*(mdvar.Pint) = (mdvar.virial->v[3] + Ek_t * 2 / unit_Ek_kT) / (::fUnit_atm_u * mmcell.h[0] * mmcell.h[1] * mmcell.h[2] * 8 * 3);
			if (mdvar.bSaveVirial) {
				(*mdvar.virial_save->out)<<"          "<<mdvar.virial->v[3] * ::unit_Ek_kT<<"  "<<Ek_t<<"  "<<mdvar.Pint<<endl;
			}
		}

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_Verlet), mdvar.job_mdcell, nMDThreads);
		// accept current NHC vars and prepair for next step
		mmcell.hdyn.nhc->accept_nextstep();

		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			if (mmcell.mm.n > 2) {
				MTOperate< MDCELL_ThreadJob<MMOL_VMD_CELL> >((void*)(&VMDCELL_CheckRotAngle), mdvar.job_mdcell, nMDThreads);
			}
			else {
				for (nm = 0; nm < mmcell.mm.n; nm++) {
					mm = mmcell.mm.m[nm];
					check_angle(*mm, *mmcell.lfmd.m[nm]);
				}
			}
		}
		iCheckAngle++; if (iCheckAngle >= nCheckAngle) iCheckAngle = 0;

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mmcell.CalibrateTorsionAngle();

		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "checking procedure takes %d ms", time.glance());
			show_infor(errmsg);
		#endif

		if (::bVirial) mdvar.virial_save->save_over();

		nloop++;
	}
	/*
	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
	if (::bVirial) mdvar.virial_save.reset();
	*/
}

#endif

void MD_PROC(MMOL_VMD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, int md_mode) {
	using namespace _spme_;

	NVT_MDVar<MMOL_VMD_CELL, BASIC_CLUSTER> var;
	
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}
	// we need to make sure the mass center of the cell is at zero

	molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell

	// *****************  over  ******************
	show_all_molecules();

	NVT_VMDVar<MMOL_VMD_CELL, BASIC_CLUSTER> mdvar;
	mdvar.UseVar(var);

	bool bStorageChain = false;
	mdvar.bStorageChain = bStorageChain;


#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

	init_BSplineFunc<2>(::bsp, ::indx_bspline);

	mdvar.init_SPME_SYS(mmcell, bPolarize);
	mdvar.init_spme_vars(cmm->xd[0] * cmm->xd[1] * cmm->xd[2]);

extern void cal_dipole(MMOL_VMD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
#endif

	int i = 0;
	
	InitMD(mmcell, cmm, bStorageChain);

	long nloop = 0;
	int nc = 0, nm = 0, ncell = 0;

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	int nMDThreads = MAX_THREADS;
	mdvar.init_job(mmcell, cmm);

	MATRIX<3> mjob[MAX_THREADS];
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];
	// make sure the spatial momentum of the whole cell is zero
	MDCELL_mass_center(mmcell);
	ZeroSpatialMomentum(mdvar.job_mdcell, sv3job, MAX_THREADS);

	mmcell.Init_MD(); //copy the initial torsion angle, speed to the mdpar

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	MD_SAVE *mdsave = mmcell.msave;
	MD_BASE *mdbase = mmcell.md_infor();
	long LOOPS = mdbase->LOOPS;
	long nLoops = 1000;

	bool bSaveVirial = false;
	mdvar.bSaveVirial = bSaveVirial;
	mdvar.init_virial_save(mmcell.msave);

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command
		NVT_MD(mmcell, cmm, bPolarize, mdvar, nLoops, nloop);
		nloop += nLoops;
	}

	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
	if (::bVirial) mdvar.virial_save->reset();
}
