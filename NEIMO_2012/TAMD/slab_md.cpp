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
#include "CMM_2d.h"
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
#include "complex.h"
#include "fftw3.h"
#include "EwaldSum.h"
#include "spme.h"
#include "spme_interact.h"
#include "spme_2d.h"
#include "spme_interact_2d.h"

#include "mdvar.h"
#include "slab_md.h"

extern BSplineFunc<2> bsp;

#include "interface.h"
using namespace _interface_constraint_;
using namespace _EwaldSum_real_;

//************************************************************************************
//    2d-periodical slab is periodical in x-y, while not periodical in z axis
//************************************************************************************

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
ITER_SAVE tsave;
#endif

void cluster_check_2d_periodic_cell(MMOL_MD_CELL &mmcell) {
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 2; i++) {
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
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(sm->r0.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
	}
	*/
}

void molecule_check_2d_periodic_cell(MMOL_MD_CELL &mmcell) {
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		m->calMassCenter(false);
		for (i = 0; i < 2; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(m->mDynKin.cm.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) m->shiftMol(ds);

		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 2; i++) {
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
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 2; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(sm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m + nm;
		for (i = 0; i < 2; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(pm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) pm->shiftMol(ds);
	}

	//cluster_check_periodic_cell(mmcell);
}

void MDCELL_molecule_2d_PeriodicCell_check(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0, imol = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		m = mdcell->mm.m + imol;
		m->calMassCenter(false);
		for (i = 0; i < 2; i++) {
			xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
			x = (float)(m->mDynKin.cm.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0) m->shiftMol(ds);

		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 2; i++) {
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
		sm = mdcell->sm.m + imol; 
		for (i = 0; i < 2; i++) {
			xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
			x = (float)(sm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0) sm->shiftMol(ds);
		/*
		if (sm->r.v[0] > mdcell->h[0] || sm->r.v[0] < -mdcell->h[0] ||
			sm->r.v[1] > mdcell->h[1] || sm->r.v[1] < -mdcell->h[1]){
				sprintf(msg, "strange: sm %d has coordinate out of cell [%f, %f, %f]", sm->r.v[0], sm->r.v[1], sm->r.v[2]);
				show_log(msg, true);
		}
		*/
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		for (i = 0; i < 2; i++) {
			xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
			x = (float)(pm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0) pm->shiftMol(ds);
	}
}

// **********************************************************************************************************************
// calculating the Hoover Force on each cluster
// **********************************************************************************************************************
extern void cal_Hoover_force(MMOLECULE *mm, double ksi);
extern void cal_Hoover_force(SMOLECULE *sm, double ksi);
extern void cal_Hoover_force(PMOLECULE *pm, double ksi);



void MMCELL_ExplicitInteract_2d(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *job, INTERACT_2d *interact) {
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
		mm = job->job->mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			bpc = (BASIC_CLUSTER*)(mm->cluster + ic);
			CLUSTER_LJ_INTERACT(*cmm, bpc, interact->var, tres);
			interact->LJ_res += tres;
			if (!bpc->eNeutral) {
				_spme_2d_::ExplicitEInteract_2d(*cmm, bpc, interact->evar, tres);
				interact->eres += tres;
			}
		}
	}
	for (i = 0; i < job->job->smIndx.n; i++) {
		imol = job->job->smIndx.m[i];
		if (imol >= job->job->mdcell->sm.n) continue;
		bpc = job->job->mdcell->sm.m[imol].c;
		CLUSTER_LJ_INTERACT(*cmm, bpc, interact->var, tres);
		interact->LJ_res += tres;
		if (!bpc->eNeutral) {
			ExplicitEInteract_2d(*cmm, bpc, interact->evar, tres);
			interact->eres += tres;
		}
	}
	for (i = 0; i < job->job->pmIndx.n; i++) {
		imol = job->job->pmIndx.m[i];
		if (imol >= job->job->mdcell->pm.n) continue;
		bpc = job->job->mdcell->pm.m[imol].c;
		CLUSTER_LJ_INTERACT(*cmm, bpc, interact->var, tres);
		interact->LJ_res += tres;
		if (!bpc->eNeutral) {
			ExplicitEInteract_2d(*cmm, bpc, interact->evar, tres);
			interact->eres += tres;
		}
	}
}

void MMCELL_CMM_CheckRelship_2d(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator) {
	int i, imol, nc;
	MMOLECULE *mm = NULL;
	BASIC_CLUSTER *pc = NULL;
	for (i = 0; i < allocator->job->mmIndx.n; i++) {
		imol = allocator->job->mmIndx.m[i];
		if (imol >= allocator->job->mdcell->mm.n) continue;
		mm = allocator->job->mdcell->mm.m + imol;
		for (nc = 0; nc < mm->nCluster; nc++) {
			_cmm_2d_::CheckClusterRelship((BASIC_CLUSTER*)(mm->cluster + nc), *(allocator->cmm)); 
		}
	}
	for (i = 0; i < allocator->job->smIndx.n; i++) {
		imol = allocator->job->smIndx.m[i];
		if (imol >= allocator->job->mdcell->sm.n) continue;
		_cmm_2d_::CheckClusterRelship(allocator->job->mdcell->sm.m[imol].c, *(allocator->cmm)); 
	}
	for (i = 0; i < allocator->job->pmIndx.n; i++) {
		imol = allocator->job->pmIndx.m[i];
		if (imol >= allocator->job->mdcell->pm.n) continue;
		_cmm_2d_::CheckClusterRelship(allocator->job->mdcell->pm.m[imol].c, *(allocator->cmm)); 
	}
}


void MMCELL_CMM_CheckImpSolvRelship_2d(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator) {
	extern float rSolvent, dSolv_max;
	//float r0_solvent = rSolvent;
	CMM_CELL3D<BASIC_CLUSTER> *cmm = allocator->cmm;
	float xd[3] = {cmm->xd[0] / cmm->bcell.nx, cmm->xd[1] / cmm->bcell.ny, cmm->xd[2] / cmm->bcell.nz};
	int dnc[3] = {0, 0, 0};
	float rcut = dSolv_max;
	
	if (::local_relax) {dnc[0] = 1; dnc[1] = 1;}
	else {
		dnc[0] = int(rcut / xd[0] + 0.5f); if (dnc[0] == 0) dnc[0] = 1; if (dnc[0] * xd[0] < rcut) dnc[0] += 1;
		dnc[1] = int(rcut / xd[1] + 0.5f); if (dnc[1] == 0) dnc[1] = 1; if (dnc[1] * xd[1] < rcut) dnc[1] += 1;
		dnc[0] += 1; dnc[1] += 1;
	}

	int i, imol, nc;
	MMOLECULE *mm = NULL;
	BASIC_CLUSTER *pc = NULL;
	for (i = 0; i < allocator->job->mmIndx.n; i++) {
		imol = allocator->job->mmIndx.m[i];
		if (imol >= allocator->job->mdcell->mm.n) continue;
		mm = allocator->job->mdcell->mm.m + imol;
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = (BASIC_CLUSTER*)(mm->cluster + nc);
			_cmm_2d_::ClusterImpSolvNeighborsInCell<BASIC_CLUSTER>(cmm, pc, dnc); 
		}
	}
	for (i = 0; i < allocator->job->smIndx.n; i++) {
		imol = allocator->job->smIndx.m[i];
		if (imol >= allocator->job->mdcell->sm.n) continue;
		pc = allocator->job->mdcell->sm.m[imol].c;
		_cmm_2d_::ClusterImpSolvNeighborsInCell<BASIC_CLUSTER>(cmm, pc, dnc); 
	}
	for (i = 0; i < allocator->job->pmIndx.n; i++) {
		imol = allocator->job->pmIndx.m[i];
		if (imol >= allocator->job->mdcell->pm.n) continue;
		pc = allocator->job->mdcell->pm.m[imol].c;
		_cmm_2d_::ClusterImpSolvNeighborsInCell<BASIC_CLUSTER>(cmm, pc, dnc); 
	}
}

extern void calTransSpatialMomentum(MMOLECULE &mm, VECTOR3 &P);
void MDCELL_SpatialMomentumZ(MDCELL_ThreadJob<MMOL_MD_CELL> *job, double *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
// do we need to check the rotation momentum ?
	VECTOR3 P, v, p;
	//VECTOR3 Pw, w, pw, pw1;
	VECTOR3 Rcm;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;

		calTransSpatialMomentum(*mm, v); // velocity of each cluster at its mass-center
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;

			p.v[0] = pc->M * pc->mcDyn->v.v[0];
			p.v[1] = pc->M * pc->mcDyn->v.v[1];
			p.v[2] = pc->M * pc->mcDyn->v.v[2];

			V3plusV3(p, P, P)
		}
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; M = sm->c->M;

		p.v[0] = M * sm->c->bkm->V0.v[3];
		p.v[1] = M * sm->c->bkm->V0.v[4];
		p.v[2] = M * sm->c->bkm->V0.v[5];

		V3plusV3(p, P, P)
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol; M = pm->c->M;

		p.v[0] = M * pm->v.v[0];
		p.v[1] = M * pm->v.v[1];  
		p.v[2] = M * pm->v.v[2];

		V3plusV3(p, P, P)
	}
	//mp->v[0] = Pw.v[0]; mp->v[1] = Pw.v[1]; mp->v[2] = Pw.v[2];
	//mp->v[3] = P.v[0]; mp->v[4] = P.v[1]; mp->v[5] = P.v[2];
	*mp = P.v[2];
}

void MDCELL_ZVelocity(MDCELL_ThreadJob<MMOL_MD_CELL> *job, int nJob) {
	ARRAY<double> mp;
	mp.set_array(nJob); memset(mp.m, 0, sizeof(double) * mp.n);
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, double >((void*)(&MDCELL_SpatialMomentumZ), job, nJob, mp.m);
	double Pz = 0;
	for (int i = 0; i < mp.n; i++) Pz += mp.m[i];
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	mdcell->Vmc.v[2] = Pz / mdcell->M;
}

void MDCELL_ZeroSpatialMomentumZ(MDCELL_ThreadJob< MMOL_MD_CELL> *job, double *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;

	double v = *mp;
	VECTOR3 P;

	if (FABS(v) < 1e-10) return;
	else show_log("zero the total spatial momentum of the cell", true);

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		mm->V0.v[5] -= v;
		for (ic = 0; ic < mm->nCluster; ic++) {
			mm->cluster[ic].bkm->V0.v[5] -= v;
			mm->cluster[ic].bkm->P0.v[5] -= mm->cluster[ic].M * v;
		}
		mm->mDynKin.V.v[5] -= v;
		mm->mDynKin.P.v[5] -= mm->mDynKin.M * v;
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		sm->c->bkm->V0.v[5] -= v;
		sm->c->bkm->P0.v[5] -= sm->c->M * v;
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		pm->v.v[2] -= v;
		pm->c->bkm->P0.v[5] -= pm->c->M * v;
	}
}

void ZeroSpatialMomentumZ(MDCELL_ThreadJob< MMOL_MD_CELL> *job, double *mp, int nJob) {
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, double >((void*)(&MDCELL_SpatialMomentumZ), job, nJob, mp);
	//VECTOR3 P, v, Pw, w;
	double Pz = 0, vz = 0;
	float M = job->job->mdcell->M;
	int i;
	for (i = 0; i < nJob; i++) Pz += mp[i];
	vz = Pz / M;
	
	for (i = 0; i < nJob; i++) mp[i] = vz;
	
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, double >((void*)(&MDCELL_ZeroSpatialMomentumZ), job, nJob, mp);
}


void MDCELL_NetForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job, VECTOR3 *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	VECTOR3 F;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			F.v[0] += pc->dyn->fc.v[3];
			F.v[1] += pc->dyn->fc.v[4];
			F.v[2] += pc->dyn->fc.v[5];
		}
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; M = sm->c->M;
		F.v[0] += sm->c->dyn->fc.v[3];
		F.v[1] += sm->c->dyn->fc.v[4];
		F.v[2] += sm->c->dyn->fc.v[5];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol; M = pm->c->M;
		F.v[0] += pm->c->dyn->fc.v[3];
		F.v[1] += pm->c->dyn->fc.v[4];
		F.v[2] += pm->c->dyn->fc.v[5];
	}
	mp->v[0] = F.v[0]; mp->v[1] = F.v[1]; mp->v[2] = F.v[2];
}

void NetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob) {
	VECTOR3 fc[MAX_THREADS];
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, VECTOR3 >((void*)(&MDCELL_NetForce), job, nJob, fc);
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	VECTOR3 F;
	float M = mdcell->M;
	int i;
	for (i = 0; i < nJob; i++) {
		F.v[0] += fc[i].v[0];
		F.v[1] += fc[i].v[1];
		F.v[2] += fc[i].v[2];
	}
	memcpy(mdcell->fc.v, F.v, 3 * SIZE_DOUBLE);
	/*
	mdcell->ac.v[0] = F.v[0] / M;
	mdcell->ac.v[1] = F.v[1] / M;
	mdcell->ac.v[2] = F.v[2] / M;
	*/
}

void MDCELL_CompensateNetForceZ(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 3> *sv3) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	double *ac = mdcell->ac.v;
	double *Vmc = mdcell->Vmc.v;

	double *str = sv3->v;
	memset(str, 0, SIZE_V3);

	double Fz, rz;
	double ita_mc = ::ita_mc;

	VECTOR3 force, torq, dr;
	// we do not need to calculate the torque to mass center of macromolecule

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			//pc->dyn->fc.v[3] -= pc->M * ac[0];
			//pc->dyn->fc.v[4] -= pc->M * ac[1];
			Fz = pc->M * ac[2];
			pc->dyn->fc.v[5] -= Fz;

			force.v[2] = Fz;
			V3PV3(pc->cm0, force, torq)
			pc->dyn->fc.v[0] -= torq.v[0];
			pc->dyn->fc.v[1] -= torq.v[1];
			pc->dyn->fc.v[2] -= torq.v[2];
			
			//rz = pc->cm0.v[2] - mdcell->rc.v[2];
			rz = pc->cm0.v[2];
			str[2] -= Fz * rz;
		}
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; M = sm->c->M;
		//sm->c->dyn->fc.v[3] -= M * ac[0];
		//sm->c->dyn->fc.v[4] -= M * ac[1];

		Fz = M * ac[2];
		sm->c->dyn->fc.v[5] -= Fz;

		//rz = sm->r.v[2] - mdcell->rc.v[2];
		rz = sm->r.v[2];
		str[2] -= Fz * rz;
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol; M = pm->c->M;
		//pm->c->dyn->fc.v[3] -= M * ac[0];
		//pm->c->dyn->fc.v[4] -= M * ac[1];

		Fz = M * ac[2];
		pm->c->dyn->fc.v[5] -= Fz;

		//rz = pm->r.v[2] - mdcell->rc.v[2];
		rz = pm->r.v[2];
		str[2] -= Fz * rz;
	}
}

void MDCELL_CompensateNetForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 3> *sv3) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	double *ac = mdcell->ac.v;

	double *str = sv3->v;
	memset(str, 0, SIZE_V3);

	VECTOR3 F, T, dr;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			F.v[0] = pc->M * ac[0];
			F.v[1] = pc->M * ac[1];
			F.v[2] = pc->M * ac[2];
			dr = pc->cm0;
			V3PV3(dr, F, T); 

			pc->dyn->fc.v[0] -= T.v[0];
			pc->dyn->fc.v[1] -= T.v[1];
			pc->dyn->fc.v[2] -= T.v[2];
			pc->dyn->fc.v[3] -= F.v[0];
			pc->dyn->fc.v[4] -= F.v[1];
			pc->dyn->fc.v[5] -= F.v[2];

			V32V3(pc->cm, dr)
			//VECT3(mdcell->rc, pc->cm0, dr)
			str[0] -= F.v[0] * dr.v[0];
			str[1] -= F.v[1] * dr.v[1];
			str[2] -= F.v[2] * dr.v[2];
		}
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; M = sm->c->M;
		F.v[0] = M * ac[0];
		F.v[1] = M * ac[1];
		F.v[2] = M * ac[2];
		sm->c->dyn->fc.v[3] -= F.v[0];
		sm->c->dyn->fc.v[4] -= F.v[1];
		sm->c->dyn->fc.v[5] -= F.v[2];

		V32V3(sm->r, dr)
		//VECT3(mdcell->rc, sm->r, dr)
		str[0] -= F.v[0] * dr.v[0];
		str[1] -= F.v[1] * dr.v[1];
		str[2] -= F.v[2] * dr.v[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol; M = pm->c->M;
		F.v[0] = M * ac[0];
		F.v[1] = M * ac[1];
		F.v[2] = M * ac[2];
		pm->c->dyn->fc.v[3] -= F.v[0];
		pm->c->dyn->fc.v[4] -= F.v[1];
		pm->c->dyn->fc.v[5] -= F.v[2];

		V32V3(pm->r, dr)
		//VECT3(mdcell->rc, pm->r, dr)
		str[0] -= F.v[0] * dr.v[0];
		str[1] -= F.v[1] * dr.v[1];
		str[2] -= F.v[2] * dr.v[2];
	}
}

void CompensateNetForceZ(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 3> *sv3) {
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 3> >((void*)(&MDCELL_CompensateNetForceZ), job, nJob, sv3);
}

void CompensateNetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 3> *sv3) {
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 3> >((void*)(&MDCELL_CompensateNetForce), job, nJob, sv3);
}
/*
void MDCELL_MassCenterZ(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 6> *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	VECTOR3 Rcm;
	double Rz = 0;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			//Rcm.v[0] += pc->M * (pc->cm.v[0] + pc->dr_cmm.v[0]);
			//Rcm.v[1] += pc->M * (pc->cm.v[1] + pc->dr_cmm.v[1]);
			//Rcm.v[2] += pc->M * (pc->cm.v[2] + pc->dr_cmm.v[2]);
			Rz += pc->M * pc->vp->v[2]; //(pc->cm.v[2] + pc->dr_cmm.v[2]);
		}
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	double *r = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; M = sm->c->M; r = sm->r.v;
		//Rcm.v[0] += M * r[0];
		//Rcm.v[1] += M * r[1];
		//Rcm.v[2] += M * r[2];
		Rz += M * r[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol; M = pm->c->M; r = pm->r.v;
		//Rcm.v[0] += M * r[0];
		//Rcm.v[1] += M * r[1];
		//Rcm.v[2] += M * r[2];
		Rz += M * r[2];
	}
	//memcpy(mp->v, Rcm.v, SIZE_V3);
	mp->v[0] = Rz;
}

void MDCELL_ApplyZSpring(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 6> *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;

	double az = mp->v[0];
	if (FABS(az) < 1e-8) return;
	VECTOR3 f, t, dr;
	memset(mp->v + 3, 0, SIZE_V3);

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			f.v[2] = pc->M * az;
			VECT3((*(pc->Op)), pc->cm, dr)
			V3PV3(dr, f, t)
			pc->dyn->fc.v[0] += t.v[0];
			pc->dyn->fc.v[1] += t.v[1];
			pc->dyn->fc.v[2] += t.v[2];
			pc->dyn->fc.v[5] += f.v[2];

			mp->v[5] += f.v[2] * pc->cm0.v[2];
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		f.v[2] = sm->c->M * az;
		sm->c->dyn->fc.v[5] += f.v[2];

		mp->v[5] += f.v[2] * sm->r.v[2];

	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		f.v[2] = pm->c->M * az;
		pm->c->dyn->fc.v[5] += f.v[2];

		mp->v[5] += f.v[2] * pm->r.v[2];
	}
}

void ApplyZSpring(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 6> *mp, int nJob, float *k) {
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 6> >((void*)(&MDCELL_MassCenterZ), job, nJob, mp);
	float M = job->job->mdcell->M;
	double Rz = 0;
	int i;
	for (i = 0; i < nJob; i++) Rz += mp[i].v[0];
	Rz /= M;
	mp[0].v[0] = (k[1] * (Rz > 0 ? -1 : 1) * Rz * Rz - k[0] * Rz) / M; // add a spring-force
	for (i = 1; i < nJob; i++) mp[i].v[0] = mp[0].v[0];  

	job->job->mdcell->rc.v[2] = Rz;
	
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 6> >((void*)(&MDCELL_ApplyZSpring), job, nJob, mp);
}
*/
void MDCELL_calHooverVirial(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 6> *sv6) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;

	double *str = sv6->v;
	memset(str, 0, SIZE_V6);
	double *f, *r;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			f = pc->dyn->f_Hoover.v + 3;
			r = pc->cm0.v;
			str[0] += f[0] * r[0];
			str[1] += f[1] * r[1];
			str[2] += f[2] * r[2];
		}
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; 
		f = sm->c->dyn->f_Hoover.v + 3;
		r = sm->r.v;
		str[0] += f[0] * r[0];
		str[1] += f[1] * r[1];
		str[2] += f[2] * r[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		f = pm->c->dyn->f_Hoover.v + 3;
		r = pm->r.v;
		str[0] += f[0] * r[0];
		str[1] += f[1] * r[1];
		str[2] += f[2] * r[2];
	}
}

void calHooverVirial(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 6> *sv6) {
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 6> >((void*)(&MDCELL_calHooverVirial), job, nJob, sv6);
}


extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
extern void MDCELL_CalMolTransKineticEnergy(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4> *Ekt);
extern void InterfaceConstraints_CMM(CMM_CELL3D<BASIC_CLUSTER> *cmm, int nThreads, double &Ut);
extern char parfname_IConstraint[256];

void SLAB_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm) {
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}

	// we need to make sure the mass center of the cell is at zero
	molecule_check_2d_periodic_cell(mmcell); // check the molecule mass center is inside the cell
	VECTOR3 mass_center;
	set_free_cell_mass_center(mmcell, mass_center);

	show_log("running molecular dynamic on slab ...", true);

	molecule_check_2d_periodic_cell(mmcell); // check the molecule mass center is inside the cell

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************
	show_all_molecules();

	if (strlen(parfname_IConstraint) > 0) {
		if (init_clusterIC_DB(2, cmm->xl[2], cmm->xr[2])) {
			sprintf(errmsg, "Interface defined in %s is enabled!", parfname_IConstraint); show_log(errmsg, true);
			apply_slab_barrier_constraint(clusterIC_db);
		}
		else if (clusterIC_db.n == 0 || clusterIC_db.n != 2) {
			sprintf(errmsg, "2 interfaces are required to define a slab. check %s !", parfname_IConstraint);
			show_log(errmsg, true); return;
		}
	}
	else {
		sprintf(errmsg, "To ensure the species are not evaporated out of FFT space, external interfacial constraint is required!");
		show_log(errmsg, true); //return;
	}

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	if (::rcut_Ewd > cmm->xd[0] / 3) ::rcut_Ewd = cmm->xd[0] / 3;
	if (::rcut_Ewd > cmm->xd[1] / 3) ::rcut_Ewd = cmm->xd[1] / 3;
	if (::rcut_Ewd > cmm->xd[2] / 3) ::rcut_Ewd = cmm->xd[2] / 3;

	bool bStorageChain = false;
	InitMD(mmcell, cmm, bStorageChain);


	Slab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	mdvar.q0 = 0; // a charge for the interactions between induced dipoles. The charge has to be set to 0, important!!

	mdvar.set_virial(::bVirial, true);
	mdvar.spme_var.spme_3d_var.bEfieldOnly = false;
	mdvar.spme_var.spme_3d_var.esv1.bEwald = true; 
	//mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(cmm->xd[0] * cmm->xd[1] * cmm->xd[2], ::rcut_Ewd);
	mdvar.spme_var.spme_3d_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	mdvar.spme_var.spme_3d_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;

	init_BSplineFunc<2>(::bsp, ::indx_bspline);
	mdvar.vspme_q.mMP = 0; // charge only
	mdvar.vspme_mu.mMP = 1; // with dipole
	mdvar.vspme_mu_induced.mMP = 1; // with dipole of course
	if (bPolarize) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		init_spme(&mmcell, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.set_BSplineFunc(&(::bsp));
		mdvar.vspme_mu.set_cell(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], true, ::bExpandSPMECellZ); // important, after init();
		mdvar.vspme_mu.init_k(); mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();

		init_spme_induced(&mmcell, mdvar.vspme_mu_induced, &(mdvar.q0)); // init the viral atoms
		mdvar.vspme_mu_induced.set_BSplineFunc(&(::bsp));
		mdvar.vspme_mu_induced.set_cell(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], true, ::bExpandSPMECellZ); // important, after init();
		mdvar.vspme_mu_induced.init_k(); mdvar.vspme_mu_induced.init_b(); mdvar.vspme_mu_induced.init_C(); mdvar.vspme_mu_induced.init_vars();

		mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(mdvar.vspme_mu.V, ::rcut_Ewd);

		if (::bExpandSPMECellZ) {
			sprintf(errmsg, "SPME: expand cell z-axis from %f to %f", cmm->xd[2], mdvar.vspme_mu.xd[2]); show_log(errmsg, true);
		}
		else {
			sprintf(errmsg, "SPME: keep cell z-axis : %f / %f", cmm->xd[2], mdvar.vspme_mu.xd[2]); show_log(errmsg, true);
		}
	}
	else if (::bDipole) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		init_spme(&mmcell, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.set_BSplineFunc(&(::bsp));
		mdvar.vspme_mu.set_cell(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], true, ::bExpandSPMECellZ); // important, after init();
		mdvar.vspme_mu.init_k(); mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();

		mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(mdvar.vspme_mu.V, ::rcut_Ewd);

		if (::bExpandSPMECellZ) {
			sprintf(errmsg, "SPME: expand cell z-axis from %f to %f", cmm->xd[2], mdvar.vspme_mu.xd[2]); show_log(errmsg, true);
		}
		else {
			sprintf(errmsg, "SPME: keep cell z-axis : %f / %f", cmm->xd[2], mdvar.vspme_mu.xd[2]); show_log(errmsg, true);
		}
	}
	else { // charge only
		mmcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		init_spme(&mmcell, mdvar.vspme_q); // init the viral atoms
		mdvar.vspme_q.set_BSplineFunc(&(::bsp));
		mdvar.vspme_q.set_cell(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], true, ::bExpandSPMECellZ); // important, after init();
		mdvar.vspme_q.init_k(); mdvar.vspme_q.init_b(); mdvar.vspme_q.init_C(); mdvar.vspme_q.init_vars();

		mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(mdvar.vspme_q.V, ::rcut_Ewd);

		if (::bExpandSPMECellZ) {
			sprintf(errmsg, "SPME: expand cell z-axis from %f to %f", cmm->xd[2], mdvar.vspme_q.xd[2]); show_log(errmsg, true);
		}
		else {
			sprintf(errmsg, "SPME: keep cell z-axis : %f / %f", cmm->xd[2], mdvar.vspme_q.xd[2]); show_log(errmsg, true);
		}
	}

	double kappa = mdvar.spme_var.spme_3d_var.esv1.kappa;
	mdvar.spme_var.K0var.bEfieldOnly = false;
	mdvar.spme_var.K0var.set(kappa, cmm->xd[0] * cmm->xd[1]);

	mdvar.spme_var.spme_3d_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	mdvar.spme_var.K0var.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);

	sprintf(errmsg, "Ewald-Sum kappa = %f, rcut = %f", kappa, ::rcut_Ewd); show_log(errmsg, true);
#endif

	int i = 0;
	int nloop = 0, nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0, U_LJ = 0, U_estat = 0;
	double Ektxx, Ektyy, Ektzz;

	double Uext = 0; // external electric field

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	int nMDThreads = MAX_THREADS;

	MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> mdcellHDynThreadVar[MAX_THREADS];
	AssignJob_MDCELL_HDYN<MMOL_MD_CELL>(&mmcell, mdcellHDynThreadVar, nMDThreads, false);

	mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::init_job(mmcell, cmm);
	sprintf(errmsg, "Job # : %d; working with molecules:", MAX_THREADS); show_log(errmsg, true);
	for (i = 0; i < MAX_THREADS; i++) {
		sprintf(errmsg, "%d -- MM %d, SM %d, PM %d", i, mdvar.mdcellJob[i].mmIndx.n, mdvar.mdcellJob[i].smIndx.n, mdvar.mdcellJob[i].pmIndx.n);
		show_log(errmsg, true);
		//for (nm = 0; nm < mdvar.mdcellJob[i].mmIndx.n; nm++) {
		//	sprintf(errmsg, "%d ", mdvar.mdcellJob[i].mmIndx.m[nm]); show_log(errmsg, false);
		//}
		//show_log("", true);
	}

	double djob[MAX_THREADS];
	MATRIX<3> mjob[MAX_THREADS];
	bool status_job[MAX_THREADS];
	//INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];

	INTERACT_2d explicit_interact_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		explicit_interact_var[i].evar.bEfieldOnly = false;
		explicit_interact_var[i].evar.bVirial = bVirial;
		explicit_interact_var[i].evar.bST_diagonal = true;
		explicit_interact_var[i].evar.esv1.bEwald = false;
	}


	// make sure the spatial momentum of the whole cell along z is zero
	if (mmcell.bHDynMC_T) MDCELL_mass_center(mmcell);
	ZeroSpatialMomentumZ(mdvar.job_mdcell, djob, MAX_THREADS);

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;
	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;
	bool bCheckClusterInCMM = true;
	bool Overlap = false;
	int nConvg = 0;
	bool bSpeedVerlet = true;

	mmcell.Init_MD(); //copy the initial torsion angle, speed to the mdpar
	if (bPolarize) mmcell.init_dipole_hist(); // has to be done after Init_MD();

	long LOOPS = mmcell.LOOPS;
	int save_indx = 0;
	MD_SAVE *mdsave = mmcell.msave;
	AsyncMDSave_VARS asyncMDSave;

	bool bSaveVirial = true;
	ITER_SAVE virial_save;
	if (::bVirial) {
		virial_save.set_file(mdsave->title, "virial");
		virial_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	ITER_SAVE Pz_save;
	bool bSavePz = true;
	if (bSavePz) {
		Pz_save.set_file(mdsave->title, "Pz");
		Pz_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	ITER_SAVE mu_save;
	bool bSave_mu = (bPolarize || ::bDipole ? true : false);
	if (bSave_mu) {
		mu_save.set_file(mdsave->title, "mu");
		mu_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	LJ_AsynVARS asynVar_LJ;

	int ncs_mmcell = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) ncs_mmcell += mmcell.mm.m[nm].nCluster;
	ncs_mmcell += mmcell.sm.n + mmcell.pm.n;

#if _CHECK_TIME_STAMP_
	tsave.set_file(mdsave->title, "t");
	//tsave.set_iter(mdsave->nsave, mdsave->max_save);
	tsave.set_iter(1, mdsave->max_save);
	long dt = 0;
#endif

	mdvar.spme_var.spme_3d_var.esv1.iSurface = 1; // 2d
	mdvar.spme_var.spme_3d_var.esv1.iSurfaceBoundary = ::iSurfaceBoundary; // normal
	mdvar.spme_var.spme_3d_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	mdvar.spme_var.K0var.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command

		iCheckCMM = 0;

		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			//case MULTI_MM_PERIDIC_CELL_MD:
				// here we have to check the periodity of the molecule position
				//molecule_check_periodic_cell(mmcell);
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_molecule_2d_PeriodicCell_check), 
					mdvar.job_mdcell, nMDThreads);
			//	break;
		}

#if _CHECK_TIME_STAMP_
tsave.one_more();
mt.start();
if (tsave.bsave()) (*tsave.out)<<nloop<<"  ";
#endif

		check_atom_rg(mmcell, nMDThreads);

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalInertiaMoment), mdvar.job_mdcell, nMDThreads);
		
		MDCellMassCenter(mdvar.job_mdcell, sv3job, nMDThreads);

		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_InitClusterForce), mdvar.job_mdcell, nMDThreads);

	// re-distribute / check clusters in CMM, use vp of each cluster only
		if (bStorageChain) {
			if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
			else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
		}
		else {
			if (nloop == 0 || bCheckClusterInCMM) {
				if (MAX_THREADS == 1) {
					mdvar.job_mdcellCMM[0].cmm = cmm; // directly allocate to global cmm
					cmm->reset_cell_chain();
					MMCELL_CMM_cluster_allocate(mdvar.job_mdcellCMM, status_job);
				}
				else {
					mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::job_mdcell_CMM_use_internal_cmm(); // set the JOB-related CMM is the job-related CMM
					MTOperate1<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, bool>(
						(void*)(&MMCELL_CMM_cluster_allocate), mdvar.job_mdcellCMM, nMDThreads, status_job);
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
				mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::job_mdcell_CMM_use_external_cmm(cmm);
				//for (i = 0; i < MAX_THREADS; i++) job_mdcellCMM[i].cmm = cmm; // reset the JOB-related CMM is the main CMM
			}

			if (bCheckClusterInCMM) {
			{
				nc = cmm->number_atoms();
				/*
				for (i = 0; i < cmm->acell.n; i++) {
					nc += cmm->acell.m[i]->length();
					//sprintf(errmsg, "%d %d", i, cmm->acell.m[i]->length()); show_log(errmsg, true);
				}
				*/
				if (nc != ncs_mmcell) {
					sprintf(errmsg, "%d : cluster in CMM -- %d, real -- %d", nloop, nc, ncs_mmcell); show_log(errmsg, true);
				}
			}
			}
			
		}
	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		dt = mt.elapse();
		(*tsave.out)<<dt<<"  ";
	#endif

		Ep = 0;
		// check relationship of each cluster, electrostatic and LJ
		if (iCheckCMM == 0) {
		//#if _CHECK_TIME_STAMP_ == 1
		//	time.start();
		//#endif
			if (bStorageChain) {
				//ParallelCheckRelshipInCells(*cmm, 0, cmm->acell.n - 1, MAX_THREADS); // setup relationship for each atom with CMM
				if (mmcell.mm.n > 0) _cmm_2d_::MMCheckRelshipInCell(mmcell.mm.m, mmcell.mm.n, *cmm, nMDThreads);
				if (mmcell.sm.n > 0) _cmm_2d_::SMCheckRelshipInCell(mmcell.sm.m, mmcell.sm.n, *cmm, nMDThreads);
				if (mmcell.pm.n > 0) _cmm_2d_::PMCheckRelshipInCell(mmcell.pm.m, mmcell.pm.n, *cmm, nMDThreads);
			}
			else {
				MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >(
					(void*)(&MMCELL_CMM_CheckRelship_2d), mdvar.job_mdcellCMM, nMDThreads);
			}
		//#if _CHECK_TIME_STAMP_ == 1
		//	sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
		//	show_infor(errmsg);
		//#endif
		}

		// explicit interaction calculation 
		/*
		MTOperate1< MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, INTERACT_2d>((void*)(&MMCELL_ExplicitInteract_2d), mdvar.job_mdcellCMM, nMDThreads, explicit_interact_var);
		Ep = 0; mdvar.virial.v[0] = 0; mdvar.virial.v[1] = 0; mdvar.virial.v[2] = 0; mdvar.virial.v[3] = 0;
		for (i = 0; i < nMDThreads; i++) {
			Ep += explicit_interact_var[i].LJ_res.U + explicit_interact_var[i].eres.U; 
			mdvar.virial.v[0] += explicit_interact_var[i].LJ_res.STxx + explicit_interact_var[i].eres.STxx;
			mdvar.virial.v[1] += explicit_interact_var[i].LJ_res.STyy + explicit_interact_var[i].eres.STyy;
			mdvar.virial.v[2] += explicit_interact_var[i].LJ_res.STzz + explicit_interact_var[i].eres.STzz;
			mdvar.virial.v[3] += explicit_interact_var[i].LJ_res.STtr + explicit_interact_var[i].eres.STtr;
		}
		*/


		// important : calculate LJ interaction first, here cluster force & torque are reset
		if (::bMTForce) {
			Asyn_LJ_Interact(&mmcell, cmm, MAX_THREADS, mdvar.MDVar_LJ::LJ_var, mdvar.MDVar_LJ::LJ_res, asynVar_LJ);
		}
		else {
			LJ_Interact(&mmcell, cmm, MAX_THREADS, mdvar.MDVar_LJ::LJ_var, mdvar.MDVar_LJ::LJ_res);
			U_LJ = mdvar.MDVar_LJ::LJ_res.U;
		}

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!mmcell.eNeutral && ::bRealEwaldSumOnly) {
			if (::bPolarize || ::bDipole) _spme_2d_::LocalEF(&mmcell, cmm, mdvar.vspme_mu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			else _spme_2d_::LocalEF(&mmcell, cmm, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
		}
		else if (!::local_interact && !mmcell.eNeutral) {
			// For 2d x-y periodical EwalsSum, we do not need the total net dipole for EwaldSum!!
			// And, we have to set net dipole to zero because real part 2d EwaldSum used 3d EwaldSum program
			// where net cell dipole is used for energy, E and F calculation!!
			//V3zero(mmcell.mu_surf.mu) V3zero(mmcell.mu_surf.mu0) V3zero(mmcell.mu_surf.qr) mmcell.mu_surf.q = 0;
			// calculat the total dipole of the whole cell
			//cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));

			mdvar.spme_var.spme_3d_var.bEfieldOnly = false;
			mdvar.spme_var.K0var.bEfieldOnly = false;
			// important : calculate Coulomb secondly, here cluster force & torque are accumulated
			if (::bPolarize) {
				if (nloop >= 3) Guess_induced_dipole(mmcell, MAX_THREADS);
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				if (::bExcludeInducedDipole) {
					if (!_spme_2d_::Polarize_SPME_Interact2(&mmcell, cmm, mdvar.vspme_mu, mdvar.vspme_mu_induced, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
						sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
					}
				}
				else if (!_spme_2d_::Polarize_SPME_Interact(&mmcell, cmm, mdvar.vspme_mu, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
				}
				Backup_induced_dipole_hist(mmcell, MAX_THREADS);
			}
			else if (::bDipole) {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				_spme_2d_::SPME_Interact(&mmcell, cmm, mdvar.vspme_mu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			}
			else {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				_spme_2d_::SPME_Interact(&mmcell, cmm, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			}

			U_estat = mdvar.spme_res.U;
		}
#endif
		
#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
		}
#endif

		if (::bMTForce) {
			asynVar_LJ.thread_var.WaitUntilThreadOver();
			U_LJ = mdvar.MDVar_LJ::LJ_res.U;
		}

		Ep = U_LJ + U_estat;
		if (::bVirial) {
			mdvar.virial.v[0] = mdvar.LJ_res.STxx + mdvar.spme_res.STxx;
			mdvar.virial.v[1] = mdvar.LJ_res.STyy + mdvar.spme_res.STyy;
			mdvar.virial.v[2] = mdvar.LJ_res.STzz + mdvar.spme_res.STzz;
			mdvar.virial.v[3] = mdvar.LJ_res.STtr + mdvar.spme_res.STtr;

			if (!virial_save.one_more()) break;
			bSaveVirial = virial_save.bsave(); 
		}
		
		//virial.v[0] = 0; virial.v[1] = 0; virial.v[2] = 0; virial.v[3] = 0; Ep = 0;

		if (bSavePz) {
			if (!Pz_save.one_more()) break;
		}

		if (bSave_mu) {
			mu_save.one_more();
			if (mu_save.bsave()) {
				*(mu_save.out)<<nloop<<"  "<<cmm->mu_surf.mu.v[0]<<"  "<<cmm->mu_surf.mu.v[1]<<"  "<<cmm->mu_surf.mu.v[2]<<endl;
			}
		}

		if (mmcell.mm.n > 0) {
			MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, TorsionVar, InteractRes>(
				(void*)(&MMCELL_MM_Torsion), mdvar.job_mdcell, nMDThreads, mdvar.torsion_var, mdvar.mires);
			Etorsion = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ep += mdvar.mires[i].U; Etorsion += mdvar.mires[i].U;
				if (mdvar.torsion_var[i].bVirial) {
					mdvar.virial.v[0] += mdvar.mires[i].STxx;
					mdvar.virial.v[1] += mdvar.mires[i].STyy;
					mdvar.virial.v[2] += mdvar.mires[i].STzz;
					mdvar.virial.v[3] += mdvar.mires[i].STtr;
				}
			}
			//sprintf(errmsg, "%d: Torsion %f kT", nloop, Etorsion); show_log(errmsg, true);
		}

		if (::bEext) {
			if (mmcell.mm.m != NULL) {
				MOperate3<MMOLECULE, float, double>((void*)(&MM_ExternalEField), mmcell.mm.m, mmcell.mm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.sm.m != NULL) {
				MOperate3<SMOLECULE, float, double>((void*)(&SM_ExternalEField), mmcell.sm.m, mmcell.sm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.pm.m != NULL) {
				MOperate3<PMOLECULE, float, double>((void*)(&PM_ExternalEField), mmcell.pm.m, mmcell.pm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
		}


		//if (mmcell.mm.m != NULL) MOperate<MMOLECULE>(MM_TOTAL_FORCE, mmcell.mm.m, mmcell.mm.n, MAX_THREADS);
		//if (mmcell.sm.m != NULL) MOperate<SMOLECULE>(SM_TOTAL_FORCE, mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
		//if (mmcell.pm.m != NULL) MOperate<PMOLECULE>(PM_TOTAL_FORCE, mmcell.pm.m, mmcell.pm.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ClusterForce), mdvar.job_mdcell, nMDThreads);

		if (::bVirial && bSaveVirial && (mmcell.mm.n > 0 || mmcell.sm.n > 0)) {
			for (i = 0; i < nMDThreads; i++) memset(sv4job[i].v, 0, 4 * SIZE_DOUBLE);
			switch (::method_Virial) {
			case CLUSTER_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >(
					(void*)(&MDCELL_VirialCorrect_RigidCluster), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			case MOLECULE_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >(
					(void*)(&MDCELL_VirialCorrect_Molecule), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			}
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] -= sv4job[i].v[0];
				mdvar.virial.v[1] -= sv4job[i].v[1];
				mdvar.virial.v[2] -= sv4job[i].v[2];
				mdvar.virial.v[3] -= sv4job[i].v[0] + sv4job[i].v[1] + sv4job[i].v[2];
			}
		}
		
		if (::bVirial) {
			if (virial_save.bsave()) {
				(*virial_save.out)<<nloop<<"  "<<mdvar.virial.v[0]<<"  "<<mdvar.virial.v[1]
					<<"  "<<mdvar.virial.v[2]<<"  "<<mdvar.virial.v[3];
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
			MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >(
				(void*)(&MMCELL_CMM_CheckImpSolvRelship_2d), mdvar.job_mdcellCMM, nMDThreads);
			MTOperate1< MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, double>(
				(void*)(&MMCELL_cluster_impsolv_interact), mdvar.job_mdcellCMM, nMDThreads, djob);
			for (i = 0; i < nMDThreads; i++) Ep += djob[i];
		}

		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ScaleClusterForce), mdvar.job_mdcell, nMDThreads);

		if (bRandomForce) gen_random_force(mmcell, mdvar.rand[0], mdvar.rand[1]);

		if (_bias_force_::biasOn) apply_point_bias_force(mmcell, mdcell_bias);

		if (!::local_interact && ::mode_NetForce > 0) {
			NetForce(mdvar.job_mdcell, nMDThreads); // calculate net force and acceleration of the whole cell
			if (::bZSpring) {
				mmcell.fc.v[2] += ::k_zspring[0] * mmcell.rc.v[2] + (mmcell.rc.v[2] > 0 ? 1 : -1) * ::k_zspring[1] * mmcell.rc.v[2] * mmcell.rc.v[2];
				MDCELL_ZVelocity(mdvar.job_mdcell, nMDThreads); // calculate z-velocity of cell
				mmcell.fc.v[2] += ::ita_mc * mmcell.Vmc.v[2] * mmcell.M;
			}
			if (::mode_NetForce == 1) {
				mmcell.ac.v[0] = mmcell.fc.v[0] / mmcell.M;
				mmcell.ac.v[1] = mmcell.fc.v[1] / mmcell.M;
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}
			else if (::mode_NetForce == 2) {
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}

			if (::mode_NetForce == 1) CompensateNetForce(mdvar.job_mdcell, nMDThreads, sv3job); // compensate the net force on all specises in the cell
			else if (::mode_NetForce == 2) CompensateNetForceZ(mdvar.job_mdcell, nMDThreads, sv3job); // compensate the net force on all specises in the cell
			/*
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] += sv6job[i].v[0];
				mdvar.virial.v[1] += sv6job[i].v[1];
				mdvar.virial.v[2] += sv6job[i].v[2];
				mdvar.virial.v[3] += (sv6job[i].v[0] + sv6job[i].v[1] + sv6job[i].v[2]);
			}
			*/
		}

	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
		}
	#endif

		if (clusterIC_db.n > 0) {
			// the force from interface, and the torque
			// however, the force will not contribute to the virial
			InterfaceConstraints_CMM(cmm, nMDThreads, djob[0]); Ep += djob[0];
		}

	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
		}
	#endif

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		// self-consistent velocity and acceleration, Nose-Hoover Chain
		bSpeedVerlet = (nloop == 0 ? false : true);
		if (!LeapFrogVerlet_SelfConsistent_SpeedAccel(&mmcell, mdvar.job_mdcell, mdcellHDynThreadVar, bSpeedVerlet, nMDThreads)) {
			sprintf(errmsg, "failure to get consistent velocity and acceleration at loop %d", nloop);
			show_log(errmsg, true);
		}

	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
		}
	#endif

		Ep_t = Ep; Ek_t = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) Ek_t += mmcell.mm.m[nm].Ek;
		for (nm = 0; nm < mmcell.sm.n; nm++) Ek_t += mmcell.sm.m[nm].Ek;
		for (nm = 0; nm < mmcell.pm.n; nm++) Ek_t += mmcell.pm.m[nm].Ek;

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);

		// save md asynchrotronly
		//if (mmcell.msave->bBuffered) asyncMDSave.WaitUntilThreadOver();
		//mmcell.MDSave(nloop, (float)Ep_t, false, true);
		//AsyncMDSave(mmcell, asyncMDSave);

		mmcell.MDSave(nloop, (float)Ep_t, false);

		//if (::bPolarize && mdsave->mKinSave.bsave()) log_dipole(&mmcell); // write the dipole of cluster into log file

		if (::bVirial && bSaveVirial && ::method_Virial == CLUSTER_VIRIAL) {
			if (mmcell.mm.n > 0 && mmcell.bHingeConstraint) {
				MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, HingeConstraintVar, InteractRes>(
					(void*)(&MDCELL_ClusterVirialCorrect_HingeConstraint), mdvar.job_mdcell, nMDThreads, 
					mdvar.hinge_constraint_var, mdvar.mires);
				for (i = 0; i < nMDThreads; i++) {
					mdvar.virial.v[0] += mdvar.mires[i].STxx;
					mdvar.virial.v[1] += mdvar.mires[i].STyy;
					mdvar.virial.v[2] += mdvar.mires[i].STzz;
					mdvar.virial.v[3] += mdvar.mires[i].STtr;
				}
			}
/*
			calHooverVirial(mdvar.job_mdcell, nMDThreads, sv6job);
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] += sv6job[i].v[0];
				mdvar.virial.v[1] += sv6job[i].v[1];
				mdvar.virial.v[2] += sv6job[i].v[2];
				mdvar.virial.v[3] += sv6job[i].v[0] + sv6job[i].v[1] + sv6job[i].v[2];
			}
*/
			MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_CalMolTransKineticEnergy), mdvar.job_mdcell, nMDThreads, sv4job);
			Ek_t = 0; Ektxx = 0; Ektyy = 0; Ektzz = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ektxx += sv4job[i].v[0]; 
				Ektyy += sv4job[i].v[1];
				Ektzz += sv4job[i].v[2];
				Ek_t += sv4job[i].v[3];
			}
			mdvar.Pint = (mdvar.virial.v[3] + Ek_t * 2 / unit_Ek_kT) / (::fUnit_atm_u * mmcell.h[0] * mmcell.h[1] * mmcell.h[2] * 8 * 3);
			if (bSaveVirial) {
				(*virial_save.out)<<"          "<<(mdvar.virial.v[0] + Ektxx * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pxx -- mN/meter
				(*virial_save.out)<<"  "<<(mdvar.virial.v[1] + Ektyy * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pyy -- mN/meter
				(*virial_save.out)<<"  "<<(mdvar.virial.v[2] + Ektzz * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"  "<<0.5 * (mdvar.virial.v[2] - (mdvar.virial.v[0] + mdvar.virial.v[1]) * 0.5 + (Ektzz - 0.5 * (Ektxx + Ektyy)) * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"  "<<0.5 * (mdvar.virial.v[2] - (mdvar.virial.v[0] + mdvar.virial.v[1]) * 0.5) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"          "<<mdvar.Pint<<"      "<<mmcell.rc.v[2]<<endl;
/*
				(*virial_save.out)<<"          "<<Ektxx;  // Pxx -- mN/meter
				(*virial_save.out)<<"  "<<Ektyy;  // Pyy -- mN/meter
				(*virial_save.out)<<"  "<<Ektzz;  // Pzz -- mN/meter
				(*virial_save.out)<<"  "<<Ektzz - 0.5 * (Ektxx + Ektyy);  // Pzz -- mN/meter
				(*virial_save.out)<<"  "<<0.5 * (mdvar.virial.v[2] - (mdvar.virial.v[0] + mdvar.virial.v[1]) * 0.5) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"          "<<mdvar.Pint<<endl;
				*/
			}
		}

		if (bSavePz && Pz_save.bsave()) {
			{
				MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, double >((void*)(&MDCELL_SpatialMomentumZ), mdvar.job_mdcell, MAX_THREADS, djob);
				double Pz = 0, vz = 0;
				float M = mmcell.M;
				int i;
				for (i = 0; i < MAX_THREADS; i++) Pz += djob[i];
				vz = Pz / M;
				*(Pz_save.out)<<nloop<<"  ";
				*(Pz_save.out)<<mmcell.fc.v[0]<<"  "<<mmcell.fc.v[1]<<"  "<<mmcell.fc.v[2]<<"  "<<mmcell.Ekt_mc<<"  "<<Pz * Pz / M * 0.5 * unit_Ek_kT<<"  "<<Pz<<"  "<<vz;
				if (::bZSpring) *(Pz_save.out)<<"  "<<mmcell.rc.v[2];
				*(Pz_save.out)<<endl;
			}
		}

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_Verlet), mdvar.job_mdcell, nMDThreads);
		// accept current NHC vars and prepair for next step
		for (i = 0; i < mmcell.hdyn.n; i++) mmcell.hdyn.m[i].nhc->accept_nextstep();
		if (mmcell.bHDynMC_T) mmcell.mcHDyn.nhc->accept_nextstep();

		#if _CHECK_TIME_STAMP_ == 1
			//mt.start();
		#endif
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			if (mmcell.mm.n > 2) {
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CheckRotAngle), mdvar.job_mdcell, nMDThreads);
			}
			else {
				for (nm = 0; nm < mmcell.mm.n; nm++) {
					mm = mmcell.mm.m + nm;
					check_angle(*mm, mmcell.lfmd.m[nm]);
				}
			}
		}
		iCheckAngle++; if (iCheckAngle >= nCheckAngle) iCheckAngle = 0;

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mmcell.CalibrateTorsionAngle();

		#if _CHECK_TIME_STAMP_ == 1
			//sprintf(::errmsg, "checking procedure takes %d ms", time.glance());
			//show_infor(errmsg);
		#endif

		if (::bVirial) virial_save.save_over();
		if (bSavePz) Pz_save.save_over();
		if (bSave_mu) mu_save.save_over();
#if _CHECK_TIME_STAMP_
		if (tsave.bsave()) (*tsave.out)<<endl;
		tsave.save_over();
#endif

		nloop++;
	}
	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
	mdsave->flush_buf();
	if (::bVirial) virial_save.reset();
	if (bSavePz) Pz_save.reset();
	if (bSave_mu) mu_save.reset();
#if _CHECK_TIME_STAMP_
	tsave.reset();
#endif
}




extern char parfname_IConstraint[256];
void ImageSlab_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm) {
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}

	// we need to make sure the mass center of the cell is at zero
	VECTOR3 mass_center;
	set_free_cell_mass_center(mmcell, mass_center);

	molecule_check_2d_periodic_cell(mmcell); // check the molecule mass center is inside the cell

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************
	show_all_molecules();

	if (strlen(parfname_IConstraint) > 0) {
		if (init_clusterIC_DB(2, cmm->xl[2], cmm->xr[2])) {
			sprintf(errmsg, "Interface defined in %s is enabled!", parfname_IConstraint); show_log(errmsg, true);
			apply_slab_barrier_constraint(clusterIC_db);
		}
	}

	if (clusterIC_db.n == 0 || clusterIC_db.n != 2) {
		sprintf(errmsg, "2 interfaces are required to define a slab. check %s !", parfname_IConstraint);
		show_log(errmsg, true); return;
	}


#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	if (::rcut_Ewd > cmm->xd[0] / 3) ::rcut_Ewd = cmm->xd[0] / 3;
	if (::rcut_Ewd > cmm->xd[1] / 3) ::rcut_Ewd = cmm->xd[1] / 3;
	if (::rcut_Ewd > cmm->xd[2] / 3) ::rcut_Ewd = cmm->xd[2] / 3;

	bool bStorageChain = false;
	InitMD(mmcell, cmm, bStorageChain);

	int nc_cmm = cmm->number_atoms();

	float f1_img = 0, z1_img = 0, f2_img = 0, z2_img = 0;
	if (!read_interface_image_constant(parfname_IConstraint, f1_img, z1_img, f2_img, z2_img)) return;

	ARRAY<ImageCluster> cimg1, cimg2;
	ARRAY<MATOM*> aimg1, aimg2;
	_spme_2d_::CreateImageCluster(mmcell, z1_img, f1_img, cimg1, aimg1, mmcell.ccluster.n);
	_spme_2d_::CreateImageCluster(mmcell, z2_img, f2_img, cimg2, aimg2, mmcell.ccluster.n * 2);
	//UpdatePos_ImageCluster(cimg1, z1, cimg2, z2, nThreads); // position was updated in CreateImageCluster already
	_spme_2d_::cmm_init_distribute_cluster(cimg1.m, cimg1.n, *cmm, false);
	_spme_2d_::cmm_init_distribute_cluster(cimg2.m, cimg2.n, *cmm, false);

	ARRAY<BASIC_CLUSTER*> clusters;
	clusters.set_array(mmcell.ccluster.n + cimg1.n + cimg2.n);
	{
		int i;
		for (i = 0; i < mmcell.ccluster.n; i++) clusters.m[i] = mmcell.ccluster.m[i];
		for (i = 0; i < cimg1.n; i++) clusters.m[mmcell.ccluster.n + i] = cimg1.m + i;
		for (i = 0; i < cimg2.n; i++) clusters.m[mmcell.ccluster.n + cimg1.n + i] = cimg2.m + i;
	}

	nc_cmm = cmm->number_atoms();

	ImageSlab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	mdvar.CMM_AtomJob_Var<BASIC_CLUSTER>::init_job(clusters);

	init_BSplineFunc<2>(::bsp, ::indx_bspline);
	mdvar.vspme_q.set_MP(0); // charge only
	mdvar.vspme_mu.set_MP(1); // with dipole
	if (bPolarize) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		_spme_2d_::UpdateDipole_ImageCluster(cimg1, f1_img, cimg2, f2_img, MAX_THREADS);

		_spme_2d_::init_spme(mmcell.catom, aimg1, aimg2, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.f1 = f1_img;
		mdvar.vspme_mu.f2 = f2_img;
		mdvar.vspme_mu.z1 = z1_img;
		mdvar.vspme_mu.z2 = z2_img;
		mdvar.vspme_mu.set_BSplineFunc(&(::bsp));
		//mdvar.vspme_mu.set_unit_cells(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]); // important, after init();
		mdvar.vspme_mu.set_cells(cmm->xl[0], cmm->xl[1], (z1_img < z2_img ? z1_img : z2_img), cmm->xd[0], cmm->xd[1], fabs(z2_img - z1_img), ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]); // important, after init();

		mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(mdvar.vspme_mu.vspme.V, ::rcut_Ewd);
	}
	else { // charge only
		mmcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		_spme_2d_::init_spme(mmcell.catom, aimg1, aimg2, mdvar.vspme_q); // init the viral atoms
		mdvar.vspme_q.f1 = f1_img;
		mdvar.vspme_q.f2 = f2_img;
		mdvar.vspme_q.z1 = z1_img;
		mdvar.vspme_q.z2 = z2_img;
		mdvar.vspme_q.set_BSplineFunc(&(::bsp));
		//mdvar.vspme_q.set_unit_cells(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]); // important, after init();
		mdvar.vspme_q.set_cells(cmm->xl[0], cmm->xl[1], (z1_img < z2_img ? z1_img : z2_img), cmm->xd[0], cmm->xd[1], fabs(z2_img - z1_img), ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]); // important, after init();

		mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(mdvar.vspme_q.vspme.V, ::rcut_Ewd);
	}
#endif

	int i = 0;
	int nloop = 0, nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0, U_LJ = 0, U_estat = 0;
	double Ektxx, Ektyy, Ektzz;

	mdvar.set_virial(::bVirial, true);

	mdvar.spme_var.spme_3d_var.bEfieldOnly = false;
	mdvar.spme_var.spme_3d_var.esv1.bEwald = true; 
	//mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(cmm->xd[0] * cmm->xd[1] * cmm->xd[2], ::rcut_Ewd);
	mdvar.spme_var.spme_3d_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	mdvar.spme_var.spme_3d_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	
	double kappa = mdvar.spme_var.spme_3d_var.esv1.kappa;
	mdvar.spme_var.K0var.bEfieldOnly = false;
	mdvar.spme_var.K0var.set(kappa, cmm->xd[0] * cmm->xd[1]);

	mdvar.spme_var.spme_3d_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	mdvar.spme_var.K0var.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);


	double Uext = 0; // external electric field

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	int nMDThreads = MAX_THREADS;

	MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> mdcellHDynThreadVar[MAX_THREADS];
	AssignJob_MDCELL_HDYN<MMOL_MD_CELL>(&mmcell, mdcellHDynThreadVar, nMDThreads, false);

	mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::init_job(mmcell, cmm);
	sprintf(errmsg, "Job # : %d; working with molecules:", MAX_THREADS); show_log(errmsg, true);
	for (i = 0; i < MAX_THREADS; i++) {
		sprintf(errmsg, "%d -- MM %d, SM %d, PM %d", i, mdvar.mdcellJob[i].mmIndx.n, mdvar.mdcellJob[i].smIndx.n, mdvar.mdcellJob[i].pmIndx.n);
		show_log(errmsg, true);
		//for (nm = 0; nm < mdvar.mdcellJob[i].mmIndx.n; nm++) {
		//	sprintf(errmsg, "%d ", mdvar.mdcellJob[i].mmIndx.m[nm]); show_log(errmsg, false);
		//}
		//show_log("", true);
	}

	double djob[MAX_THREADS];
	MATRIX<3> mjob[MAX_THREADS];
	bool status_job[MAX_THREADS];
	//INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];

	INTERACT_2d explicit_interact_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		explicit_interact_var[i].evar.bEfieldOnly = false;
		explicit_interact_var[i].evar.bVirial = bVirial;
		explicit_interact_var[i].evar.bST_diagonal = true;
		explicit_interact_var[i].evar.esv1.bEwald = false;
	}


	// make sure the spatial momentum of the whole cell is zero
	if (mmcell.bHDynMC_T) MDCELL_mass_center(mmcell);
	ZeroSpatialMomentumZ(mdvar.job_mdcell, djob, MAX_THREADS);


	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;
	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;
	bool bCheckClusterInCMM = true;
	bool Overlap = false;
	int nConvg = 0;
	bool bSpeedVerlet = true;

	mmcell.Init_MD(); //copy the initial torsion angle, speed to the mdpar
	if (bPolarize) mmcell.init_dipole_hist(); // has to be done after Init_MD();

	long LOOPS = mmcell.LOOPS;
	int save_indx = 0;
	MD_SAVE *mdsave = mmcell.msave;
	AsyncMDSave_VARS asyncMDSave;

	bool bSaveVirial = false;
	ITER_SAVE virial_save;
	if (::bVirial) {
		virial_save.set_file(mdsave->title, "virial");
		virial_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	ITER_SAVE Pz_save;
	bool bSavePz = true;
	if (bSavePz) {
		Pz_save.set_file(mdsave->title, "Pz");
		Pz_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	int ncs_mmcell = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) ncs_mmcell += mmcell.mm.m[nm].nCluster;
	ncs_mmcell += mmcell.sm.n + mmcell.pm.n;

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command

		iCheckCMM = 0;

		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			//switch(md_mode) {
			//case MULTI_MM_PERIDIC_CELL_MD:
				// here we have to check the periodity of the molecule position
				//molecule_check_periodic_cell(mmcell);
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_molecule_2d_PeriodicCell_check), 
					mdvar.job_mdcell, nMDThreads);
			//	break;
			//}
		}

		check_atom_rg(mmcell, nMDThreads);
		_spme_2d_::UpdatePos_ImageCluster(cimg1, z1_img, cimg2, z2_img, nMDThreads);

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalInertiaMoment), mdvar.job_mdcell, nMDThreads);

		MDCellMassCenter(mdvar.job_mdcell, sv3job, nMDThreads);
		
		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_InitClusterForce), mdvar.job_mdcell, nMDThreads);

	// re-distribute / check clusters in CMM, use vp of each cluster only
	#if _CHECK_TIME_STAMP_ == 1
		mt.start();
	#endif
		if (bStorageChain) {
			if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
			else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
		}
		else {
			if (nloop == 0 || bCheckClusterInCMM) {
				/*
				if (MAX_THREADS == 1) {
					mdvar.job_mdcellCMM[0].cmm = cmm; // directly allocate to global cmm
					MMCELL_CMM_cluster_allocate(mdvar.job_mdcellCMM, status_job);
				}
				else {
					mdvar.MDVar_Job::job_mdcell_CMM_use_internal_cmm(); // set the JOB-related CMM is the job-related CMM
					MTOperate1<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, bool>(
						(void*)(&MMCELL_CMM_cluster_allocate), mdvar.job_mdcellCMM, nMDThreads, status_job);
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
				mdvar.MDVar_Job::job_mdcell_CMM_use_external_cmm(cmm);
				//for (i = 0; i < MAX_THREADS; i++) job_mdcellCMM[i].cmm = cmm; // reset the JOB-related CMM is the main CMM
				*/

				if (MAX_THREADS == 1) {
					// directly allocate to global cmm
					cmm->reset_cell_chain();
					AllocateCluster_CMM_array(clusters, cmm);
				}
				else {
					mdvar.cmm_distribute_cluster_using_internal_cmm(); // set the JOB-related CMM is the job-related CMM
					MTOperate< JOB_THREAD_VARS<CMM_AtomJob<BASIC_CLUSTER> > >(
						(void*)(&AllocateCluster_CMM), mdvar.cluster_distribt_JobVar, nMDThreads);
					//for (i = 0; i < nMDThreads; i++) {
					//	nc = mdvar.job_cmm[i].number_atoms();
					//}
					if (!combine_AtomAllocatedCMM_array<BASIC_CLUSTER>(mdvar.job_cmm, nMDThreads, cmm, true, true)) {
						sprintf(errmsg, "loop %d : buffer of CMM cell leaks when combining allocated cluster", nloop);
						show_log(errmsg, true);
					}
				}
			}

			/*
			if (bCheckClusterInCMM) {
			{
				nc = cmm->number_atoms();
				//if (nc != ncs_mmcell) {
				if (nc != nc_cmm) {
					sprintf(errmsg, "%d : cluster in CMM -- %d, real -- %d", nloop, nc, ncs_mmcell); show_log(errmsg, true);
				}
			}
			}
			*/
			
		}
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "cmm_check takes %d ms", mt.glance());
		show_infor(errmsg);
	#endif

		Ep = 0;
		// check relationship of each cluster, electrostatic and LJ
		if (iCheckCMM == 0) {
		#if _CHECK_TIME_STAMP_ == 1
			mt.start();
		#endif
			if (bStorageChain) {
				if (mmcell.mm.n > 0) _cmm_2d_::MMCheckRelshipInCell(mmcell.mm.m, mmcell.mm.n, *cmm, nMDThreads);
				if (mmcell.sm.n > 0) _cmm_2d_::SMCheckRelshipInCell(mmcell.sm.m, mmcell.sm.n, *cmm, nMDThreads);
				if (mmcell.pm.n > 0) _cmm_2d_::PMCheckRelshipInCell(mmcell.pm.m, mmcell.pm.n, *cmm, nMDThreads);
				//ParallelCheckRelshipInCells(*cmm, 0, cmm->acell.n - 1, MAX_THREADS); // setup relationship for each atom with CMM
			}
			else {
				MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >(
					(void*)(&MMCELL_CMM_CheckRelship_2d), mdvar.job_mdcellCMM, nMDThreads);
			}
			_cmm_2d_::ImageClusterCheckRelshipInCell(cimg1.m, cimg1.n, *cmm, nMDThreads);
			_cmm_2d_::ImageClusterCheckRelshipInCell(cimg2.m, cimg2.n, *cmm, nMDThreads);

		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "Relationship_check takes %d ms", mt.glance());
			show_infor(errmsg);
		#endif
		}

		// explicit interaction calculation 
		/*
		MTOperate1< MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, INTERACT_2d>((void*)(&MMCELL_ExplicitInteract_2d), mdvar.job_mdcellCMM, nMDThreads, explicit_interact_var);
		Ep = 0; mdvar.virial.v[0] = 0; mdvar.virial.v[1] = 0; mdvar.virial.v[2] = 0; mdvar.virial.v[3] = 0;
		for (i = 0; i < nMDThreads; i++) {
			Ep += explicit_interact_var[i].LJ_res.U + explicit_interact_var[i].eres.U; 
			mdvar.virial.v[0] += explicit_interact_var[i].LJ_res.STxx + explicit_interact_var[i].eres.STxx;
			mdvar.virial.v[1] += explicit_interact_var[i].LJ_res.STyy + explicit_interact_var[i].eres.STyy;
			mdvar.virial.v[2] += explicit_interact_var[i].LJ_res.STzz + explicit_interact_var[i].eres.STzz;
			mdvar.virial.v[3] += explicit_interact_var[i].LJ_res.STtr + explicit_interact_var[i].eres.STtr;
		}
		*/

		// important : calculate LJ interaction first, here cluster force & torque are reset
		
		LJ_Interact(&mmcell, cmm, MAX_THREADS, mdvar.MDVar_LJ::LJ_var, mdvar.MDVar_LJ::LJ_res);
		U_LJ = mdvar.MDVar_LJ::LJ_res.U;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!::local_interact && !mmcell.eNeutral) {
			// For 2d x-y periodical EwalsSum, we do not need the total net dipole for EwaldSum!!
			// And, we have to set net dipole to zero because real part 2d EwaldSum used 3d EwaldSum program
			// where net cell dipole is used for energy, E and F calculation!!
			//V3zero(mmcell.mu_surf.mu) V3zero(mmcell.mu_surf.mu0) V3zero(mmcell.mu_surf.qr) mmcell.mu_surf.q = 0;
			// calculat the total dipole of the whole cell
			//cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));

			mdvar.spme_var.spme_3d_var.bEfieldOnly = false;
			mdvar.spme_var.K0var.bEfieldOnly = false;
			// important : calculate Coulomb secondly, here cluster force & torque are accumulated
			if (::bPolarize) {
				if (nloop > 3) Guess_induced_dipole(mmcell, MAX_THREADS);
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				if (!_spme_2d_::Polarize_SPME_Interact(&mmcell, cimg1, f1_img, cimg2, f2_img, cmm, mdvar.vspme_mu, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
				}
				Backup_induced_dipole_hist(mmcell, MAX_THREADS);
			}
			else {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				_spme_2d_::SPME_Interact(&mmcell, cimg1, cimg2, cmm, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			}

			U_estat = mdvar.spme_res.U;
		}
#endif
		
		Ep = U_LJ + U_estat;
		if (::bVirial) {
			mdvar.virial.v[0] = mdvar.LJ_res.STxx + mdvar.spme_res.STxx;
			mdvar.virial.v[1] = mdvar.LJ_res.STyy + mdvar.spme_res.STyy;
			mdvar.virial.v[2] = mdvar.LJ_res.STzz + mdvar.spme_res.STzz;
			mdvar.virial.v[3] = mdvar.LJ_res.STtr + mdvar.spme_res.STtr;

			if (!virial_save.one_more()) break;
			bSaveVirial = virial_save.bsave(); 
		}
		
		//virial.v[0] = 0; virial.v[1] = 0; virial.v[2] = 0; virial.v[3] = 0; Ep = 0;

		if (bSavePz) {
			if (!Pz_save.one_more()) break;
		}

		if (mmcell.mm.n > 0) {
			MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, TorsionVar, InteractRes>(
				(void*)(&MMCELL_MM_Torsion), mdvar.job_mdcell, nMDThreads, mdvar.torsion_var, mdvar.mires);
			Etorsion = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ep += mdvar.mires[i].U; Etorsion += mdvar.mires[i].U;
				if (mdvar.torsion_var[i].bVirial) {
					mdvar.virial.v[0] += mdvar.mires[i].STxx;
					mdvar.virial.v[1] += mdvar.mires[i].STyy;
					mdvar.virial.v[2] += mdvar.mires[i].STzz;
					mdvar.virial.v[3] += mdvar.mires[i].STtr;
				}
			}
			//sprintf(errmsg, "%d: Torsion %f kT", nloop, Etorsion); show_log(errmsg, true);
		}

		if (::bEext) {
			if (mmcell.mm.m != NULL) {
				MOperate3<MMOLECULE, float, double>((void*)(&MM_ExternalEField), mmcell.mm.m, mmcell.mm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.sm.m != NULL) {
				MOperate3<SMOLECULE, float, double>((void*)(&SM_ExternalEField), mmcell.sm.m, mmcell.sm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.pm.m != NULL) {
				MOperate3<PMOLECULE, float, double>((void*)(&PM_ExternalEField), mmcell.pm.m, mmcell.pm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
		}


		//if (mmcell.mm.m != NULL) MOperate<MMOLECULE>(MM_TOTAL_FORCE, mmcell.mm.m, mmcell.mm.n, MAX_THREADS);
		//if (mmcell.sm.m != NULL) MOperate<SMOLECULE>(SM_TOTAL_FORCE, mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
		//if (mmcell.pm.m != NULL) MOperate<PMOLECULE>(PM_TOTAL_FORCE, mmcell.pm.m, mmcell.pm.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ClusterForce), mdvar.job_mdcell, nMDThreads);

		if (::bVirial && bSaveVirial && (mmcell.mm.n > 0 || mmcell.sm.n > 0)) {
			for (i = 0; i < nMDThreads; i++) {
				sv4job[i].v[0] = 0; sv4job[i].v[1] = 0; sv4job[i].v[2] = 0; sv4job[i].v[3] = 0;
			}
			switch (::method_Virial) {
			case CLUSTER_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >(
					(void*)(&MDCELL_VirialCorrect_RigidCluster), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			case MOLECULE_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >(
					(void*)(&MDCELL_VirialCorrect_Molecule), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			}
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] -= sv4job[i].v[0];
				mdvar.virial.v[1] -= sv4job[i].v[1];
				mdvar.virial.v[2] -= sv4job[i].v[2];
				mdvar.virial.v[3] -= sv4job[i].v[0] + sv4job[i].v[1] + sv4job[i].v[2];
			}
		}
		
		if (::bVirial) {
			if (virial_save.bsave()) {
				(*virial_save.out)<<nloop<<"  "<<mdvar.virial.v[0]<<"  "<<mdvar.virial.v[1]
					<<"  "<<mdvar.virial.v[2]<<"  "<<mdvar.virial.v[3];
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
			MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >(
				(void*)(&MMCELL_CMM_CheckImpSolvRelship_2d), mdvar.job_mdcellCMM, nMDThreads);
			MTOperate1< MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, double>(
				(void*)(&MMCELL_cluster_impsolv_interact), mdvar.job_mdcellCMM, nMDThreads, djob);
			for (i = 0; i < nMDThreads; i++) Ep += djob[i];
		}

		//CellOperate<BASIC_CLUSTER>((void*)(&SCALE_FORCE_CMM), *cmm, 0, cmm->acell.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ScaleClusterForce), mdvar.job_mdcell, nMDThreads);

		if (bRandomForce) gen_random_force(mmcell, mdvar.rand[0], mdvar.rand[1]);

		if (_bias_force_::biasOn) apply_point_bias_force(mmcell, mdcell_bias);

		if (!::local_interact && ::mode_NetForce > 0) {
			NetForce(mdvar.job_mdcell, nMDThreads); // calculate net force and acceleration of the whole cell
			if (::bZSpring) {
				mmcell.fc.v[2] += ::k_zspring[0] * mmcell.rc.v[2] + (mmcell.rc.v[2] > 0 ? 1 : -1) * ::k_zspring[1] * mmcell.rc.v[2] * mmcell.rc.v[2];
				MDCELL_ZVelocity(mdvar.job_mdcell, nMDThreads); // calculate z-velocity of cell
				mmcell.fc.v[2] += ::ita_mc * mmcell.Vmc.v[2] * mmcell.M;
			}
			if (::mode_NetForce == 1) {
				mmcell.ac.v[0] = mmcell.fc.v[0] / mmcell.M;
				mmcell.ac.v[1] = mmcell.fc.v[1] / mmcell.M;
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}
			else if (::mode_NetForce == 2) {
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}

			if (::mode_NetForce == 1) CompensateNetForce(mdvar.job_mdcell, nMDThreads, sv3job); // compensate the net force on all specises in the cell
			else if (::mode_NetForce == 2) CompensateNetForceZ(mdvar.job_mdcell, nMDThreads, sv3job); // compensate the net force on all specises in the cell
			/*
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] += sv6job[i].v[0];
				mdvar.virial.v[1] += sv6job[i].v[1];
				mdvar.virial.v[2] += sv6job[i].v[2];
				mdvar.virial.v[3] += (sv6job[i].v[0] + sv6job[i].v[1] + sv6job[i].v[2]);
			}
			*/
		}

		if (clusterIC_db.n > 0) {
			// the force from interface, and the torque
			// however, the force will not contribute to the virial
			InterfaceConstraints_CMM(cmm, nMDThreads, djob[0]); Ep += djob[0];
		}


		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		// self-consistent velocity and acceleration, Nose-Hoover Chain
		bSpeedVerlet = (nloop == 0 ? false : true);
		if (!LeapFrogVerlet_SelfConsistent_SpeedAccel(&mmcell, mdvar.job_mdcell, mdcellHDynThreadVar, bSpeedVerlet, nMDThreads)) {
			sprintf(errmsg, "failure to get consistent velocity and acceleration at loop %d", nloop);
			show_log(errmsg, true);
		}

		Ep_t = Ep; Ek_t = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) Ek_t += mmcell.mm.m[nm].Ek;
		for (nm = 0; nm < mmcell.sm.n; nm++) Ek_t += mmcell.sm.m[nm].Ek;
		for (nm = 0; nm < mmcell.pm.n; nm++) Ek_t += mmcell.pm.m[nm].Ek;

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);
		
		// save md asynchrotronly
		//if (mmcell.msave->bBuffered) asyncMDSave.WaitUntilThreadOver();
		//mmcell.MDSave(nloop, (float)Ep_t, false, true);
		//AsyncMDSave(mmcell, asyncMDSave);

		mmcell.MDSave(nloop, (float)Ep_t, false);

		//if (::bPolarize && mdsave->mKinSave.bsave()) log_dipole(&mmcell); // write the dipole of cluster into log file

		if (::bVirial && bSaveVirial && ::method_Virial == CLUSTER_VIRIAL) {
			if (mmcell.mm.n > 0 && mmcell.bHingeConstraint) {
				MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, HingeConstraintVar, InteractRes>(
					(void*)(&MDCELL_ClusterVirialCorrect_HingeConstraint), mdvar.job_mdcell, nMDThreads, 
					mdvar.hinge_constraint_var, mdvar.mires);
				for (i = 0; i < nMDThreads; i++) {
					mdvar.virial.v[0] += mdvar.mires[i].STxx;
					mdvar.virial.v[1] += mdvar.mires[i].STyy;
					mdvar.virial.v[2] += mdvar.mires[i].STzz;
					mdvar.virial.v[3] += mdvar.mires[i].STtr;
				}
			}

			MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_CalMolTransKineticEnergy), mdvar.job_mdcell, nMDThreads, sv4job);
			Ek_t = 0; Ektxx = 0; Ektyy = 0; Ektzz = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ektxx += sv4job[i].v[0]; 
				Ektyy += sv4job[i].v[1];
				Ektzz += sv4job[i].v[2];
				Ek_t += sv4job[i].v[3];
			}
			mdvar.Pint = (mdvar.virial.v[3] + Ek_t * 2 / unit_Ek_kT) / (::fUnit_atm_u * mmcell.h[0] * mmcell.h[1] * mmcell.h[2] * 8 * 3);
			if (bSaveVirial) {
				(*virial_save.out)<<"          "<<(mdvar.virial.v[0] + Ektxx * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pxx -- mN/meter
				(*virial_save.out)<<"  "<<(mdvar.virial.v[1] + Ektyy * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pyy -- mN/meter
				(*virial_save.out)<<"  "<<(mdvar.virial.v[2] + Ektzz * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"  "<<0.5 * (mdvar.virial.v[2] - (mdvar.virial.v[0] + mdvar.virial.v[1]) * 0.5 + (Ektzz - 0.5 * (Ektxx + Ektyy)) * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"          "<<mdvar.Pint<<"      "<<mmcell.rc.v[2]<<endl;
			}
		}

		if (bSavePz && Pz_save.bsave()) {
			{
				MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, double >((void*)(&MDCELL_SpatialMomentumZ), mdvar.job_mdcell, MAX_THREADS, djob);
				double Pz = 0, vz = 0;
				float M = mmcell.M;
				int i;
				for (i = 0; i < MAX_THREADS; i++) Pz += djob[i];
				vz = Pz / M;
				*(Pz_save.out)<<nloop<<"  "<<mmcell.Ekt_mc<<"  "<<Pz * Pz / M * 0.5 * unit_Ek_kT<<"  "<<Pz<<"  "<<vz<<endl;
			}
		}

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_Verlet), mdvar.job_mdcell, nMDThreads);
		// accept current NHC vars and prepair for next step
		for (i = 0; i < mmcell.hdyn.n; i++) mmcell.hdyn.m[i].nhc->accept_nextstep();
		if (mmcell.bHDynMC_T) mmcell.mcHDyn.nhc->accept_nextstep();

		#if _CHECK_TIME_STAMP_ == 1
			mt.start();
		#endif
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			if (mmcell.mm.n > 2) {
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CheckRotAngle), mdvar.job_mdcell, nMDThreads);
			}
			else {
				for (nm = 0; nm < mmcell.mm.n; nm++) {
					mm = mmcell.mm.m + nm;
					check_angle(*mm, mmcell.lfmd.m[nm]);
				}
			}
		}
		iCheckAngle++; if (iCheckAngle >= nCheckAngle) iCheckAngle = 0;

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mmcell.CalibrateTorsionAngle();

		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "checking procedure takes %d ms", mt.glance());
			show_infor(errmsg);
		#endif

		if (::bVirial) virial_save.save_over();
		if (bSavePz) Pz_save.save_over();

		nloop++;
	}
	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
	mdsave->flush_buf();
	if (::bVirial) virial_save.reset();
	if (bSavePz) Pz_save.reset();
}
