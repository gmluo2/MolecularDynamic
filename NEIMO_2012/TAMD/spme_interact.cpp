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

#include "EwaldSum.h"
#include "spme.h"
#include "spme_interact.h"

using namespace _EwaldSum_real_;

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

// pair-wise database about Thole-radius between ions in Thole polarization model, s = a * (alpha_i * alpha_j)^1/6
extern PAIRWISE_DB< PAIRWISE<float> > TholeRadius_db; //



void CLUSTER_LJ_INTERACT_ATOM(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, InteractVar &var, InteractRes &res) {
	MATOM *patom = NULL;
	VECTOR3 dr, u;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, LJ_RELSHIP> LJit;

	double R, R2;
	double t1, t2, t3;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, local = true;

	double U_LJ = 0;
	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;

	PAIRWISE<EF_DPROF> *pw = NULL;
	unsigned int key;
	char aindx[2];

	res.reset();

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	//V6zero(pc->dyn->fc)

	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		V32V3(patom->rg, r0)

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		if (patom->par->bSR) aindx[0] = patom->par->indx;

		// the atom has to have short-range interaction, Lennard-Jones and/or other short-range interaction
		if (patom->par->bSR || patom->par->bLJ) {
			//ch_imag_cluster = pc->LJRship.ch;
			//while (ch_imag_cluster != NULL) {
				//pc1 = ch_imag_cluster->p->p;
			LJit.set(&(pc->LJRship)); LJit.start_iterate();
			while (!LJit.eoi()) {
				imag_cluster = LJit.current(); if (imag_cluster == NULL) {LJit.next_iterate(); continue;}
				pc1 = imag_cluster->p;
				if (pc1->cID < 0) {LJit.next_iterate(); continue;}
				//if (!imag_cluster->local) {LJit.next_iterate(); continue;}
				//if (pc1->mIndx == pc->mIndx && imag_cluster->local) {LJit.next_iterate(); continue;}
				if (pc1 == pc && imag_cluster->local) {
					//ch_imag_cluster = ch_imag_cluster->next; 
					LJit.next_iterate();
					continue;
				} // ingnore the interaction inside cluster
				drCell.v[0] = imag_cluster->nx * cell.xd[0];
				drCell.v[1] = imag_cluster->ny * cell.xd[1];
				drCell.v[2] = imag_cluster->nz * cell.xd[2];

				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					patom1 = pc1->atom + na1;
					if (!(patom1->par->bLJ || patom1->par->bSR)) continue;
					for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->rg.v[i] - drCell.v[i];
					//ABS2(dr, i, R2)
					R2 = dr.v[0] * dr.v[0]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[1] * dr.v[1]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[2] * dr.v[2]; if (R2 > r2cut_LJ) continue;
					if (R2 > r2cut_LJ) continue;
					
					if (R2 < _R2min_Interact) continue; // ignore the interaction when they are too close

					if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}
					if (bound) continue; // ignore bound interaction

					R = sqrt(R2); u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;

					V3zero(force)
					if (patom->par->bLJ) pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					else pLJ = NULL;
					if (pLJ != NULL && pLJ->epsLJ > 0.0001 && R2 < pLJ->rcut2) {
						//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
						LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, R, u, i, t1, t3, t2, force)
						//if (scale_force && t3 > force_max) {
						//	t3 = force_max / t3;
						//	force.v[0] *= t3; force.v[1] *= t3; force.v[2] *= t3;
						//}
						if (bound14 && pc != pc1) { //1-4 LJ interaction / 2
							t2 *= A14_LJ; t3 *= A14_LJ;
							for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
						}
						U_LJ += t2 * 0.5;
					}
					if (patom1->par->bSR) {
						aindx[1] = patom1->par->indx;
						key = construct_char2_key_pairwise(aindx);
						pw = sr_db.search(key);
						if (pw != NULL) {
							interpolate2(t1, t2, pw->pv, R, i, t3)
							U_LJ += t1 * 0.5;
							if (t2 > force_max && scale_force) {t2 = force_max; patom->bShocked = true;}
							force.v[0] += t2 * u.v[0]; force.v[1] += t2 * u.v[1]; force.v[2] += t2 * u.v[2];
						}
					}

					if (R < (patom->par->r_min + patom1->par->r_min)) patom->bShocked = true;

					if (var.bVirial) {
						if (var.bST_diagonal) Accumulate_Virial_Diagonal(res, force.v, dr.v);
						else Accumulate_Virial_Trace(res, t3, R);
					}
					//V3PV3(patom->r0, force, torque)
					for (i = 0; i < 3; i++) {
						//pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];
						patom->F.v[i] += force.v[i]; // force has mass unit, was converted in LJ12_6 profile
					}
					//USelf += (dr.v[0] * force.v[0] + dr.v[1] * force.v[1] + dr.v[2] * force.v[2]) * 0.5;
				}
				//ch_imag_cluster = ch_imag_cluster->next;
				LJit.next_iterate();
			}
		}
		#if _CHECK_TIME_STAMP_ == 1
		//t_LJ += time.elapse();
		#endif

	}
	// total energy in kT
	res.U = U_LJ;

	//USelf = 0;
	//char msg[256] = "\0";
	//sprintf(msg, "RF = %f", USelf); show_log(msg, true);
}


void CLUSTER_LJ_INTERACT_FATOM(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, InteractVar &var, InteractRes &res) {
	MATOM *patom = NULL;
	VECTOR3 dr, u;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, LJ_RELSHIP> LJit;

	double R, R2;
	double t1, t2, t3;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, local = true;

	double U_LJ = 0;
	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;

	PAIRWISE<EF_DPROF> *pw = NULL;
	char aindx[2];
	unsigned int key;

	double *f = NULL;

	res.reset();

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];
		V32V3(patom->rg, r0)

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		if (::bMTForce) f = patom->f.m[_IF_SR].v;
		else f = patom->F.v;

		if (patom->par->bSR) aindx[0] = patom->par->indx;

		// LJ interaction
		if (patom->par->bLJ || patom->par->bSR) {
			//ch_imag_cluster = pc->LJRship.ch;
			//while (ch_imag_cluster != NULL) {
				//pc1 = ch_imag_cluster->p->p;
			LJit.set(&(pc->LJRship)); LJit.start_iterate();
			while (!LJit.eoi()) {
				imag_cluster = LJit.current(); if (imag_cluster == NULL) {LJit.next_iterate(); continue;}
				pc1 = imag_cluster->p;
				if (pc1->cID < 0) {
					LJit.next_iterate(); continue;
				}
				//if (!imag_cluster->local) {LJit.next_iterate(); continue;}
				//if (pc1->mIndx == pc->mIndx && imag_cluster->local) {LJit.next_iterate(); continue;}
				if (pc1 == pc && imag_cluster->local) {
					//ch_imag_cluster = ch_imag_cluster->next; 
					LJit.next_iterate();
					continue;
				} // ingnore the interaction inside cluster
				drCell.v[0] = imag_cluster->nx * cell.xd[0];
				drCell.v[1] = imag_cluster->ny * cell.xd[1];
				drCell.v[2] = imag_cluster->nz * cell.xd[2];

				for (na1 = 0; na1 < pc1->fatom.n; na1++) {
					patom1 = pc1->fatom.m[na1];
					if (!(patom1->par->bLJ || patom1->par->bSR)) continue;
					for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->rg.v[i] - drCell.v[i];
					//ABS2(dr, i, R2)
					R2 = dr.v[0] * dr.v[0]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[1] * dr.v[1]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[2] * dr.v[2]; if (R2 > r2cut_LJ) continue;
					if (R2 > r2cut_LJ) continue;
					
					//if (R2 < _R2min_Interact) {patom->bShocked = true; continue;} // ignore the interaction when they are too close

					if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}
					if (bound) continue; // ignore bound interaction

					R = sqrt(R2); u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;

					V3zero(force)
					if (patom1->par->bLJ) pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					else pLJ = NULL;
					if (pLJ != NULL && pLJ->epsLJ > 0.0001 && R2 < pLJ->rcut2) {
						//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
						LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, R, u, i, t1, t3, t2, force)
						//if (scale_force && t3 > force_max) {
						//	t3 = force_max / t3;
						//	force.v[0] *= t3; force.v[1] *= t3; force.v[2] *= t3;
						//}
						if (bound14 && pc != pc1) { //1-4 LJ interaction / 2
							t2 *= A14_LJ; t3 *= A14_LJ;
							for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
						}
						U_LJ += t2 * 0.5;
					}
					if (patom1->par->bSR) {
						aindx[1] = patom1->par->indx;
						key = construct_char2_key_pairwise(aindx);
						pw = sr_db.search(key);
						if (pw != NULL) {
							interpolate2(t1, t2, pw->pv, R, i, t3)
							U_LJ += t1 * 0.5;
							if (t2 > force_max && scale_force) {t2 = force_max; patom->bShocked = true;}
							force.v[0] += t2 * u.v[0]; force.v[1] += t2 * u.v[1]; force.v[2] += t2 * u.v[2];
						}
					}

					if (R < (patom->par->r_min + patom1->par->r_min)) {
						t2 = force_max;
						force.v[0] += t2 * u.v[0]; force.v[1] += t2 * u.v[1]; force.v[2] += t2 * u.v[2];
						patom->bShocked = true;
					}

					if (var.bVirial) {
						if (var.bST_diagonal) Accumulate_Virial_Diagonal(res, force.v, dr.v);
						else Accumulate_Virial_Trace(res, t3, R);
					}
					//V3PV3(patom->r0, force, torque)
					for (i = 0; i < 3; i++) {
						//pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];
						//patom->F.v[i] += force.v[i]; // force has mass unit, was converted in LJ12_6 profile
						f[i] += force.v[i];
					}
					//USelf += (dr.v[0] * force.v[0] + dr.v[1] * force.v[1] + dr.v[2] * force.v[2]) * 0.5;
				}
				//ch_imag_cluster = ch_imag_cluster->next;
				LJit.next_iterate();
			}
		}
		#if _CHECK_TIME_STAMP_ == 1
		//t_LJ += time.elapse();
		#endif

	}
	// total energy in kT
	res.U = U_LJ;

	//USelf = 0;
	//char msg[256] = "\0";
	//sprintf(msg, "RF = %f", USelf); show_log(msg, true);
}

void CLUSTER_LJ_INTERACT(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, InteractVar &var, InteractRes &res) {
	if (pc->fatom.m == NULL) return CLUSTER_LJ_INTERACT_ATOM(cell, pc, var, res);
	else return CLUSTER_LJ_INTERACT_FATOM(cell, pc, var, res);
}

void SingleMM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, InteractVar &var, InteractRes &res) {
	int nc;
	BASIC_CLUSTER *pc = NULL;
	InteractRes tres;
	res.reset();
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = (BASIC_CLUSTER*)(mm.cluster + nc);
		CLUSTER_LJ_INTERACT(cell, pc, var, tres);
		res += tres;
	}
}

void MM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE *mm, int nMM, InteractVar &var, InteractRes &res) {
	int nm;
	InteractRes tres;
	res.reset();
	for (nm = 0; nm < nMM; nm++) {
		SingleMM_LJ(cell, mm[nm], 0, mm[nm].nCluster, var, tres);
		res += tres;
	}
}

void VMM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE **mm, int nMM, InteractVar &var, InteractRes &res) {
	int nm;
	InteractRes tres;
	res.reset();
	for (nm = 0; nm < nMM; nm++) {
		SingleMM_LJ(cell, *mm[nm], 0, mm[nm]->nCluster, var, tres);
		res += tres;
	}
}

void SM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE *sm, int nSM, InteractVar &var, InteractRes &res) {
	int nm;
	InteractRes tres;
	res.reset();
	for (nm = 0; nm < nSM; nm++) {
		CLUSTER_LJ_INTERACT(cell, sm[nm].c, var, tres);
		res += tres;
	}
}

void VSM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE **sm, int nSM, InteractVar &var, InteractRes &res) {
	int nm;
	InteractRes tres;
	res.reset();
	for (nm = 0; nm < nSM; nm++) {
		CLUSTER_LJ_INTERACT(cell, sm[nm]->c, var, tres);
		res += tres;
	}
}


void PM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, PMOLECULE *pm, int nPM, InteractVar &var, InteractRes &res) {
	int nm;
	InteractRes tres;
	res.reset();
	for (nm = 0; nm < nPM; nm++) {
		CLUSTER_LJ_INTERACT(cell, pm[nm].c, var, tres);
		res += tres;
	}
	//char msg[256] = "\0";
	//sprintf(msg, "Rij*Fij = %f", Uself);
	//show_log(msg, true);
}

void VPM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, PMOLECULE **pm, int nPM, InteractVar &var, InteractRes &res) {
	int nm;
	InteractRes tres;
	res.reset();
	for (nm = 0; nm < nPM; nm++) {
		CLUSTER_LJ_INTERACT(cell, pm[nm]->c, var, tres);
		res += tres;
	}
	//char msg[256] = "\0";
	//sprintf(msg, "Rij*Fij = %f", Uself);
	//show_log(msg, true);
}

void LJ_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, int nThreads, InteractVar &var, InteractRes &res) {
	InteractRes tres;
	res.reset();

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, InteractVar, InteractRes>((void*)(&MM_LJ), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, InteractVar, InteractRes>((void*)(&SM_LJ), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, InteractVar, InteractRes>((void*)(&PM_LJ), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var, tres);
		res += tres;
	}
}

void LJ_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, int nThreads, InteractVar &var, InteractRes &res) {
	InteractRes tres;
	res.reset();

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE*, BASIC_CLUSTER, InteractVar, InteractRes>((void*)(&MM_LJ), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE*, BASIC_CLUSTER, InteractVar, InteractRes>((void*)(&SM_LJ), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE*, BASIC_CLUSTER, InteractVar, InteractRes>((void*)(&PM_LJ), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var, tres);
		res += tres;
	}
}



#if _SYS_ == _WINDOWS_SYS_
int _LJ_AsynFunc(LPVOID *vp) {
#elif _SYS_ == _LINUX_SYS_
void* _LJ_AsynFunc(void *vp) {
#endif
	LJ_AsynVARS *par = (LJ_AsynVARS*)vp;
	par->res->reset();

	LJ_Interact(par->mdcell, par->cmm_cell, par->nThreads, *(par->var), *(par->res));
	
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->thread_var.hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void Asyn_LJ_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, int nThreads, InteractVar &var, InteractRes &res, LJ_AsynVARS &av) {
	av.thread_var.set_pars(0);
	av.mdcell = mdcell; av.cmm_cell = cmm_cell;
	av.var = &var;
	av.res = &res;
	av.nThreads = nThreads;
#if _SYS_ == _WINDOWS_SYS_
	ResetEvent(av.thread_var.hMMThread);
	av.thread_var.thread = AfxBeginThread((AFX_THREADPROC)_LJ_AsynFunc, (LPVOID)(&av), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
#elif _SYS_ == _LINUX_SYS_
	pthread_create(&(av.thread_var.thread), NULL, &_LJ_AsynFunc, (void *)(&av));
#endif
}



void CLUSTER_TOTAL_FORCE(BASIC_CLUSTER *pc) {
	int na, i, j;
	MATOM *patom = NULL;
	VECTOR3 torque;

	SVECTOR<double, 3> r0;
	double *v0 = pc->Op->v, *v1, *v = r0.v;

	V6zero(pc->dyn->fc) 
	pc->bShocked = false;
	if (pc->fatom.m == NULL) {
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			if (patom->bShocked) pc->bShocked = true;
			if (::bMTForce) {
				for (j = 0; j < patom->f.n; j++) {
					V3plusV3(patom->F, patom->f.m[j], patom->F)
				}
			}
			V3PV3(patom->r0, patom->F, torque)
			for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[3+i] += patom->F.v[i];}
		}
	}
	else {
		for (na = 0; na < pc->fatom.n; na++) {
			patom = pc->fatom.m[na];
			if (patom->bShocked) pc->bShocked = true;
			if (::bMTForce) {
				for (j = 0; j < patom->f.n; j++) {
					V3plusV3(patom->F, patom->f.m[j], patom->F)
				}
			}
			v1 = patom->r.v;
			DVECT3(v0, v1, v)
			V3PV3(r0, patom->F, torque)
			for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[3+i] += patom->F.v[i];}
		}
	}
}

void SingleMM_TOTAL_FORCE(MMOLECULE *mm, int nc1, int nc2) {
	for (int nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		CLUSTER_TOTAL_FORCE((BASIC_CLUSTER*)(mm->cluster + nc));
	}
}

void MM_TOTAL_FORCE(MMOLECULE *mm, int nMM) {
	for (int nm = 0; nm < nMM; nm++) SingleMM_TOTAL_FORCE(mm + nm, 0, mm[nm].nCluster);
}

void VMM_TOTAL_FORCE(MMOLECULE **mm, int nMM) {
	for (int nm = 0; nm < nMM; nm++) SingleMM_TOTAL_FORCE(mm[nm], 0, mm[nm]->nCluster);
}

void SM_TOTAL_FORCE(SMOLECULE *sm, int nSM) {
	for (int nm = 0; nm < nSM; nm++) {
		CLUSTER_TOTAL_FORCE(sm[nm].c);
	}
}

void VSM_TOTAL_FORCE(SMOLECULE **sm, int nSM) {
	for (int nm = 0; nm < nSM; nm++) {
		CLUSTER_TOTAL_FORCE(sm[nm]->c);
	}
}

void PM_TOTAL_FORCE(PMOLECULE *pm, int nPM) {
	for (int nm = 0; nm < nPM; nm++) {
		CLUSTER_TOTAL_FORCE(pm[nm].c);
	}
}

void VPM_TOTAL_FORCE(PMOLECULE **pm, int nPM) {
	for (int nm = 0; nm < nPM; nm++) {
		CLUSTER_TOTAL_FORCE(pm[nm]->c);
	}
}

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

namespace _EwaldSum_real_ {
// interaction on patom1 from patom2, 
// dr -- vector from patom2 to patom1
// R -- distance |dr|
// R -- |dr|^2
// E -- electric field from patom2
// F -- force from patom2 to patom1
// U -- interaction energy
void EwaldReal(EwaldSumRealVars1 &esv, MATOM *patom1, MATOM *patom2, VECTOR3& dr, double R, double R2, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a, bool bExcludeInducedDipole) {
	int i, j, k;
	double *mu1 = patom1->mu.v, *mu2 = patom2->mu.v;
	double *mu1_ind = patom1->ind_mu.v, *mu2_ind = patom2->ind_mu.v;
	double c1 = patom1->c, c2 = patom2->c;
	double t;

	U = 0; V3zero(E) V3zero(F)
	bool bCharge = (patom1->mMP == 0 && patom2->mMP == 0 ? true : false); // interaction between point charges, without multipole ?

	if (bCharge) esv.EwaldSumRealVars0::init_vars(dr, R, R2);
	else {
		if (patom1->bIMu || patom2->bIMu) esv.init_vars(dr, R, R2);
		else if (bExcludeInducedDipole) esv.init_vars_E(dr, R, R2); // no induced dipole and dipole-dipole interaction will be excluded
		else esv.init_vars(dr, R, R2);
	}

	// U
	t = c1 * c2 / R; // q1 - q2
	U = t * esv.f0;
	if (bEwaldKickoff) U -= t * a;

	if (patom1->mMP > 0) {
		if (bEwaldKickoff) { // mu1 - q2
			U += c2 * (mu1[0] * (esv.mET1.v[0] - a * esv.mT1.v[0]) 
				+ mu1[1] * (esv.mET1.v[1] - a * esv.mT1.v[1]) 
				+ mu1[2] * (esv.mET1.v[2] - a * esv.mT1.v[2]));
		}
		else U += c2 * (mu1[0] * esv.mET1.v[0] + mu1[1] * esv.mET1.v[1] + mu1[2] * esv.mET1.v[2]);

		if (patom2->mMP > 0) { // mu1 - mu2
			if (bEwaldKickoff) {
				U -= mu1[0] * ( (esv.mET2.m[0][0] - a * esv.mT2.m[0][0]) * mu2[0] 
				               + (esv.mET2.m[0][1] - a * esv.mT2.m[0][1]) * mu2[1] 
							   + (esv.mET2.m[0][2] - a * esv.mT2.m[0][2]) * mu2[2] )
					+ mu1[1] * ( (esv.mET2.m[1][0]  - a * esv.mT2.m[1][0]) * mu2[0] 
							   + (esv.mET2.m[1][1]  - a * esv.mT2.m[1][1]) * mu2[1] 
							   + (esv.mET2.m[1][2]  - a * esv.mT2.m[1][2]) * mu2[2] )
					+ mu1[2] * ( (esv.mET2.m[2][0]  - a * esv.mT2.m[2][0]) * mu2[0] 
							   + (esv.mET2.m[2][1]  - a * esv.mT2.m[2][1]) * mu2[1] 
							   + (esv.mET2.m[2][2]  - a * esv.mT2.m[2][2]) * mu2[2] );
			}
			else {
				U -= mu1[0] * ( esv.mET2.m[0][0] * mu2[0] + esv.mET2.m[0][1] * mu2[1] + esv.mET2.m[0][2] * mu2[2])
					+ mu1[1] * ( esv.mET2.m[1][0] * mu2[0] + esv.mET2.m[1][1] * mu2[1] + esv.mET2.m[1][2] * mu2[2])
					+ mu1[2] * ( esv.mET2.m[2][0] * mu2[0] + esv.mET2.m[2][1] * mu2[1] + esv.mET2.m[2][2] * mu2[2]);
			}

			if (bExcludeInducedDipole) { // excluding the interaction between induced dipole part
				if (bEwaldKickoff) {
					U += mu1_ind[0] * ( (esv.mET2.m[0][0] - a * esv.mT2.m[0][0]) * mu2_ind[0] 
					               + (esv.mET2.m[0][1] - a * esv.mT2.m[0][1]) * mu2_ind[1] 
								   + (esv.mET2.m[0][2] - a * esv.mT2.m[0][2]) * mu2_ind[2] )
						+ mu1_ind[1] * ( (esv.mET2.m[1][0]  - a * esv.mT2.m[1][0]) * mu2_ind[0] 
								   + (esv.mET2.m[1][1]  - a * esv.mT2.m[1][1]) * mu2_ind[1] 
								   + (esv.mET2.m[1][2]  - a * esv.mT2.m[1][2]) * mu2_ind[2] )
						+ mu1_ind[2] * ( (esv.mET2.m[2][0]  - a * esv.mT2.m[2][0]) * mu2_ind[0] 
								   + (esv.mET2.m[2][1]  - a * esv.mT2.m[2][1]) * mu2_ind[1] 
								   + (esv.mET2.m[2][2]  - a * esv.mT2.m[2][2]) * mu2_ind[2] );
				}
				else {
					U += mu1_ind[0] * ( esv.mET2.m[0][0] * mu2_ind[0] + esv.mET2.m[0][1] * mu2_ind[1] + esv.mET2.m[0][2] * mu2_ind[2])
						+ mu1_ind[1] * ( esv.mET2.m[1][0] * mu2_ind[0] + esv.mET2.m[1][1] * mu2_ind[1] + esv.mET2.m[1][2] * mu2_ind[2])
						+ mu1_ind[2] * ( esv.mET2.m[2][0] * mu2_ind[0] + esv.mET2.m[2][1] * mu2_ind[1] + esv.mET2.m[2][2] * mu2_ind[2]);
				}
			}
		}
	}
	if (patom2->mMP > 0) { // c1 - mu2
		if (bEwaldKickoff) {
			U -= c1 * ( mu2[0] * (esv.mET1.v[0] - a * esv.mT1.v[0]) 
				      + mu2[1] * (esv.mET1.v[1] - a * esv.mT1.v[1])
					  + mu2[2] * (esv.mET1.v[2] - a * esv.mT1.v[2]) );
		}
		else {
			U -= c1 * (mu2[0] * esv.mET1.v[0] + mu2[1] * esv.mET1.v[1] + mu2[2] * esv.mET1.v[2]);
		}
	}

	// E
	for (i = 0; i < 3; i++) {
		if (bEwaldKickoff) { // from q2
			E.v[i] -= c2 * (esv.mET1.v[i] - a * esv.mT1.v[i]);
			if (patom2->mMP > 0) { // from mu2
				for (j = 0; j < 3; j++) E.v[i] += (esv.mET2.m[i][j] - a * esv.mT2.m[i][j]) * mu2[j];
			}
		}
		else {
			E.v[i] -= c2 * esv.mET1.v[i]; // from q2
			if (patom2->mMP > 0) { // from mu2
				for (j = 0; j < 3; j++) E.v[i] += esv.mET2.m[i][j] * mu2[j];
			}
		}
	}

	// F
	for (i = 0; i < 3; i++) {
		if (bEwaldKickoff) {
			F.v[i] -= c1 * c2 * (esv.mET1.v[i] - a * esv.mT1.v[i]); // q1 - q2
			if (patom2->mMP > 0) { // q1 - mu2
				for (j = 0; j < 3; j++) F.v[i] += c1 * (esv.mET2.m[i][j] - a * esv.mT2.m[i][j]) * mu2[j];
				if (patom1->mMP > 0) { // mu1 - mu2
					for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
						F.v[i] += mu1[j] * (esv.mET3.m[i][j][k] - a * esv.mT3.m[i][j][k]) * mu2[k];
					}}

					if (bExcludeInducedDipole) { // excluding the interaction between induced dipole part
						for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
							F.v[i] -= mu1_ind[j] * (esv.mET3.m[i][j][k] - a * esv.mT3.m[i][j][k]) * mu2_ind[k];
						}}
					}
				}
			}
			if (patom1->mMP > 0) { // mu1 - q2
				for (j = 0; j < 3; j++) F.v[i] -= c2 * (esv.mET2.m[i][j] - a * esv.mT2.m[i][j]) * mu1[j];
			}
		}
		else {
			F.v[i] -= c1 * c2 * esv.mET1.v[i]; // q1 - q2  -- must be wrong
			//F.v[i] -= c1 * c2 * esv.mET1.v[i]; // q1 - q2
			if (patom2->mMP > 0) { // q1 - mu2
				for (j = 0; j < 3; j++) F.v[i] += c1 * esv.mET2.m[i][j] * mu2[j];
				if (patom1->mMP > 0) { // mu1 - mu2
					for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
						F.v[i] += mu1[j] * esv.mET3.m[i][j][k] * mu2[k];
					}}
					if (bExcludeInducedDipole) { // excluding the interaction between induced dipole part
						for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
							F.v[i] -= mu1_ind[j] * esv.mET3.m[i][j][k] * mu2_ind[k];
						}}
					}
				}
			}
			if (patom1->mMP > 0) { // mu1 - q2
				for (j = 0; j < 3; j++) F.v[i] -= c2 * esv.mET2.m[i][j] * mu1[j];
			}
		}
	}
}

// calculate the electric field only about the real part of Ewald Sum
void EwaldReal_E(EwaldSumRealVars1 &esv, MATOM *patom1, MATOM *patom2, VECTOR3& dr, double R, double R2, VECTOR3 &E, bool bEwaldKickoff, double a) {
	int i, j, k;
	double *mu1 = patom1->mu.v, *mu2 = patom2->mu.v;
	double c1 = patom1->c, c2 = patom2->c;
	double t;

	V3zero(E)
	bool bCharge = (patom1->mMP == 0 && patom2->mMP == 0 ? true : false); // interaction between point charges, without multipole ?

	if (bCharge) esv.EwaldSumRealVars0::init_vars(dr, R, R2);
	else esv.init_vars_E(dr, R, R2);

	// E
	for (i = 0; i < 3; i++) {
		if (bEwaldKickoff) { // from q2
			E.v[i] -= c2 * (esv.mET1.v[i] - a * esv.mT1.v[i]);
			if (patom2->mMP > 0) { // from mu2
				for (j = 0; j < 3; j++) E.v[i] += (esv.mET2.m[i][j] - a * esv.mT2.m[i][j]) * mu2[j];
			}
		}
		else {
			E.v[i] -= c2 * esv.mET1.v[i]; // from q2
			if (patom2->mMP > 0) { // from mu2
				for (j = 0; j < 3; j++) E.v[i] += esv.mET2.m[i][j] * mu2[j];
			}
		}
	}
}

void EwaldSum_SurfaceDipole(SURF_DIPOLE &mu_surf, EwaldSumRealVars0 *esv, InteractRes &res) {
	double Uestat = 0, STxx = 0, STyy = 0, STzz = 0;
	double *v1, *v2;
	// energy from the surface induced dipole 
	if (esv->iSurface == 0) { // cubic dipole
		switch (::iSurfaceBoundary) {
		case 0: // normal
			Uestat = 0.5 * esv->inv_3V * (mu_surf.mu.v[0] * mu_surf.mu.v[0] + mu_surf.mu.v[1] * mu_surf.mu.v[1] + mu_surf.mu.v[2] * mu_surf.mu.v[2]);
			if (bExcludeInducedDipole) Uestat -= 0.5 * esv->inv_3V * (mu_surf.mu_ind.v[0] * mu_surf.mu_ind.v[0] 
					+ mu_surf.mu_ind.v[1] * mu_surf.mu_ind.v[1] + mu_surf.mu_ind.v[2] * mu_surf.mu_ind.v[2]);

			// virial coming from the surface dipole, U(V), from receiprcal part
			v1 = mu_surf.qr.v; v2 = mu_surf.mu0.v;
			STxx = (v1[0] * v1[0] * 0.5 + 2 * v1[0] * v2[0] + 1.5 * v2[0] * v2[0]) * esv->inv_3V * ::fUnit_estat_mass;
			STyy = (v1[1] * v1[1] * 0.5 + 2 * v1[1] * v2[1] + 1.5 * v2[1] * v2[1]) * esv->inv_3V * ::fUnit_estat_mass;
			STzz = (v1[2] * v1[2] * 0.5 + 2 * v1[2] * v2[2] + 1.5 * v2[2] * v2[2]) * esv->inv_3V * ::fUnit_estat_mass;

			if (bExcludeInducedDipole) {
				v2 = mu_surf.mu_ind.v;
				STxx -= (1.5 * v2[0] * v2[0]) * esv->inv_3V * ::fUnit_estat_mass;
				STyy -= (1.5 * v2[1] * v2[1]) * esv->inv_3V * ::fUnit_estat_mass;
				STzz -= (1.5 * v2[2] * v2[2]) * esv->inv_3V * ::fUnit_estat_mass;
			}

			break;
		case 1: // disable all axis
			break;
		case 2: // disable z axis
			Uestat = 0.5 * esv->inv_3V * (mu_surf.mu.v[0] * mu_surf.mu.v[0] + mu_surf.mu.v[1] * mu_surf.mu.v[1]);
			if (bExcludeInducedDipole) Uestat -= 0.5 * esv->inv_3V * (mu_surf.mu_ind.v[0] * mu_surf.mu_ind.v[0] 
					+ mu_surf.mu_ind.v[1] * mu_surf.mu_ind.v[1]);

			// virial coming from the surface dipole, U(V), from receiprcal part
			v1 = mu_surf.qr.v; v2 = mu_surf.mu0.v;
			STxx = (v1[0] * v1[0] * 0.5 + 2 * v1[0] * v2[0] + 1.5 * v2[0] * v2[0]) * esv->inv_3V * ::fUnit_estat_mass;
			STyy = (v1[1] * v1[1] * 0.5 + 2 * v1[1] * v2[1] + 1.5 * v2[1] * v2[1]) * esv->inv_3V * ::fUnit_estat_mass;
			//STzz = (v1[2] * v1[2] * 0.5 + 2 * v1[2] * v2[2] + 1.5 * v2[2] * v2[2]) * esv->inv_3V * ::fUnit_estat_mass;

			if (bExcludeInducedDipole) {
				v2 = mu_surf.mu_ind.v;
				STxx -= (1.5 * v2[0] * v2[0]) * esv->inv_3V * ::fUnit_estat_mass;
				STyy -= (1.5 * v2[1] * v2[1]) * esv->inv_3V * ::fUnit_estat_mass;
				//STzz -= (1.5 * v2[2] * v2[2]) * esv->inv_3V * ::fUnit_estat_mass;
			}
			break;
		case 3: // disable x-y axis
			Uestat = 0.5 * esv->inv_3V * mu_surf.mu.v[2] * mu_surf.mu.v[2];
			if (bExcludeInducedDipole) Uestat -= 0.5 * esv->inv_3V * mu_surf.mu_ind.v[2] * mu_surf.mu_ind.v[2];

			// virial coming from the surface dipole, U(V), from receiprcal part
			v1 = mu_surf.qr.v; v2 = mu_surf.mu0.v;
			//STxx = (v1[0] * v1[0] * 0.5 + 2 * v1[0] * v2[0] + 1.5 * v2[0] * v2[0]) * esv->inv_3V * ::fUnit_estat_mass;
			//STyy = (v1[1] * v1[1] * 0.5 + 2 * v1[1] * v2[1] + 1.5 * v2[1] * v2[1]) * esv->inv_3V * ::fUnit_estat_mass;
			STzz = (v1[2] * v1[2] * 0.5 + 2 * v1[2] * v2[2] + 1.5 * v2[2] * v2[2]) * esv->inv_3V * ::fUnit_estat_mass;

			if (bExcludeInducedDipole) {
				v2 = mu_surf.mu_ind.v;
				//STxx -= (1.5 * v2[0] * v2[0]) * esv->inv_3V * ::fUnit_estat_mass;
				//STyy -= (1.5 * v2[1] * v2[1]) * esv->inv_3V * ::fUnit_estat_mass;
				STzz -= (1.5 * v2[2] * v2[2]) * esv->inv_3V * ::fUnit_estat_mass;
			}
			break;
		default:
			break;
		}
	}
	else if (esv->iSurface == 1) { // slab dipole
		switch (::iSurfaceBoundary) {
		case 0: // do nothing
			break;
		default: // disable surface dipole effect along x-y
				// for a slab, net dipole is included inside the receiprocal term already
				// obviously, surface dipole along z axis does NOT exist!
				// to disable the surface dipole effect along x-y, inverse surface dipole along x and y has to be included

			// disable this, do nothing

			/*
			// In another word, assuming surface dipole induced energy is U = 4PI/3V * |M|^2
			// so, here we are to remove it
			Uestat = -0.5 * esv->inv_3V * (mu_surf.mu.v[0] * mu_surf.mu.v[0] + mu_surf.mu.v[1] * mu_surf.mu.v[1]);
			if (bExcludeInducedDipole) Uestat += 0.5 * esv->inv_3V * (mu_surf.mu_ind.v[0] * mu_surf.mu_ind.v[0] 
					+ mu_surf.mu_ind.v[1] * mu_surf.mu_ind.v[1]);

			// virial coming from the surface dipole, U(V), from receiprcal part
			v1 = mu_surf.qr.v; v2 = mu_surf.mu0.v;
			STxx = -(v1[0] * v1[0] * 0.5 + 2 * v1[0] * v2[0] + 1.5 * v2[0] * v2[0]) * esv->inv_3V * ::fUnit_estat_mass;
			STyy = -(v1[1] * v1[1] * 0.5 + 2 * v1[1] * v2[1] + 1.5 * v2[1] * v2[1]) * esv->inv_3V * ::fUnit_estat_mass;
			//STzz = -(v1[2] * v1[2] * 0.5 + 2 * v1[2] * v2[2] + 1.5 * v2[2] * v2[2]) * esv->inv_3V * ::fUnit_estat_mass;

			if (bExcludeInducedDipole) {
				v2 = mu_surf.mu_ind.v;
				STxx += (1.5 * v2[0] * v2[0]) * esv->inv_3V * ::fUnit_estat_mass;
				STyy += (1.5 * v2[1] * v2[1]) * esv->inv_3V * ::fUnit_estat_mass;
				//STzz += (1.5 * v2[2] * v2[2]) * esv->inv_3V * ::fUnit_estat_mass;
			}
			*/
			break;
		}
	}

	res.U = Uestat * eUnit_estat_kT;
	res.STxx = STxx;
	res.STyy = STyy;
	res.STzz = STzz;
	res.STtr = (STxx + STyy + STzz);
}

void CLUSTER_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, SPME_VAR& var, InteractRes& res) {
	using namespace _Thole_polarize_;
	EwaldSumRealVars1 *esv = &(var.esv1);
	TholeCorrTensor *tv = &(var.tv);

	MATOM *patom = NULL;
	VECTOR3 E_real, F_real, E_Thole, F_Thole;
	VECTOR3 dr, u;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	CHAIN<BASIC_CLUSTER> *ch_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, SPME_RELSHIP> rshipIt;

	double R, R2;
	double t1, t2;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, local = true, bEwaldCut = false;
	bool bEwaldKickoff = false, bSameCluster = false;

	double alpha_Thole = 0, a_Thole = ::_Thole_polarize_::a;

	double Utotal = 0, Uestat = 0, Uext = 0, U_real = 0, U_Thole = 0;

	char aindx[2] = {0x00, 0x00};
	unsigned int key;
	PAIRWISE<float> *TholeRadius = NULL;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	double *f = NULL;

	res.reset();

	//char msg[256] = "\0";
	//sprintf(msg, "CLUSTER %d :", pc->mIndx); show_log(msg, true);

	//for (na = 0; na < pc->nAtoms; na++) {
		//patom = pc->atom + na;
	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];
		// reset E of the atom
		V3zero(patom->E)
		//	DO NOT reset F! L-J interaction is stored! 
		//if (!var.bEfieldOnly) {V3zero(patom->F)}
		if (patom->eNeutral) continue;
		Uestat = 0;
		//for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		V32V3(patom->rg, r0)

		if (!var.bEfieldOnly) {
			if (::bMTForce) {f = patom->f.m[_IF_ES].v; memset(f, 0, SIZE_V3);}
			else f = patom->F.v;
		}

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		//ch_imag_cluster = pc->spmeRship.ch;
		//while (ch_imag_cluster != NULL) {
		rshipIt.set(&(pc->spmeRship)); rshipIt.start_iterate();
		while (!rshipIt.eoi()) {
			imag_cluster = rshipIt.current(); if (imag_cluster == NULL) {rshipIt.next_iterate(); continue;}
			pc1 = imag_cluster->p;
			local = imag_cluster->local;
			drCell.v[0] = imag_cluster->nx * cell.xd[0];
			drCell.v[1] = imag_cluster->ny * cell.xd[1];
			drCell.v[2] = imag_cluster->nz * cell.xd[2];
			//for (na1 = 0; na1 < pc1->nAtoms; na1++) {
				//patom1 = pc1->atom + na1;
			for (na1 = 0; na1 < pc1->fatom.n; na1++) {
				patom1 = pc1->fatom.m[na1];
				bSameCluster = (local && pc1 == pc ? true : false);
				if (bSameCluster && na1 == na) continue; // same atom
				if (patom1->eNeutral) continue; // ignore this atom for charge interaction
				for (i = 0; i < 3; i++) {
					//dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
					dr.v[i] = r0.v[i] - patom1->rg.v[i] - drCell.v[i];
				}
				if (FABS(dr.v[0]) > rcut_Ewd || FABS(dr.v[1]) > rcut_Ewd || FABS(dr.v[2]) > rcut_Ewd) {
					bEwaldCut = true;
					continue;
				}
				V3ABS2(dr, R2)
				//if (R2 < _R2min_Interact) continue;
				bEwaldCut = (R2 > r2cut_Ewd ? true : false);
				if (bEwaldCut) continue;

				R = sqrt(R2); 

				//if (R < 3) {
				//	sprintf(msg, "C %d M %d -- C %d M %d,  CELL [%d, %d, %d] : |R| = %f", pc->mIndx, na, pc1->mIndx, na1, imag_cluster->nx, imag_cluster->ny, imag_cluster->nz, R); show_log(msg, true);
				//}

				if (pc->mIndx == pc1->mIndx && R2 < d2_14) {
					BOUND(patom, patom1, i, bound) // bound atom
					if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
					if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
					else bound14 = false;
				}
				else {bound = false; bound14 = false;}

				if (bound || bSameCluster) { // bounding atoms or they are in the same cluster
					/*
					for (i = 0; i < 3; i++) {
						t1 = patom1->c * dr.v[i] / R3;
						E_real.v[i] -= t1 * (1 - f1); // bound interaction has to be kicked off in Ewald Sum
					}
					// energy, no self interaction
					t1 = patom1->c / R;
					Uestat +=  0.5 * t1 * (f0 - 1); // bound interaction has to be kicked off in Ewald Sum
					*/
					bEwaldKickoff = true; t1 = 1;
				}
				else if (bound14 && pc != pc1) { // 1-4 interaction
					/*
					for (i = 0; i < 3; i++) {
						t1 = patom1->c * dr.v[i] / R3;
						E_real.v[i] += t1 * f1;
						E_14.v[i] += (A14_E - 1) * t1; // 1-4 interaction will to be corrected later
					}
					// energy
					t1 = patom1->c / R;
					Uestat +=  0.5 * t1 * (f0 - 1 + A14_E); // 1-4 interaction is corrected by * A14_E
					*/
					bEwaldKickoff = true; t1 = 1 - A14_E;
				}
				else {
					/*
					for (i = 0; i < 3; i++) E_real.v[i] += f1 * patom1->c * dr.v[i] / R3;
					// energy
					t1 = patom1->c / R;
					Uestat += 0.5 * t1 * f0;
					*/
					bEwaldKickoff = false; t1 = 0;
				}
				if (var.bEfieldOnly) {
					EwaldReal_E(*esv, patom, patom1, dr, R, R2, E_real, bEwaldKickoff, t1);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr_E(*tv, a_Thole, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v)) {
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
								}
								V3plusV3(E_real, E_Thole, E_real)
							}
						}
					}
				}
				else {
					EwaldReal(*esv, patom, patom1, dr, R, R2, E_real, F_real, U_real, bEwaldKickoff, t1, bExcludeInducedDipole);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr(*tv, a_Thole, patom->c, patom->mu.v, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v, F_Thole.v, U_Thole)) {
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
									F_Thole.v[0] *= t2; F_Thole.v[1] *= t2; F_Thole.v[2] *= t2;
									U_Thole *= t2;
								}
								V3plusV3(E_real, E_Thole, E_real) V3plusV3(F_real, F_Thole, F_real) U_real += U_Thole;
							}
						}
					}
				}
				patom->E.v[0] += E_real.v[0]; patom->E.v[1] += E_real.v[1]; patom->E.v[2] += E_real.v[2];
				if (!var.bEfieldOnly) {
					// the forece on atom has unit on mass unit
					F_real.v[0] *= fUnit_estat_mass; 
					F_real.v[1] *= fUnit_estat_mass; 
					F_real.v[2] *= fUnit_estat_mass;

					//patom->F.v[0] += F_real.v[0]; 
					//patom->F.v[1] += F_real.v[1]; 
					//patom->F.v[2] += F_real.v[2];
					f[0] += F_real.v[0];
					f[1] += F_real.v[1];
					f[2] += F_real.v[2];

					Uestat += 0.5 * U_real;
					if (var.bVirial) Accumulate_Virial_Diagonal(res, F_real.v, dr.v);
				}
			}
			//ch_imag_cluster = ch_imag_cluster->next;
			rshipIt.next_iterate();
		}

		//show_log("\n", true);
		// the electric field from the surface induced dipole, virial of surface dipole will be calculated later
		if (esv->iSurface == 0) { // cubic
			switch (::iSurfaceBoundary) { // normal
			case 0:
				for (i = 0; i < 3; i++) {
					t1 = esv->inv_3V * cell.mu_surf.mu.v[i];
					patom->E.v[i] -= t1;
					if (!var.bEfieldOnly) {
						//patom->F.v[i] -= patom->c * t1 * fUnit_estat_mass;
						f[i] -= patom->c * t1 * fUnit_estat_mass;
						if (bExcludeInducedDipole) {
							//patom->F.v[i] += esv->inv_3V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;
							f[i] += esv->inv_3V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;
						}
					}
				}
				break;
			case 1: // disable surface dipole in all direction
				break;
			case 2: // disable surface dipole along z axis
				for (i = 0; i < 2; i++) {
					t1 = esv->inv_3V * cell.mu_surf.mu.v[i];
					patom->E.v[i] -= t1;
					if (!var.bEfieldOnly) {
						//patom->F.v[i] -= patom->c * t1 * fUnit_estat_mass;
						f[i] -= patom->c * t1 * fUnit_estat_mass;
						if (bExcludeInducedDipole) {
							//patom->F.v[i] += esv->inv_3V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;
							f[i] += esv->inv_3V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;
						}
					}
				}
				break;
			case 3: // disable surface dipole in x-y axises
				//for (i = 0; i < 3; i++) {
					i = 2;
					t1 = esv->inv_3V * cell.mu_surf.mu.v[i];
					patom->E.v[i] -= t1;
					if (!var.bEfieldOnly) {
						//patom->F.v[i] -= patom->c * t1 * fUnit_estat_mass;
						f[i] -= patom->c * t1 * fUnit_estat_mass;
						if (bExcludeInducedDipole) {
							//patom->F.v[i] += esv->inv_3V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;
							f[i] += esv->inv_3V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;
						}
					}
				break;
			}
		}
		else if (esv->iSurface == 1) { // slab
			switch (::iSurfaceBoundary) {
			case 0: // normal
				// confirmed with explicit calculation, this boundary condition should be disabled
				/*
				for (i = 0; i < 2; i++) {
					//i = 2;
					t1 = esv->inv_V * cell.mu_surf.mu.v[i];
					patom->E.v[i] -= t1;
					if (!var.bEfieldOnly) {
						//patom->F.v[i] -= patom->c * t1 * fUnit_estat_mass;
						f[i] -= patom->c * t1 * fUnit_estat_mass;
						if (bExcludeInducedDipole) {
						//	patom->F.v[i] += esv->inv_V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;
							f[i] += esv->inv_V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;
						}
					}
				}
				*/
				break;
			default: // do not do anything
				/*
				// disable surface dipole effect along x-y
				// for a slab, net dipole is included inside the receiprocal term already
				// obviously, surface dipole along z axis does NOT exist!
				// to disable the surface dipole effect along x-y, inverse surface dipole along x and y has to be included
				
				for (i = 0; i < 2; i++) {
					t1 = esv->inv_V * cell.mu_surf.mu.v[i];
					patom->E.v[i] += t1; // inverse dipole
					if (!var.bEfieldOnly) {
						//patom->F.v[i] -= patom->c * t1 * fUnit_estat_mass;
						f[i] += patom->c * t1 * fUnit_estat_mass; // inverse dipole
						if (bExcludeInducedDipole) {
						//	patom->F.v[i] -= esv->inv_V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;
							f[i] -= esv->inv_V * patom->c * cell.mu_surf.mu_ind.v[i] * fUnit_estat_mass;  // inverse dipole
						}
					}
				}
				*/
				break;
			}
		}

		// self term of electric field
		if (patom->mMP > 0) {
			t1 = 4 * esv->kappa3 / 3 * INV_SQRT_PI;
			patom->E.v[0] += t1 * patom->mu.v[0];
			patom->E.v[1] += t1 * patom->mu.v[1];
			patom->E.v[2] += t1 * patom->mu.v[2];
		}

		if (var.bEfieldOnly) {
			// VERY IMPORTANT :
			// electric field calculation stop here
			// if more code is required for electric field calculation, the code needs to added before this line
			continue;
		}

		Uestat -= esv->kappa * INV_SQRT_PI * patom->c * patom->c; // self energy from charge
		if (patom->mMP > 0) { // self energy from dipole
			V3ABS2(patom->mu, t1)
			Uestat -= 2 * esv->kappa3 / 3 * INV_SQRT_PI * t1;
		}

		// energy and virial term from surface dipole will be calculated later, independently with function
		// void EwaldSum_SurfaceDipole(SURF_DIPOLE &mu_surf, EwaldSumRealVars1 *esv, InteractRes &res)

		Utotal += Uestat;

	#if _CHECK_TIME_STAMP_ == 1
		//t_real += time.elapse();
		//time.start();
	#endif
	}

	if (var.bEfieldOnly) return;

	// virial coming from the surface dipole will be calculated later

	// total energy in kT
	res.U = Utotal * eUnit_estat_kT + Uext;
}

void EFIELD_SPME_REAL_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, SPME_VAR& var, InteractRes& res) {
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> >it;

	res.reset();
	InteractRes tres;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			CLUSTER_EFIELD_SPME_REAL(cell, pc, var, tres);
			res += tres;
			it.next_iterate();
		}
	}
}

void SingleMM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, SPME_VAR& var, InteractRes& res) {
	int nc;
	BASIC_CLUSTER *pc = NULL;

	res.reset();
	InteractRes tres; 
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = (BASIC_CLUSTER*)(mm.cluster + nc);
		CLUSTER_EFIELD_SPME_REAL(cell, pc, var, tres);
		res += tres;
	}
}

void MM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE *mm, int nMM, SPME_VAR& var, InteractRes& res) {
	int nm, nc;
	
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < nMM; nm++) {
		//SingleMM_EFIELD_SPME_REAL(cell, mm[nm], 0, mm[nm].nCluster, U, bEfieldOnly);
		for (nc = 0; nc < mm[nm].nCluster; nc++) {
			CLUSTER_EFIELD_SPME_REAL(cell, mm[nm].cluster + nc, var, tres);
			res += tres;
		}
	}
}

void VMM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE **mm, int nMM, SPME_VAR& var, InteractRes& res) {
	int nm, nc;
	
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < nMM; nm++) {
		//SingleMM_EFIELD_SPME_REAL(cell, mm[nm], 0, mm[nm].nCluster, U, bEfieldOnly);
		for (nc = 0; nc < mm[nm]->nCluster; nc++) {
			CLUSTER_EFIELD_SPME_REAL(cell, mm[nm]->cluster + nc, var, tres);
			res += tres;
		}
	}
}

void SM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE *sm, int nSM, SPME_VAR& var, InteractRes& res) {
	int nm;
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < nSM; nm++) {
		CLUSTER_EFIELD_SPME_REAL(cell, sm[nm].c, var, tres);
		res += tres;
	}
}

void VSM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE **sm, int nSM, SPME_VAR& var, InteractRes& res) {
	int nm;
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < nSM; nm++) {
		CLUSTER_EFIELD_SPME_REAL(cell, sm[nm]->c, var, tres);
		res += tres;
	}
}

void PM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, PMOLECULE *pm, int nPM, SPME_VAR& var, InteractRes& res) {
	int nm;
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < nPM; nm++) {
		CLUSTER_EFIELD_SPME_REAL(cell, pm[nm].c, var, tres);
		res += tres;
	}
}

void VPM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, PMOLECULE **pm, int nPM, SPME_VAR& var, InteractRes& res) {
	int nm;
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < nPM; nm++) {
		CLUSTER_EFIELD_SPME_REAL(cell, pm[nm]->c, var, tres);
		res += tres;
	}
}

} // end of namespace _EwaldSum_real

namespace _spme_ {
using namespace _EwaldSum_real_;

/* multi-thread function to calculate real part of Ewald Sum for electrostatic interaction, to calculate Electric field and F on each atom */
void LocalEF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	InteractRes tres_real;

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var, tres);
		res += tres;
	}
}

/* multi-thread function to calculate real part of Ewald Sum for electrostatic interaction, to calculate Electric field and F on each atom */
void LocalEF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_evatom_::_VMUatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	InteractRes tres_real;

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var, tres);
		res += tres;
	}
}



#if _SYS_ == _WINDOWS_SYS_
int SPME_EF_Real(LPVOID *vp) {
#elif _SYS_ == _LINUX_SYS_
void* SPME_EF_Real(void *vp) {
#endif
	SPME_REAL_VARS *par = (SPME_REAL_VARS*)vp;
	par->res->reset();
	InteractRes tres;
	if (par->mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL), *(par->cmm_cell), par->mdcell->mm.m, par->mdcell->mm.n, par->nThreads, *(par->var), tres);
		*(par->res) += tres;
	}
	if (par->mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL), *(par->cmm_cell), par->mdcell->sm.m, par->mdcell->sm.n, par->nThreads, *(par->var), tres);
		*(par->res) += tres;
	}
	if (par->mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL), *(par->cmm_cell), par->mdcell->pm.m, par->mdcell->pm.n, par->nThreads, *(par->var), tres);
		*(par->res) += tres;
	}
	if (par->cmm_cell->mu_surf.q != 0) par->res->U -= 1.5 * par->var->esv1.inv_3V / (par->var->esv1.kappa * par->var->esv1.kappa) * par->cmm_cell->mu_surf.q * par->cmm_cell->mu_surf.q * eUnit_estat_kT;

#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->thread_var.hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void Asyn_SPME_EF_Real(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, int nThreads, SPME_VAR& var, InteractRes& res, SPME_REAL_VARS &svar) {
	svar.thread_var.set_pars(0);
	svar.mdcell = mdcell;
	svar.cmm_cell = cmm_cell;
	svar.var = &var;
	svar.res = &res;
	svar.nThreads = nThreads;
#if _SYS_ == _WINDOWS_SYS_
	ResetEvent(svar.thread_var.hMMThread);
	svar.thread_var.thread = AfxBeginThread((AFX_THREADPROC)SPME_EF_Real, (LPVOID)(&svar), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
#elif _SYS_ == _LINUX_SYS_
	pthread_create(&(svar.thread_var.thread), NULL, &SPME_EF_Real, (void *)(&svar));
#endif
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	SPME_REAL_VARS svar;
	InteractRes tres_real;

	if (nThreads == 1) {
		if (mdcell->mm.m != NULL) {
			MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var, tres);
			res += tres;
		}
		if (mdcell->sm.m != NULL) {
			MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var, tres);
			res += tres;
		}
		if (mdcell->pm.m != NULL) {
			MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var, tres);
			res += tres;
		}
		if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;

		//sprintf(errmsg, "Real : %f, %f", res.U, res.STtr); show_log(errmsg, true);
	}
	else Asyn_SPME_EF_Real(mdcell, cmm_cell, nThreads, var, tres_real, svar);

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
		EwaldSum_SurfaceDipole(cmm_cell->mu_surf, (EwaldSumRealVars0*)(&(var.esv1)), tres);
		res += tres;
	}
	
	//sprintf(errmsg, "After surface dipole correction, virial changes to %f", res.STtr); show_log(errmsg, true);

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VQatom> *pvcell = (_VE_CELL<_VQatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	if (nThreads > 1) {
		svar.thread_var.WaitUntilThreadOver();
		res += tres_real;
	}

	vspme.dump_EF(!var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
}

void SPME_EF(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE*, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&VMM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE*, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&VSM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE*, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&VPM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var, tres);
		res += tres;
	}
	if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;

	//sprintf(errmsg, "Real : %f, %f", res.U, res.STtr); show_log(errmsg, true);

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
		EwaldSum_SurfaceDipole(cmm_cell->mu_surf, (EwaldSumRealVars0*)(&(var.esv1)), tres);
		res += tres;
	}
	//sprintf(errmsg, "After surface dipole correction, virial changes to %f", res.STtr); show_log(errmsg, true);

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VQatom> *pvcell = (_VE_CELL<_VQatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	vspme.dump_EF(!var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
void SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res);

	if (var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}

void SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res);

	if (var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE*>(VMM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE*>(VSM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE*>(VPM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}

// with dipole
/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_evatom_::_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	SPME_REAL_VARS svar;
	InteractRes tres_real;
	if (nThreads == 1) {
		if (mdcell->mm.m != NULL) {
			MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var, tres);
			res += tres;
		}
		if (mdcell->sm.m != NULL) {
			MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var, tres);
			res += tres;
		}
		if (mdcell->pm.m != NULL) {
			MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var, tres);
			res += tres;
		}
		if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;
	}
	else Asyn_SPME_EF_Real(mdcell, cmm_cell, nThreads, var, tres_real, svar);

	// accumulate the electric field from the image part
	//vspme.reset_local_EF(); // reset the variable for E & F
	//init_normalized_coordinates(vspme, nThreads); -- vspme is initilized outside before this funciton is called
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
		EwaldSum_SurfaceDipole(cmm_cell->mu_surf, (EwaldSumRealVars0*)(&(var.esv1)), tres);
		res += tres;
	}

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VMUatom> *pvcell = (_VE_CELL<_VMUatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	if (nThreads > 1) {
		svar.thread_var.WaitUntilThreadOver();
		res += tres_real;
	}

	vspme.dump_EF(!var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
}

void SPME_EF(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_evatom_::_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE*, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&VMM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE*, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&VSM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE*, BASIC_CLUSTER, SPME_VAR, InteractRes>((void*)(&VPM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var, tres);
		res += tres;
	}
	if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;

	// accumulate the electric field from the image part
	//vspme.reset_local_EF(); // reset the variable for E & F
	//init_normalized_coordinates(vspme, nThreads); -- vspme is initilized outside before this funciton is called
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
		EwaldSum_SurfaceDipole(cmm_cell->mu_surf, (EwaldSumRealVars0*)(&(var.esv1)), tres);
		res += tres;
	}

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VMUatom> *pvcell = (_VE_CELL<_VMUatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	vspme.dump_EF(!var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
void SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_evatom_::_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	init_normalized_coordinates(vspme, nThreads); // vspme has to be initialized before SPME_EF for atom
	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res);

	if (var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}

void SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_evatom_::_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	init_normalized_coordinates(vspme, nThreads); // vspme has to be initialized before SPME_EF for atom
	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res);

	if (var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}

void dmax(double &x1, double &x2, double &res) {
	res = (x1 > x2 ? x1 : x2);
}
/*
void AtomicEFieldCorrect_SM(SMOLECULE *sm, int nSM) {
	EwaldSumRealVars1 esv1;
	esv1.bEwald = false; // direct interaction
	VECTOR3 E, dr;
	double R, R2;
	int im, ia, ia1;
	BASIC_CLUSTER *pc;
	MATOM *patom, *patom1;
	for (im = 0; im < nSM; im++) {
		pc = sm[im].c;
		for (ia = 0; ia < pc->nAtoms; ia++) {
			patom = pc->atom + ia;
			for (ia1 = 0; ia1 < pc->nAtoms; ia1++) {
				if (ia == ia1) continue;
				patom1 = pc->atom + ia1;
				VECT3(patom->r, patom1->r, dr)
				V3ABS2(dr, R2) R = sqrt(R2);
				esv1.init_vars(dr, R, R2);
				MTP_E(esv1, patom1->c, patom1->mu.v, E, false, 0);
				patom->E.v[0] += E.v[0];
				patom->E.v[1] += E.v[1];
				patom->E.v[2] += E.v[2];
			}
		}
	}
}
*/
/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
bool Polarize_SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_evatom_::_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
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
	while (bPolarize && iloop < max_loops) {
		SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, tres); // calculate the dipole only

		// add the electric field from the atom inside the cluster
		//if (mdcell->sm.n > 0) MOperate<SMOLECULE>((void*)(&AtomicEFieldCorrect_SM), mdcell->sm.m, mdcell->sm.n, nThreads);

		// polarization
		if (mdcell->mm.n > 0) MOperate<MMOLECULE>((void*)(&polarize_MM), mdcell->mm.m, mdcell->mm.n, nThreads);
		if (mdcell->sm.n > 0) MOperate<SMOLECULE>((void*)(&polarize_SM), mdcell->sm.m, mdcell->sm.n, nThreads);
		if (mdcell->pm.n > 0) MOperate<PMOLECULE>((void*)(&polarize_PM), mdcell->pm.m, mdcell->pm.n, nThreads);
		cal_dipole(*mdcell, nThreads, mdcell->mu_surf); // accumulate dipole of each atom, and the whole mdcell
		memcpy(&(cmm_cell->mu_surf), &(mdcell->mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
		dmu_max = 0; vt = 0;
		if (mdcell->ind_mu_mm.n > 0) MOperate2<MOL_POL<MMOLECULE>, double>((void*)(&MM_maxdiff_backup_polarization), mdcell->ind_mu_mm.m, mdcell->ind_mu_mm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_sm.n > 0) MOperate2<MOL_POL<SMOLECULE>, double>((void*)(&SM_maxdiff_backup_polarization), mdcell->ind_mu_sm.m, mdcell->ind_mu_sm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_pm.n > 0) MOperate2<MOL_POL<PMOLECULE>, double>((void*)(&PM_maxdiff_backup_polarization), mdcell->ind_mu_pm.m, mdcell->ind_mu_pm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;

		iloop++;
		if (dmu_max <= muConvg && iloop > 0) {status = true; break;}
	}

	//char msg[256] = "\0";
	//if (status) sprintf(msg, "Convergent polarized dipole is obtainted with %d loops", iloop); 
	//else sprintf(msg, "No convergent polarized dipole is obtained with %d loops", iloop);
	//show_log(msg, true);

	var.bEfieldOnly = bEfieldOnly;

	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res); // calculate the dipole and force

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}

	return status;
}

bool Polarize_SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_evatom_::_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
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
	while (bPolarize && iloop < max_loops) {
		SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, tres); // calculate the dipole only

		// add the electric field from the atom inside the cluster
		//if (mdcell->sm.n > 0) MOperate<SMOLECULE>((void*)(&AtomicEFieldCorrect_SM), mdcell->sm.m, mdcell->sm.n, nThreads);

		// polarization
		if (mdcell->mm.n > 0) MOperate<MMOLECULE*>((void*)(&polarize_VMM), mdcell->mm.m, mdcell->mm.n, nThreads);
		if (mdcell->sm.n > 0) MOperate<SMOLECULE*>((void*)(&polarize_VSM), mdcell->sm.m, mdcell->sm.n, nThreads);
		if (mdcell->pm.n > 0) MOperate<PMOLECULE*>((void*)(&polarize_VPM), mdcell->pm.m, mdcell->pm.n, nThreads);
		cal_dipole(*mdcell, nThreads, mdcell->mu_surf); // accumulate dipole of each atom, and the whole mdcell
		memcpy(&(cmm_cell->mu_surf), &(mdcell->mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
		dmu_max = 0; vt = 0;
		if (mdcell->ind_mu_mm.n > 0) MOperate2<MOL_POL<MMOLECULE>, double>((void*)(&MM_maxdiff_backup_polarization), mdcell->ind_mu_mm.m, mdcell->ind_mu_mm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_sm.n > 0) MOperate2<MOL_POL<SMOLECULE>, double>((void*)(&SM_maxdiff_backup_polarization), mdcell->ind_mu_sm.m, mdcell->ind_mu_sm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_pm.n > 0) MOperate2<MOL_POL<PMOLECULE>, double>((void*)(&PM_maxdiff_backup_polarization), mdcell->ind_mu_pm.m, mdcell->ind_mu_pm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;

		iloop++;
		if (dmu_max <= muConvg && iloop > 0) {status = true; break;}
	}

	//char msg[256] = "\0";
	//if (status) sprintf(msg, "Convergent polarized dipole is obtainted with %d loops", iloop); 
	//else sprintf(msg, "No convergent polarized dipole is obtained with %d loops", iloop);
	//show_log(msg, true);

	var.bEfieldOnly = bEfieldOnly;

	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res); // calculate the dipole and force

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}

	return status;
}

// polarized SPME, with the interaction between induced dipole excluded
// in some polarized model, the interaction between induced dipole is excluded, e.g. only q-induced dipole is considered.
void VSPME_F_exclude_induced(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_evatom_::_VMUatom>& vspme_induced, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	// accumulate the electric field from the image part
	//vspme_induced.reset_local_EF(); // reset the variable for E & F
	init_normalized_coordinates(vspme_induced, nThreads); //-- vspme is initilized outside before this funciton is called
	cal_Q_spme(vspme_induced, nThreads); // SPME Q matrix
	double U = vspme_induced.cal_U_recp(); //calculate the energy from the image part 
	res.U -= U;
	cal_E_recp(vspme_induced, true, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (var.bVirial) {
		cal_StressTensor(vspme_induced, nThreads);
		if (var.bST_diagonal) {
			res.STxx -= vspme_induced.STxx; res.STyy -= vspme_induced.STyy; res.STzz -= vspme_induced.STzz;
		}
		res.STtr -= vspme_induced.STtr;
	}

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VMUatom> *pvcell = (_VE_CELL<_VMUatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	//vspme_induced.dump_EF(!var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
	bool bDumpF = false; // useless
	if (nThreads < 2 || vspme_induced.va.n < N4PARALLEL) dump_F_exclude_induced_mu(vspme_induced, 0, vspme_induced.va.n - 1, bDumpF);
	else MRunSys2< VSPME<_evatom_::_VMUatom>, bool >((void*)(&dump_F_exclude_induced_mu), vspme_induced, 0, vspme_induced.va.n, nThreads, bDumpF);
}

bool Polarize_SPME_Interact2(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_evatom_::_VMUatom>& vspme, VSPME<_evatom_::_VMUatom>& vspme_induced, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
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
}

#endif // _IGNORE_ELECTROSTATIC_INTERACTION_

} // end of namespace _spme_

#define TEST 1
#if TEST == 1
using namespace _spme_;

double EwaldSum_Recip(MMOL_MD_CELL &mdcell) {
#define _NH 4
#define NH 9
#define AH_MIN 5e-3

	float AhMin = (float)1e-8;

	float V = mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8;
	double kappa, kappa3;
	double inv_3V;

	InitEwaldSumPars(V, inv_3V, kappa, kappa3);

	complex f;
	VECTOR3 h, dr;
	double inv_k2 = 0.25 / (kappa * kappa);
	double fcos0, fsin0, fcos1, fsin1, hR, hk, muh, H2, Ah;
	double t1, t2;
	double U_recip, f_recip;
	int nH = int(kappa * mdcell.h[0] * 2 + 0.5); if (nH < 1) nH = 1;
	int nc[3];
	U_recip = 0;
	int iatom, n, i;
	int m = 0, jatom;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;

	double EF = 0;

	for (n = 0; n < mdcell.sm.n; n++) {
		pc = mdcell.sm.m[n].c;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom;
			for (i = 0; i < 3; i++) {
				patom->E.v[i] = 0;
			}
		}
	}

	nH = 20;
	VECTOR3 r;
	for (nc[0] = -nH; nc[0] <= nH; nc[0]++) {for (nc[1] = -nH; nc[1] <= nH; nc[1]++) {for (nc[2] = -nH; nc[2] <= nH; nc[2]++) {
	//for (nc[0] = 0; nc[0] <= nH; nc[0]++) {for (nc[1] = 0; nc[1] <= nH; nc[1]++) {for (nc[2] = 0; nc[2] <= nH; nc[2]++) {
		if (nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
		//if ((nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2]) > nH * nH) continue;
		for (i = 0; i < 3; i++) h.v[i] = PI / mdcell.h[i] * nc[i];
		V3ABS2(h, H2)
		hk = H2 * inv_k2;
		//EXP(Ah, hk, i)
		Ah = exp(-hk);
		//if (Ah < 1e-8) continue;
		Ah *= PI2 / V / H2;

		f.x = 0; f.y = 0;
		for (n = 0; n < mdcell.mm.n; n++) {
			for (m = 0; m < mdcell.mm.m[n].nCluster; m++) {
				pc = mdcell.mm.m[n].cluster + m;
				for (iatom = 0; iatom < pc->nAtoms; iatom++) {
					patom = pc->atom + iatom;
					scale_uv3(h, patom->r, hR)
					PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
					f.x += patom->c * fcos0; f.y += patom->c * fsin0;
				}
			}
		}
		for (n = 0; n < mdcell.sm.n; n++) {
			pc = mdcell.sm.m[n].c;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom;
				for (i = 0; i < 3; i++) r.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
				scale_uv3(h, patom->r, hR)
				//PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
				fcos0 = cos(hR); fsin0 = sin(hR);
				f.x += patom->c * fcos0; f.y += patom->c * fsin0;
			}
		}
		U_recip += Ah * (f.x * f.x + f.y * f.y);

		// electric field
		for (n = 0; n < mdcell.sm.n; n++) {
			pc = mdcell.sm.m[n].c;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom;
				for (i = 0; i < 3; i++) r.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
				scale_uv3(h, patom->r, hR)
				//PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
				fcos1 = cos(hR); fsin1 = sin(hR);

				EF =  fsin1 * f.x - fcos1 * f.y;
				for (i = 0; i < 3; i++) {
					patom->E.v[i] += 2 * Ah * h.v[i] * EF;
				}
			}
		}
	}}}


	return U_recip * eUnit_estat_kT;
}

double EwaldSum_Recip_dipole(MMOL_MD_CELL &mdcell) {
#define _NH 4
#define NH 9
#define AH_MIN 5e-3

	float AhMin = (float)1e-8;

	float V = mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8;
	double kappa, kappa3;
	double inv_3V;

	InitEwaldSumPars(V, inv_3V, kappa, kappa3);

	complex f;
	VECTOR3 h, dr;
	double inv_k2 = 0.25 / (kappa * kappa);
	double fcos0, fsin0, fcos1, fsin1, hR, hk, muh, H2, Ah;
	double t1, t2;
	double U_recip, f_recip;
	int nH = int(kappa * mdcell.h[0] * 2 + 0.5); if (nH < 1) nH = 1;
	int nc[3];
	U_recip = 0;
	int iatom, n, i;
	int m = 0, jatom;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;

	double EF = 0;

	for (n = 0; n < mdcell.sm.n; n++) {
		pc = mdcell.sm.m[n].c;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom;
			for (i = 0; i < 3; i++) {
				patom->E.v[i] = 0;
			}
		}
	}

	nH = 20;
	VECTOR3 r;
	for (nc[0] = -nH; nc[0] <= nH; nc[0]++) {for (nc[1] = -nH; nc[1] <= nH; nc[1]++) {for (nc[2] = -nH; nc[2] <= nH; nc[2]++) {
	//for (nc[0] = 0; nc[0] <= nH; nc[0]++) {for (nc[1] = 0; nc[1] <= nH; nc[1]++) {for (nc[2] = 0; nc[2] <= nH; nc[2]++) {
		if (nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
		//if ((nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2]) > nH * nH) continue;
		for (i = 0; i < 3; i++) h.v[i] = PI / mdcell.h[i] * nc[i];
		V3ABS2(h, H2)
		hk = H2 * inv_k2;
		//EXP(Ah, hk, i)
		Ah = exp(-hk);
		//if (Ah < 1e-8) continue;
		Ah *= PI2 / V / H2;

		f.x = 0; f.y = 0;
		for (n = 0; n < mdcell.mm.n; n++) {
			for (m = 0; m < mdcell.mm.m[n].nCluster; m++) {
				pc = mdcell.mm.m[n].cluster + m;
				for (iatom = 0; iatom < pc->nAtoms; iatom++) {
					patom = pc->atom + iatom;
					scale_uv3(h, patom->r, hR)
					PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
					f.x += patom->c * fcos0; f.y += patom->c * fsin0;
					f.x -= (patom->mu.v[0] * h.v[0] + patom->mu.v[1] * h.v[1] + patom->mu.v[2] * h.v[2])* fsin0;
					f.y += (patom->mu.v[0] * h.v[0] + patom->mu.v[1] * h.v[1] + patom->mu.v[2] * h.v[2])* fcos0;
				}
			}
		}
		for (n = 0; n < mdcell.sm.n; n++) {
			pc = mdcell.sm.m[n].c;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom;
				for (i = 0; i < 3; i++) r.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
				scale_uv3(h, patom->r, hR)
				//PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
				fcos0 = cos(hR); fsin0 = sin(hR);
				f.x += patom->c * fcos0; f.y += patom->c * fsin0;
				f.x -= (patom->mu.v[0] * h.v[0] + patom->mu.v[1] * h.v[1] + patom->mu.v[2] * h.v[2])* fsin0;
				f.y += (patom->mu.v[0] * h.v[0] + patom->mu.v[1] * h.v[1] + patom->mu.v[2] * h.v[2])* fcos0;
			}
		}
		U_recip += Ah * (f.x * f.x + f.y * f.y);

		// electric field
		for (n = 0; n < mdcell.sm.n; n++) {
			pc = mdcell.sm.m[n].c;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom;
				for (i = 0; i < 3; i++) r.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
				scale_uv3(h, patom->r, hR)
				//PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
				fcos1 = cos(hR); fsin1 = sin(hR);

				EF =  fsin1 * f.x - fcos1 * f.y;
				EF -= patom->mu.v[0] * h.v[0] + patom->mu.v[1] * h.v[1] + patom->mu.v[2] * h.v[2]; 
				for (i = 0; i < 3; i++) {
					patom->E.v[i] += 2 * Ah * h.v[i] * EF;
				}
			}
		}
	}}}

	for (n = 0; n < mdcell.sm.n; n++) {
			pc = mdcell.sm.m[n].c;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom;
				for (i = 0; i < 3; i++) {
					patom->E.v[i] += 4 * PI / V / 3 * patom->mu.v[i];
				}
			}
		}


	return U_recip * eUnit_estat_kT;
}

extern void get_atoms(MMOL_MD_CELL &mdcell, ARRAY<MATOM*> &atom);
extern void get_bclusters(MMOL_MD_CELL &mdcell, ARRAY<BASIC_CLUSTER*>& bcluster);

void test_SPME_EwaldSum() {
	int nThreads = MAX_THREADS;

	using namespace _cmm_3d_;

	::eps_dc = 1;
	MMOL_MD_CELL mdcell;
	//mdcell.SetMol(0, 40, 0);
	mdcell.SetMol(0, 5, 0);
	CMM_CELL3D<BASIC_CLUSTER> cmm;
	::mlog.init("test.log", false);

	//float h[3] = {20, 20, 20};
	float h[3] = {6, 6, 6};
	cmm.set_cell(10, 10, 10);
	cmm.set_cell_range(true, -h[0], h[0], -h[1], h[1], -h[2], h[2]);
	mdcell.h[0] = h[0]; mdcell.h[1] = h[1]; mdcell.h[2] = h[2];  mdcell.period = true;
	cmm.set_cell_pos();
	int ncs = 10;
	CMM_array_set_storage<BASIC_CLUSTER>(&cmm, ncs);

	int imol = 0, ia;
	SMOLECULE cw;
	VECTOR3 dr, axis;

	extern bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol);
	extern bool cp(SMOLECULE *d, SMOLECULE *s, bool setup = true);
	extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
	extern void init_LJ_Pars();
	extern void construct_LJ12_6_prof(int mode);
	extern double cal_charge(MMOL_MD_CELL &mmcell);

	::bPolarize = true; // polarizable water molecule

	ConstructSimpleMolecule(&cw, "SPC");
	cw.c->cal_geometry_radius();
	dr.v[0] = -cw.r.v[0]; dr.v[1] = -cw.r.v[1]; dr.v[2] = -cw.r.v[2];
	cw.shiftMol(dr);
	::iFormat_LJ = 1;  // Amber
	init_LJ_Pars();
	construct_LJ12_6_prof(::iFormat_LJ);

	double Ep = 0;
	
	//cmm_init_subcell<BASIC_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (h[2]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		mdcell.sm.m[imol].shiftMol(dr);
		rand_vect3(axis.v);
		rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI, true);

		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
	get_atoms(mdcell, mdcell.atom);
	get_bclusters(mdcell, mdcell.bcluster);
	mdcell.check_eneutral();

	mdcell.set_atom_multipole((bPolarize ? 1 : 0));
	mdcell.setMolIndx();
	molecule_check_periodic_cell(mdcell);
	check_atom_rg(mdcell, nThreads);

	cmm_init_distribute_cluster(mdcell.sm.m, mdcell.sm.n, cmm, true); // add simple solid molecule into cmm

	cmm_check_cluster<BASIC_CLUSTER>(cmm);
	cmm_check<BASIC_CLUSTER>(cmm, MAX_THREADS);
	SMCheckRelshipInCell(mdcell.sm.m, mdcell.sm.n, cmm, MAX_THREADS);

	/*
	{
		for (int imol = 0; imol < mdcell.sm.n; imol++) {
			sprintf(errmsg, "m %d [%f, %f, %f]", imol, mdcell.sm.m[imol].r0.v[0], mdcell.sm.m[imol].r0.v[1], mdcell.sm.m[imol].r0.v[2]);
			show_log(errmsg, true);
			sprintf(errmsg, "m %d has %d of neighbors", imol, number<CMM_IMAGINE<BASIC_CLUSTER> >(mdcell.sm.m[imol].c->spmeRship.ch));
			show_log(errmsg, true);
		}
	}
	*/

	// calculat the total dipole of the whole cell
	cal_dipole(mdcell, 1, mdcell.mu_surf);

	// SPME parameters
	extern BSplineFunc<2> bsp;
	init_BSplineFunc(bsp, 6);

	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 32;
	::dim_spme[0] = 64; ::dim_spme[1] = 64; ::dim_spme[2] = 64;

	VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	vspme_q.bST = true; vspme_q.bST_diagonal = true;
	VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	vspme_mu.bST = true; vspme_mu.bST_diagonal = true;
	if (bPolarize) {
		mdcell.init_polarization_buff();
		mdcell.init_dipole_hist();
		mdcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		vspme_mu.bsp = &(bsp);
		init_spme(&mdcell, vspme_mu); // init the viral atoms
		vspme_mu.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_mu.xl[0] = cmm.xl[0]; vspme_mu.xl[1] = cmm.xl[1]; vspme_mu.xl[2] = cmm.xl[2];
		vspme_mu.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();
		vspme_mu.init_k();
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();
	}
	else { // charge only
		mdcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		vspme_q.bsp = &(bsp);
		init_spme(&mdcell, vspme_q); // init the viral atoms
		vspme_q.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_q.xl[0] = cmm.xl[0]; vspme_q.xl[1] = cmm.xl[1]; vspme_q.xl[2] = cmm.xl[2];
		vspme_q.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();
		vspme_q.init_k();
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}

	SPME_VAR spme_var;
	spme_var.esv1.bEwald = true; spme_var.esv1.init_EwaldSum(mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8, ::rcut_Ewd);
	spme_var.bVirial = true; spme_var.bST_diagonal = true;
	spme_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	spme_var.esv1.iSurface = 0; // cubic

	InteractRes ires;
	if (bPolarize) {
		if (!Polarize_SPME_Interact(&mdcell, &cmm, vspme_mu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else SPME_Interact(&mdcell, &cmm, vspme_q, nThreads, spme_var, ires);
	Ep = ires.U;

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);
	
	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-SPME-EwaldSum.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);



	if (bPolarize) Ep = vspme_mu.cal_U_recp();
	else Ep = vspme_q.cal_U_recp();
	sprintf(errmsg, "SPME Ewald Sum [reciprocal] : %f kT", Ep); show_log(errmsg, true); 


	if (bPolarize) Ep = EwaldSum_Recip_dipole(mdcell);
	else Ep = EwaldSum_Recip(mdcell);
	sprintf(errmsg, "Restrict Ewald Sum [reciprocal] : %f kT", Ep); show_log(errmsg, true); 

	//sprintf(errmsg, "Restrict Ewald Sum [reciprocal] efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	//show_log(errmsg, true);

	if (bPolarize) vspme_mu.release_plans();
	else vspme_q.release_plans();

	// explicit calculation
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;
	int i;

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		V6zero(mdcell.sm.m[imol].c->dyn->fc)
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			V3zero(mdcell.sm.m[imol].c->atom[ia].E)
			V3zero(mdcell.sm.m[imol].c->atom[ia].F)
			for (i = 0; i < mdcell.sm.m[imol].c->atom[ia].f.n; i++) {
				V3zero(mdcell.sm.m[imol].c->atom[ia].f.m[i])
			}
		}
	}

	extern void ExplicitEInteract(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, EInteractVar& evar, InteractRes &res);

	EInteractVar evar;
	evar.bVirial = true; evar.bST_diagonal = true; 
	evar.bEfieldOnly = false; evar.bTholePolarize = false; evar.tv.iTholeModel = 0;
	evar.esv1.bEwald = false; // direct interaction, not EwaldSum
	InteractRes tres;

	ires.reset();
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		ExplicitEInteract(cmm, mdcell.sm.m[imol].c, evar, tres);
		ires += tres;
	}
	Ep = ires.U;

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);

	sprintf(errmsg, "explicit calculation : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-explicit.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	
	::mlog.close();
}


#endif

