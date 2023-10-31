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

//#include "NEIMO.h"
#include "MD.h"
//#include "cg-mm.h"
//#include "cg-md.h"
//#include "MD_POLYMER.h"
#include "MD_SM.h"
#include "md_cell.h"

#include "vmd_cell.h"

#include "EwaldSum.h"
#include "spme.h"
#include "spme_interact.h"

#include "Interact_gc.h"

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
//long t_real = 0, t_recp = 0, t_emp = 0, t_cluster = 0, t_LJ = 0;
#endif

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, _EwaldSum_real_::SURF_DIPOLE &mu_surf);
extern void cal_dipole(MMOL_VMD_CELL &mmcell, int nThreads, _EwaldSum_real_::SURF_DIPOLE &mu_surf);

extern PAIRWISE_DB< PAIRWISE<EF_DPROF> > sr_db; // pair-wise short-range interaction database

// pair-wise database about Thole-radius between ions in Thole polarization model, s = a * (alpha_i * alpha_j)^1/6
extern PAIRWISE_DB< PAIRWISE<float> > TholeRadius_db; //

//#error more work on fg, force on fractal coefficient, in all functions
// Interaction of Fractal Cluster with other clusters defined in Grand-Canonical Essemble 
// see J. Chem. Phys. 102, 925 (1995), by Chaomei Lo and Bruce Palmer
void GC_CLUSTER_LJ_INTERACT_ATOM_CORRECT(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, GC_InteractVar &var, GC_InteractRes &res) {
	MATOM *patom = NULL;
	VECTOR3 dr, u;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

	double rc_gc = var.r_c(); // additional distance for distance correction

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, LJ_RELSHIP> LJit;

	double R, R2, Rc, Uc;
	double t1, t2, t3, t2_c;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, local = true;

	double U_LJ = 0, epsLJ;
	VECTOR3 force, force_c, torque;
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
				//if (!imag_cluster->local) {LJit.next_iterate(); continue;}
				//if (pc1->mIndx == pc->mIndx && imag_cluster->local) {LJit.next_iterate(); continue;}
				if (imag_cluster->local && pc1->mIndx == pc->mIndx) {
					// same molecule, we do not do correction about the interaction in same molecule
					LJit.next_iterate();
					continue;
				}
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

					if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}
					if (bound) continue; // ignore bound interaction

					R = sqrt(R2); 
					u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;

					Rc = R + rc_gc;

					V3zero(force)
					if (patom->par->bLJ) pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					else pLJ = NULL;
					if (pLJ != NULL && pLJ->epsLJ > 0.0001) {
						if (R >= _Rmin_Interact) {
							LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, R, u, i, t1, t3, t2, force)
						}
						else {
							t2 = 0; V3zero(force)
						}
						epsLJ = pLJ->epsLJ * var.g;
						LJ_V3FORCE(epsLJ, pLJ->rLJ, Rc, u, i, t1, t3, t2_c, force_c)
						t2 = t2_c - t2; // energy correction
						V3minusV3(force_c, force, force_c) // force correction
						if (bound14 && pc != pc1) { //1-4 LJ interaction / 2
							t2 *= A14_LJ; t3 *= A14_LJ;
							for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
						}
						U_LJ += t2;
					}
					if (patom1->par->bSR) {
						aindx[1] = patom1->par->indx;
						key = construct_char2_key_pairwise(aindx);
						pw = sr_db.search(key);
						if (pw != NULL) {
							interpolate2(t1, t2, pw->pv, R, i, t3)
							interpolate2(Uc, t2_c, pw->pv, Rc, i, t3)
							Uc *= var.g; t2_c *= var.g;
							U_LJ += Uc - t1;
							t2 = t2_c - t2;
							if (t2 > force_max && scale_force) t2 = force_max; \
							force.v[0] += t2 * u.v[0]; force.v[1] += t2 * u.v[1]; force.v[2] += t2 * u.v[2];
						}
					}
					if (var.bVirial) {
						if (var.bST_diagonal) Accumulate_Virial_Diagonal_full(res, force.v, dr.v);
						else Accumulate_Virial_Trace_full(res, force.v, dr.v);
					}
					//V3PV3(patom->r0, force, torque)
					for (i = 0; i < 3; i++) {
						//pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];
						patom->F.v[i] += force.v[i]; // force has mass unit, was converted in LJ12_6 profile
						patom1->F.v[i] -= force.v[i];
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


void GC_CLUSTER_LJ_INTERACT_FATOM_CORRECT(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, GC_InteractVar &var, GC_InteractRes &res) {
	MATOM *patom = NULL;
	VECTOR3 dr, u;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

	double rc_gc = var.r_c(); // additional distance for distance correction

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, LJ_RELSHIP> LJit;

	double R, R2, Rc, Uc;
	double t1, t2, t3, t2_c;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, local = true;

	double U_LJ = 0, epsLJ;
	VECTOR3 force, torque, force_c;
	LJ_PARS *pLJ = NULL;

	PAIRWISE<EF_DPROF> *pw = NULL;
	char aindx[2];
	unsigned int key;

	res.reset();

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	//V6zero(pc->dyn->fc)

	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];
		V32V3(patom->rg, r0)

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

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
				//if (!imag_cluster->local) {LJit.next_iterate(); continue;}
				//if (pc1->mIndx == pc->mIndx && imag_cluster->local) {LJit.next_iterate(); continue;}
				if (imag_cluster->local && pc1->mIndx == pc->mIndx) {
					// same molecule, we do not do correction about the interaction in same molecule
					LJit.next_iterate();
					continue;
				}
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
					
					if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}
					if (bound) continue; // ignore bound interaction

					R = sqrt(R2); 
					u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;

					Rc = R + rc_gc;

					V3zero(force)
					if (patom1->par->bLJ) pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					else pLJ = NULL;
					if (pLJ != NULL && pLJ->epsLJ > 0.0001) {
						if (R >= _Rmin_Interact) {
							LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, R, u, i, t1, t3, t2, force)
						}
						else {
							t2 = 0; V3zero(force)
						}
						epsLJ = pLJ->epsLJ * var.g;
						LJ_V3FORCE(epsLJ, pLJ->rLJ, Rc, u, i, t1, t3, t2_c, force_c)
						t2 = t2_c - t2; // energy correction
						V3minusV3(force_c, force, force_c) // force correction
						if (bound14 && pc != pc1) { //1-4 LJ interaction / 2
							t2 *= A14_LJ; t3 *= A14_LJ;
							for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
						}
						U_LJ += t2;
					}
					if (patom1->par->bSR) {
						aindx[1] = patom1->par->indx;
						key = construct_char2_key_pairwise(aindx);
						pw = sr_db.search(key);
						if (pw != NULL) {
							interpolate2(t1, t2, pw->pv, R, i, t3)
							interpolate2(Uc, t2_c, pw->pv, Rc, i, t3)
							Uc *= var.g; t2_c *= var.g;
							U_LJ += Uc - t1;
							t2 = t2_c - t2;
							if (t2 > force_max && scale_force) t2 = force_max; \
							force.v[0] += t2 * u.v[0]; force.v[1] += t2 * u.v[1]; force.v[2] += t2 * u.v[2];
						}
					}
					if (var.bVirial) {
						if (var.bST_diagonal) Accumulate_Virial_Diagonal_full(res, force.v, dr.v);
						else Accumulate_Virial_Trace_full(res, force.v, dr.v);
					}
					//V3PV3(patom->r0, force, torque)
					for (i = 0; i < 3; i++) {
						//pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];
						patom->F.v[i] += force.v[i]; // force has mass unit, was converted in LJ12_6 profile
						patom1->F.v[i] -= force.v[i]; // force has mass unit, was converted in LJ12_6 profile
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

void GC_CLUSTER_LJ_INTERACT_CORRECT(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, GC_InteractVar &var, GC_InteractRes &res) {
	if (pc->fatom.m == NULL) return GC_CLUSTER_LJ_INTERACT_ATOM_CORRECT(cell, pc, var, res);
	else return GC_CLUSTER_LJ_INTERACT_FATOM_CORRECT(cell, pc, var, res);
}


void GC_LJ_Interact_Correct(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, GC_InteractVar &var, GC_InteractRes &res) {
	GC_InteractRes tres;
	res.reset();

	int im, ic;
	MMOLECULE *mm = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	for (im = 0; im < mdcell->fmm.n; im++) {
		mm = mdcell->mm.m[mdcell->fmm.m[im]];
		for (ic = 0; ic < mm->nCluster; ic++) {
			GC_CLUSTER_LJ_INTERACT_CORRECT(*cmm_cell, (BASIC_CLUSTER*)(mm->cluster + ic), var, tres);
			res += tres;
		}
	}
	for (im = 0; im < mdcell->fsm.n; im++) {
		sm = mdcell->sm.m[mdcell->fsm.m[im]];
		GC_CLUSTER_LJ_INTERACT_CORRECT(*cmm_cell, sm->c, var, tres);
		res += tres;
	}
	for (im = 0; im < mdcell->fpm.n; im++) {
		pm = mdcell->pm.m[mdcell->fpm.m[im]];
		GC_CLUSTER_LJ_INTERACT_CORRECT(*cmm_cell, pm->c, var, tres);
		res += tres;
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
void GC_CLUSTER_EFIELD_CORRECT(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, GC_SPME_VAR& var, GC_InteractRes& res) {
	using namespace _Thole_polarize_;
	EwaldSumRealVars1 *esv = &(var.esv1);
	TholeCorrTensor *tv = &(var.tv);

	esv->bEwald = false; // correction is only the direct interaction in short-range

	MATOM *patom = NULL;
	VECTOR3 E_real, F_real, E_Thole, F_Thole;
	VECTOR3 dr, u;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

	VECTOR3 dr_c, E_real_c, F_real_c, E_Thole_c, F_Thole_c;


	VECTOR3 dr1, E1_real, E1_Thole;
	VECTOR3 dr1_c, E1_real_c, E1_Thole_c;

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	CHAIN<BASIC_CLUSTER> *ch_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, SPME_RELSHIP> rshipIt;

	double R, R2, Rc, R2c;
	double t1, t2, t2_c;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, local = true, bEwaldCut = false;
	bool bEwaldKickoff = false, bSameCluster = false;

	double alpha_Thole = 0, a_Thole = ::_Thole_polarize_::a;

	double Utotal = 0, Uestat = 0, Uext = 0, U_real = 0, U_Thole = 0, Uc = 0, Uc_Thole = 0;

	char aindx[2] = {0x00, 0x00};
	unsigned int key;
	PAIRWISE<float> *TholeRadius = NULL;

	bool bCharge;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	res.reset();

	//char msg[256] = "\0";
	//sprintf(msg, "CLUSTER %d :", pc->mIndx); show_log(msg, true);

	//for (na = 0; na < pc->nAtoms; na++) {
		//patom = pc->atom + na;
	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];
		//if (!bEfieldOnly) {V3zero(patom->F)}
		if (patom->eNeutral) continue;
		Uestat = 0;
		//for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		V32V3(patom->rg, r0)

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
			if (local && pc1->mIndx == pc->mIndx) {
				// same molecule, we do not do correction about the interaction in same molecule
				rshipIt.next_iterate();
				continue;
			}
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

				R = sqrt(R2); var.set_dist(R); Rc = R + var.r_c; R2c = Rc * Rc;
				t1 = Rc / R; dr_c.v[0] = dr.v[0] * t1; dr_c.v[1] = dr.v[1] * t1; dr_c.v[2] = dr.v[2] * t1;

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
					// no correction is required -- inside same cluster
					continue;
				}
				else if (bound14 && pc != pc1) { // 1-4 interaction
					// no correction is required -- inside the same molecule
					continue;
				}
				else {
					/*
					for (i = 0; i < 3; i++) E_real.v[i] += f1 * patom1->c * dr.v[i] / R3;
					// energy
					t1 = patom1->c / R;
					Uestat += 0.5 * t1 * f0;
					*/
					//bEwaldKickoff = false; t1 = 0;
				}
				if (var.bEfieldOnly) {
					if (R >= _Rmin_Interact) {
						MTP_E(*esv, patom1, dr, R, R2, E_real, false, 0);
						if (var.bTholePolarize && !bSameCluster) {
							aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
							TholeRadius = ::TholeRadius_db.search(key);
							if (TholeRadius != NULL) {
								a_Thole = TholeRadius->pv;
								if (::_Thole_polarize_::TholeCorr_E(*tv, a_Thole, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v)) {
									V3plusV3(E_real, E_Thole, E_real)
								}
							}
						}
					}
					else {
						V3zero(E_real)
					}

					MTP_E(*esv, patom1, dr_c, Rc, R2c, E_real_c, false, 0);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr_E(*tv, a_Thole, patom1->c, patom1->mu.v, Rc, R2c, dr_c.v, E_Thole_c.v)) {
								/*
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
								}
								*/
								V3plusV3(E_real_c, E_Thole_c, E_real_c)
							}
						}
					}

					V3minusV3(E_real_c, E_real, E_real)
				}
				else {
					if (R >= _Rmin_Interact) {
						MTP_Interact(*esv, patom, patom1, dr, R, R2, E_real, F_real, U_real, false, 0, bExcludeInducedDipole);
						if (var.bTholePolarize && !bSameCluster) {
							aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
							TholeRadius = ::TholeRadius_db.search(key);
							if (TholeRadius != NULL) {
								a_Thole = TholeRadius->pv;
								if (::_Thole_polarize_::TholeCorr(*tv, a_Thole, patom->c, patom->mu.v, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v, F_Thole.v, U_Thole)) {
									V3plusV3(E_real, E_Thole, E_real) V3plusV3(F_real, F_Thole, F_real) U_real += U_Thole;
								}
							}
						}
					}
					else {
						V3zero(E_real) V3zero(F_real) U_real = 0;
					}

					MTP_Interact(*esv, patom, patom1, dr_c, Rc, R2c, E_real_c, F_real_c, Uc, false, 0, bExcludeInducedDipole);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr(*tv, a_Thole, patom->c, patom->mu.v, patom1->c, patom1->mu.v, Rc, R2c, dr_c.v, E_Thole_c.v, F_Thole_c.v, Uc_Thole)) {
								/*
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
									F_Thole.v[0] *= t2; F_Thole.v[1] *= t2; F_Thole.v[2] *= t2;
									U_Thole *= t2;
								}
								*/
								V3plusV3(E_real_c, E_Thole_c, E_real_c) V3plusV3(F_real_c, F_Thole_c, F_real_c) Uc += Uc_Thole;
							}
						}
					}

					V3minusV3(E_real_c, E_real, E_real)
					V3minusV3(F_real_c, F_real, F_real)
					U_real = Uc - U_real;
				}

				// electric field on atom 1
				if (R >= _Rmin_Interact) {
					dr1.v[0] = -dr.v[0]; dr1.v[1] = -dr.v[1]; dr.v[2] = -dr.v[2];
					MTP_E(*esv, patom, dr1, R, R2, E1_real, false, 0);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr_E(*tv, a_Thole, patom->c, patom->mu.v, R, R2, dr1.v, E1_Thole.v)) {
								V3plusV3(E1_real, E1_Thole, E1_real)
							}
						}
					}
				}
				else {V3zero(E1_real)}

				dr1_c.v[0] = -dr_c.v[0]; dr1_c.v[1] = -dr_c.v[1]; dr1_c.v[2] = -dr_c.v[2];
				MTP_E(*esv, patom, dr1_c, Rc, R2c, E1_real_c, false, 0);
				if (var.bTholePolarize && !bSameCluster) {
					aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
					TholeRadius = ::TholeRadius_db.search(key);
					if (TholeRadius != NULL) {
						a_Thole = TholeRadius->pv;
						if (::_Thole_polarize_::TholeCorr_E(*tv, a_Thole, patom->c, patom->mu.v, Rc, R2c, dr1_c.v, E1_Thole_c.v)) {
							V3plusV3(E1_real_c, E1_Thole_c, E1_real_c)
						}
					}
				}

				V3minusV3(E1_real_c, E1_real, E1_real)
				// end of E on atom 1


				patom->E.v[0] += E_real.v[0]; patom->E.v[1] += E_real.v[1]; patom->E.v[2] += E_real.v[2];
				patom1->E.v[0] += E1_real.v[0]; patom1->E.v[1] += E1_real.v[1]; patom1->E.v[2] += E1_real.v[2];
				if (!var.bEfieldOnly) {
					// the forece on atom has unit on mass unit
					F_real.v[0] *= fUnit_estat_mass; 
					F_real.v[1] *= fUnit_estat_mass; 
					F_real.v[2] *= fUnit_estat_mass;

					patom->F.v[0] += F_real.v[0]; 
					patom->F.v[1] += F_real.v[1]; 
					patom->F.v[2] += F_real.v[2];

					patom1->F.v[0] -= F_real.v[0]; 
					patom1->F.v[1] -= F_real.v[1]; 
					patom1->F.v[2] -= F_real.v[2];

					Uestat += U_real;
					if (var.bVirial) Accumulate_Virial_Diagonal_full(res, F_real.v, dr.v);
				}
			}
			//ch_imag_cluster = ch_imag_cluster->next;
			rshipIt.next_iterate();
		}
		/*
		//show_log("\n", true);
		// the electric field from the surface induced dipole
		for (i = 0; i < 3; i++) {
			t1 = esv->inv_3V * cell.mu_surf.mu.v[i];
			patom->E.v[i] -= t1;
			if (!var.bEfieldOnly) patom->F.v[i] -= patom->c * t1 * fUnit_estat_mass;
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

		// energy from the surface induced dipole 
		switch (patom->mMP) {
		case 0:
			Uestat += 0.5 * esv->inv_3V * patom->c * (patom->r0.v[0] * cell.mu_surf.mu.v[0] + patom->r0.v[1] * cell.mu_surf.mu.v[1] + patom->r0.v[2] * cell.mu_surf.mu.v[2]);
			break;
		case 1:
			Uestat += 0.5 * esv->inv_3V * ((patom->c * patom->r0.v[0] + patom->mu.v[0]) * cell.mu_surf.mu.v[0] 
			+ (patom->c * patom->r0.v[1] + patom->mu.v[1]) * cell.mu_surf.mu.v[1] 
			+ (patom->c * patom->r0.v[2] + patom->mu.v[2]) * cell.mu_surf.mu.v[2]);
			break;
		default:
			break;
		}
		*/
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

void GC_Electrostatic_Interact_Correct(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, GC_SPME_VAR &var, GC_InteractRes &res) {
	GC_InteractRes tres;
	res.reset();

	int im, ic;
	MMOLECULE *mm = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	for (im = 0; im < mdcell->fmm.n; im++) {
		mm = mdcell->mm.m[mdcell->fmm.m[im]];
		for (ic = 0; ic < mm->nCluster; ic++) {
			GC_CLUSTER_EFIELD_CORRECT(*cmm_cell, (BASIC_CLUSTER*)(mm->cluster + ic), var, tres);
			res += tres;
		}
	}
	for (im = 0; im < mdcell->fsm.n; im++) {
		sm = mdcell->sm.m[mdcell->fsm.m[im]];
		GC_CLUSTER_EFIELD_CORRECT(*cmm_cell, sm->c, var, tres);
		res += tres;
	}
	for (im = 0; im < mdcell->fpm.n; im++) {
		pm = mdcell->pm.m[mdcell->fpm.m[im]];
		GC_CLUSTER_EFIELD_CORRECT(*cmm_cell, pm->c, var, tres);
		res += tres;
	}
}

} // end of namespace _EwaldSum_real_

#endif  // _IGNORE_ELECTROSTATIC_INTERACTION_
