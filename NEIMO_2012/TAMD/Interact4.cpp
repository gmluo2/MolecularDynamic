/*********************************************************************
 Similiar to Interaction.cpp, while working on CG_MMOLECULE and CMM
**********************************************************************/

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
#include "ZMatrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "cluster.h"

#include "cg-mm.h"

#define INTERACT_DEF 1
#define _INTERACTION 0
#include "Interaction.h"

#include "var.h"

#include "Interact4.h"

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

namespace _coarse_grain_ {
extern CGBOUND_DB cg_non_neighbor_db;
extern CGBOUND_DB cg_intramol_non_neighbor_db;

void gen_add_random_force_cgmm(CG_MMOLECULE &mm, _rand_::RAND_ARRAY<int> &randa) {
	SVECTOR<double, 3> f;
	double d = 0;
	int i;
	CG_CLUSTER *pc = NULL;
	for (i = 0; i < randa.nw; i++) {
		if (randa.a.m[i] >=  mm.nCluster) break;
		pc = mm.cluster + randa.a.m[i];
		// random force and torque on each cluster
		gauss_random_force(randf, fwhm_randf, f.v, 1);
		//random_force(randf * 1.5f, t);
		//d = 1.5f * ranf();
		
		for (i = 0; i < 3; i++) pc->dyn->fc.v[i] += f.v[i];
	}
}

// calculate the interaction energy, Electric field, gradient electric field on each atom and LJ force
void CG_EFIELD_CMM(CMM_CELL3D<CG_CLUSTER> &cell, CG_CLUSTER *pc, double &UTotal, bool &OVERLAP) {
	CGATOM *patom = NULL, *patom1 = NULL;
	VECTOR3 r0, r1, dr, u;

	//CHAIN< CMM_IMAGINE<CG_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<CG_CLUSTER> *imag_cluster = NULL;
	CG_CLUSTER *pc1 = NULL;

	RELSHIP_ITERATOR< CG_CLUSTER, LOCAL_RELSHIP> LJit;

	int iBound = -1;
	bool bound = false, bound14 = false;

	unsigned int key = 0;
	char aindx[2];
	BOUND_PROF *pBoundProf = NULL;

	double U_LJ = 0;
	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;

	double eUnit_estat_kT = U_ElectForce / kT;

	double Uestat = 0;
	UTotal = 0;

	int i, j, k, l, na, na1;

	double R, R2, R3, R4, R5, R7, R9;
	double t1, t2, t3;

	VECTOR3 dr_sm; // dr in the same molecule
	double R2_sm = 0;

	/*
	VECTOR3 E;
	DMATRIX<3> Eg;

	CHAIN<ELECTRO_MULTIPOLE> *ch_emp = NULL;
	ELECTRO_MULTIPOLE *emp = NULL;

	VECTOR3 mT1;
	DMATRIX<3> mT2, mD2;
	CMATRIX<3> mT3;
	QMATRIX<3> mT4;

	VECTOR3 *h = NULL;
	*/

	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		Uestat = 0;
		for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		/*
		V3zero(E)
		Mzero(Eg, i, j)

#if _ATOM_INDUCED_DIPOLE == 0
		goto _EOVER_EFIELD_CMM;
#endif
		// using relation-ship
		if (local_interact) ch_emp = NULL; // ignore the interaction with multipoles
		else ch_emp = pc->ERship.ch_mpl;
		while (ch_emp != NULL) {
			emp = ch_emp->p;
			//for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - emp->r0.v[i];
			//ABS2(dr, i, R2)
			VECT3(emp->r0, r0, dr)
			scale_uv3(dr, dr, R2)

			R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2;
			MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
			MTP_T3(mT3, dr, R2, R7, i, j, k, t1) 
#if _ATOM_INDUCED_DIPOLE == 1
			R9 = R7 * R2;
			MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2)
#endif
			// energy
			t1 = (emp->q == 0 ? 0 : patom->c * emp->q / R); // q-q
			for (i = 0; i < 3; i++) {
				t1 += -patom->c * mT1.v[i] * emp->mu.v[i]; // q-mu
#if _ATOM_INDUCED_DIPOLE == 1
				if (emp->q != 0) t1 += patom->ind_mu.v[i] * mT1.v[i] * emp->q; // mu-q
#endif
				for (j = 0; j < 3; j++) {
					t1 += patom->c * mT2.m[i][j] * emp->Q.m[i][j] / 3; // q-quadrupole
#if _ATOM_INDUCED_DIPOLE == 1
					t1 -= patom->ind_mu.v[i] * mT2.m[i][j] * emp->mu.v[j]; // mu-mu
					for (k = 0; k < 3; k++) t1 += patom->ind_mu.v[k] * mT3.m[i][j][k] * emp->Q.m[i][j] / 3; // mu-quadrupole
#endif
				}
			}
			Uestat += t1;

			// electric field
			for (i = 0; i < 3; i++) {
				if (emp->q != 0) E.v[i] -= emp->q * mT1.v[i];
				for (j = 0; j < 3; j++) {
					E.v[i] += emp->mu.v[j] * mT2.m[i][j];
					for (k = 0; k < 3; k++) E.v[i] -= emp->Q.m[j][k] * mT3.m[i][j][k] / 3;
				}
			}
#if _ATOM_INDUCED_DIPOLE == 1
			// electric field gradient
			for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
				if (emp->q != 0) Eg.m[i][j] -= emp->q * mT2.m[i][j];
				for (k = 0; k < 3; k++) {
					Eg.m[i][j] += emp->mu.v[k] * mT3.m[i][j][k];
					for (l = 0; l < 3; l++) Eg.m[i][j] -= emp->Q.m[k][l] * mT4.m[i][j][k][l] / 3;
				}
			}}
#endif
			ch_emp = ch_emp->next;
		}

		ch_cluster = pc->ERship.ch_cluster;
		while (ch_cluster != NULL) {
			pc1 = ch_cluster->p;
			if (pc1 == pc) {ch_cluster = ch_cluster->next; continue;} // same cluster
			BASIC_CLUSTER_IN_COARSE_GRAIN(pc1, pc, same_coarse_grain, i)
			if (same_coarse_grain) {ch_cluster = ch_cluster->next; continue;} // in the same coarse-grain cluster

			for (na1 = 0; na1 < pc1->nAtoms; na1++) {
				patom1 = pc1->atom + na1;
#if _ATOM_INDUCED_DIPOLE == 0
				if (FABS(patom1->c) < 0.00001) continue; // this atom has no charge
#endif
				for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i];
				//ABS2(dr, i, R2)
				scale_uv3(dr, dr, R2)

				R = sqrt(R2); R3 = R * R2; 
#if _ATOM_INDUCED_DIPOLE == 1
				R4 = R2 * R2; R5 = R3 * R2; R7 = R5 * R2;
				MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
				MTP_T3(mT3, dr, R2, R7, i, j, k, t1)
#endif
				if (R2 < d2_14) {
					BOUND(patom, patom1, i, bound) // bound atom
					if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
					if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
					else bound14 = false;
				}
				else {bound = false; bound14 = false;}

				if (bound || pc == pc1) { // bounding atoms or they are in the same cluster
				}
				//else if (pc->parent == pc1 || pc1->parent == pc || bound14) {
				//}
				else if (bound14) { // 1-4 interaction
					// electric field
					for (i = 0; i < 3; i++) {
						t1 = patom1->c * dr.v[i] / R3;
						E.v[i] += t1 * A14_E;
#if _ATOM_INDUCED_DIPOLE == 1
						for (j = 0; j < 3; j++) {
							E.v[i] += patom1->ind_mu.v[j] * mT2.m[i][j] * A14_E;
						}
#endif
					}
#if _ATOM_INDUCED_DIPOLE == 1
					// electric field gradient
					for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
						Eg.m[i][j] -= patom1->c * mT2.m[i][j] * A14_E;
						for (k = 0; k < 3; k++) {
							Eg.m[i][j] += patom1->ind_mu.v[k] * mT3.m[i][j][k] * A14_E;
						}
					}}
#endif
					// energy
					Uestat += patom->c * patom1->c / R * A14_E; // 1-4 interaction is corrected by * A14_E
#if _ATOM_INDUCED_DIPOLE == 1
					for (i = 0; i < 3; i++) {
						Uestat += (-patom->c * patom1->ind_mu.v[i] + patom1->c * patom->ind_mu.v[i]) * 
							mT1.v[i] * A14_E;
						for (j = 0; j < 3; j++) Uestat -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mT2.m[i][j] * A14_E;
					}
#endif
				}
				else {
					// electric field
					for (i = 0; i < 3; i++) {
						E.v[i] += patom1->c * dr.v[i] / R3;
#if _ATOM_INDUCED_DIPOLE == 1
						for (j = 0; j < 3; j++) E.v[i] -= patom1->ind_mu.v[j] * mT2.m[i][j];
#endif
					}
#if _ATOM_INDUCED_DIPOLE == 1
					// electric field gradient
					for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
						Eg.m[i][j] -= patom1->c * mT2.m[i][j];
						for (k = 0; k < 3; k++) Eg.m[i][j] += patom1->ind_mu.v[k] * mT3.m[i][j][k];
					}}
#endif
					// energy
					Uestat += patom->c * patom1->c / R;
#if _ATOM_INDUCED_DIPOLE == 1
					for (i = 0; i < 3; i++) {
						Uestat += (patom1->c * patom->ind_mu.v[i] - patom1->ind_mu.v[i] * patom->c) * mT1.v[i];
						for (j = 0; j < 3; j++) Uestat -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mT2.m[i][j];
					}
#endif
				}
			}
			ch_cluster = ch_cluster->next;
		}

_EOVER_EFIELD_CMM:
		for (i = 0; i < 3; i++) {
			patom->E.v[i] = E.v[i];
#if _ATOM_INDUCED_DIPOLE == 1
			for (j = 0; j < 3; j++) patom->Eg.m[i][j] = Eg.m[i][j];
#endif
		}
		*/

#if _CG_NEIGHBOR_INTERACT_ == 0
		// LJ interaction
		if (patom->par->epsLJ > 0.0001) {
			ch_imag_cluster = pc->LRship.ch;
			while (ch_imag_cluster != NULL) {
				pc1 = ch_imag_cluster->p->p;
				if (pc1 == pc) {ch_imag_cluster = ch_imag_cluster->next; continue;}

				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					patom1 = pc1->atom + na1;
					if (patom1->par->epsLJ < 0.0001) continue;
					for (i = 0; i < 3; i++) {
						dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i];
					}
					//ABS2(dr, i, R2)
					R2 = dr.v[0] * dr.v[0]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[1] * dr.v[1]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[2] * dr.v[2]; if (R2 > r2cut_LJ) continue;

					if (R2 < _R2min_Interact) continue; // the two atoms are too close

					if (pc->mIndx == pc1->mIndx) {
						iBound = pc->cBound(pc1);
						if (bound = (iBound == 1 ? true : false)) continue; // ignore interaction between bound clusters
						bound14 = (iBound == 3 ? true : false);
					}

					pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					if (pLJ != NULL && pLJ->epsLJ > 0.0001) {
						//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
						R = sqrt(R2); u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;
						LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, pLJ->rLJ2, R, u, i, t1, t3, t2, force)
						/*
						if (bound14) { // 1-4 interaction, LJ /= 2
							t2 *= A14_LJ;
							for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
						}
						*/
						//U_LJ += t2 * 0.5; // all clusters are in same system
						U_LJ += t2; // will be corrected later
						//V3PV3(patom->r0, force, torque)
						//for (i = 0; i < 3; i++) {pc->dyn->fc.v[i+3] += force.v[i]; pc->dyn->fc.v[i] += torque.v[i];}
						V3plusV3(pc->dyn->fc, force, pc->dyn->fc)
					}
				}
				ch_imag_cluster = ch_imag_cluster->next;
			}
		}

#elif _CG_NEIGHBOR_INTERACT_ == 1
		//ch_imag_cluster = pc->LRship.ch;
		//while (ch_imag_cluster != NULL) {
			//pc1 = ch_imag_cluster->p->p;
		LJit.set(&(pc->LRship)); LJit.start_iterate();
		while (!LJit.eoi()) {
			imag_cluster = LJit.current();  if (imag_cluster == NULL) {LJit.next_iterate(); continue;}
			pc1 = imag_cluster->p;
			if (pc1 == pc && imag_cluster->local) {
				//ch_imag_cluster = ch_imag_cluster->next; 
				LJit.next_iterate();
				continue;
			}

			for (na1 = 0; na1 < pc1->nAtoms; na1++) {
				patom1 = pc1->atom + na1;
				for (i = 0; i < 3; i++) {
					dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i];
				}
				//ABS2(dr, i, R2)
				R2 = dr.v[0] * dr.v[0]; if (R2 > r2cut_LJ) continue;
				R2 += dr.v[1] * dr.v[1]; if (R2 > r2cut_LJ) continue;
				R2 += dr.v[2] * dr.v[2]; if (R2 > r2cut_LJ) continue;

				if (R2 < _R2min_Interact) continue; // the two atoms are too close

				aindx[0] = patom->aindx; aindx[1] = patom1->aindx;
				key = construct_char2_key_pairwise(aindx);
			#if _INTRA_INTER_INTERACTS == 1
				if (pc1->mIndx == pc->mIndx) {
					VECT3(patom1->r, patom->r, dr_sm)
					V3ABS2(dr_sm, R2_sm)
					R2_sm -= R2; R2_sm = (R2_sm >= 0 ? R2_sm : -R2_sm);
					if (R2_sm > 1) pBoundProf = cg_non_neighbor_db.search(key); // not in the same molecule
					else continue; // intra-molecule non-neighbor interaction will be calculated independently
					//else pBoundProf = cg_intramol_non_neighbor_db.search(key);
				}
				else pBoundProf = cg_non_neighbor_db.search(key);
			#else
				pBoundProf = cg_non_neighbor_db.search(key);
			#endif
				if (pBoundProf != NULL) {
					R = sqrt(R2); 
					if (R < pBoundProf->efprof.xt) {
						u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;
						interpolate2(t2, t3, (pBoundProf->efprof), R, i, t1)
						U_LJ += t2;
						pc->dyn->fc.v[0] += t3 * u.v[0];
						pc->dyn->fc.v[1] += t3 * u.v[1];
						pc->dyn->fc.v[2] += t3 * u.v[2];
					}
				}
			}
			//ch_imag_cluster = ch_imag_cluster->next;
			LJit.next_iterate();
		}
#endif

		// energy
		UTotal += Uestat * eUnit_estat_kT;
	}
	// total energy in kT
	UTotal += U_LJ;
}


// calculate L-J interaction energy and force between the atoms in same molecule
void cgmm_IntraMolInteract(CG_MMOLECULE *mm, int ncf, int nct, double &V) {
#if _INTRA_INTER_INTERACTS == 0
	V = 0; return;
#elif _INTRA_INTER_INTERACTS == 1
	CG_CLUSTER *pc = NULL, *pc1 = NULL;
	CGATOM *patom = NULL, *patom1 = NULL;
	VECTOR3 r0, r1, dr, u;

	int iBound = -1;
	unsigned int key = 0;
	char aindx[2];
	BOUND_PROF *pBoundProf = NULL;

	double U_LJ = 0;
	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;

	int i, j, k, l, nc, nc1;

	double R, R2;
	double t1, t2, t3;

	for (nc = ncf; nc <= nct; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		patom = pc->atom;
		for (nc1 = 0; nc1 < pc->non_neighbor_sm.n; nc1++) {
			pc1 = pc->non_neighbor_sm.m[nc1];
			patom1 = pc1->atom;
			VECT3(patom1->r, patom->r, dr)
			V3ABS2(dr, R2)
			if (R2 < _R2min_Interact) continue; // the two atoms are too close

			aindx[0] = patom->aindx; aindx[1] = patom1->aindx;
			key = construct_char2_key_pairwise(aindx);
			
			pBoundProf = cg_intramol_non_neighbor_db.search(key);
			if (pBoundProf != NULL) {
				R = sqrt(R2); 
				if (R < pBoundProf->efprof.xt) {
					u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;
					interpolate2(t2, t3, (pBoundProf->efprof), R, i, t1)
					U_LJ += t2; // energy in efprof is in kT
					pc->dyn->fc.v[0] += t3 * u.v[0];
					pc->dyn->fc.v[1] += t3 * u.v[1];
					pc->dyn->fc.v[2] += t3 * u.v[2];
				}
			}
		}
	}
	V = U_LJ;
#endif
}


void EFIELD_CMM_SINGLE_CGMM(CMM_CELL3D<CG_CLUSTER> &cell, CG_MMOLECULE &mm, int n1, int n2, double &Utotal, bool &OVERLAP) {
	int nc, n;
	CG_CLUSTER *pc = NULL;
	double U1 = 0;
	bool overlap = false; OVERLAP = false;
	Utotal = 0;
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = mm.cluster + nc;
		CG_EFIELD_CMM(cell, pc, U1, overlap);
		if (overlap) OVERLAP = true;
		Utotal += U1;
	}
	Utotal *= 0.5; // electrostatic and LJ are pairwise interaction
}

void FreeCellCMMInteraction_SingleCGMM(CMM_CELL3D<CG_CLUSTER> &cell, CG_MMOLECULE &mm, double &V, bool &OVERLAP, bool hard_sphere, int n1, int n2) { // V in kT
	int nc, na, i, j;
	CG_CLUSTER *pc = NULL;
	CGATOM *patom = NULL;
	VECTOR3 F;

	double unit_estat_mass = U_ElectForce / U_MassAccel;
	double unit_Eext_estat_force = U_EextForce / U_ElectForce;
	float coeff = sqrt(::eps_dc);

	// calculate the electric field on each atom in these cells
	EFIELD_CMM_SINGLE_CGMM(cell, mm, n1, n2, V, OVERLAP);

	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = mm.cluster + nc;
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			if (patom->c != 0) {
				for (i = 0; i < 3; i++) {
					F.v[i] = patom->c * patom->E.v[i];
					if (bEext) {if (i == 2) F.v[i] += patom->c0 * Eext * unit_Eext_estat_force;} // external electric field
#if _ATOM_INDUCED_DIPOLE == 1
					for (j = 0; j < 3; j++) {
						F.v[i] += patom->ind_mu.v[j] * patom->Eg.m[i][j];
					}
#endif
					//pc->dyn->fc.v[3+i] += F.v[i];
					pc->dyn->fc.v[i] += F.v[i];
				}
				//V3PV3(patom->r0, F, torque)
				//for (i = 0; i < 3; i++) pc->dyn->fc.v[i] += torque.v[i];
				if (bEext) {
					V -= patom->c0 * patom->r0.v[2] * Eext * U_EextFr / kT;
				}
			}
		}

		//for (i = 0; i < 6; i++) pc->dyn->fc.v[i] *= unit_estat_mass; // change unit
		for (i = 0; i < 3; i++) pc->dyn->fc.v[i] *= unit_estat_mass; // change unit

		//if (bEext) {
		//	V -= coeff * (pc->emp->mu.v[2] + pc->emp->q * pc->emp->r0.v[2]) * Eext * U_EextFr / kT; // external electric field
		//}
	}
}

void FreeCellCMMInteraction_CGSM(CMM_CELL3D<CG_CLUSTER> &cell, CGSM* sm, int nSM, double &V, bool &OVERLAP) { // V in kT
	int nm, na, i, j;
	CG_CLUSTER *pc = NULL;
	CGATOM *patom = NULL;
	VECTOR3 F;

	double unit_estat_mass = U_ElectForce / U_MassAccel;
	double unit_Eext_estat_force = U_EextForce / U_ElectForce;
	float coeff = sqrt(::eps_dc);

	for (nm = 0; nm <= nSM; nm++) {
		pc = (CG_CLUSTER*)(sm + nm);
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			if (patom->c != 0) {
				for (i = 0; i < 3; i++) {
					F.v[i] = patom->c * patom->E.v[i];
					if (bEext) {if (i == 2) F.v[i] += patom->c0 * Eext * unit_Eext_estat_force;} // external electric field
#if _ATOM_INDUCED_DIPOLE == 1
					for (j = 0; j < 3; j++) {
						F.v[i] += patom->ind_mu.v[j] * patom->Eg.m[i][j];
					}
#endif
					pc->dyn->fc.v[i] += F.v[i];
				}
				if (bEext) {
					V -= patom->c0 * patom->r0.v[2] * Eext * U_EextFr / kT;
				}
			}
		}

		for (i = 0; i < 3; i++) pc->dyn->fc.v[i] *= unit_estat_mass; // change unit

		//if (bEext) {
		//	V -= coeff * (pc->emp->mu.v[2] + pc->emp->q * pc->emp->r0.v[2]) * Eext * U_EextFr / kT; // external electric field
		//}
	}
}

void scale_force_cgcluster(CG_MMOLECULE *mm, int nc1, int nc2) {
	double f = 0;
	VECTOR3 force;
	int nc = 0;

	CG_CLUSTER *pc = NULL;
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;

		V32V3(pc->dyn->fc, force)
		V3ABS2(force, f)
		if (f > force2_max) {
			f = force_max / sqrt(f);
			force.v[0] *= f; force.v[1] *= f; force.v[2] *= f;

			V32V3(force, pc->dyn->fc)
		}
	}
}

void scale_force_cgsm(CGSM *sm, int nSM) {
	double f = 0;
	VECTOR3 force;
	int n = 0;

	CG_CLUSTER *pc = NULL;
	for (n = 0; n < nSM; n++) {
		pc = (CG_CLUSTER*)(sm + n);

		V32V3(pc->dyn->fc, force)
		V3ABS2(force, f)
		if (f > force2_max) {
			f = force_max / sqrt(f);
			force.v[0] *= f; force.v[1] *= f; force.v[2] *= f;

			V32V3(force, pc->dyn->fc)
		}
	}
}

void CG_Brownian_force(CG_MMOLECULE *mm, int nc1, int nc2) { // f0 -- total Brownian force at mass center
	//V6zero(f0)
	int nc = 0, i;
	VECTOR3 dr, torque;
	CG_CLUSTER *pc = NULL;

	double v2 = 0, v = 0;
	double *u = NULL;
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) return;// f0;
		pc = mm->cluster + nc;
		//TIME(pc->bkm.V0, (-ita_Brownian), pc->dyn->f_Brown, i)
		u = pc->bkm->v.v; v2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]; v = sqrt(v2);
		for (i = 0; i < 3; i++) {
			pc->dyn->f_Brown.v[i] = -pc->ita_Brownian * pc->bkm->v.v[i] - pc->cd_trans * pc->bkm->v.v[i] * v;
		}
		// we do not need to calculate the total force and torque at mass center
		//V32V3(mm->cm, (*(pc->r)), dr)
		//V3PV3(dr, pc->dyn->f_Brown, torque)
		//for (i = 0; i < 3; i++) {f.v[i] += torque.v[i]; f0.v[i + 3] += pc->dyn->f_Brown.v[i];}
	}
}

// Implicit Solvation Energy & Force
void CGClusterImplicitSolvationForce(CG_CLUSTER *pc, CMM_CELL3D<CG_CLUSTER> *cell) {
	//CHAIN< CMM_IMAGINE<CG_CLUSTER> > *ch = NULL;
	CMM_IMAGINE<CG_CLUSTER> *imag_cluster = NULL;
	CG_CLUSTER *bpc1 = NULL;

	RELSHIP_ITERATOR< CG_CLUSTER, IMPSOLV_RELSHIP<CG_CLUSTER> > ImpIt;
	
	double S0 = 0, S = 0, b = 0, cs = 0, fb;
	double radius0, radius1, rw = ::rSolvent, delta_r = 0;
	VECTOR3 f, dr, r0, r1;
	double d = 0, d2, d_cut = 0, r0_w = 0, t;
	int i, n_neighbors = 0;

	IMPSOLV_F *fpar = NULL;

	CLUSTER_IMPSOLV_PARS *solv_par = NULL, *solv1_par = NULL;

	solv_par = pc->psolv; radius0 = solv_par->r0; r0_w = radius0 + rw;
	double r0_w_2 = (radius0 + rw) * (radius0 + rw), r1_w_2 = 0;

	pc->E_ImplicitSolvation = 0;
	V3zero(f) 
	if (cell->periodic_cell) {
		V3plusV3(pc->rc_solv, pc->dr_cmm, r0) // r0 = rc_solv + dr_cmm of pc
	}
	else {
		V32V3(pc->rc_solv, r0)
	}
	//ch = pc->impsolv_rship.ch;
	//n_neighbors = number< CMM_IMAGINE<CG_CLUSTER> >(ch);
	n_neighbors = pc->impsolv_rship.nrship();

	S0 = 4 * PI * r0_w * r0_w;
	cs = 1;

	if (n_neighbors == 0) return;
	else {
		fpar = new IMPSOLV_F[n_neighbors];
		i = 0;
		cs = 1;
		ImpIt.set(&(pc->impsolv_rship)); ImpIt.start_iterate();
		while (!ImpIt.eoi()) {
			imag_cluster = ImpIt.current(); if (imag_cluster == NULL) {ImpIt.next_iterate(); continue;}
			bpc1 = imag_cluster->p;
			if (bpc1 == NULL) {
				fpar[i].fb = 0; fpar[i].fbp = 0;
				fpar[i].fb1 = 0; fpar[i].fbp1 = 0;
				//ch = ch->next; 
				ImpIt.next_iterate();
				i++; continue;
			}
			solv1_par = bpc1->psolv; radius1 = solv1_par->r0;
			d_cut = radius0 + radius1 + rw + rw;
			delta_r = radius1 - radius0;
			r1_w_2 = (radius1 + rw) * (radius1 + rw);

			if (cell->periodic_cell) {
				if (imag_cluster->local && bpc1 == pc) { // same cluster on itself, ignore it
					fpar[i].fb = 0; fpar[i].fbp = 0;
					fpar[i].fb1 = 0; fpar[i].fbp1 = 0;
					//ch = ch->next; 
					ImpIt.next_iterate();
					i++; continue;
				}

				V3plusV3(bpc1->rc_solv, bpc1->dr_cmm, r1)  // r1 = r_solv + dr_cmm of bpc1
				if (!imag_cluster->local) {
					r1.v[0] -= imag_cluster->nx * cell->xd[0];
					r1.v[1] -= imag_cluster->ny * cell->xd[1];
					r1.v[2] -= imag_cluster->nz * cell->xd[2];
				}
			}
			else {
				V32V3(bpc1->rc_solv, r1)
			}
			
			VECT3(r0, r1, dr) // dr = r1 - r0
			V3ABS2(dr, d2) d = sqrt(d2); 
			if (d < 0.1) {
				fpar[i].fb = 0; fpar[i].fbp = 0;
				fpar[i].fb1 = 0; fpar[i].fbp1 = 0;
				//ch = ch->next; 
				ImpIt.next_iterate();
				i++; continue;
			}
			fpar[i].u.v[0] = dr.v[0] / d; fpar[i].u.v[1] = dr.v[1] / d; fpar[i].u.v[2] = dr.v[2] / d;

			if (d < d_cut) {
				b = PI * r0_w * (d_cut - d) * (1 + delta_r / d);
				if (d >= FABS(delta_r)) fb = -PI * r0_w * ((d_cut - d) * delta_r / d2 + (1 + delta_r / d));
				else fb = -PI * r0_w * ((d_cut - d) / delta_r + 1 + sign(delta_r));
				cs *= (1 - b / S0);
				if (cs < 0) cs = 0;

				fb /= (S0 - b);
				fpar[i].fb = fb;
			}
			else {
				fpar[i].fb = 0;
			}
			i++;
			//ch = ch->next;
			ImpIt.next_iterate();
		}
		S = S0 * cs;
		pc->E_ImplicitSolvation = solv_par->sigma * S;

		V3zero(f)
		for (i = 0; i < n_neighbors; i++) {
			fpar[i].fb *= solv_par->f_sigma * S;
			f.v[0] -= fpar[i].fb * fpar[i].u.v[0];
			f.v[1] -= fpar[i].fb * fpar[i].u.v[1];
			f.v[2] -= fpar[i].fb * fpar[i].u.v[2];
		}
	}
	if (fpar != NULL) {delete[] fpar; fpar = NULL;}
	for (i = 0; i < 3; i++) {pc->dyn->fc.v[i+3] += f.v[i];}
}


// solvent effect between cluster in the same macro-molecules
void CMM_CGImplicitSolvationForce(CMM_CELL3D<CG_CLUSTER> &cell, int n1Cell, int n2Cell) {
	CG_CLUSTER *pc = NULL;
	CMM_CELL<CG_CLUSTER> *pcell = NULL;
	int i, j, k;

	ITERATOR<CG_CLUSTER, CMM_CELL<CG_CLUSTER> > it;

	for (i = n1Cell; i <= n2Cell; i++) {
		if (i >= cell.bcell.nx) break;
		for (j = 0; j < cell.bcell.ny; j++) {for (k = 0; k < cell.bcell.nz; k++) {
			pcell = &(cell.bcell.m[i].m[j].m[k]); it.set(pcell); it.start_iterate();
			while (!it.eoi()) {
				pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
				CGClusterImplicitSolvationForce(pc, &cell);
				it.next_iterate();
			}
		}}
	}
}

void CGMM_ImplicitSolvationForce(CMM_CELL3D<CG_CLUSTER> *cell, CG_MMOLECULE *mm, int nc1, int nc2) {
	CG_CLUSTER *pc = NULL;
	for (int nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		CGClusterImplicitSolvationForce(pc, cell);
	}
}


} // end of namespace _coarse_grain_
