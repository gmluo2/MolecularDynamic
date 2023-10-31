/**************************************************************************************
    similiar to Interaction.cpp, except calculating cluster interaction force 
	between clusters by multipole-multipole calculation directly, without calculating
	the electric field from those multipole
***************************************************************************************/

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

#define INTERACT_DEF 1
#define _INTERACTION 0
#include "Interaction.h"
#include "Interaction1.h"

#include "var.h"

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

#if _DISABLE_ == 0
// calculate the interaction energy, Electric field, gradient electric field on each atom and LJ force
void EFIELD1_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, double &UTotal, bool &OVERLAP) {
	MATOM *patom = NULL;
	VECTOR3 E;
	DMATRIX<3> Eg;
	VECTOR3 dr, u;
	int i, j, k, l, na, na1;
	VECTOR3 r0, r1;
	double t1, t2, t3;

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
#if _EWALD_SUM == _MultiPole_EWALD_SUM
	CHAIN<ELECTRO_MULTIPOLE> *ch_emp = NULL;
	ELECTRO_MULTIPOLE *emp = NULL;
#endif

	RELSHIP_ITERATOR<BASIC_CLUSTER, LJ_RELSHIP> LJit;

	CHAIN<BASIC_CLUSTER> *ch_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	VECTOR3 mT1;
	DMATRIX<3> mT2, mD2;
	CMATRIX<3> mT3;
	QMATRIX<3> mT4;

	double R, R2, R3, R4, R5, R7, R9;
	VECTOR3 *h = NULL;

	bool bound = false, bound14 = false;

	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;

	double Uestat = 0, U_LJ = 0;
	UTotal = 0;

	bool bQuadrupole = false;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 1
	// ****** ignore the electric interaction between atoms *******
	goto _local_interact;
	// ************************************************************
#endif

	if (::local_relax) goto _local_interact; // do parent-child interaction only

#if _EWALD_SUM == _MultiPole_EWALD_SUM

	for (i = 0; i < 3; i++) r0.v[i] = pc->emp->r0.v[i];
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

		bQuadrupole = (R < (pc->emp->d + emp->d) * 10 ? true : false);

		if (bQuadrupole) {
			R9 = R7 * R2;
			MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2)
		}

		// energy
		t1 = (emp->q == 0 ? 0 : pc->emp->q * emp->q / R); // q-q
		for (i = 0; i < 3; i++) {
			t1 += (-pc->emp->q * emp->mu.v[i] + emp->q * pc->emp->mu.v[i]) * mT1.v[i]; // q-mu
			for (j = 0; j < 3; j++) {
				if (bQuadrupole) t1 += (pc->emp->q * emp->Q.m[i][j] + emp->q * pc->emp->Q.m[i][j]) * mT2.m[i][j] / 3; // q-quadrupole
				t1 += -pc->emp->mu.v[i] * emp->mu.v[j] * mT2.m[i][j]; // mu-mu
			}
			// ignore the quadrupole - quadrupole
		}
		Uestat += t1;

		// force
		for (i = 0; i < 3; i++) {
			force.v[i] += -pc->emp->q * emp->q * mT1.v[i];  // q-q
			for (j = 0; j < 3; j++) {
				force.v[i] += (pc->emp->q * emp->mu.v[j] - emp->q * pc->emp->mu.v[j]) * mT2.m[i][j]; //q-mu
				for (k = 0; k < 3; k++) {
					force.v[i] -= (pc->emp->q * emp->Q.m[j][k] + emp->q * pc->emp->Q.m[j][k]) * mT3.m[i][j][k] / 3; // q-Q
					force.v[i] += pc->emp->mu.v[j] * emp->mu.v[k] * mT3.m[i][j][k]; // mu-mu
					if (bQuadrupole) {
						for (l = 0; l < 3; l++) {
							force.v[i] += (emp->mu.v[j] * pc->emp->Q.m[k][l] - pc->emp->mu.v[j] * emp->Q.m[k][l]) * mT4.m[i][j][k][l] / 3; // mu-Q
						}
					}
				}
			}
		}
		ch_emp = ch_emp->next;
	}
	VECT3((*(pc->Op)), pc->emp->r0, r0)  // r0 = pc->Op ==> pc->emp.r0
	force.v[0] *= fUnit_estat_mass; force.v[1] *= fUnit_estat_mass; force.v[2] *= fUnit_estat_mass;
	V3PV3(r0, force, torque)
	for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
#endif //#if _EWALD_SUM == _MultiPole_EWALD_SUM

_local_interact:

	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		V3zero(E)
	#if _EWALD_SUM == _MultiPole_EWALD_SUM
		Mzero(Eg, i, j)

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 1
	// ****** ignore the electric interaction between atoms *******
	goto _EOVER_EFIELD_CMM;
	// ************************************************************
#endif

#if _ATOM_INDUCED_DIPOLE == 0
		if (FABS(patom->c) < 0.0001) goto _EOVER_EFIELD_CMM;
#endif
		ch_cluster = pc->ERship.ch_cluster;
		while (ch_cluster != NULL) {
			pc1 = ch_cluster->p;
			if (pc1 == pc) {ch_cluster = ch_cluster->next; continue;} // same cluster

			if (::local_relax) { // do parent-child interaction only
				if (pc1->parent != pc && pc->parent != pc1) {ch_cluster = ch_cluster->next; continue;}
			}

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
				if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
					BOUND(patom, patom1, i, bound) // bound atom
					if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
					if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
					else bound14 = false;
				}
				else {bound = false; bound14 = false;}

				if (bound || pc == pc1) { // bounding atoms or they are in the same cluster
				}
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
	#endif //_EWALD_SUM == _MultiPole_EWALD_SUM

		// LJ interaction
		if (patom->par->epsLJ > 0.0001) {
			//ch_imag_cluster = pc->LJRship.ch;
			//while (ch_imag_cluster != NULL) {
				//pc1 = ch_imag_cluster->p->p;
			LJit.set(&(pc->LJRship)); LJit.start_iterate();
			while (!LJit.eoi()) {
				imag_cluster = LJit.current(); if (imag_cluster == NULL) {LJit.next_iterate(); continue;}
				pc1 = imag_cluster->p;
				if (pc1 == pc) {
					//ch_imag_cluster = ch_imag_cluster->next; 
					LJit.next_iterate();
					continue;
				}
				if (::local_relax) { // do parent-child interaction only
					if (pc1->parent != pc && pc->parent != pc1) {
						//ch_imag_cluster = ch_imag_cluster->next; 
						LJit.next_iterate();
						continue;
					}
				}

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

					if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}
					if (bound) continue; // ignore bound interaction

					if (R2 < _R2min_Interact) continue; // the two atoms are too close

					pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					if (pLJ != NULL && pLJ->epsLJ > 0.0001) {
						if (R2 > pLJ->rcut2) continue;
						//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
						R = sqrt(R2); u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;
						LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, pLJ->rLJ2, R, u, i, t1, t3, t2, force)
						if (bound14) { // 1-4 interaction, LJ /= 2
							t2 *= A14_LJ;
							for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
						}
						//U_LJ += t2 * 0.5; // all clusters are in same system
						U_LJ += t2; // will be corrected later
						V3PV3(patom->r0, force, torque)
						for (i = 0; i < 3; i++) {pc->dyn->fc.v[i+3] += force.v[i]; pc->dyn->fc.v[i] += torque.v[i];}
					}
				}
				//ch_imag_cluster = ch_imag_cluster->next;
				LJit.next_iterate();
			}
		}
	}
	// total energy in kT
	UTotal = Uestat * eUnit_estat_kT + U_LJ;
}

void EFIELD1_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, double &Utotal, bool &OVERLAP) {
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;
	double U1 = 0;
	bool overlap = false; OVERLAP = false;
	Utotal = 0;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			EFIELD1_CMM(cell, pc, U1, overlap);
			if (overlap) OVERLAP = true;
			Utotal += U1;
			it.next_iterate();
		}
	}
	Utotal *= 0.5; // electrostatic and LJ are pairwise interaction
}

void EFIELD1_CMM_SINGLE_MM(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, double &Utotal, bool &OVERLAP) {
	int nc;
	BASIC_CLUSTER *pc = NULL;
	double U1 = 0;
	bool overlap = false; OVERLAP = false;
	Utotal = 0;
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = (BASIC_CLUSTER*)(mm.cluster + nc);
		EFIELD1_CMM(cell, pc, U1, overlap);
		if (overlap) OVERLAP = true;
		Utotal += U1;
	}
	Utotal *= 0.5; // electrostatic and LJ are pairwise interaction
}

void FreeCellCMMInteraction1_SingleMM(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, double &V, bool &OVERLAP, bool hard_sphere, int n1, int n2) { // V in kT
	int nc, na, i;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	float coeff = sqrt(::eps_dc);

	// calculate the electric field on each atom in these cells
	EFIELD1_CMM_SINGLE_MM(cell, mm, n1, n2, V, OVERLAP);

	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = (BASIC_CLUSTER*)(mm.cluster + nc);
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			if (patom->c != 0) {
				for (i = 0; i < 3; i++) {
					F.v[i] = patom->c * patom->E.v[i] * fUnit_estat_mass;
					if (bEext) {if (i == 2) F.v[i] += patom->c0 * Eext * fUnit_ext_mass;} // external electric field
#if _ATOM_INDUCED_DIPOLE == 1
					for (j = 0; j < 3; j++) {
						F.v[i] += patom->ind_mu.v[j] * patom->Eg.m[i][j] * fUnit_estat_mass;
					}
#endif
					pc->dyn->fc.v[3+i] += F.v[i];
				}
				V3PV3(patom->r0, F, torque)
				for (i = 0; i < 3; i++) pc->dyn->fc.v[i] += torque.v[i];
				if (bEext) {
					V -= patom->c0 * patom->r0.v[2] * Eext * U_EextFr / kT;
				}
			}
		}

		//if (bEext) {
			//V -= coeff * (pc->emp->mu.v[2] + pc->emp->q * pc->emp->r0.v[2]) * Eext * U_EextFr / kT; // external electric field
		//}
	}
}

#endif //_DISABLE_
