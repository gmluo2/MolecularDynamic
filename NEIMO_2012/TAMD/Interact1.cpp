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
#include "Interact2.h"
#include "var.h"

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

#if _DISABLE_ == 0
// calculate the interaction energy, Electric field, gradient electric field on each atom and LJ force
// the interaction is calculated between atoms
void EFIELD(EwaldSumRealVars1 &esv, CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, InteractVar& var, InteractRes &res) {
	MATOM *patom = NULL;
	VECTOR3 E;
	DMATRIX<3> Eg;
	VECTOR3 dr, u;
	int i, j, k, l, na, na1;
	VECTOR3 r0, r1;
	double t1, t2, t3;

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
//#if _EWALD_SUM == _MultiPole_EWALD_SUM
	//CHAIN<ELECTRO_MULTIPOLE> *ch_emp = NULL;
	//ELECTRO_MULTIPOLE *emp = NULL;
//#endif

	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;

	VECTOR3 mT1;
	DMATRIX<3> mT2, mD2;
	CMATRIX<3> mT3;
	QMATRIX<3> mT4;

	double R, R2, R3, R4, R5, R7, R9;
	int npc[3], nCell, NC = 6;
	VECTOR3 drCell;
	VECTOR3 *h = NULL;

	bool bound = false, bound14 = false, intrinsic_cell = true;

	double U_LJ = 0;
	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;

	double Uestat = 0, Uself = 0;

	V6zero(pc->dyn->fc)

	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;

		for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];

		V3zero(E) V3zero(patom->E)
		Mzero(Eg, i, j)
		// periodic cell iteration
		for (npc[0] = -NC; npc[0] <= NC; npc[0]++) {drCell.v[0] = npc[0] * cell.xd[0]; 
		for (npc[1] = -NC; npc[1] <= NC; npc[1]++) {drCell.v[1] = npc[1] * cell.xd[1]; 
		for (npc[2] = -NC; npc[2] <= NC; npc[2]++) {drCell.v[2] = npc[2] * cell.xd[2];
			intrinsic_cell = ((npc[0] == 0 && npc[1] == 0 && npc[2] == 0) ? true : false);
			// image cell which is far from the cluster
			for (nCell = 0; nCell < cell.acell.n; nCell++) {
				pcell = cell.acell.m[nCell]; it.set(pcell); it.start_iterate();
				while (!it.eoi()) {
					pc1 = it.current(); if (pc1 == NULL) {it.next_iterate(); continue;}
					if (intrinsic_cell && pc1 == pc) { // same cluster in same cell
						it.next_iterate();
						continue;
					}
					for (na1 = 0; na1 < pc1->nAtoms; na1++) {
						patom1 = pc1->atom + na1;
						for (i = 0; i < 3; i++) {
							r1.v[i] = patom1->r.v[i] + pc1->dr_cmm.v[i];
							dr.v[i] = r0.v[i] - r1.v[i] - drCell.v[i];
						}
						V3ABS2(dr, R2)

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 1
						// ****** ignore the electric interaction between atoms *******
						goto _LJ_INTERACT1;
						// ************************************************************
#endif
#if _ATOM_INDUCED_DIPOLE == 0
						if (patom->c == 0) goto _LJ_INTERACT1; // ignore the electric calculation on this atom
#endif


#if _ATOM_INDUCED_DIPOLE == 0
						if (patom1->c == 0) goto _LJ_INTERACT1; // ignore the electric calculation with this atom
#endif
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

						if (bound || (pc == pc1 && intrinsic_cell)) { // bounding atoms, or inside the same cluster
							continue;
						}
						else if (bound14 && pc != pc1) { // 1-4 interaction
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
							t1 = patom->c * patom1->c / R * A14_E; // 1-4 interaction is corrected by * A14_E
#if _ATOM_INDUCED_DIPOLE == 1
							for (i = 0; i < 3; i++) {
								t1 += (-patom->c * patom1->ind_mu.v[i] + patom1->c * patom->ind_mu.v[i]) * 
									mT1.v[i] * A14_E;
								for (j = 0; j < 3; j++) t1 -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mT2.m[i][j] * A14_E;
							}
#endif
							Uestat += t1 * 0.5;
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
							t1 = patom->c * patom1->c / R;
#if _ATOM_INDUCED_DIPOLE == 1
							for (i = 0; i < 3; i++) {
								t1 += (patom1->c * patom->ind_mu.v[i] - patom1->ind_mu.v[i] * patom->c) * mT1.v[i];
								for (j = 0; j < 3; j++) t1 -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mT2.m[i][j];
							}
#endif
							Uestat += t1 * 0.5;
							USelf -= patom->c * patom1->c / R * ::eUnit_estat_kT * 0.5;
						}

_LJ_INTERACT1:
						// LJ interaction
						//if (R2 > r2cut_LJ || patom->par->epsLJ < 0.0001) continue;

						pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
						if (pLJ != NULL && pLJ->epsLJ > 0.0001) {
							if (bound) continue; // ignore bound interaction

							//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
							R = sqrt(R2); if (R > rcut_LJ) continue;
							u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;
							LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, pLJ->rLJ2, R, u, i, t1, t3, t2, force)
							if (bound14 && pc != pc1) {
								t2 *= A14_LJ; //1-4 LJ interaction / 2
								for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
							}
							U_LJ += t2 * 0.5;
							V3PV3(patom->r0, force, torque)
							for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
							USelf -= R * t3 * 0.5;
						}
					}
					it.next_iterate();
				}
			}
		}}} // end of periodic cell iteration

		for (i = 0; i < 3; i++) {
			patom->E.v[i] = E.v[i];
#if _ATOM_INDUCED_DIPOLE == 1
			for (j = 0; j < 3; j++) patom->Eg.m[i][j] = Eg.m[i][j];
#endif
			force.v[i] = patom->c * E.v[i] * ::fUnit_estat_mass;
		}
		V3PV3(patom->r0, force, torque);
		for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
	}
	// total energy in kT
	UTotal = Uestat * eUnit_estat_kT + U_LJ;
}

// calculate the interaction energy, Electric field, gradient electric field on each atom and LJ force
void EFIELD_CMM_1(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, double &UTotal, double &USelf, bool &OVERLAP) {
	MATOM *patom = NULL;
	VECTOR3 E;
	DMATRIX<3> Eg;
	VECTOR3 dr, u;
	int i, j, k, l, na, na1;
	VECTOR3 r0, r1;
	double t1, t2, t3;

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, LJ_RELSHIP> LJit;

#if _EWALD_SUM == _MultiPole_EWALD_SUM
	//CHAIN<ELECTRO_MULTIPOLE> *ch_emp = NULL;
	//ELECTRO_MULTIPOLE *emp = NULL;
#endif

	CHAIN<BASIC_CLUSTER> *ch_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	VECTOR3 mT1;
	DMATRIX<3> mT2, mD2;
	CMATRIX<3> mT3;
	QMATRIX<3> mT4;

	double R, R2, R3, R4, R5, R7, R9;
	int npc[3], nCell, NC = 8;
	VECTOR3 drCell;
	VECTOR3 *h = NULL;

	bool bound = false, bound14 = false, intrinsic_cell = true;

	double U_LJ = 0;
	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;

	double Uestat = 0, Uself = 0;
	UTotal = 0; USelf = 0;

	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		Vzero(E, i)
		Mzero(Eg, i, j)

	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 1
		// ****** ignore the electric interaction between atoms *******
		goto _EOVER_EFIELD_CMM_1;
		// ************************************************************
	#endif

	#if _ATOM_INDUCED_DIPOLE == 0
		if (FABS(patom->c) < 0.00001) goto _EOVER_EFIELD_CMM_1;
	#endif

	#if _EWALD_SUM == _MultiPole_EWALD_SUM
		// periodic cell iteration
		for (npc[0] = -NC; npc[0] <= NC; npc[0]++) {drCell.v[0] = npc[0] * cell.xd[0]; 
		for (npc[1] = -NC; npc[1] <= NC; npc[1]++) {drCell.v[1] = npc[1] * cell.xd[1]; 
		for (npc[2] = -NC; npc[2] <= NC; npc[2]++) {drCell.v[2] = npc[2] * cell.xd[2];
		if (npc[0] > 1 || npc[0] < -1 || npc[1] > 1 || npc[1] < -1 || npc[2] > 1 || npc[2] < -1) {
			// image cell which is far from the cluster
			for (nCell = 0; nCell < cell.nCells; nCell++) {
				emp = cell.cell[nCell].emp;
				for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - emp->r0.v[i] - drCell.v[i];
				V3ABS2(dr, R2)
				R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2; R9 = R7 * R2;
				MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
				MTP_T3(mT3, dr, R2, R7, i, j, k, t1) MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2)

				// energy
				if (emp->q != 0) Uestat += patom->c * emp->q / R; // q-q
				for (i = 0; i < 3; i++) {
					Uestat += -patom->c * mT1.v[i] * emp->mu.v[i]; // q-mu
				#if _ATOM_INDUCED_DIPOLE == 1
					if (emp->q != 0) Uestat += patom->ind_mu.v[i] * mT1.v[i] * emp->q;  // mu-q
				#endif
					for (j = 0; j < 3; j++) {
						Uestat += patom->c * mT2.m[i][j] * emp->Q.m[i][j] / 3; // q-quadrupole
					#if _ATOM_INDUCED_DIPOLE == 1
						Uestat -= patom->ind_mu.v[i] * mT2.m[i][j] * emp->mu.v[j]; // mu-mu
						for (k = 0; k < 3; k++) Uestat += patom->ind_mu.v[k] * mT3.m[i][j][k] * emp->Q.m[i][j] / 3; // mu-quadrupole
					#endif
					}
				}
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
			}
		}
		else { // using relation-ship
			intrinsic_cell = (npc[0] == 0 && npc[1] == 0 && npc[2] == 0 ? true : false);
			ch_emp = pc->ERship.ch_mpl;
			while (ch_emp != NULL) {
				emp = ch_emp->p;
			
				for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - emp->r0.v[i] - drCell.v[i];
				V3ABS2(dr, R2)
				R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2; R9 = R7 * R2;
				MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
				MTP_T3(mT3, dr, R2, R7, i, j, k, t1) MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2)

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
				if (intrinsic_cell) Uself += t1 * eUnit_estat_kT;
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
				if (intrinsic_cell && pc1 == pc) { // same cluster in same cell
					ch_cluster = ch_cluster->next; continue;
				}
				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					patom1 = pc1->atom + na1;
				#if _ATOM_INDUCED_DIPOLE == 0
					if (FABS(patom1->c) < 0.00001) continue;
				#endif
					for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
					V3ABS2(dr, R2)
					R = sqrt(R2); R3 = R * R2; R4 = R2 * R2; R5 = R3 * R2; R7 = R5 * R2;
					MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
					MTP_T3(mT3, dr, R2, R7, i, j, k, t1)

					if (R2 < d2_14) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}

					if (intrinsic_cell && bound) { // bounding atoms 
					}
					else if (bound14 && pc != pc1) { // 1-4 interaction
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
						t1 = patom->c * patom1->c / R * A14_E; // 1-4 interaction is corrected by * A14_E
					#if _ATOM_INDUCED_DIPOLE == 1
						for (i = 0; i < 3; i++) {
							t1 += (-patom->c * patom1->ind_mu.v[i] + patom1->c * patom->ind_mu.v[i]) * 
								mT1.v[i] * A14_E;
							for (j = 0; j < 3; j++) t1 -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mT2.m[i][j] * A14_E;
						}
					#endif
						Uestat += t1;
						if (intrinsic_cell) Uself += t1 * eUnit_estat_kT;
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
						t1 = patom->c * patom1->c / R;
					#if _ATOM_INDUCED_DIPOLE == 1
						for (i = 0; i < 3; i++) {
							t1 += (patom1->c * patom->ind_mu.v[i] - patom1->ind_mu.v[i] * patom->c) * mT1.v[i];
							for (j = 0; j < 3; j++) t1 -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mT2.m[i][j];
						}
					#endif
						Uestat += t1;
						if (intrinsic_cell) Uself += t1 * eUnit_estat_kT;
					}
				}
				ch_cluster = ch_cluster->next;
			}
			}
		}}} // end of periodic cell iteration
	#endif // #if _EWALD_SUM == _MultiPole_EWALD_SUM

_EOVER_EFIELD_CMM_1:

		for (i = 0; i < 3; i++) {
			patom->E.v[i] = E.v[i];
#if _ATOM_INDUCED_DIPOLE == 1
			for (j = 0; j < 3; j++) patom->Eg.m[i][j] = Eg.m[i][j];
#endif
		}

		// LJ interaction
		if (patom->par->epsLJ > 0.0001) {
			//ch_imag_cluster = pc->LJRship.ch;
			LJit.set(&(pc->LJRship)); LJit.start_iterate();
			//while (ch_imag_cluster != NULL) {
			while (!LJit.eoi()) {
				imag_cluster = LJit.current(); if (imag_cluster == NULL) {LJit.next_iterate(); continue;}
				pc1 = imag_cluster->p;
				if (pc1 == pc) {
					//ch_imag_cluster = ch_imag_cluster->next; 
					LJit.next_iterate();
					continue;
				} // ingnore the interaction inside cluster
				drCell.v[0] = imag_cluster->nx * cell.xd[0];
				drCell.v[1] = imag_cluster->ny * cell.xd[1];
				drCell.v[2] = imag_cluster->nz * cell.xd[2];

				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					patom1 = pc1->atom + na1;
					if (patom1->par->epsLJ < 0.0001) continue;
					for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
					V3ABS2(dr, R2)
					if (R2 > r2cut_LJ) continue;

					pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					if (pLJ != NULL && pLJ->epsLJ > 0.0001) {
						if (R2 > pLJ->rcut2) continue;
						if (R2 < d2_14) {
							BOUND(patom, patom1, i, bound) // bound atom
							if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
							if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
							else bound14 = false;
						}
						else {bound = false; bound14 = false;}
						if (bound) continue; // ignore bound interaction

						//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
						R = sqrt(R2); u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;
						LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, pLJ->rLJ2, R, u, i, t1, t3, t2, force)
						if (bound14) {
							t2 *= A14_LJ; //1-4 LJ interaction / 2
							for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
						}
						if (imag_cluster->local) t2 *= 0.5; // will be calculated twice
						U_LJ += t2;
						V3PV3(patom->r0, force, torque)
						for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
					}
				}
				//ch_imag_cluster = ch_imag_cluster->next;
				LJit.next_iterate();
			}
		}
	}
	// total energy in kT
	UTotal = Uestat * eUnit_estat_kT + U_LJ;
	USelf = Uself;
}

void EFIELD_CMM_1(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, double &Utotal, bool &OVERLAP) {
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;
	double Uself = 0, U1 = 0, U2 = 0;
	bool overlap = false; OVERLAP = false;
	Utotal = 0;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
#define _interact 0
#if _interact == 0 // using multipole method
			EFIELD_CMM_1(cell, pc, U1, U2, overlap);
#elif _interact == 1 // using all atom-atom interactions
			EFIELD(cell, pc, U1, U2, overlap);
#endif
#undef _interact 
			if (overlap) OVERLAP = true;
			Utotal += U1; Uself += U2;
			it.next_iterate();
		}
	}
	Utotal -= Uself * 0.5;
	return;
}

void CMMInteraction(CMM_CELL3D<BASIC_CLUSTER> &cell, double &V, bool &OVERLAP, bool hard_sphere, int n1Cell, int n2Cell) { // V in kT
	int nc, na, i;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;

	float coeff = sqrt(::eps_dc);

	// calculate the electric field on each atom in these cells
	EFIELD_CMM_1(cell, n1Cell, n2Cell, V, OVERLAP);

	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
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
			it.next_iterate();
		}
		//if (bEext) {
		//	V -= coeff * (cell.cell[nc].emp->mu.v[2] + cell.cell[nc].emp->q * cell.cell[nc].emp->r0.v[2]) * Eext * U_EextFr / kT; // external electric field
		//}
	}
}

#endif // _DISABLE_

// For periodical cell, explicitly calculate the electrostatic interaction energy, Electric field, and force on each atom
// the interaction is calculated between atoms
extern PAIRWISE_DB< PAIRWISE<float> > TholeRadius_db;
void ExplicitEInteract(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, EInteractVar& evar, InteractRes &res) {
	_EwaldSum_real_::EwaldSumRealVars1 *esv = &(evar.esv1);
	esv->bEwald = false; // direct interaction
	MATOM *patom = NULL;
	VECTOR3 E, F, dr, E_Thole, F_Thole;
	VECTOR3 Ei, Fi;
	int i, j, k, l, na, na1;
	VECTOR3 r0, r1;

	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;

	double R, R2, U, Ui;
	int npc[3], nCell, NC = 20;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, intrinsic_cell = true;

	double Uestat = 0, U_Thole = 0;

	char aindx[2] = {0x00, 0x00};
	unsigned int key = 0;
	PAIRWISE<float> *TholeRadius = NULL;
	double a_Thole = 0;

	res.reset();

	//for (na = 0; na < pc->nAtoms; na++) {
		//patom = pc->atom + na;
	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];

		V32V3(patom->rg, r0)

		// periodic cell iteration
		for (npc[0] = -NC; npc[0] <= NC; npc[0]++) {drCell.v[0] = npc[0] * cell.xd[0]; 
		for (npc[1] = -NC; npc[1] <= NC; npc[1]++) {drCell.v[1] = npc[1] * cell.xd[1]; 
		for (npc[2] = -NC; npc[2] <= NC; npc[2]++) {drCell.v[2] = npc[2] * cell.xd[2];
			intrinsic_cell = ((npc[0] == 0 && npc[1] == 0 && npc[2] == 0) ? true : false);
			// image cell which is far from the cluster
			for (nCell = 0; nCell < cell.acell.n; nCell++) {
				pcell = cell.acell.m[nCell]; it.set(pcell); it.start_iterate();
				while (!it.eoi()) {
					pc1 = it.current(); if (pc1 == NULL) {it.next_iterate(); continue;}
					if (intrinsic_cell && pc1 == pc) { // same cluster in same cell
						it.next_iterate();
						continue;
					}
					//for (na1 = 0; na1 < pc1->nAtoms; na1++) {
						//patom1 = pc1->atom + na1;
					for (na1 = 0; na1 < pc1->fatom.n; na1++) {
						patom1 = pc1->fatom.m[na1];
						V32V3(patom1->rg, r1)
						for (i = 0; i < 3; i++) {
							dr.v[i] = r0.v[i] - r1.v[i] - drCell.v[i];
						}
						V3ABS2(dr, R2)
						R = sqrt(R2);

						if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
							BOUND(patom, patom1, i, bound) // bound atom
							if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
							if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
							else bound14 = false;
						}
						else {bound = false; bound14 = false;}

						if (bound || (pc == pc1 && intrinsic_cell)) { // bounding atoms, or inside the same cluster
							continue;
						}

						if (patom->mMP == 0) {
							esv->EwaldSumRealVars0::init_vars(dr, R, R2);
							MTP_Interact(*esv, patom->c, patom1->c, E, F, U, false, 0);
						}
						else if (patom->mMP == 1) {
							esv->init_vars(dr, R, R2);
							MTP_Interact(*esv, patom->c, patom->mu.v, patom1->c, patom1->mu.v, E, F, U, false, 0);

							if (esv->iExIndDipole) {
								MTP_Interact(*esv, 0, patom->ind_mu.v, 0, patom1->ind_mu.v, Ei, Fi, Ui, false, 0);
								//V3minusV3(E, Ei, E) 
								V3minusV3(F, Fi, F)
								U -= Ui;
							}
						}

						if (bound14 && pc != pc1) { // 1-4 interaction
							U *= A14_E;
							// electric field
							for (i = 0; i < 3; i++) {
								E.v[i] *= A14_E; F.v[i] *= A14_E;
							}
						}

						if (evar.bTholePolarize) {
							aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
							TholeRadius = ::TholeRadius_db.search(key);
							if (TholeRadius != NULL) {
								a_Thole = TholeRadius->pv;
								if (evar.tv.init_vars(a_Thole, R, R2, dr.v)) {
									::_Thole_polarize_::TholeCorr(evar.tv, a_Thole, patom->c, patom->mu.v, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v, F_Thole.v, U_Thole);
									if (bound14) {
										E_Thole.v[0] *= A14_E; E_Thole.v[1] *= A14_E; E_Thole.v[2] *= A14_E;
										F_Thole.v[0] *= A14_E; F_Thole.v[1] *= A14_E; F_Thole.v[2] *= A14_E;
										U_Thole *= A14_E;
									}
									V3plusV3(E, E_Thole, E) V3plusV3(F, F_Thole, F) U += U_Thole;
								}
							}
						}

						Uestat += U;
						for (i = 0; i < 3; i++) {
							patom->E.v[i] += E.v[i];
							patom->F.v[i] += F.v[i] * ::fUnit_estat_mass;
						}
						if (evar.bVirial) Accumulate_Virial_Diagonal(res, F.v, dr.v);
					}

					it.next_iterate();
				}
			}
		}}} // end of periodic cell iteration

		//V3PV3(patom->r0, force, torque);
		//for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
	}
	// total energy in kT
	res.U = Uestat * ::eUnit_estat_kT;
	if (evar.bVirial) {
		if (evar.bST_diagonal) {
			res.STxx *= fUnit_estat_mass;
			res.STyy *= fUnit_estat_mass;
			res.STzz *= fUnit_estat_mass;
		}
		res.STtr *= fUnit_estat_mass;
	}
}

void CMM_Explicit_EInteract(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, EInteractVar& evar, InteractRes &res) {
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;
	InteractRes tres;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			ExplicitEInteract(cell, pc, evar, tres);
			res += tres;
			it.next_iterate();
		}
	}
}





