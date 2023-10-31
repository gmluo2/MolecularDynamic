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
#include "Interact3.h"

#include "var.h"

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

#if _DISABLE_ == 0

#if _ATOM_INDUCED_DIPOLE == 1
// calculating electric field from atomic charge and dipole
void EFIELD_EwaldSum_Cluster(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc) {
	MATOM *patom = NULL;
	VECTOR3 E_real, E_recip, E_14;
	VECTOR3 dr;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

	CHAIN<ELECTRO_MULTIPOLE> *ch_emp = NULL;
	ELECTRO_MULTIPOLE *emp = NULL;
	CHAIN<BASIC_CLUSTER> *ch_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	VECTOR3 mT1, mET1;
	DMATRIX<3> mT2, mET2, mD2;
	CMATRIX<3> mT3, mET3;
	//QMATRIX<3> mT4, mET4;

	double R, R2, R3, R4, R5, R7;//, R9;
	double kR, kR2;
	double fexp2, f0, f1, f2, f3, f4;//, f5, f6;
	double t1, t2, Ah, fcos0, fsin0, fcos1, fsin1, hR, fss, fcc, fsc, fcs;
	int npc[3];
	VECTOR3 drCell;
	VECTOR3 *h = NULL;

	bool bound = false, bound14 = false, intrinsic_cell = false;

	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		Vzero(E_real, i) Vzero(E_recip, i) Vzero(E_14, i)
		// periodic cell iteration
		for (npc[0] = -NPC; npc[0] <= NPC; npc[0]++) {drCell.v[0] = npc[0] * cell.xd[0]; 
		for (npc[1] = -NPC; npc[1] <= NPC; npc[1]++) {drCell.v[1] = npc[1] * cell.xd[1]; 
		for (npc[2] = -NPC; npc[2] <= NPC; npc[2]++) {drCell.v[2] = npc[2] * cell.xd[2];
			intrinsic_cell = (npc[0] == 0 && npc[1] == 0 && npc[2] == 0 ? true : false);
			ch_emp = pc->ERship.ch_mpl;
			while (ch_emp != NULL) {
				emp = ch_emp->p;
			
				for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - emp->r0.v[i] - drCell.v[i];
				if (FABS(dr.v[0]) > rcut_Ewd || FABS(dr.v[1]) > rcut_Ewd || FABS(dr.v[2]) > rcut_Ewd) {
					ch_emp = ch_emp->next; continue;
				}
				V3ABS2(dr, R2)
				if (R2 > r2cut_Ewd) {ch_emp = ch_emp->next; continue;}
				R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2;// R9 = R7 * R2;
				kR = kappa * R; kR2 = kR * kR;
				//fexp2 = exp(-kR2);
				Ewald_f04(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4)
				MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
				MTP_T3(mT3, dr, R2, R7, i, j, k, t1)

				EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j) 
				EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)

				// electric field
				for (i = 0; i < 3; i++) {
					if (emp->q != 0) E_real.v[i] -= emp->q * mET1.v[i];
					for (j = 0; j < 3; j++) {
						E_real.v[i] += emp->mu.v[j] * mET2.m[i][j];
						for (k = 0; k < 3; k++) E_real.v[i] -= emp->Q.m[j][k] * mET3.m[i][j][k] / 3;
					}
				}
				ch_emp = ch_emp->next;
			}

			ch_cluster = pc->ERship.ch_cluster;
			while (ch_cluster != NULL) {
				pc1 = ch_cluster->p;
				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					if (intrinsic_cell && pc1 == pc && na1 == na) continue; // same atom
					patom1 = pc1->atom + na1;
					for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
					if (FABS(dr.v[0]) > rcut_Ewd || FABS(dr.v[1]) > rcut_Ewd || FABS(dr.v[2]) > rcut_Ewd) {
						continue;
					}
					V3ABS2(dr, R2)
					if (R2 > r2cut_Ewd) continue;
					// adding the dipole interaction inside
					R = sqrt(R2); R3 = R * R2; R4 = R2 * R2; R5 = R3 * R2; kR = kappa * R; 
					EXP2(fexp2, kR, i) Ewald_f0(f0, kR, i) Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2)
					MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) //MTP_D2(mD2, dr, i, j)

					EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j)

					if (R2 < d2_14) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}

					if (bound || (intrinsic_cell && pc == pc1)) { // bounding atoms or they are in the same cluster
						for (i = 0; i < 3; i++) {
							// bound interaction has to be kicked off in Ewald Sum
							t1 = 1 - f1;
							E_real.v[i] -= patom1->c * mT1.v[i] * t1;
							for (j = 0; j < 3; j++) E_real.v[i] -= patom1->ind_mu.v[j] * (mT2.m[i][j] - mET2.m[i][j]);
						}
					}
					else if (bound14 && pc != pc1) { // 1-4 interaction
						for (i = 0; i < 3; i++) {
							E_real.v[i] += patom1->c * mET1.v[i];
							E_14.v[i] += (A14_E - 1) * patom1->c * mT1.v[i]; // 1-4 interaction will to be corrected later
							for (j = 0; j < 3; j++) {
								E_real.v[i] += patom1->ind_mu.v[j] * mET2.m[i][j];
								E_14.v[i] += (A14_E - 1) * patom1->ind_mu.v[j] * mT2.m[i][j];
							}
						}
					}
					else {
						for (i = 0; i < 3; i++) {
							E_real.v[i] += patom1->c * mET1.v[i];
							for (j = 0; j < 3; j++) E_real.v[i] += patom1->ind_mu.v[j] * mET2.m[i][j];
						}
					}
				}
				ch_cluster = ch_cluster->next;
			}
		}}} // end of periodic cell iteration

		for (npc[0] = 0; npc[0] < NH; npc[0]++) {for (npc[1] = 0; npc[1] < NH; npc[1]++) {for (npc[2] = 0; npc[2] < NH; npc[2]++) {
			Ah = ::Ah[npc[0]][npc[1]][npc[2]]; if (Ah < AhMin) continue;
			if (npc[0] == _NH && npc[1] == _NH && npc[2] == _NH) continue;
			ch_emp = pc->ERship.ch_mpl;
			h = &(::h[npc[0]][npc[1]][npc[2]]); 
			scale_uv3((*h), r0, hR)
			PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t2)
			while (ch_emp != NULL) {
				emp = ch_emp->p;
				fcos1 = emp->fcos.m[npc[0]][npc[1]][npc[2]]; fsin1 = emp->fsin.m[npc[0]][npc[1]][npc[2]];
				fcc = fcos0 * fcos1; fss = fsin0 * fsin1; fcs = fcos0 * fsin1; fsc = fsin0 * fcos1;
				for (i = 0; i < 3; i++) {
					t1 = 0;
					if (emp->q != 0) t1 += emp->q * (fsc - fcs);
					t1 -= emp->muh.m[npc[0]][npc[1]][npc[2]] * (fcc + fss);
					t1 += emp->Qh.m[npc[0]][npc[1]][npc[2]] * (fcs - fsc) / 3;
					E_recip.v[i] += Ah * h->v[i] * t1;
				}
				ch_emp = ch_emp->next;
			}

			ch_cluster = pc->ERship.ch_cluster;
			while (ch_cluster != NULL) {
				pc1 = ch_cluster->p;
				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					patom1 = pc1->atom + na1;
					if (patom1->c == 0) continue;
					r1.v[0] = patom1->r.v[0] + pc1->dr_cmm.v[0]; r1.v[1] = patom1->r.v[1] + pc1->dr_cmm.v[1]; r1.v[2] = patom1->r.v[2] + pc1->dr_cmm.v[2];
					scale_uv3((*h), r1, hR)
					PERIOD(hR, 0, PI2, t1) COS(fcos1, t1, i, t1) SIN(fsin1, t1, i, t2)
					fcc = fcos0 * fcos1; fss = fsin0 * fsin1; fcs = fcos0 * fsin1; fsc = fsin0 * fcos1;
					t1 = patom1->c * (fsc - fcs);
					scale_uv3(patom1->ind_mu, ::h[npc[0]][npc[1]][npc[2]], t2)
					t1 -= t2 * (fcc - fss);
					E_recip.v[i] += Ah * h->v[i] * t1;
				}
				scale_uv3(patom1->ind_mu, ::h[npc[0]][npc[1]][npc[2]], t1)
				E_recip.v[i] += Ah * h->v[i] * t1;
				ch_cluster = ch_cluster->next;
			}
		}}} // end of reciprocal sum
		for (i = 0; i < 3; i++) patom->E.v[i] = E_real.v[i] + E_recip.v[i] + E_14.v[i] - inv_V * cell.mu.v[i];
	}
}

void EFIELD_EwaldSum(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell) {
	int nc;
	CHAIN<BASIC_CLUSTER> *ch = NULL;
	BASIC_CLUSTER *pc = NULL;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		ch = cell.cell[nc].ch;
		while (ch != NULL) {
			pc = ch->p;
			EFIELD_EwaldSum_Cluster(cell, pc);
			ch = ch->next;
		}
	}
}

void Polarization(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell) {
	int nc, na;
	CHAIN<BASIC_CLUSTER> *ch = NULL;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		ch = cell.cell[nc].ch;
		while (ch != NULL) {
			pc = ch->p;
			for (na = 0; na < pc->nAtoms; na++) {
				patom = pc->atom + na;
				V3PC(patom->E, patom->par->alpha, patom->ind_mu)
			}
			ch = ch->next;
		}
	}
}

/***************************************************************
     Multi-Thread EFIELD Calculation based on EWALD-SUM & 
	             CELL MULTIPOLE METHOD
*****************************************************************/

void Efield_Polarization(CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads, int ncycles) {
	int ncycle = 0;
	for (ncycle = 0; ncycle < ncycles; ncycle++) {
		// function to calculate the electro-multipoles here
		CellOperate<BASIC_CLUSTER>((void*)(&EFIELD_EwaldSum), cell, nThreads); // calculate the electric field on each atoms
		CellOperate<BASIC_CLUSTER>((void*)(&Polarization), cell, nThreads); // calculate the induced dipoles of each atoms
		calColumbMultipoles<BASIC_CLUSTER>(cell, true, nThreads); // calculate the multipole again
	}
}


/********************************************************************
   INTERACTION CALCULATION on ATOM WITH CHARGE & INDUCED DIPOLE
*********************************************************************/

// calculate the interaction energy, Electric field, gradient electric field on each atom and LJ force
void Polar_EFIELD_EwaldSum_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, double &UTotal, double &USelf, bool &OVERLAP) {
	MATOM *patom = NULL;
	VECTOR3 E_real, E_recip, E_14;
	DMATRIX<3> Eg_real, Eg_recip, Eg_14;
	VECTOR3 dr, u;
	int i, j, k, l, na, na1;
	VECTOR3 r0, r1;

	CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CHAIN<ELECTRO_MULTIPOLE> *ch_emp = NULL;
	ELECTRO_MULTIPOLE *emp = NULL;
	CHAIN<BASIC_CLUSTER> *ch_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	VECTOR3 mT1, mET1;
	DMATRIX<3> mT2, mET2, mD2;
	CMATRIX<3> mT3, mET3;
	QMATRIX<3> mT4, mET4;

	double R, R2, R3, R4, R5, R7, R9;
	double kR, kR2;
	double fexp2, f0, f1, f2, f3, f4, f5, f6;
	double t1, t2, t3, Ah, fcos0, fsin0, fcos1, fsin1, fss, fcc, fsc, fcs, hR;
	int npc[3];
	VECTOR3 drCell;
	VECTOR3 *h = NULL;
	double muh0, muh1;

	bool bound = false, bound14 = false, intrinsic_cell = true, bEwaldCut = false;

	double U_LJ = 0;
	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;
	double unit_LJ_estat_force = U_LJForce / U_ElectForce;

	double Utotal = 0, Uself = 0;
	UTotal = 0; USelf = 0;

	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		Utotal = 0; Uself = 0;
		for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		Vzero(E_real, i) Vzero(E_recip, i) Vzero(E_14, i)
		Mzero(Eg_real, i, j) Mzero(Eg_recip, i, j) Mzero(Eg_14, i, j)
		// periodic cell iteration
		for (npc[0] = -NPC; npc[0] <= NPC; npc[0]++) {drCell.v[0] = npc[0] * cell.xd[0]; 
		for (npc[1] = -NPC; npc[1] <= NPC; npc[1]++) {drCell.v[1] = npc[1] * cell.xd[1]; 
		for (npc[2] = -NPC; npc[2] <= NPC; npc[2]++) {drCell.v[2] = npc[2] * cell.xd[2];
			intrinsic_cell = (npc[0] == 0 && npc[1] == 0 && npc[2] == 0 ? true : false);
			ch_emp = pc->ERship.ch_mpl;
			while (ch_emp != NULL) {
				emp = ch_emp->p;
				for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - emp->r0.v[i] - drCell.v[i];
				if (FABS(dr.v[0]) > rcut_Ewd || FABS(dr.v[1]) > rcut_Ewd || FABS(dr.v[2]) > rcut_Ewd) {
					bEwaldCut = true;
					if (!intrinsic_cell) {ch_emp = ch_emp->next; continue;}
				}
				V3ABS2(dr, R2)
				bEwaldCut = (R2 > r2cut_Ewd ? true : false);
				if (!intrinsic_cell && bEwaldCut) {ch_emp = ch_emp->next; continue;}
				R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2; 
				MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) 
				MTP_T3(mT3, dr, R2, R7, i, j, k, t1) 

				// self energy
				if (intrinsic_cell) {
					if (emp->q != 0) Uself += patom->c * emp->q / R; // q-q
					for (i = 0; i < 3; i++) {
						Uself += -patom->c * mT1.v[i] * emp->mu.v[i];
						if (emp->q != 0) Uself += patom->ind_mu.v[i] * mT1.v[i] * emp->q;  // q-mu + mu-q
						for (j = 0; j < 3; j++) {
							Uself += patom->c * mT2.m[i][j] * emp->Q.m[i][j] / 3; // q-quadrupole
							for (k = 0; k < 3; k++) Uself += patom->ind_mu.v[k] * mT3.m[i][j][k] * emp->Q.m[i][j] / 3; // mu-quadrupole
						}
					}
					if (bEwaldCut) {ch_emp = ch_emp->next; continue;}
				}

				R9 = R7 * R2;
				MTP_D2(mD2, dr, i, j) MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2)
				kR = kappa * R; kR2 = kR * kR;
				//fexp2 = exp(-kR2);
				Ewald_f(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4, f5, f6)
				EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j) 
				EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)
				EWALD_T4(mET4, f1, f2, f3, f4, f5, f6, dr, mT2, mD2, mT3, mT4, i, j, k, l)

				// energy
				if (emp->q != 0) Utotal += patom->c * emp->q / R * f0; // q-q
				for (i = 0; i < 3; i++) {
					Utotal += -patom->c * mET1.v[i] * emp->mu.v[i];
					if (emp->q != 0) Utotal += patom->ind_mu.v[i] * mET1.v[i] * emp->q;  // q-mu + mu-q
					for (j = 0; j < 3; j++) {
						Utotal -= patom->ind_mu.v[i] * mET2.m[i][j] * emp->mu.v[j]; // mu-mu
						Utotal += patom->c * mET2.m[i][j] * emp->Q.m[i][j] / 3; // q-quadrupole
						for (k = 0; k < 3; k++) Utotal += patom->ind_mu.v[k] * mET3.m[i][j][k] * emp->Q.m[i][j] / 3; // mu-quadrupole
					}
				}

				// electric field
				for (i = 0; i < 3; i++) {
					if (emp->q != 0) E_real.v[i] -= emp->q * mET1.v[i];
					for (j = 0; j < 3; j++) {
						E_real.v[i] += emp->mu.v[j] * mET2.m[i][j];
						for (k = 0; k < 3; k++) E_real.v[i] -= emp->Q.m[j][k] * mET3.m[i][j][k] / 3;
					}
				}
				// electric field gradient
				for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
					if (emp->q != 0) Eg_real.m[i][j] -= emp->q * mET2.m[i][j];
					for (k = 0; k < 3; k++) {
						Eg_real.m[i][j] += emp->mu.v[k] * mET3.m[i][j][k];
						for (l = 0; l < 3; l++) Eg_real.m[i][j] -= emp->Q.m[k][l] * mET4.m[i][j][k][l] / 3;
					}
				}}
				ch_emp = ch_emp->next;
			}

			ch_cluster = pc->ERship.ch_cluster;
			while (ch_cluster != NULL) {
				pc1 = ch_cluster->p;
				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					if (intrinsic_cell && pc1 == pc && na1 == na) continue; // same atom
					patom1 = pc1->atom + na1;
					for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
					if (FABS(dr.v[0]) > rcut_Ewd || FABS(dr.v[1]) > rcut_Ewd || FABS(dr.v[2]) > rcut_Ewd) {
						bEwaldCut = true;
						if (!intrinsic_cell) continue;
					}
					V3ABS2(dr, R2)
					bEwaldCut = (R2 > r2cut_Ewd ? true : false);
					if (!intrinsic_cell && bEwaldCut) continue;
					R = sqrt(R2); R3 = R * R2; R4 = R2 * R2; R5 = R3 * R2; R7 = R5 * R2; kR = kappa * R; kR2 = kR * kR;
					//EXP2(fexp2, kR, i) Ewald_f0(f0, kR, i) Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2)
					Ewald_f04(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4)
					MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) //MTP_D2(mD2, dr, i, j)
					MTP_T3(mT3, dr, R2, R7, i, j, k, t1)

					EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j) 
					EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)

					if (R2 < d2_14) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}

					if (bound || (intrinsic_cell && pc == pc1)) { // bounding atoms or they are in the same cluster
						// bound interaction has to be kicked off in Ewald Sum
						// electric field
						for (i = 0; i < 3; i++) {
							E_real.v[i] += patom1->c * (mT1.v[i] - mET1.v[i]); 
							for (j = 0; j < 3; j++) E_real.v[i] += patom1->ind_mu.v[j] * (mET2.m[i][j] - mT2.m[i][j]);
						}
						// electric field gradient
						for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
							Eg_real.m[i][j] += patom1->c * (mT2.m[i][j] - mET2.m[i][j]);
							for (k = 0; k < 3; k++) Eg_real.m[i][j] += patom1->ind_mu.v[k] * (mET3.m[i][j][k] - mT3.m[i][j][k]);
						}}

						// bound interaction has to be kicked off in Ewald Sum
						// energy, no self interaction
						Utotal -= patom->c * patom1->c / R * (1 - f0); // q-q
						for (i = 0; i < 3; i++){
							Utotal += patom->c * patom1->ind_mu.v[i] * (mT1.v[i] - mET1.v[i]) +
								patom->ind_mu.v[i] * (mET1.v[i] - mT1.v[i]); // q-mu & mu-q
							for (j = 0; j < 3; j++) Utotal -= patom->ind_mu.v[i] * (mET2.m[i][j] - mT2.m[i][j]) * patom1->ind_mu.v[j]; // mu-mu
						}
					}
					else if (bound14 && pc != pc1) { // 1-4 interaction
						// electric field
						for (i = 0; i < 3; i++) {
							t1 = patom1->c * dr.v[i] / R3;
							E_real.v[i] += t1 * f1;
							E_14.v[i] += (A14_E - 1) * t1; // 1-4 interaction will to be corrected later
							for (j = 0; j < 3; j++) {
								E_real.v[i] += patom1->ind_mu.v[j] * mET2.m[i][j];
								E_14.v[i] += patom1->ind_mu.v[j] * (A14_E - 1) * mT2.m[i][j];
							}
						}
						// electric field gradient
						for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
							Eg_real.m[i][j] -= patom1->c * mET2.m[i][j];
							Eg_14.m[i][j] -= patom1->c * (A14_E - 1) * mT2.m[i][j];
							for (k = 0; k < 3; k++) {
								Eg_real.m[i][j] += patom1->ind_mu.v[k] * mET3.m[i][j][k];
								Eg_14.m[i][j] += patom1->ind_mu.v[k] * (A14_E - 1) * mT3.m[i][j][k];
							}
						}}
						// energy
						Utotal += patom->c * patom1->par->c / R * (f0 - 1 + A14_E); // 1-4 interaction is corrected by * A14_E
						for (i = 0; i < 3; i++) {
							Utotal += (-patom->c * patom1->ind_mu.v[i] + patom1->par->c * patom->ind_mu.v[i]) * 
								(mET1.v[i] + (A14_E - 1) * mT1.v[i]);
							for (j = 0; j < 3; j++) Utotal -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * (mET2.m[i][j] + (A14_E - 1) * mT2.m[i][j]);
						}
						if (intrinsic_cell) {
							Uself += patom->c * patom1->c / R * A14_E;
							for (i = 0; i < 3; i++) {
								Uself += (-patom->c * patom1->ind_mu.v[i] + patom1->c * patom->ind_mu.v[i]) * 
									A14_E * mT1.v[i];
								for (j = 0; j < 3; j++) Uself -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * A14_E * mT2.m[i][j];
							}
						}
					}
					else {
						// electric field
						for (i = 0; i < 3; i++) {
							E_real.v[i] += f1 * patom1->c * dr.v[i] / R3;
							for (j = 0; j < 3; j++) E_real.v[i] -= patom1->ind_mu.v[j] * mET2.m[i][j];
						}
						// electric field gradient
						for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
							Eg_real.m[i][j] -= patom1->c * mET2.m[i][j];
							for (k = 0; k < 3; k++) Eg_real.m[i][j] += patom1->ind_mu.v[k] * mET3.m[i][j][k];
						}}
						// energy
						Utotal += patom->c * patom1->c / R * f0;
						for (i = 0; i < 3; i++) {
							Utotal += (patom1->c * patom->ind_mu.v[i] - patom1->ind_mu.v[i] * patom->c) * mET1.v[i];
							for (j = 0; j < 3; j++) Utotal -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mET2.m[i][j];
						}
						if (intrinsic_cell) {
							Uself += patom->c * patom1->c / R;
							for (i = 0; i < 3; i++) {
								Uself += (patom1->c * patom->ind_mu.v[i] - patom1->ind_mu.v[i] * patom->c) * mT1.v[i];
								for (j = 0; j < 3; j++) Uself -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mT2.m[i][j];
							}
						}
					}
				}
				ch_cluster = ch_cluster->next;
			}
		}}} // end of periodic cell iteration

		for (npc[0] = 0; npc[0] < NH; npc[0]++) {for (npc[1] = 0; npc[1] < NH; npc[1]++) {for (npc[2] = 0; npc[2] < NH; npc[2]++) {
			Ah = ::Ah[npc[0]][npc[1]][npc[2]]; if (Ah < AhMin) continue;
			if (npc[0] == _NH && npc[1] == _NH && npc[2] == _NH) continue;
			ch_emp = pc->ERship.ch_mpl;
			h = &(::h[npc[0]][npc[1]][npc[2]]);
			scale_uv3((*h), r0, hR)
			PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t2)
			scale_uv3((*h), patom->ind_mu, muh0)
			while (ch_emp != NULL) {
				emp = ch_emp->p;
				fcos1 = emp->fcos.m[npc[0]][npc[1]][npc[2]]; fsin1 = emp->fsin.m[npc[0]][npc[1]][npc[2]];
				fss = fsin0 * fsin1; fcc = fcos0 * fcos1; fsc = fsin0 * fcos1; fcs = fcos0 * fsin1;
				// electric field
				for (i = 0; i < 3; i++) {
					t1 = 0;
					if (emp->q != 0) t1 += emp->q * (fsc - fcs);
					t1 -= emp->muh.m[npc[0]][npc[1]][npc[2]] * (fcc + fss);
					t1 += emp->Qh.m[npc[0]][npc[1]][npc[2]] * (fcs - fsc) / 3;
					E_recip.v[i] += Ah * h->v[i] * t1;
				}
				// electric field gradient
				for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
					t1 = (emp->q == 0 ? 0 : emp->q * (fcc + fss));
					t1 += emp->muh.m[npc[0]][npc[1]][npc[2]] * (fsc - fcs);
					t1 -= emp->Qh.m[npc[0]][npc[1]][npc[2]] * (fss + fcc) / 3;
					Eg_recip.m[i][j] += Ah * h->v[i] * h->v[j] * t1;
				}}
				// energy
				Utotal += Ah * patom->c * (emp->q * (fcc + fss) + emp->muh.m[npc[0]][npc[1]][npc[2]] * (fsc - fcs)
					- emp->Qh.m[npc[0]][npc[1]][npc[2]] * (fcc + fss) / 3);
				Utotal += Ah * muh0 * (emp->q * (fcs - fsc) + emp->muh.m[npc[0]][npc[1]][npc[2]] * (fcc + fss)
					+ emp->Qh.m[npc[0]][npc[1]][npc[2]] * (fsc - fcs) / 3);

				ch_emp = ch_emp->next;
			}

			ch_cluster = pc->ERship.ch_cluster;
			while (ch_cluster != NULL) {
				pc1 = ch_cluster->p;
				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					patom1 = pc1->atom + na1;
					if (patom1->c == 0) continue;
					r1.v[0] = patom1->r.v[0] + pc1->dr_cmm.v[0]; r1.v[1] = patom1->r.v[1] + pc1->dr_cmm.v[1]; r1.v[2] = patom1->r.v[2] + pc1->dr_cmm.v[2];
					scale_uv3((*h), r1, hR)
					PERIOD(hR, 0, PI2, t1) COS(fcos1, t1, i, t1) SIN(fsin1, t1, i, t2)
					fss = fsin0 * fsin1; fcc = fcos0 * fcos1; fsc = fsin0 * fcos1; fcs = fcos0 * fsin1;
					scale_uv3(patom1->ind_mu, (*h), muh1)
					// electric field
					for (i = 0; i < 3; i++) {
						t1 = patom1->c * (fsc - fcs);
						t1 -= muh1 * (fcc + fss);
						E_recip.v[i] += Ah * h->v[i] * t1;
					}
					// electric field gradient
					for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
						t1 = patom1->c * (fcc + fss);
						t1 += muh1 * (fsc - fcs);
						Eg_recip.m[i][j] += Ah * h->v[i] * h->v[j] * t1;
					}}
					// energy
					Utotal += Ah * patom->c * (patom1->par->c * (fcc + fss) + muh1 * (fsc - fcs));
					Utotal += Ah * muh0 * (patom1->c * (fcs - fsc) + muh1 * (fcc + fss));
				}
				ch_cluster = ch_cluster->next;
			}
			for (i = 0; i < 3; i++) {
				t1 = Ah * h->v[i];
				E_recip.v[i] += t1 * muh0;
				for (j = 0; j < 3; j++) Eg_recip.m[i][j] -= t1 * h->v[j] * patom->c;
			}
		}}} // end of reciprocal sum
		for (i = 0; i < 3; i++) {
			patom->E.v[i] = E_real.v[i] + E_recip.v[i] + E_14.v[i] - inv_V * cell.mu.v[i];
			for (j = 0; j < 3; j++) patom->Eg.m[i][j] = Eg_real.m[i][j] + Eg_recip.m[i][j] + Eg_14.m[i][j] - inv_V * patom->c;
		}

		// LJ interaction
		if (patom->par->epsLJ > 0.0001) {
			ch_imag_cluster = pc->LJRship.ch;
			while (ch_imag_cluster != NULL) {
				pc1 = ch_imag_cluster->p->p;
				if (pc1 == pc) {ch_imag_cluster = ch_imag_cluster->next; continue;} // ingnore the interaction inside cluster
				drCell.v[0] = ch_imag_cluster->p->nx * cell.xd[0];
				drCell.v[1] = ch_imag_cluster->p->ny * cell.xd[1];
				drCell.v[2] = ch_imag_cluster->p->nz * cell.xd[2];

				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					patom1 = pc1->atom + na1;
					if (patom1->par->epsLJ < 0.0001) continue;
					for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
					//ABS2(dr, i, R2)
					R2 = dr.v[0] * dr.v[0]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[1] * dr.v[1]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[2] * dr.v[2]; if (R2 > r2cut_LJ) continue;
					if (R2 > r2cut_LJ) continue;

					pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					if (pLJ != NULL && pLJ->epsLJ > 0.0001) {
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
						if (ch_imag_cluster->p->local) t2 *= 0.5; // will be calculated twice
						U_LJ += t2;
						for (i = 0; i < 3; i++) force.v[i] *= unit_LJ_estat_force; // change to electostatic interaction unit
						V3PV3(patom->r0, force, torque)
						for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
					}
				}
				ch_imag_cluster = ch_imag_cluster->next;
			}
		}

		// energy
		scale_uv3(r0, cell.mu, t1)
		t2 = patom->c * t1 * inv_V - 2 * kappa * INV_SQRT_PI * patom->c * patom->c;
		UTotal += Utotal + t2;
		USelf += Uself;
	}
	// total energy in kT
	UTotal = (UTotal * U_ElectFr + U_LJ * U_eps) / kT;
	USelf = USelf * U_ElectFr / kT;
}

void Polar_EFIELD_EwaldSum_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, double &Utotal, bool &OVERLAP) {
	int nc;
	CHAIN<BASIC_CLUSTER> *ch = NULL;
	BASIC_CLUSTER *pc = NULL;
	double Uself = 0, U1 = 0, U2 = 0;
	bool overlap = false; OVERLAP = false;
	Utotal = 0;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		ch = cell.cell[nc].ch;
		while (ch != NULL) {
			pc = ch->p;
			Polar_EFIELD_EwaldSum_CMM(cell, pc, U1, U2, overlap);
			if (overlap) OVERLAP = true;
			Utotal += U1; Uself += U2;
			ch = ch->next;
		}
	}
	Utotal -= Uself * 0.5;
}

void Polar_CMMEwaldSumInteraction(CMM_CELL3D<BASIC_CLUSTER> &cell, double &V, bool &OVERLAP, bool hard_sphere, int n1Cell, int n2Cell) { // V in kT
	int nc, na, i, j;
	CHAIN<BASIC_CLUSTER> *ch = NULL;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	double unit_estat_mass = U_ElectForce / U_MassAccel;
	double unit_Eext_estat_force = U_EextForce / U_ElectForce;
	float coeff = sqrt(::eps_dc);

	// calculate the electric field on each atom in these cells
	Polar_EFIELD_EwaldSum_CMM(cell, n1Cell, n2Cell, V, OVERLAP);

	for (nc = n1Cell; nc <= n2Cell; nc++) {
		ch = cell.cell[nc].ch;
		while (ch != NULL) {
			pc = ch->p;
			for (na = 0; na < pc->nAtoms; na++) {
				patom = pc->atom + na;
				if (patom->c != 0) {
					for (i = 0; i < 3; i++) {
						F.v[i] = patom->c * patom->E.v[i];
						for (j = 0; j < 3; j++) {
							F.v[i] += patom->ind_mu.v[j] * patom->Eg.m[i][j];
						}
						if (bEext) {if (i == 2) F.v[i] += patom->c0 * Eext * unit_Eext_estat_force;} // external electric field
						pc->dyn->fc.v[3+i] += F.v[i];
					}
					V3PV3(patom->r0, F, torque)
					for (i = 0; i < 3; i++) pc->dyn->fc.v[i] += torque.v[i];
					if (bEext) {
						V -= patom->c0 * patom->r0.v[2] * Eext * U_EextFr / kT;
					}
				}
			}
			for (i = 0; i < 6; i++) pc->dyn->fc.v[i] *= unit_estat_mass; // change unit
			ch = ch->next;
		}
		//if (bEext) {
		//	V -= coeff * (cell.cell[nc].emp->mu.v[2] + cell.cell[nc].emp->q * cell.cell[nc].emp->r0.v[2]) * Eext * U_EextFr / kT; // external electric field
		//}
	}
}

#endif

#endif // _DISABLE_
