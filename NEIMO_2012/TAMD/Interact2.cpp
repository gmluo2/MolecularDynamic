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

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
long t_real = 0, t_recp = 0, t_emp = 0, t_cluster = 0, t_LJ = 0;
#endif

#if _DISABLE_ == 0

void EFIELD_EwaldSum_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, double& UTotal, double& USelf, bool &OVERLAP) {
	MATOM *patom = NULL;
	VECTOR3 E_real, E_recip, E_14;
	VECTOR3 dr, u;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

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

	VECTOR3 mT1, mET1;
	DMATRIX<3> mT2, mET2, mD2;
	CMATRIX<3> mT3, mET3;
	//QMATRIX<3> mT4, mET4;

	double R, R2, R3, R4, R5, R7;//, R9;
	double kR, kR2;
	double fexp2, f0, f1, f2, f3, f4;//, f5, f6;
	double t1, t2, t3, Ah, fcos0, fsin0, fcos1, fsin1, hR, fss, fcc, fcs, fsc;
	int npc[3];
	VECTOR3 drCell;
	VECTOR3 *h = NULL;

	bool bound = false, bound14 = false, electric = true, intrinsic_cell = true, bEwaldCut = false;

	double U_LJ = 0;
	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;

	double Uestat = 0, Uself = 0;
	UTotal = 0; USelf = 0;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif


	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		Uestat = 0; Uself = 0;
		for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];

	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 1
		// ****** ignore the electric interaction between atoms *******
		goto _LJ_INTERACT;
		// ************************************************************
	#else
		if (FABS(patom->c) < 0.0001) goto _LJ_INTERACT;  // disable efield calculation
	#endif

	#if _EWALD_SUM == _MultiPole_EWALD_SUM

		if (patom->c != 0) {
			V3zero(E_real) V3zero(E_recip) V3zero(E_14)

		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif

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
					R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; 
					MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j)

					// self energy
					if (intrinsic_cell) {
						if (emp->q != 0) Uself += emp->q / R;
						for (i = 0; i < 3; i++) {
							Uself -= mT1.v[i] * emp->mu.v[i]; 
							for (j = 0; j < 3; j++) Uself += mT2.m[i][j] * emp->Q.m[i][j] / 3;
						}
						if (bEwaldCut) {ch_emp = ch_emp->next; continue;}
					}

					R7 = R5 * R2;// R9 = R7 * R2;
					MTP_T3(mT3, dr, R2, R7, i, j, k, t1)  MTP_D2(mD2, dr, i, j)
					kR = kappa * R; kR2 = kR * kR;
					//fexp2 = exp(-kR2);
					Ewald_f04(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4)
					EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j) 
					EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)

					// energy
					if (emp->q != 0) Uestat += emp->q / R * f0;
					for (i = 0; i < 3; i++) {
						Uestat -= mET1.v[i] * emp->mu.v[i];
						for (j = 0; j < 3; j++) Uestat += mET2.m[i][j] * emp->Q.m[i][j] / 3;
					}

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
						electric = (patom1->c == 0 ? false : true) ; // this atom has charge ?
						if (!electric) continue; // ignore this atom for charge interaction
						for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
						if (FABS(dr.v[0]) > rcut_Ewd || FABS(dr.v[1]) > rcut_Ewd || FABS(dr.v[2]) > rcut_Ewd) {
							bEwaldCut = true;
							if (!intrinsic_cell) continue;
						}
						V3ABS2(dr, R2)
						bEwaldCut = (R2 > r2cut_Ewd ? true : false);
						if (!intrinsic_cell && bEwaldCut) continue;
						if (electric) {
							R = sqrt(R2); R3 = R * R2; kR = kappa * R; 
							EXP2(fexp2, kR, i) Ewald_f0(f0, kR, i) Ewald_f1(f1, kR, f0, fexp2)
						}

						if (R2 < d2_14) {
							BOUND(patom, patom1, i, bound) // bound atom
							if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
							if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
							else bound14 = false;
						}
						else {bound = false; bound14 = false;}

						if (bound || (intrinsic_cell && pc == pc1)) { // bounding atoms or they are in the same cluster
							if (electric) {// electric field
								for (i = 0; i < 3; i++) {
									t1 = patom1->c * dr.v[i] / R3;
									E_real.v[i] -= t1 * (1 - f1); // bound interaction has to be kicked off in Ewald Sum
								}
								// energy, no self interaction
								Uestat -= patom1->c / R * (1 - f0); // bound interaction has to be kicked off in Ewald Sum
							}
						}
						else if (bound14 && pc != pc1) { // 1-4 interaction
							if (electric) {// electric field
								for (i = 0; i < 3; i++) {
									t1 = patom1->c * dr.v[i] / R3;
									E_real.v[i] += t1 * f1;
									E_14.v[i] += (A14_E - 1) * t1; // 1-4 interaction will to be corrected later
								}
								// energy
								Uestat += patom1->c / R * (f0 - 1 + A14_E); // 1-4 interaction is corrected by * A14_E
								if (intrinsic_cell) Uself += patom1->c / R * A14_E;
							}
						}
						else {
							if (electric) {// electric field
								for (i = 0; i < 3; i++) E_real.v[i] += f1 * patom1->c * dr.v[i] / R3;
								// energy
								Uestat += patom1->c / R * f0;
								if (intrinsic_cell) Uself += patom1->c / R;
							}
						}
					}
					ch_cluster = ch_cluster->next;
				}
			}}} // end of periodic cell iteration
			#if _CHECK_TIME_STAMP_ == 1
			t_real += time.elapse();
			//time.start();
			#endif

			for (npc[0] = 0; npc[0] < NH; npc[0]++) {for (npc[1] = 0; npc[1] < NH; npc[1]++) {for (npc[2] = 0; npc[2] < NH; npc[2]++) {
				Ah = ::Ah[npc[0]][npc[1]][npc[2]]; if (Ah < AhMin) continue;
				if (npc[0] == _NH && npc[1] == _NH && npc[2] == _NH) continue;
				h = &(::h[npc[0]][npc[1]][npc[2]]); 
				ch_emp = pc->ERship.ch_mpl;
				scale_uv3((*h), r0, hR)
				PERIOD(hR, 0, PI2, t1) fcos0 = cos(t1); fsin0 = sin(t1); //SCOS(fcos0, t1, i, t1) SSIN(fsin0, t1, i, t2)
				/*
				{
					int nemp = number(ch_emp);
					nemp = nemp;
				}
				*/
				#if _CHECK_TIME_STAMP_ == 1
				time.start();
				#endif
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
					// energy
					t1 = (emp->q == 0 ? 0 : emp->q * (fcc + fss));
					t1 += emp->muh.m[npc[0]][npc[1]][npc[2]] * (fsc - fcs)
						- emp->Qh.m[npc[0]][npc[1]][npc[2]] * (fcc + fss) / 3;
					Uestat += Ah * t1;

					ch_emp = ch_emp->next;
				}
				#if _CHECK_TIME_STAMP_ == 1
				t_emp += time.elapse();
				t_recp += time.elapse();
				#endif

				#if _CHECK_TIME_STAMP_ == 1
				time.start();
				#endif
				ch_cluster = pc->ERship.ch_cluster;
				/*
				{
					int ncluster = number(ch_cluster);
					ncluster = ncluster;
				}
				*/
				while (ch_cluster != NULL) {
					pc1 = ch_cluster->p;
					for (na1 = 0; na1 < pc1->nAtoms; na1++) {
						patom1 = pc1->atom + na1;
						if (patom1->c == 0) continue;
						r1.v[0] = patom1->r.v[0] + pc1->dr_cmm.v[0]; r1.v[1] = patom1->r.v[1] + pc1->dr_cmm.v[1]; r1.v[2] = patom1->r.v[2] + pc1->dr_cmm.v[2];
						scale_uv3((*h), r1, hR)
						PERIOD(hR, 0, PI2, t1) fcos1 = cos(t1); fsin1 = sin(t1); //SCOS(fcos1, t1, i, t1) SSIN(fsin1, t1, i, t2)
						fcc = fcos0 * fcos1; fss = fsin0 * fsin1; fcs = fcos0 * fsin1; fsc = fsin0 * fcos1;
						for (i = 0; i < 3; i++) E_recip.v[i] += Ah * h->v[i] * patom1->c * (fsc - fcs);

						// energy
						Uestat += Ah * patom1->c * (fcc + fss);
					}
					ch_cluster = ch_cluster->next;
				}
				#if _CHECK_TIME_STAMP_ == 1
				t_cluster += time.elapse();
				t_recp += time.elapse();
				#endif
			}}} // end of reciprocal sum
			#if _CHECK_TIME_STAMP_ == 1
			//t_recp += time.elapse();
			#endif
			for (i = 0; i < 3; i++) patom->E.v[i] = E_real.v[i] + E_recip.v[i] + E_14.v[i] - inv_V * cell.mu.v[i];
		}
	#endif // _EWALD_SUM == _MultiPole_EWALD_SUM

_LJ_INTERACT:

		#if _CHECK_TIME_STAMP_ == 1
		time.start();
		#endif
		// LJ interaction
		if (patom->par->epsLJ > 0.0001) {
			//ch_imag_cluster = pc->LJRship.ch;
			LJit.set(&(pc->LJRship)); LJit.start_iterate();
			//while (ch_imag_cluster != NULL) {
			while (!LJit.eoi()) {
				//pc1 = ch_imag_cluster->p->p;
				imag_cluster = LJit.current(); if (imag_cluster == NULL) {LJit.next_iterate(); continue;}
				pc1 = imag_cluster->p;
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
					if (patom1->par->epsLJ < 0.0001) continue;
					for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
					//ABS2(dr, i, R2)
					R2 = dr.v[0] * dr.v[0]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[1] * dr.v[1]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[2] * dr.v[2]; if (R2 > r2cut_LJ) continue;
					if (R2 > r2cut_LJ) continue;
					
					if (R2 < _R2min_Interact) continue; // ignore the interaction when they are too close

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
						if (bound14) { //1-4 LJ interaction / 2
							t2 *= A14_LJ; 
							for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
						}
						if (imag_cluster->local) t2 *= 0.5;  // will be calculated twice
						U_LJ += t2;
						V3PV3(patom->r0, force, torque)
						for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
					}
				}
				//ch_imag_cluster = ch_imag_cluster->next;
				LJit.next_iterate();
			}
		}
		#if _CHECK_TIME_STAMP_ == 1
		t_LJ += time.elapse();
		#endif

		// energy
	#if _EWALD_SUM == _MultiPole_EWALD_SUM
		if (!pc->eNeutral && patom->c != 0) {
			Uself = patom->c * Uself;
			Uestat -= 2 * kappa * INV_SQRT_PI * patom->c;
			Uestat = patom->c * Uestat;

			scale_uv3(r0, cell.mu, t1)
			t2 = patom->c * t1 * inv_V;
			Uestat += t2;

			UTotal += Uestat; USelf += Uself;
		}
	#endif //#if _EWALD_SUM == _MultiPole_EWALD_SUM
	}
	// total energy in kT
	UTotal = UTotal * eUnit_estat_kT + U_LJ;
	USelf = USelf * eUnit_estat_kT;
}

void EFIELD_EwaldSum_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, double &Utotal, bool &OVERLAP) {
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	double Uself = 0, U1 = 0, U2 = 0;
	bool overlap = false; OVERLAP = false;
	Utotal = 0;
	
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			EFIELD_EwaldSum_CMM(cell, pc, U1, U2, overlap);
			if (overlap) OVERLAP = true;
			Utotal += U1; Uself += U2;
			it.next_iterate();
		}
	}

	Utotal -= Uself * 0.5;
	return;
}

void CMMEwaldSumInteraction(CMM_CELL3D<BASIC_CLUSTER> &cell, double &V, bool &OVERLAP, bool hard_sphere, int n1Cell, int n2Cell) { // V in kT
	int nc, na, i;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;

	float coeff = sqrt(::eps_dc);

	#if _CHECK_TIME_STAMP_ == 1
	t_real = 0, t_recp = 0, t_emp = 0, t_cluster = 0, t_LJ = 0;
	char msg[256] = "\0";
	#endif
	// calculate the electric field on each atom in these cells
	EFIELD_EwaldSum_CMM(cell, n1Cell, n2Cell, V, OVERLAP);
	#if _CHECK_TIME_STAMP_ == 1
	sprintf(msg, "real space: %d, recp. space: %d, emp: %d, cluster: %d, LJ interact: %d", t_real, t_recp, t_emp, t_cluster, t_LJ);
	show_infor(msg);
	#endif


	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			for (na = 0; na < pc->nAtoms; na++) {
				patom = pc->atom + na;
				if (patom->c != 0) {
					for (i = 0; i < 3; i++) {
						F.v[i] = patom->c * patom->E.v[i] * fUnit_estat_mass;
						if (bEext) {if (i == 2) F.v[i] += patom->c0 * Eext * fUnit_estat_mass;} // external electric field
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

void EFIELD_EwaldSum_CMM_SINGLE_MM(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, double &Utotal, bool &OVERLAP) {
	int nc;
	BASIC_CLUSTER *pc = NULL;
	double Uself = 0, U1 = 0, U2 = 0;
	bool overlap = false; OVERLAP = false;
	Utotal = 0;
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = (BASIC_CLUSTER*)(mm.cluster + nc);
		EFIELD_EwaldSum_CMM(cell, pc, U1, U2, overlap);
		if (overlap) OVERLAP = true;
		Utotal += U1; Uself += U2;
	}
	Utotal -= Uself * 0.5;
}

void EwaldSumCMMInteraction_SingleMM(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, double &V, bool &OVERLAP, bool hard_sphere, int n1, int n2) { // V in kT
	int nc, na, i;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	// calculate the electric field on each atom in these cells
	EFIELD_EwaldSum_CMM_SINGLE_MM(cell, mm, n1, n2, V, OVERLAP);

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

void EwaldSumCMMInteraction_SM(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE *sm, int nSM, double &V, bool &OVERLAP) { // V in kT
	double Utotal = 0, Uself = 0;
	EFIELD_EwaldSum_CMM(cell, sm->c, Utotal, Uself, OVERLAP);
	V = Utotal - Uself * 0.5;
}

#endif // _DISABLE_

#define TEST 0
#if TEST == 1
void efield_test() {
	int N = 6;
	VECTOR3 r0[6]; double q[6] = {1, 2, 1, 2, 1, 2};
	r0[0] = VECTOR3(1, 0, 1);
	r0[1] = VECTOR3(1, 1, 1);
	r0[2] = VECTOR3(0, -1, 0);
	r0[3] = VECTOR3(0, 1, -1);
	r0[4] = VECTOR3(1, 1, 0);
	r0[5] = VECTOR3(-1, 0, 0);

	VECTOR3 r1[2];

	double kappa = 1;
	double q1 = 0, q2 = 1;
	VECTOR3 r;
	VECTOR3 mu1, mu2;
	DMATRIX<3> Q1, Q2;

	VECTOR3 mT1, mET1;
	DMATRIX<3> mT2, mET2, mD2;
	CMATRIX<3> mT3, mET3;
	QMATRIX<3> mT4, mET4;

	double R, R2, R3, R4, R5, R7, R9;
	//double f0, f1, f2, f3, f4, f5, f6;
	double t1, t2;

	VECTOR3 E, F;

	VECTOR3 E_dir, F_dir, r_dir;

	int i, j, k, l;
	double dis, dis2, dis3;

	q1 = 0; V3zero(mu1); Mzero(Q1, i, j);
	for (i = 0; i < N; i++) {
		q1 += q[i];
		for (j = 0; j < 3; j++) mu1.v[j] += q[i] * r0[i].v[j];
		V3ABS2(r0[i], R2)
		for (k = 0; k < 3; k++) for (l = 0; l < 3; l++) {
			Q1.m[k][l] += 1.5 * q[i] * r0[i].v[k] * r0[i].v[l];
		}
	}

	//q2 = q[0] + q[1];
	//for (i = 0; i < 3; i++) mu2.v[i] = q[0] * r1[0].v[i] + q[1] * r1[1].v[i];
	//for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) Q2.m[i][j] = 1.5 * (q[0] * r1[0].v[i] * r1[0].v[j] + q[1] * r1[1].v[i] * r1[1].v[j]);

	q2 = q[0];

	ofstream out;
	out.open("test-F.dat");
	for (R = 3; R < 20; R += 1) {
		r.v[0] = 0; r.v[1] = 0; r.v[2] = R;
		R2 = R * R; R3 = R * R2; R4 = R2 * R2; R5 = R * R4; R7 = R2 * R5; R9 = R2 * R7;
		MTP_T1(mT1, r, R3)
		MTP_T2(mT2, r, R2, R5, i, j)
		MTP_T3(mT3, r, R2, R7, i, j, k, t1)
		MTP_T4(mT4, r, R2, R4, R9, i, j, k, l, t1, t2)

		//EFIELD_MULTIPOL_REAL(E, q1, mu1, Q1, mT1, mT2, mT3, i, j, k)
		//FORCE1_REAL(F, q2, mu2, q1, mu1, mT1, mT2, mT3, i, j, k)
		//FORCE2_REAL(F, q2, mu2, q1, mu1, Q1, mT1, mT2, mT3, mT4, i, j, k, l)
		FORCE3_REAL(F, q2, mu2, Q2, q1, mu1, Q1, mT1, mT2, mT3, mT4, i, j, k, l)

		V3zero(E_dir, i)
		Vzero(F_dir, i)
		for (i = 0; i < N; i++) { 
			//j = 0;
			for (j = 0; j < 2; j++) {
				r_dir = r + r1[j] - r0[i];
				//r_dir = r - r0[i];
				V3ABS2(r_dir, dis2)
				dis = sqrt(dis2);
				dis3 = dis * dis2;
				for (k = 0; k < 3; k++) {
					//E_dir.v[k] += q[i] / dis3 * r_dir.v[k];
					F_dir.v[k] += q[j] * q[i] / dis3 * r_dir.v[k];
				}
			}
		}

		//out<<R<<"  "<<E.v[0]<<"  "<<E.v[1]<<"  "<<E.v[2]<<"  "<<E_dir.v[0]<<"  "<<E_dir.v[1]<<"  "<<E_dir.v[2]<<endl;
		out<<R<<"  "<<F.v[0]<<"  "<<F.v[1]<<"  "<<F.v[2]<<"  "<<F_dir.v[0]<<"  "<<F_dir.v[1]<<"  "<<F_dir.v[2]<<"  ";
		t1 = F.abs(); t2 = F_dir.abs();
		out<<t1<<"  "<<t2<<endl;
	}
	out.close();
}

#define NATOMS 4
class EMOL {
public:
	float c[NATOMS];
	VECTOR3 u[NATOMS]; // (induced) dipole moment of of the atom
	VECTOR3 r[NATOMS];
	VECTOR3 r0;
	float q;
	VECTOR3 mu;
	DMATRIX<3> Q;
	void rotate(VECTOR3 &axis, double omega) {
		int i = 0, k, l;
		R rmatrix; VECTOR3 u;
		RMatrix(axis, omega, rmatrix, true);
		for (i = 0; i < NATOMS; i++) {
			u = r[i];
			MpV(rmatrix, u, r[i], k, l, 3)
		}
	};
	EMOL() {
		int i = 0, j = 0;
		for (i = 0; i < NATOMS; i++) {
			c[i] = ((i / 2 * 2 == i) ? 1 : -1) * 10;
			rand_vect3(u[i].v); u[i] *= 10;
			for (j = 0; j < 3; j++) r[i].v[j] = 2 * (ranf() - 0.5); // ranf() / 2;
		}
	};
	void shift(VECTOR3 dr) {
		int i = 0;
		for (i = 0; i < 3; i++) r0.v[i] += dr.v[i];
	};
	void calMultipole() {
		int i, j, k;
		q = 0;
		Vzero(mu, i)
		Mzero(Q, i, j)
		for (i = 0; i < NATOMS; i++) {
			q += c[i];
			for (j = 0; j < 3; j++) {
				mu.v[j] += c[i] * r[i].v[j];
				mu.v[j] += u[i].v[j]; // plus induced dipole of each atom
				for (k = 0; k < 3; k++) {
					Q.m[j][k] += 1.5 * c[i] * r[i].v[j] * r[i].v[k];
				}
			}
		}
	};
};

void test_EwaldSum_EF() {
	//int NMOL = 125, NDIM = 5;
	//EMOL mol[125];
	int NMOL = 2, NDIM = 2;
	EMOL mol[8];
	float dcell = 20;
	int i, j, k, l, m, n;
	VECTOR3 dr, mu_total;
	for (i = 0; i < NDIM; i++) for (j = 0; j < NDIM; j++) for (k = 0; k < NDIM; k++) {
		dr.v[0] = i * dcell / NDIM; dr.v[1] = j * dcell / NDIM; dr.v[2] = k * dcell / NDIM;
		n = i * NDIM * NDIM + j * NDIM + k;
		mol[n].shift(dr);
		mol[n].calMultipole();
		mu_total += mol[n].mu;
	}

	int nc[3], nCell = 10, nH = 6;
	float kappa = 3.2f / dcell;
	VECTOR3 E, E_recip, Edir, F, F_recip, Fdir, dR;

	VECTOR3 mT1, mET1;
	DMATRIX<3> mT2, mET2, mD2;
	CMATRIX<3> mT3, mET3;
	QMATRIX<3> mT4, mET4;

	double R, R2, R3, R4, R5, R7, R9;
	double kappa3 = kappa * kappa * kappa;
	double kR, kR2;
	double f0, f1, f2, f3, f4, f5, f6, fexp2, Ah;
	VECTOR3 h;
	double inv_k2 = 0.25 / (kappa * kappa), inv_V = 4 * PI / (dcell * dcell * dcell);
	double fcos0, fsin0, fcos1, fsin1, hR, hk, muh, muh0, Qh, Qh0, H2;
	double t1, t2;

	nCell = 20;
	m = 0; int matom = 0, natom;

	Vzero(E, i) Vzero(F, i) Vzero(Edir, i) Vzero(Fdir, i)
	for (nc[0] = -nCell; nc[0] <= nCell; nc[0]++) {for (nc[1] = -nCell; nc[1] <= nCell; nc[1]++) {for (nc[2] = -nCell; nc[2] <= nCell; nc[2]++) {
		for (i = 0; i < 3; i++) dR.v[i] = dcell * nc[i];
		for (n = 0; n < NMOL; n++) {
			if (n < NMOL) {
			//if (n == m && nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
			//if (n == m) {
			//if (n > NMOL) {
				for (natom = 0; natom < NATOMS; natom++) {
					if (n == m && natom == matom && nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
					for (i = 0; i < 3; i++) dr.v[i] = mol[m].r0.v[i] + mol[m].r[matom].v[i] - mol[n].r0.v[i] - mol[n].r[natom].v[i] - dR.v[i];
					//dr = mol[m].r0 + mol[m].r[matom] - mol[n].r0 - mol[n].r[natom] - dR;
					V3ABS2(dr, R2)
					R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2; R9 = R7 * R2;
					kR = kappa * R; kR2 = kR * kR;
					//fexp2 = exp(-kR2);
					Ewald_f04(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4)
					MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
					MTP_T3(mT3, dr, R2, R7, i, j, k, t1)

					EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j) 
					EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)

					for (i = 0; i < 3; i++) {
						E.v[i] -= mol[n].c[natom] * mET1.v[i]; // from charge
						Edir.v[i] -= mol[n].c[natom] * mT1.v[i];
						
						F.v[i] -= mol[m].c[matom] * mET1.v[i] * mol[n].c[natom]; // charge - charge
						Fdir.v[i] -= mol[m].c[matom] * mT1.v[i] * mol[n].c[natom];
						
						for (j = 0; j < 3; j++) {
							E.v[i] += mET2.m[i][j] * mol[n].u[natom].v[j]; // from dipole
							Edir.v[i] += mT2.m[i][j] * mol[n].u[natom].v[j];
							
							F.v[i] += mol[m].c[matom] * mET2.m[i][j] * mol[n].u[natom].v[j]
								- mol[m].u[matom].v[j] * mET2.m[i][j] * mol[n].c[natom]; // charge - dipole & dipole - charge
							Fdir.v[i] += mol[m].c[matom] * mT2.m[i][j] * mol[n].u[natom].v[j]
								- mol[m].u[matom].v[j] * mT2.m[i][j] * mol[n].c[natom];
							
							for (k = 0; k < 3; k++) {
								F.v[i] += mol[m].u[matom].v[j] * mET3.m[i][j][k] * mol[n].u[natom].v[k]; // dipole - dipole
								Fdir.v[i] += mol[m].u[matom].v[j] * mT3.m[i][j][k] * mol[n].u[natom].v[k];
							}
						}
					}
				}
			}
			else {
				dr = mol[m].r0 + mol[m].r[matom] - mol[n].r0 - dR;
				V3ABS2(dr, R2)
				R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2; R9 = R7 * R2;
				kR = kappa * R; kR2 = kR * kR;

				Ewald_f(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4, f5, f6)

				MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_T3(mT3, dr, R2, R7, i, j, k, t1)
				MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2) MTP_D2(mD2, dr, i, j)

				EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j)
				EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)
				EWALD_T4(mET4, f1, f2, f3, f4, dr, mT2, mD2, mT3, mT4, i, j, k, l)

				for (i = 0; i < 3; i++) {
					for (j = 0; j < 3; j++) {
						E.v[i] += mol[n].mu.v[j] * mET2.m[i][j];
						for (k = 0; k < 3; k++) E.v[i] -= mol[n].Q.m[j][k] * mET3.m[i][j][k] / 3;
	
						Edir.v[i] += mol[n].mu.v[j] * mT2.m[i][j];
						for (k = 0; k < 3; k++) Edir.v[i] -= mol[n].Q.m[j][k] * mT3.m[i][j][k] / 3;
					}
				}
			}
		}
	}}}

	nH = int(kappa * dcell + 0.5); if (nH < 1) nH = 1;
	Vzero(E_recip, i)
	double e_recip, f_recip;
	for (nc[0] = -nH; nc[0] <= nH; nc[0]++) {for (nc[1] = -nH; nc[1] <= nH; nc[1]++) {for (nc[2] = -nH; nc[2] <= nH; nc[2]++) {
	//for (nc[0] = 0; nc[0] <= nH; nc[0]++) {for (nc[1] = -nH; nc[1] <= nH; nc[1]++) {for (nc[2] = -nH; nc[2] <= nH; nc[2]++) {
		if (nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
		if (nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2] > 2 * nH * nH) continue;
		for (i = 0; i < 3; i++) h.v[i] = PI2 / dcell * nc[i];

		dr = mol[m].r0 + mol[m].r[matom];
		scale_uv3(h, dr, hR)
		PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
		//fcos0 = cos(hR); fsin0 = sin(hR);
		scale_uv3(mol[m].u[matom], h, muh0)

		V3ABS2(h, H2)
		hk = H2 * inv_k2;
		//EXP(Ah, hk, i)
		Ah = exp(-hk);
		if (Ah < 1e-8) continue;
		Ah *= inv_V / H2;
		e_recip = 0; f_recip = 0;
		for (n = 0; n < NMOL; n++) {
			//if (n > NMOL) {
			//if (n < NMOL) {
			if (n == m) {
				for (natom = 0; natom < NATOMS; natom++) {
					for (i = 0; i < 3; i++) dr.v[i] = mol[n].r0.v[i] + mol[n].r[natom].v[i];
					scale_uv3(h, dr, hR)
					PERIOD(hR, 0, PI2, t1) COS(fcos1, t1, i, t1) SIN(fsin1, t1, i, t1)
					//fcos1 = cos(hR); fsin1 = sin(hR);
					scale_uv3(mol[n].u[natom], h, muh)

					e_recip += mol[n].c[natom] * (fsin0 * fcos1 - fcos0 * fsin1);
					e_recip -= muh * (fcos0 * fcos1 + fsin0 * fsin1);

					f_recip += mol[m].c[matom] * (mol[n].c[natom] * (fsin0 * fcos1 - fcos0 * fsin1) 
						- muh * (fcos0 * fcos1 + fsin0 * fsin1)); // q-q & q - dipole
					f_recip += muh0 * (fcos0 * fcos1 + fsin0 * fsin1) * mol[n].c[natom]; // dipole - q
					f_recip += muh0 * (fsin0 * fcos1 - fcos0 * fsin1) * muh; // dipole - dipole
				}
			}
			else {
				dr = mol[n].r0;
				scale_uv3(h, dr, hR)
				PERIOD(hR, 0, PI2, t1) COS(fcos1, t1, i, t1) SIN(fsin1, t1, i, t1)
				//fcos1 = cos(hR); fsin1 = sin(hR);
				scale_uv3(mol[n].mu, h, muh)
				scale_VMV(h, mol[n].Q, h, Qh, i, j, 3)
				e_recip -= (fcos0 * muh * fcos1 + fsin0 * muh * fsin1);
				e_recip += 1.0 / 3 * (-fsin0 * Qh * fcos1 + fcos0 * Qh * fsin1);

				f_recip += mol[m].c[matom] * (mol[n].q * (fsin0 * fcos1 - fcos0 * fsin1) 
					- muh * (fcos0 * fcos1 + fsin0 * fsin1)); // q~q ; q~mu
				f_recip += mol[m].c[matom] * (fcos0 * fsin1 - fsin0 * fcos1) * Qh / 3; // q~Q
				f_recip += muh0 * (fcos0 * fcos1 + fsin0 * fsin1) * mol[n].q; // mu~q
				f_recip += muh0 * (fsin0 * fcos1 - fcos0 * fsin1) * muh; // mu ~ mu
				f_recip -= muh0 * (fcos0 * fcos1 + fsin0 * fsin1) * Qh / 3; // mu~Q
			}
		}
		//scale_uv3(mol[m].mu, h, t1)
		e_recip += muh0;
		for (i = 0; i < 3; i++) E_recip.v[i] += Ah * h.v[i] * e_recip;
		for (i = 0; i < 3; i++) F_recip.v[i] += Ah * h.v[i] * f_recip;
	}}}

	for (i = 0; i < 3; i++) {
		t1 = 0;
		for (n = 0; n < NMOL; n++) {
			for (natom = 0; natom < NATOMS; natom++) {
				t1 += mol[n].u[natom].v[i];
				t1 += mol[n].c[natom] * (mol[n].r[natom].v[i] + mol[n].r0.v[i]);
			}
		}
		E_recip.v[i] -= inv_V / 3 * t1;
		F_recip.v[i] -= inv_V / 3 * (mol[m].c[matom] * t1 + mol[m].u[matom].v[i] * mol[m].c[matom]);
	}

	for (i = 0; i < 3; i++) {E.v[i] += E_recip.v[i]; F.v[i] += F_recip.v[i];}

	VECTOR3 Eself, dE;

	dE = Edir - E;

	/*
	for (natom = 0; natom < NATOMS; natom++) {
		if (natom == matom) continue;
		dr = mol[m].r[matom] - mol[m].r[natom];
		V3ABS2(dr, R2)
		R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2; R9 = R7 * R2;
		kR = kappa * R; kR2 = kR * kR;
		fexp2 = exp(-kR2); Ewald_f0(f0, kR, i) Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2)
		MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
		EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j)

		for (i = 0; i < 3; i++) {
			Eself.v[i] += mol[m].c[natom] * (mT1.v[i] - mET1.v[i]);
		}
	}
	*/


	VECTOR3 E_Ewald = E + Eself;

	char msg[512] = "\0", msg1[256] = "\0", msg2[256] = "\0";
	E.Show(msg1); Edir.Show(msg2);
	sprintf(msg, "Ewald E : \n %s \n\n Direct E : \n %s\n", msg1, msg2);
	ofstream out;
	out.open("test-Ewald.dat");
	out<<msg<<endl;
	out.close();
}

void test_EwaldSum() {
	int NMOL = 2, NDIM = 2;
	float q[8];
	VECTOR3 r[8];
	float dcell = 20;
	int i, j, k, l, m, n;
	VECTOR3 dr;
	float q0 = 0;
	for (i = 0; i < NDIM; i++) for (j = 0; j < NDIM; j++) for (k = 0; k < NDIM; k++) {
		n = i * NDIM * NDIM + j * NDIM + k;
		r[n].v[0] = dcell * (ranf() - 0.5);
		r[n].v[1] = dcell * (ranf() - 0.5);
		r[n].v[2] = dcell * (ranf() - 0.5);
		//r[n].v[0] = dcell * i / NDIM;
		//r[n].v[1] = dcell * j / NDIM;
		//r[n].v[2] = dcell * k / NDIM;
		q[n] = (n == 0 ? 1 : -q[n-1]);
		q0 += q[n];
	}

	int nc[3], nCell = 4, nH = 6;
	float kappa = 3 / dcell;
	VECTOR3 E, F, Edir, E_recip, Fdir, dR;
	double U, Udir, U_recip;

	VECTOR3 mT1, mET1;
	DMATRIX<3> mT2, mET2, mD2;
	CMATRIX<3> mT3, mET3;
	QMATRIX<3> mT4, mET4;

	double R, R2, R3, R4, R5, R7, R9;
	double kappa3 = kappa * kappa * kappa;
	double kR, kR2;
	double f0, f1, f2, f3, f4, f5, f6, fexp2, Ah;
	VECTOR3 h;
	double inv_k2 = 0.25 / (kappa * kappa), inv_V = 4 * PI / (dcell * dcell * dcell);
	double fcos0, fsin0, fcos1, fsin1, hR, hk, muh, Qh, H2;
	double t1, t2;

	Vzero(E, i) Vzero(F, i) Vzero(Edir, i) Vzero(Fdir, i)
	U = 0; Udir = 0; 
	nCell = 20;
	m = 0;
	for (m = 0; m < NMOL; m++) {
	for (nc[0] = -nCell; nc[0] <= nCell; nc[0]++) {for (nc[1] = -nCell; nc[1] <= nCell; nc[1]++) {for (nc[2] = -nCell; nc[2] <= nCell; nc[2]++) {
		if (nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2] > nCell * nCell) continue;
		for (i = 0; i < 3; i++) dR.v[i] = dcell * nc[i];
		for (n = 0; n < NMOL; n++) {
			if (n == m && nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
			dr = r[n] - r[m] + dR;
			V3ABS2(dr, R2)
			R = sqrt(R2); kR = kappa * R;
			Ewald_f0(f0, kR, i)
			/*
			R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2; R9 = R7 * R2;
			kR = kappa * R; kR2 = kR * kR;

			Ewald_f(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4, f5, f6)

			MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_T3(mT3, dr, R2, R7, i, j, k, t1)
			MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2) MTP_D2(mD2, dr, i, j)

			EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j)
			EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)
			EWALD_T4(mET4, f1, f2, f3, f4, dr, mT2, mD2, mT3, mT4, i, j, k, l)
			*/

			if (nc[0] == 0 && nc[1] == 0 && nc[2] == 0) {
				U += q[m] * q[n] * f0 / R / 2;
				Udir += q[m] * q[n] / R / 2;
			}
			else {
				U += q[m] * q[n] * f0 / R;
				Udir += q[m] * q[n] / R;
			}
		}
	}}}
	}

	nH = int(kappa * dcell + 0.5); if (nH < 1) nH = 1;
	Vzero(E_recip, i)
	double f_recip;
	U_recip = 0;
	m = 0;
	for (m = 0; m < NMOL; m++) {
	for (nc[0] = -nH; nc[0] <= nH; nc[0]++) {for (nc[1] = -nH; nc[1] <= nH; nc[1]++) {for (nc[2] = -nH; nc[2] <= nH; nc[2]++) {
		if (nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
		//if (nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2] > nH * nH) continue;
		for (i = 0; i < 3; i++) h.v[i] = PI2 / dcell * nc[i];
		V3ABS2(h, H2)
		hk = H2 * inv_k2;
		EXP(Ah, hk, i)
		//Ah = exp(-hk);
		if (Ah < 1e-8) continue;
		Ah *= inv_V / H2;
		scale_uv3(h, r[m], hR)
		PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
		f_recip = 0; 
		for (n = 0; n < NMOL; n++) {
			dr = r[n] - r[m];
			scale_uv3(h, dr, hR)
			PERIOD(hR, 0, PI2, t1) COS(fcos1, t1, i, t1) SIN(fsin1, t1, i, t1)
			//fcos1 = cos(hR); fsin1 = sin(hR);
			//U_recip += Ah * q[m] * q[n] * (fcos0 * fcos1 + fsin0 * fsin1);
			U_recip += Ah * q[m] * q[n] * fcos1;
			/*
			scale_uv3(mol[n].mu, h, muh)
			scale_VMV(h, mol[n].Q, h, Qh, i, j, 3)
			f_recip += -fsin0 * mol[n].q * fcos1 + fcos0 * mol[n].q * fsin1;
			f_recip -= fcos0 * muh * fcos1 + fsin0 * muh * fsin1;
			//f_recip += 2.0 / 3 * Ah * h.v[i] * (-fsin0 * Qh * fcos1 + fcos0 * Qh * fsin1);
			*/
		}
		//scale_uv(mol[0].mu, h, t1, i, 3)
		//f_recip += t1;
		//for (i = 0; i < 3; i++) E_recip.v[i] += Ah * h.v[i] * f_recip;
	}}}
	}
	U_recip /= 2;
	//for (i = 0; i < 3; i++) E_recip.v[i] += inv_V / 3 * mol[0].mu.v[i];

	//for (i = 0; i < 3; i++) E.v[i] += E_recip.v[i];

	double U_self = 0;
	for (m = 0; m < NMOL; m++) U_self -= kappa / sqrt(3.1415926) * q[m] * q[m];
	
	VECTOR3 mu0;
	for (m = 0; m < NMOL; m++) {for (i = 0; i < 3; i++) mu0.v[i] += q[m] * r[m].v[i];}
	U_self += PI / 3 / (dcell * dcell * dcell) * (mu0.v[0] * mu0.v[0] + mu0.v[1] * mu0.v[1] + mu0.v[2] * mu0.v[2]);

	double dU = Udir - U - U_recip;

	n = n;
}

void test_EwaldSum1() {
	int NMOL = 8;
	float q[8];
	VECTOR3 r[8];
	float dcell = 50;
	int i, j, k, l, m, n;
	VECTOR3 dr;
	float q0 = 0;
	for (n = 0; n < NMOL; n++) {
		r[n].v[0] = dcell * (ranf() - 0.5);
		r[n].v[1] = dcell * (ranf() - 0.5);
		r[n].v[2] = dcell* (ranf() - 0.5);
		q[n] = (n == 0 ? 1 : -q[n-1]);
		q0 += q[n];
	}

	int nc[3], nCell = 4, nH = 6;
	float kappa = 3 / dcell;
	VECTOR3 E, F, Edir, E_recip, Fdir, dR;
	double U, Udir, U_recip;

	VECTOR3 mT1, mET1;
	DMATRIX<3> mT2, mET2, mD2;
	CMATRIX<3> mT3, mET3;
	QMATRIX<3> mT4, mET4;

	double R, R2, R3, R4, R5, R7, R9;
	double kappa3 = kappa * kappa * kappa;
	double kR, kR2;
	double f0, f1, f2, f3, f4, f5, f6, fexp2, Ah;
	VECTOR3 h;
	double inv_k2 = 0.25 / (kappa * kappa), inv_V = 4 * PI / (dcell * dcell * dcell);
	double fcos0, fsin0, fcos1, fsin1, hR, hk, muh, H2;
	double t1, t2;

	Vzero(E, i) Vzero(F, i) Vzero(Edir, i) Vzero(Fdir, i)
	U = 0; Udir = 0;
	nCell = 20;
	m = 0;
	for (nc[0] = -nCell; nc[0] <= nCell; nc[0]++) {for (nc[1] = -nCell; nc[1] <= nCell; nc[1]++) {for (nc[2] = -nCell; nc[2] <= nCell; nc[2]++) {
		if (nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2] > nCell * nCell) continue;
		for (i = 0; i < 3; i++) dR.v[i] = dcell * nc[i];
		for (n = 0; n < NMOL; n++) {
			if (n == m && nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
			dr = r[m] - r[n] - dR;
			V3ABS2(dr, R2)
			R = sqrt(R2); kR = kappa * R;
			Ewald_f0(f0, kR, i)
			/*
			R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2; R9 = R7 * R2;
			kR = kappa * R; kR2 = kR * kR;

			Ewald_f(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4, f5, f6)

			MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_T3(mT3, dr, R2, R7, i, j, k, t1)
			MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2) MTP_D2(mD2, dr, i, j)

			EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j)
			EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)
			EWALD_T4(mET4, f1, f2, f3, f4, f5, f6, dr, mT2, mD2, mT3, mT4, i, j, k, l)
			*/

			U += q[m] * q[n] * f0 / R / 2;
			Udir += q[m] * q[n] / R / 2;
		}
	}}}

	nH = int(kappa * dcell + 0.5); if (nH < 1) nH = 1;
	Vzero(E_recip, i)
	double f_recip;
	U_recip = 0;
	m = 0;
	for (nc[0] = -nH; nc[0] <= nH; nc[0]++) {for (nc[1] = -nH; nc[1] <= nH; nc[1]++) {for (nc[2] = -nH; nc[2] <= nH; nc[2]++) {
		if (nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
		if (nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2] > nH * nH) continue;
		for (i = 0; i < 3; i++) h.v[i] = PI2 / dcell * nc[i];
		V3ABS2(h, H2)
		hk = H2 * inv_k2;
		EXP(Ah, hk, i)
		//Ah = exp(-hk);
		if (Ah < 1e-8) continue;
		Ah *= inv_V / H2;
		scale_uv3(h, r[m], hR)
		PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
		f_recip = 0;
		for (n = 0; n < NMOL; n++) {
			dr = r[m] - r[n];
			scale_uv3(h, dr, hR)
			PERIOD(hR, 0, PI2, t1) COS(fcos1, t1, i, t1) SIN(fsin1, t1, i, t1)
			//fcos1 = cos(hR); fsin1 = sin(hR);
			//U_recip += Ah * q[m] * q[n] * (fcos0 * fcos1 + fsin0 * fsin1);
			U_recip += Ah * q[m] * q[n] * fcos1;
			/*
			scale_uv3(mol[n].mu, h, muh)
			scale_VMV(h, mol[n].Q, h, Qh, i, j, 3)
			f_recip += -fsin0 * mol[n].q * fcos1 + fcos0 * mol[n].q * fsin1;
			f_recip -= fcos0 * muh * fcos1 + fsin0 * muh * fsin1;
			//f_recip += 2.0 / 3 * Ah * h.v[i] * (-fsin0 * Qh * fcos1 + fcos0 * Qh * fsin1);
			*/
		}
		//scale_uv(mol[0].mu, h, t1, i, 3)
		//f_recip += t1;
		//for (i = 0; i < 3; i++) E_recip.v[i] += Ah * h.v[i] * f_recip;
	}}}
	U_recip /= 2;
	//for (i = 0; i < 3; i++) E_recip.v[i] += inv_V / 3 * mol[0].mu.v[i];

	//for (i = 0; i < 3; i++) E.v[i] += E_recip.v[i];

	double U_self = 0;
	U_self -= kappa / sqrt(3.1415926) * q[m] * q[m];
	
	VECTOR3 mu0;
	for (n = 0; n < NMOL; n++) {for (i = 0; i < 3; i++) mu0.v[i] += q[n] * r[n].v[i];}
	U_self += PI / 3 / (dcell * dcell * dcell) * (q[m] * r[m].v[0] * mu0.v[0] + q[m] * r[m].v[1] * mu0.v[1] + q[m] * r[m].v[2] * mu0.v[2]);
	
	double dU = Udir - U - U_recip;

	n = n;
}

void test() {
	efield_test();
	//test_EwaldSum_EF();
	//test_EwaldSum();
	//test_EwaldSum1();
}

#endif


