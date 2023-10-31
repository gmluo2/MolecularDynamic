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

#include "EwaldSum.h"

namespace _EwaldSum_real_ {

void EwaldSumRealVars0::init_vars(VECTOR3& dr, double R, double R2) {
	int i;
	this->dr = dr.v; this->R = R; this->R2 = R2; R3 = R2 * R; 
	MTP_T1(mT1, dr, R3) 
	if (bEwald) {
		kR = kappa * R;
		EXP2(fexp2, kR, i) Ewald_f0(f0, kR, i) Ewald_f1(f1, kR, f0, fexp2)
		EWALD_T1(mET1, f1, mT1) 
	}
}

void EwaldSumRealVars1::init_vars(VECTOR3& dr, double R, double R2) {
	int i, j, k;
	double t1;
	this->dr = dr.v; this->R = R; this->R2 = R2; 
	R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2;// R9 = R7 * R2;
	MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
	MTP_T3(mT3, dr, R2, R7, i, j, k, t1)

	if (bEwald) {
		kR = kappa * R; kR2 = kR * kR;
		//fexp2 = exp(-kR2);
		Ewald_f04(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4)
		EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j) 
		EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)
	}
}

void EwaldSumRealVars1::init_vars_E(VECTOR3& dr, double R, double R2) {
	int i, j, k;
	double t1;
	this->dr = dr.v; this->R = R; this->R2 = R2; 
	R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2;// R9 = R7 * R2;
	MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
	//MTP_T3(mT3, dr, R2, R7, i, j, k, t1)

	if (bEwald) {
		kR = kappa * R; kR2 = kR * kR;
		//fexp2 = exp(-kR2);
		//Ewald_f04(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4)
		Ewald_f02(kappa3, R, R2, kR, kR2, i, fexp2, f0, f1, f2)
		EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j) 
		//EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)
	}
}

void EwaldSumRealVars2::init_vars(VECTOR3& dr, double R, double R2) {
	int i, j, k, l;
	double t1, t2;
	this->dr = dr.v; this->R = R; this->R2 = R2; 
	R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2; R9 = R7 * R2;
	MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
	MTP_T3(mT3, dr, R2, R7, i, j, k, t1)
	MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2) MTP_D2(mD2, dr, i, j)

	if (bEwald) {
		kR = kappa * R; kR2 = kR * kR;
		//fexp2 = exp(-kR2);
		Ewald_f(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4, f5, f6)
		EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j) 
		EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)
		EWALD_T4(mET4, f1, f2, f3, f4, f5, f6, dr, mT2, mD2, mT3, mT4, i, j, k, l)
	}
}

void cMTP_E(EwaldSumRealVars0 &esv, double c2, VECTOR3 &E, bool bEwaldKickoff, double a) {
	V3zero(E)
	int i;
	double ac2;
	if (bEwaldKickoff) {ac2 = a * c2;}
	if (esv.bEwald) {// ewald sum
		E.v[0] = -c2 * esv.mET1.v[0];
		E.v[1] = -c2 * esv.mET1.v[1];
		E.v[2] = -c2 * esv.mET1.v[2];

		if (bEwaldKickoff) {
			E.v[0] += ac2 * esv.mT1.v[0];
			E.v[1] += ac2 * esv.mT1.v[1];
			E.v[2] += ac2 * esv.mT1.v[2];
		}
	}
	else {// direct interaction
		E.v[0] = -c2 * esv.mT1.v[0];
		E.v[1] = -c2 * esv.mT1.v[1];
		E.v[2] = -c2 * esv.mT1.v[2];
	}
}

void MTP_Interact(EwaldSumRealVars0 &esv, double c1, double c2, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a) {
	U = 0; V3zero(E) V3zero(F)
	int i;
	double c12 = c1 * c2, Ut = c12 / esv.R; // q1 - q2
	double ac2, ac12;
	if (bEwaldKickoff) {ac2 = a * c2; ac12 = a * c12;}
	if (esv.bEwald) {// ewald sum
		U = Ut * esv.f0;
		E.v[0] = -c2 * esv.mET1.v[0];
		E.v[1] = -c2 * esv.mET1.v[1];
		E.v[2] = -c2 * esv.mET1.v[2];
		F.v[0] = -c12 * esv.mET1.v[0];
		F.v[1] = -c12 * esv.mET1.v[1];
		F.v[2] = -c12 * esv.mET1.v[2];

		if (bEwaldKickoff) {
			U -= Ut * a;
			E.v[0] += ac2 * esv.mT1.v[0];
			E.v[1] += ac2 * esv.mT1.v[1];
			E.v[2] += ac2 * esv.mT1.v[2];
			F.v[0] += ac12 * esv.mT1.v[0];
			F.v[1] += ac12 * esv.mT1.v[1];
			F.v[2] += ac12 * esv.mT1.v[2];
		}
	}
	else {// direct interaction
		U = Ut;
		E.v[0] = -c2 * esv.mT1.v[0];
		E.v[1] = -c2 * esv.mT1.v[1];
		E.v[2] = -c2 * esv.mT1.v[2];
		F.v[0] = -c12 * esv.mT1.v[0];
		F.v[1] = -c12 * esv.mT1.v[1];
		F.v[2] = -c12 * esv.mT1.v[2];
	}
	U *= 0.5;
}

void MTP_E(EwaldSumRealVars1 &esv, double c2, double *mu2, VECTOR3 &E, bool bEwaldKickoff, double a) {
	V3zero(E)
	int i, j, k;
	double ac2;
	double amu2[3];
	if (bEwaldKickoff) {
		ac2 = a * c2;
		amu2[0] = a * mu2[0]; amu2[1] = a * mu2[1]; amu2[2] = a * mu2[2];
	}

	if (esv.bEwald) {
		// E, -- q2, mu2
		for (i = 0; i < 3; i++) {
			// E from q2
			E.v[i] -= c2 * esv.mET1.v[i];
			for (j = 0; j < 3; j++) {
				// E from mu2
				E.v[i] += esv.mET2.m[i][j] * mu2[j];
			}
		}

		if (bEwaldKickoff) {
			// E, - q2, mu2
			for (i = 0; i < 3; i++) {
				// E from q2
				E.v[i] += ac2 * esv.mT1.v[i];
				for (j = 0; j < 3; j++) {
					// E from mu2
					E.v[i] -= esv.mT2.m[i][j] * amu2[j];
				}
			}
		}
	}
	else {
		// E, - q2, mu2
		for (i = 0; i < 3; i++) {
			// E from q2
			E.v[i] -= c2 * esv.mT1.v[i];
			for (j = 0; j < 3; j++) {
				// E from mu2
				E.v[i] += esv.mT2.m[i][j] * mu2[j];
			}
		}
	}
}

void MTP_Interact(EwaldSumRealVars1 &esv, double c1, double *mu1, double c2, double *mu2, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a) {
	U = 0; V3zero(E) V3zero(F)
	int i, j, k;
	double c12 = c1 * c2, Ut = c12 / esv.R; // q1 - q2
	double ac1, ac2, ac12;
	double amu1[3], amu2[3];
	if (bEwaldKickoff) {
		ac1 = a * c1; ac2 = a * c2; ac12 = a * c12;
		amu1[0] = a * mu1[0]; amu1[1] = a * mu1[1]; amu1[2] = a * mu1[2];
		amu2[0] = a * mu2[0]; amu2[1] = a * mu2[1]; amu2[2] = a * mu2[2];
	}

	if (esv.bEwald) {
		U = Ut * esv.f0; // q - q
		// mu1 - q2
		U += c2 * (mu1[0] * esv.mET1.v[0] + mu1[1] * esv.mET1.v[1] + mu1[2] * esv.mET1.v[2]);
		// mu2 - q1
		U -= c1 * (mu2[0] * esv.mET1.v[0] + mu2[1] * esv.mET1.v[1] + mu2[2] * esv.mET1.v[2]);
		// mu1 - mu2
		U -= mu1[0] * ( esv.mET2.m[0][0] * mu2[0] + esv.mET2.m[0][1] * mu2[1] + esv.mET2.m[0][2] * mu2[2])
					+ mu1[1] * ( esv.mET2.m[1][0] * mu2[0] + esv.mET2.m[1][1] * mu2[1] + esv.mET2.m[1][2] * mu2[2])
					+ mu1[2] * ( esv.mET2.m[2][0] * mu2[0] + esv.mET2.m[2][1] * mu2[1] + esv.mET2.m[2][2] * mu2[2]);

		// E & F, q1, mu1 - q2, mu2
		for (i = 0; i < 3; i++) {
			// E from q2
			E.v[i] -= c2 * esv.mET1.v[i];
			// force
			F.v[i] -= c12 * esv.mET1.v[i]; // q1 - q2
			for (j = 0; j < 3; j++) {
				// E from mu2
				E.v[i] += esv.mET2.m[i][j] * mu2[j];
				// force
				F.v[i] += c1 * esv.mET2.m[i][j] * mu2[j]; // q1 - mu2
				F.v[i] -= c2 * esv.mET2.m[i][j] * mu1[j]; // q2 - mu1
				for (k = 0; k < 3; k++) { // mu1 - mu2
					F.v[i] += mu1[j] * esv.mET3.m[i][j][k] * mu2[k];
				}
			}
		}

		if (bEwaldKickoff) {
			U -= Ut * a; // q - q
			U -= ac2 * (mu1[0] * esv.mT1.v[0] + mu1[1] * esv.mT1.v[1] + mu1[2] * esv.mT1.v[2]); // mu1 - q2
			U += ac1 * (mu2[0] * esv.mT1.v[0] + mu2[1] * esv.mT1.v[1] + mu2[2] * esv.mT1.v[2]); // q1 - mu2
			// mu1 - mu2
			U +=  amu1[0] * (esv.mT2.m[0][0] * mu2[0] + esv.mT2.m[0][1] * mu2[1] + esv.mT2.m[0][2] * mu2[2])
				+ amu1[1] * (esv.mT2.m[1][0] * mu2[0] + esv.mT2.m[1][1] * mu2[1] + esv.mT2.m[1][2] * mu2[2])
				+ amu1[2] * (esv.mT2.m[2][0] * mu2[0] + esv.mT2.m[2][1] * mu2[1] + esv.mT2.m[2][2] * mu2[2]);

			// E & F, q1, mu1 - q2, mu2
			for (i = 0; i < 3; i++) {
				// E from q2
				E.v[i] += ac2 * esv.mT1.v[i];
				// force
				F.v[i] += ac12 * esv.mT1.v[i]; // q1 - q2
				for (j = 0; j < 3; j++) {
					// E from mu2
					E.v[i] -= esv.mT2.m[i][j] * amu2[j];
					// force
					F.v[i] -= ac1 * esv.mT2.m[i][j] * mu2[j]; // q1 - mu2
					F.v[i] += ac2 * esv.mT2.m[i][j] * mu1[j]; // q2 - mu1
					for (k = 0; k < 3; k++) { // mu1 - mu2
						F.v[i] -= amu1[j] * esv.mT3.m[i][j][k] * mu2[k];
					}
				}
			}
		}
	}
	else {
		U = Ut; // q - q
		U += c2 * (mu1[0] * esv.mT1.v[0] + mu1[1] * esv.mT1.v[1] + mu1[2] * esv.mT1.v[2]); // mu1 - q2
		U -= c1 * (mu2[0] * esv.mT1.v[0] + mu2[1] * esv.mT1.v[1] + mu2[2] * esv.mT1.v[2]); // q1 - mu2
		// mu1 - mu2
		U -=  mu1[0] * (esv.mT2.m[0][0] * mu2[0] + esv.mT2.m[0][1] * mu2[1] + esv.mT2.m[0][2] * mu2[2])
			+ mu1[1] * (esv.mT2.m[1][0] * mu2[0] + esv.mT2.m[1][1] * mu2[1] + esv.mT2.m[1][2] * mu2[2])
			+ mu1[2] * (esv.mT2.m[2][0] * mu2[0] + esv.mT2.m[2][1] * mu2[1] + esv.mT2.m[2][2] * mu2[2]);

		// E & F, q1, mu1 - q2, mu2
		for (i = 0; i < 3; i++) {
			// E from q2
			E.v[i] -= c2 * esv.mT1.v[i];
			// force
			F.v[i] -= c12 * esv.mT1.v[i]; // q1 - q2
			for (j = 0; j < 3; j++) {
				// E from mu2
				E.v[i] += esv.mT2.m[i][j] * mu2[j];
				// force
				F.v[i] += c1 * esv.mT2.m[i][j] * mu2[j]; // q1 - mu2
				F.v[i] -= c2 * esv.mT2.m[i][j] * mu1[j]; // q2 - mu1
				for (k = 0; k < 3; k++) { // mu1 - mu2
					F.v[i] += mu1[j] * esv.mT3.m[i][j][k] * mu2[k];
				}
			}
		}
	}
	U *= 0.5;
}

void MTP_Interact(EwaldSumRealVars1 &esv, MATOM *patom1, MATOM *patom2, VECTOR3& dr, double R, double R2, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a, bool bExcludeInducedDipole) {
	U = 0; V3zero(E) V3zero(F)
	bool bCharge = (patom1->mMP == 0 && patom2->mMP == 0 ? true : false); // interaction between point charges, without multipole ?
	if (bCharge) {
		esv.EwaldSumRealVars0::init_vars(dr, R, R2);
		MTP_Interact((*(EwaldSumRealVars0*)(&esv)), patom1->c, patom2->c, E, F, U, bEwaldKickoff, a);
	}
	else {
		esv.init_vars(dr, R, R2);
		MTP_Interact(esv, patom1->c, patom1->mu.v, patom2->c, patom2->mu.v, E, F, U, bEwaldKickoff, a);
		if (bExcludeInducedDipole) {
			{
				VECTOR3 E1, F1;
				double U1 = 0;
				MTP_Interact(esv, 0, patom1->ind_mu.v, 0, patom2->ind_mu.v, E1, F1, U1, bEwaldKickoff, a);
				V3minusV3(F, F1, F)
				U -= U1;
			}
		}
	}
}

void MTP_E(EwaldSumRealVars1 &esv, MATOM *patom2, VECTOR3& dr, double R, double R2, VECTOR3 &E, bool bEwaldKickoff, double a) {
	V3zero(E)
	bool bCharge = (patom2->mMP == 0 ? true : false); // interaction between point charges, without multipole ?
	if (bCharge) {
		esv.EwaldSumRealVars0::init_vars(dr, R, R2);
		cMTP_E(*(EwaldSumRealVars0*)(&esv), (double)(patom2->c), E, bEwaldKickoff, a);
	}
	else {
		esv.init_vars(dr, R, R2);
		MTP_E(esv, patom2->c, patom2->mu.v, E, bEwaldKickoff, a);
	}
}

void MTP_E(EwaldSumRealVars2 &esv, double c2, double *mu2, DMATRIX<3>& Q, VECTOR3 &E, bool bEwaldKickoff, double a) {
	V3zero(E)
	int i, j, k, l;
	double ac2;
	double amu2[3];
	if (bEwaldKickoff) {
		ac2 = a * c2;
		amu2[0] = a * mu2[0]; amu2[1] = a * mu2[1]; amu2[2] = a * mu2[2];
	}

	if (esv.bEwald) {
		// E, - q2, mu2
		for (i = 0; i < 3; i++) {
			// E from q2
			E.v[i] -= c2 * esv.mET1.v[i];
			for (j = 0; j < 3; j++) {
				// E from mu2
				E.v[i] += esv.mET2.m[i][j] * mu2[j];
				// E from Q2
				for (k = 0; k < 3; k++) E.v[i] -= Q.m[j][k] * esv.mET3.m[i][j][k] / 3; 
			}
		}

		if (bEwaldKickoff) {
			// E, - q2, mu2
			for (i = 0; i < 3; i++) {
				// E from q2
				E.v[i] += ac2 * esv.mT1.v[i];
				for (j = 0; j < 3; j++) {
					// E from mu2
					E.v[i] -= esv.mT2.m[i][j] * amu2[j];
					// E from Q2
					for (k = 0; k < 3; k++) E.v[i] += a * Q.m[j][k] * esv.mT3.m[i][j][k] / 3;
				}
			}
		}
	}
	else {
		// E, - q2, mu2
		for (i = 0; i < 3; i++) {
			// E from q2
			E.v[i] -= c2 * esv.mT1.v[i];
			for (j = 0; j < 3; j++) {
				// E from mu2
				E.v[i] += esv.mT2.m[i][j] * mu2[j];
				// E from Q2
				for (k = 0; k < 3; k++) E.v[i] -= Q.m[j][k] * esv.mT3.m[i][j][k] / 3;
			}
		}
	}
}

void MTP_Interact(EwaldSumRealVars2 &esv, double c1, double *mu1, double c2, double *mu2, DMATRIX<3>& Q, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a) {
	U = 0; V3zero(E) V3zero(F)
	int i, j, k, l;
	double c12 = c1 * c2, Ut = c12 / esv.R; // q1 - q2
	double ac1, ac2, ac12;
	double amu1[3], amu2[3];
	if (bEwaldKickoff) {
		ac1 = a * c1; ac2 = a * c2; ac12 = a * c12;
		amu1[0] = a * mu1[0]; amu1[1] = a * mu1[1]; amu1[2] = a * mu1[2];
		amu2[0] = a * mu2[0]; amu2[1] = a * mu2[1]; amu2[2] = a * mu2[2];
	}

	if (esv.bEwald) {
		U = Ut * esv.f0; // q - q
		// mu1 - q2
		U += c2 * (mu1[0] * esv.mET1.v[0] + mu1[1] * esv.mET1.v[1] + mu1[2] * esv.mET1.v[2]);
		// mu2 - q1
		U -= c1 * (mu2[0] * esv.mET1.v[0] + mu2[1] * esv.mET1.v[1] + mu2[2] * esv.mET1.v[2]);
		// mu1 - mu2
		U -= mu1[0] * ( esv.mET2.m[0][0] * mu2[0] + esv.mET2.m[0][1] * mu2[1] + esv.mET2.m[0][2] * mu2[2])
					+ mu1[1] * ( esv.mET2.m[1][0] * mu2[0] + esv.mET2.m[1][1] * mu2[1] + esv.mET2.m[1][2] * mu2[2])
					+ mu1[2] * ( esv.mET2.m[2][0] * mu2[0] + esv.mET2.m[2][1] * mu2[1] + esv.mET2.m[2][2] * mu2[2]);
		for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
			U += c1 * esv.mET2.m[i][j] * Q.m[i][j] / 3; // q-quadrupole
			for (k = 0; k < 3; k++) U += mu1[i] * esv.mET3.m[i][j][k] * Q.m[j][k] / 3; // mu-quadrupole
		}}

		// E & F, q1, mu1 - q2, mu2
		for (i = 0; i < 3; i++) {
			// E from q2
			E.v[i] -= c2 * esv.mET1.v[i];
			// force
			F.v[i] -= c12 * esv.mET1.v[i]; // q1 - q2
			for (j = 0; j < 3; j++) {
				// E from mu2
				E.v[i] += esv.mET2.m[i][j] * mu2[j];
				// E from Q2
				for (k = 0; k < 3; k++) E.v[i] -= Q.m[j][k] * esv.mET3.m[i][j][k] / 3; 
				// force
				F.v[i] += c1 * esv.mET2.m[i][j] * mu2[j]; // q1 - mu2
				F.v[i] -= c2 * esv.mET2.m[i][j] * mu1[j]; // q2 - mu1
				for (k = 0; k < 3; k++) { 
					F.v[i] += mu1[j] * esv.mET3.m[i][j][k] * mu2[k]; // mu1 - mu2
					F.v[i] -= c1 * esv.mET3.m[i][j][k] * Q.m[j][k] / 3; // q1 - Q2
					for (l = 0; l < 3; l++) {
						F.v[i] -= mu1[j] * esv.mET4.m[i][j][k][l] * Q.m[k][l] / 3; // mu1 - Q2
					}
				}
			}
		}

		if (bEwaldKickoff) {
			U -= Ut * a; // q - q
			U -= ac2 * (mu1[0] * esv.mT1.v[0] + mu1[1] * esv.mT1.v[1] + mu1[2] * esv.mT1.v[2]); // mu1 - q2
			U += ac1 * (mu2[0] * esv.mT1.v[0] + mu2[1] * esv.mT1.v[1] + mu2[2] * esv.mT1.v[2]); // q1 - mu2
			// mu1 - mu2
			U +=  amu1[0] * (esv.mT2.m[0][0] * mu2[0] + esv.mT2.m[0][1] * mu2[1] + esv.mT2.m[0][2] * mu2[2])
				+ amu1[1] * (esv.mT2.m[1][0] * mu2[0] + esv.mT2.m[1][1] * mu2[1] + esv.mT2.m[1][2] * mu2[2])
				+ amu1[2] * (esv.mT2.m[2][0] * mu2[0] + esv.mT2.m[2][1] * mu2[1] + esv.mT2.m[2][2] * mu2[2]);
			for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
				U -= ac1 * esv.mT2.m[i][j] * Q.m[i][j] / 3; // q-quadrupole
				for (k = 0; k < 3; k++) U -= amu1[i] * esv.mT3.m[i][j][k] * Q.m[j][k] / 3; // mu-quadrupole
			}}

			// E & F, q1, mu1 - q2, mu2
			for (i = 0; i < 3; i++) {
				// E from q2
				E.v[i] += ac2 * esv.mT1.v[i];
				// force
				F.v[i] += ac12 * esv.mT1.v[i]; // q1 - q2
				for (j = 0; j < 3; j++) {
					// E from mu2
					E.v[i] -= esv.mT2.m[i][j] * amu2[j];
					// E from Q2
					for (k = 0; k < 3; k++) E.v[i] += a * Q.m[j][k] * esv.mT3.m[i][j][k] / 3;
					// force
					F.v[i] -= ac1 * esv.mT2.m[i][j] * mu2[j]; // q1 - mu2
					F.v[i] += ac2 * esv.mT2.m[i][j] * mu1[j]; // q2 - mu1
					for (k = 0; k < 3; k++) { // mu1 - mu2
						F.v[i] -= amu1[j] * esv.mT3.m[i][j][k] * mu2[k]; //mu1 - mu2
						F.v[i] += ac1 * esv.mT3.m[i][j][k] * Q.m[j][k] / 3; // q1 - Q2
						for (l = 0; l < 3; l++) {
							F.v[i] += amu1[j] * esv.mT4.m[i][j][k][l] * Q.m[k][l] / 3; // mu1 - Q2
						}
					}
				}
			}
		}
	}
	else {
		U = Ut; // q - q
		U += c2 * (mu1[0] * esv.mT1.v[0] + mu1[1] * esv.mT1.v[1] + mu1[2] * esv.mT1.v[2]); // mu1 - q2
		U -= c1 * (mu2[0] * esv.mT1.v[0] + mu2[1] * esv.mT1.v[1] + mu2[2] * esv.mT1.v[2]); // q1 - mu2
		// mu1 - mu2
		U -=  mu1[0] * (esv.mT2.m[0][0] * mu2[0] + esv.mT2.m[0][1] * mu2[1] + esv.mT2.m[0][2] * mu2[2])
			+ mu1[1] * (esv.mT2.m[1][0] * mu2[0] + esv.mT2.m[1][1] * mu2[1] + esv.mT2.m[1][2] * mu2[2])
			+ mu1[2] * (esv.mT2.m[2][0] * mu2[0] + esv.mT2.m[2][1] * mu2[1] + esv.mT2.m[2][2] * mu2[2]);
		for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
			U += c1 * esv.mT2.m[i][j] * Q.m[i][j] / 3; // q-quadrupole
			for (k = 0; k < 3; k++) U += mu1[i] * esv.mT3.m[i][j][k] * Q.m[j][k] / 3; // mu-quadrupole
		}}

		// E & F, q1, mu1 - q2, mu2
		for (i = 0; i < 3; i++) {
			// E from q2
			E.v[i] -= c2 * esv.mT1.v[i];
			// force
			F.v[i] -= c12 * esv.mT1.v[i]; // q1 - q2
			for (j = 0; j < 3; j++) {
				// E from mu2
				E.v[i] += esv.mT2.m[i][j] * mu2[j];
				// E from Q2
				for (k = 0; k < 3; k++) E.v[i] -= Q.m[j][k] * esv.mT3.m[i][j][k] / 3;
				// force
				F.v[i] += c1 * esv.mT2.m[i][j] * mu2[j]; // q1 - mu2
				F.v[i] -= c2 * esv.mT2.m[i][j] * mu1[j]; // q2 - mu1
				for (k = 0; k < 3; k++) { // mu1 - mu2
					F.v[i] += mu1[j] * esv.mT3.m[i][j][k] * mu2[k];
					F.v[i] -= c1 * esv.mT3.m[i][j][k] * Q.m[j][k] / 3; // q1 - Q2
					for (l = 0; l < 3; l++) {
						F.v[i] -= mu1[j] * esv.mT4.m[i][j][k][l] * Q.m[k][l] / 3; // mu1 - Q2
					}
				}
			}
		}
	}
	U *= 0.5;
}

} // end of namespace _evatom_

namespace _EwaldSum_recp_ {
using namespace _EwaldSum_real_;
double EwaldSum_Recip(MMOL_MD_CELL &mdcell) {
#define _NH 4
#define NH 9
#define AH_MIN 5e-3

	float AhMin = (float)1e-8;
	float V = mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8;
	double kappa, kappa3, inv_3V;

	InitEwaldSumPars(V, inv_3V, kappa, kappa3);

	complex f;
	VECTOR3 h, dr;
	double inv_k2 = 0.25 / (kappa * kappa);
	double fcos0, fsin0, fcos1, fsin1, hR, hk, H2, Ah;
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
					for (i = 0; i < 3; i++) r.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
					scale_uv3(h, r, hR)
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
				scale_uv3(h, r, hR)
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
				scale_uv3(h, r, hR)
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

double EwaldSum_Recip(_VE_CELL<_VQatom> &vc) {
#define _NH 4
#define NH 9
#define AH_MIN 5e-3

	float AhMin = (float)1e-8;

	float V = (vc.xd[0] * vc.xd[1] * vc.xd[2]);
	double kappa, kappa3, inv_3V;

	InitEwaldSumPars(V, inv_3V, kappa, kappa3);

	complex f;
	VECTOR3 h, dr;
	double inv_k2 = 0.25 / (kappa * kappa);
	double hR, hk, H2, Ah;
	ARRAY<double> fcos0, fsin0;
	fcos0.SetArray(vc.va.n); fsin0.SetArray(vc.va.n); 
	double t1, t2;
	double U_recip, f_recip;
	int nH = int(kappa * vc.xd[0] + 0.5); if (nH < 1) nH = 1;
	int nc[3];
	U_recip = 0;
	int iatom, i;
	_VQatom *patom = NULL;

	double EF = 0;

	for (iatom = 0; iatom < vc.va.n; iatom++) {
		patom = vc.va.m + iatom;
		V3zero(patom->E_recp)
		V3zero(patom->F_recp)
	}

	nH = 20;
	VECTOR3 r;
	float q;
	for (nc[0] = -nH; nc[0] <= nH; nc[0]++) {for (nc[1] = -nH; nc[1] <= nH; nc[1]++) {for (nc[2] = -nH; nc[2] <= nH; nc[2]++) {
	//for (nc[0] = 0; nc[0] <= nH; nc[0]++) {for (nc[1] = 0; nc[1] <= nH; nc[1]++) {for (nc[2] = 0; nc[2] <= nH; nc[2]++) {
		if (nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
		//if ((nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2]) > nH * nH) continue;
		for (i = 0; i < 3; i++) h.v[i] = PI2 / vc.xd[i] * nc[i];
		V3ABS2(h, H2)
		hk = H2 * inv_k2;
		//EXP(Ah, hk, i)
		Ah = exp(-hk);
		if (Ah < 1e-8) continue;
		Ah *= PI2 / V / H2;

		f.x = 0; f.y = 0;
		for (iatom = 0; iatom < vc.va.n; iatom++) {
			patom = vc.va.m + iatom; q = *(patom->q);
			scale_uv3(h, (*(patom->r)), hR)
			PERIOD(hR, 0, PI2, t1) COS(fcos0.m[iatom], t1, i, t1) SIN(fsin0.m[iatom], t1, i, t1)
			f.x += q * fcos0.m[iatom]; f.y += q * fsin0.m[iatom];
		}
		U_recip += Ah * (f.x * f.x + f.y * f.y);

		// electric field
		for (iatom = 0; iatom < vc.va.n; iatom++) {
			patom = vc.va.m + iatom; q = *(patom->q);
			//scale_uv3(h, (*(patom->r)), hR)
			//PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
			//fcos1 = cos(hR); fsin1 = sin(hR);

			EF =  fsin0.m[iatom] * f.x - fcos0.m[iatom] * f.y;
			for (i = 0; i < 3; i++) {
				patom->E_recp.v[i] += 2 * Ah * h.v[i] * EF;
				patom->F_recp.v[i] += q * 2 * Ah * h.v[i] * EF;
			}
		}
	}}}


	return U_recip * eUnit_estat_kT;
}


double EwaldSum_Recip(_VE_CELL<_VMUatom> &vc) {
#define _NH 4
#define NH 9
#define AH_MIN 5e-3

	float AhMin = (float)1e-8;

	float V = (vc.xd[0] * vc.xd[1] * vc.xd[2]);
	double kappa, kappa3, inv_3V;

	InitEwaldSumPars(V, inv_3V, kappa, kappa3);

	complex f;
	VECTOR3 h, dr;
	double inv_k2 = 0.25 / (kappa * kappa);
	double hR, hk, H2, Ah;
	ARRAY<double> fcos0, fsin0;
	fcos0.SetArray(vc.va.n); fsin0.SetArray(vc.va.n); 
	double t1, t2;
	double U_recip, f_recip;
	int nH = int(kappa * vc.xd[0] + 0.5); if (nH < 1) nH = 1;
	int nc[3];
	U_recip = 0;
	int iatom, i;
	_VMUatom *patom = NULL;

	double EF = 0;

	for (iatom = 0; iatom < vc.va.n; iatom++) {
		patom = vc.va.m + iatom;
		V3zero(patom->E_recp)
		V3zero(patom->F_recp)
	}

	nH = 20;
	VECTOR3 r;
	float q;
	ARRAY<double> muk;
	muk.SetArray(vc.va.n);
	for (nc[0] = -nH; nc[0] <= nH; nc[0]++) {for (nc[1] = -nH; nc[1] <= nH; nc[1]++) {for (nc[2] = -nH; nc[2] <= nH; nc[2]++) {
	//for (nc[0] = 0; nc[0] <= nH; nc[0]++) {for (nc[1] = 0; nc[1] <= nH; nc[1]++) {for (nc[2] = 0; nc[2] <= nH; nc[2]++) {
		if (nc[0] == 0 && nc[1] == 0 && nc[2] == 0) continue;
		//if ((nc[0] * nc[0] + nc[1] * nc[1] + nc[2] * nc[2]) > nH * nH) continue;
		for (i = 0; i < 3; i++) h.v[i] = PI2 / vc.xd[i] * nc[i];
		V3ABS2(h, H2)
		hk = H2 * inv_k2;
		//EXP(Ah, hk, i)
		Ah = exp(-hk);
		if (Ah < 1e-8) continue;
		Ah *= PI2 / V / H2;

		f.x = 0; f.y = 0;
		for (iatom = 0; iatom < vc.va.n; iatom++) {
			patom = vc.va.m + iatom; q = *(patom->q);
			scale_uv3(h, (*(patom->r)), hR)
			PERIOD(hR, 0, PI2, t1) COS(fcos0.m[iatom], t1, i, t1) SIN(fsin0.m[iatom], t1, i, t1)
			f.x += q * fcos0.m[iatom]; f.y += q * fsin0.m[iatom];

			// dipole term
			muk.m[iatom] = h.v[0] * patom->mu->v[0] + h.v[1] * patom->mu->v[1] + h.v[2] + patom->mu->v[2];
			f.x -= muk.m[iatom] * fsin0.m[iatom]; f.y = muk.m[iatom] * fcos0.m[iatom];
		}
		U_recip += Ah * (f.x * f.x + f.y * f.y);

		// electric field
		for (iatom = 0; iatom < vc.va.n; iatom++) {
			patom = vc.va.m + iatom; q = *(patom->q);
			//scale_uv3(h, (*(patom->r)), hR)
			//PERIOD(hR, 0, PI2, t1) COS(fcos0, t1, i, t1) SIN(fsin0, t1, i, t1)
			//fcos1 = cos(hR); fsin1 = sin(hR);

			EF =  fsin0.m[iatom] * f.x - fcos0.m[iatom] * f.y;
			EF -= muk.m[iatom];
			for (i = 0; i < 3; i++) {
				patom->E_recp.v[i] += 2 * Ah * h.v[i] * EF;
				patom->F_recp.v[i] += q * 2 * Ah * h.v[i] * EF;
			}
		}
	}}}


	return U_recip * eUnit_estat_kT;
}

} // end of namespace _EwaldSum_recp_


namespace _Thole_polarize_ {
using namespace _EwaldSum_real_;
bool TholeCorrTensor::init_vars_0(double a, double R, double R2, double *dr, bool bF) {
	// r-exponential screening function (exp(-a * r))
	double r = a * R;
	if (r > this->xcut_exp) return false;
	double fexp = 0;
	int i, j, k;
	interpolate(fexp, ::ExpProf, r, i);
	if (fexp == 0) return false;

	double r2 = r * r, R3 = R * R2, R5 = R2 * R3, R7 = R2 * R5;
	double g1 = (1 + r + 0.5 * r2) * fexp;
	double g2 = (1 + r + 0.5 * r2 + r * r2 / 6) * fexp;
	double g3 = g2 + r2 * r2 / 30 * fexp;

	f = -(1 + 0.5 * r) * fexp / R;
	
	double *v1 = mT1.v, inv = g1 / R3;
	v1[0] = dr[0] * inv;
	v1[1] = dr[1] * inv;
	v1[2] = dr[2] * inv;

	inv = -3 * g2 / R5;
	double inv1 = g1 / R3;
	v1 = mT2.m[0]; 
	v1[0] = dr[0] * dr[0] * inv + inv1;
	v1[1] = dr[0] * dr[1] * inv;
	v1[2] = dr[0] * dr[2] * inv;

	v1 = mT2.m[1]; 
	v1[0] = dr[1] * dr[0] * inv;
	v1[1] = dr[1] * dr[1] * inv + inv1;
	v1[2] = dr[1] * dr[2] * inv;

	v1 = mT2.m[2]; 
	v1[0] = dr[2] * dr[0] * inv;
	v1[1] = dr[2] * dr[1] * inv;
	v1[2] = dr[2] * dr[2] * inv + inv1;

	if (!bF) return true;

	double t = 0;
	inv = g3 * 15 / R7; inv1 = -g2 * 3 / R5;
	for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
		mT3.m[i][j][k] = inv * dr[i] * dr[j] * dr[k];
		t = 0;
		if (i == j) t += dr[k];
		if (j == k) t += dr[i];
		if (i == k) t += dr[j];
		if (t != 0) mT3.m[i][j][k] += inv1 * t;
	}}}
	return true;
};

bool TholeCorrTensor::init_vars_1(double a, double R, double R2, double *dr, bool bF) {
	// r3-exponential screening function (exp(-a * r^3))
	double R3 = R * R2, R5 = R2 * R3, R7 = R2 * R5, r3 = a * R3, r6 = r3 * r3;
	if (r3 > xcut_exp) return false;

	int i, j, k;
	double fexp = 0, gamma = 0;
	interpolate(fexp, ::ExpProf, r3, i);
	if (fexp == 0) return false;
	interpolate(gamma, ::IncGammaProf_2_3, r3, i);

	double g1 = fexp;
	double g2 = (1 + r3) * fexp;
	double g3 = (1 + r3 + 0.6 * r6) * fexp;

	f = (pow(r3, 1./3) * gamma - fexp) / R;
	
	double *v1 = mT1.v, inv = g1 / R3;
	v1[0] = dr[0] * inv;
	v1[1] = dr[1] * inv;
	v1[2] = dr[2] * inv;

	inv = -3 * g2 / R5;
	double inv1 = g1 / R3;
	v1 = mT2.m[0]; 
	v1[0] = dr[0] * dr[0] * inv + inv1;
	v1[1] = dr[0] * dr[1] * inv;
	v1[2] = dr[0] * dr[2] * inv;

	v1 = mT2.m[1]; 
	v1[0] = dr[1] * dr[0] * inv;
	v1[1] = dr[1] * dr[1] * inv + inv1;
	v1[2] = dr[1] * dr[2] * inv;

	v1 = mT2.m[2]; 
	v1[0] = dr[2] * dr[0] * inv;
	v1[1] = dr[2] * dr[1] * inv;
	v1[2] = dr[2] * dr[2] * inv + inv1;

	if (!bF) return true;

	double t = 0;
	inv = g3 * 15 / R7; inv1 = -g2 * 3 / R5;
	for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
		mT3.m[i][j][k] = inv * dr[i] * dr[j] * dr[k];
		t = 0;
		if (i == j) t += dr[k];
		if (j == k) t += dr[i];
		if (i == k) t += dr[j];
		if (t != 0) mT3.m[i][j][k] += inv1 * t;
	}}}
	return true;
};

// correction of electric field from atom2 [q2, mu2] to atom1 [q1, mu1], with dr = r1 - r2; // 2==>1
// Correction: [in format of paper by Andres Aguado et al, J. Chem. Phys., 119, 7471 (2003)]
// U = f * q1 * q2;
// E & F are from 2 to 1, with dr = r1 - r2 (from r2 to 1)
// E = -q2 * mT1 + mT2 * mu2;  
// F = -q1 * q2 * mT1 + q1 * mT2 * mu2 - q2 * mT2 * mu1
bool TholeCorr_E(TholeCorrTensor &tt, double a, double q2, double *mu2, double R, double R2, double *dr, double *E) {
	DV3zero(E);
	if (!tt.init_vars(a, R, R2, dr, false)) return false;
	double *m = tt.mT1.v;
	E[0] = -q2 * m[0];
	E[1] = -q2 * m[1];
	E[2] = -q2 * m[2];

	m = tt.mT2.m[0]; E[0] += m[0] * mu2[0] + m[1] * mu2[1] + m[2] * mu2[2];
	m = tt.mT2.m[1]; E[1] += m[0] * mu2[0] + m[1] * mu2[1] + m[2] * mu2[2];
	m = tt.mT2.m[2]; E[2] += m[0] * mu2[0] + m[1] * mu2[1] + m[2] * mu2[2];
	return true;
}

// correction of electric field, force and U from atom2 [q2, mu2] to atom1 [q1, mu1], with dr = r1 - r2; // 2==>1
// Correction: [in format of paper by Andres Aguado et al, J. Chem. Phys., 119, 7471 (2003)]
// U = f * q1 * q2;
// E & F are from 2 to 1, with dr = r1 - r2 (from r2 to 1)
// E = -q2 * mT1 + mT2 * mu2;  
// F = -q1 * q2 * mT1 + q1 * mT2 * mu2 - q2 * mT2 * mu1
bool TholeCorr(TholeCorrTensor &tt, double a, double q1, double *mu1, double q2, double *mu2, double R, double R2, double *dr, double *E, double *F, double &U) {
	DV3zero(E) DV3zero(F) U = 0;
	if (!tt.init_vars(a, R, R2, dr, true)) return false;
	double *m = tt.mT1.v;
	double q12 = q1 * q2;
	U = q12 * tt.f;
	U -= q1 * (m[0] * mu2[0] + m[1] * mu2[1] + m[2] * mu2[2]);
	U += q2 * (m[0] * mu1[0] + m[1] * mu1[1] + m[2] * mu1[2]);

	U -= mu1[0] * (tt.mT2.m[0][0] * mu2[0] + tt.mT2.m[0][1] * mu2[1] + tt.mT2.m[0][2] * mu2[2])
		+ mu1[1] * (tt.mT2.m[1][0] * mu2[0] + tt.mT2.m[1][1] * mu2[1] + tt.mT2.m[1][2] * mu2[2])
		+ mu1[2] * (tt.mT2.m[2][0] * mu2[0] + tt.mT2.m[2][1] * mu2[1] + tt.mT2.m[2][2] * mu2[2]);

	m = tt.mT1.v;
	E[0] = -q2 * m[0]; F[0] = -q12 * m[0];
	E[1] = -q2 * m[1]; F[1] = -q12 * m[1];
	E[2] = -q2 * m[2]; F[2] = -q12 * m[2];

	double t;
	m = tt.mT2.m[0]; t = m[0] * mu2[0] + m[1] * mu2[1] + m[2] * mu2[2]; E[0] += t; F[0] += q1 * t;
	m = tt.mT2.m[1]; t = m[0] * mu2[0] + m[1] * mu2[1] + m[2] * mu2[2]; E[1] += t; F[1] += q1 * t;
	m = tt.mT2.m[2]; t = m[0] * mu2[0] + m[1] * mu2[1] + m[2] * mu2[2]; E[2] += t; F[2] += q1 * t;

	m = tt.mT2.m[0]; F[0] -= q2 * (m[0] * mu1[0] + m[1] * mu1[1] + m[2] * mu1[2]);
	m = tt.mT2.m[1]; F[1] -= q2 * (m[0] * mu1[0] + m[1] * mu1[1] + m[2] * mu1[2]);
	m = tt.mT2.m[2]; F[2] -= q2 * (m[0] * mu1[0] + m[1] * mu1[1] + m[2] * mu1[2]);

	int i, j, k;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
			F[i] += mu1[j] * tt.mT3.m[i][j][k] * mu2[k];
		}}
	}
	return true;
}

void initTholeRadius(ATOM_COMM_PARS_DATABASE &atompar_db, PAIRWISE_DB< PAIRWISE<float> > &TholeRadius_db) {
	extern short iTholeModel;

	char msg[256] = "\0";
	TholeRadius_db.releaseChain(true);

	CHAIN<ATOM_COMM_PARS> *ch1, *ch2;
	ATOM_COMM_PARS *apar1, *apar2;

	char aindx[2] = {0x00, 0x00};
	unsigned int key = 0;
	PAIRWISE<float> *pw = NULL;
	float radius = 0;

	ch1 = atompar_db.ch;
	while (ch1 != NULL) {
		apar1 = ch1->p;
		if (apar1->alpha <= 0) {ch1 = ch1->next; continue;}
		ch2 = atompar_db.ch;
		while (ch2 != NULL) {
			apar2 = ch2->p;
			if (apar2->alpha <= 0) {ch2 = ch2->next; continue;}
			aindx[0] = apar1->indx; aindx[1] = apar2->indx;
			key = construct_char2_key_pairwise(aindx);
			pw = TholeRadius_db.search_chain(key);
			if (pw == NULL) {
				radius = pow(double(apar1->alpha * apar2->alpha), 1./6);
				pw = new PAIRWISE<float>;
				pw->aindx[0] = apar1->indx; pw->aindx[1] = apar2->indx; pw->key = key;
				if (iTholeModel == 0) pw->pv = ::_Thole_polarize_::a / radius;
				else if (iTholeModel == 1) pw->pv = ::_Thole_polarize_::a / (radius * radius * radius);
				TholeRadius_db.attach(pw);
				pw = NULL;
				sprintf(msg, "Thole-radius [%s - %s] : %f Angs.", apar1->atom, apar2->atom, radius);
				show_log(msg, true);
			}

			ch2 = ch2->next;
		}
		ch1 = ch1->next;
	}
	TholeRadius_db.Chain2Array();
	show_log("", true);
	return;
}

} // end of namespace _Thole_Polarize_
