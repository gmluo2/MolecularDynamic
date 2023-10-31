#include "project.h"

#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
/*
#include "show.h"
#include "ranlib.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"
*/

#ifdef __GPU__

//#include "atom_sr.h"
#include "gpu_vector.h"
#include "gpu_esr.h"

#include "gpu_sr.h"

//extern double erfc(double x);

namespace _gpu_md_ {
	namespace _gpu_ewald_ {
	// following is copied from EwaldSum.cpp

	__device__ void EwaldSumRealVars0::init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2) {
		int i;
		//this->dr = dr.v; 
		this->R = R; this->R2 = R2; R3 = R2 * R; 
		MTP_T1(mT1, dr, R3) 
		if (bEwald) {
			kR = kappa * R;
			fexp2 = exp(-kR * kR); Ewald_f0(f0, kR); Ewald_f1(f1, kR, f0, fexp2)
			EWALD_T1(mET1, f1, mT1) 
		}
	}

	__device__ void EwaldSumRealVars0::init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2, DOUBLE f_kR, DOUBLE f_exp2) {
		int i;
		//this->dr = dr.v;
		this->R = R; this->R2 = R2; R3 = R2 * R;
		MTP_T1(mT1, dr, R3)
		if (bEwald) {
			kR = kappa * R;
			//fexp2 = exp(-kR * kR); Ewald_f0(f0, kR);
			fexp2 = f_exp2; f0 = f_kR;
			Ewald_f1(f1, kR, f0, fexp2)
			EWALD_T1(mET1, f1, mT1)
		}
	}


	__device__ void EwaldSumRealVars1::init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2) {
		int i, j, k;
		DOUBLE t1;
		//this->dr = dr.v; 
		this->R = R; this->R2 = R2; 
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

	__device__ void EwaldSumRealVars1::init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2, DOUBLE f_kR, DOUBLE f_exp2) {
		int i, j, k;
		DOUBLE t1;
		//this->dr = dr.v;
		this->R = R; this->R2 = R2;
		R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2;// R9 = R7 * R2;
		MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
		MTP_T3(mT3, dr, R2, R7, i, j, k, t1)

		if (bEwald) {
			kR = kappa * R; kR2 = kR * kR;
			//fexp2 = exp(-kR2);
			fexp2 = f_exp2; f0 = f_kR;
			_Ewald_f04_simplify(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4)
			EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j)
			EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)
		}
	}

	__device__ void EwaldSumRealVars1::init_vars_E(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2) {
		int i, j, k;
		DOUBLE t1;
		//this->dr = dr.v; 
		this->R = R; this->R2 = R2; 
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

	__device__ void EwaldSumRealVars1::init_vars_E(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2, DOUBLE f_kR, DOUBLE f_exp2) {
		int i, j, k;
		DOUBLE t1;
		//this->dr = dr.v;
		this->R = R; this->R2 = R2;
		R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2;// R9 = R7 * R2;
		MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
		//MTP_T3(mT3, dr, R2, R7, i, j, k, t1)

		if (bEwald) {
			kR = kappa * R; kR2 = kR * kR;
			//fexp2 = exp(-kR2);
			//Ewald_f04(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4)
			f0 = f_kR; fexp2 = f_exp2;
			_Ewald_f02_simplify(kappa3, R, R2, kR, kR2, i, fexp2, f0, f1, f2)
			EWALD_T1(mET1, f1, mT1) EWALD_T2(mET2, f1, f2, mT2, mD2, i, j)
			//EWALD_T3(mET3, f1, f2, f3, f4, dr, mT2, mD2, mT3, i, j, k)
		}
	}

	__device__ void EwaldSumRealVars2::init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2) {
		int i, j, k, l;
		DOUBLE t1, t2;
		//this->dr = dr.v; 
		this->R = R; this->R2 = R2; 
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



	// the code below is copied / modified from atom_sr_interact.cpp

	// interaction on patom1 from patom2, 
	// dr -- vector from patom2 to patom1
	// R -- distance |dr|
	// R -- |dr|^2
	// E -- electric field from patom2
	// F -- force from patom2 to patom1
	// U -- interaction energy
	__device__ void EwaldReal_Q(EwaldSumRealVars1 &esv, DOUBLE c1, DOUBLE c2, DOUBLE *dr, DOUBLE R, DOUBLE R2, /*double *E, */DOUBLE *F, DOUBLE &U) {
		double t, c12 = c1 * c2;
		esv.EwaldSumRealVars0::init_vars(dr, R, R2);

		// U
		t = c12 / R; // q1 - q2
		U = t * esv.f0;

		// E
		/*
		E.v[0] = -c2 * esv.mET1.v[0]; // from q2
		E.v[1] = -c2 * esv.mET1.v[1]; // from q2
		E.v[2] = -c2 * esv.mET1.v[2]; // from q2
		*/

		// F
		F[0] = -c12 * esv.mET1.v[0]; // q1 - q2
		F[1] = -c12 * esv.mET1.v[1]; // q1 - q2
		F[2] = -c12 * esv.mET1.v[2]; // q1 - q2
	}

	__device__ void EInteract_Q(_gpu_ewald_::EwaldSumRealVars1 &esv, DOUBLE c1, DOUBLE c2, DOUBLE *dr, DOUBLE R, DOUBLE R2, /*DOUBLE *E, */DOUBLE *F, DOUBLE &U, bool init_esv) {
		DOUBLE c12 = c1 * c2;

		if (init_esv) esv.EwaldSumRealVars0::init_vars(dr, R, R2);

		// U
		U = c12 / R; // q1 - q2
		//U = t * esv.f0;

		// E
		/*
		E.v[0] = -c2 * esv.mT1.v[0]; // from q2
		E.v[1] = -c2 * esv.mT1.v[1]; // from q2
		E.v[2] = -c2 * esv.mT1.v[2]; // from q2
		*/

		// F
		F[0] = -c12 * esv.mT1.v[0]; // q1 - q2
		F[1] = -c12 * esv.mT1.v[1]; // q1 - q2
		F[2] = -c12 * esv.mT1.v[2]; // q1 - q2
	}

	__device__ void EwaldReal_MU(EwaldSumRealVars1 &esv, DOUBLE c1, DOUBLE *mu1, DOUBLE c2, DOUBLE *mu2, DOUBLE *dr, DOUBLE R, DOUBLE R2, /*double *E,*/ DOUBLE *F, DOUBLE &U, char bExcludeDipole) {
		int i, j, k;
		DOUBLE t;

		if (bExcludeDipole) esv.init_vars_E(dr, R, R2); // no induced dipole and dipole-dipole interaction will be excluded
		else esv.init_vars(dr, R, R2);

		// U
		t = c1 * c2 / R; // q1 - q2
		U = t * esv.f0;

		U += c2 * (mu1[0] * esv.mET1.v[0] + mu1[1] * esv.mET1.v[1] + mu1[2] * esv.mET1.v[2]);
		U -= c1 * (mu2[0] * esv.mET1.v[0] + mu2[1] * esv.mET1.v[1] + mu2[2] * esv.mET1.v[2]);

		if (bExcludeDipole == 0) {
			U -= mu1[0] * ( esv.mET2.m[0][0] * mu2[0] + esv.mET2.m[0][1] * mu2[1] + esv.mET2.m[0][2] * mu2[2])
				+ mu1[1] * ( esv.mET2.m[1][0] * mu2[0] + esv.mET2.m[1][1] * mu2[1] + esv.mET2.m[1][2] * mu2[2])
				+ mu1[2] * ( esv.mET2.m[2][0] * mu2[0] + esv.mET2.m[2][1] * mu2[1] + esv.mET2.m[2][2] * mu2[2]);
		}

		// E
		/*
		for (i = 0; i < 3; i++) {
			E.v[i] -= c2 * esv.mET1.v[i]; // from q2
			for (j = 0; j < 3; j++) E.v[i] += esv.mET2.m[i][j] * mu2[j];
		}
		*/

		// F
		for (i = 0; i < 3; i++) {
			F[i] = -c1 * c2 * esv.mET1.v[i]; // q1 - q2
			for (j = 0; j < 3; j++) {
				F[i] += c1 * esv.mET2.m[i][j] * mu2[j];
				F[i] -= c2 * esv.mET2.m[i][j] * mu1[j];
			}
			if (!bExcludeDipole) { // mu1 - mu2
				for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
					F[i] += mu1[j] * esv.mET3.m[i][j][k] * mu2[k];
				}}
			}
		}
	}

	__device__ void EInteract_MU(EwaldSumRealVars1 &esv, DOUBLE c1, DOUBLE *mu1, DOUBLE c2, DOUBLE *mu2, DOUBLE *dr, DOUBLE R, DOUBLE R2, /*DOUBLE *E,*/ DOUBLE *F, DOUBLE &U, char bExcludeDipole, bool init_esv) {
		int i, j, k;

		if (init_esv) {
			if (bExcludeDipole) esv.init_vars_E(dr, R, R2); // no induced dipole and dipole-dipole interaction will be excluded
			else esv.init_vars(dr, R, R2);
		}

		// U
		U = c1 * c2 / R; // q1 - q2
		//U = t * esv.f0;

		U += c2 * (mu1[0] * esv.mT1.v[0] + mu1[1] * esv.mT1.v[1] + mu1[2] * esv.mT1.v[2]);
		U -= c1 * (mu2[0] * esv.mT1.v[0] + mu2[1] * esv.mT1.v[1] + mu2[2] * esv.mT1.v[2]);

		if (bExcludeDipole == 0) {
			U -= mu1[0] * ( esv.mT2.m[0][0] * mu2[0] + esv.mT2.m[0][1] * mu2[1] + esv.mT2.m[0][2] * mu2[2])
				+ mu1[1] * ( esv.mT2.m[1][0] * mu2[0] + esv.mT2.m[1][1] * mu2[1] + esv.mT2.m[1][2] * mu2[2])
				+ mu1[2] * ( esv.mT2.m[2][0] * mu2[0] + esv.mT2.m[2][1] * mu2[1] + esv.mT2.m[2][2] * mu2[2]);
		}

		// E
		/*
		for (i = 0; i < 3; i++) {
			E.v[i] -= c2 * esv.mET1.v[i]; // from q2
			for (j = 0; j < 3; j++) E.v[i] += esv.mET2.m[i][j] * mu2[j];
		}
		*/

		// F
		for (i = 0; i < 3; i++) {
			F[i] = -c1 * c2 * esv.mT1.v[i]; // q1 - q2
			for (j = 0; j < 3; j++) {
				F[i] += c1 * esv.mT2.m[i][j] * mu2[j];
				F[i] -= c2 * esv.mT2.m[i][j] * mu1[j];
			}
			if (!bExcludeDipole) { // mu1 - mu2
				for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
					F[i] += mu1[j] * esv.mT3.m[i][j][k] * mu2[k];
				}}
			}
		}
	}


	// calculate the electric field only about the real part of Ewald Sum
	__device__ void EwaldReal_E_Q(EwaldSumRealVars1 &esv, DOUBLE c2, DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE *E) {
		esv.EwaldSumRealVars0::init_vars(dr, R, R2);

		// E
		E[0] = -c2 * esv.mET1.v[0]; // from q2
		E[1] = -c2 * esv.mET1.v[1];
		E[2] = -c2 * esv.mET1.v[2];
	}

	__device__ void EInteract_E_Q(EwaldSumRealVars1 &esv, DOUBLE c2, DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE *E, bool init_esv) {
		if (init_esv) esv.EwaldSumRealVars0::init_vars(dr, R, R2);

		// E
		E[0] = -c2 * esv.mT1.v[0]; // from q2
		E[1] = -c2 * esv.mT1.v[1];
		E[2] = -c2 * esv.mT1.v[2];
	}

	__device__ void EwaldReal_E_MU(EwaldSumRealVars1 &esv, DOUBLE c2, DOUBLE *mu2, DOUBLE* dr, DOUBLE R, DOUBLE R2, DOUBLE *E) {
		int i, j;
		esv.init_vars_E(dr, R, R2);

		// E
		for (i = 0; i < 3; i++) {
			E[i] = -c2 * esv.mET1.v[i]; // from q2
			for (j = 0; j < 3; j++) E[i] += esv.mET2.m[i][j] * mu2[j];
		}
	}

	__device__ void EInteract_E_MU(EwaldSumRealVars1 &esv, DOUBLE c2, DOUBLE *mu2, DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE *E, bool init_esv) {
		int i, j;
		if (init_esv) esv.init_vars_E(dr, R, R2);

		// E
		for (i = 0; i < 3; i++) {
			E[i] = -c2 * esv.mT1.v[i]; // from q2
			for (j = 0; j < 3; j++) E[i] += esv.mT2.m[i][j] * mu2[j];
		}
	}

	} // namespace _gpu_ewald_

	__global__ void EwaldSum_Real(GPU_EAtCELL *acell) {
		// works on 1 dimensional grid and one dimensional blocks
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int Na = acell->Na;
		int nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;
		
		int i, ia1, ia, ia2, cIndx1, cIndx2, i3a1, i3a2;
		int icx, icy, icz, ncx, ncy, ncz, dcx, dcy, dcz;
		int Ncx, Ncy, Nc;

		DOUBLE R, R2;
		DOUBLE t1;
		DOUBLE xd[3] = {acell->cw[0] * acell->Nx, acell->cw[1] * acell->Ny, acell->cw[2] * acell->Nz};
		DOUBLE doff[3];

		DOUBLE U_estat = 0;
		DOUBLE dr[3];//, u[3];
		DOUBLE c1, c2, mu1[3], mu2[3], E[3], F[3], Et[3], Ft[3], Ut;
		
		DOUBLE U, Uestat = 0;

		int Nm, ic, ica1, ica2;

		int hx = acell->Nx / 2, hy = acell->Ny, hz = acell->Nz / 2;
		int *cubic;

		char bExcludeDipole = (char)(acell->esv.iExIndDipole);
		int cal_mode = acell->bCalE;
		if (cal_mode <= 0 || cal_mode >= 3) return; // with charge only, so far choice 1 -- electric field; 2 -- force and energy
		if (acell->bDipole) cal_mode += 10; // with dipole

		_gpu_ewald_::EwaldSumRealVars1 esv;
		//esv.cp_EwaldSum_vars(*(EwaldSumRealVars0*)(&acell->esv));
		esv.bEwald = true;
		esv.init(true, acell->esv.iSurface, acell->esv.iSurfaceBoundary, acell->esv.iExIndDipole, acell->esv.kappa, acell->esv.V);
		
		// direct interaction
		_gpu_ewald_::EwaldSumRealVars1 esv0;
		//esv.cp_EwaldSum_vars(*(EwaldSumRealVars0*)(&acell->esv));
		esv0.init(false, acell->esv.iSurface, acell->esv.iSurfaceBoundary, acell->esv.iExIndDipole, acell->esv.kappa, acell->esv.V);

		DOUBLE rcut = acell->rcut, r2cut = rcut * rcut, RMAX = rcut * 10, R2MAX = RMAX * RMAX;
		DOUBLE cw_max = acell->cw[0];
		if (acell->cw[1] > cw_max) cw_max = acell->cw[1];
		if (acell->cw[2] > cw_max) cw_max = acell->cw[2];
		DOUBLE rcut_m = rcut + cw_max * 2, r2cut_m = rcut_m * rcut_m;
		DOUBLE d2x, d2y, d2z, d2xy, d2;

		bool bIgnore = true;
		
		dcx = int(rcut / acell->cw[0]) + 1;//2;
		dcy = int(rcut / acell->cw[1]) + 1;
		dcz = int(rcut / acell->cw[2]) + 1;

		int drIdx, srIdx, nsr, n_srnb = acell->n_srnb;
		int *isr = NULL;
		DOUBLE *drsr = NULL;

		bool bInitializedSRNeighbor = ((acell->srnb_dev != NULL && acell->bSRInitialized) ? true : false);

		for (i = 0; i < nEachThread; i++) {
			ia1 = tid * nEachThread + i; i3a1 = x_idx(ia1);
			if (ia1 >= Na) break;

			c1 = acell->c_dev[ia1]; mu1[0] = acell->mu_dev[i3a1]; mu1[1] = acell->mu_dev[i3a1 + 1]; mu1[2] = acell->mu_dev[i3a1 + 2];
			cIndx1 = acell->cIndx_dev[ia1];
			
			if (acell->srnb_dev != NULL) {
				srIdx = acell->srIndx(ia1);
				drIdx = acell->drIndx(ia1);
				isr = acell->srnb_dev + srIdx;
				drsr = acell->dr_dev + drIdx;
			}

			Uestat = 0;
			if (bInitializedSRNeighbor) {
				Nc = isr[0];
				for (ia = 0; ia < Nc; ia++) {
					ia2 = isr[ia + 1];
//					if (ia2 >= Na) continue;
//					if (ia2 == ia1) continue;
					i3a2 = x_idx(ia2);
					cIndx2 = acell->cIndx_dev[ia2];
					//if (acell->aID_dev[ia2] < 0) continue;
					c2 = acell->c_dev[ia2]; mu2[0] = acell->mu_dev[i3a2]; mu2[1] = acell->mu_dev[i3a2 + 1]; mu2[2] = acell->mu_dev[i3a2 + 2];
					// drsr is from center atom to ia2, while dr want that from ia2 to center atom
					dr[0] = -drsr[0]; dr[1] = -drsr[1]; dr[2] = -drsr[2]; R = drsr[3];
					//if (R2 > r2cut) continue;
					R2 = R * R;

					if (R > rcut) continue;

					// calculation of the interaction
					// important, we are assuming rcut >> size of cluster
					switch (cal_mode) {
					case 1: // electric field only
						EwaldReal_E_Q(esv, c2, dr, R, R2, E);
						//if (acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) { // excluding that from the atoms inside the cluster
						if (cIndx1 == cIndx2) {
							esv.bEwald = false; EInteract_E_Q(esv, c2, dr, R, R2, Et, false);
							//EInteract_E_Q(esv0, c2, dr, R, R2, E, true);
							E[0] -= Et[0]; E[1] -= Et[1]; E[2] -= Et[2];
							esv.bEwald = true;
						}
						acell->E_dev[i3a1] += E[0]; acell->E_dev[i3a1 + 1] += E[1]; acell->E_dev[i3a1 + 2] += E[2];
						break;
					case 2: // force and U
						EwaldReal_Q(esv, c1, c2, dr, R, R2, F, U);
						//if (blocal && acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) { // excluding that from the atoms inside the cluster
						if (cIndx1 == cIndx2) {
							esv.bEwald = false; EInteract_Q(esv, c1, c2, dr, R, R2, Ft, Ut, false);
							//EInteract_Q(esv0, c1, c2, dr, R, R2, Ft, Ut, true);
							F[0] -= Ft[0]; F[1] -= Ft[1]; F[2] -= Ft[2];
							U -= Ut;
							esv.bEwald = true;
						}
						acell->Fe_dev[i3a1] += F[0]; acell->Fe_dev[i3a1 + 1] += F[1]; acell->Fe_dev[i3a1 + 2] += F[2];
						Uestat += U * 0.5;
						// virial
						acell->ve_dev[i3a1] += F[0] * dr[0] * 0.5;
						acell->ve_dev[i3a1 + 1] += F[1] * dr[1] * 0.5;
						acell->ve_dev[i3a1 + 2] += F[2] * dr[2] * 0.5;
						break;
					case 11:
						EwaldReal_E_MU(esv, c2, mu2, dr, R, R2, E);
						//if (blocal && acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) { // excluding that from the atoms inside the cluster
						if (cIndx1 == cIndx2) {
							esv.bEwald = false; EInteract_E_MU(esv, c2, mu2, dr, R, R2, Et, false);
							//EInteract_E_MU(esv0, c2, mu2, dr, R, R2, Et, true);
							E[0] -= Et[0]; E[1] -= Et[1]; E[2] -= Et[2];
							esv.bEwald = true;
						}
						acell->E_dev[i3a1] += E[0]; acell->E_dev[i3a1 + 1] += E[1]; acell->E_dev[i3a1 + 2] += E[2];
						break;
					case 12: // force and U
						EwaldReal_MU(esv, c1, mu1, c2, mu2, dr, R, R2, F, U, (char)bExcludeDipole);
						//if (blocal && acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) { // excluding that from the atoms inside the cluster
						if (cIndx1 == cIndx2) {
							esv.bEwald = false; EInteract_MU(esv, c1, mu1, c2, mu2, dr, R, R2, Ft, Ut, (char)bExcludeDipole, false);
							//EInteract_MU(esv0, c1, mu1, c2, mu2, dr, R, R2, Ft, Ut, (char)bExcludeDipole, true);
							F[0] -= Ft[0]; F[1] -= Ft[1]; F[2] -= Ft[2];
							U -= Ut;
							esv.bEwald = true;
						}
						acell->Fe_dev[i3a1] += F[0]; acell->Fe_dev[i3a1 + 1] += F[1]; acell->Fe_dev[i3a1 + 2] += F[2];
						Uestat += U * 0.5;
						// virial
						acell->ve_dev[i3a1] += F[0] * dr[0] * 0.5;
						acell->ve_dev[i3a1 + 1] += F[1] * dr[1] * 0.5;
						acell->ve_dev[i3a1 + 2] += F[2] * dr[2] * 0.5;
						break;
					default:
						break;
					}
				}
			}
			else {
				for (icx = -dcx; icx <= dcx; icx++) {
					ncx = acell->ir_dev[i3a1] + icx;
					if (ncx < 0) {ncx += acell->Nx; doff[0] = -xd[0];}
					else if (ncx >= acell->Nx) {ncx -= acell->Nx; doff[0] = xd[0];}
					else doff[0] = 0;
					//Ncx = ncx * acell->Nyz;

					d2x = icx * acell->cw[0]; d2x *= d2x;
			
					for (icy = -dcy; icy <= dcy; icy++) {
						ncy = acell->ir_dev[i3a1 + 1] + icy;
						if (ncy < 0) {ncy += acell->Ny; doff[1] = -xd[1];}
						else if (ncy >= acell->Ny) {ncy -= acell->Ny; doff[1] = xd[1];}
						else doff[1] = 0;
						//Ncy = ncy * acell->Nz;

						d2y = icy * acell->cw[1]; d2y *= d2y;
						d2xy = d2x + d2y;
						if (d2xy > r2cut_m) continue;

						for (icz = -dcz; icz <= dcz; icz++) {
							ncz = acell->ir_dev[i3a1 + 2] + icz;
							if (ncz < 0) {ncz += acell->Nz; doff[2] = -xd[2];}
							else if (ncz >= acell->Nz) {ncz -= acell->Nz; doff[2] = xd[2];}
							else doff[2] = 0;

							d2z = icz * acell->cw[2]; d2z *= d2z;
							d2 = d2xy + d2z;
							if (d2 > r2cut_m) continue;
						
							//Nc = (Ncx + Ncy + ncz) * (acell->NC + 1);
							//for (ia = 0; ia < acell->cmm_dev[Nc]; ia++) {
							//ia2 = acell->cmm_dev[Nc + ia + 1];
							
							cubic = acell->device_cubic(ncx, ncy, ncz);
							Nc = cubic[0];
							for (ia = 0; ia < Nc; ia++) {
								ia2 = cubic[ia + 1];
#ifdef __IGNORE_THIS_LINE__
								if (ia2 >= Na) continue;
								if (ia2 == ia1) continue;
#else
								if (ia2 >= Na) continue;
								if (ia2 == ia1) continue;
								bIgnore = false;
#endif
								i3a2 = x_idx(ia2);
								cIndx2 = acell->cIndx_dev[ia2];
								//if (acell->aID_dev[ia2] < 0) continue;
								c2 = acell->c_dev[ia2]; mu2[0] = acell->mu_dev[i3a2]; mu2[1] = acell->mu_dev[i3a2 + 1]; mu2[2] = acell->mu_dev[i3a2 + 2];
								dr[0] = acell->r_dev[i3a1] - acell->r_dev[i3a2] - doff[0];
								dr[1] = acell->r_dev[i3a1 + 1] - acell->r_dev[i3a2 + 1] - doff[1];
								dr[2] = acell->r_dev[i3a1 + 2] - acell->r_dev[i3a2 + 2] - doff[2];
								R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
								if (R2 > r2cut) continue;
								R = GPUSQRT(R2);

								// if the atoms are too close, too big induced dipole could be generated
								// so we have a minimum distance here, although this is an artificial effect
								//if (R < 0.5 && cIndx1 != cIndx2) {
								//	dr[0] *= 0.5 / R; dr[1] *= 0.5 / R; dr[2] *= 0.5 / R; R = 0.5; R2 = R * R;
								//}
//#ifdef __IGNORE_THIS_LINE__
								//if (R > rcut) continue;
							//u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;
//#else
//							if (ia2 >= Na || ia2 == ia1 || R > rcut) bIgnore = true;
//							if (bIgnore) { // issuing a position further than rcut, the direction can be any direction
//								R = RMAX; R2 = R2MAX; dr[0] = R; dr[1] = 0; dr[2] = 0;
//							}
//#endif

								// calculation of the interaction
								// important, we are assuming rcut >> size of cluster
								switch (cal_mode) {
								case 1: // electric field only
									EwaldReal_E_Q(esv, c2, dr, R, R2, E);
									//if (acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) { // excluding that from the atoms inside the cluster
									if (cIndx1 == cIndx2) {
										esv.bEwald = false; EInteract_E_Q(esv, c2, dr, R, R2, Et, false);
										//EInteract_E_Q(esv0, c2, dr, R, R2, E, true);
										E[0] -= Et[0]; E[1] -= Et[1]; E[2] -= Et[2];
										esv.bEwald = true;
									}
									acell->E_dev[i3a1] += E[0]; acell->E_dev[i3a1 + 1] += E[1]; acell->E_dev[i3a1 + 2] += E[2];
									break;
								case 2: // force and U
									EwaldReal_Q(esv, c1, c2, dr, R, R2, F, U);
									//if (blocal && acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) { // excluding that from the atoms inside the cluster
									if (cIndx1 == cIndx2) {
										esv.bEwald = false; EInteract_Q(esv, c1, c2, dr, R, R2, Ft, Ut, false);
										//EInteract_Q(esv0, c1, c2, dr, R, R2, Ft, Ut, true);
										F[0] -= Ft[0]; F[1] -= Ft[1]; F[2] -= Ft[2];
										U -= Ut;
										esv.bEwald = true;
									}
									acell->Fe_dev[i3a1] += F[0]; acell->Fe_dev[i3a1 + 1] += F[1]; acell->Fe_dev[i3a1 + 2] += F[2];
									Uestat += U * 0.5;
									// virial
									acell->ve_dev[i3a1] += F[0] * dr[0] * 0.5;
									acell->ve_dev[i3a1 + 1] += F[1] * dr[1] * 0.5;
									acell->ve_dev[i3a1 + 2] += F[2] * dr[2] * 0.5;
									break;
								case 11:
									EwaldReal_E_MU(esv, c2, mu2, dr, R, R2, E);
									//if (blocal && acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) { // excluding that from the atoms inside the cluster
									if (cIndx1 == cIndx2) {
										esv.bEwald = false; EInteract_E_MU(esv, c2, mu2, dr, R, R2, Et, false);
										//EInteract_E_MU(esv0, c2, mu2, dr, R, R2, Et, true);
										E[0] -= Et[0]; E[1] -= Et[1]; E[2] -= Et[2];
										esv.bEwald = true;
									}
									acell->E_dev[i3a1] += E[0]; acell->E_dev[i3a1 + 1] += E[1]; acell->E_dev[i3a1 + 2] += E[2];
									break;
								case 12: // force and U
									EwaldReal_MU(esv, c1, mu1, c2, mu2, dr, R, R2, F, U, (char)bExcludeDipole);
									//if (blocal && acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) { // excluding that from the atoms inside the cluster
									if (cIndx1 == cIndx2) {
										esv.bEwald = false; EInteract_MU(esv, c1, mu1, c2, mu2, dr, R, R2, Ft, Ut, (char)bExcludeDipole, false);
										//EInteract_MU(esv0, c1, mu1, c2, mu2, dr, R, R2, Ft, Ut, (char)bExcludeDipole, true);
										F[0] -= Ft[0]; F[1] -= Ft[1]; F[2] -= Ft[2];
										U -= Ut;
										esv.bEwald = true;
									}
									acell->Fe_dev[i3a1] += F[0]; acell->Fe_dev[i3a1 + 1] += F[1]; acell->Fe_dev[i3a1 + 2] += F[2];
									Uestat += U * 0.5;
									// virial
									acell->ve_dev[i3a1] += F[0] * dr[0] * 0.5;
									acell->ve_dev[i3a1 + 1] += F[1] * dr[1] * 0.5;
									acell->ve_dev[i3a1 + 2] += F[2] * dr[2] * 0.5;
									break;
								default:
									break;
								}
							}
						}
					}
				}
			}

			switch (cal_mode) {
			case 1: // q, electric field
				// no correction
				break;
			case 2:
				Uestat -= esv.kappa * INV_SQRT_PI * c1 * c1; // self energy from charge
				break;
			case 11: // electric field with dipole
				// self term of electric field from dipole
				t1 = 4 * esv.kappa3 / 3 * INV_SQRT_PI;
				acell->E_dev[i3a1] += t1 * mu1[0];
				acell->E_dev[i3a1 + 1] += t1 * mu1[1];
				acell->E_dev[i3a1 + 2] += t1 * mu1[2];
				break;
			case 12:
				Uestat -= esv.kappa * INV_SQRT_PI * c1 * c1; // self energy from charge
				// self energy from dipole
				Uestat -= 2 * esv.kappa3 / 3 * INV_SQRT_PI * (mu1[0] * mu1[0] + mu1[1] * mu1[1] + mu1[2] * mu1[2]);
				break;
			default:
				break;
			}

			// the electric field from the surface induced dipole, virial of surface dipole will be calculated later
			if (esv.iSurface == 0) { // cubic
				switch (esv.iSurfaceBoundary) { // normal
				case 0:
					if (cal_mode == 1 || cal_mode == 11) { // electric field only
						E[0] = -esv.inv_3V * acell->mu_surf[0];
						E[1] = -esv.inv_3V * acell->mu_surf[1];
						E[2] = -esv.inv_3V * acell->mu_surf[2];

						acell->E_dev[i3a1] += E[0];
						acell->E_dev[i3a1 + 1] += E[1];
						acell->E_dev[i3a1 + 2] += E[2];
					}
					else if (cal_mode == 2 || cal_mode == 12) {// force only. energy and virial will be calculated somewhere else
						if (esv.iExIndDipole == 0) {
							F[0] = -c1 * esv.inv_3V * acell->mu_surf[0];// * fUnit_estat_mass;
							F[1] = -c1 * esv.inv_3V * acell->mu_surf[1];// * fUnit_estat_mass;
							F[2] = -c1 * esv.inv_3V * acell->mu_surf[2];// * fUnit_estat_mass;
						}
						else if (esv.iExIndDipole == 1) {
							F[0] = -c1 * esv.inv_3V * (acell->mu_surf[0] - acell->ind_mu_surf[0]);// * fUnit_estat_mass;
							F[1] = -c1 * esv.inv_3V * (acell->mu_surf[1] - acell->ind_mu_surf[1]);// * fUnit_estat_mass;
							F[2] = -c1 * esv.inv_3V * (acell->mu_surf[2] - acell->ind_mu_surf[2]);// * fUnit_estat_mass;
						}

						acell->Fe_dev[i3a1] += F[0];
						acell->Fe_dev[i3a1 + 1] += F[1];
						acell->Fe_dev[i3a1 + 2] += F[2];
					}
					break;
				case 1: // disable surface dipole in all direction
					break;
				case 2: // disable surface dipole along z axis
					if (cal_mode == 1 || cal_mode == 11) { // electric field only
						E[0] = -esv.inv_3V * acell->mu_surf[0];
						E[1] = -esv.inv_3V * acell->mu_surf[1];
						//E[2] = -esv->inv_3V * acell->mu_surf[2];

						acell->E_dev[i3a1] += E[0];
						acell->E_dev[i3a1 + 1] += E[1];
						//acell->E_dev[i3a1 + 2] += E[2];
					}
					else if (cal_mode == 2 || cal_mode == 12) {// force only. energy and virial will be calculated somewhere else
						if (esv.iExIndDipole == 0) {
							F[0] = -c1 * esv.inv_3V * acell->mu_surf[0];// * fUnit_estat_mass;
							F[1] = -c1 * esv.inv_3V * acell->mu_surf[1];// * fUnit_estat_mass;
							//F[2] = -c1 * esv->inv_3V * acell->mu_surf[2];// * fUnit_estat_mass;
						}
						else if (esv.iExIndDipole == 1) {
							F[0] = -c1 * esv.inv_3V * (acell->mu_surf[0] - acell->ind_mu_surf[0]);// * fUnit_estat_mass;
							F[1] = -c1 * esv.inv_3V * (acell->mu_surf[1] - acell->ind_mu_surf[1]);// * fUnit_estat_mass;
							F[2] = -c1 * esv.inv_3V * (acell->mu_surf[2] - acell->ind_mu_surf[2]);// * fUnit_estat_mass;
						}

						acell->Fe_dev[i3a1] += F[0];
						acell->Fe_dev[i3a1 + 1] += F[1];
						//acell->Fe_dev[i3a1 + 2] += F[2];
					}
					break;
				default:
					break;
				}
			}
			//else if (esv.iSurface == 1) { // slab
				//switch (esv.iSurfaceBoundary) {
				//case 0: // normal
					// confirmed with explicit calculation, this boundary condition should be disabled
				//	break;
				//default: // do not do anything
				//	break;
				//}
			//}

			if (cal_mode == 2 || cal_mode == 12) acell->U_estat_dev[ia1] += Uestat;
		}
		
		__threadfence(); // its result can be seen in all other threads in this stream
		__syncthreads();
		

#ifdef _TO_CLUSTER_CALCULATED_
		// was it calculated ?
		// interactions inside cluster
		Nm = acell->Nc;
		nEachThread = Nm / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Nm) nEachThread += 1;

		//esv.bEwald = false;

		for (i = 0; i < nEachThread; i++) {
			ic = tid * nEachThread + i;
			if (ic >= acell->Nc) break;
			Nm = ic * (acell->Nca + 1);
			if (acell->cluster[Nm] <= 1) continue; // has only one atom
			for (ica1 = 0; ica1 < acell->cluster[Nm]; ica1++) {
				ia1 = acell->cluster[Nm + ica1 + 1];
				i3a1 = x_idx(ia1);
				c1 = acell->c_dev[ia1]; mu1[0] = acell->mu_dev[i3a1]; mu1[1] = acell->mu_dev[i3a1 + 1]; mu1[2] = acell->mu_dev[i3a1 + 2];
				Uestat = 0;
				for (ica2 = ica1 + 1; ica2 < acell->cluster[Nm]; ica2++) {
					ia2 = acell->cluster[Nm + ica2 + 1]; i3a2 = x_idx(ia2);
					c2 = acell->c_dev[ia2]; mu2[0] = acell->mu_dev[i3a2]; mu2[1] = acell->mu_dev[i3a2 + 1]; mu2[2] = acell->mu_dev[i3a2 + 2];
					// in this program, cluster is considered as unit in periodical coordination checking
					// so, we do not need to check close-neighbor cell
					dr[0] = acell->r_dev[i3a1] - acell->r_dev[i3a2];
					dr[1] = acell->r_dev[i3a1 + 1] - acell->r_dev[i3a2 + 1];
					dr[2] = acell->r_dev[i3a1 + 2] - acell->r_dev[i3a2 + 2];
					R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]; if (R2 < 0.2) R2 = 0.2;
					R = SQRT(R2);
					//if (R > acell->rcut) {continue; atomicAdd(acell->ir_map, 1);}
					//u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

					// calculate electrostatic interaction, excluding the interactions between atoms inside the same cluster
					switch (cal_mode) {
					case 1: // electric field only
						EInteract_E_Q(esv0, c2, dr, R, R2, E, true);
						acell->E_dev[i3a1] -= E[0]; acell->E_dev[i3a1 + 1] -= E[1]; acell->E_dev[i3a1 + 2] -= E[2]; // excluding
						acell->E_dev[i3a2] += E[0]; acell->E_dev[i3a2 + 1] += E[1]; acell->E_dev[i3a2 + 2] += E[2]; // excluding
						break;
					case 2: // force and U
						EInteract_Q(esv0, c1, c2, dr, R, R2, F, U, true);
						acell->Fe_dev[i3a1] -= F[0]; acell->Fe_dev[i3a1 + 1] -= F[1]; acell->Fe_dev[i3a1 + 2] -= F[2]; // excluding
						acell->Fe_dev[i3a2] += F[0]; acell->Fe_dev[i3a2 + 1] += F[1]; acell->Fe_dev[i3a2 + 2] += F[2]; // excluding
						Uestat -= U;
						// virial -- fully add to ia1
						acell->ve_dev[i3a1] -= F[0] * dr[0];
						acell->ve_dev[i3a1 + 1] -= F[1] * dr[1];
						acell->ve_dev[i3a1 + 2] -= F[2] * dr[2];
						break;
					case 11: 
						EInteract_E_MU(esv0, c2, mu2, dr, R, R2, E, true);
						acell->E_dev[i3a1] -= E[0]; acell->E_dev[i3a1 + 1] -= E[1]; acell->E_dev[i3a1 + 2] -= E[2];
						acell->E_dev[i3a2] += E[0]; acell->E_dev[i3a2 + 1] += E[1]; acell->E_dev[i3a2 + 2] += E[2];
						break;
					case 12: // force and U
						EInteract_MU(esv0, c1, mu1, c2, mu2, dr, R, R2, F, U, (char)bExcludeDipole, true);
						acell->Fe_dev[i3a1] -= F[0]; acell->Fe_dev[i3a1 + 1] -= F[1]; acell->Fe_dev[i3a1 + 2] -= F[2];
						acell->Fe_dev[i3a2] += F[0]; acell->Fe_dev[i3a2 + 1] += F[1]; acell->Fe_dev[i3a2 + 2] += F[2];
						Uestat -= U;
						// virial -- fully add to ia1
						acell->ve_dev[i3a1] -= F[0] * dr[0];
						acell->ve_dev[i3a1 + 1] -= F[1] * dr[1];
						acell->ve_dev[i3a1 + 2] -= F[2] * dr[2];
						break;
					default:
						break;
					}
				}
				acell->U_estat_dev[ia1] += Uestat; // Uestat is excluded energy, fully add to ia1
			}
		}

		__syncthreads();
		//if (threadIdx.x == 0) { // 0 thread in its block
			__threadfence(); // its result can be seen in all other thread, in all stream
		//}
#endif


		// should we calculate SR of bound and bound14 ? for those of macromolecule
		int ix1, iy1, iz1, ix2, iy2, iz2, dx, dy, dz;

		//tid = threadIdx.x + blockIdx.x * blockDim.x;
		//Na = acell->Na;
		DOUBLE rmax_bound = 3, Rmax_boundLength;
		int Nb = acell->Nb, ib, i3b;
		Rmax_boundLength = 3 * 2; // bound include 1-2 and 1-3 bound, latter could be twice of closest bound
		//esv.bEwald = false;
		if (Nb > 0) {
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;
			
			for (i = 0; i < nEachThread; i++) {
				ib = tid * nEachThread + i; 
				if (ib > acell->Nb) break;
				Nm = ib * 2; i3b = x_idx(ib);
				ia1 = acell->bound_dev[Nm]; ia2 = acell->bound_dev[Nm + 1]; i3a1 = x_idx(ia1); i3a2 = x_idx(ia2);
				c1 = acell->c_dev[ia1]; c2 = acell->c_dev[ia2];
				if (acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) continue; // they are in the same cluster, which was ignored in short-range interaction calculation
		
				ix1 = acell->ir_dev[i3a1]; iy1 = acell->ir_dev[i3a1 + 1]; iz1 = acell->ir_dev[i3a1 + 2];
				ix2 = acell->ir_dev[i3a2]; iy2 = acell->ir_dev[i3a2 + 1]; iz2 = acell->ir_dev[i3a2 + 2];
				dx = ix2 - ix1; dy = iy2 - iy1; dz = iz2 - iz1;
				if (dx > hx) doff[0] = -xd[0]; else if (dx < -hx) doff[0] = xd[0]; else doff[0] = 0;
				if (dy > hy) doff[1] = -xd[1]; else if (dy < -hy) doff[1] = xd[1]; else doff[1] = 0;
				if (dz > hz) doff[2] = -xd[2]; else if (dz < -hz) doff[2] = xd[2]; else doff[2] = 0;

				dr[0] = acell->r_dev[i3a2] - acell->r_dev[i3a1] + doff[0];
				dr[1] = acell->r_dev[i3a2 + 1] - acell->r_dev[i3a1 + 1] + doff[1];
				dr[2] = acell->r_dev[i3a2 + 2] - acell->r_dev[i3a1 + 2] + doff[2];
				R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]; if (R2 < 0.2) R2 = 0.2;
				R = GPUSQRT(R2);
				//if (R > rcut) continue;
				//u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

				U = 0;

				// calculate electrostatic interaction
				switch (cal_mode) {
				case 1: // electric field only
					//MTP_E_Q(esv, c2, dr, R, R2, E);
					EInteract_E_Q(esv0, c2, dr, R, R2, E, true);
					acell->E_b_dev[i3b] = E[0]; acell->E_b_dev[i3b + 1] = E[1]; acell->E_b_dev[i3b + 2] = E[2];
					break;
				case 2: // force and U
					//MTP_Interact_Q(esv, c1, c2, dr, R, R2, F, U);
					EInteract_Q(esv0, c1, c2, dr, R, R2, F, U, true);
					acell->Fe_b_dev[i3b] = F[0]; acell->Fe_b_dev[i3b + 1] = F[1]; acell->Fe_b_dev[i3b + 2] = F[2];
					acell->Ub_estat_dev[ib] = U;
					// virial
					acell->ve_b_dev[i3b] = F[0] * dr[0];
					acell->ve_b_dev[i3b + 1] = F[1] * dr[1];
					acell->ve_b_dev[i3b + 2] = F[2] * dr[2];
					break;
				case 11: 
					//MTP_E_MU(esv, c2, mu2, dr, R, R2, E);
					EInteract_E_MU(esv0, c2, mu2, dr, R, R2, E, true);
					acell->E_b_dev[i3b] = E[0]; acell->E_b_dev[i3b + 1] = E[1]; acell->E_b_dev[i3b + 2] = E[2];
					break;
				case 12: // force and U
					//MTP_Interact_MU(esv, c1, mu1, c2, mu2, dr, R, R2, F, U, bExcludeDipole);
					EInteract_MU(esv0, c1, mu1, c2, mu2, dr, R, R2, F, U, bExcludeDipole, true);
					acell->Fe_b_dev[i3b] = F[0]; acell->Fe_b_dev[i3b + 1] = F[1]; acell->Fe_b_dev[i3b + 2] = F[2];
					acell->Ub_estat_dev[ib] = U;
					// virial
					acell->ve_b_dev[i3b] = F[0] * dr[0];
					acell->ve_b_dev[i3b + 1] = F[1] * dr[1];
					acell->ve_b_dev[i3b + 2] = F[2] * dr[2];
					break;
				default:
					break;
				}
			}
		}
		// bound14
		Nb = acell->Nb14;
		Rmax_boundLength = 3 * 3; // bound 1-4
		if (Nb > 0) {
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;
			
			for (i = 0; i < nEachThread; i++) {
				ib = tid * nEachThread + i; 
				if (ib > acell->Nb14) break;
				Nm = ib * 2; i3b = x_idx(ib);
				ia1 = acell->bound_14_dev[Nm]; ia2 = acell->bound_14_dev[Nm + 1]; i3a1 = x_idx(ia1); i3a2 = x_idx(ia2);
				c1 = acell->c_dev[ia1]; c2 = acell->c_dev[ia2];
				if (acell->cIndx_dev[ia1] == acell->cIndx_dev[ia2]) continue; // they are in the same cluster, which was ignored in short-range interaction calculation
		
				ix1 = acell->ir_dev[i3a1]; iy1 = acell->ir_dev[i3a1 + 1]; iz1 = acell->ir_dev[i3a1 + 2];
				ix2 = acell->ir_dev[i3a2]; iy2 = acell->ir_dev[i3a2 + 1]; iz2 = acell->ir_dev[i3a2 + 2];
				dx = ix2 - ix1; dy = iy2 - iy1; dz = iz2 - iz1;
				if (dx > hx) doff[0] = -xd[0]; else if (dx < -hx) doff[0] = xd[0]; else doff[0] = 0;
				if (dy > hy) doff[1] = -xd[1]; else if (dy < -hy) doff[1] = xd[1]; else doff[1] = 0;
				if (dz > hz) doff[2] = -xd[2]; else if (dz < -hz) doff[2] = xd[2]; else doff[2] = 0;

				dr[0] = acell->r_dev[i3a2] - acell->r_dev[i3a1] + doff[0];
				dr[1] = acell->r_dev[i3a2 + 1] - acell->r_dev[i3a1 + 1] + doff[1];
				dr[2] = acell->r_dev[i3a2 + 2] - acell->r_dev[i3a1 + 2] + doff[2];
				R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]; if (R2 < 0.2) R2 = 0.2;
				R = GPUSQRT(R2);
				//if (R > rcut) continue;
				//u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

				U = 0;

				// calculate electrostatic interaction
				switch (cal_mode) {
				case 1: // electric field only
					//MTP_E_Q(esv, c2, dr, R, R2, E);
					EInteract_E_Q(esv0, c2, dr, R, R2, E, true);
					acell->E_14_dev[i3b] = E[0]; acell->E_14_dev[i3b + 1] = E[1]; acell->E_14_dev[i3b + 2] = E[2];
					break;
				case 2: // force and U
					//MTP_Interact_Q(esv, c1, c2, dr, R, R2, F, U);
					EInteract_Q(esv0, c1, c2, dr, R, R2, F, U, true);
					acell->Fe_14_dev[i3b] = F[0]; acell->Fe_14_dev[i3b + 1] = F[1]; acell->Fe_14_dev[i3b + 2] = F[2];
					acell->U14_estat_dev[ib] = U;
					// virial
					acell->ve_14_dev[i3b] = F[0] * dr[0];
					acell->ve_14_dev[i3b + 1] = F[1] * dr[1];
					acell->ve_14_dev[i3b + 2] = F[2] * dr[2];
					break;
				case 11: 
					//MTP_E_MU(esv, c2, mu2, dr, R, R2, E);
					EInteract_E_MU(esv0, c2, mu2, dr, R, R2, E, true);
					acell->E_14_dev[i3b] = E[0]; acell->E_14_dev[i3b + 1] = E[1]; acell->E_14_dev[i3b + 2] = E[2];
					break;
				case 12: // force and U
					//MTP_Interact_MU(esv, c1, mu1, c2, mu2, dr, R, R2, F, U, bExcludeDipole);
					EInteract_MU(esv0, c1, mu1, c2, mu2, dr, R, R2, F, U, bExcludeDipole, true);
					acell->Fe_14_dev[i3b] = F[0]; acell->Fe_14_dev[i3b + 1] = F[1]; acell->Fe_14_dev[i3b + 2] = F[2];
					acell->U14_estat_dev[ib] = U;
					// virial
					acell->ve_14_dev[i3b] = F[0] * dr[0];
					acell->ve_14_dev[i3b + 1] = F[1] * dr[1];
					acell->ve_14_dev[i3b + 2] = F[2] * dr[2];
					break;
				default:
					break;
				}
			}
		}
	}



	__global__ void collect_eforce_res(GPU_EAtCELL *acell) {
		// works on 1 dimensional grid and one dimensional blocks
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int Na = acell->Na;
		int nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;

		int i, j, ia1, ia2, ia, ib, iba1, iba2;
		int i3a1, i3a2;

		int Nb, Nm;

		if (acell->Na_b > 0) {
			Nb = acell->Na_b;
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;

			for (i = 0; i < nEachThread; i++) {
				ia = tid * nEachThread + i;
				if (ia >= Nb) break;
				Nm = ia * (acell->Nm_b + 2);
				ia1 = acell->b_relship_dev[Nm];
				i3a1 = x_idx(ia1);
				for (j = 0; j < acell->b_relship_dev[Nm + 1]; j++) {
					ib = acell->b_relship_dev[Nm + 2 + j];
					i3a2 = x_idx(ib);
					iba1 = acell->bound_dev[2 * ib];
					iba2 = acell->bound_dev[2 * ib + 1];
					if (iba1 == ia1) { // int this bound, ia1 is the first atom
						acell->Fe_dev[i3a1] -= acell->Fe_b_dev[i3a2];
						acell->Fe_dev[i3a1 + 1] -= acell->Fe_b_dev[i3a2 + 1];
						acell->Fe_dev[i3a1 + 2] -= acell->Fe_b_dev[i3a2 + 2];

						// collecting energy and virial term at this case only
						acell->U_estat_dev[ia1] -= acell->Ub_estat_dev[ib];
						acell->ve_dev[i3a1] -= acell->ve_b_dev[i3a2];
						acell->ve_dev[i3a1 + 1] -= acell->ve_b_dev[i3a2 + 1];
						acell->ve_dev[i3a1 + 2] -= acell->ve_b_dev[i3a2 + 2];
					}
					else if (iba2 == ia1) { // ia1 is the second atom of the bound
						acell->Fe_dev[i3a1] += acell->Fe_b_dev[i3a2];
						acell->Fe_dev[i3a1 + 1] += acell->Fe_b_dev[i3a2 + 1];
						acell->Fe_dev[i3a1 + 2] += acell->Fe_b_dev[i3a2 + 2];
					}
				}
			}
			
			__threadfence(); // its result can be seen in all other thread, in all stream
			__syncthreads();
		}

		if (acell->Na_b14 > 0) {
			Nb = acell->Na_b14;
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;

			for (i = 0; i < nEachThread; i++) {
				ia = tid * nEachThread + i;
				if (ia >= Nb) break;
				Nm = ia * (acell->Nm_b14 + 2);
				ia1 = acell->b14_relship_dev[Nm];
				i3a1 = x_idx(ia1);
				for (j = 0; j < acell->b14_relship_dev[Nm + 1]; j++) {
					ib = acell->b14_relship_dev[Nm + 2 + j];
					i3a2 = x_idx(ib);
					iba1 = acell->bound_14_dev[2 * ib];
					iba2 = acell->bound_14_dev[2 * ib + 1];
					if (iba1 == ia1) { // int this bound, ia1 is the first atom
						acell->Fe_dev[i3a1] -= acell->Fe_14_dev[i3a2] * acell->c14;
						acell->Fe_dev[i3a1 + 1] -= acell->Fe_14_dev[i3a2 + 1] * acell->c14;
						acell->Fe_dev[i3a1 + 2] -= acell->Fe_14_dev[i3a2 + 2] * acell->c14;

						// collecting energy and virial term at this case only
						acell->U_estat_dev[ia1] -= acell->U14_estat_dev[ib] * acell->c14;
						acell->ve_dev[i3a1] -= acell->ve_14_dev[i3a2] * acell->c14;
						acell->ve_dev[i3a1 + 1] -= acell->ve_14_dev[i3a2 + 1] * acell->c14;
						acell->ve_dev[i3a1 + 2] -= acell->ve_14_dev[i3a2 + 2] * acell->c14;
					}
					else if (iba2 == ia1) { // ia1 is the second atom of the bound
						acell->Fe_dev[i3a1] += acell->Fe_14_dev[i3a2] * acell->c14;
						acell->Fe_dev[i3a1 + 1] += acell->Fe_14_dev[i3a2 + 1] * acell->c14;
						acell->Fe_dev[i3a1 + 2] += acell->Fe_14_dev[i3a2 + 2] * acell->c14;
					}
				}
			}
			
			__threadfence(); // its result can be seen in all other thread, in all stream
			__syncthreads();
		}
		
		Na = acell->Na;
		nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;

		for (i = 0; i < nEachThread; i++) {
			ia1 = tid * nEachThread + i; i3a1 = x_idx(ia1);
			if (ia1 >= Na) break;
			//memcpy(acell->Fe_map + i3a1, acell->Fe_dev + i3a1, SIZE_V3);
			acell->Fe_map[i3a1] = acell->Fe_dev[i3a1];
			acell->Fe_map[i3a1 + 1] = acell->Fe_dev[i3a1 + 1];
			acell->Fe_map[i3a1 + 2] = acell->Fe_dev[i3a1 + 2];
			//memcpy(acell->ve_map + i3a1, acell->ve_dev + i3a1, SIZE_V3);
			acell->ve_map[i3a1] = acell->ve_dev[i3a1];
			acell->ve_map[i3a1 + 1] = acell->ve_dev[i3a1 + 1];
			acell->ve_map[i3a1 + 2] = acell->ve_dev[i3a1 + 2];
			acell->U_estat_map[ia1] = acell->U_estat_dev[ia1];
		}
	}




	__global__ void collect_efield_res(GPU_EAtCELL *acell) {
		// works on 1 dimensional grid and one dimensional blocks
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int Na = acell->Na;
		int nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;

		int i, j, ia1, ia2, ia, ib, iba1, iba2;
		int i3a1, i3a2;

		int Nb, Nm;

		if (acell->Na_b > 0) {
			Nb = acell->Na_b;
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;

			for (i = 0; i < nEachThread; i++) {
				ia = tid * nEachThread + i;
				if (ia >= Nb) break;
				Nm = ia * (acell->Nm_b + 2);
				ia1 = acell->b_relship_dev[Nm];
				i3a1 = x_idx(ia1);
				for (j = 0; j < acell->b_relship_dev[Nm + 1]; j++) {
					ib = acell->b_relship_dev[Nm + 2 + j];
					i3a2 = x_idx(ib);
					iba1 = acell->bound_dev[2 * ib];
					iba2 = acell->bound_dev[2 * ib + 1];
					if (iba1 == ia1) { // int this bound, ia1 is the first atom
						acell->E_dev[i3a1] -= acell->E_b_dev[i3a2];
						acell->E_dev[i3a1 + 1] -= acell->E_b_dev[i3a2 + 1];
						acell->E_dev[i3a1 + 2] -= acell->E_b_dev[i3a2 + 2];
					}
					else if (iba2 == ia1) { // ia1 is the second atom of the bound
						acell->E_dev[i3a1] += acell->E_b_dev[i3a2];
						acell->E_dev[i3a1 + 1] += acell->E_b_dev[i3a2 + 1];
						acell->E_dev[i3a1 + 2] += acell->E_b_dev[i3a2 + 2];
					}
				}
			}
			
			__threadfence(); // its result can be seen in all other thread, in all stream
			__syncthreads();
		}

		if (acell->Na_b14 > 0) {
			Nb = acell->Na_b14;
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;

			for (i = 0; i < nEachThread; i++) {
				ia = tid * nEachThread + i;
				if (ia >= Nb) break;
				Nm = ia * (acell->Nm_b14 + 2);
				ia1 = acell->b14_relship_dev[Nm];
				i3a1 = x_idx(ia1);
				for (j = 0; j < acell->b14_relship_dev[Nm + 1]; j++) {
					ib = acell->b14_relship_dev[Nm + 2 + j];
					i3a2 = x_idx(ib);
					iba1 = acell->bound_14_dev[2 * ib];
					iba2 = acell->bound_14_dev[2 * ib + 1];
					if (iba1 == ia1) { // int this bound, ia1 is the first atom
						acell->E_dev[i3a1] -= acell->E_14_dev[i3a2] * acell->c14;
						acell->E_dev[i3a1 + 1] -= acell->E_14_dev[i3a2 + 1] * acell->c14;
						acell->E_dev[i3a1 + 2] -= acell->E_14_dev[i3a2 + 2] * acell->c14;
					}
					else if (iba2 == ia1) { // ia1 is the second atom of the bound
						acell->E_dev[i3a1] += acell->E_14_dev[i3a2] * acell->c14;
						acell->E_dev[i3a1 + 1] += acell->E_14_dev[i3a2 + 1] * acell->c14;
						acell->E_dev[i3a1 + 2] += acell->E_14_dev[i3a2 + 2] * acell->c14;
					}
				}
			}
			
			__threadfence(); // its result can be seen in all other thread, in all stream
			__syncthreads();
		}

		Na = acell->Na;
		nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;

		for (i = 0; i < nEachThread; i++) {
			ia1 = tid * nEachThread + i; i3a1 = x_idx(ia1);
			if (ia1 >= Na) break;
			//memcpy(acell->E_map + i3a1, acell->E_dev + i3a1, SIZE_V3);
			acell->E_map[i3a1] = acell->E_dev[i3a1];
			acell->E_map[i3a1 + 1] = acell->E_dev[i3a1 + 1];
			acell->E_map[i3a1 + 2] = acell->E_dev[i3a1 + 2];
		}
		
		__syncthreads();
	}


#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#define printf(f, ...) ((void)(f, __VA_ARGS__), 0)
#endif

	__global__ void test(GPU_EAtCELL *acell) {
		printf("%d atoms\n", acell->Na);
	}
	
	__global__ void test1(int *f) {
		atomicAdd(f, 1);
	}
	
	__global__ void reset_esr_neighbors_kernel(GPU_EAtCELL *acell) {
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int Na = acell->Na;
		int nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;

		if (acell->srnb_dev == NULL) return;
		int srIdx;
		int *isr = NULL;
		int i, ia;
		for (i = 0; i < nEachThread; i++) {
			ia = tid * nEachThread;
			if (ia >= Na) break;
			srIdx = acell->srIndx(ia);
			isr = acell->srnb_dev + srIdx;
			isr[0] = 0;
		}

		if (tid == 0) acell->bSRInitialized = false;
	}

	__host__ void reset_esr_neighbors(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev, cudaStream_t &stream) {
		int dimBlock = 128, ngrid = 1024;
		_gpu_dstruct_::_gpu_util_::gpu_grid(acell_host->Na, dimBlock, ngrid);
		reset_esr_neighbors_kernel<<< ngrid, dimBlock, 0, stream >>>(acell_dev);
	}


	__host__ void GPU_EwaldSum_Real_3d(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev, cudaStream_t &stream) {
		int nEachBlock = 128;
		int nGRID = acell_host->Na / nEachBlock;
		if (nGRID < 1) nGRID = 1;
		if (nGRID * nEachBlock < acell_host->Na) nGRID += 1;
		
		EwaldSum_Real<<< nGRID, nEachBlock, 0, stream >>>(acell_dev);
		//test<<< nGRID, nEachBlock >>>(acell_dev);
		cudaStreamSynchronize(stream);
		
		//cudaDeviceSynchronize();
		if (acell_host->bCalE == 2) collect_eforce_res<<<nGRID, nEachBlock, 0, stream>>>(acell_dev);
		else if (acell_host->bCalE == 1)collect_efield_res<<<nGRID, nEachBlock, 0, stream>>>(acell_dev);
		//cudaDeviceSynchronize();
	};
	
	extern __global__ void EwaldSum_2D_K0(GPU_EAtCELL *cell);
	
	__host__ void GPU_EwaldSumReal_2d(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev, cudaStream_t &stream) {
		int nEachBlock = 128;
		int nGRID = acell_host->Na / nEachBlock;
		if (nGRID < 1) nGRID = 1;
		if (nGRID * nEachBlock < acell_host->Na) nGRID += 1;
		
		//test<<< nGRID, nEachBlock >>>(acell_dev);
		
		EwaldSum_Real<<< nGRID, nEachBlock, 0, stream >>>(acell_dev);
		cudaStreamSynchronize(stream);
		EwaldSum_2D_K0<<< nGRID, nEachBlock, 0, stream >>>(acell_dev);
		cudaStreamSynchronize(stream);
		
		//cudaDeviceSynchronize();
		if (acell_host->bCalE == 2) collect_eforce_res<<<nGRID, nEachBlock, 0, stream>>>(acell_dev);
		else if (acell_host->bCalE == 1)collect_efield_res<<<nGRID, nEachBlock, 0, stream>>>(acell_dev);
		//cudaDeviceSynchronize();
	};

} // namespace _gpu_md_

#endif // __GPU__
