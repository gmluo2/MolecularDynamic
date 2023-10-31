namespace _gpu_md_ {
	using namespace _gpu_dstruct_;
	namespace _gpu_ewald_ {
		
#ifndef _KAPPA
#define _KAPPA     3.2
#endif

#ifndef INV_SQRT_PI
#define INV_SQRT_PI  0.5641895835477563
#endif

#define T1(r, u, i, dis3) for (i = 0; i < u.n; i++) u.v[i] = -r.v[i] / dis3;

#define T2(r, mt, i, j, dis, dis2, dis5, temp) for (i = 0; i < mt.n; i++) {for (j = 0; j < mt.n; j++) {\
		temp = 3 * r.v[i] * r.v[j]; if (i == j) temp -= dis2; mt.m[i][j] = temp / dis5; \
	}}

#define T3(r, cm, i, j, k, dis, dis2, dis7, temp) \
	for (i = 0; i < cm.n; i++) {for (j = 0; j < cm.n; j++) {for(k = 0; k < cm.n; k++) { \
		cm.m[i][j][k] = 15 * r.v[i] * r.v[j] * r.v[k]; temp = 0.0; \
		if (i == j) temp += r.v[k]; if (i == k) temp += r.v[j]; if (j == k) temp += r.v[i]; \
		cm.m[i][j][k] -= 3 * temp * dis2; \
		cm.m[i][j][k] /= -(dis7); \
	}}}

#define T4(r, qm, i, j, k, l, dis, dis2, dis4, dis9, temp) \
	for (i = 0; i < qm.n; i++) {for (j = 0; j < qm.n; j++) {for(k = 0; k < qm.n; k++) {for (l = 0; l < qm.n; l++) { \
		qm.m[i][j][k][l] = 105 * r.v[i] * r.v[j] * r.v[k] * r.v[l]; temp = 0; \
		if (i == j) temp += r.v[k] * r.v[l]; if (i == k) temp += r.v[j] * r.v[l]; \
		if (i == l) temp += r.v[j] * r.v[k]; if (j == k) temp += r.v[i] * r.v[l]; \
		if (j == l) temp += r.v[i] * r.v[k]; if (k == l) temp += r.v[i] * r.v[j]; \
		if (temp != 0) {temp *= 15 * dis2; qm.m[i][j][k][l] -= temp;} \
		temp = 0; \
		if (i == j && k == l) temp += 1; \
		if (i == k && j == l) temp += 1; \
		if (i == l && j == k) temp += 1; \
		if (temp != 0) qm.m[i][j][k][l] += 3 * temp * dis4; \
		qm.m[i][j][k][l] /= dis9; \
	}}}}

#define MTP_T1(u, r, R3) u.v[0] = -r.v[0] / R3; u.v[1] = -r.v[1] / R3; u.v[2] = -r.v[2] / R3;
#define MTP_T2(m2, r, R2, R5, i, j) for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {\
	m2.m[i][j] = 3 * r.v[i] * r.v[j]; if (i == j) m2.m[i][j] -= R2; m2.m[i][j] /= R5; \
	}}

#define MTP_T3(m3, r, R2, R7, i, j, k, t) t = 3 * R2; for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {\
	m3.m[i][j][k] = 15 * r.v[i] * r.v[j] * r.v[k]; \
	if (i == j) m3.m[i][j][k] -= t * r.v[k]; if (i == k) m3.m[i][j][k] -= t * r.v[j]; \
	if (j == k) m3.m[i][j][k] -= t * r.v[i]; m3.m[i][j][k] /= -R7; \
	}}}

#define MTP_T4(m4, r, R2, R4, R9, i, j, k, l, t1, t2) t1 = 15 * R2; t2 = 3 * R4; \
	for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {for (l = 0; l < 3; l++) {\
	m4.m[i][j][k][l] = 105 * r.v[i] * r.v[j] * r.v[k] * r.v[l]; \
	if (i == j) m4.m[i][j][k][l] -= t1 * r.v[k] * r.v[l]; \
	if (i == k) m4.m[i][j][k][l] -= t1 * r.v[j] * r.v[l]; \
	if (i == l) m4.m[i][j][k][l] -= t1 * r.v[j] * r.v[k]; \
	if (j == k) m4.m[i][j][k][l] -= t1 * r.v[i] * r.v[l]; \
	if (j == l) m4.m[i][j][k][l] -= t1 * r.v[i] * r.v[k]; \
	if (k == l) m4.m[i][j][k][l] -= t1 * r.v[i] * r.v[j]; \
	if (i == j && k == l) m4.m[i][j][k][l] += t2; \
	if (i == k && j == l) m4.m[i][j][k][l] += t2; \
	if (i == l && j == k) m4.m[i][j][k][l] += t2; \
	m4.m[i][j][k][l] /= R9; \
	}}}}

#define MTP_D2(m2, r, i, j) for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {\
	m2.m[i][j] = r.v[i] * r.v[j]; \
	}}
#define EXP2(fexp2, kR) fexp2 = GPUEXP(-kR * kR);
#define Ewald_f0(f, kR) f = GPUERFC(kR);
#define Ewald_f1(f, kR, f0, fexp2) f = f0 + 2 * kR * fexp2 * INV_SQRT_PI;
#define Ewald_f2(f, kappa3, R2, fexp2) f = 4 * kappa3 * INV_SQRT_PI * fexp2 / R2;
#define Ewald_f3(f, R3, f2) f = R3 * f2;
#define Ewald_f4(f, kappa3, kR2, R4, fexp2) f = 8 * kappa3 * INV_SQRT_PI * (kR2 + 1) / R4 * fexp2;
#define Ewald_f5(f, R, R3, f2, f4) f = R3 * f4 - 3 * R * f2;
#define Ewald_f6(f, R2, kR2, R4, f2, f4) f = 2 * (kR2 + 1) / R2 * f4 + 4 / R4 * f2;
#define Ewald_f(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4, f5, f6) \
	EXP2(fexp2, kR) \
	Ewald_f0(f0, kR) Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2) \
	Ewald_f3(f3, R3, f2) Ewald_f4(f4, kappa3, kR2, R4, fexp2) \
	Ewald_f5(f5, R, R3, f2, f4) Ewald_f6(f6, R2, kR2, R4, f2, f4)

#define Ewald_f04(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4) \
	EXP2(fexp2, kR) \
	Ewald_f0(f0, kR) Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2) \
	Ewald_f3(f3, R3, f2) Ewald_f4(f4, kappa3, kR2, R4, fexp2)

#define Ewald_f02(kappa3, R, R2, kR, kR2, i, fexp2, f0, f1, f2) \
	EXP2(fexp2, kR) \
	Ewald_f0(f0, kR) Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2) 


#define _Ewald_f_simplify(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4, f5, f6) \
	Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2) \
	Ewald_f3(f3, R3, f2) Ewald_f4(f4, kappa3, kR2, R4, fexp2) \
	Ewald_f5(f5, R, R3, f2, f4) Ewald_f6(f6, R2, kR2, R4, f2, f4)

#define _Ewald_f04_simplify(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4) \
	Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2) \
	Ewald_f3(f3, R3, f2) Ewald_f4(f4, kappa3, kR2, R4, fexp2)

#define _Ewald_f02_simplify(kappa3, R, R2, kR, kR2, i, fexp2, f0, f1, f2) \
	Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2)


#define EWALD_T1(u, f1, uT1) u.v[0] = f1 * uT1.v[0]; u.v[1] = f1 * uT1.v[1]; u.v[2] = f1 * uT1.v[2];

#define EWALD_T2(m2, f1, f2, mT2, mD2, i, j) for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {\
	m2.m[i][j] = f1 * mT2.m[i][j] + f2 * mD2.m[i][j]; \
	}}

#define EWALD_T3(m3, f1, f2, f3, f4, r, mT2, mD2, mT3, i, j, k) for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {\
	m3.m[i][j][k] = f1 * mT3.m[i][j][k] - f3 * r.v[k] * mT2.m[i][j] - f4 * r.v[k] * mD2.m[i][j]; \
	if (j == k) m3.m[i][j][k] += f2 * r.v[i]; if (i == k) m3.m[i][j][k] += f2 * r.v[j]; \
	}}}

#define EWALD_T4(m4, f1, f2, f3, f4, f5, f6, r, mT2, mD2, mT3, mT4, i, j, k, l) \
	for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) { for (l = 0; l < 3; l++) { \
	m4.m[i][j][k][l] = f1 * mT4.m[i][j][k][l] - f3 * (r.v[l] * mT3.m[i][j][k] + r.v[k] * mT3.m[i][j][l]) + f5 * mD2.m[k][l] * mT2.m[i][j] + f6 * mD2.m[k][l] * mD2.m[i][j]; \
	if (i == l && j == k) m4.m[i][j][k][l] += f2; \
	if (i == k && j == l) m4.m[i][j][k][l] += f2; \
	if (k == l) m4.m[i][j][k][l] -= f3 * mT2.m[i][j]; \
	if (k == l) m4.m[i][j][k][l] -= f4 * mD2.m[i][j]; \
	if (j == l) m4.m[i][j][k][l] -= f4 * mD2.m[i][k]; \
	if (j == k) m4.m[i][j][k][l] -= f4 * mD2.m[i][l]; \
	if (i == l) m4.m[i][j][k][l] -= f4 * mD2.m[j][k]; \
	if (i == k) m4.m[i][j][k][l] -= f4 * mD2.m[j][l]; \
	}}}}

#define GPU_VECTOR3 VECTOR3<DOUBLE>
		// the variable which is used for Ewald-Sum for real term calculation
		// important: this class hs to be the exactly same as the one defined in EwaldSumVar.h
		// because the memory will be copied from host to GPU
		class EwaldSumRealVars0 { // for charge - charge interaction
		public:
			GPU_VECTOR3 mT1, mET1;
			//double* dr;
			DOUBLE R, R2, R3;
			DOUBLE kR;
			DOUBLE fexp2, f0, f1;

			bool bEwald;
			DOUBLE kappa, kappa3;
			short iSurface; // 0 -- cubic use 4 * PI / 3V,  1 -- slab use 4PI/V
			short iSurfaceBoundary; // cubic or 2d cell 0/1
			short iExIndDipole; // excluding induced dipole -- 0/1
			DOUBLE inv_3V; // 4 * PI / (3 * V)
			DOUBLE inv_V; // 4 * PI / V

			__device__ void init(bool bEwald, int iSurface, int iSurfaceBoundary, int iExIndDipole, DOUBLE kappa, DOUBLE V) {
				this->bEwald = bEwald;
				this->iSurface = (short)iSurface;
				this->iSurfaceBoundary = (short)iSurfaceBoundary;
				this->iExIndDipole = (short)iExIndDipole;
				this->kappa = kappa; this->kappa3 = kappa * kappa * kappa;
				if (bEwald || V > 1) {
					this->inv_V = 4 * PI / V;
					this->inv_3V = this->inv_V / 3;
				}
			};

			__device__ EwaldSumRealVars0() {bEwald = true; inv_3V = 0; iSurface = 0; inv_V = 0; iExIndDipole = 0;};
			__device__ void init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2);
			__device__ void init_vars(DOUBLE *dr, DOUBLE R, DOUBLE R2) {GPU_VECTOR3 dR(dr); init_vars(dR, R, R2);};

			__device__ void init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2, DOUBLE f_kR, DOUBLE f_exp2);
			__device__ void init_vars(DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE f_kR, DOUBLE f_exp2) {GPU_VECTOR3 dR(dr); init_vars(dR, R, R2, f_kR, f_exp2);};

			__device__ void init_EwaldSum(float V, float rcut_Ewd) {
				inv_V = 4 * PI / V;
				inv_3V = inv_V / 3;
				kappa = _KAPPA / rcut_Ewd;
				kappa3 = kappa * kappa * kappa;
			};
			__device__ void cp_EwaldSum_vars(EwaldSumRealVars0 &var) {
				bEwald = var.bEwald;
				iSurfaceBoundary = var.iSurfaceBoundary;
				iSurface = var.iSurface;
				iExIndDipole = var.iExIndDipole;
				inv_3V = var.inv_3V;
				inv_V = var.inv_V;
				kappa = var.kappa;
				kappa3 = var.kappa3;
			};
		};

		class EwaldSumRealVars1 : public EwaldSumRealVars0 {
		public:
			// more variable for dipole interaction
			MATRIX<DOUBLE, 3> mT2, mET2, mD2;
			CMATRIX<DOUBLE, 3> mT3, mET3;
			//QMATRIX<DOUBLE, 3> mT4, mET4;

			DOUBLE R4, R5, R7;//, R9;
			DOUBLE kR2;
	
			DOUBLE f2, f3, f4;//, f5, f6;

			__device__ EwaldSumRealVars1() {};
			__device__ void init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2);
			__device__ void init_vars(DOUBLE *dr, DOUBLE R, DOUBLE R2) {GPU_VECTOR3 dR(dr); init_vars(dR, R, R2);};

			__device__ void init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2, DOUBLE f_kR, DOUBLE f_exp2);
			__device__ void init_vars(DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE f_kR, DOUBLE f_exp2) {GPU_VECTOR3 dR(dr); init_vars(dR, R, R2, f_kR, f_exp2);};

			__device__ void init_vars_E(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2); // for E calculation only
			__device__ void init_vars_E(DOUBLE *dr, DOUBLE R, DOUBLE R2) {GPU_VECTOR3 dR(dr); init_vars_E(dR, R, R2);}; // for E calculation only

			__device__ void init_vars_E(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2, DOUBLE f_kR, DOUBLE f_exp2); // for E calculation only
			__device__ void init_vars_E(DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE f_kR, DOUBLE f_exp2) {GPU_VECTOR3 dR(dr); init_vars_E(dR, R, R2, f_kR, f_exp2);}; // for E calculation only
		};

		class EwaldSumRealVars2 : public EwaldSumRealVars1 {
		public:
			// more variable for dipole-quardupole interaction
			QMATRIX<DOUBLE, 3> mT4, mET4;

			DOUBLE R9, f5, f6;

			__device__ EwaldSumRealVars2() {};
			__device__ void init_vars(GPU_VECTOR3 &dr, DOUBLE R, DOUBLE R2);
		};

		__device__ void EwaldReal_Q(EwaldSumRealVars1 &esv, DOUBLE c1, DOUBLE c2, DOUBLE *dr, DOUBLE R, DOUBLE R2, /*double *E, */DOUBLE *F, DOUBLE &U);
		__device__ void EInteract_Q(_gpu_ewald_::EwaldSumRealVars1 &esv, DOUBLE c1, DOUBLE c2, DOUBLE *dr, DOUBLE R, DOUBLE R2, /*double *E, */DOUBLE *F, DOUBLE &U, bool init_esv);
		__device__ void EwaldReal_MU(EwaldSumRealVars1 &esv, DOUBLE c1, DOUBLE *mu1, DOUBLE c2, DOUBLE *mu2, DOUBLE *dr, DOUBLE R, DOUBLE R2, /*double *E,*/ DOUBLE *F, DOUBLE &U, char bExcludeDipole);
		__device__ void EInteract_MU(EwaldSumRealVars1 &esv, DOUBLE c1, DOUBLE *mu1, DOUBLE c2, DOUBLE *mu2, DOUBLE *dr, DOUBLE R, DOUBLE R2, /*double *E,*/ DOUBLE *F, DOUBLE &U, char bExcludeDipole, bool init_esv);
		__device__ void EwaldReal_E_Q(EwaldSumRealVars1 &esv, DOUBLE c2, DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE *E);
		__device__ void EInteract_E_Q(EwaldSumRealVars1 &esv, DOUBLE c2, DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE *E, bool init_esv);
		__device__ void EwaldReal_E_MU(EwaldSumRealVars1 &esv, DOUBLE c2, DOUBLE *mu2, DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE *E);
		__device__ void EInteract_E_MU(EwaldSumRealVars1 &esv, DOUBLE c2, DOUBLE *mu2, DOUBLE *dr, DOUBLE R, DOUBLE R2, DOUBLE *E, bool init_esv);






		//*******************************************************************************************************************
		//  ********   functions for 2d Ewald-Sum, with k=0, which actually belongs to the receipt part of 2d-Ewald-Sum *****
		//*******************************************************************************************************************
		class EwaldSumRecpK0Var {
		public:
			bool bEfieldOnly;
			bool bVirial;
			short iExIndDipole; // excluding induced dipole -- 0/1
			double kappa, A;

			// explicit parameters to calculate the reciprocal term with k = 0
			double PI_over_A, PI2_over_A;  // PI/A, 2PI / A
			double W; // dz * erf(kappa * dz) + 1 / (kappa * Sqrt(PI)) * exp(-(kappa * dz)^2)
			double dz, dz_abs, az, az2, az_abs, ferf, fexp;

			double inverse_kappa_sqrt_pi;  // 1 / (kappa * sqrt(PI))
			double KAPPA2_over_SQRT_PI, kappa2;  // 2 * kappa/sqrt(PI), kappa^2

			__device__ EwaldSumRecpK0Var() {bEfieldOnly = false; bVirial = true; iExIndDipole = 0;};

			__device__ void set(double kappa, double A) {
				this->kappa = kappa; this->A = A; 
				PI_over_A = PI / A; PI2_over_A = PI_over_A * 2;
				inverse_kappa_sqrt_pi = 1 / kappa * INV_SQRT_PI;

				kappa2 = kappa * kappa;
				KAPPA2_over_SQRT_PI = 2 * kappa * INV_SQRT_PI;
			};
			__device__ void init(bool bEfieldOnly, int iExIndDipole, double kappa, double A) {
				this->bEfieldOnly = bEfieldOnly;
				this->iExIndDipole = (short)iExIndDipole;
				this->set(kappa, A);
				// we always calculate virial
			};
			__device__ void init_var(double dz) {
				int i;
				this->dz = dz; dz_abs = FABS(dz);
				az = kappa * dz; az_abs = FABS(az); az2 = az * az;
				ferf = erf(az); 
				W = dz * ferf;
				if (az2 < 16) {
					fexp = exp(-az2);//EXP2(fexp, az_abs, i)
					W += inverse_kappa_sqrt_pi * fexp;
				}
				else fexp = 0;
			};
			__device__ void operator = (EwaldSumRecpK0Var &var1) {
				this->bEfieldOnly = var1.bEfieldOnly;
				this->bVirial = var1.bVirial;
				this->iExIndDipole = var1.iExIndDipole;
				this->set(var1.kappa, var1.A);
				// init_var() is a dynamic part of the variables, 
				// the related variables will not be copied by operator =
			};
		};

		__device__ void EwaldSumRecpKzero(EwaldSumRecpK0Var &var, DOUBLE q1, DOUBLE q2, DOUBLE &Ez, DOUBLE &Fz, DOUBLE &U);
		__device__ void EwaldSumRecpKzero(EwaldSumRecpK0Var &var, DOUBLE q1, DOUBLE mu1z, DOUBLE q2, DOUBLE mu2z, DOUBLE &Ez, DOUBLE &Fz, DOUBLE &U);

	
	} // namespace _gpu_ewald_

	//__global__ void EwaldSum_Real(GPU_EAtCELL *acell);
	//__global__ void collect_efield_res(GPU_EAtCELL *acell);
	//__global__ void collect_eforce_res(GPU_EAtCELL *acell);
} // namespace _gpu_md_
