// macro defined for interaction calculation

// we can disable exp profile by define 
// #define EXP2(f, r) f = exp(-r * r)
#define Ewald_f0(f, kR, i) ERFC(f, kR, i)
#define Ewald_f1(f, kR, f0, fexp2) f = f0 + 2 * kR * fexp2 * INV_SQRT_PI;
#define Ewald_f2(f, kappa3, R2, fexp2) f = 4 * kappa3 * INV_SQRT_PI * fexp2 / R2;
#define Ewald_f3(f, R3, f2) f = R3 * f2;
#define Ewald_f4(f, kappa3, kR2, R4, fexp2) f = 8 * kappa3 * INV_SQRT_PI * (kR2 + 1) / R4 * fexp2;
#define Ewald_f5(f, R, R3, f2, f4) f = R3 * f4 - 3 * R * f2;
#define Ewald_f6(f, R2, kR2, R4, f2, f4) f = 2 * (kR2 + 1) / R2 * f4 + 4 / R4 * f2;
#define Ewald_f(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4, f5, f6) \
	EXP2(fexp2, kR, i) \
	Ewald_f0(f0, kR, i) Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2) \
	Ewald_f3(f3, R3, f2) Ewald_f4(f4, kappa3, kR2, R4, fexp2) \
	Ewald_f5(f5, R, R3, f2, f4) Ewald_f6(f6, R2, kR2, R4, f2, f4)

#define Ewald_f04(kappa3, R, R2, R3, R4, kR, kR2, i, fexp2, f0, f1, f2, f3, f4) \
	EXP2(fexp2, kR, i) \
	Ewald_f0(f0, kR, i) Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2) \
	Ewald_f3(f3, R3, f2) Ewald_f4(f4, kappa3, kR2, R4, fexp2)

#define Ewald_f02(kappa3, R, R2, kR, kR2, i, fexp2, f0, f1, f2) \
	EXP2(fexp2, kR, i) \
	Ewald_f0(f0, kR, i) Ewald_f1(f1, kR, f0, fexp2) Ewald_f2(f2, kappa3, R2, fexp2) 

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

#define EWALD_EFIELD_Q_REAL(E, q, mEwaldT1) E.v[0] = -q * mEwaldT1.v[0]; E.v[1] = -q * mEwaldT1.v[1]; E.v[2] = -q * mEwaldT1.v[2];

#define EWALD_EFIELD_DIPOLE_REAL(E, mu, mEwaldT2, i, j) for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {\
	E.v[i] = mEwaldT2.m[i][j] * mu.v[j]; \
	}

#define EWALD_EFIELD_QUAD_REAL(E, Q, mEwaldT3, i, j, k) for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {\
	E.v[i] = -mEwaldT3.m[i][j][k] * Q.m[j][k] / 3; \
	}}}

#define EWALD_EFIELD_MULTIPOL_REAL(E, q, mu, Q, mEwaldT1, mEwaldT2, mEwaldT3, i, j, k) for (i = 0; i < 3; i++) {\
	E.v[i] = -q * mEwaldT1.v[i]; \
	for (j = 0; j < 3; j++) {E.v[i] += mEwaldT2.m[i][j] * mu.v[j]; \
		for (k = 0; k < 3; k++) E.v[i] -= mEwaldT3.m[i][j][k] * Q.m[j][k] / 3; \
	}}

#define EWALD_FORCE1_REAL(F, q1, mu1, q2, mu2, mEwaldT1, mEwaldT2, mEwaldT3, i, j, k) \
	for (i = 0; i < 3; i++) { F.v[i] = -q1 * q2 * mEwaldT1.v[i]; \
		for (j = 0; j < 3; j++) {F.v[i] += mEwaldT2.m[i][j] * (q1 * mu2.v[j] - q2 * mu1.v[j]); \
			for (k = 0; k < 3; k++) {\
				F.v[i] += mEwaldT3.m[i][j][k] * mu1.v[j] * mu2.v[k]; \
			} \
	}}

#define EWALD_FORCE2_REAL(F, q1, mu1, q2, mu2, Q2, mEwaldT1, mEwaldT2, mEwaldT3, mEwaldT4, i, j, k, l) \
	for (i = 0; i < 3; i++) { F.v[i] = -q1 * q2 * mEwaldT1.v[i]; \
		for (j = 0; j < 3; j++) {F.v[i] += mEwaldT2.m[i][j] * (q1 * mu2.v[j] - q2 * mu1.v[j]); \
			for (k = 0; k < 3; k++) {\
				F.v[i] += mEwaldT3.m[i][j][k] * (mu1.v[j] * mu2.v[k] - q1 * Q2.m[j][k] / 3); \
				for (l = 0; l < 3; l++) {\
					F.v[i] -= mEwaldT4.m[i][j][k][l] * mu1.v[j] * Q2.m[k][l] / 3; \
				} \
	}}}




#define EFIELD_MULTIPOL_REAL(E, q, mu, Q, mT1, mT2, mT3, i, j, k) for (i = 0; i < 3; i++) {\
	E.v[i] = -q * mT1.v[i]; \
	for (j = 0; j < 3; j++) {E.v[i] += mT2.m[i][j] * mu.v[j]; \
		for (k = 0; k < 3; k++) E.v[i] -= mT3.m[i][j][k] * Q.m[j][k] / 3; \
	}}

#define FORCE1_REAL(F, q1, mu1, q2, mu2, mT1, mT2, mT3, i, j, k) \
	for (i = 0; i < 3; i++) { F.v[i] = -q1 * q2 * mT1.v[i]; \
		for (j = 0; j < 3; j++) { \
			F.v[i] += mT2.m[i][j] * (q1 * mu2.v[j] - q2 * mu1.v[j]); \
			for (k = 0; k < 3; k++) {\
				F.v[i] += mT3.m[i][j][k] * mu1.v[j] * mu2.v[k]; \
			} \
		}\
	}

#define FORCE2_REAL(F, q1, mu1, q2, mu2, Q2, mT1, mT2, mT3, mT4, i, j, k, l) \
	for (i = 0; i < 3; i++) { F.v[i] = -q1 * q2 * mT1.v[i]; \
		for (j = 0; j < 3; j++) {F.v[i] += mT2.m[i][j] * (q1 * mu2.v[j] - q2 * mu1.v[j]); \
			for (k = 0; k < 3; k++) {\
				F.v[i] += mT3.m[i][j][k] * (mu1.v[j] * mu2.v[k] - q1 * Q2.m[j][k] / 3); \
				for (l = 0; l < 3; l++) {\
					F.v[i] -= mT4.m[i][j][k][l] * mu1.v[j] * Q2.m[k][l] / 3; \
				} \
	}}}


#define FORCE3_REAL(F, q1, mu1, Q1, q2, mu2, Q2, mT1, mT2, mT3, mT4, i, j, k, l) \
	for (i = 0; i < 3; i++) { F.v[i] = -q1 * q2 * mT1.v[i]; \
		for (j = 0; j < 3; j++) {F.v[i] += mT2.m[i][j] * (q1 * mu2.v[j] - q2 * mu1.v[j]); \
			for (k = 0; k < 3; k++) {\
				F.v[i] += mT3.m[i][j][k] * (mu1.v[j] * mu2.v[k] - q1 * Q2.m[j][k] / 3 - q2 * Q1.m[j][k] / 3); \
				for (l = 0; l < 3; l++) {\
					F.v[i] -= mT4.m[i][j][k][l] * (mu1.v[j] * Q2.m[k][l] - Q1.m[k][l] * mu2.v[j]) / 3 ; \
				} \
	}}}

void test();

#if _DISABLE_ == 0
void CMMEwaldSumInteraction(CMM_CELL3D<BASIC_CLUSTER> &cell, double &V, bool &OVERLAP, bool hard_sphere, int n1Cell, int n2Cell);

void EwaldSumCMMInteraction_SingleMM(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, double &V, bool &OVERLAP, bool hard_sphere, int n1, int n2);

void EwaldSumCMMInteraction_SM(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE *sm, int nSM, double &V, bool &OVERLAP); // V in kT
#endif // _DISABLE_
