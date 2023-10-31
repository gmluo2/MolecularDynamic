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
#define EXP2(fexp2, kR) fexp2 = exp(-kR * kR);
#define Ewald_f0(f, kR) f = erfc(kR);
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
