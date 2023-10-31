#define _CHECK_TIME_STAMP_  0

#define UNKNOWN_BOUND       0
#define SINGLE_BOUND        1
#define DOUBLE_BOUND        2
#define TRIPLE_BOUND        3

 // method form Ewald-Sum
#define _EWALD_SUM   1   // type needs to be defined, 0 : 1

#define _MultiPole_EWALD_SUM  0
#define _SPME_EWALD_SUM  1
// end of definition of Ewald-Sum

// Length of Nose-Hoover Chain
#define MNHC  3
// end of definition of length of Nose-Hoover Chain

#define MD_NVT   0
#define MD_NPT   1

// MD constraint
#define MD_NO_CONSTRAINT      0
#define MD_PMF_CONSTRAINT     1
#define MD_SLAB               2
#define MD_IMAGE_SLAB         3


#define _LOCAL_SPEED_HOOVER        0
#define _GLOBAL_SPEED_HOOVER       1

#define MIN_ROT   1e-10

#define NCLUSTERS_PARALLEL        100
#define MAX_THREADS               4

// the macro-molecule which can be considered as big-macromolecules 
#define NCLUSTERS_BIG_MM          20

// check cluster in the cell of CMM for every some MD steps
#define N_MDSTEP_CHECK_CMM        5

// sometime we want to make a judgement of whether run lots of simple operations in parallel
#define N4PARALLEL             80

// virial calculation
#define CLUSTER_VIRIAL      0
#define MOLECULE_VIRIAL     1

// coarse-grained MD, use FENE or given neighbor-interaction from files?
// 0 -- FENE
// 1 -- given neighbor-interaction from files
#define _CG_NEIGHBOR_INTERACT_    1

#define RELEASE_ARRAY(p) if (p != NULL) {delete[] p; p = NULL;}

#define same_sign(a, b) (((a > 0 && b > 0) || (a < 0 && b < 0)) ? true : false)
#define sign(a) (a > 0 ? 1 : -1)

#define FABS(x) (x > 0 ? x : -(x))

#define PERIOD_RANGE(x, xl, xr, period) while (x > xr) x -= period; while (x < xl) x += period;

#define PERIOD_CHECK(x, xl, xr, period, delta) delta = 0; \
	if (x > xr) {while (x > xr) {x -= period; delta -= period;}} \
	else if (x < xl) {while (x < xl) {x += period; delta += period;}}

//#define PERIOD_CHECK(x, xl, xr, period, delta) if (x <= xr && x >= xl) delta = 0; else {\
//	delta = int((x - xl) / period); if (delta < 0) delta -= 1; delta = (delta == 0 ? 0 : -delta * period); \
//	delta = (delta + x < xl ? delta + period : delta);\
//	}

// memory for char, int, float, double, check them with sizeof(char), sizeof(int), sizeof(float), sizeof(double)
#define SIZE_CHAR    1
#define SIZE_INT     4
#define SIZE_FLOAT   4
#define SIZE_DOUBLE  8
// 3 * double
#define SIZE_V3      24   
// 6 * double
#define SIZE_V6      48   
// 3 x 3 * double
#define SIZE_M3      72  
// 6 x 6 * double
#define SIZE_M6      288  

// matrix element
#define me(m, ni, nj, i, j) m[ni + i][nj + j]

#define V2V(u1, u2, i, n) for (i = 0; i < n; i++) u2.v[i] = u1.v[i];

#define INV_V3(u) u.v[0] = -u.v[0]; u.v[1] = -u.v[1]; u.v[2] = -u.v[2];
#define INV_V6(u) u.v[0] = -u.v[0]; u.v[1] = -u.v[1]; u.v[2] = -u.v[2]; u.v[3] = -u.v[3]; u.v[4] = -u.v[4]; u.v[5] = -u.v[5];

//#define V32V3(u1, u2) u2.v[0] = u1.v[0]; u2.v[1] = u1.v[1]; u2.v[2] = u1.v[2];
//#define V62V6(u1, u2) u2.v[0] = u1.v[0]; u2.v[1] = u1.v[1]; u2.v[2] = u1.v[2]; u2.v[3] = u1.v[3]; u2.v[4] = u1.v[4]; u2.v[5] = u1.v[5];

#define V32V3(u1, u2) memcpy(u2.v, u1.v, SIZE_V3);
#define DV32V3(u1, u2) memcpy(u2, u1, SIZE_V3);
#define V62V6(u1, u2) memcpy(u2.v, u1.v, SIZE_V6);
#define DV62V6(u1, u2) memcpy(u2, u1, SIZE_V6);

#define V3PC(u1, a, u2) u2.v[0] = u1.v[0] * a; u2.v[1] = u1.v[1] * a; u2.v[2] = u1.v[2] * a;

#define M2M(m1, m2, i, j, n) for (i = 0; i < n; i++) for (j = 0; j < n; j++) m2.m[i][j] = m1.m[i][j];

//#define Vzero(u, i) for (i = 0; i < u.n; i++) u.v[i] = 0;
//#define V3zero(u) u.v[0] = 0; u.v[1] = 0; u.v[2] = 0;
//#define V6zero(u) u.v[0] = 0; u.v[1] = 0; u.v[2] = 0; u.v[3] = 0; u.v[4] = 0; u.v[5] = 0;

#define Vzero(u, i) memset(u.v, 0, u.n * sizeof(u.v));
#define V3zero(u) memset(u.v, 0, SIZE_V3);
#define V6zero(u) memset(u.v, 0, SIZE_V6);
#define DV3zero(u) memset(u, 0, SIZE_V3);

#define Vabs(u1, u2, i) for (i = 0; i < u1.n; i++) u2.v[i] = (u1.v[i] > 0 ? u1.v[i] : -(u1.v[i]));

#define VplusV(v1, v2, res, i, n) for (i = 0; i < n; i++) res.v[i] = v1.v[i] + v2.v[i];

#define VminusV(v1, v2, res, i, n) for (i = 0; i < n; i++) res.v[i] = v1.v[i] - v2.v[i];

#define V3plusV3(v1, v2, res) res.v[0] = v1.v[0] + v2.v[0]; res.v[1] = v1.v[1] + v2.v[1]; res.v[2] = v1.v[2] + v2.v[2];
#define V3minusV3(v1, v2, res) res.v[0] = v1.v[0] - v2.v[0]; res.v[1] = v1.v[1] - v2.v[1]; res.v[2] = v1.v[2] - v2.v[2];

#define V6plusV6(v1, v2, res) res.v[0] = v1.v[0] + v2.v[0]; res.v[1] = v1.v[1] + v2.v[1]; res.v[2] = v1.v[2] + v2.v[2]; \
	res.v[3] = v1.v[3] + v2.v[3]; res.v[4] = v1.v[4] + v2.v[4]; res.v[5] = v1.v[5] + v2.v[5];
#define V6minusV6(v1, v2, res) res.v[0] = v1.v[0] - v2.v[0]; res.v[1] = v1.v[1] - v2.v[1]; res.v[2] = v1.v[2] - v2.v[2]; \
	res.v[3] = v1.v[3] - v2.v[3]; res.v[4] = v1.v[4] - v2.v[4]; res.v[5] = v1.v[5] - v2.v[5];

#define V6plus(v1, v2, res) res[0] = v1[0] + v2[0]; res[1] = v1[1] + v2[1]; res[2] = v1[2] + v2[2]; \
	res[3] = v1[3] + v2[3]; res[4] = v1[4] + v2[4]; res[5] = v1[5] + v2[5];
#define V6minus(v1, v2, res) res[0] = v1[0] - v2[0]; res[1] = v1[1] - v2[1]; res[2] = v1[2] - v2[2]; \
	res[3] = v1[3] - v2[3]; res[4] = v1[4] - v2[4]; res[5] = v1[5] - v2[5];

#define V3plus(v1, v2, res) res[0] = v1[0] + v2[0]; res[1] = v1[1] + v2[1]; res[2] = v1[2] + v2[2];
#define V3minus(v1, v2, res) res[0] = v1[0] - v2[0]; res[1] = v1[1] - v2[1]; res[2] = v1[2] - v2[2];


#define EachAxisInRange(u, i, vmin, vmax, status) status = true; \
	for (i = 0; i < u.n; i++) {if (u.v[i] > vmax || u.v[i] < vmin) {status = false; break;}}

#define ABS2(u, i, dis) dis = 0; for (i = 0; i < u.n; i++) dis += u.v[i] * u.v[i];
#define V3ABS2(u, dis) dis = u.v[0] * u.v[0] + u.v[1] * u.v[1] + u.v[2] * u.v[2];
#define DV3ABS2(u, dis) dis = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];

#define ABS(u, i, dis) dis = 0; for (i = 0; i < u.n; i++) dis += u.v[i] * u.v[i]; dis = sqrt(dis);

#define DIV(u1, c, u2, i) for (i = 0; i < u1.n; i++) u2.v[i] = u1.v[i] / c;
#define DIV3(u1, c, u2) u2.v[0] = u1.v[0] / c; u2.v[1] = u1.v[1] / c; u2.v[2] = u1.v[2] / c;
#define DIV6(u1, c, u2) u2.v[0] = u1.v[0] / c; u2.v[1] = u1.v[1] / c; u2.v[2] = u1.v[2] / c; u2.v[3] = u1.v[3] / c; u2.v[4] = u1.v[4] / c; u2.v[5] = u1.v[5] / c;

#define TIME(u1, c, u2, i) for (i = 0; i < u1.n; i++) u2.v[i] = u1.v[i] * c;
#define TIME3(u1, c, u2) u2.v[0] = u1.v[0] * c; u2.v[1] = u1.v[1] * c; u2.v[2] = u1.v[2] * c;
#define TIME6(u1, c, u2) u2.v[0] = u1.v[0] * c; u2.v[1] = u1.v[1] * c; u2.v[2] = u1.v[2] * c; u2.v[3] = u1.v[3] * c; u2.v[4] = u1.v[4] * c; u2.v[5] = u1.v[5] * c;

#define MAX(u, i, max) for (i = 0; i < u.n; i++) {if (i == 0) max = u.v[i]; else max = (FABS(u.v[i]) > FABS(max) ? u.v[i] : max);}

/*
#define Mzero(m1, i, j) for (i = 0; i < m1.n; i++) for (j = 0; j < m1.n; j++) m1.m[i][j] = 0;

#define M3zero(m1) \
	m1.m[0][0] = 0; m1.m[0][1] = 0; m1.m[0][2] = 0;\
	m1.m[1][0] = 0; m1.m[1][1] = 0; m1.m[1][2] = 0;\
	m1.m[2][0] = 0; m1.m[2][1] = 0; m1.m[2][2] = 0;\

#define M6zero(m1) \
	m1.m[0][0] = 0; m1.m[0][1] = 0; m1.m[0][2] = 0; m1.m[0][3] = 0; m1.m[0][4] = 0; m1.m[0][5] = 0; \
	m1.m[1][0] = 0; m1.m[1][1] = 0; m1.m[1][2] = 0; m1.m[1][3] = 0; m1.m[1][4] = 0; m1.m[1][5] = 0; \
	m1.m[2][0] = 0; m1.m[2][1] = 0; m1.m[2][2] = 0; m1.m[2][3] = 0; m1.m[2][4] = 0; m1.m[2][5] = 0; \
	m1.m[3][0] = 0; m1.m[3][1] = 0; m1.m[3][2] = 0; m1.m[3][3] = 0; m1.m[3][4] = 0; m1.m[3][5] = 0; \
	m1.m[4][0] = 0; m1.m[4][1] = 0; m1.m[4][2] = 0; m1.m[4][3] = 0; m1.m[4][4] = 0; m1.m[4][5] = 0; \
	m1.m[5][0] = 0; m1.m[5][1] = 0; m1.m[5][2] = 0; m1.m[5][3] = 0; m1.m[5][4] = 0; m1.m[5][5] = 0;

#define IMatrix(m1, i, j) for (i = 0; i < m1.n; i++) for (j = 0; j < m1.n; j++) {\
	if (i == j) m1.m[i][j] = 1.0; else m1.m[i][j] = 0.0;}
*/
#define M3zero(m1) memset(m1.m[0], 0, SIZE_M3);
#define M6zero(m1) memset(m1.m[0], 0, SIZE_M6);
#define I3Matrix(m1) memset(m1.m[0], 0, SIZE_M3); m1.m[0][0] = 1; m1.m[1][1] = 1; m1.m[2][2] = 1; 
#define I6Matrix(m1) memset(m1.m[0], 0, SIZE_M6); m1.m[0][0] = 1; m1.m[1][1] = 1; m1.m[2][2] = 1; m1.m[3][3] = 1; m1.m[4][4] = 1; m1.m[5][5] = 1;

#define V3PV3(u1, u2, res)  res.v[0] = u1.v[1] * u2.v[2] - u1.v[2] * u2.v[1]; \
		res.v[1] = u1.v[2] * u2.v[0] - u1.v[0] * u2.v[2]; \
		res.v[2] = u1.v[0] * u2.v[1] - u1.v[1] * u2.v[0];

#define DV3PV3(u1, u2, res)  res[0] = u1[1] * u2[2] - u1[2] * u2[1]; \
		res[1] = u1[2] * u2[0] - u1[0] * u2[2]; \
		res[2] = u1[0] * u2[1] - u1[1] * u2[0];

#define VECT(u1, u2, r, i) for (i = 0; i < r.n; i++) r.v[i] = u2.v[i] - u1.v[i];

#define VECT3(u1, u2, r) r.v[0] = u2.v[0] - u1.v[0]; r.v[1] = u2.v[1] - u1.v[1]; r.v[2] = u2.v[2] - u1.v[2];
#define DVECT3(u1, u2, r) r[0] = u2[0] - u1[0]; r[1] = u2[1] - u1[1]; r[2] = u2[2] - u1[2];

#define u2Matrix(u, m, ni, nj) \
	me(m, ni, nj, 0, 0) = 0; me(m, ni, nj, 0, 1) = -u.v[2]; me(m, ni, nj, 0, 2) = u.v[1]; \
	me(m, ni, nj, 1, 0) = u.v[2]; me(m, ni, nj, 1, 1) = 0; me(m, ni, nj, 1, 2) = -u.v[0]; \
	me(m, ni, nj, 2, 0) = -u.v[1]; me(m, ni, nj, 2, 1) = u.v[0]; me(m, ni, nj, 2, 2) = 0; 

#define u2Matrix3(u, m) \
	m[0][0] = 0; m[0][1] = -u.v[2]; m[0][2] = u.v[1]; \
	m[1][0] = u.v[2]; m[1][1] = 0; m[1][2] = -u.v[0]; \
	m[2][0] = -u.v[1]; m[2][1] = u.v[0]; m[2][2] = 0; 
/*
#define V2wv(w, u, V) \
	w.v[0] = V.v[0]; w.v[1] = V.v[1]; w.v[2] = V.v[2]; \
	u.v[0] = V.v[3]; u.v[1] = V.v[4]; u.v[2] = V.v[5];

#define wv2V(w, u, V) \
	V.v[0] = w.v[0]; V.v[1] = w.v[1]; V.v[2] = w.v[2]; \
	V.v[3] = u.v[0]; V.v[4] = u.v[1]; V.v[5] = u.v[2];
*/
#define V2wv(w, u, V) memcpy(w.v, V.v, SIZE_V3); memcpy(u.v, V.v + 3, SIZE_V3);

#define wv2V(w, u, V) memcpy(V.v, w.v, SIZE_V3); memcpy(V.v + 3, u.v, SIZE_V3);

#define uv2M(u1, u2, M, i, j, n) for (i = 0; i < n; i++) { \
	for (j = 0; j < n; j++) M.m[i][j] = u1.v[i] * u2.v[j];}

#define uv2TAO6(u1, u2, M) M.m[0][0] = 1 - u1.v[0] * u2.v[0]; M.m[0][1] = -u1.v[0] * u2.v[1]; M.m[0][2] = -u1.v[0] * u2.v[2]; \
	M.m[0][3] = -u1.v[0] * u2.v[3]; M.m[0][4] = -u1.v[0] * u2.v[4]; M.m[0][5] = -u1.v[0] * u2.v[5]; \
	M.m[1][0] = -u1.v[1] * u2.v[0]; M.m[1][1] = 1 - u1.v[1] * u2.v[1]; M.m[1][2] = -u1.v[1] * u2.v[2]; \
	M.m[1][3] = -u1.v[1] * u2.v[3]; M.m[1][4] = -u1.v[1] * u2.v[4]; M.m[1][5] = -u1.v[1] * u2.v[5]; \
	M.m[2][0] = -u1.v[2] * u2.v[0]; M.m[2][1] = -u1.v[2] * u2.v[1]; M.m[2][2] = 1 - u1.v[2] * u2.v[2]; \
	M.m[2][3] = -u1.v[2] * u2.v[3]; M.m[2][4] = -u1.v[2] * u2.v[4]; M.m[2][5] = -u1.v[2] * u2.v[5]; \
	M.m[3][0] = -u1.v[3] * u2.v[0]; M.m[3][1] = -u1.v[3] * u2.v[1]; M.m[3][2] = -u1.v[3] * u2.v[2]; \
	M.m[3][3] = 1 - u1.v[3] * u2.v[3]; M.m[3][4] = -u1.v[3] * u2.v[4]; M.m[3][5] = -u1.v[3] * u2.v[5]; \
	M.m[4][0] = -u1.v[4] * u2.v[0]; M.m[4][1] = -u1.v[4] * u2.v[1]; M.m[4][2] = -u1.v[4] * u2.v[2]; \
	M.m[4][3] = -u1.v[4] * u2.v[3]; M.m[4][4] = 1 - u1.v[4] * u2.v[4]; M.m[4][5] = -u1.v[4] * u2.v[5]; \
	M.m[5][0] = -u1.v[5] * u2.v[0]; M.m[5][1] = -u1.v[5] * u2.v[1]; M.m[5][2] = -u1.v[5] * u2.v[2]; \
	M.m[5][3] = -u1.v[5] * u2.v[3]; M.m[5][4] = -u1.v[5] * u2.v[4]; M.m[5][5] = 1 - u1.v[5] * u2.v[5]; 

#define MpM(m1, m2, res, i, j, k, n) for (i = 0; i < n; i++) {for (j = 0; j < n; j++) {\
	res.m[i][j] = 0; for (k = 0; k < n; k++) res.m[i][j] += m1.m[i][k] * m2.m[k][j];}}

#define MpV(m1, u, res, i, j, n) for (i = 0; i < n; i++) { res.v[i] = 0; \
	for (j = 0; j < n; j++) res.v[i] += m1.m[i][j] * u.v[j];}
	
#define VTpM(u, m1, res, i, j, n) for (i = 0; i < n; i++) {\
	res.v[i] = 0; for (j = 0; j < n; j++) res.v[i] += u.v[j] * m1.m[j][i];}

#define MpV3(m1, u, res) \
	res.v[0] = m1.m[0][0] * u.v[0] + m1.m[0][1] * u.v[1] + m1.m[0][2] * u.v[2]; \
	res.v[1] = m1.m[1][0] * u.v[0] + m1.m[1][1] * u.v[1] + m1.m[1][2] * u.v[2]; \
	res.v[2] = m1.m[2][0] * u.v[0] + m1.m[2][1] * u.v[1] + m1.m[2][2] * u.v[2];

#define MpV6(m1, u, res) \
	res.v[0] = m1.m[0][0] * u.v[0] + m1.m[0][1] * u.v[1] + m1.m[0][2] * u.v[2] + m1.m[0][3] * u.v[3] + m1.m[0][4] * u.v[4] + m1.m[0][5] * u.v[5]; \
	res.v[1] = m1.m[1][0] * u.v[0] + m1.m[1][1] * u.v[1] + m1.m[1][2] * u.v[2] + m1.m[1][3] * u.v[3] + m1.m[1][4] * u.v[4] + m1.m[1][5] * u.v[5]; \
	res.v[2] = m1.m[2][0] * u.v[0] + m1.m[2][1] * u.v[1] + m1.m[2][2] * u.v[2] + m1.m[2][3] * u.v[3] + m1.m[2][4] * u.v[4] + m1.m[2][5] * u.v[5]; \
	res.v[3] = m1.m[3][0] * u.v[0] + m1.m[3][1] * u.v[1] + m1.m[3][2] * u.v[2] + m1.m[3][3] * u.v[3] + m1.m[3][4] * u.v[4] + m1.m[3][5] * u.v[5]; \
	res.v[4] = m1.m[4][0] * u.v[0] + m1.m[4][1] * u.v[1] + m1.m[4][2] * u.v[2] + m1.m[4][3] * u.v[3] + m1.m[4][4] * u.v[4] + m1.m[4][5] * u.v[5]; \
	res.v[5] = m1.m[5][0] * u.v[0] + m1.m[5][1] * u.v[1] + m1.m[5][2] * u.v[2] + m1.m[5][3] * u.v[3] + m1.m[5][4] * u.v[4] + m1.m[5][5] * u.v[5]; 

#define scale_uv(u1, u2, res, i, n) res = 0; \
	for (i = 0; i < n; i++) res += u1.v[i] * u2.v[i];

#define scale_uv3(u1, u2, res) res = u1.v[0] * u2.v[0] + u1.v[1] * u2.v[1] + u1.v[2] * u2.v[2];

#define scale_uv6(u1, u2, res) res = u1.v[0] * u2.v[0] + u1.v[1] * u2.v[1] + u1.v[2] * u2.v[2] \
	+ u1.v[3] * u2.v[3] + u1.v[4] * u2.v[4] + u1.v[5] * u2.v[5];

#define SWAP(a, b, temp) temp = a; a = b; b = temp;

/*
#define M6T2M6(m1, m2, i, j, t) for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {\
	me(m2.m, 0, 0, i, j) = me(m1.m, 0, 0, i, j); me(m2.m, 3, 3, i, j) = me(m1.m, 3, 3, i, j); \
	me(m2.m, 0, 3, i, j) = me(m1.m, 3, 0, i, j); me(m2.m, 3, 0, i, j) = me(m1.m, 0, 3, i, j);}
*/

#define M6T2M6(m1, m2, i, j) for (i = 0; i < 6; i++) for (j = 0; j < 6; j++) {\
	m2.m[i][j] = m1.m[j][i];}

#define Matrix2Vect(m, ni, nj, v) v[0] = me(m, ni, nj, 2, 1); v[1] = me(m, ni, nj, 0, 2);\
	v[2] = me(m, ni, nj, 1, 0);

#define PHAI(phai, r, i, j) for (i = 0; i < phai.n; i++) for (j = 0; j < phai.n; j++) phai.m[i][j] = (i == j ? 1 : 0); \
	u2Matrix(r, phai.m, 0, 3)

#define MplusM(m1, m2, res, i, j) for (i = 0; i < res.n; i++) for (j = 0; j < res.n; j++) res.m[i][j] = m1.m[i][j] + m2.m[i][j];

#define M6plusM6(m1, m2, res) res.m[0][0] = m1.m[0][0] + m2.m[0][0]; res.m[0][1] = m1.m[0][1] + m2.m[0][1]; res.m[0][2] = m1.m[0][2] + m2.m[0][2]; \
	res.m[0][3] = m1.m[0][3] + m2.m[0][3]; res.m[0][4] = m1.m[0][4] + m2.m[0][4]; res.m[0][5] = m1.m[0][5] + m2.m[0][5]; \
	res.m[1][0] = m1.m[1][0] + m2.m[1][0]; res.m[1][1] = m1.m[1][1] + m2.m[1][1]; res.m[1][2] = m1.m[1][2] + m2.m[1][2]; \
	res.m[1][3] = m1.m[1][3] + m2.m[1][3]; res.m[1][4] = m1.m[1][4] + m2.m[1][4]; res.m[1][5] = m1.m[1][5] + m2.m[1][5]; \
	res.m[2][0] = m1.m[2][0] + m2.m[2][0]; res.m[2][1] = m1.m[2][1] + m2.m[2][1]; res.m[2][2] = m1.m[2][2] + m2.m[2][2]; \
	res.m[2][3] = m1.m[2][3] + m2.m[2][3]; res.m[2][4] = m1.m[2][4] + m2.m[2][4]; res.m[2][5] = m1.m[2][5] + m2.m[2][5]; \
	res.m[3][0] = m1.m[3][0] + m2.m[3][0]; res.m[3][1] = m1.m[3][1] + m2.m[3][1]; res.m[3][2] = m1.m[3][2] + m2.m[3][2]; \
	res.m[3][3] = m1.m[3][3] + m2.m[3][3]; res.m[3][4] = m1.m[3][4] + m2.m[3][4]; res.m[3][5] = m1.m[3][5] + m2.m[3][5]; \
	res.m[4][0] = m1.m[4][0] + m2.m[4][0]; res.m[4][1] = m1.m[4][1] + m2.m[4][1]; res.m[4][2] = m1.m[4][2] + m2.m[4][2]; \
	res.m[4][3] = m1.m[4][3] + m2.m[4][3]; res.m[4][4] = m1.m[4][4] + m2.m[4][4]; res.m[4][5] = m1.m[4][5] + m2.m[4][5]; \
	res.m[5][0] = m1.m[5][0] + m2.m[5][0]; res.m[5][1] = m1.m[5][1] + m2.m[5][1]; res.m[5][2] = m1.m[5][2] + m2.m[5][2]; \
	res.m[5][3] = m1.m[5][3] + m2.m[5][3]; res.m[5][4] = m1.m[5][4] + m2.m[5][4]; res.m[5][5] = m1.m[5][5] + m2.m[5][5]; 

#define MminusM(m1, m2, res, i, j) for (i = 0; i < res.n; i++) for (j = 0; j < res.n; j++) res.m[i][j] = m1.m[i][j] - m2.m[i][j];

#define scale_MM(m1, m2, res, i, j) res = 0; for (i = 0; i < m1.n; i++) for (j = 0; j < m1.n; j++) res += m1.m[i][j] * m2.m[i][j];

#define scale_VMV(u1, m1, u2, res, i, j, n) res = 0; for (i = 0; i < n; i++) for (j = 0; j < n; j++) res += u1.v[i] * m1.m[i][j] * u2.v[j];

#define scale_CMCM(cm1, cm2, res, i, j, k) res = 0; for (i = 0; i < cm1.n; i++) for (j = 0; j < cm1.n; j++) for (k = 0; k < cm1.n; k++) res += cm1.m[i][j][k] * cm2.m[i][j][k];

#define scale_QMQM(qm1, qm2, res, i, j, k, l) res = 0;  for (i = 0; i < qm1.n; i++) for (j = 0; j < qm1.n; j++) for (k = 0; k < qm1.n; k++) for (l = 0; l < qm1.n; l++) res += qm1.m[i][j][k][l] * cm2.m[i][j][k][l];

#define Vect_MV(m1, u, res, i, j) for (i = 0; i < res.n; i++) {res.v[i] = 0; \
		for (j = 0; j < u.n; j++) res.v[i] += m1.m[i][j] * u.v[j]; \
	}

#define Vect_CMM(cm, m1, res, i, j, k) for (i = 0; i < res.n; i++) {res.v[i] = 0; \
		for (j = 0; j < m1.n; j++) for (k = 0; k < m1.n; k++) res.v[i] += cm.m[i][j][k] * m1.m[j][k]; \
	}

#define Vect_QMCM(qm, cm, res, i, j, k, l) for (i = 0; i < res.n; i++) {res.v[i] = 0; \
		for (j = 0; j < cm.n; j++) for (k = 0; k < cm.n; k++) for (l = 0; l < cm.n; l++) res.v[i] += qm.m[i][j][k][l] * cm.m[j][k][l]; \
	}


/*
#define v3cp(dest, source) dest[0] = source[0]; dest[1] = source[1]; dest[2] = source[2];
#define v6cp(dest, source) dest[0] = source[0]; dest[1] = source[1]; dest[2] = source[2]; \
	dest[3] = source[3]; dest[4] = source[4]; dest[5] = source[5];
*/

#define v3cp(dest, source) memcpy(dest, source, SIZE_V3); 
#define v6cp(dest, source) memcpy(dest, source, SIZE_V6); 


#define dmatrix2array(dv, dm, i, j, k) k = 0; for (i = 0; i < dm.n; i++) for (j = 0; j < dm.n; j++) {\
	dv[k] = dm.m[i][j]; k++;}

#define array2dmatrix(dm, dv, i, j, k) k = 0; for (i = 0; i < dm.n; i++) for (j = 0; j < dm.n; j++) {\
	dm.m[i][j] = dv[k]; k++;}

#define cmatrix2array(dv, cm, i, j, k, l) l = 0; for (i = 0; i < cm.n; i++) for (j = 0; j < cm.n; j++) for (k = 0; k < cm.n; k++) {\
	dv[l] = cm.m[i][j][k]; l++;}

#define array2cmatrix(cm, dv, i, j, k, l) l = 0; for (i = 0; i < cm.n; i++) for (j = 0; j < cm.n; j++) for (k = 0; k < cm.n; k++) {\
	cm.m[i][j][k] = dv[l]; l++;}


#define VplusV3(v1, v2, v3) v3[0] = v1[0] + v2[0]; v3[1] = v1[1] + v2[1]; v3[2] = v1[2] + v2[2];

#define TRUNCATE_RANGE(x, x1, x2) if (x < x1) x = x1; else if (x > x2) x = x2;
