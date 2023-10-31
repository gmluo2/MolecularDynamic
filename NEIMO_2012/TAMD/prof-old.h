#define INV_SQRT_PI  0.5641895835477563

double erf(double x);
double erfc(double x);

template <class T> class _PROFILE {
public:
	T *x, *y;
	int npts;
	T dx;
	_PROFILE() {x = NULL; y = NULL; npts = 0;};
	void release() {
		if (x != NULL) {delete[] x; x = NULL;}
		if (y != NULL) {delete[] y; y = NULL;}
		npts = 0;
		dx = 0;
	};
	~_PROFILE() {release();};
	void set(int npts, T dx, T x1 = 0) {
		release(); this->npts = npts; this->dx = dx;
		if (npts <= 0) return;
		x = new T[npts]; y = new T[npts];
		int i = 0; x[0] = x1; y[0] = 0;
		for (i = 1; i < npts; i++) {x[i] = x[i-1] + dx; y[i] = 0;}
	};
};

class DPROF : public _PROFILE<double> {
public:
	DPROF() {};
	~DPROF() {};
};

#define interpolate(f, prof, X, i) if (X <= prof.x[0]) f = prof.y[0]; else if (X >= prof.x[prof.npts-1]) f = prof.y[prof.npts-1]; else { \
	i = (int)((X - prof.x[0]) / prof.dx + 0.5); \
	if (i <= 0) i = 0; else if (i >= prof.npts - 1) i = prof.npts-2; \
	if (X == prof.x[i]) f = prof.y[i]; \
	else if (X > prof.x[i]) f = prof.y[i] + (prof.y[i+1] - prof.y[i]) / prof.dx * (X - prof.x[i]); \
	else f = prof.y[i] + (prof.y[i-1] - prof.y[i]) / prof.dx * (X - prof.x[i]);}

template <int dim, class T> class MEMORY_2D {
public:
	int N;
	T m[dim][dim];
	MEMORY_2D() {
		int i = 0, j = 0;
		N = dim;
		for (i = 0; i < dim; i++) for (j = 0; j < dim; j++) m[i][j] = 0;
	};
};

template <int dim, int dim_2d, class T> class MEMORY_3D {
public:
	int N;
	MEMORY_2D<dim_2d, T>* mm[dim];
	MEMORY_3D() {
		int i = 0;
		N = dim;
		for (i = 0; i < N; i++) mm[i] = new MEMORY_2D<dim_2d, T>;
	};
	~MEMORY_3D() {
		int i = 0;
		for (i = 0; i < N; i++) delete mm[i];
	};
};

#define mm3d(m3d, i, j, k) m3d.mm[i]->m[j][k]


// erfc function profile -- erfc_prof
// exp function profile -- exp_prof
#define ERFC(f, r, i) interpolate(f, ErfcProf, r, i)
#define EXP2(f, r, i) interpolate(f, Exp2Prof, r, i)
#define EXP(f, r, i) interpolate(f, ExpProf, r, i)
#define PERIOD(x, r1, rd, r) r = x - int((x - r1) / rd) * rd; if (r < r1) r += rd;
#define SIN(f, x, i, t) t = x; while (t < 0) t += PI2; while (t > PI2) t -= PI2; interpolate(f, SinProf, t, i)
#define COS(f, x, i, t) t = x; while (t < 0) t += PI2; while (t > PI2) t -= PI2; interpolate(f, CosProf, t, i)

#define SSIN(f, x, i, t) t = x; interpolate(f, SinProf, t, i)
#define SCOS(f, x, i, t) t = x; interpolate(f, CosProf, t, i)

extern DPROF ErfcProf, Exp2Prof, ExpProf, SinProf, CosProf;

template <class T, int N> class _MPROFILE {
public:
	T *x, *y[N];
	int npts, n;
	T dx;
	_MPROFILE() {
		x = NULL; npts = 0; n = N;
		for (int i = 0; i < N; i++) y[i] = NULL;
	};
	void release() {
		if (x != NULL) {delete[] x; x = NULL;}
		for (int i = 0; i < n; i++) {
			if (y[i] != NULL) {delete[] y[i]; y[i] = NULL;}
		}
		npts = 0;
		dx = 0;
	};
	~_MPROFILE() {release();};
	void set(int npts, T dx, T x1 = 0) {
		release();
		int i = 0, ip = 0; 
		release(); this->npts = npts; this->dx = dx;
		if (npts <= 0) return;
		x = new T[npts]; x[0] = x1;
		for (ip = 0; ip < n; ip++) {y[ip] = new T[npts]; y[ip][0] = 0;}
		for (i = 1; i < npts; i++) {
			x[i] = x[i-1] + dx; 
			for (ip = 0; ip < n; ip++) y[ip][i] = 0;
		}
	};
};

class INTERACT_PROF : public _MPROFILE<double, 2> { // y[0] -- energy, y[1] -- force
public:
	INTERACT_PROF() {};
	~INTERACT_PROF() {};
};

#define interact_interpolate(E, f, prof, r, i, t) if (r <= prof.x[0]) {E = prof.y[0][0]; f = prof.y[1][0];} \
	else if (r >= prof.x[prof.npts-1]) {E = prof.y[0][prof.npts-1]; f = prof.y[1][prof.npts-1];} \
	else {\
		i = (int)((r - prof.x[0]) / prof.dx + 0.5); \
		if (i <= 0) i = 0; else if (i >= prof.npts - 1) i = prof.npts-2; \
		if (r == prof.x[i]) {E = prof.y[0][i]; f = prof.y[1][i];} \
		else if (r > prof.x[i]) {\
			t = (r - prof.x[i]) / prof.dx; \
			E = prof.y[0][i] + (prof.y[0][i+1] - prof.y[0][i]) * t; \
			f = prof.y[0][i] + (prof.y[1][i+1] - prof.y[1][i]) * t; \
		}\
		else {\
			t = (r - prof.x[i]) / prof.dx; \
			E = prof.y[0][i] + (prof.y[0][i-1] - prof.y[0][i]) * t;\
			f = prof.y[1][i] + (prof.y[1][i-1] - prof.y[1][i]) * t;\
		}\
	}

