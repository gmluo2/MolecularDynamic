#define INV_SQRT_PI  0.5641895835477563

double erf(double x);
double erfc(double x);


template <class T> class _PROFILE_1D_ {
public:
	T *v;
	int N;
	T *t; // the first and the last point in the curve
	_PROFILE_1D_() {v = NULL; N = 0; t = NULL;};
	void release() {
		if (v != NULL) {delete[] v; v = NULL;}
		N = 0; t = NULL;
	};
	void set(int npts) {
		release(); if (npts <= 0) return;
		N = npts;
		v = new T[N];
		memset(v, 0, sizeof(T) * npts);
		t = v + N - 1;
	};
	void update_tail() {
		if (v == NULL) {t = NULL;}
		else {t = v + N - 1;}
	};
	void set_v_zero() {
		memset(v, 0, N * sizeof(T));
	};
	void normalize(T v) {
		for (int i = 0; i < N; i++) this->v[i] /= v;
	};
	~_PROFILE_1D_(){release();};
	void operator = (_PROFILE_1D_<T> &prof) {
		if (this->N != prof.N) this->set(prof.N);
		memcpy(this->v, prof.v, sizeof(T) * prof.N);
	};
	void operator *= (T v) {
		for (int i = 0; i < this->N; i++) this->v[i] *= v;
	};
	void operator /= (T v) {
		for (int i = 0; i < this->N; i++) this->v[i] /= v;
	};
};

template <class Tx, class Ty> class _PROFILE_ {
public:
	_PROFILE_1D_<Tx> x;
	_PROFILE_1D_<Ty> y;
	void release() {
		if (x.v != NULL) x.release();
		if (y.v != NULL) y.release();
	};
	void set(int npts) {
		if (x.N != npts) x.set(npts);
		if (y.N != npts) y.set(npts);
	};
	void set_y_zero() {
		memset(y.v, 0, y.N * sizeof(Ty));
	};
	~_PROFILE_(){release();};
	void save(ofstream &out) {
		if (x.v == NULL || y.v == NULL) return;
		for (int i = 0; i < x.N; i++) out<<x.v[i]<<"  "<<y.v[i]<<endl;
	};
	bool save(char *ofname) {
		ofstream out;
		out.open(ofname);
		if (!out.is_open()) return false;
		save(out); out.close(); return true;
	};
};

template <class Tx, class Ty, class Ts> class _PROFILE_EB_ {
public:
	_PROFILE_1D_<Tx> x;
	_PROFILE_1D_<Ty> y;
	_PROFILE_1D_<Ts> s; // error bar
	void release() {
		if (x.v != NULL) x.release();
		if (y.v != NULL) y.release();
		if (s.v != NULL) s.release();
	};
	void set(int npts) {
		if (x.N != npts) x.set(npts);
		if (y.N != npts) y.set(npts);
		if (s.N != npts) s.set(npts);
	};
	void set_ys_zero() {
		memset(y.v, 0, y.N * sizeof(Ty));
		memset(s.v, 0, s.N * sizeof(Ts));
	};
	~_PROFILE_EB_(){release();};
	void save(ofstream &out) {
		if (x.v == NULL || y.v == NULL || s.v == NULL) return;
		for (int i = 0; i < x.N; i++) out<<x.v[i]<<"  "<<y.v[i]<<"  "<<s.v[i]<<endl;
	};
	bool save(char *ofname) {
		ofstream out;
		out.open(ofname);
		if (!out.is_open()) return false;
		save(out); out.close(); return true;
	};
};

template <class Tx, class Ty> class _PROFILE_1D {
public:
	Tx *x;
	Ty *y;
	int N;
	_PROFILE_1D() {x = NULL; y = NULL; N = 0;};
	void release() {
		if (x != NULL) {delete[] x; x = NULL;}
		if (y != NULL) {delete[] y; y = NULL;}
		N = 0;
	};
	void set(int npts) {
		release(); if (npts <= 0) return;
		N = npts;
		x = new Tx[N]; y = new Ty[N];
		//for (int i = 0; i < N; i++) y[i] = 0;
		memset(x, 0, sizeof(Tx) * npts);
		memset(y, 0, sizeof(Ty) * npts);
	};
	void set_y_zero() {
		memset(y, 0, N * sizeof(Ty));
	};
	void normalize_y(Ty v) {
		for (int i = 0; i < N; i++) y[i] /= v;
	};
	~_PROFILE_1D(){release();};
	void save(ofstream &out) {
		if (x == NULL || y == NULL) return;
		for (int i = 0; i < N; i++) out<<x[i]<<"  "<<y[i]<<endl;
	};
	bool save(char *ofname) {
		ofstream out;
		out.open(ofname);
		if (!out.is_open()) return false;
		save(out); out.close(); return true;
	};
};

template <class T> class _PROFILE : public _PROFILE_1D<T, T> {
public:
	T dx;
	_PROFILE() {};
	void release() {
		_PROFILE_1D<T, T>::release();
		dx = 0;
	};
	~_PROFILE() {release();};
	void set(int npts, T dx, T x1 = 0) {
		_PROFILE_1D<T, T>::set(npts); this->dx = dx;
		if (npts <= 0) return;
		int i = 0; _PROFILE_1D<T, T>::x[0] = x1; 
		for (i = 1; i < npts; i++) {_PROFILE_1D<T, T>::x[i] = _PROFILE_1D<T, T>::x[i-1] + dx;}
	};
};

class DPROF : public _PROFILE<double> {
public:
	DPROF() {};
	~DPROF() {};
};

#define interpolate(f, prof, X, i) if (X <= prof.x[0]) f = prof.y[0]; else if (X >= prof.x[prof.N-1]) f = prof.y[prof.N-1]; else { \
	i = (int)((X - prof.x[0]) / prof.dx + 0.5); \
	if (i <= 0) i = 0; else if (i >= prof.N - 1) i = prof.N-2; \
	if (X == prof.x[i]) f = prof.y[i]; \
	else if (X > prof.x[i]) f = prof.y[i] + (prof.y[i+1] - prof.y[i]) / prof.dx * (X - prof.x[i]); \
	else f = prof.y[i] - (prof.y[i-1] - prof.y[i]) / prof.dx * (X - prof.x[i]);}


template <class Tx, class Ty, int N> class _MPROFILE {
public:
	Tx *x;
	Ty* y[N];
	int n, npts;
	Tx dx;
	_MPROFILE() {
		x = NULL; n = N; npts = 0;
		for (int i = 0; i < N; i++) y[i] = NULL;
	};
	void release() {
		if (x != NULL) {delete[] x; x = NULL;}
		for (int i = 0; i < N; i++) {
			if (y[i] != NULL) {delete[] y[i]; y[i] = NULL;}
		}
		npts = 0; dx = 0;
	};
	~_MPROFILE() {release();};
	void set(int npts, Tx dx, Tx x1 = 0) {
		release(); this->npts = npts; this->dx = dx;
		if (npts <= 0) return;
		x = new Tx[npts]; //memset(x, 0, sizeof(Tx) * npts);
		int i = 0, iy = 0;
		for (i = 0; i < npts; i++) x[i] = x1 + dx * i;
		for (iy = 0; iy < N; iy++) {
			y[iy] = new Ty[npts];
			memset(y[iy], 0, sizeof(Ty) * npts);
		}
		/*
		for (i = 0; i < npts; i++) {
			x[i] = x1 + i * dx; 
			for (iy = 0; iy < N; iy++) y[iy][i] = 0;
		}
		*/
	};
	void set_y_zero(int ny) {
		int i;
		if (ny > 0) {
			if (ny < N) memset(y[ny], 0, sizeof(Ty) * npts);
			return;
		}
		else {
			for (i = 0; i < N; i++) memset(y[i], 0, sizeof(Ty) * npts);
		}
	};
};

class EF_DPROF : public _MPROFILE<double, double, 2> {
public:
	double xf, xt;
	EF_DPROF() {xf = 0; xt = 0;};
	~EF_DPROF() {};
	void set(int npts, double dx, double x1 = 0) {
		_MPROFILE<double, double, 2>::set(npts, dx, x1);
		if (npts > 0) {xf = _MPROFILE<double, double, 2>::x[0]; xt = _MPROFILE<double, double, 2>::x[npts-1];}
	};
};

#define interpolate2(E, f, efprof, r, i, t) \
	if (r <= efprof.x[0]) {E = efprof.y[0][0]; f = efprof.y[1][0];}  \
	else if (r >= efprof.xt) {E = efprof.y[0][efprof.npts-1]; f = efprof.y[1][efprof.npts-1];} \
	else { \
		i = int((r - efprof.x[0]) / efprof.dx + 0.5); \
		if (i <= 0) i = 0; else if (i >= efprof.npts - 1) i = efprof.npts-2; \
		if (r == efprof.x[i]) {E = efprof.y[0][i]; f = efprof.y[1][i];} \
		else if (r > efprof.x[i]) { \
			t = (r - efprof.x[i]) / efprof.dx; \
			E = efprof.y[0][i] + (efprof.y[0][i+1] - efprof.y[0][i]) * t; \
			f = efprof.y[1][i] + (efprof.y[1][i+1] - efprof.y[1][i]) * t; \
		} \
		else { \
			t = (r - efprof.x[i]) / efprof.dx; \
			E = efprof.y[0][i] - (efprof.y[0][i-1] - efprof.y[0][i]) * t; \
			f = efprof.y[1][i] - (efprof.y[1][i-1] - efprof.y[1][i]) * t; \
		} \
	} 


#define interpolate1(f, efprof, r, i, iy) \
	if (r <= efprof.x[0]) {f = efprof.y[iy][0];}  \
	else if (r >= efprof.xt) {f = efprof.y[iy][efprof.npts-1];} \
	else { \
		i = int((r - efprof.x[0]) / efprof.dx + 0.5); \
		if (i <= 0) i = 0; else if (i >= efprof.npts - 1) i = efprof.npts-2; \
		if (r == efprof.x[i]) {f = efprof.y[iy][i];} \
		else if (r > efprof.x[i]) { \
			t = (r - efprof.x[i]) / efprof.dx; \
			f = efprof.y[iy][i] + (efprof.y[iy][i+1] - efprof.y[iy][i]) * t; \
		} \
		else { \
			t = (r - efprof.x[i]) / efprof.dx; \
			f = efprof.y[iy][i] - (efprof.y[iy][i-1] - efprof.y[iy][i]) * t; \
		} \
	} 

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
#define ERF(f, r, i) if (r >= 0) {interpolate(f, ErfProf, r, i)} else {interpolate(f, ErfProf, (-r), i) f = -f;}
#define ERFC(f, r, i) if (r >= 0) {interpolate(f, ErfcProf, r, i)} else {interpolate(f, ErfcProf, (-r), i) f = 2 - f;}
#define EXP2(f, r, i) if (r >= 0) {interpolate(f, Exp2Prof, r, i)} else {interpolate(f, Exp2Prof, (-r), i)}
#define EXP(f, r, i) interpolate(f, ExpProf, r, i)
#define PERIOD(x, r1, rd, r) r = x - int((x - r1) / rd) * rd; if (r < r1) r += rd;
#define SIN(f, x, i, t) t = x; while (t < 0) t += PI2; while (t > PI2) t -= PI2; interpolate(f, SinProf, t, i)
#define COS(f, x, i, t) t = x; while (t < 0) t += PI2; while (t > PI2) t -= PI2; interpolate(f, CosProf, t, i)

#define SSIN(f, x, i, t) t = x; interpolate(f, SinProf, t, i)
#define SCOS(f, x, i, t) t = x; interpolate(f, CosProf, t, i)

extern DPROF ErfProf, ErfcProf, Exp2Prof, ExpProf, SinProf, CosProf, IncGammaProf_2_3;
extern double ranf();

template <class Tx, class Ty> Tx xrand(_PROFILE_1D<Tx, Ty> &distrbt) {
	int n = 0;
	if (distrbt.N <= 0) return 0;
	else if (distrbt.N == 1) return distrbt.x[0];
	
	while (1) {
		n = int(distrbt.N * ranf());
		if (n >= distrbt.N) continue;
		if (ranf() < distrbt.y[n]) return distrbt.x[n];
	}
};
