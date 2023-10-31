class BSplineFuncs {
private:
	int nprof, npts; // number of profile, npts for each profile
public:
	double dx;
	double xmax;
	int N; // index of the B-spline, from 2 on
	ARRAY<double> x;
	FLEX_ASYM_MATRIX<double> bsp; // B-spline functions, from 2 on
	FLEX_ASYM_MATRIX<double> dbsp; //1st order differential of B-spline functions, from 2 on

	BSplineFuncs() {dx = 0.01; xmax = 0; N = 0;};
	void release() {x.release(); bsp.release(); dbsp.release();};
	~BSplineFuncs() {release();};
	void init(int indx);
	int get_pt(double x) {
		if (x < 0) return -1;
		else if (x > xmax) return npts + 1;
		return int(x / dx + 0.1);
	};
	double bsp_interpolate(int indx, int npt, double x) {
		if (npt < 0 || npt >= npts) return 0;
		int n = indx - 2;
		if (n < 0 || indx > N) return 0;
		double y;
		double *f = bsp.m[n].m;
		if (x == this->x.m[npt]) return f[npt];
		if (npt == npts - 1) y = f[npt] - (f[npt-1] - f[npt]) / dx * (x - this->x.m[npt]);
		else y = f[npt] + (f[npt+1] - f[npt]) / dx * (x - this->x.m[npt]);
		return y;
	};
	double dbsp_interpolate(int indx, int npt, double x) {
		if (npt < 0 || npt >= npts) return 0;
		int n = indx - 2;
		if (n < 0 || indx > N) return 0;
		double y;
		double *f = dbsp.m[n].m;
		if (x == this->x.m[npt]) return f[npt];
		if (npt == npts - 1) y = f[npt] - (f[npt-1] - f[npt]) / dx * (x - this->x.m[npt]);
		else y = f[npt] + (f[npt+1] - f[npt]) / dx * (x - this->x.m[npt]);
		return y;
	};
};

// nth B-Spline function in bsp, indx-th derivation
double diff(BSplineFuncs &bsps, int indx, int nth, double x);

template <int nD> class BSplineFunc : public _MPROFILE<double, double, nD + 1> {
public:
	short indx; // B-Spline function index, start from 2
	short ndiff; // nth order of derivation is required, start from 1st ...
	double xmax;
	int mpt; // the point with x = 1
	BSplineFunc() {indx = 0; xmax = 0; ndiff = nD;};
	~BSplineFunc() {};
	int get_pt(double xv) {
		if (xv < 0) return -1;
		else if (xv > xmax) return _MPROFILE<double, double, nD + 1>::npts + 1;
		return int(xv / _MPROFILE<double, double, nD + 1>::dx + 0.1);
	};
	double bsp_interpolate(int npt, double xv) {
		if (npt < 0 || npt >= _MPROFILE<double, double, nD + 1>::npts) return 0;
		double yv;
		double *f = _MPROFILE<double, double, nD + 1>::y[0];
		if (xv == _MPROFILE<double, double, nD + 1>::x[npt]) return f[npt];
		if (npt == _MPROFILE<double, double, nD + 1>::npts - 1) yv = f[npt] - (f[npt-1] - f[npt]) / _MPROFILE<double, double, nD + 1>::dx * (xv - _MPROFILE<double, double, nD + 1>::x[npt]);
		else yv = f[npt] + (f[npt+1] - f[npt]) / _MPROFILE<double, double, nD + 1>::dx * (xv - _MPROFILE<double, double, nD + 1>::x[npt]);
		return yv;
	};
	double dbsp_interpolate(int npt, double xv, int nth_diff = 1) {
		if (npt < 0 || npt >= _MPROFILE<double, double, nD + 1>::npts || nth_diff > ndiff) return 0;
		double yv;
		double *f = _MPROFILE<double, double, nD + 1>::y[nth_diff];
		if (xv == _MPROFILE<double, double, nD + 1>::x[npt]) return f[npt];
		if (npt == _MPROFILE<double, double, nD + 1>::npts - 1) yv = f[npt] - (f[npt-1] - f[npt]) / _MPROFILE<double, double, nD + 1>::dx * (xv - _MPROFILE<double, double, nD + 1>::x[npt]);
		else yv = f[npt] + (f[npt+1] - f[npt]) / _MPROFILE<double, double, nD + 1>::dx * (xv - _MPROFILE<double, double, nD + 1>::x[npt]);
		return yv;
	};
};

template <int nD> void transfer_func(BSplineFuncs &bsps, int indx, BSplineFunc<nD> &bsp) {
//void transfer_func(BSplineFuncs &bsps, int indx, BSplineFunc &bsp) {
	if (indx < 2 || indx > bsps.N) return;
	bsp.xmax = indx; bsp.dx = bsps.dx;
	int npts = int(double(indx) / bsp.dx + 0.1) + 1;
	int n = indx - 2;
	int ipt, idiff;
	bsp.set(npts, bsps.dx, 0);
	bsp.indx = indx;
	for (ipt = 0; ipt < bsp.npts; ipt++) {
		bsp.x[ipt] = bsps.x.m[ipt];
		bsp.y[0][ipt] = bsps.bsp.m[n].m[ipt];
		bsp.y[1][ipt] = bsps.dbsp.m[n].m[ipt];
	}
	bsp.mpt = int(1.0 / bsp.dx + 0.1);
	if (nD <= 1) return;
	for (idiff = 2; idiff <= nD; idiff++) {
		for (ipt = 0; ipt < bsp.npts; ipt++) {
			bsp.y[idiff][ipt] = diff(bsps, indx, idiff, bsp.x[ipt]);
		}
	}
};

template <int nD> void init_BSplineFunc(BSplineFunc<nD> &bsp, int indx) {
	if (bsp.indx == indx && bsp.y[0] != NULL) return;
	BSplineFuncs bspf;
	bspf.init(indx);
	transfer_func<nD>(bspf, indx, bsp);
	/*
	int i, j;
	ofstream out;
	out.open("B-Spline.dat");
	for (i = 0; i < bsp.npts; i++) {
		out<<bsp.x[i];
		for (j = 0; j <= nD; j++) out<<"  "<<bsp.y[j][i];
		out<<endl;
	}	
	out.close();
	*/
};

class AtomBSP {
public:
	ARRAY<double> buff;
	ARRAY<double> Mn0, Mn1, Mn2;
	ARRAY<double> dMn0, dMn1, dMn2;
	ARRAY<double> d2Mn0, d2Mn1, d2Mn2;

	/*
	void init(int bsp_indx) {
		Mn0.set_array(bsp_indx);
		Mn1.set_array(bsp_indx);
		Mn2.set_array(bsp_indx);
		dMn0.set_array(bsp_indx);
		dMn1.set_array(bsp_indx);
		dMn2.set_array(bsp_indx);

		d2Mn0.set_array(bsp_indx);
		d2Mn1.set_array(bsp_indx);
		d2Mn2.set_array(bsp_indx);
	};
	*/
	void init(int bsp_indx) {
		int ipos = 0;
		buff.set_array(bsp_indx * 9);
		Mn0.n = bsp_indx; Mn0.m = buff.m; ipos += bsp_indx;
		Mn1.n = bsp_indx; Mn1.m = buff.m + ipos; ipos += bsp_indx;
		Mn2.n = bsp_indx; Mn2.m = buff.m + ipos; ipos += bsp_indx;

		dMn0.n = bsp_indx; dMn0.m = buff.m + ipos; ipos += bsp_indx;
		dMn1.n = bsp_indx; dMn1.m = buff.m + ipos; ipos += bsp_indx;
		dMn2.n = bsp_indx; dMn2.m = buff.m + ipos; ipos += bsp_indx;

		d2Mn0.n = bsp_indx; d2Mn0.m = buff.m + ipos; ipos += bsp_indx;
		d2Mn1.n = bsp_indx; d2Mn1.m = buff.m + ipos; ipos += bsp_indx;
		d2Mn2.n = bsp_indx; d2Mn2.m = buff.m + ipos; ipos += bsp_indx;
	};
	void release() {
		Mn0.m = NULL; Mn0.n = 0; Mn1.m = NULL; Mn1.n = 0; Mn2.m = NULL; Mn2.n = 0;
		dMn0.m = NULL; dMn0.n = 0; dMn1.m = NULL; dMn1.n = 0; dMn2.m = NULL; dMn2.n = 0;
		d2Mn0.m = NULL; d2Mn0.n = 0; d2Mn1.m = NULL; d2Mn1.n = 0; d2Mn2.m = NULL; d2Mn2.n = 0;
		buff.release();
	};
	AtomBSP() {};
	~AtomBSP() {release();};
};

class __AtomBSP_BUFF {
public:
	ARRAY<AtomBSP> absp;
	bool bReady_absp;
	void init_absp(int natoms, int bsp_indx) {
		int i;
		absp.set_array(natoms);
		for (i = 0; i < natoms; i++) absp.m[i].init(bsp_indx);
	};
};

inline void cp_absp_buff(__AtomBSP_BUFF *s, __AtomBSP_BUFF *d) {
	int i;
	int n = s->absp.m[0].buff.n * SIZE_DOUBLE;
	//int n = s->absp.m[0].Mn0.n * SIZE_DOUBLE;
	for (i = 0; i < s->absp.n; i++) {
		memcpy(d->absp.m[i].buff.m, s->absp.m[i].buff.m, n);
		/*
		memcpy(d->absp.m[i].Mn0.m, s->absp.m[i].Mn0.m, n);
		memcpy(d->absp.m[i].Mn1.m, s->absp.m[i].Mn1.m, n);
		memcpy(d->absp.m[i].Mn2.m, s->absp.m[i].Mn2.m, n);

		memcpy(d->absp.m[i].dMn0.m, s->absp.m[i].dMn0.m, n);
		memcpy(d->absp.m[i].dMn1.m, s->absp.m[i].dMn1.m, n);
		memcpy(d->absp.m[i].dMn2.m, s->absp.m[i].dMn2.m, n);

		memcpy(d->absp.m[i].d2Mn0.m, s->absp.m[i].d2Mn0.m, n);
		memcpy(d->absp.m[i].d2Mn1.m, s->absp.m[i].d2Mn1.m, n);
		memcpy(d->absp.m[i].d2Mn2.m, s->absp.m[i].d2Mn2.m, n);
		*/
	}
}


namespace _spme_ {
using namespace _evatom_;
using namespace _EwaldSum_real_;
using namespace _EwaldSum_recp_;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

template <class atom> class __VSPME :  public _VE_CELL<atom>, public __AtomBSP_BUFF {
public:
	float md[3]; // full mesh-size
};

template <class atom> class _VSPME : public __VSPME<atom> {
public:
	ARRAY<double> ftd;
	ARRAY<fftw_complex> ftc;
	fftw_plan ft_plan;
	fftw_plan ift_plan;

public:
	int N; // B-spine 
	BSplineFunc<2> *bsp; // has 2nd order derivation
	int Nx, Ny, Nz, Nyz, Nv, nhx, nhy, nhz;
	//float V;
	ARRAY<double> k[3];
	ARRAY<double> b[3]; // x, y, z
	FLEX_ASYM_CMATRIX<double> C;

	// Stress Tensor calculation
	// This is a cubic cell. So, the off-diagonal elements of the tensor are 0
	// We calculate the diagonal elements, xx, yy and zz
	bool bST; // calculate stress tensor ? 
	bool bST_diagonal; // ist == 0 -- xx, yy and zz, ist == 1 -- trace of st tensor only
	FLEX_ASYM_CMATRIX<double> Cxx, Cyy, Czz, Ctr;
	double STxx, STyy, STzz, STtr;

	_VSPME() {
		N = 0; bsp = NULL;
		Nx = 0; Ny = 0; Nz = 0;  
		ft_plan = 0; ift_plan = 0;
		nhx = 0; nhy = 0; nhz = 0;

		bST = false; bST_diagonal = false; // trace of stress tensor only
		STxx = 0; STyy = 0; STzz = 0; STtr = 0;
	};
	void init_fft_plans() {
		if (ft_plan != 0 || ift_plan != 0) {release_plans();}
		if (MAX_THREADS > 1 && Nv > 4096) { // 16^3
			fftw_init_threads(); // enable multi-thread
			fftw_plan_with_nthreads(MAX_THREADS);
		}
		ft_plan = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, ftd.m, ftc.m, FFTW_ESTIMATE);
		ift_plan = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, ftc.m, ftd.m, FFTW_ESTIMATE);
	};
	void release_plans() {
		if (ft_plan != 0) {fftw_destroy_plan(ft_plan); ft_plan = 0;}
		if (ift_plan != 0) {fftw_destroy_plan(ift_plan); ift_plan = 0;}
	};
	~_VSPME() {release_plans(); bsp = NULL; N = 0;};

	void init(int bsp_indx, int nx, int ny, int nz) {
		N = bsp_indx;
		k[0].SetArray(nx); k[1].SetArray(ny); k[2].SetArray(nz);
		b[0].SetArray(nx); b[1].SetArray(ny); b[2].SetArray(nz);
		this->Nx = nx; this->Ny = ny; this->Nz = nz;
		nhx = Nx / 2; nhy = Ny / 2; nhz = Nz / 2;

		C.set(nx, ny, nz);
		if (bST) {
			if (bST_diagonal) {
				Cxx.set(nx, ny, nz); Cyy.set(nx, ny, nz); Czz.set(nx, ny, nz);
			}
			else Ctr.set(nx, ny, nz);
		}

		Nyz = Ny * Nz; Nv = Nx * Nyz;
		ftd.SetArray(Nv); ftc.SetArray(Nv);

		init_fft_plans();
	};
	void cell_dim(float dx, float dy, float dz) {
		this->xd[0] = dx; this->xd[1] = dy; this->xd[2] = dz;
		this->md[0] = this->xd[0] / this->Nx; 
		this->md[1] = this->xd[1] / this->Ny;
		this->md[2] = this->xd[2] / this->Nz;

		this->V = dx * dy * dz;
		this->inv_V = 4 * PI / this->V;
	};
	void init_k(/*float dx, float dy, float dz*/) {
		int i;
		double dk = PI2 / this->xd[0]; //dx;
		for (i = 0; i < k[0].n; i++) k[0].m[i] = dk * (i <= nhx ? i : i - Nx);
		dk = PI2 / this->xd[1]; //dy;
		for (i = 0; i < k[1].n; i++) k[1].m[i] = dk * (i <= nhy ? i : i - Ny);
		dk = PI2 / this->xd[2]; //dz;
		for (i = 0; i < k[2].n; i++) k[2].m[i] = dk * (i <= nhz ? i : i - Nz);
	};
	void init_b() {
		int i, k, n;
		double phi = 0, Mn;
		complex vc;
		for (i = 0; i < Nx; i++) {
			vc.x = 0; vc.y = 0; 
			for (k = 0; k <= N - 2; k++) {
				//phi = PI2 * i * k / Nx;
				phi = PI2 * (i <= nhx ? i : i - Nx) * k / Nx;
				n = bsp->mpt * (k + 1);
				Mn = bsp->y[0][n];
				vc.x += Mn * cos(phi); vc.y += Mn * sin(phi);
			}
			//phi = PI2 * (N - 1) * i / Nx;
			//b[0].m[i] = complex(cos(phi), sin(phi)) / vc;
			b[0].m[i] = 1.0 / (vc.x * vc.x + vc.y * vc.y);
		}

		for (i = 0; i < Ny; i++) {
			vc.x = 0; vc.y = 0; 
			for (k = 0; k <= N - 2; k++) {
				//phi = PI2 * i * k / Ny;
				phi = PI2 * (i <= nhy ? i : i - Ny) * k / Ny;
				n = bsp->mpt * (k + 1);
				Mn = bsp->y[0][n];
				vc.x += Mn * cos(phi); vc.y += Mn * sin(phi);
			}
			//phi = PI2 * (N - 1) * i / Ny;
			//b[1].m[i] = complex(cos(phi), sin(phi)) / b[1].m[i];
			b[1].m[i] = 1.0 / (vc.x * vc.x + vc.y * vc.y);
		}

		for (i = 0; i < Nz; i++) {
			vc.x = 0; vc.y = 0; 
			for (k = 0; k <= N - 2; k++) {
				//phi = PI2 * i * k / Nz;
				phi = PI2 * (i <= nhz ? i : i - Nz) * k / Nz;
				n = bsp->mpt * (k + 1);
				Mn = bsp->y[0][n];
				vc.x += Mn * cos(phi); vc.y += Mn * sin(phi);
			}
			//phi = PI2 * (N - 1) * i / Nx;
			//b[2].m[i] = complex(cos(phi), sin(phi)) / b[2].m[i];
			b[2].m[i] = 1.0 / (vc.x * vc.x + vc.y * vc.y);
		}
	};
	void init_C() {
		// important: the summation in reciprocal space is shifted to the center of FFT space so that we can use FFT
		int ni, nj, nk;
		double kx2, ky2, kz2, k2;
		__VSPME<atom>::kappa = _KAPPA / ::rcut_Ewd;
		double a2 = 4 * __VSPME<atom>::kappa * __VSPME<atom>::kappa;
		double f = 0;
		for (ni = 0; ni < C.nx; ni++) {
			kx2 = k[0].m[ni] * k[0].m[ni];
			for (nj = 0; nj < C.ny; nj++) {
				ky2 = k[1].m[nj] * k[1].m[nj];
				for (nk = 0; nk < C.nz; nk++) {
					kz2 = k[2].m[nk] * k[2].m[nk];
					k2 = kx2 + ky2 + kz2;
					if (ni == 0 && nj == 0 && nk == 0) {
						C.m[ni].m[nj].m[nk] = 0;
						if (bST) {
							if (bST_diagonal) {
								Cxx.m[ni].m[nj].m[nk] = 0;
								Cyy.m[ni].m[nj].m[nk] = 0;
								Czz.m[ni].m[nj].m[nk] = 0;
							}
							else Ctr.m[ni].m[nj].m[nk] = 0;
						}
					}
					else {
						f = exp(-(k2) / a2) / k2;
						if (f < 1e-8) f = 0;
						C.m[ni].m[nj].m[nk] = f;
						if (bST) {
							if (bST_diagonal) {
								Cxx.m[ni].m[nj].m[nk] = f * (1 - 2 * (1 + k2 / a2) / k2 * kx2);
								Cyy.m[ni].m[nj].m[nk] = f * (1 - 2 * (1 + k2 / a2) / k2 * ky2);
								Czz.m[ni].m[nj].m[nk] = f * (1 - 2 * (1 + k2 / a2) / k2 * kz2);
							}
							else Ctr.m[ni].m[nj].m[nk] = f * (3 - 2 * (1 + k2 / a2));
						}
					}
				}
			}
		}
	};
};

template <class atom> class VSPME : public _VSPME<atom> {
public:
	short mMP; // 0 -- charge only, 1 -- with dipole
	FLEX_ASYM_CMATRIX<double> Q, BCQ; // C[i][j][k] = B[i][j][k] * real C[i][j][k] 
	FLEX_ASYM_CMATRIX<complex> FQ, FBCQ;

	// for stress tensor
	//FLEX_ASYM_CMATRIX<double> BTQ;
	FLEX_ASYM_CMATRIX<complex> FBTQ;

	// buffers used for multi-thread Q calculation
	FLEX_ASYM_CMATRIX<double> Qbuf[MAX_THREADS];

	void init_vars() {
		int i;
		Q.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);
		BCQ.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);
		FQ.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);
		FBCQ.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);
		if (_VSPME<atom>::bST) {
			//BTQ.set(Nx, Ny, Nz); 
			FBTQ.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);
		}

		for (i = 0; i < MAX_THREADS; i+=1) Qbuf[i].set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);

		_VSPME<atom>::init_absp(_VSPME<atom>::va.n, _VSPME<atom>::bsp->indx);
	};

	double cal_U_recp(bool bU = true);
	double cal_ST_recp_0(FLEX_ASYM_CMATRIX<double> &CT);
	void cal_StressTensor_recp_0() {
		if (!_VSPME<atom>::bST) return;
		if (_VSPME<atom>::bST_diagonal) {
			_VSPME<atom>::STxx = cal_ST_recp_0(_VSPME<atom>::Cxx);
			_VSPME<atom>::STyy = cal_ST_recp_0(_VSPME<atom>::Cyy);
			_VSPME<atom>::STzz = cal_ST_recp_0(_VSPME<atom>::Czz);
			_VSPME<atom>::STtr = _VSPME<atom>::STxx + _VSPME<atom>::STyy + _VSPME<atom>::STzz;
		}
		else _VSPME<atom>::STtr = cal_ST_recp_0(_VSPME<atom>::Ctr);
	};
	
	VSPME() {mMP = 0;};
};

template <class atom> double VSPME<atom>::cal_U_recp(bool bU) {
	// Q ==> F[Q]
	int i, j, k, n;
	int inx, iny, inz;

	// Q ==>ftd
	n = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		//nx = i * Nyz;
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			//ny = j * Nz;
			/*
			for (k = 0; k < _VSPME<atom>::Nz; k++) {
				//ftd.m[nx + ny + k] = Q.m[i].m[j].m[k];
				_VSPME<atom>::ftd.m[n] = Q.m[i].m[j].m[k]; n++;
			}
			*/
			memcpy(_VSPME<atom>::ftd.m + n, Q.m[i].m[j].m, _VSPME<atom>::Nz * SIZE_DOUBLE);
			n += _VSPME<atom>::Nz;
		}
	}
	fftw_execute(_VSPME<atom>::ft_plan); // ftd ==FFT==> ftc
	// ftc ==> FQ
	n = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			/*
			for (k = 0; k <= _VSPME<atom>::nhz; k++) {
				FQ.m[i].m[j].m[k].x = _VSPME<atom>::ftc.m[n][0];
				FQ.m[i].m[j].m[k].y = _VSPME<atom>::ftc.m[n][1];
				n++;
			}
			*/
			memcpy(FQ.m[i].m[j].m, _VSPME<atom>::ftc.m + n, (_VSPME<atom>::nhz + 1) * sizeof(complex));
			n += _VSPME<atom>::nhz + 1;
		}
	}
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		inx = (i == 0 ? 0 : _VSPME<atom>::Nx - i); // the inverse index 
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			iny = (j == 0 ? 0 : (_VSPME<atom>::Ny - j));
			for (k = _VSPME<atom>::nhz + 1; k < _VSPME<atom>::Nz; k++) {
				inz = _VSPME<atom>::Nz - k; // k > hz
				FQ.m[i].m[j].m[k].x = FQ.m[inx].m[iny].m[inz].x;
				FQ.m[i].m[j].m[k].y = -FQ.m[inx].m[iny].m[inz].y;
			}
		}
	}
	/*
	{
		char msg[256] = "\0", title[50] = "Q", fname[256] = "\0";
		double v = 0;
		int i, j, k;
		ofstream *out = NULL;
		for (i = 0; i < Nx; i++) {
		sprintf(fname, "%s-%d.dat", title, i);
		out = new ofstream;
		out->open(fname);
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) {
				v = Q.m[i].m[j].m[k];
				sprintf(msg, "%d %d  %f ", j, k, v);
				(*out)<<msg<<endl;
			}
		}
		out->close();
		delete out; out = NULL;
		}
	}
	
	{
		char msg[256] = "\0", title[50] = "recp", fname[256] = "\0";
		double v = 0;
		int i, j, k;
		ofstream *out = NULL;
		for (i = 0; i < Nx; i++) {
		sprintf(fname, "%s-%d.dat", title, i);
		out = new ofstream;
		out->open(fname);
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) {
				v = sqrt(FQ.m[i].m[j].m[k].x * FQ.m[i].m[j].m[k].x + FQ.m[i].m[j].m[k].y * FQ.m[i].m[j].m[k].y);
				sprintf(msg, "%d %d  %f ", j, k, v);
				(*out)<<msg<<endl;
			}
		}
		out->close();
		delete out; out = NULL;
		}
	}
	*/


	// Energy
	double vc, vb, vbc;
	double E = 0;
	complex *fq = NULL;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {for (j = 0; j < _VSPME<atom>::Ny; j++) {for (k = 0; k < _VSPME<atom>::Nz; k++) {
		if (_VSPME<atom>::C.m[i].m[j].m[k] == 0) continue;
		vc = _VSPME<atom>::C.m[i].m[j].m[k]; vb = _VSPME<atom>::b[0].m[i] * _VSPME<atom>::b[1].m[j] * _VSPME<atom>::b[2].m[k]; fq = &(FQ.m[i].m[j].m[k]);
		vbc = vc * vb;
		FBCQ.m[i].m[j].m[k].x = vbc * fq->x; FBCQ.m[i].m[j].m[k].y = vbc * fq->y;
	}}}

	// FBCQ ==> ftc
	n = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			/*
			for (k = 0; k <= _VSPME<atom>::nhz; k++) {
				_VSPME<atom>::ftc.m[n][0] = FBCQ.m[i].m[j].m[k].x;
				_VSPME<atom>::ftc.m[n][1] = FBCQ.m[i].m[j].m[k].y;
				n++;
			}
			*/
			memcpy(_VSPME<atom>::ftc.m + n, FBCQ.m[i].m[j].m, (_VSPME<atom>::nhz + 1) * sizeof(complex));
			n += _VSPME<atom>::nhz + 1;
		}
	}
	fftw_execute(_VSPME<atom>::ift_plan); // ftc == inverse ft ==> ftd
	//ftd ==> BCQ, here BCQ is inverse FT of FBCQ
	n = 0; E = 0;
	if (bU) {
		for (i = 0; i < _VSPME<atom>::Nx; i++) {
			for (j = 0; j < _VSPME<atom>::Ny; j++) {
				for (k = 0; k < _VSPME<atom>::Nz; k++) {
					E += Q.m[i].m[j].m[k] * _VSPME<atom>::ftd.m[n];
					BCQ.m[i].m[j].m[k] = _VSPME<atom>::ftd.m[n]; n++;
				}
			}
		}
	}
	else {
		for (i = 0; i < _VSPME<atom>::Nx; i++) {
			for (j = 0; j < _VSPME<atom>::Ny; j++) {
				/*
				for (k = 0; k < _VSPME<atom>::Nz; k++) {
					E += Q.m[i].m[j].m[k] * _VSPME<atom>::ftd.m[n];
					BCQ.m[i].m[j].m[k] = _VSPME<atom>::ftd.m[n]; n++;
				}
				*/
				memcpy(BCQ.m[i].m[j].m, _VSPME<atom>::ftd.m + n, _VSPME<atom>::Nz * SIZE_DOUBLE);
				n += _VSPME<atom>::Nz;
			}
		}
	}

	E *= eUnit_estat_kT * PI2 / _VSPME<atom>::V;
	/*
	{
		char msg[256] = "\0", title[50] = "Fcbq", fname[256] = "\0";
		double v = 0;
		int i, j, k;
		ofstream *out = NULL;
		for (i = 0; i < Nx; i++) {
		sprintf(fname, "%s-%d.dat", title, i);
		out = new ofstream;
		out->open(fname);
		for (j = 0; j < Ny; j++) {
			for (k = 0; k < Nz; k++) {
				v = Q.m[i].m[j].m[k];
				sprintf(msg, "%d %d  %f ", j, k, v);
				(*out)<<msg<<endl;
			}
		}
		out->close();
		delete out; out = NULL;
		}
	}
	*/

	return E;
};

//calculate stress sensor with FFTW method
#define _ST_METHOD_ 1
#if _ST_METHOD_ == 0
template <class atom> double VSPME<atom>::cal_ST_recp_0(FLEX_ASYM_CMATRIX<double> &CST) {
	if (!_VSPME<atom>::bST) return 0;
	int i, j, k, n;
	int inx, iny, inz;

	double vc, vb, vbc;
	complex *fq = NULL;
	double ST = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {for (j = 0; j < _VSPME<atom>::Ny; j++) {for (k = 0; k < _VSPME<atom>::Nz; k++) {
		if (CST.m[i].m[j].m[k] == 0) continue;
		vc = CST.m[i].m[j].m[k]; vb = _VSPME<atom>::b[0].m[i] * _VSPME<atom>::b[1].m[j] * _VSPME<atom>::b[2].m[k]; fq = &(FQ.m[i].m[j].m[k]);
		vbc = vc * vb;
		FBTQ.m[i].m[j].m[k].x = vbc * fq->x; FBTQ.m[i].m[j].m[k].y = vbc * fq->y;
	}}}

	// FBTQ ==> ftc
	n = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			/*
			for (k = 0; k <= _VSPME<atom>::nhz; k++) {
				_VSPME<atom>::ftc.m[n][0] = FBTQ.m[i].m[j].m[k].x;
				_VSPME<atom>::ftc.m[n][1] = FBTQ.m[i].m[j].m[k].y;
				n++;
			}
			*/
			memcpy(_VSPME<atom>::ftc.m + n, FBTQ.m[i].m[j].m, (_VSPME<atom>::nhz + 1) * sizeof(complex));
			n += _VSPME<atom>::nhz + 1;
		}
	}
	fftw_execute(_VSPME<atom>::ift_plan); // ftc == inverse ft ==> ftd
	//ftd ==> BTQ, here BTQ is inverse FFT of FBCQ
	n = 0; ST = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			for (k = 0; k < _VSPME<atom>::Nz; k++) {
				ST += Q.m[i].m[j].m[k] * _VSPME<atom>::ftd.m[n];
				//BCQ.m[i].m[j].m[k] = ftd.m[n]; 
				n++;
			}
		}
	}

	ST *= PI2 / _VSPME<atom>::V * ::fUnit_estat_mass;

	return ST;
};
#else
// calculate directly from reciprocal space
template <class atom> double VSPME<atom>::cal_ST_recp_0(FLEX_ASYM_CMATRIX<double> &CST) {
	if (!_VSPME<atom>::bST) return 0;
	int i, j, k, n;
	int inx, iny, inz;

	double vc, vb, vbc;
	complex *fq = NULL;
	double ST = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {for (j = 0; j < _VSPME<atom>::Ny; j++) {for (k = 0; k < _VSPME<atom>::Nz; k++) {
		if (CST.m[i].m[j].m[k] == 0) continue;
		vc = CST.m[i].m[j].m[k]; vb = _VSPME<atom>::b[0].m[i] * _VSPME<atom>::b[1].m[j] * _VSPME<atom>::b[2].m[k]; fq = &(FQ.m[i].m[j].m[k]);
		vbc = vc * vb;
		//FBTQ.m[i].m[j].m[k].x = vbc * fq->x; FBTQ.m[i].m[j].m[k].y = vbc * fq->y;
		ST += vbc * (fq->x * fq->x + fq->y * fq->y);
	}}}

	ST *= PI2 / _VSPME<atom>::V * ::fUnit_estat_mass;

	return ST;

	// FBTQ ==> ftc
	n = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			/*
			for (k = 0; k <= _VSPME<atom>::nhz; k++) {
				_VSPME<atom>::ftc.m[n][0] = FBTQ.m[i].m[j].m[k].x;
				_VSPME<atom>::ftc.m[n][1] = FBTQ.m[i].m[j].m[k].y;
				n++;
			}
			*/
			memcpy(_VSPME<atom>::ftc.m + n, FBTQ.m[i].m[j].m, (_VSPME<atom>::nhz + 1) * sizeof(complex));
			n += _VSPME<atom>::nhz + 1;
		}
	}
	fftw_execute(_VSPME<atom>::ift_plan); // ftc == inverse ft ==> ftd
	//ftd ==> BTQ, here BTQ is inverse FFT of FBCQ
	n = 0; ST = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			for (k = 0; k < _VSPME<atom>::Nz; k++) {
				ST += Q.m[i].m[j].m[k] * _VSPME<atom>::ftd.m[n];
				//BCQ.m[i].m[j].m[k] = ftd.m[n]; 
				n++;
			}
		}
	}

	ST *= PI2 / _VSPME<atom>::V * ::fUnit_estat_mass;

	return ST;
};
#endif

template <class atom> void init_normalized_coordinates_tmpl(_VSPME<atom> *vspme, int na1, int na2) {
	int iatom;
	atom *patom = NULL;
	for (iatom = na1; iatom <= na2; iatom++) {
		if (iatom >= vspme->va.n) break;
		patom = vspme->va.m + iatom;
		patom->u.v[0] = (patom->r->v[0] - vspme->xl[0]) / vspme->md[0];
		patom->u.v[1] = (patom->r->v[1] - vspme->xl[1]) / vspme->md[1];
		patom->u.v[2] = (patom->r->v[2] - vspme->xl[2]) / vspme->md[2];

		if (patom->u.v[0] < 0) patom->u.v[0] += vspme->Nx;
		else if (patom->u.v[0] >= vspme->Nx) patom->u.v[0] -= vspme->Nx;

		if (patom->u.v[1] < 0) patom->u.v[1] += vspme->Ny;
		else if (patom->u.v[1] >= vspme->Ny) patom->u.v[1] -= vspme->Ny;

		if (patom->u.v[2] < 0) patom->u.v[2] += vspme->Nz;
		else if (patom->u.v[2] >= vspme->Nz) patom->u.v[2] -= vspme->Nz;
	}
};

inline void init_normalized_coordinates_q(VSPME<_VQatom> *vspme, int na1, int na2) {
	_VSPME<_VQatom>* _vspme = (_VSPME<_VQatom> *)vspme;
	init_normalized_coordinates_tmpl<_VQatom>(_vspme, na1, na2);
};

inline void init_normalized_coordinates(VSPME<_VQatom> &vspme, int nThreads) {
	MacroMoleculeOperate< VSPME<_VQatom> >((void*)(&init_normalized_coordinates_q), vspme, 0, vspme.va.n - 1, nThreads);
	vspme.bReady_absp = false;
}

inline void init_normalized_coordinates_mu(VSPME<_VMUatom> *vspme, int na1, int na2) {
	_VSPME<_VMUatom>* _vspme = (_VSPME<_VMUatom>*)vspme;
	init_normalized_coordinates_tmpl< _VMUatom >(_vspme, na1, na2);
}

inline void init_normalized_coordinates(VSPME<_VMUatom> &vspme, int nThreads) {
	MacroMoleculeOperate< VSPME<_VMUatom> >((void*)(&init_normalized_coordinates_mu), vspme, 0, vspme.va.n - 1, nThreads);
	/*
	int ia;
	char msg[256] = "\0";
	for (ia = 0; ia < vspme.va.n; ia++) {
		if (vspme.va.m[ia].u.v[0] >= vspme.Nx || vspme.va.m[ia].u.v[0] < 0 ||
			vspme.va.m[ia].u.v[1] >= vspme.Ny || vspme.va.m[ia].u.v[1] < 0 ||
			vspme.va.m[ia].u.v[2] >= vspme.Nz || vspme.va.m[ia].u.v[2] < 0) {
				sprintf(msg, "atom %d has coordinates [%f, %f, %f], normalized coordinates [%f, %f, %f]", ia,
					vspme.va.m[ia].r->v[0], vspme.va.m[ia].r->v[1], vspme.va.m[ia].r->v[2],
					vspme.va.m[ia].u.v[0], vspme.va.m[ia].u.v[1], vspme.va.m[ia].u.v[2]);
				show_log(msg, true);
		}
	}
	*/
	vspme.bReady_absp = false;
}

/***************************************************************************************
	CALCULATE Q, Q is accumulated in the calculation, without resetting
	And, mesh dimension is assumed to be bigger than the order of the B-spline function
	In another word, the spline points is within the nearest neighbor cell only
***************************************************************************************/

class RECP_SPME_RES {
public:
	short mMP;
	FLEX_ASYM_CMATRIX<double> *Q;
	void *f_MU_accumulate;
	void *f_mu_accumulate;

	RECP_SPME_RES() {mMP = 0; f_MU_accumulate = NULL; f_mu_accumulate = NULL; Q = NULL;};

	void operator = (RECP_SPME_RES &var) {
		this->mMP = var.mMP;
		// the dimension will be checked automatically in operator = for FLEX_ASYM_CMATRIX
		//this->Q = var.Q;
		this->f_MU_accumulate = var.f_MU_accumulate;
		this->f_mu_accumulate = var.f_mu_accumulate;
	};
};

template <class atom, class VSPME, class RES> void MU_accumulate_tmpl(atom *patom, VSPME* vspme, double* Mn, double* dMn, RES& res, int kx, int ky, int kz) {
	if (patom->mMP == 0 || patom->mu == NULL) return;

	double *mu = patom->mu->v;
	float *md = vspme->md;

	// accumulating
	res.Q->m[kx].m[ky].m[kz] += mu[0] * dMn[0] * Mn[1] * Mn[2] / md[0] + mu[1] * Mn[0] * dMn[1] * Mn[2] / md[1] + mu[2] * Mn[0] * Mn[1] * dMn[2] / md[2];
};

inline void MU_accumulate(_VMUatom *atom, VSPME<_VMUatom>* vspme, double *Mn, double* dMn, RECP_SPME_RES& res, int kx, int ky, int kz) {
	return MU_accumulate_tmpl<_VMUatom, VSPME<_VMUatom>, RECP_SPME_RES>(atom, vspme, Mn, dMn, res, kx, ky, kz);
};

template <class atom, class VSPME, class VAR, class RES> void VSPME_Q_accum_tmpl(VSPME *vspme, int na1, int na2, RES& res) {
	typedef void (*FUNC)(atom*, VSPME*, double*, double*, RES&, int, int, int);

	float q;
	int i, j, k, nx, ny, nz, kx, ky, kz;
	int iatom, ipt;
	atom *patom = NULL;
	AtomBSP *pabsp = NULL;
	int ib;
	//VECTOR3 r;
	double Mn[3], dMn[3], ux, uy, uz, dux, duy, duz;
	BSplineFunc<2> *bsp = vspme->bsp;
	int nb = bsp->indx + 1;
	for (iatom = na1; iatom <= na2; iatom++) {
		if (iatom >= vspme->va.n) break;
		patom = vspme->va.m + iatom;
		q = *(patom->q);
		if (patom->eNeutral) continue;
		pabsp = vspme->absp.m + iatom;
		/*
		r.v[0] = patom->r->v[0] - vspme->xl[0]; // make 0 is relative to the corner of the cell
		if (r.v[0] < 0) r.v[0] += vspme->xd[0];
		else if (r.v[0] >= vspme->xd[0]) r.v[0] -= vspme->xd[0];
		ux = r.v[0] / vspme->md[0]; nx = int(ux);
		//if (nx < 0) nx = 0;
		//else if (nx >= vspme->Nx) nx = vspme->Nx - 1;
		*/
		ux = patom->u.v[0]; nx = int(ux);

		/*
		r.v[1] = patom->r->v[1] - vspme->xl[1]; // make 0 is relative to the corner of the cell
		if (r.v[1] < 0) r.v[1] += vspme->xd[1];
		else if (r.v[1] >= vspme->xd[1]) r.v[1] -= vspme->xd[1];
		uy = r.v[1] / vspme->md[1]; ny = int(uy); 
		//if (ny < 0) ny = 0;
		//else if (ny >= vspme->Ny) ny = vspme->Ny - 1;
		*/
		uy = patom->u.v[1]; ny = int(uy);

		/*
		r.v[2] = patom->r->v[2] - vspme->xl[2]; // make 0 is relative to the corner of the cell
		if (r.v[2] < 0) r.v[2] += vspme->xd[2];
		else if (r.v[2] >= vspme->xd[2]) r.v[2] -= vspme->xd[2];
		uz = r.v[2] / vspme->md[2]; nz = int(uz); 
		//if (nz < 0) nz = 0;
		//else if (nz >= vspme->Nz) nz = vspme->Nz - 1;
		*/
		uz = patom->u.v[2]; nz = int(uz);
		
		for (i = -nb; i <= nb; i++) {
			kx = nx + i; 
			dux = ux - kx; if (dux < 0) break;
			ib = int(dux); if (ib >= bsp->indx) continue;
			ipt = bsp->get_pt(dux); Mn[0] = bsp->bsp_interpolate(ipt, dux); pabsp->Mn0.m[ib] = Mn[0]; if (Mn[0] == 0) continue;
			dMn[0] = bsp->dbsp_interpolate(ipt, dux); pabsp->dMn0.m[ib] = dMn[0];
			if (patom->mMP > 0) {pabsp->d2Mn0.m[ib] = bsp->dbsp_interpolate(ipt, dux, 2);}
			if (kx < 0) kx += vspme->Nx; else if (kx >= vspme->Nx) kx -= vspme->Nx;
			for (j = -nb; j <= nb; j++) {
				ky = ny + j; 
				duy = uy - ky; if (duy < 0) break;
				ib = int(duy); if (ib >= bsp->indx) continue;
				ipt = bsp->get_pt(duy); Mn[1] = bsp->bsp_interpolate(ipt, duy); pabsp->Mn1.m[ib] = Mn[1]; if (Mn[1] == 0) continue;
				dMn[1] = bsp->dbsp_interpolate(ipt, duy); pabsp->dMn1.m[ib] = dMn[1];
				if (patom->mMP > 0) {pabsp->d2Mn1.m[ib] = bsp->dbsp_interpolate(ipt, duy, 2);}
				if (ky < 0) ky += vspme->Ny; else if (ky >= vspme->Ny) ky -= vspme->Ny;
				for (k = -nb; k <= nb; k++) {
					kz = nz + k; 
					duz = uz - kz; if (duz < 0) break;
					ib = int(duz); if (ib >= bsp->indx) continue;
					ipt = bsp->get_pt(duz); Mn[2] = bsp->bsp_interpolate(ipt, duz); pabsp->Mn2.m[ib] = Mn[2]; if (Mn[2] == 0) continue;
					dMn[2] = bsp->dbsp_interpolate(ipt, duz); pabsp->dMn2.m[ib] = dMn[2];
					if (patom->mMP > 0) {pabsp->d2Mn2.m[ib] = bsp->dbsp_interpolate(ipt, duz, 2);}
					if (kz < 0) kz += vspme->Nz; else if (kz >= vspme->Nz) kz -= vspme->Nz;

					//vspme->Q->m[kx].m[ky].m[kz] += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
					res.Q->m[kx].m[ky].m[kz] += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
					if (patom->mMP > 0) {
						((FUNC)res.f_mu_accumulate)(patom, vspme, Mn, dMn, res, kx, ky, kz);
					}
				}
			}
		}
	}
};


// assuming Mn[3], dMn[3] and d2Mn[3] for each atom were calculated
template <class atom, class VSPME, class VAR, class RES> void VSPME_Q_accum_tmpl_absp(VSPME *vspme, int na1, int na2, RES& res) {
	typedef void (*FUNC)(atom*, VSPME*, double*, double*, RES&, int, int, int);

	float q;
	int i, j, k, nx, ny, nz, kx, ky, kz;
	int iatom;//, ipt;
	atom *patom = NULL;
	AtomBSP *pabsp = NULL;
	int ib;
	//VECTOR3 r;
	double Mn[3], dMn[3], ux, uy, uz, dux, duy, duz;
	BSplineFunc<2> *bsp = vspme->bsp;
	int nb = bsp->indx + 1;
	for (iatom = na1; iatom <= na2; iatom++) {
		if (iatom >= vspme->va.n) break;
		patom = vspme->va.m + iatom;
		q = *(patom->q);
		if (patom->eNeutral) continue;
		pabsp = vspme->absp.m + iatom;
		/*
		r.v[0] = patom->r->v[0] - vspme->xl[0]; // make 0 is relative to the corner of the cell
		if (r.v[0] < 0) r.v[0] += vspme->xd[0];
		else if (r.v[0] >= vspme->xd[0]) r.v[0] -= vspme->xd[0];
		ux = r.v[0] / vspme->md[0]; nx = int(ux);
		//if (nx < 0) nx = 0;
		//else if (nx >= vspme->Nx) nx = vspme->Nx - 1;
		*/
		ux = patom->u.v[0]; nx = int(ux);

		/*
		r.v[1] = patom->r->v[1] - vspme->xl[1]; // make 0 is relative to the corner of the cell
		if (r.v[1] < 0) r.v[1] += vspme->xd[1];
		else if (r.v[1] >= vspme->xd[1]) r.v[1] -= vspme->xd[1];
		uy = r.v[1] / vspme->md[1]; ny = int(uy); 
		//if (ny < 0) ny = 0;
		//else if (ny >= vspme->Ny) ny = vspme->Ny - 1;
		*/
		uy = patom->u.v[1]; ny = int(uy);

		/*
		r.v[2] = patom->r->v[2] - vspme->xl[2]; // make 0 is relative to the corner of the cell
		if (r.v[2] < 0) r.v[2] += vspme->xd[2];
		else if (r.v[2] >= vspme->xd[2]) r.v[2] -= vspme->xd[2];
		uz = r.v[2] / vspme->md[2]; nz = int(uz); 
		//if (nz < 0) nz = 0;
		//else if (nz >= vspme->Nz) nz = vspme->Nz - 1;
		*/
		uz = patom->u.v[2]; nz = int(uz);
		
		for (i = -nb; i <= nb; i++) {
			kx = nx + i; 
			dux = ux - kx; if (dux < 0) break;
			ib = int(dux); if (ib >= bsp->indx) continue;
			//ipt = bsp->get_pt(dux); Mn[0] = bsp->bsp_interpolate(ipt, dux); pabsp->Mn0.m[ib] = Mn[0]; if (Mn[0] == 0) continue;
			//dMn[0] = bsp->dbsp_interpolate(ipt, dux); pabsp->dMn0.m[ib] = dMn[0];
			//if (patom->mMP > 0) {pabsp->d2Mn0.m[ib] = bsp->dbsp_interpolate(ipt, dux, 2);}
			Mn[0] = pabsp->Mn0.m[ib]; if (patom->mMP > 0) dMn[0] = pabsp->dMn0.m[ib];
			if (kx < 0) kx += vspme->Nx; else if (kx >= vspme->Nx) kx -= vspme->Nx;
			for (j = -nb; j <= nb; j++) {
				ky = ny + j; 
				duy = uy - ky; if (duy < 0) break;
				ib = int(duy); if (ib >= bsp->indx) continue;
				//ipt = bsp->get_pt(duy); Mn[1] = bsp->bsp_interpolate(ipt, duy); pabsp->Mn1.m[ib] = Mn[1]; if (Mn[1] == 0) continue;
				//dMn[1] = bsp->dbsp_interpolate(ipt, duy); pabsp->dMn1.m[ib] = dMn[1];
				//if (patom->mMP > 0) {pabsp->d2Mn1.m[ib] = bsp->dbsp_interpolate(ipt, duy, 2);}
				Mn[1] = pabsp->Mn1.m[ib]; if (patom->mMP > 0) dMn[1] = pabsp->dMn1.m[ib];
				if (ky < 0) ky += vspme->Ny; else if (ky >= vspme->Ny) ky -= vspme->Ny;
				for (k = -nb; k <= nb; k++) {
					kz = nz + k; 
					duz = uz - kz; if (duz < 0) break;
					ib = int(duz); if (ib >= bsp->indx) continue;
					//ipt = bsp->get_pt(duz); Mn[2] = bsp->bsp_interpolate(ipt, duz); pabsp->Mn2.m[ib] = Mn[2]; if (Mn[2] == 0) continue;
					//dMn[2] = bsp->dbsp_interpolate(ipt, duz); pabsp->dMn2.m[ib] = dMn[2];
					//if (patom->mMP > 0) {pabsp->d2Mn2.m[ib] = bsp->dbsp_interpolate(ipt, duz, 2);}
					Mn[2] = pabsp->Mn2.m[ib]; if (patom->mMP > 0) dMn[2] = pabsp->dMn2.m[ib];
					if (kz < 0) kz += vspme->Nz; else if (kz >= vspme->Nz) kz -= vspme->Nz;

					//vspme->Q->m[kx].m[ky].m[kz] += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
					res.Q->m[kx].m[ky].m[kz] += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
					if (patom->mMP > 0) {
						((FUNC)res.f_mu_accumulate)(patom, vspme, Mn, dMn, res, kx, ky, kz);
					}
				}
			}
		}
	}
};



inline void VSPME_Q_accum_q(VSPME<_VQatom> *vspme, int na1, int na2, RECP_SPME_RES& res) {
	if (vspme->bReady_absp) return VSPME_Q_accum_tmpl_absp<_VQatom, VSPME<_VQatom>, RECP_SPME_RES >(vspme, na1, na2, res);
	else return VSPME_Q_accum_tmpl<_VQatom, VSPME<_VQatom>, RECP_SPME_RES >(vspme, na1, na2, res);
};

inline void VSPME_Q_accum_mu(VSPME<_VMUatom> *vspme, int na1, int na2, RECP_SPME_RES& res) {
	if (vspme->bReady_absp) return VSPME_Q_accum_tmpl_absp<_VMUatom, VSPME<_VMUatom>, RECP_SPME_RES >(vspme, na1, na2, res);
	else return VSPME_Q_accum_tmpl<_VMUatom, VSPME<_VMUatom>, RECP_SPME_RES >(vspme, na1, na2, res);
};

void cal_Q_spme(VSPME<_VQatom>& vspme, int nThreads);

void cal_Q_spme(VSPME<_VMUatom>& vspme, int nThreads);


/***************************************************************************************
	CALCULATE E_recp -- E_recp = F_recp / q
	Mesh dimension is assumed to be bigger than the order of the B-spline function
	In another word, the spline points is within the nearest neighbor cell only
***************************************************************************************/

class ERECP_SPME_VAR {
public:
	FLEX_ASYM_CMATRIX<double> *BCQ; // here Q is inverse FT of F[BCQ] 

	bool bE, bF;
	void* fFunc;

	ERECP_SPME_VAR() {BCQ = NULL; bE = true; bF = false; fFunc = NULL;};
	void operator = (ERECP_SPME_VAR &var) {
		BCQ = var.BCQ; bE = var.bE; bF = var.bF; fFunc = var.fFunc;
	};
};

template <class atom, class VSPME, class EVAR> void F_spme_accumulate_tmpl(atom *patom, VSPME *vspme, double* Mn, double* dMn, double* d2Mn, EVAR& var, int kx, int ky, int kz) {
	double tv_q[3] = {0, 0, 0}, tv_mu[3] = {0, 0, 0};
	float *md = vspme->md;
	double *mu = patom->mu->v;
	float q = *(patom->q);
	if (*(patom->q) != 0) {
		tv_q[0] = q * dMn[0] / md[0] * Mn[1] * Mn[2];
		tv_q[1] = q * Mn[0] * dMn[1] / md[1] * Mn[2];
		tv_q[2] = q * Mn[0] * Mn[1] * dMn[2] / md[2];
	}
	if (patom->mMP > 0 && mu != NULL) {
		tv_mu[0] = mu[0] * d2Mn[0] * Mn[1] * Mn[2] / (md[0] * md[0]) + mu[1] * dMn[0] * dMn[1] * Mn[2] / (md[0] * md[1]) + mu[2] * dMn[0] * Mn[1] * dMn[2] / (md[0] * md[2]);
		tv_mu[1] = mu[0] * dMn[0] * dMn[1] * Mn[2] / (md[0] * md[1]) + mu[1] * Mn[0] * d2Mn[1] * Mn[2] / (md[1] * md[1]) + mu[2] * Mn[0] * dMn[1] * dMn[2] / (md[1] * md[2]);
		tv_mu[2] = mu[0] * dMn[0] * Mn[1] * dMn[2] / (md[0] * md[2]) + mu[1] * Mn[0] * dMn[1] * dMn[2] / (md[1] * md[2]) + mu[2] * Mn[0] * Mn[1] * d2Mn[2] / (md[2] * md[2]);
	}
	patom->F_recp.v[0] -= vspme->inv_V * (tv_q[0] + tv_mu[0]) * var.BCQ->m[kx].m[ky].m[kz]; 
	patom->F_recp.v[1] -= vspme->inv_V * (tv_q[1] + tv_mu[1]) * var.BCQ->m[kx].m[ky].m[kz];
	patom->F_recp.v[2] -= vspme->inv_V * (tv_q[2] + tv_mu[2]) * var.BCQ->m[kx].m[ky].m[kz];
};

inline void F_spme_accumulate_mu(_VMUatom *patom, VSPME<_VMUatom> *vspme, double* Mn, double* dMn, double* d2Mn, ERECP_SPME_VAR& var, int kx, int ky, int kz) {
	return F_spme_accumulate_tmpl<_VMUatom, VSPME<_VMUatom>, ERECP_SPME_VAR>(patom, vspme, Mn, dMn, d2Mn, var, kx, ky, kz);
};

template <class atom, class VSPME, class EVAR> void E_recp_tmpl(VSPME *vspme, int na1, int na2, EVAR &var) { // E of atom is accumulated ,without initilized
	typedef void (*FUNC)(atom*, VSPME*, double*, double*, double*, EVAR&, int, int, int);

	int i, j, k, nx, ny, nz, kx, ky, kz;
	int iatom;//, ipt;
	atom *patom = NULL;
	AtomBSP *pabsp = NULL;
	int ib;
	//VECTOR3 r;
	double Mn[3], dMn[3], d2Mn[3];
	double ux, uy, uz, dux, duy, duz;
	BSplineFunc<2> *bsp = vspme->bsp;
	int nb = bsp->indx + 1;

	float *md = vspme->md;

	for (iatom = na1; iatom <= na2; iatom++) {
		if (iatom >= vspme->va.n) break;
		patom = vspme->va.m + iatom;
		pabsp = vspme->absp.m + iatom;
		V3zero(patom->E_recp)
		if (var.bF) {V3zero(patom->F_recp)}
		//if (patom->c == 0) continue;
		/*
		r.v[0] = patom->r->v[0] - vspme->xl[0]; // make 0 is relative to the corner of the cell
		if (r.v[0] > vspme->xd[0]) r.v[0] -= vspme->xd[0]; 
		else if (r.v[0] < 0) r.v[0] += vspme->xd[0];
		ux = r.v[0] / vspme->md[0]; nx = int(ux);
		// if (nx >= vspme->Nx) nx = vspme->Nx - 1;
		*/
		ux = patom->u.v[0]; nx = int(ux);

		/*
		r.v[1] = patom->r->v[1] - vspme->xl[1]; // make 0 is relative to the corner of the cell
		if (r.v[1] > vspme->xd[1]) r.v[1] -= vspme->xd[1];
		else if (r.v[1] < 0) r.v[1] += vspme->xd[1];
		uy = r.v[1] / vspme->md[1]; ny = int(uy); 
		//if (ny >= vspme->Ny) nx = vspme->Ny - 1;
		*/
		uy = patom->u.v[1]; ny = int(uy);

		/*
		r.v[2] = patom->r->v[2] - vspme->xl[2]; // make 0 is relative to the corner of the cell
		if (r.v[2] > vspme->xd[2]) r.v[2] -= vspme->xd[2];
		else if (r.v[2] < 0) r.v[2] += vspme->xd[2];
		uz = r.v[2] / vspme->md[2]; nz = int(uz); 
		//if (nz >= vspme->Nz) nx = vspme->Nz - 1;
		*/
		uz = patom->u.v[2]; nz = int(uz);

		for (i = -nb; i <= nb; i++) {
			kx = nx + i; dux = ux - kx; if (dux < 0 || dux >= bsp->indx) continue;
			//ipt = bsp->get_pt(dux); 
			//Mn[0] = bsp->bsp_interpolate(ipt, dux);
			//dMn[0] = bsp->dbsp_interpolate(ipt, dux);
			//if (var.bF) d2Mn[0] = bsp->dbsp_interpolate(ipt, dux, 2);
			ib = int(dux); Mn[0] = pabsp->Mn0.m[ib]; dMn[0] = pabsp->dMn0.m[ib];
			if (var.bF) d2Mn[0] = pabsp->d2Mn0.m[ib];
			if (kx < 0) kx += vspme->Nx; else if (kx >= vspme->Nx) kx -= vspme->Nx;

			//kx = (kx == 0 ? 0 : var.N[0] - kx);

			for (j = -nb; j <= nb; j++) {
				ky = ny + j; duy = uy - ky; if (duy < 0 || duy >= bsp->indx) continue;
				//ipt = bsp->get_pt(duy); 
				//Mn[1] = bsp->bsp_interpolate(ipt, duy);
				//dMn[1] = bsp->dbsp_interpolate(ipt, duy);
				//if (var.bF) d2Mn[1] = bsp->dbsp_interpolate(ipt, duy, 2);
				ib = int(duy); Mn[1] = pabsp->Mn1.m[ib]; dMn[1] = pabsp->dMn1.m[ib];
				if (var.bF) d2Mn[1] = pabsp->d2Mn1.m[ib];
				if (ky < 0) ky += vspme->Ny; else if (ky >= vspme->Ny) ky -= vspme->Ny;

				//ky = (ky == 0 ? 0 : vspme->Ny - ky);

				for (k = -nb; k <= nb; k++) {
					kz = nz + k; duz = uz - kz;  if (duz < 0 || duz >= bsp->indx) continue;
					//ipt = bsp->get_pt(duz); 
					//Mn[2] = bsp->bsp_interpolate(ipt, duz);
					//dMn[2] = bsp->dbsp_interpolate(ipt, duz);
					//if (var.bF) d2Mn[2] = bsp->dbsp_interpolate(ipt, duz, 2);
					ib = int(duz); Mn[2] = pabsp->Mn2.m[ib]; dMn[2] = pabsp->dMn2.m[ib];
					if (var.bF) d2Mn[2] = pabsp->d2Mn2.m[ib];
					if (kz < 0) kz += vspme->Nz; else if (kz >= vspme->Nz) kz -= vspme->Nz;

					//kz = (kz == 0 ? 0 : vspme->Nz - kz);

					if (var.bE) {// accumulating E
						patom->E_recp.v[0] -= vspme->inv_V * dMn[0] / md[0] * Mn[1] * Mn[2] * var.BCQ->m[kx].m[ky].m[kz]; // accumulating
						patom->E_recp.v[1] -= vspme->inv_V * Mn[0] * dMn[1] / md[1] * Mn[2] * var.BCQ->m[kx].m[ky].m[kz]; // accumulating
						patom->E_recp.v[2] -= vspme->inv_V * Mn[0] * Mn[1] * dMn[2] / md[2] * var.BCQ->m[kx].m[ky].m[kz]; // accumulating
					}

					if (var.bF) {// accumulating F
						((FUNC)var.fFunc)(patom, vspme, Mn, dMn, d2Mn, var, kx, ky, kz);
					}
				}
			}
		}
	}
};

inline void E_recp_q(VSPME<_VQatom> *vspme, int na1, int na2, ERECP_SPME_VAR &var) {
	return E_recp_tmpl<_VQatom, VSPME<_VQatom>, ERECP_SPME_VAR>(vspme, na1, na2, var);
};

inline void E_recp_mu(VSPME<_VMUatom> *vspme, int na1, int na2, ERECP_SPME_VAR &var) {
	return E_recp_tmpl<_VMUatom, VSPME<_VMUatom>, ERECP_SPME_VAR>(vspme, na1, na2, var);
};

void cal_E_recp(VSPME<_VQatom>& vspme, bool bF, int nThreads);
void cal_E_recp(VSPME<_VMUatom>& vspme, bool bF, int nThreads);

void cal_F_recp(VSPME<_VMUatom>& vspme, int nThreads);

/*
template <class atom> void init_spme_tmpl(MMOL_MD_CELL *mdcell, VSPME<atom> &vspme, void *init_vatom_func) {
	typedef void (*FUNC)(atom*, MATOM*);
	FUNC setup_vatom = (FUNC)init_vatom_func;

	int na = 0;
	int imol = 0, ic, iatom;
	for (imol = 0; imol < mdcell->mm.n; imol++) {
		for (ic = 0; ic < mdcell->mm.m[imol].nCluster; ic++) {
			if (mdcell->mm.m[imol].cluster[ic].eNeutral) continue;
			na += mdcell->mm.m[imol].cluster[ic].nAtoms;
		}
	}
	for (imol = 0; imol < mdcell->sm.n; imol++) {
		if (mdcell->sm.m[imol].c->eNeutral) continue;
		na += mdcell->sm.m[imol].c->nAtoms;
	}
	for (imol = 0; imol < mdcell->pm.n; imol++) {
		if (mdcell->pm.m[imol].c->eNeutral) continue;
		na += 1;
	}

	vspme.va.SetArray(na);
	na = 0;
	MATOM *patom = NULL;
	for (imol = 0; imol < mdcell->mm.n; imol++) {
		for (ic = 0; ic < mdcell->mm.m[imol].nCluster; ic++) {
			if (mdcell->mm.m[imol].cluster[ic].eNeutral) continue;
			for (iatom = 0; iatom < mdcell->mm.m[imol].cluster[ic].nAtoms; iatom++) {
				patom = mdcell->mm.m[imol].cluster[ic].atom + iatom;
				setup_vatom(vspme.va.m + na, patom);
				na++;
			}
		}
	}
	for (imol = 0; imol < mdcell->sm.n; imol++) {
		if (mdcell->sm.m[imol].c->eNeutral) continue;
		for (iatom = 0; iatom < mdcell->sm.m[imol].c->nAtoms; iatom++) {
			patom = mdcell->sm.m[imol].c->atom + iatom;
			setup_vatom(vspme.va.m + na, patom);
			na++;
		}
	}
	for (imol = 0; imol < mdcell->pm.n; imol++) {
		if (mdcell->pm.m[imol].c->eNeutral) continue;
		patom = mdcell->pm.m[imol].c->atom;
		setup_vatom(vspme.va.m + na, patom);
		na++;
	}
	vspme.nCluster = vspme.va.n;
};

template <class atom> void init_spme_induced_tmpl(MMOL_MD_CELL *mdcell, VSPME<atom> &vspme, void *init_vatom_func, float *q) {
	typedef void (*FUNC)(atom*, MATOM*, float*);
	FUNC setup_vatom = (FUNC)init_vatom_func;
	*q = 0;

	int na = 0;
	int imol = 0, ic, iatom;
	BASIC_CLUSTER *pc;
	for (imol = 0; imol < mdcell->mm.n; imol++) {
		for (ic = 0; ic < mdcell->mm.m[imol].nCluster; ic++) {
			if (mdcell->mm.m[imol].cluster[ic].eNeutral) continue;
			pc = mdcell->mm.m[imol].cluster + ic;
			na += mdcell->mm.m[imol].cluster[ic].nAtoms;
		}
	}
	for (imol = 0; imol < mdcell->sm.n; imol++) {
		if (mdcell->sm.m[imol].c->eNeutral) continue;
		na += mdcell->sm.m[imol].c->nAtoms;
	}
	for (imol = 0; imol < mdcell->pm.n; imol++) {
		if (mdcell->pm.m[imol].c->eNeutral) continue;
		na += 1;
	}

	vspme.va.SetArray(na);
	na = 0;
	MATOM *patom = NULL;
	for (imol = 0; imol < mdcell->mm.n; imol++) {
		for (ic = 0; ic < mdcell->mm.m[imol].nCluster; ic++) {
			if (mdcell->mm.m[imol].cluster[ic].eNeutral) continue;
			for (iatom = 0; iatom < mdcell->mm.m[imol].cluster[ic].nAtoms; iatom++) {
				patom = mdcell->mm.m[imol].cluster[ic].atom + iatom;
				setup_vatom(vspme.va.m + na, patom, q);
				na++;
			}
		}
	}
	for (imol = 0; imol < mdcell->sm.n; imol++) {
		if (mdcell->sm.m[imol].c->eNeutral) continue;
		for (iatom = 0; iatom < mdcell->sm.m[imol].c->nAtoms; iatom++) {
			patom = mdcell->sm.m[imol].c->atom + iatom;
			setup_vatom(vspme.va.m + na, patom, q);
			na++;
		}
	}
	for (imol = 0; imol < mdcell->pm.n; imol++) {
		if (mdcell->pm.m[imol].c->eNeutral) continue;
		patom = mdcell->pm.m[imol].c->atom;
		setup_vatom(vspme.va.m + na, patom, q);
		na++;
	}
	vspme.nCluster = vspme.va.n;
};
*/

template <class atom> void init_spme_tmpl(MMOL_MD_CELL *mdcell, VSPME<atom> &vspme, void *init_vatom_func) {
	typedef void (*FUNC)(atom*, MATOM*);
	FUNC setup_vatom = (FUNC)init_vatom_func;

	int na = 0;
	int iatom;

	vspme.va.SetArray(mdcell->catom.n);
	na = 0;
	MATOM *patom = NULL;
	for (iatom = 0; iatom < mdcell->catom.n; iatom++) {
		patom = mdcell->catom.m[iatom];
		setup_vatom(vspme.va.m + na, patom);
		na++;
	}

	vspme.nCluster = vspme.va.n;
};

template <class atom> void init_spme_induced_tmpl(MMOL_MD_CELL *mdcell, VSPME<atom> &vspme, void *init_vatom_func, float *q) {
	typedef void (*FUNC)(atom*, MATOM*, float*);
	FUNC setup_vatom = (FUNC)init_vatom_func;
	*q = 0;

	int na = 0;
	int iatom;
	
	vspme.va.SetArray(mdcell->catom.n);
	na = 0;
	MATOM *patom = NULL;
	for (iatom = 0; iatom < mdcell->catom.n; iatom++) {
		patom = mdcell->catom.m[iatom];
		setup_vatom(vspme.va.m + na, patom, q);
		na++;
	}

	vspme.nCluster = vspme.va.n;
};




void init_spme(MMOL_MD_CELL *mdcell, VSPME<_VQatom> &vspme);
void init_spme(MMOL_MD_CELL *mdcell, VSPME<_VMUatom> &vspme);
void init_spme_induced(MMOL_MD_CELL *mdcell, VSPME<_VMUatom> &vspme, float *q);

void dump_F_exclude_induced_mu(VSPME<_VMUatom> &vspme, int n1, int n2, bool &status, bool bCleanBuff = true);
void dumpJob_F_exclude_induced_mu(JobIndx *job, VSPME<_VMUatom> *vspme);

void dump_EF_Q_Job(JobIndx *job, VSPME<_VQatom> *vspme, bool *bE, bool *bF);
void dump_EF_MU_Job(JobIndx *job, VSPME<_VMUatom> *vspme, bool *bE, bool *bF);


#ifdef _VMD_CELL_

template <class atom> void init_spme_tmpl_vcell(MMOL_VMD_CELL *mdcell, VSPME<atom> &vspme, void *init_vatom_func) {
	typedef void (*FUNC)(atom*, MATOM*);
	FUNC setup_vatom = (FUNC)init_vatom_func;

	int na = 0;
	int imol = 0, ic, iatom;
	for (imol = 0; imol < mdcell->mm.n; imol++) {
		for (ic = 0; ic < mdcell->mm.m[imol]->nCluster; ic++) {
			if (mdcell->mm.m[imol]->cluster[ic].eNeutral) continue;
			na += mdcell->mm.m[imol]->cluster[ic].nAtoms;
		}
	}
	for (imol = 0; imol < mdcell->sm.n; imol++) {
		if (mdcell->sm.m[imol]->c->eNeutral) continue;
		na += mdcell->sm.m[imol]->c->nAtoms;
	}
	for (imol = 0; imol < mdcell->pm.n; imol++) {
		if (mdcell->pm.m[imol]->c->eNeutral) continue;
		na += 1;
	}

	vspme.va.SetArray(na);
	na = 0;
	MATOM *patom = NULL;
	for (imol = 0; imol < mdcell->mm.n; imol++) {
		for (ic = 0; ic < mdcell->mm.m[imol]->nCluster; ic++) {
			if (mdcell->mm.m[imol]->cluster[ic].eNeutral) continue;
			for (iatom = 0; iatom < mdcell->mm.m[imol]->cluster[ic].nAtoms; iatom++) {
				patom = mdcell->mm.m[imol]->cluster[ic].atom + iatom;
				setup_vatom(vspme.va.m + na, patom);
				na++;
			}
		}
	}
	for (imol = 0; imol < mdcell->sm.n; imol++) {
		if (mdcell->sm.m[imol]->c->eNeutral) continue;
		for (iatom = 0; iatom < mdcell->sm.m[imol]->c->nAtoms; iatom++) {
			patom = mdcell->sm.m[imol]->c->atom + iatom;
			setup_vatom(vspme.va.m + na, patom);
			na++;
		}
	}
	for (imol = 0; imol < mdcell->pm.n; imol++) {
		if (mdcell->pm.m[imol]->c->eNeutral) continue;
		patom = mdcell->pm.m[imol]->c->atom;
		setup_vatom(vspme.va.m + na, patom);
		na++;
	}
	vspme.nCluster = vspme.va.n;
};

void init_spme(MMOL_VMD_CELL *mdcell, VSPME<_VQatom> &vspme);
void init_spme(MMOL_VMD_CELL *mdcell, VSPME<_VMUatom> &vspme);
#endif

void cal_StressTensor(VSPME<_VMUatom>& vspme, int nThreads);

#endif

} // end of namespace _spme_
