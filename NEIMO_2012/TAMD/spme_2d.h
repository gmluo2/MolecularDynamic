extern void get_atoms(MMOL_MD_CELL &mdcell, ARRAY<MATOM*> &atom);

namespace _spme_2d_ {
using namespace _evatom_;
using namespace _EwaldSum_real_;
using namespace _EwaldSum_recp_;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

template <class atom> class __VSPME :  public _VE_CELL<atom>, public __AtomBSP_BUFF  {
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
		if (ft_plan != 0 || ift_plan != 0) {
			release_plans();
		}
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

	void init(int bsp_indx, int nx, int ny, int nz, bool bInitFFTW) {
		N = bsp_indx;
		k[0].SetArray(nx); k[1].SetArray(ny); k[2].SetArray(nz);
		b[0].SetArray(nx); b[1].SetArray(ny); b[2].SetArray(nz);
		this->Nx = nx; this->Ny = ny; this->Nz = nz;
		nhx = Nx / 2; nhy = Ny / 2; nhz = Nz / 2;

		C.set(nx, ny, nz);
		//if (bST) {
		//	if (bST_diagonal) {
				Cxx.set(nx, ny, nz); Cyy.set(nx, ny, nz); Czz.set(nx, ny, nz); Ctr.set(nx, ny, nz);
		//	}
		//	else Ctr.set(nx, ny, nz);
		//}

		Nyz = Ny * Nz; Nv = Nx * Nyz;
		ftd.SetArray(Nv); ftc.SetArray(Nv);

		if (bInitFFTW) init_fft_plans();
	};
	void set_BSplineFunc(BSplineFunc<2> *bsp) {
		this->bsp = bsp; this->N = bsp->indx;
	};
	void cell_dim(float dx, float dy, float dz) {
		this->xd[0] = dx; this->xd[1] = dy; this->xd[2] = dz;
		this->md[0] = this->xd[0] / this->Nx; 
		this->md[1] = this->xd[1] / this->Ny;
		this->md[2] = this->xd[2] / this->Nz;

		this->V = dx * dy * dz;
		this->inv_V = 4 * PI / this->V;
	};

	// very important: for x-y 2d slab system, it is important to make sure z axis has enough grid,
	// because z axis is not periodical, and we have to make sure enough grids for spline function
	// So, we resize the cell, and redimension the knotes, while keeping md[3] as expected
	void set_cell(float xl, float yl, float zl, float xd, float yd, float zd, int Nx, int Ny, int Nz, bool bInitFFTW, bool bExpandSPMECellZ = true) {
		int nBSpline = this->bsp->indx;
		float mdz = zd / Nz;
		// expand the dimension of cell, along z axis
		int Next = Nz / 2 + nBSpline; // it is better to expand z length to double
		if (!bExpandSPMECellZ) Next = 0;
		Nz += Next * 2;
		zd += Next * 2 * mdz;
		zl -= Next * mdz;
		this->xl[0] = xl; this->xl[1] = yl; this->xl[2] = zl;
		init(nBSpline, Nx, Ny, Nz, bInitFFTW); // buffer will keep the same if the dimension does not change
		cell_dim(xd, yd, zd);
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
					if (ni == 0 && nj == 0/* && nk == 0*/) {
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
								Czz.m[ni].m[nj].m[nk] = f * k[2].m[nk]; 
							}
							else {
								show_log("error: diagonal elements of stress sensor have to be calculated separately!", true);
								Ctr.m[ni].m[nj].m[nk] = 0; //f * (3 - 2 * (1 + k2 / a2)); // this is WRONG !
							}
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
	FLEX_ASYM_CMATRIX<double> TQ; // defined in J. Chem. Phys., 115, 4457, (2001), Eq. (29) on page 4460 
	FLEX_ASYM_CMATRIX<complex> FTQ;
	FLEX_ASYM_CMATRIX<complex> FBTQ;

	// buffers used for multi-thread Q calculation
	FLEX_ASYM_CMATRIX<double> Qbuf[MAX_THREADS];

	void init_vars() {
		int i;
		Q.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);
		BCQ.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);
		FQ.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);

		TQ.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);
		FTQ.set(_VSPME<atom>::Nx, _VSPME<atom>::Ny, _VSPME<atom>::Nz);

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
	double cal_ST_recp_z(FLEX_ASYM_CMATRIX<double> &CT);
	void cal_FTQ();
	void cal_StressTensor_recp_0() {
		if (!_VSPME<atom>::bST) return;
		if (_VSPME<atom>::bST_diagonal) {
			_VSPME<atom>::STxx = cal_ST_recp_0(_VSPME<atom>::Cxx);
			_VSPME<atom>::STyy = cal_ST_recp_0(_VSPME<atom>::Cyy);

			_VSPME<atom>::STzz = cal_ST_recp_z(_VSPME<atom>::Czz);
			_VSPME<atom>::STtr = _VSPME<atom>::STxx + _VSPME<atom>::STyy + _VSPME<atom>::STzz;
		}
		else {
			show_log("error: diagonal elements of stress sensor have to be calculated separately!", true);
			_VSPME<atom>::STtr = 0; //cal_ST_recp_0(_VSPME<atom>::Ctr);
		}
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
			//ny = j * _VSPME<atom>::Nz;
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
				if (Q.m[i].m[j].m[k] != 0) {
					ST += Q.m[i].m[j].m[k] * _VSPME<atom>::ftd.m[n];
					//BCQ.m[i].m[j].m[k] = _VSPME<atom>::ftd.m[n]; 
				}
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
	//int inx, iny, inz;

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
				if (Q.m[i].m[j].m[k] != 0) {
					ST += Q.m[i].m[j].m[k] * _VSPME<atom>::ftd.m[n];
					//BCQ.m[i].m[j].m[k] = _VSPME<atom>::ftd.m[n]; 
				}
				n++;
			}
		}
	}

	ST *= PI2 / _VSPME<atom>::V * ::fUnit_estat_mass;

	return ST;
};
#endif

template <class atom> void VSPME<atom>::cal_FTQ() {
	// TQ ==> F[TQ]
	int i, j, k, n;
	int inx, iny, inz;

	// TQ ==>ftd
	n = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		//nx = i * Nyz;
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			//ny = j * _VSPME<atom>::Nz;
			/*
			for (k = 0; k < _VSPME<atom>::Nz; k++) {
				//ftd.m[nx + ny + k] = Q.m[i].m[j].m[k];
				_VSPME<atom>::ftd.m[n] = Q.m[i].m[j].m[k]; n++;
			}
			*/
			memcpy(_VSPME<atom>::ftd.m + n, TQ.m[i].m[j].m, _VSPME<atom>::Nz * SIZE_DOUBLE);
			n += _VSPME<atom>::Nz;
		}
	}
	fftw_execute(_VSPME<atom>::ft_plan); // ftd ==FFT==> ftc
	// ftc ==> FTQ
	n = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			/*
			for (k = 0; k <= _VSPME<atom>::nhz; k++) {
				FTQ.m[i].m[j].m[k].x = _VSPME<atom>::ftc.m[n][0];
				FTQ.m[i].m[j].m[k].y = _VSPME<atom>::ftc.m[n][1];
				n++;
			}
			*/
			memcpy(FTQ.m[i].m[j].m, _VSPME<atom>::ftc.m + n, (_VSPME<atom>::nhz + 1) * sizeof(complex));
			n += _VSPME<atom>::nhz + 1;
		}
	}
	for (i = 0; i < _VSPME<atom>::Nx; i++) {
		inx = (i == 0 ? 0 : _VSPME<atom>::Nx - i); // the inverse index 
		for (j = 0; j < _VSPME<atom>::Ny; j++) {
			iny = (j == 0 ? 0 : (_VSPME<atom>::Ny - j));
			for (k = _VSPME<atom>::nhz + 1; k < _VSPME<atom>::Nz; k++) {
				inz = _VSPME<atom>::Nz - k; // k > hz
				FTQ.m[i].m[j].m[k].x = FTQ.m[inx].m[iny].m[inz].x;
				FTQ.m[i].m[j].m[k].y = -FTQ.m[inx].m[iny].m[inz].y;
			}
		}
	}
};

template <class atom> double VSPME<atom>::cal_ST_recp_z(FLEX_ASYM_CMATRIX<double> &CST) {
	if (!_VSPME<atom>::bST) return 0;
	int i, j, k;
	//int inx, iny, inz;

	double vc, vb, vbc;
	complex *fq = NULL, *ftq = NULL;
	double ST = 0, st = 0;
	for (i = 0; i < _VSPME<atom>::Nx; i++) {for (j = 0; j < _VSPME<atom>::Ny; j++) {for (k = 0; k < _VSPME<atom>::Nz; k++) {
		if (CST.m[i].m[j].m[k] == 0) continue;
		vc = CST.m[i].m[j].m[k]; vb = _VSPME<atom>::b[0].m[i] * _VSPME<atom>::b[1].m[j] * _VSPME<atom>::b[2].m[k]; 
		fq = &(FQ.m[i].m[j].m[k]); ftq = &(FTQ.m[i].m[j].m[k]);
		vbc = vc * vb;
		st = 2 * (-fq->x * ftq->y + ftq->x * fq->y);
		ST += st * vbc;
	}}}

	ST *= PI2 / _VSPME<atom>::V * ::fUnit_estat_mass;

	return ST;
};


template <class atom> void init_normalized_coordinates_tmpl(_VSPME<atom> *vspme, int na1, int na2) {
	char msg[256] = "\0";
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

		// z axis is not periodical in 2d slab
		if (patom->u.v[2] < 0) {
			sprintf(msg, "ERROR: SPME CELL dimension in z axis is not big enough, with normalized z : %f", patom->u.v[2]); show_log(msg, true);
			patom->u.v[2] = 0;
		}
		else if (patom->u.v[2] >= vspme->Nz) {
			sprintf(msg, "ERROR: SPME CELL dimension in z axis is not big enough, with normalized z : %f", patom->u.v[2]); show_log(msg, true);
			patom->u.v[2] = 0;
		}
		//if (patom->u.v[2] < 0) patom->u.v[2] += vspme->Nz;
		//else if (patom->u.v[2] >= vspme->Nz) patom->u.v[2] -= vspme->Nz;
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
	char msg[256] = "\0";

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

		ux = patom->u.v[0]; nx = int(ux);

		uy = patom->u.v[1]; ny = int(uy);

		uz = patom->u.v[2]; nz = int(uz);
		
		for (i = -nb; i <= nb; i++) {
			kx = nx + i; 
			dux = ux - kx; if (dux < 0) break;
			ipt = bsp->get_pt(dux); Mn[0] = bsp->bsp_interpolate(ipt, dux); if (Mn[0] == 0) continue;
			ib = int(dux); pabsp->Mn0.m[ib] = Mn[0]; 
			dMn[0] = bsp->dbsp_interpolate(ipt, dux); pabsp->dMn0.m[ib] = dMn[0]; 
			if (patom->mMP > 0) {pabsp->d2Mn0.m[ib] = bsp->dbsp_interpolate(ipt, dux, 2);}
			if (kx < 0) kx += vspme->Nx; else if (kx >= vspme->Nx) kx -= vspme->Nx;
			for (j = -nb; j <= nb; j++) {
				ky = ny + j; 
				duy = uy - ky; if (duy < 0) break;
				ipt = bsp->get_pt(duy); Mn[1] = bsp->bsp_interpolate(ipt, duy); if (Mn[1] == 0) continue;
				ib = int(duy); pabsp->Mn1.m[ib] = Mn[1]; 
				dMn[1] = bsp->dbsp_interpolate(ipt, duy); pabsp->dMn1.m[ib] = dMn[1]; 
				if (patom->mMP > 0) {pabsp->d2Mn1.m[ib] = bsp->dbsp_interpolate(ipt, duy, 2);}
				if (ky < 0) ky += vspme->Ny; else if (ky >= vspme->Ny) ky -= vspme->Ny;
				for (k = -nb; k <= nb; k++) {
					kz = nz + k; 
					if (kz < 0 || kz >= vspme->Nz) {
						sprintf(msg, "ERROR: grid along z is not enough, atom %d -- uz = %f, kz = %d, range [0, %d]", iatom, uz, kz, vspme->Nz); show_log(msg, true);
						continue; // z axis is not periodical
					}
					duz = uz - kz; if (duz < 0) break;
					ipt = bsp->get_pt(duz); Mn[2] = bsp->bsp_interpolate(ipt, duz); if (Mn[2] == 0) continue;
					ib = int(duz); pabsp->Mn2.m[ib] = Mn[2]; 
					dMn[2] = bsp->dbsp_interpolate(ipt, duz); pabsp->dMn2.m[ib] = dMn[2]; 
					if (patom->mMP > 0) {pabsp->d2Mn2.m[ib] = bsp->dbsp_interpolate(ipt, duz, 2);}
					//if (kz < 0) kz += vspme->Nz; else if (kz >= vspme->Nz) kz -= vspme->Nz;

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

// assuming Mn and dMn are calculated already
template <class atom, class VSPME, class VAR, class RES> void VSPME_Q_accum_tmpl_absp(VSPME *vspme, int na1, int na2, RES& res) {
	typedef void (*FUNC)(atom*, VSPME*, double*, double*, RES&, int, int, int);
	char msg[256] = "\0";

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

		ux = patom->u.v[0]; nx = int(ux);

		uy = patom->u.v[1]; ny = int(uy);

		uz = patom->u.v[2]; nz = int(uz);
		
		for (i = -nb; i <= nb; i++) {
			kx = nx + i; 
			dux = ux - kx; if (dux < 0) break;
			//ipt = bsp->get_pt(dux); Mn[0] = bsp->bsp_interpolate(ipt, dux); if (Mn[0] == 0) continue;
			//ib = int(dux); pabsp->Mn0.m[ib] = Mn[0];
			//if (patom->mMP > 0) {dMn[0] = bsp->dbsp_interpolate(ipt, dux); pabsp->dMn0.m[ib] = dMn[0];}
			ib = int(dux); if (ib >= bsp->indx) continue;
			Mn[0] = pabsp->Mn0.m[ib]; if (Mn[0] == 0) continue; if (patom->mMP > 0) dMn[0] = pabsp->dMn0.m[ib];
			if (kx < 0) kx += vspme->Nx; else if (kx >= vspme->Nx) kx -= vspme->Nx;
			for (j = -nb; j <= nb; j++) {
				ky = ny + j; 
				duy = uy - ky; if (duy < 0) break;
				//ipt = bsp->get_pt(duy); Mn[1] = bsp->bsp_interpolate(ipt, duy); if (Mn[1] == 0) continue;
				//ib = int(duy); pabsp->Mn1.m[ib] = Mn[1];
				//if (patom->mMP > 0) {dMn[1] = bsp->dbsp_interpolate(ipt, duy); pabsp->dMn1.m[ib] = dMn[1];}
				ib = int(duy); if (ib >= bsp->indx) continue;
				Mn[1] = pabsp->Mn1.m[ib]; if (Mn[1] == 0) continue; if (patom->mMP > 0) dMn[1] = pabsp->dMn1.m[ib];
				if (ky < 0) ky += vspme->Ny; else if (ky >= vspme->Ny) ky -= vspme->Ny;
				for (k = -nb; k <= nb; k++) {
					kz = nz + k; 
					if (kz < 0 || kz >= vspme->Nz) {
						sprintf(msg, "ERROR: grid along z is not enough, atom %d -- uz = %f, kz = %d, range [0, %d]", iatom, uz, kz, vspme->Nz); show_log(msg, true);
						continue; // z axis is not periodical
					}
					duz = uz - kz; if (duz < 0) break;
					//ipt = bsp->get_pt(duz); Mn[2] = bsp->bsp_interpolate(ipt, duz); if (Mn[2] == 0) continue;
					//ib = int(duz); pabsp->Mn2.m[ib] = Mn[2];
					//if (patom->mMP > 0) {dMn[2] = bsp->dbsp_interpolate(ipt, duz); pabsp->Mn2.m[ib] = Mn[2];}
					ib = int(duz); if (ib >= bsp->indx) continue;
					Mn[2] = pabsp->Mn2.m[ib]; if (Mn[2] == 0) continue; if (patom->mMP > 0) dMn[2] = pabsp->dMn2.m[ib];
					//if (kz < 0) kz += vspme->Nz; else if (kz >= vspme->Nz) kz -= vspme->Nz;

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


/************************************************************************************************
  In J. Chem. Phys., 115, 4457 (2001), Eq. (27, 28, 29) for calculation of stress sensors
  ST_xz, ST_yz and ST_zz, T = SUM_i[q_i * z_i * exp(i * h * z_i)] is required

  for case with dipole, T = SUM_i[(q_i + p_i * delta_i) * z_i * exp(i * h * z_i)] is required
  T_i = (q_i * z_i + p_i * z_i * delta_i + pz_i) * exp(i * h * z_i)
*************************************************************************************************/
template <class atom, class VSPME, class RES> void TMU_accumulate_tmpl(atom *patom, VSPME* vspme, double* Mn, double* dMn, RES& res, int kx, int ky, int kz) {
	if (patom->mMP == 0 || patom->mu == NULL) return;

	double *mu = patom->mu->v;
	float *md = vspme->md;

	double z = patom->u.v[2] * vspme->md[2];

	// accumulating
	res.Q->m[kx].m[ky].m[kz] += z * (mu[0] * dMn[0] * Mn[1] * Mn[2] / md[0] + mu[1] * Mn[0] * dMn[1] * Mn[2] / md[1] + mu[2] * Mn[0] * Mn[1] * dMn[2] / md[2]);
	res.Q->m[kx].m[ky].m[kz] += mu[2] * Mn[0] * Mn[1] * Mn[2]; // accumulating
};

inline void TMU_accumulate(_VMUatom *atom, VSPME<_VMUatom>* vspme, double *Mn, double* dMn, RECP_SPME_RES& res, int kx, int ky, int kz) {
	return TMU_accumulate_tmpl<_VMUatom, VSPME<_VMUatom>, RECP_SPME_RES>(atom, vspme, Mn, dMn, res, kx, ky, kz);
};

template <class atom, class VSPME, class VAR, class RES> void VSPME_TQ_accum_tmpl(VSPME *vspme, int na1, int na2, RES& res) {
	typedef void (*FUNC)(atom*, VSPME*, double*, double*, RES&, int, int, int);
	char msg[256] = "\0";

	float q;
	int i, j, k, nx, ny, nz, kx, ky, kz;
	int iatom, ipt;
	atom *patom = NULL;
	AtomBSP *pabsp = NULL;
	int ib;
	//VECTOR3 r;
	double Mn[3], dMn[3], ux, uy, uz, dux, duy, duz, qz;
	BSplineFunc<2> *bsp = vspme->bsp;
	int nb = bsp->indx + 1;
	for (iatom = na1; iatom <= na2; iatom++) {
		if (iatom >= vspme->va.n) break;
		patom = vspme->va.m + iatom;
		q = *(patom->q);
		if (patom->eNeutral) continue;
		pabsp = vspme->absp.m + iatom;

		ux = patom->u.v[0]; nx = int(ux);

		uy = patom->u.v[1]; ny = int(uy);

		uz = patom->u.v[2]; nz = int(uz);

		qz = q * patom->u.v[2] * vspme->md[2];
		
		for (i = -nb; i <= nb; i++) {
			kx = nx + i; 
			dux = ux - kx; if (dux < 0) break;
			ipt = bsp->get_pt(dux); Mn[0] = bsp->bsp_interpolate(ipt, dux); if (Mn[0] == 0) continue;
			ib = int(dux); pabsp->Mn0.m[ib] = Mn[0];
			dMn[0] = bsp->dbsp_interpolate(ipt, dux); pabsp->dMn0.m[ib] = dMn[0]; 
			if (patom->mMP > 0) {pabsp->d2Mn0.m[ib] = bsp->dbsp_interpolate(ipt, dux, 2);}
			if (kx < 0) kx += vspme->Nx; else if (kx >= vspme->Nx) kx -= vspme->Nx;
			for (j = -nb; j <= nb; j++) {
				ky = ny + j; 
				duy = uy - ky; if (duy < 0) break;
				ipt = bsp->get_pt(duy); Mn[1] = bsp->bsp_interpolate(ipt, duy); if (Mn[1] == 0) continue;
				ib = int(duy); pabsp->Mn1.m[ib] = Mn[1];
				dMn[1] = bsp->dbsp_interpolate(ipt, duy); pabsp->dMn1.m[ib] = dMn[1]; 
				if (patom->mMP > 0) {pabsp->d2Mn1.m[ib] = bsp->dbsp_interpolate(ipt, duy, 2);}
				if (ky < 0) ky += vspme->Ny; else if (ky >= vspme->Ny) ky -= vspme->Ny;
				for (k = -nb; k <= nb; k++) {
					kz = nz + k; 
					if (kz < 0 || kz >= vspme->Nz) {
						sprintf(msg, "ERROR: grid along z is not enough, atom %d -- uz = %f, kz = %d, range [0, %d]", iatom, uz, kz, vspme->Nz); show_log(msg, true);
						continue; // z axis is not periodical
					}
					duz = uz - kz; if (duz < 0) break;
					ipt = bsp->get_pt(duz); Mn[2] = bsp->bsp_interpolate(ipt, duz); if (Mn[2] == 0) continue;
					ib = int(duz); pabsp->Mn2.m[ib] = Mn[2];
					dMn[2] = bsp->dbsp_interpolate(ipt, duz); pabsp->dMn2.m[ib] = dMn[2]; 
					if (patom->mMP > 0) {pabsp->d2Mn2.m[ib] = bsp->dbsp_interpolate(ipt, duz, 2);}
					//if (kz < 0) kz += vspme->Nz; else if (kz >= vspme->Nz) kz -= vspme->Nz;

					//vspme->Q->m[kx].m[ky].m[kz] += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
					res.Q->m[kx].m[ky].m[kz] += qz * Mn[0] * Mn[1] * Mn[2]; // accumulating
					if (patom->mMP > 0) {
						((FUNC)res.f_mu_accumulate)(patom, vspme, Mn, dMn, res, kx, ky, kz);
					}
				}
			}
		}
	}
};


// assuming atomic bsp buffer, Mn[3] and dMn[3] were calculated
template <class atom, class VSPME, class VAR, class RES> void VSPME_TQ_accum_tmpl_absp(VSPME *vspme, int na1, int na2, RES& res) {
	typedef void (*FUNC)(atom*, VSPME*, double*, double*, RES&, int, int, int);
	char msg[256] = "\0";

	float q;
	int i, j, k, nx, ny, nz, kx, ky, kz;
	int iatom;//, ipt;
	atom *patom = NULL;
	AtomBSP *pabsp = NULL;
	int ib;
	//VECTOR3 r;
	double Mn[3], dMn[3], ux, uy, uz, dux, duy, duz, qz;
	BSplineFunc<2> *bsp = vspme->bsp;
	int nb = bsp->indx + 1;
	for (iatom = na1; iatom <= na2; iatom++) {
		if (iatom >= vspme->va.n) break;
		patom = vspme->va.m + iatom;
		q = *(patom->q);
		if (patom->eNeutral) continue;
		pabsp = vspme->absp.m + iatom;

		ux = patom->u.v[0]; nx = int(ux);

		uy = patom->u.v[1]; ny = int(uy);

		uz = patom->u.v[2]; nz = int(uz);

		qz = q * patom->u.v[2] * vspme->md[2];
		
		for (i = -nb; i <= nb; i++) {
			kx = nx + i; 
			dux = ux - kx; if (dux < 0) break;
			//ipt = bsp->get_pt(dux); Mn[0] = bsp->bsp_interpolate(ipt, dux); if (Mn[0] == 0) continue;
			//ib = int(dux); pabsp->Mn0.m[ib] = Mn[0];
			//if (patom->mMP > 0) {dMn[0] = bsp->dbsp_interpolate(ipt, dux); pabsp->dMn0.m[ib] = dMn[0];}
			ib = int(dux); if (ib >= bsp->indx) continue;
			Mn[0] = pabsp->Mn0.m[ib]; if (Mn[0] == 0) continue; if (patom->mMP > 0) dMn[0] = pabsp->dMn0.m[ib];
			if (kx < 0) kx += vspme->Nx; else if (kx >= vspme->Nx) kx -= vspme->Nx;
			for (j = -nb; j <= nb; j++) {
				ky = ny + j; 
				duy = uy - ky; if (duy < 0) break;
				//ipt = bsp->get_pt(duy); Mn[1] = bsp->bsp_interpolate(ipt, duy); if (Mn[1] == 0) continue;
				//ib = int(duy); pabsp->Mn1.m[ib] = Mn[1];
				//if (patom->mMP > 0) {dMn[1] = bsp->dbsp_interpolate(ipt, duy); pabsp->dMn1.m[ib] = dMn[1];}
				ib = int(duy); if (ib >= bsp->indx) continue;
				Mn[1] = pabsp->Mn1.m[ib]; if (Mn[1] == 0) continue; if (patom->mMP > 0) dMn[1] = pabsp->dMn1.m[ib];
				if (ky < 0) ky += vspme->Ny; else if (ky >= vspme->Ny) ky -= vspme->Ny;
				for (k = -nb; k <= nb; k++) {
					kz = nz + k; 
					if (kz < 0 || kz >= vspme->Nz) {
						sprintf(msg, "ERROR: grid along z is not enough, atom %d -- uz = %f, kz = %d, range [0, %d]", iatom, uz, kz, vspme->Nz); show_log(msg, true);
						continue; // z axis is not periodical
					}
					duz = uz - kz; if (duz < 0) break;
					//ipt = bsp->get_pt(duz); Mn[2] = bsp->bsp_interpolate(ipt, duz); if (Mn[2] == 0) continue;
					//ib = int(duz); pabsp->Mn2.m[ib] = Mn[2];
					//if (patom->mMP > 0) {dMn[2] = bsp->dbsp_interpolate(ipt, duz); pabsp->dMn2.m[ib] = dMn[2];}
					ib = int(duz); if (ib >= bsp->indx) continue;
					Mn[2] = pabsp->Mn2.m[ib]; if (Mn[2] == 0) continue; if (patom->mMP > 0) dMn[2] = pabsp->dMn2.m[ib];
					//if (kz < 0) kz += vspme->Nz; else if (kz >= vspme->Nz) kz -= vspme->Nz;

					//vspme->Q->m[kx].m[ky].m[kz] += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
					res.Q->m[kx].m[ky].m[kz] += qz * Mn[0] * Mn[1] * Mn[2]; // accumulating
					if (patom->mMP > 0) {
						((FUNC)res.f_mu_accumulate)(patom, vspme, Mn, dMn, res, kx, ky, kz);
					}
				}
			}
		}
	}
};


inline void VSPME_TQ_accum_q(VSPME<_VQatom> *vspme, int na1, int na2, RECP_SPME_RES& res) {
	if (vspme->bReady_absp) return VSPME_TQ_accum_tmpl_absp<_VQatom, VSPME<_VQatom>, RECP_SPME_RES >(vspme, na1, na2, res);
	else return VSPME_TQ_accum_tmpl<_VQatom, VSPME<_VQatom>, RECP_SPME_RES >(vspme, na1, na2, res);
};

inline void VSPME_TQ_accum_mu(VSPME<_VMUatom> *vspme, int na1, int na2, RECP_SPME_RES& res) {
	if (vspme->bReady_absp) return VSPME_TQ_accum_tmpl_absp<_VMUatom, VSPME<_VMUatom>, RECP_SPME_RES >(vspme, na1, na2, res);
	else return VSPME_TQ_accum_tmpl<_VMUatom, VSPME<_VMUatom>, RECP_SPME_RES >(vspme, na1, na2, res);
};

void cal_TQ_spme(VSPME<_VQatom>& vspme, int nThreads);

void cal_TQ_spme(VSPME<_VMUatom>& vspme, int nThreads);




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

	char msg[256] = "\0";

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

	double bcq = 0;

	for (iatom = na1; iatom <= na2; iatom++) {
		if (iatom >= vspme->va.n) break;
		patom = vspme->va.m + iatom;
		pabsp = vspme->absp.m + iatom;
		V3zero(patom->E_recp)
		if (var.bF) {V3zero(patom->F_recp)}
		//if (patom->c == 0) continue;
		ux = patom->u.v[0]; nx = int(ux);

		uy = patom->u.v[1]; ny = int(uy);

		uz = patom->u.v[2]; nz = int(uz);

		for (i = -nb; i <= nb; i++) {
			kx = nx + i; dux = ux - kx; if (dux < 0 || dux >= bsp->indx) continue;
			//ipt = bsp->get_pt(dux); 
			//Mn[0] = bsp->bsp_interpolate(ipt, dux);
			//dMn[0] = bsp->dbsp_interpolate(ipt, dux);
			//if (var.bF) d2Mn[0] = bsp->dbsp_interpolate(ipt, dux, 2);
			ib = int(dux); Mn[0] = pabsp->Mn0.m[ib]; dMn[0] = pabsp->dMn0.m[ib]; d2Mn[0] = pabsp->d2Mn0.m[ib];
			if (kx < 0) kx += vspme->Nx; else if (kx >= vspme->Nx) kx -= vspme->Nx;

			//kx = (kx == 0 ? 0 : var.N[0] - kx);

			for (j = -nb; j <= nb; j++) {
				ky = ny + j; duy = uy - ky; if (duy < 0 || duy >= bsp->indx) continue;
				//ipt = bsp->get_pt(duy); 
				//Mn[1] = bsp->bsp_interpolate(ipt, duy);
				//dMn[1] = bsp->dbsp_interpolate(ipt, duy);
				//if (var.bF) d2Mn[1] = bsp->dbsp_interpolate(ipt, duy, 2);
				ib = int(duy); Mn[1] = pabsp->Mn1.m[ib]; dMn[1] = pabsp->dMn1.m[ib]; d2Mn[1] = pabsp->d2Mn1.m[ib];
				if (ky < 0) ky += vspme->Ny; else if (ky >= vspme->Ny) ky -= vspme->Ny;

				//ky = (ky == 0 ? 0 : vspme->Ny - ky);

				for (k = -nb; k <= nb; k++) {
					kz = nz + k; 

					if (kz < 0 || kz >= vspme->Nz) {
						sprintf(msg, "ERROR: grid along z is not enough, atom %d -- uz = %f, kz = %d, range [0, %d]", iatom, uz, kz, vspme->Nz); show_log(msg, true);
						continue; // z axis is not periodical
					}
					bcq = var.BCQ->m[kx].m[ky].m[kz];
					if (FABS(bcq) < 1e-8) continue;

					duz = uz - kz;  if (duz < 0 || duz >= bsp->indx) continue;
					//ipt = bsp->get_pt(duz); 
					//Mn[2] = bsp->bsp_interpolate(ipt, duz);
					//dMn[2] = bsp->dbsp_interpolate(ipt, duz);
					//if (var.bF) d2Mn[2] = bsp->dbsp_interpolate(ipt, duz, 2);
					ib = int(duz); Mn[2] = pabsp->Mn2.m[ib]; dMn[2] = pabsp->dMn2.m[ib]; d2Mn[2] = pabsp->d2Mn2.m[ib];
					if (kz < 0) kz += vspme->Nz; else if (kz >= vspme->Nz) kz -= vspme->Nz;

					//kz = (kz == 0 ? 0 : vspme->Nz - kz);
					if (var.bE) {// accumulating E
						patom->E_recp.v[0] -= vspme->inv_V * dMn[0] / md[0] * Mn[1] * Mn[2] * bcq; //var.BCQ->m[kx].m[ky].m[kz]; // accumulating
						patom->E_recp.v[1] -= vspme->inv_V * Mn[0] * dMn[1] / md[1] * Mn[2] * bcq; //var.BCQ->m[kx].m[ky].m[kz]; // accumulating
						patom->E_recp.v[2] -= vspme->inv_V * Mn[0] * Mn[1] * dMn[2] / md[2] * bcq; //var.BCQ->m[kx].m[ky].m[kz]; // accumulating
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


void cal_E_recp(VSPME<_VQatom>& vspme, int na1, int na2, bool bF, int nThreads);
void cal_E_recp(VSPME<_VMUatom>& vspme, int na1, int na2, bool bF, int nThreads);
void cal_F_recp(VSPME<_VMUatom>& vspme, int na1, int na2, int nThreads);

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
				// E & F for K0 term
				vspme.va.m[na].aE.set_array(1);
				vspme.va.m[na].aF.set_array(1);
				na++;
			}
		}
	}
	for (imol = 0; imol < mdcell->sm.n; imol++) {
		if (mdcell->sm.m[imol].c->eNeutral) continue;
		for (iatom = 0; iatom < mdcell->sm.m[imol].c->nAtoms; iatom++) {
			patom = mdcell->sm.m[imol].c->atom + iatom;
			setup_vatom(vspme.va.m + na, patom);
			// E & F for K0 term
			vspme.va.m[na].aE.set_array(1);
			vspme.va.m[na].aF.set_array(1);
			na++;
		}
	}
	for (imol = 0; imol < mdcell->pm.n; imol++) {
		if (mdcell->pm.m[imol].c->eNeutral) continue;
		patom = mdcell->pm.m[imol].c->atom;
		setup_vatom(vspme.va.m + na, patom);
		// E & F for K0 term
		vspme.va.m[na].aE.set_array(1);
		vspme.va.m[na].aF.set_array(1);
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
				// E & F for K0 term
				vspme.va.m[na].aE.set_array(1);
				vspme.va.m[na].aF.set_array(1);
				na++;
			}
		}
	}
	for (imol = 0; imol < mdcell->sm.n; imol++) {
		if (mdcell->sm.m[imol].c->eNeutral) continue;
		for (iatom = 0; iatom < mdcell->sm.m[imol].c->nAtoms; iatom++) {
			patom = mdcell->sm.m[imol].c->atom + iatom;
			setup_vatom(vspme.va.m + na, patom, q);
			// E & F for K0 term
			vspme.va.m[na].aE.set_array(1);
			vspme.va.m[na].aF.set_array(1);
			na++;
		}
	}
	for (imol = 0; imol < mdcell->pm.n; imol++) {
		if (mdcell->pm.m[imol].c->eNeutral) continue;
		patom = mdcell->pm.m[imol].c->atom;
		setup_vatom(vspme.va.m + na, patom, q);
		// E & F for K0 term
		vspme.va.m[na].aE.set_array(1);
		vspme.va.m[na].aF.set_array(1);
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
		// E & F for K0 term
		vspme.va.m[na].aE.set_array(1);
		vspme.va.m[na].aF.set_array(1);
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
		// E & F for K0 term
		vspme.va.m[na].aE.set_array(1);
		vspme.va.m[na].aF.set_array(1);
		na++;
	}
	vspme.nCluster = vspme.va.n;
};


void Slab_dump_EF_Q(VSPME<_VQatom> &vspme, int n1, int n2, bool& bF, bool bCleanBuff = true);
void Slab_dump_EF_MU(VSPME<_VMUatom> &vspme, int n1, int n2, bool& bF, bool bCleanBuff = true);
void Slab_dump_F_exclude_induced_mu(VSPME<_VMUatom> &vspme, int n1, int n2, bool &status, bool bCleanBuff = true);

void Slab_dumpJob_F_exclude_induced_mu(JobIndx *job, VSPME<_VMUatom> *vspme);

void Slab_dump_EF_Q_Job(JobIndx *job, VSPME<_VQatom> *vspme, bool *bE, bool *bF);
void Slab_dump_EF_MU_Job(JobIndx *job, VSPME<_VMUatom> *vspme, bool *bE, bool *bF);

void init_spme(MMOL_MD_CELL *mdcell, VSPME<_VQatom> &vspme);
void init_spme(MMOL_MD_CELL *mdcell, VSPME<_VMUatom> &vspme);
void init_spme_induced(MMOL_MD_CELL *mdcell, VSPME<_VMUatom> &vspme, float* q);

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
				// E & F for K0 term
				vspme.va.m[na].aE.set_array(1);
				vspme.va.m[na].aF.set_array(1);
				na++;
			}
		}
	}
	for (imol = 0; imol < mdcell->sm.n; imol++) {
		if (mdcell->sm.m[imol]->c->eNeutral) continue;
		for (iatom = 0; iatom < mdcell->sm.m[imol]->c->nAtoms; iatom++) {
			patom = mdcell->sm.m[imol]->c->atom + iatom;
			setup_vatom(vspme.va.m + na, patom);
			// E & F for K0 term
			vspme.va.m[na].aE.set_array(1);
			vspme.va.m[na].aF.set_array(1);
			na++;
		}
	}
	for (imol = 0; imol < mdcell->pm.n; imol++) {
		if (mdcell->pm.m[imol]->c->eNeutral) continue;
		patom = mdcell->pm.m[imol]->c->atom;
		setup_vatom(vspme.va.m + na, patom);
		// E & F for K0 term
		vspme.va.m[na].aE.set_array(1);
		vspme.va.m[na].aF.set_array(1);
		na++;
	}
	vspme.nCluster = vspme.va.n;
};

void init_spme(MMOL_VMD_CELL *mdcell, VSPME<_VQatom> &vspme);
void init_spme(MMOL_VMD_CELL *mdcell, VSPME<_VMUatom> &vspme);

#endif // of _VMD_CELL

void cal_StressTensor(VSPME<_VMUatom>& vspme, int nThreads);



//*******************************************************************
//     2-d Slab with Image Charge / Atom by One Interface only
//*******************************************************************
template <class atom> class VSPME_ImageInterface : public VSPME<atom> {
public:
	// in the slab with image charge / atom, ARRAY<atom> has atoms [0 ~ n_real - 1] are real atom
	// atoms [n_real ~ n_real + n_img - 1] are image atoms of the interface,
	// we also assume that va.n = n_real + n_img
	// Important: index of image atom is the same as that of original atoms, e.g., original atom with index i
	//            has its image atoms with indexes n_real + i
	int n_real, n_img;

	void set_atoms(int n_real, int n_img) {
		this->n_real = n_real; this->n_img = n_img;
		VSPME<atom>::va.set_array(n_real + n_img);
		VSPME<atom>::nCluster = VSPME<atom>::va.n;
	};
	VSPME_ImageInterface() {n_real = 0; n_img = 0;};
	~VSPME_ImageInterface() {};
};

//**************************************************************************************************************
//     VSPME_ImageSlab for 2-d Slab with Image Charge / Atom by two and bottom interfaces
//  VSPME calculate interaction between all the atoms. So, a slab with image charges needs 
//  to exclude the interactions between image charges(atoms), and the interactions between 
//  real atoms with their image charges(atoms) respectively with same x-y, non-periodical image-cell
//***************************************************************************************************************


//***************************************************************************************************************
//     for a slab with two interfaces, the electrostatic interaction energy is 
//  U = 0.5 * [Qr*Qr + Qr*Qimg1 + Qr*Qimg2], where Qr is the real charges defined in the slab,
//  Qimg1 and Qimg2 are image charges by interfaces 1 and 2 respectively.
//
//  For a system A with [Qr + Qimg1], U_A = 0.5 * [Qr * Qr + Qimg1 * Qimg1] + Qr * Qimg1 
//  So, Qr * Qimg1 = U_A - 0.5 * Qr * Qr - 0.5 * Qimg1 * Qimg1 == U_A - U_r - U_img1,
//  U_r and U_img1 are the energy of system with real charges, image charges 1 only respectively.
//
//  For a system B with [Qr + Qimg2], U_B = 0.5 * [Qr * Qr + Qimg2 * Qimg2] + Qr * Qimg2
//  So, Qr * Qimg2 = U_B - 0.5 * Qr * Qr - 0.5 * Qimg2 * Qimg2 == U_A - U_r - U_img2,
//  U_img2 is the energy of system with image charges 2 only.
//
//  In summary, energy of a slab with two interfaces (image charges), the electrostatic interaction
//  energy is: U = U_r + 0.5 * (U_A - U_r - U_img1) + 0.5 * (U_B - U_r - U_img2)
//               = 0.5 * (U_A + U_B - U_img1 - U_img2)

//  Because image1 and image2 are the mirrors of the charge in slab by two different interfaces respectively, 
//  from the aspect of energy, U_img1 == U_img2 in case the two interface have same dielectric constants. So, 
//     U = 0.5 * (U_A + U_B) - U_img1. Also, it is easy to obtain, U_img1 = f1 * f1 * U_r, f1 = Q_img / Q_real

//  About the force on real atoms (charges) in the slab can be derived from U. For atom i, the force on it 
//         F_i = -dU_i / dr_i = 0.5 * (dU_A_i / dr_i + dU_B_i / dr_i)
//  It is easy to find out that F_i is not relating to the systems with image charges only, img1 and img2.

//  Using Ewald Sum, the real parts and self energy of U_A + U_B, can be calculated together with same kappa, 
//  in a bigger system, real charge + image1 + image2. The interaction with image1 and image2 are counted in half,
//  while the interaction between the real charges are counted fully.

//  The reciprocal parts, k!= 0 & k=0, have to be calculated from systems A and B independently.


//********************************************************************
//   !!!! Important Modification on 02/01/2012  !!!!
//  energy of a slab with two interfaces (image charges), the electrostatic interaction
//  energy is: U = U_r + (U_A - U_r - U_img1) + (U_B - U_r - U_img2)
//               = U_A + U_B - U_r - U_img1 - U_img2
// The force calculation above is also not correct, because the coordination of image charge is associated with real charge.
// Correct way to calculate force is:
//    electrostatic interaction energy between real charges A and B with one interface is:
//      U_AB = U_A_B + U_A_Ai + U_A_Bi + U_Ai_B + U_B_Bi, where Ai and Bi are image charge of A and B respectively
//    because r_Ai = r_x_A - r_z_A, the force on A is:
//   f_A = f_A_B + f_A_Ai + f_A_Bi + f_Ai_B_T + f_Ai_A_T, 
// here f_A_ri = f_A_B + f_A_Ai + f_A_Bi is the calculated force on A within system (real + imagine)
// f_Ai_B_T is imagine of f_Ai_B, f_x_Ai_B_T = f_Ai_B and f_z_Ai_B_T = -f_z_Ai_B
// f_Ai_A_T is imagine of f_Ai_A, f_x_Ai_A_T = f_Ai_A and f_z_Ai_A_T = -f_z_Ai_A
// and, f_Ai_ri = f_Ai_A + f_Ai_B + f_Ai_Bi is the calculated force on A within system (real + imagine)
// so, f_Ai_A + f_Ai_B = f_Ai_ri - f_Ai_Bi, where f_Ai_Bi can be calculated within imagine system
// this means the final force on A should be calculated:
// f_A = f_A_ri + (f_Ai_ri - f_Ai_i)_T, here _T means imagine transform of the force on Ai, keeping xy but inverse of z

// with two interfaces, f_A = f_A_ri1 + (f_Ai_ri1 - f_Ai_i1)_T + f_A_ri2 + (f_Ai_ri2 - f_Ai_i2)_T - f_A_B_r

// same expression has to be done for E
//********************************************************************
//****************************************************************************************************************

template <class atom> class VSPME_ImageSlab {
public:
	VSPME_ImageInterface<atom> Ivspme1, Ivspme2;
	VSPME<atom> vspme; // for atoms in the slab

	VSPME<atom> img1_vspme;  // for image atoms by interface 1
	VSPME<atom> img2_vspme;  // for image atoms by interface 2
	double f1, f2; // image force coefficients for interfaces 1 and 2
	float z1, z2; // z values of interfaces 1 and 2

	VSPME_ImageSlab() {f1 = 1; f2 = 1; z1 = 0; z2 = 0;};
	void set_BSplineFunc(BSplineFunc<2> *bsp) {
		Ivspme1.set_BSplineFunc(bsp); 
		Ivspme2.set_BSplineFunc(bsp);
		img1_vspme.set_BSplineFunc(bsp);
		img2_vspme.set_BSplineFunc(bsp);
		vspme.set_BSplineFunc(bsp);
		
	};
	void set_MP(int mp) {
		Ivspme1.mMP = mp; Ivspme2.mMP = mp; 
		img1_vspme.mMP = mp; 
		img2_vspme.mMP = mp;
		vspme.mMP = mp; 
	};
	void set_cells(float xl, float yl, float zl, float xd, float yd, float zd, int Nx, int Ny, int Nz) {
		// this function must be called after set_BSplineFunc()
		float dz, z_max, z0;
		if (this->z1 < 0) { // z1 is bottom interface
			z_max = zd + zl;
			z0 = this->z1 * 2 - z_max;
			dz = z_max - z0;
		}
		else { // z1 is up interface
			z0 = zl;
			z_max = this->z1 * 2 - z0;
			dz = z_max - z0;
		}
		Ivspme1.set_cell(xl, yl, z0, xd, yd, dz, Nx, Ny, Nz + Nz, false);
		Ivspme1.init_k(); Ivspme1.init_b(); Ivspme1.init_C(); Ivspme1.init_vars();

		if (this->z2 < 0) { // z2 is bottom interface
			z_max = zd + zl;
			z0 = this->z2 * 2 - z_max;
			dz = z_max - z0;
		}
		else { // z2 is up interface
			z0 = zl;
			z_max = this->z2 * 2 - z0;
			dz = z_max - z0;
		}
		Ivspme2.set_cell(xl, yl, z0, xd, yd, dz, Nx, Ny, Nz + Nz, false);
		Ivspme2.init_k(); Ivspme2.init_b(); Ivspme2.init_C(); Ivspme2.init_vars();

		z_max = zd + zl;
		z0 = this->z1 * 2 - z_max;
		z_max = this->z1 * 2 - zl;
		dz = z_max - z0;

		img1_vspme.set_cell(xl, yl, z0, xd, yd, dz, Nx, Ny, Nz, false);
		img1_vspme.init_k(); img1_vspme.init_b(); img1_vspme.init_C(); img1_vspme.init_vars();

		z_max = zd + zl;
		z0 = this->z2 * 2 - z_max;
		z_max = this->z2 * 2 - zl;
		dz = z_max - z0;

		img2_vspme.set_cell(xl, yl, z0, xd, yd, dz, Nx, Ny, Nz, false);
		img2_vspme.init_k(); img2_vspme.init_b(); img2_vspme.init_C(); img2_vspme.init_vars();

		vspme.set_cell(xl, yl, zl, xd, yd, zd, Nx, Ny, Nz, false);
		vspme.init_k(); vspme.init_b(); vspme.init_C(); vspme.init_vars();
	};
	void set_unit_cells(float xl, float yl, float zl, float xd, float yd, float zd, int Nx, int Ny, int Nz) {
		// this function must be called after set_BSplineFunc()
		Ivspme1.set_cell(xl, yl, zl, xd, yd, zd, Nx, Ny, Nz, false);
		Ivspme1.init_k(); Ivspme1.init_b(); Ivspme1.init_C(); Ivspme1.init_vars();

		Ivspme2.set_cell(xl, yl, zl, xd, yd, zd, Nx, Ny, Nz, false);
		Ivspme2.init_k(); Ivspme2.init_b(); Ivspme2.init_C(); Ivspme2.init_vars();

		img1_vspme.set_cell(xl, yl, zl, xd, yd, zd, Nx, Ny, Nz, false);
		img1_vspme.init_k(); img1_vspme.init_b(); img1_vspme.init_C(); img1_vspme.init_vars();

		img2_vspme.set_cell(xl, yl, zl, xd, yd, zd, Nx, Ny, Nz, false);
		img2_vspme.init_k(); img2_vspme.init_b(); img2_vspme.init_C(); img2_vspme.init_vars();

		vspme.set_cell(xl, yl, zl, xd, yd, zd, Nx, Ny, Nz, false);
		vspme.init_k(); vspme.init_b(); vspme.init_C(); vspme.init_vars();
	};
	void release_fftw_plan() {
		Ivspme1.release_plans(); Ivspme2.release_plans(); vspme.release_plans();
		img1_vspme.release_plans();
		img2_vspme.release_plans();
	};
	// E & F calculation has to be done for both real and imagine atoms, and we assume they are done already
	//f_A = f_A_ri1 + (f_Ai_ri1 - f_Ai_i1)_T + f_A_ri2 + (f_Ai_ri2 - f_Ai_i2)_T - f_A_B_r
	void dump_EF(bool bF = true, bool bCleanBuff = true) {
		int i, i1, i2, j;
		double F, v;
		for (i = 0; i < Ivspme1.n_real; i++) {
			i1 = i + Ivspme1.n_real; i2 = i + Ivspme2.n_real;
			Ivspme1.va.m[i].E->v[0] += Ivspme1.va.m[i].E_recp.v[0] + Ivspme2.va.m[i].E_recp.v[0] - vspme.va.m[i].E_recp.v[0];
			Ivspme1.va.m[i].E->v[0] += Ivspme1.va.m[i1].E_recp.v[0] + Ivspme2.va.m[i2].E_recp.v[0]; // on image charge
			Ivspme1.va.m[i].E->v[0] -= img1_vspme.va.m[i].E_recp.v[0] + img2_vspme.va.m[i].E_recp.v[0];
			v = 0;
			for (j = 0; j < Ivspme1.va.m[i].aE.n; j++) {
				v += Ivspme1.va.m[i].aE.m[j].v[0] + Ivspme2.va.m[i].aE.m[j].v[0] - vspme.va.m[i].aE.m[j].v[0];
				v += Ivspme1.va.m[i1].aE.m[j].v[0] + Ivspme2.va.m[i2].aE.m[j].v[0]; // on image charge
				v -= img1_vspme.va.m[i].aE.m[j].v[0] + img2_vspme.va.m[i].aE.m[j].v[0];
			}
			Ivspme1.va.m[i].E->v[0] += v;

			Ivspme1.va.m[i].E->v[1] += Ivspme1.va.m[i].E_recp.v[1] + Ivspme2.va.m[i].E_recp.v[1] - vspme.va.m[i].E_recp.v[1];
			Ivspme1.va.m[i].E->v[1] += Ivspme1.va.m[i1].E_recp.v[1] + Ivspme2.va.m[i2].E_recp.v[1]; // on image charge
			Ivspme1.va.m[i].E->v[1] -= img1_vspme.va.m[i].E_recp.v[1] + img2_vspme.va.m[i].E_recp.v[1];
			v = 0;
			for (j = 0; j < Ivspme1.va.m[i].aE.n; j++) {
				v += Ivspme1.va.m[i].aE.m[j].v[1] + Ivspme2.va.m[i].aE.m[j].v[1] - vspme.va.m[i].aE.m[j].v[1];
				v += Ivspme1.va.m[i1].aE.m[j].v[1] + Ivspme2.va.m[i2].aE.m[j].v[1]; // on image charge
				v -= img1_vspme.va.m[i].aE.m[j].v[1] + img2_vspme.va.m[i].aE.m[j].v[1];
			}
			Ivspme1.va.m[i].E->v[1] += v;

			Ivspme1.va.m[i].E->v[2] += Ivspme1.va.m[i].E_recp.v[2] + Ivspme2.va.m[i].E_recp.v[2] - vspme.va.m[i].E_recp.v[2];
			Ivspme1.va.m[i].E->v[2] -= Ivspme1.va.m[i1].E_recp.v[2] + Ivspme2.va.m[i2].E_recp.v[2]; // on image charge
			Ivspme1.va.m[i].E->v[2] += img1_vspme.va.m[i].E_recp.v[2] + img2_vspme.va.m[i].E_recp.v[2];
			v = 0;
			for (j = 0; j < Ivspme1.va.m[i].aE.n; j++) {
				v += Ivspme1.va.m[i].aE.m[j].v[2] + Ivspme2.va.m[i].aE.m[j].v[2] - vspme.va.m[i].aE.m[j].v[2];
				v -= Ivspme1.va.m[i1].aE.m[j].v[2] + Ivspme2.va.m[i2].aE.m[j].v[2]; // on image charge
				v += img1_vspme.va.m[i].aE.m[j].v[2] + img2_vspme.va.m[i].aE.m[j].v[2];
			}
			Ivspme1.va.m[i].E->v[2] += v;

			if (bF) {
				F = Ivspme1.va.m[i].F_recp.v[0] + Ivspme2.va.m[i].F_recp.v[0] - vspme.va.m[i].F_recp.v[0];
				F += Ivspme1.va.m[i1].F_recp.v[0] + Ivspme2.va.m[i2].F_recp.v[0]; // on image charge
				F -= img1_vspme.va.m[i].F_recp.v[0] + img2_vspme.va.m[i].F_recp.v[0];
				for (j = 0; j < Ivspme1.va.m[i].aF.n; j++) {
					F += Ivspme1.va.m[i].aF.m[j].v[0] + Ivspme2.va.m[i].aF.m[j].v[0] - vspme.va.m[i].aF.m[j].v[0];
					F += Ivspme1.va.m[i1].aF.m[j].v[0] + Ivspme2.va.m[i2].aF.m[j].v[0]; // on image charge
					F -= img1_vspme.va.m[i].aF.m[j].v[0] + img2_vspme.va.m[i].aF.m[j].v[0];
				}
				Ivspme1.va.m[i].F->v[0] += F * fUnit_estat_mass;

				F = Ivspme1.va.m[i].F_recp.v[1] + Ivspme2.va.m[i].F_recp.v[1] - vspme.va.m[i].F_recp.v[1];
				F += Ivspme1.va.m[i1].F_recp.v[1] + Ivspme2.va.m[i2].F_recp.v[1]; // on image charge
				F -= img1_vspme.va.m[i].F_recp.v[1] + img2_vspme.va.m[i].F_recp.v[1];
				for (j = 0; j < Ivspme1.va.m[i].aF.n; j++) {
					F += Ivspme1.va.m[i].aF.m[j].v[1] + Ivspme2.va.m[i].aF.m[j].v[1] - vspme.va.m[i].aF.m[j].v[1];
					F += Ivspme1.va.m[i1].aF.m[j].v[1] + Ivspme2.va.m[i2].aF.m[j].v[1]; // on image charge
					F -= img1_vspme.va.m[i].aF.m[j].v[1] + img2_vspme.va.m[i].aF.m[j].v[1];
				}
				Ivspme1.va.m[i].F->v[1] += F * fUnit_estat_mass;

				F = Ivspme1.va.m[i].F_recp.v[2] + Ivspme2.va.m[i].F_recp.v[2] - vspme.va.m[i].F_recp.v[2];
				F -= Ivspme1.va.m[i1].F_recp.v[2] + Ivspme2.va.m[i2].F_recp.v[2]; // on image charge
				F += img1_vspme.va.m[i].F_recp.v[2] + img2_vspme.va.m[i].F_recp.v[2];
				for (j = 0; j < Ivspme1.va.m[i].aF.n; j++) {
					F += Ivspme1.va.m[i].aF.m[j].v[2] + Ivspme2.va.m[i].aF.m[j].v[2] - vspme.va.m[i].aF.m[j].v[2];
					F -= Ivspme1.va.m[i1].aF.m[j].v[2] + Ivspme2.va.m[i2].aF.m[j].v[2]; // on image charge
					F += img1_vspme.va.m[i].aF.m[j].v[2] + img2_vspme.va.m[i].aF.m[j].v[2];
				}
				Ivspme1.va.m[i].F->v[2] += F * fUnit_estat_mass;
			}
			if (bCleanBuff) {
				Ivspme1.va.m[i].reset_local_EF(bF);
				Ivspme2.va.m[i].reset_local_EF(bF);
				img1_vspme.va.m[i].reset_local_EF(bF);
				img2_vspme.va.m[i].reset_local_EF(bF);
				vspme.va.m[i].reset_local_EF(bF);
			}
		}
	};
	void set_virial(bool bST_diagonal) {
		Ivspme1.bST = true;
		Ivspme1.bST_diagonal = bST_diagonal;
		Ivspme2.bST = true;
		Ivspme2.bST_diagonal = bST_diagonal;
		img1_vspme.bST = true;
		img1_vspme.bST_diagonal = bST_diagonal;
		img2_vspme.bST = true;
		img2_vspme.bST_diagonal = bST_diagonal;
		vspme.bST = true;
		vspme.bST_diagonal = bST_diagonal;
	};
};

void ImageSlab_dump_EF_Q(VSPME_ImageSlab<_VQatom> &vspme, int n1, int n2, bool& bF, bool bCleanBuff = true);
void ImageSlab_dump_EF_MU(VSPME_ImageSlab<_VMUatom> &vspme, int n1, int n2, bool& bF, bool bCleanBuff = true);

template <class atom> void init_spme_tmpl_ImageSlab(ARRAY<MATOM*> &atom_r, ARRAY<MATOM*> &img1, ARRAY<MATOM*> &img2, VSPME_ImageSlab<atom> &vspme_slab, void *init_vatom_func) {
	typedef void (*FUNC)(atom*, MATOM*);
	FUNC setup_vatom = (FUNC)init_vatom_func;

	int ia, iatom;
	
	MATOM *patom = NULL;
	// main vspme has all the atoms
	vspme_slab.Ivspme1.set_atoms(atom_r.n, img1.n);
	vspme_slab.Ivspme2.set_atoms(atom_r.n, img2.n);

	for (iatom = 0; iatom < atom_r.n; iatom++) {
		patom = atom_r.m[iatom];
		setup_vatom(vspme_slab.Ivspme1.va.m + iatom, patom);
		setup_vatom(vspme_slab.Ivspme2.va.m + iatom, patom);
		// E & F for K0 term
		vspme_slab.Ivspme1.va.m[iatom].aE.set_array(1);
		vspme_slab.Ivspme1.va.m[iatom].aF.set_array(1);
		vspme_slab.Ivspme2.va.m[iatom].aE.set_array(1);
		vspme_slab.Ivspme2.va.m[iatom].aF.set_array(1);
	}
	for (iatom = 0; iatom < img1.n; iatom++) {
		patom = img1.m[iatom];
		setup_vatom(vspme_slab.Ivspme1.va.m + atom_r.n + iatom, patom);
		// E & F for K0 term
		vspme_slab.Ivspme1.va.m[atom_r.n + iatom].aE.set_array(1);
		vspme_slab.Ivspme1.va.m[atom_r.n + iatom].aF.set_array(1);
	}
	for (iatom = 0; iatom < img2.n; iatom++) {
		patom = img2.m[iatom];
		setup_vatom(vspme_slab.Ivspme2.va.m + atom_r.n + iatom, patom);
		// E & F for K0 term
		vspme_slab.Ivspme2.va.m[atom_r.n + iatom].aE.set_array(1);
		vspme_slab.Ivspme2.va.m[atom_r.n + iatom].aF.set_array(1);
	}

	// img1_vspme, img2_vspme
	vspme_slab.img1_vspme.set_atoms(img1.n);
	vspme_slab.img2_vspme.set_atoms(img2.n);

	for (iatom = 0; iatom < img1.n; iatom++) {
		patom = img1.m[iatom];
		setup_vatom(vspme_slab.img1_vspme.va.m + iatom, patom);
		// E & F for K0 term
		vspme_slab.img1_vspme.va.m[iatom].aE.set_array(1);
		vspme_slab.img1_vspme.va.m[iatom].aF.set_array(1);
	}
	for (iatom = 0; iatom < img2.n; iatom++) {
		patom = img2.m[iatom];
		setup_vatom(vspme_slab.img2_vspme.va.m + iatom, patom);
		// E & F for K0 term
		vspme_slab.img2_vspme.va.m[iatom].aE.set_array(1);
		vspme_slab.img2_vspme.va.m[iatom].aF.set_array(1);
	}

	// vspme
	vspme_slab.vspme.set_atoms(atom_r.n);

	for (iatom = 0; iatom < atom_r.n; iatom++) {
		patom = atom_r.m[iatom];
		setup_vatom(vspme_slab.vspme.va.m + iatom, patom);
		// E & F for K0 term
		vspme_slab.vspme.va.m[iatom].aE.set_array(1);
		vspme_slab.vspme.va.m[iatom].aF.set_array(1);
	}
};



void init_spme(ARRAY<MATOM*> &atomr_r, ARRAY<MATOM*> &img1, ARRAY<MATOM*> &img2, VSPME_ImageSlab<_VQatom> &vspme);
void init_spme(ARRAY<MATOM*> &atomr_r, ARRAY<MATOM*> &img1, ARRAY<MATOM*> &img2, VSPME_ImageSlab<_VMUatom> &vspme);

#endif // _IGNORE_ELECTROSTATIC_INTERACTION_
} // end of namespace _spme_
