namespace _gpu_md_ {
// following is copied from vspme.h
	namespace _spme_ {
#ifdef __GPU__

	class _MPROFILE {
	public:
		double *x;
		double *y0, *y1, *y2;
		int npts;
		double dx;

		__host__ __device__ void _init() {x = NULL; y0 = NULL; y1 = NULL; y2 = NULL; npts = 0; dx = 0.1;};
		__device__ _MPROFILE() {
			x = NULL; npts = 0;
			y0 = NULL; y1 = NULL; y2 = NULL;
		};
		__host__ void set_profiles_host(int n) {
			release_profiles_host(); this->npts = n; if (n == 0) return;
			cudaMalloc(&x, n * sizeof(double));
			cudaMalloc(&y0, n * sizeof(double));
			cudaMalloc(&y1, n * sizeof(double));
			cudaMalloc(&y2, n * sizeof(double));
		};
		__host__ void release_profiles_host() {
			if (x != NULL) {
				//delete[] x; x = NULL;
				cudaFree(x); x = NULL;
			}
			if (y0 != NULL) {
				//delete[] y0; y0 = NULL;
				cudaFree(y0); x = NULL;
			}
			if (y1 != NULL) {
				//delete[] y1; y1 = NULL;
				cudaFree(y1); x = NULL;
			}
			if (y2 != NULL) {
				//delete[] y2; y2 = NULL;
				cudaFree(y2); x = NULL;
			}
		};
		__device__ void release() {
			x = NULL; npts = 0;
			y0 = NULL; y1 = NULL; y2 = NULL;
			npts = 0; dx = 0;
		};
		__device__ ~_MPROFILE() {release();};
		__device__ void set(int npts, double dx, double x1 = 0) {
			release(); this->npts = npts; this->dx = dx;
			if (npts <= 0) return;
			x = new double[npts]; memset(x, 0, sizeof(double) * npts);
			int i = 0;
			for (i = 0; i < npts; i++) x[i] = x1 + dx * i;
			y0 = new double[npts]; memset(y0, 0, sizeof(double) * npts);
			y1 = new double[npts]; memset(y1, 0, sizeof(double) * npts);
			y2 = new double[npts]; memset(y2, 0, sizeof(double) * npts);
		};
		__device__ void set_y_zero(int ny) {
			int i;
			if (ny == 0) memset(y0, 0, sizeof(double) * npts);
			else if (ny == 1) memset(y1, 0, sizeof(double) * npts);
			else if (ny == 2) memset(y2, 0, sizeof(double) * npts);
			else {
				memset(y0, 0, sizeof(double) * npts);
				memset(y1, 0, sizeof(double) * npts);
				memset(y2, 0, sizeof(double) * npts);
			}
		};
		__device__ __inline__ double* y(int ndiff) {
			switch (ndiff) {
			case 0:
				return y0;
			case 1:
				return y1;
			case 2:
				return y2;
			default:
				return NULL;
			}
		};
	};

	class BSplineFunc : public _MPROFILE {
	public:
		short indx; // B-Spline function index, start from 2
		short ndiff; // nth order of derivation is required, start from 1st ...
		double xmax;
		int mpt; // the point with x = 1

		__host__ __device__ void _init() {_MPROFILE::_init(); indx = 0; xmax = 0; ndiff = 2;};
		__host__ __device__ void set_range(int indx, int npts) {
			this->indx = indx; ndiff = 2; xmax = indx;
			_MPROFILE::npts = npts; _MPROFILE::dx = xmax / (npts - 1);
			this->mpt = int(1.0 / _MPROFILE::dx + 0.1); // the point with x = 1
		};

		__device__ BSplineFunc() {indx = 0; xmax = 0; ndiff = 2;};
		__device__ ~BSplineFunc() {};
		__device__ int get_pt(double xv) {
			if (xv < 0) return -1;
			else if (xv > xmax) return _MPROFILE::npts + 1;
			return int(xv / _MPROFILE::dx + 0.1);
		};
		__device__ double bsp_interpolate(int npt, double xv) {
			if (npt < 0 || npt >= _MPROFILE::npts) return 0;
			double yv;
			double *f = _MPROFILE::y0;
			if (xv == _MPROFILE::x[npt]) return f[npt];
			if (npt == _MPROFILE::npts - 1) {
				yv = f[npt] - (f[npt-1] - f[npt]) / _MPROFILE::dx * (xv - _MPROFILE::x[npt]);
			}
			else yv = f[npt] + (f[npt+1] - f[npt]) / _MPROFILE::dx * (xv - _MPROFILE::x[npt]);
			return yv;
		};
		__device__ double dbsp_interpolate(int npt, double xv, int nth_diff = 1) {
			if (npt < 0 || npt >= _MPROFILE::npts || nth_diff > ndiff) return 0;
			double yv;
			double *f = _MPROFILE::y(nth_diff);
			if (xv == _MPROFILE::x[npt]) return f[npt];
			if (npt == _MPROFILE::npts - 1) {
				yv = f[npt] - (f[npt-1] - f[npt]) / _MPROFILE::dx * (xv - _MPROFILE::x[npt]);
			}
			else yv = f[npt] + (f[npt+1] - f[npt]) / _MPROFILE::dx * (xv - _MPROFILE::x[npt]);
			return yv;
		};
	};


	// require q, mu, u (normalized coordinate in SPME)
	// mMP -- multipole ? (0, 1, 2, ...) short
	// E_recp, F_recp
	// aE, aF -- additional E & F
	// reset local EF, on E_recp, aE, F_recp, aF
	// pointer to E, F in host memory, for output

	class _VE_CELL {
	public:
		// a set of double array with dimension nBuff, dbuff_host -- pointer at host, dbuff_dev -- device pointer
		GPU_DOUBLE *dbuff_host;
		GPU_DOUBLE *dbuff_dev;
		int ndbuff;

	public:
		int Na;
		double *r, *q_host, *mu_host, *E_host, *F_host; // mapped host memory for atomic coordination, dipole, E and F, dimension -- Na * 3
		_gpu_dstruct_::ARRAY<double> q; // dimension -- Na
		_gpu_dstruct_::ARRAY<double> u, mu, E, F; // normalized coordination, dipole, E and F, with dimension -- Na * 3
		_gpu_dstruct_::ARRAY<int> ic; // the cubic cell, with dimension -- Na * 3
		float xl[3], xd[3]; // xl -- left bottom corner; xd -- width of the cell
		float V; // volumn
		double inv_V; //invV = 4 * PI / V; !!!!!!  NOT 4 * PI / (3 * V)
		double kappa;

		__host__ __device__ void _init() {
			Na = 0; r = NULL; q_host = NULL; mu_host = NULL; E_host = NULL; F_host = NULL;
			ic._init(); q._init(); u._init(); mu._init(); E._init(); F._init();
		};
		__device__ void reset_local_E() {memset(E.m, 0, E.n * sizeof(double));};
		__device__ void reset_local_F() {memset(F.m, 0, F.n * sizeof(double));};
		/*
		__device__ void dump_E(bool bCleanBuff) {
			memcpy(E_host, E.m, E.n * sizeof(double));
			if (bCleanBuff) memset(E.m, 0, E.n * sizeof(double));
		};
		__device__ void dump_F(bool bCleanBuff) {
			memcpy(F_host, F.m, F.n * sizeof(double));
			if (bCleanBuff) memset(F.m, 0, F.n * sizeof(double));
		};
		*/

		__device__ void set_atoms(int natoms) {
			ic.set_array(natoms * 3);
			q.set_array(natoms); memset(q.m, 0, q.n * sizeof(double));
			E.set_array(natoms * 3); F.set_array(natoms * 3);
			u.set_array(natoms * 3); mu.set_array(natoms * 3); memset(mu.m, 0, mu.n * sizeof(double));
			this->Na = natoms;
		};
		__host__ void set_atoms_host(int natoms) {
			ic.set_array_host(natoms * 3);
			q.set_array_host(natoms); 
			E.set_array_host(natoms * 3); F.set_array_host(natoms * 3);
			u.set_array_host(natoms * 3); mu.set_array_host(natoms * 3); 
			this->Na = natoms;
		};
		__host__ void release_atoms_host() {
			ic.release_host(); q.release_host(); E.release_host(); F.release_host(); u.release_host(); mu.release_host();
		};
		__device__ void release() {
			ic.release(); q.release(); E.release(); F.release(); u.release(); mu.release();
		};
		__device__ _VE_CELL() {
			dbuff_host = NULL; dbuff_dev = NULL; ndbuff = 0;
			E_host = NULL; F_host = NULL;
			memset(xl, 0, 3 * sizeof(float)); memset(xd, 0, 3 * sizeof(float));
			V = 0; inv_V = 0; kappa = 0;
		};
		__device__ ~_VE_CELL() {release();};
	};

	class __VSPME :  public _VE_CELL {
	public:
		float md[3]; // full mesh-size
		__device__ void release() {
			_VE_CELL::release();
		};
		__host__ __device__ void _init() {
			_VE_CELL::_init();
		};

	};

	class _VSPME : public __VSPME {
	public:
		// host-mapped memory for CPU-FFTW
		short bUseHostFFT;
		double *ftd_host; 
		cuDoubleComplex *ftc_host;

		//_gpu_dstruct_::ARRAY<double> ftd_h; // mapping to host
		//_gpu_dstruct_::ARRAY<cuDoubleComplex> ftc_h; // mapping to host

		//ARRAY<fftw_complex> ftc;
		//fftw_plan ft_plan;
		//fftw_plan ift_plan;

		_gpu_dstruct_::ARRAY<DOUBLE> ftd;
		_gpu_dstruct_::ARRAY<COMPLEX> ftc;
		// cufftHandle will be used in host, cuFFT offer host API function only
		cufftHandle ft_plan;
		cufftHandle ift_plan;

	public:
		//double *host_map; // swap some record back to host, temporary use only

		int N; // B-spine 
		BSplineFunc bsp; // has 2nd order derivation
		int Nx, Ny, Nz, Nyz, Nv, nhx, nhy, nhz;
		//float V;
		_gpu_dstruct_::ARRAY<DOUBLE> k[3];
		_gpu_dstruct_::ARRAY<DOUBLE> b[3]; // x, y, z
		_gpu_dstruct_::FLEX_CMATRIX<DOUBLE> C;

		int Nca; // maximum number of atom inside a cubic of CELL
		_gpu_dstruct_::ARRAY<int> cmm; // the atoms located in each cubic of meshed CELL, with dimension (Nca + 1) * Nx * Ny * nz
		
		__host__ __device__ void _init() {
			bUseHostFFT = 0; ftd_host = NULL; ftc_host = NULL;
			ftd._init(); ftc._init(); //ftd_h._init(); ftc_h._init();
			k[0]._init(); k[1]._init(); k[2]._init();
			b[0]._init(); b[1]._init(); b[2]._init();
			C._init(); cmm._init(); Nca = 5;
			Cxx._init(); Cyy._init(); Czz._init(); Ctr._init();
			__VSPME::_init();
			absp._init();
		};
		__device__ __inline__ int* cmm_cubic(int ix, int iy, int iz) {
			int ic = ix * Nyz + iy * Nz + iz;
			return cmm.m + ic * (Nca + 1);
		};
		__device__ __inline__ void reset_cmm_cubic() {memset(cmm.m, 0, cmm.n * sizeof(int));};



		// each atom has dMn[3], d2Mn[3] for each atom, Mn, dMn, d2Mn have dimension N 
		// iatom, (ndiff, dim) = (Mn[0], Mn[1], Mn[2], dMn[0], dMn[1], dMn[2], d2Mn[0], d2Mn[1], d2Mn[2]), iu
		_gpu_dstruct_::ARRAY<double> absp; // dimension Na * N * 9, Mn[3], 
		__device__ __inline__ double absp_Mn(int iatom, int ndiff, int dim, int iu) {
			int ioff = iatom * 9 * N + (ndiff * 3 + dim) * N + iu;
			return absp.m[ioff];
		};
		__device__ __inline__ double* absp_eMn(int iatom, int ndiff, int dim, int iu) {
			int ioff = iatom * 9 * N + (ndiff * 3 + dim) * N + iu;
			return absp.m + ioff;
		};
		__device__ __inline__ double* absp_M(int iatom, int ndiff, int dim) {
			int ioff = iatom * 9 * N + (ndiff * 3 + dim) * N;
			return absp.m + ioff;
		};

		// Stress Tensor calculation
		// This is a cubic cell. So, the off-diagonal elements of the tensor are 0
		// We calculate the diagonal elements, xx, yy and zz
		_gpu_dstruct_::FLEX_CMATRIX<DOUBLE> Cxx, Cyy, Czz, Ctr;
		DOUBLE STxx, STyy, STzz, STtr;

		__host__ __device__ _VSPME() {
			N = 0;
			Nx = 0; Ny = 0; Nz = 0;  
			ft_plan = 0; ift_plan = 0;
			nhx = 0; nhy = 0; nhz = 0;

			STxx = 0; STyy = 0; STzz = 0;; STtr = 0;

			Nca = 5;
		};

		__device__ void release() {
			ftd.release(); ftc.release(); //ftd_h.release(); ftc_h.release();
			k[0].release(); k[1].release(); k[2].release();
			b[0].release(); b[1].release(); b[2].release();
			C.release(); cmm.release();
			Cxx.release(); Cyy.release(); Czz.release(); Ctr.release();
			bsp.release(); absp.release();
			__VSPME::release();
		};

		// cufftHandle will be used in host, cuFFT offer host API function only
		__host__ void init_fft_plans() {
			if (ft_plan != 0 || ift_plan != 0) {release_plans();}
			cufftPlan3d(&ft_plan, Nx, Ny, Nz, CUFFT_PLAN_D2Z);
			cufftSetCompatibilityMode(ft_plan, CUFFT_COMPATIBILITY_NATIVE);
			//cufftSetCompatibilityMode(ft_plan, CUFFT_COMPATIBILITY_FFTW_ALL);
			cufftPlan3d(&ift_plan, Nx, Ny, Nz, CUFFT_PLAN_Z2D);
			cufftSetCompatibilityMode(ift_plan, CUFFT_COMPATIBILITY_NATIVE);
			//cufftSetCompatibilityMode(ift_plan, CUFFT_COMPATIBILITY_FFTW_ALL);
		};
		__host__ void release_plans() {
			if (ft_plan != 0) {
				cufftDestroy(ft_plan); //fftw_destroy_plan(ft_plan); 
				ft_plan = 0;
			}
			if (ift_plan != 0) {
				cufftDestroy(ift_plan); //fftw_destroy_plan(ift_plan); 
				ift_plan = 0;
			}
		};

		__host__ __device__ ~_VSPME() {
			//release_plans(); 
			N = 0;
		};

		__host__ bool init_buff_host(int bsp_indx, int nx, int ny, int nz) {
			N = bsp_indx;
			k[0].set_array_host(nx); k[1].set_array_host(ny); k[2].set_array_host(nz);
			b[0].set_array_host(nx); b[1].set_array_host(ny); b[2].set_array_host(nz);
			this->Nx = nx; this->Ny = ny; this->Nz = nz;
			nhx = Nx / 2; nhy = Ny / 2; nhz = Nz / 2;

			bool status = true;

			C.set_host(nx, ny, nz); 
			Cxx.set_host(nx, ny, nz); Cyy.set_host(nx, ny, nz); Czz.set_host(nx, ny, nz);
			Ctr.set_host(nx, ny, nz); 

			Nyz = Ny * Nz; Nv = Nx * Nyz;
			//if (ftd_host == NULL) ftd_h.set_array_host(Nv);
			//else {ftd_h.m = ftd_host; ftd_h.n = Nv;}
			//if (ftc_host == NULL) ftc_h.set_array_host(Nv);
			//else {ftc_h.m = ftc_host; ftc_h.n = Nv;}

			ftd.set_array_host(Nv); ftc.set_array_host(Nv);

			if (!cmm.set_array_host((Nca + 1) * Nv)) status = false;

			if (!absp.set_array_host(Na * 9 * N)) status = false;
			return status;
		};
		__host__ void release_buff_host() {
			k[0].release_host(); k[1].release_host(); k[2].release_host();
			b[0].release_host(); b[1].release_host(); b[2].release_host();
			C.release_host(); 
			Cxx.release_host(); Cyy.release_host(); Czz.release_host();
			Ctr.release_host(); 

			ftd.release_host(); ftc.release_host();

			cmm.release_host();

			this->bsp.release_profiles_host();

			absp.release_host();
		};

		__host__ void cell_dim(float dx, float dy, float dz) {
			this->xd[0] = dx; this->xd[1] = dy; this->xd[2] = dz;
			this->md[0] = this->xd[0] / this->Nx; 
			this->md[1] = this->xd[1] / this->Ny;
			this->md[2] = this->xd[2] / this->Nz;

			this->V = dx * dy * dz;
			this->inv_V = 4 * PI / this->V;
		};

		__device__ void init_k(/*float dx, float dy, float dz*/) {
			int i;
			DOUBLE dk = PI2 / this->xd[0]; //dx;
			for (i = 0; i < k[0].n; i++) k[0].m[i] = dk * (i <= nhx ? i : i - Nx);
			dk = PI2 / this->xd[1]; //dy;
			for (i = 0; i < k[1].n; i++) k[1].m[i] = dk * (i <= nhy ? i : i - Ny);
			dk = PI2 / this->xd[2]; //dz;
			for (i = 0; i < k[2].n; i++) k[2].m[i] = dk * (i <= nhz ? i : i - Nz);
		};
		__device__ void init_b() {
			int i, k, n;
			double phi = 0, Mn;
			_gpu_dstruct_::complex vc;
			for (i = 0; i < Nx; i++) {
				vc.x = 0; vc.y = 0; 
				for (k = 0; k <= N - 2; k++) {
					//phi = PI2 * i * k / Nx;
					phi = PI2 * (i <= nhx ? i : i - Nx) * k / Nx;
					n = bsp.mpt * (k + 1);
					Mn = bsp.y(0)[n];
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
					n = bsp.mpt * (k + 1);
					Mn = bsp.y(0)[n];
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
					n = bsp.mpt * (k + 1);
					Mn = bsp.y(0)[n];
					vc.x += Mn * cos(phi); vc.y += Mn * sin(phi);
				}
				//phi = PI2 * (N - 1) * i / Nx;
				//b[2].m[i] = complex(cos(phi), sin(phi)) / b[2].m[i];
				b[2].m[i] = 1.0 / (vc.x * vc.x + vc.y * vc.y);
			}
		};
		__device__ void init_C() {
			// important: the summation in reciprocal space is shifted to the center of FFT space so that we can use FFT
			int ni, nj, nk;
			DOUBLE kx2, ky2, kz2, k2;
			//__VSPME<atom>::kappa = _KAPPA / ::rcut_Ewd;
			DOUBLE a2 = 4 * __VSPME::kappa * __VSPME::kappa;
			DOUBLE f = 0;
			for (ni = 0; ni < C.nx; ni++) {
				kx2 = k[0].m[ni] * k[0].m[ni];
				for (nj = 0; nj < C.ny; nj++) {
					ky2 = k[1].m[nj] * k[1].m[nj];
					for (nk = 0; nk < C.nz; nk++) {
						kz2 = k[2].m[nk] * k[2].m[nk];
						k2 = kx2 + ky2 + kz2;
						if (ni == 0 && nj == 0 && nk == 0) {
							*(C.e(ni, nj, nk)) = 0; //C.m[ni].m[nj].m[nk] = 0;
							*(Cxx.e(ni, nj, nk)) = 0; //Cxx.m[ni].m[nj].m[nk] = 0;
							*(Cyy.e(ni, nj, nk)) = 0; //Cyy.m[ni].m[nj].m[nk] = 0;
							*(Czz.e(ni, nj, nk)) = 0; //Czz.m[ni].m[nj].m[nk] = 0;
							*(Ctr.e(ni, nj, nk)) = 0; //Ctr.m[ni].m[nj].m[nk] = 0;
						}
						else {
							f = exp(-(k2) / a2) / k2;
							if (f < 1e-8) f = 0;
							*(C.e(ni, nj, nk)) = f; //C.m[ni].m[nj].m[nk] = f;
							//Cxx.m[ni].m[nj].m[nk] = f * (1 - 2 * (1 + k2 / a2) / k2 * kx2);
							*(Cxx.e(ni, nj, nk)) = f * (1 - 2 * (1 + k2 / a2) / k2 * kx2);
							//Cyy.m[ni].m[nj].m[nk] = f * (1 - 2 * (1 + k2 / a2) / k2 * ky2);
							*(Cyy.e(ni, nj, nk)) = f * (1 - 2 * (1 + k2 / a2) / k2 * ky2);
							//Czz.m[ni].m[nj].m[nk] = f * (1 - 2 * (1 + k2 / a2) / k2 * kz2);
							*(Czz.e(ni, nj, nk)) = f * (1 - 2 * (1 + k2 / a2) / k2 * kz2);
							//Ctr.m[ni].m[nj].m[nk] = f * (3 - 2 * (1 + k2 / a2));
							*(Ctr.e(ni, nj, nk)) = f * (3 - 2 * (1 + k2 / a2));
						}
					}
				}
			}
		};
	};

	class VSPME : public _VSPME {
	public:
		float *status_host; // host mapped memory, to check the history of code, which can be view on host

		short mMP; // 0 -- charge only, 1 -- with dipole
		_gpu_dstruct_::FLEX_CMATRIX<DOUBLE> Q, BCQ; // C[i][j][k] = B[i][j][k] * real C[i][j][k]
		_gpu_dstruct_::FLEX_CMATRIX<COMPLEX> FQ, FBCQ;

		_gpu_dstruct_::FLEX_CMATRIX<DOUBLE> Qbuf[Nflows];

		// for stress tensor
		//_gpu_dstruct_::FLEX_ASYM_CMATRIX<float> BTQ;
		_gpu_dstruct_::FLEX_CMATRIX<COMPLEX> FBTQ;

		__host__ bool init_buffer_host() {
			int i;
			bool status = true;
			for (i = 0; i < Nflows; i++) Qbuf[i].set_host(_VSPME::Nx, _VSPME::Ny, _VSPME::Nz);
			Q.set_host(_VSPME::Nx, _VSPME::Ny, _VSPME::Nz);
			BCQ.set_host(_VSPME::Nx, _VSPME::Ny, _VSPME::Nz);
			FQ.set_host(_VSPME::Nx, _VSPME::Ny, _VSPME::Nz);
			FBCQ.set_host(_VSPME::Nx, _VSPME::Ny, _VSPME::Nz);
			FBTQ.set_host(_VSPME::Nx, _VSPME::Ny, _VSPME::Nz);
			
			if (!_VSPME::init_buff_host(_VSPME::N, _VSPME::Nx, _VSPME::Ny, _VSPME::Nz)) status = false;

			_VSPME::set_atoms_host(_VSPME::Na);
			return status;
		};

		__host__ void release_host() {
			int i;
			for (i = 0; i < Nflows; i++) Qbuf[i].release_host();
			Q.release_host();
			BCQ.release_host();
			FQ.release_host();
			FBCQ.release_host();
			//if (_VSPME::bST) {
				//BTQ.set(Nx, Ny, Nz); 
			FBTQ.release_host();
			//}
			
			_VSPME::release_buff_host();

			_VSPME::release_atoms_host();
		};

		__host__ __device__ void _init() {
			int i;
			for (i = 0; i < Nflows; i++) Qbuf[i]._init();
			Q._init(); BCQ._init(); FQ._init(); FBCQ._init(); FBTQ._init(); _VSPME::_init();
		};
		__device__ void release() {
			int i;
			for (i = 0; i < Nflows; i++) Qbuf[i].release();
			Q.release(); BCQ.release(); FQ.release(); FBCQ.release(); FBTQ.release();
			_VSPME::release();
		};

		__device__ VSPME() {mMP = 0;};
		__device__ ~VSPME() {release();};
	};


	// host functions for calculate SPME receiprocal terms
	// status_host is host mapped memory, so that GPU can write it and host can see it
	__host__ bool setup_vspme_host(VSPME *vspme_host, VSPME **vspme_dev, bool b_init, float *status_host, int bsp_indx, double *bsp_x, double *bsp_y0, double *bsp_y1, double *bsp_y2, int bsp_npts, float *xl, float *xd, double kappa, int Nx, int Ny, int Nz, int nBSP, short mMP, double *r_hostmap, double *q_hostmap, double *mu_hostmap, double *E_hostmap, double *F_hostmap, int Na);
	__host__ void release_device_spme(VSPME *vspme_host, VSPME *vspme_dev); // release the device memory of vspme_dev, the release vspme_dev object
	__host__ void update_dipole_device_spme(VSPME *vspme_host); // host mapped memory for dipole is updated already, transfer it to the the device local memory
	__host__ void update_q_device_spme(VSPME *vspme_host); //  host mapped memory for charge is updated already, transfer it to the device local memory
	
	// to update atomic coordination, just update the host mapped memory of r_host, and then call
	__host__ void init_normalized_coordinates(VSPME *vspme_host, VSPME *vspme_dev);
	// then call this function to calculate Q of meshed cell
	__host__ void cal_Q_spme(VSPME *vspme_host, VSPME *vspme_dev, bool bQ, bool bMU);
	// then call this function to run FT / IFT transform of Q
	__host__ void vspme_FT_IFT(VSPME *vspme_host, VSPME *vspme_dev);
	
	// then we can calculate E, F, U and stress tensors
	__host__ double cal_U_recp(VSPME *vspme_host, VSPME *vspme_dev);
	__host__ void cal_ST_recp(VSPME *vspme_host, VSPME *vspme_dev, double& STxx, double& STyy, double& STzz);
	__host__ void cal_STtr_recp(VSPME *vspme_host, VSPME *vspme_dev, double& STtr);
	// when E and F are calculated, the values will be updated to the host mapped memory
	__host__ void cal_EF_spme(VSPME *vspme_host, VSPME *vspme_dev, bool bQ, bool bMU, bool bE, bool bF);



	// to use CPU-FFTW, these three functions are required between CPU FFT and IFFT, to run on GPU
	__host__ void vspme_FT_prior(VSPME *vspme_host, VSPME *vspme_dev);
	__host__ void vspme_FT_post(VSPME *vspme_host, VSPME *vspme_dev);
	__host__ void vspme_IFT_post(VSPME *vspme_host, VSPME *vspme_dev);

#endif //__GPU__

} // end of namespace _spme_
} // end of namespace _gpu_md_

