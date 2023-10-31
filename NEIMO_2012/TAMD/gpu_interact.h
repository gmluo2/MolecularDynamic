#ifdef __GPU__
namespace _gpu_md_ {
	// host class to wrap the _gpu_md_::_spme_::VSPME
	class AtomGPUSPME {
	public:
		double *r, *mu, *E, *F; // coordination, dipole, E & F result
		float *q; // charge
		AtomGPUSPME() {r = NULL; q = NULL; mu = NULL; E = NULL; F = NULL;};
		~AtomGPUSPME() {r = NULL; q = NULL; mu = NULL; E = NULL; F = NULL;};
	};

	class GPU_SPME {
	public:
		// multi-thread CPU FFTW could be faster than GPU-cuFFT
		ARRAY<double> ftd;
		ARRAY<fftw_complex> ftc;
		fftw_plan ft_plan;
		fftw_plan ift_plan;

		_gpu_dstruct_::_gpu_util_::HostArray<GPU_DOUBLE> dbuf;

		void init_fft_plans(int Nx, int Ny, int Nz) {
			if (ft_plan != 0 || ift_plan != 0) {release_plans();}
			int Nv = Nx * Ny * Nz;
			ftd.use_GPU_MapMem(true); ftd.set_array(Nx * Ny * Nz);
			ftc.use_GPU_MapMem(true); ftc.set_array(Nx * Ny * Nz);

			if (MAX_THREADS > 1 && Nv > 512) { // 8^3
				fftw_init_threads(); // enable multi-thread
				fftw_plan_with_nthreads(MAX_THREADS);
			}
			ft_plan = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, ftd.m, ftc.m, FFTW_ESTIMATE);
			ift_plan = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, ftc.m, ftd.m, FFTW_ESTIMATE);
		};
		void release_plans() {
			if (ft_plan != 0) {fftw_destroy_plan(ft_plan); ft_plan = 0;}
			if (ift_plan != 0) {fftw_destroy_plan(ift_plan); ift_plan = 0;}
			ftc.release(); ftd.release();
		};

	public:
		ARRAY<float> status; // to check the status inside GPU, host memory mapped to GPU device
		ARRAY<AtomGPUSPME> atom; // on host
		ARRAY<double> r, q, mu, E, F; // on host, mapped to the device

		GPU_SPME() {
			ft_plan = 0; ift_plan = 0;
			status.use_GPU_MapMem(true); status.set_array(20);
		};
		void init_memory(int natoms) {
			//q.use_GPU_MapMem(true); -- charge does not need to use mapped memory
			r.use_GPU_MapMem(true); mu.use_GPU_MapMem(true); E.use_GPU_MapMem(true); F.use_GPU_MapMem(true);
			q.set_array(natoms); r.set_array(natoms * 3); mu.set_array(natoms * 3); E.set_array(natoms * 3); F.set_array(natoms * 3);
		};
	};

	class GPU_SPME_3d : public GPU_SPME {
	public:
		_gpu_md_::_spme_::VSPME vspme_host; // vspme on host, a copy of the one on device
		_gpu_md_::_spme_::VSPME *vspme_dev; // vspme on device

		GPU_SPME_3d() {
			vspme_dev = NULL;
		};
		~GPU_SPME_3d() {
		  if (vspme_dev != NULL) {
		    _gpu_md_::_spme_::release_device_spme(&vspme_host, vspme_dev);
		    cudaFree(vspme_dev); vspme_dev = NULL;
		  }
		};
	};

	void update_atomic_coordinate(GPU_SPME *gspme, int nThreads);
	void update_atomic_dipole(GPU_SPME *gspme, int nThreads);

	void dump_atomic_E(GPU_SPME *gspme, int nThreads);
	void dump_atomic_F(GPU_SPME *gspme, int nThreads);

	void dump_subtract_atom_E(JobIndx *job, GPU_SPME *gspme);
	void dump_subtract_atom_F(JobIndx *job, GPU_SPME *gspme);
	void dump_subtract_atomic_E(GPU_SPME *gspme, int nThreads);
	void dump_subtract_atomic_F(GPU_SPME *gspme, int nThreads);

	void collect_efield_sr_job(JobIndx* job, _atomic_cmm_::EAtCELL *cell, GPU_SPME* vspme);
	void collect_eforce_sr_job(JobIndx* job, _atomic_cmm_::EAtCELL *cell, GPU_SPME* vspme);
	void reset_global_efield(JobIndx *job, GPU_SPME* vspme);
	void reset_global_force(JobIndx *job, GPU_SPME* vspme);

	void init_atom(ARRAY<MATOM*> &matom, ARRAY<AtomGPUSPME> &atom);


	bool setup_gpu_spme(_gpu_md_::GPU_SPME_3d *spme_gpu, BSplineFunc<2> *bsp, float* xl, float *xd, int *dim_spme, double kappa, int natoms, int mMP);

}// _gpu_md_

namespace _atomic_cmm_ {
	void gpu_reset_esr_neighbors(EAtCELL *acell);

	void SPME_EF(EAtCELL &cell, _gpu_md_::GPU_SPME_3d &vspme, int nThreads, SPME_VAR& var, InteractRes& res, bool bInitCoordination = true);
	bool Polarize_SPME_Interact(MMOL_MD_CELL &mdcell, EAtCELL &cell, _gpu_md_::GPU_SPME_3d &vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);
	bool Polarize_SPME_Interact2(MMOL_MD_CELL &mdcell, EAtCELL &cell, _gpu_md_::GPU_SPME_3d &vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);
};

#endif //__GPU__
