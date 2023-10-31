#ifdef __GPU__
namespace _gpu_md_ {
	// host class to wrap the _gpu_md_::_spme_::VSPME
	class GPU_SPME_2d : public GPU_SPME {
	public:
		_gpu_md_::_spme_::VSPME_2d vspme_host; // vspme on host, a copy of the one on device
		_gpu_md_::_spme_::VSPME_2d *vspme_dev; // vspme on device

		GPU_SPME_2d() {
			vspme_dev = NULL;
		};
		~GPU_SPME_2d() {
		  if (vspme_dev != NULL) {
		    _gpu_md_::_spme_::release_device_spme(&vspme_host, vspme_dev);
		    cudaFree(vspme_dev); vspme_dev = NULL;
		  }
		};
	};

	bool setup_gpu_spme(_gpu_md_::GPU_SPME_2d *spme_gpu, BSplineFunc<2> *bsp, bool bAutoExpandZ, float* xl, float *xd, int *dim_spme, double kappa, int natoms, int mMP);

}// _gpu_md_

namespace _atomic_cmm_ {
	void SPME_EF(EAtCELL &cell, _gpu_md_::GPU_SPME_2d &vspme, int nThreads, _spme_2d_::SPME_VAR& var, InteractRes& res, bool bInitCoordination = true);
	bool Polarize_SPME_Interact(MMOL_MD_CELL &mdcell, EAtCELL &cell, _gpu_md_::GPU_SPME_2d &vspme, int nThreads, bool bPolarize, _spme_2d_::SPME_VAR &var, InteractRes &res);
	bool Polarize_SPME_Interact2(MMOL_MD_CELL &mdcell, EAtCELL &cell, _gpu_md_::GPU_SPME_2d &vspme, int nThreads, bool bPolarize, _spme_2d_::SPME_VAR &var, InteractRes &res);
};

#endif //__GPU__
