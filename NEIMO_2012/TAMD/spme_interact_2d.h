namespace _EwaldSum_2d_ {
	class EwaldSumRecpK0Var {
	public:
		bool bEfieldOnly;
		bool bVirial;
		short iExIndDipole; // excluding induced dipole -- 0/1
		double kappa, A;

		// explicit parameters to calculate the reciprocal term with k = 0
		double PI_over_A, PI2_over_A;  // PI/A, 2PI / A
		double W; // dz * erf(kappa * dz) + 1 / (kappa * Sqrt(PI)) * exp(-(kappa * dz)^2)
		double dz, dz_abs, az, az2, az_abs, ferf, fexp;

		double inverse_kappa_sqrt_pi;  // 1 / (kappa * sqrt(PI))
		double KAPPA2_over_SQRT_PI, kappa2;  // 2 * kappa/sqrt(PI), kappa^2

		EwaldSumRecpK0Var() {bEfieldOnly = false; bVirial = true; iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);};

		void set(double kappa, double A) {
			this->kappa = kappa; this->A = A; 
			PI_over_A = PI / A; PI2_over_A = PI_over_A * 2;
			inverse_kappa_sqrt_pi = 1 / kappa * INV_SQRT_PI;

			kappa2 = kappa * kappa;
			KAPPA2_over_SQRT_PI = 2 * kappa * INV_SQRT_PI;
		};
		void init_var(double dz) {
			int i;
			this->dz = dz; dz_abs = FABS(dz);
			az = kappa * dz; az_abs = FABS(az); az2 = az * az;
			ERF(ferf, az, i) 
			W = dz * ferf;
			if (az2 < 16) {
				fexp = exp(-az2);//EXP2(fexp, az_abs, i)
				W += inverse_kappa_sqrt_pi * fexp;
			}
			else fexp = 0;
		};
		void operator = (EwaldSumRecpK0Var &var1) {
			this->bEfieldOnly = var1.bEfieldOnly;
			this->bVirial = var1.bVirial;
			this->iExIndDipole = var1.iExIndDipole;
			this->set(var1.kappa, var1.A);
			// init_var() is a dynamic part of the variables, 
			// the related variables will not be copied by operator =
		};
	};

	void EwaldSumRecpKzero(EwaldSumRecpK0Var &var, double q1, double q2, double &Ez, double &Fz, double &U);
	void EwaldSumRecpKzero(EwaldSumRecpK0Var &var, double q1, double mu1z, double q2, double mu2z, double &Ez, double &Fz, double &U);

	void EwaldSum_2D_K0_Q(ARRAY<_evatom_::_VQatom>& sys, _evatom_::_VQatom *atom, int natoms, EwaldSumRecpK0Var& var, InteractRes& res);
	void EwaldSum_2D_K0_MU(ARRAY<_evatom_::_VMUatom>& sys, _evatom_::_VMUatom *atom, int natoms, EwaldSumRecpK0Var& var, InteractRes& res);
	void EwaldSum_2D_K0_VQ(ARRAY<_evatom_::_VQatom*>& sys, _evatom_::_VQatom *atom, int natoms, EwaldSumRecpK0Var& var, InteractRes& res);
	void EwaldSum_2D_K0_VMU(ARRAY<_evatom_::_VMUatom*>& sys, _evatom_::_VMUatom *atom, int natoms, EwaldSumRecpK0Var& var, InteractRes& res);

}; // end of namespace _EwaldSum_2d_

namespace _spme_2d_ {
	class SPME_VAR : public InteractVar {
	public:
		_EwaldSum_real_::SPME_VAR spme_3d_var;  // EwaldSum for real space part in 3d and 2d have same procedures to calculate
		_EwaldSum_2d_::EwaldSumRecpK0Var K0var;

		SPME_VAR() {};
		void set_virial(bool bVirial, bool bST_diagonal) {
			spme_3d_var.bVirial = bVirial;
			spme_3d_var.bST_diagonal = bST_diagonal;
			K0var.bVirial = bVirial;
			InteractVar::bVirial = bVirial;
			InteractVar::bST_diagonal = bST_diagonal;
		};
		void operator = (SPME_VAR &var) {
			(*(InteractVar*)this) = (*(InteractVar*)(&var));
			this->spme_3d_var = var.spme_3d_var;
			this->K0var = var.K0var;
		};
	};

	void EwaldSum_2D_K0(_VSPME<_evatom_::_VQatom> &vspme, int nThreads, _EwaldSum_2d_::EwaldSumRecpK0Var &var, InteractRes &res);
	void EwaldSum_2D_K0(_VSPME<_evatom_::_VMUatom> &vspme, int nThreads, _EwaldSum_2d_::EwaldSumRecpK0Var &var, InteractRes &res);

	class EwaldSum_2D_K0_VARS_Q {
	public:
		_VSPME<_evatom_::_VQatom> *vspme;
		_EwaldSum_2d_::EwaldSumRecpK0Var* var;
		InteractRes* res;
		int nThreads;
		_THREAD_VARS thread_var;
	};

	void Asyn_EwaldSum_2D_K0_Q(_VSPME<_evatom_::_VQatom> &vspme, int nThreads, _EwaldSum_2d_::EwaldSumRecpK0Var &var, InteractRes &res, EwaldSum_2D_K0_VARS_Q &av);

	class EwaldSum_2D_K0_VARS_MU {
	public:
		_VSPME<_evatom_::_VMUatom> *vspme;
		_EwaldSum_2d_::EwaldSumRecpK0Var* var;
		InteractRes* res;
		int nThreads;
		_THREAD_VARS thread_var;
	};

	void Asyn_EwaldSum_2D_K0_MU(_VSPME<_evatom_::_VMUatom> &vspme, int nThreads, _EwaldSum_2d_::EwaldSumRecpK0Var &var, InteractRes &res, EwaldSum_2D_K0_VARS_MU &av);

	void LocalEF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res);
	void LocalEF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res);

	void SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res);
//	void SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res);
	void SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res);
//	void SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res);

	bool Polarize_SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);
//	bool Polarize_SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);
	bool Polarize_SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);
//	bool Polarize_SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);
	bool Polarize_SPME_Interact2(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, VSPME<_VMUatom>& vspme_induced, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);
	
	class E2dInteractVar : public InteractVar {
	public:
		bool bEfieldOnly;
		_EwaldSum_real_::EwaldSumRealVars1 esv1;

		E2dInteractVar() {
			bEfieldOnly = false; 
		};
		void operator = (E2dInteractVar &var) {
			this->bEfieldOnly = var.bEfieldOnly;
			(*(InteractVar*)this) = (*(InteractVar*)(&var));
			esv1.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&esv1));
		};
	};

	void SPME_Interact(MMOL_MD_CELL *mdcell, ARRAY<ImageCluster> &cimg1, ARRAY<ImageCluster> &cimg2, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME_ImageSlab<_VQatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res);
	void SPME_Interact(MMOL_MD_CELL *mdcell, ARRAY<ImageCluster> &img1, float f1, ARRAY<ImageCluster> &img2, float f2, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME_ImageSlab<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res);
	bool Polarize_SPME_Interact(MMOL_MD_CELL *mdcell, ARRAY<ImageCluster> &img1, float f1, ARRAY<ImageCluster> &img2, float f2, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME_ImageSlab<_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);

	void ExplicitEInteract_2d(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, E2dInteractVar& evar, InteractRes &res);
	void ExplicitEInteract_2d_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, E2dInteractVar& evar, InteractRes &res);
	void ExplicitEInteract_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOL_MD_CELL &mdcell, ARRAY<MATOM*> &img1, ARRAY<MATOM*> &img2, E2dInteractVar& evar, InteractRes &res);


	void CreateImageCluster(MMOL_MD_CELL &mdcell, float z, float f, ARRAY<ImageCluster> &c_img, ARRAY<MATOM*> &a_img, int mIndx_start);
	void UpdatePos_ImageCluster(ARRAY<ImageCluster> &img1, float z1, ARRAY<ImageCluster> &img2, float z2, int nThreads);
	void UpdateDipole_ImageCluster(ARRAY<ImageCluster> &img1, float f1, ARRAY<ImageCluster> &img2, float f2, int nThreads);
	void cmm_init_distribute_cluster(ImageCluster* c, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain);

} // end of namespace _EwaldSum_2d_
