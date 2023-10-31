

class MDVar_SPME {
public:
	_EwaldSum_real_::SPME_VAR spme_var;
	_spme_::VSPME<_evatom_::_VQatom> vspme_q;
	_spme_::VSPME<_evatom_::_VMUatom> vspme_mu;

	float q0; // a charge for the interactions between induced dipoles. The charge has to be set to 0, important!!
	_spme_::VSPME<_evatom_::_VMUatom> vspme_mu_induced;
	InteractRes spme_res;
};

class MDVar_LJ {
public:
	InteractVar LJ_var; 
	InteractRes LJ_res;
};

class MDVar_Torsion {
public:
	TorsionVar torsion_var[MAX_THREADS];
	HingeConstraintVar hinge_constraint_var[MAX_THREADS];
	InteractRes mires[MAX_THREADS];
};

class MDVar_ExplicitInteract {
public:
	INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation
};

template <class MD_CELL, class CMM> class MDVar_Job {
public:
	MDCELL_Job<MD_CELL> mdcellJob[MAX_THREADS];

	CMM job_cmm[MAX_THREADS];

	MDCELL_ThreadJob<MD_CELL> job_mdcell[MAX_THREADS];
	MDCELL_CMM_ThreadJob<MD_CELL, CMM> job_mdcellCMM[MAX_THREADS];

	void init_job(MD_CELL &mdcell, CMM *cmm) { // ncs is the number of atom initilized for each subcell in CMM
		// this function has to be called after cmm has been distributed of BASIC_CLUSTER
		int i, ncs = 0;
		AssignJob_MDCELL<MD_CELL>(&mdcell, this->mdcellJob, MAX_THREADS);
		if (cmm != NULL) ncs = cmm->acell.m[0]->ba.n;
		for (i = 0; i < MAX_THREADS; i++) {
			if (MAX_THREADS > 1) {
				if (ncs > 0) { // ncs > 0 means cmm use array for atoms, while not chain
					CMM_array_construct_for_storage<BASIC_CLUSTER>(this->job_cmm + i, cmm, ncs);
				}
			}

			job_mdcell[i].set_pars(i, mdcellJob + i);
			job_mdcellCMM[i].set_pars(i, mdcellJob + i);
			job_mdcellCMM[i].cmm = cmm;
			//job_mdcellCMM[i].cmm = job_cmm + i;
		}
	};
	void job_mdcell_CMM_use_internal_cmm() {
		int i;
		for (i = 0; i < MAX_THREADS; i++) job_mdcellCMM[i].cmm = job_cmm + i;
	};
	void job_mdcell_CMM_use_external_cmm(CMM* cmm) {
		int i;
		for (i = 0; i < MAX_THREADS; i++) job_mdcellCMM[i].cmm = cmm;
	}
};

class MDVar_Virial {
public:
	SVECTOR<double, 4> virial;

	bool bSaveVirial;
	ITER_SAVE virial_save;
	MDVar_Virial() {bSaveVirial = true; memset(virial.v, 0, 4 * SIZE_DOUBLE);};
};

template <int N> class Random {
public:
	_rand_::RAND_ARRAY<int> rand[N];
};


template <class MD_CELL, class BASIC_CLUSTER> class NVT_MDVar : public MDVar_LJ, public MDVar_SPME, 
public MDVar_Torsion, public MDVar_Job<MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, public MDVar_Virial, public Random<2> {
public:
	double Pint;

	NVT_MDVar() {Pint = 0;};

	void set_virial(bool bVirial, bool bST_diagonal) {
		MDVar_LJ::LJ_var.bVirial = bVirial; MDVar_LJ::LJ_var.bST_diagonal = bST_diagonal;
		MDVar_SPME::spme_var.bVirial = bVirial; MDVar_SPME::spme_var.bST_diagonal = bST_diagonal;
		MDVar_SPME::vspme_q.bST = bVirial;
		MDVar_SPME::vspme_q.bST_diagonal = bST_diagonal;
		MDVar_SPME::vspme_mu.bST = bVirial;
		MDVar_SPME::vspme_mu.bST_diagonal = bST_diagonal;
		MDVar_SPME::vspme_mu_induced.bST = bVirial;
		MDVar_SPME::vspme_mu_induced.bST_diagonal = bST_diagonal;
		int i;
		for (i = 0; i < MAX_THREADS; i++) {
			MDVar_Torsion::torsion_var[i].bVirial = bVirial;
			MDVar_Torsion::torsion_var[i].bST_diagonal = bST_diagonal;
			MDVar_Torsion::hinge_constraint_var[i].bVirial = bVirial;
			MDVar_Torsion::hinge_constraint_var[i].bST_diagonal = bST_diagonal;
		}
	};
};
