class MDVar_SPME2d {
public:
	_spme_2d_::SPME_VAR spme_var;
	_spme_2d_::VSPME<_evatom_::_VQatom> vspme_q;
	_spme_2d_::VSPME<_evatom_::_VMUatom> vspme_mu;

	float q0; // q = 0, using for the charge in vspme_mu_induced
	_spme_2d_::VSPME<_evatom_::_VMUatom> vspme_mu_induced;
	InteractRes spme_res;
};


template <class MD_CELL, class BASIC_CLUSTER> class Slab_MDVar : public MDVar_LJ, public MDVar_SPME2d, 
public MDVar_Torsion, public MDVar_Job<MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, public MDVar_Virial, public Random<2> { //without image charge
public:
	double Pint;

	Slab_MDVar() {Pint = 0; };
	void set_virial(bool bVirial, bool bST_diagonal) {
		MDVar_LJ::LJ_var.bVirial = bVirial; MDVar_LJ::LJ_var.bST_diagonal = bST_diagonal;
		MDVar_SPME2d::spme_var.set_virial(bVirial, bST_diagonal);
		//MDVar_SPME2d::spme_var.spme_3d_var.bVirial = bVirial; 
		//MDVar_SPME2d::spme_var.spme_3d_var.bST_diagonal = bST_diagonal;
		//MDVar_SPME2d::spme_var.K0var.bVirial = bVirial;
		MDVar_SPME2d::vspme_q.bST = bVirial;
		MDVar_SPME2d::vspme_q.bST_diagonal = bST_diagonal;
		MDVar_SPME2d::vspme_mu.bST = bVirial;
		MDVar_SPME2d::vspme_mu.bST_diagonal = bST_diagonal;
		MDVar_SPME2d::vspme_mu_induced.bST = bVirial;
		MDVar_SPME2d::vspme_mu_induced.bST_diagonal = bST_diagonal;
		int i;
		for (i = 0; i < MAX_THREADS; i++) {
			MDVar_Torsion::torsion_var[i].bVirial = bVirial;
			MDVar_Torsion::torsion_var[i].bST_diagonal = bST_diagonal;
			MDVar_Torsion::hinge_constraint_var[i].bVirial = bVirial;
			MDVar_Torsion::hinge_constraint_var[i].bST_diagonal = bST_diagonal;
		}
	};
};

class INTERACT_2d { // using for explicit interaction calculation
public:
	InteractVar var;
	_spme_2d_::E2dInteractVar evar;
	InteractRes LJ_res, eres;

	INTERACT_2d() {};
	~INTERACT_2d() {};
};


class MDVar_SPME_ImageSlab {
public:
	_spme_2d_::SPME_VAR spme_var;
	_spme_2d_::VSPME_ImageSlab<_evatom_::_VQatom> vspme_q;
	_spme_2d_::VSPME_ImageSlab<_evatom_::_VMUatom> vspme_mu;
	InteractRes spme_res;
};

template <class atom> class CMM_AtomJob_Var {
public:
	CMM_AtomJob<atom> cmm_atom_Job[MAX_THREADS];
	JOB_THREAD_VARS< CMM_AtomJob<atom> > cluster_distribt_JobVar[MAX_THREADS];

	CMM_AtomJob_Var() {
		//for (int i = 0; i < MAX_THREADS; i++) cluster_distribt_JobVar[i].job = cmm_atom_Job + i;
	};

	void init_job(ARRAY<atom*> &atoms) {
		AssignAtomJob_CMM<atom>(atoms, cmm_atom_Job, MAX_THREADS);
		for (int i = 0; i < MAX_THREADS; i++) cluster_distribt_JobVar[i].set_pars(i, cmm_atom_Job + i);
	};
};

template <class MD_CELL, class BASIC_CLUSTER> class ImageSlab_MDVar : public MDVar_LJ, public MDVar_SPME_ImageSlab, 
public MDVar_Torsion, public MDVar_Job<MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, public MDVar_Virial, public Random<2>,
public CMM_AtomJob_Var<BASIC_CLUSTER> { //without image charge
public:
	double Pint;

	ImageSlab_MDVar() {Pint = 0; };
	void set_virial(bool bVirial, bool bST_diagonal) {
		MDVar_LJ::LJ_var.bVirial = bVirial; MDVar_LJ::LJ_var.bST_diagonal = bST_diagonal;
		MDVar_SPME_ImageSlab::vspme_q.set_virial(bST_diagonal);
		MDVar_SPME_ImageSlab::vspme_mu.set_virial(bST_diagonal);
		MDVar_SPME_ImageSlab::spme_var.set_virial(bVirial, bST_diagonal);

		int i;
		for (i = 0; i < MAX_THREADS; i++) {
			MDVar_Torsion::torsion_var[i].bVirial = bVirial;
			MDVar_Torsion::torsion_var[i].bST_diagonal = bST_diagonal;
			MDVar_Torsion::hinge_constraint_var[i].bVirial = bVirial;
			MDVar_Torsion::hinge_constraint_var[i].bST_diagonal = bST_diagonal;
		}
	};
	void cmm_distribute_cluster_using_internal_cmm() {
		for (int i = 0; i < MAX_THREADS; i++) CMM_AtomJob_Var<BASIC_CLUSTER>::cmm_atom_Job[i].cmm =  MDVar_Job<MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::job_cmm + i;
	};
	void cmm_distribute_cluster_using_external_cmm(CMM_CELL3D<BASIC_CLUSTER>* cmm) {
		int i;
		for (i = 0; i < MAX_THREADS; i++) CMM_AtomJob_Var<BASIC_CLUSTER>::cmm_atom_Job[i].cmm = cmm;
	}
};
