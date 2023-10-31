template <class MD_CELL, class BASIC_CLUSTER> class GC_MDVar : public MDVar_LJ, public MDVar_SPME, 
public MDVar_Torsion, public MDVar_Job<MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, public MDVar_Virial, public Random<2> {
public:
	double Pint;

	GC_MDVar() {Pint = 0;};
};

template <class MD_CELL, class BASIC_CLUSTER> class GC_MDVar1 : public GC_MDVar<MD_CELL, BASIC_CLUSTER> {
public:
	GC_InteractVar gc_InteractVar;
	_EwaldSum_real_::GC_SPME_VAR gc_EVar;
	GC_InteractRes gc_LJRes, gc_ERes;

	GC_MDVar1() {};
};


template <class MD_CELL, class BASIC_CLUSTER> class NVT_VMDVar {
public:
	InteractVar *LJ_var; 
	
	_EwaldSum_real_::SPME_VAR *spme_var;

	_spme_::VSPME<_evatom_::_VQatom> *vspme_q;
	_spme_::VSPME<_evatom_::_VMUatom> *vspme_mu;

	InteractRes *LJ_res, *spme_res;

	TorsionVar *torsion_var;  //[MAX_THREADS];
	HingeConstraintVar *hinge_constraint_var;  //[MAX_THREADS];
	InteractRes *mires;  //[MAX_THREADS];

	MDCELL_Job<MD_CELL> *mdcellJob;  //[MAX_THREADS];

	CMM_CELL3D<BASIC_CLUSTER> *job_cmm;  //[MAX_THREADS];

	MDCELL_ThreadJob<MD_CELL> *job_mdcell;  //[MAX_THREADS];
	MDCELL_CMM_ThreadJob<MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *job_mdcellCMM;  //[MAX_THREADS];

	//INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation

	_rand_::RAND_ARRAY<int> *rand0, *rand1;

	double *Pint;
	SVECTOR<double, 4> *virial;

	bool bSaveVirial;
	ITER_SAVE *virial_save;

	bool bStorageChain; // use chain to store the cluster in CMM ? false -- use array

	NVT_VMDVar() {
		bStorageChain = false;
		Pint = 0; 
		bSaveVirial = false;

		LJ_var = NULL;
		spme_var = NULL;

		vspme_q = NULL; vspme_mu = NULL;

		LJ_res = NULL; spme_res = NULL;

		torsion_var = NULL;
		hinge_constraint_var = NULL;
		mires = NULL;

		mdcellJob = NULL; job_cmm = NULL; job_mdcell = NULL; job_mdcellCMM = NULL;

		rand0 = NULL; rand1 = NULL;
		virial = NULL; virial_save = NULL;
	};
		

	void init_torsion_vars() {
		int i;
		for (i = 0; i < MAX_THREADS; i++) {
			torsion_var[i].bVirial = ::bVirial; torsion_var[i].bST_diagonal = true;
		}
	};
	void init_hinge_constraint_vars() {
		int i;
		for (i = 0; i < MAX_THREADS; i++) {
			hinge_constraint_var[i].bVirial = ::bVirial; hinge_constraint_var[i].bST_diagonal = true;
		}
	};

	void init_virial() {
		*(Pint) = 0;
		memset(virial->v, 0, 4 * SIZE_DOUBLE);
	};
	void init_virial_save(MD_SAVE *mdsave) {
		if (::bVirial) {
			virial_save->set_file(mdsave->title, "virial");
			virial_save->set_iter(mdsave->nsave, mdsave->max_save);
		}
	};

	void init_LJ_vars() {
		LJ_var->bVirial = ::bVirial; LJ_var->bST_diagonal = true;
	};
	void init_spme_vars(float V_CELL) {
		spme_var->bVirial = ::bVirial; spme_var->bST_diagonal = true; spme_var->bEfieldOnly = false;
		spme_var->esv1.bEwald = true; spme_var->esv1.init_EwaldSum(V_CELL, ::rcut_Ewd);
		spme_var->bTholePolarize = ::_Thole_polarize_::bTholePolarize;
		spme_var->tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	};

	void init_job(MD_CELL &mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm) {
		int i;
		AssignJob_MDCELL<MD_CELL>(&mdcell, mdcellJob, MAX_THREADS);
		int ncs = cmm->acell.m[0]->ba.n;
		for (i = 0; i < MAX_THREADS; i++) {
			if (MAX_THREADS > 1) {
				if (ncs > 0) {// use array while not chain
					CMM_array_construct_for_storage<BASIC_CLUSTER>(job_cmm + i, cmm, ncs);
				}
			}

			job_mdcell[i].set_pars(i, mdcellJob + i);
			job_mdcellCMM[i].set_pars(i, mdcellJob + i);
			job_mdcellCMM[i].cmm = cmm;
		}

		sprintf(errmsg, "Number of Job # : %d; working with molecules:", MAX_THREADS); show_log(errmsg, true);
		for (i = 0; i < MAX_THREADS; i++) {
			sprintf(errmsg, "%d -- MM %d, SM %d, PM %d", i, mdcellJob[i].mmIndx.n, mdcellJob[i].smIndx.n, mdcellJob[i].pmIndx.n);
			show_log(errmsg, true);
			//for (nm = 0; nm < mdcellJob[i].mmIndx.n; nm++) {
			//	sprintf(errmsg, "%d ", mdcellJob[i].mmIndx.m[nm]); show_log(errmsg, false);
			//}
			//show_log("", true);
		}
	};

	void init_SPME_SYS(MD_CELL &mdcell, bool bPolarize) {
		extern BSplineFunc<2> bsp;

		vspme_q->mMP = 0; // charge only
		vspme_q->bST = ::bVirial; vspme_q->bST_diagonal = true;
		vspme_mu->mMP = 1; // with dipole
		vspme_mu->bST = ::bVirial; vspme_mu->bST_diagonal = true;

		if (bPolarize) {
			mdcell.init_polarization_buff();
			//mdcell.init_dipole_hist();
			show_log("atomic dipole history is not initialized. Polarization is not working!", true);
			mdcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

			vspme_mu->bsp = &(::bsp);
			init_spme(&mdcell, *vspme_mu); // init the viral atoms
			vspme_mu->init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
			vspme_mu->xl[0] = -mdcell.h[0]; vspme_mu->xl[1] = -mdcell.h[1]; vspme_mu->xl[2] = -mdcell.h[2];
			vspme_mu->cell_dim(mdcell.h[0] * 2, mdcell.h[1] * 2, mdcell.h[2] * 2); // important, after init();
			vspme_mu->init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
			vspme_mu->init_b(); vspme_mu->init_C(); vspme_mu->init_vars();
		}
		else { // charge only
			mdcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

			vspme_q->bsp = &(::bsp);
			init_spme(&mdcell, *vspme_q); // init the viral atoms
			vspme_q->init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
			vspme_q->xl[0] = -mdcell.h[0]; vspme_q->xl[1] = -mdcell.h[1]; vspme_q->xl[2] = -mdcell.h[2];
			vspme_q->cell_dim(mdcell.h[0] * 2, mdcell.h[1] * 2, mdcell.h[2] * 2); // important, after init();
			vspme_q->init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
			vspme_q->init_b(); vspme_q->init_C(); vspme_q->init_vars();
		}
	};

	void UseVar(NVT_MDVar<MD_CELL, BASIC_CLUSTER> &vb) {
		this->LJ_var = &(vb.LJ_var);
		this->spme_var = &(vb.spme_var);
		this->LJ_res = &(vb.LJ_res); this->spme_res = &(vb.spme_res);
		this->vspme_q = &(vb.vspme_q); this->vspme_mu = &(vb.vspme_mu);

		this->torsion_var = vb.torsion_var;
		this->hinge_constraint_var = vb.hinge_constraint_var;
		this->mires = vb.mires;
		this->mdcellJob = vb.mdcellJob; 
		this->job_cmm = vb.job_cmm; 
		this->job_mdcell = vb.job_mdcell; 
		this->job_mdcellCMM = vb.job_mdcellCMM;
		
		this->rand0 = &(vb.rand[0]); this->rand1 = &(vb.rand[1]);
		this->virial = &(vb.virial); this->virial_save = &(vb.virial_save);

		this->Pint = &(vb.Pint);
		this->bSaveVirial = vb.bSaveVirial;
	};

};



template <class MD_CELL, class BASIC_CLUSTER> class GC_VMDVar : public NVT_VMDVar<MD_CELL, BASIC_CLUSTER> {
public:
	GC_InteractVar *gc_InteractVar;
	_EwaldSum_real_::GC_SPME_VAR *gc_EVar;
	GC_InteractRes *gc_LJRes, *gc_ERes;

	GC_VMDVar() {gc_InteractVar = NULL; gc_EVar = NULL; gc_LJRes = NULL; gc_ERes = NULL;};
	void UseVarBuffer(GC_MDVar<MD_CELL, BASIC_CLUSTER> &vb) {
		NVT_VMDVar<MD_CELL, BASIC_CLUSTER>::UseVar(*((NVT_MDVar<MD_CELL, BASIC_CLUSTER>*)(&vb)));
		this->gc_InteractVar = &(vb.gc_InteractVar);
		this->gc_EVar = &(vb.gc_EVar);
		this->gc_LJRes = &(vb.gc_LJRes);
		this->gc_ERes = &(vb.gc_ERes);
	};
};
