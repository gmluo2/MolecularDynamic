class _MD_CELL_DIM {
public:
	float h[3]; // half width / height / depth of the MD cell in x, y, z directions
	bool period;

	void set_dim(float dx, float dy, float dz) { h[0] = dx; h[1] = dy; h[2] = dz;};

	_EwaldSum_real_::SURF_DIPOLE mu_surf;

	_MD_CELL_DIM() {period = true; h[0] = 0; h[1] = 0; h[2] = 0;};
};

class MDCELL_HDYN : public HOOVER_DYN {
public:
	ARRAY<int> mmIndx, smIndx, pmIndx;
	bool mol_coupled() {
		if (mmIndx.n > 0) return true;
		if (smIndx.n > 0) return true;
		if (pmIndx.n > 0) return true;
		return false;
	};
};

template <unsigned short N> class MM_VELOCITY {
public:
	ARRAY< VECTOR<N> > cv;
};

template <unsigned short N> class MDCELL_DYN_BUFFER {
public:
	ARRAY< MM_VELOCITY<N> > mmv1, mmv2;
	ARRAY<double> mm_dvmax, mm_dwmax;
	ARRAY<bool> mm_status;
	ARRAY<bool> mm_sc_status; // speed convergent ?

	ARRAY<double> sm_dvmax, sm_dwmax;
	ARRAY<bool> sm_status;
	ARRAY<bool> sm_sc_status; // speed convergent ?

	ARRAY<double> pm_dvmax;
	ARRAY<bool> pm_status;
	ARRAY<bool> pm_sc_status; // speed convergent ?

	void reset_verlet_status(bool status);
	void reset_sc_status(bool status);
};

template <unsigned short N> void MDCELL_DYN_BUFFER<N>::reset_verlet_status(bool status) {
	int i;
	for (i = 0; i < mm_status.n; i++) mm_status.m[i] = status;
	for (i = 0; i < sm_status.n; i++) sm_status.m[i] = status;
	for (i = 0; i < pm_status.n; i++) pm_status.m[i] = status;
};

template <unsigned short N> void MDCELL_DYN_BUFFER<N>::reset_sc_status(bool status) {
	int i;
	for (i = 0; i < mm_status.n; i++) mm_sc_status.m[i] = status;
	for (i = 0; i < sm_status.n; i++) sm_sc_status.m[i] = status;
	for (i = 0; i < pm_status.n; i++) pm_sc_status.m[i] = status;
};

template <class MOL> class MOL_POL {
public:
	MOL *m;
	ARRAY<VECTOR3> mu;
	MOL_POL() {m = NULL;};
	~MOL_POL() {m = NULL;};
};

class MMOL_MD_CELL : public _MD_CELL_DIM {
public:
	ARRAY<MMOLECULE> mm;
	ARRAY<LFVERLET_MM_MD> lfmd;

	ARRAY<SMOLECULE> sm;
	ARRAY<LFVERLET_SM_MD> lfsmd;

	ARRAY<PMOLECULE> pm;
	ARRAY<LFVERLET_PM_MD> lfpmd;

	ARRAY<MMOLECULE*> bmm; // big macro-molecules
	ARRAY<MMOLECULE*> smm; // small macro-molecules

	ARRAY<MATOM*> atom;
	ARRAY<BASIC_CLUSTER*> bcluster; // all the clusters in the MMOL_MD_CELL

	ARRAY<MATOM*> catom; // charged atoms in the ccluster
	ARRAY<BASIC_CLUSTER*> ccluster; // clusters with charge (non-neutral) in the MMOL_MD_CELL

	ARRAY<MATOM*> pol_atom; // polarizable atom
	ARRAY<DIPOLE_HIST> mu_hist;

	PARS_DATABASE< COMM_PARS<50> > mTypeIndx_db; // database of molecule type-index, which is CELL-dependent
	//PARS_DATABASE< COMM_PARS<50> > cTypeIndx_db; // database of cluster type-index, which is CELL-dependent
	// it is defined in global variable, cluster_db;

	MD_SAVE *msave;
	long LOOPS;

	int NF;
	float Ethermal;
	ARRAY<MDCELL_HDYN> hdyn; 
	MDCELL_DYN_BUFFER<6> mvb;

	ARRAY< MOL_POL<MMOLECULE> > ind_mu_mm; // for the polarized dipole
	ARRAY< MOL_POL<SMOLECULE> > ind_mu_sm; // for the polarized dipole
	ARRAY< MOL_POL<PMOLECULE> > ind_mu_pm; // for the polarized dipole

	// Monitor the spatial momentum of the whole cell
	bool bHDynMC_T; // applying NHC on mass center of the whole cell, translation only ?
	float M; // total mass the cell
	VECTOR3 rc; // mass center
	VECTOR3 Vmc; // mass-center velocity of the whole cell, the rotation of the whole cell
	HOOVER_DYN mcHDyn;
	float Ethermal_mct; // mass center translation and rotational thermal kinetic energy
	char ix_mct, iy_mct, iz_mct; // the way to calculate mass center
	float Mthermal_mc; // in Nose-Hoover chain, temperature of mass center kinetic energy has to be normalized to dT/TB
	                    // normally dT/TB = Ek / Ethermal - 1. However, for mass center, its kinetic energy can be set to 0
	                    // In this case, dT / TB = (Ek - Ethermal) / kT / Mthermal_mct
	                    // effectively it changes Nose-Hoover constant tao to sqrt(Mthermal_mct) * tao

	float Ekt_mc; // translation and rotation energy of the whole cell

	VECTOR3 fc, ac; // total force of all species in the cell, and the acceleration of the whole cell

	bool eNeutral;
	bool bHingeConstraint;

	void setup_mcHDyn(float tao, float dt, float max_ksi) {
		if (bHDynMC_T) {
			if (mcHDyn.nhc == NULL) mcHDyn.nhc = new NHC<MNHC>;
			mcHDyn.set_vars(tao, dt, max_ksi);
		}
	};
	void release_mcHDyn() {
		if (mcHDyn.nhc != NULL) {delete mcHDyn.nhc; mcHDyn.nhc = NULL;}
	};

	MMOL_MD_CELL() {
		NF = 0; Ethermal = 0; bHDynMC_T = false; eNeutral = false;
		Ekt_mc = 0; msave = NULL;
		bHingeConstraint = false;
		Ethermal_mct = 0; Mthermal_mc = 1; ix_mct = 1; iy_mct = 1; iz_mct = 1;
	};
	void check_eneutral();
	void SetMol(int n_mm, int n_sm, int n_pm = 0) {
		mm.SetArray(n_mm); lfmd.SetArray(n_mm); 
		sm.SetArray(n_sm); lfsmd.SetArray(n_sm);
		pm.SetArray(n_pm); lfpmd.SetArray(n_pm);
	};
	void SetSave(MD_SAVE *msave) {
		int i = 0; 
		for (i = 0; i < lfmd.n; i++) lfmd.m[i].mmsave.mdsave = msave;
		for (i = 0; i < lfsmd.n; i++) lfsmd.m[i].smsave.mdsave = msave;
		for (i = 0; i < lfpmd.n; i++) lfpmd.m[i].pmsave.mdsave = msave;
		this->msave = msave;
	};
	void setMD(long nloops) {
		int i;
		for (i = 0; i < lfmd.n; i++) lfmd.m[i].LOOPS = nloops;
		for (i = 0; i < lfsmd.n; i++) lfsmd.m[i].LOOPS = nloops;
		for (i = 0; i < lfpmd.n; i++) lfpmd.m[i].LOOPS = nloops;
		this->LOOPS = nloops;
	};
	void set_md_timestep(float dt) {
		int i;
		for (i = 0; i < lfmd.n; i++) lfmd.m[i].set_dt(dt);
		for (i = 0; i < lfsmd.n; i++) lfsmd.m[i].set_dt(dt);
		for (i = 0; i < lfpmd.n; i++) lfpmd.m[i].set_dt(dt);
	};
	bool MDSave(long nloop, float Ep, bool bEachMolDyn = true, bool bAsync = false);
	void SetMDMemory() {
		int i = 0; for (i = 0; i < lfmd.n; i++) {lfmd.m[i].resetClusters(mm.m[i].nCluster);}
	};
	int max_ncluster() {
		int i = 0, nmax = 0;
		for (i = 0; i < mm.n; i++) if (nmax < mm.m[i].nCluster) nmax = mm.m[i].nCluster;
		return nmax;
	};
	void release() {
		int i;
		for (i = 0; i < bmm.n; i++) bmm.m[i] = NULL;
		for (i = 0; i < smm.n; i++) smm.m[i] = NULL;
		bmm.release(); smm.release();
		mm.release(); lfmd.release();
		sm.release(); lfsmd.release();
		pm.release(); lfpmd.release();
		mTypeIndx_db.release();
	};
	
	void release_hdyn() {
		int i;
		if (hdyn.m != NULL) {
			for (i = 0; i < hdyn.n; i++) {
				if (hdyn.m[i].nhc != NULL) {delete hdyn.m[i].nhc; hdyn.m[i].nhc = NULL;}
			}
			hdyn.release();
		}
	};
	void set_hdyn(int nHDyn) {
		int i;
		release_hdyn(); if (nHDyn <= 0) return;
		hdyn.SetArray(nHDyn);
		for (i = 0; i < hdyn.n; i++) {
			hdyn.m[i].nhc = new NHC<MNHC>;
			hdyn.m[i].hID = -2;
		}
	};
	HOOVER_DYN* get_HDyn(short hID) {
		HOOVER_DYN *pdyn = NULL;
		int i;
		for (i = 0; i < hdyn.n; i++) {if (hdyn.m[i].hID == hID) {pdyn = hdyn.m + i; break;}}
		if (pdyn == NULL) {
			for (i = 0; i < hdyn.n; i++) {if (hdyn.m[i].hID == -1) {pdyn = hdyn.m + i; break;}}
		}
		return pdyn;
	};
	void init_nhc_mols();  // to group the molecules into each NHC
	void check_freedom();  // freedoms for each NHC

	~MMOL_MD_CELL() {
		release(); release_hdyn(); release_mcHDyn();
	};
	void setMolIndx() {
		int n = 0, i = 0;
		COMM_PARS<50> *mpar = NULL;
		this->mTypeIndx_db.release(); //this->cTypeIndx_db.release();
		for (i = 0; i < this->mm.n; i++) {
			mpar = this->mTypeIndx_db.attach(this->mm.m[i].mol); if (mpar != NULL) mpar->uid = mpar->indx;
			if (mpar == NULL) mpar = this->mTypeIndx_db.search_par(this->mm.m[i].mol);
			if (mpar != NULL) this->mm.m[i].mTypeIndx = mpar->uid;
/*
			for (ic = 0; ic < mm.m[i].nCluster; ic++) {
				mpar = this->cTypeIndx_db.attach(this->mm.m[i].cluster[ic].cID); 
				if (mpar == NULL) this->cTypeIndx_db.search_par_uid(this->mm.m[i].cluster[ic].cID);
				if (mpar != NULL) this->mm.m[i].cluster[ic].cTypeIndx = mpar->indx;
			}
*/
			this->mm.m[i].setMolIndx(n); n++;
		}
		for (i = 0; i < this->sm.n; i++) {
			mpar = this->mTypeIndx_db.attach(this->sm.m[i].mol); if (mpar != NULL) mpar->uid = mpar->indx;
			if (mpar == NULL) mpar = this->mTypeIndx_db.search_par(this->sm.m[i].mol);
			if (mpar != NULL) this->sm.m[i].mTypeIndx = mpar->uid;
/*
			mpar = this->cTypeIndx_db.attach(this->sm.m[i].c->cID); 
			if (mpar == NULL) this->cTypeIndx_db.search_par_uid(this->sm.m[i].c->cID);
			if (mpar != NULL) this->sm.m[i].c->cTypeIndx = mpar->indx;
*/
			this->sm.m[i].mIndx = n; this->sm.m[i].c->mIndx = n; n++;
		}
		for (i = 0; i < this->pm.n; i++) {
			mpar = this->mTypeIndx_db.attach(this->pm.m[i].mol); if (mpar != NULL) mpar->uid = mpar->indx;
			if (mpar == NULL) mpar = this->mTypeIndx_db.search_par(this->pm.m[i].mol);
			if (mpar != NULL) this->pm.m[i].mTypeIndx = mpar->uid;
/*
			mpar = this->cTypeIndx_db.attach(this->pm.m[i].c->cID); 
			if (mpar == NULL) this->cTypeIndx_db.search_par_uid(this->pm.m[i].c->cID);
			if (mpar != NULL) this->pm.m[i].c->cTypeIndx = mpar->indx;
*/
			this->pm.m[i].mIndx = n; this->pm.m[i].c->mIndx = n; n++;
		}
	};
	void Init_MD() {
		int i = 0;
		for (i = 0; i < mm.n; i++) lfmd.m[i].Init(mm.m[i]);
		for (i = 0; i < sm.n; i++) lfsmd.m[i].Init(sm.m[i]);
		for (i = 0; i < pm.n; i++) lfpmd.m[i].Init(pm.m[i]);
	};
	void CalibrateTorsionAngle() {
		int i = 0;
		for (i = 0; i < mm.n; i++) lfmd.m[i].CalibrateTorsionAngle(mm.m[i]);
	};
	void init_big_small_mm() {
		int i, nbm = 0, nsm = 0, ibm = 0, ism = 0;
		for (i = 0; i < bmm.n; i++) bmm.m[i] = NULL;
		for (i = 0; i < smm.n; i++) smm.m[i] = NULL;
		bmm.release(); smm.release();
		for (i = 0; i < mm.n; i++) {
			if (mm.m[i].nCluster >= NCLUSTERS_BIG_MM) nbm++;
			else nsm++;
		}
		bmm.SetArray(nbm); smm.SetArray(nsm);
		ibm = 0; ism = 0;
		for (i = 0; i < mm.n; i++) {
			if (mm.m[i].nCluster >= NCLUSTERS_BIG_MM) {
				bmm.m[ibm] = mm.m + i; ibm++;
			}
			else {
				smm.m[ism] = mm.m + i; ism++;
			}
		}
	};
	void init_mdcell_velocity_memory();
	void init_polarization_buff();
	void set_atom_multipole(int mMP);
	void TotalMass() {
		int ic, im = 0;
		M = 0;
		for (im = 0; im < mm.n; im++) {
			for (ic = 0; ic < mm.m[im].nCluster; ic++) M += mm.m[im].cluster[ic].M;
		}
		for (im = 0; im < sm.n; im++) M += sm.m[im].c->M;
		for (im = 0; im < pm.n; im++) M += pm.m[im].c->M;
	};
	void init_groups() {
		extern void get_atoms(MMOL_MD_CELL &mdcell, ARRAY<MATOM*> &atom);
		extern void get_bclusters(MMOL_MD_CELL &mdcell, ARRAY<BASIC_CLUSTER*> &bcluster);
		get_atoms(*this, atom);
		get_bclusters(*this, bcluster);
	};

	void init_dipole_hist();
};

int nAtoms(MMOL_MD_CELL &mdcell);

void Backup_induced_dipole_hist(MMOL_MD_CELL &mdcell, int nMDThreads);
void Guess_induced_dipole(MMOL_MD_CELL &mdcell, int nMDThreads);

void MM_maxdiff_backup_polarization(MOL_POL<MMOLECULE> *mp, int nmol, double &dmu_max);
void SM_maxdiff_backup_polarization(MOL_POL<SMOLECULE> *mp, int nmol, double &dmu_max);
void PM_maxdiff_backup_polarization(MOL_POL<PMOLECULE> *mp, int nmol, double &dmu_max);

class MD_THREAD_PAR : public _THREAD_VARS {
public:
	bool* bHooverDyn;
	int* nconvg_loops;
	bool* res;
	double* Ek, *Ep;

	int md_loop;
	int n; // number of SM

	char msg_show[5120], msg_log[5120];

	void reset() {
		if (bHooverDyn != NULL) {delete[] bHooverDyn; bHooverDyn = NULL;}
		if (nconvg_loops != NULL) {delete[] nconvg_loops; nconvg_loops = NULL;}
		if (res != NULL) {delete[] res; res = NULL;}
		if (Ek != NULL) {delete[] Ek; Ek = NULL;}
		if (Ep != NULL) {delete[] Ep; Ep = NULL;}
		strcpy(msg_show, "\0"); strcpy(msg_log, "\0");
	};
	void set(int nM, int iThread) {
		reset();
		_THREAD_VARS::set_pars(iThread);
		this->md_loop = 0; this->iThread = iThread;
		this->n = nM; if (nM <= 0) return;
		bHooverDyn = new bool[nM]; res = new bool[nM]; nconvg_loops = new int[nM];
		Ek = new double[nM]; this->Ep = new double[nM];
		int i = 0; for (i = 0; i < nM; i++) {Ek[i] = 0; this->Ep[i] = 0; bHooverDyn[i] = true; nconvg_loops[i] = 0;}
	};
	void set_loop(int md_loop, bool HooverDyn, double Ep) {
		this->md_loop = md_loop;
		int i = 0; for (i = 0; i < n; i++) {bHooverDyn[i] = HooverDyn; this->Ep[i] = Ep;}
	};
	MD_THREAD_PAR() {n = 0; bHooverDyn = NULL; nconvg_loops = NULL; res = NULL; Ek = NULL; Ep = NULL; strcpy(msg_show, "\0"); strcpy(msg_log, "\0");};
	~MD_THREAD_PAR() {reset();};
};

class SM_MD_THREAD_PAR : public MD_THREAD_PAR {
public:
	SMOLECULE* sm;
	LFVERLET_SM_MD *lfsmdpar;

	SM_MD_THREAD_PAR() {sm = NULL; lfsmdpar = NULL;};
	void set(SMOLECULE *sm, LFVERLET_SM_MD* lfsmdpar, int nSM, int iThread) {
		this->sm = sm; this->lfsmdpar = lfsmdpar;
		if (nSM <= 0) return;
		MD_THREAD_PAR::set(nSM, iThread);
	};
	~SM_MD_THREAD_PAR() {MD_THREAD_PAR::reset();};
};

class MM_MD_THREAD_PAR : public MD_THREAD_PAR {
public:
	MMOLECULE *mm;
	LFVERLET_MM_MD *lfmdpar;

	VECTOR<6> *vect1, *vect2;
	double *dv1, *dv2;
	bool local_relax;

	MM_MD_THREAD_PAR() {mm = NULL; lfmdpar = NULL; vect1 = NULL; vect2 = NULL; dv1 = NULL; dv2 = NULL; local_relax = false;};
	void set(MMOLECULE *mm, LFVERLET_MM_MD* lfmdpar, bool local_relax, int nMM, int iThread) {
		this->mm = mm; this->lfmdpar = lfmdpar; this->local_relax = local_relax;
		if (nMM <= 0) return;
		MD_THREAD_PAR::set(nMM, iThread);
	};
	void reset_buffer_memory() {
		RELEASE_ARRAY(vect1); RELEASE_ARRAY(vect2); RELEASE_ARRAY(dv1); RELEASE_ARRAY(dv2);
	};
	void reset() {MD_THREAD_PAR::reset(); reset_buffer_memory();};
	void set_buffer_memory(int nset) {
		reset_buffer_memory();
		this->vect1 = new VECTOR<6>[nset]; this->vect2 = new VECTOR<6>[nset];
		this->dv1 = new double[nset]; this->dv2 = new double[nset];
	};
	~MM_MD_THREAD_PAR() {reset();};
};

#if _SYS_ == _WINDOWS_SYS_
int SM_LFVERLET_MD_thread_func(LPVOID *par);
int MM_LFVERLET_MD_thread_func(LPVOID *par);
#elif _SYS_ == _LINUX_SYS_
void* SM_LFVERLET_MD_thread_func(void *par);
void* MM_LFVERLET_MD_thread_func(void *par);
#endif

void MDCELL_mass_center(MMOL_MD_CELL &mmcell);
void set_free_cell_mass_center(MMOL_MD_CELL &mmcell, VECTOR3 &O);
void cluster_check_periodic_cell(MMOL_MD_CELL &mmcell);
void molecule_check_periodic_cell(MMOL_MD_CELL &mmcell);

void InitClusterForceInCMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell);

// calculate NEIMO_MOMENT_MATRIX for multi-macromolecules
void MultiMM_NEIMO_MOMENT(MMOLECULE *mm, int nmm);
void MultiMMP_NEIMO_MOMENT(MMOLECULE** mm, int nmm);



void MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, int md_mode);


/***********************************************************************************************
	MULTI-THREAD Operator on MDCELL, OVER ALL MOLECULES
***********************************************************************************************/

template <class MDCELL> class MDCELL_Job {
public:
	MDCELL *mdcell;
	ARRAY<int> mmIndx, smIndx, pmIndx;

	MDCELL_Job() {
		this->mdcell = NULL;
	};
	void set_pars(MDCELL* mdcell) {
		this->mdcell = mdcell;
	};
	~MDCELL_Job() {
		this->mdcell = NULL;
	};
};

template <class MDCELL> void AssignJob_MDCELL(MDCELL *mdcell, MDCELL_Job<MDCELL> *job, int nJobs) {
	//MDCELL_THREAD_VARS<MDCELL> *thread_var = new MDCELL_THREAD_VARS<MDCELL>[nThreads];
	int i = 0, imol;
	int n1, n2, nstep, n;
	for (i = 0; i < nJobs; i++) {
		job[i].mmIndx.release(); 
		job[i].smIndx.release(); 
		job[i].pmIndx.release(); 
		job[i].set_pars(mdcell);
	}
	if (mdcell->mm.n > 0) {
		n1 = 0; nstep = mdcell->mm.n / nJobs;
		if (nstep == 0 || (nstep * nJobs + 1) < mdcell->mm.n) nstep += 1;
		for (i = 0; i < nJobs; i++) {
			n2 = n1 + nstep - 1; if (i == nJobs - 1 || n2 >= mdcell->mm.n) n2 = mdcell->mm.n - 1;
			n = n2 - n1 + 1; if (n <= 0) break;
			job[i].mmIndx.SetArray(n);
			for (imol = 0; imol < n; imol++) job[i].mmIndx.m[imol] = n1 + imol;
			n1 = n2 + 1;
		}
	}
	if (mdcell->sm.n > 0) {
		n1 = 0; nstep = mdcell->sm.n / nJobs;
		if (nstep == 0 || (nstep * nJobs + 1) < mdcell->sm.n) nstep += 1;
		for (i = 0; i < nJobs; i++) {
			n2 = n1 + nstep - 1; if (i == nJobs - 1 || n2 >= mdcell->sm.n) n2 = mdcell->sm.n - 1;
			n = n2 - n1 + 1; if (n <= 0) break;
			job[i].smIndx.SetArray(n);
			for (imol = 0; imol < n; imol++) job[i].smIndx.m[imol] = n1 + imol;
			n1 = n2 + 1;
		}
	}
	if (mdcell->pm.n > 0) {
		n1 = 0; nstep = mdcell->pm.n / nJobs;
		if (nstep == 0 || (nstep * nJobs + 1) < mdcell->pm.n) nstep += 1;
		for (i = 0; i < nJobs; i++) {
			n2 = n1 + nstep - 1; if (i == nJobs - 1 || n2 >= mdcell->pm.n) n2 = mdcell->pm.n - 1;
			n = n2 - n1 + 1; if (n <= 0) break;
			job[i].pmIndx.SetArray(n);
			for (imol = 0; imol < n; imol++) job[i].pmIndx.m[imol] = n1 + imol;
			n1 = n2 + 1;
		}
	}
};

/***********************************************************************************************
	MULTI-THREAD Operator on MDCELL_HDYN in MDCELL, OVER ALL MOLECULES
***********************************************************************************************/

template <class MDCELL> class MDCELL_HDYN_THREAD_VARS : public _THREAD_VARS {
public:
	MDCELL *mdcell;
	ARRAY<int> hdynIndx;

	void* func;

	MDCELL_HDYN_THREAD_VARS() {
		this->mdcell = NULL; func = NULL;
	};
	void set_pars(int iThread, MDCELL* mdcell) {
		this->mdcell = mdcell;
		_THREAD_VARS::set_pars(iThread);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~MDCELL_HDYN_THREAD_VARS() {
		this->mdcell = NULL; func = NULL;
	};
};


#if _SYS_ == _WINDOWS_SYS_
template <class MDCELL> int MDCELL_HDYN_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MDCELL> void* MDCELL_HDYN_Thread(void *void_par) { // thread-function
#endif
	MDCELL_HDYN_THREAD_VARS<MDCELL> *par = (MDCELL_HDYN_THREAD_VARS<MDCELL>*)void_par;
	typedef void (*MDCELL_HDYN_OPERATOR)(MDCELL*, ARRAY<int> &);
	((MDCELL_HDYN_OPERATOR)(par->func))(par->mdcell, par->hdynIndx);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class MDCELL> void AssignJob_MDCELL_HDYN(MDCELL *mdcell, MDCELL_HDYN_THREAD_VARS<MDCELL> *thread_var, int nThreads, bool bExcludeLast) {
	//MDCELL_HDYN_THREAD_VARS<MDCELL> *thread_var = new MDCELL_HDYN_THREAD_VARS<MDCELL>[nThreads];
	int iThread = 0, ihdyn;
	int n1, n2, nstep, n;
	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_var[iThread].hdynIndx.release(); 
		thread_var[iThread].set_pars(iThread, mdcell);
	}
	int nHdyn = (bExcludeLast ? mdcell->hdyn.n - 1 : mdcell->hdyn.n);
	if (nHdyn <= 0) return;

	n1 = 0; nstep = nHdyn / nThreads;
	if (nstep == 0 || (nstep * nThreads + 1) < nHdyn) nstep += 1;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = mdcell->hdyn.n - 1;
		n2 = n1 + nstep - 1; 
		if (bExcludeLast) {
			if (n2 >= mdcell->hdyn.n - 1) n2 = mdcell->hdyn.n - 2;
		}
		else {
			if (n2 >= mdcell->hdyn.n) n2 = mdcell->hdyn.n - 1;
		}
		n = n2 - n1 + 1; if (n <= 0) break;
		thread_var[iThread].hdynIndx.SetArray(n);
		for (ihdyn = 0; ihdyn < n; ihdyn++) thread_var[iThread].hdynIndx.m[ihdyn] = n1 + ihdyn;
		n1 = n2 + 1;
	}
};

template <class MDCELL> void MDCELL_HDyn_MultiThreadOperate(void *func, MDCELL_HDYN_THREAD_VARS<MDCELL> *thread_var, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*MDCELL_HDYN_OPERATOR)(MDCELL*, ARRAY<int> &);
		((MDCELL_HDYN_OPERATOR)(func))(thread_var->mdcell, thread_var->hdynIndx);
		return;
	}

	int iThread;
	for (iThread = 0; iThread < nThreads; iThread++) thread_var[iThread].set_func(func);
#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_var[iThread].hMMThread);
		thread_var[iThread].thread = AfxBeginThread((AFX_THREADPROC)MDCELL_HDYN_Thread<MDCELL>, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_var[iThread].hMMThread, INFINITE);
		//CloseHandle(thread_var[iThread].hMMThread); thread_var[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_var[iThread].thread), NULL, &MDCELL_HDYN_Thread<MDCELL>, (void *)(thread_var + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_var[iThread].thread != 0) pthread_join(thread_var[iThread].thread, NULL); // wait the thread over
	}
#endif
};

template <class MDCELL> class MDCELL_ThreadJob : public JOB_THREAD_VARS< MDCELL_Job<MDCELL> > {
public:
	MDCELL_ThreadJob() {};
	~MDCELL_ThreadJob() {};
};

template <class MDCELL, class CMM> class MDCELL_CMM_ThreadJob : public JOB_THREAD_VARS< MDCELL_Job<MDCELL> > {
public:
	CMM *cmm;
	
	MDCELL_CMM_ThreadJob() {cmm = NULL;};
	~MDCELL_CMM_ThreadJob() {};
};

void MDCELL_molecule_PeriodicCell_check(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
void MDCELL_SpeedVerlet(MDCELL_ThreadJob<MMOL_MD_CELL> *job, bool *status); // with confinement on total spatial momentum
void MDCELL_SpeedVerlet0(MDCELL_ThreadJob<MMOL_MD_CELL> *job, bool *status); // without confinement on total spatial momentum 
void MDCELL_Verlet(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
void MDCELL_CalMolKineticEnergy(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
void MDCELL_CalMolNetForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job); //for macromol only
void MDCELL_ClusterRealForceBackup(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
void MDCELL_CalInertiaMoment(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
void MDCELL_InitClusterForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
void MDCELL_ClusterForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
void MDCELL_ScaleClusterForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
void MDCELL_CheckRotAngle(MDCELL_ThreadJob<MMOL_MD_CELL> *job);

void MDCELL_HDyn(MMOL_MD_CELL* mdcell, ARRAY<int> &hdynIndx);
void MDCELL_CheckHDynCovergent(MMOL_MD_CELL* mdcell, ARRAY<int> &hdynIndx);
void MDCELL_AcceptCurrentDynForNextStep(MMOL_MD_CELL* mdcell, ARRAY<int> &hdynIndx);

bool LeapFrogVerlet_SelfConsistent_SpeedAccel(MMOL_MD_CELL *mdcell, MDCELL_ThreadJob<MMOL_MD_CELL> *job_mdcell, MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> *hdyn_thread_var, bool bSpeedVerlet, int nThreads);
// without confinement on spatial momentum
bool LeapFrogVerlet0_SelfConsistent_SpeedAccel(MMOL_MD_CELL *mdcell, MDCELL_ThreadJob<MMOL_MD_CELL> *job_mdcell, MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> *hdyn_thread_var, bool bSpeedVerlet, int nThreads);


void check_atom_rg(MMOL_MD_CELL &mmcell, int nThreads);
void update_atom_rg(MMOL_MD_CELL &mmcell);

//void MM_Torsion(MMOLECULE *mm, int nMM, double& Etorsion);


class INTERACT { // using for explicit interaction calculation
public:
	InteractVar var;
	EInteractVar evar;
	InteractRes LJ_res, eres;

	INTERACT() {};
	~INTERACT() {};
};

void MMCELL_CMM_cluster_allocate(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator, bool *status);

void MMCELL_CMM_CheckRelship(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator);
void MMCELL_CMM_CheckImpSolvRelship(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator);



void MMCELL_cluster_impsolv_interact(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *mdcell_impsolv, double *v);
void MMCELL_MM_Torsion(MDCELL_ThreadJob<MMOL_MD_CELL> *job, TorsionVar *var, InteractRes *ires);

void gen_random_force(MMOL_MD_CELL &mdcell, _rand_::RAND_ARRAY<int> &randa_m, _rand_::RAND_ARRAY<int> &randa_c);
void InitMD(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, bool bStorageChain);

void MMCELL_ExplicitInteract(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *job, INTERACT *interact);

// virial correction for rigid molecules
void MDCELL_VirialCorrect_RigidCluster(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4> *virial);
void MDCELL_VirialCorrect_Molecule(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4>* virial);
void MDCELL_ClusterVirialCorrect_HingeConstraint(MDCELL_ThreadJob<MMOL_MD_CELL> *job, HingeConstraintVar &var, InteractRes &ires);

void MDCellMassCenter(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob);

void setup_mol_kin(MMOL_MD_CELL &mdcell);

void MDCELL_CalCellInertiaTensor(MDCELL_ThreadJob<MMOL_MD_CELL> *job, MATRIX<3> *mI);

class AsyncMDSave_VARS : public _THREAD_VARS {
public:
	MMOL_MD_CELL *mdcell;
	AsyncMDSave_VARS() {mdcell = NULL;};
};
void AsyncMDSave(MMOL_MD_CELL &md_cell, AsyncMDSave_VARS &av);


namespace _bias_force_ {
	// use bias force to a position with constant force, for a molecule
	class BIAS_POINT {
	public :
		ARRAY<int> mtype; // type of molecule, 0 -- macromolecule, 1 -- SM, 2 -- PM
		ARRAY<int> mindx; // molecule index in the MD_CELL
		ARRAY<VECTOR3> org;
		
		// force is the force on unit mass, or actually the acceleration
		// to ensure each particle is moving in same speed in diffussion / dislocation
		ARRAY<float> force;

		void init_buffer(int N) {
			mtype.set_array(N);
			mindx.set_array(N);
			org.set_array(N);
			force.set_array(N);
		}
		void set_bias(int ibias, int mtype, int mindx, float force) {
			if (ibias < 0 || ibias >= this->mtype.n) return;
			this->mtype.m[ibias] = mtype;
			this->mindx.m[ibias] = mindx;
			this->force.m[ibias] = force;
		}
	};
}

void apply_point_bias_force(MMOL_MD_CELL&mdcell, _bias_force_::BIAS_POINT &mdcell_bias);

extern _bias_force_::BIAS_POINT mdcell_bias;




