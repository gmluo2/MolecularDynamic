#define _VMD_CELL_ 1

class FLEX_MMOL_CELL : public _MD_CELL_DIM {
public:
	// molecules are in use
	FLEX_ARRAY<MMOLECULE> in_mm;
	FLEX_ARRAY<SMOLECULE> in_sm;
	FLEX_ARRAY<PMOLECULE> in_pm;

	// molecules are not in use, that can be used in future
	FLEX_ARRAY<MMOLECULE> un_mm;
	FLEX_ARRAY<SMOLECULE> un_sm;
	FLEX_ARRAY<PMOLECULE> un_pm;

	void release_inuse() {
		in_mm.release(true); in_sm.release(true); in_pm.release(true);
	};
	void release_unuse() {
		un_mm.release(true); un_sm.release(true); un_pm.release(true);
	};
	MMOLECULE* use_new_mm() {
		MMOLECULE* m = NULL;
		if (un_mm.a.n > 0) {
			m = un_mm.a.m[0];
			un_mm.remove(0);
		}
		else m = new MMOLECULE;
		in_mm.attach(m);
		return m;
	};
	void attach_inuse_mm(MMOLECULE* m, int n) {
		in_mm.attach(m, n);
	};
	void unuse_mm(MMOLECULE* m) { // put m to un_mm
		int indx = in_mm.indx(m);
		if (indx > 0) in_mm.remove(indx);
		un_mm.attach(m);
	};
	SMOLECULE* use_new_sm() {
		SMOLECULE* m = NULL;
		if (un_sm.a.n > 0) {
			m = un_sm.a.m[0];
			un_sm.remove(0);
		}
		else m = new SMOLECULE;
		in_sm.attach(m);
		return m;
	};
	void attach_inuse_sm(SMOLECULE* m, int n) {
		in_sm.attach(m, n);
	};
	void unuse_sm(SMOLECULE* m) { // put m to un_mm
		int indx = in_sm.indx(m);
		if (indx > 0) in_sm.remove(indx);
		un_sm.attach(m);
	};
	PMOLECULE* use_new_pm() {
		PMOLECULE* m = NULL;
		if (un_pm.a.n > 0) {
			m = un_pm.a.m[0];
			un_pm.remove(0);
		}
		else m = new PMOLECULE;
		in_pm.attach(m);
		return m;
	};
	void attach_inuse_pm(PMOLECULE* m, int n) {
		in_pm.attach(m, n);
	};
	void unuse_pm(PMOLECULE* m) { // put m to un_mm
		int indx = in_pm.indx(m);
		if (indx > 0) in_pm.remove(indx);
		un_pm.attach(m);
	};
	~FLEX_MMOL_CELL() {release_inuse(); release_unuse();};
};

class MMOL_VMD_CELL : public _MD_CELL_DIM {
public:
	ARRAY<MMOLECULE*> mm;
	ARRAY<LFVERLET_MM_MD*> lfmd;

	ARRAY<SMOLECULE*> sm;
	ARRAY<LFVERLET_SM_MD*> lfsmd;

	ARRAY<PMOLECULE*> pm;
	ARRAY<LFVERLET_PM_MD*> lfpmd;

	// index of the molecule using for fractal molecule, MM, SM, PM, which are defined in mm, sm, and pm
	ARRAY<int> fmm;
	ARRAY<int> fsm;
	ARRAY<int> fpm;

	int NF;
	float Ethermal;
	HOOVER_DYN hdyn; 

	float M; // total mass
	VECTOR3 rc; // mass center
	VECTOR3 Vmc, Wmc; // mass-center velocity of the whole cell, the rotation of the whole cell
	MATRIX<3> I, invI; // inertia moment of the whole cell, and its inverse matrix

	MDCELL_DYN_BUFFER<6> mvb;

	ARRAY< MOL_POL<MMOLECULE> > ind_mu_mm; // for the polarized dipole
	ARRAY< MOL_POL<SMOLECULE> > ind_mu_sm; // for the polarized dipole
	ARRAY< MOL_POL<PMOLECULE> > ind_mu_pm; // for the polarized dipole

	bool eNeutral;
	bool bHingeConstraint;

	MD_SAVE *msave;
	long LOOPS;

	MMOL_VMD_CELL() {
		NF = 0; Ethermal = 0; eNeutral = false; M = 0; msave = NULL; bHingeConstraint = false;
	};
	void check_eneutral() {
		eNeutral = true;
		int i, ic;
		for (i = 0; i < mm.n; i++) { for (ic = 0; ic < mm.m[i]->nCluster; ic++) {
			mm.m[i]->cluster[ic].check_eneutral();
		}}
		for (i = 0; i < sm.n; i++) { 
			sm.m[i]->c->check_eneutral();
		}
		for (i = 0; i < pm.n; i++) { 
			pm.m[i]->c->check_eneutral();
		}

		for (i = 0; i < mm.n; i++) { for (ic = 0; ic < mm.m[i]->nCluster; ic++) {
			if (!mm.m[i]->cluster[ic].eNeutral) {eNeutral = false; return;}
		}}
		for (i = 0; i < sm.n; i++) { 
			if (!sm.m[i]->c->eNeutral) {eNeutral = false; return;}
		}
		for (i = 0; i < pm.n; i++) { 
			if (!pm.m[i]->c->eNeutral) {eNeutral = false; return;}
		}
	};
	void SetMol(int n_mm, int n_sm, int n_pm = 0) {
		mm.SetArray(n_mm); lfmd.SetArray(n_mm); 
		sm.SetArray(n_sm); lfsmd.SetArray(n_sm);
		pm.SetArray(n_pm); lfpmd.SetArray(n_pm);
	};
	void SetSave(MD_SAVE *msave) {
		int i = 0; 
		for (i = 0; i < lfmd.n; i++) lfmd.m[i]->mmsave.mdsave = msave;
		for (i = 0; i < lfsmd.n; i++) lfsmd.m[i]->smsave.mdsave = msave;
		for (i = 0; i < lfpmd.n; i++) lfpmd.m[i]->pmsave.mdsave = msave;
		this->msave = msave;
	};
	void setMD(long nloops) {
		int i;
		for (i = 0; i < lfmd.n; i++) lfmd.m[i]->LOOPS = nloops;
		for (i = 0; i < lfsmd.n; i++) lfsmd.m[i]->LOOPS = nloops;
		for (i = 0; i < lfpmd.n; i++) lfpmd.m[i]->LOOPS = nloops;
		this->LOOPS = nloops;
	};
	void set_md_timestep(float dt) {
		int i;
		for (i = 0; i < lfmd.n; i++) lfmd.m[i]->set_dt(dt);
		for (i = 0; i < lfsmd.n; i++) lfsmd.m[i]->set_dt(dt);
		for (i = 0; i < lfpmd.n; i++) lfpmd.m[i]->set_dt(dt);
	};
	bool MDSave(long nloop, float Ep, bool bEachMolDyn = true);
	void SetMDMemory() {
		int i = 0; for (i = 0; i < lfmd.n; i++) {lfmd.m[i]->resetClusters(mm.m[i]->nCluster);}
	};
	int max_ncluster() {
		int i = 0, nmax = 0;
		for (i = 0; i < mm.n; i++) if (nmax < mm.m[i]->nCluster) nmax = mm.m[i]->nCluster;
		return nmax;
	};
	void release() {
		mm.release(); lfmd.release();
		sm.release(); lfsmd.release();
		pm.release(); lfpmd.release();
	};
	
	void release_hdyn(bool release_nhc = false) {
		if (hdyn.nhc != NULL) {
			if (release_nhc) delete hdyn.nhc; 
			hdyn.nhc = NULL;
		}
	};
	void create_hdyn() {
		release_hdyn(); 
		hdyn.nhc = new NHC<MNHC>;
	};
	void check_freedom();  // freedoms for each NHC

	~MMOL_VMD_CELL() {
		release(); release_hdyn(true);
	};
	void setMolIndx() {
		int n = 0, i = 0;
		for (i = 0; i < this->mm.n; i++) {this->mm.m[i]->setMolIndx(n); n++;}
		for (i = 0; i < this->sm.n; i++) {this->sm.m[i]->c->mIndx = n; n++;}
		for (i = 0; i < this->pm.n; i++) {this->pm.m[i]->c->mIndx = n; n++;}
	};
	void Init_MD() {
		int i = 0;
		for (i = 0; i < mm.n; i++) lfmd.m[i]->Init(*(mm.m[i]));
		for (i = 0; i < sm.n; i++) lfsmd.m[i]->Init(*(sm.m[i]));
		for (i = 0; i < pm.n; i++) lfpmd.m[i]->Init(*(pm.m[i]));
	};
	void CalibrateTorsionAngle() {
		int i = 0;
		for (i = 0; i < mm.n; i++) lfmd.m[i]->CalibrateTorsionAngle(*(mm.m[i]));
	};
	void init_mdcell_velocity_memory();
	void init_polarization_buff();
	void set_atom_multipole(int mMP);
	void TotalMass() {
		int ic, im = 0;
		M = 0;
		for (im = 0; im < mm.n; im++) {
			for (ic = 0; ic < mm.m[im]->nCluster; ic++) M += mm.m[im]->cluster[ic].M;
		}
		for (im = 0; im < sm.n; im++) M += sm.m[im]->c->M;
		for (im = 0; im < pm.n; im++) M += pm.m[im]->c->M;
	};
	MD_BASE* md_infor() {
		if (lfpmd.m != NULL) return (MD_BASE*)(lfpmd.m[0]);
		if (lfsmd.m != NULL) return (MD_BASE*)(lfsmd.m[0]);
		if (lfmd.m != NULL) return (MD_BASE*)(lfmd.m[0]);
		return NULL;
	};
};

int nAtoms(MMOL_VMD_CELL &mdcell);

void molecule_check_periodic_cell(MMOL_VMD_CELL &mmcell);
void VMDCELL_molecule_PeriodicCell_check(MDCELL_ThreadJob<MMOL_VMD_CELL> *job);
