namespace _coarse_grain_ {

class CGMM_MD_CELL : public _MD_CELL_DIM {
public:
	ARRAY<CG_MMOLECULE> mm;
	ARRAY<LFVERLET_CGMM_MD> lfmd;

	ARRAY<CGSM> sm;
	ARRAY<LFVERLET_CGSM_MD> lfsmd;

	ARRAY<CGSM> pm;
	ARRAY<LFVERLET_CGSM_MD> lfpmd;

	int NF; // total freedom
	float Ethermal;
	ARRAY<MDCELL_HDYN> hdyn; 
	MDCELL_DYN_BUFFER<3> mvb;

	// useless for CGMM, but required for template function
	bool bIndependentFrameHDyn;
	ARRAY<MDCELL_HDYN> mmHdyn; 

	CGMM_MD_CELL() {};
	void SetMol(int n_mm, int n_sm, int n_pm = 0) {
		mm.SetArray(n_mm); lfmd.SetArray(n_mm); 
		sm.SetArray(n_sm); lfsmd.SetArray(n_sm);
		pm.SetArray(n_pm); lfpmd.SetArray(n_pm);
	};
	void SetSave(MD_SAVE *msave) {
		int i = 0; 
		for (i = 0; i < lfmd.n; i++) lfmd.m[i].mmsave.mdsave = msave;
		for (i = 0; i < lfsmd.n; i++) lfsmd.m[i].smsave.mdsave = msave;
		for (i = 0; i < lfpmd.n; i++) lfpmd.m[i].smsave.mdsave = msave;
	};
	void SetMD(long nloops) {
		int i = 0; for (i = 0; i < lfmd.n; i++) lfmd.m[i].LOOPS = nloops;
		for (i = 0; i < lfsmd.n; i++) lfsmd.m[i].LOOPS = nloops;
		for (i = 0; i < lfpmd.n; i++) lfpmd.m[i].LOOPS = nloops;
	};
	void set_md_timestep(float dt) {
		int i;
		for (i = 0; i < lfmd.n; i++) lfmd.m[i].set_dt(dt);
		for (i = 0; i < lfsmd.n; i++) lfsmd.m[i].set_dt(dt);
		for (i = 0; i < lfpmd.n; i++) lfpmd.m[i].set_dt(dt);
	};
	bool MDSave(long nloop, float Ep);
	void SetMDMemory() {
		int i = 0; for (i = 0; i < lfmd.n; i++) {lfmd.m[i].resetClusters(mm.m[i].nCluster);}
	};
	int max_ncluster() {
		int i = 0, nmax = 0;
		for (i = 0; i < mm.n; i++) if (nmax < mm.m[i].nCluster) nmax = mm.m[i].nCluster;
		return nmax;
	};
	void release() {
		mm.SetArray(0); lfmd.SetArray(0); 
		sm.SetArray(0); lfsmd.SetArray(0);
		pm.SetArray(0); lfpmd.SetArray(0);
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
	void set_HDyn(int nHDyn) {
		int i;
		release_hdyn(); if (nHDyn <= 0) return;
		hdyn.SetArray(nHDyn);
		for (i = 0; i < hdyn.n; i++) {
			hdyn.m[i].nhc = new NHC<MNHC>;
			hdyn.m[i].hID = -2;
		}
	};
	HOOVER_DYN* get_HDyn(short mType) {
		HOOVER_DYN *pdyn = NULL;
		int i;
		for (i = 0; i < hdyn.n; i++) {if (hdyn.m[i].hID == mType) {pdyn = hdyn.m + i; break;}}
		if (pdyn == NULL) {
			for (i = 0; i < hdyn.n; i++) {if (hdyn.m[i].hID == -1) {pdyn = hdyn.m + i; break;}}
		}
		return pdyn;
	};
	/*
	void SetHooverDyn(int ihdyn, int hID, float dt, float tao, float max_ksi) {
		int i = 0; 
		HOOVER_DYN *pdyn = NULL;
		if (ihdyn < 0 || ihdyn >= hdyn.n) return;
		pdyn = hdyn.m + ihdyn;
		pdyn->hID = hID;
		hdyn.m[i].set_vars(tao, dt, max_ksi);
	};
	*/
	void init_nhc_mols();  // to group the molecules into each NHC
	void check_freedom();  // freedoms for each NHC


	~CGMM_MD_CELL() {release(); release_hdyn();};
	void setMolIndx() {
		int n = 0, i = 0;
		for (i = 0; i < this->pm.n; i++) {this->pm.m[i].mIndx = n; n++;}
		for (i = 0; i < this->sm.n; i++) {this->sm.m[i].mIndx = n; n++;}
		for (i = 0; i < this->mm.n; i++) {this->mm.m[i].setMolIndx(n); n++;}
	};
	void Init() {
		int i = 0;
		for (i = 0; i < mm.n; i++) lfmd.m[i].Init(mm.m[i]);
		for (i = 0; i < sm.n; i++) lfsmd.m[i].Init(sm.m[i]);
		for (i = 0; i < pm.n; i++) lfpmd.m[i].Init(pm.m[i]);
	};
};

class CGMD_THREAD_PAR : public _THREAD_VARS {
public:
	int md_loop, nM;
	double *Ek;

	char msg_show[5120], msg_log[5120];

	void reset() {
		RELEASE_ARRAY(Ek)
		strcpy(msg_show, "\0"); strcpy(msg_log, "\0");
	};
	void set(int nM, int iThread) {
		reset();
		_THREAD_VARS::set_pars(iThread);
		this->nM = nM; this->md_loop = 0;
		Ek = new double[nM];
		for (int i = 0; i < nM; i++) Ek[i] = 0;
	};
	void set_loop(int md_loop) {
		this->md_loop = md_loop;
	};
	CGMD_THREAD_PAR() {
		strcpy(msg_show, "\0"); strcpy(msg_log, "\0");
		Ek = NULL;
	};
	~CGMD_THREAD_PAR() {reset();};
};


class CGMM_MD_THREAD_PAR : public CGMD_THREAD_PAR {
public:
	CG_MMOLECULE *mm;
	LFVERLET_CGMM_MD *lfmdpar;
	VECTOR3 *vbf;

	CGMM_MD_THREAD_PAR() {mm = NULL; lfmdpar = NULL; vbf = NULL;};
	void set(CG_MMOLECULE *mm, LFVERLET_CGMM_MD* lfmdpar, int nMM, int iThread) {
		this->mm = mm; this->lfmdpar = lfmdpar; 
		CGMD_THREAD_PAR::set(nMM, iThread);
	};
	void reset_buffer_memory() {RELEASE_ARRAY(vbf)};
	void reset() {
		CGMD_THREAD_PAR::reset(); reset_buffer_memory();
	};
	void set_buffer_memory(int n) {
		reset_buffer_memory(); this->vbf = new VECTOR3[n];
	};
	~CGMM_MD_THREAD_PAR() {reset();};
};

class CGSM_MD_THREAD_PAR : public CGMD_THREAD_PAR {
public:
	CGSM* sm;
	LFVERLET_CGSM_MD *lfsmdpar;

	CGSM_MD_THREAD_PAR() {sm = NULL; lfsmdpar = NULL;};
	void set(CGSM *sm, LFVERLET_CGSM_MD* lfsmdpar, int nSM, int iThread) {
		this->sm = sm; this->lfsmdpar = lfsmdpar;
		CGMD_THREAD_PAR::set(nSM, iThread);
	};
	~CGSM_MD_THREAD_PAR() {CGMD_THREAD_PAR::reset();};
};

void set_free_cell_mass_center(CGMM_MD_CELL &mmcell, VECTOR3 &O);
void cluster_check_periodic_cell(CGMM_MD_CELL &mmcell);
void molecule_check_periodic_cell(CGMM_MD_CELL &mmcell);
void CGSM_MULTI_THREAD_MD(CGSM_MD_THREAD_PAR* sm_mdpar, int nThreads);
void CGMM_MULTI_THREAD_MD(CGMM_MD_THREAD_PAR* mm_mdpar, int nThreads);
void InitCGClusterForceInCMM(CMM_CELL3D<CG_CLUSTER> &cell, int n1Cell, int n2Cell);

void MD_PROC(CGMM_MD_CELL &mmcell, CMM_CELL3D<CG_CLUSTER> *cmm, int md_mode);

} // end of namespace _coarse_grain_ 

