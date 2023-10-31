namespace _coarse_grain_ {

// following classes are for macro-molecule
class CGMM_MD_SAVE {
public:
	MD_SAVE *mdsave;
	bool save_kinetic_pars(CG_MMOLECULE &cgmm, long nloop, bool check_status, char* additional_title, bool save_md_title);
	CGMM_MD_SAVE() {mdsave = NULL;};
	~CGMM_MD_SAVE() {mdsave = NULL;};
};

bool save_xyz_struct(CG_MMOLECULE &mm, char *fname);

class CGMM_MD_BASE : public MD_BASE {
public:
	int nCluster;

	VECTOR3 *v;

	CGMM_MD_SAVE mmsave;

	CGMM_MD_BASE() {
		nCluster = 0; v = NULL;
	};
	~CGMM_MD_BASE() {
		if (v != NULL) {delete[] v; v = NULL;}
	};
	void set_dt(float dt) {this->dt = dt;};
	void resetClusters(int n) {
		int i = 0;
		if (v != NULL) {delete[] v; v = NULL;}
		if (n <= 0) nCluster = 0;
		else {
			nCluster = n;
			v = new VECTOR3[n];
		}
	};
};

// leapfrog Verlet algorithm & self-consistent
class LFVERLET_CGMM_MD : public CGMM_MD_BASE {
public:
	VECTOR3 *v_ph, *v_mh; // speed of +1/2 and -1/2 of current position
	HOOVER_DYN *hdyn;
	
	LFVERLET_CGMM_MD() {v_ph = NULL; v_mh = NULL; hdyn = NULL;};
	void resetClusters(int n) {
		CGMM_MD_BASE::resetClusters(n);
		if (v_ph != NULL) {delete[] v_ph; v_ph = NULL;}
		if (v_mh != NULL) {delete[] v_mh; v_mh = NULL;}
		if (n <= 0) return;
		else {
			v_ph = new VECTOR3[n];
			v_mh = new VECTOR3[n];
		}
	};
	~LFVERLET_CGMM_MD() {
		if (v_ph != NULL) {delete[] v_ph; v_ph = NULL;}
		if (v_mh != NULL) {delete[] v_mh; v_mh = NULL;}
		hdyn = NULL;
	};
	bool save_mol_dyn_pars(CG_MMOLECULE &mm, long nloop, bool check_status, char* additional_title, bool save_md_title);
	void Init(CG_MMOLECULE &mm) {
		int i = 0;
		for (i = 0; i < this->nCluster; i++) {
			if (mm.cluster[i].bkm != NULL) CGMM_MD_BASE::v[i] = mm.cluster[i].bkm->v;
		}
	};
	void set_HDyn(HOOVER_DYN *hdyn, HOOVER_DYN *iHDyn) {
		this->hdyn = hdyn;
	};
};

void Speed_Kinetic(LFVERLET_CGMM_MD &mdpar, CG_MMOLECULE &mm);
void Speed_Verlet(LFVERLET_CGMM_MD &mdpar, CG_MMOLECULE &mm, bool reset_mol_speed);
void Verlet(LFVERLET_CGMM_MD &mdpar, CG_MMOLECULE &mm);

class CG_MD_VARS {
public:
	float dt, tao, max_ksi;
	int nCMMRefresh;
	CG_MD_VARS() {dt = 1; tao = 6; max_ksi = 0.1f; nCMMRefresh = 10;};
};



/**********************************************************************************
						CGSM MD leapfrog Verlet algorithm
**********************************************************************************/
class CGSM_MD_SAVE {
public:
	MD_SAVE *mdsave;
	bool save_kinetic_pars(CGSM &cgsm, long nloop, bool check_status, char* additional_title, bool save_md_title);
	CGSM_MD_SAVE() {mdsave = NULL;};
	~CGSM_MD_SAVE() {mdsave = NULL;};
};

bool save_xyz_struct(CGSM &sm, char *fname);

class CGSM_MD_BASE : public MD_BASE {
public:
	VECTOR3 v;

	CGSM_MD_SAVE smsave;

	CGSM_MD_BASE() {};
	~CGSM_MD_BASE() {};
	void set_dt(float dt) {this->dt = dt;};
};

// leapfrog Verlet algorithm & self-consistent
class LFVERLET_CGSM_MD : public CGSM_MD_BASE {
public:
	VECTOR3 v_ph, v_mh; // speed of +1/2 and -1/2 of current position
	HOOVER_DYN *hdyn;
	
	LFVERLET_CGSM_MD() {hdyn = NULL;};
	~LFVERLET_CGSM_MD() {hdyn = NULL;};
	bool save_mol_dyn_pars(CGSM &sm, long nloop, bool check_status, char* additional_title, bool save_md_title);
	void Init(CGSM &sm) {
		if (sm.bkm == NULL) return;
		V32V3(sm.bkm->v, v)
	};
	void set_HDyn(HOOVER_DYN *hdyn) {
		this->hdyn = hdyn;
	};
};

void Speed_Verlet(LFVERLET_CGSM_MD &mdpar, CGSM &sm, bool reset_mol_speed);
void Speed_Kinetic(LFVERLET_CGSM_MD &mdpar, CGSM &sm);
void Verlet(LFVERLET_CGSM_MD &mdpar, CGSM &sm);

} // end of namespace _coarse_grain_ 
