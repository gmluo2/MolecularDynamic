class ITER_SAVE {
protected:
	int icount, nsave;
	char fname[256];
	char title[20], ext[20];
	int ifile;
	int isave, max_save;
public:
	ostringstream sbuf;
	//char buffer[51200];
	//ostream sbuf;
	ofstream* out;
	ITER_SAVE() {
		icount = 0; nsave = 0; isave = 0; max_save = 0; out = NULL;
		strcpy(fname, "\0"); strcpy(title, "\0"); strcpy(ext, "\0");
		ifile = 0;
		//strcpy(buffer, "\0"); sbuf.rdbuf()->pubsetbuf(buffer, 51200);
	};
	void reset(bool bResetCount = true) {
		if (bResetCount) {icount = 0; ifile = 0; isave = 0;}
		if (out != NULL) {if (out->is_open()) out->close(); delete out; out = NULL;}
	};
	int save_indx() {return isave + ifile * max_save;};
	int file_indx() {return ifile;};
	~ITER_SAVE() {reset();};
	void set_iter(int nsave, int max_save) {this->nsave = nsave; this->max_save = max_save;};
	void set_file(char *ftitle, char *fext) {strcpy(title, ftitle); strcpy(ext, fext);};
	bool one_more() {
		char errmsg[256] = "\0";
		icount++; if (icount >= nsave) icount = 0;
		if (icount != 0) return true;
		if (out == NULL) {
			sprintf(fname, "%s-%d.%s", title, ifile, ext);
			out = new ofstream; 
			if (isave == 0) out->open(fname); // ios_base::out | ios_base::trunc);
			else out->open(fname, ios_base::out | ios_base::app);
			if (!out->is_open()) {
				delete out; out->close();
				sprintf(errmsg, "can not open file %s", fname);
				show_msg(errmsg); return false;
			}
		}
		return true;
	};
	bool bsave() {return (icount == 0 ? true : false);};
	void save_over() {
		if (icount != 0) return;
		isave++;
		if (isave >= max_save) {reset(false); isave = 0; ifile++;}
		//out->close(); delete out; out = NULL;
	};
	void init_buf() {
		sbuf.clear();
	};
	void flush_buf(bool bEnforce = false) { // used only when ostringstream is used
		if (out == NULL || !out->is_open()) {
			if (bEnforce) {
				sbuf.clear(); sbuf.str("");
			}
			return;
		}
		if (bEnforce || sbuf.tellp() > 25600) {
			sbuf.clear(); /*(*out)<<sbuf.buf;*/ (*out)<<sbuf.str();
			(*out).flush(); sbuf.clear(); 
			sbuf.str(""); sbuf.clear();
		} 
	};
};

class MD_SAVE {
public:
	ITER_SAVE mESave, mKinSave;
	ITER_SAVE mDynSave; // for dynamic parameters
	ITER_SAVE muSave; // for dipole, polarized

	char title[20];
	int nsave; //save every nsave
	int max_save; // maximum number of loop saved in each file

	bool bBuffered;

	MD_SAVE() {strcpy(title, "\0"); nsave = 0; max_save = 0; bBuffered = false;};
	~MD_SAVE() {};
	void set_save(int nsave, int max_save) {
		this->nsave = nsave; this->max_save = max_save;
		mESave.set_iter(nsave, max_save);
		mKinSave.set_iter(nsave, max_save);
		mDynSave.set_iter(nsave, max_save);
		muSave.set_iter(nsave, max_save);
	};
	void set_title(char *title) {
		strcpy(this->title, title);
		mESave.set_file(title, "E");
		mKinSave.set_file(title, "kin");
		mDynSave.set_file(title, "dyn");
		muSave.set_file(title, "amu");
	};
	bool save_energy(float Ek, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync = false);

	void init_buf() {
		mESave.init_buf();
		mKinSave.init_buf();
		mDynSave.init_buf();
		muSave.init_buf();
	};
	void flush_buf(bool bEnforce = true) {
		mESave.flush_buf(bEnforce);
		mKinSave.flush_buf(bEnforce);
		mDynSave.flush_buf(bEnforce);
		muSave.flush_buf(bEnforce);
		if (bEnforce) bBuffered = false;
	};
	void close() {
		mESave.reset();
		mKinSave.reset();
		mDynSave.reset();
		muSave.reset();
	};
};

class HOOVER_DYN {
public:
	// Nose-Hoover parameter -- tao_s
	/*
	double tao_s, inv_tao_ss; // inv_tao_ss = 1 / tao_s^2
	double ksi_Hoover, vksi_Hoover, max_ksi_Hoover; // Hoover parameter
	*/
	short hID;
	NHC<MNHC> *nhc;
	bool hdyn_status;

	HOOVER_DYN() {
		/*
		tao_s = 1;
		inv_tao_ss = 1 / (tao_s * tao_s);
		ksi_Hoover = 0; vksi_Hoover = 0;
		max_ksi_Hoover = 0.01;
		*/
		hID = -1;
		nhc = NULL;
	};
	void set_vars(float tao, float dt, float max_ksi) {
		/*
		this->tao_s = tao; inv_tao_ss = 1 / (tao_s * tao_s);
		*/
		if (nhc != NULL) nhc->set_vars(tao, dt, max_ksi);
	};
	void release(bool release_nhc) {
		if (release_nhc) {
			if (nhc != NULL) {delete nhc; nhc = NULL;}
		}
		else nhc = NULL;
	};
	~HOOVER_DYN() {release(false);};
};

class MM_HOOVER_DYN {
public:
	// for macro-molecule
	double ksi, vksi; 

	HOOVER_DYN *iHDyn; // for internal movement

	MM_HOOVER_DYN() {
		// for macro-molecule
		ksi = 0; vksi = 0;
		iHDyn = NULL;
	};
	~MM_HOOVER_DYN() {
		iHDyn = NULL;
	};
	void set_HDyn(HOOVER_DYN *ftHDyn, HOOVER_DYN *frHDyn, HOOVER_DYN *iHDyn) {
		this->iHDyn = iHDyn;
	};
	void fetch_ksi() {
		if (iHDyn != NULL) {
			ksi = iHDyn->nhc->ksi.v[0]; vksi = iHDyn->nhc->vksi.v[0];
		}
	};
};

// following classes are for macro-molecule
class MM_MD_SAVE {
public:
	MD_SAVE *mdsave;
	bool save_kinetic_pars(MMOLECULE &mm, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync = false);
	bool save_atomic_dipole(MMOLECULE &mm, long nloop, bool check_status, char *additional_title, bool save_md_title, bool bAsync = false);
	bool save_energy(MMOLECULE &mm, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync = false);
	MM_MD_SAVE() {mdsave = NULL;};
	~MM_MD_SAVE() {mdsave = NULL;};
};

bool save_xyz_struct(MMOLECULE &mm, char *fname);

class MD_BASE {
public:
	long LOOPS;
	double dt; // in fs
	MD_BASE() {LOOPS = 0; dt = 1;};
	void set_dt(float dt) {this->dt = dt;};
	~MD_BASE() {};
};

class MM_MD_BASE : public MD_BASE {
public:
	int nCluster;

	double *t, *tp;
	//VECTOR<6> *P; // momentum pf each cluster
	
	VECTOR<6> v0; // base cluster
	VECTOR<6> Iw; // total momentum respect to base cluster

	VECTOR<6> V0; // speed of mass center
	VECTOR<6> P0; // total momentum respect to mass center

	MM_MD_SAVE mmsave;

	MM_MD_BASE() {
		nCluster = 0;
		t = NULL; tp = NULL; //P = NULL;
	};
	~MM_MD_BASE() {
		if (t != NULL) {delete[] t; t = NULL;}
		if (tp != NULL) {delete[] tp; tp = NULL;}
		//if (P != NULL) {delete[] P; P = NULL;}
	};
	void resetClusters(int n) {
		int i = 0;
		if (t != NULL) {delete[] t; t = NULL;}
		if (tp != NULL) {delete[] tp; tp = NULL;}
		//if (P != NULL) {delete[] P; P = NULL;}
		if (n <= 0) nCluster = 0;
		else {
			nCluster = n;
			t = new double[n];
			tp = new double[n];
			//P = new VECTOR<6>[n];
		}
		for (i = 0; i < n; i++) {t[i] = 0; tp[i] = 0;}
	};
};

// leapfrog Verlet algorithm & self-consistent
class LFVERLET_MM_MD : public MM_HOOVER_DYN, public MM_MD_BASE {
public:
	double *tp_ph, *tp_mh; // speed of +1/2 and -1/2 of current position
	//VECTOR<6> *p_ph, *p_mh; // momentum of each cluster, at +1/2 and -1/2 of current position
	VECTOR<6> v0_ph, v0_mh; // base cluster
	VECTOR<6> Iw_mh, Iw_ph; // spatial mometum respect to the base

	VECTOR<6> P0_ph, P0_mh; // total momentun respect to the mass center
	VECTOR<6> V0_ph, V0_mh; // speed of mass center
	
	LFVERLET_MM_MD() {
		tp_ph = NULL; tp_mh = NULL; 
		//p_ph = NULL; p_mh = NULL;
	};
	void resetClusters(int n) {
		MM_MD_BASE::resetClusters(n);
		int i = 0;
		if (tp_ph != NULL) {delete[] tp_ph; tp_ph = NULL;}
		if (tp_mh != NULL) {delete[] tp_mh; tp_mh = NULL;}
		//if (p_ph != NULL) {delete[] p_ph; p_ph = NULL;}
		//if (p_mh != NULL) {delete[] p_mh; p_mh = NULL;}
		if (n <= 0) return;
		else {
			tp_ph = new double[n];
			tp_mh = new double[n];
			//p_ph = new VECTOR<6>[n];
			//p_mh = new VECTOR<6>[n];
		}
		for (i = 0; i < n; i++) {
			tp_ph[i] = 0; tp_mh[i] = 0;
		}
	};
	~LFVERLET_MM_MD() {
		if (tp_ph != NULL) {delete[] tp_ph; tp_ph = NULL;}
		if (tp_mh != NULL) {delete[] tp_mh; tp_mh = NULL;}
		//if (p_ph != NULL) {delete[] p_ph; p_ph = NULL;}
		//if (p_mh != NULL) {delete[] p_mh; p_mh = NULL;}
	};
	bool save_mol_dyn_pars(MMOLECULE &mm, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync = false);
	void Init(MMOLECULE &mm) {
		int i = 0;
		for (i = 0; i < this->nCluster; i++) {
			MM_MD_BASE::t[i] = mm.cluster[i].km.ta;
			MM_MD_BASE::tp[i] = mm.cluster[i].km.tp;
			//MM_MD_BASE::P[i] = mm.cluster[i].bkm->P0;
		}
		for (i = 0; i < 6; i++) {
			MM_MD_BASE::v0.v[i] = mm.V0.v[i]; // base velocity
			MM_MD_BASE::Iw.v[i] = mm.Iw.v[i]; // total momentum respect to the base cluster
			MM_MD_BASE::V0.v[i] = mm.mDynKin.V.v[i]; // velocity at mass center 
			MM_MD_BASE::P0.v[i] = mm.mDynKin.P.v[i]; // total momentum at mass center

			v0_mh.v[i] = v0.v[i]; v0_ph.v[i] = v0.v[i];
			Iw_mh.v[i] = Iw.v[i]; Iw_ph.v[i] = Iw.v[i];
			V0_mh.v[i] = V0.v[i]; V0_ph.v[i] = V0.v[i];
			P0_mh.v[i] = P0.v[i]; P0_ph.v[i] = P0.v[i];
		}
	};
	void CalibrateTorsionAngle(MMOLECULE &mm) {
		int i = 0;
		CLUSTER *pc = NULL;
		BASIC_CLUSTER *parent = NULL;
		for (i = 0; i < this->nCluster; i++) {
			pc = mm.cluster + i;
			parent = (BASIC_CLUSTER*)pc->parent;
			if (parent == NULL) continue;
			pc->km.ta = torsion_angle(parent, (BASIC_CLUSTER*)pc);
			MM_MD_BASE::t[i] = mm.cluster[i].km.ta;
		}
	};
};

double KineticEnergy_Foward(MMOLECULE &mm, LFVERLET_MM_MD &mdpar);

// torsion angle could be not in range [-PI, PI], reset the value
void check_angle(MMOLECULE &mm, LFVERLET_MM_MD &mdpar);

void Speed_Kinetic(LFVERLET_MM_MD &mdpar, MMOLECULE &mm, bool check_freebase);
void Speed_Verlet(LFVERLET_MM_MD &mdpar, MMOLECULE &mm);
void Verlet(LFVERLET_MM_MD &mdpar, MMOLECULE &mm);
void Verlet_Relax(LFVERLET_MM_MD &mdpar, MMOLECULE &mm, double dta_max);

// Verlet function is splitted to two function
// this function gives the movement parameters, and change mdpar only, but does not move the macro-molecule
void VerletStep(LFVERLET_MM_MD &mdpar, MMOLECULE &mm, VECTOR<6> &dv_base, double *dta);
// this function move the macro-molecule with given movement parameters, and mdpar
void VerletMove(LFVERLET_MM_MD &mdpar, VECTOR<6> &dv_base, double *dta, MMOLECULE &mm);

// following classes are for simple solid molecule
class SM_MD_SAVE{
public:
	MD_SAVE *mdsave;
public:
	SM_MD_SAVE() {mdsave = NULL;};
	~SM_MD_SAVE() {mdsave = NULL;};
	bool save_kinetic_pars(SMOLECULE &m, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync = false);
	bool save_atomic_dipole(SMOLECULE &mm, long nloop, bool check_status, char *additional_title, bool save_md_title, bool bAsync = false);
};

class SM_MD_BASE {
public:
	long LOOPS;
	double dt; // in fs

	VECTOR<6> v0; // speed
	VECTOR<3> p0; // angular momentum

	SM_MD_SAVE smsave;

	SM_MD_BASE() {LOOPS = 1000; dt = 0.1f;};
	~SM_MD_BASE() {};
	void set_dt(float dt) {this->dt = dt;};
};

// leapfrog Verlet algorithm & self-consistent
class LFVERLET_SM_MD : public SM_MD_BASE {
public:
	VECTOR<6> v0_ph, v0_mh; // position of base cluster -- translation, and rotation
	// we have keep tracking momentum of molecule, because leapfrog algorithm is working on momentum, while not speed
	// Since mass does not change with time, translation speed can work on leapfrog algorithm. However, this is NOT true 
	// for angular momentum, because inertial momentum varies with orientation
	// angular momentum of base cluster -- NOT the rotation speed, because inertial momentum of body is changing
	VECTOR<3> p0_ph, p0_mh; 

	HOOVER_DYN *hdyn;

	LFVERLET_SM_MD() {hdyn = NULL;};
	~LFVERLET_SM_MD() {hdyn = NULL;};
	bool save_mol_dyn_pars(SMOLECULE &m, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync = false);
	void Init(SMOLECULE &sm) {
		int i = 0;
		for (i = 0; i < 6; i++) {
			SM_MD_BASE::v0.v[i] = sm.c->bkm->V0.v[i];
			v0_mh.v[i] = v0.v[i]; v0_ph.v[i] = v0.v[i];
		}
		for (i = 0; i < 3; i++) {
			SM_MD_BASE::p0.v[i] = sm.c->bkm->P0.v[i];
			p0_mh.v[i] = p0.v[i]; p0_ph.v[i] = p0.v[i];
		}
	};
	void set_HDyn(HOOVER_DYN *hdyn) {
		this->hdyn = hdyn;
	};
};

void Speed_Kinetic(LFVERLET_SM_MD &mdpar, SMOLECULE &m);
void Speed_Verlet(LFVERLET_SM_MD &mdpar, SMOLECULE &m);
void Verlet(LFVERLET_SM_MD &mdpar, SMOLECULE &m);


// following classes are for point molecule
class PM_MD_SAVE{
public:
	MD_SAVE *mdsave;
public:
	PM_MD_SAVE() {mdsave = NULL;};
	~PM_MD_SAVE() {mdsave = NULL;};
	bool save_kinetic_pars(PMOLECULE &m, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync = false);
	bool save_atomic_dipole(PMOLECULE &m, long nloop, bool check_status, char *additional_title, bool save_md_title, bool bAsync = false);
};

class PM_MD_BASE {
public:
	long LOOPS;
	double dt; // in fs

	// PM does not need to use momentum, since no rotation exist
	VECTOR3 v0;

	PM_MD_SAVE pmsave;

	PM_MD_BASE() {LOOPS = 1000; dt = 0.1f;};
	~PM_MD_BASE() {};
	void set_dt(float dt) {this->dt = dt;};
};

// leapfrog Verlet algorithm & self-consistent
class LFVERLET_PM_MD : public PM_MD_BASE {
public:
	VECTOR3 v0_ph, v0_mh; // position of base cluster
	HOOVER_DYN *hdyn;

	LFVERLET_PM_MD() {hdyn = NULL;};
	~LFVERLET_PM_MD() {hdyn = NULL;};
	bool save_mol_dyn_pars(PMOLECULE &m, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync = false);
	void Init(PMOLECULE &pm) {
		int i = 0;
		for (i = 0; i < 3; i++) {
			PM_MD_BASE::v0.v[i] = pm.v.v[i];
			v0_mh.v[i] = v0.v[i]; v0_ph.v[i] = v0.v[i];
		}
	};
	void set_HDyn(HOOVER_DYN *hdyn) {
		this->hdyn = hdyn;
	};
};

void Speed_Kinetic(LFVERLET_PM_MD &mdpar, PMOLECULE &m);
void Speed_Verlet(LFVERLET_PM_MD &mdpar, PMOLECULE &m, bool reset_mol_speed);
void Verlet(LFVERLET_PM_MD &mdpar, PMOLECULE &m);


