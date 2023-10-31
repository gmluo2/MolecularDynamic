namespace _coarse_grain_ {
class CGATOM_COMM_PARS : public COMM_PARS<5> { // all these parameter are site indepedent
public:
	char pdb_alias[5]; // atom name saved for PDB file, corase-grained atom can not be recognized by PDB
	float m; //mass
	float epsLJ, rLJ; // L-J interaction
	float alpha; // polarizility
	float rcut; // longest distance for local interaction
	CGATOM_COMM_PARS() {
		strcpy(pdb_alias, "\0");
		indx = 0; uid = 0;
		strcpy(atom, "\0");
		m = 0; alpha = 0;
		epsLJ = 0; rLJ = 1; rcut = 1;
	};
};

bool get_atompar(CGATOM_COMM_PARS *atompar);
bool get_cgatom_pdb_alias(CGATOM_COMM_PARS *atompar, char *fname, char *title);
bool get_cgatom_uid(char *cgatom, int &indx);
//bool get_cgatom_from_alias(char *fname, char *alias, char *cgatom);

class CGATOM_COMM_PARS_DATABASE : public PARS_DATABASE<CGATOM_COMM_PARS> {
public:
	CGATOM_COMM_PARS_DATABASE(bool (*get_atompar_func)(CGATOM_COMM_PARS*)):PARS_DATABASE<CGATOM_COMM_PARS>(get_atompar_func) {};
	~CGATOM_COMM_PARS_DATABASE() {release();};
	bool get_pdb_alias(char *fname, char *title) {
		CHAIN<CGATOM_COMM_PARS> *pch = this->ch;
		bool status = true;
		while (pch != NULL) {
			status = status & get_cgatom_pdb_alias(pch->p, fname, title);
			pch = pch->next;
		}
		return status;
	};
};


// beads of the polymer chain are connected by
// finitely extensible nonlinear (FENE) springs, -1/2*k*R0^2ln[1-(r/R0)^2], R0 is the bound-length
class FENE_BOUND { // U = -0.5 * vb * blen * blen * ln[1 - (r/blen)^2]
public:
	float vb; // kJ/mol
	float blen; // bound length in Angstrom

	FENE_BOUND() {vb = 0; blen = 1;};
	~FENE_BOUND() {};
	void construct_prof(EF_DPROF &efprof);
};

class BOUND_PROF : public BOUND_BASE<char, 2> {
public:
	EF_DPROF efprof;
	BOUND_PROF() {};
	~BOUND_PROF(){efprof.release();};
};

class CGBOUND_DB : public _BOUND_DB<BOUND_PROF, char, 2> {
public:
	CGBOUND_DB() {};
	~CGBOUND_DB() {};
};

class CGBOUND_DIPOLE_DISTRBT : public BOUND_BASE<char, 2> {
public:
	_PROFILE_1D<float, float> mu, theta;
	CGBOUND_DIPOLE_DISTRBT() {};
	void normalize_distrbt() {
		int i = 0;
		float vmax = 0;
		for (i = 0; i < mu.N; i++) 	vmax = (vmax > mu.y[i] ? vmax : mu.y[i]);
		if (vmax > 0) mu.normalize_y(vmax);
		vmax = 0; for (i = 0; i < theta.N; i++) vmax = (vmax > theta.y[i] ? vmax : theta.y[i]);
		if (vmax > 0) theta.normalize_y(vmax);
		
	};
	~CGBOUND_DIPOLE_DISTRBT() {mu.release(); theta.release();};
};

class CGBOUND_DIPOLE_DISTRBT_DB : public _BOUND_DB<CGBOUND_DIPOLE_DISTRBT, char, 2> {
public:
	CGBOUND_DIPOLE_DISTRBT_DB() {};
	~CGBOUND_DIPOLE_DISTRBT_DB() {};
};


class CGATOM : public _ATOM<CGATOM> {
public:
	float c; // charge -- normalized with 1/sqrt(eps), is not a common parameter from AMBER force-field parameter
	float c0; // charge 
	CGATOM_COMM_PARS *par;

	IMPSOLV_PARS *psolv;

	// ****** parameters for electrostatic field ******
	VECTOR3 mu; // dipole moment
#if _ATOM_INDUCED_DIPOLE
	VECTOR3 ind_mu; // induced dipole moment
	DMATRIX<3> Eg; // gradient electric field 
#endif
	VECTOR3 E; // electric field from charge, multipoles of other atoms
	// ****** end of electrostatic field ******

	ARRAY<VECTOR3> f;

	CGATOM() {
		par = NULL;
		c = 0; c0 = 0;
		psolv = NULL;
	};
	~CGATOM() {
		par = NULL; _ATOM<CGATOM>::reset_bounds();
		psolv = NULL;
	};
	/*
	bool is_bound(CGATOM* p) {
		return _ATOM::is_bound((_ATOM*)p);
	};
	int ibound(CGATOM* p) {
		return _ATOM::ibound((_ATOM*)p);
	};
	*/
	void rotate(R &r) {};
};

class KINETIC {
public:
	VECTOR3 v, alpha; // speed and acceleration speed
	double Ek; // kinetic energy
	bool hoover; // to speed up /down
};

class DYNAMIC {
public:
	VECTOR3 fc; // real force acting on the cluster
	bool overlap;
	VECTOR3 f_Brown;
	VECTOR3 f; // the total force, fc + f_Brown

	DYNAMIC() {overlap = false;};
};

class CG_CLUSTER;

class LOCAL_RELSHIP : public IMAG_RELSHIP<CG_CLUSTER> {
public:
	LOCAL_RELSHIP() {};
	~LOCAL_RELSHIP() {};
};

// IMPORTANT -- assuming that each cluster has only one cgatom only
class CG_CLUSTER : public _BASIC_CLUSTER<CGATOM> {
public:
	KINETIC *bkm;
	DYNAMIC *dyn;

	VECTOR3 *r; // the cordinate of cluster, pointing to atom[0].r
	float *M; // mass of the cluster

	float d; // the size of the cluster -- used for Brownian force calculation
	float ita_Brownian, cd_trans; // drag force of continum solvent model
	
	short cID; // ID of the type of cluster 
	short mIndx; // the molecule index in the whole system, where this cluster is
	//short cgIndx; // index of coarse-grained cluster where this cluster is in
	//short monIndx; // monomer index inside the macro-molecule

	bool eNeutral;
	//ELECTRO_MULTIPOLE *emp;  // ignore the multipole of the CLUSTER
	//CMM_RELSHIP ERship;
	LOCAL_RELSHIP LRship;

	//ARRAY<CG_CLUSTER*> c12;
	//ARRAY<CG_CLUSTER*> c13;
	//ARRAY<CG_CLUSTER*> c14;

	ARRAY< ARRAY<CG_CLUSTER*> > neighbor;
	ARRAY< ARRAY<BOUND_PROF*> > cgb;

	ARRAY<CG_CLUSTER*> non_neighbor_sm; // non-neighbors in the same molecule

	//TORSION *pTors;

	bool highlight; // highlight this cluster in VMM, for view only

	CLUSTER_IMPSOLV_PARS *psolv;
	VECTOR3 rc_solv; // geometrical center of cluster, in real and local cordinates
	IMPSOLV_RELSHIP<CG_CLUSTER> impsolv_rship; // neighbor cluster in chain for solvation consideration
	double E_ImplicitSolvation;

	VECTOR3 rg_cmm;
	VECTOR3 *vp0; // vp = vp0 + dr_cmm;
	
	void update_vp() {
		if (this->vp0 == NULL) return;
		if (bCMMShift) {V3plusV3((*vp0), dr_cmm, (*vp))}
		else {V32V3((*vp0), (*vp))}
	};

	CG_CLUSTER() {
		eNeutral = true; cID = 0; mIndx = 0;
		highlight = false;
		bkm = NULL; dyn = NULL; r = NULL; M = 0;
		//pTors = NULL; emp = NULL;
		ita_Brownian = 0; cd_trans = 0;

		psolv = NULL;
		E_ImplicitSolvation = 0;
	};
	
	void setup_kin_dyn() {
		if (bkm == NULL) bkm = new KINETIC;
		if (dyn == NULL) dyn = new DYNAMIC;
		int i;
		if (::bMTForce) {
			for (i = 0; i < nAtoms; i++) atom[i].f.set_array(_NFS);
		}
	};
	void reset_kin_dyn() {
		if (bkm != NULL) {delete bkm; bkm = NULL;}
		if (dyn != NULL) {delete dyn; dyn = NULL;}
	};
	/*
	void setup_emp() {
		if (!this->eNeutral && emp == NULL) emp = new ELECTRO_MULTIPOLE;
	};
	*/
	void reset() {
		_BASIC_CLUSTER<CGATOM>::reset(); 
		reset_kin_dyn();
		//if (emp != NULL) {delete emp; emp = NULL;}

		psolv = NULL;
		impsolv_rship.release_relationship();
	};
	~CG_CLUSTER() {
		reset();
		//ERship.release_relationship();
		LRship.release_relationship();
	};

	bool linkPars(CGATOM_COMM_PARS_DATABASE &cgatompar_db, char *errmsg) {
		int i = 0, j = 0;
		CGATOM_COMM_PARS *atompar = NULL;
		for (i = 0; i < nAtoms; i++) {
			if (atom[i].par != NULL) continue;
			atom[i].par = cgatompar_db.search_par_indx((char)atom[i].aindx);
			if (atom[i].par == NULL) {
				sprintf(errmsg, "atom parameter of #%d is not defined in program", int(atom[i].aindx));
				return false;
			}
		}
#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 1
		this->eNeutral = true;
#else
		check_neutral();
		//this->setup_emp();
#endif
		return true;
	};

	//void calElectroMultipole(bool bEwaldSum);
	void check_neutral() {
		int na; eNeutral = true;
		for (na = 0; na < nAtoms; na++) {
			if (atom[na].c != 0) {eNeutral = false; return;}
		}
	};
	//float cal_geometry_radius();
	void setup_cordinate() {
		if (atom == NULL) {this->r = NULL; this->vp0 = NULL; this->M = NULL;}
		else {this->r = &(atom[0].r); this->vp0 = this->r; this->M = &(atom[0].par->m);}
		this->vp = &rg_cmm;
	};
	void setup_neighbors(int n_neighbors);
	int cBound(CG_CLUSTER* pc) {
		if (pc == this) return 0;
		int ib, nc = 0;
		for (ib = 0; ib < neighbor.n; ib++) {
			for (nc = 0; nc < neighbor.m[ib].n; nc++) {
				if (neighbor.m[ib].m[nc] == pc) return ib + 1;
			}
		}
		return -1; // unknown
	};
	void calBrownianPars();
};

class CGBOUND : public _BOUND_<CG_CLUSTER> {
public:
	CGBOUND_DIPOLE_DISTRBT *dt;
	VECTOR3 mu; // dipole moment
	// polar cordinates of dipole relative to the bound direction
	// these values will be generated randomly
	float vmu, theta; 
	VECTOR3 randv; // a random-vector, mu = rotate bound director, with axis (bound x  randv)

	CGBOUND(){dt = NULL; vmu = 0; theta = 0;};
	~CGBOUND() {dt = NULL;};
};

void cp_cgatom(CGATOM *d, CGATOM *s);

bool cp(CG_CLUSTER *d, CG_CLUSTER *s);

bool ConnectClusters(CG_CLUSTER *c1, int iHinge1, CG_CLUSTER *c2, int iHinge2, char bound_type = SINGLE_BOUND);

class CGMM_KIN_DYN {
public:
	float Ekt, Ekr; // translation kinetic energy, and rotation kinetic energy of frame at whole molecular mass center
	VECTOR<6> P; // first three -- translation spatial moment, the 2nd three -- angular spatial moment relative to mass center
	MATRIX<3> M, invM; // inertia moment of the whole molecule, and its inverse
	VECTOR3 v, w; // velocity of mass center, angular velocity relative to mass center

	VECTOR<6> f; // total force and torque, relative to mass center
};

class CG_MMOLECULE {
public:
	short mIndx; // the index of the molecule in the whole system
	short mType; // the type of the molecule
	char mol[20];
	int nCluster;
	CG_CLUSTER *cluster;

	CGBOUND *bound;
	int nBound;

	float M; // total mass
	VECTOR3 cm; // mass center

	float Ek; // kinetic and potential energy

	CGMM_KIN_DYN *mKD;

	CG_MMOLECULE() {
		mIndx = 0; mType = 0; strcpy(mol, "\n");
		cluster = NULL; nCluster = 0; Ek = 0; mKD = NULL; bound = NULL; nBound = 0;
	};
	void reset() {
		if (cluster != NULL) {delete[] cluster; cluster = NULL;}; nCluster = 0;
		if (mKD != NULL) {delete mKD; mKD = NULL;}
		if (bound != NULL) {delete[] bound; bound = NULL;} nBound = 0;
	};
	~CG_MMOLECULE() {reset();};
	bool xyz_save(ofstream &out);
	void setMolIndx(int indx) {this->mIndx = indx; for (int i = 0; i < nCluster; i++) cluster[i].mIndx = indx;};

	void shiftMol(VECTOR3 &dr);

	void calMassCenter();
	void SetMassCenter(double *v, bool cal_mass_center);
	//void calTotalForce();

	void setup_kin_dyn() {
		if (mKD == NULL) mKD = new CGMM_KIN_DYN;
		for (int nc = 0; nc < nCluster; nc++) cluster[nc].setup_kin_dyn();
	};
	void setup_cluster_neighbors(int n_neighbors) {
		for (int nc = 0; nc < nCluster; nc++) cluster[nc].setup_neighbors(n_neighbors);
	};
	void calBrownianPars() {
		for (int nc = 0; nc < nCluster; nc++) cluster[nc].calBrownianPars();
	};
	void setup_bounds();
	void bound_dipole();
	void check_non_neighbor_sm(float r2cut);
};

void calDipoleEffect(CG_MMOLECULE *mm, int nc1, int nc2, double &V);

// return energy in kT
void calCGClustersKineticEnergy(CG_MMOLECULE *mm, int n1, int n2, double &Ek);
// return energy in kT
void ParallelCalKineticEnergy(CG_MMOLECULE &mm);

void calInertiaMoment(CG_MMOLECULE *mm);
void calFrameSpatialMomentum(CG_MMOLECULE *mm); // in parallel

void calFrameEnergy(CG_MMOLECULE *mm);

// calculate acceleration
void calClusterDynamic(CG_MMOLECULE *mm, int n1, int n2);
// calculate the acceleration of each cluster in mm
void ParallelCalDynamic(CG_MMOLECULE &mm);

//bool construct_coarse_grained_mm(MMOLECULE &mm, CG_MMOLECULE &cgmm);

bool atom_indx(CG_MMOLECULE &cgmm, CGATOM* pcgatom, int &nc, int &na, int &indx);

void attach_FENE_bound_db(CG_MMOLECULE &cgmm, CGBOUND_DB &cgbound_db);  // FENE
void attach_neighbor_bound_db_from_file(CG_MMOLECULE &cgmm, CGBOUND_DB &cgbound_db); // FROM FILE
void attach_bound_db(CG_MMOLECULE &cgmm, CGBOUND_DB &cgbound_db, int cgbound_type);
void cgbound_interact(CG_MMOLECULE *cgmm, int nc1, int nc2, double& V);

void attach_non_neighbor_db(CGATOM_COMM_PARS_DATABASE &cgatompar_db, CGBOUND_DB &cg_non_neighbor_db);
void attach_intramol_non_neighbor_db(CGATOM_COMM_PARS_DATABASE &cgatompar_db, CGBOUND_DB &cg_non_neighbor_db);

class CGSM : public CG_CLUSTER {
public:
	short mIndx;
	short mType;
	char mol[20];

	float Ek; //kinetic energy
	CGSM() {Ek = 0; mIndx = 0; mType = -1; strcpy(mol, "\n");};
	~CGSM() {};
	bool xyz_save(ofstream &out) {
		out<<atom[0].par->atom<<" "<<r->v[0]<<" "<<r->v[1]<<" "<<r->v[2]<<endl;
		return true;
	};
	void shiftMol(VECTOR3 &dr) {this->shift(dr);};
};

double calCGSMKineticEnergy(CGSM &sm);
void calCGSMDynamic(CGSM &sm);

bool cp(CG_MMOLECULE *d, CG_MMOLECULE *s);
bool cp(CGSM *d, CGSM *s);

bool ConstructCGSM(CGSM* cgsm, char *mol);
bool construct_defined_cgsm(char *unit, CGSM* m);

} // end of namespace _coarse_grain_ 


namespace _coarse_grain_ {
	extern CGATOM_COMM_PARS_DATABASE cgatompar_db;
	extern char cgatom_par_fname[256], cgatom_index_fname[256];
}
