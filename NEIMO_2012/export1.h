using namespace _coarse_grain_;

class GENERAL_MOLECULE {
public:
	short mol_type;
	MMOLECULE* mm;
	CG_MMOLECULE *cgmm;
	SMOLECULE *sm;
	PMOLECULE *pm;
	CGSM *cgsm;
	GENERAL_MOLECULE() {mol_type = 0; mm = NULL; cgmm = NULL; sm = NULL; cgsm = NULL; pm = NULL;};
	void set_pmol(PMOLECULE *pm) {
		mol_type = 4; this->sm = NULL; this->mm = NULL; this->cgmm = NULL; this->cgsm = NULL; this->pm = pm;
	};
	void set_smol(SMOLECULE *sm) {
		mol_type = 0; this->sm = sm; this->mm = NULL; this->cgmm = NULL; this->cgsm = NULL; this->pm = NULL;
	};
	void set_mmol(MMOLECULE *mm) {
		mol_type = 1; this->mm = mm; this->sm = NULL; this->cgmm = NULL; this->cgsm = NULL; this->pm = NULL;
	};
	void set_cgmm(CG_MMOLECULE *cgmm) {
		mol_type = 2; this->cgmm = cgmm; this->sm = NULL; this->mm = NULL; this->cgsm = NULL; this->pm = NULL;
	};
	void set_cgsm(CGSM *cgsm) {
		mol_type = 3; this->cgsm = cgsm; this->mm = NULL; this->cgmm = NULL; this->sm = NULL; this->pm = NULL;
	};
};

class PICK_SCATT_STRUCT {
public:
	ARRAY<int> monomer_uid; // specific type of monomer ? 

	int nMol, nMonomer, nCluster, natom;
	void reset() {nMol = -1; nMonomer = -1; nCluster = -1; natom = -1;};
	// cluster iteration
	void set(int nMol, MMOLECULE &mm, int nCluster, int natom) {
		if (nCluster >= mm.nCluster) {reset(); return;}
		if (natom >= mm.cluster[nCluster].nAtoms) {nCluster++; natom = 0;}
		if (nCluster >= mm.nCluster) {reset(); return;}
		this->nMol = nMol;
		this->nCluster = nCluster; this->natom = natom;
	};
	// monomer iteration
	void set(int nMol, MMOLECULE &mm, int nMonomer, int nCluster, int natom) {
		if (nMonomer >= mm.nMonomer) {reset(); return;}
		if (natom >= mm.monomer[nMonomer].u[nCluster]->nAtoms) {nCluster++; natom = 0;}
		if (nCluster >= mm.monomer[nMonomer].nUnit) {nMonomer++; nCluster = 0; natom = 0;}
		if (nMonomer >= mm.nMonomer) {reset(); return;}
		this->nMol = nMol; this->nMonomer = nMonomer;
		this->nCluster = nCluster; this->natom = natom;
	};
	// for coarse-grained molecules, we use natom for iteration / index
	// in coarse-grained macromolecules, each monomer/cluster has only one atom inside
	void set(int nMol, CG_MMOLECULE &mm, int natom) {
		if (natom >= mm.nCluster) {reset(); return;}
		this->nMol = nMol;
		this->natom = natom;
	};
	bool exist_uid(int uid) {
		if (monomer_uid.m == NULL) return true; // if non-uid is specified, all uid is defined
		else return monomer_uid.exist(uid);
	};
	PICK_SCATT_STRUCT() {monomer_uid.SetArray(0); reset();};
};

bool get_molecule_upto_job(int nMol, GENERAL_MOLECULE& gm);

class READ_MD_PROC_FILE {
public:
	MD_SAVE msave;
	int nfile, nsaved; // currently opened file #, & the MD proc
	ifstream *in;

	bool bNPT;
	ifstream *PVin; // in NPT simulatioin, we need to update volumn of the system
	float xd[3]; // the original dimension of the periodical cell

	READ_MD_PROC_FILE() {in = NULL; nfile = -1; nsaved = -1; bNPT = false; PVin = NULL;};
	void reset_PV_opened_file() {
		if (PVin != NULL) {
			if (PVin->is_open()) PVin->close(); 
			delete PVin; PVin = NULL;
		}
	};
	void reset_opened_file() {
		if (in != NULL) {
			if (in->is_open()) in->close(); 
			delete in; in = NULL;
		}
		if (bNPT) reset_PV_opened_file();
		nfile = -1; nsaved = -1;
	};
	bool init(char *mdpar_fname);
	~READ_MD_PROC_FILE() {reset_opened_file();};
	bool search_proc(int isaved, char* buffer); // buffer saves the title -- first line
};

class CGATOM_SCATT_PARS : public COMM_PARS<5> {
public:
	float scf;
	CGATOM_SCATT_PARS() {scf = 0;};
	bool get_scatt(char *fname, int SOURCE);
};
/*
class CGATOM_SCF_DATABASE : public PARS_DATABASE<CGATOM_SCATT_PARS> {
public:
	CGATOM_SCF_DATABASE(bool (*get_atompar_func)(CGATOM_SCATT_PARS*)):PARS_DATABASE<CGATOM_SCATT_PARS>(get_atompar_func) {};
	~CGATOM_SCF_DATABASE() {release();};
	bool get_scatt_factor(char *fname, int source) {
		CHAIN<CGATOM_SCATT_PARS> *pch = this->ch;
		bool status = true;
		while (pch != NULL) {
			status = status & pch->p->get_scatt(fname, source);
			pch = pch->next;
		}
		return status;
	};
};
*/
// to get specific type of monomer
int get_monomer_total_scatt_atom_number(ARRAY<int>& monomer_uid, bool all_atom);
int get_monomer_scatt_struct(bool all_atom, float *fscatt, float *x, float *y, float *z, int MAX, int nMol, int nstart);

