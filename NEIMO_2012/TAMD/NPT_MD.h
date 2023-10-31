#define RF_NVAR  6

// dynamic of barostat for isotropic system
class ISOTROPIC_PDYN {
public:
	float W; // mass of the barostat
	double eps, v; // v is barostat momentum p / W, eps the Lagranian position in space relative to p in momentum, eps = 1/d * ln(V / V0)
	double alpha; // the acceleration of eps

	double Pint; // the internal pressure
	double fc; // the real force on p, dV(Pint - Pext) + 1/N * SUM(p^2/m) of MTK

	double V0, V; // volumn: initial and current
	double Ek; // kinetic energy of the barostat bath
};

void InitVelocity(ISOTROPIC_PDYN &pdyn, float Ek);

inline double calKineticEnergy(ISOTROPIC_PDYN &pdyn) {
	pdyn.Ek = pdyn.W * pdyn.v * pdyn.v * 0.5; 
	return pdyn.Ek;
};

class ISOTROPIC_P_MD_SAVE{
public:
	MD_SAVE *mdsave;
public:
	ISOTROPIC_P_MD_SAVE() {mdsave = NULL;};
	~ISOTROPIC_P_MD_SAVE() {mdsave = NULL;};
	bool save_kinetic_pars(ISOTROPIC_PDYN &pdyn, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync = false);
};

class ISOTROPIC_P_MD_BASE {
public:
	long LOOPS;
	float dt; // in fs

	double v0;
	double eps0;

	ISOTROPIC_P_MD_SAVE psave;

	ISOTROPIC_P_MD_BASE() {LOOPS = 1000; dt = 0;};
	~ISOTROPIC_P_MD_BASE() {};
	void set_dt(float dt) {this->dt = dt;};
};

class LFVERLET_ISOTROPIC_P_MD : public ISOTROPIC_P_MD_BASE {
public:
	double v_mh, v_ph;
	HOOVER_DYN *hdyn;

	double e_eps;

	LFVERLET_ISOTROPIC_P_MD() {hdyn = NULL;};
	~LFVERLET_ISOTROPIC_P_MD() {hdyn = NULL;};
	bool save_dyn_pars(ISOTROPIC_PDYN &pdyn, long nloop, bool check_status, char* additional_title, bool save_md_title);
	void Init(ISOTROPIC_PDYN &p) {
		ISOTROPIC_P_MD_BASE::v0 = p.v;
		v_mh = p.v; v_ph = p.v;
	};
	void set_HDyn(HOOVER_DYN *hdyn) {
		this->hdyn = hdyn;
	};
	void cal_e() { // use v_ph
		double d = v_ph * dt;
		e_eps = exp(d);
	};
};

template <class MDCELL, class CMM> class NPT {
public:
	MDCELL *mdcell;
	CMM *cmm;
	float V0; // the initial volumn of the cell

	float Pext; // in atmosphere

	int NF_pressure; // the translation part only which has contribution to the internal pressure

	float Ekt, RF; // translation energy Ek, SUM[r * f] of the cell
	SVECTOR<double, 3> virial;

	NPT() {
		mdcell = NULL; cmm = NULL; Pext = 1; V0 = 0; NF_pressure = 0; Ekt = 0; RF = 0;
		virial.v[0] = 0; virial.v[1] = 0; virial.v[2] = 0;
	};
	~NPT() {mdcell = NULL; cmm = NULL;};
	/*
	void release_vhdyn(bool bNHC) {
		if (vhdyn == NULL) return;
		if (bNHC && vhdyn->nhc != NULL) {delete vhdyn->nhc; vhdyn->nhc = NULL;};
		delete vhdyn; vhdyn = NULL;
	};
	*/
};

template <class MDCELL, class CMM> class ISOTROPIC_NPT : public NPT<MDCELL, CMM> {
public:
	ISOTROPIC_PDYN pdyn;
	LFVERLET_ISOTROPIC_P_MD pmdpar;

	HOOVER_DYN phdyn; // for PDyn

	ISOTROPIC_NPT() {pmdpar.hdyn = &phdyn;};
	void setup_pressure_NHC() {
		if (phdyn.nhc == NULL) phdyn.nhc = new NHC<MNHC>;
	};
	void release_pressure_NHC() {
		if (phdyn.nhc != NULL) {delete phdyn.nhc; phdyn.nhc = NULL;}
	}
	~ISOTROPIC_NPT() {release_pressure_NHC(); pmdpar.hdyn = NULL;};
};

class ISOTROPIC_NPT_CELL : public ISOTROPIC_NPT<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > {
public:
	ISOTROPIC_NPT_CELL() {};
	~ISOTROPIC_NPT_CELL() {};
};

int check_NF_pressure(MMOL_MD_CELL *mdcell);

void MDCELL_Ekt(MDCELL_ThreadJob<ISOTROPIC_NPT_CELL> *job, double *Ekt);

void MDCELL_PressureForce(MDCELL_ThreadJob<ISOTROPIC_NPT_CELL> *job, double *v, int nJob);

void cal_dyn(ISOTROPIC_PDYN &pdyn, float ksi);
void P_NHC(ISOTROPIC_PDYN &pdyn, HOOVER_DYN &phdyn);

void NPT_Speed_Kinetic(LFVERLET_ISOTROPIC_P_MD &mdpar, ISOTROPIC_PDYN &pdyn);
void NPT_Speed_Verlet(LFVERLET_ISOTROPIC_P_MD &mdpar, ISOTROPIC_PDYN &pdyn);
void NPT_Verlet(LFVERLET_ISOTROPIC_P_MD &mdpar, ISOTROPIC_PDYN &pdyn);

void NPT_Speed_Verlet(LFVERLET_MM_MD &mdpar, MMOLECULE &mm);
void NPT_Speed_Kinetic(LFVERLET_MM_MD &mdpar, MMOLECULE &mm);
void NPT_Verlet(LFVERLET_MM_MD &mdpar, MMOLECULE &mm);

void NPT_Speed_Verlet(LFVERLET_SM_MD &mdpar, MMOLECULE &sm);
void NPT_Speed_Kinetic(LFVERLET_SM_MD &mdpar, SMOLECULE &sm);
void NPT_Verlet(LFVERLET_SM_MD &mdpar, SMOLECULE &sm);

void NPT_Speed_Verlet(LFVERLET_PM_MD &mdpar, PMOLECULE &pm);
void NPT_Speed_Kinetic(LFVERLET_PM_MD &mdpar, PMOLECULE &pm);
void NPT_Verlet(LFVERLET_PM_MD &mdpar, PMOLECULE &pm);

void NPT_MD_PROC(ISOTROPIC_NPT_CELL &npt_cell);
