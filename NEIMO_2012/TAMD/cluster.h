/*
#define _HEAD     1
#define _TAIL     -1
#define _IN_CHAIN    0
*/
#define CH_BLENGTH   1.13f
#define CH1_BLENGTH  1.09f
#define CO_BLENGTH   1.43f

/*****************************
    BASIC_CLUSTER INDEX
******************************/
#define CINDX_SPC             0
#define CINDX_SPCE            1
#define CINDX_SPC_OH          2
#define CINDX_SPC_O           3

#define CINDX_BENZENE         10
#define CINDX_BENZENE1        11
#define CINDX_BENZENE12       12
#define CINDX_BENZENE13       13
#define CINDX_BENZENE14       14

#define CINDX_CH4             20
#define CINDX_CH3             21
#define CINDX_CH2             22
#define CINDX_CH              23

#define CINDX_C2H2            30
#define CINDX_C2H             31

#define CINDX_PEO_CH2         41
#define CINDX_PEO_O           42
#define CINDX_PEO_OH          43

#define CINDX_PPO_O           51
#define CINDX_PPO_OH          52


// to test the functions of Z-Matrix Operation
int SPC_O(BASIC_CLUSTER &c);
int SPC_OH(BASIC_CLUSTER &c);
int DoubleBoundAtom(BASIC_CLUSTER &c, char aindx, float bl, float ba, short cindx);

int SPCwater(BASIC_CLUSTER &c);
int SPCEwater(BASIC_CLUSTER &c);
// Constructing a benzene ring
int Benzene(BASIC_CLUSTER &c);

int Benzene1(CLUSTER &c);  // benzene-  C6H5- with a parent

int Benzene14(CLUSTER &c); // -benzene-  -C6H4- with hing at C1 & C4

int Benzene12(CLUSTER &c); // -benzene-  -C6H4- with hing at C1 & C2

int Benzene13(CLUSTER &c); // -benzene-  -C6H4- with hing at C1 & C3

// Constructing CH4
int CH4(CLUSTER &c);

int CH3(CLUSTER &c); // CH3-

int CH2(CLUSTER &c); //-CH2-

int CH(CLUSTER &c);  // -CH-- ; one parent, 2 children

int C2H2(CLUSTER &c);  // -CH=CH-

int C2H(CLUSTER &c);   // -CH=C-- ; one parent, with 2 children

int CH2_template(CLUSTER &c, char Cindx, char Hindx, float bl, short cindx); //-CH2-

int isoprene_monomer(MONOMER &m); //-CH2-C2H-[CH3]CH2-

int isoprene1_monomer(MONOMER &m); //CH3-C2H-[CH3]CH2-

int isoprene2_monomer(MONOMER &m); //-CH2-C2H-[CH3]CH3

int styrene_monomer(MONOMER &m);  //-CH-[C6H5]CH2-

int styrene1_monomer(MONOMER &m); //CH2-[C6H5]CH2-

int styrene2_monomer(MONOMER &m); //-CH-[C6H5]CH3

int ethylene_oxide_monomer(MONOMER &m); //-CH2-CH2-O-

int ethylene_oxide1_monomer(MONOMER &m); //-CH2-CH2-OH

int propylene_oxide_monomer(MONOMER &m); //-CH2-CH(CH3)-O-

int propylene_oxide1_monomer(MONOMER &m); //-CH2-CH(CH3)-OH

int propylene_oxide2_monomer(MONOMER &m); //CH3-CH(CH3)-O-

// copy the clusters in monomer to cluster array;
// and setup the same relationship between the clusters in monomer
bool copy_monomer(CLUSTER* cluster, MONOMER* monomer);

// twist the torsion angle of cluster along its parent hinge
void tweak(MMOLECULE &mm, bool bsave = false);

// unit is composed with cluster c in order, start with head and end with tail
// using C-C as hing bound to connect the clusters
void SingleChainPolymer(MMOLECULE &mm, MONOMER *head, MONOMER *m, int n_monomer, int *nperiod, float *ta, MONOMER *tail, float ta_tail, int &nbase, int nCG = 1);

bool get_period(char **mol, char *period);
bool get_period_number(char **str, int &n);
bool get_cluster_unit(char **str, char *unit);

bool defined_simple_molecule(char *mol);
bool defined_point_molecule(char *mol);
bool defined_cluster(char *unit);
bool defined_monomer(char *unit);

bool construct_cluster(BASIC_CLUSTER *pc, char *cname); // from structure database file
bool construct_defined_simple_molecule(char *mol, BASIC_CLUSTER* p);
bool construct_defined_point_molecule(char *mol, BASIC_CLUSTER* p);
bool construct_defined_cluster(char *unit, CLUSTER* p);
int construct_defined_monomer(char *monomer, MONOMER &m);
bool construct_monomer(char *monomer_name, MONOMER *m);

int nMonomerFromPolymerName(char *polymer);
bool ParseSingleChainPolymer(char *polymer, MONOMER **m, int &n_monomer, int **nperiod);
bool ConstructSingleChainPolymerMonomerChain(MMOLECULE &mm, char *polymer, float *ta, float ta_tail, int &nbase, int nCG = 1);

bool ConstructSingleChainPolymerFromParamFile(MMOLECULE &mm, char *fname, bool save = true);

bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol);
bool ConstructPointMolecule(PMOLECULE* pm, char *mol);

bool get_solvation_energy(CLUSTER_IMPSOLV_PARS *par);
bool readImplicitSolvationPars(CLUSTER_IMPSOLV_PARS_DATABASE &sdb, float &rSolvent);
void calImplicitSolvationArea(CLUSTER_IMPSOLV_PARS_DATABASE &sdb, float rSolvent);

bool readDeuterizedMonomer(char *fname, char *monomer);
void monomer_deuterization(MONOMER *m);

class MONOMER_LINK_RELSHIP {
public:
	int mIndx1, iHinge1; // the hinge indx of this monomer
	int mIndx2, iHinge2; // the hinge indx of another monomer
};

class MONOMER_GROUP_HINGE {
public:
	int mIndx, iHinge; // monomer, cluster, iHinge in cluster
};

class MONOMER_GROUP {
public:
	char name[256];
	ARRAY<MONOMER> monomer;
	ARRAY<MONOMER_LINK_RELSHIP> mLink; // the relationship indexed inside the group
	ARRAY<MONOMER_GROUP_HINGE> mHinge;
	MONOMER_GROUP() {strcpy(name, "\0");};
	void reset() {
		monomer.release(); mLink.release(); mHinge.release();
	};
	~MONOMER_GROUP() {reset();};
};

bool ReadMonomersFromFile(char *fname, char *mname, MONOMER_GROUP &mgroup);

bool ReadMacroMol(char *fname, MMOLECULE &mm, bool save = true);

bool get_atompar(ATOM_COMM_PARS *atompar);
