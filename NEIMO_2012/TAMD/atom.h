
template <int N> class COMM_PARS { // all these parameter are site indepedent
public:
	short indx; // atom index -- type in NEIMO 
	short uid; // type of atom, in VMM
	char atom[N];
	COMM_PARS() {
		indx = 0; uid = 0;
		strcpy(atom, "\0");
	};
};

// for different type of atom
class ATOM_COMM_PARS : public COMM_PARS<5> { // all these parameter are site indepedent
public:
	float m; //mass
	float alpha; // polarizility
	bool bLJ;
	float epsLJ, rLJ; // L-J interaction
	float r_min; //minimum radius, when the distance closer than the minimum distance, atom is shocked

	bool bSR; // explicit short-range interaction, other than Lennard-Jones

	ATOM_COMM_PARS() {
		m = 0; alpha = 0; r_min = 0.5f;
		bLJ = true; epsLJ = 0; rLJ = 1;
		bSR = false;
	};
};

// TPAR has to be derived class of COMM_PARS, like ATOM_COMM_PARS
template <class TPAR> class PARS_DATABASE {
	bool (*get_atompar)(TPAR*);
public:
	CHAIN<TPAR> *ch;
	ARRAY<TPAR*> par;

	PARS_DATABASE() {ch = NULL; get_atompar = NULL;};
	PARS_DATABASE(bool (*get_atompar_func)(TPAR*)) {ch = NULL; get_atompar = get_atompar_func;};
	void release() {
		int i;
		for (i = 0; i < par.n; i++) par.m[i] = NULL;
		release_chain<TPAR>(&ch, true);
	};
	~PARS_DATABASE() {release();};

	TPAR* attach(char *atom) {
		CHAIN<TPAR> *pch = ch;
		TPAR *par = NULL;
		while (pch != NULL) {
			par = pch->p;
			if (strcmp(par->atom, atom) == 0) return par; // this par is in the chain
			else if (pch->next == NULL) break;
			else pch = pch->next;
		};
		par = new TPAR; strcpy(par->atom, atom);
		if (get_atompar != NULL && !(*get_atompar)(par)) {delete par; return NULL;}
		if (pch == NULL) {
			ch = new CHAIN<TPAR>; ch->p = par; par->indx = 0;
		}
		else {
			pch->next = new CHAIN<TPAR>; pch->next->p = par; par->indx = pch->p->indx + 1;
		}
		return par;
	};
	TPAR* attach(short uid) {
		CHAIN<TPAR> *pch = ch;
		TPAR *par = NULL;
		while (pch != NULL) {
			par = pch->p;
			if (par->uid == uid) return par; // this par is in the chain
			else if (pch->next == NULL) break;
			else pch = pch->next;
		};
		par = new TPAR; par->uid = uid;
		if (get_atompar != NULL && !(*get_atompar)(par)) {delete par; return NULL;}
		if (pch == NULL) {
			ch = new CHAIN<TPAR>; ch->p = par; par->indx = 0;
		}
		else {
			pch->next = new CHAIN<TPAR>; pch->next->p = par; par->indx = pch->p->indx + 1;
		}
		return par;
	};
	TPAR* search_par(char *atom) {
		CHAIN<TPAR> *pch = ch;
		TPAR *par = NULL;
		while (pch != NULL) {
			par = pch->p;
			if (strcmp(par->atom, atom) == 0) return par;
			else pch = pch->next;
		}
		return NULL;
	};
	TPAR* search_par_indx(char indx) {
		CHAIN<TPAR> *pch = ch;
		TPAR *par = NULL;
		while (pch != NULL) {
			par = pch->p;
			if (par->indx == indx) return par;
			else pch = pch->next;
		}
		return NULL;
	};
	TPAR* search_par_uid(short uid) {
		CHAIN<TPAR> *pch = ch;
		TPAR *par = NULL;
		while (pch != NULL) {
			par = pch->p;
			if (par->uid == uid) return par;
			else pch = pch->next;
		}
		return NULL;
	};
	void chain2array() {
		chain_array(&ch, par);
	};
};

class ATOM_COMM_PARS_DATABASE : public PARS_DATABASE<ATOM_COMM_PARS> {
public:
	ATOM_COMM_PARS_DATABASE(bool (*get_atompar_func)(ATOM_COMM_PARS*)):PARS_DATABASE<ATOM_COMM_PARS>(get_atompar_func) {};
	~ATOM_COMM_PARS_DATABASE() {release();};
};


/*********************************************
***** database for bound-type parameters *****
*********************************************/
template <class INDX_TYPE, int N> class BOUND_BASE {
public:
	INDX_TYPE aindx[N]; // atom index -- type in NEIMO
	unsigned int key; // key from 
};

inline unsigned int construct_char2_key_pairwise(char *aindx) {
	char a1_indx, a2_indx;
	if (aindx[0] > aindx[1]) {a1_indx = aindx[0]; a2_indx = aindx[1];}
	else {a1_indx = aindx[1]; a2_indx = aindx[0];}
	unsigned int key = a1_indx;
	key = key<<4;
	key += a2_indx;
	return key;
};

inline unsigned int construct_char2_key_non_pairwise(char *aindx) {
	unsigned int key = aindx[0];
	key = key<<4;
	key += aindx[1];
	return key;
};

inline unsigned int construct_char3_key_12_pairwise(char *aindx) {
	char a1_indx, a2_indx;
	if (aindx[0] > aindx[1]) {a1_indx = aindx[0]; a2_indx = aindx[1];}
	else {a1_indx = aindx[1]; a2_indx = aindx[0];}
	unsigned int key = a1_indx; key = key<<4;
	key += a2_indx; key = key<<4;
	key += aindx[2];
	return key;
};

inline unsigned int construct_char3_key_non_pairwise(char *aindx) {
	unsigned int key = aindx[0]; key = key<<4;
	key += aindx[1]; key = key<<4;
	key += aindx[2];
	return key;
};

template <class BPAR, class INDX_TYPE, int N> class _BOUND_DB { // BPAR is derived class from BOUND_BASE
public:
	CHAIN<BPAR> *bch; // database in chain
	ARRAY<BPAR*> ba; // database in array, easir to search
	_BOUND_DB() {bch = NULL;};
	void releaseChain(bool release_context) {if (bch != NULL) release_chain<BPAR>(&bch, release_context);};
	void releaseArray(bool release_context) {
		for (int i = 0; i < ba.n; i++) {
			if (release_context) {if (ba.m[i] != NULL) {delete ba.m[i]; ba.m[i] = NULL;}}
			else ba.m[i] = NULL;
		}
		ba.release();
	};
	~_BOUND_DB() {releaseChain(true); releaseArray(true);}
	void Chain2Array() {
		int n = number<BPAR>(bch);
		ba.release();
		if (n == 0) return;
		CHAIN<BPAR> *pch = bch;
		ba.SetArray(n);
		n = 0;
		while (pch != NULL) {
			ba.m[n] = pch->p; pch->p = NULL;
			pch = pch->next; n++;
		}
		releaseChain(false);
	};
	BPAR* search_chain(unsigned int key) {
		CHAIN<BPAR> *pch = bch;
		while (pch != NULL) {
			if (pch->p->key == key) return pch->p;
			pch = pch->next;
		}
		return NULL;
	};
	BPAR* search(unsigned int key) {
		for (int i = 0; i < ba.n; i++) {if (ba.m[i]->key == key) return ba.m[i];}
		return NULL;
	};
	void attach(BPAR* bpar) {
		if (bch == NULL) {bch = new CHAIN<BPAR>; bch->p = bpar;}
		else bch->attach2tail(bpar);
	};
};


// pairwise interaction (or other pairwise operator)
// the profile related to this pairwise operator
template <class OPERATOR> class PAIRWISE : public BOUND_BASE<char, 2> {
public:
	OPERATOR pv;
	PAIRWISE() {};
	~PAIRWISE(){};
};

template <class PAIRWISE> class PAIRWISE_DB : public _BOUND_DB<PAIRWISE, char, 2> {
public:
	PAIRWISE_DB() {};
	~PAIRWISE_DB() {};
};


// parameters for exponential potential:
// V_Exp = A * exp(-B * r)
#define _SR_EXP_     0  // type of short-range interaction
#define _N_EXP_      2
#define _EXP_A_      0
#define _EXP_B_      1

// parameters for Buckingham exponential-6 (Buck6) potential:
// V_Buck6 = A * exp(-B * r) - C / r^6
#define _SR_BUCK6_   1  // type of short-range interaction
#define _N_BUCK6_    3
#define _BUCK6_A_    0
#define _BUCK6_B_    1
#define _BUCK6_C_    2

// parameters for Tosi-Fumi (TF) potential:
// V = A * exp(-B * r) - C6 / r^6 - C8 / r^8 - C10 / r^10
#define _SR_TF_   2  // type of short-range interaction
#define _N_TF_    5
#define _TF_A_    0
#define _TF_B_    1
#define _TF_C6_   2
#define _TF_C8_   3
#define _TF_C10_  4

// parameters for Lennard-Jones 12-6, with eps and sigma
// V = 4 * eps * [(sigma/r)^12 - (sigma/r)^6]
#define _SR_LJ_      3
#define _N_LJ_       2
#define _LJ_EPS_     0
#define _LJ_SIGMA_   1

// parameters for alternative format of Lennard-Jones 12-6, with eps and sigma
// V = 4 * eps * [(sigma/r)^12 - 2 * (sigma/r)^6]
#define _SR_ALJ_      4
#define _N_ALJ_       2
#define _ALJ_EPS_     0
#define _ALJ_SIGMA_   1

// parameters for Lennard-Jones 12-6, with A and B
// V = A /r^12 - B / r^6
#define _SR_GLJ_      5
#define _N_GLJ_       2
#define _GLJ_EPS_     0
#define _GLJ_SIGMA_   1


#define _SR_FILE_     6


// dipole moment on atom ? 0 / 1
#define _ATOM_DIPOLE   0
// induced dipole moment on atom ? 0 / 1
#define _ATOM_INDUCED_DIPOLE   0

#define A14_E      0.8333333333
#define A14_LJ     0.5

#define _KAPPA     3.2
#if _EWALD_SUM == _MultiPole_EWALD_SUM
// ****** Multipole EWALD ELECTROSTATIC SUM ******
#define _ALPHA_EWALD    0.8   // rcut = _ALPHA * L (size of cell)
#define _NH         4     // recipracal space summation, {NH, NH, NH}, NH ~ 3.2 * L / rcut
#define  NH         9     // 2 * NH + 1 
#define  NPC        1     // period cell repeatation

#define AH_MIN      5e-3  // 1e-9  // minumn cut for Ah -- AhMin
// ****** end of Multipole Ewald Sum ******
#endif

// ************* implicit solvation parameters **************
#define _IMPLICIT_SOLVENT_     0

#if _IMPLICIT_SOLVENT_ == 0
#define _ATOMIC_IMPLICIT_SOLVATION_    0
#define _CLUSTER_IMPLICIT_SOLVATION_   0
#define _CG_CLUSTER_IMPLICIT_SOLVATION_  0
#else
#define _ATOMIC_IMPLICIT_SOLVATION_    0
#define _CLUSTER_IMPLICIT_SOLVATION_   1
#define _CG_CLUSTER_IMPLICIT_SOLVATION_  0
#endif
// ************* end of implicit solvation parameters **************
