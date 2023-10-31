#ifndef _NO_ADDITIONAL_HEAD_MM_
#include "def.h"
#include "dstruct.h"
#include "prof.h"
#include "atom.h"

#include "amber_pars.h"
#endif

#ifndef _NFS
	#define _NFS    2
#endif

extern bool bMTForce;

class DIHEDRAL_DATABASE {
public:
	CHAIN<TORSION_PARS> *ch_dh;
	ARRAY<TORSION_PARS*> dha;
	TORSION_PARS* get_torsion_ch(short a1_indx, short a2_indx, short a3_indx, short a4_indx) {
		CHAIN<TORSION_PARS> *pch = ch_dh;
		while (pch != NULL && !(pch->p->aindx[0] == a1_indx && pch->p->aindx[1] == a2_indx && pch->p->aindx[2] == a3_indx && pch->p->aindx[3] == a4_indx)) pch = pch->next;
		if (pch == NULL) return NULL;
		else return pch->p;
	};
	TORSION_PARS* get_torsion(short a1_indx, short a2_indx, short a3_indx, short a4_indx) {
		TORSION_PARS *tpar = NULL;
		for (int i = 0; i < dha.n; i++) {
			tpar = dha.m[i];
			if (tpar->aindx[0] == a1_indx && tpar->aindx[1] == a2_indx && tpar->aindx[2] == a3_indx && tpar->aindx[3] == a4_indx) return tpar;
		}
		return NULL;
	};
	void attach(TORSION_PARS *pTorsion) {
		CHAIN<TORSION_PARS> *pch = ch_dh;
		if (ch_dh == NULL) {
			ch_dh = new CHAIN<TORSION_PARS>;
			ch_dh->p = pTorsion;
		}
		else {
			while (pch->next != NULL) pch = pch->next;
			pch->next = new CHAIN<TORSION_PARS>; pch->next->p = pTorsion;
		}
	};
	void release() {
		int i;
		if (ch_dh == NULL) release_chain<TORSION_PARS>(&ch_dh, true);
		if (dha.m != NULL) {
			for (i = 0; i < dha.n; i++) {
				if (dha.m[i] != NULL) {delete dha.m[i]; dha.m[i] = NULL;}
			}
			dha.release();
		}
	};
	void convert_array() {
		chain_array<TORSION_PARS>(&ch_dh, dha);
	};
	DIHEDRAL_DATABASE() {release();};
};

extern char errmsg[2560];

class LJ_PARS { // L-J interactions parameters
public:
	float epsLJ, rLJ, rLJ2; // rLJ2 = rLJ^2
	float rcut, rcut2; // cut-off distance
	LJ_PARS() {
		epsLJ = 0; rLJ = 100; rLJ2 = 10000;
		rcut = 1; rcut2 = 1;
	};
	void set_pars(float eps, float r0, float rc) {
		this->epsLJ = eps;
		rLJ = r0; rLJ2 = r0 * r0; rcut = rc; rcut2 = rc * rc;
	};
};

// implicit solvation parameter
class IMPSOLV_PARS {
public:
	// sigma -- solvation energy in kT, f_sigma -- force due to sigma, with unit kT/U_InertiaMoment
	// r0 is the radius -- monomer / atom ?
	float sigma, f_sigma, r0, S0; // S0 = 4 * PI * (r0 + rw)^2
};

#define NORMAL_BOUND  0x00  
#define PI_BOUND      0x01

template <class DERIVED_ATOM> class _ATOM {
public:
	short nb;
	//_ATOM **ba; // bound atoms
	DERIVED_ATOM **ba; // bound atoms
	char *bv; // bound valance
	char *btype; // bound type

	short aindx; // atom index in the atompar database -- used to find out the type of atom
	//char name[4]; // alias
	// r is the external cordinate, r0 is the inertial acordinate
	// r0 is relative to Om- of parent
	VECTOR3 r, r0;
	VECTOR3 rg; // the real position in the simulation system

	_ATOM() {
		nb = 0; ba = NULL; bv = NULL; btype = NULL;
		aindx = '\0';
		// strcpy(name, "\0"); 
	};
	void reset_bounds() {
		int i = 0;
		if (ba != NULL) {
			//for (i = 0; i < nb; i++) {ba[i] = NULL; bv[i] = 0x00; btype[i] = 0x00;}
			//memset(ba, 0, nb * sizeof(_ATOM*)); 
			memset(ba, 0, nb * sizeof(DERIVED_ATOM*)); 
			memset(bv, 0, nb * sizeof(char)); memset(btype, 0, nb * sizeof(char));
			delete[] ba; delete[] bv; delete[] btype;
			ba = NULL; bv = NULL; btype = NULL;
		}
		nb = 0;
	};
	void set_bounds(int n) {
		reset_bounds();
		int i = 0;
		if (n > 0) {
			nb = n;
			ba = new DERIVED_ATOM*[n]; //new _ATOM*[n];
			bv = new char[n]; btype = new char[n];
			//for (i = 0; i < nb; i++) {ba[i] = NULL; bv[i] = 0x00; btype[i] = 0x00;}
			//memset(ba, 0, nb * sizeof(_ATOM*)); 
			memset(ba, 0, nb * sizeof(DERIVED_ATOM*)); 
			memset(bv, 0, nb * sizeof(char)); memset(btype, 0, nb * sizeof(char));
		}
	};
	~_ATOM() {
		reset_bounds();
	};
	//bool is_bound(_ATOM* p) {
	bool is_bound(DERIVED_ATOM* p) {
		int i = 0; 
		for (i = 0; i < nb; i++) if (ba[i] == p) return true;
		return false;
	};
	//int ibound(_ATOM* p) {
	int ibound(DERIVED_ATOM* p) {
		int i = 0; 
		for (i = 0; i < nb; i++) if (ba[i] == p) return i;
		return -1;
	};
	int ibound(char bound) {
		int i = 0;
		for (i = 0; i < nb; i++) {
			if (bv[i] == bound) return i;
			//if (btype[i] == bound) return i;
		}
		return -1;
	};
	int get_free_bound() {
		int i = 0;
		for (i = 0; i < nb; i++) if (ba[i] == NULL) return i;
		return -1;
	};
	virtual void rotate1(R &r) {};
};

template <class T> class _BOUND_ {
public:
	T* a[2];
	_BOUND_() {a[0] = NULL; a[1] = NULL;};
	~_BOUND_() {a[0] = NULL; a[1] = NULL;};
};

template <class T, class BOUND> bool exist_bound(T* a1, T* a2, bool pairwise, BOUND *bound, int nBound) {
	for (int i = 0; i < nBound; i++) {
		if (bound[i].a[0] == a1 && bound[i].a[1] == a2) return true;
		if (pairwise) {if (bound[i].a[1] == a1 && bound[i].a[0] == a2) return true;}
	}
	return false;
};

#define BOUND(a1, a2, i, status) status = false; for (i = 0; i < a1->nb; i++) {\
	if (a1->ba[i] == a2) {status = true; break;}}

#define BOUND_SAME_ATOM(a1, a2, i, j, status) status = false; \
	for (i = 0; i < a1->nb; i++) {\
		for (j = 0; j < a2->nb; j++) {\
			if (a1->ba[i] == a2->ba[j]) {status = true; break;} \
		}\
		if (status) break;\
	}

#define BOUND_CLOSE_ATOMS(a1, a2, i, j, k, status) status = false; \
	for (i = 0; i < a1->nb; i++) {\
		for (j = 0; j < a2->nb; j++) {\
			BOUND(a1->ba[i], a2->ba[j], k, status) \
			if (status) break; \
		}\
		if (status) break;\
	}

//bool bound_connect(MATOM* p1, short ibound1, MATOM* p2, short ibound2, char bound_valence, char bound_type) {
// atom are derived class of _ATOM
template <class atom> bool bound_connect(atom* p1, short ibound1, atom* p2, short ibound2, char bound_valence = SINGLE_BOUND, char bound_type = NORMAL_BOUND) {
	p1->ba[ibound1] = p2; p1->bv[ibound1] = bound_valence; p1->btype[ibound1] = bound_type; 
	p2->ba[ibound2] = p1; p2->bv[ibound2] = bound_valence; p2->btype[ibound2] = bound_type;
	return true;
};

//bool bound_connect(MATOM* p1, MATOM*p2, char bound_valence, char bound_type) {
// atom are derived class of _ATOM
template <class atom> bool bound_connect(atom* p1, atom* p2, char bound_valence = SINGLE_BOUND, char bound_type = NORMAL_BOUND) {
	int i, j;
	for (i = 0; i < p1->nb; i++) {if (p1->ba[i] == NULL) break;}
	if (i >= p1->nb) {
		sprintf(errmsg, "atom 1 has no free bound to connect another atom");
		return false;
	}
	for (j = 0; j < p2->nb; j++) {if (p2->ba[j] == NULL) break;}
	if (j >= p2->nb) {
		sprintf(errmsg, "atom 2 has no free bound to connect another atom");
		return false;
	}
	return bound_connect<atom>(p1, i, p2, j, bound_valence, bound_type);
};


class DIPOLE_HIST {
public:
	SVECTOR<double, 3> mu[3];

	DIPOLE_HIST() {
		memset(mu[0].v, 0, SIZE_V3); memset(mu[1].v, 0, SIZE_V3); memset(mu[2].v, 0, SIZE_V3);
	};
};

class MATOM : public _ATOM<MATOM> {
public:
	//MATOM **ba; // bound atoms
	bool eNeutral;

	float c; // charge -- normalized with 1/sqrt(eps), is not a common parameter from AMBER force-field parameter
	float c0; // charge 
	ATOM_COMM_PARS *par;

	//IMPSOLV_PARS *psolv;  -- implicit solvation model for each atom

	// ****** parameters for electrostatic field ******
	short mMP; // multipole indx
	bool bIMu; // with intrinsic dipole
	VECTOR3 intr_mu; //intrinsic dipole
	VECTOR3 ind_mu; // induced dipole moment
	VECTOR3 mu; // total dipole

	DIPOLE_HIST* mu_hist;

	VECTOR3 E; // electric field on this atom
	// ****** end of electrostatic field ******

	VECTOR3 F; // the force on this atom

	ARRAY<VECTOR3> f;

	int aIndx; // indexed in the whole system, actually the index in atom of MMOL_MD_CELL

	bool bShocked;

	MATOM() {
		par = NULL;
		eNeutral = true;
		c = 0; c0 = 0;
		mMP = 0; // charge only
		bIMu = false;
		mu_hist = NULL;
		//psolv = NULL;
		bShocked = false;
	};
	~MATOM() {
		par = NULL; _ATOM<MATOM>::reset_bounds(); mu_hist = NULL;
		//psolv = NULL;
	};
	void rotate1(R &r) { // rotate something of atom
		// e.g. intrinsic dipole
		VECTOR3 v;
		if (mMP > 0 && bIMu) {
			MpV3(r, intr_mu, v)
			V32V3(v, intr_mu)
		}
	};
	/*
	bool is_bound(MATOM* p) {
		return _ATOM::is_bound((_ATOM*)p);
	};
	int ibound(MATOM* p) {
		return _ATOM::ibound((_ATOM*)p);
	};
	*/
};

void cp(MATOM *d, MATOM *s);


class CLUSTER_IMPSOLV_PARS : public IMPSOLV_PARS { // all these parameter are site indepedent
public:
	short id; // index in current NEIMO run
	short cID; // cluster Index -- defined in CLUSTER.dat 
	CLUSTER_IMPSOLV_PARS() {
		id = 0;
		cID = 0;
	};
};

class CLUSTER_IMPSOLV_PARS_DATABASE {
public:
	CHAIN<CLUSTER_IMPSOLV_PARS> *ch;

	CLUSTER_IMPSOLV_PARS_DATABASE() {ch = NULL;};
	void release() {release_chain<CLUSTER_IMPSOLV_PARS>(&ch, true);};
	~CLUSTER_IMPSOLV_PARS_DATABASE() {release();};

	CLUSTER_IMPSOLV_PARS* attach(short cID);
	CLUSTER_IMPSOLV_PARS* search_par(short cID) {
		CHAIN<CLUSTER_IMPSOLV_PARS> *pch = ch;
		CLUSTER_IMPSOLV_PARS *par = NULL;
		while (pch != NULL) {
			par = pch->p;
			if (par->cID == cID) return par;
			else pch = pch->next;
		}
		return NULL;
	};
};

template <class T> class HINGE {
public:
	short int indx, ibound; // hinge related atom index and bound index in the cluster

	VECTOR3 *O, *Olink, vHinge; // adjoint point and the hinge vector
	T *plink; // the cluster connecting to the hinge
	short indx_link; // the index of the hinge in connected cluster
	
	HINGE() {
		indx = -1; ibound = -1;
		plink = NULL; indx_link = -1; O = NULL; Olink = NULL;
	};
	~HINGE() {plink = NULL; O = NULL; Olink = NULL;};
};

class BASIC_CLUSTER;

// ************* these were defined in CMM.h, using for CMM multipole Ewald-Sum calculation ***************
#if _EWALD_SUM == _MultiPole_EWALD_SUM
class ELECTRO_MULTIPOLE {
public:
	float q; // total charge
	VECTOR3 r0; // charge center position
	VECTOR3 mu; // dipole
	DMATRIX<3> Q; // quadrupole
	//CMATRIX<3> O; // Octopole
	// for Ewald Sum
	CMATRIX<NH> fsin, fcos;
	CMATRIX<NH> muh, Qh;
	
	float d; // size of the multipole

	ELECTRO_MULTIPOLE() {q = 0; d = 0.1f;};
	bool neutral() {
		if (q != 0 || mu.v[0] != 0 || mu.v[1] != 0 || mu.v[2] != 0) return false;
		int i, j;
		for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) if (Q.m[i][j] != 0) return false;
		return true;
	};
};

class CMM_RELSHIP {
private:
	CHAIN<ELECTRO_MULTIPOLE> *tail_mpl;
	CHAIN<BASIC_CLUSTER> *tail_cluster;
public:
	CHAIN<ELECTRO_MULTIPOLE> *ch_mpl;
	CHAIN<BASIC_CLUSTER> *ch_cluster;
	CMM_RELSHIP() {ch_mpl = NULL; tail_mpl = NULL; ch_cluster = NULL; tail_cluster = NULL;};
	void release_relationship() {
		if (ch_mpl != NULL) release_chain<ELECTRO_MULTIPOLE>(&ch_mpl, false); tail_mpl = NULL;
		if (ch_cluster != NULL) release_chain<BASIC_CLUSTER>(&ch_cluster, false); tail_cluster = NULL;
	};
	~CMM_RELSHIP() {release_relationship();};
	void attach(ELECTRO_MULTIPOLE *emp) {
		CHAIN<ELECTRO_MULTIPOLE> *pch = new CHAIN<ELECTRO_MULTIPOLE>; pch->p = emp;
		if (ch_mpl == NULL) {ch_mpl = pch; tail_mpl = pch;}
		else {tail_mpl->next = pch; tail_mpl = pch;}//ch_mpl->attach2tail(emp);
	};
	void attach(BASIC_CLUSTER *pc) {
		CHAIN<BASIC_CLUSTER> *pch = new CHAIN<BASIC_CLUSTER>; pch->p = pc;
		if (ch_cluster == NULL) {ch_cluster = pch; tail_cluster = pch;}
		else {tail_cluster->next = pch; tail_cluster = pch;} //ch_cluster->attach2tail(pc);
	};
};

#endif

template <class T> class CMM_IMAGINE {
public:
	T *p;
	int nx, ny, nz; // imagine cell index
	bool local;
	CMM_IMAGINE() {p = NULL; nx = 0; ny = 0; nz = 0; local = true;};
	CMM_IMAGINE(T *p, int nx, int ny, int nz) {
		this->p = p; this->nx = nx; this->ny = ny; this->nz = nz;
		local = ((nx == 0 && ny == 0 && nz == 0) ? true : false);
	};
	~CMM_IMAGINE() {p = NULL; nx = 0; ny = 0; nz = 0;};
};

template <class T> class FLEX_BUFF {
public:
	CHAIN<T> *ch;
	CHAIN<T> *tail;
	//CHAIN<T> *cur_ch; -- for iteration

	ARRAY<T*> ba; // relationship array
	int nbs; // number of relationships used
	//int icur; // current relationship in array -- for iteration
	void expand_array(float c) { // c is the relative increment
		int n = ba.n, dn = int(n * c);
		if (dn < 10) dn = 10;
		T** bf = ba.m; ba.m = NULL;
		ba.SetArray(n + dn);
		/*
		for (i = 0; i < n; i++) ba.m[i] = bf[i];
		for (i = n; i < ba.n; i++) ba.m[i] = NULL;
		*/
		memcpy(ba.m, bf, n * sizeof(T*));
		delete[] bf;
	};
	void require_memory(int n) { // make sure enough memory is avaliable -- n
		int i, dn = ba.n - nbs, n0 = ba.n;
		if (dn >= n) return;
		T** bf = ba.m; 
		ba.m = NULL; ba.n = 0;
		ba.SetArray(nbs + n);
		memcpy(ba.m, bf, n0 * sizeof(T*));
		/*
		for (i = 0; i < n; i++) {
			ba.m[i].p = bf[i].p;
		}
		*/
		delete[] bf;
	};
public:
	bool bChain;
	FLEX_BUFF() {
		ch = NULL; tail = NULL; nbs = 0; bChain = true; // default -- use chain
		//cur_ch = NULL; icur = 0; 
	};
	void release(bool reset = false) {
	#if _SYS_ == _WINDOWS_SYS_
		if (ch != NULL) ::release_chain<T>(&ch, false);
	#elif _SYS_ == _LINUX_SYS_
		if (ch != NULL) ::release_chain<T>(&ch, false);
	#endif
		tail = NULL;
		//int i;
		if (reset) {
			//for (i = 0; i < ba.n; i++) ba.m[i] = NULL;
			memset(ba.m, 0, ba.n * sizeof(T*));
		}
		nbs = 0;
		//cur_ch = NULL; icur = -1; 
	};
	void use_array(int nmax) {
		bChain = false; ba.SetArray(nmax); nbs = 0;
		//for (int i = 0; i < nmax; i++) ba.m[i] = NULL;
		memset(ba.m, 0, ba.n * sizeof(T*));
	};
	CHAIN<T> *get_chain_head() {return (bChain? ch : NULL);};
	T** get_head() {return (bChain ? NULL : ba.m);}
	/*
	void restart_iterate() {
		cur_ch = ch; icur = 0;
	};
	void next_iterate() {
		if (bChain) {
			if (cur_ch == NULL) return;
			else cur_ch = cur_ch->next;
		}
		else {
			if (icur < nbs) icur++;
			if (icur >= nbs) icur = -1;
		}
	};
	T* current() {
		if (bChain) {
			if (cur_ch == NULL) return NULL;
			else return cur_ch->p;
		}
		else {
			if (icur < 0 || icur >= nbs) return NULL;
			else return ba.m[icur];
		}
	};
	*/
	~FLEX_BUFF() {release();};
	bool attach(T *p) {
		CHAIN<T> *pch = NULL;
		if (bChain) {
			pch = new CHAIN<T>; pch->p = p;
			if (ch == NULL) {ch = pch; tail = pch;}
			else {tail->next = pch; tail = pch;} //ch->attach2tail(p);
			nbs++; return true;
		}
		else {
			if (nbs >= ba.n) {// expand the array about 20%
				//show_log("relationship is expanded", true); return false;
				expand_array(0.2);
			}
			ba.m[nbs] = p; nbs++; return true;
		}
	};
	bool detach(T *p) {
		bool btail, status = true;
		int indx = 0;
		if (bChain) {
			btail = (tail != NULL && tail->p == p ? true : false);
			status = ::detach_from_chain<T>(&ch, p, false);
			if (btail) {if (ch == NULL) tail = NULL; else {tail = ch; while (tail->next != NULL) tail = tail->next;}}
			if (status) nbs--; return status;
		}
		else {
			for (indx = 0; indx < nbs; indx++) {if (ba.m[indx] == p) break;}
			if (indx < nbs) {
				nbs--; ba.m[indx] = ba.m[nbs]; return true;
			} // replace indx with the last ca
			else return false;
		}
	};
	int length() {
		if (bChain) return number<T>(ch);
		else return nbs;
		//return nbs;
	};
};

template <class T, class BUFF> class ITERATOR {
protected:
	BUFF *buff;
	CHAIN<T> *cur_ch;
	int icur;
public:
	ITERATOR() {buff = NULL; cur_ch = NULL; icur = -1;};
	ITERATOR(BUFF *buff) {this->buff = buff; cur_ch = NULL; icur = -1;};
	void set(BUFF* buff) {this->buff = buff;};
	void start_iterate() {
		if (buff->bChain) cur_ch = buff->get_chain_head(); 
		else icur = 0;
	};
	void next_iterate() {
		if (buff->bChain) {
			if (cur_ch == NULL) return;
			else cur_ch = cur_ch->next;
		}
		else {
			if (icur < buff->length()) icur++;
			if (icur >= buff->length()) icur = -1;
		}
	};
	T* current() {
		if (buff->bChain) {
			if (cur_ch == NULL) return NULL;
			else return cur_ch->p;
		}
		else {
			if (icur < 0 || icur >= buff->length()) return NULL;
			else return *(buff->get_head() + icur);
		}
	};
	bool eoi() {
		if (buff->bChain) return (cur_ch == NULL ? true : false);
		else return ((icur < 0 || icur >= buff->length()) ? true : false);
	};
};

template <class T> class IMAG_RELSHIP {
public:
	CHAIN< CMM_IMAGINE<T> > *ch;
	CHAIN< CMM_IMAGINE<T> > *tail;
	//CHAIN< CMM_IMAGINE<T> > *cur_ch;

	ARRAY< CMM_IMAGINE<T> > ba; // relationship array
	int nbs; // number of relationships used
	//int icur; // current relationship in array

	void expand_array(float c) { // c is the relative increment
		int i, n = ba.n, dn = int(n * c);
		if (dn < 10) dn = 10;
		CMM_IMAGINE<T>* bf = ba.m; 
		ba.m = NULL; ba.n = 0;
		ba.SetArray(n + dn);
		/*
		for (i = 0; i < n; i++) {
			ba.m[i].p = bf[i].p;
			ba.m[i].nx = bf[i].nx; ba.m[i].ny = bf[i].ny; ba.m[i].nz = bf[i].nz;
			ba.m[i].local = bf[i].local;
		}
		*/
		memcpy(ba.m, bf, n * sizeof(CMM_IMAGINE<T>));
		delete[] bf;
	};
	void require_memory(int n) { // make sure enough memory is avaliable -- n
		int i, dn = ba.n - nbs, n0 = ba.n;
		if (dn >= n) return;
		CMM_IMAGINE<T>* bf = ba.m; 
		ba.m = NULL; ba.n = 0;
		ba.SetArray(nbs + n);
		memcpy(ba.m, bf, n0 * sizeof(CMM_IMAGINE<T>));
		/*
		for (i = 0; i < n; i++) {
			ba.m[i].p = bf[i].p;
			ba.m[i].nx = bf[i].nx; ba.m[i].ny = bf[i].ny; ba.m[i].nz = bf[i].nz;
			ba.m[i].local = bf[i].local;
		}
		*/
		delete[] bf;
	};
public:
	bool bChain;
	IMAG_RELSHIP() {
		ch = NULL; tail = NULL; nbs = 0; bChain = true; // default -- use chain
		//cur_ch = NULL; icur = 0; 
	};
	void release_relationship(bool reset = false) {
	#if _SYS_ == _WINDOWS_SYS_
		if (ch != NULL) ::release_chain< CMM_IMAGINE<T> >(&ch, true);
	#elif _SYS_ == _LINUX_SYS_
		if (ch != NULL) ::release_chain< CMM_IMAGINE<T> >(&ch, true);
	#endif
		tail = NULL; 
		int i;
		if (reset) {
			for (i = 0; i < ba.n; i++) ba.m[i].p = NULL;
		}
		nbs = 0;
		//cur_ch = NULL; icur = -1; 
	};
	void use_array(int nmax) {
		bChain = false; ba.SetArray(nmax); nbs = 0;
	};
	CHAIN<CMM_IMAGINE<T> > *get_chain_head() {return (bChain ? ch : NULL);};
	CMM_IMAGINE<T>* get_head() {return (bChain ? NULL : ba.m);}
	/*
	void restart_iterate() {
		cur_ch = ch; icur = 0;
	};
	void next_iterate() {
		if (bChain) {if (cur_ch == NULL) return;}
		else {
			if (icur < nrs) icur++;
			if (icur >= nrs) icur = -1;
		}
	};
	CMM_IMAGINE<T>* current_relationship() {
		if (bChain) {
			if (cur_ch == NULL) return NULL;
			else return cur_ch->p;
		}
		else {
			if (icur < 0 || icur >= nrs) return NULL;
			else return ra.m + icur;
		}
	};
	*/
	~IMAG_RELSHIP() {release_relationship(); ba.release();};
	bool attach(T *pc, int nx, int ny, int nz) {
		CMM_IMAGINE<T> *p = NULL;
		CHAIN< CMM_IMAGINE<T> > *pch = NULL;
		bool status = true;
		if (bChain) {
			p = new CMM_IMAGINE<T>(pc, nx, ny, nz);
			pch = new CHAIN< CMM_IMAGINE<T> >; pch->p = p;
			if (ch == NULL) {ch = pch; tail = pch;}
			else {tail->next = pch; tail = pch;} //ch->attach2tail(p);
			nbs++;
		}
		else {
			if (nbs >= ba.n) {//we have to expand the array buffer
				//show_log("relationship is leaking", true); status = false;
				expand_array(0.2); // 20%
			}
			ba.m[nbs].p = pc; ba.m[nbs].nx = nx; ba.m[nbs].ny = ny; ba.m[nbs].nz = nz;
			ba.m[nbs].local = ((nx == 0 && ny == 0 && nz == 0) ? true : false);
			nbs++;
		}
		return status;
	};
	CMM_IMAGINE<T> *get_rship(T* pc) {
		CHAIN< CMM_IMAGINE<T> > *pch = NULL;
		int ich = 0;
		if (bChain) {
			pch = ch;
			while (pch != NULL) {
				if (pch->p->p == pc) return pch->p;
				pch = pch->next;
			}
			return NULL;
		}
		else {
			for (ich = 0; ich < nbs; ich++) {
				if (ba.m[ich].p == pc) return ba.m + ich;
			}
			return NULL;
		}
	};
	int nrship() {
		if (bChain) return number< CMM_IMAGINE<T> >(ch);
		else return nbs;
		//return nbs;
	};
	int length() {
		if (bChain) return number< CMM_IMAGINE<T> >(ch);
		else return nbs;
		//return nbs;
	};
};

template <class T, class RELSHIP> class RELSHIP_ITERATOR {
protected:
	RELSHIP *rship;
	CHAIN< CMM_IMAGINE<T> > *cur_ch;
	int icur;
public:
	RELSHIP_ITERATOR() {rship = NULL; cur_ch = NULL; icur = -1;};
	RELSHIP_ITERATOR(RELSHIP *rship) {this->rship = rship; cur_ch = NULL; icur = -1;};
	void set(RELSHIP* rship) {this->rship = rship;};
	void start_iterate() {
		if (rship->bChain) cur_ch = rship->get_chain_head(); 
		else icur = 0;
	};
	void next_iterate() {
		if (rship->bChain) {
			if (cur_ch == NULL) return;
			else cur_ch = cur_ch->next;
		}
		else {
			if (icur < rship->length()) icur++;
			if (icur >= rship->length()) icur = -1;
		}
	};
	CMM_IMAGINE<T>* current() {
		if (rship->bChain) {
			if (cur_ch == NULL) return NULL;
			else return cur_ch->p;
		}
		else {
			if (icur < 0 || icur >= rship->length()) return NULL;
			else return rship->get_head() + icur;
		}
	};
	bool eoi() {
		if (rship->bChain) return (cur_ch == NULL ? true : false);
		else return ((icur < 0 || icur >= rship->length()) ? true : false);
	};
};


/*
template <class T> CMM_IMAGINE<T>* in_imagine_rship(IMAG_RELSHIP<T> *rship, T *pc) {
	CHAIN< CMM_IMAGINE<T> > *ch = rship_chain;
	while (ch != NULL) {
		if (ch->p->p == pc) return ch->p;
		ch = ch->next;
	}
	return NULL;
};
*/

class LJ_RELSHIP : public IMAG_RELSHIP<BASIC_CLUSTER> {
public:
	LJ_RELSHIP() {};
	~LJ_RELSHIP() {};
};

class SPME_RELSHIP : public IMAG_RELSHIP<BASIC_CLUSTER> {
public:
	SPME_RELSHIP() {};
	~SPME_RELSHIP() {};
};

template <class CLUSTER> class IMPSOLV_RELSHIP : public IMAG_RELSHIP<CLUSTER> {
public:
	IMPSOLV_RELSHIP() {};
	~IMPSOLV_RELSHIP() {};
};

// **************************************** end *******************************************
//******************************************************************************************
//    FOLLOWING WAY TO DESCRIBE THE TORSION BETWEEN CLUSTER (OR ON HINGE) IS WRONG
//******************************************************************************************
/*
class TORSION {
public:
	short pID, cID; // cluster indexes of parent and child clusters, of this torsion hinge
	EF_DPROF tProf; // torsion profile of this hinge
	//DPROF torqueProf; // torque profile of this hinge
	TORSION() {pID = -1; cID = -1;};
	TORSION(short pID, short cID) {this->pID = pID; this->cID = cID;};
	void release() {tProf.release();};
	~TORSION() {release();};
	void get_torsion(double x, double &E, double &torque) {
		int i = 0;
		double t;
		PERIOD_RANGE(x, -PI, PI, PI2)
		//interpolate(E, tProf, x, i)
		//interpolate(torque, torqueProf, x, i)
		interpolate2(E, torque, tProf, x, i, t)
	};
};

class TORSION_DATABASE {
public:
	CHAIN<TORSION> *ch_torsion;
	TORSION* get_torsion(short pID, short cID) {
		CHAIN<TORSION> *pch = ch_torsion;
		while (pch != NULL && !(pch->p->cID == cID && pch->p->pID == pID)) {
			pch = pch->next;
		}
		if (pch == NULL) return NULL;
		else return pch->p;
	};
	void attach(TORSION *pTorsion) {
		CHAIN<TORSION> *pch = ch_torsion;
		if (ch_torsion == NULL) {
			ch_torsion = new CHAIN<TORSION>;
			ch_torsion->p = pTorsion;
		}
		else {
			while (pch->next != NULL) pch = pch->next;
			pch->next = new CHAIN<TORSION>; pch->next->p = pTorsion;
		}
	};
	void release() {release_chain<TORSION>(&ch_torsion, true);};
	TORSION_DATABASE() {ch_torsion = NULL;};
};
*/

// TATOM has to be derived class from _ATOM
template <class TATOM> class _BASIC_CLUSTER {
public:
	int nAtoms;
	TATOM *atom;

	int nHinges;
	HINGE<_BASIC_CLUSTER<TATOM> > *hinge;
	VECTOR3 *Op; // the position where parent hinge connecting to
	VECTOR3 *vcm; // mass center
	_BASIC_CLUSTER<TATOM> *parent;

	short n2Parent; // the hinge connecting to parent

	bool bCMMShift; // dr_cmm is non-zero ?
	VECTOR3 dr_cmm; // periodical shift in the CMM

	VECTOR3 *vp; // point to a vector for cluster position using in CMM cell 

	short cTypeIndx; // index of the type of cluster in the whole system, used to find out the type of cluster
	                 // actually, this number is indexed in global variable cluster_db

	_BASIC_CLUSTER() {
		n2Parent = -1; nAtoms = 0; atom = NULL; nHinges = 0; cTypeIndx = 0;
		hinge = NULL; Op = NULL; parent = NULL; vp = NULL;
		bCMMShift = false; vcm = NULL;
	};
	void resetAtoms() {
		int i = 0;
		if (atom != NULL) {delete[] atom; atom = NULL;}; nAtoms = 0;
	};
	void setup_hinge(int iHinge, int iatom, short ibound) {
		hinge[iHinge].indx = iatom; hinge[iHinge].ibound = ibound;
		hinge[iHinge].O = &(atom[iatom].r);
	};
	void set_parent_hinge(int iHinge) {this->n2Parent = iHinge;};
	
	void SetupParent(int natom_default = 0) {
		if (n2Parent >= 0 && hinge[n2Parent].Olink != NULL) {
			Op = hinge[n2Parent].Olink;
			parent = hinge[n2Parent].plink;
		}
		else {Op = &(atom[natom_default].r); parent = NULL;}
	};
	void disable_parent(int natom_default = 0) {
		n2Parent = -1; parent = NULL;
		Op = &(atom[natom_default].r);
	};
	
	void reset_hinges() {if (hinge != NULL) {delete[] hinge; hinge = NULL;} nHinges = 0;};
	void reset() {
		resetAtoms(); reset_hinges(); Op = NULL; parent = NULL; 
	};
	~_BASIC_CLUSTER() {
		reset(); Op = NULL; parent = NULL; 
	};

	void setAtoms(int nAtoms) {
		resetAtoms();
		this->nAtoms = nAtoms;
		if (nAtoms > 0) atom = new TATOM[nAtoms];
	};
	void set_hinges(int n) {
		reset_hinges();
		if (n > 0) {
			nHinges = n;
			hinge = new HINGE<_BASIC_CLUSTER<TATOM> >[n];
		}
	};
	int iAtom(TATOM *p) {
		int i = 0;
		for (i = 0; i < nAtoms; i++) if (atom + i == p) return i;
		return -1;
	};

	bool setup_link_relationship(int iHinge, _BASIC_CLUSTER<TATOM> *c, int jHinge, char bound_valence = SINGLE_BOUND) {
		hinge[iHinge].plink = c; c->hinge[jHinge].plink = (_BASIC_CLUSTER<TATOM> *)this;
		hinge[iHinge].indx_link = jHinge; c->hinge[jHinge].indx_link = iHinge;
		hinge[iHinge].Olink = c->hinge[jHinge].O;
		c->hinge[jHinge].Olink = hinge[iHinge].O;
		
		return bound_connect<TATOM>(&(atom[hinge[iHinge].indx]), hinge[iHinge].ibound, c->atom + c->hinge[jHinge].indx, c->hinge[jHinge].ibound, bound_valence, NORMAL_BOUND);
	};
	/*
	bool setup_child_relationship(int iHinge, CLUSTER *c, char bound_type = SINGLE_BOUND) {
		if (c->n2Parent < 0) return false; // c does not have parent hinge
		return setup_link_relationship(iHinge, c, c->n2Parent, bound_type);
	};
	*/
	void shift(VECTOR3 &d) {
		int i = 0, j;
		VECTOR3 dr = d;
		for (i = 0; i < nAtoms; i++) {for (j = 0; j < 3; j++) atom[i].r.v[j] += dr.v[j];}
	};
	void rotate(R &r) { // given rotation matrix
		int i = 0;
		VECTOR3 v;
		for (i = 0; i < nAtoms; i++) {
			v = atom[i].r; // atom r
			atom[i].r = r * v;
			atom[i].rotate(r);
		}
		// hinges
		for (i = 0; i < nHinges; i++) {
			v = hinge[i].vHinge; // hinge vector
			hinge[i].vHinge = r * v;
		}
	};
	void rotate(VECTOR3 &r1, VECTOR3 &r2, double omega, bool arc) { // r1 is base point, omega in degree ?
		if (FABS(omega) < MIN_ROT) return;
		int i;
		VECTOR3 r01, r02;
		for (i = 0; i < 3; i++) {r01.v[i] = r1.v[i]; r02.v[i] = r2.v[i];}
		VECTOR3 v0;
		for (i = 0; i < 3; i++) v0.v[i]= -(r01.v[i]);
		VECTOR3 axis = r02 - r01;
		shift(v0);
		rotate(axis, omega, arc);
		shift(r01);
	};
	void align(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &r0, VECTOR3 &u) { // align r1==>r0 & r12 || u
		VECTOR3 R0 = r0, R1 = r1;

		VECTOR3 v1 = r2 - r1;
		VECTOR3 dr = r1 * (-1);
		shift(dr); // move r1 to 0

		double cos = 0;
		int i = 0;
		scale_uv3(v1, u, cos) // c = v1 * u
		VECTOR3 nu;

		extern void non_parallel(VECTOR3 &v, VECTOR3 &nv);

		// rotate r12 ==> u
		VECTOR3 axis;
		double omega = 0;
		cp(v1, u, axis);
		if (axis.abs() > 1e-7){
			omega = angle(v1, u); // in arc degree
		}
		else { // the two axises are parallel
			if (cos < 0) { // inverse direction
				non_parallel(u, nu); // nu is an axis which is not parallel to u
				cp(u, nu, axis);
				omega = PI;
			}
			else omega = 0;
		}

		R rmatrix;
		RMatrix(axis, omega, rmatrix, true);
		rotate(rmatrix);

		// shift r1 ==> r0
		dr = R0;
		shift(dr);
	};
};

template <class __ATOM> bool cp_basic_cluster(_BASIC_CLUSTER<__ATOM> *d, _BASIC_CLUSTER<__ATOM> *s, void* cp_atom) {
typedef void (*CP_ATOM)(__ATOM*, __ATOM*);
	d->reset();
	if (s->nAtoms <= 0) return true;
	d->nAtoms = s->nAtoms;
	d->atom = new __ATOM[d->nAtoms];
	int i = 0, j = 0, n = 0, nb = 0;
	for (i = 0; i < d->nAtoms; i++) (*(CP_ATOM)cp_atom)(d->atom + i, s->atom + i);
	// forming the bound 
	for (i = 0; i < d->nAtoms; i++) {
		for (j = 0; j < d->atom[i].nb; j++) {
			if (s->atom[i].ba[j] != NULL) {
				n = s->iAtom((__ATOM*)(s->atom[i].ba[j]));
				if (n < 0) { // this bound could be in other cluster
					//sprintf(errmsg, "ERROR : Program is not consistent. Can not find atom defined in cluster");
					//return false;
				}
				else {
					nb = s->atom[n].ibound(s->atom + i);
					if (nb < 0) {
						sprintf(errmsg, "ERROR: Program is not consisten to find out the bound connection in defined cluster");
						return false;
					}
					if (!bound_connect<__ATOM>(d->atom + i, j, d->atom + n, nb, s->atom[i].bv[j], s->atom[i].btype[j])) return false;
				}
			}
		}
	}
	d->set_hinges(s->nHinges);

	for (i = 0; i < s->nHinges; i++) {
		d->hinge[i].vHinge = s->hinge[i].vHinge;
		if (s->hinge[i].indx >= 0) d->setup_hinge(i, s->hinge[i].indx, s->hinge[i].ibound);
		// the connection can not be copied
	}
	d->n2Parent = s->n2Parent;
	d->SetupParent();
	d->cTypeIndx = s->cTypeIndx;

	return true;
};

// looking for all children of the cluster in cch
// CLUSTER is a derived cluster from _BASIC_CLUSTER
template <class CLUSTER> void child_chain(CHAIN<CLUSTER> *cch, CHAIN<CLUSTER> **child) {
	if (child == NULL) return;
	if (*child != NULL) release_chain<CLUSTER>(child, false);

	CHAIN<CLUSTER> *tail = NULL, *pch = NULL, *mch = NULL, *pmch = NULL;
	CLUSTER *pc = NULL;
	int iHinge = 0;

	mch = cch; pmch = mch; *child = NULL; tail = NULL;
	while (pmch != NULL) {
		pc = pmch->p;
		for (iHinge = 0; iHinge < pc->nHinges; iHinge++) {
			if (pc->hinge[iHinge].plink == NULL) continue;
			pch = new CHAIN<CLUSTER>; pch->p = (CLUSTER*)(pc->hinge[iHinge].plink);
			if (*child == NULL) *child = pch;
			else tail->next = pch;
			tail = pch;
		}
		pmch = pmch->next;
	}
};

// looking for nth grand-children of a cluster 
// CLUSTER is a derived cluster from _BASIC_CLUSTER
template <class CLUSTER> void grandchildren_chain(CLUSTER *pcluster, int nth_grand, CHAIN<CLUSTER> **gchild) {
	CHAIN<CLUSTER> *child = NULL, *parent = NULL;
	if (gchild == NULL) return;
	if (*gchild != NULL) release_chain<CLUSTER>(gchild, false);

	CHAIN<CLUSTER> *tail = NULL, *pch = NULL, *mch = NULL, *pmch = NULL;
	CLUSTER *pc = NULL, *pc1 = NULL;
	int iHinge = 0, ith_grand = 0;

	parent = NULL; mch = NULL; child = new CHAIN<CLUSTER>; child->p = pcluster;
	while (ith_grand < nth_grand) {
		release_chain<CLUSTER>(&parent, false);
		parent = mch; mch = child; pmch = mch; child = NULL; tail = NULL;
		while (pmch != NULL) {
			pc = pmch->p;
			for (iHinge = 0; iHinge < pc->nHinges; iHinge++) {
				if (pc->hinge[iHinge].plink == NULL) continue;
				pc1 = (CLUSTER*)(pc->hinge[iHinge].plink);
				if (in_chain<CLUSTER>(parent, pc1)) continue; // this cluster is in parent
				pch = new CHAIN<CLUSTER>; pch->p = pc1;
				if (child == NULL) child = pch;
				else tail->next = pch;
				tail = pch;
			}
			pmch = pmch->next;
		}
		ith_grand++;
	}
	*gchild = child; child = NULL;
	release_chain<CLUSTER>(&parent, false); release_chain<CLUSTER>(&mch, false);
};

// looking for all children and grand-children of a cluster, within ngrands
// CLUSTER is a derived cluster from _BASIC_CLUSTER
template <class CLUSTER> void children_chain(CLUSTER *pcluster, int ngrands, CHAIN<CLUSTER> **children) {
	CHAIN<CLUSTER> *child = NULL, *parent = NULL;
	if (children == NULL) return;
	if (*children != NULL) release_chain<CLUSTER>(children, false);
	CHAIN<CLUSTER> *children_tail = NULL;

	CHAIN<CLUSTER> *tail = NULL, *pch = NULL, *mch = NULL, *pmch = NULL;
	CLUSTER *pc = NULL, *pc1 = NULL;
	int iHinge = 0, ith_grand = 0;

	parent = NULL; mch = NULL; child = new CHAIN<CLUSTER>; child->p = pcluster;
	while (ith_grand < ngrands) {
		release_chain<CLUSTER>(&parent, false);
		parent = mch; mch = child; pmch = mch; child = NULL; tail = NULL;
		while (pmch != NULL) {
			pc = pmch->p;
			for (iHinge = 0; iHinge < pc->nHinges; iHinge++) {
				if (pc->hinge[iHinge].plink == NULL) continue;
				pc1 = (CLUSTER*)(pc->hinge[iHinge].plink);
				if (in_chain<CLUSTER>(parent, pc1)) continue; // this cluster is in parent
				pch = new CHAIN<CLUSTER>; pch->p = pc1;
				if (child == NULL) child = pch;
				else tail->next = pch;
				tail = pch;

				// add pc1 to *children
				pch = new CHAIN<CLUSTER>; pch->p = pc1;
				if (*children == NULL) *children = pch;
				else children_tail->next = pch;
				children_tail = pch;
			}
			pmch = pmch->next;
		}
		ith_grand++;
	}
	release_chain<CLUSTER>(&parent, false); 
	release_chain<CLUSTER>(&mch, false);
	release_chain<CLUSTER>(&child, false); 
};

class IMPSOLV_F {
public:
	float fb, fb1, fbp, fbp1;
	VECTOR3 u;
	IMPSOLV_F() {fb = 0; fb1 = 0; fbp = 0; fbp1 = 0;};
};

class BASIC_KINEMATIC {
public:
	VECTOR<6> P0; // momentum, specifically because of angular momentum
	VECTOR<6> V0; // velocity of cluster in the cordinate of base cluster, or say experimental cordinate
	VECTOR<6> alpha; // acceleration
};

class BASIC_DYN {
public:
	VECTOR<6> fc; // real torque and force acting on the cluster
	bool overlap;
	VECTOR<6> f_Hoover; // Hoover force
	VECTOR<6> f_Brown; // Brownian force

	VECTOR<6> f; // total force on the cluster, which will be used for NEIMO calculation

	BASIC_DYN() {overlap = false;};
};

class BASIC_CLUSTER : public _BASIC_CLUSTER<MATOM> {
public:
	bool deuterized; // is this basic_cluster deuterized ?

	BASIC_KINEMATIC *bkm;
	BASIC_DYN *dyn;

	float d; // the size of the cluster -- rough diameter or dimension, using in Ewald_Sum calculation
	float M; // total mass of cluster
	
	bool eNeutral; // no charge in side at all -- all atoms are neutral
	short cID; // ID/UID of the type of cluster, defined in the database
	char mType; // type of molecule, MM, SM, PM -- 0, 1, 2
	short mIndx; // the molecule index in the whole system, where this cluster is
	int cmIndx; // the cluster index in the macro-molecule
	short cgIndx; // index of coarse-grained cluster where this cluster is in
	short monIndx; // monomer index inside the macro-molecule

#if _EWALD_SUM == _MultiPole_EWALD_SUM
	ELECTRO_MULTIPOLE *emp;
	CMM_RELSHIP ERship; // multipole moment method with CMM
#endif

	// virtual atoms, used for mass/inertia calculation, force calculation ...
	ARRAY<MATOM*> gatom; // reorganize the atom, involving in this cluster
	ARRAY<MATOM*> matom; // the atom will be involved for inertia momentum calculation for this cluster, has non-zero mass
	ARRAY<MATOM*> fatom; // the atom will be involved for force calculation (L-J and/or electrostat) for this cluster
	//***************  end  ****************


	SPME_RELSHIP spmeRship; // SPME atomic interaction, this structure is same as LJ_RELSHIP
	LJ_RELSHIP LJRship;

	//TORSION *pTors;

	bool highlight; // highlight this cluster in VMM, for view only
	char indx_img; // working as image cluster? 0 -- not image cluster, > 0 -- index the region of image, decided by system

	CLUSTER_IMPSOLV_PARS *psolv;
	VECTOR3 rc_solv, rc0_solv; // geometrical center of cluster, in real and local cordinates
	IMPSOLV_RELSHIP<BASIC_CLUSTER> impsolv_rship; // neighbor cluster in chain for solvation consideration
	double E_ImplicitSolvation;
	float SASA; // total solvent accessible surfacial area of the cluster
	ARRAY<IMPSOLV_F> impsolv_fpar; // to collect the structural / force / energy parameters, relating to the calculation of implicit solvenation 

	float r_solve; // exact size of the cluster -- using for force from solvation energy calculation
	float ita_Brownian, cd_trans, cd_rot; // drag force of continum solvent model

	bool bShocked;

	BASIC_CLUSTER() {
		cgIndx = -1; // does not belong to any coarse-grained cluster
		d = 1; eNeutral = true; cID = 0; mType = 0; cmIndx = 0; mIndx = 0; cgIndx = 0; monIndx = 0;
		highlight = false; indx_img = 0; bShocked = false;
		deuterized = false; //pTors = NULL;
		bkm = NULL; dyn = NULL; 
#if _EWALD_SUM == _MultiPole_EWALD_SUM
		emp = NULL;
#endif
		// implicit-solvation variables
		psolv = NULL;
		E_ImplicitSolvation = 0; SASA = 0;

		// Brownian 
		r_solve = 0; ita_Brownian = 0; cd_trans = 0; cd_rot = 0; 
	};
	
	void setup_kin_dyn() {
		if (bkm == NULL) bkm = new BASIC_KINEMATIC();
		if (dyn == NULL) dyn = new BASIC_DYN();
		int i;
		if (::bMTForce) {
			for (i = 0; i < nAtoms; i++) atom[i].f.set_array(_NFS);
		}
	};
	void reset_kin_dyn() {
		if (bkm != NULL) {delete bkm; bkm = NULL;}
		if (dyn != NULL) {delete dyn; dyn = NULL;}
	};
	void setup_emp() {
#if _EWALD_SUM == _MultiPole_EWALD_SUM
		if (!this->eNeutral && emp == NULL) emp = new ELECTRO_MULTIPOLE;
#endif
	};
	void reset() {
		_BASIC_CLUSTER<MATOM>::reset(); 
		reset_kin_dyn();
#if _EWALD_SUM == _MultiPole_EWALD_SUM
		if (emp != NULL) {delete emp; emp = NULL;}
#endif
		// implicit-solvation variables
		psolv = NULL;
		impsolv_rship.release_relationship();
		impsolv_fpar.release();
	};
	~BASIC_CLUSTER() {
		reset();
#if _EWALD_SUM == _MultiPole_EWALD_SUM
		ERship.release_relationship();
#endif
		LJRship.release_relationship();
	};

	bool linkPars(ATOM_COMM_PARS_DATABASE &atompar_db, char *errmsg) {
		int i = 0, j = 0;
		ATOM_COMM_PARS *atompar = NULL;
		for (i = 0; i < nAtoms; i++) {
			if (atom[i].par != NULL) continue;
			atom[i].par = atompar_db.search_par_indx((char)atom[i].aindx);
			if (atom[i].par == NULL) {
				sprintf(errmsg, "atom parameter of #%d is not defined in program", int(atom[i].aindx));
				return false;
			}
		}
#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 1
		this->eNeutral = true;
#else
		check_eneutral();
	#if _EWALD_SUM == _MultiPole_EWALD_SUM
		this->setup_emp();
	#endif
#endif
		return true;
	};

	void calElectroMultipole(bool bEwaldSum);
	void check_eneutral() {
		int na; eNeutral = true;
		for (na = 0; na < nAtoms; na++) {
			if (!atom[na].eNeutral) {eNeutral = false; return;}
		}
	};
	float cal_geometry_radius();
	void calSolvationRadius();
	void show_BrownianForcePars(CHAIN<short> **ch_id);

	void setup_cluster_matom();
	void setup_cluster_fatom(); // be sure the eneutral is checked
};

// CLUSTER ROTATE ALONG ITS PARENT -- using for parallel rotation
class CLUSTER_ROT {
public:
	double dta;
	R R0, Rg; // R0 -- from local rotation along parent hinge, R is total rotation matrix in whole macro-molecule
	VECTOR3 dr, drg; // dr is translation vector in local rotation matrix; drg is that in global rotation matrix

	CLUSTER_ROT() {dta = 0;};
};

class MOL_BASE_MOVEMENT {
public:
	bool bRot; // true -- rotate along drot; false -- rotate axises to given zaxis and xaxis
	VECTOR3 dr, drot, zaxis, xaxis;
};

class CLUSTER_MASS_CENTER_DYN { // using for NPT
public:
	MATRIX<3> I;
	VECTOR3 w, v;
};

class HINGE_CONSTRAINT {
public:
	MATOM *patom_10, *patom_11, *patom_12;
	MATOM *patom_00, *patom_01, *patom_02;
	//      11     12
	//       \    /
	//        \  /
	//         10
	//         |
	//         |  <== HINGE
	//         |
	//         00
	//        /  \
	//       /    \
	//      01     02

	VECTOR3 fc1[3], fc0[3]; // the constraint force on the 6 atoms

	HINGE_CONSTRAINT() {patom_10 = NULL; patom_11 = NULL; patom_12 = NULL; patom_00 = NULL; patom_01 = NULL; patom_02 = NULL;};
};

template <class TATOM, class TORSION_PAR> class DIHEDRAL {
public:
	TATOM* dh[4];
	TORSION_PAR *tpar;
	DIHEDRAL() {dh[0] = NULL; dh[1] = NULL; dh[2] = NULL; dh[3] = NULL; tpar = NULL;};
	~DIHEDRAL() {dh[0] = NULL; dh[1] = NULL; dh[2] = NULL; dh[3] = NULL; tpar = NULL;};
};


// kinematic parameters for cluster
class KINEMATIC {
public:
	// phai & H are Jacobian matrixs of this cluster to its parent
	// In this program, the cordinate frames in all clusters have same direction
	// The only difference is their original points are different (parent Om[])
	// So, H is a constant matrix, H[i][i] = 1. So, we do not have specific variable for H
	MATRIX<6> phai, phai_T; // phai_T = phai*
	VECTOR<6> H;
	double ta, tp, tpp; // torsion angle in arc-degree, its speed and accelerating speed
	KINEMATIC() {ta = 0; tp = 0; tpp = 0;};
};

// matrixes for Recursive solvation of Inverse Newton-Euler Equations
class INVERSE_NE {
public:
	MATRIX<3> I; // inertia matrix
	MATRIX<6> M; //spatial inertia
	//VECTOR<6> f; // torque and force of in this cluster's inertial cordinate
	VECTOR<6> a, b; // Coriolis acceleration term and gyroscopic spatial force
	MATRIX<6> P, Pp;
	VECTOR<6> G;
	VECTOR<6> zp;
	double D;
	double T; // explicit torque along the hinge
	double gamma;

	// using for calculation of free base acceleration 
	MATRIX<6> tao; // (I - H* x G*) * phai_T(parent. this)
	VECTOR<6> beta;
	// more variables
	MATRIX<6>  m0; // m0 = phai(base=>cluster) * M
	//MATRIX<6> mA; // mA = I - H* x G*

	// using for force calculation induced from the Hoover correction on tp
	VECTOR<6> dalpha_H;

	// using for force translation to mass center
	MATRIX<6> phai_cm2me, phai_T_cm2me; // phai from mass center of macromolecule to this cluster
	MATRIX<6> phai_base2me, phai_T_base2me; // phai from mass center of macromolecule to this cluster

	INVERSE_NE() {T = 0; D = 0; gamma = 0;};
};

// INVERSE_NE for base cluster
// because base cluster is free, H  = I (unit vector to any direction) while not Hinge !
// So, D = P; G = P * H * D^-1 = I; tao = I - G * H = 0; Pp = 0;
// gamma is Matrix, gamma = P^-1 * (P * a + b +/- fc)
class BASE_INVERSE_NE {
public:
	MATRIX<6> invD;
};

// atom cluster defined by A. JAIN
class CLUSTER : public BASIC_CLUSTER {
public:
	float Iw; // inertia moment along its parent hinge
	VECTOR3 cm0; //center of mass relative to hinge point (Op), of inertial cordinate
	VECTOR3 cm; //mass center in global cordinate

	// inertial cordinates of this cluster transfer to parent's inertial cordinates
	// it is actually the inertial vector of Om [to this cluster] in parent's inertial cordinate
	VECTOR3 l;
	VECTOR3 up; // the hinge / axis connecting to parent, up is the unit vector of vpHinge

	KINEMATIC km;
	INVERSE_NE *invNE;

	CLUSTER_ROT rot;

	float Ek;
	bool fast_speed;

	CLUSTER_MASS_CENTER_DYN *mcDyn;
	VECTOR<6> *fphinge; // force on parent hinge
	HINGE_CONSTRAINT *hinge_constraint;

	ARRAY< DIHEDRAL<MATOM, TORSION_PARS> > dh; // 4-atom dihedral relative to the parent hinge

	VECTOR3 rg_cmm; // position in CMM
	VECTOR3 *vp0; // position of cluster in MMOLECULE, vp = vp0 + dr_cmm;

	CLUSTER(int nAtoms = 0) {
		int i = 0;
		this->nAtoms = nAtoms;
		atom = NULL;
		if (nAtoms > 0) atom = new MATOM[nAtoms];

		M = 0; Iw = 0; 
		vp0 = &cm; vcm = &cm; vp = &rg_cmm; Ek = 0; fast_speed = false;
		invNE = NULL; mcDyn = NULL;
		fphinge = NULL; hinge_constraint = NULL;
	};

	void update_vp() {
		if (bCMMShift) {V3plusV3((*vp0), dr_cmm, (*vp))}
		else {V32V3((*vp0), (*vp))}
	};
	
	~CLUSTER() {
		reset();
		if (invNE != NULL) {delete invNE; invNE = NULL;}
		if (mcDyn != NULL) {delete mcDyn; mcDyn = NULL;}
		if (fphinge != NULL) {delete fphinge; fphinge = NULL;}
		if (hinge_constraint != NULL) {delete hinge_constraint; hinge_constraint = NULL;}
		dh.release();
	};
	
	void shift(VECTOR3 &d);
	void rotate(R &r); // given rotation matrix
	void shift_rotate(VECTOR3 &dv, VECTOR3 &v0, R &r);
	void rotate(VECTOR3 &axis, double omega, bool arc); // omega in degree ?
	void rotate(VECTOR3 &r1, VECTOR3 &r2, double omega, bool arc); // r1 is base point, omega in degree ?
	void align(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &r0, VECTOR3 &u); // align r1==>r0 & r12 || u

	bool xyz_save(ofstream &out, int indx = -1);

	void InertialCord();
	void MassCenter();
	void calJacobianMatrix();
	void InertiaTensor();
	void ab();

	void setup_hinge_dihedral();

	void setup_solvation_neighbors(int N);
	void setup_kin_dyn() {
		BASIC_CLUSTER::setup_kin_dyn();
		if (invNE == NULL) invNE = new INVERSE_NE;
	};
	void setup_mc_dyn() {
		if (mcDyn == NULL) mcDyn = new CLUSTER_MASS_CENTER_DYN;
	};
	void cal_cluster_inertia_tensor_masscenter();
};

bool setup_hinge_constraint(CLUSTER *pc);

void setup_solvation_neighbors(CLUSTER* rpc, float dis_cut);

// c1 & c2 could not has parent-child relationship, but linked, torsion angle is viewing from c2 to c1
double torsion_angle(BASIC_CLUSTER *c1, BASIC_CLUSTER *c2); // torsion angle is viewing from c2 to c1
double IdeaTorsionAngleFromSymetry(MATOM *patom);
double IdeaTorsionAngleFromSymetry(BASIC_CLUSTER *c1, BASIC_CLUSTER *c2);
// c1 & c2 could not has parent-child relationship, but linked, torsion angle is viewing from c2 to c1
//void calTorsionProfile(BASIC_CLUSTER *c1, BASIC_CLUSTER *c2);

#define PARENT(pc) (pc->n2Parent >= 0 && pc->hinge[pc->n2Parent].plink != NULL)

#define NO_CHILD(pc) (pc->nHinges == 1 && pc->n2Parent == 0)

#define PARENT_CLUSTER(pc) ((CLUSTER*)(pc->hinge[n2Parent].plink))

#define CHILD_CLUSTER(child, pc, i, j, n) n = 0; \
	for (j = 0; j < pc->nHinges; j++) {if (n == i) break; if (j != pc->n2Parent) n++;} \
	if (j < pc->nHinges) child = (CLUSTER*)(pc->hinge[j].plink); else child = NULL;

// origin of the cordinate locating at the cluster pc
#define CLUSTER_ORIGIN(pc, i) (pc->Op->v[i])

bool cp(BASIC_CLUSTER *d, BASIC_CLUSTER *s);
bool cp(CLUSTER *d, CLUSTER *s);
void shift_cluster_and_children(CLUSTER *c, VECTOR3 &d);
void tweak_cluster_and_children(CLUSTER *c, VECTOR3 &v0, R &r);
void tweak_cluster_and_children(CLUSTER *c, VECTOR3 &v0, VECTOR3 &axis, double omega, bool arc);
void shift_tweak_cluster_and_children(CLUSTER *c, VECTOR3 &dv, VECTOR3 &v0, R &r);
void shift_tweak_cluster_and_children(CLUSTER *c, VECTOR3 &dv, VECTOR3 &v0, VECTOR3 &axis, double omega, bool arc);

bool ConnectClusters(CLUSTER *c1, int iHinge1, CLUSTER *c2, int iHinge2, char bound_type = UNKNOWN_BOUND); // connect c2 [child] to c1 [parent] by aligning vpHinge in c2 parallel to vcHinge in c1

void setup_parent_child_relationship(CLUSTER *base, bool linked = false);
void setup_parent_child_relationship(CLUSTER *parent, CLUSTER *child);

bool check_parent_child_setup(CLUSTER *cluster, int nc);

// from the aspect of mass center
class MOL_DYN_KIN {
public:
	VECTOR<3> r0; // ceneter of mass in experimental cordinate
	VECTOR<3> cm; // center of mass in inertial cordinate
	MATRIX<6> I, invI;  // total inertia matrix relative to center of mass
	VECTOR<6> V, alpha; // 6-dim. of velocity and acceleration
	VECTOR<6> P; // spatial moment of molecule at center of mass
	VECTOR<6> fc; // real force and torque on mass center
	VECTOR<6> b; // Gyroscopic Spatial Force [force and torque] on mass center
	VECTOR<6> f_Hoover; // force and torque due to Hoover correction on motion
	VECTOR<6> f_Brown;
	VECTOR<6> f; // f_Hoover + fc + f_Hoover
	float M; // total mass


	// variables for free base acceleratioin calculation
	MATRIX<6> M0, invM0; //M0 += phai(base=>cluster) * M * (I - H* x G*) * phai_T(base=>cluster) for all clusters
};

template <class T> class _GROUP {
public:
	int uid; // the type of this GROUP
	int nUnit;
	T** u;
	_GROUP() {nUnit = 0; u = NULL; uid = 0;};
	void reset() {
		if (u != NULL) {
			for (int i = 0; i < nUnit; i++) u[i] = NULL;
			delete[] u; u = NULL;
		}
		nUnit = 0;
	};
	void set(int n) {
		reset(); if (n <= 0) return;
		nUnit = n;
		u = new T*[nUnit];
		for (n = 0; n < nUnit; n++) u[n] = NULL;
	};
	~_GROUP() {reset();};
};

// macromolecule is a combination of chained clusters
// total spatial momentum includes two part:
// 1.  considering the whole macromolecule as a unit rigid body (or body-frame, experimental frame), there is a momentum.
//      we can use the position of base cluster as reference, or the mass center
// 2.  each cluster inside the macromolecule are movable. The rotations inside the macromolecule along their own axis, or hinges.
//     IMPORTANT: even using the position of base cluster as reference frame, base cluster has its own Non-Zero velocity (translation and rotation).
//     In another word, the velocity of mass center is NOT a siimple transportation of base cluster's velocity to the mass center, but also the 
//     the velocity of base cluster relative the reference frame
// 3.  from the view of each cluster, the cluster itself with the child or grand-child clusters can be considered as a combination
//      of rigid unit, and the rotations inside

// The Newton Law at mass center: d P(t) / dt = F, the real force and torque relative to the mass center, including the Hoover force, but NO gyroscopic force

class MMOLECULE {
public:
	short mIndx; // the index of the molecule in the whole system, each molecule has an unique number, starting from 0
	short mType; // the type of molecule, 0 -- MM, 1 -- SM, 2 -- PM
	short mTypeIndx; // uid of the molecule in the CELL, because molecule is not uniquely defined
	char mol[250];
	int nCluster;
	CLUSTER *cluster;
	CLUSTER *base;
	CLUSTER **tip;
	int ntips;

	TREE<CLUSTER> ctree;
	ARRAY< ARRAY<int> > indx_tree;

	VECTOR3 zaxis, xaxis; // x and z axises directions rotated from originaly seted-up

	VECTOR3 r0; // the origin of base cluster's cordinate relative to the molecule's center mass

	VECTOR<6> V0, alpha0; // velocity of base cluster or external speed
	MATRIX<6> M0, invM0; // inertia matrix of the whole cluster relative to the base cluster, summation over all cluster for phai(base=>cluster) * M * phai_T(base=>cluster)
	VECTOR<6> Iw, f, b; // total spatial moment, force and gyroscopic spatial force, relative to the base cluster

	MOL_DYN_KIN mDynKin; // of mass center

	BASE_INVERSE_NE binvNE; // for free base cluster specifically

	MOL_BASE_MOVEMENT mbm;

	bool free_base; // to indicate the base cluster is free
	float Ek, Ekt, Ektx, Ekty, Ektz, Ep; // kinetic and potential energy
	float Ek_internal, Ek_frame, Ek_frame_trans, Ek_frame_rot;

	_GROUP<CLUSTER> *monomer;
	int nMonomer;

	MMOLECULE() {
		mIndx = 0; mType = -1; mTypeIndx = -1;
		cluster = NULL; nCluster = 0; strcpy(mol, "\0"); 
		base = NULL; ntips = 0;
		tip = NULL;
		free_base = true;
		zaxis.v[0] = 0; zaxis.v[1] = 0; zaxis.v[2] = 1; // z-axis : [0, 0, 1]
		xaxis.v[0] = 1; xaxis.v[1] = 0; xaxis.v[2] = 0; // x-axis : [1, 0, 0]
		Ek = 0; Ep = 0;
		nMonomer = 0; monomer = NULL;
	};
	void InitMolFrameAxis() {
		zaxis.v[0] = 0; zaxis.v[1] = 0; zaxis.v[2] = 1; // z-axis : [0, 0, 1]
		xaxis.v[0] = 1; xaxis.v[1] = 0; xaxis.v[2] = 0; // x-axis : [1, 0, 0]
	};
	void reset_tips() {
		int i = 0;
		if (tip != NULL) {
			for (i = 0; i < ntips; i++) tip[i] = NULL;
			delete[] tip; tip = NULL;
		}
		ntips = 0;
	};
	void reset() {
		base = NULL;
		reset_tips();
		ctree.release_tree(false);
		if (cluster != NULL) {delete[] cluster; cluster = NULL;}; nCluster = 0;
		if (monomer != NULL) {delete[] monomer; monomer = NULL;} nMonomer = 0;
		InitMolFrameAxis();
	};
	~MMOLECULE() {reset();};
	bool xyz_save(ofstream &out);
	void construct_tree();
	void checkMolecule(bool tree);
	bool linkPars(ATOM_COMM_PARS_DATABASE &atompar_db, char *errmsg) {
		int i = 0;
		for (i = 0; i < nCluster; i++) if (!cluster[i].linkPars(atompar_db, errmsg)) return false;
		return true;
	};
	void shift(CLUSTER *c, VECTOR3 &d); // shift the cluster and its children with given shift
	void shiftMol(VECTOR3 &d); // shift the whole molecule
	void tw(CLUSTER *c, VECTOR3 &v0, R &r); // tweak cluster and its children with given shift and rotation matrix
	void twTA(int n, double domega, bool arc); // twist nth cluster along its vpHinge with domega
	void twMol(VECTOR3 &vr0, VECTOR3 &axis, double domega, bool arc); // twist molecule with axis and r0 with domega
	void twMol(VECTOR3 &vr0, VECTOR3 &drot, bool arc); // twist molecule with dw

	void shift_tw_mol(VECTOR3 &dv, VECTOR3 &v0, VECTOR3 &axis, double domega, bool arc);
	void shift_tw_mol(VECTOR3 &dv, VECTOR3 &v0, VECTOR3 &drot, bool arc); // twist molecule with dw

	void resetDynPars();

	void CalMolInertiaMoment();
	void CalGyroscopicSpatialForce();
	void SetBaseForce();
	void CalFreeBaseVelocity(VECTOR<6> &Iw0_base);
	void CalFreeBaseVelocity_GivenMassCenterMomentum();
	void CalMassCenterVelocity_MassCenterMomentum(); // was named SetBaseKinetic
	void calMassCenter(bool cal_cluster_mc);
	void SetMassCenter(double *v, bool cal_mc, bool cal_cluster_mc);
	void SetBasePosZero();

	void setMolIndx(int indx);
	//void SetupHingeTorsion();
	//void calTorsion(double &E);
	void calTorsionAngle(); //calculate torsion angle with parent-child relationship setup
	void highlight(bool status) {
		for (int n = 0; n < nCluster; n++) cluster[n].highlight = status;
	};

	void calTotalForce();

	void set_base_move2(VECTOR3& r2, VECTOR3& zaxis2, VECTOR3& xaxis2) {
		mbm.bRot = false;
		VECT3((*(base->Op)), r2, mbm.dr)
		V32V3(zaxis2, mbm.zaxis) V32V3(xaxis2, mbm.xaxis)
	};
	void Init_MolFrameRot2(VECTOR3& dr, VECTOR3& zaxis2, VECTOR3& xaxis2);
	void set_base_move(VECTOR3& dr, VECTOR3& drot) {
		mbm.bRot = true;
		V32V3(dr, mbm.dr) V32V3(drot, mbm.drot)
	};
	void Init_MolFrameRot(VECTOR3& dr, VECTOR3& drot);
	// function using for multi-thread cluster-rotations
	void setup_base_movement(VECTOR3 &dv, VECTOR3 &v0, VECTOR3 &drot, bool arc);

	void calClusterSolvationRadius() {
		for (int nc = 0; nc < nCluster; nc++) cluster[nc].calSolvationRadius();
	};

	void setup_Solvation_neighbors(int N);
	void setup_Solvation_neighbors(float dis_cut);
	void setup_kin_dyn() {
		for (int nc = 0; nc < nCluster; nc++) cluster[nc].setup_kin_dyn();
	};
	void setup_emp() {
		for (int nc = 0; nc < nCluster; nc++) cluster[nc].setup_emp();
	};
	void setup_mcDyn() {
		for (int nc = 0; nc < nCluster; nc++) cluster[nc].setup_mc_dyn();
	};
	void show_cluster_BrownianForcePars(CHAIN<short> **ch_id);
	void init_hinge_force() {
		for (int ic = 0; ic < nCluster; ic++) {
			if (cluster[ic].fphinge == NULL) cluster[ic].fphinge = new VECTOR<6>;
		}
	};
	bool setup_hinge_constraint();
	void setup_hinge_dihedral() {
		for (int ic = 0; ic < nCluster; ic++) cluster[ic].setup_hinge_dihedral();
	};
	void setup_cluster_matom();
	void setup_cluster_fatom() { // be sure the eneutral is checked
		for (int ic = 0; ic < nCluster; ic++) cluster[ic].setup_cluster_fatom();
	};
};

//void calTorsion(MMOLECULE*mm, int nc1, int nc2, double& E);

bool cp(MMOLECULE *m2, MMOLECULE* m1, bool steup_tree = true);

// function using for multi-thread cluster-rotations
void setup_LocalRotMatrix(MMOLECULE *mm, int nc1, int nc2);
void setup_all_LocalRotMatrix(MMOLECULE &mm);
void setup_GlobalRotMatrix(MMOLECULE &mm);
void tweak_clusters(MMOLECULE *mm, int nc1, int nc2);
void tweak_all_clusters(MMOLECULE &mm);
void MacroMolMove(MMOLECULE *mm);

class RELATIONSHIP {
public:
	int n1, iHinge1, n2, iHinge2;
	char bound_type;
	RELATIONSHIP() {n1 = 0; n2 = 0; iHinge1 = 0; iHinge2 = 0; bound_type = SINGLE_BOUND;};
	void set_relationship(int n1, int iHinge1, int n2, int iHinge2, char bound_type = SINGLE_BOUND) {
		this->n1 = n1; this->iHinge1 = iHinge1;
		this->n2 = n2; this->iHinge2 = iHinge2;
		this->bound_type = bound_type;
	};
};

class MONOMER { // monomer having one parent and one child
public:
	int uid; // the type of this monomer
	int nCluster; // number of cluster in each monomer
	char name[20];
	CLUSTER *cluster;
	int nHinges;
	int *nc, *iHinge; // the cluster connect to parent monomer

	int nrship;
	RELATIONSHIP *rship;

	int nbase; // if this monomer has no parent, which cluster inside will be the parent

	MONOMER() {
		nCluster = 0; cluster = NULL; strcpy(name, "\0");
		rship = NULL; nrship = 0;
		nc = NULL; iHinge = NULL;
		nbase = 0; uid = 0;
	};
	void setHinges(int n) {
		if (nc != NULL) {delete[] nc; nc = NULL;}
		if (iHinge != NULL) {delete[] iHinge; iHinge = NULL;}
		nHinges = n;
		if (n <= 0) return;
		nc = new int[n]; iHinge = new int[n];
		for (int i = 0; i < n; i++) {nc[i] = 0; iHinge[i] = 0;}
	};
	void reset(int n, int nrship) {
		if (cluster != NULL) {delete[] cluster; cluster = NULL;}
		if (rship != NULL) {delete[] rship; rship = NULL;}
		nCluster = n;
		if (n > 0) cluster = new CLUSTER[n];
		this->nrship = nrship;
		if (nrship > 0) rship = new RELATIONSHIP[nrship];
	};
	~MONOMER() {reset(0, 0); setHinges(0);};
};

// rotate all the clusters related with pc, except the one linked with parent_hinge, with given rotation matrix, original point -- r0, and after the rotation, shift dr
void series_shift_rotate(CLUSTER *pc, int parent_hinge, VECTOR3 &dr, VECTOR3 &r0, R &rmatrix);
// rotate pc + all clusters associated with pc except the one linked with parent_hinge, along parent_hinge with domega
void series_rotate(CLUSTER *pc, int parent_hinge, double domega, bool arc);

void calClusterAlpha(MMOLECULE &mm);
void calClusterVelocity0(MMOLECULE &mm, bool zeroBaseVelocity = false);

double calKineticEnergy(MMOLECULE &mm); // return energy in kT of whole molecule
double calFrameKineticEnergy(MMOLECULE &mm); // return energy in kT of the center mass translation and rotation

// return energy in kT
// calculate kinetic energy of the cluster with its children + grand-children
double calClusterKineticEnergy(CLUSTER *pc);

// rotation of cluster is initilized based on the velocity of base cluster
double InitMolClusterVelocity(MMOLECULE &mm);

// Spatial Moment is relative to the base cluster's cordinate
void calSpatialMomentum(MMOLECULE &mm);

// Spatial Acceleration Moment is relative to the base cluster's cordinate
void calSpatialAccel(MMOLECULE &mm, VECTOR<6> &Ialpha);

// Spatial Moment is relative to molecular Mass Center
void calSpatialMomentum0(MMOLECULE &mm, bool bCalBaseMomentum = true);

// Spatial Acceleration Moment is relative to molecular Mass Center
void calSpatialAccel0(MMOLECULE &mm, VECTOR<6> &Ialpha);

double ScaleKineticEnergy(MMOLECULE &mm, double Ek0, bool reset_total_SpacialMoment = true);

void ClusterMaxVelocity(MMOLECULE &mm, double &vmax, double &wmax);

bool atom_indx(MMOLECULE &mm, MATOM* patom, int &nc, int &na, int &indx);


// return energy in kT
void calClustersKineticEnergy(MMOLECULE *mm, int n1, int n2, double &Ek);
// return energy in kT
double ParallelCalKineticEnergy(MMOLECULE &mm);
// Spatial Moment is relative to the base cluster's cordinate
void ParallelCalSpatialMomentum(MMOLECULE &mm);

// Spatial Acceleration Moment is relative to the base cluster's cordinate
void ParallelCalSpatialAccel(MMOLECULE &mm, VECTOR<6> &Ialpha);

void CalClusterElectroMultipole(MMOLECULE *mm, int nc1, int nc2, bool bEwaldSum);



//**********************************************************************************
//**  MULTI-THREAD FUNCTIONS FOR MACROMOL OPERATION WITHOUT ADDITIONAL PARAMETERS **
//**********************************************************************************

template <class MACROMOL> class _MACROMOL_THREAD_VARS : public _THREAD_VARS {
public:
	MACROMOL* mm;

	_MACROMOL_THREAD_VARS() {
		this->mm = NULL;
	};
	void set_pars(int iThread, MACROMOL* m) {
		this->mm = m;
		_THREAD_VARS::set_pars(iThread);
	};
	~_MACROMOL_THREAD_VARS() {};
};

template <class MACROMOL> class MACROMOL_THREAD_VARS : public _MACROMOL_THREAD_VARS<MACROMOL> {
public:
	typedef void (*MACROMOL_OPERATOR)(MACROMOL*, int, int);

	int n1, n2;// iThread;

	void* func;

	MACROMOL_THREAD_VARS() {
		this->mm = NULL; func = NULL; n1 = 0; n2 = 0;
	};
	void set_pars(int iThread, MACROMOL* m, int n1, int n2) {
		this->n1 = n1; this->n2 = n2;
		_MACROMOL_THREAD_VARS<MACROMOL>::set_pars(iThread, m);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~MACROMOL_THREAD_VARS() {};
};

template <class MACROMOL, class JOB> void assignMMJob(MACROMOL *mm, JOB *job, int nJobs) {
	int nJob = nJobs;
	int i = 0, nstep = mm->nCluster / nJob;
	int n1 = 0, n2 = 0;
	if (nstep < 10 && nJob > 1) {
		nJob = mm->nCluster / 10; 
		if (nJob < 1) nJob = 1;
		else nstep = 10;
	}
	if (nJob * nstep < mm->nCluster) nstep += 1;
	for (i = 0; i < nJob; i++) {
		if (n1 >= mm->nCluster) {job[i].n1 = mm->nCluster; job[i].n2 = mm->nCluster; continue;}

		if (i == nJob - 1) n2 = mm->nCluster - 1;
		else n2 = n1 + nstep - 1;

		if (n2 >= mm->nCluster) n2 = mm->nCluster - 1;

		//job[i].set_pars(i, mm, n1, n2);
		job[i].mm = mm;
		job[i].n1 = n1;
		job[i].n2 = n2;
		n1 += nstep;
	}
};

#if _SYS_ == _WINDOWS_SYS_
template <class MACROMOL> int MacroMol_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MACROMOL> void* MacroMol_Thread(void *void_par) { // thread-function
#endif
	MACROMOL_THREAD_VARS<MACROMOL> *par = (MACROMOL_THREAD_VARS<MACROMOL>*)void_par;
	//((MACROMOL_THREAD_VARS<MACROMOL>::MACROMOL_OPERATOR)(par->func))(par->mm, par->n1, par->n2);
	typedef void (*MACROMOL_OPERATOR)(MACROMOL*, int, int);
	((MACROMOL_OPERATOR)(par->func))(par->mm, par->n1, par->n2);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class MACROMOL> void MacroMoleculeOperate(void* func, MACROMOL &mm, int nc1, int nc2, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*MACROMOL_OPERATOR)(MACROMOL*, int, int);
		((MACROMOL_OPERATOR)(func))(&mm, nc1, nc2);
		return;
	}

//void MacroMoleculeOperate(MMOLECULE &mm, int nc1, int nc2, MM_PVOID func, int nThreads) {
	MACROMOL_THREAD_VARS<MACROMOL> *thread_vars = new MACROMOL_THREAD_VARS<MACROMOL>[nThreads];
	int iThread = 0;
	if (nc2 >= mm.nCluster) nc2 = mm.nCluster - 1;
	int n1 = nc1, n2 = 0, nstep = (nc2 - nc1 + 1) / nThreads;
	if (nstep < 10 && nThreads > 1) {
		nThreads = (nc2 - nc1 + 1) / 10; 
		if (nThreads < 1) nThreads = 1;
		else nstep = 10;
	}
	if (nThreads * nstep < nc2 - nc1 + 1) nstep += 1;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = nc2;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &mm, n1, n2);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MacroMol_Thread<MACROMOL>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MacroMol_Thread<MACROMOL>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].mm = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
};


//**********************************************************************************
//********  MULTI-THREAD FUNCTIONS FOR MACROMOL OPERATION WITH A PARAMETERS ********
//********  WITHOUT RETURNING RESULT                                        ********
//**********************************************************************************

template <class MACROMOL, class TVAR> class MACROMOL_THREAD_VARS_1 : public _MACROMOL_THREAD_VARS<MACROMOL> {
public:
	typedef void (*MACROMOL_OPERATOR1)(MACROMOL*, int, int, TVAR&);

	int n1, n2;
	TVAR v;

	void* func;

	MACROMOL_THREAD_VARS_1() {
		this->mm = NULL; func = NULL; n1 = 0; n2 = 0;
	};
	void set_pars(int iThread, MACROMOL* m, int n1, int n2, TVAR& v) {
		this->n1 = n1; this->n2 = n2; this->v = v;
		_MACROMOL_THREAD_VARS<MACROMOL>::set_pars(iThread, m);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~MACROMOL_THREAD_VARS_1() {};
};


#if _SYS_ == _WINDOWS_SYS_
template <class MACROMOL, class TVAR> int MACROMOL_Thread_1(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MACROMOL, class TVAR> void* MACROMOL_Thread_1(void *void_par) { // thread-function
#endif
	MACROMOL_THREAD_VARS_1<MACROMOL, TVAR> *par = (MACROMOL_THREAD_VARS_1<MACROMOL, TVAR>*)void_par;
	//((MACROMOL_THREAD_VARS_1<MACROMOL, TVAR>::MACROMOL_OPERATOR1)(par->func))(par->mm, par->n1, par->n2, par->v);
	typedef void (*MACROMOL_OPERATOR)(MACROMOL*, int, int, TVAR&);
	((MACROMOL_OPERATOR)(par->func))(par->mm, par->n1, par->n2, par->v);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class MACROMOL, class TVAR> void MacroMoleculeOperate1(void* func, MACROMOL &mm, int nc1, int nc2, TVAR& v, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*MACROMOL_OPERATOR)(MACROMOL*, int, int, TVAR&);
		((MACROMOL_OPERATOR)(func))(&mm, nc1, nc2, v);
		return;
	}

//void MacroMoleculeOperate(MMOLECULE &mm, int nc1, int nc2, MM_VECTOR6 func, int nThreads, VECTOR<6> &res) {
	MACROMOL_THREAD_VARS_1<MACROMOL, TVAR> *thread_vars = new MACROMOL_THREAD_VARS_1<MACROMOL, TVAR>[nThreads];
	int iThread = 0;
	if (nc2 >= mm.nCluster) nc2 = mm.nCluster - 1;
	int n1 = nc1, n2 = 0, nstep = (nc2 - nc1 + 1) / nThreads;
	if (nstep < 10 && nThreads > 1) {
		nThreads = (nc2 - nc1 + 1) / 10; 
		if (nThreads < 1) nThreads = 1;
		else nstep = 10;
	}
	if (nThreads * nstep < nc2 - nc1 + 1) nstep += 1;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = nc2;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &mm, n1, n2, v);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MACROMOL_Thread_1<MACROMOL, TVAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MACROMOL_Thread_1<MACROMOL, TVAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	delete[] thread_vars; thread_vars = NULL;
	return;
};



//**********************************************************************************
//********  MULTI-THREAD FUNCTIONS FOR MACROMOL OPERATION WITH A PARAMETER  ********
//********  WITHOUT RETURNING RESULT; PARAMETER FOR EACH THREAD WERE GIVEN  ********
//**********************************************************************************

template <class MACROMOL, class TVAR> class MACROMOL_THREAD_VARS_1_MV : public _MACROMOL_THREAD_VARS<MACROMOL> {
public:
	typedef void (*MACROMOL_OPERATOR1)(MACROMOL*, int, int, TVAR&);

	int n1, n2;
	TVAR *v;

	void* func;

	MACROMOL_THREAD_VARS_1_MV() {
		this->mm = NULL; func = NULL; n1 = 0; n2 = 0; v = NULL;
	};
	void set_pars(int iThread, MACROMOL* m, int n1, int n2, TVAR* v) {
		this->n1 = n1; this->n2 = n2; this->v = v;
		_MACROMOL_THREAD_VARS<MACROMOL>::set_pars(iThread, m);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~MACROMOL_THREAD_VARS_1_MV() {};
};


#if _SYS_ == _WINDOWS_SYS_
template <class MACROMOL, class TVAR> int MACROMOL_Thread_1_MV(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MACROMOL, class TVAR> void* MACROMOL_Thread_1_MV(void *void_par) { // thread-function
#endif
	MACROMOL_THREAD_VARS_1_MV<MACROMOL, TVAR> *par = (MACROMOL_THREAD_VARS_1_MV<MACROMOL, TVAR>*)void_par;
	//((MACROMOL_THREAD_VARS_1<MACROMOL, TVAR>::MACROMOL_OPERATOR1)(par->func))(par->mm, par->n1, par->n2, *(par->v));
	typedef void (*MACROMOL_OPERATOR)(MACROMOL*, int, int, TVAR&);
	((MACROMOL_OPERATOR)(par->func))(par->mm, par->n1, par->n2, *(par->v));
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class MACROMOL, class TVAR> void MacroMoleculeOperate1_mv(void* func, MACROMOL &mm, int nc1, int nc2, TVAR* v, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*MACROMOL_OPERATOR)(MACROMOL*, int, int, TVAR&);
		((MACROMOL_OPERATOR)(func))(&mm, nc1, nc2, v[0]);
		return;
	}

//void MacroMoleculeOperate(MMOLECULE &mm, int nc1, int nc2, MM_VECTOR6 func, int nThreads, VECTOR<6> &res) {
	MACROMOL_THREAD_VARS_1_MV<MACROMOL, TVAR> *thread_vars = new MACROMOL_THREAD_VARS_1_MV<MACROMOL, TVAR>[nThreads];
	int iThread = 0;
	if (nc2 >= mm.nCluster) nc2 = mm.nCluster - 1;
	int n1 = nc1, n2 = 0, nstep = (nc2 - nc1 + 1) / nThreads;
	if (nstep < 10 && nThreads > 1) {
		nThreads = (nc2 - nc1 + 1) / 10; 
		if (nThreads < 1) nThreads = 1;
		else nstep = 10;
	}
	if (nThreads * nstep < nc2 - nc1 + 1) nstep += 1;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = nc2;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &mm, n1, n2, v + iThread);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MACROMOL_Thread_1_MV<MACROMOL, TVAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MACROMOL_Thread_1_MV<MACROMOL, TVAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	delete[] thread_vars; thread_vars = NULL;
	return;
};

//**********************************************************************************
//**  MULTI-THREAD FUNCTIONS FOR MACROMOL OPERATION WITHOUT ADDITIONAL PARAMETERS **
//**  RETURN A RESULT                                                             **
//**********************************************************************************

template <class MACROMOL, class TVAR> class MACROMOL_THREAD_VARS_2 : public _MACROMOL_THREAD_VARS<MACROMOL> {
public:
	typedef void (*MACROMOL_OPERATOR2)(MACROMOL*, int, int, TVAR&);

	int n1, n2;
	TVAR result;

	void* func;

	MACROMOL_THREAD_VARS_2() {
		this->mm = NULL; func = NULL; n1 = 0; n2 = 0;
	};
	void set_pars(int iThread, MACROMOL* m, int n1, int n2) {
		this->n1 = n1; this->n2 = n2;
		_MACROMOL_THREAD_VARS<MACROMOL>::set_pars(iThread, m);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~MACROMOL_THREAD_VARS_2() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class MACROMOL, class TVAR> int MACROMOL_Thread_2(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MACROMOL, class TVAR> void* MACROMOL_Thread_2(void *void_par) { // thread-function
#endif
	MACROMOL_THREAD_VARS_2<MACROMOL, TVAR> *par = (MACROMOL_THREAD_VARS_2<MACROMOL, TVAR>*)void_par;
	//((MACROMOL_THREAD_VARS_2<MACROMOL, TVAR>::MACROMOL_OPERATOR2)(par->func))(par->mm, par->n1, par->n2, par->result);
	typedef void (*MACROMOL_OPERATOR)(MACROMOL*, int, int, TVAR&);
	((MACROMOL_OPERATOR)(par->func))(par->mm, par->n1, par->n2, par->result);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

template <class MACROMOL, class TVAR> void MacroMoleculeOperate2(void *func, MACROMOL &mm, int nc1, int nc2, int nThreads, TVAR &res, bool init_res = false, void *res_operator = NULL) {
	if (nThreads <= 1) {
		typedef void (*MACROMOL_OPERATOR)(MACROMOL*, int, int, TVAR&);
		res = (TVAR)0;
		((MACROMOL_OPERATOR)(func))(&mm, nc1, nc2, res);
		return;
	}

	res = (TVAR)0;
	MACROMOL_THREAD_VARS_2<MACROMOL, TVAR> *thread_vars = new MACROMOL_THREAD_VARS_2<MACROMOL, TVAR>[nThreads];
	int iThread = 0;
	if (nc2 >= mm.nCluster) nc2 = mm.nCluster - 1;
	int n1 = nc1, n2 = 0, nstep = (nc2 - nc1 + 1) / nThreads;
	if (nstep < 10 && nThreads > 1) {
		nThreads = (nc2 - nc1 + 1) / 10; 
		if (nThreads < 1) nThreads = 1;
		else nstep = 10;
	}
	if (nThreads * nstep < nc2 - nc1 + 1) nstep += 1;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = nc2;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &mm, n1, n2);
		thread_vars[iThread].set_func(func);
		if (init_res) thread_vars[iThread].result = res;
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MACROMOL_Thread_2<MACROMOL, TVAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MACROMOL_Thread_2<MACROMOL, TVAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	typedef void (*RES_OPERATOR)(TVAR&, TVAR&, TVAR&);
	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].mm = NULL;
		if (res_operator == NULL) res += thread_vars[iThread].result;
		else {
			((RES_OPERATOR)res_operator)(res, thread_vars[iThread].result, res);
		}
	}
	delete[] thread_vars; thread_vars = NULL;
	//return res;
};


//*****************************************************
//**  MULTI-THREAD FUNCTIONS FOR MACROMOL OPERATION  **
//**  WITH A PARAMETERS, RETURN A RESULT             **
//*****************************************************

template <class MACROMOL, class VVAR, class RESVAR> class MACROMOL_THREAD_VARS_3 : public _MACROMOL_THREAD_VARS<MACROMOL> {
public:
	typedef void (*MACROMOL_OPERATOR3)(MACROMOL*, int, int, VVAR&, RESVAR&);

	int n1, n2;
	VVAR v;
	RESVAR result;

	void* func;

	MACROMOL_THREAD_VARS_3() {
		this->mm = NULL; func = NULL; n1 = 0; n2 = 0;
	};
	void set_pars(int iThread, MACROMOL* m, int n1, int n2, VVAR& v) {
		this->n1 = n1; this->n2 = n2; this->v = v;
		_MACROMOL_THREAD_VARS<MACROMOL>::set_pars(iThread, m);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~MACROMOL_THREAD_VARS_3() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class MACROMOL, class VVAR, class RESVAR> int MACROMOL_Thread_3(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MACROMOL, class VVAR, class RESVAR> void* MACROMOL_Thread_3(void *void_par) { // thread-function
#endif
	MACROMOL_THREAD_VARS_3<MACROMOL, VVAR, RESVAR> *par = (MACROMOL_THREAD_VARS_3<MACROMOL, VVAR, RESVAR>*)void_par;
	//((MACROMOL_THREAD_VARS_3<MACROMOL, VVAR, RESVAR>::MACROMOL_OPERATOR3)(par->func))(par->mm, par->n1, par->n2, par->v, par->result);
	typedef void (*MACROMOL_OPERATOR3)(MACROMOL*, int, int, VVAR&, RESVAR&);
	((MACROMOL_OPERATOR3)(par->func))(par->mm, par->n1, par->n2, par->v, par->result);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class MACROMOL, class VVAR, class RESVAR> void MacroMoleculeOperate3(void *func, MACROMOL &mm, int nc1, int nc2, VVAR& v, int nThreads, RESVAR& res, bool init_res = false, void *res_operator = NULL) {
	if (nThreads <= 1) {
		typedef void (*MACROMOL_OPERATOR3)(MACROMOL*, int, int, VVAR&, RESVAR&);
		res = (RESVAR)0;
		((MACROMOL_OPERATOR3)(func))(&mm, nc1, nc2, v, res);
		return;
	}

	res = (RESVAR)0;
	MACROMOL_THREAD_VARS_3<MACROMOL, VVAR, RESVAR> *thread_vars = new MACROMOL_THREAD_VARS_3<MACROMOL, VVAR, RESVAR>[nThreads];
	int iThread = 0;
	if (nc2 >= mm.nCluster) nc2 = mm.nCluster - 1;
	int n1 = nc1, n2 = 0, nstep = (nc2 - nc1 + 1) / nThreads;
	if (nstep < 10 && nThreads > 1) {
		nThreads = (nc2 - nc1 + 1) / 10; 
		if (nThreads < 1) nThreads = 1;
		else nstep = 10;
	}
	if (nThreads * nstep < nc2 - nc1 + 1) nstep += 1;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = nc2;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &mm, n1, n2, v);
		thread_vars[iThread].set_func(func);
		if (init_res) thread_vars[iThread].result = res;
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MACROMOL_Thread_3<MACROMOL, VVAR, RESVAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MACROMOL_Thread_3<MACROMOL, VVAR, RESVAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	typedef void (*RES_OPERATOR)(RESVAR&, RESVAR&, RESVAR&);
	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].mm = NULL;
		if (res_operator == NULL) res += thread_vars[iThread].result;
		else {
			((RES_OPERATOR)res_operator)(res, thread_vars[iThread].result, res);
		}
	}
	delete[] thread_vars; thread_vars = NULL;
	//return res;
};

class FUNC_VAR {
public:
	void *func;
	FUNC_VAR() {func = NULL;};
	~FUNC_VAR() {func = NULL;};
	void operator = (FUNC_VAR &fvar) {
		this->func = fvar.func;
	};
};


void polarize_MM(MMOLECULE *mm, int n);
void polarize_VMM(MMOLECULE **mm, int n);

//**********************************************************************************
//**  Generized Multi-Thread Function, with the job assigned explicively **
//**********************************************************************************

template <class JOB> class JOB_THREAD_VARS : public FUNC_VAR, public _THREAD_VARS {
public:
	JOB *job;

	JOB_THREAD_VARS() {
		job = NULL;
	};
	void set_pars(int iThread, JOB* job) {
		this->job = job;
		_THREAD_VARS::set_pars(iThread);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~JOB_THREAD_VARS() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class JOB_VAR> int General_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class JOB_VAR> void* General_Thread(void *void_par) { // thread-function
#endif
	JOB_VAR *par = (JOB_VAR*)void_par;
	typedef void (*JOB_OPERATOR)(JOB_VAR *);
	((JOB_OPERATOR)(par->func))(par);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class JOB_VAR> void MTOperate(void* func, JOB_VAR *job, int nJob) {
	if (nJob == 1) {
		typedef void (*JOB_OPERATOR)(JOB_VAR *);
		((JOB_OPERATOR)(func))(job);
		return;
	}

	int i;
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nJob; i++) {
		ResetEvent(job[i].hMMThread);
		job[i].func = func;
		job[i].thread = AfxBeginThread((AFX_THREADPROC)General_Thread<JOB_VAR>, (LPVOID)(job + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nJob; i++) {
		WaitForSingleObject(job[i].hMMThread, INFINITE);
		//CloseHandle(job[iThread].hMMThread); job[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nJob; i++) {
		job[i].func = func;
		pthread_create(&(job[i].thread), NULL, &General_Thread<JOB_VAR>, (void *)(job + i));
	}
	for (i = 0; i < nJob; i++) {
		if (job[i].thread != 0) pthread_join(job[i].thread, NULL); // wait the thread over
	}
#endif
};


template <class JOB, class VAR1, class VAR2> class _GENERAL_JOB_VAR {
public:
	JOB *job;
	VAR1* v1;
	VAR2* v2;
	_GENERAL_JOB_VAR() {job = NULL; v1 = NULL; v2 = NULL;};
	~_GENERAL_JOB_VAR() {job = NULL; v1 = NULL; v2 = NULL;};
};

#if _SYS_ == _WINDOWS_SYS_
template <class JOB, class VAR> int General_Thread_1(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class JOB, class VAR> void *General_Thread_1(void *void_par) { // thread-function
#endif
	_GENERAL_JOB_VAR<JOB, VAR, void> *par = (_GENERAL_JOB_VAR<JOB, VAR, void>*)void_par;
	typedef void (*JOB_OPERATOR)(JOB*, VAR*);
	((JOB_OPERATOR)(par->job->func))(par->job, par->v1);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->job->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class JOB, class VAR> void MTOperate1(void* func, JOB *job, int nJob, VAR* var, bool bSingleVar = false) {
	if (nJob == 1) {
		typedef void (*JOB_OPERATOR)(JOB*, VAR*);
		((JOB_OPERATOR)(func))(job, var);
		return;
	}

	int i;
	_GENERAL_JOB_VAR<JOB, VAR, void> *gjob = new _GENERAL_JOB_VAR<JOB, VAR, void>[nJob];
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nJob; i++) {
		ResetEvent(job[i].hMMThread);
		job[i].func = func;
		gjob[i].job = job + i;
		if (bSingleVar) gjob[i].v1 = var;
		else gjob[i].v1 = var + i;
		gjob[i].v2 = NULL;
		job[i].thread = AfxBeginThread((AFX_THREADPROC)General_Thread_1<JOB, VAR>, (LPVOID)(gjob + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nJob; i++) {
		WaitForSingleObject(job[i].hMMThread, INFINITE);
		//CloseHandle(job[iThread].hMMThread); job[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nJob; i++) {
		job[i].func = func;
		gjob[i].job = job + i;
		if (bSingleVar) gjob[i].v1 = var;
		else gjob[i].v1 = var + i;
		gjob[i].v2 = NULL;
		pthread_create(&(job[i].thread), NULL, &General_Thread_1<JOB, VAR>, (void *)(gjob + i));
	}
	for (i = 0; i < nJob; i++) {
		if (job[i].thread != 0) pthread_join(job[i].thread, NULL); // wait the thread over
	}
#endif
	delete[] gjob;
};




#if _SYS_ == _WINDOWS_SYS_
template <class JOB, class VAR1, class VAR2> int General_Thread_2(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class JOB, class VAR1, class VAR2> void *General_Thread_2(void *void_par) { // thread-function
#endif
	_GENERAL_JOB_VAR<JOB, VAR1, VAR2> *par = (_GENERAL_JOB_VAR<JOB, VAR1, VAR2>*)void_par;
	typedef void (*JOB_OPERATOR)(JOB*, VAR1*, VAR2*);
	((JOB_OPERATOR)(par->job->func))(par->job, par->v1, par->v2);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->job->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class JOB, class VAR1, class VAR2> void MTOperate2(void* func, JOB *job, int nJob, VAR1* v1, VAR2* v2, bool bSingleVar = false) {
	if (nJob == 1) {
		typedef void (*JOB_OPERATOR)(JOB*, VAR1*, VAR2*);
		((JOB_OPERATOR)(func))(job, v1, v2);
		return;
	}

	int i;
	_GENERAL_JOB_VAR<JOB, VAR1, VAR2> *gjob = new _GENERAL_JOB_VAR<JOB, VAR1, VAR2>[nJob];
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nJob; i++) {
		ResetEvent(job[i].hMMThread);
		job[i].func = func;
		gjob[i].job = job + i;
		if (bSingleVar) {
			gjob[i].v1 = v1;
			gjob[i].v2 = v2;
		}
		else {
			gjob[i].v1 = v1 + i;
			gjob[i].v2 = v2 + i;
		}
		job[i].thread = AfxBeginThread((AFX_THREADPROC)General_Thread_2<JOB, VAR1, VAR2>, (LPVOID)(gjob + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nJob; i++) {
		WaitForSingleObject(job[i].hMMThread, INFINITE);
		//CloseHandle(job[iThread].hMMThread); job[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nJob; i++) {
		job[i].func = func;
		gjob[i].job = job + i;
		if (bSingleVar) {
			gjob[i].v1 = v1;
			gjob[i].v2 = v2;
		}
		else {
			gjob[i].v1 = v1 + i;
			gjob[i].v2 = v2 + i;
		}
		pthread_create(&(job[i].thread), NULL, &General_Thread_2<JOB, VAR1, VAR2>, (void *)(gjob + i));
	}
	for (i = 0; i < nJob; i++) {
		if (job[i].thread != 0) pthread_join(job[i].thread, NULL); // wait the thread over
	}
#endif
	delete[] gjob;
};


template <class JOB, class VAR1, class VAR2, class VAR3> class _GENERAL_JOB_VAR3 {
public:
	JOB *job;
	VAR1* v1;
	VAR2* v2;
	VAR3* v3;
	_GENERAL_JOB_VAR3() {job = NULL; v1 = NULL; v2 = NULL; v3 = NULL;};
	~_GENERAL_JOB_VAR3() {job = NULL; v1 = NULL; v2 = NULL; v3 = NULL;};
};

#if _SYS_ == _WINDOWS_SYS_
template <class JOB, class VAR1, class VAR2, class VAR3> int General_Thread_3(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class JOB, class VAR1, class VAR2, class VAR3> void *General_Thread_3(void *void_par) { // thread-function
#endif
	_GENERAL_JOB_VAR3<JOB, VAR1, VAR2, VAR3> *par = (_GENERAL_JOB_VAR3<JOB, VAR1, VAR2, VAR3>*)void_par;
	typedef void (*JOB_OPERATOR)(JOB*, VAR1*, VAR2*, VAR3*);
	((JOB_OPERATOR)(par->job->func))(par->job, par->v1, par->v2, par->v3);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->job->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class JOB, class VAR1, class VAR2, class VAR3> void MTOperate3(void* func, JOB *job, int nJob, VAR1* var1, VAR2* var2, VAR3* var3, bool bSingleVar = false) {
	if (nJob == 1) {
		typedef void (*JOB_OPERATOR)(JOB*, VAR1*, VAR2*, VAR3*);
		((JOB_OPERATOR)(func))(job, var1, var2, var3);
		return;
	}

	int i;
	_GENERAL_JOB_VAR3<JOB, VAR1, VAR2, VAR3> *gjob = new _GENERAL_JOB_VAR3<JOB, VAR1, VAR2, VAR3>[nJob];
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nJob; i++) {
		ResetEvent(job[i].hMMThread);
		job[i].func = func;
		gjob[i].job = job + i;
		if (bSingleVar) {
			gjob[i].v1 = var1;
			gjob[i].v2 = var2;
			gjob[i].v3 = var3;
		}
		else {
			gjob[i].v1 = var1 + i;
			gjob[i].v2 = var2 + i;
			gjob[i].v3 = var3 + i;
		}
		
		job[i].thread = AfxBeginThread((AFX_THREADPROC)General_Thread_3<JOB, VAR1, VAR2, VAR3>, (LPVOID)(gjob + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nJob; i++) {
		WaitForSingleObject(job[i].hMMThread, INFINITE);
		//CloseHandle(job[iThread].hMMThread); job[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nJob; i++) {
		job[i].func = func;
		gjob[i].job = job + i;
		if (bSingleVar) {
			gjob[i].v1 = var1;
			gjob[i].v2 = var2;
			gjob[i].v3 = var3;
		}
		else {
			gjob[i].v1 = var1 + i;
			gjob[i].v2 = var2 + i;
			gjob[i].v3 = var3 + i;
		}
		pthread_create(&(job[i].thread), NULL, &General_Thread_3<JOB, VAR1, VAR2, VAR3>, (void *)(gjob + i));
	}
	for (i = 0; i < nJob; i++) {
		if (job[i].thread != 0) pthread_join(job[i].thread, NULL); // wait the thread over
	}
#endif
	delete[] gjob;
};



class JobIndx : public JOB_THREAD_VARS<JobIndx> {
public:
	int n1, n2;
	JobIndx() {n1 = 0; n2 = -1;};
};

void assignJobIndx(JobIndx *job, int nThreads, int n1, int n2);
void assignJobIndx1(JobIndx *job, int &nThreads, int n1, int n2, int nStep = 100, int nMaxThreads = MAX_THREADS);



void initTorsionPars(MMOLECULE* mm);


class ImageCluster : public BASIC_CLUSTER {
public:
	VECTOR3 r; // normally we use center of mass as the position
	BASIC_CLUSTER *s; // the source of the cluster
	ImageCluster() {vp = &r; vcm = &r; s = NULL;};
	~ImageCluster() {s = NULL;};
};
