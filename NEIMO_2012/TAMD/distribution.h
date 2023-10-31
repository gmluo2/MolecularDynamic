namespace _distribution_ {

#define DISTRIBT_INDX(dstb, xv, n) n = int((xv - dstb.x.v[0]) / dstb.dx + 0.5);
#define LInRange(x, xf, xt) (x >= xf && x < xt ? true : false)

template <class Tx, class Ty> class DISTRIBT {
public:
	_PROFILE_1D_<Tx> x;  // x axis
	_PROFILE_1D_<Ty> y;  // y axis
	_PROFILE_1D_<long> sn; // sampled number
	_PROFILE_1D_<float> eb; // error bar of the sample
	Tx dx, xf, xt;
	Ty vnormal; // this value is used for normalization
	DISTRIBT() {xf = 0; xt = 0; dx = 0; vnormal = 0;};
	void set_distribt_x(Tx from, Tx to, int npts) {
		xf = from; xt = to;
		dx = (to - from) / npts;
		x.set(npts + 1);
		for (int i = 0; i < x.N; i++) x.v[i] = from + dx * i;
		
		if (y.N != x.N) y.set(npts + 1);
		else y.set_v_zero();
		
		if (sn.N != x.N) sn.set(npts + 1);
		else sn.set_v_zero();
		
		if (eb.N != x.N) eb.set(npts + 1);
		else eb.set_v_zero();
	};
	bool sample(Tx x, bool skip_exclude = false) {
		int n = int((x - xf) / dx + 0.5);
		if (n < 0) {
			if (skip_exclude) return false;
			else n = 0;
		}
		else if (n > this->y.N - 1) {
			if (skip_exclude) return false;
			else n = this->y.N - 1;
		}
		this->y.v[n] += 1;
		this->sn.v[n] += 1;
		return true;
	};
	bool sample(Tx x, float w, bool skip_exclude = false) {
		int n = int((x - xf) / dx + 0.5);
		if (n < 0) {
			if (skip_exclude) return false;
			else n = 0;
		}
		else if (n > y.N - 1) {
			if (skip_exclude) return false;
			else n = y.N - 1;
		}
		this->y.v[n] += w;
		this->sn.v[n] += 1;
		return true;
	};
	bool add_sample(int ipt, bool skip_exclude = false) {
		if (ipt < 0 || ipt >= x.N) {
			if (skip_exclude) return false;
			else ipt = (ipt < 0 ? 0 : x.N - 1);
		}
		this->y.v[ipt] += 1;
		this->sn.v[ipt] += 1;
		return true;
	};
	bool add_sample(int ipt, float w, bool skip_exclude = false) {
		if (ipt < 0 || ipt >= x.N) {
			if (skip_exclude) return false;
			else ipt = (ipt < 0 ? 0 : x.N - 1);
		}
		this->y.v[ipt] += w;
		this->sn.v[ipt] += 1;
		return true;
	};
	void eb_raw() {
		double t = 0;
		int i;
		for (i = 0; i < y.N; i++) {
			if (sn.v[i] == 0) eb.v[i] = 0;
			else eb.v[i] = sqrt(double(sn.v[i])) * y.v[i] / sn.v[i];
		}
	};
	void reset_distribt() {
		y.set_v_zero(); sn.set_v_zero(); eb.set_v_zero();
	};
	void release() {
		x.release(); y.release(); eb.release(); sn.release();
	};
	~DISTRIBT() {};
	void save_x_y_eb(ofstream &out) {
		for (int i = 0; i < x.N; i++) out<<x.v[i]<<"  "<<y.v[i]<<"  "<<eb.v[i]<<endl;
	};
	void save_x_y(ofstream &out) {
		for (int i = 0; i < x.N; i++) out<<x.v[i]<<"  "<<y.v[i]<<endl;
	};
	void save_x_y_eb_s(ofstream &out) {
		for (int i = 0; i < x.N; i++) out<<x.v[i]<<"  "<<y.v[i]<<"  "<<eb.v[i]<<"  "<<sn.v[i]<<endl;
	};
};

template <class Tx_res, class Ty_res, class Tx, class Ty> void raw_distribt_combine(DISTRIBT<Tx_res, Ty_res> &res, DISTRIBT<Tx, Ty> &dt) {
	if (dt.x.N != res.x.N) return;
	int i;
	for (i = 0; i < dt.x.N; i++) {
		res.y.v[i] += dt.y.v[i];
		res.sn.v[i] += dt.sn.v[i];
	}
}

template <class Tx, class Ty> class NDISTRIBT {
public:
	DISTRIBT<Tx, Ty> *dt;
	int ndt;

	void reset() {
		if (dt != NULL) {delete[] dt; dt = NULL;}
		ndt = 0;
	};
	void set(int nDT) {
		reset(); if (nDT <= 0) return;
		dt = new DISTRIBT<Tx, Ty>[nDT]; this->ndt = nDT;
	};
	void reset_distribt() {
		for (int i = 0; i < ndt; i++) dt[i].reset_distribt(); // all y axises are set to 0
	};
	NDISTRIBT() {dt = NULL; ndt = 0;};
	~NDISTRIBT() {reset();};
};

inline void setup_ndt(NDISTRIBT<float, long> &dest, NDISTRIBT<float, float> &source) {
	dest.set(source.ndt);
	for (int i = 0; i < dest.ndt; i++) dest.dt[i].set_distribt_x(source.dt[i].xf, source.dt[i].xt, source.dt[i].x.N - 1);
};

inline void setup_ndt(NDISTRIBT<float, long> &dest, NDISTRIBT<float, long> &source) {
	dest.set(source.ndt);
	for (int i = 0; i < dest.ndt; i++) dest.dt[i].set_distribt_x(source.dt[i].xf, source.dt[i].xt, source.dt[i].x.N - 1);
};

inline void setup_ndt(NDISTRIBT<float, float> &dest, NDISTRIBT<float, float> &source) {
	dest.set(source.ndt);
	for (int i = 0; i < dest.ndt; i++) dest.dt[i].set_distribt_x(source.dt[i].xf, source.dt[i].xt, source.dt[i].x.N - 1);
};

}  // end of namespace _distribution_

namespace _cluster_distribution_ {
using namespace _distribution_;
/*****************************************************
********  DIPOLE DISTRIBUTION OF MONOMER      ********
******************************************************/

#define N_MU_DT      7
#define MU_ABS       0
#define MU_X         1
#define MU_Y         2
#define MU_Z         3
#define MU_XY        4
#define MU_THETA     5
#define MU_OMEGA     6

void monomer_dipole_distribution(MMOLECULE &mm, bool group, int monomer_uid, bool limit_range, float mu1, float mu2, DISTRIBT<float, long> *dt);

typedef void (*DISTRIBT_DIPOLE_PVOID)(MMOLECULE&, bool group, int, bool, float, float, DISTRIBT<float, long>*);

class DISTRIBT_MU_THREAD_VARS : public _MACROMOL_THREAD_VARS<MMOLECULE> {
public:
	bool monomer_group; // calculate dipole of a group of monomer ?
	int monomer_uid;
	bool limit_range;
	float mu1, mu2; // range of mu in statistic
	DISTRIBT<float, long> *dt;

	DISTRIBT_DIPOLE_PVOID func;

	DISTRIBT_MU_THREAD_VARS() {
		this->mm = NULL; func = NULL; monomer_uid = 0; monomer_group = false;
		this->limit_range = false; mu1 = 0; mu2 = 0;
	};
	void set_pars(MMOLECULE* m, bool group, int monomer_uid, DISTRIBT<float, long> *dt, int iThread) {
		this->monomer_group = group; this->monomer_uid = monomer_uid; this->dt = dt;
		_MACROMOL_THREAD_VARS<MMOLECULE>::set_pars(iThread, m);
	};
	void set_range(bool limit, float mu1, float mu2) {
		this->limit_range = limit; this->mu1 = mu1; this->mu2 = mu2;
	};
	void set_func(DISTRIBT_DIPOLE_PVOID func) {this->func = (DISTRIBT_DIPOLE_PVOID)func;};
	~DISTRIBT_MU_THREAD_VARS() {dt = NULL;};
};

bool monomer_dipole_distribution_stat(char *mdpar_fname, bool group, int monomer_uid, int nf, int nt, int nstep, bool limit_range, float r1, float r2, DISTRIBT<float, float> *dt, int &nconfs);
extern "C" bool dipole_distribution(char *parfname);


/*****************************************************
********       ORIENTATION OF MONOMER         ********
******************************************************/

struct MONOMER_ATOM_INDX {
	int nc, na;
};

struct VECT_IN_MONOMER {
	MONOMER_ATOM_INDX atom[2];
};

#define N_ORT_DT     2
#define ORT_THETA    0
#define ORT_OMEGA    1

void monomer_orientation_distribution(MMOLECULE &mm, int monomer_uid, VECT_IN_MONOMER*, DISTRIBT<float, long> *dt);

typedef void (*DISTRIBT_ORT_PVOID)(MMOLECULE&, int, VECT_IN_MONOMER*, DISTRIBT<float, long>*);

class DISTRIBT_ORT_THREAD_VARS : public _MACROMOL_THREAD_VARS<MMOLECULE> {
public:
	int monomer_uid;
	VECT_IN_MONOMER mv;
	DISTRIBT<float, long> *dt;

	DISTRIBT_ORT_PVOID func;

	DISTRIBT_ORT_THREAD_VARS() {
		func = NULL; monomer_uid = 0;
	};
	void set_pars(MMOLECULE* m, int monomer_uid, VECT_IN_MONOMER &mv, DISTRIBT<float, long> *dt, int iThread) {
		this->monomer_uid = monomer_uid; this->dt = dt;
		this->mv.atom[0].nc = mv.atom[0].nc; this->mv.atom[0].na = mv.atom[0].na;
		this->mv.atom[1].nc = mv.atom[1].nc; this->mv.atom[1].na = mv.atom[1].na;
		_MACROMOL_THREAD_VARS<MMOLECULE>::set_pars(iThread, m);
	};
	void set_func(DISTRIBT_ORT_PVOID func) {this->func = (DISTRIBT_ORT_PVOID)func;};
	~DISTRIBT_ORT_THREAD_VARS() {dt = NULL;};
};

bool monomer_orientation_distribution_stat(char *mdpar_fname, int monomer_uid, VECT_IN_MONOMER* mv, int nf, int nt, int nstep, DISTRIBT<float, float> *dt, int &nconfs);
extern "C" bool orientation_distribution(char *parfname);


/*****************************************************************
**** RELATIVE ORIENTATION BETWEEN TWO CLUSTER IN SAME MONOMER ****
******************************************************************/

// relative orientation between two directions (given with atom's cordinates)
#define N_ORT2_DT    1
#define ORT2_THETA   0
#define ORT2_OMEGA   1 // two vector does not have omega angle

void cluster_relative_orientation_distribution(MMOLECULE &mm, int monomer_uid, VECT_IN_MONOMER* mv, DISTRIBT<float, long> *dt);

typedef void (*DISTRIBT_RELORT_PVOID)(MMOLECULE&, int, VECT_IN_MONOMER*, DISTRIBT<float, long>*);

class DISTRIBT_RELORT_THREAD_VARS : public _MACROMOL_THREAD_VARS<MMOLECULE> {
public:
	int monomer_uid;
	VECT_IN_MONOMER mv[2];
	DISTRIBT<float, long> *dt;

	DISTRIBT_RELORT_PVOID func;

	DISTRIBT_RELORT_THREAD_VARS() {
		func = NULL; monomer_uid = 0;
	};
	void set_pars(MMOLECULE* m, int monomer_uid, VECT_IN_MONOMER *mv, DISTRIBT<float, long> *dt, int iThread) {
		this->monomer_uid = monomer_uid; this->dt = dt;
		this->mv[0].atom[0].nc = mv[0].atom[0].nc; this->mv[0].atom[0].na = mv[0].atom[0].na;
		this->mv[0].atom[1].nc = mv[0].atom[1].nc; this->mv[0].atom[1].na = mv[0].atom[1].na;
		this->mv[1].atom[0].nc = mv[1].atom[0].nc; this->mv[1].atom[0].na = mv[1].atom[0].na;
		this->mv[1].atom[1].nc = mv[1].atom[1].nc; this->mv[1].atom[1].na = mv[1].atom[1].na;
		_MACROMOL_THREAD_VARS<MMOLECULE>::set_pars(iThread, m);
	};
	void set_func(DISTRIBT_RELORT_PVOID func) {this->func = (DISTRIBT_RELORT_PVOID)func;};
	~DISTRIBT_RELORT_THREAD_VARS() {dt = NULL;};
};

bool cluster_relative_orientation_distribution_stat(char *mdpar_fname, int monomer_uid, VECT_IN_MONOMER *mv, int nf, int nt, int nstep, DISTRIBT<float, float> *dt, int &nconfs);
extern "C" bool relative_orientation_distribution(char *parfname);


/*********************************************************************
*********  DIPOLE DIRECTION RELATIVE TO MONOMER ORIENTATION  *********
**********************************************************************/

// relative orientation between two directions (given with atom's cordinates)
#define N_LOCAL_MU_DT    1
#define LOCAL_MU_THETA   0
#define LOCAL_MU_OMEGA   1   // is omega reasonable ?

void dipole_local_orientation_distribution(MMOLECULE &mm, int monomer_uid, VECT_IN_MONOMER* mv, bool limit, float r1, float r2, DISTRIBT<float, long> *dt);

typedef void (*DISTRIBT_LOCAL_DIPOLE_PVOID)(MMOLECULE&, int, VECT_IN_MONOMER*, bool, float, float, DISTRIBT<float, long>*);

class DISTRIBT_LOCAL_DIPOLE_THREAD_VARS : public _MACROMOL_THREAD_VARS<MMOLECULE> {
public:
	int monomer_uid;
	VECT_IN_MONOMER mv;
	DISTRIBT<float, long> *dt;

	bool limit;
	float r1, r2;

	DISTRIBT_LOCAL_DIPOLE_PVOID func;

	DISTRIBT_LOCAL_DIPOLE_THREAD_VARS() {
		func = NULL; monomer_uid = 0; dt = NULL; limit = false; r1 = 0; r2 = 0;
	};
	void set_pars(MMOLECULE* m, int monomer_uid, VECT_IN_MONOMER &mv, DISTRIBT<float, long> *dt, int iThread) {
		this->monomer_uid = monomer_uid; this->dt = dt;
		this->mv.atom[0].nc = mv.atom[0].nc; this->mv.atom[0].na = mv.atom[0].na;
		this->mv.atom[1].nc = mv.atom[1].nc; this->mv.atom[1].na = mv.atom[1].na;
		_MACROMOL_THREAD_VARS<MMOLECULE>::set_pars(iThread, m);
	};
	void set_limit(bool limit, float r1, float r2) {
		this->limit = limit; this->r1 = r1; this->r2 = r2;
	};
	void set_func(DISTRIBT_LOCAL_DIPOLE_PVOID func) {this->func = (DISTRIBT_LOCAL_DIPOLE_PVOID)func;};
	~DISTRIBT_LOCAL_DIPOLE_THREAD_VARS() {dt = NULL;};
};

bool dipole_local_orientation_distribution_stat(char *mdpar_fname, int monomer_uid, VECT_IN_MONOMER *mv, bool limit, float r1, int nf, int nt, int nstep, DISTRIBT<float, float> *dt, int &nconfs);
extern "C" bool dipole_local_orientation(char *parfname);



/*********************************************************************
***************       TORSION ANGLE DISTRIBUTION       ***************
**********************************************************************/

template <class Tx, class Ty> class TA_DISTRIBT {
public:
	int monomer_uid, cindx;
	char title[256];
	//DISTRIBT<float, long> dt;
	DISTRIBT<Tx, Ty> dt;

	TA_DISTRIBT() {monomer_uid = 0; cindx = 0; strcpy(title, "\0");};
	~TA_DISTRIBT() {dt.release();};
};

void cluster_torsion_angle_distribution(MMOLECULE &mm, TA_DISTRIBT<float, long> *distb, int ndistb);

typedef void (*DISTRIBT_CLUSTER_TORSION_ANGLE_PVOID)(MMOLECULE&, TA_DISTRIBT<float, long>*, int);

class DISTRIBT_CLUSTER_TA_THREAD_VARS : public _MACROMOL_THREAD_VARS<MMOLECULE> {
public:
	TA_DISTRIBT<float, long> *distb;
	int ndistb;

	DISTRIBT_CLUSTER_TORSION_ANGLE_PVOID func;

	DISTRIBT_CLUSTER_TA_THREAD_VARS() {
		func = NULL; distb = NULL; ndistb = 0;
	};
	void set_pars(MMOLECULE* m, TA_DISTRIBT<float, long> *distb, int ndistb, int iThread) {
		this->distb = distb; this->ndistb = ndistb;
		_MACROMOL_THREAD_VARS<MMOLECULE>::set_pars(iThread, m);
	};
	void set_func(DISTRIBT_CLUSTER_TORSION_ANGLE_PVOID func) {this->func = (DISTRIBT_CLUSTER_TORSION_ANGLE_PVOID)func;};
	~DISTRIBT_CLUSTER_TA_THREAD_VARS() {distb = NULL; ndistb = 0;};
};

bool cluster_torsion_angle_distribution_stat(char *mdpar_fname, int nf, int nt, int nstep, TA_DISTRIBT<float, float> *distb, int ndistb, int &nconfs);



/*****************************************************
********  SPATIAL DISTRIBUTION OF CLUSTER       ********
******************************************************/

#define R_X         1
#define R_Y         2
#define R_Z         3

void cluster_spatial_distribution(ARRAY<BASIC_CLUSTER*> &bc, ARRAY<int>& c_uid, int axis, bool bLimitRange, DISTRIBT<float, long> *dt);

typedef void (*SPATIAL_DISTRIBT_PVOID)(ARRAY<BASIC_CLUSTER*>&, ARRAY<int>&, int, bool, DISTRIBT<float, long>*);

class SPATIAL_DISTRIBT_THREAD_VARS : public _THREAD_VARS {
public:
	ARRAY<BASIC_CLUSTER*> *bc;
	ARRAY<int> c_uid;
	int axis;
	bool bLimitRange;
	DISTRIBT<float, long> *dt;

	SPATIAL_DISTRIBT_PVOID func;

	SPATIAL_DISTRIBT_THREAD_VARS() {
		this->bc = NULL; this->dt = NULL; func = NULL; axis = R_Z; this->bLimitRange = false; 
	};
	void set_pars(ARRAY<BASIC_CLUSTER*>* bc, ARRAY<int>& c_uid, int axis, DISTRIBT<float, long> *dt, int iThread) {
		this->bc = bc; this->c_uid.set_array(c_uid.n); memcpy(this->c_uid.m, c_uid.m, c_uid.n * sizeof(int)); 
		this->axis = axis; this->dt = dt;
		_THREAD_VARS::set_pars(iThread);
	};
	void set_func(SPATIAL_DISTRIBT_PVOID func) {this->func = (SPATIAL_DISTRIBT_PVOID)func;};
	~SPATIAL_DISTRIBT_THREAD_VARS() {bc = NULL; dt = NULL;};
};

bool cluster_spatial_distribution_stat(char *mdpar_fname, int nf, int nt, int nstep, ARRAY<short> &c_uid, int axis, bool bLimitRange, DISTRIBT<float, float> *dt, int &nconfs);

} // end of namespace _cluster_distribution

extern "C" bool cluster_torsion_angle(char *parfname);
extern "C" bool cluster_spatial_distribution(char *parfname);

bool get_monomer_uid(char *monomer, int &uid);
bool get_cluster_uid(char *cluster, int &uid);
bool read_given_id(char *fname, char* title, ARRAY<short>& id_cluster);
bool read_given_id(char *fname, char* title, ARRAY<short>& id_cluster, bool (*get_uid)(char*, int&), char *end_str);
