// in interaction calculation, we do not calculate the direct interaction
// when the distance of two atoms are closer than _Rmin_Interact
// for electrostatic interaction, real part of electrostatic interaction are ignored when R < _Rmin_Interact
#define _Rmin_Interact   0.2
#define _R2min_Interact  0.04


//#if INTERACT_DEF == 1

// in amber force field, LJ interaction is defined as eps * [(rLJ / r)^-12 - 2 * (r_LJ / r)^-6]
// when r/rLJ = 0.72, (rLJ/r)^2 = 1.929, (rLJ/r)^12 - 2 * (rLJ/r)^6 = 37.2
// when r/rLJ = 0.7, (rLJ/r)^2 = 2.0408, (rLJ/r)^12 - 2 * (rLJ/r)^6 = 55.28

#define _LJ_V3FORCE(eps, rLJ2, r, dis2, t, f, V) \
	t = rLJ2 / dis2; \
	t = t * t * t; \
	V = eps * (t * t - t - t) * ELJ_unit; \
	t = 12 * eps * (t * t - t) / dis2 * FLJ_unit; \
	f.v[0] = t * r.v[0]; f.v[1] = t * r.v[1]; f.v[2] = t * r.v[2];  

#define SR_FORCE_MAX   1

// disable the force limitation -- if (f.v[0] > 1e4) f.v[0] = 1e4; 
#define LJ_V3FORCE(eps, rLJ, dis, u, i, t1, t2, V, f) \
	t1 = dis / rLJ; interpolate2(V, f.v[0], LJ12_6, t1, i, t2) \
	V *= eps; t2 = f.v[0] * eps / rLJ; \
	/*if (t2 > SR_FORCE_MAX) t2 = SR_FORCE_MAX; 	else if (t2 > force_max && scale_force) t2 = force_max; */\
	f.v[0] = t2 * u.v[0]; f.v[1] = t2 * u.v[1]; f.v[2] = t2 * u.v[2];  

#define LJ_V(eps, rLJ, dis, i, t, V) \
	t = dis / rLJ; interpolate1(V, LJ12_6, t, i, 0) \
	V *= eps; 

//#endif

#define T1(r, u, i, dis3) for (i = 0; i < u.n; i++) u.v[i] = -r.v[i] / dis3;

#define T2(r, mt, i, j, dis, dis2, dis5, temp) for (i = 0; i < mt.n; i++) {for (j = 0; j < mt.n; j++) {\
		temp = 3 * r.v[i] * r.v[j]; if (i == j) temp -= dis2; mt.m[i][j] = temp / dis5; \
	}}

#define T3(r, cm, i, j, k, dis, dis2, dis7, temp) \
	for (i = 0; i < cm.n; i++) {for (j = 0; j < cm.n; j++) {for(k = 0; k < cm.n; k++) { \
		cm.m[i][j][k] = 15 * r.v[i] * r.v[j] * r.v[k]; temp = 0.0; \
		if (i == j) temp += r.v[k]; if (i == k) temp += r.v[j]; if (j == k) temp += r.v[i]; \
		cm.m[i][j][k] -= 3 * temp * dis2; \
		cm.m[i][j][k] /= -(dis7); \
	}}}

#define T4(r, qm, i, j, k, l, dis, dis2, dis4, dis9, temp) \
	for (i = 0; i < qm.n; i++) {for (j = 0; j < qm.n; j++) {for(k = 0; k < qm.n; k++) {for (l = 0; l < qm.n; l++) { \
		qm.m[i][j][k][l] = 105 * r.v[i] * r.v[j] * r.v[k] * r.v[l]; temp = 0; \
		if (i == j) temp += r.v[k] * r.v[l]; if (i == k) temp += r.v[j] * r.v[l]; \
		if (i == l) temp += r.v[j] * r.v[k]; if (j == k) temp += r.v[i] * r.v[l]; \
		if (j == l) temp += r.v[i] * r.v[k]; if (k == l) temp += r.v[i] * r.v[j]; \
		if (temp != 0) {temp *= 15 * dis2; qm.m[i][j][k][l] -= temp;} \
		temp = 0; \
		if (i == j && k == l) temp += 1; \
		if (i == k && j == l) temp += 1; \
		if (i == l && j == k) temp += 1; \
		if (temp != 0) qm.m[i][j][k][l] += 3 * temp * dis4; \
		qm.m[i][j][k][l] /= dis9; \
	}}}}

#define MTP_T1(u, r, R3) u.v[0] = -r.v[0] / R3; u.v[1] = -r.v[1] / R3; u.v[2] = -r.v[2] / R3;
#define MTP_T2(m2, r, R2, R5, i, j) for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {\
	m2.m[i][j] = 3 * r.v[i] * r.v[j]; if (i == j) m2.m[i][j] -= R2; m2.m[i][j] /= R5; \
	}}

#define MTP_T3(m3, r, R2, R7, i, j, k, t) t = 3 * R2; for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {\
	m3.m[i][j][k] = 15 * r.v[i] * r.v[j] * r.v[k]; \
	if (i == j) m3.m[i][j][k] -= t * r.v[k]; if (i == k) m3.m[i][j][k] -= t * r.v[j]; \
	if (j == k) m3.m[i][j][k] -= t * r.v[i]; m3.m[i][j][k] /= -R7; \
	}}}

#define MTP_T4(m4, r, R2, R4, R9, i, j, k, l, t1, t2) t1 = 15 * R2; t2 = 3 * R4; \
	for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {for (l = 0; l < 3; l++) {\
	m4.m[i][j][k][l] = 105 * r.v[i] * r.v[j] * r.v[k] * r.v[l]; \
	if (i == j) m4.m[i][j][k][l] -= t1 * r.v[k] * r.v[l]; \
	if (i == k) m4.m[i][j][k][l] -= t1 * r.v[j] * r.v[l]; \
	if (i == l) m4.m[i][j][k][l] -= t1 * r.v[j] * r.v[k]; \
	if (j == k) m4.m[i][j][k][l] -= t1 * r.v[i] * r.v[l]; \
	if (j == l) m4.m[i][j][k][l] -= t1 * r.v[i] * r.v[k]; \
	if (k == l) m4.m[i][j][k][l] -= t1 * r.v[i] * r.v[j]; \
	if (i == j && k == l) m4.m[i][j][k][l] += t2; \
	if (i == k && j == l) m4.m[i][j][k][l] += t2; \
	if (i == l && j == k) m4.m[i][j][k][l] += t2; \
	m4.m[i][j][k][l] /= R9; \
	}}}}

#define MTP_D2(m2, r, i, j) for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {\
	m2.m[i][j] = r.v[i] * r.v[j]; \
	}}

#define _PAIRWISE  0
#define _CMM       1

#define SOLVENT_GAUSS_FORCE(eps, ur, d, f, t, i) t = (d - dSolvent) / sigmaSolvent; \
	if (t >= 4) {f.v[0] = 0; f.v[1] = 0; f.v[2] = 0;} \
	else {\
		if (t < 0) f.v[0] = 1 - 2 * t / sigmaSolvent; \
		else if (t <= 0.01) f.v[0] = 1; else {EXP2(f.v[0], t, i)} \
		f.v[0] *= eps; f.v[1] = f.v[0] * ur.v[1]; f.v[2] = f.v[0] * ur.v[2]; f.v[0] *= ur.v[0]; \
	}

class InteractVar {
public:
	bool bVirial; // calculate virial / stress tensor ?
	bool bST_diagonal; // calculate the diagonal elements, xx, yy and zz, or the trace only ?
	InteractVar() {bVirial = false; bST_diagonal = false;};
	InteractVar(double v) {};
	void operator = (const InteractVar &var) {
		this->bVirial = var.bVirial;
		this->bST_diagonal = var.bST_diagonal;
	};
};

class EInteractVar : public InteractVar {
public:
	bool bEfieldOnly;
	_EwaldSum_real_::EwaldSumRealVars1 esv1;
	bool bTholePolarize;
	_Thole_polarize_::TholeCorrTensor tv;

	EInteractVar() {
		bEfieldOnly = false; 
		bTholePolarize = false;
	};
	void operator = (const EInteractVar &var) {
		this->bEfieldOnly = var.bEfieldOnly;
		(*(InteractVar*)this) = (*(InteractVar*)(&var));
		esv1.cp_EwaldSum_vars(*((_EwaldSum_real_::EwaldSumRealVars0*)(&esv1)));
		this->bTholePolarize = var.bTholePolarize;
		this->tv.iTholeModel = var.tv.iTholeModel;
		this->tv.xcut_exp = var.tv.xcut_exp;
	};
};

class InteractRes {
public:
	double U; // potential
	double STxx, STyy, STzz, STtr; // stress tensor's diagonal elements, and the trace

	void reset(){U = 0; STxx = 0; STyy = 0; STzz = 0; STtr = 0;};
	InteractRes() {reset();};
	InteractRes(double v) {U = v; STxx = v; STyy = v; STzz = v; STtr = v;};
	void operator = (const InteractRes &var) {
		U = var.U; STxx = var.STxx; STyy = var.STyy; STzz = var.STzz; STtr = var.STtr;
	};
	void operator = (double v) {
		U = v; STxx = v; STyy = v; STzz = v; STtr = v;
	};
	InteractRes operator + (InteractRes &var) {
		InteractRes res;
		res.U = U + var.U;
		res.STxx = STxx + var.STxx; 
		res.STyy = STyy + var.STyy; 
		res.STzz = STzz + var.STzz; 
		res.STtr = STtr + var.STtr;
		return res;
	};
	InteractRes operator - (InteractRes &var) {
		InteractRes res;
		res.U = U - var.U;
		res.STxx = STxx - var.STxx; 
		res.STyy = STyy - var.STyy; 
		res.STzz = STzz - var.STzz; 
		res.STtr = STtr - var.STtr;
		return res;
	};
	void operator += (InteractRes &var) {
		U += var.U;
		STxx += var.STxx; 
		STyy += var.STyy; 
		STzz += var.STzz; 
		STtr += var.STtr;
	};
	void operator -= (InteractRes &var) {
		U -= var.U;
		STxx -= var.STxx; 
		STyy -= var.STyy; 
		STzz -= var.STzz; 
		STtr -= var.STtr;
	};
};

inline void Accumulate_Virial_Diagonal(InteractRes &res, double *f, double *dr) {
	double v[3] = {0.5 * f[0] * dr[0], 0.5 * f[1] * dr[1], 0.5 * f[2] * dr[2]};
	res.STxx += v[0]; res.STyy += v[1]; res.STzz += v[2];
	res.STtr += v[0] + v[1] + v[2];
};

inline void Accumulate_Virial_Diagonal_full(InteractRes &res, double *f, double *dr) {
	double v[3] = {f[0] * dr[0], f[1] * dr[1], f[2] * dr[2]};
	res.STxx += v[0]; res.STyy += v[1]; res.STzz += v[2];
	res.STtr += v[0] + v[1] + v[2];
};

inline void Accumulate_Virial_Trace(InteractRes &res, double *f, double *dr) {
	double v[3] = {0.5 * f[0] * dr[0], 0.5 * f[1] * dr[1], 0.5 * f[2] * dr[2]};
	res.STtr += v[0] + v[1] + v[2];
};

inline void Accumulate_Virial_Trace(InteractRes &res, double f, double dr) {
	res.STtr += f * dr * 0.5;
};

inline void Accumulate_Virial_Trace_full(InteractRes &res, double *f, double *dr) {
	double v[3] = {f[0] * dr[0], f[1] * dr[1], f[2] * dr[2]};
	res.STtr += v[0] + v[1] + v[2];
};

inline void Accumulate_Virial_Trace_full(InteractRes &res, double f, double dr) {
	res.STtr += f * dr;
};

/***************************************************************
     Multi-Thread Interaction Calculation based on 
	             CELL MULTIPOLE METHOD
*****************************************************************/

template <class CMM_CLUSTER, class VAR, class RES> class CMM_THREAD_VARS : public _THREAD_VARS {
public:
	CMM_CELL3D<CMM_CLUSTER> *cell;
	int n1Cell, n2Cell; // the interaction of the cell index in the range [n1Cell, n2Cell]
	VAR var;
	RES res;

	void* InteractFunc;

	CMM_THREAD_VARS() {
		cell = NULL; InteractFunc = NULL;
	};
	void set_pars(int i, CMM_CELL3D<CMM_CLUSTER> *c, int n1, int n2, VAR &var) {
		cell = c; n1Cell = n1; n2Cell = n2; this->var = var;
		_THREAD_VARS::set_pars(i);
	};
	void set_func(void* func) {this->InteractFunc = (void*)func;};
	~CMM_THREAD_VARS() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class CMM_CLUSTER, class VAR, class RES> int ClusterInteractionThread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class CMM_CLUSTER, class VAR, class RES> void *ClusterInteractionThread(void *void_par) { // thread-function
#endif
	typedef void (*CLUSTER_INTERACTION_FUNC)(CMM_CELL3D<CMM_CLUSTER>&, int, int, VAR&, RES&);
	CMM_THREAD_VARS<CMM_CLUSTER, VAR, RES> *par = (CMM_THREAD_VARS<CMM_CLUSTER, VAR, RES>*)void_par;
	((CLUSTER_INTERACTION_FUNC)(par->InteractFunc))(*(par->cell), par->n1Cell, par->n2Cell, par->var, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};
// MULTI-THREAD function to calculate the interaction
template <class CMM_CLUSTER, class VAR, class RES> void ClusterInteraction(void* func, CMM_CELL3D<CMM_CLUSTER> &cell, int n1Cell, int n2Cell, VAR& var, RES& res, int nThreads) { // V in kT
	if (nThreads <= 1) {
		typedef void (*CLUSTER_INTERACTION_FUNC)(CMM_CELL3D<CMM_CLUSTER>&, int, int, VAR&, RES&);
		res = (RES)0;
		((CLUSTER_INTERACTION_FUNC)(func))(cell, n1Cell, n2Cell, var, res);
		return;
	}

	if (n2Cell >= cell.acell.n) n2Cell = cell.acell.n - 1;
	int nc = 0;
	CMM_THREAD_VARS<CMM_CLUSTER, VAR, RES> *thread_vars = new CMM_THREAD_VARS<CMM_CLUSTER, VAR, RES>[nThreads];
	int iThread = 0;
	int n1 = 0, n2 = 0;
	int nstep = (n2Cell - n1Cell + 1) / nThreads; n1 = n1Cell;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = n2Cell;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &cell, n1, n2, var);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)ClusterInteractionThread<CMM_CLUSTER, VAR, RES>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &ClusterInteractionThread<CMM_CLUSTER, VAR, RES>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	res = (RES)0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		res += thread_vars[iThread].res;
		thread_vars[iThread].cell = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
};

/***************************************************************
     Multi-Thread Interaction Calculation of Cluster 
	          in MacroMolecule, with CMM
*****************************************************************/

template <class MM, class CMM_CLUSTER, class VAR, class RES> class MM_CMM_THREAD_VARS : public _THREAD_VARS {
public:
	CMM_CELL3D<CMM_CLUSTER> *cell;
	MM *mm;
	VAR var;
	RES res;
	int n1, n2; // the interaction of the cluster index in the range [n1, n2]

	void* InteractFunc;

	MM_CMM_THREAD_VARS() {
		cell = NULL; InteractFunc = NULL;
	};
	void set_pars(int i, CMM_CELL3D<CMM_CLUSTER> *c, MM* mm, int n1, int n2, VAR &var) {
		cell = c; this->mm = mm; this->n1 = n1; this->n2 = n2; this->var = var;
		_THREAD_VARS::set_pars(i);
	};
	void set_func(void* func) {this->InteractFunc = (void*)func;};
	~MM_CMM_THREAD_VARS() {};
};


#if _SYS_ == _WINDOWS_SYS_
template <class MM, class CMM_CLUSTER, class VAR, class RES> int MMClusterInteractionThread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MM, class CMM_CLUSTER, class VAR, class RES> void *MMClusterInteractionThread(void *void_par) { // thread-function
#endif
	typedef void (*MM_INTERACTION_FUNC)(CMM_CELL3D<CMM_CLUSTER>&, MM&, int, int, VAR&, RES&);
	MM_CMM_THREAD_VARS<MM, CMM_CLUSTER, VAR, RES> *par = (MM_CMM_THREAD_VARS<MM, CMM_CLUSTER, VAR, RES>*)void_par;
	((MM_INTERACTION_FUNC)(par->InteractFunc))(*(par->cell), *(par->mm), par->n1, par->n2, par->var, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

// MULTI-THREAD function to calculate the interaction
template <class MM, class CMM_CLUSTER, class VAR, class RES> void MMClusterInteraction(void* func, CMM_CELL3D<CMM_CLUSTER> &cell, MM &mm, int nc1, int nc2, int nThreads, VAR &var, RES &res) { // V in kT
	if (nThreads <= 1) {
		typedef void (*MM_INTERACTION_FUNC)(CMM_CELL3D<CMM_CLUSTER>&, MM&, int, int, VAR&, RES&);
		res = (RES)0;
		((MM_INTERACTION_FUNC)(func))(cell, mm, nc1, nc2, var, res);
		return;
	}

	if (nc2 >= mm.nCluster) nc2 = mm.nCluster - 1;
	MM_CMM_THREAD_VARS<MM, CMM_CLUSTER, VAR, RES> *thread_vars = new MM_CMM_THREAD_VARS<MM, CMM_CLUSTER, VAR, RES>[nThreads];
	int iThread = 0;
	int n1 = 0, n2 = 0;
	int nstep = (nc2 - nc1 + 1) / nThreads; n1 = nc1;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = nc2;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &cell, &mm, n1, n2, var);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MMClusterInteractionThread<MM, CMM_CLUSTER, VAR, RES>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MMClusterInteractionThread<MM, CMM_CLUSTER, VAR, RES>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	res = (RES)0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		res += thread_vars[iThread].res;
		thread_vars[iThread].cell = NULL; thread_vars[iThread].mm = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};

/****************************************************************************
     Multi-Thread Interaction Calculation of Small Molecule, with CMM
*****************************************************************************/

template <class M, class CMM_CLUSTER, class VAR, class RES> class M_CMM_THREAD_VARS : public _THREAD_VARS {
public:
	CMM_CELL3D<CMM_CLUSTER> *cell;
	M *m;
	int nM;
	VAR var;
	RES res;

	void* InteractFunc;

	M_CMM_THREAD_VARS() {
		cell = NULL; InteractFunc = NULL;
	};
	void set_pars(int i, CMM_CELL3D<CMM_CLUSTER> *c, M* m, int nM, VAR &var) {
		cell = c; this->m = m; this->nM = nM; this->var = var;
		_THREAD_VARS::set_pars(i);
	};
	void set_func(void* func) {this->InteractFunc = (void*)func;};
	~M_CMM_THREAD_VARS() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class M, class CMM_CLUSTER, class VAR, class RES> int MInteractionThread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class M, class CMM_CLUSTER, class VAR, class RES> void *MInteractionThread(void *void_par) { // thread-function
#endif
	typedef void (*M_INTERACTION_FUNC)(CMM_CELL3D<CMM_CLUSTER>&, M*, int, VAR&, RES&);
	M_CMM_THREAD_VARS<M, CMM_CLUSTER, VAR, RES> *par = (M_CMM_THREAD_VARS<M, CMM_CLUSTER, VAR, RES>*)void_par;
	((M_INTERACTION_FUNC)(par->InteractFunc))(*(par->cell), par->m, par->nM, par->var, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class M, class CMM_CLUSTER, class VAR, class RES> void MultiMolInteraction(void* func, CMM_CELL3D<CMM_CLUSTER> &cell, M *m, int nM, int nThreads, VAR& var, RES& res) {
	if (nThreads <= 1) {
		typedef void (*M_INTERACTION_FUNC)(CMM_CELL3D<CMM_CLUSTER>&, M*, int, VAR&, RES&);
		res = (RES)0;
		((M_INTERACTION_FUNC)(func))(cell, m, nM, var, res);
		return;
	}

	M_CMM_THREAD_VARS<M, CMM_CLUSTER, VAR, RES> *thread_vars = new M_CMM_THREAD_VARS<M, CMM_CLUSTER, VAR, RES>[nThreads];
	int iThread = 0;
	int nstart = 0, nm = 0;
	int nstep = nM / nThreads;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) nm = nM - nstart;
		else nm = nstep;

		thread_vars[iThread].set_pars(iThread, &cell, m + nstart, nm, var);
		thread_vars[iThread].set_func(func);
		nstart += nm;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MInteractionThread<M, CMM_CLUSTER, VAR, RES>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MInteractionThread<M, CMM_CLUSTER, VAR, RES>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	res = (RES)0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		res += thread_vars[iThread].res;
		thread_vars[iThread].cell = NULL; thread_vars[iThread].m = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};




/****************************************************************************
     Multi-Thread Interaction Calculation of atom, within a system
*****************************************************************************/

template <class SYS, class ATOM, class VAR, class RES> class SYS_THREAD_VARS : public _THREAD_VARS {
public:
	SYS *sys;
	ATOM *m;
	int nM;
	VAR var;
	RES res;

	void* func;

	SYS_THREAD_VARS() {
		sys = NULL; func = NULL; m = NULL; nM = 0;
	};
	void set_pars(int i, SYS *sys, ATOM* m, int nM, VAR &var) {
		this->sys = sys; this->m = m; this->nM = nM; this->var = var;
		_THREAD_VARS::set_pars(i);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~SYS_THREAD_VARS() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class SYS, class ATOM, class VAR, class RES> int MRunSysThread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class SYS, class ATOM, class VAR, class RES> void *MRunSysThread(void *void_par) { // thread-function
#endif
	typedef void (*MRUN_SYS_FUNC)(SYS&, ATOM*, int, VAR&, RES&);
	SYS_THREAD_VARS<SYS, ATOM, VAR, RES> *par = (SYS_THREAD_VARS<SYS, ATOM, VAR, RES>*)void_par;
	((MRUN_SYS_FUNC)(par->func))(*(par->sys), par->m, par->nM, par->var, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class SYS, class ATOM, class VAR, class RES> void MRunSys(void* func, SYS &sys, ATOM *m, int nM, int nThreads, VAR& var, RES& res) {
	if (nThreads <= 1) {
		typedef void (*MRUN_SYS_FUNC)(SYS&, ATOM*, int, VAR&, RES&);
		res = (RES)0;
		((MRUN_SYS_FUNC)(func))(sys, m, nM, var, res);
		return;
	}

	SYS_THREAD_VARS<SYS, ATOM, VAR, RES> *thread_vars = new SYS_THREAD_VARS<SYS, ATOM, VAR, RES>[nThreads];
	int iThread = 0;
	int nstart = 0, nm = 0;
	int nstep = nM / nThreads;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) nm = nM - nstart;
		else nm = nstep;

		thread_vars[iThread].set_pars(iThread, &sys, m + nstart, nm, var);
		thread_vars[iThread].set_func(func);
		nstart += nm;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MRunSysThread<SYS, ATOM, VAR, RES>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MRunSysThread<SYS, ATOM, VAR, RES>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	res = (RES)0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		res += thread_vars[iThread].res;
		thread_vars[iThread].sys = NULL; thread_vars[iThread].m = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};




template <class SYS, class VAR, class RES> class SYS1_THREAD_VARS : public _THREAD_VARS {
public:
	SYS *sys;
	int n1, n2;
	VAR var;
	RES res;

	void* func;

	SYS1_THREAD_VARS() {
		sys = NULL; func = NULL; n1 = 0; n2 = 0;
	};
	void set_pars(int i, SYS *sys, int n1, int n2, VAR &var) {
		this->sys = sys; this->n1 = n1; this->n2 = n2; this->var = var;
		_THREAD_VARS::set_pars(i);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~SYS1_THREAD_VARS() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class SYS, class VAR, class RES> int MRunSysThread_1(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class SYS, class VAR, class RES> void *MRunSysThread_1(void *void_par) { // thread-function
#endif
	typedef void (*MRUN_SYS_FUNC)(SYS&, int, int, VAR&, RES&);
	SYS1_THREAD_VARS<SYS, VAR, RES> *par = (SYS1_THREAD_VARS<SYS, VAR, RES>*)void_par;
	((MRUN_SYS_FUNC)(par->func))(*(par->sys), par->n1, par->n2, par->var, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class SYS, class VAR, class RES> void MRunSys1(void* func, SYS &sys, int n1, int n2, int nThreads, VAR& var, RES& res) {
	if (nThreads <= 1) {
		typedef void (*MRUN_SYS_FUNC)(SYS&, int, int, VAR&, RES&);
		res = (RES)0;
		((MRUN_SYS_FUNC)(func))(sys, n1, n2, var, res);
		return;
	}

	int N = n2 - n1 + 1;

	SYS1_THREAD_VARS<SYS, VAR, RES> *thread_vars = new SYS1_THREAD_VARS<SYS, VAR, RES>[nThreads];
	int iThread = 0;
	int nstart = n1, nm = 0, nstop;
	int nstep = N / nThreads;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) {nm = N - nstart; nstop = n2;}
		else {nm = nstep; nstop = nstart + nm - 1;}

		thread_vars[iThread].set_pars(iThread, &sys, nstart, nstop, var);
		thread_vars[iThread].set_func(func);
		nstart += nm;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MRunSysThread_1<SYS, VAR, RES>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MRunSysThread_1<SYS, VAR, RES>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	res = (RES)0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		res += thread_vars[iThread].res;
		thread_vars[iThread].sys = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};


template <class SYS, class VAR> class SYS2_THREAD_VARS : public _THREAD_VARS {
public:
	SYS *sys;
	int n1, n2;
	VAR var;

	void* func;

	SYS2_THREAD_VARS() {
		sys = NULL; func = NULL; n1 = 0; n2 = 0;
	};
	void set_pars(int i, SYS *sys, int n1, int n2, VAR &var) {
		this->sys = sys; this->n1 = n1; this->n2 = n2; this->var = var;
		_THREAD_VARS::set_pars(i);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~SYS2_THREAD_VARS() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class SYS, class VAR> int MRunSysThread_2(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class SYS, class VAR> void *MRunSysThread_2(void *void_par) { // thread-function
#endif
	typedef void (*MRUN_SYS_FUNC)(SYS&, int, int, VAR&);
	SYS2_THREAD_VARS<SYS, VAR> *par = (SYS2_THREAD_VARS<SYS, VAR>*)void_par;
	((MRUN_SYS_FUNC)(par->func))(*(par->sys), par->n1, par->n2, par->var);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class SYS, class VAR> void MRunSys2(void* func, SYS &sys, int n1, int n2, int nThreads, VAR& var) {
	if (nThreads <= 1) {
		typedef void (*MRUN_SYS_FUNC)(SYS&, int, int, VAR&);
		((MRUN_SYS_FUNC)(func))(sys, n1, n2, var);
		return;
	}

	int N = n2 - n1 + 1;

	SYS2_THREAD_VARS<SYS, VAR> *thread_vars = new SYS2_THREAD_VARS<SYS, VAR>[nThreads];
	int iThread = 0;
	int nstart = n1, nm = 0, nstop;
	int nstep = N / nThreads;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) {nm = N - nstart; nstop = n2;}
		else {nm = nstep; nstop = nstart + nm - 1;}

		thread_vars[iThread].set_pars(iThread, &sys, nstart, nstop, var);
		thread_vars[iThread].set_func(func);
		nstart += nm;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MRunSysThread_2<SYS, VAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MRunSysThread_2<SYS, VAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].sys = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};


template <class SYS, class RES> class SYS3_THREAD_VARS : public _THREAD_VARS {
public:
	SYS *sys;
	int n1, n2;
	RES res;

	void* func;

	SYS3_THREAD_VARS() {
		sys = NULL; func = NULL; n1 = 0; n2 = 0;
	};
	void set_pars(int i, SYS *sys, int n1, int n2) {
		this->sys = sys; this->n1 = n1; this->n2 = n2;
		_THREAD_VARS::set_pars(i);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~SYS3_THREAD_VARS() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class SYS, class RES> int MRunSysThread_3(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class SYS, class RES> void *MRunSysThread_3(void *void_par) { // thread-function
#endif
	typedef void (*MRUN_SYS_FUNC)(SYS&, int, int, RES&);
	SYS3_THREAD_VARS<SYS, RES> *par = (SYS3_THREAD_VARS<SYS, RES>*)void_par;
	((MRUN_SYS_FUNC)(par->func))(*(par->sys), par->n1, par->n2, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class SYS, class RES> void MRunSys3(void* func, SYS &sys, int n1, int n2, int nThreads, RES& res) {
	if (nThreads <= 1) {
		typedef void (*MRUN_SYS_FUNC)(SYS&, int, int, RES&);
		res = (RES)0;
		((MRUN_SYS_FUNC)(func))(sys, n1, n2, res);
		return;
	}

	int N = n2 - n1 + 1;

	SYS3_THREAD_VARS<SYS, RES> *thread_vars = new SYS3_THREAD_VARS<SYS, RES>[nThreads];
	int iThread = 0;
	int nstart = n1, nm = 0, nstop;
	int nstep = N / nThreads;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) {nm = N - nstart; nstop = n2;}
		else {nm = nstep; nstop = nstart + nm - 1;}

		thread_vars[iThread].set_pars(iThread, &sys, nstart, nstop);
		thread_vars[iThread].set_func(func);
		nstart += nm;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MRunSysThread_3<SYS, RES>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MRunSysThread_3<SYS, RES>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	res = (RES)0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		res += thread_vars[iThread].res;
		thread_vars[iThread].sys = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};



void gauss_random_force(float f0, float sigma, double *f, float ctheta);

// force on the molecular mass center
void CalNetForce(MMOLECULE &mm, int nc1, int nc2);
// scale force on each cluster
void mm_scale_force(MMOLECULE *mm, int nc1, int nc2);
void scale_force_mmolecule(MMOLECULE &mm);
void scale_force_cluster(BASIC_CLUSTER *pc);
void SCALE_FORCE_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell);


// ****** FUNCTIONS FOR IMPLICIT_SOLVATION
// calcualte SASA of cluster first, then calculate the force = d(U1+U2)/dr, U1 & U2 are the solvation energies of the two related clusters
void Cluster_SASA(BASIC_CLUSTER *pc, CMM_CELL3D<BASIC_CLUSTER> *cell);
void MM_ImpSolvSASA(CMM_CELL3D<BASIC_CLUSTER> *cell, MMOLECULE *mm, int nc1, int nc2);
void ClusterImpSolvForce(BASIC_CLUSTER *pc);
void MM_ImpSolvForce(MMOLECULE *mm, int nc1, int nc2);

// calcualte the force = d(U)/dr, U is the solvation energies of the cluster
// IS THIS A RIGHT WAY TO CALCULATE THE INTERACTION ????
void ClusterImplicitSolvationForce(BASIC_CLUSTER *pc, CMM_CELL3D<BASIC_CLUSTER> *cell); // solvent effect between cluster in the same macro-molecules
// solvent effect between cluster in the periodical CMM CELL
void CMM_ImplicitSolvationForce(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell);
// solvent effect for the clusters in the macro-molecule
void MM_ImplicitSolvationForce(CMM_CELL3D<BASIC_CLUSTER> *cell, MMOLECULE *mm, int nc1, int nc2);

// this function is to calculate the implicit solvation energy on pc, with pc1 at a specific position
// important : this is not pairwise interaction! 
double ClusterImplicitSolvationEnergy(BASIC_CLUSTER *pc, BASIC_CLUSTER *pc1, CMM_CELL3D<BASIC_CLUSTER> *cell);


#if _DISABLE == 0
// interactions calculated with CELL MULTIPOLE METHOD
void FreeCellCMMInteraction_SingleMM(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, double &V, bool &OVERLAP, bool hard_sphere, int n1, int n2); // V in kT
// interactions between clusters in MacroMolecule
void FreeCellClusterInteraction(MMOLECULE &mm, InteractVar& var, InteractRes &res); // U in kT

// interactions calculated with CELL MULTIPOLE METHOD
void FreeCellCMMClusterInteraction(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, InteractVar& var, InteractRes &res); // V in kT
#endif // _DISABLE_


void distribute_MassCenter_real_force(MMOLECULE &mm, VECTOR<6> &fc, int nThreads);
void random_torque_MassCenter(MMOLECULE &mm, float f0);

void Brownian_force(MMOLECULE *mm, int nc1, int nc2, VECTOR<6> &f0);
void SM_Brownian_force(SMOLECULE *sm);
void PM_Brownian_force(PMOLECULE *pm);

extern EF_DPROF LJ12_6;
void construct_LJ12_6_prof(int mode); // 0 -- Standard; 1 -- Alternative -- Amber


void EInteract_SINGLE_MM(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, EInteractVar &evar, InteractRes &res);

class TorsionVar : public InteractVar {
public:
	DIHEDRAL_INTERACT dh;
	TorsionVar() {};
	~TorsionVar() {};
};
void TorsionInteract(MMOLECULE*mm, int nc1, int nc2, TorsionVar &var, InteractRes& ires);

class HingeConstraintVar : public InteractVar {
public:
	HINGE_CONSTRAINT_FORCE hc;
	HingeConstraintVar() {};
	~HingeConstraintVar() {};
};
bool HingeConstraintForce(MMOLECULE &mm, HingeConstraintVar &var, InteractRes& ires);


// external EField
void ExternalEField(BASIC_CLUSTER *bpc, float &Eext, double &Uext);
void CLUSTER_ExternalEField(CLUSTER *cluster, int ncs, float &Eext, double &Uext);
void MM_ExternalEField(MMOLECULE *mm, int nM, float &Eext, double &Uext);
void SM_ExternalEField(SMOLECULE *sm, int nM, float &Eext, double &Uext);
void PM_ExternalEField(PMOLECULE *pm, int nM, float &Eext, double &Uext);

void VMM_ExternalEField(MMOLECULE **mm, int nM, float &Eext, double &Uext);
void VSM_ExternalEField(SMOLECULE **sm, int nM, float &Eext, double &Uext);
void VPM_ExternalEField(PMOLECULE **pm, int nM, float &Eext, double &Uext);

namespace _rand_ {
class GaussRand2 {
public:
	int irand;
	double v[2];

	GaussRand2() {irand = 0; v[0] = 0; v[1] = 0;};
};

double GaussianRand1();
double GaussianRand2();
void test_GaussRand(float width);

template <class T> class RAND_ARRAY {
public:
	ARRAY<T> a, buf;
	int nw;

	void init_elements(int n) {
		if (a.n < n) a.SetArray(n);
		if (buf.n < n) buf.SetArray(n);
		nw = n;
	};
	void random_order() {
		int i, ipt;
		T t;
		int n = nw;
		while (n > 1) {
			ipt = int(ranf() * n);
			if (ipt >= n) ipt = n - 1;
			a.m[nw - n] = buf.m[ipt];
			t = buf.m[ipt]; buf.m[ipt] = buf.m[n - 1]; buf.m[n - 1] = t;
			n -= 1;
		}
		a.m[nw - 1] = buf.m[0];
	};
};

void init_array(RAND_ARRAY<int> &a);

void init_random_array(RAND_ARRAY<int> &rand, int n);

} // end of namespace _rand_

// add random force on each atom

void gen_add_random_force(MMOLECULE &mm, int ncs, float f0, float f_sigma, _rand_::RAND_ARRAY<int> &randa);
void gen_add_random_force(MMOLECULE &mm, int ncs, float f0, float f_sigma, float torq, float t_sigma, _rand_::RAND_ARRAY<int> &randa);



// when calculating force on atom with multi-thread technique, short-range force and electrostatic force
// can be calculated in parallel. In this case, forces should be stored in different memory.
// In MATOM, an array is used, ARRAY< SVECTOR<double, 3> > f
#ifndef _NFS
	#define _NFS    2
#endif

#define _IF_SR  0  // force used for short-range
#define _IF_ES  1  // force for electrostatic

