// simple solid molecule -- has only one cluster
// parameter for Dynamic
class SMOLECULE {
public:
	BASIC_CLUSTER *c; // cluster's Op is located at mass center, VERY IMPORTANT!!!
	VECTOR3 xaxis, zaxis; // local cordinates of the solid molecule, with origin at mass center

	VECTOR3 r; // mass center position in experimental cordinates, also this is Op of the cluster
	MATRIX<3> I, invI; // spatial moment and its inverse

	MOL_BASE_MOVEMENT mbm;

	float Ek, Ekt, Ektx, Ekty, Ektz;
	float Ep; // the total interaction energy with other molecules in whole system

	short mIndx; // the index of molecule in the whole system, each molecule has an unique number, starting from 0
	short mType; // the type of the molecule, 0 -- MM, 1 -- SM, 2 -- PM
	short mTypeIndx; // uid of the molecule in the CELL, because molecule is not uniquely defined
	char mol[20]; 

	bool bFixed;

	SMOLECULE() {
		xaxis.v[0] = 1; zaxis.v[2] = 1; Ek = 0; Ep = 0; 
		c = NULL; bFixed = false;
		mIndx = 0; mType = -1; mTypeIndx = -1; strcpy(mol, "\n");
	};
	void reset() {
		if (c != NULL) {c->reset(); delete c; c = NULL;}
	};
	void setup_kin_dyn() {
		if (c == NULL) return; c->setup_kin_dyn();
	};
	void setup_emp() {
		if (c == NULL) return; c->setup_emp();
	};
	~SMOLECULE() {reset();};
	void check_atom_cordinate(); // set cluster's mass center at mass center
	void calInertiaMoment();
	void shiftMol(VECTOR3 &ds) {
		if (c == NULL) return;
		int i = 0, na = 0;
		for (i = 0; i < 3; i++) {this->r.v[i] += ds.v[i];}
		for (na = 0; na < c->nAtoms; na++) {
			for (i = 0; i < 3; i++) c->atom[na].r.v[i] += ds.v[i];
		}
	};
	void set_base_move2(VECTOR3& r2, VECTOR3& zaxis2, VECTOR3& xaxis2) {
		mbm.bRot = false;
		VECT3(r, r2, mbm.dr)
		V32V3(zaxis2, mbm.zaxis) V32V3(xaxis2, mbm.xaxis)
	};
	void set_base_move(VECTOR3& dr, VECTOR3& drot) {
		mbm.bRot = true;
		V32V3(dr, mbm.dr) V32V3(drot, mbm.drot)
	};
	void calTotalForce() {
		double *fc = c->dyn->fc.v, *fh = c->dyn->f_Hoover.v, *fb = c->dyn->f_Brown.v, *f = c->dyn->f.v;
		for (int i = 0; i < 6; i++) f[i] = fc[i] + fh[i] + fb[i];
	};
	void updateSpatialMomentum() {
		MpV3(I, c->bkm->V0, c->bkm->P0)
		c->bkm->P0.v[3] = c->M * c->bkm->V0.v[3];
		c->bkm->P0.v[4] = c->M * c->bkm->V0.v[4];
		c->bkm->P0.v[5] = c->M * c->bkm->V0.v[5];
	};
	void updateRotationVelocity() {
		MpV3(invI, c->bkm->P0, c->bkm->V0)  // P ==> V0
		// do not need to update translation velocity
	};
};

bool cp(SMOLECULE *d, SMOLECULE *s, bool setup = true);

void rotate_SM(SMOLECULE& m, R &r);
void rotate_SM(SMOLECULE &m, VECTOR3 &axis, double omega, bool arc);
void rotate_SM(SMOLECULE &m, VECTOR3 &drot, bool arc);

void SMolMove(SMOLECULE &sm);

double calKineticEnergy(SMOLECULE &m);
void InitVelocity(SMOLECULE &m);
double ScaleKineticEnergy(SMOLECULE &m, double Ek0);

void calAccel(SMOLECULE &m, float ksi);
void calAccel_f(SMOLECULE &m);

bool atom_indx(SMOLECULE &sm, MATOM* patom, int &na, int &indx);

class PMOLECULE { // point molecule / ion
public:
	BASIC_CLUSTER *c; // cluster's mass center is located at mass center
	VECTOR3 r, v, alpha; 
	float q; // charge
	float Ek, Ekt, Ektx, Ekty, Ektz;

	short mIndx; // the index of molecule in the whole system, each molecule has an unique number, starting from 0
	short mType; // the type of the molecule, 0 -- MM, 1 -- SM, 2 -- PM
	short mTypeIndx; // uid of the molecule in the CELL, because molecule is not uniquely defined
	char mol[20];

	bool bFixed;

	PMOLECULE() {
		c = NULL; q = 0; Ek = 0; bFixed = false;
		mIndx = 0; mType = -1; mTypeIndx = -1; strcpy(mol, "\n");
	};
	inline void shiftMol(VECTOR3& dr) {
		V3plusV3(r, dr, r) 
		if (c != NULL) c->shift(dr);
	};
	inline void check_coordinate() {
		if (c == NULL) return;
		V3zero(c->atom[0].r0)
		c->atom[0].r = r; c->atom[0].rg = r;
		c->vp = &(r); c->Op = &(r); c->vcm = &r;
	};
	void setup_dyn() {
		if (c->bkm == NULL) c->bkm = new BASIC_KINEMATIC();
		if (c->dyn == NULL) c->dyn = new BASIC_DYN();
		int i;
		if (::bMTForce) {
			for (i = 0; i < c->nAtoms; i++) c->atom[i].f.set_array(_NFS);
		}
	};
	void reset() {if (c != NULL) delete c; c = NULL;};
	~PMOLECULE() {if (c != NULL) delete c; c = NULL;};
	void calTotalForce() {
		double *fc = c->dyn->fc.v, *fh = c->dyn->f_Hoover.v, *fb = c->dyn->f_Brown.v, *f = c->dyn->f.v;
		for (int i = 0; i < 6; i++) f[i] = fc[i] + fh[i] + fb[i];
	};
	void updateSpatialMomentum() {
		c->bkm->P0.v[3] = c->M * v.v[0];
		c->bkm->P0.v[4] = c->M * v.v[1];
		c->bkm->P0.v[5] = c->M * v.v[2];
	};
};

bool cp(PMOLECULE *d, PMOLECULE *s, bool setup = true);

double calKineticEnergy(PMOLECULE &m);
void InitVelocity(PMOLECULE &m);
double ScaleKineticEnergy(PMOLECULE &m, double Ek0);
void calAccel(PMOLECULE &m, float ksi);
void calAccel_f(PMOLECULE &m);

void CalSMElectroMultipole(SMOLECULE *sm, int nSM, bool bEwaldSum);

void polarize_SM(SMOLECULE *sm, int n);
void polarize_PM(PMOLECULE *pm, int n);

void polarize_VSM(SMOLECULE **sm, int n);
void polarize_VPM(PMOLECULE **pm, int n);

template <class T> class M_THREAD_VARS : public _THREAD_VARS {
public:
	T *m;
	int n;
	void* func;

	M_THREAD_VARS() {
		m = NULL; n = 0; func = NULL;
	};
	void set_pars(int iThread, T* m, int n) {
		this->m = m; this->n = n;
		_THREAD_VARS::set_pars(iThread);
	};
	void set_func(void* func) {this->func = func;};
	~M_THREAD_VARS() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class T> int MThread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class T> void* MThread(void *void_par) { // thread-function
#endif
	typedef void (*M_PVOID)(T*, int);
	M_THREAD_VARS<T> *par = (M_THREAD_VARS<T>*)void_par;
	((M_PVOID)(par->func))(par->m, par->n);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class T> void MOperate(void* func, T *m, int n, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*M_PVOID)(T*, int);
		((M_PVOID)(func))(m, n);
		return;
	}

	int nstart = 0, nstep = n / nThreads; if (nstep < 1) {nstep = 1; nThreads = n;}

	M_THREAD_VARS<T> *thread_vars = new M_THREAD_VARS<T>[nThreads];
	int iThread = 0, nmol = 0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) nmol = n - nstart;
		else nmol = nstep;

		thread_vars[iThread].set_pars(iThread, m + nstart, nmol);
		thread_vars[iThread].set_func(func);
		nstart += nmol;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MThread<T>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MThread<T>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].m = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};

template <class T, class TVAR> class M_THREAD_VARS1 : public _THREAD_VARS {
public:
	T *m;
	int n;
	TVAR var;
	void* func;

	M_THREAD_VARS1() {
		m = NULL; n = 0; func = NULL;
	};
	void set_pars(int iThread, T* m, int n, TVAR& var) {
		this->m = m; this->n = n; this->var = var;
		_THREAD_VARS::set_pars(iThread);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~M_THREAD_VARS1() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class T, class TVAR> int MThread1(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class T, class TVAR> void* MThread1(void *void_par) { // thread-function
#endif
	typedef void (*M_VOIDV)(T*, int, TVAR&);
	M_THREAD_VARS1<T, TVAR> *par = (M_THREAD_VARS1<T, TVAR>*)void_par;
	((M_VOIDV)(par->func))(par->m, par->n, par->var);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class T, class TVAR> void MOperate1(void *func, T *m, int n, TVAR &var, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*M_VOIDV)(T*, int, TVAR&);
		((M_VOIDV)(func))(m, n, var);
		return;
	}

	int nstart = 0, nstep = n / nThreads; if (nstep < 1) {nstep = 1; nThreads = n;}

	M_THREAD_VARS1<T, TVAR> *thread_vars = new M_THREAD_VARS1<T, TVAR>[nThreads];
	int iThread = 0, nmol = 0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) nmol = n - nstart;
		else nmol = nstep;

		thread_vars[iThread].set_pars(iThread, m + nstart, nmol, var);
		thread_vars[iThread].set_func(func);
		nstart += nmol;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MThread1<T, TVAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MThread1<T, TVAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].m = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};


//*****************************************************************************************
//**  MULTI-THREAD FUNCTIONS FOR SIMPLE MOLECULE OPERATION WITHOUT ADDITIONAL PARAMETERS **
//**  RETURN A RESULT                                                                    **
//*****************************************************************************************

template <class T, class TVAR> class M_THREAD_VARS2 : public _THREAD_VARS {
public:
	T *m;
	int n;
	TVAR res;
	void* func;

	M_THREAD_VARS2() {
		m = NULL; n = 0; func = NULL;
	};
	void set_pars(int iThread, T* m, int n) {
		this->m = m; this->n = n;
		_THREAD_VARS::set_pars(iThread);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~M_THREAD_VARS2() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class T, class TVAR> int MThread2(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class T, class TVAR> void* MThread2(void *void_par) { // thread-function
#endif
	typedef void (*M_VOIDV)(T*, int, TVAR&);
	M_THREAD_VARS2<T, TVAR> *par = (M_THREAD_VARS2<T, TVAR>*)void_par;
	((M_VOIDV)(par->func))(par->m, par->n, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class T, class TVAR> void MOperate2(void *func, T *m, int n, int nThreads, TVAR &res, bool init_res = false, void *res_operator = NULL) {
	if (nThreads <= 1) {
		typedef void (*M_VOIDV)(T*, int, TVAR&);
		((M_VOIDV)(func))(m, n, res);
		return;
	}

	res = (TVAR)0;
	int nstart = 0, nstep = n / nThreads; if (nstep < 1) {nstep = 1; nThreads = n;}
	M_THREAD_VARS2<T, TVAR> *thread_vars = new M_THREAD_VARS2<T, TVAR>[nThreads];
	int iThread = 0, nmol = 0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) nmol = n - nstart;
		else nmol = nstep;

		thread_vars[iThread].set_pars(iThread, m + nstart, nmol);
		thread_vars[iThread].set_func(func);
		if (init_res) thread_vars[iThread].res = res;
		nstart += nmol;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MThread2<T, TVAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MThread2<T, TVAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	typedef void (*RES_OPERATOR)(TVAR&, TVAR&, TVAR&);
	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].m = NULL;
		if (res_operator == NULL) res += thread_vars[iThread].res;
		else {
			((RES_OPERATOR)res_operator)(res, thread_vars[iThread].res, res);
		}
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};

//*****************************************************
//**  MULTI-THREAD FUNCTIONS FOR MOLECULE OPERATION  **
//**  WITH A PARAMETERS, RETURN A RESULT             **
//*****************************************************
template <class T, class VVAR, class RESVAR> class M_THREAD_VARS3 : public _THREAD_VARS {
public:
	T *m;
	int n;
	VVAR v;
	RESVAR res;
	void* func;

	M_THREAD_VARS3() {
		m = NULL; n = 0; func = NULL;
	};
	void set_pars(int iThread, T* m, int n, VVAR &v) {
		this->m = m; this->n = n; this->v = v;
		_THREAD_VARS::set_pars(iThread);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~M_THREAD_VARS3() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class T, class VVAR, class RESVAR> int MThread3(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class T, class VVAR, class RESVAR> void *MThread3(void *void_par) { // thread-function
#endif
	M_THREAD_VARS3<T, VVAR, RESVAR> *par = (M_THREAD_VARS3<T, VVAR, RESVAR>*)void_par;
	typedef void (*M_VOIDV)(T*, int, VVAR&, RESVAR&);
	((M_VOIDV)(par->func))(par->m, par->n, par->v, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class T, class VVAR, class RESVAR> void MOperate3(void *func, T *m, int n, VVAR &var, int nThreads, RESVAR &res, bool init_res = false, void *res_operator = NULL) {
	if (nThreads <= 1) {
		typedef void (*M_VOIDV)(T*, int, VVAR&, RESVAR&);
		((M_VOIDV)(func))(m, n, var, res);
		return;
	}

	res = (RESVAR)0;
	int nstart = 0, nstep = n / nThreads; if (nstep < 1) {nstep = 1; nThreads = n;}
	M_THREAD_VARS3<T, VVAR, RESVAR> *thread_vars = new M_THREAD_VARS3<T, VVAR, RESVAR>[nThreads];
	int iThread = 0, nmol = 0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) nmol = n - nstart;
		else nmol = nstep;

		thread_vars[iThread].set_pars(iThread, m + nstart, nmol, var);
		thread_vars[iThread].set_func(func);
		if (init_res) thread_vars[iThread].res = res;
		nstart += nmol;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MThread3<T, VVAR, RESVAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MThread3<T, VVAR, RESVAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	typedef void (*RES_OPERATOR)(RESVAR&, RESVAR&, RESVAR&);
	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].m = NULL;
		if (res_operator == NULL) res += thread_vars[iThread].res;
		else {
			((RES_OPERATOR)res_operator)(res, thread_vars[iThread].res, res);
		}
	}
	delete[] thread_vars; thread_vars = NULL;
	return;
};
