class FMD_THREAD_PAR : public _THREAD_VARS {
public:
	double ksi;
	bool* res;
	bool local_relax;
	float Ek_total;

	int md_loop, nM;

	char msg_show[5120], msg_log[5120];

	void reset() {
		RELEASE_ARRAY(res) nM = 0;
		strcpy(msg_show, "\0"); strcpy(msg_log, "\0");
	};
	void set(int nM, int iThread) {
		reset();
		_THREAD_VARS::set_pars(iThread);
		this->iThread = iThread;
		if (nM <= 0) return;
		this->nM = nM;
		res = new bool[nM];
	};
	void set_loop(int md_loop) {this->md_loop = md_loop;};
	FMD_THREAD_PAR() {res = NULL; ksi = 0; strcpy(msg_show, "\0"); strcpy(msg_log, "\0"); local_relax = false; md_loop = 0;};
	~FMD_THREAD_PAR() {reset();};
};

class PM_FMD_THREAD_PAR : public FMD_THREAD_PAR {
public:
	PMOLECULE* pm;
	LFVERLET_PM_MD *lfpmdpar;
	double *dv;

	PM_FMD_THREAD_PAR() {pm = NULL; lfpmdpar = NULL; dv = NULL;};
	void set(PMOLECULE *pm, LFVERLET_PM_MD* lfpmdpar, int nM, int iThread) {
		this->pm = pm; this->lfpmdpar = lfpmdpar;
		RELEASE_ARRAY(dv);
		if (nM <= 0) return;
		FMD_THREAD_PAR::set(nM, iThread);
		dv = new double[nM];
	};
	~PM_FMD_THREAD_PAR() {FMD_THREAD_PAR::reset(); RELEASE_ARRAY(dv);};
};

class SM_FMD_THREAD_PAR : public FMD_THREAD_PAR {
public:
	SMOLECULE* sm;
	LFVERLET_SM_MD *lfsmdpar;
	double *dv, *dw;

	SM_FMD_THREAD_PAR() {sm = NULL; lfsmdpar = NULL; dv = NULL; dw = NULL;};
	void set(SMOLECULE *sm, LFVERLET_SM_MD* lfsmdpar, int nSM, int iThread) {
		this->sm = sm; this->lfsmdpar = lfsmdpar;
		RELEASE_ARRAY(dv) RELEASE_ARRAY(dw)
		if (nSM <= 0) return;
		FMD_THREAD_PAR::set(nSM, iThread);
		dv = new double[nSM]; dw = new double[nSM];
	};
	~SM_FMD_THREAD_PAR() {FMD_THREAD_PAR::reset(); RELEASE_ARRAY(dv) RELEASE_ARRAY(dw)};
};

class MM_FMD_THREAD_PAR : public FMD_THREAD_PAR {
public:
	MMOLECULE *mm;
	LFVERLET_MM_MD *lfmdpar;

	VECTOR<6> *vect1, *vect2;
	double* dv, *dw;

	MM_FMD_THREAD_PAR() {mm = NULL; lfmdpar = NULL; vect1 = NULL; vect2 = NULL; dv = NULL; dw = NULL;};
	void set(MMOLECULE *mm, LFVERLET_MM_MD* lfmdpar, bool local_relax, int nMM, int iThread) {
		this->mm = mm; this->lfmdpar = lfmdpar; this->local_relax = local_relax;
		if (nMM <= 0) return;
		FMD_THREAD_PAR::set(nMM, iThread);
	};
	void reset_buffer_memory() {
		RELEASE_ARRAY(vect1); RELEASE_ARRAY(vect2); RELEASE_ARRAY(dv); RELEASE_ARRAY(dw);
	};
	void reset() {FMD_THREAD_PAR::reset(); reset_buffer_memory();};
	void set_buffer_memory(int nset) {
		reset_buffer_memory();
		this->vect1 = new VECTOR<6>[nset]; this->vect2 = new VECTOR<6>[nset];
		this->dv = new double[nset]; this->dw = new double[nset];
	};
	~MM_FMD_THREAD_PAR() {reset();};
};

#if _SYS_ == _WINDOWS_SYS_
int PM_SpeedVerlet_FMD_thread_func(LPVOID *par);
int SM_SpeedVerlet_FMD_thread_func(LPVOID *par);
int MM_SpeedVerlet_FMD_thread_func(LPVOID *par);
#elif _SYS_ == _LINUX_SYS_
void *PM_SpeedVerlet_FMD_thread_func(LPVOID *par);
void *SM_SpeedVerlet_FMD_thread_func(void *par);
void *MM_SpeedVerlet_FMD_thread_func(void *par);
#endif

void FMD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm);
