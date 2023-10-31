int MM_mpi_process(int imol);
int SM_mpi_process(int imol);
void mpi_multi_mol_md_proc(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, int md_mode);


class MM_EF_VAR {
public:
	bool checkCMM, bCheckClusterInCMM;
	MM_EF_VAR() {checkCMM = false; bCheckClusterInCMM = false;};
	void operator =(MM_EF_VAR& v) {  // important funciton
		this->checkCMM = v.checkCMM;
		this->bCheckClusterInCMM = v.bCheckClusterInCMM;
	};
};

class MM_EF_RES {
public:
	double Ep;
	bool OVERLAP;
	MM_EF_RES() {Ep = 0;};
	void initialize() {Ep = 0;};  // important funciton
	void accumulate(MM_EF_RES &r) {  // important funciton
		this->Ep += r.Ep;
		if (OVERLAP) this->OVERLAP = true;
	};
};



class MM_IMPSOLV_VAR {
public:
	bool bCheckClusterInCMM;
	MM_IMPSOLV_VAR() {bCheckClusterInCMM = false;};
	void operator =(MM_IMPSOLV_VAR& v) {  // important funciton
		this->bCheckClusterInCMM = v.bCheckClusterInCMM;
	};
};

class MM_IMPSOLV_RES {
public:
	int res;
	MM_IMPSOLV_RES() {res = 0;};
	void initialize() {res = 0;};  // important funciton
	void accumulate(MM_IMPSOLV_RES &r) {  // important funciton
		this->res += r.res;
	};
};

template <class MACROMOL> class MACROMOL_RUN_BKGD : public _MACROMOL_THREAD_VARS<MACROMOL> {
public:
	void *func;
	MACROMOL_RUN_BKGD() {func = NULL;};
	~MACROMOL_RUN_BKGD() {func = NULL;};
};

#if _SYS_ == _WINDOWS_SYS_
template <class MACROMOL> int MacroMol_RunBkgd_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MACROMOL> void* MacroMol_RunBkgd_Thread(void *void_par) { // thread-function
#endif
	MACROMOL_RUN_BKGD<MACROMOL> *par = (MACROMOL_RUN_BKGD<MACROMOL>*)void_par;
	typedef void (*MACROMOL_OPERATOR)(MACROMOL*);
	((MACROMOL_OPERATOR)(par->func))(par->mm);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class MACROMOL> void MacroMol_Run_bkgd(MACROMOL_RUN_BKGD<MACROMOL> &mpar) {
	if (mpar.mm == NULL) return;
	mpar._THREAD_VARS::set_pars(1);
#if _SYS_ == _WINDOWS_SYS_
	ResetEvent(mpar.hMMThread);
	mpar.thread = AfxBeginThread((AFX_THREADPROC)MacroMol_RunBkgd_Thread<MACROMOL>, (LPVOID)(&mpar), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
#elif _SYS_ == _LINUX_SYS_
	pthread_create(&(mpar.thread), NULL, &MacroMol_RunBkgd_Thread<MACROMOL>, (void *)(&mpar));
#endif
};





template <class MACROMOL, class atom> class CELL_MACROMOL_RUN_BKGD : public _MACROMOL_THREAD_VARS<MACROMOL> {
public:
	CMM_CELL3D<atom> *cell;
	void *func;
	CELL_MACROMOL_RUN_BKGD() {cell = NULL; func = NULL;};
	~CELL_MACROMOL_RUN_BKGD() {cell = NULL; func = NULL;};
};

#if _SYS_ == _WINDOWS_SYS_
template <class MACROMOL, class atom> int CellMacroMol_RunBkgd_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MACROMOL, class atom> void* CellMacroMol_RunBkgd_Thread(void *void_par) { // thread-function
#endif
	CELL_MACROMOL_RUN_BKGD<MACROMOL, atom> *par = (CELL_MACROMOL_RUN_BKGD<MACROMOL, atom>*)void_par;
	typedef void (*CELL_MACROMOL_OPERATOR)(CMM_CELL3D<atom>*, MACROMOL*);
	((CELL_MACROMOL_OPERATOR)(par->func))(par->cell, par->mm);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class MACROMOL, class atom> void CellMacroMol_Run_bkgd(CELL_MACROMOL_RUN_BKGD<MACROMOL, atom> &mpar) {
	if (mpar.cell == NULL || mpar.func == NULL || mpar.mm == NULL) return;
	mpar._THREAD_VARS::set_pars(1);
#if _SYS_ == _WINDOWS_SYS_
	ResetEvent(mpar.hMMThread);
	mpar.thread = AfxBeginThread((AFX_THREADPROC)CellMacroMol_RunBkgd_Thread<MACROMOL, atom>, (LPVOID)(&mpar), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
#elif _SYS_ == _LINUX_SYS_
	pthread_create(&(mpar.thread), NULL, &CellMacroMol_RunBkgd_Thread<MACROMOL, atom>, (void *)(&mpar));
#endif
};

