//#ifndef MAX_THREADS
//#define MAX_THREADS 4
//#endif

class _THREAD_VARS {
public:
	int iThread;

#if _SYS_ == _WINDOWS_SYS_
	HANDLE hMMThread;
	CWinThread *thread;
#elif _SYS_ == _LINUX_SYS_
	pthread_t thread;
#endif

	_THREAD_VARS() {
		iThread = 0;
	#if _SYS_ == _WINDOWS_SYS_
		hMMThread = NULL; thread = NULL;
	#elif _SYS_ == _LINUX_SYS_
		thread = 0;
	#endif
	};
	void set_pars(int iThread) {
		this->iThread = iThread;
	#if _SYS_ == _WINDOWS_SYS_
		if (hMMThread == NULL) hMMThread = CreateEvent(NULL, FALSE, FALSE, NULL);
		else ResetEvent(hMMThread);
		thread = NULL;
	#elif _SYS_ == _LINUX_SYS_
		thread = 0;
	#endif
	};
	~_THREAD_VARS() {
	#if _SYS_ == _WINDOWS_SYS_
		if (hMMThread != NULL) {CloseHandle(hMMThread); hMMThread = NULL;}
		thread = NULL;
	#elif _SYS_ == _LINUX_SYS_
		thread = 0;
	#endif
	};
	void WaitUntilThreadOver() {
	#if _SYS_ == _WINDOWS_SYS_
		if (hMMThread == NULL) return;
		WaitForSingleObject(hMMThread, INFINITE);
	#elif _SYS_ == _LINUX_SYS_
		if (thread == 0) return;
		pthread_join(thread, NULL); // wait the thread over
	#endif
	};
};


//**********************************************************************************
//**            a general class to run a function in a new thread                 **
//**********************************************************************************

template <class VAR> class AsynRun;

#if _SYS_ == _WINDOWS_SYS_
template <class VAR> int asynThread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class VAR> void* asynThread(void *void_par) { // thread-function
#endif
	AsynRun<VAR> *par = (AsynRun<VAR>*)void_par;
	typedef void (*FUNC)(void*);
	((FUNC)(par->func))(par->v);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};


template <class VAR> class AsynRun : public _THREAD_VARS {
public:
	VAR *v;
	void *func;

	AsynRun() {v = NULL; func = NULL;};
	bool set_par(VAR* v) {
		if (this->thread != NULL) return false; // the background thread is busy
		this->v = v;
		_THREAD_VARS::set_pars(0);
		return true;
	};
	bool is_busy() {
		return (this->thread == NULL ? false : true);
	};

	bool start_func(void* func) {
		if (this->thread == NULL) {
			this->func = func;
#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(this->hMMThread);
			this->thread = AfxBeginThread((AFX_THREADPROC)asynThread<VAR>, (LPVOID)(this), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(this->thread), NULL, &asynThread<VAR>, (void*)(this));
#endif
			return true;
		}
		else return false;
	};

	bool wait() {
#if _SYS_ == _WINDOWS_SYS_ 
		if (this->hMMThread == NULL) return false;
		WaitForSingleObject(this->hMMThread, INFINITE);
#elif _SYS_ == _LINUX_SYS_
		if (this->thread == NULL) return false;
		pthread_join(this->thread, NULL); // wait the thread over
#endif
		this->thread = NULL;
		return true;
	};

	~AsynRun() { // if the background thread is still running, the thread will be lost
#if _SYS_ == _WINDOWS_SYS_ 
		CloseHandle(this->hMMThread); this->hMMThread = NULL;
#elif _SYS_ == _LINUX_SYS_
		this->thread = NULL;
#endif
	};
};


//**********************************************************************************
//**  MULTI-THREAD FUNCTIONS FOR MACROMOL OPERATION WITHOUT ADDITIONAL PARAMETERS **
//**********************************************************************************

template <class CM, class T> class CM3OPERATE_THREAD_VARS : public _THREAD_VARS {
public:
	CM *m1, *m2, *mres;
	int n1, n2;
	void* func;

	CM3OPERATE_THREAD_VARS() {
		m1 = NULL; m2 = NULL; mres = NULL;
	};
	void set_pars(int iThread, CM *m1, CM *m2, CM *mres, int n1, int n2) {
		this->m1 = m1; this->m2 = m2; this->mres = mres; this->n1 = n1; this->n2 = n2;
		_THREAD_VARS::set_pars(iThread);
	};
	void set_func(void* func) {this->func = (void*)func;};
	~CM3OPERATE_THREAD_VARS() {};
};


#if _SYS_ == _WINDOWS_SYS_
template <class CM, class T> int CM3Operate_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class CM, class T> void* CM3Operate_Thread(void *void_par) { // thread-function
#endif
	CM3OPERATE_THREAD_VARS<CM, T> *par = (CM3OPERATE_THREAD_VARS<CM, T>*)void_par;
	typedef void (*CM3_OPERATOR)(CM*, CM*, CM*, int, int);
	((CM3_OPERATOR)(par->func))(par->m1, par->m2, par->mres, par->n1, par->n2);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class CM, class T> void CM3Operate(void* func, CM &m1, CM &m2, CM &mres, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*CM3_OPERATOR)(CM*, CM*, CM*, int, int);
		((CM3_OPERATOR)(func))(&m1, &m2, &mres, 0, m1.n);
		return;
	}

	CM3OPERATE_THREAD_VARS<CM, T> *thread_vars = new CM3OPERATE_THREAD_VARS<CM, T>[nThreads];
	int iThread = 0;
	int n1 = 0, n2 = 0, nstep = m1.n / nThreads;
	for (iThread = 0; iThread < nThreads; iThread++) {
		nstep = (m1.n - n1) / (nThreads - iThread);
		if (nstep * (nThreads - iThread) < (m1.n - n1)) nstep += 1;
		if (iThread == nThreads - 1) n2 = m1.n - 1;
		else n2 = n1 + nstep - 1;

		if (n2 >= m1.n) n2 = m1.n - 1;

		thread_vars[iThread].set_pars(iThread, &m1, &m2, &mres, n1, n2);
		thread_vars[iThread].set_func(func);
		n1 += n2 + 1;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)CM3Operate_Thread<CM, T>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &CM3Operate_Thread<CM, T>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		//thread_vars[iThread].mm = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
};



#ifndef __GPU__
template <class T> class ARRAY {
public:
	int n;
	T *m;

	ARRAY() {
		n = 0; m = NULL; 
	};

	void release() {
		if (m != NULL) {delete[] m; m = NULL;} n = 0;
	};

	bool SetArray(int n) {
		if (this->n == n) return true;
		release(); if (n <= 0) return true;
		m = new T[n]; this->n = n; return true;
	};

	bool set_array(int n) {return SetArray(n);};
	bool exist(T v) {
		for (int i = 0; i < n; i++) {
			if (m[i] == v) return true;
		}
		return false;
	};
	int indx(T v) {
		for (int i = 0; i < n; i++) {
			if (m[i] == v) return i;
		}
		return -1;
	};
	void release_pointer() {
		m = NULL; n = 0;
	};
	void expand(int nIncrease) {
		int n1 = n;
		int n2 = n + nIncrease;
		T *bf = m;
		release_pointer();
		set_array(n2);
		if (bf == NULL) return;
		//for (i = 0; i < n1; i++) m[i] = bf[i];
		memcpy(m, bf, n1 * sizeof(T));
		delete[] bf;
	};
	void remove_one(int indx) {
		if (n == 0 || indx < 0) indx += n;
		if (indx < 0 || indx >= n) return;
		T* bf = m;
		int n2 = n - 1;
		if (n2 <= 0) {release(); return;}
		else {
			release_pointer();
			set_array(n2);
			if (indx > 0) memcpy(m, bf, indx * sizeof(T));
			if (indx < n) memcpy(m + indx, bf + indx + 1, (n - indx - 1) * sizeof(T));
			delete bf; bf = NULL;
			return;
		}
	};
	~ARRAY() {release();};
};

// T is the type without pointer inside
template <class T> class ARRAY_FLEX {
public:
	int n, nm; // n is the dimension used, nm is the memory size
	T *m;

	ARRAY_FLEX() {
		n = 0; m = NULL; nm = 0;
	};
	void release(bool release_memory = false) {
		if (release_memory) {
			if (m != NULL) {delete[] m; m = NULL;} n = 0; nm = 0;
		}
		else n = 0;
	};

	void SetArray(int nm) {
		if (this->nm == nm) return;
		release(true); if (nm <= 0) return;
		m = new T[nm]; this->nm = nm; n = 0;
	};

	void set_array(int nm) {SetArray(nm);};
	void release_pointer() {
		m = NULL; nm = 0; n = 0;
	};

	void expand(int nIncrease) {
		int n1 = nm;
		int n2 = nm + nIncrease;
		T *bf = m;
		release_pointer();
		set_array(n2);
		if (bf == NULL) return;
		memcpy(m, bf, n1 * sizeof(T));
		delete[] bf;
	};
	void attach(T v) {
		if (n == nm) expand(5);
		memcpy(m + n, &v, sizeof(T)); 
		n++;
	};
	~ARRAY_FLEX() {release(true);};
};

#endif



#ifdef __GPU__
template <class T> class ARRAY {
public:
	int n;
	T *m;

	bool bGPU_MapMem;
	T *m_dev; // device pointer to the mapped memory, looking from GPU kernel

	ARRAY() {
		n = 0; m = NULL; 
		bGPU_MapMem = false;
		m_dev = NULL;
	};

	ARRAY(bool bGPU_MapMem) {
		n = 0; m = NULL; 
		this->bGPU_MapMem = bGPU_MapMem;
		m_dev = NULL;
	};

	void release() {
		if (this->bGPU_MapMem) {
			if (m_dev != NULL) {
				//cudaFree(m_dev); m_dev = NULL; // can not free this from device, although it is mapped to device
				m_dev = NULL;
			}
			if (m != NULL) {
				cudaFreeHost(m); m = NULL;
			}
			n = 0;
		}
		else if (m != NULL) {delete[] m; m = NULL;} 
		n = 0;
	};

	void use_GPU_MapMem(bool bGPU_MapMem = false) {
		int n = this->n; T* bf = this->m;
		bool bOldMode = this->bGPU_MapMem;
		if (this->bGPU_MapMem != bGPU_MapMem && m != NULL) {
			release_pointer();
			this->bGPU_MapMem = bGPU_MapMem;
			SetArray(n);
			memcpy(this->m, bf, n * sizeof(T));
			if (bOldMode) cudaFreeHost(bf);
			else delete[] bf; 
			bf = NULL;
		}
		else this->bGPU_MapMem = bGPU_MapMem;
	};
	bool SetArray(int n) {
		extern void show_msg(char* msg, bool endline = true);
		if (this->n == n) return true;
		release(); if (n <= 0) return true;

		if (this->bGPU_MapMem) {
			if (cudaHostAlloc((void**)(&m), n * sizeof(T), cudaHostAllocMapped) != cudaSuccess) {
				show_msg("failure to create GPU-mapped host-memory"); m = NULL; n = 0; return false;
			}
			if (cudaHostGetDevicePointer((void**)&m_dev, m, 0) != cudaSuccess) {
				show_msg("failure to get device pointer to mapped memory"); release(); return false;
			}
			this->n = n;
			return true;
		}
		else {
			m = new T[n]; this->n = n; return true;
		}
	};

	bool set_array(int n) {return SetArray(n);};
	bool exist(T v) {
		for (int i = 0; i < n; i++) {
			if (m[i] == v) return true;
		}
		return false;
	};
	int indx(T v) {
		for (int i = 0; i < n; i++) {
			if (m[i] == v) return i;
		}
		return -1;
	};
	void release_pointer() {
		m = NULL; n = 0;
		m_dev = NULL;
	};
	void expand(int nIncrease) {
		int n1 = n;
		int n2 = n + nIncrease;
		T *bf = m;
		release_pointer();
		set_array(n2);
		if (bf == NULL) return;
		//for (i = 0; i < n1; i++) m[i] = bf[i];
		memcpy(m, bf, n1 * sizeof(T));
		if (this->bGPU_MapMem) cudaFreeHost(bf);
		else delete[] bf;
	};
	void remove_one(int indx) {
		if (n == 0 || indx < 0) indx += n;
		if (indx < 0 || indx >= n) return;
		T* bf = m;
		int n2 = n - 1;
		if (n2 <= 0) {release(); return;}
		else {
			release_pointer();
			set_array(n2);
			if (indx > 0) memcpy(m, bf, indx * sizeof(T));
			if (indx < n) memcpy(m + indx, bf + indx + 1, (n - indx - 1) * sizeof(T));
			delete bf; bf = NULL;
			return;
		}
	};
	~ARRAY() {release();};
};

// T is the type without pointer inside
template <class T> class ARRAY_FLEX {
public:
	int n, nm; // n is the dimension used, nm is the memory size
	T *m;
	bool bGPU_MapMem;
	T* m_dev; // device pointer to the mapped memory, looking from GPU kernel

	ARRAY_FLEX() {
		n = 0; m = NULL; nm = 0;
		bGPU_MapMem = false;
		m_dev = NULL;
	};
	void release(bool release_memory = false) {
		if (this->bGPU_MapMem) {
			if (release_memory) {
				if (m_dev != NULL) {
					//cudaFree(m_dev); m_dev = NULL; // can not free this from device, although it is mapped to device
					m_dev = NULL;
				}
				if (m != NULL) {
					cudaFreeHost(m); m = NULL; 
				}
				n = 0; nm = 0;
			}
			else n = 0;
			return;
		}
		else if (release_memory) {
			if (m != NULL) {delete[] m; m = NULL;} n = 0; nm = 0;
		}
		else n = 0;
	};

	void use_GPU_MapMem(bool bGPU_MapMem = false) {
		int nm = this->nm; T* bf = this->m;
		bool bOldMode = this->bGPU_MapMem;
		if (this->bGPU_MapMem != bGPU_MapMem && m != NULL) {
			release_pointer();
			this->bGPU_MapMem = bGPU_MapMem;
			SetArray(nm);
			memcpy(this->m, bf, nm * sizeof(T));
			if (bOldMode) cudaFreeHost(bf);
			else delete[] bf; 
			bf = NULL;
		}
		else this->bGPU_MapMem = bGPU_MapMem;
	};
	bool SetArray(int nm) {
		extern void show_msg(char* msg, bool endline = true);
		if (this->nm == nm) return true;
		release(true); if (nm <= 0) return true;
		if (bGPU_MapMem) {
			if (cudaHostAlloc((void**)(&m), n * sizeof(T), cudaHostAllocMapped) != cudaSuccess) {
				show_msg("failure to create GPU-mapped host-memory"); m = NULL; nm = 0; return false;
			}
			if (cudaHostGetDevicePointer((void**)&m_dev, m, 0) != cudaSuccess) {
				show_msg("failure to get device pointer to mapped memory"); release(true); return false;
			}
			this->nm = nm;
			return true;
		}
		else {
			m = new T[nm]; this->nm = nm; n = 0;
		}
		return true;
	};

	void set_array(int nm) {SetArray(nm);};
	void release_pointer() {
		m = NULL; nm = 0; n = 0;
		m_dev = NULL;
	};

	void expand(int nIncrease) {
		int n1 = nm;
		int n2 = nm + nIncrease;
		T *bf = m;
		release_pointer();
		set_array(n2);
		if (bf == NULL) return;
		memcpy(m, bf, n1 * sizeof(T));
		if (this->bGPU_MapMem) cudaFreeHost(bf);
		else delete[] bf;
	};
	void attach(T v) {
		if (n == nm) expand(5);
		memcpy(m + n, &v, sizeof(T)); 
		n++;
	};
	~ARRAY_FLEX() {release(true);};
};
#endif



template <class T> class FLEX_ARRAY {
public:
	ARRAY<T*> a;

	T* get(int indx) {
		if (indx < 0) indx += a.n;
		if (indx < 0) return NULL;
		else if (indx < a.n) return a.m[indx];
		else return NULL;
	};
	void attach(T* v) {
		a.expand(1);
		a.m[a.n - 1] = v;
	};
	T** add(int n) {
		int n0 = a.n;
		a.expand(n);
		T* bf = new T[n];
		for (int i = 0; i < n; i++) a.m[n0 + i] = bf + i;
		return a.m + n0;
	};
	T** attach(T* v, int n) {
		int n0 = a.n;
		a.expand(n);
		for (int i = 0; i < n; i++) a.m[n0 + i] = v + i;
		return a.m + n0;
	};
	int indx(T* v) {
		for (int i = 0; i < a.n; i++) {
			if (a.m[i] == v) return i;
		}
		return -1;
	};
	void remove(int indx) {
		return a.remove_one(indx);
	};
	void del(T* v) {
		int vindx = indx(v);
		if (vindx < 0) return;
		remove(vindx); 
	};
	void release(bool release_context) {
		int i;
		if (release_context) {
			for (i = 0; i < a.n; i++) {delete a.m[i]; a.m[i] = NULL;}
		}
		else {
			for (i = 0; i < a.n; i++) a.m[i] = NULL;
		}
		a.release();
	};
	~FLEX_ARRAY() {release(false);};
};


template <class T> class FLEX_MATRIX {
	ARRAY<T> mb; // whole memory block of the matrix
	int nx, ny;
public:
	int n;
	ARRAY<T> *m;
	FLEX_MATRIX() {m = NULL; n = 0; nx = 0; ny = 0;};
	FLEX_MATRIX(T v) {m = NULL; n = 0; nx = 0; ny = 0;};
	void use_GPU_MapMem(bool bGPU_MapMem) {
		if (mb.bGPU_MapMem == bGPU_MapMem) return;
		if (mb.m == NULL) {mb.use_GPU_MapMem(bGPU_MapMem); return;}
		else {
			if (m != NULL) release_pointer();
			mb.use_GPU_MapMem(bGPU_MapMem);
			if (m != NULL) relocate_pointer();
		}
	};
	void release() {
		int i;
		if (m != NULL) {
			for (i = 0; i < n; i++) m[i].release_pointer();
			delete[] m; m = NULL; n = 0;
		}
		mb.release();
	};
	void release_pointer() {
		int i;
		if (m != NULL) {
			for (i = 0; i < n; i++) m[i].release_pointer();
		}
	};
	void release_mblock() {
		mb.release();
	};
	void release_mblock_pointer() {
		mb.m = NULL; mb.n = 0; this->nx = 0; this->ny = 0;
	};
	void set_mblock(T* m, int n) {
		mb.m = m; mb.n = n;
	};
	void set_dim(int nx, int ny) {this->nx = nx; this->ny = ny;};
	void relocate_pointer() {
		int i;
		if (mb.n == NULL || mb.n != nx * ny || m == NULL) return;
		for (i = 0; i < n; i++) {m[i].m = mb.m + i * ny; m[i].n = ny;}
	};
	T* mblock() {return mb.m;};
	int mem_size() {return mb.n;};
	int dim_x() {return nx;};
	int dim_y() {return ny;};
	void set(int n) {
		if (this->n == n) return;
		release(); this->n = n; if (n <= 0) {n = 0; return;}
		mb.set_array(n * n); nx = n; ny = n;
		m = new ARRAY<T>[n];
		int i = 0;
		for (i = 0; i < n; i++) {
			//m[i].SetArray(n);
			m[i].m = mb.m + i * n; m[i].n = n;
		}
	};
	~FLEX_MATRIX() {release();};
	bool check_dim(FLEX_MATRIX<T> &mv) {
		if (n == mv.n && m != NULL) return true;
		release(); if (mv.m == NULL) return true;
		set(mv.n); return false;
	};
	void operator = (FLEX_MATRIX<T> &fm) {
		//int i, j;
		if (this->m == NULL) { // not initialized
			//for (i = 0; i < n; i++) memset(m[i].m, 0, m[i].n * sizeof(T));
			set(fm.n);
		}
		if (fm.m == NULL) return;
		//for (i = 0; i < n; i++) {
		//	memcpy(m[i].m, fm.m[i].m, m[i].n * sizeof(T));
		//}
		memcpy(this->mb.m, fm.mblock(), this->mb.n * sizeof(T)); return;
	};
	void operator = (const int v) {
		if (mb.m == NULL) return;
		memset(this->mb.m, v, this->mb.n * sizeof(T)); return;
	};
	void operator += (FLEX_MATRIX<T> &fm) {
		int i;//, j;
		if (this->mb.n != fm.mem_size()) {this->set(fm.n); memset(this->mb.m, 0, this->mb.n * sizeof(T));}
		//for (i = 0; i < n; i++) for (j = 0; j < n; j++) m[i].m[j] += fm.m[i].m[j];
		T* fmb = fm.mblock(), mb = this->mb.m;
		for (i = 0; i < mb.n; i++) mb[i] += fmb[i];
	};
	void operator -= (FLEX_MATRIX<T> &fm) {
		int i;//, j;
		if (!check_dim(fm)) (*this) = 0;
		if (this->mb.n != fm.mem_size()) {this->set(fm.n); memset(this->mb.m, 0, this->mb.n * sizeof(T));}
		//for (i = 0; i < n; i++) for (j = 0; j < n; j++) m[i].m[j] += fm.m[i].m[j];
		T* fmb = fm.mblock(), mb = this->mb.m;
		for (i = 0; i < mb.n; i++) mb[i] -= fmb[i];
	};
};

template <class T> class FLEX_ASYM_MATRIX {
	ARRAY<T> mb; // whole memory block of the matrix
	int nx, ny; // dimension
public:
	int n;
	ARRAY<T> *m;
	FLEX_ASYM_MATRIX() {m = NULL; n = 0; nx = 0; ny = 0;};
	FLEX_ASYM_MATRIX(T v) {m = NULL; n = 0; nx = 0; ny = 0;};
	void use_GPU_MapMem(bool bGPU_MapMem) {
		if (mb.bGPU_MapMem == bGPU_MapMem) return;
		if (mb.m == NULL) {mb.use_GPU_MapMem(bGPU_MapMem); return;}
		else {
			if (m != NULL) release_pointer();
			mb.use_GPU_MapMem(bGPU_MapMem);
			if (m != NULL) relocate_pointer();
		}
	};
	void release() {
		int i;
		if (m != NULL) {
			for (i = 0; i < n; i++) m[i].release_pointer();
			delete[] m; m = NULL; n = 0; nx = 0; ny = 0;
		}
		mb.release();
	};
	void release_pointer() {
		int i;
		if (m != NULL) {
			for (i = 0; i < n; i++) m[i].release_pointer();
		}
	};
	void release_mblock() {
		mb.release();
	};
	void release_mblock_pointer() {
		mb.m = NULL; mb.n = 0; this->nx = 0; this->ny = 0;
	};
	void set_mblock(T* m, int n) {
		mb.m = m; mb.n = n;
	};
	void set_dim(int nx, int ny) {this->nx = nx; this->ny = ny;};
	void relocate_pointer() {
		int i;
		if (mb.m == NULL || mb.n != nx * ny || m == NULL) return;
		for (i = 0; i < n; i++) {m[i].m = mb.m + i * ny; m[i].n = ny;}
	};
	T* mblock() {return mb.m;};
	int mem_size() {return mb.n;};
	int dim_x() {return nx;};
	int dim_y() {return ny;};
	void set(int n1, int n2) {
		if (this->nx == n1 && this->ny == n2) return;
		release(); this->nx = n1; this->ny = n2; this->n = n1; if (n1 <= 0) return;
		mb.set_array(n1 * n2);
		m = new ARRAY<T>[n1];
		int i = 0;
		for (i = 0; i < n1; i++) {
			//m[i].SetArray(n);
			m[i].m = mb.m + i * n2; m[i].n = n2;
		}
	};
	~FLEX_ASYM_MATRIX() {release();};
	bool check_dim(FLEX_ASYM_MATRIX<T> &mv) {
		if (n == mv.n && m != NULL) return true;
		release(); if (mv.m == NULL) return false;
		set(mv.dim_x(), mv.dim_y()); return false;
	};
	void operator = (FLEX_ASYM_MATRIX<T> &fm) {
		//int i, j;
		check_dim(fm);
		memcpy(mb.m, fm.mblock(), mb.n * sizeof(T));
	};
	void operator += (FLEX_ASYM_MATRIX<T> &fm) {
		int i;
		if (!check_dim(fm)) memset(mb.m, 0, mb.n * sizeof(T));
		T* fmb = fm.mblock(), *mmb = mb.m;
		for (i = 0; i < mb.n; i++) mmb[i] += fmb[i];
	};
	void operator -= (FLEX_ASYM_MATRIX<T> &fm) {
		int i;
		if (!check_dim(fm)) memset(mb.m, 0, mb.n * sizeof(T));
		T* fmb = fm.mblock(), *mmb = mb.m;
		for (i = 0; i < mb.n; i++) mmb[i] -= fmb[i];
	};
};



template <class M, class T> void cmplus_thread(M *m1, M *m2, M *mres, int n1, int n2) {
	int i, j, k;
	int nj = m1->m[0].n, nk = m1->m[0].m[0].n;
	T *v1, *v2, *vres;
	for (i = n1; i <= n2; i++) {
		for (j = 0; j < nj; j++) {
			v1 = m1->m[i].m[j].m; v2 = m2->m[i].m[j].m; vres = mres->m[i].m[j].m;
			for (k = 0; k < nk; k++) {
				vres[k] = v1[k] + v2[k];
			}
		}
	}
};


#if _SYS_ == _WINDOWS_SYS_
template <class CM, class T> int CM3plus_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class CM, class T> void* CM3plus_Thread(void *void_par) { // thread-function
#endif
	CM3OPERATE_THREAD_VARS<CM, T> *par = (CM3OPERATE_THREAD_VARS<CM, T>*)void_par;
	//typedef void (*CM3_OPERATOR)(CM*, CM*, CM*, int, int);
	//((CM3_OPERATOR)(par->func))(par->m1, par->m2, par->mres, par->n1, par->n2);
	cmplus_thread<CM, T>(par->m1, par->m2, par->mres, par->n1, par->n2);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class CM, class T> void CM3_plus(CM &m1, CM &m2, CM &mres, int nThreads) {
	if (nThreads <= 1) {
		//typedef void (*CM3_OPERATOR)(CM*, CM*, CM*, int, int);
		//((CM3_OPERATOR)(func))(&m1, &m2, &mres, 0, m1.n - 1);
		cmplus_thread<CM, T>(&m1, &m2, &mres, 0, m1.n - 1);
		return;
	}

	CM3OPERATE_THREAD_VARS<CM, T> *thread_vars = new CM3OPERATE_THREAD_VARS<CM, T>[nThreads];
	int iThread = 0;
	int n1 = 0, n2 = 0, nstep = m1.n / nThreads;
	//char msg[256] = "\0";
	//extern void show_log(char *msg, bool endline);
	for (iThread = 0; iThread < nThreads; iThread++) {
		nstep = (m1.n - n1) / (nThreads - iThread);
		if (nstep * (nThreads - iThread) < (m1.n - n1)) nstep += 1;
		if (iThread == nThreads - 1) n2 = m1.n - 1;
		else n2 = n1 + nstep - 1;

		if (n2 >= m1.n) n2 = m1.n - 1;

		thread_vars[iThread].set_pars(iThread, &m1, &m2, &mres, n1, n2);
		//sprintf(msg, "%d: [%d - %d], ", iThread, n1, n2); show_log(msg, false);
		//thread_vars[iThread].set_func(func);
		n1 = n2 + 1;
	}
	//show_log(" ", true);

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)CM3plus_Thread<CM, T>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &CM3plus_Thread<CM, T>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		//thread_vars[iThread].mm = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
};


template <class T> class FLEX_CMATRIX { // flexible cubic matrix
	ARRAY<T> mb; // memory block
	int nx, ny, nz;
public:
	int n;
	FLEX_MATRIX<T> *m;
	FLEX_CMATRIX() {m = NULL; n = 0; nx = 0; ny = 0; nz = 0;};
	FLEX_CMATRIX(T v) {m = NULL; n = 0; nx = 0; ny = 0; nz = 0;};
	void use_GPU_MapMem(bool bGPU_MapMem) {
		if (mb.bGPU_MapMem == bGPU_MapMem) return;
		if (mb.m == NULL) {mb.use_GPU_MapMem(bGPU_MapMem); return;}
		else {
			if (m != NULL) release_pointer();
			mb.use_GPU_MapMem(bGPU_MapMem);
			if (m != NULL) relocate_pointer();
		}
	};
	void release() {
		int i;
		if (m == NULL) {mb.release(); nx = 0; ny = 0; nz = 0; n = 0; return;}
		for (i = 0; i < n; i++) {
			m[i].release_pointer(); m[i].release_mblock_pointer(); m[i].set_dim(0, 0);
			m[i].release();
		}
		delete[] m; m = NULL;
		mb.release();
	};
	void release_pointer() {
		int i;
		if (m == NULL) {mb.release(); nx = 0; ny = 0; nz = 0; n = 0; return;}
		for (i = 0; i < n; i++) {
			m[i].release_pointer(); m[i].release_mblock_pointer(); m[i].set_dim(0, 0);
		}
	};
	void release_mblock() {mb.release();};
	void release_mblock_pointer() {mb.m = NULL; mb.n = 0; this->nx = 0; this->ny = 0; this->nz = 0;};
	void relocate_pointer() {
		int i;
		if (m == NULL || mb.m == NULL) return;
		this->n = this->nx;
		for (i = 0; i < n; i++) {
			m[i].set_mblock(mb.m + i * ny * nz, ny * nz); m[i].set_dim(ny, nz); m[i].relocate_pointer(); 
		}
	};
	T* mblock() {return mb.m;};
	int mem_size() {return mb.n;};
	int dim_x() {return nx;};
	int dim_y() {return ny;};
	int dim_z() {return nz;};

	void set_dim(int nx, int ny, int nz) {this->nx = nx; this->ny = ny; this->nz = nz;};
	
	void set(int n) {
		if (this->n == n) return;
		release(); this->n = n; if (n <= 0) {n = 0; return;}
		mb.set_array(n * n * n); nx = n; ny = n; nz = n;
		m = new FLEX_MATRIX<T>[n]; 
		int i = 0, j;
		for (i = 0; i < n; i++) {
			m[i].set_mblock(mb.m + i * ny * nz, ny * nz); m[i].set_dim(ny, nz);
			m[i].m = new ARRAY<T>[ny]; m[i].n = ny;
			for (j = 0; j < ny; j++) m[i].m[j].n = nz;
			m[i].relocate_pointer();
		}
	};
	~FLEX_CMATRIX() {release();};
	bool check_dim(FLEX_CMATRIX<T> &mv) {
		if (n == mv.n && m != NULL) return true;
		release(); if (mv.m == NULL) return true;
		set(mv.n); return false;
	};
	void operator = (FLEX_CMATRIX<T> &fm) {
		if (fm.m == NULL) {memset(mb.m, 0, mb.n * sizeof(T)); return;}
		if (!check_dim(fm)) return;
		memcpy(mb.m, fm.mblock(), mb.n * sizeof(T));
	};
	void operator += (FLEX_CMATRIX<T> &fm) {
		int i;
		T* fmb = fm.mblock();
		if (!check_dim(fm)) memset(mb.m, 0, mb.n * sizeof(T));
		if (n < 5 * MAX_THREADS) {
			//for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++) m[i].m[j].m[k] += fm.m[i].m[j].m[k];
			for (i = 0; i < mb.n; i++) mb.m[i] += fmb[i];
		}
		else {
			CM3_plus< FLEX_CMATRIX<T>, T> ((*this), fm, (*this), MAX_THREADS);
		}
	};
	void operator -= (FLEX_CMATRIX<T> &fm) {
		int i;
		if (!check_dim(fm)) memset(mb.m, 0, mb.n * sizeof(T));
		T* fmb = fm.mblock();
		for (i = 0; i < mb.n; i++) mb.m[i] -= fmb[i];
	};
};

template <class T> class FLEX_ASYM_CMATRIX { // flexible asymmetrical cubic matrix
	ARRAY<T> mb; // memory block
public:
	int nx, ny, nz;
	int n; // same as nx
	FLEX_ASYM_MATRIX<T> *m;
	FLEX_ASYM_CMATRIX() {
		m = NULL; nx = 0; ny = 0; nz = 0; n = 0;
	};
	FLEX_ASYM_CMATRIX(T v) {m = NULL; nx = 0; ny = 0; nz = 0; n = 0;};
	void use_GPU_MapMem(bool bGPU_MapMem) {
		if (mb.bGPU_MapMem == bGPU_MapMem) return;
		if (mb.m == NULL) {mb.use_GPU_MapMem(bGPU_MapMem); return;}
		else {
			if (m != NULL) release_pointer();
			mb.use_GPU_MapMem(bGPU_MapMem);
			if (m != NULL) relocate_pointer();
		}
	};
	void release() {
		int i;
		if (m == NULL) {mb.release(); nx = 0; ny = 0; nz = 0; n = 0; return;}
		for (i = 0; i < n; i++) {
			m[i].release_pointer(); m[i].release_mblock_pointer(); m[i].set_dim(0, 0);
			m[i].release();
		}
		delete[] m; m = NULL;
		mb.release();
	};
	void release_pointer() {
		int i;
		if (m == NULL) {mb.release(); nx = 0; ny = 0; nz = 0; n = 0; return;}
		for (i = 0; i < n; i++) {
			m[i].release_pointer(); m[i].release_mblock_pointer(); m[i].set_dim(0, 0);
		}
	};
	void release_mblock() {mb.release();};
	void release_mblock_pointer() {mb.m = NULL; mb.n = 0; this->nx = 0; this->ny = 0; this->nz = 0;};
	void relocate_pointer() {
		int i;
		if (m == NULL || mb.m == NULL) return;
		this->n = this->nx;
		for (i = 0; i < n; i++) {
			m[i].set_mblock(mb.m + i * ny * nz, ny * nz); m[i].set_dim(ny, nz); m[i].relocate_pointer(); 
		}
	};
	T* mblock() {return mb.m;};
	int mem_size() {return mb.n;};
	int dim_x() {return nx;};
	int dim_y() {return ny;};
	int dim_z() {return nz;};

	void set_dim(int nx, int ny, int nz) {this->nx = nx; this->ny = ny; this->nz = nz;};
	
	void set(int nx, int ny, int nz) {
		if (mb.m != NULL && this->nx == nx && this->ny == ny && this->nz == nz) return;
		release(); this->n = nx; if (this->n <= 0) {this->n = 0; return;}
		mb.set_array(nx * ny * nz);
		m = new FLEX_ASYM_MATRIX<T>[n]; this->nx = nx; this->ny = ny; this->nz = nz;
		int i = 0, j;
		for (i = 0; i < n; i++) {
			m[i].set_mblock(mb.m + i * ny * nz, ny * nz); m[i].set_dim(ny, nz);
			m[i].m = new ARRAY<T>[ny]; m[i].n = ny;
			for (j = 0; j < ny; j++) m[i].m[j].n = nz;
			m[i].relocate_pointer();
		}
	};

	~FLEX_ASYM_CMATRIX() {release();};

	bool check_dim(FLEX_ASYM_CMATRIX<T> &mv) {
		if (n == mv.n && m != NULL) return true;
		release(); if (mv.m == NULL) return true;
		set(mv.nx, mv.ny, mv.nz); return false;
	};
	void operator = (FLEX_ASYM_CMATRIX<T> &fm) {
		if (fm.m == NULL) {memset(mb.m, 0, mb.n * sizeof(T)); return;}
		if (!check_dim(fm)) return;
		memcpy(mb.m, fm.mblock(), mb.n * sizeof(T));
	};
	void operator += (FLEX_ASYM_CMATRIX<T> &fm) {
		int i;
		T* fmb = fm.mblock();
		if (!check_dim(fm)) memset(mb.m, 0, mb.n * sizeof(T));
		if (n < 5 * MAX_THREADS) {
			//for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++) m[i].m[j].m[k] += fm.m[i].m[j].m[k];
			for (i = 0; i < mb.n; i++) mb.m[i] += fmb[i];
		}
		else {
			CM3_plus< FLEX_ASYM_CMATRIX<T>, T> ((*this), fm, (*this), MAX_THREADS);
		}
	};
	void operator -= (FLEX_ASYM_CMATRIX<T> &fm) {
		int i;
		if (!check_dim(fm)) memset(mb.m, 0, mb.n * sizeof(T));
		T* fmb = fm.mblock();
		for (i = 0; i < mb.n; i++) mb.m[i] -= fmb[i];
	};
};

template <class DATA> class BUFFER {
public:
	int n_size, n_block, nbs; // total size of buffer, number of blocks, maximum size of block
	ARRAY< ARRAY<DATA> > d;
	BUFFER(int nbs) {
		this->nbs = nbs; n_block = 0; n_size = 0;
	};
	void release() {
		for (int i = 0; i < d.n; i++) {
			d.m[i].release();
		}
		d.release();
		n_block = 0; n_size = 0;
	};
	~BUFFER() {release();};
	void init(int n_size) {
		if (d.n > 0) release();
		this->n_size = n_size;
		n_block = n_size / nbs;
		if (n_block * nbs < n_size) n_block+=1;
		d.set_array(n_block);
		int n = n_size, i;
		for (i = 0; i < n_block; i++) {
			if (n > nbs) d.m[i].set_array(nbs);
			else d.m[i].set_array(n);
			n -= d.m[i].n;
		}
	};
	DATA* get_data(int ipt) {
		if (ipt >= n_size) return NULL;
		int iblock = ipt / nbs, ia = ipt - iblock * nbs;
		return d.m[iblock].m + ia;
	};
};

