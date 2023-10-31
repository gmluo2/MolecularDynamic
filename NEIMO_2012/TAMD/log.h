class SerialObject {
public:
#if _SYS_ == _WINDOWS_SYS_
	HANDLE hMutex;
#elif _SYS_ == _LINUX_SYS_
	pthread_mutex_t mutex;
#endif

	SerialObject() {
	#if _SYS_ == _WINDOWS_SYS_
		hMutex = NULL;
		hMutex = CreateMutex(NULL, FALSE, NULL);
	#elif _SYS_ == _LINUX_SYS_
		//mutex = PTHREAD_MUTEX_INITIALIZER;
		pthread_mutex_init(&mutex, NULL);
	#endif
	};
	void release_mutex() {
	#if _SYS_ == _WINDOWS_SYS_
		if (hMutex != NULL) {CloseHandle(hMutex); hMutex = NULL;}
	#elif _SYS_ == _LINUX_SYS_
		//mutex = PTHREAD_MUTEX_INITIALIZER;
	#endif
	};
	~SerialObject() {
		release_mutex();
	};
	// before use the serial object, please wait for the object avaliable
	void wait() {
	#if _SYS_ == _WINDOWS_SYS_
		WaitForSingleObject(hMutex, INFINITE);
	#elif _SYS_ == _LINUX_SYS_
		pthread_mutex_lock(&mutex);
	#endif
	};
	// after using the serial object, please call leave to release the object
	void leave() {
	#if _SYS_ == _WINDOWS_SYS_
		ReleaseMutex(hMutex);
	#elif _SYS_ == _LINUX_SYS_
		pthread_mutex_unlock(&mutex);
	#endif
	};
};

#ifndef of_pos
#define of_pos streampos
#endif

inline of_pos current_ofpos(ofstream &out) {
	return out.tellp();
}

inline void set_ofpos(ofstream &out, of_pos &pos) {
	of_pos cpos = out.tellp();
	out.seekp(pos - cpos, ios_base::cur);
}

class LOG : public SerialObject {
public:
	ofstream *out;
	char fname[256];
	bool bKeepFileClose;
	
	of_pos fpos;
	bool bReplaceLog;

	int nMaxLines, nLines;

	LOG() {out = NULL; strcpy(fname, "\0"); this->bKeepFileClose = true; this->bReplaceLog = false; fpos = 0; nLines = 0; nMaxLines = -1;};
	~LOG() {close();};
	bool init(char *fname, bool bKeepFileClose = true, bool bReplaceLog = false) {
		strcpy(this->fname, fname);
		out = new ofstream;
		out->open(fname, ios_base::out | ios_base::trunc);
		bool status = out->is_open(); 
		if (bKeepFileClose) {
			out->close();
			delete out; out = NULL;
		}
		fpos = 0;
		this->bKeepFileClose = bKeepFileClose;
		this->bReplaceLog = bReplaceLog;
		return status;
	};
	void close(bool wait_event = true) {
		if (wait_event) wait(); // wait for the object for avaliable
		if (out != NULL) {
			if (out->is_open()) {out->flush(); out->close();}; 
			delete out; out = NULL;
		}
		if (wait_event) leave(); // after use it, release the usage of the object
		fpos = 0;
	};
	void record_pos() {
		if (out != NULL) fpos = current_ofpos(*out);
		else fpos = 0;
	};
	void go_end() {
		if (out != NULL) out->seekp(0, ios_base::end);
	};
	void go_begin() {
		if (out != NULL) out->seekp(0, ios_base::beg);
	};
	void show(char *msg, bool endline = true, bool onscreen = false, bool bReplaceLog = true) {
		wait(); // wait for the object for avaliable
		bool status = false;
		if (strlen(fname) > 0) {
			if (out == NULL) {
				out = new ofstream;
				out->open(fname, ios_base::out | ios_base::app);
			}
			status = out->is_open();
			if (status) {
				if (bReplaceLog && this->bReplaceLog && nLines == 0) set_ofpos(*out, fpos);
				(*out)<<msg; if (endline) (*out)<<endl;
			}
			if (bKeepFileClose) close(false);
		}
		else {
#if _SYS_ == _DOS_SYS_ || _SYS_ == _LINUX_SYS_
			cout<<msg;
			if (endline) cout<<endl;
			cout.flush();
#elif _SYS_ == _WINDOWS_SYS_
			//AfxMessageBox(LPCTSTR(CString(msg)));
#endif
		}
		if (endline) {
			nLines++;
			if (nLines == nMaxLines) nLines = 0;
		}
		leave(); // after use it, release the usage of the object
	};
};


