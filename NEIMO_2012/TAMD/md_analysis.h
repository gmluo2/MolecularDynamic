namespace _md_analysis_ {

template <class VAR, class RES> class AnalysisThread : public _THREAD_VARS {
public:
	VAR *var;
	RES *res;
	void *func;
};

#if _SYS_ == _WINDOWS_SYS_
template <class VAR, class RES> int analysis_thread(LPVOID *void_par) // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class VAR, class RES> void* analysis_thread(void *void_par) // thread-function
#endif
{
	typedef void (*FUNC) (VAR*, RES*);
	AnalysisThread<VAR, RES> *par = (AnalysisThread<VAR, RES>*)void_par;
	((FUNC)(par->func))(par->var, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}


template <class MD_CELL, class VAR, class RES> 
bool analysis_single_md_multithread(MD_CELL &md_cell, char *mdpar_fname, int nf, int nt, int nstep, 
						   void *read_molstruct, void *post_init_mdcell, void *update_cell, VAR *var, RES *res, int nThreads, 
						   void *init_var, void *update_var, void *analysis_func, int &nconfs) 
{
	char msg[256] = "\0", msg1[256] = "\0", buffer[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
	int MD_MODE = 0;
	if (sscanf(msg, "%d", &MD_MODE) != 1 || (MD_MODE != MULTI_MM_FREE_CELL_MD && MD_MODE != MULTI_MM_PERIDIC_CELL_MD
		&& MD_MODE != MULTI_CGMM_FREE_CELL_MD && MD_MODE != MULTI_CGMM_PERIDIC_CELL_MD)) { // MultiMol ?
		sprintf(msg, "md is not working on multi-molecule");
		show_msg(msg); return false;
	}

	set_md_mode(MD_MODE);
	if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) set_job(GET_MULTI_MM_MD_PROC);
	else if (MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) set_job(GET_MULTI_CGMM_MD_PROC);

	typedef void (*POST_INIT_MDCELL)(MDCELL&);

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	if (post_init_mdcell != NULL) ((POST_INIT_MDCELL)post_init_mdcell)(mdcell);
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int iconf = 0, nconf = 0, iThread, nthread = 0;

	typedef void (*INIT_VAR)(MD_CELL&, VAR&);
	ARRAY< AnalysisThread<VAR, RES> > thread_var; thread_var.SetArray(nThreads);
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (init_var != NULL) {
			((INIT_VAR)init_var)(md_cell, var[iThread]);
		}

		thread_var.m[iThread].set_pars(iThread);
		thread_var.m[iThread].var = var + iThread;
		thread_var.m[iThread].res = res + iThread;
		thread_var.m[iThread].func = analysis_func;
	}

	typedef void (*UPDATE_VAR)(MD_CELL&, VAR&);
	typedef bool (*READ_Moltruct)(int);
	typedef void (*UPDATE_CELL)(MD_CELL&);

	sprintf(msg, "analyzing md config. in %s : ", mdpar_fname); show_log(msg, true);
	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		
		nthread = 0; 
		for (iThread = 0; iThread < nThreads; iThread++) {
			//if (!serial_multi_mol_md_proc(iconf)) break;
			if (!((READ_Moltruct)read_molstruct)(iconf)) break;
			sprintf(msg, " %d, ", iconf); show_log(msg, false);

			if (update_cell != NULL) ((UPDATE_CELL)update_cell)(md_cell);
			if (update_var != NULL) ((UPDATE_VAR)update_var)(md_cell, var[iThread]);

		#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(thread_var.m[iThread].hMMThread);
			AfxBeginThread((AFX_THREADPROC)analysis_thread<VAR, RES>, (LPVOID)(thread_var.m + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
		#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(thread_var.m[iThread].thread), NULL, &analysis_thread<VAR, RES>, (void *)(thread_var.m + iThread));
		#endif

			nthread++; iconf += nstep; nconf++;
		}
		if (nthread > 0) show_log("are read, ", false);

		if (nthread == 0) break; // no thread to wait
		#if _SYS_ == _WINDOWS_SYS_ 
		for (iThread = 0; iThread < nthread; iThread++) {
			WaitForSingleObject(thread_var.m[iThread].hMMThread, INFINITE);
		}
		#elif _SYS_ == _LINUX_SYS_
		for (iThread = 0; iThread < nthread; iThread++) {
			if (thread_var.m[iThread].thread != 0) pthread_join(thread_var.m[iThread].thread, NULL); // wait the thread over
		}
		#endif
		show_log("analyzed.", true);
	}

	clear_read_mdproc();

	nconfs = nconf;
	
	return true;
}

/*********************************************************************************
	MULTI-THREAD FUNCTION TO DO RDF-ANALYSIS OF MULTI-MOLECULE WITH SEVERAL MD SIMULATION
**********************************************************************************/
template <class MD_CELL, class VAR, class RES> 
bool analysis_multi_md_multithread(char *parfname, MD_CELL &md_cell, VAR *var, RES *res,
						    void *read_molstruct, void *post_init_mdcell, void *update_cell, 
							void *init_var_template, void *init_res_template, /* from the template of var and res to initialize a new var and res */
							void* init_var, void *update_var,  /* from the MDCELL to initialize the var, and update the var */
							void *accumulate_res, void *handle_res,
							void *analysis_func, int &nconfs) 
{
	VAR tvar[MAX_THREADS];
	RES tres[MAX_THREADS];
	int i;
	
	typedef void (*INIT_VAR)(VAR&, VAR&); // init_var_template(VAR &dest, VAR &source) -- using source as a template to initialize the destination
	typedef void (*INIT_RES)(RES&, RES&); // init_res_template(RES &dest, RES &source) -- using source as a template to initialize the destination
	for (i = 0; i < MAX_THREADS; i++) {
		if (init_var_template != NULL) ((INIT_VAR)init_var_template)(tvar[i], *var);
		if (init_res_template != NULL) ((INIT_RES)init_res_template)(tres[i], *res);
	}

	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	if (!read_item(parfname, "[OUTPUT]", buffer)) {sprintf(msg, "can not find [OUTPUT] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%s", title) != 1) {sprintf(msg, "can not get monomer to analyze, and the output filename after [OUTPUT]"); show_msg(msg); return false;}

	int nMD = 1;
	if (!read_item(parfname, "[MD_CONFIGS]", buffer)) {sprintf(msg, "can not find [MD_CONFIGS] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%d", &nMD) != 1 || nMD <= 0) {sprintf(msg, "nMD has to be >= 1, see %s", buffer); show_msg(msg); return false;}
	
	int nf = 0, nt = 0, nstep = 1;
	char mdpar_fname[256] = "\0";
	int nconf = 0, n, MD_MODE = 0;
	bool status = true;

	nconfs = 0;

	ifstream in_cfg;
	in_cfg.open(parfname);
	search(in_cfg, "[MD_CONFIGS]", buffer);
	in_cfg.getline(buffer, 250); // ignore nMD
	
	n = 0; nconfs = 0;
	while (n < nMD) {
		in_cfg.getline(buffer, 250);
		if (sscanf(buffer, "%s %d %d %d", mdpar_fname, &nf, &nt, &nstep) != 4 || nstep < 1) break;

		if (!read_item(mdpar_fname, "[MD_CELL]", msg)) {
			sprintf(msg, "[MD_CELL] is not defined in %s, UNKNOWN MD parameter file", mdpar_fname);
			show_msg(msg, true); break;
		}
		if (sscanf(msg, "%d", &MD_MODE) != 1) { // MultiMol ?
			sprintf(msg, "md is not working on multi-molecule");
			show_msg(msg); break;
		}

		if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD || MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			status = analysis_single_md_multithread<MD_CELL, VAR, RES >(md_cell, mdpar_fname, nf, nt, nstep,  
				read_molstruct, post_init_mdcell, update_cell, tvar, tres, MAX_THREADS, init_var, update_var, analysis_func, nconf);
			if (!status) return false;
		}
		else {
			sprintf(msg, "MD in %s, is not multi-mm", mdpar_fname); show_msg(msg, true);
		}

		if (status) nconfs += nconf;

		n++;
	}
	in_cfg.close();

	typedef void (*RES_ACCUMULATE)(RES&, RES&); // accumulate_res(RES& dest, RES& source) -- accumulate source to destination
	for (i = 0; i < MAX_THREADS; i++) ((RES_ACCUMULATE)accumulate_res)(*res, tres[i]);
	sprintf(msg, "totally %d configs. are analyzed", nconfs); show_log(msg, true);

	typedef void (*HANDLE_RES)(RES&);
	if (handle_res != NULL) ((HANDLE_RES)handle_res)(*res);
	return true;
};



// single-thread function

template <class MD_CELL, class VAR, class RES> 
bool analysis_single_md_singlethread(MD_CELL &md_cell, char *mdpar_fname, int nf, int nt, int nstep, 
						   VAR* var, RES* res, void *read_molstruct, void *post_init_mdcell, void *update_cell, void* init_var, void *update_var, 
						   void *analysis_func, int &nconfs) 
{
	char msg[256] = "\0", msg1[256] = "\0", buffer[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
	int MD_MODE = 0;
	if (sscanf(msg, "%d", &MD_MODE) != 1 || (MD_MODE != MULTI_MM_FREE_CELL_MD && MD_MODE != MULTI_MM_PERIDIC_CELL_MD
		&& MD_MODE != MULTI_CGMM_FREE_CELL_MD && MD_MODE != MULTI_CGMM_PERIDIC_CELL_MD)) { // MultiMol ?
		sprintf(msg, "md is not working on multi-molecule");
		show_msg(msg); return false;
	}

	set_md_mode(MD_MODE);
	if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) set_job(GET_MULTI_MM_MD_PROC);
	else if (MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) set_job(GET_MULTI_CGMM_MD_PROC);

	typedef void (*POST_INIT_MDCELL)(MDCELL&);

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	if (post_init_mdcell != NULL) ((POST_INIT_MDCELL)post_init_mdcell)(mdcell);
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int iconf = 0, nconf = 0;

	typedef void (*INIT_VAR)(MD_CELL&, VAR&);
	typedef void (*UPDATE_VAR)(MD_CELL&, VAR&);
	typedef bool (*READ_Moltruct)(int);
	typedef void (*UPDATE_CELL)(MD_CELL&);
	//typedef void (*UPDATE_CELL)(MD_CELL&, CELL&);

	((INIT_VAR)init_var)(md_cell, *var);

	typedef void (*FUNC) (VAR&, RES&);

	sprintf(msg, "analyzing md config. in %s : ", mdpar_fname);
	show_log(msg, true);
	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		
		//if (!serial_multi_mol_md_proc(iconf)) break;
		if (!((READ_Moltruct)read_molstruct)(iconf)) break;
		sprintf(msg, " %d is read ", iconf); show_log(msg, false);

		if (update_cell != NULL) ((UPDATE_CELL)update_cell)(md_cell);
		if (update_vcell != NULL) ((UPDATE_VAR)update_var)(md_cell, *var);

		((FUNC)(analysis_func))(var, res);

		iconf += nstep; nconf++;

		show_log(", analyzed.", true);
	}

	clear_read_mdproc();

	nconfs = nconf;
	
	return true;
}

/*********************************************************************************
	SINGLE-THREAD FUNCTION TO DO RDF-ANALYSIS OF MULTI-MOLECULE WITH SEVERAL MD SIMULATION
**********************************************************************************/
template <class MD_CELL, class VAR, class RES> 
bool analysis_multi_md_singlethread(char *parfname, MD_CELL &md_cell, VAR *var, RES *res, 
						    void *read_molstruct, void *update_cell, void* init_var, void *update_var, void *init_res_template,
							void *analysis_func, void *accumulate_func, void *handle_res, int &nconfs) 
{
	int i;

	RES tres;
	typedef void (*INIT_RES)(RES&, RES&); // init_res_template(RES &dest, RES& source); // initialize destination with source as template
	((INIT_RES)init_res_template)(tres, res);

	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	if (!read_item(parfname, "[OUTPUT]", buffer)) {sprintf(msg, "can not find [OUTPUT] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%s", title) != 1) {sprintf(msg, "can not get monomer to analyze, and the output filename after [OUTPUT]"); show_msg(msg); return false;}

	int nMD = 1;
	if (!read_item(parfname, "[MD_CONFIGS]", buffer)) {sprintf(msg, "can not find [MD_CONFIGS] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%d", &nMD) != 1 || nMD <= 0) {sprintf(msg, "nMD has to be >= 1, see %s", buffer); show_msg(msg); return false;}
	
	int nf = 0, nt = 0, nstep = 1;
	char mdpar_fname[256] = "\0";
	int nconf = 0, n, MD_MODE = 0;
	bool status = true;

	ifstream in_cfg;
	in_cfg.open(parfname);
	search(in_cfg, "[MD_CONFIGS]", buffer);
	in_cfg.getline(buffer, 250); // ignore nMD

	typedef void (*ACCUMULATE_RES)(RES&, RES&);
	
	n = 0; nconfs = 0;
	while (n < nMD) {
		in_cfg.getline(buffer, 250);
		if (sscanf(buffer, "%s %d %d %d", mdpar_fname, &nf, &nt, &nstep) != 4 || nstep < 1) break;

		if (!read_item(mdpar_fname, "[MD_CELL]", msg)) {
			sprintf(msg, "[MD_CELL] is not defined in %s, UNKNOWN MD parameter file", mdpar_fname);
			show_msg(msg, true); break;
		}
		if (sscanf(msg, "%d", &MD_MODE) != 1) { // MultiMol ?
			sprintf(msg, "md is not working on multi-molecule");
			show_msg(msg); break;
		}

		if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD || MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			status = analysis_single_md_singlethread<MD_CELL, VAR, RES>(md_cell, mdpar_fname, nf, nt, nstep,  
				var, &tres, read_molstruct, update_cell, init_var, update_var, analysis_func, nconf);
			if (!status) return false;
		}
		else {
			sprintf(msg, "MD in %s, is not multi-mm", mdpar_fname); show_msg(msg, true);
		}

		n++;
		if (status) {
			((ACCUMULATE_RES)accumulate_res)(res, tres);
			nconfs += nconf;
		}
	}
	in_cfg.close();

	sprintf(msg, "totally %d configs. are analyzed", nconfs); show_log(msg, true);

	typedef void (*HANDLE_RES)(RES& res);
	if (handle_res != NULL) ((HANDLE_RES)handle_res)(res);

	return true;
};


}  // end of namesapce _md_analysis_


