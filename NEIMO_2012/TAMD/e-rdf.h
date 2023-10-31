extern void init_ImplicitSolvation_Pars();
extern void init_LJ_Pars();
extern void construct_LJ12_6_prof();
extern double ClusterImplicitSolvationEnergy(BASIC_CLUSTER *pc, BASIC_CLUSTER *pc1, CMM_CELL3D<BASIC_CLUSTER> *cell);

namespace _erdf_ {
using namespace _distribution_;
using namespace _MM_rdf_;

template <class MM, class CLUSTER, class CMM_CLUSTER> class NEIGHBOR_ERDF_THREAD_VAR : public NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> {
public:
	NDISTRIBT<float, float> eLJ_dt, eImpSolv_dt, edt; // e-distribt
	//NDISTRIBT<float, long> sdt; // distribt of number of samples

	CMM_CELL3D<CMM_CLUSTER> *cmm;
	void *efunc; // the function to calculate the interaction 

	NEIGHBOR_ERDF_THREAD_VAR() {cmm = NULL; efunc = NULL;};
	~NEIGHBOR_ERDF_THREAD_VAR() {cmm = NULL; efunc = NULL;};
};

class NEIGHBOR_ERDF_VAR : public NEIGHBOR_RDF_VAR {
public:
	NDISTRIBT<float, float> eLJ_dt, eImpSolv_dt, edt; // e-distribt
	//NDISTRIBT<float, long> sdt; // distribt of number of samples

	NEIGHBOR_ERDF_VAR() {};
	~NEIGHBOR_ERDF_VAR() {};
};



template <class CMM_CLUSTER> bool construct_cmm(char *mdpar_fname, CMM_CELL3D<CMM_CLUSTER> &cmm) {
	//float size_cluster = 0;//, rLJ_cut = 0, rshell = 0;
	float dx, dy, dz;
	char buffer[256] = "\0", md_cell_fname[256] = "\0";
	int nx = 0, ny = 0, nz = 0;

	if (!read_msg(mdpar_fname, buffer, "[CLUSTER_SIZE]", true) || sscanf(buffer, "%f", &size_cluster) != 1) size_cluster = 2;

	::min_cmm_size = size_cluster * 2;

	cmm.nCheckCluster = 1;

	if (!read_msg(mdpar_fname, buffer, "[MD_CELL_TEMPLATE]", true) || sscanf(buffer, "%s", md_cell_fname) != 1) {
		sprintf(errmsg, "can not get [MD_CELL_TEMPLATE] for cell size in %s", mdpar_fname); return false;
	}
	if (!ReadDefinedCellSize(md_cell_fname, dx, dy, dz)) return false;
	dx /= 2; dy /= 2; dz /= 2;
	nx = int(2 * dx / min_cmm_size); nx = (nx < 1 ? 1 : nx);
	ny = int(2 * dy / min_cmm_size); ny = (ny < 1 ? 1 : ny);
	nz = int(2 * dz / min_cmm_size); nz = (nz < 1 ? 1 : nz);
	cmm.set_cell(nx, ny, nz);
	if (::MD_MODE == SINGLE_MM_FREE_CELL_MD) {
		cmm.set_cell_range(false, -dx, dx, -dy, dy, -dz, dz);
		cmm.set_cell_pos();
	}
	else if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
		cmm.set_cell_range(true, -dx, dx, -dy, dy, -dz, dz);
		cmm.set_cell_pos();
	}
	else if (::MD_MODE == MULTI_MM_FREE_CELL_MD) {
		cmm.set_cell_range(false, -dx, dx, -dy, dy, -dz, dz);
		cmm.set_cell_pos();
	}
	else return false;
	return true;
}

/****************************************************************************
	FUNCTION TO DO ERDF-ANALYSIS OF SINGLE MOLECULAR CONFORMATION
	WITH GIVEN MD SIMULATION, WITH MULTI-THREAD

	this function is modified from rdf_distribt_single_mm in rdf-distribt.h
	by including CMM and interaction calculation, with single-processor
****************************************************************************/

// RDF_VAR has DISTRIBT<> and ndt, for rdf parameter and result
// VM_VAR has the parameter to construct MM_VM
// RDF_THREAD_VAR is the rdf calculation of each thread
// read_molstruct is to read molecular structure from saved MD procedure
// calPosition is to calculate the virtual monomer's position
// vm_constructor is to construct MM_VAR from VM_VAR
// erdf_thread_var_func is to construct ERDF_THREAD_VAR from MM_VM and RDF_VAR
// erdf_func is to calculate the erdf
// efunc is the function to calculate the interaction, used by erdf_func
template <class MM, class CLUSTER, class CMM_CLUSTER, class VM_VAR, class RDF_VAR, class RDF_THREAD_VAR> 
bool erdf_distribt_single_mm(char *mdpar_fname, int nf, int nt, int nstep, 
							void* read_molstruct, void *calPosition, 
							VM_VAR &vm_var, void* vm_constructor, 
							RDF_VAR *rdf_var, int nrdf, 
							void *rdf_thread_var_func, void *erdf_func, void *efunc, int &nconfs, void *accumulate_rdf) 
{
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", mol_template[256] = "\0", buffer[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
	if (sscanf(msg, "%d", &(::MD_MODE)) != 1 || (::MD_MODE != SINGLE_MM_FREE_CELL_MD && ::MD_MODE != SINGLE_CGMM_FREE_CELL_MD)) {
		sprintf(msg, "md is not working on single molecule"); show_msg(msg); return false;
	}
	if (!read_item(mdpar_fname, "[MOLECULE_TEMPLATE]", msg)) return false;
	int single_mol = 0;
	if (sscanf(msg, "%d %s", &single_mol, mol_template) != 2 || single_mol != 1) {
		sprintf(msg, "md is not working on single molecule, or macromolecule template filename is not given");
		show_msg(msg); return false;
	}
	READ_MD_PROC_FILE rmdproc;
	rmdproc.reset_opened_file();
	bool status = rmdproc.init(mdpar_fname) && read_MD_SAVE_Pars(mdpar_fname, rmdproc.msave); //init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	//set_md_mode(MD_MODE);
	//if (MD_MODE == SINGLE_MM_FREE_CELL_MD) set_job(GET_SINGLE_MM_MD_PROC);
	//else if (MD_MODE == SINGLE_CGMM_FREE_CELL_MD) set_job(GET_SINGLE_CGMM_MD_PROC);

	int nThreads = 0, iconf = 0, nconf = 0, iThread, i, j;

#define THREADS  4
	MM mm[THREADS];
	CMM_CELL3D<CMM_CLUSTER> cmm[THREADS];
	MM_VM<MM, CLUSTER> vm[THREADS];
	ARRAY<RDF_THREAD_VAR> rdf_thread_var[THREADS];
	SINGLE_MM_RDF_THREAD_VAR<RDF_THREAD_VAR> thread_var[THREADS];

	typedef void (*RDFVAR_INIT)(RDF_VAR*, RDF_THREAD_VAR*);

	for (i = 0; i < THREADS; i++) {
		status = ConstructSingleChainPolymerFromParamFile(mm[i], mol_template);
		if (!status) {
			sprintf(msg, "failure to construct the single-chain polymer with template %s", mol_template); show_msg(msg); return false;
		}
		else {
			mm[i].calMassCenter(true);
			// atompar_db was setup during constructing the molecules;
			// setup LJ parameters of all atoms
			if (i == 0) {
				init_LJ_Pars(); construct_LJ12_6_prof(::iFormat_LJ);
				if (bImplicitSolvent) init_ImplicitSolvation_Pars();
			}
			mm[i].calClusterSolvationRadius();
		}
		vm_var.m = mm + i;
		if (status) status = construct_virtual_monomers< VIRTUAL_MONOMER<CLUSTER>, VM_VAR >(vm_constructor, &vm_var, &(vm[i].vm), vm[i].nVM);
		if (!status) {
			sprintf(msg, "failure to initilize the virtual monomer"); show_msg(msg); return false;
		}
		vm[i].m = mm + i;
		rdf_thread_var[i].SetArray(nrdf);
		for (j = 0; j < rdf_thread_var[i].n; j++) {
			rdf_thread_var[i].m[j].cmm = cmm + i;
			rdf_thread_var[i].m[j].efunc = efunc;
		}

		thread_var[i].set_func(erdf_func);
		thread_var[i].nrdf = nrdf;
		thread_var[i].rdf_thread_var = rdf_thread_var[i].m;
		for (j = 0; j < nrdf; j++) {
			rdf_thread_var[i].m[j].vm = vm + i;
			((RDFVAR_INIT)rdf_thread_var_func)(rdf_var + j, rdf_thread_var[i].m + j);
		}
	}

	typedef bool (*READ_MOLSTRUCT)(ifstream&, MM&);
	typedef void (*VM_POS)(MM_VM<MM, CLUSTER> &vm);

	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		nThreads = 0;
		strcpy(msg1, "\0");
		for (i = 0; i < THREADS; i++) {
			if (nt > 0 && iconf > nt) break;
			if (!rmdproc.search_proc(iconf, buffer)) break;
			//ReadMoleculeStruct(*(rmdproc.in), mm[i]);
			if (!((READ_MOLSTRUCT)read_molstruct)(*(rmdproc.in), mm[i])) break;
			sprintf(msg, " %d ", iconf); strcat(msg1, msg);

			if (nconf == 0) {
				if (!construct_cmm<CMM_CLUSTER>(mdpar_fname, cmm[i])) return false;
				cmm_init_distribute_cluster(mm[i], cmm[i], true); // set base cluster at the center cell
				cmm_init_subcell<CMM_CLUSTER>(cmm[i]); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
			}

			//for (j = 0; j < vm[i].nVM; j++) calPosition(vm[i].vm + j); //, int method = 0)
			((VM_POS)calPosition)(vm[i]);

			if (nconf == 0) cmm_check_cluster<CMM_CLUSTER>(cmm[i]);
			else cmm_check<CMM_CLUSTER>(cmm[i], MAX_THREADS);
			MMCheckRelshipInCell(mm[i], 0, mm[i].nCluster, cmm[i], MAX_THREADS);
		#if _CLUSTER_IMPLICIT_SOLVATION_ == 1
			MacroMolCellOperate<MM, CMM_CLUSTER>((void*)(&MMClusterImpSolvNeighborsInCell), cmm[i], mm[i], 0, mm[i].nCluster - 1, MAX_THREADS);
		#endif
			iThread = i;
			thread_var[iThread].setThread(iThread);
		#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(thread_var[iThread].hMMThread);
			AfxBeginThread((AFX_THREADPROC)rdf_thread_single_mm<RDF_THREAD_VAR>, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
		#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(thread_var[iThread].thread), NULL, &rdf_thread_single_mm<RDF_THREAD_VAR>, (void *)(thread_var + iThread));
		#endif

			iconf += nstep;
			nThreads++; nconf++;
		}
		if (nThreads == 0) break;
		sprintf(msg, "%sin %s are read", msg1, mdpar_fname); show_log(msg, false);
		
	#if _SYS_ == _WINDOWS_SYS_ 
		for (iThread = 0; iThread < nThreads; iThread++) {
			WaitForSingleObject(thread_var[iThread].hMMThread, INFINITE);
		}
	#elif _SYS_ == _LINUX_SYS_
		for (iThread = 0; iThread < nThreads; iThread++) {
			if (thread_var[iThread].thread != 0) pthread_join(thread_var[iThread].thread, NULL); // wait the thread over
		}
	#endif
		sprintf(msg, ", and analyzed"); show_log(msg, true);
	}
	rmdproc.reset_opened_file();

	typedef void (*ACCUMULATE)(RDF_VAR*, RDF_THREAD_VAR*);
	
	for (int irdf = 0; irdf < nrdf; irdf++) {
		for (iThread = 0; iThread < THREADS; iThread++) 
			((ACCUMULATE)accumulate_rdf)(rdf_var + irdf, rdf_thread_var[iThread].m + irdf);
	}
	nconfs = nconf;

	return true;
#undef THREADS
};


/****************************************************************************
	FUNCTION TO DO RDF-ANALYSIS OF MULTI-MOLECULAR CONFORMATION
	WITH GIVEN MD SIMULATION, WITH MULTI-THREAD

	this function is modified from rdf_distribt_multi_mm in rdf-distribt.h
	by including CMM and interaction calculation, with single-processor
****************************************************************************/

// md_cell has all the macromolecules, because we use construct_macromolecule(**) to construct
//        multi-mm system as a global variable, md_cell has to be a given global variable
//        relating to construct_macromolecule(**)

// RDF_VAR has DISTRIBT<> and ndt, for rdf parameter and result
// VM_VAR has the parameter to construct MM_VM
// RDF_THREAD_VAR is the rdf calculation of each thread
// read_molstruct is to read molecular structure from saved MD procedure
// calPosition is to calculate the virtual monomer's position
// vm_constructor is to construct MM_VAR from VM_VAR
// rdf_thread_var_func is to construct RDF_THREAD_VAR from MM_VM and RDF_VAR
// rdf_func is to calculate the rdf
// efunc is the function to calculate the interaction, used by erdf_func

template <class MD_CELL, class MM, class CLUSTER, class CMM_CLUSTER, class VM_VAR, class RDF_VAR, class RDF_THREAD_VAR> 
bool erdf_distribt_multi_mm(MD_CELL &md_cell, char *mdpar_fname, int nf, int nt, int nstep, 
						   ARRAY<short> &mol_id, void *read_molstruct, void *calPosition, 
						   VM_VAR &vm_var, void* vm_constructor, 
						   RDF_VAR *rdf_var, int nrdf, 
						   void *rdf_thread_var_func, void *erdf_func, void *efunc, int &nconfs, void *accumulate_rdf) 
{
	nconfs = 0;
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

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int nThreads = 0, iconf = 0, nconf = 0, iThread, imol = 0, i, j;

	init_LJ_Pars(); construct_LJ12_6_prof(::iFormat_LJ);
	if (bImplicitSolvent) init_ImplicitSolvation_Pars();
	for (imol = 0; imol < md_cell.mm.n; imol++) {
		md_cell.mm.m[imol].calClusterSolvationRadius();
	}

	typedef void (*RDFVAR_INIT)(RDF_VAR*, RDF_THREAD_VAR*);

	MULTI_MM_RDF_THREAD_VAR<MD_CELL, MM_VM<MM, CLUSTER>, RDF_THREAD_VAR> thread_var[MAX_THREADS];
	CMM_CELL3D<CMM_CLUSTER> cmm;
	ARRAY< MM_VM<MM, CLUSTER> > vm;
	FLEX_ASYM_MATRIX<RDF_THREAD_VAR> mrdf_thread_var;
	vm.SetArray(md_cell.mm.n);
	mrdf_thread_var.set(md_cell.mm.n, nrdf);
	for (i = 0; i < vm.n; i++) {
		vm_var.m = md_cell.mm.m + i;
		if (!construct_virtual_monomers<VIRTUAL_MONOMER<CLUSTER>, VM_VAR>(vm_constructor, &vm_var, &(vm.m[i].vm), vm.m[i].nVM)) {
			sprintf(msg, "failure to construct virtual monomer of molecule %d", i);
			show_msg(msg); return false;
		}
		vm.m[i].m = md_cell.mm.m + i;
		for (j = 0; j < nrdf; j++) {
			mrdf_thread_var.m[i].m[j].vm = vm.m + i;
			mrdf_thread_var.m[i].m[j].efunc = efunc;
			mrdf_thread_var.m[i].m[j].cmm = &cmm;
			((RDFVAR_INIT)rdf_thread_var_func)(rdf_var + j, mrdf_thread_var.m[i].m + j);
		}
	}

	typedef void (*VM_POS)(ARRAY< MM_VM<MM, CLUSTER> > &vm);
	typedef bool (*READ_MMStruct)(int);

	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		nThreads = 0;
		//if (!serial_multi_mol_md_proc(iconf)) break;
		if (!((READ_MMStruct)read_molstruct)(iconf)) break;
		sprintf(msg, " %d in %s are read", iconf, mdpar_fname); show_log(msg, true);

		((VM_POS)calPosition)(vm);

		if (nconf == 0) {
			if (!construct_cmm<CMM_CLUSTER>(mdpar_fname, cmm)) return false;
			for (imol = 0; imol < md_cell.mm.n; imol++) {
				cmm_init_distribute_cluster(md_cell.mm.m[imol], cmm, (imol == 0 ? true : false)); // set base cluster at the center cell
			}
			cmm_init_subcell<CMM_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
		}
		if (nconf == 0) cmm_check_cluster<CMM_CLUSTER>(cmm);
		else cmm_check<CMM_CLUSTER>(cmm, MAX_THREADS);

		iThread = 0; nThreads = 0; show_log("analyzing molecules ", false);
		for (imol = 0; imol < md_cell.mm.n; imol++) { // working on macro-molecules only
			if (mol_id.m != NULL && !mol_id.exist(imol)) continue;
			thread_var[iThread].md_cell = &md_cell;
			thread_var[iThread].vm_array = &vm;
			thread_var[iThread].rdf_thread_var = mrdf_thread_var.m[imol].m;
			thread_var[iThread].nrdf = mrdf_thread_var.m[imol].n;
			thread_var[iThread].rdf_func = erdf_func;
			thread_var[iThread].setThread(iThread);

			MMCheckRelshipInCell(md_cell.mm.m[imol], 0, md_cell.mm.m[imol].nCluster, cmm, MAX_THREADS);
		#if _CLUSTER_IMPLICIT_SOLVATION_ == 1
			MacroMolCellOperate<MM, CMM_CLUSTER>((void*)(&MMClusterImpSolvNeighborsInCell), cmm, md_cell.mm.m[imol], 0, md_cell.mm.m[imol].nCluster - 1, MAX_THREADS);
		#endif

		#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(thread_var[iThread].hMMThread);
			AfxBeginThread((AFX_THREADPROC)rdf_thread_multi_mm<MD_CELL, MM_VM<MM, CLUSTER>, RDF_THREAD_VAR>, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
		#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(thread_var[iThread].thread), NULL, &rdf_thread_multi_mm<MD_CELL, MM_VM<MM, CLUSTER>, RDF_THREAD_VAR>, (void *)(thread_var + iThread));
		#endif

			sprintf(msg, "%d ", imol); show_log(msg, false);
			iThread++; nThreads++;


			if (iThread == MAX_THREADS || imol == md_cell.mm.n - 1) { // reaching maximum threads or all molecules are run
				if (nThreads == 0) break; // no thread to wait
			#if _SYS_ == _WINDOWS_SYS_ 
				for (iThread = 0; iThread < nThreads; iThread++) {
					WaitForSingleObject(thread_var[iThread].hMMThread, INFINITE);
				}
			#elif _SYS_ == _LINUX_SYS_
				for (iThread = 0; iThread < nThreads; iThread++) {
					if (thread_var[iThread].thread != 0) pthread_join(thread_var[iThread].thread, NULL); // wait the thread over
				}
			#endif
				iThread = 0; nThreads = 0;
			}
		}
		show_log("", true);

		iconf += nstep;
		nconf++;
	}

	clear_read_mdproc();

	typedef void (*ACCUMULATE)(RDF_VAR*, RDF_THREAD_VAR*);
	
	for (int irdf = 0; irdf < nrdf; irdf++) {
		for (iThread = 0; iThread < mrdf_thread_var.n; iThread++) 
			((ACCUMULATE)accumulate_rdf)(rdf_var + irdf, mrdf_thread_var.m[iThread].m + irdf);
	}
	if (mol_id.m == NULL) nconfs = nconf * ::mm_md_cell.mm.n;
	else nconfs = nconf * mol_id.n;
	return true;
}


/*********************************************************************************
	FUNCTION TO DO RDF-ANALYSIS OF SINGLE-MOLECULE WITH SEVERAL MD SIMULATION

	this functions are copied exactly from rdf-distrbt.h, 
	except the function parameters
**********************************************************************************/

template <class MM, class CLUSTER, class CMM_CLUSTER, class VM_VAR, class RDF_VAR, class RDF_THREAD_VAR> 
bool erdf_distribt_single_mm_multi_md(char *parfname, void *read_molstruct, void *calPosition, 
									 VM_VAR &vm_var, void* vm_constructor, 
									 RDF_VAR* rdf_var, int nrdf, 
									 void *erdf_thread_var_func, void *erdf_func, void *efunc,
									 void *accumulate_rdf, void *normalize_rdf, void *save_rdf) 
{
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	if (!read_item(parfname, "[OUTPUT]", buffer)) {sprintf(msg, "can not find [OUTPUT] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%s", title) != 1) {sprintf(msg, "can not get monomer to analyze, and the output filename after [OUTPUT]"); show_msg(msg); return false;}

	int nMD = 1;
	if (!read_item(parfname, "[MD_CONFIGS]", buffer)) {sprintf(msg, "can not find [MD_CONFIGS] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%d", &nMD) != 1 || nMD <= 0) {sprintf(msg, "nMD has to be >= 1, see %s", buffer); show_msg(msg); return false;}
	
	int nf = 0, nt = 0, nstep = 1;
	char mdpar_fname[256] = "\0";
	int nconfs = 0, nconf = 0, n, MD_MODE = 0;
	bool status = true;

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

		if (MD_MODE == SINGLE_MM_FREE_CELL_MD || MD_MODE == SINGLE_CGMM_FREE_CELL_MD) {
			status = erdf_distribt_single_mm<MM, CLUSTER, CMM_CLUSTER, VM_VAR, RDF_VAR, RDF_THREAD_VAR>(mdpar_fname, nf, nt, nstep, 
			read_molstruct, calPosition, vm_var, vm_constructor, rdf_var, nrdf, erdf_thread_var_func, erdf_func, efunc, nconf, accumulate_rdf);
			if (!status) return false;
		}
		else {
			sprintf(msg, "MD in %s is not for single mm", mdpar_fname); show_msg(msg, true); nconf = 0;
		}

		nconfs += nconf;
		n++;
	}
	in_cfg.close();

	typedef void (*NORMALIZE)(RDF_VAR*);
	typedef void (*SAVE_RDF)(RDF_VAR*, char*, char*);
	for (n = 0; n < nrdf; n++) {
		((NORMALIZE)normalize_rdf)(rdf_var + n);
		sprintf(ofname, "%s-%drdf", title, n);
		sprintf(msg, "%drdf ", n);
		((SAVE_RDF)save_rdf)(rdf_var + n, ofname, msg);
	}

	return true;
};

/*********************************************************************************
	FUNCTION TO DO RDF-ANALYSIS OF MULTI-MOLECULE WITH SEVERAL MD SIMULATION

	this functions are copied exactly from rdf-distrbt.h, 
	except the function parameters
**********************************************************************************/

template <class MD_CELL, class MM, class CLUSTER, class CMM_CLUSTER, class VM_VAR, class RDF_VAR, class RDF_THREAD_VAR> 
bool erdf_distribt_multi_mm_multi_md(char *parfname, MD_CELL &md_cell, 
									void *read_molstruct, void *calPosition, 
									VM_VAR &vm_var, void* vm_constructor, 
									RDF_VAR* rdf_var, int nrdf, 
									void *erdf_thread_var_func, void *erdf_func, void *efunc,
									void *accumulate_rdf, void *normalize_rdf, void *save_rdf) 
{
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	if (!read_item(parfname, "[OUTPUT]", buffer)) {sprintf(msg, "can not find [OUTPUT] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%s", title) != 1) {sprintf(msg, "can not get monomer to analyze, and the output filename after [OUTPUT]"); show_msg(msg); return false;}

	int nMD = 1;
	if (!read_item(parfname, "[MD_CONFIGS]", buffer)) {sprintf(msg, "can not find [MD_CONFIGS] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%d", &nMD) != 1 || nMD <= 0) {sprintf(msg, "nMD has to be >= 1, see %s", buffer); show_msg(msg); return false;}
	
	int nf = 0, nt = 0, nstep = 1;
	char mdpar_fname[256] = "\0";
	int nconfs = 0, nconf = 0, n, MD_MODE = 0;
	bool status = true;

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

		ARRAY<short> mol_id;

		if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD || MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			read_given_id(parfname, "[MOLECULES]", mol_id);
			status = erdf_distribt_multi_mm<MD_CELL, MM, CLUSTER, CMM_CLUSTER, VM_VAR, RDF_VAR, RDF_THREAD_VAR>(md_cell, mdpar_fname, nf, nt, nstep, mol_id, 
				read_molstruct, calPosition, vm_var, vm_constructor, rdf_var, nrdf, erdf_thread_var_func, erdf_func, efunc, nconf, accumulate_rdf);
			if (!status) return false;
		}
		else {
			sprintf(msg, "MD in %s, is not multi-mm", mdpar_fname); show_msg(msg, true); nconf = 0;
		}

		nconfs += nconf;
		n++;
	}
	in_cfg.close();

	typedef void (*NORMALIZE)(RDF_VAR*);
	typedef void (*SAVE_RDF)(RDF_VAR*, char*, char*);

	for (n = 0; n < nrdf; n++) {
		((NORMALIZE)normalize_rdf)(rdf_var + n);
		sprintf(ofname, "%s-%drdf", title, n);
		sprintf(msg, "%drdf ", n);
		((SAVE_RDF)save_rdf)(rdf_var + n, ofname, msg);
	}

	return true;
};

/********************************************************************
	RDF FUNCTION BETWEEN NEIGHBORS / NON-NEIGHBORS
********************************************************************/


template <class MM, class CLUSTER, class CMM_CLUSTER> void single_vm_neighbor_erdf_func1(NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *erdf_thread_var) {
	int idt = 0, im, im1, npt = 0;
	NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *var;
	cNeighbor< VIRTUAL_MONOMER<CLUSTER> > *neighbor = NULL;
	_GROUP<CLUSTER> *pgroup = NULL;
	VIRTUAL_MONOMER<CLUSTER> *vm = NULL, *vm1 = NULL;
	VECTOR3 r0, r1, dr;
	double R = 0, U_LJ = 0, U_ImpSolv = 0;

	typedef void (*EFUNC)(VIRTUAL_MONOMER<CLUSTER>*, VIRTUAL_MONOMER<CLUSTER>*, CMM_CELL3D<CMM_CLUSTER>*, double&, double&);

	var = erdf_thread_var;
	for (im = 0; im < var->neighbor.n; im++) {
		vm = var->vm->vm + im;
		if (var->cm_uid > 0 && vm->g.uid != var->cm_uid) continue;
		V32V3(vm->r, r0)
		neighbor = var->neighbor.neighbor + im;
		for (im1 = 0; im1 < neighbor->neighbor.n; im1++) {
			pgroup = &(neighbor->neighbor.m[im1]->g);
			if (var->nm_uid > 0 && pgroup->uid != var->nm_uid) continue;
			vm1 = neighbor->neighbor.m[im1];
			V32V3(vm1->r, r1)
			VECT3(r0, r1, dr)
			scale_uv3(dr, dr, R) R = sqrt(R);

			((EFUNC)(var->efunc))(vm, vm1, var->cmm, U_LJ, U_ImpSolv);

			for (idt = 0; idt < var->edt.ndt; idt++) {
				DISTRIBT_INDX(var->edt.dt[idt], R, npt) 
				if (LInRange(npt, 0, var->edt.dt[idt].x.N)) {
					var->eLJ_dt.dt[idt].add_sample(npt, U_LJ, true);
					var->eImpSolv_dt.dt[idt].add_sample(npt, U_ImpSolv, true);
					var->edt.dt[idt].add_sample(npt, U_LJ + U_ImpSolv, true);
					//var->sdt.dt[idt].add_sample(npt, true);

					//var->sdt.dt[idt].vnormal += 1; // total number of samples
					var->eLJ_dt.dt[idt].vnormal += 1; // total number of samples
					var->eImpSolv_dt.dt[idt].vnormal += 1; // total number of samples
					var->edt.dt[idt].vnormal += 1; // total number of samples
				}
			}
		}
	}
}

template <class MM, class CLUSTER, class CMM_CLUSTER> void single_vm_neighbor_erdf_func(NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *erdf_thread_var, int nrdf) {
	int irdf = 0;
	NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *var;

	for (irdf = 0; irdf < nrdf; irdf++) {
		var = erdf_thread_var + irdf;
		single_vm_neighbor_erdf_func1<MM, CLUSTER, CMM_CLUSTER>(var);
	}
}


// checking the neighbor vm distribution of multi-mm system, vmarray is not useful
template <class MD_CELL, class MM, class CLUSTER, class CMM_CLUSTER> void multi_vm_neighbor_erdf_func(MD_CELL *md_cell, ARRAY< MM_VM<MM, CLUSTER> > *vmarray, NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *erdf_thread_var, int nrdf) {
	single_vm_neighbor_erdf_func<MM, CLUSTER, CMM_CLUSTER>(erdf_thread_var, nrdf);
}

template <class MM, class CLUSTER, class CMM_CLUSTER> void single_vm_non_neighbor_erdf_func1(NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *erdf_thread_var) {
	int idt = 0, im, im1, npt = 0;
	NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *var;
	cNeighbor< VIRTUAL_MONOMER<CLUSTER> > *neighbor = NULL;
	_GROUP<CLUSTER> *pgroup = NULL;
	VIRTUAL_MONOMER<CLUSTER> *vm = NULL, *vm1 = NULL;
	VECTOR3 r0, r1, dr;
	double R = 0, U_LJ = 0, U_ImpSolv = 0;

	typedef void (*EFUNC)(VIRTUAL_MONOMER<CLUSTER>*, VIRTUAL_MONOMER<CLUSTER>*, CMM_CELL3D<CMM_CLUSTER>*, double&, double&);

	var = erdf_thread_var;
	for (im = 0; im < var->vm->nVM; im++) {
		vm = var->vm->vm + im;
		if (var->cm_uid > 0 && vm->g.uid != var->cm_uid) continue;
		V32V3(vm->r, r0)
		neighbor = var->neighbor.neighbor + im;
		for (im1 = 0; im1 < var->vm->nVM; im1++) {
			if (im == im1) continue;
			vm1 = var->vm->vm + im1;
			pgroup = &(vm1->g);
			if (var->nm_uid > 0 && pgroup->uid != var->nm_uid) continue;
			if (neighbor->neighbor.exist(vm1)) continue;
				
			V32V3(vm1->r, r1)
			VECT3(r0, r1, dr)
			scale_uv3(dr, dr, R) R = sqrt(R);

			((EFUNC)(var->efunc))(vm, vm1, var->cmm, U_LJ, U_ImpSolv);

			for (idt = 0; idt < var->edt.ndt; idt++) {
				DISTRIBT_INDX(var->edt.dt[idt], R, npt) 
				if (LInRange(npt, 0, var->edt.dt[idt].x.N)) {
					var->eLJ_dt.dt[idt].add_sample(npt, U_LJ, true);
					var->eImpSolv_dt.dt[idt].add_sample(npt, U_ImpSolv, true);
					var->edt.dt[idt].add_sample(npt, U_LJ + U_ImpSolv, true);
					//var->sdt.dt[idt].add_sample(npt, true);

					//var->sdt.dt[idt].vnormal += 1; // total number of samples
					var->eLJ_dt.dt[idt].vnormal += 1; // total number of samples
					var->eImpSolv_dt.dt[idt].vnormal += 1; // total number of samples
					var->edt.dt[idt].vnormal += 1; // total number of samples
				}
			}
		}
	}
}

template <class MM, class CLUSTER, class CMM_CLUSTER> void single_vm_non_neighbor_erdf_func(NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *erdf_thread_var, int nrdf) {
	int irdf = 0;
	NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *var;

	for (irdf = 0; irdf < nrdf; irdf++) {
		var = erdf_thread_var + irdf;
		single_vm_non_neighbor_erdf_func1<MM, CLUSTER, CMM_CLUSTER>(var);
	}
}


template <class MD_CELL, class MM, class CLUSTER, class CMM_CLUSTER> void multi_vm_non_neighbor_erdf_func1(MD_CELL *md_cell, ARRAY< MM_VM<MM, CLUSTER> >* vmarray, NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *erdf_thread_var) {
	int idt = 0, im, im1, npt = 0, imol = 0;
	NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *var;
	MM_VM<MM, CLUSTER> *mvm = NULL, *mvm1 = NULL;
	bool same_mol = false;

	cNeighbor< VIRTUAL_MONOMER<CLUSTER> > *neighbor = NULL;
	_GROUP<CLUSTER> *pgroup = NULL;
	VIRTUAL_MONOMER<CLUSTER> *vm = NULL, *vm1 = NULL;
	VECTOR3 r0, r1, dr;
	double R = 0, U_LJ = 0, U_ImpSolv = 0;

	typedef void (*EFUNC)(VIRTUAL_MONOMER<CLUSTER>*, VIRTUAL_MONOMER<CLUSTER>*, CMM_CELL3D<CMM_CLUSTER>*, double&, double&);

	var = erdf_thread_var;
	mvm = var->vm;
	for (im = 0; im < mvm->nVM; im++) {
		vm = mvm->vm + im;
		if (var->cm_uid > 0 && vm->g.uid != var->cm_uid) continue;
		V32V3(vm->r, r0)
		neighbor = var->neighbor.neighbor + im;

		for (imol = 0; imol < vmarray->n; imol++) {
			mvm1 = vmarray->m + imol;
			same_mol = (mvm == mvm1 ? true : false);
			for (im1 = 0; im1 < mvm1->nVM; im1++) {
				vm1 = mvm1->vm + im1;
				pgroup = &(vm1->g);
				if (var->nm_uid > 0 && pgroup->uid != var->nm_uid) continue;
				if (same_mol) {
					if (im == im1) continue;
					if (neighbor->neighbor.exist(vm1)) continue;
				}
				
				V32V3(vm1->r, r1)
				VECT3(r0, r1, dr)
				scale_uv3(dr, dr, R) R = sqrt(R);

				((EFUNC)(var->efunc))(vm, vm1, var->cmm, U_LJ, U_ImpSolv);

				for (idt = 0; idt < var->edt.ndt; idt++) {
					DISTRIBT_INDX(var->edt.dt[idt], R, npt) 
					if (LInRange(npt, 0, var->edt.dt[idt].x.N)) {
						var->eLJ_dt.dt[idt].add_sample(npt, U_LJ, true);
						var->eImpSolv_dt.dt[idt].add_sample(npt, U_ImpSolv, true);
						var->edt.dt[idt].add_sample(npt, U_LJ + U_ImpSolv, true);
						//var->sdt.dt[idt].add_sample(npt, true);

						//var->sdt.dt[idt].vnormal += 1; // total number of samples
						var->eLJ_dt.dt[idt].vnormal += 1; // total number of samples
						var->eImpSolv_dt.dt[idt].vnormal += 1; // total number of samples
						var->edt.dt[idt].vnormal += 1; // total number of samples
					}
				}
			}
		}
	}
}

template <class MD_CELL, class MM, class CLUSTER, class CMM_CLUSTER> void multi_vm_non_neighbor_erdf_func(MD_CELL *md_cell, ARRAY< MM_VM<MM, CLUSTER> >* vmarray, NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *erdf_thread_var, int nrdf) {
	int irdf = 0;
	NEIGHBOR_ERDF_THREAD_VAR<MM, CLUSTER, CMM_CLUSTER> *var;

	for (irdf = 0; irdf < nrdf; irdf++) {
		var = erdf_thread_var + irdf;
		multi_vm_non_neighbor_erdf_func1<MD_CELL, MM, CLUSTER, CMM_CLUSTER>(md_cell, vmarray, var);
	}
}


} // end of namespace _erdf_
