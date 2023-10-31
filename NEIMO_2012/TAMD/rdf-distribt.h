namespace _MM_rdf_ {
using namespace _distribution_;

template <class T> class cNeighbor {
public:
	ARRAY<T*> neighbor;
	int nth;

	cNeighbor() {nth = 0;};
	~cNeighbor() {
		for (int i = 0; i < neighbor.n; i++) neighbor.m[i] = NULL;
		neighbor.release();
	};
};

template <class T> class cNeighbors {
public:
	cNeighbor<T> *neighbor;
	int n;

	cNeighbors() {neighbor = NULL; n = 0;};
	void release() {
		if (neighbor == NULL) return;
		delete[] neighbor; neighbor = NULL; n = 0;
	};
	~cNeighbors() {release();};
	void set_elements(int n, int nth) {
		release(); if (n <= 0) return;
		neighbor = new cNeighbor<T>[n]; this->n = n;
		for (int i = 0; i < n; i++) neighbor[i].nth = nth;
	};
};

// cneighbor is an array, with m->nCluster of elements
void init_cluster_neighbor(MMOLECULE *m, cNeighbor<CLUSTER> *cneighbor);

// DERIVED_CLASS is derived from VIRTUAL_MONOMER
template <class DERIVED_CLASS> class _VIRTUAL_CLUSTER { 
public:
	int nHinges;
	HINGE<DERIVED_CLASS> *hinge;
	_VIRTUAL_CLUSTER() {hinge = NULL; nHinges = 0;};
	void reset_hinges() {
		if (hinge != NULL) {delete[] hinge; hinge = NULL;}
		nHinges = 0;
	};
	~_VIRTUAL_CLUSTER() {reset_hinges();};
	void set_hinges(int nHinges) {
		reset_hinges(); if (nHinges <= 0) return;
		hinge = new HINGE<DERIVED_CLASS>[nHinges]; this->nHinges = nHinges;
	};
	bool setup_link_relationship(int iHinge, DERIVED_CLASS *c, int jHinge) {
		hinge[iHinge].plink = c; c->hinge[jHinge].plink = (DERIVED_CLASS*)this;
		hinge[iHinge].indx_link = jHinge; c->hinge[jHinge].indx_link = iHinge;
		hinge[iHinge].Olink = c->hinge[jHinge].O;
		c->hinge[jHinge].Olink = hinge[iHinge].O;
		
		return true;
	};
};

template <class T> class VIRTUAL_MONOMER : public _VIRTUAL_CLUSTER< VIRTUAL_MONOMER<T> > {
public:
	_GROUP<T> g;
	VECTOR3 r;
	VIRTUAL_MONOMER() {};
	void reset_hinges() {
		_VIRTUAL_CLUSTER< VIRTUAL_MONOMER<T> >::reset_hinges();
		//if (indx != NULL) {delete[] indx; indx = NULL;}
		//if (iHinge != NULL) {delete[] iHinge; iHinge = NULL;}
	};
	void set_hinges(int nHinges) {
		_VIRTUAL_CLUSTER< VIRTUAL_MONOMER<T> >::set_hinges(nHinges);
		//if (nHinges > 0) {
		//	indx = new int[nHinges]; iHinge = new int[nHinges];
		//}
	};
	~VIRTUAL_MONOMER() {
		g.reset(); reset_hinges();
	};
};


#define _USER_MONOMER_UID   500000

template <class T> void construct_virtual_monomer(CHAIN<T> *ch, VIRTUAL_MONOMER<T> *monomer) {
	monomer->g.reset(); if (ch == NULL) return;
	CHAIN<T> *pch = ch;
	T* p = NULL;
	int n = number<T>(ch);
	monomer->g.set(n);
	n = 0;
	while (pch != NULL) {monomer->g.u[n] = pch->p; n++; pch = pch->next;}
};

class HINGE_INDX {
public:
	int indx, iHinge;
};

template <class T> void virtual_monomer_setup_hinges(VIRTUAL_MONOMER<T> *monomer) { // if T has hinge also
	if (monomer->g.nUnit == 0) return;

	CHAIN<HINGE_INDX> *mhinge_ch = NULL, *tail = NULL, *pch = NULL;
	T *pc = NULL, *pc1 = NULL;
	int i, j, iHinge = 0, nHinges = 0;
	bool free_hinge = false;
	for (i = 0; i < monomer->g.nUnit; i++) {
		pc = monomer->g.u[i];
		for (iHinge = 0; iHinge < pc->nHinges; iHinge++) {
			free_hinge = true;
			pc1 = (T*)(pc->hinge[iHinge].plink);
			for (j = 0; j < monomer->g.nUnit; j++) {
				if (j == i) continue;
				if (pc1 == monomer->g.u[j]) {free_hinge = false; break;}
			}
			if (free_hinge) {
				pch = new CHAIN<HINGE_INDX>; pch->p = new HINGE_INDX; pch->p->indx = i; pch->p->iHinge = iHinge;
				if (mhinge_ch == NULL) mhinge_ch = pch;
				else tail->next = pch;
				tail = pch;
			}
		}
	}
	nHinges = number<HINGE_INDX>(mhinge_ch);
	if (nHinges == 0) return;
	monomer->set_hinges(nHinges);
	pch = mhinge_ch; i = 0;
	while (pch != NULL) {
		monomer->hinge[i].indx = pch->p->indx; monomer->hinge[i].ibound = pch->p->iHinge;
		pch = pch->next; i++;
	}
	release_chain<HINGE_INDX>(&mhinge_ch, true);
};

template <class T> void connect_virtual_monomers(VIRTUAL_MONOMER<T> *monomer, int nMonomer) { // if T has hinge also
	char msg[256] = "\0";
	VIRTUAL_MONOMER<T> *pm = NULL;
	T *pc = NULL, *pc1 = NULL, *pc2 = NULL, *pc3 = NULL;
	int i, j, iHinge = 0, jHinge = 0;
	bool hinge_match;
	for (i = 0; i < nMonomer; i++) {
		pm = monomer + i;
		for (iHinge = 0; iHinge < pm->nHinges; iHinge++) {
			hinge_match = false;
			pc = pm->g.u[pm->hinge[iHinge].indx];
			pc1 = (T*)(pc->hinge[pm->hinge[iHinge].ibound].plink);
			for (j = 0; j < nMonomer; j++) {
				for (jHinge = 0; jHinge < monomer[j].nHinges; jHinge++) {
					pc2 = (T*)(monomer[j].g.u[monomer[j].hinge[jHinge].indx]);
					pc3 = (T*)(pc2->hinge[monomer[j].hinge[jHinge].ibound].plink);
					if (pc2 == pc1 && pc == pc3) {
						hinge_match = true; break;
					}
				}
				if (hinge_match) break;
			}
			if (hinge_match) pm->setup_link_relationship(iHinge, monomer + j, jHinge);
			else {
				sprintf(msg, "virtual monomer %d has free hinge %d", i, iHinge);
				show_log(msg, true);
			}
		}
	}
};

bool construct_monomer_single_chain_homopolymer(MMOLECULE *m, int m_uid, int nm_cg, int vm_uid, VIRTUAL_MONOMER<CLUSTER> **cgmonomer, int& ncgMonomers);

template <class MM> class VM_VAR {
public:
	MM *m;
	int nBlocks; // number of blocks of the copolymer
	int *muid, *nm_vm, *vmuid;
	int ncgm;
	VM_VAR() {m = NULL; nBlocks = 0; muid = NULL; nm_vm = NULL; vmuid = NULL;};
	void reset() {
		if (muid != NULL) {delete[] muid; muid = NULL;}
		if (nm_vm != NULL) {delete[] nm_vm; nm_vm = NULL;}
		if (vmuid != NULL) {delete[] vmuid; vmuid = NULL;}
		ncgm = 0;
	};
	void set_vm(int ncgm) {
		reset(); if (ncgm < 1) return;
		muid = new int[ncgm]; nm_vm = new int[ncgm]; vmuid = new int[ncgm]; this->ncgm = ncgm;
	};
	~VM_VAR() {reset();};
};

bool construct_monomer_single_chain_block_copolymer(VM_VAR<MMOLECULE>* vm_var, VIRTUAL_MONOMER<CLUSTER> **cgmonomer, int& ncgMonomers);


// mneighbor is an array, with m->nMonomer of elements
template <class MMOLECULE, class CLUSTER> void init_virtual_monomer_neighbor(VIRTUAL_MONOMER<CLUSTER> *m, int nm, cNeighbor< VIRTUAL_MONOMER<CLUSTER> > *mneighbor) {
	VIRTUAL_MONOMER<CLUSTER> *pm = NULL;
	CHAIN< VIRTUAL_MONOMER<CLUSTER> > *ch = NULL;
	for (int i = 0; i < nm; i++) {
		grandchildren_chain< VIRTUAL_MONOMER<CLUSTER> >(m + i, mneighbor[i].nth, &ch);
		chain_array< VIRTUAL_MONOMER<CLUSTER> >(&ch, mneighbor[i].neighbor);
		release_chain< VIRTUAL_MONOMER<CLUSTER> >(&ch, false);
	}
}

// mneighbor is an array, with m->nMonomer of elements
// the neighbor includes all the neighbors within nth_neighbor
template <class MMOLECULE, class CLUSTER> void init_virtual_monomer_all_neighbors(VIRTUAL_MONOMER<CLUSTER> *m, int nm, cNeighbor< VIRTUAL_MONOMER<CLUSTER> > *mneighbor) {
	VIRTUAL_MONOMER<CLUSTER> *pm = NULL;
	CHAIN< VIRTUAL_MONOMER<CLUSTER> > *ch = NULL;
	for (int i = 0; i < nm; i++) {
		children_chain< VIRTUAL_MONOMER<CLUSTER> >(m + i, mneighbor[i].nth, &ch);
		chain_array< VIRTUAL_MONOMER<CLUSTER> >(&ch, mneighbor[i].neighbor);
		release_chain< VIRTUAL_MONOMER<CLUSTER> >(&ch, false);
	}
}

template <class VIRTUAL_MONOMER, class VM_VAR> bool construct_virtual_monomers(void *func, VM_VAR *var, VIRTUAL_MONOMER **cgmonomer, int& ncgMonomers) {
	typedef bool (*CONSTRUCTOR)(VM_VAR*, VIRTUAL_MONOMER**, int&);
	return ((CONSTRUCTOR)(func))(var, cgmonomer, ncgMonomers);
};

template <class MM, class CLUSTER> class MM_VM {
public:
	MM *m;
	VIRTUAL_MONOMER<CLUSTER> *vm;
	int nVM;

	MM_VM() {m = NULL; vm = NULL; nVM = 0;};
	void reset() {
		if (vm != NULL) {delete[] vm; vm = NULL;} 
		nVM = 0;
	};
	~MM_VM() {reset();};
	void set_vm(int nVM) {
		reset(); if (nVM < 1) return;
		vm = new VIRTUAL_MONOMER<CLUSTER>[nVM]; this->nVM = nVM;
	};
};

void calPosition(VIRTUAL_MONOMER<CLUSTER> *vm, int method = 0);

template <class MM, class CLUSTER> class NEIGHBOR_RDF_THREAD_VAR {
public:
	int cm_uid, nm_uid; // center monomer's uid, neibhbor-monomer's uid
	int nth_neighbor;

	MM_VM<MM, CLUSTER> *vm;
	cNeighbors< VIRTUAL_MONOMER<CLUSTER> > neighbor;

	NDISTRIBT<float, long> dt;

	NEIGHBOR_RDF_THREAD_VAR() {cm_uid = 0; nm_uid = 0; nth_neighbor = 0; vm = NULL;};
	~NEIGHBOR_RDF_THREAD_VAR() {neighbor.release();};
};

class NEIGHBOR_RDF_VAR {
public:
	int cm_uid, nm_uid; // center monomer's uid, neibhbor-monomer's uid
	int nth_neighbor;

	NDISTRIBT<float, float> dt;

	NEIGHBOR_RDF_VAR() {};
	~NEIGHBOR_RDF_VAR() {};
};

template <class RDF_THREAD_VAR> class SINGLE_MM_RDF_THREAD_VAR : public _THREAD_VARS {
public:
	RDF_THREAD_VAR *rdf_thread_var;
	int nrdf;
	void *rdf_func; // the function to do rdf-analysis
	SINGLE_MM_RDF_THREAD_VAR() {rdf_thread_var = NULL; rdf_func = NULL; nrdf = 0;};
	~SINGLE_MM_RDF_THREAD_VAR() {rdf_func = NULL; rdf_thread_var = NULL; nrdf = 0;};
	void set_func(void* rdf_func) {this->rdf_func = (void*)rdf_func;};
	void setThread(int iThread) {
		_THREAD_VARS::set_pars(iThread);
	};
};

template <class MMOLECULE, class CLUSTER> void init_vm_neighbor_rdf_thread_var(NEIGHBOR_RDF_VAR* rdf_var, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER>* rdf_thread_var) {
	setup_ndt(rdf_thread_var->dt, rdf_var->dt);
	rdf_thread_var->cm_uid = rdf_var->cm_uid;
	rdf_thread_var->nm_uid = rdf_var->nm_uid;
	rdf_thread_var->nth_neighbor = rdf_var->nth_neighbor;
	rdf_thread_var->neighbor.set_elements(rdf_thread_var->vm->nVM, rdf_var->nth_neighbor);
	init_virtual_monomer_neighbor<MMOLECULE, CLUSTER>(rdf_thread_var->vm->vm, rdf_thread_var->vm->nVM, rdf_thread_var->neighbor.neighbor);
}

template <class MMOLECULE, class CLUSTER> void init_vm_all_neighbors_rdf_thread_var(NEIGHBOR_RDF_VAR* rdf_var, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER>* rdf_thread_var) {
	setup_ndt(rdf_thread_var->dt, rdf_var->dt);
	rdf_thread_var->cm_uid = rdf_var->cm_uid;
	rdf_thread_var->nm_uid = rdf_var->nm_uid;
	rdf_thread_var->nth_neighbor = rdf_var->nth_neighbor;
	rdf_thread_var->neighbor.set_elements(rdf_thread_var->vm->nVM, rdf_var->nth_neighbor);
	init_virtual_monomer_all_neighbors<MMOLECULE, CLUSTER>(rdf_thread_var->vm->vm, rdf_thread_var->vm->nVM, rdf_thread_var->neighbor.neighbor);
}

#if _SYS_ == _WINDOWS_SYS_
template <class RDF_THREAD_VAR> int rdf_thread_single_mm(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class RDF_THREAD_VAR> void* rdf_thread_single_mm(void *void_par) { // thread-function
#endif
	typedef void (*RDF_FUNC) (RDF_THREAD_VAR*, int);
	SINGLE_MM_RDF_THREAD_VAR<RDF_THREAD_VAR> *par = (SINGLE_MM_RDF_THREAD_VAR<RDF_THREAD_VAR>*)void_par;
	((RDF_FUNC)(par->rdf_func))(par->rdf_thread_var, par->nrdf);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

/****************************************************************************
	FUNCTION TO DO RDF-ANALYSIS OF SINGLE MOLECULAR CONFORMATION
	WITH GIVEN MD SIMULATION, WITH MULTI-THREAD
****************************************************************************/

// RDF_VAR has DISTRIBT<> and ndt, for rdf parameter and result
// VM_VAR has the parameter to construct MM_VM
// RDF_THREAD_VAR is the rdf calculation of each thread
// read_molstruct is to read molecular structure from saved MD procedure
// calPosition is to calculate the virtual monomer's position
// vm_constructor is to construct MM_VAR from VM_VAR
// rdf_thread_var_func is to construct RDF_THREAD_VAR from MM_VM and RDF_VAR
// rdf_func is to calculate the rdf
template <class MM, class CLUSTER, class VM_VAR, class RDF_VAR, class RDF_THREAD_VAR> 
bool rdf_distribt_single_mm(char *mdpar_fname, int nf, int nt, int nstep, 
							void* read_molstruct, void *calPosition, 
							VM_VAR &vm_var, void* vm_constructor, 
							RDF_VAR *rdf_var, int nrdf, 
							void *rdf_thread_var_func, void *rdf_func, int &nconfs, void* accumulate_rdf) 
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

	MM mm[MAX_THREADS];
	MM_VM<MM, CLUSTER> vm[MAX_THREADS];
	ARRAY<RDF_THREAD_VAR> rdf_thread_var[MAX_THREADS];
	SINGLE_MM_RDF_THREAD_VAR<RDF_THREAD_VAR> thread_var[MAX_THREADS];

	typedef void (*RDFVAR_INIT)(RDF_VAR*, RDF_THREAD_VAR*);

	for (i = 0; i < MAX_THREADS; i++) {
		if (status) status = ConstructSingleChainPolymerFromParamFile(mm[i], mol_template);
		vm_var.m = mm + i;
		if (status) status = construct_virtual_monomers< VIRTUAL_MONOMER<CLUSTER>, VM_VAR >(vm_constructor, &vm_var, &(vm[i].vm), vm[i].nVM);
		if (!status) {
			sprintf(msg, "failure to initilize the virtual monomer"); show_msg(msg); return false;
		}
		vm[i].m = mm + i;
		rdf_thread_var[i].SetArray(nrdf);

		thread_var[i].set_func(rdf_func);
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
		for (i = 0; i < MAX_THREADS; i++) {
			if (nt > 0 && iconf > nt) break;
			if (!rmdproc.search_proc(iconf, buffer)) break;
			//ReadMoleculeStruct(*(rmdproc.in), mm[i]);
			if (!((READ_MOLSTRUCT)read_molstruct)(*(rmdproc.in), mm[i])) break;
			sprintf(msg, " %d ", iconf); strcat(msg1, msg);

			//for (j = 0; j < vm[i].nVM; j++) calPosition(vm[i].vm + j); //, int method = 0)
			((VM_POS)calPosition)(vm[i]);

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
		for (iThread = 0; iThread < MAX_THREADS; iThread++) 
			((ACCUMULATE)accumulate_rdf)(rdf_var + irdf, rdf_thread_var[iThread].m + irdf);
	}
	nconfs = nconf;

	return true;
};

template <class MD_CELL, class MM_VM, class RDF_THREAD_VAR> class MULTI_MM_RDF_THREAD_VAR : public _THREAD_VARS {
public:
	ARRAY<MM_VM> *vm_array;
	RDF_THREAD_VAR *rdf_thread_var;
	int nrdf;

	MD_CELL *md_cell;
	void *rdf_func; // function to do rdf-analysis

	MULTI_MM_RDF_THREAD_VAR() {md_cell = NULL; rdf_thread_var = NULL; rdf_func = NULL; nrdf = 0; vm_array = NULL;};
	~MULTI_MM_RDF_THREAD_VAR() {md_cell = NULL; rdf_thread_var = NULL; rdf_func = NULL;};
	void set_func(void* rdf_func) {this->rdf_func = (void*)rdf_func;};
	void setThread(int iThread) {
		_THREAD_VARS::set_pars(iThread);
	};
};

#if _SYS_ == _WINDOWS_SYS_
template <class MD_CELL, class MM_VM, class RDF_THREAD_VAR> int rdf_thread_multi_mm(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MD_CELL, class MM_VM, class RDF_THREAD_VAR> void* rdf_thread_multi_mm(void *void_par) { // thread-function
#endif
	typedef void (*RDF_FUNC) (MD_CELL*, ARRAY<MM_VM>*, RDF_THREAD_VAR*, int);
	MULTI_MM_RDF_THREAD_VAR<MD_CELL, MM_VM, RDF_THREAD_VAR> *par = (MULTI_MM_RDF_THREAD_VAR<MD_CELL, MM_VM, RDF_THREAD_VAR>*)void_par;
	((RDF_FUNC)(par->rdf_func))(par->md_cell, par->vm_array, par->rdf_thread_var, par->nrdf);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

/****************************************************************************
	FUNCTION TO DO RDF-ANALYSIS OF MULTI-MOLECULAR CONFORMATION
	WITH GIVEN MD SIMULATION, WITH MULTI-THREAD
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

template <class MD_CELL, class MM, class CLUSTER, class VM_VAR, class RDF_VAR, class RDF_THREAD_VAR> 
bool rdf_distribt_multi_mm(MD_CELL &md_cell, char *mdpar_fname, int nf, int nt, int nstep, 
						   ARRAY<short> &mol_id, void *read_molstruct, void *calPosition, 
						   VM_VAR &vm_var, void* vm_constructor, 
						   RDF_VAR *rdf_var, int nrdf, 
						   void *rdf_thread_var_func, void *rdf_func, int &nconfs, void *accumulate_rdf) 
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

	typedef void (*RDFVAR_INIT)(RDF_VAR*, RDF_THREAD_VAR*);

	MULTI_MM_RDF_THREAD_VAR<MD_CELL, MM_VM<MM, CLUSTER>, RDF_THREAD_VAR> thread_var[MAX_THREADS];
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

		iThread = 0; nThreads = 0; show_log("analyzing molecules ", false);
		for (imol = 0; imol < md_cell.mm.n; imol++) { // working on macro-molecules only
			if (mol_id.m != NULL && !mol_id.exist(imol)) continue;
			thread_var[iThread].md_cell = &md_cell;
			thread_var[iThread].vm_array = &vm;
			thread_var[iThread].rdf_thread_var = mrdf_thread_var.m[imol].m;
			thread_var[iThread].nrdf = mrdf_thread_var.m[imol].n;
			thread_var[iThread].rdf_func = rdf_func;
			thread_var[iThread].setThread(iThread);

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
**********************************************************************************/

template <class MM, class CLUSTER, class VM_VAR, class RDF_VAR, class RDF_THREAD_VAR> 
bool rdf_distribt_single_mm_multi_md(char *parfname, void *read_molstruct, void *calPosition, 
									 VM_VAR &vm_var, void* vm_constructor, 
									 RDF_VAR* rdf_var, int nrdf, 
									 void *rdf_thread_var_func, void *rdf_func, 
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
			status = rdf_distribt_single_mm<MM, CLUSTER, VM_VAR, RDF_VAR, RDF_THREAD_VAR>(mdpar_fname, nf, nt, nstep, 
			read_molstruct, calPosition, vm_var, vm_constructor, rdf_var, nrdf, rdf_thread_var_func, rdf_func, nconf, accumulate_rdf);
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
**********************************************************************************/

template <class MD_CELL, class MM, class CLUSTER, class VM_VAR, class RDF_VAR, class RDF_THREAD_VAR> 
bool rdf_distribt_multi_mm_multi_md(char *parfname, MD_CELL &md_cell, 
									void *read_molstruct, void *calPosition, 
									VM_VAR &vm_var, void* vm_constructor, 
									RDF_VAR* rdf_var, int nrdf, 
									void *rdf_thread_var_func, void *rdf_func, 
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
			status = rdf_distribt_multi_mm<MD_CELL, MM, CLUSTER, VM_VAR, RDF_VAR, RDF_THREAD_VAR>(md_cell, mdpar_fname, nf, nt, nstep, mol_id, 
				read_molstruct, calPosition, vm_var, vm_constructor, rdf_var, nrdf, rdf_thread_var_func, rdf_func, nconf, accumulate_rdf);
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
template <class MM, class CLUSTER, class Tx, class Ty> void single_vm_neighbor_rdf_func1(NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, NDISTRIBT<Tx, Ty> &distrbt) {
	int idt = 0, im, im1, npt = 0;
	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *var;
	cNeighbor< VIRTUAL_MONOMER<CLUSTER> > *neighbor = NULL;
	_GROUP<CLUSTER> *pgroup = NULL;
	VIRTUAL_MONOMER<CLUSTER> *vm = NULL, *vm1 = NULL;
	VECTOR3 r0, r1, dr;
	double R = 0;

	var = rdf_thread_var;
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

			for (idt = 0; idt < distrbt.ndt; idt++) {
				if (distrbt.dt[idt].sample(R, true)) distrbt.dt[idt].vnormal += 1; // total number of samples
			}
		}
	}
}

template <class MM, class CLUSTER> void single_vm_neighbor_rdf_func(NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, int nrdf) {
	int irdf = 0;
	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *var;

	for (irdf = 0; irdf < nrdf; irdf++) {
		var = rdf_thread_var + irdf;
		single_vm_neighbor_rdf_func1<MM, CLUSTER, float, long>(var, var->dt);
	}
}

// checking the neighbor vm distribution of multi-mm system, vmarray is not useful
template <class MD_CELL, class MM, class CLUSTER> void multi_vm_neighbor_rdf_func(MD_CELL *md_cell, ARRAY< MM_VM<MM, CLUSTER> > *vmarray, NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, int nrdf) {
	single_vm_neighbor_rdf_func<MM, CLUSTER>(rdf_thread_var, nrdf);
}

template <class MM, class CLUSTER, class Tx, class Ty> void single_vm_non_neighbor_rdf_func1(NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, NDISTRIBT<Tx, Ty> &distrbt) {
	int idt = 0, im, im1, npt = 0;
	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *var;
	cNeighbor< VIRTUAL_MONOMER<CLUSTER> > *neighbor = NULL;
	_GROUP<CLUSTER> *pgroup = NULL;
	VIRTUAL_MONOMER<CLUSTER> *vm = NULL, *vm1 = NULL;
	VECTOR3 r0, r1, dr;
	double R = 0;

	var = rdf_thread_var;
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

			for (idt = 0; idt < var->dt.ndt; idt++) {
				if (distrbt.dt[idt].sample(R, true)) distrbt.dt[idt].vnormal += 1; // total number of samples
			}
		}
	}
}

template <class MM, class CLUSTER> void single_vm_non_neighbor_rdf_func(NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, int nrdf) {
	int irdf = 0;
	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *var;

	for (irdf = 0; irdf < nrdf; irdf++) {
		var = rdf_thread_var + irdf;
		single_vm_non_neighbor_rdf_func1<MM, CLUSTER, float, long>(var, var->dt);
	}
}

// checking the neighbor vm distribution of multi-mm system, vmarray is not useful
template <class MD_CELL, class MM, class CLUSTER> void multi_vm_non_neighbor_in_sm_rdf_func(MD_CELL *md_cell, ARRAY< MM_VM<MM, CLUSTER> > *vmarray, NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, int nrdf) {
	single_vm_non_neighbor_rdf_func<MM, CLUSTER>(rdf_thread_var, nrdf);
}

template <class MD_CELL, class MM, class CLUSTER, class Tx, class Ty> void multi_vm_non_neighbor_rdf_func1(MD_CELL *md_cell, ARRAY< MM_VM<MM, CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, NDISTRIBT<Tx, Ty> &distrbt) {
	int idt = 0, im, im1, npt = 0, imol = 0, i;
	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *var;
	MM_VM<MM, CLUSTER> *mvm = NULL, *mvm1 = NULL;
	bool same_mol = false;

	cNeighbor< VIRTUAL_MONOMER<CLUSTER> > *neighbor = NULL;
	_GROUP<CLUSTER> *pgroup = NULL;
	VIRTUAL_MONOMER<CLUSTER> *vm = NULL, *vm1 = NULL;
	VECTOR3 r0, r1, dr;
	double R = 0;

	var = rdf_thread_var;
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
				if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD || ::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
					for (i = 0; i < 3; i++) {
						if (dr.v[i] > md_cell->h[i]) dr.v[i] -= md_cell->h[i] + md_cell->h[i];
						else if (dr.v[i] < -md_cell->h[i]) dr.v[i] += md_cell->h[i] + md_cell->h[i];
					}
				}
				scale_uv3(dr, dr, R) R = sqrt(R);

				for (idt = 0; idt < var->dt.ndt; idt++) {
					if (distrbt.dt[idt].sample(R, true)) distrbt.dt[idt].vnormal += 1; // total number of samples
				}
			}
		}
	}
}

template <class MD_CELL, class MM, class CLUSTER> void multi_vm_non_neighbor_rdf_func(MD_CELL *md_cell, ARRAY< MM_VM<MM, CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, int nrdf) {
	int irdf = 0;
	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *var;

	for (irdf = 0; irdf < nrdf; irdf++) {
		var = rdf_thread_var + irdf;
		multi_vm_non_neighbor_rdf_func1<MD_CELL, MM, CLUSTER, float, long>(md_cell, vmarray, var, var->dt);
	}
}

template <class MD_CELL, class MM, class CLUSTER, class Tx, class Ty> void multi_vm_nsm_rdf_func1(MD_CELL *md_cell, ARRAY< MM_VM<MM, CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, NDISTRIBT<Tx, Ty> &distrbt) {
	int idt = 0, im, im1, npt = 0, imol = 0, i;
	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *var;
	MM_VM<MM, CLUSTER> *mvm = NULL, *mvm1 = NULL;
	bool same_mol = false;

	_GROUP<CLUSTER> *pgroup = NULL;
	VIRTUAL_MONOMER<CLUSTER> *vm = NULL, *vm1 = NULL;
	VECTOR3 r0, r1, dr;
	double R = 0;

	var = rdf_thread_var;
	mvm = var->vm;
	for (im = 0; im < mvm->nVM; im++) {
		vm = mvm->vm + im; 
		if (var->cm_uid > 0 && vm->g.uid != var->cm_uid) continue;
		V32V3(vm->r, r0)

		for (imol = 0; imol < vmarray->n; imol++) {
			mvm1 = vmarray->m + imol;
			same_mol = (mvm == mvm1 ? true : false);
			if (same_mol) continue; // in different molecule
			for (im1 = 0; im1 < mvm1->nVM; im1++) {
				vm1 = mvm1->vm + im1;
				pgroup = &(vm1->g);
				if (var->nm_uid > 0 && pgroup->uid != var->nm_uid) continue;
				
				V32V3(vm1->r, r1)
				VECT3(r0, r1, dr)
				if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD || ::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
					for (i = 0; i < 3; i++) {
						if (dr.v[i] > md_cell->h[i]) dr.v[i] -= md_cell->h[i] + md_cell->h[i];
						else if (dr.v[i] < -md_cell->h[i]) dr.v[i] += md_cell->h[i] + md_cell->h[i];
					}
				}
				scale_uv3(dr, dr, R) R = sqrt(R);

				for (idt = 0; idt < var->dt.ndt; idt++) {
					if (distrbt.dt[idt].sample(R, true)) distrbt.dt[idt].vnormal += 1; // total number of samples
				}
			}
		}
	}
}

template <class MD_CELL, class MM, class CLUSTER> void multi_vm_nsm_rdf_func(MD_CELL *md_cell, ARRAY< MM_VM<MM, CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_thread_var, int nrdf) {
	int irdf = 0;
	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *var;

	for (irdf = 0; irdf < nrdf; irdf++) {
		var = rdf_thread_var + irdf;
		multi_vm_nsm_rdf_func1<MD_CELL, MM, CLUSTER, float, long>(md_cell, vmarray, var, var->dt);
	}
}

/************************************************************************************
	CALTULATING : <S_m(r)> -- RDF,  <S_m(r) * S_n(r)>
************************************************************************************/

class MULTI_NEIGHBOR_RDF_VAR {
public:
	// important: we suppose that all the distribution has the same x-axis, range and npts
	// we have to seperate neighbor-rdf and non-neighbor-rdf 
	//         because they will be initialized differently
	NEIGHBOR_RDF_VAR *rdf_var;  // neighbor rdf
	int n_rdf;
	NEIGHBOR_RDF_VAR *nn_rdf_var;  // non-neighbor rdf
	int n_nn_rdf;

	NDISTRIBT<float, float> crdf; // rdf correlation, <S(r) * S'(r)>

	MULTI_NEIGHBOR_RDF_VAR() {rdf_var = NULL; n_rdf = 0; nn_rdf_var = NULL; n_nn_rdf = 0;};
	void reset() {
		if (rdf_var == NULL) {delete[] rdf_var; rdf_var = NULL;} n_rdf = 0;
		if (nn_rdf_var == NULL) {delete[] nn_rdf_var; nn_rdf_var = NULL;} n_nn_rdf = 0;
		crdf.reset();
	};
	~MULTI_NEIGHBOR_RDF_VAR() {reset();};
	void set_neighbor_rdf(int nrdf) {
		if (rdf_var == NULL) {delete[] rdf_var; rdf_var = NULL;} n_rdf = 0;
		if (nrdf <= 0) return;
		rdf_var = new NEIGHBOR_RDF_VAR[nrdf];
		this->n_rdf = nrdf;
	};
	void set_non_neighbor_rdf(int nrdf) {
		if (nn_rdf_var == NULL) {delete[] nn_rdf_var; nn_rdf_var = NULL;} n_nn_rdf = 0;
		if (nrdf <= 0) return;
		nn_rdf_var = new NEIGHBOR_RDF_VAR[nrdf];
		this->n_nn_rdf = nrdf;
	};
	void setup_correlated_distribts() {
		int n = n_rdf + n_nn_rdf;
		if (n <= 0) return;
		crdf.set(n * n);
		for (n = 0; n < crdf.ndt; n++) {
			if (rdf_var != NULL) crdf.dt[n].set_distribt_x(rdf_var[0].dt.dt[0].xf, rdf_var[0].dt.dt[0].xt, rdf_var[0].dt.dt[0].x.N - 1);
			else if (nn_rdf_var != NULL) crdf.dt[n].set_distribt_x(nn_rdf_var[0].dt.dt[0].xf, nn_rdf_var[0].dt.dt[0].xt, nn_rdf_var[0].dt.dt[0].x.N - 1);;
		}
	};
};

template <class MM, class CLUSTER> class NEIGHBOR_RDF_CORRELATION_THREAD_VAR {
public:
	int cm_uid, nm_uid; // center monomer's uid, neibhbor-monomer's uid
	int nth_neighbor;

	MM_VM<MM, CLUSTER> *vm;

	// important: we suppose that all the distribution has the same x-axis, range and npts
	// we have to seperate neighbor-rdf and non-neighbor-rdf 
	//         because they will be initialized differently
	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *rdf_var; // neighbor rdf
	int n_rdf;

	NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> *nn_rdf_var; // non-neighbor rdf
	int n_nn_rdf;

	NDISTRIBT<float, long> crdf; // rdf correlation, <S(r) * S'(r)>

	NEIGHBOR_RDF_CORRELATION_THREAD_VAR() {
		vm = NULL;
		rdf_var = NULL; n_rdf = 0;
		nn_rdf_var = NULL; n_nn_rdf = 0;
	};
	void reset() {
		if (rdf_var == NULL) {delete[] rdf_var; rdf_var = NULL;} n_rdf = 0;
		if (nn_rdf_var == NULL) {delete[] nn_rdf_var; nn_rdf_var = NULL;} n_nn_rdf = 0;
		crdf.reset();
	};
	~NEIGHBOR_RDF_CORRELATION_THREAD_VAR() {reset();};
	void set_neighbor_rdf(int nrdf) {
		if (rdf_var == NULL) {delete[] rdf_var; rdf_var = NULL;} n_rdf = 0;
		if (nrdf <= 0) return;
		rdf_var = new NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER>[nrdf];
		this->n_rdf = nrdf;
	};
	void set_non_neighbor_rdf(int nrdf) {
		if (nn_rdf_var == NULL) {delete[] nn_rdf_var; nn_rdf_var = NULL;} n_nn_rdf = 0;
		if (nrdf <= 0) return;
		nn_rdf_var = new NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER>[nrdf];
		this->n_nn_rdf = nrdf;
	};
	void setup_correlated_distribts() {
		int n = n_rdf + n_nn_rdf;
		if (n <= 0) return;
		crdf.set(n * n);
		for (n = 0; n < crdf.ndt; n++) {
			if (rdf_var != NULL) crdf.dt[n].set_distribt_x(rdf_var[0].dt.dt[0].xf, rdf_var[0].dt.dt[0].xt, rdf_var[0].dt.dt[0].x.N - 1);
			else if (nn_rdf_var != NULL) crdf.dt[n].set_distribt_x(nn_rdf_var[0].dt.dt[0].xf, nn_rdf_var[0].dt.dt[0].xt, nn_rdf_var[0].dt.dt[0].x.N - 1);;
		}
	};
};

template <class MM, class CLUSTER> void init_vm_rdf_correlation_thread_var(MULTI_NEIGHBOR_RDF_VAR* mrdf_var, NEIGHBOR_RDF_CORRELATION_THREAD_VAR<MM, CLUSTER>* crdf_thread_var) {
	int i = 0;
	crdf_thread_var->reset();
	crdf_thread_var->set_neighbor_rdf(mrdf_var->n_rdf);
	for (i = 0; i < mrdf_var->n_rdf; i++) {
		setup_ndt(crdf_thread_var->rdf_var[i].dt, mrdf_var->rdf_var[i].dt);
		crdf_thread_var->rdf_var[i].cm_uid = mrdf_var->rdf_var[i].cm_uid;
		crdf_thread_var->rdf_var[i].nm_uid = mrdf_var->rdf_var[i].nm_uid;
		crdf_thread_var->rdf_var[i].vm = crdf_thread_var->vm;
		crdf_thread_var->rdf_var[i].nth_neighbor = mrdf_var->rdf_var[i].nth_neighbor;
		crdf_thread_var->rdf_var[i].neighbor.set_elements(crdf_thread_var->rdf_var[i].vm->nVM, mrdf_var->rdf_var[i].nth_neighbor);
		init_virtual_monomer_neighbor<MM, CLUSTER>(crdf_thread_var->vm->vm, crdf_thread_var->vm->nVM, crdf_thread_var->rdf_var[i].neighbor.neighbor);
	}
	crdf_thread_var->set_non_neighbor_rdf(mrdf_var->n_nn_rdf);
	for (i = 0; i < mrdf_var->n_nn_rdf; i++) {
		setup_ndt(crdf_thread_var->nn_rdf_var[i].dt, mrdf_var->nn_rdf_var[i].dt);
		crdf_thread_var->nn_rdf_var[i].cm_uid = mrdf_var->nn_rdf_var[i].cm_uid;
		crdf_thread_var->nn_rdf_var[i].nm_uid = mrdf_var->nn_rdf_var[i].nm_uid;
		crdf_thread_var->nn_rdf_var[i].vm = crdf_thread_var->vm;
		crdf_thread_var->nn_rdf_var[i].nth_neighbor = mrdf_var->nn_rdf_var[i].nth_neighbor;
		crdf_thread_var->nn_rdf_var[i].neighbor.set_elements(crdf_thread_var->nn_rdf_var[i].vm->nVM, mrdf_var->nn_rdf_var[i].nth_neighbor);
		init_virtual_monomer_all_neighbors<MM, CLUSTER>(crdf_thread_var->vm->vm, crdf_thread_var->vm->nVM, crdf_thread_var->nn_rdf_var[i].neighbor.neighbor);
	}
	crdf_thread_var->setup_correlated_distribts();
};

void single_mm_cal_Position(MM_VM<MMOLECULE, CLUSTER>& mvm);
void multi_mm_cal_Position_NonCell(ARRAY< MM_VM<MMOLECULE, CLUSTER> >& vmarray);
void multi_mm_cal_Position_InCell(ARRAY< MM_VM<MMOLECULE, CLUSTER> >& vmarray);

bool read_vm_var(char *parfname, VM_VAR<MMOLECULE> &vm_var);

bool read_neighbor_rdf_pars(char *parfname, NEIGHBOR_RDF_VAR** rdf_var, int &nrdf);
bool read_non_neighbor_rdf_pars(char *parfname, NEIGHBOR_RDF_VAR** rdf_var, int &nrdf);

void raw_accumulate_rdf(NEIGHBOR_RDF_VAR *rdf, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER>* nrdf);

void normalize_rdf(NEIGHBOR_RDF_VAR *rdf);
void normalize_mrdf(MULTI_NEIGHBOR_RDF_VAR *mrdf);

void save_rdf(NEIGHBOR_RDF_VAR *rdf_var, char *title, char *msg_title);
void save_mrdf(MULTI_NEIGHBOR_RDF_VAR *mrdf_var, char *title, char *msg_title);

}  // end of namespace _MM_rdf_





//***************************************************************************
//        CORRELATION OF ATOMIC COORDINATES
//***************************************************************************
namespace coordinate_correlation {
using namespace _distribution_;

class VATOM {
public:
	int uid;
	VECTOR3 *r;
	float w; // weight 
};

class VATOM1 {
public:
	int uid;
	SVECTOR<double, 3> r;
	float w;
};

template <class VATOM> class VATOM_CELL {
public:
	ARRAY<VATOM> va;
	float xd[3];
	bool period;

	~VATOM_CELL() {
		va.release();
	};
};

class CORREL_DISTRIBT {
public:
	int cuid, ruid;
	DISTRIBT<float, double> dt;
	int ncatom;
	bool bRadialDistribt;
	CORREL_DISTRIBT() {
		ncatom = 0; cuid = 0; ruid = 0;
		bRadialDistribt = true;
	};
	void operator = (CORREL_DISTRIBT &c1) {
		this->cuid = c1.cuid;
		this->ruid = c1.ruid;
		this->ncatom = c1.ncatom;
		this->bRadialDistribt = c1.bRadialDistribt;
	};
};


class CORREL_2D_DISTRIBT : public CORREL_DISTRIBT {
public:
	class Zregion {
	public:
		float zf, zt;
	};
	Zregion rg1, rg2;
	CORREL_2D_DISTRIBT() {};
	~CORREL_2D_DISTRIBT() {};
	bool inZrg1(double z) {
		if (z < rg1.zf || z > rg1.zt) return false;
		else return true;
	};
	bool inZrg2(double z) {
		if (z < rg2.zf || z > rg2.zt) return false;
		else return true;
	};
	void operator = (CORREL_2D_DISTRIBT &c1) {
		this->cuid = c1.cuid;
		this->ruid = c1.ruid;
		this->ncatom = c1.ncatom;
		this->bRadialDistribt = c1.bRadialDistribt;
		this->rg1.zf = c1.rg1.zf;
		this->rg1.zt = c1.rg1.zt;
		this->rg2.zf = c1.rg2.zf;
		this->rg2.zt = c1.rg2.zt;
	};
};


template <class VCELL, class COR_DISTRIBT> class CORRELATION_THREAD_VAR {
public:
	VCELL cell;
	ARRAY< COR_DISTRIBT > dt;
	int nconf;
};

template <class CORREL_DISTRIBT> class CORREL_THREAD : public _THREAD_VARS {
public:
	CORRELATION_THREAD_VAR< VATOM_CELL<VATOM1>, CORREL_DISTRIBT > *var;
	void *func;
};

#if _SYS_ == _WINDOWS_SYS_
template <class CORREL_DISTRIBT> int correlation_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class CORREL_DISTRIBT> void* correlation_thread(void *void_par) { // thread-function
#endif
	typedef void (*FUNC)(VATOM_CELL<VATOM1>&, ARRAY< CORREL_DISTRIBT >&);
	CORREL_THREAD<CORREL_DISTRIBT> *par = (CORREL_THREAD<CORREL_DISTRIBT>*)void_par;
	((FUNC)(par->func))(par->var->cell, par->var->dt);
	par->var->nconf++;
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

/*
#if _SYS_ == _WINDOWS_SYS_
extern int correlation_2d_thread(LPVOID *void_par); // thread-function
#elif _SYS_ == _LINUX_SYS_
extern void* correlation_2d_thread(void *void_par); // thread-function
#endif
*/

template <class MD_CELL, class VCELL, class CORREL_DISTRIBT> 
bool correlation_single_md_multithread(MD_CELL &md_cell, char *mdpar_fname, int nf, int nt, int nstep, 
						   void *read_molstruct, void *update_cell, void* vcell_constructor, void *update_vcell, 
						   CORRELATION_THREAD_VAR<VCELL, CORREL_DISTRIBT> *corr_thread_var, int nThreads, void *correlation_func) 
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

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int iconf = 0, nconf = 0, iThread, nthread = 0;

	typedef void (*VCELL_CONSTRUCTOR)(MD_CELL&, VCELL&);

	ARRAY<CORREL_THREAD<CORREL_DISTRIBT> > thread_var; thread_var.SetArray(nThreads);
	for (iThread = 0; iThread < nThreads; iThread++) {
		((VCELL_CONSTRUCTOR)vcell_constructor)(md_cell, corr_thread_var[iThread].cell);

		thread_var.m[iThread].set_pars(iThread);
		thread_var.m[iThread].var = corr_thread_var + iThread;
		thread_var.m[iThread].func = correlation_func;
	}

	typedef void (*UPDATE_VCELL)(MD_CELL&, VCELL&);
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
			if (update_vcell != NULL) ((UPDATE_VCELL)update_vcell)(md_cell, corr_thread_var[iThread].cell);

		#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(thread_var.m[iThread].hMMThread);
			AfxBeginThread((AFX_THREADPROC)coordinate_correlation::correlation_thread<CORREL_DISTRIBT>, (LPVOID)(thread_var.m + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
		#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(thread_var.m[iThread].thread), NULL, &correlation_thread<CORREL_DISTRIBT>, (void *)(thread_var.m + iThread));
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
	
	return true;
}

/*********************************************************************************
	MULTI-THREAD FUNCTION TO DO RDF-ANALYSIS OF MULTI-MOLECULE WITH SEVERAL MD SIMULATION
**********************************************************************************/
template <class MD_CELL, class VATOM_CELL, class CORREL_DISTRIBT> 
bool correlation_multi_md_multithread(char *parfname, MD_CELL &md_cell, ARRAY<CORREL_DISTRIBT> &dt,
						    void *read_molstruct, void *update_cell, void* vcell_constructor, void *update_vcell, 
							void *correlation_func) 
{
	CORRELATION_THREAD_VAR< VATOM_CELL, CORREL_DISTRIBT > corr_thread_var[MAX_THREADS];
	int i, j, k;
	for (i = 0; i < dt.n; i++) {
		dt.m[i].dt.reset_distribt(); // reset the distribution
	}
	for (i = 0; i < MAX_THREADS; i++) {
		corr_thread_var[i].dt.SetArray(dt.n);
		corr_thread_var[i].nconf = 0;
		for (j = 0; j < dt.n; j++) {
			//corr_thread_var[i].dt.m[j].cuid = dt.m[j].cuid;
			//corr_thread_var[i].dt.m[j].ruid = dt.m[j].ruid;
			//corr_thread_var[i].dt.m[j].ncatom = 0;
			//corr_thread_var[i].dt.m[j].bRadialDistribt = dt.m[j].bRadialDistribt;
			corr_thread_var[i].dt.m[j] = dt.m[j];

			corr_thread_var[i].dt.m[j].dt.set_distribt_x(dt.m[j].dt.xf, dt.m[j].dt.xt, dt.m[j].dt.x.N - 1);
			corr_thread_var[i].dt.m[j].dt.dx = dt.m[j].dt.dx;
		}
	}

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

		ARRAY<int> mol_id;

		if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD || MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			status = correlation_single_md_multithread<MD_CELL, VATOM_CELL, CORREL_DISTRIBT >(md_cell, mdpar_fname, nf, nt, nstep,  
				read_molstruct, update_cell, vcell_constructor, update_vcell, corr_thread_var, MAX_THREADS, correlation_func);
			if (!status) return false;
		}
		else {
			sprintf(msg, "MD in %s, is not multi-mm", mdpar_fname); show_msg(msg, true);
		}

		n++;
	}
	in_cfg.close();

	for (i = 0; i < dt.n; i++) {
		for (j = 0; j < MAX_THREADS; j++) {
			raw_distribt_combine<float, double, float, double>(dt.m[i].dt, corr_thread_var[j].dt.m[i].dt);
			dt.m[i].dt.eb_raw();
		}
	}
	nconfs = 0;
	for (i = 0; i < MAX_THREADS; i++) nconfs += corr_thread_var[i].nconf;
	sprintf(msg, "totally %d configs. are analyzed", nconfs); show_log(msg, true);

	for (i = 0; i < dt.n; i++) {
		dt.m[i].ncatom = 0;
		for (j = 0; j < MAX_THREADS; j++) dt.m[i].ncatom += corr_thread_var[j].dt.m[i].ncatom;
	}
	double v = 0;
	ofstream *out = NULL;
	for (i = 0; i < dt.n; i++) {
		dt.m[i].dt.vnormal = dt.m[i].ncatom;
		dt.m[i].dt.y /= dt.m[i].ncatom;  // normalize the number of center atom accumulated
		dt.m[i].dt.eb /= dt.m[i].ncatom;
		if (dt.m[i].bRadialDistribt) {
			for (j = 0; j < dt.m[i].dt.x.N; j++) {
				if (dt.m[i].dt.x.v[j] == 0) continue;
				v = 4 * PI * (dt.m[i].dt.x.v[j] * dt.m[i].dt.x.v[j]) * dt.m[i].dt.dx;
				dt.m[i].dt.y.v[j] /= v; dt.m[i].dt.eb.v[j] /= v;
			}
		}
		sprintf(ofname, "%s-corr-%d.dat", title, i);
		out = new ofstream;
		out->open(ofname);
		if (out->is_open()) {
			dt.m[i].dt.save_x_y_eb(*out); out->close();
		}
		delete out; out = NULL;
	}

	return true;
};



// single-thread function

template <class MD_CELL, class VCELL, class CORREL_DISTRIBT> 
bool correlation_single_md_singlethread(MD_CELL &md_cell, char *mdpar_fname, int nf, int nt, int nstep, 
						   void *read_molstruct, void *update_cell, void* vcell_constructor, void *update_vcell, 
						   CORRELATION_THREAD_VAR<VCELL, CORREL_DISTRIBT> &corr_thread_var, void *correlation_func) 
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

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int iconf = 0, nconf = 0;

	typedef void (*VCELL_CONSTRUCTOR)(MD_CELL&, VCELL&);

	((VCELL_CONSTRUCTOR)vcell_constructor)(md_cell, corr_thread_var.cell);

	typedef void (*UPDATE_VCELL)(MD_CELL&, VCELL&);
	typedef bool (*READ_Moltruct)(int);
	typedef void (*UPDATE_CELL)(MD_CELL&);
	typedef void (*UPDATE_VCELL)(MD_CELL&, VCELL&);

	typedef void (*CORR_FUNC) (VCELL&, ARRAY< CORREL_DISTRIBT >&);

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
		if (update_vcell != NULL) ((UPDATE_VCELL)update_vcell)(md_cell, corr_thread_var.cell);

		((CORR_FUNC)(correlation_func))(corr_thread_var.cell, corr_thread_var.dt);
		corr_thread_var.nconf++;

		iconf += nstep; nconf++;

		show_log(", analyzed.", true);
	}

	clear_read_mdproc();
	
	return true;
}

/*********************************************************************************
	SINGLE-THREAD FUNCTION TO DO RDF-ANALYSIS OF MULTI-MOLECULE WITH SEVERAL MD SIMULATION
**********************************************************************************/
template <class MD_CELL, class VATOM_CELL, class CORREL_DISTRIBT> 
bool correlation_multi_md_singlethread(char *parfname, MD_CELL &md_cell, ARRAY<CORREL_DISTRIBT> &dt,
						    void *read_molstruct, void *update_cell, void* vcell_constructor, void *update_vcell, 
							void *correlation_func) 
{
	CORRELATION_THREAD_VAR<VATOM_CELL, CORREL_DISTRIBT> corr_thread_var;
	int i, j;
	for (i = 0; i < dt.n; i++) {
		dt.m[i].dt.reset_distribt(); // reset the distribution
	}

	corr_thread_var.dt.SetArray(dt.n);
	corr_thread_var.nconf = 0;
	for (j = 0; j < dt.n; j++) {
		//corr_thread_var.dt.m[j].cuid = dt.m[j].cuid;
		//corr_thread_var.dt.m[j].ruid = dt.m[j].ruid;
		//corr_thread_var.dt.m[j].ncatom = 0;
		//corr_thread_var.dt.m[j].bRadialDistribt = dt.m[j].bRadialDistribt;
		corr_thread_var.dt.m[j] = dt.m[j];

		corr_thread_var.dt.m[j].dt.set_distribt_x(dt.m[j].dt.xf, dt.m[j].dt.xt, dt.m[j].dt.x.N - 1);
		corr_thread_var.dt.m[j].dt.dx = dt.m[j].dt.dx;
	}

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

		if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD || MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			status = correlation_single_md_singlethread<MD_CELL, VATOM_CELL, CORREL_DISTRIBT >(md_cell, mdpar_fname, nf, nt, nstep,  
				read_molstruct, update_cell, vcell_constructor, update_vcell, corr_thread_var, correlation_func);
			if (!status) return false;
		}
		else {
			sprintf(msg, "MD in %s, is not multi-mm", mdpar_fname); show_msg(msg, true);
		}

		n++;
	}
	in_cfg.close();

	for (i = 0; i < dt.n; i++) {
		raw_distribt_combine<float, double, float, double>(dt.m[i].dt, corr_thread_var.dt.m[i].dt);
		dt.m[i].ncatom = corr_thread_var.dt.m[i].ncatom;
	}
	nconfs = corr_thread_var.nconf;
	sprintf(msg, "totally %d configs. are analyzed", nconfs); show_log(msg, true);

	double v;
	ofstream *out = NULL;
	for (i = 0; i < dt.n; i++) {
		dt.m[i].dt.vnormal = dt.m[i].ncatom;
		dt.m[i].dt.eb_raw();
		dt.m[i].dt.y /= dt.m[i].ncatom;  // normalize the number of center atom accumulated
		dt.m[i].dt.eb /= dt.m[i].ncatom;  // normalize the number of center atom accumulated
		if (dt.m[i].bRadialDistribt) {
			for (j = 0; j < dt.m[i].dt.x.N; j++) {
				if (dt.m[i].dt.x.v[j] == 0) continue;
				v = 4 * PI * (dt.m[i].dt.x.v[j] * dt.m[i].dt.x.v[j]) * dt.m[i].dt.dx;
				dt.m[i].dt.y.v[j] /= v; dt.m[i].dt.eb.v[j] /= v;
			}
		}
		sprintf(ofname, "%s-corr-%d.dat", title, i);
		out = new ofstream;
		out->open(ofname);
		if (out->is_open()) {
			dt.m[i].dt.save_x_y_eb(*out); out->close();
		}
		delete out; out = NULL;
	}

	return true;
};


#define WEIGHT_NUMBER    0
#define WEIGHT_MASS      1
#define WEIGHT_INDX      2
}  // end of namesapce coordinate_correlation


