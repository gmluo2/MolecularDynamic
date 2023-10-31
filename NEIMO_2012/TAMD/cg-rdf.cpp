#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "ranlib.h"

#include "def.h"
#include "show.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"

#include "ZMatrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
//#include "cluster.h"

#include "var.h"

#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"
#include "Interaction.h"
#include "Interact1.h"
#include "md_cell.h"
#include "cg-md-cell.h"

#define _READ_MM_MD_ 1
#include "read.h"

#include "ReadMDProc.h"
#include "mol-cell.h"
#include "show_md.h"
#if _SYS_ == _WINDOWS_SYS_ 
#include "..\export1.h"
#elif _SYS_ == _LINUX_SYS_
#include "export1.h"
#endif

extern "C"  bool construct_macromolecule(char *fname);
extern "C" bool init_read_mdproc(char *mdpar_fname);
extern "C" bool serial_md_proc(long isave);
extern "C" void clear_read_mdproc();
extern "C" void set_md_mode(int mode);
extern "C" void set_job(int job);
extern char struct_fname[256];
extern MMOL_MD_CELL mm_md_cell;
extern int MD_MODE;
extern READ_MD_PROC_FILE rmdproc;

extern bool get_monomer_uid(char *monomer, int &uid);

extern bool read_given_id(char *fname, char* title, ARRAY<short>& id_cluster);


#include "distribution.h"
#include "rdf-distribt.h"

#include "cg-cluster.h"
#include "cg-rdf.h"

using namespace _coarse_grain_;
extern CGMM_MD_CELL cgmm_md_cell;

namespace _CGMM_rdf_ {
using namespace _coarse_grain_;

using namespace _distribution_;
using namespace _MM_rdf_;

// cneighbor is an array, with m->nCluster of elements
void init_cluster_neighbor(CG_MMOLECULE *m, cNeighbor<CG_CLUSTER> *cneighbor) {
	CG_CLUSTER *pc = NULL;
	CHAIN<CG_CLUSTER> *ch = NULL;
	for (int nc = 0; nc < m->nCluster; nc++) {
		pc = m->cluster + nc;
		grandchildren_chain<CG_CLUSTER>(pc, cneighbor[nc].nth, &ch);
		chain_array<CG_CLUSTER>(&ch, cneighbor[nc].neighbor);
		release_chain<CG_CLUSTER>(&ch, false);
	}
}

#define cgcluster_uid(cgcluster) cgcluster.atom[0].par->uid
#define pcgcluster_uid(cgcluster) cgcluster->atom[0].par->uid

bool construct_monomer_single_chain_homopolymer(CG_MMOLECULE *m, int m_uid, int nm_cg, int vm_uid, VIRTUAL_MONOMER<CG_CLUSTER> **cgmonomer, int& ncgMonomers) {
	char msg[256] = "\0";
	CHAIN<CG_CLUSTER> *ch = NULL, *tail = NULL, *pch = NULL;
	int nc = 0, nm = 0, im = 0, icg = 0;
	int nms = 0;
	if (*cgmonomer != NULL) {delete[] *cgmonomer; *cgmonomer = NULL;}
	for (nm = 0; nm < m->nCluster; nm++) {
		if (cgcluster_uid(m->cluster[nm]) == m_uid) nms++;
	}
	ncgMonomers = nms / nm_cg;
	if (ncgMonomers == 0) {
		sprintf(msg, "macromolecule has only %d monomers, while each virtual_monomer needs %d monomers", nms, nm_cg);
		show_log(msg, true); return false;
	}
	else if (nms > ncgMonomers * nm_cg) {
		sprintf(msg, "macromolecule has only %d monomers, each virtual_monomer needs %d monomers. So, some monomer will be not used.", nms, nm_cg);
		show_log(msg, true);
	}
	*cgmonomer = new VIRTUAL_MONOMER<CG_CLUSTER>[ncgMonomers];

	for (nm = 0; nm < m->nCluster; nm++) {
		if (cgcluster_uid(m->cluster[nm]) != m_uid) continue;
		pch = new CHAIN<CG_CLUSTER>; pch->p = m->cluster + nm;
		if (ch == NULL) ch = pch;
		else tail->next = pch;
		tail = pch;
		
		im += 1;
		if (im == nm_cg) {
			construct_virtual_monomer<CG_CLUSTER>(ch, *cgmonomer + icg);
			virtual_monomer_setup_hinges<CG_CLUSTER>(*cgmonomer + icg);
			(*cgmonomer)[icg].g.uid = vm_uid;
			im = 0; icg += 1;
			release_chain<CG_CLUSTER>(&ch, false); tail = NULL;
		}
	}

	connect_virtual_monomers<CG_CLUSTER>(*cgmonomer, ncgMonomers);
	return true;
}

void init_cgneighbor_rdf_thread_var(NEIGHBOR_RDF_VAR* rdf_var, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER>* rdf_thread_var) {
	return init_vm_neighbor_rdf_thread_var<CG_MMOLECULE, CG_CLUSTER>(rdf_var, rdf_thread_var);
}

void init_all_cgneighbors_rdf_thread_var(NEIGHBOR_RDF_VAR* rdf_var, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER>* rdf_thread_var) {
	return init_vm_all_neighbors_rdf_thread_var<CG_MMOLECULE, CG_CLUSTER>(rdf_var, rdf_thread_var);
}

bool construct_monomer_single_chain_block_cgcopolymer(VM_VAR<CG_MMOLECULE>* vm_var, VIRTUAL_MONOMER<CG_CLUSTER> **cgmonomer, int& ncgMonomers) {
	char msg[256] = "\0";
	CHAIN<CG_CLUSTER> *ch = NULL, *tail = NULL, *pch = NULL;
	int nc = 0, nm = 0, im = 0, icg = 0, iBlock = 0;
	int *nms = new int[vm_var->nBlocks];
	int *ncgm = new int[vm_var->nBlocks];
	if (*cgmonomer != NULL) {delete[] *cgmonomer; *cgmonomer = NULL;}
	for (iBlock = 0; iBlock < vm_var->nBlocks; iBlock++) {nms[iBlock] = 0; ncgm[iBlock] = 0;}

#define RELEASE if (nms != NULL) {delete[] nms; nms = NULL;} \
	if (ncgm != NULL) {delete[] ncgm; ncgm = NULL;}

	CG_MMOLECULE *m = vm_var->m;
	for (nm = 0; nm < m->nCluster; nm++) {
		for (iBlock = 0; iBlock < vm_var->nBlocks; iBlock++) {
			if (cgcluster_uid(m->cluster[nm]) == vm_var->muid[iBlock]) {nms[iBlock]++; break;}
		}
	}
	ncgMonomers = 0;
	for (iBlock = 0; iBlock < vm_var->nBlocks; iBlock++) {
		ncgm[iBlock] = nms[iBlock] / vm_var->nm_vm[iBlock];
		ncgMonomers += ncgm[iBlock];
		if (ncgm[iBlock] == 0) {
			sprintf(msg, "macromolecule has only %d monomers of block %d, while each virtual_monomer needs %d monomers", nms[iBlock], iBlock, vm_var->nm_vm[iBlock]);
			show_log(msg, true); RELEASE return false;
		}
		else if (nms[iBlock] > ncgm[iBlock] * vm_var->nm_vm[iBlock]) {
			sprintf(msg, "macromolecule has only %d monomers of block %d, each virtual_monomer needs %d monomers. So, some monomer will be not used.", nms[iBlock], iBlock, vm_var->nm_vm[iBlock]);
			show_log(msg, true);
		}
	}
	RELEASE
	if (ncgMonomers == 0) return false;
	*cgmonomer = new VIRTUAL_MONOMER<CG_CLUSTER>[ncgMonomers];

	bool status = false;
	int muid = 0;
	for (nm = 0; nm < m->nCluster; nm++) {
		status = false;
		for (iBlock = 0; iBlock < vm_var->nBlocks; iBlock++) {
			if (cgcluster_uid(m->cluster[nm]) == vm_var->muid[iBlock]) {status = true; break;}
		}
		if (!status) continue;
		if (ch == NULL) {
			pch = new CHAIN<CG_CLUSTER>; pch->p = m->cluster + nm;
			ch = pch;
			muid = cgcluster_uid(m->cluster[nm]);
			tail = pch;
		}
		else if (cgcluster_uid(m->cluster[nm]) != muid) { // this is a different type of monomer
			// this virtual_monomer can not be constructed
			release_chain<CG_CLUSTER>(&ch, false); tail = NULL;
			status = false;
		}
		else {
			pch = new CHAIN<CG_CLUSTER>; pch->p = m->cluster + nm;
			tail->next = pch;
			tail = pch;
		}

		if (!status) continue; // this virtual_monomer can not be constructed
		
		im += 1;
		if (im == vm_var->nm_vm[iBlock]) {
			construct_virtual_monomer<CG_CLUSTER>(ch, *cgmonomer + icg);
			virtual_monomer_setup_hinges<CG_CLUSTER>(*cgmonomer + icg);
			(*cgmonomer)[icg].g.uid = vm_var->vmuid[iBlock];
			im = 0; icg += 1;
			release_chain<CG_CLUSTER>(&ch, false); tail = NULL;
		}
	}

	connect_virtual_monomers<CG_CLUSTER>(*cgmonomer, ncgMonomers);
	return true;
#undef RELEASE
}

void CGcalPosition(VIRTUAL_MONOMER<CG_CLUSTER> *vm, int method) {
	int nc;
	float M = 0, m = 0;
	if (method == 0 || method == 1) { // mass center
		V3zero(vm->r)
		if (vm->g.nUnit == 1) {
			V32V3(vm->g.u[0]->atom[0].r, vm->r)
		}
		else {
			for (nc = 0; nc < vm->g.nUnit; nc++) {
				m = *(vm->g.u[nc]->M);
				vm->r.v[0] += vm->g.u[nc]->r->v[0] * m; 
				vm->r.v[1] += vm->g.u[nc]->r->v[1] * m; 
				vm->r.v[2] += vm->g.u[nc]->r->v[2] * m;
				M += m;
			}
			vm->r.v[0] /= M; vm->r.v[1] /= M; vm->r.v[2] /= M;
		}
	}
	if (method == 1) {
		// use the first cluser's dr_cmm for the CMM-shift
		vm->r.v[0] += vm->g.u[0]->dr_cmm.v[0];
		vm->r.v[1] += vm->g.u[0]->dr_cmm.v[1];
		vm->r.v[2] += vm->g.u[0]->dr_cmm.v[2];
	}
}

void single_cgmm_cal_Position(MM_VM<CG_MMOLECULE, CG_CLUSTER>& mvm) {
	for (int i = 0; i < mvm.nVM; i++) CGcalPosition(mvm.vm + i, 0);
}

void multi_cgmm_cal_Position_NonCell(ARRAY< MM_VM<CG_MMOLECULE, CG_CLUSTER> >& vmarray) { // position not in cell
	int i, j;
	for (i = 0; i < vmarray.n; i++) {
		for (j = 0; j < vmarray.m[i].nVM; j++) CGcalPosition(vmarray.m[i].vm + j, 0);
	}
}

void multi_cgmm_cal_Position_InCell(ARRAY< MM_VM<CG_MMOLECULE, CG_CLUSTER> >& vmarray) { // position in cell
	int i, j;
	for (i = 0; i < vmarray.n; i++) {
		for (j = 0; j < vmarray.m[i].nVM; j++) CGcalPosition(vmarray.m[i].vm + j, 1);
	}
}

void single_cgmm_neighbor_rdf_func(NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> *rdf_thread_var, int nrdf) {
	return single_vm_neighbor_rdf_func<CG_MMOLECULE, CG_CLUSTER>(rdf_thread_var, nrdf);
}

void multi_cgmm_neighbor_rdf_func(CGMM_MD_CELL *md_cell, ARRAY< MM_VM<CG_MMOLECULE, CG_CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> *rdf_thread_var, int nrdf) {
	return multi_vm_neighbor_rdf_func<CGMM_MD_CELL, CG_MMOLECULE, CG_CLUSTER>(md_cell, vmarray, rdf_thread_var, nrdf);
}

void single_cgmm_non_neighbor_rdf_func(NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> *rdf_thread_var, int nrdf) {
	return single_vm_non_neighbor_rdf_func<CG_MMOLECULE, CG_CLUSTER>(rdf_thread_var, nrdf);
}

void multi_cgmm_non_neighbor_rdf_func(CGMM_MD_CELL *md_cell, ARRAY< MM_VM<CG_MMOLECULE, CG_CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> *rdf_thread_var, int nrdf) {
	return multi_vm_non_neighbor_rdf_func<CGMM_MD_CELL, CG_MMOLECULE, CG_CLUSTER>(md_cell, vmarray, rdf_thread_var, nrdf);
}

void multi_cgmm_non_neighbor_in_sm_rdf_func(CGMM_MD_CELL *md_cell, ARRAY< MM_VM<CG_MMOLECULE, CG_CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> *rdf_thread_var, int nrdf) {
	return multi_vm_non_neighbor_in_sm_rdf_func<CGMM_MD_CELL, CG_MMOLECULE, CG_CLUSTER>(md_cell, vmarray, rdf_thread_var, nrdf);
}

void multi_cgmm_nsm_rdf_func(CGMM_MD_CELL *md_cell, ARRAY< MM_VM<CG_MMOLECULE, CG_CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> *rdf_thread_var, int nrdf) {
	return multi_vm_nsm_rdf_func<CGMM_MD_CELL, CG_MMOLECULE, CG_CLUSTER>(md_cell, vmarray, rdf_thread_var, nrdf);
}

void raw_accumulate_rdf(NEIGHBOR_RDF_VAR *rdf, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER>* nrdf) {
	int irdf, i;
	for (irdf = 0; irdf < rdf->dt.ndt; irdf++) {
		raw_distribt_combine<float, float, float, long>(rdf->dt.dt[irdf], nrdf->dt.dt[irdf]);
		rdf->dt.dt[irdf].vnormal += nrdf->dt.dt[irdf].vnormal;
	}
}

bool read_vm_var(char *parfname, VM_VAR<CG_MMOLECULE> &vm_var) {
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	int i, nBlockCopolymer = 0;
	if (!search(in, "[BLOCKS]", buffer)) {
		sprintf(msg, "[BLOCKS] is not defined in %s", parfname); show_log(msg, true); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nBlockCopolymer) != 1) {
		sprintf(msg, "blocks of copolymer is not give after [BLOCKS] in %s", parfname); show_log(msg, true); in.close(); return false;
	}
	else if (nBlockCopolymer < 1) {
		sprintf(msg, "blocks of copolymer has to be more than 1", parfname); show_log(msg, true); in.close(); return false;
	}
	vm_var.nBlocks = nBlockCopolymer;

	rewind(in);
	int ncgm = 0;
	if (!search(in, "[MONOMER-GROUP]", buffer)) {sprintf(msg, "can not find [MONOMER-GROUP] in %s", parfname); show_msg(msg, true); in.close(); return false;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &ncgm) != 1 || ncgm <= 0) {
		sprintf(msg, "number of coarse-grained  monomer has to be > 0, after [MONOMER-GROUP]"); show_msg(msg, true);
		in.close(); return true;
	}
	vm_var.set_vm(ncgm);
	char monomer[256] = "\0";
	for (i = 0; i < vm_var.ncgm; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s %d %d", monomer, vm_var.nm_vm + i, vm_var.vmuid + i) != 3) {
			sprintf(msg, "number of coarse-grained  monomer has to be > 0, after [MONOMER-GROUP]"); show_msg(msg, true);
			in.close(); return false;
		}
		else if (!get_cgatom_uid(monomer, vm_var.muid[i])) {
			in.close(); return false;
		}
		vm_var.vmuid[i] += _USER_MONOMER_UID;
	}

	in.close(); return true;
}

} // end of namespace _CGMM_rdf_

extern bool ReadCGMMStruct(ifstream &in, CG_MMOLECULE &mm);
extern bool serial_multi_mol_md_proc(long isave);
extern bool serial_multi_cgmm_md_proc(long isave);

extern "C" bool single_cgmm_neighbor_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _CGMM_rdf_;
	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<CG_MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!_MM_rdf_::read_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_single_mm_multi_md<CG_MMOLECULE, CG_CLUSTER, VM_VAR<CG_MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> >(
		parfname, (void*)(&ReadCGMMStruct), (void*)(&single_cgmm_cal_Position), vm_var, (void*)(&construct_monomer_single_chain_block_cgcopolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_cgneighbor_rdf_thread_var), (void*)(&single_cgmm_neighbor_rdf_func),
		(void*)(&_CGMM_rdf_::raw_accumulate_rdf), (void*)(&_CGMM_rdf_::normalize_rdf), (void*)(&_CGMM_rdf_::save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_cgmm_neighbor_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _CGMM_rdf_;
	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<CG_MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_multi_mm_multi_md<CGMM_MD_CELL, CG_MMOLECULE, CG_CLUSTER, VM_VAR<CG_MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> >(
		parfname, ::cgmm_md_cell, (void*)(&serial_multi_cgmm_md_proc), (void*)(&multi_cgmm_cal_Position_NonCell), vm_var, (void*)(&construct_monomer_single_chain_block_cgcopolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_cgneighbor_rdf_thread_var), (void*)(&multi_cgmm_neighbor_rdf_func),
		(void*)(&_CGMM_rdf_::raw_accumulate_rdf), (void*)(&_CGMM_rdf_::normalize_rdf), (void*)(&_CGMM_rdf_::save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}


extern "C" bool single_cgmm_non_neighbor_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _CGMM_rdf_;
	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<CG_MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_single_mm_multi_md<CG_MMOLECULE, CG_CLUSTER, VM_VAR<CG_MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> >(
		parfname, (void*)(&ReadCGMMStruct), (void*)(&single_cgmm_cal_Position), vm_var, (void*)(&construct_monomer_single_chain_block_cgcopolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_all_cgneighbors_rdf_thread_var), (void*)(&single_cgmm_non_neighbor_rdf_func),
		(void*)(&_CGMM_rdf_::raw_accumulate_rdf), (void*)(&_CGMM_rdf_::normalize_rdf), (void*)(&_CGMM_rdf_::save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_cgmm_non_neighbor_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _CGMM_rdf_;
	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<CG_MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_multi_mm_multi_md<CGMM_MD_CELL, CG_MMOLECULE, CG_CLUSTER, VM_VAR<CG_MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> >(
		parfname, ::cgmm_md_cell, (void*)(&serial_multi_cgmm_md_proc), (void*)(&multi_cgmm_cal_Position_InCell), vm_var, (void*)(&construct_monomer_single_chain_block_cgcopolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_all_cgneighbors_rdf_thread_var), (void*)(&multi_cgmm_non_neighbor_rdf_func),
		(void*)(&_CGMM_rdf_::raw_accumulate_rdf), (void*)(&_CGMM_rdf_::normalize_rdf), (void*)(&_CGMM_rdf_::save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_cgmm_non_neighbor_sm_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _CGMM_rdf_;
	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<CG_MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_multi_mm_multi_md<CGMM_MD_CELL, CG_MMOLECULE, CG_CLUSTER, VM_VAR<CG_MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> >(
		parfname, ::cgmm_md_cell, (void*)(&serial_multi_cgmm_md_proc), (void*)(&multi_cgmm_cal_Position_InCell), vm_var, (void*)(&construct_monomer_single_chain_block_cgcopolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_all_cgneighbors_rdf_thread_var), (void*)(&multi_cgmm_non_neighbor_in_sm_rdf_func),
		(void*)(&_CGMM_rdf_::raw_accumulate_rdf), (void*)(&_CGMM_rdf_::normalize_rdf), (void*)(&_CGMM_rdf_::save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_cgmm_nsm_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _CGMM_rdf_;
	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<CG_MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_multi_mm_multi_md<CGMM_MD_CELL, CG_MMOLECULE, CG_CLUSTER, VM_VAR<CG_MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> >(
		parfname, ::cgmm_md_cell, (void*)(&serial_multi_cgmm_md_proc), (void*)(&multi_cgmm_cal_Position_InCell), vm_var, (void*)(&construct_monomer_single_chain_block_cgcopolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_all_cgneighbors_rdf_thread_var), (void*)(&multi_cgmm_nsm_rdf_func),
		(void*)(&_CGMM_rdf_::raw_accumulate_rdf), (void*)(&_CGMM_rdf_::normalize_rdf), (void*)(&_CGMM_rdf_::save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}



/**************************************************************************************
	FUNCTION FOR RDF S(r), and THEIR CORRELATION <S(r) * S(r)>
**************************************************************************************/
namespace _CGMM_rdf_ {
void init_cgneighbor_rdf_correlation_thread_var(MULTI_NEIGHBOR_RDF_VAR *mrdf_var, NEIGHBOR_RDF_CORRELATION_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> *crdf_var) {
	init_vm_rdf_correlation_thread_var<CG_MMOLECULE, CG_CLUSTER>(mrdf_var, crdf_var);
}

void raw_accumulate_mrdf(MULTI_NEIGHBOR_RDF_VAR *mrdf, NEIGHBOR_RDF_CORRELATION_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER>* nrdf) {
	int irdf, n, i;
	for (irdf = 0; irdf < mrdf->n_rdf; irdf++) {
		for (n = 0; n < mrdf->rdf_var[irdf].dt.ndt; n++) {
			raw_distribt_combine<float, float, float, long>(mrdf->rdf_var[irdf].dt.dt[n], nrdf->rdf_var[irdf].dt.dt[n]);
			mrdf->rdf_var[irdf].dt.dt[n].vnormal += nrdf->rdf_var[irdf].dt.dt[n].vnormal;
		}
	}
	for (irdf = 0; irdf < mrdf->n_nn_rdf; irdf++) {
		for (n = 0; n < mrdf->nn_rdf_var[irdf].dt.ndt; n++) {
			raw_distribt_combine<float, float, float, long>(mrdf->nn_rdf_var[irdf].dt.dt[n], nrdf->nn_rdf_var[irdf].dt.dt[n]);
			mrdf->nn_rdf_var[irdf].dt.dt[n].vnormal += nrdf->nn_rdf_var[irdf].dt.dt[n].vnormal;
		}
	}
	for (n = 0; n < mrdf->crdf.ndt; n++) {
		raw_distribt_combine<float, float, float, long>(mrdf->crdf.dt[n], nrdf->crdf.dt[n]);
		mrdf->crdf.dt[n].vnormal += nrdf->crdf.dt[n].vnormal;
	}
}

void single_cgmm_neighbor_crdf_func(NEIGHBOR_RDF_CORRELATION_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> *rdf_thread_var, int nrdf) {
	int irdf, n, i, j, idt, jdt;
	NEIGHBOR_RDF_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> *rdf_var = NULL;
	NDISTRIBT<float, long> *dt = NULL, *rdf = NULL, *rdf1 = NULL, *rdf2 = NULL;
	
	int ndt = 0;
	for (irdf = 0; irdf < nrdf; irdf++) {
		ndt = rdf_thread_var[irdf].n_rdf + rdf_thread_var[irdf].n_nn_rdf;
		dt = new NDISTRIBT<float, long>[ndt];
		idt = 0;
		for (n = 0; n < rdf_thread_var[irdf].n_rdf; n++) {
			setup_ndt(dt[idt], rdf_thread_var[irdf].rdf_var[n].dt); // all rdfs in rdf_thread_var have same range and npts
			single_vm_neighbor_rdf_func1<CG_MMOLECULE, CG_CLUSTER, float, long>(rdf_thread_var[irdf].rdf_var + n, dt[idt]);
			for (i = 0; i < dt[idt].ndt; i++) raw_distribt_combine<float, long, float, long>(rdf_thread_var[irdf].rdf_var[n].dt.dt[i], dt[idt].dt[i]);
			idt++;
		}
		for (n = 0; n < rdf_thread_var[irdf].n_nn_rdf; n++) {
			setup_ndt(dt[idt], rdf_thread_var[irdf].nn_rdf_var[n].dt); // all rdfs in rdf_thread_var have same range and npts
			single_vm_non_neighbor_rdf_func1<CG_MMOLECULE, CG_CLUSTER, float, long>(rdf_thread_var[irdf].nn_rdf_var + n, dt[idt]);
			for (i = 0; i < dt[idt].ndt; i++) raw_distribt_combine<float, long, float, long>(rdf_thread_var[irdf].nn_rdf_var[n].dt.dt[i], dt[idt].dt[i]);
			idt++;
		}

		n = 0;
		for (idt = 0; idt < ndt; idt++) {
			rdf1 = dt + idt;
			for (jdt = 0; jdt < ndt; jdt++) {
				rdf2 = dt + jdt;
				for (j = 0; j < rdf1->dt[0].y.N; j++) {
					rdf_thread_var[irdf].crdf.dt[n].y.v[j] += rdf1->dt[0].y.v[j] * rdf2->dt[0].y.v[j];
					rdf_thread_var[irdf].crdf.dt[n].sn.v[j] += (int)(sqrt((double)(rdf1->dt[0].sn.v[j] * rdf2->dt[0].sn.v[j])));
				}
				n++;
			}
		}
		delete[] dt; dt = NULL;
	}
}

bool read_neighbor_crdf_pars(char *parfname, MULTI_NEIGHBOR_RDF_VAR& mrdf_var) {
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";
	float xf = 0, xt = 0, dxf = 0, dxt = 0;
	int npts = 0, i;
	int nrdf = 0;

	NEIGHBOR_RDF_VAR *rdf_var = NULL;

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	if (!search(in, "[NEIGHBOR-RDF]", buffer)) {
		sprintf(msg, "can not find [NEIGHBOR-RDF] in %s", parfname); show_msg(msg); 
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nrdf) == 1 && nrdf >= 0) {
		mrdf_var.set_neighbor_rdf(nrdf);
	}
	else {
		sprintf(msg, "number of distribution has to be >= 0. see %s", buffer); 
		show_msg(msg, true); in.close(); return false;
	}
	for (i = 0; i < nrdf; i++) {
		rdf_var = mrdf_var.rdf_var + i;
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d %d %f %f %d", &((*rdf_var).cm_uid), &((*rdf_var).nm_uid), &((*rdf_var).nth_neighbor), &xf, &xt, &npts) != 6 || (*rdf_var).nth_neighbor < 1 || xt <= xf || npts <= 0) {
			sprintf(msg, "NEIGHBOR-RDF format : [center monomer uid] [neighbor-monomer uid] [nth-neighbor] [from] [to] [npts], see %s", buffer); show_msg(msg); 
			in.close(); return false;
		}
		(*rdf_var).dt.set(1);
		(*rdf_var).cm_uid += _USER_MONOMER_UID;
		(*rdf_var).nm_uid += _USER_MONOMER_UID;
		(*rdf_var).dt.dt[0].set_distribt_x(xf, xt, npts);

		// check the range and the npts are the same for all distribution
		if (i != 0) {
			dxf = xf - mrdf_var.rdf_var[0].dt.dt[0].xf; dxf = FABS(dxf);
			dxt = xt - mrdf_var.rdf_var[0].dt.dt[0].xt; dxt = FABS(dxt);
			if (mrdf_var.rdf_var[i].dt.dt[0].x.N != mrdf_var.rdf_var[0].dt.dt[0].x.N ||
				dxf >= mrdf_var.rdf_var[0].dt.dt[0].dx || dxt >= mrdf_var.rdf_var[0].dt.dt[0].dx) {
				sprintf(msg, "RDF range [from] [to] [npts] has to be the same for all RDFs"); show_msg(msg); 
				in.close(); return false;
			}
		}
	}

	rewind(in);
	if (!search(in, "[NON-NEIGHBOR-RDF]", buffer)) {
		sprintf(msg, "can not find [NEIGHBOR-RDF] in %s", parfname); show_msg(msg); 
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nrdf) == 1 && nrdf >= 0) {
		mrdf_var.set_non_neighbor_rdf(nrdf);
	}
	else {
		sprintf(msg, "number of distribution has to be >= 0. see %s", buffer); 
		show_msg(msg, true); in.close(); return false;
	}
	for (i = 0; i < nrdf; i++) {
		rdf_var = mrdf_var.nn_rdf_var + i;
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d %d %f %f %d", &((*rdf_var).cm_uid), &((*rdf_var).nm_uid), &((*rdf_var).nth_neighbor), &xf, &xt, &npts) != 6 || (*rdf_var).nth_neighbor < 1 || xt <= xf || npts <= 0) {
			sprintf(msg, "NEIGHBOR-RDF format : [center monomer uid] [neighbor-monomer uid] [nth-neighbor] [from] [to] [npts], see %s", buffer); show_msg(msg); 
			in.close(); return false;
		}
		(*rdf_var).dt.set(1);
		(*rdf_var).cm_uid += _USER_MONOMER_UID;
		(*rdf_var).nm_uid += _USER_MONOMER_UID;
		(*rdf_var).dt.dt[0].set_distribt_x(xf, xt, npts);

		// check the range and the npts are the same for all distribution
		if (mrdf_var.n_rdf > 0) {
			dxf = xf - mrdf_var.rdf_var[0].dt.dt[0].xf; dxf = FABS(dxf);
			dxt = xt - mrdf_var.rdf_var[0].dt.dt[0].xt; dxt = FABS(dxt);
			if (mrdf_var.nn_rdf_var[i].dt.dt[0].x.N != mrdf_var.rdf_var[0].dt.dt[0].x.N ||
				dxf >= mrdf_var.rdf_var[0].dt.dt[0].dx || dxt >= mrdf_var.rdf_var[0].dt.dt[0].dx) {
				sprintf(msg, "RDF range [from] [to] [npts] has to be the same for all RDFs"); show_msg(msg); 
				in.close(); return false;
			}
		}
		else if (i != 0) {
			dxf = xf - mrdf_var.nn_rdf_var[0].dt.dt[0].xf; dxf = FABS(dxf);
			dxt = xt - mrdf_var.nn_rdf_var[0].dt.dt[0].xt; dxt = FABS(dxt);
			if (mrdf_var.nn_rdf_var[i].dt.dt[0].x.N != mrdf_var.nn_rdf_var[0].dt.dt[0].x.N ||
				dxf >= mrdf_var.nn_rdf_var[0].dt.dt[0].dx || dxt >= mrdf_var.nn_rdf_var[0].dt.dt[0].dx) {
				sprintf(msg, "RDF range [from] [to] [npts] has to be the same for all RDFs"); show_msg(msg); 
				in.close(); return false;
			}
		}
	}

	if (mrdf_var.n_rdf + mrdf_var.n_nn_rdf <= 0) {
		sprintf(msg, "total number of distribution has to be > 0. see %s", buffer); 
		show_msg(msg, true); in.close(); return false;
	}
	else mrdf_var.setup_correlated_distribts();

	in.close(); return true;
}

} // end of namespace _CGMM_rdf_

extern "C" bool single_cgmm_neighbor_crdf_single_chain_block_copolymer(char *parfname) {
	using namespace _CGMM_rdf_;
	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<CG_MMOLECULE> vm_var;
	MULTI_NEIGHBOR_RDF_VAR mrdf_var;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_neighbor_crdf_pars(parfname, mrdf_var)) return false;

	rdf_distribt_single_mm_multi_md<CG_MMOLECULE, CG_CLUSTER, VM_VAR<CG_MMOLECULE>, MULTI_NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_CORRELATION_THREAD_VAR<CG_MMOLECULE, CG_CLUSTER> >(
		parfname, (void*)(&ReadCGMMStruct), (void*)(&single_cgmm_cal_Position), vm_var, (void*)(&construct_monomer_single_chain_block_cgcopolymer), 
		&mrdf_var, 1, (void*)(&init_cgneighbor_rdf_correlation_thread_var), (void*)(&single_cgmm_neighbor_crdf_func),
		(void*)(&_CGMM_rdf_::raw_accumulate_rdf), (void*)(&_MM_rdf_::normalize_mrdf), (void*)(&_MM_rdf_::save_mrdf));

	mlog.close();

	return true;
}

