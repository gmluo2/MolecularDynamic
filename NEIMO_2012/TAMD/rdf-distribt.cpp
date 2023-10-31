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
#include "cluster.h"

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
extern bool serial_single_mol_md_proc(long isave);
extern bool serial_multi_mol_md_proc(long isave);
extern "C" bool serial_md_proc(long isave);
extern "C" void clear_read_mdproc();
extern "C" void set_md_mode(int mode);
extern "C" void set_job(int job);
extern void molecule_check_periodic_cell(MMOL_MD_CELL &mmcell);
extern char struct_fname[256];
extern MMOL_MD_CELL mm_md_cell;
extern int MD_MODE;
extern READ_MD_PROC_FILE rmdproc;

extern CGMM_MD_CELL cgmm_md_cell;
extern void molecule_check_periodic_cell(CGMM_MD_CELL &mmcell);

extern bool get_monomer_uid(char *monomer, int &uid);

extern bool read_given_id(char *fname, char* title, ARRAY<short>& id_cluster);

extern bool get_atom_index(char *atom, int &aindx);

#include "distribution.h"
#include "rdf-distribt.h"

namespace _MM_rdf_ {
using namespace _distribution_;

// cneighbor is an array, with m->nCluster of elements
void init_cluster_neighbor(MMOLECULE *m, cNeighbor<CLUSTER> *cneighbor) {
	CLUSTER *pc = NULL;
	CHAIN<CLUSTER> *ch = NULL;
	for (int nc = 0; nc < m->nCluster; nc++) {
		pc = m->cluster + nc;
		grandchildren_chain<CLUSTER>(pc, cneighbor[nc].nth, &ch);
		chain_array<CLUSTER>(&ch, cneighbor[nc].neighbor);
		release_chain<CLUSTER>(&ch, false);
	}
}

bool construct_monomer_single_chain_homopolymer(MMOLECULE *m, int m_uid, int nm_cg, int vm_uid, VIRTUAL_MONOMER<CLUSTER> **cgmonomer, int& ncgMonomers) {
	char msg[256] = "\0";
	CHAIN<CLUSTER> *ch = NULL, *tail = NULL, *pch = NULL;
	int nc = 0, nm = 0, im = 0, icg = 0;
	int nms = 0;
	if (*cgmonomer != NULL) {delete[] *cgmonomer; *cgmonomer = NULL;}
	for (nm = 0; nm < m->nMonomer; nm++) {
		if (m->monomer[nm].uid == m_uid) nms++;
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
	*cgmonomer = new VIRTUAL_MONOMER<CLUSTER>[ncgMonomers];

	for (nm = 0; nm < m->nMonomer; nm++) {
		if (m->monomer[nm].uid != m_uid) continue;
		for (nc = 0; nc < m->monomer[nm].nUnit; nm++) {
			pch = new CHAIN<CLUSTER>; pch->p = m->monomer[nm].u[nc];
			if (ch == NULL) ch = pch;
			else tail->next = pch;
			tail = pch;
		}
		
		im += 1;
		if (im == nm_cg) {
			construct_virtual_monomer<CLUSTER>(ch, *cgmonomer + icg);
			virtual_monomer_setup_hinges<CLUSTER>(*cgmonomer + icg);
			(*cgmonomer)[icg].g.uid = vm_uid;
			im = 0; icg += 1;
			release_chain<CLUSTER>(&ch, false); tail = NULL;
		}
	}

	connect_virtual_monomers<CLUSTER>(*cgmonomer, ncgMonomers);
	return true;
}

void init_neighbor_rdf_thread_var(NEIGHBOR_RDF_VAR* rdf_var, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER>* rdf_thread_var) {
	return init_vm_neighbor_rdf_thread_var<MMOLECULE, CLUSTER>(rdf_var, rdf_thread_var);
}

void init_all_neighbors_rdf_thread_var(NEIGHBOR_RDF_VAR* rdf_var, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER>* rdf_thread_var) {
	return init_vm_all_neighbors_rdf_thread_var<MMOLECULE, CLUSTER>(rdf_var, rdf_thread_var);
}

bool construct_monomer_single_chain_block_copolymer(VM_VAR<MMOLECULE>* vm_var, VIRTUAL_MONOMER<CLUSTER> **cgmonomer, int& ncgMonomers) {
	char msg[256] = "\0";
	CHAIN<CLUSTER> *ch = NULL, *tail = NULL, *pch = NULL;
	int nc = 0, nm = 0, im = 0, icg = 0, iBlock = 0;
	int *nms = new int[vm_var->nBlocks];
	int *ncgm = new int[vm_var->nBlocks];
	if (*cgmonomer != NULL) {delete[] *cgmonomer; *cgmonomer = NULL;}
	for (iBlock = 0; iBlock < vm_var->nBlocks; iBlock++) {nms[iBlock] = 0; ncgm[iBlock] = 0;}

#define RELEASE if (nms != NULL) {delete[] nms; nms = NULL;} \
	if (ncgm != NULL) {delete[] ncgm; ncgm = NULL;}

	MMOLECULE *m = vm_var->m;
	for (nm = 0; nm < m->nMonomer; nm++) {
		for (iBlock = 0; iBlock < vm_var->nBlocks; iBlock++) {
			if (m->monomer[nm].uid == vm_var->muid[iBlock]) {nms[iBlock]++; break;}
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
	*cgmonomer = new VIRTUAL_MONOMER<CLUSTER>[ncgMonomers];

	bool status = false;
	int muid = 0;
	for (nm = 0; nm < m->nMonomer; nm++) {
		status = false;
		for (iBlock = 0; iBlock < vm_var->nBlocks; iBlock++) {
			if (m->monomer[nm].uid == vm_var->muid[iBlock]) {status = true; break;}
		}
		if (!status) continue;
		for (nc = 0; nc < m->monomer[nm].nUnit; nc++) {
			if (ch == NULL) {
				pch = new CHAIN<CLUSTER>; pch->p = m->monomer[nm].u[nc];
				ch = pch;
				muid = m->monomer[nm].uid;
			}
			else if (m->monomer[nm].uid != muid) { // this is a different type of monomer
				// this virtual_monomer can not be constructed
				release_chain<CLUSTER>(&ch, false); tail = NULL;
				status = false; break;
			}
			else {
				pch = new CHAIN<CLUSTER>; pch->p = m->monomer[nm].u[nc];
				tail->next = pch;
			}
			tail = pch;
		}
		if (!status) continue; // this virtual_monomer can not be constructed
		
		im += 1;
		if (im == vm_var->nm_vm[iBlock]) {
			construct_virtual_monomer<CLUSTER>(ch, *cgmonomer + icg);
			virtual_monomer_setup_hinges<CLUSTER>(*cgmonomer + icg);
			(*cgmonomer)[icg].g.uid = vm_var->vmuid[iBlock];
			im = 0; icg += 1;
			release_chain<CLUSTER>(&ch, false); tail = NULL;
		}
	}

	connect_virtual_monomers<CLUSTER>(*cgmonomer, ncgMonomers);
	return true;
#undef RELEASE
}

void calPosition(VIRTUAL_MONOMER<CLUSTER> *vm, int method) {
	int nc;
	float M = 0, m = 0;
	if (method == 0 || method == 1) { // mass center
		V3zero(vm->r)
		for (nc = 0; nc < vm->g.nUnit; nc++) {
			vm->g.u[nc]->MassCenter();
			m = vm->g.u[nc]->M;
			vm->r.v[0] += vm->g.u[nc]->cm.v[0] * m; 
			vm->r.v[1] += vm->g.u[nc]->cm.v[1] * m; 
			vm->r.v[2] += vm->g.u[nc]->cm.v[2] * m;
			M += m;
		}
		vm->r.v[0] /= M; vm->r.v[1] /= M; vm->r.v[2] /= M;
	}
	if (method == 1) {
		// use the first cluser's dr_cmm for the CMM-shift
		vm->r.v[0] += vm->g.u[0]->dr_cmm.v[0];
		vm->r.v[1] += vm->g.u[0]->dr_cmm.v[1];
		vm->r.v[2] += vm->g.u[0]->dr_cmm.v[2];
	}
}

void single_mm_cal_Position(MM_VM<MMOLECULE, CLUSTER>& mvm) {
	for (int i = 0; i < mvm.nVM; i++) calPosition(mvm.vm + i, 0);
}

void multi_mm_cal_Position_NonCell(ARRAY< MM_VM<MMOLECULE, CLUSTER> >& vmarray) { // not position in cell
	int i, j;
	for (i = 0; i < vmarray.n; i++) {
		for (j = 0; j < vmarray.m[i].nVM; j++) calPosition(vmarray.m[i].vm + j, 0);
	}
}

void multi_mm_cal_Position_InCell(ARRAY< MM_VM<MMOLECULE, CLUSTER> >& vmarray) { // position in cell
	int i, j;
	for (i = 0; i < vmarray.n; i++) {
		for (j = 0; j < vmarray.m[i].nVM; j++) calPosition(vmarray.m[i].vm + j, 1);
	}
}

void single_mm_neighbor_rdf_func(NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> *rdf_thread_var, int nrdf) {
	return single_vm_neighbor_rdf_func<MMOLECULE, CLUSTER>(rdf_thread_var, nrdf);
}

void multi_mm_neighbor_rdf_func(MMOL_MD_CELL *md_cell, ARRAY< MM_VM<MMOLECULE, CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> *rdf_thread_var, int nrdf) {
	return multi_vm_neighbor_rdf_func<MMOL_MD_CELL, MMOLECULE, CLUSTER>(md_cell, vmarray, rdf_thread_var, nrdf);
}

void single_mm_non_neighbor_rdf_func(NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> *rdf_thread_var, int nrdf) {
	return single_vm_non_neighbor_rdf_func<MMOLECULE, CLUSTER>(rdf_thread_var, nrdf);
}

void multi_mm_non_neighbor_rdf_func(MMOL_MD_CELL *md_cell, ARRAY< MM_VM<MMOLECULE, CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> *rdf_thread_var, int nrdf) {
	return multi_vm_non_neighbor_rdf_func<MMOL_MD_CELL, MMOLECULE, CLUSTER>(md_cell, vmarray, rdf_thread_var, nrdf);
}

void multi_mm_non_neighbor_in_sm_rdf_func(MMOL_MD_CELL *md_cell, ARRAY< MM_VM<MMOLECULE, CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> *rdf_thread_var, int nrdf) {
	return multi_vm_non_neighbor_in_sm_rdf_func<MMOL_MD_CELL, MMOLECULE, CLUSTER>(md_cell, vmarray, rdf_thread_var, nrdf);
}

void multi_mm_nsm_rdf_func(MMOL_MD_CELL *md_cell, ARRAY< MM_VM<MMOLECULE, CLUSTER> >* vmarray, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> *rdf_thread_var, int nrdf) {
	return multi_vm_nsm_rdf_func<MMOL_MD_CELL, MMOLECULE, CLUSTER>(md_cell, vmarray, rdf_thread_var, nrdf);
}

void raw_accumulate_rdf(NEIGHBOR_RDF_VAR *rdf, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER>* nrdf) {
	int irdf, i;
	for (irdf = 0; irdf < rdf->dt.ndt; irdf++) {
		raw_distribt_combine<float, float, float, long>(rdf->dt.dt[irdf], nrdf->dt.dt[irdf]);
		rdf->dt.dt[irdf].vnormal += nrdf->dt.dt[irdf].vnormal;
	}
}

void normalize_rdf(NEIGHBOR_RDF_VAR *rdf) {
	int irdf, i;
	for (irdf = 0; irdf < rdf->dt.ndt; irdf++) {
		rdf->dt.dt[irdf].eb_raw();
		rdf->dt.dt[irdf].y /= rdf->dt.dt[irdf].vnormal;
		rdf->dt.dt[irdf].eb /= rdf->dt.dt[irdf].vnormal;
	}
}

void normalize_mrdf(MULTI_NEIGHBOR_RDF_VAR *mrdf) {
	int irdf, n, i;
	for (irdf = 0; irdf < mrdf->n_rdf; irdf++) {
		for (n = 0; n < mrdf->rdf_var[irdf].dt.ndt; n++) {
			mrdf->rdf_var[irdf].dt.dt[n].eb_raw();
			mrdf->rdf_var[irdf].dt.dt[n].y /= mrdf->rdf_var[irdf].dt.dt[n].vnormal;
			mrdf->rdf_var[irdf].dt.dt[n].eb /= mrdf->rdf_var[irdf].dt.dt[n].vnormal;
		}
	}
	for (irdf = 0; irdf < mrdf->n_nn_rdf; irdf++) {
		for (n = 0; n < mrdf->nn_rdf_var[irdf].dt.ndt; n++) {
			mrdf->nn_rdf_var[irdf].dt.dt[n].eb_raw();
			mrdf->nn_rdf_var[irdf].dt.dt[n].y /= mrdf->nn_rdf_var[irdf].dt.dt[n].vnormal;
			mrdf->nn_rdf_var[irdf].dt.dt[n].eb /= mrdf->nn_rdf_var[irdf].dt.dt[n].vnormal;
		}
	}
	for (n = 0; n < mrdf->crdf.ndt; n++) {
		mrdf->crdf.dt[n].eb_raw();
		mrdf->crdf.dt[n].y /= mrdf->crdf.dt[n].vnormal;
		mrdf->crdf.dt[n].eb /= mrdf->crdf.dt[n].vnormal;
	}
}

void save_rdf(NEIGHBOR_RDF_VAR *rdf_var, char *title, char *msg_title) {
	char ofname[256] = "\0", msg[256] = "\0";
	ofstream *out = NULL;

	int i, n;
	for (i = 0; i < rdf_var->dt.ndt; i++) {
		sprintf(ofname, "%s-%d.dat", title, i);
		out = new ofstream;
		out->open(ofname);
		rdf_var->dt.dt[i].save_x_y_eb(*out);
		/*
		for (n = 0; n < rdf_var->dt.dt[i].N; n++) {
			(*out)<<rdf_var->dt.dt[i].x[n]<<" "<<rdf_var->dt.dt[i].y[n] / rdf_var->dt.dt[i].dx<<endl;
		}
		*/
		out->close(); delete out; out = NULL;
		sprintf(msg, "%s # %d ==> %s", msg_title, i, ofname);
		show_log(msg, true);
	}
	return;
}

void save_mrdf(MULTI_NEIGHBOR_RDF_VAR *mrdf_var, char *title, char *msg_title) {
	char ofname[256] = "\0", msg[256] = "\0";
	NEIGHBOR_RDF_VAR *rdf_var = NULL;

	int irdf, i;
	ofstream *out = NULL;
	for (irdf = 0; irdf < mrdf_var->n_rdf; irdf++) {
		sprintf(ofname, "%s-%dnrdf", title, irdf);
		sprintf(msg, "%s neighbor-rdf #%d", msg_title, irdf);
		rdf_var = mrdf_var->rdf_var + irdf;
		save_rdf(rdf_var, ofname, msg);
	}
	for (irdf = 0; irdf < mrdf_var->n_nn_rdf; irdf++) {
		sprintf(ofname, "%s-%dnnrdf", title, irdf);
		sprintf(msg, "%s non-neighbor-rdf #%d", msg_title, irdf);
		rdf_var = mrdf_var->nn_rdf_var + irdf;
		save_rdf(rdf_var, ofname, msg);
	}
	for (irdf = 0; irdf < mrdf_var->crdf.ndt; irdf++) {
		sprintf(ofname, "%s-%dcrdf.dat", title, irdf);
		out = new ofstream; out->open(ofname);
		/*
		for (i = 0; i < mrdf_var->crdf.dt[irdf].N; i++) {
			(*out)<<mrdf_var->crdf.dt[irdf].x[i]<<" "<<mrdf_var->crdf.dt[irdf].y[i]<<endl;
		}
		*/
		mrdf_var->crdf.dt[irdf].save_x_y_eb(*out);
		out->close(); delete out; out = NULL;
		sprintf(msg, "%s rdf-correlation #%d ==> %s", msg_title, irdf, ofname);
		show_log(msg, true);
	}
	return;
}

bool read_vm_var(char *parfname, VM_VAR<MMOLECULE> &vm_var) {
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
		else if (!get_monomer_uid(monomer, vm_var.muid[i])) {
			in.close(); return false;
		}
		vm_var.vmuid[i] += _USER_MONOMER_UID;
	}

	in.close(); return true;
}

bool read_neighbor_rdf_pars(char *parfname, NEIGHBOR_RDF_VAR** rdf_var, int &nrdf) {
	if (*rdf_var != NULL) {delete[] *rdf_var; *rdf_var = NULL;} 
	nrdf = 0;

#define RELEASE_RDF if (*rdf_var != NULL) {delete[] *rdf_var; *rdf_var = NULL;} nrdf = 0;

	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";
	float xf = 0, xt = 0;
	int npts = 0, i;

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	if (!search(in, "[NEIGHBOR-RDF]", buffer)) {
		sprintf(msg, "can not find [NEIGHBOR-RDF] in %s", parfname); show_msg(msg); 
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nrdf) == 1 && nrdf > 0) {
		*rdf_var = new NEIGHBOR_RDF_VAR[nrdf];
	}
	else {
		sprintf(msg, "number of distribution has to be > 0. see %s", buffer); 
		show_msg(msg, true); in.close(); return false;
	}
	for (i = 0; i < nrdf; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d %d %f %f %d", &((*rdf_var)[i].cm_uid), &((*rdf_var)[i].nm_uid), &((*rdf_var)[i].nth_neighbor), &xf, &xt, &npts) != 6 || (*rdf_var)[i].nth_neighbor < 1 || xt <= xf || npts <= 0) {
			sprintf(msg, "NEIGHBOR-RDF format : [center monomer uid] [neighbor-monomer uid] [nth-neighbor] [from] [to] [npts], see %s", buffer); show_msg(msg); 
			RELEASE_RDF in.close(); return false;
		}
		(*rdf_var)[i].dt.set(1);
		(*rdf_var)[i].cm_uid += _USER_MONOMER_UID;
		(*rdf_var)[i].nm_uid += _USER_MONOMER_UID;
		(*rdf_var)[i].dt.dt[0].set_distribt_x(xf, xt, npts);
	}
	in.close(); return true;
#undef RELEASE_RDF
}

bool read_non_neighbor_rdf_pars(char *parfname, NEIGHBOR_RDF_VAR** rdf_var, int &nrdf) {
	if (*rdf_var != NULL) {delete[] *rdf_var; *rdf_var = NULL;} 
	nrdf = 0;

#define RELEASE_RDF if (*rdf_var != NULL) {delete[] *rdf_var; *rdf_var = NULL;} nrdf = 0;

	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";
	float xf = 0, xt = 0;
	int npts = 0, i;

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	if (!search(in, "[NON-NEIGHBOR-RDF]", buffer)) {
		sprintf(msg, "can not find [NEIGHBOR-RDF] in %s", parfname); show_msg(msg); 
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nrdf) == 1 && nrdf > 0) {
		*rdf_var = new NEIGHBOR_RDF_VAR[nrdf];
	}
	else {
		sprintf(msg, "number of distribution has to be > 0. see %s", buffer); 
		show_msg(msg, true); in.close(); return false;
	}
	for (i = 0; i < nrdf; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d %d %f %f %d", &((*rdf_var)[i].cm_uid), &((*rdf_var)[i].nm_uid), &((*rdf_var)[i].nth_neighbor), &xf, &xt, &npts) != 6 || (*rdf_var)[i].nth_neighbor < 1 || xt <= xf || npts <= 0) {
			sprintf(msg, "NEIGHBOR-RDF format : [center monomer uid] [neighbor-monomer uid] [nth-neighbor] [from] [to] [npts], see %s", buffer); show_msg(msg); 
			RELEASE_RDF in.close(); return false;
		}
		(*rdf_var)[i].dt.set(1);
		(*rdf_var)[i].cm_uid += _USER_MONOMER_UID;
		(*rdf_var)[i].nm_uid += _USER_MONOMER_UID;
		(*rdf_var)[i].dt.dt[0].set_distribt_x(xf, xt, npts);
	}
	in.close(); return true;
#undef RELEASE_RDF
}

} // end of namespace _MM_rdf_

extern bool ReadMoleculeStruct(ifstream &in, MMOLECULE &mm);
extern bool serial_multi_mol_md_proc(long isave);

extern "C" bool single_mm_neighbor_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;

	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_single_mm_multi_md<MMOLECULE, CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> >(
		parfname, (void*)(&ReadMoleculeStruct), (void*)(&_MM_rdf_::single_mm_cal_Position), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_neighbor_rdf_thread_var), (void*)(&single_mm_neighbor_rdf_func), (void*)(&raw_accumulate_rdf), (void*)(&normalize_rdf), (void*)(&save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_mm_neighbor_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;

	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_multi_mm_multi_md<MMOL_MD_CELL, MMOLECULE, CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> >(
		parfname, ::mm_md_cell, (void*)(&serial_multi_mol_md_proc), (void*)(&multi_mm_cal_Position_NonCell), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_neighbor_rdf_thread_var), (void*)(&multi_mm_neighbor_rdf_func), (void*)(&raw_accumulate_rdf), (void*)(&normalize_rdf), (void*)(&save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}


extern "C" bool single_mm_non_neighbor_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;

	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_single_mm_multi_md<MMOLECULE, CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> >(
		parfname, (void*)(&ReadMoleculeStruct), (void*)(&single_mm_cal_Position), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_all_neighbors_rdf_thread_var), (void*)(&single_mm_non_neighbor_rdf_func), (void*)(&raw_accumulate_rdf), (void*)(&normalize_rdf), (void*)(&save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_mm_non_neighbor_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;

	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_multi_mm_multi_md<MMOL_MD_CELL, MMOLECULE, CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> >(
		parfname, ::mm_md_cell, (void*)(&serial_multi_mol_md_proc), (void*)(&multi_mm_cal_Position_InCell), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_all_neighbors_rdf_thread_var), (void*)(&multi_mm_non_neighbor_rdf_func), (void*)(&raw_accumulate_rdf), (void*)(&normalize_rdf), (void*)(&save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_mm_non_neighbor_sm_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;

	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_multi_mm_multi_md<MMOL_MD_CELL, MMOLECULE, CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> >(
		parfname, ::mm_md_cell, (void*)(&serial_multi_mol_md_proc), (void*)(&multi_mm_cal_Position_InCell), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_all_neighbors_rdf_thread_var), (void*)(&multi_mm_non_neighbor_in_sm_rdf_func), (void*)(&raw_accumulate_rdf), (void*)(&normalize_rdf), (void*)(&save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_mm_nsm_rdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;

	::COMMAND_MD_STOP = false;
	mlog.init("rdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_RDF_VAR *neighbor_rdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_rdf_pars(parfname, &neighbor_rdf_var, nrdf) || neighbor_rdf_var == NULL) return false;

	rdf_distribt_multi_mm_multi_md<MMOL_MD_CELL, MMOLECULE, CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_RDF_VAR, NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER> >(
		parfname, ::mm_md_cell, (void*)(&serial_multi_mol_md_proc), (void*)(&multi_mm_cal_Position_InCell), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_rdf_var, nrdf, (void*)(&init_all_neighbors_rdf_thread_var), (void*)(&multi_mm_nsm_rdf_func), (void*)(&raw_accumulate_rdf), (void*)(&normalize_rdf), (void*)(&save_rdf));

	mlog.close();

	if (neighbor_rdf_var != NULL) {delete[] neighbor_rdf_var; neighbor_rdf_var = NULL;} nrdf = 0;
	return true;
}



//***************************************************************************
//        CORRELATION OF ATOMIC COORDINATES
//***************************************************************************
namespace coordinate_correlation {
int weight_mode = WEIGHT_NUMBER;
/*
#if _SYS_ == _WINDOWS_SYS_
int correlation_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* correlation_thread(void *void_par) { // thread-function
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
*/

void correlation(VATOM_CELL<VATOM> &cell, CORREL_DISTRIBT &cdt) {
	VATOM *patom, *patom1;
	int i, j;
	int n[3];
	float dr0[3], dr[3], dcell[3], d;
	float rmax = cdt.dt.x.v[cdt.dt.x.N - 1], r2max = rmax * rmax, w;
	int ncatom = 0;
	for (i = 0; i < cell.va.n; i++) {
		patom = cell.va.m + i;
		if (cdt.cuid >= 0 && patom->uid != cdt.cuid) continue;
		ncatom++;
		for (j = 0; j < cell.va.n; j++) {
			patom1 = cell.va.m + j;
			if (cdt.ruid >= 0 && patom1->uid != cdt.ruid) continue;
			DVECT3(patom1->r->v, patom->r->v, dr0)
			w = patom->w * patom1->w;

			if (cell.period) {
				for (n[0] = -1; n[0] <= 1; n[0]++) {dcell[0] = n[0] * cell.xd[0]; dr[0] = dr0[0] + dcell[0];
				for (n[1] = -1; n[1] <= 1; n[1]++) {dcell[1] = n[1] * cell.xd[1]; dr[1] = dr0[1] + dcell[1];
				for (n[2] = -1; n[2] <= 1; n[2]++) {dcell[2] = n[2] * cell.xd[2]; dr[2] = dr0[2] + dcell[2];
					DV3ABS2(dr, d) if (d < 1e-6) continue; d = sqrt(d);
					cdt.dt.sample(d, w, true);
				}}}
			}
			else {
				DV3ABS2(dr0, d) if (d < 1e-7) continue; d = sqrt(d);
				cdt.dt.sample(d, w, true);
			}
		}
	}

	cdt.ncatom += ncatom;
}

void correlation1(VATOM_CELL<VATOM1> &cell, CORREL_DISTRIBT &cdt) {
	VATOM1 *patom, *patom1;
	int i, j;
	int n[3];
	float dr0[3], dr[3], dcell[3], dr2[3], dr2_12, d;
	float rmax = cdt.dt.x.v[cdt.dt.x.N - 1], r2max = rmax  * rmax, w;
	int ncatom = 0;
	for (i = 0; i < cell.va.n; i++) {
		patom = cell.va.m + i;
		if (cdt.cuid >= 0 && patom->uid != cdt.cuid) continue;
		ncatom++;
		for (j = 0; j < cell.va.n; j++) {
			patom1 = cell.va.m + j;
			if (cdt.ruid >= 0 && patom1->uid != cdt.ruid) continue;
			DVECT3(patom1->r.v, patom->r.v, dr0)
			w = patom->w * patom1->w;

			if (cell.period) {
				for (n[0] = -1; n[0] <= 1; n[0]++) {dcell[0] = n[0] * cell.xd[0]; dr[0] = dr0[0] + dcell[0]; dr2[0] = dr[0] * dr[0]; if (dr2[0] > r2max) continue;
				for (n[1] = -1; n[1] <= 1; n[1]++) {dcell[1] = n[1] * cell.xd[1]; dr[1] = dr0[1] + dcell[1]; dr2[1] = dr[1] * dr[1]; if (dr2[1] > r2max) continue; dr2_12 = dr2[0] + dr2[1]; if (dr2_12 > r2max) continue;
				for (n[2] = -1; n[2] <= 1; n[2]++) {dcell[2] = n[2] * cell.xd[2]; dr[2] = dr0[2] + dcell[2];
					d = dr2_12 + dr[2] * dr[2]; if (d < 1e-7) continue;
					d = sqrt(d); cdt.dt.sample(d, w, true);
				}}}
			}
			else {
				DV3ABS2(dr0, d) if (d < 1e-7) continue; d = sqrt(d);
				cdt.dt.sample(d, w, true);
			}
		}
	}

	cdt.ncatom += ncatom;
}

void correlation1_ndt(VATOM_CELL<VATOM1> &cell, ARRAY< CORREL_DISTRIBT >& ndt) {
	for (int i = 0; i < ndt.n; i++) correlation1(cell, ndt.m[i]);
}

void correlation_ndt(VATOM_CELL<VATOM> &cell, ARRAY< CORREL_DISTRIBT >& ndt) {
	for (int i = 0; i < ndt.n; i++) correlation(cell, ndt.m[i]);
}

float get_atom_weight(MATOM* patom) {
	float w = 0;
	switch (weight_mode) {
	case WEIGHT_NUMBER: // 0
		w = 1;
		break;
	case WEIGHT_MASS: // 1
		w = patom->par->m;
		break;
	case WEIGHT_INDX: // 2
		w = patom->par->uid;
		break;
	default:
		w = 1;
	}
	return w;
}

void set_vatom(MATOM *patom, VATOM *va) {
	va->uid = patom->par->uid;
	va->r = &(patom->rg);
	va->w = get_atom_weight(patom);
}

void set_vatom(MATOM *patom, VATOM1 *va) {
	va->uid = patom->par->uid;
	memcpy(va->r.v, patom->rg.v, SIZE_V3);
	va->w = get_atom_weight(patom);
}

void construct_vcell_1(MMOL_MD_CELL &mdcell, VATOM_CELL<VATOM1> &vcell) {
	int na = nAtoms(mdcell);
	vcell.va.SetArray(na);
	int imol = 0, ic = 0, ia = 0, n = 0;
	BASIC_CLUSTER *pc = NULL;
	n = 0;
	for (imol = 0; imol < mdcell.mm.n; imol++) {
		for (ic = 0; ic < mdcell.mm.m[imol].nCluster; ic++) {
			pc = (BASIC_CLUSTER*)(mdcell.mm.m[imol].cluster + ic);
			for (ia = 0; ia < pc->nAtoms; ia++) {set_vatom(pc->atom + ia, vcell.va.m + n); n++;}
		}
	}
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		pc = mdcell.sm.m[imol].c;
		for (ia = 0; ia < pc->nAtoms; ia++) {set_vatom(pc->atom + ia, vcell.va.m + n); n++;}
	}
	for (imol = 0; imol < mdcell.pm.n; imol++) {
		pc = mdcell.pm.m[imol].c;
		for (ia = 0; ia < pc->nAtoms; ia++) {set_vatom(pc->atom + ia, vcell.va.m + n); n++;}
	}
	vcell.xd[0] = mdcell.h[0] * 2;
	vcell.xd[1] = mdcell.h[1] * 2;
	vcell.xd[2] = mdcell.h[2] * 2;
	vcell.period = mdcell.period;
}

void construct_vcell(MMOL_MD_CELL &mdcell, VATOM_CELL<VATOM> &vcell) {
	int na = nAtoms(mdcell);
	vcell.va.SetArray(na);
	int imol = 0, ic = 0, ia = 0, n = 0;
	BASIC_CLUSTER *pc = NULL;
	n = 0;
	for (imol = 0; imol < mdcell.mm.n; imol++) {
		for (ic = 0; ic < mdcell.mm.m[imol].nCluster; ic++) {
			pc = (BASIC_CLUSTER*)(mdcell.mm.m[imol].cluster + ic);
			for (ia = 0; ia < pc->nAtoms; ia++) {set_vatom(pc->atom + ia, vcell.va.m + n); n++;}
		}
	}
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		pc = mdcell.sm.m[imol].c;
		for (ia = 0; ia < pc->nAtoms; ia++) {set_vatom(pc->atom + ia, vcell.va.m + n); n++;}
	}
	for (imol = 0; imol < mdcell.pm.n; imol++) {
		pc = mdcell.pm.m[imol].c;
		for (ia = 0; ia < pc->nAtoms; ia++) {set_vatom(pc->atom + ia, vcell.va.m + n); n++;}
	}
	vcell.xd[0] = mdcell.h[0] * 2;
	vcell.xd[1] = mdcell.h[1] * 2;
	vcell.xd[2] = mdcell.h[2] * 2;
	vcell.period = mdcell.period;
}

void update_vcell_1(MMOL_MD_CELL &mdcell, VATOM_CELL<VATOM1> &vcell) {
	int imol = 0, ic = 0, ia = 0, n = 0;
	BASIC_CLUSTER *pc = NULL;
	n = 0;
	for (imol = 0; imol < mdcell.mm.n; imol++) {
		for (ic = 0; ic < mdcell.mm.m[imol].nCluster; ic++) {
			pc = (BASIC_CLUSTER*)(mdcell.mm.m[imol].cluster + ic);
			for (ia = 0; ia < pc->nAtoms; ia++) {memcpy(vcell.va.m[n].r.v, pc->atom[ia].rg.v, SIZE_V3); n++;}
		}
	}
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		pc = mdcell.sm.m[imol].c;
		for (ia = 0; ia < pc->nAtoms; ia++) {memcpy(vcell.va.m[n].r.v, pc->atom[ia].rg.v, SIZE_V3); n++;}
	}
	for (imol = 0; imol < mdcell.pm.n; imol++) {
		pc = mdcell.pm.m[imol].c;
		for (ia = 0; ia < pc->nAtoms; ia++) {memcpy(vcell.va.m[n].r.v, pc->atom[ia].rg.v, SIZE_V3); n++;}
	}

	vcell.xd[0] = mdcell.h[0] * 2;
	vcell.xd[1] = mdcell.h[1] * 2;
	vcell.xd[2] = mdcell.h[2] * 2;
}

void update_vcell(MMOL_MD_CELL &mdcell, VATOM_CELL<VATOM1> &vcell) {
	// do not need to update the atomic position
	vcell.xd[0] = mdcell.h[0] * 2;
	vcell.xd[1] = mdcell.h[1] * 2;
	vcell.xd[2] = mdcell.h[2] * 2;
}

bool read_correlation_pars(char *parfname, ARRAY<CORREL_DISTRIBT> &cdt) {
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";
	char atom1[50] = "\0", atom2[50] = "\0", *bf = NULL;
	float xf = 0, xt = 0;
	int npts = 0, i;
	int ndt = 0;
	int iRadialDistribt = 0;
	int nitem = 0;

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	if (!search(in, "[ATOMIC_CORRELATION]", buffer)) {
		sprintf(msg, "can not find [ATOMIC_CORRELATION] in %s", parfname); show_msg(msg); 
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &ndt) == 1 && ndt > 0) {
		cdt.SetArray(ndt);
	}
	else {
		sprintf(msg, "number of correlations has to be > 0. see %s, from %s", buffer, parfname); 
		show_msg(msg, true); in.close(); return false;
	}
	for (i = 0; i < ndt; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d", &(cdt.m[i].cuid), &(cdt.m[i].ruid)) != 2) {
			if (sscanf(buffer, "%s %s", atom1, atom2) != 2) {
				sprintf(msg, "ATOMIC CORRELATION format : [center-atom uid] [surrounded-atom uid] [from] [to] [npts] [Radial Distribt ? 1/0]; [z1 from] -- [z1 to], [z2 from] -- [z2 to]. see :"); 
				show_msg(msg);  show_msg(buffer, true);
				cdt.release(); in.close(); return false;
			}
			else {
				if (!get_atom_index(atom1, cdt.m[i].cuid)) {
					sprintf(msg, "unknown atom: %s. see :", atom1); 
					show_msg(msg);  show_msg(buffer, true);
					cdt.release(); in.close(); return false;
				}
				if (!get_atom_index(atom2, cdt.m[i].ruid)) {
					sprintf(msg, "unknown atom: %s. see :", atom2); 
					show_msg(msg);  show_msg(buffer, true);
					cdt.release(); in.close(); return false;
				}
			}
		}
		bf = buffer; NextSection(&bf, seperator); NextSection(&bf, seperator);
		nitem = sscanf(bf, "%f %f %d %d", &xf, &xt, &npts, &iRadialDistribt);
		if (nitem == 3) {
			iRadialDistribt = 0; // default is partitial distribution
			sprintf(msg, "use default choice -- partitial distribution, not radial distribution for correlation between atom [%d - %d]", cdt.m[i].cuid, cdt.m[i].ruid);
			show_log(msg, true);
		}
		if (nitem < 3 || xt <= xf || npts <= 0) {
			sprintf(msg, "ATOMIC CORRELATION format : [center-atom uid] [surrounded-atom uid] [from] [to] [npts] [Radial Distribt ? 1/0], see %s", buffer); show_msg(msg); 
			cdt.release(); in.close(); return false;
		}
		cdt.m[i].dt.set_distribt_x(xf, xt, npts);
		cdt.m[i].dt.reset_distribt();
		cdt.m[i].dt.vnormal = 0;
		cdt.m[i].bRadialDistribt = (iRadialDistribt > 0 ? true : false); 
	}

	in.close(); return true;
}








//********************************************************************************************************
//     X-Y 2d CORRELATION OF ATOMIC COORDINATES, WITH Z IN REGIONS [z1, z2], [z3, z4] RESPECTIVELY
//********************************************************************************************************
/*
#if _SYS_ == _WINDOWS_SYS_
int correlation_2d_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* correlation_2d_thread(void *void_par) { // thread-function
#endif
	typedef void (*FUNC)(VATOM_CELL<VATOM1>&, ARRAY< CORREL_2D_DISTRIBT >&);
	CORREL_2D_THREAD *par = (CORREL_2D_THREAD*)void_par;
	((FUNC)(par->func))(par->var->cell, par->var->dt);
	par->var->nconf++;
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};
*/

void correlation_2d(VATOM_CELL<VATOM> &cell, CORREL_2D_DISTRIBT &cdt) {
	VATOM *patom, *patom1;
	int i, j;
	int n[3];
	float dr0[3], dr[3], dcell[3], d;
	float rmax = cdt.dt.x.v[cdt.dt.x.N - 1], r2max = rmax * rmax, w;
	int ncatom = 0;
	for (i = 0; i < cell.va.n; i++) {
		patom = cell.va.m + i;
		if (!cdt.inZrg1(patom->r->v[2])) continue;
		if (cdt.cuid >= 0 && patom->uid != cdt.cuid) continue;
		ncatom++;
		for (j = 0; j < cell.va.n; j++) {
			patom1 = cell.va.m + j;
			if (!cdt.inZrg2(patom1->r->v[2])) continue;
			if (cdt.ruid >= 0 && patom1->uid != cdt.ruid) continue;
			DVECT3(patom1->r->v, patom->r->v, dr0)
			w = patom->w * patom1->w;

			if (cell.period) {
				for (n[0] = -1; n[0] <= 1; n[0]++) {dcell[0] = n[0] * cell.xd[0]; dr[0] = dr0[0] + dcell[0];
				for (n[1] = -1; n[1] <= 1; n[1]++) {dcell[1] = n[1] * cell.xd[1]; dr[1] = dr0[1] + dcell[1];
				//for (n[2] = -1; n[2] <= 1; n[2]++) {dcell[2] = n[2] * cell.xd[2]; dr[2] = dr0[2] + dcell[2];
				dr[2] = dr0[2]; d = dr[0] * dr[0] + dr[1] * dr[1];
					if (d < 1e-6) continue; d = sqrt(d);
					cdt.dt.sample(d, w, true);
				}}
			}
			else {
				d = dr0[0] * dr0[0] + dr0[1] * dr0[1];
				if (d < 1e-7) continue; d = sqrt(d);
				cdt.dt.sample(d, w, true);
			}
		}
	}

	cdt.ncatom += ncatom;
}

void correlation1_2d(VATOM_CELL<VATOM1> &cell, CORREL_2D_DISTRIBT &cdt) {
	VATOM1 *patom, *patom1;
	int i, j;
	int n[3];
	float dr0[3], dr[3], dcell[3], dr2[3], dr2_12, d;
	float rmax = cdt.dt.x.v[cdt.dt.x.N - 1], r2max = rmax  * rmax, w;
	int ncatom = 0;
	for (i = 0; i < cell.va.n; i++) {
		patom = cell.va.m + i;
		if (!cdt.inZrg1(patom->r.v[2])) continue;
		if (cdt.cuid >= 0 && patom->uid != cdt.cuid) continue;
		ncatom++;
		for (j = 0; j < cell.va.n; j++) {
			patom1 = cell.va.m + j;
			if (!cdt.inZrg2(patom1->r.v[2])) continue;
			if (cdt.ruid >= 0 && patom1->uid != cdt.ruid) continue;
			DVECT3(patom1->r.v, patom->r.v, dr0)
			w = patom->w * patom1->w;

			if (cell.period) {
				for (n[0] = -1; n[0] <= 1; n[0]++) {dcell[0] = n[0] * cell.xd[0]; dr[0] = dr0[0] + dcell[0]; dr2[0] = dr[0] * dr[0]; if (dr2[0] > r2max) continue;
				for (n[1] = -1; n[1] <= 1; n[1]++) {dcell[1] = n[1] * cell.xd[1]; dr[1] = dr0[1] + dcell[1]; dr2[1] = dr[1] * dr[1]; if (dr2[1] > r2max) continue; dr2_12 = dr2[0] + dr2[1]; if (dr2_12 > r2max) continue;
				//for (n[2] = -1; n[2] <= 1; n[2]++) {dcell[2] = n[2] * cell.xd[2]; dr[2] = dr0[2] + dcell[2];
				dr[2] = dr0[2]; 
					d = dr2_12; if (d < 1e-7) continue;
					d = sqrt(d); cdt.dt.sample(d, w, true);
				}}
			}
			else {
				d = dr0[0] * dr0[0] + dr0[1] * dr0[1];
				if (d < 1e-7) continue; d = sqrt(d);
				cdt.dt.sample(d, w, true);
			}
		}
	}

	cdt.ncatom += ncatom;
}

void correlation1_2d_ndt(VATOM_CELL<VATOM1> &cell, ARRAY< CORREL_2D_DISTRIBT >& ndt) {
	for (int i = 0; i < ndt.n; i++) correlation1_2d(cell, ndt.m[i]);
}

void correlation_2d_ndt(VATOM_CELL<VATOM> &cell, ARRAY< CORREL_2D_DISTRIBT >& ndt) {
	for (int i = 0; i < ndt.n; i++) correlation_2d(cell, ndt.m[i]);
}


bool read_2d_correlation_pars(char *parfname, ARRAY<CORREL_2D_DISTRIBT> &cdt) {
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";
	char atom1[50] = "\0", atom2[50] = "\0", *bf = NULL;
	float xf = 0, xt = 0;
	int npts = 0, i;
	int ndt = 0;
	int iRadialDistribt = 0;
	int nitem = 0;

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	if (!search(in, "[ATOMIC_CORRELATION]", buffer)) {
		sprintf(msg, "can not find [ATOMIC_CORRELATION] in %s", parfname); show_msg(msg); 
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &ndt) == 1 && ndt > 0) {
		cdt.SetArray(ndt);
	}
	else {
		sprintf(msg, "number of correlations has to be > 0. see %s, from %s", buffer, parfname); 
		show_msg(msg, true); in.close(); return false;
	}
	for (i = 0; i < ndt; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d", &(cdt.m[i].cuid), &(cdt.m[i].ruid)) != 2) {
			if (sscanf(buffer, "%s %s", atom1, atom2) != 2) {
				sprintf(msg, "ATOMIC CORRELATION format : [center-atom uid] [surrounded-atom uid] [from] [to] [npts] [Radial Distribt ? 1/0]; [z1 from] -- [z1 to], [z2 from] -- [z2 to]. see :"); 
				show_msg(msg);  show_msg(buffer, true);
				cdt.release(); in.close(); return false;
			}
			else {
				if (!get_atom_index(atom1, cdt.m[i].cuid)) {
					sprintf(msg, "unknown atom: %s. see :", atom1); 
					show_msg(msg);  show_msg(buffer, true);
					cdt.release(); in.close(); return false;
				}
				if (!get_atom_index(atom2, cdt.m[i].ruid)) {
					sprintf(msg, "unknown atom: %s. see :", atom2); 
					show_msg(msg);  show_msg(buffer, true);
					cdt.release(); in.close(); return false;
				}
			}
		}
		bf = buffer; NextSection(&bf, seperator); NextSection(&bf, seperator);
		nitem = sscanf(bf, "%f %f %d %d ; %f -- %f , %f -- %f", &xf, &xt, &npts, &iRadialDistribt,
			&(cdt.m[i].rg1.zf), &(cdt.m[i].rg1.zt), &(cdt.m[i].rg2.zf), &(cdt.m[i].rg2.zt));
		if (nitem < 8 || xt <= xf || npts <= 0) {
			sprintf(msg, "ATOMIC CORRELATION format : [center-atom uid] [surrounded-atom uid] [from] [to] [npts] [Radial Distribt ? 1/0]; [z1 from] -- [z1 to], [z2 from] -- [z2 to]. see :"); 
			show_msg(msg);  show_msg(buffer, true);
			cdt.release(); in.close(); return false;
		}
		cdt.m[i].dt.set_distribt_x(xf, xt, npts);
		cdt.m[i].dt.reset_distribt();
		cdt.m[i].dt.vnormal = 0;
		cdt.m[i].bRadialDistribt = (iRadialDistribt > 0 ? true : false); 
	}

	in.close(); return true;
}

} // namespace coordinate_correlation

extern "C" bool atomic_correlation_single_thread(char *parfname) {
	using namespace coordinate_correlation;
	::COMMAND_MD_STOP = false;
	mlog.init("correlation.log", false, true);

	ARRAY<CORREL_DISTRIBT> cdt;
	bool bRadial = true;

	if (!read_correlation_pars(parfname, cdt)) return false;

	bool status = correlation_multi_md_singlethread<MMOL_MD_CELL, VATOM_CELL<VATOM>, CORREL_DISTRIBT >(
		parfname, ::mm_md_cell, cdt, (void*)(&serial_multi_mol_md_proc), (void*)(&update_atom_rg), (void*)(&construct_vcell), (void*)(&update_vcell), 
		(void*)(&correlation_ndt));

	mlog.close();

	return true;
}

extern "C" bool atomic_correlation_multi_thread(char *parfname) {
	using namespace coordinate_correlation;
	::COMMAND_MD_STOP = false;
	mlog.init("correlation.log", false, true);

	ARRAY<CORREL_DISTRIBT> cdt;

	if (!read_correlation_pars(parfname, cdt)) return false;

	correlation_multi_md_multithread<MMOL_MD_CELL, VATOM_CELL<VATOM1>, CORREL_DISTRIBT >(
		parfname, ::mm_md_cell, cdt, (void*)(&serial_multi_mol_md_proc), (void*)(&update_atom_rg), (void*)(&construct_vcell_1), (void*)(&update_vcell_1), 
		(void*)(&correlation1_ndt));

	mlog.close();

	return true;
}




extern "C" bool atomic_2d_correlation_single_thread(char *parfname) {
	using namespace coordinate_correlation;
	::COMMAND_MD_STOP = false;
	mlog.init("correlation.log", false, true);

	ARRAY<CORREL_2D_DISTRIBT> cdt;
	bool bRadial = true;

	if (!read_2d_correlation_pars(parfname, cdt)) return false;

	bool status = correlation_multi_md_singlethread<MMOL_MD_CELL, VATOM_CELL<VATOM>, CORREL_2D_DISTRIBT >(
		parfname, ::mm_md_cell, cdt, (void*)(&serial_multi_mol_md_proc), (void*)(&update_atom_rg), (void*)(&construct_vcell), (void*)(&update_vcell), 
		(void*)(&correlation_2d_ndt));

	mlog.close();

	return true;
}

extern "C" bool atomic_2d_correlation_multi_thread(char *parfname) {
	using namespace coordinate_correlation;
	::COMMAND_MD_STOP = false;
	mlog.init("correlation.log", false, true);

	ARRAY<CORREL_2D_DISTRIBT> cdt;

	if (!read_2d_correlation_pars(parfname, cdt)) return false;

	correlation_multi_md_multithread<MMOL_MD_CELL, VATOM_CELL<VATOM1>, CORREL_2D_DISTRIBT >(
		parfname, ::mm_md_cell, cdt, (void*)(&serial_multi_mol_md_proc), (void*)(&update_atom_rg), (void*)(&construct_vcell_1), (void*)(&update_vcell_1), 
		(void*)(&correlation1_2d_ndt));

	mlog.close();

	return true;
}









namespace coordinate_correlation {
//********************************************************************************************************
//     X-Y-Z 3d CORRELATION OF ATOMIC COORDINATES, WITH Z IN REGIONS [z1, z2], [z3, z4] RESPECTIVELY
//********************************************************************************************************
/*
#if _SYS_ == _WINDOWS_SYS_
int correlation0_2d_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* correlation0_2d_thread(void *void_par) { // thread-function
#endif
	typedef void (*FUNC)(VATOM_CELL<VATOM1>&, ARRAY< CORREL_2D_DISTRIBT >&);
	CORREL_2D_THREAD *par = (CORREL_2D_THREAD*)void_par;
	((FUNC)(par->func))(par->var->cell, par->var->dt);
	par->var->nconf++;
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};
*/

void correlation_3d_2d(VATOM_CELL<VATOM> &cell, CORREL_2D_DISTRIBT &cdt) {
	VATOM *patom, *patom1;
	int i, j;
	int n[3];
	float dr0[3], dr[3], dcell[3], d;
	float rmax = cdt.dt.x.v[cdt.dt.x.N - 1], r2max = rmax * rmax, w;
	int ncatom = 0;
	for (i = 0; i < cell.va.n; i++) {
		patom = cell.va.m + i;
		if (!cdt.inZrg1(patom->r->v[2])) continue;
		if (cdt.cuid >= 0 && patom->uid != cdt.cuid) continue;
		ncatom++;
		for (j = 0; j < cell.va.n; j++) {
			patom1 = cell.va.m + j;
			if (!cdt.inZrg2(patom1->r->v[2])) continue;
			if (cdt.ruid >= 0 && patom1->uid != cdt.ruid) continue;
			DVECT3(patom1->r->v, patom->r->v, dr0)
			w = patom->w * patom1->w;

			if (cell.period) {
				for (n[0] = -1; n[0] <= 1; n[0]++) {dcell[0] = n[0] * cell.xd[0]; dr[0] = dr0[0] + dcell[0];
				for (n[1] = -1; n[1] <= 1; n[1]++) {dcell[1] = n[1] * cell.xd[1]; dr[1] = dr0[1] + dcell[1];
				//for (n[2] = -1; n[2] <= 1; n[2]++) {dcell[2] = n[2] * cell.xd[2]; dr[2] = dr0[2] + dcell[2];
				dr[2] = dr0[2]; d = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
					if (d < 1e-6) continue; d = sqrt(d);
					cdt.dt.sample(d, w, true);
				}}
			}
			else {
				d = dr0[0] * dr0[0] + dr0[1] * dr0[1] + dr0[2] * dr0[2];
				if (d < 1e-7) continue; d = sqrt(d);
				cdt.dt.sample(d, w, true);
			}
		}
	}

	cdt.ncatom += ncatom;
}

void correlation1_3d_2d(VATOM_CELL<VATOM1> &cell, CORREL_2D_DISTRIBT &cdt) {
	VATOM1 *patom, *patom1;
	int i, j;
	int n[3];
	float dr0[3], dr[3], dcell[3], dr2[3], dr2_12, d;
	float rmax = cdt.dt.x.v[cdt.dt.x.N - 1], r2max = rmax  * rmax, w;
	int ncatom = 0;
	for (i = 0; i < cell.va.n; i++) {
		patom = cell.va.m + i;
		if (!cdt.inZrg1(patom->r.v[2])) continue;
		if (cdt.cuid >= 0 && patom->uid != cdt.cuid) continue;
		ncatom++;
		for (j = 0; j < cell.va.n; j++) {
			patom1 = cell.va.m + j;
			if (!cdt.inZrg2(patom1->r.v[2])) continue;
			if (cdt.ruid >= 0 && patom1->uid != cdt.ruid) continue;
			DVECT3(patom1->r.v, patom->r.v, dr0)
			w = patom->w * patom1->w;

			if (cell.period) {
				for (n[0] = -1; n[0] <= 1; n[0]++) {dcell[0] = n[0] * cell.xd[0]; dr[0] = dr0[0] + dcell[0]; dr2[0] = dr[0] * dr[0]; if (dr2[0] > r2max) continue;
				for (n[1] = -1; n[1] <= 1; n[1]++) {dcell[1] = n[1] * cell.xd[1]; dr[1] = dr0[1] + dcell[1]; dr2[1] = dr[1] * dr[1]; if (dr2[1] > r2max) continue; dr2_12 = dr2[0] + dr2[1]; if (dr2_12 > r2max) continue;
				//for (n[2] = -1; n[2] <= 1; n[2]++) {dcell[2] = n[2] * cell.xd[2]; dr[2] = dr0[2] + dcell[2];
				dr[2] = dr0[2]; dr2[2] = dr[2] * dr[2];
				dr2_12 += dr2[2];
					d = dr2_12; if (d < 1e-7) continue;
					d = sqrt(d); cdt.dt.sample(d, w, true);
				}}
			}
			else {
				d = dr0[0] * dr0[0] + dr0[1] * dr0[1] + dr0[2] * dr0[2];
				if (d < 1e-7) continue; d = sqrt(d);
				cdt.dt.sample(d, w, true);
			}
		}
	}

	cdt.ncatom += ncatom;
}

void correlation1_3d_2d_ndt(VATOM_CELL<VATOM1> &cell, ARRAY< CORREL_2D_DISTRIBT >& ndt) {
	for (int i = 0; i < ndt.n; i++) correlation1_3d_2d(cell, ndt.m[i]);
}

void correlation_3d_2d_ndt(VATOM_CELL<VATOM> &cell, ARRAY< CORREL_2D_DISTRIBT >& ndt) {
	for (int i = 0; i < ndt.n; i++) correlation_3d_2d(cell, ndt.m[i]);
}


bool read_3d_2d_correlation_pars(char *parfname, ARRAY<CORREL_2D_DISTRIBT> &cdt) {
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";
	char atom1[50] = "\0", atom2[50] = "\0", *bf = NULL;
	float xf = 0, xt = 0;
	int npts = 0, i;
	int ndt = 0;
	int iRadialDistribt = 0;
	int nitem = 0;

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	if (!search(in, "[ATOMIC_CORRELATION]", buffer)) {
		sprintf(msg, "can not find [ATOMIC_CORRELATION] in %s", parfname); show_msg(msg); 
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &ndt) == 1 && ndt > 0) {
		cdt.SetArray(ndt);
	}
	else {
		sprintf(msg, "number of correlations has to be > 0. see %s, from %s", buffer, parfname); 
		show_msg(msg, true); in.close(); return false;
	}
	for (i = 0; i < ndt; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d", &(cdt.m[i].cuid), &(cdt.m[i].ruid)) != 2) {
			if (sscanf(buffer, "%s %s", atom1, atom2) != 2) {
				sprintf(msg, "ATOMIC CORRELATION format : [center-atom uid] [surrounded-atom uid] [from] [to] [npts] [Radial Distribt ? 1/0]; [z1 from] -- [z1 to], [z2 from] -- [z2 to]. see :"); 
				show_msg(msg);  show_msg(buffer, true);
				cdt.release(); in.close(); return false;
			}
			else {
				if (!get_atom_index(atom1, cdt.m[i].cuid)) {
					sprintf(msg, "unknown atom: %s. see :", atom1); 
					show_msg(msg);  show_msg(buffer, true);
					cdt.release(); in.close(); return false;
				}
				if (!get_atom_index(atom2, cdt.m[i].ruid)) {
					sprintf(msg, "unknown atom: %s. see :", atom2); 
					show_msg(msg);  show_msg(buffer, true);
					cdt.release(); in.close(); return false;
				}
			}
		}
		bf = buffer; NextSection(&bf, seperator); NextSection(&bf, seperator);
		nitem = sscanf(bf, "%f %f %d %d ; %f -- %f , %f -- %f", &xf, &xt, &npts, &iRadialDistribt,
			&(cdt.m[i].rg1.zf), &(cdt.m[i].rg1.zt), &(cdt.m[i].rg2.zf), &(cdt.m[i].rg2.zt));
		if (nitem < 8 || xt <= xf || npts <= 0) {
			sprintf(msg, "ATOMIC CORRELATION format : [center-atom uid] [surrounded-atom uid] [from] [to] [npts] [Radial Distribt ? 1/0]; [z1 from] -- [z1 to], [z2 from] -- [z2 to]. see :"); 
			show_msg(msg);  show_msg(buffer, true);
			cdt.release(); in.close(); return false;
		}
		cdt.m[i].dt.set_distribt_x(xf, xt, npts);
		cdt.m[i].dt.reset_distribt();
		cdt.m[i].dt.vnormal = 0;
		cdt.m[i].bRadialDistribt = (iRadialDistribt > 0 ? true : false); 
	}

	in.close(); return true;
}

} // namespace coordinate_correlation



extern "C" bool atomic_2d_correlation_3d_single_thread(char *parfname) {
	using namespace coordinate_correlation;
	::COMMAND_MD_STOP = false;
	mlog.init("correlation.log", false, true);

	ARRAY<CORREL_2D_DISTRIBT> cdt;
	bool bRadial = true;

	if (!read_2d_correlation_pars(parfname, cdt)) return false;

	bool status = correlation_multi_md_singlethread<MMOL_MD_CELL, VATOM_CELL<VATOM>, CORREL_2D_DISTRIBT >(
		parfname, ::mm_md_cell, cdt, (void*)(&serial_multi_mol_md_proc), (void*)(&update_atom_rg), (void*)(&construct_vcell), (void*)(&update_vcell), 
		(void*)(&correlation_3d_2d_ndt));

	mlog.close();

	return true;
}

extern "C" bool atomic_2d_correlation_3d_multi_thread(char *parfname) {
	using namespace coordinate_correlation;
	::COMMAND_MD_STOP = false;
	mlog.init("correlation.log", false, true);

	ARRAY<CORREL_2D_DISTRIBT> cdt;

	if (!read_2d_correlation_pars(parfname, cdt)) return false;

	correlation_multi_md_multithread<MMOL_MD_CELL, VATOM_CELL<VATOM1>, CORREL_2D_DISTRIBT >(
		parfname, ::mm_md_cell, cdt, (void*)(&serial_multi_mol_md_proc), (void*)(&update_atom_rg), (void*)(&construct_vcell_1), (void*)(&update_vcell_1), 
		(void*)(&correlation1_3d_2d_ndt));

	mlog.close();

	return true;
}

