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

#include "MM.h"
#include "Mol.h"
#include "EwaldSumVar.h"
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


extern bool read_given_id(char *fname, char* title, ARRAY<short>& id_cluster);

#include "distribution.h"
#include "rdf-distribt.h"
#include "dipole-distrbt.h"

namespace _dipole_distribution_ {
using namespace _distribution_;
using namespace _MM_rdf_;

void init_dipole_distrbt_thread_var(DIPOLE_DISTRIBT_VAR* distrbt_var, DIPOLE_DISTRBT_THREAD_VAR<MMOLECULE, CLUSTER>* distrbt_thread_var) {
	init_vm_neighbor_rdf_thread_var<MMOLECULE, CLUSTER>(distrbt_var, distrbt_thread_var);
	
	int i = 0;
	distrbt_thread_var->mu_dt.set(distrbt_var->mu_dt.ndt);
	distrbt_thread_var->theta_dt.set(distrbt_var->theta_dt.ndt);
	for (i = 0; i < distrbt_thread_var->mu_dt.ndt; i++) {
		distrbt_thread_var->mu_dt.dt[i].set_distribt_x(distrbt_var->mu_dt.dt[i].xf, distrbt_var->mu_dt.dt[i].xt, distrbt_var->mu_dt.dt[i].x.N - 1);
	}
	for (i = 0; i < distrbt_thread_var->theta_dt.ndt; i++) {
		distrbt_thread_var->theta_dt.dt[i].set_distribt_x(distrbt_var->theta_dt.dt[i].xf, distrbt_var->theta_dt.dt[i].xt, distrbt_var->theta_dt.dt[i].x.N - 1);
	}
}

void cal_dipole(VIRTUAL_MONOMER<CLUSTER> *vm, VECTOR3 &mu, float &q) {
	q = 0; V3zero(mu)
	int nc, na;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	for (nc = 0; nc < vm->g.nUnit; nc++) {
		pc = vm->g.u[nc];
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			mu.v[0] += patom->c0 * (patom->r.v[0] - vm->r.v[0]);
			mu.v[1] += patom->c0 * (patom->r.v[1] - vm->r.v[1]);
			mu.v[2] += patom->c0 * (patom->r.v[2] - vm->r.v[2]);
			if (patom->mMP > 0) {
				mu.v[0] += patom->mu.v[0];
				mu.v[1] += patom->mu.v[1];
				mu.v[2] += patom->mu.v[2];
			}
			q += patom->c0;
		}
	}
}

void single_mm_dipole_distrbt(DIPOLE_DISTRBT_THREAD_VAR<MMOLECULE, CLUSTER> *dipole_distrbt_thread_var, NDISTRIBT<float, long> &mu_dt, NDISTRIBT<float, long> &theta_dt) {
	DIPOLE_DISTRBT_THREAD_VAR<MMOLECULE, CLUSTER> *var = dipole_distrbt_thread_var;
	float mu_v = 0, theta = 0, q1 = 0, q2 = 0;
	VECTOR3 dr, mu, mu1, mu2;
	VIRTUAL_MONOMER<CLUSTER> *vm1 = NULL, *vm2 = NULL;
	int idt, ivm1 = 0, ivm2 = 0, npt = 0;

	for (ivm1 = 0; ivm1 < var->vm->nVM; ivm1++) {
		vm1 = var->vm->vm + ivm1;
		if (vm1->g.uid != var->cm_uid) continue;
		cal_dipole(vm1, mu1, q1);
		ivm2 = ivm1 + 1; if (ivm2 >= var->vm->nVM) break;
		vm2 = vm1 + 1;
		if (vm2->g.uid != var->nm_uid) continue;
		VECT3(vm1->r, vm2->r, dr)
		cal_dipole(vm2, mu2, q2);
		mu.v[0] = (mu1.v[0] + mu2.v[0] + q1 * vm1->r.v[0] + q2 * vm2->r.v[0]) / 2;
		mu.v[1] = (mu1.v[1] + mu2.v[0] + q1 * vm1->r.v[1] + q2 * vm2->r.v[1]) / 2;
		mu.v[2] = (mu1.v[2] + mu2.v[0] + q1 * vm1->r.v[2] + q2 * vm2->r.v[2]) / 2;
		mu_v = mu.abs();
		theta = angle(mu, dr);
		for (idt = 0; idt < mu_dt.ndt; idt++) {
			if (mu_dt.dt[idt].sample(mu_v, true)) mu_dt.dt[idt].vnormal += 1; // total number of samples
		}
		for (idt = 0; idt < theta_dt.ndt; idt++) {
			if (theta_dt.dt[idt].sample(theta, true)) theta_dt.dt[idt].vnormal += 1; // total number of samples
		}
	}
}

void single_mm_dipole_distrbt_func(DIPOLE_DISTRBT_THREAD_VAR<MMOLECULE, CLUSTER> *dipole_distrbt_thread_var, int nrdf) {
	int irdf;
	DIPOLE_DISTRBT_THREAD_VAR<MMOLECULE, CLUSTER> *var = NULL;
	for (irdf = 0; irdf < nrdf; irdf++) {
		var = dipole_distrbt_thread_var + irdf;
		single_mm_dipole_distrbt(var, var->mu_dt, var->theta_dt);
	}
}

void multi_mm_dipole_distrbt_func(MMOL_MD_CELL *md_cell, ARRAY< MM_VM<MMOLECULE, CLUSTER> >* vmarray, DIPOLE_DISTRBT_THREAD_VAR<MMOLECULE, CLUSTER> *dipole_distrbt_thread_var, int nrdf) {
	return single_mm_dipole_distrbt_func(dipole_distrbt_thread_var, nrdf);
}

void raw_accumulate_rdf(DIPOLE_DISTRIBT_VAR *distrbt_var, DIPOLE_DISTRBT_THREAD_VAR<MMOLECULE, CLUSTER>* distrbt_thread_var) {
	int irdf, i;
	for (irdf = 0; irdf < distrbt_var->mu_dt.ndt; irdf++) {
		raw_distribt_combine<float, float, float, long>(distrbt_var->mu_dt.dt[irdf], distrbt_thread_var->mu_dt.dt[irdf]);
		distrbt_var->mu_dt.dt[irdf].vnormal += distrbt_thread_var->mu_dt.dt[irdf].vnormal;
	}
	for (irdf = 0; irdf < distrbt_var->theta_dt.ndt; irdf++) {
		raw_distribt_combine<float, float, float, long>(distrbt_var->theta_dt.dt[irdf], distrbt_thread_var->theta_dt.dt[irdf]);
		distrbt_var->theta_dt.dt[irdf].vnormal += distrbt_thread_var->theta_dt.dt[irdf].vnormal;
	}
}

void normalize_rdf(DIPOLE_DISTRIBT_VAR *rdf) {
	int irdf, i;
	for (irdf = 0; irdf < rdf->mu_dt.ndt; irdf++) {
		rdf->mu_dt.dt[irdf].eb_raw();
		rdf->mu_dt.dt[irdf].y /= rdf->mu_dt.dt[irdf].vnormal;
		rdf->mu_dt.dt[irdf].eb /= rdf->mu_dt.dt[irdf].vnormal;
	}
	for (irdf = 0; irdf < rdf->theta_dt.ndt; irdf++) {
		rdf->theta_dt.dt[irdf].eb_raw();
		rdf->theta_dt.dt[irdf].y /= rdf->theta_dt.dt[irdf].vnormal;
		rdf->theta_dt.dt[irdf].eb /= rdf->theta_dt.dt[irdf].vnormal;
	}
}

void save_rdf(DIPOLE_DISTRIBT_VAR *rdf_var, char *title, char *msg_title) {
	char ofname[256] = "\0", msg[256] = "\0";
	ofstream *out = NULL;

	int i, n;
	for (i = 0; i < rdf_var->mu_dt.ndt; i++) {
		sprintf(ofname, "%s-%d-mu.dat", title, i);
		out = new ofstream;
		out->open(ofname);
		rdf_var->mu_dt.dt[i].y /= rdf_var->mu_dt.dt[i].dx;
		rdf_var->mu_dt.dt[i].eb /= rdf_var->mu_dt.dt[i].dx;
		rdf_var->mu_dt.dt[i].save_x_y_eb(*out);
		out->close(); delete out; out = NULL;
		sprintf(msg, "%s # %d -- |mu| distrbt ==> %s", msg_title, i, ofname);
		show_log(msg, true);
	}
	for (i = 0; i < rdf_var->theta_dt.ndt; i++) {
		sprintf(ofname, "%s-%d-theta.dat", title, i);
		out = new ofstream;
		out->open(ofname);
		rdf_var->theta_dt.dt[i].dx *= fRadian2Angle;
		rdf_var->theta_dt.dt[i].xf *= fRadian2Angle;
		rdf_var->theta_dt.dt[i].xt *= fRadian2Angle;
		rdf_var->theta_dt.dt[i].x *= fRadian2Angle;
		rdf_var->theta_dt.dt[i].y /= rdf_var->theta_dt.dt[i].dx;
		rdf_var->theta_dt.dt[i].eb /= rdf_var->theta_dt.dt[i].dx;
		rdf_var->theta_dt.dt[i].save_x_y_eb(*out);
		out->close(); delete out; out = NULL;
		sprintf(msg, "%s # %d -- theta distrbt ==> %s", msg_title, i, ofname);
		show_log(msg, true);
	}
	return;
}

bool read_dipole_distrbt_pars(char *parfname, DIPOLE_DISTRIBT_VAR** distrbt_var, int &ndistrbt) {
	if (*distrbt_var != NULL) {delete[] *distrbt_var; *distrbt_var = NULL;} 
	ndistrbt = 0;

#define RELEASE_DISTRBT if (*distrbt_var != NULL) {delete[] *distrbt_var; *distrbt_var = NULL;} ndistrbt = 0;

	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";
	float xf = 0, xt = 0;
	int npts = 0, npts_theta = 0, i;

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	if (!search(in, "[DIPOLE-DISTRBT]", buffer)) {
		sprintf(msg, "can not find [DIPOLE-DISTRBT] in %s", parfname); show_msg(msg); 
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &ndistrbt) == 1 && ndistrbt > 0) {
		*distrbt_var = new DIPOLE_DISTRIBT_VAR[ndistrbt];
	}
	else {
		sprintf(msg, "number of distribution has to be > 0. see %s", buffer); 
		show_msg(msg, true); in.close(); return false;
	}
	for (i = 0; i < ndistrbt; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d %f %f %d %d", &((*distrbt_var)[i].cm_uid), &((*distrbt_var)[i].nm_uid), &xf, &xt, &npts, &npts_theta) != 6 || xt <= xf || npts <= 0 || npts_theta < 0) {
			sprintf(msg, "DIPOLE-DISTRBT format : [center monomer uid] [neighbor-monomer uid] [|mu| -- from] [|mu| -- to] [|mu| -- npts] [theta -- npts], see %s", buffer); show_msg(msg); 
			RELEASE_DISTRBT in.close(); return false;
		}
		(*distrbt_var)[i].nth_neighbor = 0;
		(*distrbt_var)[i].mu_dt.set(1);
		(*distrbt_var)[i].theta_dt.set(1);
		(*distrbt_var)[i].cm_uid += _USER_MONOMER_UID;
		(*distrbt_var)[i].nm_uid += _USER_MONOMER_UID;
		(*distrbt_var)[i].mu_dt.dt[0].set_distribt_x(xf, xt, npts);
		(*distrbt_var)[i].theta_dt.dt[0].set_distribt_x(0, PI, npts_theta);
	}
	in.close(); return true;
#undef RELEASE_DISTRBT
}

} // end of namespace _dipole_distribution_

extern bool ReadMoleculeStruct(ifstream &in, MMOLECULE &mm);
extern bool serial_multi_mol_md_proc(long isave);

extern "C" bool single_mm_dipole_distrbt_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;
	using namespace _dipole_distribution_;

	//extern void single_mm_cal_Position(MM_VM<MMOLECULE, CLUSTER>& mvm);

	::COMMAND_MD_STOP = false;
	mlog.init("dipole_distrbt.log", false);
	VM_VAR<MMOLECULE> vm_var;
	DIPOLE_DISTRIBT_VAR *dipole_distrbt_var = NULL;
	int ndistrbt = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_dipole_distrbt_pars(parfname, &dipole_distrbt_var, ndistrbt) || dipole_distrbt_var == NULL) return false;

	rdf_distribt_single_mm_multi_md<MMOLECULE, CLUSTER, VM_VAR<MMOLECULE>, DIPOLE_DISTRIBT_VAR, DIPOLE_DISTRBT_THREAD_VAR<MMOLECULE, CLUSTER> >(
		parfname, (void*)(&ReadMoleculeStruct), (void*)(&_MM_rdf_::single_mm_cal_Position), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		dipole_distrbt_var, ndistrbt, (void*)(&init_dipole_distrbt_thread_var), (void*)(&single_mm_dipole_distrbt_func),
		(void*)(&_dipole_distribution_::raw_accumulate_rdf), (void*)(&_dipole_distribution_::normalize_rdf), (void*)(&_dipole_distribution_::save_rdf));

	mlog.close();

	if (dipole_distrbt_var != NULL) {delete[] dipole_distrbt_var; dipole_distrbt_var = NULL;} ndistrbt = 0;
	return true;
}

extern "C" bool multi_mm_dipole_distrbt_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;
	using namespace _dipole_distribution_;

	::COMMAND_MD_STOP = false;
	mlog.init("dipole_distrbt.log", false);
	VM_VAR<MMOLECULE> vm_var;
	DIPOLE_DISTRIBT_VAR *dipole_distrbt_var = NULL;
	int ndistrbt = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_dipole_distrbt_pars(parfname, &dipole_distrbt_var, ndistrbt) || dipole_distrbt_var == NULL) return false;

	rdf_distribt_multi_mm_multi_md<MMOL_MD_CELL, MMOLECULE, CLUSTER, VM_VAR<MMOLECULE>, DIPOLE_DISTRIBT_VAR, DIPOLE_DISTRBT_THREAD_VAR<MMOLECULE, CLUSTER> >(
		parfname, ::mm_md_cell, (void*)(&serial_multi_mol_md_proc), (void*)(&multi_mm_cal_Position_NonCell), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		dipole_distrbt_var, ndistrbt, (void*)(&init_dipole_distrbt_thread_var), (void*)(&multi_mm_dipole_distrbt_func),
		(void*)(&_dipole_distribution_::raw_accumulate_rdf), (void*)(&_dipole_distribution_::normalize_rdf), (void*)(&_dipole_distribution_::save_rdf));

	mlog.close();

	if (dipole_distrbt_var != NULL) {delete[] dipole_distrbt_var; dipole_distrbt_var = NULL;} ndistrbt = 0;
	return true;
}
