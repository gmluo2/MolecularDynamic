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

#include "distribution.h"

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

using namespace _coarse_grain_;

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

bool get_monomer_uid(char *monomer, int &uid) {
	char title[250] = "\0", buffer[250] = "\0";
	sprintf(title, "[%s]", monomer);
	ifstream in;
	in.open(struct_fname); // this file should exist, because the macromolecule is constructed with this file
	if (!in.is_open()) {
		sprintf(buffer, "can not open %s to get monomer structure", struct_fname); show_msg(buffer); return false;
	}
	if (!search(in, title, buffer)) {
		sprintf(buffer, "can not find monomer %s", monomer); show_msg(buffer); 
		in.close(); return false;
	}
	if (!search(in, "UID", buffer, "[") || sscanf(buffer, "%s %d", title, &uid) != 2) {
		sprintf(buffer, "do not find UID of monomer %s", monomer); show_msg(buffer); 
		in.close(); return false;
	}
	in.close();
	return true;
}

bool get_cluster_uid(char *cluster, int &uid) {
	char title[250] = "\0", buffer[250] = "\0";
	sprintf(title, "[%s]", cluster);
	ifstream in;
	in.open(struct_fname); // this file should exist, because the macromolecule is constructed with this file
	if (!in.is_open()) {
		sprintf(buffer, "can not open %s to get cluster structure", struct_fname); show_msg(buffer); return false;
	}
	if (!search(in, title, buffer)) {
		sprintf(buffer, "can not find cluster %s", cluster); show_msg(buffer); 
		in.close(); return false;
	}
	if (!search(in, "INDEX", buffer, "[")) {
		sprintf(buffer, "do not find INDEX of cluster %s", cluster); show_msg(buffer); 
		in.close(); return false;
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d", &uid) != 1) {
			sprintf(buffer, "index of cluster [%s] is not given after INDEX", cluster); show_msg(buffer); 
			in.close(); return false;
		}
	}
	in.close();
	return true;
}

bool read_given_id(char *fname, char* title, ARRAY<short>& id_cluster) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return false;

	char buffer[256] = "\0", *bf = NULL;
	char seperator[5] = ",; \t";
	int uid = 0, *p_uid = NULL, n = 0;
	CHAIN<int> *mch = NULL, *tail = NULL;
	if (!search(in, title, buffer)) {in.close(); return false;}
	while (!in.eof()) {
		in.getline(buffer, 250); bf = buffer;
		n = 0;
		while (strlen(bf) > 0) {
			if (sscanf(bf, "%d", &uid) != 1) break;
			p_uid = new int; *p_uid = uid;
			if (mch == NULL) {
				mch = new CHAIN<int>; tail = mch; mch->p = p_uid;
			}
			else {
				tail->next = new CHAIN<int>; tail = tail->next; tail->p = p_uid;
			}
			n++;
			if (!NextSection(&bf, seperator)) break;
		}
		if (n == 0) break;
	}
	in.close();
	n = number<int>(mch);
	if (n == 0) return false;
	id_cluster.SetArray(n);
	tail = mch;
	for (n = 0; n < id_cluster.n; n++) {
		if (tail == NULL) {id_cluster.SetArray(0); return false;}
		id_cluster.m[n] = *(tail->p);
		tail = tail->next;
	}
	release_chain<int>(&mch, true); tail = NULL;
	return true;
}

bool read_given_id(char *fname, char* title, ARRAY<short>& id_cluster, bool (*get_uid)(char*, int&), char *end_str) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return false;

	int nlen_end = 0;
	if (end_str != NULL) nlen_end = strlen(end_str);

	char buffer[256] = "\0", *bf = NULL;
	char seperator[5] = ",; \t";
	char unit[256] = "\0", start[256] = "\0";
	int uid = 0, *p_uid = NULL, n = 0;
	CHAIN<int> *mch = NULL, *tail = NULL;
	if (!search(in, title, buffer)) {in.close(); return false;}
	while (!in.eof()) {
		in.getline(buffer, 250); bf = buffer;
		n = 0;
		while (strlen(bf) > 0) {
			if (sscanf(bf, "%s", unit) != 1) break;
			if (nlen_end > 0) {
				pickup_str(unit, 0, nlen_end, start);
				if (strcmp(start, end_str) == 0) break;
			}
			if (!get_uid(unit, uid)) break;
			p_uid = new int; *p_uid = uid;
			if (mch == NULL) {
				mch = new CHAIN<int>; tail = mch; mch->p = p_uid;
			}
			else {
				tail->next = new CHAIN<int>; tail = tail->next; tail->p = p_uid;
			}
			n++;
			if (!NextSection(&bf, seperator)) break;
		}
		if (n == 0) break;
	}
	in.close();
	n = number<int>(mch);
	if (n == 0) return false;
	id_cluster.SetArray(n);
	tail = mch;
	for (n = 0; n < id_cluster.n; n++) {
		if (tail == NULL) {id_cluster.SetArray(0); return false;}
		id_cluster.m[n] = *(tail->p);
		tail = tail->next;
	}
	release_chain<int>(&mch, true); tail = NULL;
	return true;
}

namespace _cluster_distribution_ {

/*****************************************************
********  DIPOLE DISTRIBUTION OF MONOMER      ********
******************************************************/

void monomer_dipole_distribution(MMOLECULE &mm, bool group, int monomer_uid, bool limit_range, float mu1, float mu2, DISTRIBT<float, long> *dt) {
	int na = 0, nc = 0, n_monomer = 0, npt = 0;
	_GROUP<CLUSTER> *monomer = NULL;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;

	VECTOR3 mu, mu_total, mu0;
	float c = 0, c_total = 0, x = 0, mu_abs = 0;

	if (group) { // group of monomer
		for (n_monomer = 0; n_monomer < mm.nMonomer; n_monomer++) {
			monomer = mm.monomer + n_monomer;
			if (monomer_uid > 0 && monomer->uid != monomer_uid) continue;
			for (nc = 0; nc < monomer->nUnit; nc++) {
				pc = monomer->u[nc];
				for (na = 0; na < pc->nAtoms; na++) {
					patom = pc->atom + na;
					mu.v[0] += patom->c0 * patom->r.v[0];
					mu.v[1] += patom->c0 * patom->r.v[1];
					mu.v[2] += patom->c0 * patom->r.v[2];
				}
			}
		}
		mu_abs = mu.abs();

		if (limit_range) {
			if (mu_abs > mu2 || mu_abs < mu1) return;
		}

		//DISTRIBT_INDX(dt[MU_ABS], mu_abs, npt) 
		//if (LInRange(npt, 0, dt[MU_ABS].npts)) dt[MU_ABS].y[npt]++;
		dt[MU_ABS].sample(mu_abs, true);
		//DISTRIBT_INDX(dt[MU_X], mu.v[0], npt)
		//if (LInRange(npt, 0, dt[MU_X].npts)) dt[MU_X].y[npt]++;
		dt[MU_X].sample(mu.v[0], true);
		//DISTRIBT_INDX(dt[MU_Y], mu.v[1], npt)
		//if (LInRange(npt, 0, dt[MU_Y].npts)) dt[MU_Y].y[npt]++;
		dt[MU_Y].sample(mu.v[1], true);
		//DISTRIBT_INDX(dt[MU_Z], mu.v[2], npt)
		//if (LInRange(npt, 0, dt[MU_Z].npts)) dt[MU_Z].y[npt]++;
		dt[MU_Z].sample(mu.v[2], true);
		x = sqrt(mu.v[0] * mu.v[0] + mu.v[1] * mu.v[1]);
		//DISTRIBT_INDX(dt[MU_XY], x, npt)
		//if (LInRange(npt, 0, dt[MU_XY].npts)) dt[MU_XY].y[npt]++;
		dt[MU_XY].sample(x, true);
		x = mu.theta();
		//DISTRIBT_INDX(dt[MU_THETA], x, npt)
		//if (LInRange(npt, 0, dt[MU_THETA].npts)) dt[MU_THETA].y[npt]++;
		dt[MU_THETA].sample(x, true);
		x = mu.phi();
		//DISTRIBT_INDX(dt[MU_OMEGA], x, npt)
		//if (LInRange(npt, 0, dt[MU_OMEGA].npts)) dt[MU_OMEGA].y[npt]++;
		dt[MU_OMEGA].sample(x, true);
	}
	else {
		for (n_monomer = 0; n_monomer < mm.nMonomer; n_monomer++) {
			monomer = mm.monomer + n_monomer;
			if (monomer_uid > 0 && monomer->uid != monomer_uid) continue;
			V3zero(mu) c = 0;
			for (nc = 0; nc < monomer->nUnit; nc++) {
				pc = monomer->u[nc];
				for (na = 0; na < pc->nAtoms; na++) {
					patom = pc->atom + na;
					c += patom->c0;
					TIME3(patom->r, patom->c0, mu0)
					V3plusV3(mu, mu0, mu)
				}
			}
			if (limit_range) {
				mu_abs = mu.abs();
				if (mu_abs > mu2 || mu_abs < mu1) continue;
			}
			V3plusV3(mu_total, mu, mu_total)
			c_total += c;
		
			//DISTRIBT_INDX(dt[MU_ABS], mu_abs, npt)
			//if (LInRange(npt, 0, dt[MU_ABS].N)) dt[MU_ABS].y[npt]++;
			dt[MU_ABS].sample(mu_abs, true);
			//DISTRIBT_INDX(dt[MU_X], mu.v[0], npt)
			//if (LInRange(npt, 0, dt[MU_X].N)) dt[MU_X].y[npt]++;
			dt[MU_X].sample(mu.v[0], true);
			//DISTRIBT_INDX(dt[MU_Y], mu.v[1], npt)
			//if (LInRange(npt, 0, dt[MU_Y].N)) dt[MU_Y].y[npt]++;
			dt[MU_Y].sample(mu.v[1], true);
			//DISTRIBT_INDX(dt[MU_Z], mu.v[2], npt)
			//if (LInRange(npt, 0, dt[MU_Z].N)) dt[MU_Z].y[npt]++;
			dt[MU_Z].sample(mu.v[2], true);
			x = sqrt(mu.v[0] * mu.v[0] + mu.v[1] * mu.v[1]);
			//DISTRIBT_INDX(dt[MU_XY], x, npt)
			//if (LInRange(npt, 0, dt[MU_XY].N)) dt[MU_XY].y[npt]++;
			dt[MU_XY].sample(x, true);
			x = mu.theta();
			//DISTRIBT_INDX(dt[MU_THETA], x, npt)
			//if (LInRange(npt, 0, dt[MU_THETA].N)) dt[MU_THETA].y[npt]++;
			dt[MU_THETA].sample(x, true);
			x = mu.phi();
			//DISTRIBT_INDX(dt[MU_OMEGA], x, npt)
			//if (LInRange(npt, 0, dt[MU_OMEGA].N)) dt[MU_OMEGA].y[npt]++;
			dt[MU_OMEGA].sample(x, true);
		}
	}
}

#if _SYS_ == _WINDOWS_SYS_
int monomer_dipole_distribution_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* monomer_dipole_distribution_thread(void *void_par) { // thread-function
#endif
	DISTRIBT_MU_THREAD_VARS *par = (DISTRIBT_MU_THREAD_VARS*)void_par;
	(par->func)(*(par->mm), par->monomer_group, par->monomer_uid, par->limit_range, par->mu1, par->mu2, par->dt);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

bool monomer_dipole_distribution_stat(char *mdpar_fname, bool group, int monomer_uid, int nf, int nt, int nstep, bool limit_range, float r1, float r2, DISTRIBT<float, float> *dt, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", mol_template[256] = "\0", buffer[256] = "\0";
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

	int nThreads = 0, iconf = 0, nconf = 0, iThread, i, j, k;

	MMOLECULE mm[MAX_THREADS];
	DISTRIBT<float, long> ndt[MAX_THREADS][N_MU_DT];
	DISTRIBT_MU_THREAD_VARS thread_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		if (!ConstructSingleChainPolymerFromParamFile(mm[i], mol_template)) return false;
		for (j = 0; j < N_MU_DT; j++) ndt[i][j].set_distribt_x(dt[j].xf, dt[j].xt, dt[j].x.N - 1);
		thread_var[i].set_func(monomer_dipole_distribution);
		thread_var[i].set_pars(mm + i, group, monomer_uid, ndt[i], i);
		thread_var[i].set_range(limit_range, r1, r2);
	}

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
			ReadMoleculeStruct(*(rmdproc.in), mm[i]);
			sprintf(msg, " %d ", iconf); strcat(msg1, msg);

			iThread = i;
		#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(thread_var[iThread].hMMThread);
			AfxBeginThread((AFX_THREADPROC)monomer_dipole_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
		#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(thread_var[iThread].thread), NULL, &monomer_dipole_distribution_thread, (void *)(thread_var + iThread));
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
	
	for (i = 0; i < N_MU_DT; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(dt[i], ndt[k][i]);
	}
	nconfs = nconf;

	return true;
}


bool MultiMol_monomer_dipole_distribution_stat(char *mdpar_fname, ARRAY<short> &mol_id, bool group, int monomer_uid, int nf, int nt, int nstep, bool limit_range, float r1, float r2, DISTRIBT<float, float> *dt, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", mol_template[256] = "\0", buffer[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
	int MD_MODE = 0;
	if (sscanf(msg, "%d", &MD_MODE) != 1 || (MD_MODE != MULTI_MM_FREE_CELL_MD && MD_MODE != MULTI_MM_PERIDIC_CELL_MD)) { // MultiMol ?
		sprintf(msg, "md is not working on multi-molecule");
		show_msg(msg); return false;
	}

	set_md_mode(MD_MODE);
	set_job(GET_MULTI_MM_MD_PROC);

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int nThreads = 0, iThread, imol = 0, i, j, k;
	long iconf = 0, nconf = 0;

	DISTRIBT<float, long> ndt[MAX_THREADS][N_MU_DT];
	DISTRIBT_MU_THREAD_VARS thread_var[MAX_THREADS];

	for (i = 0; i < MAX_THREADS; i++) {
		for (j = 0; j < N_MU_DT; j++) ndt[i][j].set_distribt_x(dt[j].xf, dt[j].xt, dt[j].x.N - 1);
		thread_var[i].set_func(monomer_dipole_distribution);
		thread_var[i].set_pars(NULL, group, monomer_uid, ndt[i], i);
		thread_var[i].set_range(limit_range, r1, r2);
	}

	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		nThreads = 0;
		if (!serial_multi_mol_md_proc(iconf)) break;
		sprintf(msg, " %d in %s are read", iconf, mdpar_fname); show_log(msg, true);

		iThread = 0; nThreads = 0; show_log("analyzing molecules ", false);
		for (imol = 0; imol < ::mm_md_cell.mm.n; imol++) { // working on macro-molecules only
			if (mol_id.n > 0 && !mol_id.exist(imol)) {
				if (imol != ::mm_md_cell.mm.n - 1) continue; // this molecule is ignored
			}
			else {
				thread_var[iThread].mm = ::mm_md_cell.mm.m + imol;
			#if _SYS_ == _WINDOWS_SYS_ 
				ResetEvent(thread_var[iThread].hMMThread);
				AfxBeginThread((AFX_THREADPROC)monomer_dipole_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
			#elif _SYS_ == _LINUX_SYS_
				pthread_create(&(thread_var[iThread].thread), NULL, &monomer_dipole_distribution_thread, (void *)(thread_var + iThread));
			#endif
				sprintf(msg, "%d ", imol); show_log(msg, false);
				iThread++; nThreads++;
			}

			if (iThread == MAX_THREADS || imol == ::mm_md_cell.mm.n - 1) { // reaching maximum threads or all molecules are run
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
	
	for (i = 0; i < N_MU_DT; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(dt[i], ndt[k][i]);
	}
	nconfs = nconf;

	return true;
}




/*****************************************************
********  SPATIAL DISTRIBUTION OF CLUSTER     ********
******************************************************/

void cluster_spatial_distribution(ARRAY<BASIC_CLUSTER*> &bc, ARRAY<short>& c_uid, int axis, bool bLimitRange, DISTRIBT<float, long> *dt) {
	int ic = 0, npt = 0;
	BASIC_CLUSTER *pc = NULL;

	double vpos;

	for (ic = 0; ic < bc.n; ic++) {
		pc = bc.m[ic];
		if (!c_uid.exist(pc->cID)) continue;
		switch (axis) {
		case R_X:
			vpos = pc->vp->v[0];
			break;
		case R_Y:
			vpos = pc->vp->v[1];
			break;
		case R_Z:
			vpos = pc->vp->v[2];
			break;
		default:
			return;
		}
		dt[0].sample(vpos, bLimitRange);
	}
}

#if _SYS_ == _WINDOWS_SYS_
int cluster_spatial_distribution_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* cluster_spatial_distribution_thread(void *void_par) { // thread-function
#endif
	SPATIAL_DISTRIBT_THREAD_VARS *par = (SPATIAL_DISTRIBT_THREAD_VARS*)void_par;
	(par->func)(*(par->bc), par->c_uid, par->axis, par->bLimitRange, par->dt);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

bool cluster_spatial_distribution_stat(char *mdpar_fname, int nf, int nt, int nstep, ARRAY<short> &c_uid, int axis, bool bLimitRange, DISTRIBT<float, float> *dt, int &nconfs) {
	char msg[256] = "\0", msg1[256] = "\0", buffer[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
	int MD_MODE = 0;
	if (sscanf(msg, "%d", &MD_MODE) != 1 || (MD_MODE != MULTI_MM_PERIDIC_CELL_MD)) { // MultiMol ?
		sprintf(msg, "md is not working on multi-molecule");
		show_msg(msg); return false;
	}

	set_md_mode(MD_MODE);
	set_job(GET_MULTI_MM_MD_PROC);

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	::mm_md_cell.init_groups();
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	long iconf = 0, nconf = 0;
	DISTRIBT<float, long> ndt;
	ndt.set_distribt_x(dt[0].xf, dt[0].xt, dt[0].x.N - 1);

	iconf = nf;
	nconf = 0;
	::mlog.bReplaceLog = true;
	::mlog.record_pos();
	while (1) {	
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		if (!serial_multi_mol_md_proc(iconf)) break;
		sprintf(msg, " %d in %s are read", iconf, mdpar_fname); show_log1(msg, false, true);

		cluster_spatial_distribution(::mm_md_cell.bcluster, c_uid, axis, bLimitRange, &ndt);
		show_log1("and analyzed", true, false);

		iconf += nstep;
		nconf++;
	}
	clear_read_mdproc();
	
	raw_distribt_combine<float, float, float, long>(dt[0], ndt);
	nconfs = nconf;
	return true;
}


} // end of namespace _cluster_distribution_

extern "C" bool dipole_distribution(char *parfname) {
	using namespace _cluster_distribution_;
	::COMMAND_MD_STOP = false;
	mlog.init("dipole.log", false);
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	float xf = 0, xt = 0;
	int npts = 0, mode = 0;

	DISTRIBT<float, float> dt[N_MU_DT];
	
	if (!search(in, "[DISTRIBUTION-RANGE]", buffer)) {sprintf(msg, "can not find [DISTRIBUTION-RANGE] in %s", parfname); show_msg(msg); return false;}
	// |mu|
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %d", &xf, &xt, &npts) != 3 || xt <= xf || npts <= 0) {
		sprintf(msg, "RANGE format : [from] [to] [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[MU_ABS].set_distribt_x(xf, xt, npts);
	// mu_x
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %d", &xf, &xt, &npts) != 3 || xt <= xf || npts <= 0) {
		sprintf(msg, "RANGE format : [from] [to] [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[MU_X].set_distribt_x(xf, xt, npts);
	// mu_y
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %d", &xf, &xt, &npts) != 3 || xt <= xf || npts <= 0) {
		sprintf(msg, "RANGE format : [from] [to] [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[MU_Y].set_distribt_x(xf, xt, npts);
	// mu_z
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %d", &xf, &xt, &npts) != 3 || xt <= xf || npts <= 0) {
		sprintf(msg, "RANGE format : [from] [to] [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[MU_Z].set_distribt_x(xf, xt, npts);
	// mu_xy
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %d", &xf, &xt, &npts) != 3 || xt <= xf || npts <= 0) {
		sprintf(msg, "RANGE format : [from] [to] [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[MU_XY].set_distribt_x(xf, xt, npts);
	// mu_theta
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &npts) != 1 || npts <= 0) {
		sprintf(msg, "theta : [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[MU_THETA].set_distribt_x(0, fPI, npts);
	// mu_omega
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &npts) != 1 || npts <= 0) {
		sprintf(msg, "omega : [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[MU_OMEGA].set_distribt_x(0, fPI + fPI, npts);
	in.close();

	bool group = false;
	int igroup = 0;
	if (!read_item(parfname, "[MONOMER_GROUP]", buffer)) {
		sprintf(msg, "[MONOMER_GROUP] is not defined in %s", parfname); show_log(msg, true); 
		sprintf(msg, "dipole moment of each monomer is statistical analyzed", parfname); show_log(msg, true); 
	}
	else {
		if (sscanf(buffer, "%d", &igroup) == 1 && igroup >= 1) group = true;
		else group = false;
		if (group) sprintf(msg, "statistical analysis of the dipole moment of monomer group");
		else sprintf(msg, "statistical analysis of the dipole moment of each monomer"); 
		show_log(msg, true); 
	}

	/*
	bool limit_range = false;
	int nlimit = 0;
	float r1 = 0, r2 = 0;
	if (!read_item(parfname, "[LIMIT]", buffer)) {sprintf(msg, "can not find [LIMIT] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%d %f %f", &nlimit, &r1, &r2) != 3 || r1 < r2) {
		sprintf(msg, "LIMIT format : [disable/enable -- 0/1] [from] [to], see %s", buffer); show_msg(msg); return false;
	}
	limit_range = (nlimit > 0 ? true : false);
	*/

	char monomer[256] = "\0";
	int monomer_uid = 0;
	if (!read_item(parfname, "[OUTPUT]", buffer)) {sprintf(msg, "can not find [OUTPUT] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%s %s", monomer, title) != 2) {sprintf(msg, "can not get monomer to analyze, and the output filename after [OUTPUT]"); show_msg(msg); return false;}
	if (strcmp(monomer, "NULL") == 0 || strcmp(monomer, "null") == 0) monomer_uid = -1; // the whole molecule
	else if (!get_monomer_uid(monomer, monomer_uid)) {
		sprintf(msg, "unknown monomer %s", monomer); show_msg(msg); return false;
	}

	ARRAY<short> mol_id;
	read_given_id(parfname, "[MOLECULES]", mol_id);

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

		if (MD_MODE == SINGLE_MM_FREE_CELL_MD) status = monomer_dipole_distribution_stat(mdpar_fname, group, monomer_uid, nf, nt, nstep, true, dt[MU_ABS].xf, dt[MU_ABS].xt, dt, nconf);
		else if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) 
			status = MultiMol_monomer_dipole_distribution_stat(mdpar_fname, mol_id, group, monomer_uid, nf, nt, nstep, true, dt[MU_ABS].xf, dt[MU_ABS].xt, dt, nconf);
		else {
			sprintf(msg, "UNKNOWN MD in %s", mdpar_fname); show_msg(msg, true); nconf = 0;
		}
		nconfs += nconf;
		n++;
	}
	in_cfg.close();

	ofstream *out = NULL;

	sprintf(ofname, "%s-mu-abs.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[MU_ABS].eb_raw();
	dt[MU_ABS].y /= nconfs;
	dt[MU_ABS].eb /= nconfs;
	dt[MU_ABS].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(msg, "|mu| ==> %s", ofname);

	sprintf(ofname, "%s-mu-x.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[MU_X].eb_raw();
	dt[MU_X].y /= nconfs;
	dt[MU_X].eb /= nconfs;
	dt[MU_X].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(buffer, "\nmu_x ==> %s", ofname); strcat(msg, buffer);

	sprintf(ofname, "%s-mu-y.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[MU_Y].eb_raw();
	dt[MU_Y].y /= nconfs;
	dt[MU_Y].eb /= nconfs;
	dt[MU_Y].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(buffer, "\nmu_y ==> %s", ofname);  strcat(msg, buffer);

	sprintf(ofname, "%s-mu-z.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[MU_Z].eb_raw();
	dt[MU_Z].y /= nconfs;
	dt[MU_Z].eb /= nconfs;
	dt[MU_Z].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(buffer, "\nmu_z ==> %s", ofname);  strcat(msg, buffer);

	sprintf(ofname, "%s-mu-xy.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[MU_XY].eb_raw();
	dt[MU_XY].y /= nconfs;
	dt[MU_XY].eb /= nconfs;
	dt[MU_XY].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(buffer, "\nmu_xy ==> %s", ofname); strcat(msg, buffer);

	sprintf(ofname, "%s-mu-theta.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[MU_THETA].eb_raw();
	dt[MU_THETA].y /= nconfs;
	dt[MU_THETA].eb /= nconfs;
	dt[MU_THETA].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(buffer, "\nmu-theta ==> %s", ofname); strcat(msg, buffer);

	sprintf(ofname, "%s-mu-omega.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[MU_OMEGA].eb_raw();
	dt[MU_OMEGA].y /= nconfs;
	dt[MU_OMEGA].eb /= nconfs;
	dt[MU_OMEGA].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(buffer, "\nmu-omega ==> %s", ofname); strcat(msg, buffer);
	show_msg(msg);

	mlog.close();
	return true;
}


extern "C" bool cluster_spatial_distribution(char *parfname) {
	using namespace _cluster_distribution_;
	::COMMAND_MD_STOP = false;
	mlog.init("spatial_distrbt.log", false, true);
	mlog.nMaxLines = 20;
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	float xf = 0, xt = 0;
	int npts = 0, mode = 0, axis = 0;

	DISTRIBT<float, float> dt;

	if (!search(in, "[AXIS]", buffer)) {
		sprintf(msg, "can not find [AXIS] in %s", parfname); show_msg(msg); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%s", msg) != 1) {
		sprintf(msg, "AXIS : x/y/z, see %s", buffer); show_msg(msg); in.close(); return false;
	}
	Upper(msg);
	if (strcmp(msg, "X") == 0) axis = R_X;
	else if (strcmp(msg, "Y") == 0) axis = R_Y;
	else if (strcmp(msg, "Z") == 0) axis = R_Z;
	else {
		sprintf(msg, "AXIS : x/y/z. Unknown : %s", buffer); show_msg(msg); in.close(); return false;
	}
	
	if (!search(in, "[DISTRIBUTION-RANGE]", buffer)) {sprintf(msg, "can not find [DISTRIBUTION-RANGE] in %s", parfname); show_msg(msg); return false;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %d", &xf, &xt, &npts) != 3 || xt <= xf || npts <= 0) {
		sprintf(msg, "RANGE format : [from] [to] [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt.set_distribt_x(xf, xt, npts);
	in.close();

	bool bLimitRange = true;
	/*
	int nlimit = 0;
	float r1 = 0, r2 = 0;
	if (!read_item(parfname, "[LIMIT]", buffer)) {sprintf(msg, "can not find [LIMIT] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%d %f %f", &nlimit, &r1, &r2) != 3 || r1 < r2) {
		sprintf(msg, "LIMIT format : [disable/enable -- 0/1] [from] [to], see %s", buffer); show_msg(msg); return false;
	}
	limit_range = (nlimit > 0 ? true : false);
	*/

	if (!read_item(parfname, "[OUTPUT]", buffer)) {sprintf(msg, "can not find [OUTPUT] in %s", parfname); show_msg(msg); return false;}
	if (sscanf(buffer, "%s", ofname) != 1) {sprintf(msg, "can not get the output filename after [OUTPUT]"); show_msg(msg); return false;}

	ARRAY<short> c_uid;
	//read_given_id(parfname, "[CLUSTERS]", c_uid);
	read_given_id(parfname, "[CLUSTERS]", c_uid, get_cluster_uid, "/");

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
		mlog.go_end();
		mlog.record_pos();

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

		if (MD_MODE == MULTI_MM_PERIDIC_CELL_MD) 
			status = cluster_spatial_distribution_stat(mdpar_fname, nf, nt, nstep, c_uid, axis, bLimitRange, &dt, nconf);
		else {
			sprintf(msg, "only periodical MD is acceptable. check %s", mdpar_fname); show_msg(msg, true); nconf = 0;
		}
		nconfs += nconf;
		n++;
	}
	in_cfg.close();

	ofstream *out = NULL;

	out = new ofstream;
	out->open(ofname);
	dt.eb_raw();
	dt.y /= nconfs;
	dt.eb /= nconfs;
	dt.save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(msg, "spatial distribution ==> %s", ofname);
	show_msg(msg);

	mlog.close();
	return true;
}

namespace _cluster_distribution_ {
/*****************************************************
********       ORIENTATION OF MONOMER         ********
******************************************************/

void monomer_orientation_distribution(MMOLECULE &mm, int monomer_uid, VECT_IN_MONOMER* mv, DISTRIBT<float, long> *dt) {
	int na = 0, n_monomer = 0, npt = 0;
	_GROUP<CLUSTER> *monomer = NULL;
	MATOM *patom1 = NULL, *patom2 = NULL;

	VECTOR3 v_monomer;
	float theta = 0, omega = 0;

	for (n_monomer = 0; n_monomer < mm.nMonomer; n_monomer++) {
		monomer = mm.monomer + n_monomer;
		if (monomer->uid != monomer_uid) continue;
		if (monomer->nUnit <= mv->atom[0].nc || monomer->nUnit <= mv->atom[1].nc) continue; // this is not supposed happen, but to make program safe, we check it here

		patom1 = monomer->u[mv->atom[0].nc]->atom + mv->atom[0].na;
		patom2 = monomer->u[mv->atom[1].nc]->atom + mv->atom[1].na;
		VECT3(patom1->r, patom2->r, v_monomer) // patom1 -> patom2
		theta = v_monomer.theta(); omega = v_monomer.phi();
		
		//DISTRIBT_INDX(dt[ORT_THETA], theta, npt)
		//if (LInRange(npt, 0, dt[ORT_THETA].N)) dt[ORT_THETA].y[npt]++;
		dt[ORT_THETA].sample(theta, true);
		//DISTRIBT_INDX(dt[ORT_OMEGA], omega, npt) 
		//if (LInRange(npt, 0, dt[ORT_OMEGA].N)) dt[ORT_OMEGA].y[npt]++;
		dt[ORT_OMEGA].sample(omega, true);
	}
}

#if _SYS_ == _WINDOWS_SYS_
int monomer_orientation_distribution_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* monomer_orientation_distribution_thread(void *void_par) { // thread-function
#endif
	DISTRIBT_ORT_THREAD_VARS *par = (DISTRIBT_ORT_THREAD_VARS*)void_par;
	(par->func)(*(par->mm), par->monomer_uid, &(par->mv), par->dt);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

bool monomer_orientation_distribution_stat(char *mdpar_fname, int monomer_uid, VECT_IN_MONOMER *mv, int nf, int nt, int nstep, DISTRIBT<float, float> *dt, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", mol_template[256] = "\0", buffer[256] = "\0";
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

	int nThreads = 0, iconf = 0, nconf = 0, iThread, i, j, k;

	MMOLECULE mm[MAX_THREADS];
	DISTRIBT<float, long> ndt[MAX_THREADS][N_ORT_DT];
	DISTRIBT_ORT_THREAD_VARS thread_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		if (!ConstructSingleChainPolymerFromParamFile(mm[i], mol_template)) return false;
		for (j = 0; j < N_ORT_DT; j++) ndt[i][j].set_distribt_x(dt[j].xf, dt[j].xt, dt[j].x.N - 1);
		thread_var[i].set_func(monomer_orientation_distribution);
		thread_var[i].set_pars(mm + i, monomer_uid, *mv, ndt[i], i);
	}

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
			ReadMoleculeStruct(*(rmdproc.in), mm[i]);
			sprintf(msg, " %d ", iconf); strcat(msg1, msg);
			
			iThread = i;
		#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(thread_var[iThread].hMMThread);
			AfxBeginThread((AFX_THREADPROC)monomer_orientation_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
		#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(thread_var[iThread].thread), NULL, &monomer_orientation_distribution_thread, (void *)(thread_var + iThread));
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
	
	for (i = 0; i < N_ORT_DT; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(dt[i], ndt[k][i]);
	}
	nconfs = nconf;
	return true;
}

bool MultiMol_monomer_orientation_distribution_stat(char *mdpar_fname, ARRAY<short> &mol_id, int monomer_uid, VECT_IN_MONOMER *mv, int nf, int nt, int nstep, DISTRIBT<float, float> *dt, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", buffer[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
	int MD_MODE = 0;
	if (sscanf(msg, "%d", &MD_MODE) != 1 || (MD_MODE != MULTI_MM_FREE_CELL_MD && MD_MODE != MULTI_MM_PERIDIC_CELL_MD)) { // MultiMol ?
		sprintf(msg, "md is not working on multi-molecule");
		show_msg(msg); return false;
	}

	set_md_mode(MD_MODE);
	set_job(GET_MULTI_MM_MD_PROC);

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int nThreads = 0, iThread, imol = 0, i, j, k;
	long iconf = 0, nconf = 0;

	DISTRIBT<float, long> ndt[MAX_THREADS][N_ORT_DT];
	DISTRIBT_ORT_THREAD_VARS thread_var[MAX_THREADS];

	for (i = 0; i < MAX_THREADS; i++) {
		for (j = 0; j < N_ORT_DT; j++) ndt[i][j].set_distribt_x(dt[j].xf, dt[j].xt, dt[j].x.N - 1);
		thread_var[i].set_func(monomer_orientation_distribution);
		thread_var[i].set_pars(NULL, monomer_uid, *mv, ndt[i], i);
	}

	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		nThreads = 0;
		if (!serial_multi_mol_md_proc(iconf)) break;
		sprintf(msg, " %d in %s are read", iconf, mdpar_fname); show_log(msg, true);

		iThread = 0; nThreads = 0; show_log("analyzing molecules ", false);
		for (imol = 0; imol < ::mm_md_cell.mm.n; imol++) { // working on macro-molecules only
			if (mol_id.n > 0 && !mol_id.exist(imol)) {
				if (imol != ::mm_md_cell.mm.n - 1) continue; // this molecule is ignored
			}
			else {
				thread_var[iThread].mm = ::mm_md_cell.mm.m + imol;
			#if _SYS_ == _WINDOWS_SYS_ 
				ResetEvent(thread_var[iThread].hMMThread);
				AfxBeginThread((AFX_THREADPROC)monomer_orientation_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
			#elif _SYS_ == _LINUX_SYS_
				pthread_create(&(thread_var[iThread].thread), NULL, &monomer_orientation_distribution_thread, (void *)(thread_var + iThread));
			#endif

				sprintf(msg, "%d ", imol); show_log(msg, false);
				iThread++; nThreads++;
			}

			if (iThread == MAX_THREADS || imol == ::mm_md_cell.mm.n - 1) { // reaching maximum threads or all molecules are run
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
	
	for (i = 0; i < N_ORT_DT; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(dt[i], ndt[k][i]);
	}
	nconfs = nconf;
	return true;
}

} // end of namespace _cluster_distribution_

extern "C" bool orientation_distribution(char *parfname) {
	using namespace _cluster_distribution_;
	::COMMAND_MD_STOP = false;
	mlog.init("orientation.log", false);
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	int npts = 0;
	int monomer_uid = 0;
	VECT_IN_MONOMER mv;
	char monomer[256] = "\0";

	DISTRIBT<float, float> dt[N_ORT_DT];

	if (!search(in, "[MONOMER]", buffer)) {sprintf(msg, "can not find [MONOMER] in %s", parfname); show_msg(msg); in.close(); return false;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%s", monomer) != 1) {
		sprintf(msg, "monomer is not given after [MONOMER] in %s", parfname); show_msg(msg); in.close(); return false;
	}
	if (!get_monomer_uid(monomer, monomer_uid)) {
		sprintf(msg, "unknown monomer %s", monomer); show_msg(msg); in.close(); return false;
	}
	if (!search(in, "[VECTOR]", buffer)) {sprintf(msg, "can not find [MONOMER] in %s", parfname); show_msg(msg); in.close(); return false;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d %d %d", &(mv.atom[0].nc), &(mv.atom[0].na), &(mv.atom[1].nc), &(mv.atom[1].na)) != 4) {
		sprintf(msg, "vector definition : [cluster1] [atom1] [cluster2] [atom2], see %s", buffer);
		show_msg(msg); in.close(); return false;
	}
	
	if (!search(in, "[DISTRIBUTION-RANGE]", buffer)) {sprintf(msg, "can not find [DISTRIBUTION-RANGE] in %s", parfname); show_msg(msg); in.close(); return false;}
	// theta
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &npts) != 1 || npts <= 0) {
		sprintf(msg, "RANGE format : [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[ORT_THETA].set_distribt_x(0, fPI, npts);
	// omega
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &npts) != 1 || npts <= 0) {
		sprintf(msg, "RANGE format : [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[ORT_OMEGA].set_distribt_x(0, fPI2, npts);
	in.close();

	ARRAY<short> mol_id;
	read_given_id(parfname, "[MOLECULES]", mol_id);

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

		if (MD_MODE == SINGLE_MM_FREE_CELL_MD) status = monomer_orientation_distribution_stat(mdpar_fname, monomer_uid, &mv, nf, nt, nstep, dt, nconf);
		else if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) 
			status = MultiMol_monomer_orientation_distribution_stat(mdpar_fname, mol_id, monomer_uid, &mv, nf, nt, nstep, dt, nconf);
		else {
			sprintf(msg, "UNKNOWN MD in %s", mdpar_fname); show_msg(msg, true); nconf = 0;
		}

		nconfs += nconf;
		n++;
	}
	in_cfg.close();

	ofstream *out = NULL;

	sprintf(ofname, "%s-theta.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[ORT_THETA].eb_raw();
	dt[ORT_THETA].y /= nconfs;
	dt[ORT_THETA].eb /= nconfs;
	dt[ORT_THETA].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(msg, "theta ==> %s", ofname);

	sprintf(ofname, "%s-omega.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[ORT_OMEGA].eb_raw();
	dt[ORT_OMEGA].y /= nconfs;
	dt[ORT_OMEGA].eb /= nconfs;
	dt[ORT_OMEGA].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(buffer, "\nomega ==> %s", ofname); strcat(msg, buffer);

	show_msg(msg);

	mlog.close();
	return true;
}

namespace _cluster_distribution_ {
using namespace _distribution_;
/*****************************************************************
**** RELATIVE ORIENTATION BETWEEN TWO CLUSTER IN SAME MONOMER ****
******************************************************************/

void cluster_relative_orientation_distribution(MMOLECULE &mm, int monomer_uid, VECT_IN_MONOMER* mv, DISTRIBT<float, long> *dt) {
	int na = 0, n_monomer = 0, npt = 0;
	_GROUP<CLUSTER> *monomer = NULL;
	MATOM *patom1 = NULL, *patom2 = NULL;

	VECTOR3 v_monomer[2], vr;
	float theta = 0, omega = 0;

	for (n_monomer = 0; n_monomer < mm.nMonomer; n_monomer++) {
		monomer = mm.monomer + n_monomer;
		if (monomer->uid != monomer_uid) continue;
		if (monomer->nUnit <= mv->atom[0].nc || monomer->nUnit <= mv->atom[1].nc) continue; // this is not supposed happen, but to make program safe, we check it here

		patom1 = monomer->u[mv[0].atom[0].nc]->atom + mv[0].atom[0].na;
		patom2 = monomer->u[mv[0].atom[1].nc]->atom + mv[0].atom[1].na;
		VECT3(patom1->r, patom2->r, v_monomer[0]) // patom1 -> patom2

		patom1 = monomer->u[mv[1].atom[0].nc]->atom + mv[1].atom[0].na;
		patom2 = monomer->u[mv[1].atom[1].nc]->atom + mv[1].atom[1].na;
		VECT3(patom1->r, patom2->r, v_monomer[1]) // patom1 -> patom2

		//VECT3(v_monomer[0], v_monomer[1], vr) // 0 ==> 1
		//theta = vr.theta(); omega = vr.phi();
		theta = angle(v_monomer[0], v_monomer[1]);
		
		//DISTRIBT_INDX(dt[ORT2_THETA], theta, npt)
		//if (LInRange(npt, 0, dt[ORT2_THETA].N)) dt[ORT2_THETA].y[npt]++;
		dt[ORT2_THETA].sample(theta, true);
		//DISTRIBT_INDX(dt[ORT2_OMEGA], omega, npt)
		//if (LInRange(npt, 0, dt[ORT2_OMEGA].N)) dt[ORT2_OMEGA].y[npt]++;
	}
}

#if _SYS_ == _WINDOWS_SYS_
int cluster_relative_orientation_distribution_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* cluster_relative_orientation_distribution_thread(void *void_par) { // thread-function
#endif
	DISTRIBT_RELORT_THREAD_VARS *par = (DISTRIBT_RELORT_THREAD_VARS*)void_par;
	par->func(*(par->mm), par->monomer_uid, par->mv, par->dt);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

bool cluster_relative_orientation_distribution_stat(char *mdpar_fname, int monomer_uid, VECT_IN_MONOMER *mv, int nf, int nt, int nstep, DISTRIBT<float, float> *dt, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", mol_template[256] = "\0", buffer[256] = "\0";
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

	int nThreads = 0, iconf = 0, nconf = 0, iThread, i, j, k;

	MMOLECULE mm[MAX_THREADS];
	DISTRIBT<float, long> ndt[MAX_THREADS][N_ORT2_DT];
	DISTRIBT_RELORT_THREAD_VARS thread_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		if (!ConstructSingleChainPolymerFromParamFile(mm[i], mol_template)) return false;
		for (j = 0; j < N_ORT2_DT; j++) ndt[i][j].set_distribt_x(dt[j].xf, dt[j].xt, dt[j].x.N - 1);
		thread_var[i].set_func(cluster_relative_orientation_distribution);
		thread_var[i].set_pars(mm + i, monomer_uid, mv, ndt[i], i);
	}

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
			ReadMoleculeStruct(*(rmdproc.in), mm[i]);
			sprintf(msg, " %d ", iconf); strcat(msg1, msg);

			iThread = i;			
		#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(thread_var[iThread].hMMThread);
			AfxBeginThread((AFX_THREADPROC)cluster_relative_orientation_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
		#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(thread_var[iThread].thread), NULL, &cluster_relative_orientation_distribution_thread, (void *)(thread_var + iThread));
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
	
	for (i = 0; i < N_ORT2_DT; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(dt[i], ndt[k][i]);
	}
	nconfs = nconf;
	return true;
}

bool MultiMol_cluster_relative_orientation_distribution_stat(char *mdpar_fname, ARRAY<short> &mol_id, int monomer_uid, VECT_IN_MONOMER *mv, int nf, int nt, int nstep, DISTRIBT<float, float> *dt, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", buffer[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
	int MD_MODE = 0;
	if (sscanf(msg, "%d", &MD_MODE) != 1 || (MD_MODE != MULTI_MM_FREE_CELL_MD && MD_MODE != MULTI_MM_PERIDIC_CELL_MD)) { // MultiMol ?
		sprintf(msg, "md is not working on multi-molecule");
		show_msg(msg); return false;
	}

	set_md_mode(MD_MODE);
	set_job(GET_MULTI_MM_MD_PROC);

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int nThreads = 0, iThread, imol = 0, i, j, k;
	long iconf = 0, nconf = 0;

	DISTRIBT<float, long> ndt[MAX_THREADS][N_ORT2_DT];
	DISTRIBT_RELORT_THREAD_VARS thread_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		for (j = 0; j < N_ORT2_DT; j++) ndt[i][j].set_distribt_x(dt[j].xf, dt[j].xt, dt[j].x.N - 1);
		thread_var[i].set_func(cluster_relative_orientation_distribution);
		thread_var[i].set_pars(NULL, monomer_uid, mv, ndt[i], i);
	}

	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		nThreads = 0;
		if (!serial_multi_mol_md_proc(iconf)) break;
		sprintf(msg, " %d in %s are read", iconf, mdpar_fname); show_log(msg, true);

		iThread = 0; nThreads = 0; show_log("analyzing molecules ", false);
		for (imol = 0; imol < ::mm_md_cell.mm.n; imol++) { // working on macro-molecules only
			if (mol_id.n > 0 && !mol_id.exist(imol)) {
				if (imol != ::mm_md_cell.mm.n - 1) continue; // this molecule is ignored
			}
			else {
				thread_var[iThread].mm = ::mm_md_cell.mm.m + imol;
			#if _SYS_ == _WINDOWS_SYS_ 
				ResetEvent(thread_var[iThread].hMMThread);
				AfxBeginThread((AFX_THREADPROC)cluster_relative_orientation_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
			#elif _SYS_ == _LINUX_SYS_
				pthread_create(&(thread_var[iThread].thread), NULL, &cluster_relative_orientation_distribution_thread, (void *)(thread_var + iThread));
			#endif

				sprintf(msg, "%d ", imol); show_log(msg, false);
				iThread++; nThreads++;
			}

			if (iThread == MAX_THREADS || imol == ::mm_md_cell.mm.n - 1) { // reaching maximum threads or all molecules are run
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

	for (i = 0; i < N_ORT2_DT; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(dt[i], ndt[k][i]);
	}
	nconfs = nconf;
	return true;
}

} // end of namespace _dipole_distribution

extern "C" bool relative_orientation_distribution(char *parfname) {
	using namespace _cluster_distribution_;
	::COMMAND_MD_STOP = false;
	mlog.init("rel-ort.log", false);
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	int npts = 0;
	int monomer_uid = 0;
	VECT_IN_MONOMER mv[2];
	char monomer[256] = "\0";

	DISTRIBT<float, float> dt[N_ORT2_DT];

	if (!search(in, "[MONOMER]", buffer)) {sprintf(msg, "can not find [MONOMER] in %s", parfname); show_msg(msg); in.close(); return false;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%s", monomer) != 1) {
		sprintf(msg, "monomer is not given"); show_msg(msg); in.close(); return false;
	}
	if (!get_monomer_uid(monomer, monomer_uid)) {
		sprintf(msg, "unknown monomer %s", monomer); show_msg(msg); in.close(); return false;
	}

	if (!search(in, "[VECTOR]", buffer)) {
		sprintf(msg, "vectors are not given"); show_msg(msg); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d %d %d", &(mv[0].atom[0].nc), &(mv[0].atom[0].na), &(mv[0].atom[1].nc), &(mv[0].atom[1].na)) != 4) {
		sprintf(msg, "vector definition : [cluster1] [atom1] [cluster2] [atom2], see %s", buffer);
		show_msg(msg); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d %d %d", &(mv[1].atom[0].nc), &(mv[1].atom[0].na), &(mv[1].atom[1].nc), &(mv[1].atom[1].na)) != 4) {
		sprintf(msg, "vector definition : [cluster1] [atom1] [cluster2] [atom2], see %s", buffer);
		show_msg(msg); in.close(); return false;
	}
	
	if (!search(in, "[DISTRIBUTION-RANGE]", buffer)) {sprintf(msg, "can not find [DISTRIBUTION-RANGE] in %s", parfname); show_msg(msg); in.close(); return false;}
	// theta
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &npts) != 1 || npts <= 0) {
		sprintf(msg, "RANGE format : [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[ORT2_THETA].set_distribt_x(0, fPI, npts);
	// omega
	/*
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &npts) != 1 || npts <= 0) {
		sprintf(msg, "RANGE format : [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[ORT2_OMEGA].set_distribt_x(0, fPI2, npts);
	*/
	in.close();

	ARRAY<short> mol_id;
	read_given_id(parfname, "[MOLECULES]", mol_id);

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

		if (MD_MODE == SINGLE_MM_FREE_CELL_MD) status = cluster_relative_orientation_distribution_stat(mdpar_fname, monomer_uid, mv, nf, nt, nstep, dt, nconf);
		else if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) 
			status = MultiMol_cluster_relative_orientation_distribution_stat(mdpar_fname, mol_id, monomer_uid, mv, nf, nt, nstep, dt, nconf);
		else {
			sprintf(msg, "UNKNOWN MD in %s", mdpar_fname); show_msg(msg, true); nconf = 0;
		}

		nconfs += nconf;
		n++;
	}
	in_cfg.close();

	ofstream *out = NULL;

	sprintf(ofname, "%s-theta.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[ORT2_THETA].eb_raw();
	dt[ORT2_THETA].y /= nconfs;
	dt[ORT2_THETA].eb /= nconfs;
	dt[ORT2_THETA].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(msg, "theta ==> %s", ofname);

	/*
	sprintf(ofname, "%s-omega.dat", title);
	out = new ofstream;
	out->open(ofname);
	for (n = 0; n < dt[ORT2_OMEGA].N; n++) {
		dt[ORT2_OMEGA].y[n] /= nconfs;
		(*out)<<dt[ORT2_OMEGA].x[n] / fPI * 180 <<"  "<<dt[ORT2_OMEGA].y[n]<<endl;
	}
	out->close(); delete out; out = NULL;
	sprintf(buffer, "\nomega ==> %s", ofname); strcat(msg, buffer);
	*/

	show_msg(msg);

	mlog.close();
	return true;
}


namespace _cluster_distribution_ {

/*********************************************************************
*********  DIPOLE DIRECTION RELATIVE TO MONOMER ORIENTATION  *********
**********************************************************************/

void dipole_local_orientation_distribution(MMOLECULE &mm, int monomer_uid, VECT_IN_MONOMER* mv, bool limit, float r1, float r2, DISTRIBT<float, long> *dt) {
	int nc = 0, na = 0, n_monomer = 0, npt = 0, i;
	_GROUP<CLUSTER> *monomer = NULL;
	CLUSTER *pc = NULL;
	MATOM *patom1 = NULL, *patom2 = NULL;

	VECTOR3 v_monomer, mu;
	float theta = 0, omega = 0, mu_abs = 0;

	for (n_monomer = 0; n_monomer < mm.nMonomer; n_monomer++) {
		monomer = mm.monomer + n_monomer;
		if (monomer->uid != monomer_uid) continue;
		if (monomer->nUnit <= mv->atom[0].nc || monomer->nUnit <= mv->atom[1].nc) continue; // this is not supposed happen, but to make program safe, we check it here

		patom1 = monomer->u[mv->atom[0].nc]->atom + mv->atom[0].na;
		patom2 = monomer->u[mv->atom[1].nc]->atom + mv->atom[1].na;
		VECT3(patom1->r, patom2->r, v_monomer) // patom1 -> patom2

		V3zero(mu)
		for (nc = 0; nc < monomer->nUnit; nc++) {
			pc = monomer->u[nc];
			for (na = 0; na < pc->nAtoms; na++) {
				for (i = 0; i < 3; i++) mu.v[i] += pc->atom[na].c0 * pc->atom[na].r.v[i];
			}
		}
		mu_abs = mu.abs();
		if (limit && (mu_abs > r2 || mu_abs < r1)) continue;

		if (mu_abs < 1e-6) continue; // |mu| is too small, 0 ?
		theta = angle(mu, v_monomer);
		
		//DISTRIBT_INDX(dt[LOCAL_MU_THETA], theta, npt) 
		//if (LInRange(npt, 0, dt[LOCAL_MU_THETA].N)) dt[LOCAL_MU_THETA].y[npt]++;
		dt[LOCAL_MU_THETA].sample(theta, true);
		//DISTRIBT_INDX(dt[ORT2_OMEGA], omega, npt) 
		//if (LInRange(npt, 0, dt[ORT2_OMEGA].N)) dt[ORT2_OMEGA].y[npt]++;
	}
}

#if _SYS_ == _WINDOWS_SYS_
int dipole_local_orientation_distribution_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* dipole_local_orientation_distribution_thread(void *void_par) { // thread-function
#endif
	DISTRIBT_LOCAL_DIPOLE_THREAD_VARS *par = (DISTRIBT_LOCAL_DIPOLE_THREAD_VARS*)void_par;
	par->func(*(par->mm), par->monomer_uid, &(par->mv), par->limit, par->r1, par->r2, par->dt);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

bool dipole_local_orientation_distribution_stat(char *mdpar_fname, int monomer_uid, VECT_IN_MONOMER *mv, bool limit, float r1, float r2, int nf, int nt, int nstep, DISTRIBT<float, float> *dt, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", mol_template[256] = "\0", buffer[256] = "\0";
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

	int nThreads = 0, iconf = 0, nconf = 0, iThread, i, j, k;

	MMOLECULE mm[MAX_THREADS];
	DISTRIBT<float, long> ndt[MAX_THREADS][N_LOCAL_MU_DT];
	DISTRIBT_LOCAL_DIPOLE_THREAD_VARS thread_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		if (!ConstructSingleChainPolymerFromParamFile(mm[i], mol_template)) return false;
		for (j = 0; j < N_LOCAL_MU_DT; j++) ndt[i][j].set_distribt_x(dt[j].xf, dt[j].xt, dt[j].x.N - 1);
		thread_var[i].set_func(dipole_local_orientation_distribution);
		thread_var[i].set_pars(mm + i, monomer_uid, *mv, ndt[i], i);
		thread_var[i].set_limit(limit, r1, r2);
	}

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
			ReadMoleculeStruct(*(rmdproc.in), mm[i]);
			sprintf(msg, " %d ", iconf); strcat(msg1, msg);

			iThread = i;
		#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(thread_var[iThread].hMMThread);
			AfxBeginThread((AFX_THREADPROC)dipole_local_orientation_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
		#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(thread_var[iThread].thread), NULL, &dipole_local_orientation_distribution_thread, (void *)(thread_var + iThread));
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
	
	for (i = 0; i < N_LOCAL_MU_DT; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(dt[i], ndt[k][i]);
	}
	nconfs = nconf;
	return true;
}

bool MultiMol_dipole_local_orientation_distribution_stat(char *mdpar_fname, ARRAY<short> &mol_id, int monomer_uid, VECT_IN_MONOMER *mv, bool limit, float r1, float r2, int nf, int nt, int nstep, DISTRIBT<float, float> *dt, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", buffer[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
	int MD_MODE = 0;
	if (sscanf(msg, "%d", &MD_MODE) != 1 || (MD_MODE != MULTI_MM_FREE_CELL_MD && MD_MODE != MULTI_MM_PERIDIC_CELL_MD)) { // MultiMol ?
		sprintf(msg, "md is not working on multi-molecule");
		show_msg(msg); return false;
	}

	set_md_mode(MD_MODE);
	set_job(GET_MULTI_MM_MD_PROC);

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int nThreads = 0, iThread, imol = 0, i, j, k;
	long iconf = 0, nconf = 0;

	DISTRIBT<float, long> ndt[MAX_THREADS][N_LOCAL_MU_DT];
	DISTRIBT_LOCAL_DIPOLE_THREAD_VARS thread_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		for (j = 0; j < N_LOCAL_MU_DT; j++) ndt[i][j].set_distribt_x(dt[j].xf, dt[j].xt, dt[j].x.N - 1);
		thread_var[i].set_func(dipole_local_orientation_distribution);
		thread_var[i].set_pars(NULL, monomer_uid, *mv, ndt[i], i);
		thread_var[i].set_limit(limit, r1, r2);
	}

	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		nThreads = 0;
		if (!serial_multi_mol_md_proc(iconf)) break;
		sprintf(msg, " %d in %s are read", iconf, mdpar_fname); show_log(msg, true);

		iThread = 0; nThreads = 0; show_log("analyzing molecules ", false);
		for (imol = 0; imol < ::mm_md_cell.mm.n; imol++) { // working on macro-molecules only
			if (mol_id.n > 0 && !mol_id.exist(imol)) {
				if (imol != ::mm_md_cell.mm.n - 1) continue; // this molecule is ignored
			}
			else {
				thread_var[iThread].mm = ::mm_md_cell.mm.m + imol;
			#if _SYS_ == _WINDOWS_SYS_ 
				ResetEvent(thread_var[iThread].hMMThread);
				AfxBeginThread((AFX_THREADPROC)dipole_local_orientation_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
			#elif _SYS_ == _LINUX_SYS_
				pthread_create(&(thread_var[iThread].thread), NULL, &dipole_local_orientation_distribution_thread, (void *)(thread_var + iThread));
			#endif

				sprintf(msg, "%d ", imol); show_log(msg, false);
				iThread++; nThreads++;
			}

			if (iThread == MAX_THREADS || imol == ::mm_md_cell.mm.n - 1) { // reaching maximum threads or all molecules are run
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

	for (i = 0; i < N_LOCAL_MU_DT; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(dt[i], ndt[k][i]);
	}
	nconfs = nconf;
	return true;
}
} // end of namespace _cluster_distribution_

extern "C" bool dipole_local_orientation(char *parfname) {
	using namespace _cluster_distribution_;
	::COMMAND_MD_STOP = false;
	mlog.init("dipole-local-ort.log", false);
	char buffer[256] = "\0", msg[2560] = "\0", title[250] = "\0", ofname[256] = "\0";

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	int npts = 0;
	int monomer_uid = 0;
	VECT_IN_MONOMER mv;
	char monomer[256] = "\0";

	DISTRIBT<float, float> dt[N_LOCAL_MU_DT];

	if (!search(in, "[MONOMER]", buffer)) {sprintf(msg, "can not find [MONOMER] in %s", parfname); show_msg(msg); in.close(); return false;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%s", monomer) != 1) {
		sprintf(msg, "monomer is not given"); show_msg(msg); in.close(); return false;
	}
	if (!get_monomer_uid(monomer, monomer_uid)) {
		sprintf(msg, "unknown monomer %s", monomer); show_msg(msg); in.close(); return false;
	}

	if (!search(in, "[VECTOR]", buffer)) {
		sprintf(msg, "vectors are not given"); show_msg(msg); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d %d %d", &(mv.atom[0].nc), &(mv.atom[0].na), &(mv.atom[1].nc), &(mv.atom[1].na)) != 4) {
		sprintf(msg, "vector definition : [cluster1] [atom1] [cluster2] [atom2], see %s", buffer);
		show_msg(msg); in.close(); return false;
	}

	bool limit = false;
	float r1 = 0, r2 = 0;
	int nlimit = 0;
	if (!search(in, "[LIMIT]", buffer))  {sprintf(msg, "can not find [LIMIT] in %s", parfname); show_msg(msg); in.close(); return false;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nlimit) != 1) {
		sprintf(msg, "can not get limit [0/1] [from] [to] in %s", parfname); 
		show_msg(msg); in.close(); return false;
	}
	limit = (nlimit > 0 ? true : false);
	if (limit && sscanf(buffer, "%d %f %f", &nlimit, &r1, &r2) != 3) {
		sprintf(msg, "can not get limit [0/1] [from] [to] in %s", parfname); 
		show_msg(msg); in.close(); return false;
	}
	
	if (!search(in, "[DISTRIBUTION-RANGE]", buffer)) {sprintf(msg, "can not find [DISTRIBUTION-RANGE] in %s", parfname); show_msg(msg); in.close(); return false;}
	// theta
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &npts) != 1 || npts <= 0) {
		sprintf(msg, "RANGE format : [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[LOCAL_MU_THETA].set_distribt_x(0, fPI, npts);
	// omega
	/*
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &npts) != 1 || npts <= 0) {
		sprintf(msg, "RANGE format : [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	dt[ORT2_OMEGA].set_distribt_x(0, fPI2, npts);
	*/
	in.close();

	ARRAY<short> mol_id;
	read_given_id(parfname, "[MOLECULES]", mol_id);

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

		if (MD_MODE == SINGLE_MM_FREE_CELL_MD) status = dipole_local_orientation_distribution_stat(mdpar_fname, monomer_uid, &mv, limit, r1, r2, nf, nt, nstep, dt, nconf);
		else if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) 
			status = MultiMol_dipole_local_orientation_distribution_stat(mdpar_fname, mol_id, monomer_uid, &mv, limit, r1, r2, nf, nt, nstep, dt, nconf);
		else {
			sprintf(msg, "UNKNOWN MD in %s", mdpar_fname); show_msg(msg, true); nconf = 0;
		}

		nconfs += nconf;
		n++;
	}
	in_cfg.close();

	ofstream *out = NULL;

	sprintf(ofname, "%s-theta.dat", title);
	out = new ofstream;
	out->open(ofname);
	dt[LOCAL_MU_THETA].eb_raw();
	dt[LOCAL_MU_THETA].y /= nconfs;
	dt[LOCAL_MU_THETA].eb /= nconfs;
	dt[LOCAL_MU_THETA].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;
	sprintf(msg, "theta ==> %s", ofname);

	/*
	sprintf(ofname, "%s-omega.dat", title);
	out = new ofstream;
	out->open(ofname);
	for (n = 0; n < dt[ORT_OMEGA].N; n++) {
		dt[ORT_OMEGA].y[n] /= nconfs;
		(*out)<<dt[ORT_OMEGA].x[n] / fPI * 180 <<"  "<<dt[ORT_OMEGA].y[n]<<endl;
	}
	out->close(); delete out; out = NULL;
	sprintf(buffer, "\nomega ==> %s", ofname); strcat(msg, buffer);
	*/

	show_msg(msg);

	mlog.close();
	return true;
}

namespace _cluster_distribution_ {

/*********************************************************************
***************       TORSION ANGLE DISTRIBUTION       ***************
**********************************************************************/

void cluster_torsion_angle_distribution(MMOLECULE &mm, TA_DISTRIBT<float, long> *distb, int ndistb) {
	int nc = 0, na = 0, n_monomer = 0, npt = 0, idistb;
	_GROUP<CLUSTER> *monomer = NULL;
	CLUSTER *pc = NULL;
	double ta = 0;

	for (n_monomer = 0; n_monomer < mm.nMonomer; n_monomer++) {
		monomer = mm.monomer + n_monomer;
		for (idistb = 0; idistb < ndistb; idistb++) {
			if (monomer->uid != distb[idistb].monomer_uid) continue;
			if (distb[idistb].cindx >= monomer->nUnit) continue; // does not exist in monomer
			pc = monomer->u[distb[idistb].cindx];
			//if (pc->cID != cluster_uid) continue;
			ta = pc->km.ta;
			PERIOD_RANGE(ta, -PI, PI, PI2)
			//DISTRIBT_INDX(distb[idistb].dt, ta, npt) 
			//if (LInRange(npt, 0, distb[idistb].dt.N)) distb[idistb].dt.y[npt]++;
			distb[idistb].dt.sample(ta, true);
		}
	}
}

#if _SYS_ == _WINDOWS_SYS_
int cluster_torsion_angle_distribution_thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* cluster_torsion_angle_distribution_thread(void *void_par) { // thread-function
#endif
	DISTRIBT_CLUSTER_TA_THREAD_VARS *par = (DISTRIBT_CLUSTER_TA_THREAD_VARS*)void_par;
	par->func(*(par->mm), par->distb, par->ndistb);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

bool cluster_torsion_angle_distribution_stat(char *mdpar_fname, int nf, int nt, int nstep, TA_DISTRIBT<float, float> *distb, int ndistb, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", mol_template[256] = "\0", buffer[256] = "\0";
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

	int nThreads = 0, iconf = 0, nconf = 0, iThread, i, j, k;

	MMOLECULE mm[MAX_THREADS];
	FLEX_ASYM_MATRIX< TA_DISTRIBT<float, long> > mdt;
	mdt.set(MAX_THREADS, ndistb);
	DISTRIBT_CLUSTER_TA_THREAD_VARS thread_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		if (!ConstructSingleChainPolymerFromParamFile(mm[i], mol_template)) return false;
		for (j = 0; j < ndistb; j++) {
			mdt.m[i].m[j].monomer_uid = distb[j].monomer_uid;
			mdt.m[i].m[j].cindx = distb[j].cindx;
			mdt.m[i].m[j].dt.set_distribt_x(distb[j].dt.xf, distb[j].dt.xt, distb[j].dt.x.N - 1);
		}
		thread_var[i].set_func(cluster_torsion_angle_distribution);
		thread_var[i].set_pars(mm + i, mdt.m[i].m, ndistb, i);
	}

	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		nThreads = 0;
		strcpy(msg1, "\0");
		for (iThread = 0; iThread < MAX_THREADS; iThread++) {
			if (nt > 0 && iconf > nt) break;
			if (!rmdproc.search_proc(iconf, buffer)) break;
			//ReadMoleculeStruct(*(rmdproc.in), mm[iThread]);
			ReadMoleculeStructInfo(*(rmdproc.in), mm[iThread]);
			sprintf(msg, " %d ", iconf); strcat(msg1, msg);

		#if _SYS_ == _WINDOWS_SYS_ 
			ResetEvent(thread_var[iThread].hMMThread);
			thread_var[iThread].thread = AfxBeginThread((AFX_THREADPROC)cluster_torsion_angle_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
		#elif _SYS_ == _LINUX_SYS_
			pthread_create(&(thread_var[iThread].thread), NULL, &cluster_torsion_angle_distribution_thread, (void *)(thread_var + iThread));
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
	
	for (i = 0; i < ndistb; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(distb[i].dt, mdt.m[k].m[i].dt);
	}
	nconfs = nconf;
	return true;
}

bool MultiMol_cluster_torsion_angle_distribution_stat(char *mdpar_fname, ARRAY<short> &mol_id, int nf, int nt, int nstep, TA_DISTRIBT<float, float> *distb, int ndistb, int &nconfs) {
	nconfs = 0;
	char msg[256] = "\0", msg1[256] = "\0", buffer[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
	int MD_MODE = 0;
	if (sscanf(msg, "%d", &MD_MODE) != 1 || (MD_MODE != MULTI_MM_FREE_CELL_MD && MD_MODE != MULTI_MM_PERIDIC_CELL_MD)) { // MultiMol ?
		sprintf(msg, "md is not working on multi-molecule");
		show_msg(msg); return false;
	}

	set_md_mode(MD_MODE);
	set_job(GET_MULTI_MM_MD_PROC);

	bool status = construct_macromolecule(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
	status = init_read_mdproc(mdpar_fname);
	if (!status) {sprintf(msg, "failure to initilize reading md procedure"); show_msg(msg); return false;}

	int nThreads = 0, iThread, imol = 0, i, j, k;
	long iconf = 0, nconf = 0;

	FLEX_ASYM_MATRIX< TA_DISTRIBT<float, long> > mdt;
	mdt.set(MAX_THREADS, ndistb);
	DISTRIBT_CLUSTER_TA_THREAD_VARS thread_var[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		for (j = 0; j < ndistb; j++) {
			mdt.m[i].m[j].monomer_uid = distb[j].monomer_uid;
			mdt.m[i].m[j].cindx = distb[j].cindx;
			mdt.m[i].m[j].dt.set_distribt_x(distb[j].dt.xf, distb[j].dt.xt, distb[j].dt.x.N - 1);
		}
		thread_var[i].set_func(cluster_torsion_angle_distribution);
		thread_var[i].set_pars(NULL, mdt.m[i].m, ndistb, i);
	}

	iconf = nf;
	nconf = 0;
	while (1) {
		if (::COMMAND_MD_STOP) break;
		if (nt > 0 && iconf > nt) break;
		nThreads = 0;
		if (!serial_multi_mol_md_proc(iconf)) break;
		sprintf(msg, " %d in %s are read", iconf, mdpar_fname); show_log(msg, true);

		iThread = 0; nThreads = 0; show_log("analyzing molecules ", false);
		for (imol = 0; imol < ::mm_md_cell.mm.n; imol++) { // working on macro-molecules only
			if (mol_id.n > 0 && !mol_id.exist(imol)) {
				if (imol != ::mm_md_cell.mm.n - 1) continue; // this molecule is ignored
			}
			else {
				thread_var[iThread].mm = ::mm_md_cell.mm.m + imol;
				if (!search(*(::rmdproc.in), "MOL", imol, buffer)) break;
				//ReadMoleculeStruct(*(rmdproc.in), mm[iThread]);
				ReadMoleculeStructInfo(*(::rmdproc.in), ::mm_md_cell.mm.m[imol]);

			#if _SYS_ == _WINDOWS_SYS_ 
				ResetEvent(thread_var[iThread].hMMThread);
				thread_var[iThread].thread = AfxBeginThread((AFX_THREADPROC)cluster_torsion_angle_distribution_thread, (LPVOID)(thread_var + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
			#elif _SYS_ == _LINUX_SYS_
				pthread_create(&(thread_var[iThread].thread), NULL, &cluster_torsion_angle_distribution_thread, (void *)(thread_var + iThread));
			#endif

				sprintf(msg, "%d ", imol); show_log(msg, false);
				iThread++; nThreads++;
			}

			if (iThread == MAX_THREADS || imol == ::mm_md_cell.mm.n - 1) { // reaching maximum threads or all molecules are run
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
	
	for (i = 0; i < ndistb; i++) {
		for (k = 0; k < MAX_THREADS; k++) raw_distribt_combine<float, float, float, long>(distb[i].dt, mdt.m[k].m[i].dt);
	}
	nconfs = nconf;
	return true;
}
} // end of namespace _cluster_distribution_

extern "C" bool cluster_torsion_angle(char *parfname) {
	using namespace _cluster_distribution_;
	::COMMAND_MD_STOP = false;
	mlog.init("ta.log", false);
	char buffer[256] = "\0", msg[2560] = "\0", ofname[256] = "\0";

	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {sprintf(msg, "can not open file %s", parfname); show_log(msg, true); return false;}

	int npts = 90;
	int ncluster = 0, i;
	char monomer[256] = "\0";

	ARRAY< TA_DISTRIBT<float, float> > dt;

	if (!search(in, "[CLUSTER]", buffer)) {sprintf(msg, "can not find [CLUSTER] in %s", parfname); show_msg(msg); in.close(); return false;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &ncluster) != 1) {
		sprintf(msg, "number of clusters is not given"); show_msg(msg); in.close(); return false;
	}
	if (ncluster <= 0) {
		sprintf(msg, "number of clusters is %d", ncluster); show_msg(msg); in.close(); return true;
	}
	dt.SetArray(ncluster);
	for (i = 0; i < ncluster; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s %d %s", monomer, &(dt.m[i].cindx), dt.m[i].title) != 3) {
			sprintf(msg, "cluster is not given in format [monomer, cluster id in monomer, output-filename's title] : %s", buffer); 
			show_msg(msg); in.close(); return false;
		}
		if (!get_monomer_uid(monomer, dt.m[i].monomer_uid)) {
			sprintf(msg, "unknown monomer %s", monomer); show_msg(msg); in.close(); return false;
		}
	}
	
	if (!search(in, "[DISTRIBUTION-RANGE]", buffer)) {sprintf(msg, "can not find [DISTRIBUTION-RANGE] in %s", parfname); show_msg(msg); in.close(); return false;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &npts) != 1 || npts <= 0) {
		sprintf(msg, "RANGE format : [npts], see %s", buffer); show_msg(msg); in.close(); return false;
	}
	for (i = 0; i < dt.n; i++) {
		dt.m[i].dt.set_distribt_x(-fPI, fPI, npts); // torsion angle in NEIMO is [-PI, PI]
	}

	in.close();

	ARRAY<short> mol_id;
	read_given_id(parfname, "[MOLECULES]", mol_id);

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

		if (MD_MODE == SINGLE_MM_FREE_CELL_MD) status = cluster_torsion_angle_distribution_stat(mdpar_fname, nf, nt, nstep, dt.m, dt.n, nconf);
		else if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) 
			status = MultiMol_cluster_torsion_angle_distribution_stat(mdpar_fname, mol_id, nf, nt, nstep, dt.m, dt.n, nconf);
		else {
			sprintf(msg, "UNKNOWN MD in %s", mdpar_fname); show_msg(msg, true); nconf = 0;
		}

		nconfs += nconf;
		n++;
	}
	in_cfg.close();

	ofstream *out = NULL;

	for (i = 0; i < dt.n; i++) {
		sprintf(ofname, "%s-ta.dat", dt.m[i].title);
		out = new ofstream;
		out->open(ofname);
		if (!out->is_open()) {
			sprintf(msg, "can not open file %s", ofname); show_msg(msg); out->close();
			delete out; out = NULL; continue;
		}
		else {
			dt.m[i].dt.eb_raw();
			dt.m[i].dt.y /= nconfs;
			dt.m[i].dt.eb /= nconfs;
			dt.m[i].dt.save_x_y_eb(*out);
			out->close(); delete out; out = NULL;
		}
		sprintf(msg, "torsion angle distribution ==> %s", ofname);
		show_log(msg, true);
	}
	show_msg("torsion angle distribution is calculated");

	mlog.close();
	return true;
}

bool read_given_monomer(char *fname, char* title, ARRAY<short>& id_monomer) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return false;

	char buffer[256] = "\0", *bf = NULL, monomer[250] = "\0";
	char seperator[5] = ",; \t";
	int uid = 0, n = 0;
	short *p_uid = NULL;
	CHAIN<short> *mch = NULL, *tail = NULL;
	if (!search(in, title, buffer)) {in.close(); return false;}
	while (!in.eof()) {
		in.getline(buffer, 250); bf = buffer;
		n = 0;
		while (strlen(bf) > 0) {
			if (sscanf(bf, "%s", monomer) != 1) break;
			if (!get_monomer_uid(monomer, uid)) break;
			p_uid = new short; *p_uid = (short)uid;
			if (mch == NULL) {
				mch = new CHAIN<short>; tail = mch; mch->p = p_uid;
			}
			else {
				tail->next = new CHAIN<short>; tail = tail->next; tail->p = p_uid;
			}
			n++;
			if (!NextSection(&bf, seperator)) break;
		}
		if (n == 0) break;
	}
	in.close();
	n = number<short>(mch);
	if (n == 0) return false;
	id_monomer.SetArray(n);
	tail = mch;
	for (n = 0; n < id_monomer.n; n++) {
		if (tail == NULL) {id_monomer.SetArray(0); return false;}
		id_monomer.m[n] = *(tail->p);
		tail = tail->next;
	}
	release_chain<short>(&mch, true); tail = NULL;
	return true;
}

bool read_given_cgatoms(char *fname, char* title, ARRAY<short>& id_cgatom) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return false;

	char buffer[256] = "\0", *bf = NULL, cgatom[250] = "\0";
	char seperator[5] = ",; \t";
	int uid = 0, n = 0;
	short *p_uid = NULL;
	CHAIN<short> *mch = NULL, *tail = NULL;
	if (!search(in, title, buffer)) {in.close(); return false;}
	while (!in.eof()) {
		in.getline(buffer, 250); bf = buffer;
		n = 0;
		while (strlen(bf) > 0) {
			if (sscanf(bf, "%s", cgatom) != 1) break;
			if (!get_cgatom_uid(cgatom, uid)) break;
			p_uid = new short; *p_uid = (short)uid;
			if (mch == NULL) {
				mch = new CHAIN<short>; tail = mch; mch->p = p_uid;
			}
			else {
				tail->next = new CHAIN<short>; tail = tail->next; tail->p = p_uid;
			}
			n++;
			if (!NextSection(&bf, seperator)) break;
		}
		if (n == 0) break;
	}
	in.close();
	n = number<short>(mch);
	if (n == 0) return false;
	id_cgatom.SetArray(n);
	tail = mch;
	for (n = 0; n < id_cgatom.n; n++) {
		if (tail == NULL) {id_cgatom.SetArray(0); return false;}
		id_cgatom.m[n] = *(tail->p);
		tail = tail->next;
	}
	release_chain<short>(&mch, true); tail = NULL;
	return true;
}

bool SaveConfigXYZ(char *spar_fname) {
	char msg[256] = "\0", buffer[256] = "\0", ofname[256] = "\0";
	if (!read_item(spar_fname, "[OUTPUT_FILE]", buffer) || sscanf(buffer, "%s", ofname) != 1) {
		sprintf(msg, "failure to get output filename from %s", spar_fname); show_log(msg, true); return false;
	}

	ofstream out;
	out.open(ofname);
	if (!out.is_open()) {
		sprintf(msg, "can not open file %s", ofname); show_log(msg, true); return false;
	}

	ARRAY<short> id_mm, id_monomer, id_cluster, id_cgatom;
	if (::MD_MODE == SINGLE_CGMM_FREE_CELL_MD || ::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD || ::MD_MODE == MULTI_CGMM_FREE_CELL_MD) {
		read_given_id(spar_fname, "[MOLECULES]", id_mm);
		read_given_cgatoms(spar_fname, "[CGATOMS]", id_cgatom);
	}
	else {
		read_given_id(spar_fname, "[MOLECULES]", id_mm);
		read_given_monomer(spar_fname, "[MONOMERS]", id_monomer);
		read_given_id(spar_fname, "[CLUSTER_UID]", id_cluster);
	}

	int dn[3] = {0, 0, 0}, i, j, k;
	double dv[3] = {0, 0, 0};
	if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD || ::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
		if (!read_item(spar_fname, "[CELL_EXTENSION]", buffer) || sscanf(buffer, "%d %d %d", dn, dn+1, dn+2) != 3) {
			dn[0] = 0; dn[1] = 0; dn[2] = 0; dv[0] = 0; dv[1] = 0; dv[2] = 0;
		}
	}
	
	GENERAL_MOLECULE gm;
	int imol = 0;
	for (i = -dn[0]; i <= dn[0]; i++) {for (j = -dn[1]; j <= dn[1]; j++) {for (k = -dn[2]; k <= dn[2]; k++) {
		if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
			dv[0] = i * ::mm_md_cell.h[0] * 2;
			dv[1] = j * ::mm_md_cell.h[1] * 2;
			dv[2] = k * ::mm_md_cell.h[2] * 2;
		}
		else if (::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			dv[0] = i * ::cgmm_md_cell.h[0] * 2;
			dv[1] = j * ::cgmm_md_cell.h[1] * 2;
			dv[2] = k * ::cgmm_md_cell.h[2] * 2;
		}
		imol = 0;
		while (1) {
			if (!get_molecule_upto_job(imol, gm)) break;
			if (gm.mol_type == 1) { // macro-molecule
				if (gm.mm == NULL) break;
				MMSaveXYZ(&out, gm.mm, id_monomer, id_cluster, true, dv); // save the cordinates in CMM, + dr_cmm
			}
			else if (gm.mol_type == 0) { // simple molecule
				if (gm.sm == NULL) break;
				SMSaveXYZ(&out, gm.sm, id_cluster);
			}
			else if (gm.mol_type == 4) { // point molecule
				if (gm.pm == NULL) break;
				PMSaveXYZ(&out, gm.pm, id_cluster);
			}
			else if (gm.mol_type == 2) { // coarse-grained macro-molecule
				if (gm.cgmm == NULL) break;
				CGMMSaveXYZ(&out, gm.cgmm, id_cgatom, true, dv); // save the cordinates in CMM, + dr_cmm
			}
			else if (gm.mol_type == 3) { // coarse-grained simple molecule
				if (gm.cgsm == NULL) break;
				CGSMSaveXYZ(&out, gm.cgsm, id_cgatom);
			}
			imol += 1;
		}
	}}}
	out.close();
	sprintf(msg, "xyz-struct is saved in file %s", ofname); show_log(msg, true);
	return true;
}

extern "C" bool MultiSaveConfigXYZ(char *parfname) {
	char buffer[256] = "\0", mdpar_fname[256] = "\0", spar_fname[256] = "\0";
	char msg[256] = "\0";
	int iConfig = 0;
	mlog.init("xyz.log", false);

	if (!read_item(parfname, "[MDPAR_FILE]", buffer) || sscanf(buffer, "%s", mdpar_fname) != 1) {
		sprintf(msg, "mdpar filename is not given after [MDPAR_FILE] in %s", parfname); show_log(msg, true); 
		::mlog.close(); return false;
	}

	int MD_MODE = 0;
	if (!read_item(mdpar_fname, "[MD_CELL]", buffer) || sscanf(buffer, "%d", &MD_MODE) != 1) {
		sprintf(msg, "failure to get the type of MD from %s", mdpar_fname); show_log(msg, true);
		::mlog.close(); return false;
	}
	set_md_mode(MD_MODE);
	if (MD_MODE == SINGLE_MM_FREE_CELL_MD) set_job(GET_SINGLE_MM_MD_PROC);
	else if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) 
		set_job(GET_MULTI_MM_MD_PROC);
	else if (MD_MODE == SINGLE_CGMM_FREE_CELL_MD) set_job(GET_SINGLE_CGMM_MD_PROC);
	else if (MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) 
		set_job(GET_MULTI_CGMM_MD_PROC);
	else {
		sprintf(msg, "the type of MD in %s is UNKNOWN", mdpar_fname); show_log(msg, true);
		::mlog.close(); return false;
	}
	
	ifstream in;
	in.open(parfname);
	if (!in.is_open()) return false;
	if (!search(in, "[XYZ_STRUCT]", buffer)) {
		sprintf(msg, "[XYZ_STRUCT] is not defined in %s", parfname); show_log(msg, true);
		::mlog.close(); in.close(); return false;
	}

	if (!construct_macromolecule(mdpar_fname)) {
		sprintf(msg, "failure to initilize molecule structure from %s", mdpar_fname); show_log(msg, true);
		::mlog.close(); in.close(); return false;
	}
	if (!init_read_mdproc(mdpar_fname)) {
		sprintf(msg, "failure to initilize serial-reading-md-procedure from %s", mdpar_fname); show_log(msg, true);
		::mlog.close(); in.close(); return false;
	}

	if (MD_MODE == SINGLE_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_FREE_CELL_MD ||
		MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			if (!cgatompar_db.get_pdb_alias(parfname, "[CGATOM_PDB_ALIAS]")) return false;
	}

	while (!in.eof()) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %s", &iConfig, spar_fname) != 2) break;
		if (!serial_md_proc(iConfig)) break;
		if (MD_MODE == MULTI_MM_PERIDIC_CELL_MD) molecule_check_periodic_cell(::mm_md_cell);
		else if (MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) _coarse_grain_::molecule_check_periodic_cell(::cgmm_md_cell);
		if (!SaveConfigXYZ(spar_fname)) break;
	}
	in.close(); ::mlog.close();
	return true;
}



// save the configuration, with extented periodical cell, keeping the interested cluster and its neighbors (within a radius)
bool read_interested_cluster(char *parfname, char *title, ARRAY<short> &cid, ARRAY<float> &iradius, ARRAY<float> &oradius, ARRAY<int> &mindx) {
	char msg[256] = "\0", buffer[256] = "\0", cname[256] = "\0", *bf = NULL;
	int n = 0, i, nt;
	ifstream in;
	in.open(parfname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open file %s", parfname); show_msg(msg, true); return false;
	}
	if (!search(in, title, buffer)) {
		sprintf(msg, "%s is not defined in file %s", title, parfname); show_msg(msg, true); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &n) != 1 || n <= 0) {
		sprintf(msg, "no interested cluster is defined after %s", title); show_msg(msg); in.close(); return false;
	}
	cid.SetArray(n); iradius.SetArray(n); oradius.SetArray(n); mindx.SetArray(n);
	for (i = 0; i < n; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d", &nt) != 1) {
			if (sscanf(buffer, "%s", cname) != 1) {
				sprintf(msg, "cluster is not given : %s", buffer); show_log(msg, true); in.close(); return false;
			}
			else if (!get_cluster_uid(cname, nt)) {
				sprintf(msg, "unknown cluster : %s", cname); show_log(msg, true); in.close(); return false;
			}
		}
		bf = buffer; NextSection(&bf, seperator);
		if (sscanf(bf, "%f %f", iradius.m + i, oradius.m + i) != 2) {
			sprintf(msg, "can not get definition of %dth interested cluster [with format: cluster [name/UID]  inner-radius  outter-radius], see : %s", i, buffer); show_msg(msg); in.close(); return false;
		}
		else cid.m[i] = nt;
		NextSection(&bf); NextSection(&bf);
		if (bf == NULL || sscanf(bf, "%d", &nt) == 1) mindx.m[i] = nt; // nth of this type of cluster
		else mindx.m[i] = -1; // not specified
	}
	in.close(); return true;
}

bool init_interested_cluster(MMOL_MD_CELL &mdcell, ARRAY<short> &cid, ARRAY<int> &mindx, ARRAY<BASIC_CLUSTER*> &icluster) {
	char msg[256] = "\0";
	int im, ic, type_indx;
	BASIC_CLUSTER *pc = NULL;
	int nc = 0;
	for (im = 0; im < mdcell.mm.n; im++) {
		for (ic = 0; ic < mdcell.mm.m[im].nCluster; ic++) {
			pc = (BASIC_CLUSTER*)(mdcell.mm.m[im].cluster + ic);
			if (cid.exist(pc->cID)) {
				type_indx = cid.indx(pc->cID);
				if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im) continue;
				nc++;
			}
		}
	}
	for (im = 0; im < mdcell.sm.n; im++) {
		pc = mdcell.sm.m[im].c;
		if (cid.exist(pc->cID)) {
			type_indx = cid.indx(pc->cID);
			if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im + mdcell.mm.n) continue;
			nc++;
		}
	}
	for (im = 0; im < mdcell.pm.n; im++) {
		pc = mdcell.pm.m[im].c;
		if (cid.exist(pc->cID)) {
			type_indx = cid.indx(pc->cID);
			if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im + mdcell.mm.n + mdcell.sm.n) continue;
			nc++;
		}
	}
	icluster.SetArray(nc);
	if (nc == 0) {sprintf(msg, "no interested cluster exists"); show_msg(msg); return false;}
	nc = 0;
	for (im = 0; im < mdcell.mm.n; im++) {
		for (ic = 0; ic < mdcell.mm.m[im].nCluster; ic++) {
			pc = (BASIC_CLUSTER*)(mdcell.mm.m[im].cluster + ic);
			if (cid.exist(pc->cID)) {
				type_indx = cid.indx(pc->cID);
				if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im) continue;
				icluster.m[nc] = pc; nc++;
			}
		}
	}
	for (im = 0; im < mdcell.sm.n; im++) {
		pc = mdcell.sm.m[im].c;
		if (cid.exist(pc->cID)) {
			type_indx = cid.indx(pc->cID);
			if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im + mdcell.mm.n) continue;
			icluster.m[nc] = pc; nc++;
		}
	}
	for (im = 0; im < mdcell.pm.n; im++) {
		pc = mdcell.pm.m[im].c;
		if (cid.exist(pc->cID)) {
			type_indx = cid.indx(pc->cID);
			if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im + mdcell.mm.n + mdcell.sm.n) continue;
			icluster.m[nc] = pc; nc++;
		}
	}
	return true;
}

bool init_interested_cluster(CGMM_MD_CELL &mdcell, ARRAY<short> &cid, ARRAY<int> &mindx, ARRAY<CG_CLUSTER*> &icluster) {
	char msg[256] = "\0";
	int im, ic, type_indx;
	CG_CLUSTER *pc = NULL;
	int nc = 0;
	for (im = 0; im < mdcell.mm.n; im++) {
		for (ic = 0; ic < mdcell.mm.m[im].nCluster; ic++) {
			pc = (CG_CLUSTER*)(mdcell.mm.m[im].cluster + ic);
			if (cid.exist(pc->cID)) {
				type_indx = cid.indx(pc->cID);
				if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im) continue;
				nc++;
			}
		}
	}
	for (im = 0; im < mdcell.sm.n; im++) {
		pc = (CG_CLUSTER*)(mdcell.sm.m + im);
		if (cid.exist(pc->cID)) {
			type_indx = cid.indx(pc->cID);
			if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im + mdcell.mm.n) continue;
			nc++;
		}
	}
	for (im = 0; im < mdcell.pm.n; im++) {
		pc = (CG_CLUSTER*)(mdcell.pm.m + im);
		if (cid.exist(pc->cID)) {
			type_indx = cid.indx(pc->cID);
			if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im + mdcell.mm.n + mdcell.sm.n) continue;
			nc++;
		}
	}
	icluster.SetArray(nc);
	if (nc == 0) {sprintf(msg, "no interested cluster exists"); show_msg(msg); return false;}
	nc = 0;
	for (im = 0; im < mdcell.mm.n; im++) {
		for (ic = 0; ic < mdcell.mm.m[im].nCluster; ic++) {
			pc = (CG_CLUSTER*)(mdcell.mm.m[im].cluster + ic);
			if (cid.exist(pc->cID)) {
				type_indx = cid.indx(pc->cID);
				if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im) continue;
				icluster.m[nc] = pc; nc++;
			}
		}
	}
	for (im = 0; im < mdcell.sm.n; im++) {
		pc = (CG_CLUSTER*)(mdcell.sm.m + im);
		if (cid.exist(pc->cID)) {
			type_indx = cid.indx(pc->cID);
			if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im + mdcell.mm.n) continue;
			icluster.m[nc] = pc; nc++;
		}
	}
	for (im = 0; im < mdcell.pm.n; im++) {
		pc = (CG_CLUSTER*)(mdcell.pm.m + im);
		if (cid.exist(pc->cID)) {
			type_indx = cid.indx(pc->cID);
			if (mindx.m[type_indx] > 0 && mindx.m[type_indx] != im + mdcell.mm.n + mdcell.sm.n) continue;
			icluster.m[nc] = pc; nc++;
		}
	}
	return true;
}

bool within_interested_range(double *v, double *vc, float ir2, float or2) {
	double dr[3] = {v[0] - vc[0], v[1] - vc[1], v[2] - vc[2]};
	double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
	return (d2 >= ir2 && d2 <= or2 ? true : false);
}

bool XYZsave_cluster_within_interested_range(ofstream &out, BASIC_CLUSTER *pc, bool periodic_cell, int *ncell, double *xd, ARRAY<BASIC_CLUSTER*> &icluster, ARRAY<short> &cid, ARRAY<float> &ir2, ARRAY<float> &or2) {
	double r[3], rc[3], dr[3], d2, drcell[3];
	int i, j, k;
	int l, m, n;
	int ic;
	BASIC_CLUSTER *pc1;

	//if (cid.exist(pc->cID)) {
	if (icluster.exist(pc)) {
		/*
		if (periodic_cell) {
			for (i = -ncell[0]; i <= ncell[0]; i++) {drcell[0] = i * xd[0];
			for (j = -ncell[1]; j <= ncell[1]; j++) {drcell[1] = j * xd[1];
			for (k = -ncell[2]; k <= ncell[2]; k++) {drcell[2] = k * xd[2];
				ClusterSaveXYZ(&out, pc, true, drcell);
			}}}
		}
		else ClusterSaveXYZ(&out, pc, false, NULL);
		*/
		 ClusterSaveXYZ(&out, pc, false, NULL);
		return true;
	}


	if (periodic_cell) {
		for (i = -ncell[0] -1; i <= ncell[0] + 1; i++) {drcell[0] = i * xd[0]; r[0] = pc->vp->v[0] + drcell[0];
		for (j = -ncell[1] -1; j <= ncell[1] + 1; j++) {drcell[1] = j * xd[1]; r[1] = pc->vp->v[1] + drcell[1];
		for (k = -ncell[2] -1; k <= ncell[2] + 1; k++) {drcell[2] = k * xd[2]; r[2] = pc->vp->v[2] + drcell[2];
			
			for (ic = 0; ic < icluster.n; ic++) {
				pc1 = icluster.m[ic];
				rc[0] = pc1->vp->v[0];
				rc[1] = pc1->vp->v[1];
				rc[2] = pc1->vp->v[2];
				/*
				for (l = -ncell[0]; l <= ncell[0]; l++) {rc[0] = pc1->vp->v[0] + l * xd[0];
				for (m = -ncell[1]; m <= ncell[1]; m++) {rc[1] = pc1->vp->v[1] + m * xd[1];
				for (n = -ncell[2]; n <= ncell[2]; n++) {rc[2] = pc1->vp->v[2] + n * xd[2];
				*/
				if (within_interested_range(r, rc, ir2.m[ic], or2.m[ic])) {
					ClusterSaveXYZ(&out, pc, true, drcell); goto _out;
				}
				//}}}
			}
_out:       continue;
		}}}
	}
	else {
		for (ic = 0; ic < icluster.n; ic++) {
			pc1 = icluster.m[ic];
			if (within_interested_range(pc->vp->v, pc1->vp->v, ir2.m[ic], or2.m[ic])) {
				ClusterSaveXYZ(&out, pc, false, NULL); break;
			}
		}
	}
	return true;
}

bool XYZsave_cluster_within_interested_range(ofstream &out, CG_CLUSTER *pc, bool periodic_cell, int *ncell, double *xd, ARRAY<CG_CLUSTER*> &icluster, ARRAY<short> &cid, ARRAY<float> &ir2, ARRAY<float> &or2) {
	double r[3], rc[3], dr[3], d2, drcell[3];
	int i, j, k;
	int l, m, n;
	int ic;
	CG_CLUSTER *pc1;

	//if (cid.exist(pc->cID)) {
	if (icluster.exist(pc)) {
		/*
		if (periodic_cell) {
			for (i = -ncell[0]; i <= ncell[0]; i++) {drcell[0] = i * xd[0];
			for (j = -ncell[1]; i <= ncell[1]; i++) {drcell[1] = j * xd[1];
			for (k = -ncell[2]; i <= ncell[2]; i++) {drcell[2] = k * xd[2];
				CGClusterSaveXYZ(&out, pc, true, drcell);
			}}}
		}
		else CGClusterSaveXYZ(&out, pc, false, NULL);
		*/
		CGClusterSaveXYZ(&out, pc, false, NULL);
		return true;
	}


	if (periodic_cell) {
		for (i = -ncell[0] -1; i <= ncell[0] + 1; i++) {drcell[0] = i * xd[0]; r[0] = pc->vp->v[0] + drcell[0];
		for (j = -ncell[1] -1; i <= ncell[1] + 1; i++) {drcell[1] = j * xd[1]; r[1] = pc->vp->v[1] + drcell[1];
		for (k = -ncell[2] -1; i <= ncell[2] + 1; i++) {drcell[2] = k * xd[2]; r[2] = pc->vp->v[2] + drcell[2];
			
			for (ic = 0; ic < icluster.n; ic++) {
				pc1 = icluster.m[ic];
				rc[0] = pc1->vp->v[0];
				rc[1] = pc1->vp->v[1];
				rc[2] = pc1->vp->v[2];
				/*
				for (l = -ncell[0]; l <= ncell[0]; l++) {rc[0] = pc1->vp->v[0] + l * xd[0];
				for (m = -ncell[1]; m <= ncell[1]; m++) {rc[1] = pc1->vp->v[1] + m * xd[1];
				for (n = -ncell[2]; n <= ncell[2]; n++) {rc[2] = pc1->vp->v[2]] + n * xd[2];
				*/

				if (within_interested_range(r, rc, ir2.m[ic], or2.m[ic])) {
					CGClusterSaveXYZ(&out, pc, true, drcell); goto _out;
				}
				//}}}
			}
_out:       continue;
		}}}
	}
	else {
		for (ic = 0; ic < icluster.n; ic++) {
			pc1 = icluster.m[ic];
			if (within_interested_range(pc->vp->v, pc1->vp->v, ir2.m[ic], or2.m[ic])) {
				CGClusterSaveXYZ(&out, pc, false, NULL); break;
			}
		}
	}
	return true;
}

bool SaveConfigXYZ_interested_cluster(char *spar_fname) {
	char msg[256] = "\0", buffer[256] = "\0", ofname[256] = "\0";
	if (!read_item(spar_fname, "[OUTPUT_FILE]", buffer) || sscanf(buffer, "%s", ofname) != 1) {
		sprintf(msg, "failure to get output filename from %s", spar_fname); show_log(msg, true); return false;
	}

	ofstream out;
	out.open(ofname);
	if (!out.is_open()) {
		sprintf(msg, "can not open file %s", ofname); show_log(msg, true); return false;
	}

	ARRAY<short> id_cluster;
	ARRAY<float> ir1, ir2, or1, or2;
	ARRAY<int> mindx;
	ARRAY<BASIC_CLUSTER*> icluster;
	ARRAY<CG_CLUSTER*> icgcluster;
	if (!read_interested_cluster(spar_fname, "[INTERESTED_CLUSTER]", id_cluster, ir1, or1, mindx)) return false;
	if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD || ::MD_MODE == SINGLE_MM_FREE_CELL_MD || ::MD_MODE == MULTI_MM_FREE_CELL_MD) {
		if (!init_interested_cluster(::mm_md_cell, id_cluster, mindx, icluster)) return false;
	}
	else if (::MD_MODE == SINGLE_CGMM_FREE_CELL_MD || ::MD_MODE == MULTI_CGMM_FREE_CELL_MD || ::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
		if (!init_interested_cluster(::cgmm_md_cell, id_cluster, mindx, icgcluster)) return false;
	}

	int dn[3] = {0, 0, 0}, i;
	int ic;
	double dv[3] = {0, 0, 0};
	if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD || ::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
		if (!read_item(spar_fname, "[CELL_EXTENSION]", buffer) || sscanf(buffer, "%d %d %d", dn, dn+1, dn+2) != 3) {
			dn[0] = 0; dn[1] = 0; dn[2] = 0; dv[0] = 0; dv[1] = 0; dv[2] = 0;
		}
	}

	ir2.SetArray(ir1.n); or2.SetArray(or1.n);
	for (i = 0; i < ir1.n; i++) {or2.m[i] = or1.m[i] * or1.m[i]; ir2.m[i] = ir1.m[i] * ir1.m[i];}
	
	GENERAL_MOLECULE gm;
	int imol = 0;
	double r[3];
	if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
		dv[0] = ::mm_md_cell.h[0] * 2;
		dv[1] = ::mm_md_cell.h[1] * 2;
		dv[2] = ::mm_md_cell.h[2] * 2;
	}
	else if (::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
		dv[0] = ::cgmm_md_cell.h[0] * 2;
		dv[1] = ::cgmm_md_cell.h[1] * 2;
		dv[2] = ::cgmm_md_cell.h[2] * 2;
	}
	imol = 0;
	while (1) {
		if (!get_molecule_upto_job(imol, gm)) break;
		if (gm.mol_type == 1) { // macro-molecule
			if (gm.mm == NULL) break;
			for (ic = 0; ic < gm.mm->nCluster; ic++) {
				XYZsave_cluster_within_interested_range(out, (BASIC_CLUSTER*)(gm.mm->cluster + ic), ::mm_md_cell.period, dn, dv, icluster, id_cluster, ir2, or2);
			}
		}
		else if (gm.mol_type == 0) { // simple molecule
			if (gm.sm == NULL) break;
			XYZsave_cluster_within_interested_range(out, gm.sm->c, ::mm_md_cell.period, dn, dv, icluster, id_cluster, ir2, or2);
		}
		else if (gm.mol_type == 4) { // point molecule
			if (gm.pm == NULL) break;
			XYZsave_cluster_within_interested_range(out, gm.pm->c, ::mm_md_cell.period, dn, dv, icluster, id_cluster, ir2, or2);
		}
		else if (gm.mol_type == 2) { // coarse-grained macro-molecule
			if (gm.cgmm == NULL) break;
			for (ic = 0; ic < gm.cgmm->nCluster; ic++) {
				XYZsave_cluster_within_interested_range(out, gm.cgmm->cluster + ic, ::cgmm_md_cell.period, dn, dv, icgcluster, id_cluster, ir2, or2);
			}
		}
		else if (gm.mol_type == 3) { // coarse-grained simple molecule
			if (gm.cgsm == NULL) break;
			XYZsave_cluster_within_interested_range(out, (CG_CLUSTER*)(gm.cgsm), ::cgmm_md_cell.period, dn, dv, icgcluster, id_cluster, ir2, or2);
		}
		imol += 1;
	}
	out.close();
	sprintf(msg, "xyz-struct is saved in file %s", ofname); show_log(msg, true);
	return true;
}

extern "C" bool MultiSaveConfigXYZ_interested_cluster(char *parfname) {
	char buffer[256] = "\0", mdpar_fname[256] = "\0", spar_fname[256] = "\0";
	char msg[256] = "\0";
	int iConfig = 0;
	mlog.init("xyz.log", false);

	if (!read_item(parfname, "[MDPAR_FILE]", buffer) || sscanf(buffer, "%s", mdpar_fname) != 1) {
		sprintf(msg, "mdpar filename is not given after [MDPAR_FILE] in %s", parfname); show_log(msg, true); 
		::mlog.close(); return false;
	}

	int MD_MODE = 0;
	if (!read_item(mdpar_fname, "[MD_CELL]", buffer) || sscanf(buffer, "%d", &MD_MODE) != 1) {
		sprintf(msg, "failure to get the type of MD from %s", mdpar_fname); show_log(msg, true);
		::mlog.close(); return false;
	}
	set_md_mode(MD_MODE);
	if (MD_MODE == SINGLE_MM_FREE_CELL_MD) set_job(GET_SINGLE_MM_MD_PROC);
	else if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) 
		set_job(GET_MULTI_MM_MD_PROC);
	else if (MD_MODE == SINGLE_CGMM_FREE_CELL_MD) set_job(GET_SINGLE_CGMM_MD_PROC);
	else if (MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) 
		set_job(GET_MULTI_CGMM_MD_PROC);
	else {
		sprintf(msg, "the type of MD in %s is UNKNOWN", mdpar_fname); show_log(msg, true);
		::mlog.close(); return false;
	}
	
	ifstream in;
	in.open(parfname);
	if (!in.is_open()) return false;
	if (!search(in, "[XYZ_STRUCT]", buffer)) {
		sprintf(msg, "[XYZ_STRUCT] is not defined in %s", parfname); show_log(msg, true);
		::mlog.close(); in.close(); return false;
	}

	if (!construct_macromolecule(mdpar_fname)) {
		sprintf(msg, "failure to initilize molecule structure from %s", mdpar_fname); show_log(msg, true);
		::mlog.close(); in.close(); return false;
	}
	if (!init_read_mdproc(mdpar_fname)) {
		sprintf(msg, "failure to initilize serial-reading-md-procedure from %s", mdpar_fname); show_log(msg, true);
		::mlog.close(); in.close(); return false;
	}

	if (MD_MODE == SINGLE_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_FREE_CELL_MD ||
		MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			if (!cgatompar_db.get_pdb_alias(parfname, "[CGATOM_PDB_ALIAS]")) return false;
	}

	while (!in.eof()) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %s", &iConfig, spar_fname) != 2) break;
		if (!serial_md_proc(iConfig)) break;
		if (MD_MODE == MULTI_MM_PERIDIC_CELL_MD) molecule_check_periodic_cell(::mm_md_cell);
		else if (MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) _coarse_grain_::molecule_check_periodic_cell(::cgmm_md_cell);
		if (!SaveConfigXYZ_interested_cluster(spar_fname)) break;
	}
	in.close(); ::mlog.close();
	return true;
}
