#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
#include "show.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"

#include "ZMatrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "cg-mm.h"
#include "Mol.h"
#include "CMM.h"
#include "cluster.h"
#include "MD.h"
#include "Interaction.h"
#include "var.h"

#include "MD_POLYMER.h"
#include "test.h"

using namespace _coarse_grain_;

extern void rewind(ifstream& in);
// maximum distance for L-J interaction cut-off
float rLJ_cut_max() {
	float rcut = LJPARS.m[0].m[0].rcut;
	int i = 0, j = 0;
	for (i = 0; i < LJPARS.n; i++) {
		for (j = 0; j < LJPARS.n; j++) {
			if (rcut < LJPARS.m[i].m[j].rcut) rcut = LJPARS.m[i].m[j].rcut;
	}}
	return rcut;
}

void init_LJ_Pars() {
// eps -- kJ/mol; r0 in Angs., rcut
	float rcut = 15;
	int n1 = 0, n2 = 0;
	float rLJ = 0, epsLJ = 0;
	n1 = number<ATOM_COMM_PARS>(atompar_db.ch);
	LJPARS.set(n1);
	CHAIN<ATOM_COMM_PARS> *ch1 = NULL, *ch2 = NULL;

	ch1 = atompar_db.ch; n1 = 0; 
	ch2 = NULL; n2 = 0;
	while (ch1 != NULL) {
		ch2 = atompar_db.ch; n2 = 0;
		if (ch1->p->epsLJ > 0.00001f) {
			while (ch2 != NULL) {
				if (ch2->p->epsLJ > 0.00001f) {
					//epsLJ = 0; // disable LJ interaction
					epsLJ = (float)sqrt(ch1->p->epsLJ * ch2->p->epsLJ);
					rLJ = (ch1->p->rLJ + ch2->p->rLJ) / 2;
					if (::iFormat_LJ == 0) {// 4 * U0 * [(sigma/r)^12 - (sigma/r)^6], rm = 1.122462 * sigma
						rcut = ::rcut_LJ_norm * rLJ * 1.122462f;
					}
					else if (::iFormat_LJ == 1) {// U0 * [(rm/r)^12 - 2 * (rm/r)^6]
						rcut = ::rcut_LJ_norm * rLJ;
					}
					LJPARS.m[n1].m[n2].set_pars(epsLJ, rLJ, rcut);
				}
				ch2 = ch2->next; n2++;
			}
		}
		ch1 = ch1->next; n1++;
	}

	::rcut_LJ = rLJ_cut_max();
	::r2cut_LJ = ::rcut_LJ * ::rcut_LJ;
}

extern CGBOUND_DB cgbound_db;
void init_LJ_Pars_CG() {
// eps -- kJ/mol; r0 in Angs., rcut
	float rcut = 15;
	int n1 = 0, n2 = 0;
	float rLJ = 0, epsLJ = 0;
	n1 = number<CGATOM_COMM_PARS>(cgatompar_db.ch);
	LJPARS.set(n1);
	CHAIN<CGATOM_COMM_PARS> *ch1 = NULL, *ch2 = NULL;

	ch1 = cgatompar_db.ch; n1 = 0; 
	ch2 = NULL; n2 = 0;
	while (ch1 != NULL) {
		ch2 = cgatompar_db.ch; n2 = 0;
		while (ch2 != NULL) {
			if (ch1->p->epsLJ > 0.00001f && ch2->p->epsLJ > 0.00001f) {
				//epsLJ = 0; // disable LJ interaction
				epsLJ = (float)sqrt(ch1->p->epsLJ * ch2->p->epsLJ);
				rLJ = (ch1->p->rLJ + ch2->p->rLJ) / 2;
				rcut = 3.5f * rLJ;
				LJPARS.m[n1].m[n2].set_pars(epsLJ, rLJ, rcut);
			}
			ch2 = ch2->next; n2++;
		}
		ch1 = ch1->next; n1++;
	}

	::rcut_LJ = rLJ_cut_max();
	::r2cut_LJ = ::rcut_LJ * ::rcut_LJ;
}

void init_ImplicitSolvation_Pars() {
	char msg[256] = "\0";
// eps -- kJ/mol; r0 in Angs., rcut
	readImplicitSolvationPars(ESolv_db, ::rSolvent);
	calImplicitSolvationArea(ESolv_db, ::rSolvent);
	CHAIN<CLUSTER_IMPSOLV_PARS> *ch = ESolv_db.ch;
	CLUSTER_IMPSOLV_PARS *psolv = NULL;
	sprintf(msg, "--------------  Implicit Solvation Energy --------------", ::rSolvent); show_log(msg, true);
	sprintf(msg, "Probe [solvent] radius : %8.2f", ::rSolvent); show_log(msg, true);
	sprintf(msg, "--------------------------------------------------------"); show_log(msg, true);
	float rmax = 0;
	while (ch != NULL) {
		psolv = ch->p;
		if (psolv->r0 > rmax) rmax = psolv->r0;

		sprintf(msg, "cluster id : %d", psolv->cID); show_log(msg, true);
		sprintf(msg, "Solv. Energy of unit area [kT] : %6.3E", psolv->sigma); show_log(msg, true);
		sprintf(msg, "Surface area : %8.2f", psolv->S0); show_log(msg, true);
		sprintf(msg, "Total solvation energy [kT] : %8.3f", psolv->S0 * psolv->sigma); show_log(msg, true); 
		sprintf(msg, "--------------------------------------------------------"); show_log(msg, true);

		ch = ch->next;
	}
	::dSolv_max = rmax + rmax + ::rSolvent * 2;
	//::dSolv_max += rmax; // need this ?
}

void MD(MMOLECULE &mm, LFVERLET_MM_MD &mdpar, CMM_CELL3D<BASIC_CLUSTER> *cmm, int interact_method, int nNormalMD) {
	mdpar.resetClusters(mm.nCluster);
	if (interact_method == _CMM) {
		cmm_init_distribute_cluster(mm, *cmm, true); // set base cluster at the center cell
		cmm_init_subcell<BASIC_CLUSTER>(*cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
	}
	if (mm.free_base) MD_POLYMER(mm, mdpar, cmm, interact_method, nNormalMD);
	else FIXEDBASE_MD_LPOLYMER(mm, mdpar, cmm, interact_method);
}

extern void setall(long iseed1,long iseed2);
extern void DefCommands();
extern void DelCommands();
extern void clear_memory();
extern void construct_profiles();
//extern void test();

extern char struct_fname[256], force_field_fname[256], atom_index_fname[256], sr_fname[256];
extern char cluster_solvation_energy_fname[256];
extern char monomer_infor_fname[256];
extern char cgatom_scatt_fname[256];
extern bool search(ifstream &in, char *title, char *buffer);

extern char prof_path[256];
void append_msg(char *msg, char *add_msg, bool endline = true) {
	strcat(msg, add_msg); 
	if (endline) strcat(msg, "\n");
}

extern "C" bool init_env() {
	// some parameter filename
	char msg[2560] = "\0", errmsg[256] = "\0";
	bool status = true;
	char buffer[256] = "\0";
	ifstream in;
	in.open("NEIMO.INI");
	if (!in.is_open()) {sprintf(errmsg, "can not open NEIMO.INI for database filename"); status = false; return status;}
	if (!search(in, "[FORCE_FIELD]", buffer)) {
		rewind(in);
		if (!search(in, "[AMBER_FORCE_FIELD]", buffer)) {
			sprintf(errmsg, "can not get amber force field item"); append_msg(msg, errmsg, true); status = false;
		}
	}
	
	in.getline(buffer, 250);
	if (sscanf(buffer, "%s", force_field_fname) != 1) {
		sprintf(errmsg, "can get force field filename"); append_msg(msg, errmsg, true); status = false;
	}

	int n = 0;
	rewind(in);
	if (!search(in, "[LENNARD_JONES]", buffer)) {
		sprintf(errmsg, "can not get item about Lennard-Jones"); append_msg(msg, errmsg, true); status = false;
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d", &n) != 1) {
			sprintf(errmsg, "can get Lennard-Jones interaction format"); append_msg(msg, errmsg, true); status = false;
		}
		else iFormat_LJ = n;
	}

	rewind(in);
	if (!search(in, "[SHORT_RANGE_INTERACTION]", buffer)) {
		strcpy(sr_fname, "\0");
		sprintf(errmsg, "additional short-range interaction is disabled"); append_msg(msg, errmsg, true);  
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", sr_fname) != 1) {
			sprintf(errmsg, "filename for short-range interaction can not be obtained following [SHORT_RANGE_INTERACTION]"); append_msg(msg, errmsg, true); status = false;
		}
	}

	rewind(in);
	if (!search(in, "[STRUCT]", buffer)) {
		sprintf(errmsg, "can not get cluster struct definition item"); append_msg(msg, errmsg, true); status = false;
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", struct_fname) != 1) {
			sprintf(errmsg, "can get cluster structure definition filename"); append_msg(msg, errmsg, true); status = false;
		}
	}

	rewind(in);

	if (!search(in, "[ATOM_INDEX]", buffer)) {
		sprintf(errmsg, "can not get atom definition for view"); append_msg(msg, errmsg, true); status = false;
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", atom_index_fname) != 1) {
			sprintf(errmsg, "can not get atom definition filename"); append_msg(msg, errmsg, true); status = false;
		}
	}

	rewind(in);

	if (!search(in, "[SOLVATION_ENERGY]", buffer)) {
		sprintf(errmsg, "can not get cluster-solvation-energy definition for view"); append_msg(msg, errmsg, true); status = false;
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", cluster_solvation_energy_fname) != 1) {
			sprintf(errmsg, "can not get filename for solvation energy"); append_msg(msg, errmsg, true); status = false;
		}
	}

	rewind(in);

	if (!search(in, "[MONOMER_INFOR]", buffer)) {
		sprintf(errmsg, "no additional monomer information"); append_msg(msg, errmsg, true); 
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", monomer_infor_fname) != 1) {
			strcpy(monomer_infor_fname, "\0");
			sprintf(errmsg, "no additional monomer information"); append_msg(msg, errmsg, true); 
		}
	}

	rewind(in);
	if (!search(in, "[CG_FORCE_FIELD]", buffer)) {
		sprintf(errmsg, "force-field of coarse-grained atom is not given"); append_msg(msg, errmsg, true); 
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", cgatom_par_fname) != 1) {
			strcpy(cgatom_par_fname, "\0");
			sprintf(errmsg, "filename of the force-field of coarse-grained atom is not given"); append_msg(msg, errmsg, true); 
		}
	}

	rewind(in);
	if (!search(in, "[CGATOM_INDEX]", buffer)) {
		sprintf(errmsg, "index of coarse-grained atom is not given"); append_msg(msg, errmsg, true); 
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", cgatom_index_fname) != 1) {
			strcpy(cgatom_index_fname, "\0");
			sprintf(errmsg, "filename of the index of coarse-grained atom is not given"); append_msg(msg, errmsg, true); 
		}
	}

	rewind(in);
	if (!search(in, "[CGATOM_SCATT_FACTOR]", buffer)) {
		sprintf(errmsg, "[CG_SCATT_FACTOR] -- filename for scattering factor of coarse-grained atom is not given"); append_msg(msg, errmsg, true); 
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", cgatom_scatt_fname) != 1) {
			strcpy(cgatom_scatt_fname, "\0");
			sprintf(errmsg, "filename of scattering factor of coarse-grained atom is not given"); append_msg(msg, errmsg, true); 
		}
	}

	rewind(in);
	if (!search(in, "[SPECIAL_FUNCTION_PATH]", buffer)) strcpy(prof_path, "\0");
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", prof_path) != 1) strcpy(prof_path, "\0");
		if (strlen(prof_path) > 0) {
			{
				int n = strlen(prof_path);
#if _SYS_ == _WINDOWS_SYS_
				if (prof_path[n - 1] != '\\') {prof_path[n] = '\\'; prof_path[n+1] = '\0';}
#elif _SYS_ == _LINUX_SYS_
				if (prof_path[n - 1] != '/') {prof_path[n] = '/'; prof_path[n+1] = '\0';}
#endif
			}
		}
	}
	
	in.close();

	if (strlen(msg) > 0) show_msg(msg);

#ifdef __GPU__
	//cudaSetDevice(_CudaDevice_);
	//cudaDeviceReset();
	//show_log("use one GPU device only!", true);
#endif

#if _MPI_ == 1
	long sd1 = (::mpi_id + 1) * 4897654;
	long sd2 = (::mpi_id + 1) * 6746156;
	setall(sd1, sd2); // random number
#else
	setall(45649741, 546550948); // random number
	//setall(546550948, 546550948); // random number
#endif
	DefCommands(); // command
	// normally, disable test function
	//test();

	return status;
}

#include "time.h"
#include "ranlib.h"
/*
void reset_random_seeds() {
	TIME t;
	long tsec, tms;
	t.get_second(tsec, tms);
	tsec = tsec - tsec / 1000 * 1000;
	int i, s1, s2;
	for (i = 0; i < tsec; i++) s1 = ranf();
	for (i = 0; i < tms; i++) s2 = ranf();
	setall(long(546550948 * s1) + 45649741, long(546550948 * s2) + 68756456);
}
*/

#include "time.h"

void reset_random_seeds() {
	long tsec, tns;
#if _SYS_ == _WINDOWS_SYS_
	_timeb time;
	_ftime_s(&time);
	tsec = time.time;
	tns = time.millitm;
#elif _SYS_ == _LINUX_SYS_
	//get_clocktime(tsec, tns);
	struct timeb time;
	ftime(&time);
	tsec = time.time;
	tns = time.millitm;
#endif
	int i, s1, s2;
	long seed1, seed2;
	int nr1 = tsec, nr2 = tns;
	while (nr1 > 20480) nr1 /= 100;
	while (nr2 > 20480) nr2 /= 100;
	for (i = 0; i < nr1; i++) s1 = ranf();
	for (i = 0; i < nr2; i++) s2 = ranf();

	if (tns > 1000000) tns -= tns / 100000 * 100000;

	seed1 = long(546948l * s1) + tns;
	seed2 = long(283959l * s2) + tns;

	while (seed1 < 34597590l) seed1 += tns;
	while (seed2 < 34534590l) seed2 += tns;

	setall(seed1, seed2);
}


extern "C" void clear_env() {
	DelCommands();
	clear_memory();
}
