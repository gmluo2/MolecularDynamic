#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
#include "show.h"
#include "ranlib.h"
#include "vector.h"
#include "prof.h"
#include "dstruct.h"
#include "atom.h"

//#include "bound.h"
//#include "ZMatrix.h"
//#include "MM.h"
//#include "var.h"

//#include "Matrix.h"

#include "sr.h"

char sr_fname[256] = "\0";
extern PAIRWISE_DB< PAIRWISE<EF_DPROF> > sr_db; // pair-wise short-range interaction database
extern bool is_bound(char *buffer, char *atom1, char *atom2, char **pars);

extern void pickup_str(char *source, int nstart, int len, char *dest);
extern bool NextSection(char** original, char* ch);

extern bool read_profile(char *fname, int nIgnoreLines, ARRAY<float*> &va, ARRAY<int> &icol, int& npts);

extern float fUnit_LJ_mass;
extern float eUnit_LJ_kT;
extern float kT, U_MassAccel;

bool find_pair_pars(ifstream &in, char *atom1, char *atom2, char *msg, char *endline) {
	char buffer[256] = "\0", *pars = NULL, *bf = NULL;
	char buff[256] = "\0";
	int n = 0;
	if (endline != NULL) n = (int)strlen(endline);
	while (!in.eof()) {
		in.getline(buffer, 250);
		if (n > 0) {
			bf = buffer; if (!NextSection(&bf, " \t")) continue;
			pickup_str(bf, 0, n, buff);
			if (strcmp(buff, endline) == 0) return false;
		}
		if (is_bound(buffer, atom1, atom2, &pars)) {
			strcpy(msg, pars); return true;
		}
	}
	return false;
}

extern void rewind(ifstream &in);
extern bool NextSection(char** original, char* ch);
extern char seperator[6];
extern int read_float(char *buffer, float *v, int max_v, char* format, char* seperator);
extern bool search(ifstream &in, char *title, char *buffer);

void add_exponential_potential(EF_DPROF &ef_prof, double A, double B) {
	int i;
	double v = 0;
	for (i = 0; i < ef_prof.npts; i++) {
		v = A * exp(-B * ef_prof.x[i]);
		ef_prof.y[0][i] += v * eUnit_LJ_kT;
		ef_prof.y[1][i] += B * v * fUnit_LJ_mass;
	}
}

void add_Buck6_potential(EF_DPROF &ef_prof, double A, double B, double C) {
	int i;
	double vexp = 0, r6;
	for (i = 1; i < ef_prof.npts; i++) {
		r6 = pow(ef_prof.x[i], 6);
		vexp = A * exp(-B * ef_prof.x[i]);
		ef_prof.y[0][i] += (vexp - C / r6) * eUnit_LJ_kT;
		ef_prof.y[1][i] += (B * vexp - C * 6 / (r6 * ef_prof.x[i])) * fUnit_LJ_mass;
	}
	// 0th point has x = 0
	ef_prof.y[0][0] = ef_prof.y[0][1];
	ef_prof.y[1][0] = ef_prof.y[1][1];
}

void add_TF_potential(EF_DPROF &ef_prof, double A, double B, double C6, double C8, double C10) {
	int i;
	double vexp = 0, r2, r6, r8, r10;
	for (i = 1; i < ef_prof.npts; i++) {
		r2 = ef_prof.x[i] * ef_prof.x[i];
		r6 = pow(ef_prof.x[i], 6); r8 = r6 * r2; r10 = r8 * r2;
		vexp = A * exp(-B * ef_prof.x[i]);
		ef_prof.y[0][i] += (vexp - C6 / r6 - C8 / r8 - C10 / r10) * eUnit_LJ_kT;
		ef_prof.y[1][i] += (B * vexp - C6 * 6 / (r6 * ef_prof.x[i]) - C8 * 8 / (r8 * ef_prof.x[i]) - C10 * 10 / (r10 * ef_prof.x[i])) * fUnit_LJ_mass;
	}
	// 0th point has x = 0
	ef_prof.y[0][0] = ef_prof.y[0][1];
	ef_prof.y[1][0] = ef_prof.y[1][1];
}

void add_LJ12_6_potential(EF_DPROF &ef_prof, double eps, double sigma) {
	int i;
	double r6, r12, s, r, v;
	eps *= 4;
	for (i = 1; i < ef_prof.npts; i++) {
		r = ef_prof.x[i]; s = sigma / r;
		r6 = pow(s, 6); r12 = r6 * r6;
		v = eps * (r12 - r6);
		ef_prof.y[0][i] += v * eUnit_LJ_kT;
		ef_prof.y[1][i] += (eps * (12 * r12 / r - 6 * r6 / r)) * fUnit_LJ_mass;
	}
	// 0th point has x = 0
	ef_prof.y[0][0] = ef_prof.y[0][1];
	ef_prof.y[1][0] = ef_prof.y[1][1];
}

void add_ALJ12_6_potential(EF_DPROF &ef_prof, double eps, double sigma) {
	int i;
	double r6, r12, s, r, v;
	eps *= 4;
	for (i = 1; i < ef_prof.npts; i++) {
		r = ef_prof.x[i]; s = sigma / r;
		r6 = pow(s, 6); r12 = r6 * r6;
		v = eps * (r12 - 2 * r6);
		ef_prof.y[0][i] += v * eUnit_LJ_kT;
		ef_prof.y[1][i] += (12 * eps * (r12 - r6) / r) * fUnit_LJ_mass;
	}
	// 0th point has x = 0
	ef_prof.y[0][0] = ef_prof.y[0][1];
	ef_prof.y[1][0] = ef_prof.y[1][1];
}

void add_GLJ12_6_potential(EF_DPROF &ef_prof, double A, double B) {
	int i;
	double r6, r12, r, v;
	for (i = 1; i < ef_prof.npts; i++) {
		r = ef_prof.x[i]; 
		r6 = pow(r, 6); r12 = r6 * r6;
		v = A / r12 - B / r6;
		ef_prof.y[0][i] += v * eUnit_LJ_kT;
		ef_prof.y[1][i] += ((12 * A / r12 - 6 * B / r6) / r) * fUnit_LJ_mass;
	}
	// 0th point has x = 0
	ef_prof.y[0][0] = ef_prof.y[0][1];
	ef_prof.y[1][0] = ef_prof.y[1][1];
}

bool construct_short_range_interaction_prof(ifstream &in, char *atom1, char *atom2, PAIRWISE<EF_DPROF> &sr) {
	char buffer[256] = "\0", buff[256] = "\0", msg[256] = "\0", fname[256] = "\0", *bf = NULL;
	int type = 0;
	bool status = false, prof_status = false;
	ARRAY<float> var;
	while (!in.eof()) {
		status = find_pair_pars(in, atom1, atom2, buffer, "[");
		// the energy unit here is kJ/mol
		if (!status) break;
		if (sscanf(buffer, "%d", &type) != 1) {
			sprintf(msg, "UNKNOWN type of short-range interaction, see : %s-%s %s", atom1, atom2, buffer); show_log(msg, true);
			continue;
		}
		else {
			bf = buffer;
			if (!NextSection(&bf, seperator)) {
				sprintf(msg, "UNKNOWN format for short-range interaction, see : %s-%s %s", atom1, atom2, buffer); show_log(msg, true);
				continue;
			}
			switch (type) {
			case _SR_EXP_ :
				var.SetArray(_N_EXP_);
				if (read_float(bf, var.m, var.n, "%f", seperator) != var.n) {
					sprintf(msg, "variables for exponential pairwise potential [A, B]is not given, see : %s-%s %s", atom1, atom2, buffer);
					show_log(msg, true); continue;
				}
				else {add_exponential_potential(sr.pv, var.m[0], var.m[1]); prof_status = true;}
				break;
			case _SR_BUCK6_ :
				var.SetArray(_N_BUCK6_);
				if (read_float(bf, var.m, var.n, "%f", seperator) != var.n) {
					sprintf(msg, "variables for Buck6 pairwise potential [A, B, C]is not given, see : %s-%s %s", atom1, atom2, buffer);
					show_log(msg, true); continue;
				}
				else {add_Buck6_potential(sr.pv, var.m[0], var.m[1], var.m[2]); prof_status = true;}
				break;
			case _SR_TF_ :
				var.SetArray(_N_TF_);
				if (read_float(bf, var.m, var.n, "%f", seperator) != var.n) {
					sprintf(msg, "variables for Tosi-Fumi pairwise potential [A, B, C6, C8, C10]is not given, see : %s-%s %s", atom1, atom2, buffer);
					show_log(msg, true); continue;
				}
				else {add_TF_potential(sr.pv, var.m[0], var.m[1], var.m[2], var.m[3], var.m[4]); prof_status = true;}
				break;
			case _SR_LJ_ :
				var.SetArray(_N_LJ_);
				if (read_float(bf, var.m, var.n, "%f", seperator) != var.n) {
					sprintf(msg, "variables for Lennard-Jones pairwise potential [eps, sigma] is not given, see : %s-%s %s", atom1, atom2, buffer);
					show_log(msg, true); continue;
				}
				else {add_LJ12_6_potential(sr.pv, var.m[0], var.m[1]); prof_status = true;}
				break;
			case _SR_ALJ_ :
				var.SetArray(_N_ALJ_);
				if (read_float(bf, var.m, var.n, "%f", seperator) != var.n) {
					sprintf(msg, "variables for alternative format of Lennard-Jones pairwise potential [eps, sigma] is not given, see : %s-%s %s", atom1, atom2, buffer);
					show_log(msg, true); continue;
				}
				else {add_ALJ12_6_potential(sr.pv, var.m[0], var.m[1]); prof_status = true;}
				break;
			case _SR_GLJ_ :
				var.SetArray(_N_GLJ_);
				if (read_float(bf, var.m, var.n, "%f", seperator) != var.n) {
					sprintf(msg, "variables for general Lennard-Jones pairwise potential [A, B] is not given, see : %s-%s %s", atom1, atom2, buffer);
					show_log(msg, true); continue;
				}
				else {add_GLJ12_6_potential(sr.pv, var.m[0], var.m[1]); prof_status = true;}
				break;
			case _SR_FILE_:
				if (sscanf(bf, "%s", fname) != 1 || strlen(fname) == 0) {
					sprintf(msg, "filename for short-range interaction [A, B] is not given, see : %s-%s %s", atom1, atom2, buffer);
					show_log(msg, true); continue;
				}
				{
					ARRAY<float*> fp;
					float *x = NULL, *y1 = NULL, *y2 = NULL;
					fp.set_array(3); fp.m[0] = NULL; fp.m[1] = NULL; fp.m[2] = NULL;
					ARRAY<int> col;
					col.set_array(3); col.m[0] = 0; col.m[1] = 1; col.m[2] = 2;
					int npts = 0, i, ti;
					double t, funit = kT / U_MassAccel;
					EF_DPROF prof;
					if (!read_profile(fname, 0, fp, col, npts) || npts == 0) {
						sprintf(msg, "failure to read short-range interfaction from file %s for %s - %s: %s", fname, atom1, atom2, buffer);
						show_log(msg, true);
					}
					else {
						x = fp.m[0]; y1 = fp.m[1]; y2 = fp.m[2];
						prof.set(npts, x[1] - x[0], 0);
						for (i = 0; i < npts; i++) {
							prof.x[i] = x[i]; 
							prof.y[0][i] = y1[i]; 
							prof.y[1][i] = y2[i] * funit;
						}
						for (i = 0; i < sr.pv.npts; i++) {
							interpolate2(sr.pv.y[0][i], sr.pv.y[1][i], prof, sr.pv.x[i], ti, t)
						}
						for (i = 0; i < sr.pv.npts; i++) {
							sr.pv.y[0][i] -= sr.pv.y[0][sr.pv.npts - 5];
							sr.pv.y[1][i] -= sr.pv.y[1][sr.pv.npts - 5];
						}
						sr.pv.y[0][0] = sr.pv.y[0][1];
						sr.pv.y[1][0] = sr.pv.y[1][1];
						sr.pv.y[0][sr.pv.npts-4] = 0;
						sr.pv.y[1][sr.pv.npts-4] = 0;
						sr.pv.y[0][sr.pv.npts-3] = 0;
						sr.pv.y[1][sr.pv.npts-3] = 0;
						sr.pv.y[0][sr.pv.npts-2] = 0;
						sr.pv.y[1][sr.pv.npts-2] = 0;
						sr.pv.y[0][sr.pv.npts-1] = 0;
						sr.pv.y[1][sr.pv.npts-1] = 0;

						for (i = 0; i < fp.n; i++) {if (fp.m[i] != NULL) {delete[] fp.m[i]; fp.m[i] = NULL;}};
						prof_status = true;
					}
				}
				break;
			default:
				sprintf(msg, "UNKNOWN type of short-range interaction %d, see : %s-%s %s", atom1, atom2, buffer); show_log(msg, true);
				break;
			}
		}
	}
	return prof_status;
}

bool initPairwiseShortRangeInteract(char *par_fname, ATOM_COMM_PARS_DATABASE &atompar_db, PAIRWISE_DB< PAIRWISE<EF_DPROF> > &sr_db, float &rcut, float &rmin) {
	char buffer[256] = "\0", buff[256] = "\0", msg[256] = "\0";
	float v1, v2;
	sr_db.releaseChain(true);
	ifstream in;
	in.open(par_fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open %s for short-range interaction.\nShort-range interaction is disabled", par_fname); show_log(msg, true); return false;
	}
	if (!search(in, "[R_RANGE]", buffer)) {
		sprintf(msg, "[R_RANGE] can not be found in %s", par_fname); show_log(msg, true); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f", &v1, &v2) != 2) {
		sprintf(msg, "can not get cut distance and the minimum distance to flat the potential, following [R_RANGE], see :"); show_log(msg, true);
		show_log(buffer, true); in.close(); return false;
	}
	if (v1 < 0) {
		sprintf(msg, "the cut distance can not be less than 0, rcut = %f", v1); show_log(msg, true); 
		show_log("distable the short-range pairwise interaction", true); in.close(); return false;
	}
	if (v1 <= v2) {
		sprintf(msg, "the parameters about range is given with cut-distance first, while the 2nd one is the flatting distance in close-distance");
		show_log(msg, true);
		sprintf(msg, "However, the two parameters is switched with cutting distance %f Angs, while flatting close-distance %f Angs", v2, v1);
		show_log(msg, true);
		rcut = v2; v2 = v1; v1 = rcut;
	}
	if (v2 < 0) v2 = 0;
	rcut = v1; rmin = v2;
	float dr = 0.02f;
	int npt_add = 2;
	int npts = int(rcut / dr + 0.5) + npt_add;
	int i, nmin = int(rmin / dr + 0.5);

	CHAIN<ATOM_COMM_PARS> *ch1, *ch2;
	ATOM_COMM_PARS *apar1, *apar2;

	char aindx[2] = {0x00, 0x00};
	unsigned int key = 0;
	PAIRWISE<EF_DPROF> *pw = NULL;
	bool status;

	ch1 = atompar_db.ch;
	while (ch1 != NULL) {
		apar1 = ch1->p;
		ch2 = atompar_db.ch;
		while (ch2 != NULL) {
			apar2 = ch2->p;
			aindx[0] = apar1->indx; aindx[1] = apar2->indx;
			key = construct_char2_key_pairwise(aindx);
			pw = sr_db.search_chain(key);
			if (pw != NULL) {ch2 = ch2->next; continue;}

			rewind(in);
			if (!search(in, "[SHORT_RANGE_POTENTIAL]", buffer)) {
				sprintf(msg, "[SHORT_RANGE_POTENTIAL] is not defined in %s", par_fname); show_log(msg, true); in.close(); 
				sr_db.releaseArray(true); sr_db.releaseChain(true); return true;
			}
			pw = new PAIRWISE<EF_DPROF>;
			pw->pv.set(npts, dr, 0);
			memset(pw->pv.y[0], 0, pw->pv.npts * SIZE_DOUBLE);
			memset(pw->pv.y[1], 0, pw->pv.npts * SIZE_DOUBLE);
			status = construct_short_range_interaction_prof(in, apar1->atom, apar2->atom, *pw);
			if (status) {
				for (i = 0; i < nmin; i++) {
					pw->pv.y[0][i] = pw->pv.y[0][nmin];
					pw->pv.y[1][i] = pw->pv.y[1][nmin];
				}
				pw->aindx[0] = apar1->indx; pw->aindx[1] = apar2->indx;
				pw->key = construct_char2_key_pairwise(pw->aindx);
				for (i = 1; i <= npt_add; i++) {
					pw->pv.y[0][pw->pv.npts - i] = 0;
					pw->pv.y[1][pw->pv.npts - i] = 0;
				}
				sr_db.attach(pw);
				apar1->bSR = true;
				apar2->bSR = true;

				{
					ofstream out;
					sprintf(buffer, "%s-%s-sr.dat", apar1->atom, apar2->atom);
					out.open(buffer);
					for (i = 0; i < pw->pv.npts; i++) out<<pw->pv.x[i]<<"  "<<pw->pv.y[0][i]<<"  "<<pw->pv.y[1][i]<<endl;
					out.close();
					sprintf(msg, "short-range potential [%s - %s] is saved in file %s", apar1->atom, apar2->atom, buffer);
					show_log(msg, true);
				}
			}
			else delete pw;
			pw = NULL;

			ch2 = ch2->next;
		}
		ch1 = ch1->next;
	}
	sr_db.Chain2Array();
	show_log("", true);
	return true;
}

extern bool search(ifstream &in, char *title, char *buffer, char *stop_line);
extern char force_field_fname[256];
bool update_atomic_pars(char *par_fname, ATOM_COMM_PARS_DATABASE &atompar_db) {
	char msg[256] = "\0", buffer[256] = "\0", buff[256] = "\0";
	float v = 0;
	int n = 0;
	ifstream in;
	in.open(par_fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open file %s", par_fname); show_log(msg, true); return false;
	}
	CHAIN<ATOM_COMM_PARS> *ch;
	ATOM_COMM_PARS *apar;

	rewind(in);
	if (search(in, "[DISABLE_PREVIOUS_LJ]", buffer)) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d", &n) == 1 && n == 1) { // disable the Lennard-Jones interaction of all the atom
			ch = atompar_db.ch;
			while (ch != NULL) {
				apar = ch->p;
				apar->bLJ = false;
				ch = ch->next;
			}
			sprintf(msg, "Lennard-Jones potential from force-field %s are disabled", ::force_field_fname); show_log(msg, true);
		}
	}

	rewind(in);
	if (search(in, "[POLARIZABILITY]", buffer)) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d", &n) == 1 && n == 1) { // update the polarizability
			ch = atompar_db.ch;
			while (ch != NULL) {
				apar = ch->p;
				rewind(in);
				if (!search(in, "[POLARIZABILITY]", buffer)) break;
				if (search(in, apar->atom, buffer, "[") && sscanf(buffer, "%s %f", buff, &v) == 2) {
					apar->alpha = v;
					sprintf(msg, "polarizability of %s is updated to %f", apar->atom, v); show_log(msg, true);
				}
				ch = ch->next;
			}
		}
	}
	in.close(); return true;
}
