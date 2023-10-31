#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "ranlib.h"

#include "def.h"
#include "show.h"
#include "dstruct.h"
#include "vector.h"

#include "prof.h"
#include "read.h"

#include "pmf.h"

extern float kT, U_MassAccel;

extern bool read_chars(char **buffer, char *str, char seperator);

bool get_bound(char *buffer, char *atom1, char *atom2, char **pars) {
	char *bf = buffer;
	char atom[2][50] = {"\0", "\0"};
	if (read_chars(&bf, atom[0], '-') && read_chars(&bf, atom[1], '-')) {
		strcpy(atom1, atom[0]); strcpy(atom2, atom[1]);
		if (strlen(bf) > 2) {
			if (bf[0] == '-' || bf[1] == '-') return false;
		}
		*pars = bf; return true;
	}
	return false;
}

namespace _pmf_constraint_ {
void PMF_CONSTRAINT::construct_constraint_efprofile() {
	float fUnit_kT_mass = ::kT / ::U_MassAccel;

	float dx = 0.02f;
	int i, n = int(dr_max / dx + 0.5f);
	int npts = 2 * n + 1;
	this->dr_max = dx * n;
	ef_prof.release();
	ef_prof.set(npts, dx, -(this->dr_max));
	for (i = 0; i < ef_prof.npts; i++) {
		ef_prof.x[i] = (i - n) * dx;
		if (cModel == PMF_CONSTRAINT_0 || cModel == PMF_CONSTRAINT_2) {
			ef_prof.y[0][i] = 0.5 * v_f.m[KC] * ef_prof.x[i] * ef_prof.x[i];
			ef_prof.y[1][i] = -v_f.m[KC] * ef_prof.x[i] * fUnit_kT_mass;
		}
		else if (cModel == PMF_CONSTRAINT_1 || cModel == PMF_CONSTRAINT_3) {
			ef_prof.y[0][i] = 0.5 * v_f.m[KC] * ef_prof.x[i] * ef_prof.x[i] + v_f.m[A] * ef_prof.x[i];
			ef_prof.y[1][i] = (-v_f.m[KC] * ef_prof.x[i] + v_f.m[A]) * fUnit_kT_mass;
		}
		else {
			ef_prof.y[0][i] = 0;
			ef_prof.y[1][i] = 0;
		}
	}
	/*
	ofstream out;
	out.open("PMF_EF.dat");
	for (i = 0; i < ef_prof.npts; i++) out<<ef_prof.x[i]<<"  "<<ef_prof.y[0][i]<<"  "<<ef_prof.y[1][i]<<endl;
	out.close();
	*/
}

void PMF_CONSTRAINT::ef(float R, double &U, double &f) {
	int n, i;
	double t;

	double dR = R - r0;
	if (dR <= -dr_max) {U = ef_prof.y[0][0]; f = ef_prof.y[1][0];}
	else if (dR >= dr_max) {
		n = ef_prof.npts - 1;
		U = ef_prof.y[0][n]; f = ef_prof.y[1][n];
	}
	else {
		interpolate2(U, f, ef_prof, dR, i, t)
	}
	return;
}

bool read_pmf_constraint(char *par_fname, PMF_CONSTRAINT &cpmf) {
	char buffer[256] = "\0", msg[256] = "\0";

	int model = 0;
	if (!read_msg(par_fname, buffer, "[PMF_CONSTRAINT_MODEL]", true)) return false;
	if (sscanf(buffer, "%d", &model) != 1) {
		sprintf(msg, "failure to read PMF constraint model, see: %s", buffer); show_log(msg, true);
		return false;
	}
	cpmf.set_model(model);

	char atom1[50] = "\0", atom2[50] = "\0";
	char *par = NULL;
	int n1 = 0, n2 = 0;
	float v[5];
	float dr_max = 0, r0 = 0;

	if (!read_msg(par_fname, buffer, "[PMF_CONSTRAINT]", true)) return false;
	if (model == PMF_CONSTRAINT_0 || model == PMF_CONSTRAINT_1 || model == PMF_CONSTRAINT_2 || model == PMF_CONSTRAINT_3) { // constraint between two PMs
		if (get_bound(buffer, atom1, atom2, &par)) {
			if (sscanf(atom1, "%d", &n1) == 1 && sscanf(atom2, "%d", &n2) == 1) {
				if (model == PMF_CONSTRAINT_0 || model == PMF_CONSTRAINT_2) {
					if (sscanf(par, "%f %f %f", &r0, &dr_max, v) == 3) {
						cpmf.v_int.m[0] = n1; cpmf.v_int.m[1] = n2; 
						cpmf.v_f.m[KC] = v[0]; cpmf.r0 = r0; cpmf.dr_max = dr_max;
					}
					else {
						sprintf(msg, "ERROR: unknown format for PMF constraint, see : %s", buffer); show_log(msg, true);
						return false;
					}
				}
				else if (model == PMF_CONSTRAINT_1 || model == PMF_CONSTRAINT_3) {
					if (sscanf(par, "%f %f %f %f", &r0, &dr_max, v, v + 1) == 4) {
						cpmf.v_int.m[0] = n1; cpmf.v_int.m[1] = n2; 
						cpmf.v_f.m[KC] = v[0]; cpmf.v_f.m[A] = v[1];
						cpmf.r0 = r0; cpmf.dr_max = dr_max;
					}
					else {
						sprintf(msg, "ERROR: unknown format for PMF constraint, see : %s", buffer); show_log(msg, true);
						return false;
					}
				}
				else {
					sprintf(msg, "ERROR: unknown PMF constraint model, %d", model); show_log(msg, true);
					return false;
				}
			}
			else {
				sprintf(msg, "PM indx for constraint between PMs are not defined as number, see : %s", buffer); show_log(msg, true);
				return false;
			}
		}
		else {
			sprintf(msg, "ERROR: unknown PMF constraint definition [PM_indx1] - [PM_indx2]  [variables]), see : %s", buffer); show_log(msg, true);
			return false;
		}
	}
	else {
		sprintf(msg, "ERROR: Unknown contraint model for PMF : %d", model); show_log(msg, true);
		return false;
	}

	return true;
}

} // end of namespace _pmf_constraint_ 

extern LOG mlog;
extern "C" bool cal_biased_potential(char *var) {
	char fname[256] = "\0", ofname[256] = "\0", par_fname[256] = "\0", msg[256] = "\0", buffer[256] = "\0";
	mlog.init("pmf.log");

	using namespace _pmf_constraint_;
	if (sscanf(var, "%s %s %s", fname, par_fname, ofname) != 3) {
		sprintf(msg, "variable should be [fname] [PMF par_fname] [output_fname], see : %s", var); 
		show_msg(msg, true); mlog.close(); return false;
	}
	PMF_CONSTRAINT cpmf;
	if (!read_pmf_constraint(par_fname, cpmf)) {mlog.close(); return false;}
	cpmf.construct_constraint_efprofile();
	
	ifstream in;
	ofstream out;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open file %s", fname); show_log(msg, true); mlog.close(); return false;
	}
	out.open(ofname);
	if (!out.is_open()) {
		sprintf(msg, "failure to open file %s", ofname); show_log(msg, true); mlog.close(); return false;
	}
	float x = 0;
	double U = 0, f = 0;
	while (!in.eof()) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f", &x) == 1) {
			cpmf.ef(x, U, f);
			out<<buffer<<"  "<<U<<endl;
		}
	}
	in.close(); out.close();
	mlog.close();
	return true;
}
