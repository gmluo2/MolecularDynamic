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

#include "atom.h"
#define _INTERFACE_INIT_
#include "interface.h"

extern float kT, U_MassAccel;
extern PARS_DATABASE< COMM_PARS<50> > cluster_db;

extern char parfname_IConstraint[256];

namespace _interface_constraint_ {

void INTERFACE_CONSTRAINT::construct_efprofile_erf(float U0, float z0, float sigma) {
	float fUnit_kT_mass = ::kT / ::U_MassAccel;

	float dz_max = sigma * 4;
	//this->z0 = z0; this->z1 = z0 - dz_max; this->z2 = z0 + dz_max;
	float z1 = z0 - dz_max, z2 = z0 + dz_max;

	float Uh = U0 / 2;
	double y, x, xabs, cf = 2 * INV_SQRT_PI * fUnit_kT_mass;
	int it, sign;

	float dz = 0.02f;
	int i, n = int(dz_max / dz + 0.5f);
	dz_max = n * dz;
	int npts = 2 * n + 1;
	ef_prof.release();
	ef_prof.set(npts, dz, z1);
	for (i = 0; i < ef_prof.npts; i++) {
		ef_prof.x[i] = z0 + (i - n) * dz ;
		x = dz * (i - n) / sigma; xabs = (x < 0 ? -x : x); sign = (x < 0 ? -1 : 1);
		ERFC(y, xabs, it)
			ef_prof.y[0][i] = (x < 0 ? 0.5 * Uh * y : 0.5 * Uh * (2 - y));
		EXP2(y, xabs, it)
		ef_prof.y[1][i] = -0.5 * Uh * y * cf;
	}
	ef_prof.y[1][0] = 0; ef_prof.y[1][npts-1] = 0;
	/*
	ofstream out;
	out.open("INTERFACE_EF.dat");
	for (i = 0; i < ef_prof.npts; i++) out<<ef_prof.x[i]<<"  "<<ef_prof.y[0][i]<<"  "<<ef_prof.y[1][i]<<endl;
	out.close();
	*/
}

void INTERFACE_CONSTRAINT::ef(float z, double &U, double &f) {
	int n, i;
	double t;

	//if (z <= ef_prof.xf || z >= ef_prof.xt) {U = 0; f = 0;}
	if (z <= ef_prof.xf) {U = ef_prof.y[0][0]; f = ef_prof.y[1][0];}
	else if (z >= ef_prof.xt) {U = ef_prof.y[0][ef_prof.npts-1]; f = ef_prof.y[1][ef_prof.npts-1];}
	else {
		interpolate2(U, f, ef_prof, z, i, t)
	}
	return;
}

bool read_interface_barrier(Barrier *b, int indx) {
	char buffer[256] = "\0", msg[256] = "\0", fname[256] = "\0", cname[256] = "\0";

	float z0 = 0, zf = 0, zt = 0;

	if (strlen(parfname_IConstraint) == 0) return false;

	ifstream in;
	in.open(parfname_IConstraint);
	if (!in.is_open()) {
		sprintf(msg, "failure to open file %s", parfname_IConstraint); show_log(msg, true); return false;
	}
	int model = 0;

	char title[256] = "\0";
	sprintf(title, "[INTERFACE_%d]", indx);
	if (!search(in, title, buffer)) {
		sprintf(msg, "%s is not defined in %s", parfname_IConstraint); show_log(msg, true);
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f", &z0, &zf, &zt) != 3) {
		sprintf(msg, "failure to obtain interface region in %s, see %s", parfname_IConstraint, buffer); show_log(msg, true);
		in.close(); return false;
	}
	b->z0 = z0; 
	if (zf < zt) {b->zf = zf; b->zt = zt;}
	else {b->zf = zt; b->zt = zf;}
	in.close();
	return true;
}

bool read_interface_constraints(IConstraint *ic) {
	char buffer[256] = "\0", msg[256] = "\0", fname[256] = "\0", cname[256] = "\0";

	float fUnit_kT_mass = ::kT / ::U_MassAccel;
	float z0 = 0, zf = 0, zt = 0;

	if (strlen(parfname_IConstraint) == 0) return false;

	ifstream in;
	in.open(parfname_IConstraint);
	if (!in.is_open()) {
		sprintf(msg, "failure to open file %s", parfname_IConstraint); show_log(msg, true); ic->ef_prof.release(); return false;
	}
	int model = 0;

	char title[256] = "\0";
	sprintf(title, "[INTERFACE_%d]", ic->I_indx);
	if (!search(in, title, buffer)) {
		sprintf(msg, "%s is not defined in %s", parfname_IConstraint); show_log(msg, true);
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f", &z0, &zf, &zt) != 3) {
		sprintf(msg, "failure to obtain interface region in %s, see %s", parfname_IConstraint, buffer); show_log(msg, true);
		in.close(); return false;
	}
	ic->z0 = z0;
	ic->zf = (zf < zt ? zf : zt);
	ic->zt = (zf > zt ? zf : zt);

	ARRAY<float*> ef;
	ARRAY<int> col;
	ef.set_array(3); col.set_array(3);
	int npts = 0, ipt;
	float v[4];
	char seperator[5] = " ,;\t";
	char *bf = NULL;
	while (!in.eof()) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", cname) == 1) {
			if (sscanf(buffer, "%s %d", cname, &model) == 2) {
				bf = buffer; NextSection(&bf, seperator); NextSection(&bf, seperator);
				if (bf == NULL) {
					sprintf(msg, "ERROR: unknown format for interface constraint [cluster_name, model[0/1], [0 -- z0, sigma, dG] / [1 -- filename], see : %s", buffer); show_log(msg, true);
					in.close(); ic->ef_prof.release(); return false;
				}
			}
			else if (cname[0] == '[' || cname[0] == '#' || cname[0] == '/') {
				sprintf(msg, "ERROR: %d interface constraint of %s is not defined in %s", ic->I_indx, ic->atom, parfname_IConstraint);
				show_log(msg, true); in.close(); ic->ef_prof.release(); return false;
			}
			else continue;
		}
		else continue;

		if (strcmp(cname, ic->atom) != 0) continue;

		if (model == 0) { 
			strcpy(ic->title, cname);
			if (read_float(bf, v, 3, "%f", seperator) == 3) { // format: z0, sigma, dG
				ic->construct_efprofile_erf(v[2], v[0], v[1]);
			}
			else {
				sprintf(msg, "ERROR: unknown format for interface constraint [cluster_name, model[0], z0, sigma, dG], see : %s", buffer); show_log(msg, true);
				in.close(); ic->ef_prof.release(); return false;
			}
		}
		else if (model == 1) {
			strcpy(ic->title, cname);
			if (sscanf(bf, "%s", fname) == 1 && strlen(fname) > 0) {
				//ic->z0 = z0;
				if (read_profile(fname, 0, ef, col, npts) && npts > 1) {
					ic->ef_prof.set(npts, 0, 0);
					for (ipt = 0; ipt < npts; ipt++) {
						ic->ef_prof.x[ipt] = ef.m[0][ipt];
						ic->ef_prof.y[0][ipt] = ef.m[1][ipt];
						ic->ef_prof.y[1][ipt] = ef.m[2][ipt] * fUnit_kT_mass;
						ic->ef_prof.xf = (ef.m[0][0] < ef.m[0][npts-1] ? ef.m[0][0] : ef.m[0][npts-1]);
						ic->ef_prof.xt = (ef.m[0][0] < ef.m[0][npts-1] ? ef.m[0][npts-1] : ef.m[0][0]);
						ic->ef_prof.dx = ef.m[0][1] - ef.m[0][0];
					}
					delete[] ef.m[0]; delete[] ef.m[1]; delete[] ef.m[2];
					ef.m[0] = NULL; ef.m[1] = NULL; ef.m[2] = NULL;
					in.close(); return true;
				}
				else {
					sprintf(msg, "failure to read interface constraint of %s for interface %d from %s", cname, ic->I_indx, fname); show_log(msg, true);
					in.close(); ic->ef_prof.release(); return false;
				}
			}
			else {
				sprintf(msg, "ERROR: unknown format for interface constraint, see : %s", buffer); show_log(msg, true);
				in.close(); ic->ef_prof.release(); return false;
			}
		}
		else {
			sprintf(msg, "ERROR: disable interface model, %d", model); show_log(msg, true);
			in.close(); ic->ef_prof.release(); return false;
		}

		{
			// shift minimum potential of profile to 0
			int k;
			float ymin = ic->ef_prof.y[0][0];
			if (ymin > ic->ef_prof.y[0][ic->ef_prof.npts - 1]) ymin = ic->ef_prof.y[0][ic->ef_prof.npts - 1];
			if (ymin != 0) for (k = 0; k < ic->ef_prof.npts; k++) ic->ef_prof.y[0][k] -= ymin;
		}

		{
			char ofname[256] = "\0";
			sprintf(ofname, "interface_%d_%s.dat", ic->I_indx, cname);
			ofstream out;
			out.open(ofname);
			for (ipt = 0; ipt < ic->ef_prof.npts; ipt++) {
				out<<ic->ef_prof.x[ipt]<<"  "<<ic->ef_prof.y[0][ipt]
				<<"  "<<ic->ef_prof.y[1][ipt] / fUnit_kT_mass<<endl;
			}
			out.flush(); out.close();
		}
		in.close(); return true;
	}

	in.close();
	return false;
}

bool init_clusterIC_DB(int nInterface, float z1, float z2) {
	clusterIC_db.release();
	int ncs = number< COMM_PARS<50> >(cluster_db.ch);
	if (ncs == 0) return true;
	clusterIC_db.set(ncs, nInterface);
	int i, j;
	COMM_PARS<50> *p = NULL;
	bool status = true;
	char msg[256] = "\0";
	for (i = 0; i < ncs; i++) {
		p = cluster_db.search_par_indx(i);
		for (j = 0; j < nInterface; j++) {
			clusterIC_db.m[i].m[j].uid = p->uid;
			strcpy(clusterIC_db.m[i].m[j].atom, p->atom);
			clusterIC_db.m[i].m[j].indx = p->indx;
			clusterIC_db.m[i].m[j].I_indx = j;
			status = status & read_interface_constraints(clusterIC_db.m[i].m + j);
			if (clusterIC_db.m[i].m[j].z0 >= z2 || clusterIC_db.m[i].m[j].z0 <= z1) {
				sprintf(msg, "error on interface constraint definition : constraint %d for interface %d  has position at %f, out of the range of cell size [%f, %f]", i, j, clusterIC_db.m[i].m[j].z0, z1, z2);
				show_msg(msg); status = false;
			}
		}
	}
	for (i = 0; i < nInterface; i++) {
		if (i >= 2) break;
		read_interface_barrier(::_interface_constraint_::barrier + i, i);
	}
	return status;
}

bool read_interface_image_constant(char *par_fname, float &f1_img, float &z1_img, float &f2_img, float &z2_img) {
	char buffer[256] = "\0", msg[256] = "\0", fname[256] = "\0";
	ifstream in;
	in.open(par_fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open file %s", par_fname); show_log(msg, true); return false;
	}
	if (!search(in, "[CHARGE_IMAGE_CONSTANT]", buffer)) {
		sprintf(msg, "[CHARGE_IMAGE_CONSTANT] is not defined in %s", par_fname);
		show_log(msg, true); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f", &f1_img, &z1_img) != 2) {
		show_log("failure to get charge image constants of the 1st interfaces of slab, see: ", false);
		show_log(buffer, true); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f", &f2_img, &z2_img) != 2) {
		show_log("failure to get charge image constants of the 1st interfaces of slab, see: ", false);
		show_log(buffer, true); in.close(); return false;
	}
	in.close();
	return true;
}

void get_slab_region(FLEX_ASYM_MATRIX<IConstraint> &clusterIC_db, float &z1, float &z2) {
	int ic, i;
	if (clusterIC_db.m == NULL) return;
	float z;
	z1 = 0; z2 = 0;
	for (i = 0; i < clusterIC_db.n; i++) {
		for (ic = 0; ic < clusterIC_db.m[i].n; ic++) {
			z = clusterIC_db.m[i].m[ic].z0;
			if (z1 > z) z1 = z;
			if (z2 < z) z2 = z;
		}
	}
}

void apply_slab_barrier_constraint(FLEX_ASYM_MATRIX<IConstraint> &clusterIC_db) {
	int ic, i, ib, ipt;
	char msg[2560] = "\0";
	if (clusterIC_db.m == NULL) return;
	float z1 = 0, z2 = 0;
	get_slab_region(clusterIC_db, z1, z2);
	IConstraint *pic = NULL;
	float t;
	if (z1 > z2) {SWAP(z1, z2, t)}
	for (i = 0; i < clusterIC_db.n; i++) {
		for (ic = 0; ic < clusterIC_db.m[i].n; ic++) {
			pic = clusterIC_db.m[i].m + ic;

			{
				float fmax = 0;
				ib = 0;
				for (ipt = 0; ipt < pic->ef_prof.npts; ipt++) {
					if (fmax < FABS(pic->ef_prof.y[1][ipt])) {
						ib = ipt;
						fmax = FABS(pic->ef_prof.y[1][ipt]);
					}
				}
			}

			if (pic->ef_prof.xf < z1) { // left side of ef_prof is out of z1
				if (pic->ef_prof.xt > z1) { // across the left edge
					//ib = int((z1 - pic->ef_prof.x[0]) / pic->ef_prof.dx);
					for (ipt = 0; ipt < ib; ipt++) pic->ef_prof.y[1][ipt] = pic->ef_prof.y[1][ib];

					if (pic->ef_prof.xt < z2) { // right side is within z2
					}
					else {
						sprintf(msg, "constraint [%d, %d] is NOT fully inside the slab [%f, %f]. Its region: %f - %f ", i, ic, z1, z2, pic->ef_prof.xf, pic->ef_prof.xt);
						show_log(msg, true);
					}
				}
				else {
					//ib = int((z1 - pic->ef_prof.x[0]) / pic->ef_prof.dx);
					for (ipt = 0; ipt < ib; ipt++) pic->ef_prof.y[1][ipt] = pic->ef_prof.y[1][ib];

					sprintf(msg, "constraint [%d, %d] is out of the slab [%f, %f] at left side. Its region: %f - %f ", i, ic, z1, z2, pic->ef_prof.xf, pic->ef_prof.xt);
					show_log(msg, true);
				}
			}
			else if (pic->ef_prof.xt > z2) { // right side of ef_prof is out of z2
				if (pic->ef_prof.xf < z2) { // across the right edge
					//ib = int((z2 - pic->ef_prof.x[0]) / pic->ef_prof.dx);
					for (ipt = ib + 1; ipt < pic->ef_prof.npts; ipt++) pic->ef_prof.y[1][ipt] = pic->ef_prof.y[1][ib];

					if (pic->ef_prof.xf > z1) { // left side is within z1
					}
					else {
						sprintf(msg, "constraint [%d, %d] is NOT fully inside the slab [%f, %f]. Its region: %f - %f ", i, ic, z1, z2, pic->ef_prof.xf, pic->ef_prof.xt);
						show_log(msg, true);
					}
				}
				else {
					//ib = int((z2 - pic->ef_prof.x[0]) / pic->ef_prof.dx);
					for (ipt = ib + 1; ipt < pic->ef_prof.npts; ipt++) pic->ef_prof.y[1][ipt] = pic->ef_prof.y[1][ib];

					sprintf(msg, "constraint [%d, %d] is out of the right side of the slab [%f, %f]. Its region: %f - %f ", i, ic, z1, z2, pic->ef_prof.xf, pic->ef_prof.xt);
					show_log(msg, true);
				}
			}
			else {
				sprintf(msg, "constraint [%d, %d] is within the slab [%f, %f]. Its region: %f - %f ", i, ic, z1, z2, pic->ef_prof.xf, pic->ef_prof.xt);
				show_log(msg, true);

				if ((pic->z0 - z1) < (z2 - pic->z0)) { // closer to left edge
					//ib = int((z1 - pic->ef_prof.x[0]) / pic->ef_prof.dx);
					for (ipt = 0; ipt < ib; ipt++) pic->ef_prof.y[1][ipt] = pic->ef_prof.y[1][ib];
				}
				else {
					//ib = int((z2 - pic->ef_prof.x[0]) / pic->ef_prof.dx);
					for (ipt = ib + 1; ipt < pic->ef_prof.npts; ipt++) pic->ef_prof.y[1][ipt] = pic->ef_prof.y[1][ib];
				}
			}
			{
				char ofname[256] = "\0";
				sprintf(ofname, "interface_%d_%s.dat", pic->I_indx, pic->title);
				ofstream out;
				out.open(ofname);
				float fUnit_kT_mass = ::kT / ::U_MassAccel;
				for (ipt = 0; ipt < pic->ef_prof.npts; ipt++) {
					out<<pic->ef_prof.x[ipt]<<"  "<<pic->ef_prof.y[0][ipt]
					<<"  "<<pic->ef_prof.y[1][ipt] / fUnit_kT_mass<<endl;
				}
				out.flush(); out.close();
			}
		}
	}
}

} // end of namespace _interface_constraint_ 
