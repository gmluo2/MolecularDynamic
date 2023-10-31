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
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "MD.h"
#include "CMM.h"
#include "var.h"

#include "Interaction.h"
#include "Interact1.h"

#include "md_cell.h"
#include "ReadMDProc.h"
#include "cluster.h"
#include "read.h"
#include "mol-cell.h"

#include "cg-mm.h"
#include "cg-md.h"
#include "cg-md-cell.h"

extern MMOL_MD_CELL mm_md_cell;
extern LIST<MMOLECULE> *mlist;

/*
MMOLECULE* copy_mol_cell(MMOLECULE* m, int iMol) {
	if (iMol >= mm_md_cell.mm.n) return NULL;
	MMOLECULE *m2 = mm_md_cell.mm.m + iMol;
	m2->reset();
	cp(m2, m, true); // molecule will be copied, the tree will be setup
	return m2;
}
*/




/*
bool ReadDefinedSimpleMolecule(ifstream &in, SMOLECULE& sm, CHAIN<VECTOR3> **pos, int& nMol) {
	sm.reset();
	if (pos != NULL) release_chain<VECTOR3>(pos, true);
	nMol = 0;

	char buffer[256] = "\0", item[200] = "\0", mol_template_fname[250] = "\0";
	int iMol = 0;
	float dx = 0, dy = 0, dz = 0;
	VECTOR3 *ds = NULL;

	strcpy(item, "[SIMPLE-MOLECULE]");
	if (!search(in, item, buffer)) {
		sprintf(errmsg, "Can not find %s", item);
		return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%s %d", mol_template_fname, &nMol) != 2) {
		sprintf(errmsg, "can not get molecule template filename and the number of the molecule from %s", buffer);
		return false;
	}
	if (nMol <= 0) return true;

	strcpy(item, "[SHIFT]");
	if (!search(in, item, buffer)) {sprintf(errmsg, "Can not fine %s", item); return false;}
	iMol = 0;
	while(iMol < nMol) {
		if (in.eof()) {
			sprintf(errmsg, "get %d of molecule-relative-position only, less than %d", iMol, nMol);
			release_chain<VECTOR3>(pos, true);
			return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f %f %f", &dx, &dy, &dz) != 3) {
			sprintf(errmsg, "can not get molecule-relative-position from : %s", buffer);
			release_chain<VECTOR3>(pos, true);
			return false;
		}
		ds = new VECTOR3(dx, dy, dz);
		if (iMol == 0) {*pos = new CHAIN<VECTOR3>; (*pos)->p = ds;}
		else (*pos)->attach2tail(ds);
		iMol++;
	}
	if (!ConstructSimpleMolecule(sm, mol_template_fname)) {
		sprintf(buffer, "Can not construct molecule with parameter file %s", mol_template_fname);
		strcat(errmsg, buffer);
		release_chain<VECTOR3>(pos, true);
		sm.reset();
		return false;
	}
	return true;
}
*/

// minumn distance to the atoms in MMOLECULE
/*
bool Overlap(SMOLECULE &sample, VECTOR3& r0, float radius, ARRAY<MMOLECULE> &mm, ARRAY<CHAIN<VECTOR3>*>& mm_pos, int nMM) {
	double dis2 = 0, mdis2 = radius * radius + 1e-4, mdis2_min = radius * radius;
	if (sample.c == NULL || sample.c->nAtoms == 0 || mm.m == NULL || nMM <= 0) return false;
	MMOLECULE *pm = NULL;
	CLUSTER *pc = NULL;
	CHAIN<VECTOR3>* pos = NULL;
	VECTOR3 dr, *ds_m = NULL;
	int nm = 0, nc = 0, na = 0, i = 0, k = 0;
	bool ignore = false;
	for (nm = 0; nm < nMM; nm++) {
		pm = mm.m + nm;
		pos = mm_pos.m[nm];
		while (pos != NULL) {
			ds_m = pos->p;
			if (ds_m == NULL) break;
			for (nc = 0; nc < pm->nCluster; nc++) {
				pc = pm->cluster + nc;
				for (na = 0; na < pc->nAtoms; na++) {
					for (k = 0; k < sample.c->nAtoms; k++) {
						dis2 = 0; ignore = false;
						for (i = 0; i < 3; i++) {
							dr.v[i] = pc->atom[na].r.v[i] + ds_m->v[i] - sample.c->atom[k].r.v[i] - r0.v[i];
							if (FABS(dr.v[i]) > radius) {ignore = true; break;}
							dis2 += dr.v[i] * dr.v[i];
						}
						if (ignore) continue;
						else mdis2 = (dis2 < mdis2 ? dis2 : mdis2);
						if (mdis2 < mdis2_min) return true;
					}
				}
			}
			pos = pos->next;
		}
	}
	return false;
}

// minumn distance to the atoms in MMOLECULE
bool Overlap(SMOLECULE &sample, VECTOR3& r0, float radius, ARRAY<SMOLECULE> &sm, ARRAY<CHAIN<VECTOR3>*>& sm_pos, int nSM) {
	double dis2 = 0, mdis2 = radius * radius + 1e-4, mdis2_min = radius * radius;
	if (sm.m == NULL || nSM <= 0) return false;
	BASIC_CLUSTER *pc = NULL;
	SMOLECULE *pm = NULL;
	CHAIN<VECTOR3>* pos = NULL;

	VECTOR3 dr, *ds_m = NULL;
	int nm = 0, na = 0, i = 0, k = 0;
	bool ignore = false;
	for (nm = 0; nm < nSM; nm++) {
		pm = sm.m + nm;
		pos = sm_pos.m[nm];
		pc = pm->c;
		while (pos != NULL) {
			ds_m = pos->p;
			if (ds_m == NULL) break;
			for (na = 0; na < pc->nAtoms; na++) {
				for (k = 0; k < sample.c->nAtoms; k++) {
					dis2 = 0; ignore = false;
					for (i = 0; i < 3; i++) {
						dr.v[i] = pc->atom[na].r.v[i] + ds_m->v[i] - sample.c->atom[k].r.v[i] - r0.v[i];
						if (FABS(dr.v[i]) > radius) {ignore = true; break;}
						dis2 += dr.v[i] * dr.v[i];
					}
					if (ignore) continue;
					else mdis2 = (dis2 < mdis2 ? dis2 : mdis2);
					if (mdis2 < mdis2_min) return true;
				}
			}
			pos = pos->next;
		}
	}
	return false;
}
*/

bool GenSolventMolPos(VECTOR3& u, int nx, int ny, int nz, CHAIN<VECTOR3> **pos_chain, int& nSolv) {
	CHAIN<VECTOR3> *mpos = NULL;
	if (*pos_chain != NULL) release_chain<VECTOR3>(pos_chain, true);

	VECTOR3 *v = NULL;
	int i, j, k;
	nSolv = 0;
	float vc[3] = {-float(nx - 1) / 2, -float(ny - 1) / 2, -float(nz - 1) / 2};
	for (i = 0; i < nx; i++) {for (j = 0; j < ny; j++) {for (k = 0; k < nz; k++) {
		v = new VECTOR3((vc[0] + i) * u.v[0], (vc[1] + j) * u.v[1], (vc[2] + k) * u.v[2]);
		v->v[0] += (2 * ranf() - 1) * 0.5 * u.v[0];
		v->v[1] += (2 * ranf() - 1) * 0.5 * u.v[1];
		v->v[2] += (2 * ranf() - 1) * 0.5 * u.v[2];
		if (mpos == NULL) {mpos = new CHAIN<VECTOR3>; mpos->p = v;}
		else mpos->attach2tail(v);
		nSolv++;
	}}}
	*pos_chain = mpos;
	return true;
};

bool ReadDefinedCellSize(char *fname, float &x, float &y, float &z) {
	char title[50] = "\0", buffer[256] = "\0", buff[50] = "\0";
	int ndim = 0, mol_type;
	float ux, uy, uz;

	ifstream in;

	in.open(fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", fname); show_msg(errmsg); return false;}

	strcpy(title, "[SOLVENT]");
	search(in, title, buffer);
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %s %d", &mol_type, buff, &ndim) != 3) {
		sprintf(errmsg, "can not get number of solvent molecule in each axis from %s", buffer);
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f", &ux, &uy, &uz) != 3) {
		sprintf(errmsg, "Can not get [a, b, c] of cubic cell from %s", buffer);
		in.close(); return false;
	}
	x = float(ndim) * ux;
	y = float(ndim) * uy;
	z = float(ndim) * uz;
	in.close(); return true;
}

