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

#include "MM.h"
#include "Mol.h"
#include "MD.h"

#include "read.h"
#include "var.h"
#include "cluster.h"

extern void pickup_str(char *source, int nstart, int len, char *dest);
bool search(ifstream &in, char *title, char *buffer) {
	char line[512] = "\0";
	char buff[256] = "\0";
	strcpy(buffer, "\0"); 
	while (!in.eof()) {
		in.getline(line, 510);
		sscanf(line, "%s", buff);
		if (strcmp(buff, title) == 0) {strcpy(buffer, line); return true;}
	}
	//if (in.eof()) return false;
	return false;
}

bool search(ifstream &in, char *title, char *buffer, char* stop_line) {
	int nstop = 0;
	if (stop_line != NULL) nstop = strlen(stop_line);

	char line[512] = "\0";
	char buff[256] = "\0";
	char *bf = NULL;
	strcpy(buffer, "\0"); 
	while (!in.eof()) {
		in.getline(line, 510);
		sscanf(line, "%s", buff);
		if (strcmp(buff, title) == 0) {strcpy(buffer, line); return true;}
		else if (stop_line != NULL) {
			bf = line; if (!NextSection(&bf, " \t")) continue; 
			pickup_str(bf, 0, nstop, buff);
			if (strcmp(buff, stop_line) == 0) return false;
		}
	}
	//if (in.eof()) return false;
	return false;
}

bool search(ifstream &in, char *title, int indx, char *buffer) {
	char line[512] = "\0";
	char buff[256] = "\0";
	int n = indx - 1;
	strcpy(buffer, "\0"); 
	while (!in.eof()) {
		in.getline(line, 510);
		if (sscanf(line, "%s %d", buff, &n) == 2 &&
			(strcmp(buff, title) == 0 && n == indx)) {strcpy(buffer, line); return true;}
	}
	//if (in.eof()) return false;
	return false;
}

bool search(ifstream &in, char *title, int indx, char *buffer, char *stop_line) {
	int nstop = 0;
	if (stop_line != NULL) nstop = strlen(stop_line);

	char line[512] = "\0";
	char buff[256] = "\0";
	char *bf = NULL;
	int n = indx - 1;
	strcpy(buffer, "\0"); 
	while (!in.eof()) {
		in.getline(line, 510);
		if (sscanf(line, "%s %d", buff, &n) == 2 &&
			(strcmp(buff, title) == 0 && n == indx)) {strcpy(buffer, line); return true;}
		else if (stop_line != NULL) {
			bf = line; if (!NextSection(&bf, " \t")) continue; 
			pickup_str(bf, 0, nstop, buff);
			if (strcmp(buff, stop_line) == 0) return false;
		}
	}
	//if (in.eof()) return false;
	return false;
}

bool read_item(char *fname, char *title, char *item) {
	char buffer[256] = "\0";
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return false;
	if (!search(in, title, buffer)) {in.close(); return false;}
	in.getline(item, 250);
	in.close();
	return true;
}
/*
bool ReadMoleculeStruct(ifstream &in, MMOLECULE &mm) {
	char buffer[256] = "\0";
	int nc = 0, i;
	float ta = 0, dta = 0;
	float rbase[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	VECTOR3 dr, axis, raxis;
	double domega = 0;

	VECTOR3 ZeroV;

	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f", rbase, rbase + 1, rbase + 2, rbase + 3 , rbase + 4, rbase + 5, rbase + 6, rbase + 7, rbase + 8) != 9) {
		return false;
	}
	for (i = 0; i < 3; i++) {dr.v[i] = (double)rbase[i]; dr.v[i] -= mm.base->Op->v[i];}
	mm.shiftMol(dr);
	for (i = 0; i < 3; i++) axis.v[i] = rbase[i + 3];
	cp(mm.zaxis, axis, raxis);
	domega = angle(mm.zaxis, axis);
	mm.twMol((*(mm.base->Op)), raxis, domega, true);

	for (i = 0; i < 3; i++) axis.v[i] = rbase[i + 6];
	cp(mm.xaxis, axis, raxis);
	domega = angle(mm.xaxis, axis);
	mm.twMol((*(mm.base->Op)), raxis, domega, true);

	//in.getline(buffer, 250); // ignore the speed of base
	//in.getline(buffer, 250); // ignore the acceleration of base

	if (mm.nCluster < MAX_THREADS) {
		while (!in.eof()) {
			in.getline(buffer, 250);
			if (sscanf(buffer, "%f %f %f %f %f %f", rbase, rbase+1, rbase+2, rbase+3, rbase+4, rbase+5) >= 5)
				continue; //this line could be speed and acceleration of base cluster, ignore it
			if (sscanf(buffer, "%d %f", &nc, &ta) != 2) break;
			dta = ta - mm.cluster[nc].km.ta;
			if (dta != 0) {
				mm.twTA(nc, dta, true);
				mm.cluster[nc].km.ta = ta;
			}
			if (nc == mm.nCluster - 1) break;
		}
	}
	else {
	//*******************************************************
	//			rotation all clusters parallelly
	//		NOTE: lines started with // is commented
	//*******************************************************
		// move all clusters 

		mm.setup_base_movement(ZeroV, *(mm.base->Op), ZeroV, true);

		while (!in.eof()) {
			in.getline(buffer, 250);
			if (sscanf(buffer, "%f %f %f %f %f %f", rbase, rbase+1, rbase+2, rbase+3, rbase+4, rbase+5) >= 5)
				continue; //this line could be speed and acceleration of base cluster, ignore it
			if (sscanf(buffer, "%d %f", &nc, &ta) != 2) break;
			if (nc < 0 || nc >= mm.nCluster) break;
			dta = ta - mm.cluster[nc].km.ta;
			mm.cluster[nc].rot.dta = dta;
			if (nc == mm.nCluster - 1) break;
		}

		setup_all_LocalRotMatrix(mm);
		setup_GlobalRotMatrix(mm);
		tweak_all_clusters(mm);

		for (nc = 0; nc < mm.nCluster; nc++) {
			if (mm.free_base && mm.cluster + nc == mm.base) continue;
			mm.cluster[nc].km.ta += mm.cluster[nc].rot.dta;
			//ta = torsion_angle(mm.cluster[nc].parent, mm.cluster + nc); // check torsion angle
		}
	}
	return true;
}
*/

bool ReadMoleculeStruct(ifstream &in, MMOLECULE &mm, bool full = false) {
	char buffer[256] = "\0";
	int nc = 0, i;
	float ta = 0, dta = 0;
	float rbase[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	VECTOR3 dr, zaxis, xaxis, raxis;
	double domega = 0;

	VECTOR3 ZeroV;

	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f", rbase, rbase + 1, rbase + 2, rbase + 3 , rbase + 4, rbase + 5, rbase + 6, rbase + 7, rbase + 8) != 9) {
		return false;
	}
	for (i = 0; i < 3; i++) dr.v[i] = (double)rbase[i]; 
	for (i = 0; i < 3; i++) zaxis.v[i] = rbase[i + 3];
	for (i = 0; i < 3; i++) xaxis.v[i] = rbase[i + 6];

	mm.set_base_move2(dr, zaxis, xaxis);

	in.getline(buffer, 250); // ignore the speed of base
	if (full) {
		if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f", rbase, rbase + 1, rbase + 2, rbase + 3 , rbase + 4, rbase + 5) != 6) {
			show_log("Unknown format for macromolecule's speed", true);
			return false;
		}
		else {
			mm.V0.v[0] = rbase[0]; mm.V0.v[1] = rbase[1]; mm.V0.v[2] = rbase[2]; mm.V0.v[3] = rbase[3]; mm.V0.v[4] = rbase[4]; mm.V0.v[5] = rbase[5];
		}
	}
	in.getline(buffer, 250); // ignore the acceleration of base
	if (full) {
		if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f", rbase, rbase + 1, rbase + 2, rbase + 3 , rbase + 4, rbase + 5) != 6) {
			show_log("Unknown format for macromolecule's acceleration", true);
			return false;
		}
		else {
			mm.alpha0.v[0] = rbase[0]; mm.alpha0.v[1] = rbase[1]; mm.alpha0.v[2] = rbase[2]; mm.alpha0.v[3] = rbase[3]; mm.alpha0.v[4] = rbase[4]; mm.alpha0.v[5] = rbase[5];
		}
	}

	while (!in.eof()) {
		in.getline(buffer, 250);
		//if (sscanf(buffer, "%f %f %f %f %f %f", rbase, rbase+1, rbase+2, rbase+3, rbase+4, rbase+5) >= 5)
		//	continue; //this line could be speed and acceleration of base cluster, ignore it
		if (sscanf(buffer, "%d %f", &nc, &ta) != 2) break;
		if (nc < 0 || nc >= mm.nCluster) break;
		dta = ta - mm.cluster[nc].km.ta;
		mm.cluster[nc].rot.dta = dta;
		if (full) {
			if (sscanf(buffer, "%d %f %f %f", &nc, &ta, rbase, rbase + 1) == 4) {
				mm.cluster[nc].km.tp = rbase[0];
				mm.cluster[nc].km.tpp = rbase[1];
			}
		}
		if (nc == mm.nCluster - 1) break;
	}

	MacroMolMove(&mm);
	return true;
}

bool ReadMoleculeStruct(ifstream &in, MMOLECULE &mm) {
	return ReadMoleculeStruct(in, mm, false);
}

bool ReadMoleculeStructInfo(ifstream &in, MMOLECULE &mm) { // do not rotate the cluster // do not rotate the cluster actually, but reset the torsion angle value
	char buffer[256] = "\0";
	int nc = 0, i;
	float ta = 0, dta = 0;
	float rbase[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	VECTOR3 dr, axis, raxis;
	double domega = 0;

	VECTOR3 ZeroV;

	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f", rbase, rbase + 1, rbase + 2, rbase + 3 , rbase + 4, rbase + 5, rbase + 6, rbase + 7, rbase + 8) != 9) {
		return false;
	}
	for (i = 0; i < 3; i++) {dr.v[i] = (double)rbase[i]; dr.v[i] -= mm.base->Op->v[i];}
	//mm.shiftMol(dr);
	for (i = 0; i < 3; i++) axis.v[i] = rbase[i + 3];
	cp(mm.zaxis, axis, raxis);
	domega = angle(mm.zaxis, axis);
	//mm.twMol((*(mm.base->Op)), raxis, domega, true);

	for (i = 0; i < 3; i++) axis.v[i] = rbase[i + 6];
	cp(mm.xaxis, axis, raxis);
	domega = angle(mm.xaxis, axis);
	//mm.twMol((*(mm.base->Op)), raxis, domega, true);

	//in.getline(buffer, 250); // ignore the speed of base
	//in.getline(buffer, 250); // ignore the acceleration of base

	while (!in.eof()) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f %f %f %f %f %f", rbase, rbase+1, rbase+2, rbase+3, rbase+4, rbase+5) >= 5)
			continue; //this line could be speed and acceleration of base cluster, ignore it
		if (sscanf(buffer, "%d %f", &nc, &ta) != 2) break;
		
		mm.cluster[nc].km.ta = ta;
		/*
		dta = ta - mm.cluster[nc].km.ta;
		if (dta != 0) {
			mm.twTA(nc, dta, true);
			mm.cluster[nc].km.ta = ta;
		}
		*/
		if (nc == mm.nCluster - 1) break;
	}
	return true;
}

bool ConstructMoleculeStruct(MMOLECULE *m, char *fname) {
	if (m == NULL) return false;
	m->reset();
	if (!ConstructSingleChainPolymerFromParamFile(*m, fname, false)) {
		sprintf(errmsg, "failure to construct molecule with the parameter in file %s", fname);
		show_msg(errmsg); return false;
	}

	ifstream in;
	in.open(fname);
	if (!in.is_open()) {sprintf(errmsg, "Can not open file %s", fname); show_msg(errmsg); return false;}
	char buffer[256] = "\0";
	int nCluster = 0, i = 0;
	float ta = 0;
	float v[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	VECTOR3 ds, raxis, axis;
	double domega = 0;
	if (search(in, "[MOL_STRUCT]", buffer)) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f", v, v+1, v+2, v+3, v+4, v+5, v+6, v+7, v+8) == 9) {
			for (i = 0; i < 3; i++) {ds.v[i] = v[i] - m->base->Op->v[i];}
			m->shiftMol(ds);

			for (i = 0; i < 3; i++) axis.v[i] = v[i + 3];
			cp(m->zaxis, axis, raxis);
			domega = angle(m->zaxis, axis);
			m->twMol((*(m->base->Op)), raxis, domega, true);

			for (i = 0; i < 3; i++) axis.v[i] = v[i + 6];
			cp(m->xaxis, axis, raxis);
			domega = angle(m->xaxis, axis);
			m->twMol((*(m->base->Op)), raxis, domega, true);
		}
		while (!in.eof()) {
			in.getline(buffer, 250);
			if (sscanf(buffer, "%f %f %f %f %f %f", v, v+1, v+2, v+3, v+4, v+5) >= 5)
				continue; //this line could be speed and acceleration of base cluster, ignore it
			if (strlen(buffer) > 0 && sscanf(buffer, "%d %f", &nCluster, &ta) == 2) {
				if (nCluster >= 0 && nCluster < m->nCluster) {
					m->twTA(nCluster, ta - m->cluster[nCluster].km.ta, true);
					m->cluster[nCluster].km.ta = ta;
				}
				else break;
			}
			else break;
		}
	}
	in.close();
	return true;
}

bool ReadSolidMoleculeStruct(ifstream &in, SMOLECULE &m, bool full = false) {
	char buffer[256] = "\0";
	int nc = 0, i;
	float ta = 0, dta = 0;
	float rbase[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	VECTOR3 dr, axis, raxis;
	double domega = 0;

	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f", rbase, rbase + 1, rbase + 2, rbase + 3 , rbase + 4, rbase + 5, rbase + 6, rbase + 7, rbase + 8) != 9) {
		return false;
	}
	for (i = 0; i < 3; i++) dr.v[i] = rbase[i] - m.r.v[i];
	m.shiftMol(dr);
	for (i = 0; i < 3; i++) axis.v[i] = rbase[i + 3];
	cp(m.zaxis, axis, raxis);
	domega = angle(m.zaxis, axis);
	rotate_SM(m, raxis, domega, true);

	for (i = 0; i < 3; i++) axis.v[i] = rbase[i + 6];
	cp(m.xaxis, axis, raxis);
	domega = angle(m.xaxis, axis);
	rotate_SM(m, raxis, domega, true);

	if (full) {
		in.getline(buffer, 250); // ignore the speed of base
		if (sscanf(buffer, "%f %f %f %f %f %f", rbase, rbase + 1, rbase + 2, rbase + 3 , rbase + 4, rbase + 5) == 6) {
			m.c->bkm->V0.v[0] = rbase[0]; m.c->bkm->V0.v[1] = rbase[1]; m.c->bkm->V0.v[2] = rbase[2];
			m.c->bkm->V0.v[3] = rbase[3]; m.c->bkm->V0.v[4] = rbase[4]; m.c->bkm->V0.v[5] = rbase[5];
		}
		in.getline(buffer, 250); // ignore the acceleration of base
		if (sscanf(buffer, "%f %f %f %f %f %f", rbase, rbase + 1, rbase + 2, rbase + 3 , rbase + 4, rbase + 5) == 6) {
			m.c->bkm->alpha.v[0] = rbase[0]; m.c->bkm->alpha.v[1] = rbase[1]; m.c->bkm->alpha.v[2] = rbase[2];
			m.c->bkm->alpha.v[3] = rbase[3]; m.c->bkm->alpha.v[4] = rbase[4]; m.c->bkm->alpha.v[5] = rbase[5];
		}
	}

	return true;
}

bool ReadSolidMoleculeStruct(ifstream &in, SMOLECULE &m) {
	return ReadSolidMoleculeStruct(in, m, false);
}

bool ReadPointMoleculeStruct(ifstream &in, PMOLECULE &m, bool full = false) {
	char buffer[256] = "\0";
	float rbase[3] = {0, 0, 0};
	VECTOR3 dr;
	int i;

	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f", rbase, rbase + 1, rbase + 2) != 3) {
		return false;
	}
	for (i = 0; i < 3; i++) dr.v[i] = rbase[i] - m.r.v[i];
	m.shiftMol(dr);
	if (full) {
		in.getline(buffer, 250); // ignore the speed of base
		if (sscanf(buffer, "%f %f %f", rbase, rbase + 1, rbase + 2) == 3) {
			m.v.v[0] = rbase[0]; m.v.v[1] = rbase[1]; m.v.v[2] = rbase[2];
		}
		in.getline(buffer, 250); // ignore the acceleration of base
		if (sscanf(buffer, "%f %f %f", rbase, rbase + 1, rbase + 2) == 3) {
			m.alpha.v[0] = rbase[0]; m.alpha.v[1] = rbase[1]; m.alpha.v[2] = rbase[2];
		}
	}

	return true;
}

bool ReadPointMoleculeStruct(ifstream &in, PMOLECULE &m) {
	return ReadPointMoleculeStruct(in, m, false);
}

#include "cg-mm.h"
#include "cg-cluster.h"
using namespace _coarse_grain_;

bool ConstructCGMMStruct(CG_MMOLECULE *m, char *fname) {
	if (m == NULL) return false;
	m->reset();
	if (!ConstructSingleChainPolymerFromParamFile(*m, fname, false)) {
		sprintf(errmsg, "failure to construct molecule with the parameter in file %s", fname);
		show_msg(errmsg); return false;
	}

	ifstream in;
	in.open(fname);
	if (!in.is_open()) {sprintf(errmsg, "Can not open file %s", fname); show_msg(errmsg); return false;}
	char buffer[256] = "\0";
	int nCluster = 0, ic = 0;
	float v[6] = {0, 0, 0, 0, 0, 0};
	if (search(in, "[MOL_STRUCT]", buffer)) {
		while (!in.eof()) {
			in.getline(buffer, 250);
			if (sscanf(buffer, "%d %f %f %f", &ic, v, v+1, v+2) < 4) 
				continue; //this line is not right
			else {
				if (ic >= m->nCluster) break;
				m->cluster[ic].r->v[0] = v[0]; 
				m->cluster[ic].r->v[1] = v[1];
				m->cluster[ic].r->v[2] = v[2];
				if (bImplicitSolvent) {
					m->cluster[ic].rc_solv.v[0] = v[0];
					m->cluster[ic].rc_solv.v[1] = v[1];
					m->cluster[ic].rc_solv.v[2] = v[2];
				}
			}
		}
	}
	in.close();
	return true;
}

bool ReadCGMMStruct(ifstream &in, CG_MMOLECULE &mm) {
	char buffer[256] = "\0";
	int nCluster = 0, ic = 0;
	float v[3] = {0, 0, 0};

	while (!in.eof()) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %f %f %f", &ic, v, v+1, v+2) < 4) break; //this line is not right
		else {
			if (ic < 0 || ic >= mm.nCluster) break;
			mm.cluster[ic].r->v[0] = v[0]; 
			mm.cluster[ic].r->v[1] = v[1];
			mm.cluster[ic].r->v[2] = v[2];
			if (bImplicitSolvent) {
				mm.cluster[ic].rc_solv.v[0] = v[0];
				mm.cluster[ic].rc_solv.v[1] = v[1];
				mm.cluster[ic].rc_solv.v[2] = v[2];
			}
			if (ic == mm.nCluster - 1) break;
		}
	}
	return true;
}

bool ReadCGSMStruct(ifstream &in, CGSM &sm) {
	char buffer[256] = "\0";
	float v[3] = {0, 0, 0};

	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f", v, v+1, v+2) < 3) return false; //this line is not right
	else {
		sm.r->v[0] = v[0]; 
		sm.r->v[1] = v[1];
		sm.r->v[2] = v[2];
		if (bImplicitSolvent) {
			sm.rc_solv.v[0] = v[0];
			sm.rc_solv.v[1] = v[1];
			sm.rc_solv.v[2] = v[2];
		}
		return true;
	}
}


