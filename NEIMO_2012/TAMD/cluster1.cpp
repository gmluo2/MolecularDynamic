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

extern bool search(ifstream &in, char *title, char *buffer);
extern bool search(ifstream &in, char *title, int indx, char *buffer);
extern bool NextSection(char** original, char ch);

char struct_fname[256] = "\0", force_field_fname[256] = "\0", atom_index_fname[256] = "\0";
char cluster_solvation_energy_fname[256] = "\0";
char monomer_infor_fname[256] = "\0";

ATOM_COMM_PARS_DATABASE atompar_db(&(get_atompar));
PARS_DATABASE< COMM_PARS<50> > cluster_db(NULL);
bool get_atompar(ATOM_COMM_PARS *atompar) {
	ifstream in;
	char buffer[256] = "\0", buff[256] = "\0";

	in.open(force_field_fname);
	if (!in.is_open()) {
		sprintf(errmsg, "can not find force field parameter file %s", force_field_fname); return false;
	}
	if (search(in, atompar->atom, buffer)) {
		if (sscanf(buffer, "%s %f %f", buff, &(atompar->m), &(atompar->alpha)) != 3) {
			sprintf(errmsg, "parameter wrong in amber force field parameter file [atom mass alpha] : %s", buffer);
			in.close(); return false;
		}
	}
	else {
		sprintf(errmsg, "can not find atom : %s", atompar->atom);
		in.close(); return false;
	}

	// Lennard-Jones interaction parameter
	if (!search(in, "MOD4", buffer)) {
		sprintf(errmsg, "can not get Lennard-Jones interaction parameters in %s", force_field_fname); 
		atompar->epsLJ = 0; atompar->rLJ = 1;
		//in.close(); return false;
	}
	if (search(in, atompar->atom, buffer)) {
		sscanf(buffer, "%s %f %f", buff, &(atompar->rLJ), &(atompar->epsLJ));
		atompar->r_min = atompar->rLJ * 0.5; // minimum radius of atom, can not be smaller than rLJ * 0.6
		atompar->rLJ *= 2; //in file is the radius, not distance
		atompar->epsLJ *= 4.176f; // in file the unit is kCal/mol, in this program we want kJ/mol
	}
	else {
		sprintf(errmsg, "can not get Lenard-Jones interaction parameters of %s in %s", atompar->atom, force_field_fname); 
		atompar->epsLJ = 0; atompar->rLJ = 1;
		//in.close(); return false;
	}
	in.close();

	ifstream in1;
	int n = 0;
	in1.open(atom_index_fname);
	if (!in1.is_open()) {
		sprintf(errmsg, "can not find atom index paramete file %s", atom_index_fname); return false;
	}
	if (search(in1, atompar->atom, buffer)) {
		sscanf(buffer, "%s %d", buff, &n);
		atompar->uid = n;
	}
	else {
		sprintf(errmsg, "%s atom index is not defined", atompar->atom); in1.close(); return false;
	}
	in1.close();

	return true;
}

bool get_atom_index(char *atom, int &aindx) {
	char buffer[256] = "\0", buff[250] = "\0";
	ifstream in;
	int n = 0;
	in.open(atom_index_fname);
	if (!in.is_open()) {
		sprintf(errmsg, "can not find atom index paramete file %s", atom_index_fname); return false;
	}
	if (search(in, atom, buffer)) {
		sscanf(buffer, "%s %d", buff, &n);
		aindx = n;
	}
	else {
		sprintf(errmsg, "%s atom index is not defined", atom); in.close(); return false;
	}
	in.close();
	return true;
}

extern bool search(ifstream &in, char *title, char *buffer, char* stop_line);

template <class CLUSTER> bool get_charge(ifstream &in, CLUSTER *pc) {
	// charge
	char buffer[256] = "\0", *buff = NULL;
	file_pos fpos = current_fpos(in);
	if (!search(in, "CHARGE", buffer, "[")) {set_fpos(in, fpos); return false;}
	int i = 0;
	float c = 0, c_norm = 1 / sqrt(eps_dc);
	while (i < pc->nAtoms) {
		in.getline(buffer, 250);
		buff = buffer;
		if (sscanf(buffer, "%f", &c) != 1) {
			sprintf(errmsg, "charge of %dth atom is not defined", i); 
			set_fpos(in, fpos); return false;
		}
		pc->atom[i].c0 = c; // the charge
		pc->atom[i].c = c * c_norm; // normalize the charge value by / sqrt(eps_dc)
		pc->atom[i].eNeutral = (c == 0 ? true : false);
		i++;
		while (i < pc->nAtoms && NextSection(&buff) && buff != NULL && sscanf(buff, "%f", &c) == 1) {
			pc->atom[i].c0 = c; // the charge
			pc->atom[i].c = c * c_norm;  // normalize the charge value by / sqrt(eps_dc)
			pc->atom[i].eNeutral = (c == 0 ? true : false);
			i++;
		}
	}
	set_fpos(in, fpos);
	return true;
}

template <class CLUSTER> bool get_dipole(ifstream &in, CLUSTER *pc) {
	// dipole
	char buffer[256] = "\0", *buff = NULL;
	file_pos fpos = current_fpos(in);
	if (!search(in, "DIPOLE", buffer, "[")) {set_fpos(in, fpos); return false;}
	int i = 0;
	float mu[3] = {0, 0, 0}, c_norm = 1 / sqrt(eps_dc);
	bool eneutral = false;
	while (i < pc->nAtoms) {
		in.getline(buffer, 250);
		buff = buffer;
		if (sscanf(buffer, "%f %f %f", &(mu[0]), &(mu[1]), &(mu[2])) != 3) {
			sprintf(errmsg, "dipole of %dth atom is not defined", i); 
			set_fpos(in, fpos); return false;
		}
		eneutral = true;
		pc->atom[i].intr_mu.v[0] = mu[0]; // dipole
		pc->atom[i].intr_mu.v[1] = mu[1]; // dipole
		pc->atom[i].intr_mu.v[2] = mu[2]; // dipole
		//pc->atom[i].c = c * c_norm; // normalize the charge value by / sqrt(eps_dc)
		if (mu[0] != 0 || mu[1] != 0 || mu[2] != 0) {
			if (!::bDipole) bDipole = true;
			eneutral = false;
		}
		if (!eneutral) {
			if (pc->atom[i].mMP < 1) pc->atom[i].mMP = 1;
			pc->atom[i].bIMu = true;
		}
		if (pc->atom[i].eNeutral && !eneutral) pc->atom[i].eNeutral = false;
		i++;
	}
	set_fpos(in, fpos);
	return true;
}


bool get_index(ifstream &in, short &index) {
	char buffer[256] = "\0";
	int indx = 0;
	file_pos fpos = current_fpos(in);
	if (!search(in, "INDEX", buffer, "[")) {
		sprintf(errmsg, "can not get cluster index"); set_fpos(in, fpos); return false;
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d", &indx) != 1) {
			sprintf(errmsg, "can not get cluster UNIT INDEX"); set_fpos(in, fpos); return false;
		}
		else index = indx;
	}
	set_fpos(in, fpos);
	return true;
}

bool defined_simple_molecule(char *unit) {
	char buffer[256] = "\0", title[256] = "\0";
	int nflag = 0;

	ifstream in;
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", struct_fname); return false;}
	sprintf(title, "[%s]", unit);
	while (1) {
		if (in.eof()) {sprintf(errmsg, "can not find cluster %s", unit); in.close(); return false;}
		if (!search(in, title, buffer)) {sprintf(errmsg, "can not find cluster %s", unit); in.close(); return false;}
		if (sscanf(buffer, "%s %d", title, &nflag) == 2 && nflag == 0) break; // simple cluster
	}
	in.close(); return true;
}

bool defined_point_molecule(char *unit) {
	char buffer[256] = "\0", title[256] = "\0";
	int nflag = 0;

	ifstream in;
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", struct_fname); return false;}
	sprintf(title, "[%s]", unit);
	while (1) {
		if (in.eof()) {sprintf(errmsg, "can not find cluster %s", unit); in.close(); return false;}
		if (!search(in, title, buffer)) {sprintf(errmsg, "can not find cluster %s", unit); in.close(); return false;}
		if (sscanf(buffer, "%s %d", title, &nflag) == 2 && nflag == 0) break; // simple cluster
	}
	in.close(); return true;
}

bool defined_cluster(char *unit) {
	char buffer[256] = "\0", title[256] = "\0";
	int nflag = 0;

	ifstream in;
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", struct_fname); return false;}
	sprintf(title, "[%s]", unit);
	while (1) {
		if (in.eof()) {sprintf(errmsg, "can not find cluster %s", unit); in.close(); return false;}
		if (!search(in, title, buffer)) {sprintf(errmsg, "can not find cluster %s", unit); in.close(); return false;}
		if (sscanf(buffer, "%s %d", title, &nflag) == 2 && nflag == 0) break; // simple cluster
	}
	in.close(); return true;
}

bool defined_monomer(char *unit) {
	char buffer[256] = "\0", title[256] = "\0";
	int nflag = 0;

	ifstream in;
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", struct_fname); return false;}
	sprintf(title, "[%s]", unit);
	while (1) {
		if (in.eof()) {sprintf(errmsg, "can not find cluster %s", unit); in.close(); return false;}
		if (!search(in, title, buffer)) {sprintf(errmsg, "can not find cluster %s", unit); in.close(); return false;}
		if (sscanf(buffer, "%s %d", title, &nflag) == 2 && nflag == 1) break; // simple cluster
	}
	in.close(); return true;
}

extern bool NextSection(char **original, char ch);

/**************************************************
	Z-Matrix definition:
	[atom0] //this position is [0, 0, 0]
	[atom1] [bound length] [bound atom indx] [bound valence -- single, double, triple]
	[atom2] [bound length] [bound atom indx] [bound angle] [angle atom indx] [bound valence -- single, double, triple]
	[atom3] [bound length] [bound atom indx] [bound angle] [angle atom indx] [torsion angle] [torsion atom indx] [bound valence -- single, double, triple]
	[atom4] [bound length] [bound atom indx] [bound angle] [angle atom indx] [torsion angle] [torsion atom indx] [bound valence -- single, double, triple]
***************************************************/
bool read_cluster_zmatrix(char *cname, BASIC_CLUSTER *pc) {
class ZMATRIX_ATOM {
public:
	char atom[5]; // the atom name for force field parameter
	char infor[250];
	VECTOR3 r;
	int ba, bound_valence, bound_type;

	ZMATRIX_ATOM() {strcpy(atom, "\0"); strcpy(infor, "\0"); ba = 0; bound_type = NORMAL_BOUND;};
	bool get_infor(char *buffer) {
		char *bf = buffer;
		if (sscanf(buffer, "%s", atom) != 1) return false;
		if (NextSection(&bf, ' ')) strcpy(infor, bf);
		else strcpy(infor, "\0");
		return true;
	};
};

class BOUND {
public:
	int na1, na2, bound_valence, bound_type; //the two atoms
};

	bool status = true;
	char title[256] = "\0", buffer[256] = "\0", *buff = NULL;
	int natoms = 0, nHinges = 0, nbounds = 0, i, j, k;
	int nflag = 0;

	pc->reset();

	ifstream in;
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", struct_fname); show_msg(errmsg); return false;}
	sprintf(title, "[%s]", cname);
	while (1) {
		if (in.eof()) {sprintf(errmsg, "can not find cluster %s", cname); show_msg(errmsg); in.close(); return false;}
		if (!search(in, title, buffer)) {sprintf(errmsg, "can not find cluster %s", cname); show_msg(errmsg); in.close(); return false;}
		if (sscanf(buffer, "%s %d", title, &nflag) == 2 && nflag == 0) break; // simple cluster
		else {sprintf(errmsg, "UNKNOWN format : %s\n has to be %s 0", buffer, title); show_msg(errmsg); in.close(); return false;}
	}

	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d %d", &natoms, &nHinges, &nbounds) != 3) {
		sprintf(errmsg, "can not get numbers of atom, additional bounds & hinges from %s after [%s]", buffer, cname); show_msg(errmsg); in.close(); return false;
	}
	if (natoms <= 0) {sprintf(errmsg, "atom number can not be < 1"); show_msg(errmsg); in.close(); return false;}
	if (nHinges < 0) nHinges = 0;
	if (nbounds < 0) nbounds = 0;

	sprintf(title, "ZMATRIX");
	if (!search(in, title, buffer)) {sprintf(errmsg, "can not find %s of %s", title, cname); show_msg(errmsg); in.close(); return false;}

	pc->reset(); pc->setAtoms(natoms); pc->set_hinges(nHinges);

	ZMATRIX_ATOM *zatom = new ZMATRIX_ATOM[natoms + nHinges];
	int ba = 0, aa = 0, ta = 0, nRead = 0;
	float blen = 0, bangl = 0, tors_angle = 0;
	ATOM_COMM_PARS *atompar = NULL;
	VECTOR3 vect;

	for (i = 0; i < natoms + nHinges; i++) {
		in.getline(buffer, 250);
		// get cordinates
		if (!zatom[i].get_infor(buffer)) {
			sprintf(errmsg, "failure to get %d-th zmatrix atom from %s", i, buffer); show_msg(errmsg); 
			delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
		}
		if (i == 0) {}
		else if (i == 1) {
			if ((nRead = sscanf(zatom[i].infor, "%f %d %d %d", &blen, &ba, &(zatom[i].bound_valence), &(zatom[i].bound_type))) >= 3) {
				vect = VECTOR3(0, 0, blen); ZMatrix(zatom[0].r, vect, zatom[i].r);
				if (nRead < 4) zatom[i].bound_type = NORMAL_BOUND;
			}
			else {
				sprintf(errmsg, "failure to get 2nd zmatrix atom from %s", buffer); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
			}
		}
		else if (i == 2) {
			if ((nRead = sscanf(zatom[i].infor, "%f %d %f %d %d %d", &blen, &ba, &bangl, &aa, &(zatom[i].bound_valence), &(zatom[i].bound_type))) >= 5) {
				if (aa >= i || ba >= i) {
					sprintf(errmsg, "bounding atom and angle atom have to be defined for atom %d", i); show_msg(errmsg); 
					delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
				}
				vect = VECTOR3(0, 1, 0); ZMatrix(zatom[ba].r, zatom[aa].r, vect, blen, bangl, zatom[i].r);
				if (nRead < 6) zatom[i].bound_type = NORMAL_BOUND;
			}
			else {
				sprintf(errmsg, "failure to get 3rd zmatrix atom from %s", buffer); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
			}
		}
		else {
			if ((nRead = sscanf(zatom[i].infor, "%f %d %f %d %f %d %d %d", &blen, &ba, &bangl, &aa, &tors_angle, &ta, &(zatom[i].bound_valence), &(zatom[i].bound_type))) >= 7) {
				if (ba >= i || aa >= i || ta >= i) {
					sprintf(errmsg, "bounding atom, angle atom and torsion atom have to be defined for atom %d", i); show_msg(errmsg); 
					delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
				}
				ZMatrix(zatom[ba].r, zatom[aa].r, zatom[ta].r, blen, bangl, tors_angle, zatom[i].r);
				if (nRead < 8) zatom[i].bound_type = NORMAL_BOUND;
			}
			else {
				sprintf(errmsg, "failure to get %d-th zmatrix atom from %s", i, buffer); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
			}
		}

		if (i != 0) zatom[i].ba = ba;

		if (i < natoms) { //this is atom
			//add atom to the database and the cluster
			atompar = atompar_db.attach(zatom[i].atom);
			if (atompar == NULL) {
				show_msg(errmsg); 
				delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
			}
			V32V3(zatom[i].r, pc->atom[i].r)
			pc->atom[i].par = atompar;
			pc->atom[i].aindx = atompar->indx;
		}
		else { // this is hinge
			// do nothing here
		}
	}

	// additional bounds
	BOUND *bound = NULL;
	if (nbounds > 0) {
		strcpy(title, "BOUND");
		if (!search(in, title, buffer)) {
			show_msg(errmsg); delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
		}
		bound = new BOUND[nbounds];
		for (i = 0; i < nbounds; i++) {
			in.getline(buffer, 250);
			if ((nRead = sscanf(buffer, "%d %d %d %d", &(bound[i].na1), &(bound[i].na2), &(bound[i].bound_valence), &(bound[i].bound_type))) < 3) {
				sprintf(errmsg, "can not get additional bound %d", i); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
			else if (nRead < 4) bound[i].bound_type = NORMAL_BOUND;
		}
	}

	// setup bounds & HINGEs for each atoms in cluster
	int nb = 0;
	if (pc->nAtoms > 1 || nHinges > 0) {
		for (i = 0; i < pc->nAtoms; i++) {
			nb = 0;
			for (j = 1; j < natoms + nHinges; j++) { // ZMATRIX_ATOM LINE 0 has NO bound atom defined
				if (j == i || zatom[j].ba == i) nb++;
			}
			for (j = 0; j < nbounds; j++) {
				if (bound[j].na1 == i || bound[j].na2 == i) nb++;
			}
			pc->atom[i].set_bounds(nb);
		}
		for (i = 1; i < natoms; i++) {
			if (!bound_connect<MATOM>(pc->atom + i, pc->atom + zatom[i].ba, zatom[i].bound_valence, zatom[i].bound_type)) {
				sprintf(errmsg, "failure to setup bound between atoms %d & %d in %s", i, zatom[i].ba, cname); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
		}
		for (i = 0; i < nbounds; i++) {
			if (!bound_connect<MATOM>(pc->atom + bound[i].na1, pc->atom + bound[i].na2, bound[i].bound_valence, bound[i].bound_type)) {
				sprintf(errmsg, "failure to setup bound between atoms %d & %d in %s", bound[i].na1, bound[i].na2, cname); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
		}

		// setup hinges
		for (i = 0; i < nHinges; i++) {
			j = i + natoms;
			nb = pc->atom[zatom[j].ba].get_free_bound();
			if (nb < 0) {
				sprintf(errmsg, "failure to setup hinge %d on atom %d in %s -- no free bound", i, zatom[j].ba, cname); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
			// how many hinges was on this bounding atom of i hinge ?
			// there could be more than 1 hinge on the bounding atom
			for (k = natoms; k < j; k++) {
				if (zatom[k].ba == zatom[j].ba) nb++;
			}
			if (nb >= pc->atom[zatom[j].ba].nb) {
				sprintf(errmsg, "failure to setup hinge %d on atom %d in %s -- no free bound", i, zatom[j].ba, cname); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
			pc->setup_hinge(i, zatom[j].ba, nb);
			pc->hinge[i].vHinge = zatom[j].r - zatom[zatom[j].ba].r;
		}
	}

	// clean memory -- zatoms and bounds here
	delete[] zatom; zatom = NULL;
	if (bound != NULL) delete[] bound; bound = NULL;

	//size
	if (!search(in, "SIZE", buffer)) {
		sprintf(errmsg, "can not get cluster size of %s", cname); show_msg(errmsg); in.close(); pc->reset(); return false;
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f", &(pc->d)) != 1) {
			sprintf(errmsg, "can not get cluster size of %s", cname); show_msg(errmsg); in.close(); pc->reset(); return false;
		}
	}

	// charge & cluster index
	COMM_PARS<50> *cpar = NULL;
	if (!get_charge(in, pc) || !get_index(in, pc->cID)) {
		in.close(); pc->reset(); return false;
	}
	else {
		get_dipole(in, pc);
		cpar = cluster_db.attach(pc->cID);
		if (cpar != NULL) pc->cTypeIndx = cpar->indx;
		if (cpar != NULL && strlen(cpar->atom) == 0) strcpy(cpar->atom, cname); 
	}
	in.close();

	pc->psolv = ESolv_db.attach(pc->cID);
	//pc->cal_geometry_radius();
	return true;
}

bool read_cluster_xyz(char *cname, BASIC_CLUSTER *pc) {
class XYZ_ATOM {
public:
	char atom[5]; // the atom name for force field parameter
	char infor[250];
	float r[3];

	XYZ_ATOM() {strcpy(atom, "\0"); strcpy(infor, "\0"); r[0] = 0; r[1] = 0; r[2] = 0;};
	bool get_infor(char *buffer) {
		char *bf = buffer;
		if (sscanf(buffer, "%s", atom) != 1) return false;
		if (NextSection(&bf, ' ')) strcpy(infor, bf);
		else strcpy(infor, "\0");
		if (sscanf(infor, "%f %f %f", r, r+1, r+2) == 3) return true;
		return false;
	};
	int inlist(int natoms, char *atom_name) {
		int i;
		for (i = 0; i < natoms; i++) {
			if (strcmp((this + i)->atom, atom_name) == 0) return i;
		}
		return -1;
	};
};


class BOUND {
public:
	int na1, na2, bound_valence, bound_type; //the two atoms
};

	bool status = true;
	char title[256] = "\0", buffer[256] = "\0", *buff = NULL;
	int natoms = 0, nHinges = 0, nbounds = 0, i, j, k;
	int nflag = 0;

	pc->reset();

	ifstream in;
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", struct_fname); show_msg(errmsg); return false;}
	sprintf(title, "[%s]", cname);
	while (1) {
		if (in.eof()) {sprintf(errmsg, "can not find cluster %s", cname); show_msg(errmsg); in.close(); return false;}
		if (!search(in, title, buffer)) {sprintf(errmsg, "can not find cluster %s", cname); show_msg(errmsg); in.close(); return false;}
		if (sscanf(buffer, "%s %d", title, &nflag) == 2 && nflag == 0) break; // simple cluster
		else {sprintf(errmsg, "UNKNOWN format : %s\n has to be %s 0", buffer, title); show_msg(errmsg); in.close(); return false;}
	}

	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d %d", &natoms, &nHinges, &nbounds) != 3) {
		sprintf(errmsg, "can not get numbers of atom, additional bounds & hinges from %s after [%s]", buffer, cname); show_msg(errmsg); in.close(); return false;
	}
	if (natoms <= 0) {sprintf(errmsg, "atom number can not be < 1"); show_msg(errmsg); in.close(); return false;}
	if (nHinges < 0) nHinges = 0;
	if (nbounds < 0) nbounds = 0;

	sprintf(title, "XYZ");
	if (!search(in, title, buffer)) {sprintf(errmsg, "can not find %s of %s", title, cname); show_msg(errmsg); in.close(); return false;}

	pc->reset(); pc->setAtoms(natoms); pc->set_hinges(nHinges);

	XYZ_ATOM *atom = new XYZ_ATOM[natoms + nHinges];
	int ba = 0, aa = 0, ta = 0, nRead = 0;
	float blen = 0, bangl = 0, tors_angle = 0;
	ATOM_COMM_PARS *atompar = NULL;
	VECTOR3 vect;

	for (i = 0; i < natoms + nHinges; i++) {
		in.getline(buffer, 250);
		// get cordinates
		if (!atom[i].get_infor(buffer)) {
			sprintf(errmsg, "failure to get %d-th atom from %s", i, buffer); show_msg(errmsg); 
			delete[] atom; atom = NULL; in.close(); pc->reset(); return false;
		}

		if (i < natoms) { //this is atom
			//add atom to the database and the cluster
			atompar = atompar_db.attach(atom[i].atom);
			if (atompar == NULL) {
				show_msg(errmsg); 
				delete[] atom; atom = NULL; in.close(); pc->reset(); return false;
			}
			pc->atom[i].r.v[0] = atom[i].r[0]; pc->atom[i].r.v[1] = atom[i].r[1]; pc->atom[i].r.v[2] = atom[i].r[2];
			pc->atom[i].par = atompar;
			pc->atom[i].aindx = atompar->indx;
		}
		else { // this is hinge
			// do nothing here
		}
	}

	// bounds -- all bounds are defined here 
	BOUND *bound = NULL;
	if (nbounds > 0) {
		strcpy(title, "BOUND");
		if (!search(in, title, buffer)) {
			show_msg(errmsg); delete[] atom; atom = NULL; in.close(); pc->reset(); return false;
		}
		bound = new BOUND[nbounds];
		for (i = 0; i < nbounds; i++) {
			in.getline(buffer, 250);
			if ((nRead = sscanf(buffer, "%d %d %d %d", &(bound[i].na1), &(bound[i].na2), &(bound[i].bound_valence), &(bound[i].bound_type))) < 3) {
				sprintf(errmsg, "can not get additional bound %d", i); show_msg(errmsg); 
				delete[] atom; atom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
			else if (nRead < 4) bound[i].bound_type = NORMAL_BOUND;
		}
	}

	// setup bounds & HINGEs for each atoms in cluster
	int nb = 0, ib = 0, ih = 0;
	if (pc->nAtoms > 1 || nHinges > 0) {
		for (i = 0; i < pc->nAtoms; i++) {
			nb = 0;
			for (j = 0; j < nbounds; j++) {
				if (bound[j].na1 == i || bound[j].na2 == i) nb++;
			}
			pc->atom[i].set_bounds(nb);
		}
		for (i = 0; i < nbounds; i++) { // non-hinge bound
			if (atom[natoms].inlist(nHinges, atom[bound[i].na1].atom) >= 0 || atom[natoms].inlist(nHinges, atom[bound[i].na2].atom) >= 0) continue; // hinge atom
			if (!bound_connect<MATOM>(pc->atom + bound[i].na1, pc->atom + bound[i].na2, bound[i].bound_valence, bound[i].bound_type)) {
				sprintf(errmsg, "failure to setup bound between atoms %d & %d in %s", bound[i].na1, bound[i].na2, cname); show_msg(errmsg); 
				delete[] atom; atom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
		}

		// setup hinges
		ih = 0;
		for (i = 0; i < nbounds; i++) {
			// ib is for bounding atom (real atom), j is the hinge atom
			if (atom[natoms].inlist(nHinges, atom[bound[i].na1].atom) >= 0) { // na1 is the hinge atom
				ib = bound[i].na2; j = bound[i].na1;
			}
			else if (atom[natoms].inlist(nHinges, atom[bound[i].na2].atom) >= 0) { // na2 is the hinge atom
				ib = bound[i].na1; j = bound[i].na2;
			}
			else continue; // this bound is not a hinge

			nb = pc->atom[ib].get_free_bound();
			if (nb < 0) {
				sprintf(errmsg, "failure to setup hinge %d (bound %d) on atom %d in %s -- no free bound", ih, i, ib, cname); show_msg(errmsg); 
				delete[] atom; atom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
			// how many hinges was on this bounding atom of i hinge ?
			// there could be more than 1 hinge on the bounding atom
			for (k = 0; k < i; k++) {
				if (atom[natoms].inlist(nHinges, atom[bound[k].na1].atom) >= 0 && bound[k].na2 == ib) nb++;
				else if (atom[natoms].inlist(nHinges, atom[bound[k].na2].atom) >= 0 && bound[k].na1 == ib) nb++;
			}
			if (nb >= pc->atom[ib].nb) {
				sprintf(errmsg, "failure to setup hinge %d (bound %d) on atom %d in %s -- no free bound", ih, i, ib, cname); show_msg(errmsg); 
				delete[] atom; atom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
			pc->setup_hinge(ih, ib, nb);
			pc->hinge[ih].vHinge.v[0] = atom[ib].r[0] - atom[j].r[0];
			pc->hinge[ih].vHinge.v[1] = atom[ib].r[1] - atom[j].r[1];
			pc->hinge[ih].vHinge.v[2] = atom[ib].r[2] - atom[j].r[2];

			ih++;
		}
	}

	// clean memory -- zatoms and bounds here
	delete[] atom; atom = NULL;
	if (bound != NULL) delete[] bound; bound = NULL;

	//size
	if (!search(in, "SIZE", buffer)) {
		sprintf(errmsg, "can not get cluster size of %s", cname); show_msg(errmsg); in.close(); pc->reset(); return false;
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f", &(pc->d)) != 1) {
			sprintf(errmsg, "can not get cluster size of %s", cname); show_msg(errmsg); in.close(); pc->reset(); return false;
		}
	}

	// charge & cluster index
	COMM_PARS<50> *cpar = NULL;
	if (!get_charge(in, pc) || !get_index(in, pc->cID)) {
		in.close(); pc->reset(); return false;
	}
	else {
		get_dipole(in, pc);
		cpar = cluster_db.attach(pc->cID);
		if (cpar != NULL) pc->cTypeIndx = cpar->indx;
		if (cpar != NULL && strlen(cpar->atom) == 0) strcpy(cpar->atom, cname);
	}
	in.close();

	pc->psolv = ESolv_db.attach(pc->cID);
	//pc->cal_geometry_radius();
	return true;
}

bool read_cluster(char *cname, BASIC_CLUSTER *pc) {
	char title[256] = "\0", buffer[256] = "\0";
	int nflag = 0;
	ifstream in;
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", struct_fname); show_msg(errmsg); return false;}
	sprintf(title, "[%s]", cname);
	while (1) {
		if (in.eof()) {sprintf(errmsg, "can not find cluster %s", cname); show_msg(errmsg); in.close(); return false;}
		if (!search(in, title, buffer)) {sprintf(errmsg, "can not find cluster %s", cname); show_msg(errmsg); in.close(); return false;}
		if (sscanf(buffer, "%s %d", title, &nflag) == 2 && nflag == 0) break; // simple cluster
		else {sprintf(errmsg, "UNKNOWN format : %s\n has to be %s 0", buffer, title); show_msg(errmsg); in.close(); return false;}
	}

	in.getline(buffer, 250);// ignore one line
	char format[250] = "\0";
	while (1) {
		if (in.eof()) {
			sprintf(errmsg, "structure of cluster %s is not defined after %s in file %s", cname, title, struct_fname); in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", format) == 0) continue;
		else break;
	}
	in.close();
	if (strcmp(format, "ZMATRIX") == 0) return read_cluster_zmatrix(cname, pc);
	else if (strcmp(format, "XYZ") == 0) return read_cluster_xyz(cname, pc);
	else {
		sprintf(errmsg, "in file %s, 2nd line after %s is not the atom-format definition, ZMATRIX or XYZ", struct_fname, title); show_msg(errmsg);
		return false;
	}
}

bool read_monomer(char *monomer, MONOMER *m) {
	char buffer[256] = "\0", buff[256] = "\0", title[100] = "\0", cname[250] = "\0", *bf = NULL;
	int nUnits = 0, nRelShips = 0, i, nflag = 0;
	ifstream *in_cluster = NULL;
	BASIC_CLUSTER *pc = NULL, *mc = NULL;

	COMM_PARS<50> *cpar = NULL;
	
	ifstream in;
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", struct_fname); show_msg(errmsg); return false;}
	sprintf(title, "[%s]", monomer);
	strcpy(m->name, monomer);
	while (1) {
		if (in.eof()) {sprintf(errmsg, "can not find cluster %s", monomer); show_msg(errmsg); in.close(); return false;}
		if (!search(in, title, buffer)) {sprintf(errmsg, "can not find cluster %s", monomer); show_msg(errmsg); in.close(); return false;}
		if (sscanf(buffer, "%s %d", title, &nflag) == 2 && nflag == 1) break; // monomer
		else {sprintf(errmsg, "UNKNOWN format : %s\n has to be %s 1", buffer, title); show_msg(errmsg); in.close(); return false;}
	}

	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d", &nUnits, &nRelShips) != 2) {
		sprintf(errmsg, "can not get number of units & relation-ships in the monomer from %s", buffer); show_msg(errmsg); in.close(); return false;
	}
	if (nUnits <= 0 || nRelShips < 0) {
		sprintf(errmsg, "units can not be less than 1, relationships can not be less than 0"); show_msg(errmsg); in.close(); return false;
	}
	m->reset(nUnits, nRelShips);

	for (i = 0; i < nUnits; i++) {
		sprintf(title, "#%d", i);
		if (!search(in, title, buffer)) {
			sprintf(errmsg, "can not find #%d cluster", i); m->reset(0, 0); show_msg(errmsg); in.close(); return false;
		}
		if (sscanf(buffer, "%s %s %d", title, cname, &nflag) != 3) {
			sprintf(errmsg, "UNKNOWN format of cluster, [# cluster  flag] : %s", buffer); show_msg(errmsg); in.close(); m->reset(0, 0); return false;
		}
		pc = m->cluster + i;
		if (!read_cluster(cname, pc)) { show_msg(errmsg); in.close(); m->reset(0, 0); return false;}
		if (nflag == 1) { // reset charge and cluster index
			if (!get_charge(in, pc) || !get_index(in, pc->cID)) { show_msg(errmsg); in.close(); m->reset(0, 0); return false;}
			else {
				get_dipole(in, pc);
				cpar = cluster_db.attach(pc->cID);
				if (cpar != NULL) pc->cTypeIndx = cpar->indx;
				if (cpar != NULL && strlen(cpar->atom) == 0) strcpy(cpar->atom, cname);
			}
			// if index is changed, solvation will change also
			pc->psolv = ESolv_db.attach(pc->cID);
		}
	}

	strcpy(title, "CONNECT");
	int nc1 = 0, nc2 = 0, iHinge1 = 0, iHinge2 = 0;
	float torsion_angle = 0;
	for (i = 0; i < nRelShips; i++) {
		if (!search(in, title, buffer)) {
			sprintf(errmsg, "can not find connection #%d", i); show_msg(errmsg); in.close(); m->reset(0, 0); return false;
		}
		bf = buffer;
		NextSection(&bf, ' ');
		if (sscanf(bf, "%d %d %d %d", &nc1, &nc2, &iHinge1, &iHinge2) == 4) {
			if (!ConnectClusters(m->cluster + nc1, iHinge1, m->cluster + nc2, iHinge2, SINGLE_BOUND)) {
				sprintf(errmsg, "failure to connect cluster #%d to #%d", nc1, nc2); show_msg(errmsg); in.close();
				m->reset(0, 0); return false;
			}
			torsion_angle = IdeaTorsionAngleFromSymetry((BASIC_CLUSTER*)(m->cluster + nc1), (BASIC_CLUSTER*)(m->cluster + nc2));
			series_rotate(m->cluster + nc2, iHinge2, torsion_angle - m->cluster[nc2].km.ta, true);
			m->rship[i].set_relationship(nc1, iHinge1, nc2, iHinge2);
		}
		else {
			sprintf(errmsg, "UNKNOWN connection format [CONNECT cluster1 cluster2 hinge1, hinge2] : %s", buffer); show_msg(errmsg); in.close();
			m->reset(0, 0); return false;
		}
	}

	strcpy(title, "FREE-HINGES");
	int nFreeHinges = 0;
	if (!search(in, title, buffer)) {
		sprintf(errmsg, "can not find free hinge for monomer %s", m->name); show_msg(errmsg); in.close(); m->reset(0, 0); return false;
	}
	if (sscanf(buffer, "%s %d", title, &nFreeHinges) != 2) {
		sprintf(errmsg, "UNKNOWN FORMAT : %s", buffer); show_msg(errmsg); m->reset(0, 0); in.close(); return false;
	}
	m->setHinges(nFreeHinges);
	for (i = 0; i < nFreeHinges; i++) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d", &nc1, &iHinge1) == 2) {
			if (m->cluster[nc1].hinge[iHinge1].plink != NULL) {
				sprintf(errmsg, "hinge %d of unit #%d is not free", iHinge1, nc1); show_msg(errmsg); m->reset(0, 0); in.close(); return false;
			}
			m->nc[i] = nc1; m->iHinge[i] = iHinge1;
		}
		else {
			sprintf(errmsg, "UNKNOWN FORMAT for free hinge definition: %s", buffer); show_msg(errmsg); in.close(); m->reset(0, 0); return false;
		}
	}

	if (nFreeHinges <= 1) {
		strcpy(title, "BASE");
		if (!search(in, title, buffer) || sscanf(buffer, "%s %d", title, &nc1) != 2) {
			sprintf(errmsg, "can not get base cluster from %s", buffer); m->reset(0, 0); show_msg(errmsg); in.close(); return false;
		}
		m->nbase = nc1;
	}

	int uid = 0;
	strcpy(title, "UID");
	if (!search(in, title, buffer) || sscanf(buffer, "%s %d", title, &uid) != 2) {
		sprintf(errmsg, "can not get UID of this monomer from %s", buffer); m->reset(0, 0); show_msg(errmsg); in.close(); return false;
	}
	m->uid = uid;

	return true;
}

bool construct_defined_simple_molecule(char *mol, BASIC_CLUSTER* p) {
	if (p == NULL) return false;
	if (!read_cluster(mol, p)) return false;
	return true;
}

bool construct_defined_point_molecule(char *mol, BASIC_CLUSTER* p) {
	if (p == NULL) return false;
	if (!read_cluster(mol, p)) return false;
	if (p->nAtoms > 1) {
		sprintf(errmsg, "the defined cluster %s has more than 1 atom. This is not a point-molecule!", mol);
		show_msg(errmsg); p->reset(); return false;
	}
	return true;
}

bool construct_defined_cluster(char *unit, CLUSTER* p) {
	if (p == NULL) return false;
	if (!read_cluster(unit, p)) return false;
	return true;
}

int construct_defined_monomer(char *monomer, MONOMER &m) {
	if (!read_monomer(monomer, &m)) return false;
	else return true;
}

bool get_implicit_solvation_pars(CLUSTER_IMPSOLV_PARS *par) {
	par->sigma = 0; par->f_sigma = 0;
	ifstream in;
	in.open(cluster_solvation_energy_fname);
	if (!in.is_open()) {
		sprintf(errmsg, "can not find cluster's solvation parameters parameter file %s", cluster_solvation_energy_fname);
		mlog.show(errmsg);
		return true;
	}
	float e = 0, r = 0;
	int indx = 0;
	char buffer[256] = "\0", title[250] = "\0";
	if (!search(in, "CLUSTER", par->cID, buffer)) {
		sprintf(errmsg, "cluster solvation parameters of cluster [ID: %d] is not defined in %s", par->cID, cluster_solvation_energy_fname);
		mlog.show(errmsg);
		in.close();
		return true;
	}
	in.close();
	if (sscanf(buffer, "%s %d %f %f", title, &indx, &e, &r) != 4) {
		sprintf(errmsg, "cluster solvation parameters of cluster [ID: %d] is not defined in %s", par->cID, cluster_solvation_energy_fname);
		mlog.show(errmsg);
		return true;
	}
	else {
		par->sigma = e; 
		par->f_sigma = e * ::kT / U_InertiaMoment;
		par->r0 = r;
	}
	return true;
}

bool readSolventSize(float &rSolvent) {
	rSolvent = 1;
	ifstream in;
	in.open(cluster_solvation_energy_fname);
	if (!in.is_open()) {
		sprintf(errmsg, "can not find cluster's implicit solvation parameter file %s", cluster_solvation_energy_fname);
		mlog.show(errmsg);
		return true;
	}

	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	float v1 = 0, v2 = 0;

	strcpy(item, "SOLVENT");
	if (!search(in, item, buffer)) {
		sprintf(errmsg, "can not find solvent radius from parameter file %s", cluster_solvation_energy_fname);
		mlog.show(errmsg);
		in.close(); return true;
	}
	if (sscanf(buffer, "%s %f", item, &v1) != 2) {
		sprintf(errmsg, "can not find solvent radius from parameter file %s", cluster_solvation_energy_fname);
		mlog.show(errmsg);
		in.close(); return true;
	}
	else {
		rSolvent = v1; //exclusion distance
	}
	in.close();
	return true;
}

bool readImplicitSolvationPars(CLUSTER_IMPSOLV_PARS_DATABASE &sdb, float &rSolvent) {
	bool status = readSolventSize(rSolvent);
	if (!status) return true;
	CHAIN<CLUSTER_IMPSOLV_PARS> *pch = sdb.ch;
	CLUSTER_IMPSOLV_PARS *par = NULL;
	while (pch != NULL) {
		par = pch->p;
		status = status & get_implicit_solvation_pars(par);
		pch = pch->next;
	}
	return status;
}

void calImplicitSolvationArea(CLUSTER_IMPSOLV_PARS_DATABASE &sdb, float rSolvent) {
	CHAIN<CLUSTER_IMPSOLV_PARS> *pch = sdb.ch;
	CLUSTER_IMPSOLV_PARS *par = NULL;
	float r0w = 0;
	while (pch != NULL) {
		par = pch->p;
		r0w = par->r0 + rSolvent;
		par->S0 = 4 * PI * r0w * r0w;
		pch = pch->next;
	}
}

bool readDeuterizedMonomer(char *fname, char *monomer) {
	char buffer[256] = "\0", name[250] = "\0";
	bool status = false;
	int nindx1 = 0, nindx2 = 0, nitems = 0;

	if (strlen(fname) == 0) return false;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return false; // file does not exist, so, monomer is not deuterized
	if (!search(in, "[MONOMER-DEUTERIZE]", buffer)) {in.close(); return false;}
	status = false;
	while (!status && !in.eof()) {
		in.getline(buffer, 250);
		nitems = sscanf(buffer, "%s %d %d", name, &nindx1, &nindx2);
		if (nitems == 0) {status = false; break;} // end of definition, important
		if (strcmp(monomer, name) == 0) {
			switch(nitems) {
			case 1: // only monomer is given, all this kind of monomer is deuterized
				status = true; break;
			default:
				status = true; break;
			}
		}
		if (status) break;
	}
	in.close();
	return status;
}

void monomer_deuterization(MONOMER *m) {
	bool status = readDeuterizedMonomer(monomer_infor_fname, m->name);
	int n = 0;
	for (n = 0; n < m->nCluster; n++) m->cluster->deuterized = status;
	return;
}




/*********************************************************************************
**********                    READ _BASIC_CLUSTER                      **********
*********************************************************************************/

/**************************************************
	Z-Matrix definition:
	[atom0] //this position is [0, 0, 0]
	[atom1] [bound length] [bound atom indx] [bound valence -- single, double, triple]
	[atom2] [bound length] [bound atom indx] [bound angle] [angle atom indx] [bound valence -- single, double, triple]
	[atom3] [bound length] [bound atom indx] [bound angle] [angle atom indx] [torsion angle] [torsion atom indx] [bound valence -- single, double, triple]
	[atom4] [bound length] [bound atom indx] [bound angle] [angle atom indx] [torsion angle] [torsion atom indx] [bound valence -- single, double, triple]
***************************************************/
template <class CLUSTER_ATOM> bool read_basic_cluster(char *cname, _BASIC_CLUSTER<CLUSTER_ATOM> *pc) {
class ZMATRIX_ATOM {
public:
	char atom[5]; // the atom name for force field parameter
	char infor[250];
	VECTOR3 r;
	int ba, bound_valence, bound_type;

	ZMATRIX_ATOM() {strcpy(atom, "\0"); strcpy(infor, "\0"); ba = 0; bound_type = NORMAL_BOUND;};
	bool get_infor(char *buffer) {
		char *bf = buffer;
		if (sscanf(buffer, "%s", atom) != 1) return false;
		if (NextSection(&bf, ' ')) strcpy(infor, bf);
		else strcpy(infor, "\0");
		return true;
	};
};

class BOUND {
public:
	int na1, na2, bound_valence, bound_type; //the two atoms
};

	bool status = true;
	char title[256] = "\0", buffer[256] = "\0", *buff = NULL;
	int natoms = 0, nHinges = 0, nbounds = 0, i, j, k;
	int nflag = 0;

	pc->reset();

	ifstream in;
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", struct_fname); show_msg(errmsg); return false;}
	sprintf(title, "[%s]", cname);
	while (1) {
		if (in.eof()) {sprintf(errmsg, "can not find cluster %s", cname); show_msg(errmsg); in.close(); return false;}
		if (!search(in, title, buffer)) {sprintf(errmsg, "can not find cluster %s", cname); show_msg(errmsg); in.close(); return false;}
		if (sscanf(buffer, "%s %d", title, &nflag) == 2 && nflag == 0) break; // simple cluster
		else {sprintf(errmsg, "UNKNOWN format : %s\n has to be %s 0", buffer, title); show_msg(errmsg); in.close(); return false;}
	}

	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d %d", &natoms, &nHinges, &nbounds) != 3) {
		sprintf(errmsg, "can not get numbers of atom, additional bounds & hinges from %s after [%s]", buffer, cname); show_msg(errmsg); in.close(); return false;
	}
	if (natoms <= 0) {sprintf(errmsg, "atom number can not be < 1"); show_msg(errmsg); in.close(); return false;}
	if (nHinges < 0) nHinges = 0;
	if (nbounds < 0) nbounds = 0;

	sprintf(title, "ZMATRIX");
	if (!search(in, title, buffer)) {sprintf(errmsg, "can not find %s of %s", title, cname); show_msg(errmsg); in.close(); return false;}

	pc->reset(); pc->setAtoms(natoms); pc->set_hinges(nHinges);

	ZMATRIX_ATOM *zatom = new ZMATRIX_ATOM[natoms + nHinges];
	int ba = 0, aa = 0, ta = 0, nRead = 0;
	float blen = 0, bangl = 0, tors_angle = 0;
	VECTOR3 vect;

	for (i = 0; i < natoms + nHinges; i++) {
		in.getline(buffer, 250);
		// get cordinates
		if (!zatom[i].get_infor(buffer)) {
			sprintf(errmsg, "failure to get %d-th zmatrix atom from %s", i, buffer); show_msg(errmsg); 
			delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
		}
		if (i == 0) {}
		else if (i == 1) {
			if ((nRead = sscanf(zatom[i].infor, "%f %d %d %d", &blen, &ba, &(zatom[i].bound_valence), &(zatom[i].bound_type))) >= 3) {
				vect = VECTOR3(0, 0, blen); ZMatrix(zatom[0].r, vect, zatom[i].r);
				if (nRead < 4) zatom[i].bound_type = NORMAL_BOUND;
			}
			else {
				sprintf(errmsg, "failure to get 2nd zmatrix atom from %s", buffer); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
			}
		}
		else if (i == 2) {
			if ((nRead = sscanf(zatom[i].infor, "%f %d %f %d %d %d", &blen, &ba, &bangl, &aa, &(zatom[i].bound_valence), &(zatom[i].bound_type))) >= 5) {
				if (aa >= i || ba >= i) {
					sprintf(errmsg, "bounding atom and angle atom have to be defined for atom %d", i); show_msg(errmsg); 
					delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
				}
				vect = VECTOR3(0, 1, 0); ZMatrix(zatom[ba].r, zatom[aa].r, vect, blen, bangl, zatom[i].r);
				if (nRead < 6) zatom[i].bound_type = NORMAL_BOUND;
			}
			else {
				sprintf(errmsg, "failure to get 3rd zmatrix atom from %s", buffer); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
			}
		}
		else {
			if ((nRead = sscanf(zatom[i].infor, "%f %d %f %d %f %d %d %d", &blen, &ba, &bangl, &aa, &tors_angle, &ta, &(zatom[i].bound_valence), &(zatom[i].bound_type))) >= 7) {
				if (ba >= i || aa >= i || ta >= i) {
					sprintf(errmsg, "bounding atom, angle atom and torsion atom have to be defined for atom %d", i); show_msg(errmsg); 
					delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
				}
				ZMatrix(zatom[ba].r, zatom[aa].r, zatom[ta].r, blen, bangl, tors_angle, zatom[i].r);
				if (nRead < 8) zatom[i].bound_type = NORMAL_BOUND;
			}
			else {
				sprintf(errmsg, "failure to get %dth zmatrix atom from %s", i, buffer); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
			}
		}

		if (i != 0) zatom[i].ba = ba;

		if (i < natoms) { //this is atom
			V32V3(zatom[i].r, pc->atom[i].r)
			//pc->atom[i].aindx = atompar->aindx;
		}
		else { // this is hinge
			// do nothing here
		}
	}

	// additional bounds
	BOUND *bound = NULL;
	if (nbounds > 0) {
		strcpy(title, "BOUND");
		if (!search(in, title, buffer)) {
			show_msg(errmsg); delete[] zatom; zatom = NULL; in.close(); pc->reset(); return false;
		}
		bound = new BOUND[nbounds];
		for (i = 0; i < nbounds; i++) {
			in.getline(buffer, 250);
			if ((nRead = sscanf(buffer, "%d %d %d %d", &(bound[i].na1), &(bound[i].na2), &(bound[i].bound_valence), &(bound[i].bound_type))) < 3) {
				sprintf(errmsg, "can not get additional bound %d", i); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
			else if (nRead < 4) bound[i].bound_type = NORMAL_BOUND;
		}
	}

	// setup bounds & HINGEs for each atoms in cluster
	int nb = 0;
	if (pc->nAtoms > 1 || nHinges > 0) {
		for (i = 0; i < pc->nAtoms; i++) {
			nb = 0;
			for (j = 1; j < natoms + nHinges; j++) { // ZMATRIX_ATOM LINE 0 has NO bound atom defined
				if (j == i || zatom[j].ba == i) nb++;
			}
			for (j = 0; j < nbounds; j++) {
				if (bound[j].na1 == i || bound[j].na2 == i) nb++;
			}
			pc->atom[i].set_bounds(nb);
		}
		for (i = 1; i < natoms; i++) {
			if (!bound_connect<CLUSTER_ATOM>(pc->atom + i, pc->atom + zatom[i].ba, zatom[i].bound_valence, zatom[i].bound_type)) {
				sprintf(errmsg, "failure to setup bound between atoms %d & %d in %s", i, zatom[i].ba, cname); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
		}
		for (i = 0; i < nbounds; i++) {
			if (!bound_connect<CLUSTER_ATOM>(pc->atom + bound[i].na1, pc->atom + bound[i].na2, bound[i].bound_valence, bound[i].bound_type)) {
				sprintf(errmsg, "failure to setup bound between atoms %d & %d in %s", bound[i].na1, bound[i].na2, cname); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
		}

		// setup hinges
		for (i = 0; i < nHinges; i++) {
			j = i + natoms;
			nb = pc->atom[zatom[j].ba].get_free_bound();
			if (nb < 0) {
				sprintf(errmsg, "failure to setup hinge %d on atom %d in %s -- no free bound", i, zatom[j].ba, cname); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
			// how many hinges was on this bounding atom of i hinge ?
			// there could be more than 1 hinge on the bounding atom
			for (k = natoms; k < j; k++) {
				if (zatom[k].ba == zatom[j].ba) nb++;
			}
			if (nb >= pc->atom[zatom[j].ba].nb) {
				sprintf(errmsg, "failure to setup hinge %d on atom %d in %s -- no free bound", i, zatom[j].ba, cname); show_msg(errmsg); 
				delete[] zatom; zatom = NULL; if (bound != NULL) delete[] bound; bound = NULL; in.close(); pc->reset(); return false;
			}
			pc->setup_hinge(i, zatom[j].ba, nb);
			pc->hinge[i].vHinge = zatom[j].r - zatom[zatom[j].ba].r;
		}
	}

	// clean memory -- zatoms and bounds here
	delete[] zatom; zatom = NULL;
	if (bound != NULL) delete[] bound; bound = NULL;

	// cluster index
	COMM_PARS<50> *cpar = NULL;
	short uid = 0;
	if (get_index(in, uid)) {
		cpar = cluster_db.attach(uid);
		if (cpar != NULL) pc->cTypeIndx = cpar->indx;
		if (cpar != NULL && strlen(cpar->atom) == 0) strcpy(cpar->atom, cname);
	}

	in.close();
	return true;
}

/*********************************************************************************
**********   READ COARSE-GRAINED CLUSTER, WITH ONE ATOMS INSIDE ONLY,   **********
**********   HINGE / BOUND ARE NOT CONFINED IN ANGLE AND LENGTH ******************
*********************************************************************************/
#include "cg-mm.h"

namespace _coarse_grain_ {
CGATOM_COMM_PARS_DATABASE cgatompar_db(&(get_atompar));
char cgatom_par_fname[256] = "\0", cgatom_index_fname[256] = "\0";

bool get_cgatomname_from_alias(char *fname, char *alias, char *cgatom) {
	char msg[256] = "\0", buffer[256] = "\0";
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open file, %s, looking for cgatom name for %s", fname, alias);
		show_msg(msg); return false;
	}
	char alias_name[256] = "\0", atom[256] = "\0";
	if (!search(in, "ALIAS", buffer)) {in.close(); return false;}
	while (!in.eof()) {
		in.getline(msg, 250);
		if (sscanf(msg, "%s %s", alias_name, atom) == 2) {
			if (strcmp(alias_name, alias) == 0) {
				strcpy(cgatom, atom); in.close(); return true;
			}
		}
		else break; // the format is not for alias anymore, stop here, eg. space line
	}
	in.close(); return false;
}

bool get_cgatom_indx(char *cgatom, int &indx) {
	char msg[256] = "\0", buffer[256] = "\0", buff[256] = "\0";
	ifstream in;
	in.open(cgatom_index_fname);
	if (!in.is_open()) {
		sprintf(errmsg, "can not find atom index paramete file %s", cgatom_index_fname); 
		show_msg(errmsg); return false;
	}
	if (search(in, cgatom, buffer)) {
		sscanf(buffer, "%s %d", buff, &indx);
		in.close(); return true;
	}
	else {
		sprintf(errmsg, "%s atom index is not defined", cgatom); show_msg(errmsg);
		in.close(); return false;
	}
}
bool get_atompar(CGATOM_COMM_PARS *cgatom_par) {
	char msg[256] = "\0", buffer[256] = "\0", buff[256] = "\0";
	ifstream in;
	in.open(cgatom_par_fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open cgatom par. file, %s, looking for cgatom's information of %s", cgatom_par_fname, cgatom_par->atom);
		show_msg(msg); return false;
	}
	if (!search(in, "CGATOM", buffer)) {
		sprintf(msg, "CGATOM is not defined in %s", cgatom_par_fname); show_msg(msg); 
		in.close(); return false;
	}
	if (search(in, cgatom_par->atom, buffer)) {
		if (sscanf(buffer, "%s %f %f", buff, &(cgatom_par->m), &(cgatom_par->alpha)) != 3) {
			sprintf(msg, "failure to get atom %s's infor. from %s", cgatom_par->atom, buffer); show_msg(msg); 
			in.close(); return false;
		}
	}
	else {
		sprintf(msg, "can not get cgatom %s's infor. in %s", cgatom_par->atom, cgatom_par_fname); show_msg(msg); 
		in.close(); return false;
	}

	// Lennard-Jones interaction parameter
	if (!search(in, "MOD4", buffer)) {
		sprintf(errmsg, "can not get Lennard-Jones interaction parameters in %s", cgatom_par_fname); 
		show_msg(errmsg); 
		cgatom_par->epsLJ = 0; cgatom_par->rLJ = 1;
		//in.close(); return false;
	}
	if (search(in, cgatom_par->atom, buffer)) {
		sscanf(buffer, "%s %f %f", buff, &(cgatom_par->rLJ), &(cgatom_par->epsLJ));
		cgatom_par->rLJ *= 2; //in file is the radius, not distance
		cgatom_par->epsLJ *= 4.176f; // in file the unit is kCal/mol, in this program we want kJ/mol
		cgatom_par->rcut = cgatom_par->rLJ * 4; // rcut will be used for LocalRelationship estimation
	}
	else {
		sprintf(errmsg, "can not get Lenard-Jones interaction parameters of %s in %s", cgatom_par->atom, cgatom_par_fname); 
		show_msg(errmsg); 
		cgatom_par->epsLJ = 0; cgatom_par->rLJ = 1; cgatom_par->rcut = cgatom_par->rLJ * 4;
		//in.close(); return false;
	}

	in.close();

	int n = 0;
	if (!get_cgatom_indx(cgatom_par->atom, n)) return false;
	cgatom_par->uid = n;
	return true;
}

bool read_cgcluster(char *cname, CG_CLUSTER *pc) {
	bool status = true;
	char title[256] = "\0", buffer[256] = "\0", *buff = NULL;
	int nflag = 0;
	sprintf(title, "[%s]", cname);

	pc->reset();

	if (!read_basic_cluster<CGATOM>(cname, (_BASIC_CLUSTER<CGATOM>*)pc)) return false;
	pc->setup_cordinate();
/*
	ifstream in;
	in.open(struct_fname);
	search(in, title, buffer);
	//size
	if (!search(in, "SIZE", buffer)) {
		sprintf(errmsg, "can not get cluster size of %s", cname); show_msg(errmsg); in.close(); pc->reset(); return false;
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f", &(pc->d)) != 1) {
			sprintf(errmsg, "can not get cluster size of %s", cname); show_msg(errmsg); in.close(); pc->reset(); return false;
		}
	}

	// charge & cluster index
	if (!get_charge<CG_CLUSTER>(in, pc) || !get_index<CG_CLUSTER>(in, pc)) {
		in.close(); return false;
	}
	else {
		get_dipole(in, pc);
	}

	pc->psolv = ESolv_db.attach(pc->cID);
	//pc->cal_geometry_radius();
	in.close();
	*/

	return true;
}

bool read_cgcluster_infor(char *fname, CG_CLUSTER *pc) {
	bool status = true;
	char title[256] = "\0", buffer[256] = "\0", *buff = NULL;
	char cname[10] = "\0"; strcpy(cname, pc->atom[0].par->atom);
	int nflag = 0;
	//sprintf(title, "[%s]", cname);

	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(buffer, "can not open file %s to read coarse-grained cluster information", fname);
		show_msg(buffer); return false;
	}
	/*
	if (!search(in, title, buffer)) {
		sprintf(buffer, "coarse-grained cluster information of %s is not given in %s", cname, fname);
		show_msg(buffer); in.close(); return false;
	}
	*/

	// to be finished

	//size
	if (!search(in, "BROWNIAN", buffer) || !search(in, cname, buffer)) {
		sprintf(errmsg, "can not get Brownian cluster-radius of %s", cname); show_msg(errmsg); in.close(); return false;
	}
	else {
		//in.getline(buffer, 250);
		if (sscanf(buffer, "%s %f", title, &(pc->d)) != 2) {
			sprintf(errmsg, "can not get Brownian cluster-radius of %s", cname); show_msg(errmsg); in.close(); return false;
		}
	}

	in.close();
	return true;
}

bool construct_defined_cgcluster(char *unit, int nb, CG_CLUSTER* p) {
	if (p == NULL) return false;
	char uunit[50] = "\0";
	switch (nb) {
	case 0: // no bound
		strcpy(uunit, "CG0"); break;
	case 1: // 1 bound, eg. head or tail
		strcpy(uunit, "CG1"); break;
	case 2: // 2 bounds, eg. head or tail
		strcpy(uunit, "CG2"); break;
	case 3: // 3 bounds, eg. head or tail
		strcpy(uunit, "CG3"); break;
	case 4: // 4 bounds, eg. head or tail
		strcpy(uunit, "CG4"); break;
	default: // not defined in cluster.dat
		return false;
	}

	if (!read_cgcluster(uunit, p)) return false;
	if (!get_cgatomname_from_alias(cgatom_par_fname, unit, uunit)) strcpy(uunit, unit);
	//strcpy(p->atom->par->atom, uunit);
	CGATOM_COMM_PARS *cgatompar = cgatompar_db.attach(uunit);
	if (cgatompar == NULL) {p->reset(); return false;}
	p->atom[0].par = cgatompar;
	p->atom[0].aindx = p->atom[0].par->indx;
	p->cID = p->atom[0].par->uid; // we use the cgatom indx as the cluster index

	// cluster index
	COMM_PARS<50> *cpar = NULL;
	cpar = cluster_db.attach(p->cID);
	if (cpar != NULL) p->cTypeIndx = cpar->indx;
	if (cpar != NULL && strlen(cpar->atom) == 0) strcpy(cpar->atom, uunit);

	p->psolv = ESolv_db.attach(p->cID);

	read_cgcluster_infor(cgatom_par_fname, p); // Brownian radius
	return true;
}

bool construct_defined_cgsm(char *unit, CGSM* m) {
	if (m == NULL) return false;
	CG_CLUSTER *pc = (CG_CLUSTER*)m;
	char uunit[50] = "\0";
	strcpy(uunit, "CG0"); // CGSM has no bound

	if (!read_cgcluster(uunit, pc)) return false;
	if (!get_cgatomname_from_alias(cgatom_par_fname, unit, uunit)) strcpy(uunit, unit);
	//strcpy(p->atom->par->atom, uunit);
	CGATOM_COMM_PARS *cgatompar = cgatompar_db.attach(uunit);
	if (cgatompar == NULL) {pc->reset(); return false;}
	pc->atom[0].par = cgatompar;
	pc->atom[0].aindx = pc->atom[0].par->indx;
	pc->cID = pc->atom[0].par->uid; // we use the cgatom indx as the cluster index

	COMM_PARS<50> *cpar = NULL;
	cpar = cluster_db.attach(pc->cID);
	if (cpar != NULL) pc->cTypeIndx = cpar->indx;
	if (cpar != NULL && strlen(cpar->atom) == 0) strcpy(cpar->atom, uunit);

	read_cgcluster_infor(cgatom_par_fname, pc);
	return true;
}

bool ConstructCGSM(CGSM* cgsm, char *mol) {
	return construct_defined_cgsm(mol, cgsm);
}

bool get_cgatom_uid(char *cgatom, int &indx) {
	char atom[50] = "\0";
	if (!get_cgatomname_from_alias(cgatom_par_fname, cgatom, atom)) strcpy(atom, cgatom);
	return get_cgatom_indx(atom, indx);
}

} // end of namespace _coarse_grain_
