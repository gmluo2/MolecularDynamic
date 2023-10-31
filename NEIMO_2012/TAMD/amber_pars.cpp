#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>

using namespace std;
#include "def.h"
#include "dstruct.h"
#include "amber_pars.h"

extern bool NextSection(char** original, char ch);  //default ch = ' '
extern bool pass(char** original, char ch);

bool read_chars(char **buffer, char *str, char seperator) {
	int i = 0;
	char *buff = *buffer;
	if (!pass(&buff, ' ') || !pass(&buff, seperator) || !pass(&buff, ' ')) return false;
	while (strlen(buff) > 0 && *buff != seperator && *buff != ' ') {
		str[i] = *buff; i++; buff = buff + 1;
	}
	str[i] = '\0';
	*buffer = buff;
	if (strlen(str) > 0) return true;
	else return false;
}

bool is_torsion_specific_atom(char *atom) {
	if (strlen(atom) == 0) return false;
	if (strcmp(atom, "X") == 0 || strcmp(atom, "x") == 0 || strcmp(atom, "*") == 0) return false;
	else return true;
}

bool is_torsion(char *buffer, char *atom1, char *atom2, char *atom3, char *atom4, char **pars) {
	char *buff = buffer, atom[4][50];
	bool status = true;
	if (read_chars(&buff, atom[0], '-') && read_chars(&buff, atom[1], '-') && read_chars(&buff, atom[2], '-') && read_chars(&buff, atom[3], '-')) {
		if (strcmp(atom1, atom[0]) == 0 && strcmp(atom2, atom[1]) == 0 && strcmp(atom3, atom[2]) == 0 && strcmp(atom4, atom[3]) == 0) status = true;
		else if (strcmp(atom1, atom[3]) == 0 && strcmp(atom2, atom[2]) == 0 && strcmp(atom3, atom[1]) == 0 && strcmp(atom4, atom[0]) == 0) status = true;
		else status = false;
		if (status) *pars = buff;
		return status;
	}
	else return false;
}

bool find_torsion_pars(ifstream &in, bool improper, char *atom1, char *atom2, char *atom3, char *atom4, float &v0, float &phase, float &pn, char *errmsg) {
	char buffer[256] = "\0", *buff = NULL;
	float idivf = 1, pk = 0;

	while (!in.eof()) {
		in.getline(buffer, 250);
		if (is_torsion(buffer, atom1, atom2, atom3, atom4, &buff)) {
			if (improper) {
				if (sscanf(buff, "%f %f %f", &pk, &phase, &pn) == 3) {
					v0 = pk; return true;
				}
				else {sprintf(errmsg, "UNKNOWN torsion format : %s", buffer); return false;}
			}
			else {
				if (sscanf(buff, "%f %f %f %f", &idivf, &pk, &phase, &pn) != 4) {
					sprintf(errmsg, "UNKNOWN torsion format : %s", buffer); return false;
				}
				else {v0 = pk / idivf; return true;}
			}
		}
	}
	return false;
}

bool find_all_torsion_pars(char *fname, bool improper, char *atom1, char *atom2, char *atom3, char *atom4, TORSION_PARS &tpar, char *errmsg) {
	tpar.release_tail();
	float v0 = 0, phase = 0, pn = 0;
	bool status = false;
	TORSION_PARS *pTorPar = NULL;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", fname); return false;}
	status = find_torsion_pars(in, improper, atom1, atom2, atom3, atom4, v0, phase, pn, errmsg);
	if (!status) {in.close(); return status;}
	
	tpar.v0 = v0; tpar.phase = phase; tpar.pn = pn;
	pTorPar = &tpar;

	while (pn < -0.01f && find_torsion_pars(in, improper, atom1, atom2, atom3, atom4, v0, phase, pn, errmsg)) {
		pTorPar->next = new TORSION_PARS();
		pTorPar->next->v0 = v0; pTorPar->next->phase = phase; pTorPar->next->pn = pn;
		pTorPar = pTorPar->next;
	}
	in.close(); return true;
}

void show_all_torsion_pars(char *fname, char *atom1, char *atom2, char *atom3, char *atom4, bool improper) {
	ifstream in;
	bool status = true;
	float c = 0, phase = 0;
	float pn = 0;
	char errmsg[256] = "\0", buffer[256] = "\0";
	in.open(fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", fname); return;}
	while (!in.eof()) {
		status = find_torsion_pars(in, improper, atom1, atom2, atom3, atom4, c, phase, pn, errmsg);
		if (status) {
			sprintf(buffer, "%s-%s-%s-%s : %f  %f  %f", atom1, atom2, atom3, atom4, c, phase, pn);
			cout<<buffer<<endl; cout.flush();
		}
		else cout<<errmsg<<endl;
	}
	in.close();
	return;
}

bool is_bound(char *buffer, char *atom1, char *atom2, char **pars) {
	char *bf = buffer;
	char atom[2][50] = {"\0", "\0"};
	if (read_chars(&bf, atom[0], '-') && read_chars(&bf, atom[1], '-') 
		&& (
			(strcmp(atom[0], atom1) == 0 && strcmp(atom[1], atom2) == 0) ||
			(strcmp(atom[0], atom2) == 0 && strcmp(atom[1], atom1) == 0)
			)
		) {
		if (strlen(bf) > 2) {
			if (bf[0] == '-' || bf[1] == '-') return false;
		}
		*pars = bf; return true;
	}
	return false;
}

bool find_bound_pars(char *fname, char *atom1, char *atom2, float &blen, char *errmsg) {
	char buffer[256] = "\0", *pars = NULL;
	ifstream in;
	float f = 0;
	in.open(fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", fname); return false;}

	while (!in.eof()) {
		in.getline(buffer, 250);
		if (is_bound(buffer, atom1, atom2, &pars)) {
			in.close();
			if (sscanf(pars, "%f %f", &f, &blen) == 2) {in.close(); return true;}
			else {sprintf(errmsg, "UNKNOWN bound format : %s", buffer); in.close(); return false;}
		}
	}
	sprintf(errmsg, "can not find bound %s-%s", atom1, atom2); 
	in.close(); return false;
}

bool is_bound_angle(char *buffer, char *atom1, char *atom2, char *atom3, char **pars) {
	char *bf = buffer;
	char atom[3][50] = {"\0", "\0", "\0"};
	if (read_chars(&bf, atom[0], '-') && read_chars(&bf, atom[1], '-') && read_chars(&bf, atom[2], '-') 
		&& strcmp(atom[1], atom2) == 0 && (
			(strcmp(atom[0], atom1) == 0 && strcmp(atom[2], atom3) == 0) 
			|| (strcmp(atom[0], atom3) == 0 && strcmp(atom[2], atom1) == 0)
			)
		) {
		if (strlen(bf) > 2) {
			if (bf[0] == '-' || bf[1] == '-') return false;
		}
		*pars = bf; return true;
	}
	return false;
}

bool find_bound_angle_pars(char *fname, char *atom1, char *atom2, char *atom3, float &angle, char *errmsg) {
	char buffer[256] = "\0", *pars = NULL;
	float f = 0;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", fname); return false;}

	while (!in.eof()) {
		in.getline(buffer, 250);
		if (is_bound_angle(buffer, atom1, atom2, atom3, &pars)) {
			if (sscanf(pars, "%f %f", &f, &angle) == 2) {in.close(); return true;}
			else {sprintf(errmsg, "UNKNOWN bound format : %s", buffer); in.close(); return false;}
		}
	}
	sprintf(errmsg, "can not find bound angle %s-%s-%s", atom1, atom2, atom3); in.close(); return false;
}

bool general_find_torsion_pars(char *fname, bool improper, char *atom1, char *atom2, char *atom3, char *atom4, TORSION_PARS &tpar, char *errmsg) {
	bool status = true;
	status = find_all_torsion_pars(fname, improper, atom1, atom2, atom3, atom4, tpar, errmsg);
	if (status) return true;
	else {
		status = find_all_torsion_pars(fname, improper, "X", atom2, atom3, atom4, tpar, errmsg);
		if (status) return true;
		else {
			status = find_all_torsion_pars(fname, improper, atom1, atom2, atom3, "X", tpar, errmsg);
			if (status) return true;
			else {
				status = find_all_torsion_pars(fname, improper, "X", atom2, atom3, "X", tpar, errmsg);
				if (!status) sprintf(errmsg, "can not find %s-%s-%s-%s torsion parameters", atom1, atom2, atom3, atom4); 
				return status;
			}
		}
	}
}

