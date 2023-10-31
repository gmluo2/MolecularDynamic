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

extern char monomer_infor_fname[256];

/*
bool defined_simple_molecule(char *mol) {
	if (strcmp(mol, "BENZENE") == 0) return true;
	else if (strcmp(mol, "SPC_WATER") == 0) return true;
	else if (strcmp(mol, "SPCE_WATER") == 0) return true;
	else return false;
}

bool defined_cluster(char *unit) {
	if (strcmp(unit, "CH3") == 0) return true;
	else if (strcmp(unit, "CH2") == 0) return true;
	else if (strcmp(unit, "C2H2") == 0) return true;
	else if (strcmp(unit, "CH") == 0) return true;
	else if (strcmp(unit, "C2H") == 0) return true;
	else if (strcmp(unit, "BENZENE1") == 0) return true;
	else if (strcmp(unit, "BENZENE12") == 0) return true;
	else if (strcmp(unit, "BENZENE13") == 0) return true;
	else if (strcmp(unit, "BENZENE14") == 0) return true;
	else if (strcmp(unit, "SPC_O") == 0) return true;
	else if (strcmp(unit, "SPC_OH") == 0) return true;
	return false;
}

bool defined_monomer(char *unit) {
	if (strcmp(unit, "ISOPRENE") == 0) return true;
	else if (strcmp(unit, "ISOPRENE1") == 0) return true;
	else if (strcmp(unit, "ISOPRENE2") == 0) return true;
	else if (strcmp(unit, "STYRENE") == 0) return true;
	else if (strcmp(unit, "STYRENE1") == 0) return true;
	else if (strcmp(unit, "STYRENE2") == 0) return true;
	else if (strcmp(unit, "OXYETHYLENE") == 0 || strcmp(unit, "EO") == 0) return true;
	else if (strcmp(unit, "OXYETHYLENE1") == 0 || strcmp(unit, "EO1") == 0) return true;
	else if (strcmp(unit, "OXYPROPYLENE") == 0 || strcmp(unit, "PO") == 0) return true;
	else if (strcmp(unit, "OXYPROPYLENE1") == 0 || strcmp(unit, "PO1") == 0) return true;
	else if (strcmp(unit, "OXYPROPYLENE2") == 0 || strcmp(unit, "PO2") == 0) return true;
	else return false;
}

int SPCwater(BASIC_CLUSTER &c) {
	float bl = 1.0f, ba = 109.47f / 2; // rOH and <HOH> angle
	c.set_hinges(0);
	c.setAtoms(3);
	float x = bl * sin(ba * fAngle2Radian);
	float y = bl * cos(ba * fAngle2Radian);

	c.atom[0].aindx = SPC_WATER_O; //strcpy(c.atom[0].name, "O");
	c.atom[0].r.v[0] = 0; c.atom[0].r.v[1] = 0; c.atom[0].r.v[2] = 0;
	c.atom[0].set_bounds(2); // 2 bound atoms

	c.atom[1].aindx = SPC_WATER_H; //strcpy(c.atom[1].name, "H");
	c.atom[1].r.v[0] = x; c.atom[1].r.v[1] = y; c.atom[1].r.v[2] = 0;
	c.atom[1].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom, c.atom + 1, SINGLE_BOUND);

	c.atom[2].aindx = SPC_WATER_H; //strcpy(c.atom[2].name, "H");
	c.atom[2].r.v[0] = -x; c.atom[2].r.v[1] = y; c.atom[2].r.v[2] = 0;
	c.atom[2].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom, c.atom + 2, SINGLE_BOUND);

	c.d = 2 * bl * cos(ba * Angle2Radian) + 1; // H has diameter 1 Angs.
	c.cID = CINDX_SPC;

	return 3;
}

int SPCEwater(BASIC_CLUSTER &c) {
	float bl = 1.0f, ba = 109.47f / 2; // rOH and <HOH> angle
	c.set_hinges(0);
	c.setAtoms(3);
	float x = bl * sin(ba * fAngle2Radian);
	float y = bl * cos(ba * fAngle2Radian);

	c.atom[0].aindx = E_SPC_WATER_O; //strcpy(c.atom[0].name, "O");
	c.atom[0].r.v[0] = 0; c.atom[0].r.v[1] = 0; c.atom[0].r.v[2] = 0;
	c.atom[0].set_bounds(2); // 2 bound atoms

	c.atom[1].aindx = E_SPC_WATER_H; //strcpy(c.atom[1].name, "H");
	c.atom[1].r.v[0] = x; c.atom[1].r.v[1] = y; c.atom[1].r.v[2] = 0;
	c.atom[1].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom, c.atom + 1, SINGLE_BOUND);

	c.atom[2].aindx = E_SPC_WATER_H; //strcpy(c.atom[2].name, "H");
	c.atom[2].r.v[0] = -x; c.atom[2].r.v[1] = y; c.atom[2].r.v[2] = 0;
	c.atom[2].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom, c.atom + 2, SINGLE_BOUND);

	c.d = 2 * bl * cos(ba * Angle2Radian) + 1; // H has diameter 1 Angs.
	c.cID = CINDX_SPCE;

	return 3;
}

int SPC_OH(BASIC_CLUSTER &c) {
	float bl = 1.0f, ba = 109.47f / 2; // rOH and <HOH> angle
	c.set_hinges(1);
	c.setAtoms(2);
	float x = bl * sin(ba * fAngle2Radian);
	float y = bl * cos(ba * fAngle2Radian);

	c.atom[0].aindx = SPC_WATER_O; //strcpy(c.atom[0].name, "O");
	c.atom[0].r.v[0] = 0; c.atom[0].r.v[1] = 0; c.atom[0].r.v[2] = 0;
	c.atom[0].set_bounds(2); // 2 bound atoms

	c.atom[1].aindx = SPC_WATER_H; //strcpy(c.atom[1].name, "H");
	c.atom[1].r.v[0] = x; c.atom[1].r.v[1] = y; c.atom[1].r.v[2] = 0;
	c.atom[1].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom, 1, c.atom + 1, 0, SINGLE_BOUND);

	c.hinge[0].vHinge.v[0] = -x; c.hinge[0].vHinge.v[1] = y; c.hinge[0].vHinge.v[2] = 0;
	c.setup_hinge(0, 0, 0);
	c.hinge[0].vHinge -= (*(c.hinge[0].O));
	float HingeLength = 1.2f;
	c.hinge[0].vHinge.stretch(HingeLength);

	c.d = bl + 1.2; // H has diameter 1 Angs., O, 1.4 ?
	c.cID = CINDX_SPC_OH;
	return 2;
}

int OH_template(BASIC_CLUSTER &c, char Oindx, char Hindx, short cindx) {
	float bl = 1.0f, ba = 109.47f / 2; // rOH and <HOH> angle
	c.set_hinges(1);
	c.setAtoms(2);
	float x = bl * sin(ba * fAngle2Radian);
	float y = bl * cos(ba * fAngle2Radian);

	c.atom[0].aindx = Oindx;
	c.atom[0].r.v[0] = 0; c.atom[0].r.v[1] = 0; c.atom[0].r.v[2] = 0;
	c.atom[0].set_bounds(2); // 2 bound atoms

	c.atom[1].aindx = Hindx;
	c.atom[1].r.v[0] = x; c.atom[1].r.v[1] = y; c.atom[1].r.v[2] = 0;
	c.atom[1].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom, 1, c.atom + 1, 0, SINGLE_BOUND);

	c.setup_hinge(0, 0, 0);
	c.hinge[0].vHinge.v[0] = -x; c.hinge[0].vHinge.v[1] = y; c.hinge[0].vHinge.v[2] = 0;
	c.hinge[0].vHinge -= (*(c.hinge[0].O));
	//float HingeLength = 1.2f;
	//c.hinge[0].vHinge.stretch(HingeLength);

	c.d = bl + 1.2; // H has diameter 1 Angs., O, 1.4 ?
	c.cID = cindx;
	return 2;
}

int SPC_O(BASIC_CLUSTER &c) {
	float bl = 1.0f, ba = 109.47f / 2; // rOH and <HOH> angle
	c.set_hinges(2);
	c.setAtoms(1);
	float x = bl * sin(ba * fAngle2Radian);
	float y = bl * cos(ba * fAngle2Radian);

	c.atom[0].aindx = SPC_WATER_O; //strcpy(c.atom[0].name, "O");
	c.atom[0].r.v[0] = 0; c.atom[0].r.v[1] = 0; c.atom[0].r.v[2] = 0;
	c.atom[0].set_bounds(2); // 2 bound atoms

	c.hinge[0].vHinge.v[0] = x; c.hinge[0].vHinge.v[1] = y; c.hinge[0].vHinge.v[2] = 0;
	c.hinge[1].vHinge.v[0] = -x; c.hinge[1].vHinge.v[1] = y; c.hinge[1].vHinge.v[2] = 0;
	
	float HingeLength = 1.2f;
	c.setup_hinge(0, 0, 0); c.setup_hinge(1, 0, 1);
	c.hinge[0].vHinge -= (*(c.hinge[0].O)); c.hinge[1].vHinge -= (*(c.hinge[1].O));
	c.hinge[0].vHinge.stretch(HingeLength); c.hinge[1].vHinge.stretch(HingeLength);

	c.d = 1.4; // ?
	c.cID = CINDX_SPC_O;
	return 1;
}

int DoubleBoundAtom(BASIC_CLUSTER &c, char aindx, float bl, float ba, short cindx) {
	c.set_hinges(2);
	c.setAtoms(1);
	ba /= 2;
	float x = bl * sin(ba * fAngle2Radian);
	float y = bl * cos(ba * fAngle2Radian);

	c.atom[0].aindx = aindx;
	c.atom[0].r.v[0] = 0; c.atom[0].r.v[1] = 0; c.atom[0].r.v[2] = 0;
	c.atom[0].set_bounds(2); // 2 bound atoms

	c.setup_hinge(0, 0, 0); c.setup_hinge(1, 0, 1);
	c.hinge[0].vHinge.v[0] = x; c.hinge[0].vHinge.v[1] = y; c.hinge[0].vHinge.v[2] = 0;
	c.hinge[1].vHinge.v[0] = -x; c.hinge[1].vHinge.v[1] = y; c.hinge[1].vHinge.v[2] = 0;
	
	//float HingeLength = 1.2f;
	//c.hinge[0].vHinge -= (*(c.hinge[0].O)); c.hinge[1].vHinge -= (*(c.hinge[1].O));
	//c.hinge[0].vHinge.stretch(HingeLength); c.hinge[1].vHinge.stretch(HingeLength);

	c.d = bl + 1.4; //?
	c.cID = cindx;
	return 1;
}

// to test the functions of Z-Matrix Operation
// Constructing a benzene ring
int Benzene(BASIC_CLUSTER &c) {
	float bl1 = 1.54f, ba1 = 120;
	float bl2 = 1.0f, ba2 = 120;

	c.set_hinges(0);

	int nAtoms = 12;
	c.setAtoms(nAtoms); //C6H6
	c.atom[0].aindx = BENZENE_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(3); // 3 bound atoms

	c.atom[1].aindx = BENZENE_C; //strcpy(c.atom[1].name, "C");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl1), c.atom[1].r);
	c.atom[1].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom, c.atom + 1, SINGLE_BOUND);

	c.atom[2].aindx = BENZENE_C; //strcpy(c.atom[2].name, "C");
	ZMatrix(c.atom[1].r, c.atom[0].r, VECTOR3(0, 1, 0), bl1, ba1, c.atom[2].r);
	c.atom[2].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 1, c.atom + 2, DOUBLE_BOUND);

	c.atom[3].aindx = BENZENE_C; // strcpy(c.atom[3].name, "C");
	ZMatrix(c.atom[2].r, c.atom[1].r, c.atom[0].r, bl1, ba1, 0, c.atom[3].r);
	c.atom[3].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 2, c.atom + 3);

	c.atom[4].aindx = BENZENE_C; //strcpy(c.atom[4].name, "C");
	ZMatrix(c.atom[3].r, c.atom[2].r, c.atom[1].r, bl1, ba1, 0, c.atom[4].r);
	c.atom[4].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 3, c.atom + 4, DOUBLE_BOUND);

	c.atom[5].aindx = BENZENE_C; //strcpy(c.atom[5].name, "C");
	ZMatrix(c.atom[4].r, c.atom[3].r, c.atom[2].r, bl1, ba1, 0, c.atom[5].r);
	c.atom[5].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 4, c.atom + 5);

	bound_connect(c.atom + 5, c.atom, DOUBLE_BOUND);

	c.atom[6].aindx = BENZENE_H; //strcpy(c.atom[6].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba2, 180, c.atom[6].r);
	c.atom[6].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 6, c.atom);

	c.atom[7].aindx = BENZENE_H; //strcpy(c.atom[7].name, "H");
	ZMatrix(c.atom[1].r, c.atom[2].r, c.atom[3].r, bl2, ba2, 180, c.atom[7].r);
	c.atom[7].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 7, c.atom + 1);

	c.atom[8].aindx = BENZENE_H; //strcpy(c.atom[8].name, "H");
	ZMatrix(c.atom[2].r, c.atom[3].r, c.atom[4].r, bl2, ba2, 180, c.atom[8].r);
	c.atom[8].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 8, c.atom + 2);

	c.atom[9].aindx = BENZENE_H; //strcpy(c.atom[9].name, "H");
	ZMatrix(c.atom[3].r, c.atom[4].r, c.atom[5].r, bl2, ba2, 180, c.atom[9].r);
	c.atom[9].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 9, c.atom + 3);

	c.atom[10].aindx = BENZENE_H; //strcpy(c.atom[10].name, "H");
	ZMatrix(c.atom[4].r, c.atom[5].r, c.atom[0].r, bl2, ba2, 180, c.atom[10].r);
	c.atom[10].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 10, c.atom + 4);

	c.atom[11].aindx = BENZENE_H; //strcpy(c.atom[11].name, "H");
	ZMatrix(c.atom[5].r, c.atom[0].r, c.atom[1].r, bl2, ba2, 180, c.atom[11].r);
	c.atom[11].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 11, c.atom + 5);

	c.n2Parent = -1; // no atoms to connect to parent

	c.d = 5 + 1; // H diameter 1 ?

	//c.align(c.atom[0].r, c.atom[3].r, VECTOR3(0, 0, -bl1), VECTOR3(0, 0, 1));
	c.cID = CINDX_BENZENE;
	return nAtoms; // number of atoms
}

// Constructing a benzene ring with one open link to parent
int Benzene1(CLUSTER &c) {
	float bl1 = 1.54f, ba1 = 120;
	float bl2 = 1.0f, ba2 = 120;

	c.set_hinges(1);

	int nAtoms = 11;
	c.setAtoms(nAtoms); //C6H6
	c.atom[0].aindx = BENZENE_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(3); // 3 bound atoms

	c.atom[1].aindx = BENZENE_C; //strcpy(c.atom[1].name, "C");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl1), c.atom[1].r);
	c.atom[1].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom, c.atom + 1);

	c.atom[2].aindx = BENZENE_C; //strcpy(c.atom[2].name, "C");
	ZMatrix(c.atom[1].r, c.atom[0].r, VECTOR3(0, 1, 0), bl1, ba1, c.atom[2].r);
	c.atom[2].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 1, c.atom + 2, DOUBLE_BOUND);

	c.atom[3].aindx = BENZENE_C; //strcpy(c.atom[3].name, "C");
	ZMatrix(c.atom[2].r, c.atom[1].r, c.atom[0].r, bl1, ba1, 0, c.atom[3].r);
	c.atom[3].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 2, c.atom + 3);

	c.atom[4].aindx = BENZENE_C; //strcpy(c.atom[4].name, "C");
	ZMatrix(c.atom[3].r, c.atom[2].r, c.atom[1].r, bl1, ba1, 0, c.atom[4].r);
	c.atom[4].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 3, c.atom + 4, DOUBLE_BOUND);

	c.atom[5].aindx = BENZENE_C; //strcpy(c.atom[5].name, "C");
	ZMatrix(c.atom[4].r, c.atom[3].r, c.atom[2].r, bl1, ba1, 0, c.atom[5].r);
	c.atom[5].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 4, c.atom + 5);

	bound_connect(c.atom + 5, c.atom, DOUBLE_BOUND);

	c.atom[6].aindx = BENZENE_H; //strcpy(c.atom[6].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba2, 180, c.atom[6].r);
	c.atom[6].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 6, c.atom);

	c.atom[7].aindx = BENZENE_H; //strcpy(c.atom[7].name, "H");
	ZMatrix(c.atom[1].r, c.atom[2].r, c.atom[3].r, bl2, ba2, 180, c.atom[7].r);
	c.atom[7].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 7, c.atom + 1);

	c.atom[8].aindx = BENZENE_H; //strcpy(c.atom[8].name, "H");
	ZMatrix(c.atom[2].r, c.atom[3].r, c.atom[4].r, bl2, ba2, 180, c.atom[8].r);
	c.atom[8].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 8, c.atom + 2);

	c.atom[9].aindx = BENZENE_H; //strcpy(c.atom[9].name, "H");
	ZMatrix(c.atom[3].r, c.atom[4].r, c.atom[5].r, bl2, ba2, 180, c.atom[9].r);
	c.atom[9].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 9, c.atom + 3);

	c.atom[10].aindx = BENZENE_H; //strcpy(c.atom[10].name, "H");
	ZMatrix(c.atom[4].r, c.atom[5].r, c.atom[0].r, bl2, ba2, 180, c.atom[10].r);
	c.atom[10].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 10, c.atom + 4);

	ZMatrix(c.atom[5].r, c.atom[0].r, c.atom[1].r, bl2, ba2, 180, c.hinge[0].vHinge);

	c.setup_hinge(0, 5, 2);
	c.hinge[0].vHinge -= (*(c.hinge[0].O));
	c.hinge[0].vHinge.stretch(HingeLength);

	c.align(c.hinge[0].vHinge, c.atom[5].r, VECTOR3(0, 0, 0), VECTOR3(0, 0, 1));

	c.d = 5 + 1; // H diameter 1 ?
	c.cID = CINDX_BENZENE1;

	return nAtoms; // number of atoms
}

int Benzene14(CLUSTER &c) { // -benzene- -C6H4- with hing at C1 & C4
	float bl1 = 1.54f, ba1 = 120;
	float bl2 = 1.0f, ba2 = 120;

	c.set_hinges(2); // 2 hinges

	int nAtoms = 10;
	c.setAtoms(nAtoms); //C6H4
	c.atom[0].aindx = BENZENE_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(3); // 3 bound atoms

	c.atom[1].aindx = BENZENE_C; //strcpy(c.atom[1].name, "C");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl1), c.atom[1].r);
	c.atom[1].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 1, c.atom);

	c.atom[2].aindx = BENZENE_C; //strcpy(c.atom[2].name, "C");
	ZMatrix(c.atom[1].r, c.atom[0].r, VECTOR3(0, 1, 0), bl1, ba1, c.atom[2].r);
	c.atom[2].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 2, c.atom + 1, DOUBLE_BOUND);

	c.atom[3].aindx = BENZENE_C; //strcpy(c.atom[3].name, "C");
	ZMatrix(c.atom[2].r, c.atom[1].r, c.atom[0].r, bl1, ba1, 0, c.atom[3].r);
	c.atom[3].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 3, c.atom + 2);

	c.atom[4].aindx = BENZENE_C; //strcpy(c.atom[4].name, "C");
	ZMatrix(c.atom[3].r, c.atom[2].r, c.atom[1].r, bl1, ba1, 0, c.atom[4].r);
	c.atom[4].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 4, c.atom + 3, DOUBLE_BOUND);

	c.atom[5].aindx = BENZENE_C; //strcpy(c.atom[5].name, "C");
	ZMatrix(c.atom[4].r, c.atom[3].r, c.atom[2].r, bl1, ba1, 0, c.atom[5].r);
	c.atom[5].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 5, c.atom + 4);
	
	bound_connect(c.atom + 5, c.atom, DOUBLE_BOUND);

	// Open bound for parent
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba2, 180, c.hinge[0].vHinge);
	//strcpy(c.atom[6].name, "H");
	//ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba2, 180, c.atom[6].r);

	c.atom[6].aindx = BENZENE_H; //strcpy(c.atom[6].name, "H");
	ZMatrix(c.atom[1].r, c.atom[2].r, c.atom[3].r, bl2, ba2, 180, c.atom[6].r);
	c.atom[6].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 6, c.atom + 1);

	c.atom[7].aindx = BENZENE_H; //strcpy(c.atom[7].name, "H");
	ZMatrix(c.atom[2].r, c.atom[3].r, c.atom[4].r, bl2, ba2, 180, c.atom[7].r);
	c.atom[7].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 7, c.atom + 2);

	//Open bound for child
	ZMatrix(c.atom[3].r, c.atom[4].r, c.atom[5].r, bl2, ba2, 180, c.hinge[1].vHinge); //actually here we get the hing end point
	//strcpy(c.atom[8].name, "H");
	//ZMatrix(c.atom[3].r, c.atom[4].r, c.atom[5].r, bl2, ba2, 180, c.atom[8].r);

	c.atom[8].aindx = BENZENE_H; //strcpy(c.atom[8].name, "H");
	ZMatrix(c.atom[4].r, c.atom[5].r, c.atom[0].r, bl2, ba2, 180, c.atom[8].r);
	c.atom[8].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 8, c.atom + 4);

	c.atom[9].aindx = BENZENE_H; //strcpy(c.atom[9].name, "H");
	ZMatrix(c.atom[5].r, c.atom[0].r, c.atom[1].r, bl2, ba2, 180, c.atom[9].r);
	c.atom[9].set_bounds(1); // 1 bound atoms
	bound_connect(c.atom + 9, c.atom + 5);
	
	c.setup_hinge(0, 0, 2);
	c.hinge[0].vHinge -= (*(c.hinge[0].O));
	c.setup_hinge(1, 3, 2);
	c.hinge[1].vHinge -= (*(c.hinge[1].O));

	c.hinge[0].vHinge.stretch(HingeLength);
	c.hinge[1].vHinge.stretch(HingeLength);

	c.align(c.atom[0].r, c.atom[3].r, VECTOR3(0, 0, 0), VECTOR3(0, 0, 1));

	c.d = 5 + 1; // H diameter 1 ?
	c.cID = CINDX_BENZENE14;

	return nAtoms; // number of atoms
}

int Benzene13(CLUSTER &c) { // -benzene- -C6H4- with hing at C1 & C3
	float bl1 = 1.54f, ba1 = 120;
	float bl2 = 1.0f, ba2 = 120;

	c.set_hinges(2); // 2 hinges

	int nAtoms = 10;
	c.setAtoms(nAtoms); //C6H4
	c.atom[0].aindx = BENZENE_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(3); // 3 bound atoms

	c.atom[1].aindx = BENZENE_C; //strcpy(c.atom[1].name, "C");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl1), c.atom[1].r);
	c.atom[1].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 1, c.atom + 0);

	c.atom[2].aindx = BENZENE_C; //strcpy(c.atom[2].name, "C");
	ZMatrix(c.atom[1].r, c.atom[0].r, VECTOR3(0, 1, 0), bl1, ba1, c.atom[2].r);
	c.atom[2].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 2, c.atom + 1, DOUBLE_BOUND);

	c.atom[3].aindx = BENZENE_C; //strcpy(c.atom[3].name, "C");
	ZMatrix(c.atom[2].r, c.atom[1].r, c.atom[0].r, bl1, ba1, 0, c.atom[3].r);
	c.atom[3].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 3, c.atom + 2);

	c.atom[4].aindx = BENZENE_C; //strcpy(c.atom[4].name, "C");
	ZMatrix(c.atom[3].r, c.atom[2].r, c.atom[1].r, bl1, ba1, 0, c.atom[4].r);
	c.atom[4].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 4, c.atom + 3, DOUBLE_BOUND);

	c.atom[5].aindx = BENZENE_C; //strcpy(c.atom[5].name, "C");
	ZMatrix(c.atom[4].r, c.atom[3].r, c.atom[2].r, bl1, ba1, 0, c.atom[5].r);
	c.atom[5].set_bounds(3); // 3 bound atoms
	bound_connect(c.atom + 5, c.atom + 4);

	bound_connect(c.atom + 5, c.atom, DOUBLE_BOUND);

	// Open bound for parent
	c.setup_hinge(0, 0, 2);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba2, 180, c.hinge[0].vHinge);
	c.hinge[0].vHinge -= c.atom[0].r;
	c.hinge[0].vHinge.stretch(HingeLength);
	//strcpy(c.atom[6].name, "H");
	//ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba2, 180, c.atom[6].r);

	c.atom[6].aindx = BENZENE_H; //strcpy(c.atom[6].name, "H");
	ZMatrix(c.atom[1].r, c.atom[2].r, c.atom[3].r, bl2, ba2, 180, c.atom[6].r);
	c.atom[6].set_bounds(1);
	bound_connect(c.atom + 6, c.atom + 1);

	//Open bound for child
	c.setup_hinge(1, 2, 2);
	ZMatrix(c.atom[2].r, c.atom[3].r, c.atom[4].r, bl2, ba2, 180, c.hinge[1].vHinge); //actually here we get the hing end point
	c.hinge[1].vHinge -= c.atom[2].r;
	c.hinge[1].vHinge.stretch(HingeLength);
	//strcpy(c.atom[7].name, "H");
	//ZMatrix(c.atom[2].r, c.atom[3].r, c.atom[4].r, bl2, ba2, 180, c.atom[7].r);

	c.atom[7].aindx = BENZENE_H; //strcpy(c.atom[7].name, "H");
	ZMatrix(c.atom[3].r, c.atom[4].r, c.atom[5].r, bl2, ba2, 180, c.atom[7].r);
	c.atom[7].set_bounds(1);
	bound_connect(c.atom + 7, c.atom + 3);

	c.atom[8].aindx = BENZENE_H; //strcpy(c.atom[8].name, "H");
	ZMatrix(c.atom[4].r, c.atom[5].r, c.atom[0].r, bl2, ba2, 180, c.atom[8].r);
	c.atom[8].set_bounds(1);
	bound_connect(c.atom + 8, c.atom + 4);

	c.atom[9].aindx = BENZENE_H; //strcpy(c.atom[9].name, "H");
	ZMatrix(c.atom[5].r, c.atom[0].r, c.atom[1].r, bl2, ba2, 180, c.atom[9].r);
	c.atom[9].set_bounds(1);
	bound_connect(c.atom + 9, c.atom + 5);

	c.align(c.atom[0].r, c.atom[2].r, VECTOR3(0, 0, 0), VECTOR3(0, 0, 1));

	c.d = 5 + 1; // H diameter 1 ?
	c.cID = CINDX_BENZENE13;

	return nAtoms; // number of atoms
}

int Benzene12(CLUSTER &c) { // -benzene- -C6H4- with hing at C1 & C2
	float bl1 = 1.54f, ba1 = 120;
	float bl2 = 1.0f, ba2 = 120;

	c.set_hinges(2); // 2 hinges

	int nAtoms = 10;
	c.setAtoms(nAtoms); //C6H4
	c.atom[0].aindx = BENZENE_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(3);

	c.atom[1].aindx = BENZENE_C; //strcpy(c.atom[1].name, "C");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl1), c.atom[1].r);
	c.atom[1].set_bounds(3);
	bound_connect(c.atom + 1, c.atom + 0);

	c.atom[2].aindx = BENZENE_C; //strcpy(c.atom[2].name, "C");
	ZMatrix(c.atom[1].r, c.atom[0].r, VECTOR3(0, 1, 0), bl1, ba1, c.atom[2].r);
	c.atom[2].set_bounds(3);
	bound_connect(c.atom + 2, c.atom + 1, DOUBLE_BOUND);

	c.atom[3].aindx = BENZENE_C; //strcpy(c.atom[3].name, "C");
	ZMatrix(c.atom[2].r, c.atom[1].r, c.atom[0].r, bl1, ba1, 0, c.atom[3].r);
	c.atom[3].set_bounds(3);
	bound_connect(c.atom + 3, c.atom + 2);

	c.atom[4].aindx = BENZENE_C; //strcpy(c.atom[4].name, "C");
	ZMatrix(c.atom[3].r, c.atom[2].r, c.atom[1].r, bl1, ba1, 0, c.atom[4].r);
	c.atom[4].set_bounds(3);
	bound_connect(c.atom + 4, c.atom + 3, DOUBLE_BOUND);

	c.atom[5].aindx = BENZENE_C; //strcpy(c.atom[5].name, "C");
	ZMatrix(c.atom[4].r, c.atom[3].r, c.atom[2].r, bl1, ba1, 0, c.atom[5].r);
	c.atom[5].set_bounds(3);
	bound_connect(c.atom + 5, c.atom + 4);

	bound_connect(c.atom + 5, c.atom, DOUBLE_BOUND);

	// Open bound for parent
	c.setup_hinge(0, 0, 2);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba2, 180, c.hinge[0].vHinge);  //actually here we get the hing end point
	c.hinge[0].vHinge -= c.atom[0].r;
	c.hinge[0].vHinge.stretch(HingeLength);
	//strcpy(c.atom[6].name, "H");
	//ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba2, 180, c.atom[6].r);

	//Open bound for child
	c.setup_hinge(1, 1, 2);
	ZMatrix(c.atom[1].r, c.atom[2].r, c.atom[3].r, bl2, ba2, 180, c.hinge[1].vHinge); //actually here we get the hing end point
	c.hinge[1].vHinge -= c.atom[1].r;
	c.hinge[1].vHinge.stretch(HingeLength);
	//strcpy(c.atom[6].name, "H");
	//ZMatrix(c.atom[1].r, c.atom[2].r, c.atom[3].r, bl2, ba2, 180, c.atom[6].r);

	c.atom[6].aindx = BENZENE_H; //strcpy(c.atom[6].name, "H");
	ZMatrix(c.atom[2].r, c.atom[3].r, c.atom[4].r, bl2, ba2, 180, c.atom[6].r);
	c.atom[6].set_bounds(1);
	bound_connect(c.atom + 6, c.atom + 2);

	c.atom[7].aindx = BENZENE_H; //strcpy(c.atom[7].name, "H");
	ZMatrix(c.atom[3].r, c.atom[4].r, c.atom[5].r, bl2, ba2, 180, c.atom[7].r);
	c.atom[7].set_bounds(1);
	bound_connect(c.atom + 7, c.atom + 3);

	c.atom[8].aindx = BENZENE_H; //strcpy(c.atom[8].name, "H");
	ZMatrix(c.atom[4].r, c.atom[5].r, c.atom[0].r, bl2, ba2, 180, c.atom[8].r);
	c.atom[8].set_bounds(1);
	bound_connect(c.atom + 8, c.atom + 4);

	c.atom[9].aindx = BENZENE_H; //strcpy(c.atom[9].name, "H");
	ZMatrix(c.atom[5].r, c.atom[0].r, c.atom[1].r, bl2, ba2, 180, c.atom[9].r);
	c.atom[9].set_bounds(1);
	bound_connect(c.atom + 9, c.atom + 5);

	c.align(c.atom[0].r, c.atom[1].r, VECTOR3(0, 0, 0), VECTOR3(0, 0, 1));

	c.d = 5 + 1; // H diameter 1 ?
	c.cID = CINDX_BENZENE12;

	return nAtoms; // number of atoms
}

// Constructing CH4
int CH4(CLUSTER &c) {
	float bl1 = 1.52f; // C-C
	float bl2 = 1.13f; // C-H
	float ba = 109.5f; // bound angle, C-C-C or H-C-H
	float ta = 120;    // torsion angle of CH4 structure

	c.set_hinges(0); // no hinge

	int nAtoms = 5;
	c.setAtoms(nAtoms);
	c.atom[0].aindx = SB_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(4);

	c.atom[1].aindx = SB_H; //strcpy(c.atom[1].name, "H");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl2), c.atom[1].r);
	c.atom[1].set_bounds(1);
	bound_connect(c.atom + 1, c.atom + 0);

	c.atom[2].aindx = SB_H; //strcpy(c.atom[2].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, VECTOR3(0, 1, 0), bl2, ba, c.atom[2].r);
	c.atom[2].set_bounds(1);
	bound_connect(c.atom + 2, c.atom + 0);

	c.atom[3].aindx = SB_H; //strcpy(c.atom[3].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba, ta, c.atom[3].r);
	c.atom[3].set_bounds(1);
	bound_connect(c.atom + 3, c.atom + 0);

	c.atom[4].aindx = SB_H; //strcpy(c.atom[4].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba, -ta, c.atom[4].r);
	c.atom[4].set_bounds(1);
	bound_connect(c.atom + 4, c.atom + 0);

	//c.align(c.atom[3].r, c.atom[0].r, VECTOR3(0, 0, 0), VECTOR3(0, 0, 1));

	c.d = 2.6 + 1; // H diameter 1 ? 2.6 is overestimated ?
	c.cID = CINDX_CH4;

	return nAtoms;
}

int CH3(CLUSTER &c) { // CH3-
	float bl1 = 1.52f; // C-C
	float bl2 = 1.13f; // C-H
	float ba = 109.5f; // bound angle, C-C-C or H-C-H
	float ta = 120;    // torsion angle of CH4 structure

	c.set_hinges(1); // 1 hinge

	int nAtoms = 4;
	c.setAtoms(nAtoms);
	c.atom[0].aindx = SB_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(4);

	c.atom[1].aindx = SB_H; //strcpy(c.atom[1].name, "H");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl2), c.atom[1].r);
	c.atom[1].set_bounds(1);
	bound_connect(c.atom + 1, c.atom + 0);

	c.atom[2].aindx = SB_H; //strcpy(c.atom[2].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, VECTOR3(0, 1, 0), bl2, ba, c.atom[2].r);
	c.atom[2].set_bounds(1);
	bound_connect(c.atom + 2, c.atom + 0);

	c.atom[3].aindx = SB_H; //strcpy(c.atom[3].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba, ta, c.atom[3].r);
	c.atom[3].set_bounds(1);
	bound_connect(c.atom + 3, c.atom + 0);

	c.setup_hinge(0, 0, 3);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba, -ta, c.hinge[0].vHinge);
	c.hinge[0].vHinge -= c.atom[0].r;
	c.hinge[0].vHinge.stretch(HingeLength);

	c.d = 2.6 + 1; // H diameter 1 ? 2.6 is overestimated ?
	c.cID = CINDX_CH3;

	return nAtoms;
}

int CH2(CLUSTER &c) { //-CH2-
	float bl1 = 1.52f; // C-C
	float bl2 = 1.13f; // C-H
	float ba = 109.5f; // bound angle, C-C-C or H-C-H
	float ta = 120;    // torsion angle of CH4 structure

	c.set_hinges(2); // 2 hinges

	int nAtoms = 3;
	c.setAtoms(3);
	c.atom[0].aindx = SB_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(4);

	c.atom[1].aindx = SB_H; //strcpy(c.atom[1].name, "H");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl2), c.atom[1].r);
	c.atom[1].set_bounds(1);
	bound_connect(c.atom + 1, c.atom + 0);

	c.atom[2].aindx = SB_H; //strcpy(c.atom[2].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, VECTOR3(0, 1, 0), bl2, ba, c.atom[2].r);
	c.atom[2].set_bounds(1);
	bound_connect(c.atom + 2, c.atom + 0);

	c.setup_hinge(0, 0, 2);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba, ta, c.hinge[0].vHinge);
	c.hinge[0].vHinge -= c.atom[0].r;
	c.hinge[0].vHinge.stretch(HingeLength);

	c.setup_hinge(1, 0, 3);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl2, ba, -ta, c.hinge[1].vHinge);
	c.hinge[1].vHinge -= c.atom[0].r;
	c.hinge[1].vHinge.stretch(HingeLength);

	c.d = bl2 + 1; // H diameter 1 ? 
	c.cID = CINDX_CH2;

	return nAtoms;
}

int CH2_template(CLUSTER &c, char Cindx, char Hindx, float bl, short cindx) { //-CH2-
	float ba = 109.5f; // bound angle, C-C-C or H-C-H
	float ta = 120;    // torsion angle of CH4 structure

	c.set_hinges(2); // 2 hinges

	int nAtoms = 3;
	c.setAtoms(3);
	c.atom[0].aindx = Cindx;
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(4);

	c.atom[1].aindx = Hindx;
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl), c.atom[1].r);
	c.atom[1].set_bounds(1);
	bound_connect(c.atom + 1, c.atom + 0);

	c.atom[2].aindx = Hindx;
	ZMatrix(c.atom[0].r, c.atom[1].r, VECTOR3(0, 1, 0), bl, ba, c.atom[2].r);
	c.atom[2].set_bounds(1);
	bound_connect(c.atom + 2, c.atom + 0);

	c.setup_hinge(0, 0, 2);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl, ba, ta, c.hinge[0].vHinge);
	c.hinge[0].vHinge -= c.atom[0].r;
	//c.hinge[0].vHinge.stretch(HingeLength);

	c.setup_hinge(1, 0, 3);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl, ba, -ta, c.hinge[1].vHinge);
	c.hinge[1].vHinge -= c.atom[0].r;
	//c.hinge[1].vHinge.stretch(HingeLength);

	c.d = bl + 1; // H diameter 1 ? overestimated ?
	c.cID = cindx;

	return nAtoms;
}

int CH(CLUSTER &c) { // -CH-- ; one parent, 2 children
	float bl1 = 1.52f; // C-C
	float bl2 = 1.13f; // C-H
	float ba = 109.5f; // bound angle, C-C-C or H-C-H
	float ta = 120;    // torsion angle of CH4 structure

	c.set_hinges(3); // 3 hinges

	int nAtoms = 2;
	c.setAtoms(2);
	c.atom[0].aindx = SB_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(4);

	c.atom[1].aindx = SB_H; //strcpy(c.atom[1].name, "H");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl2), c.atom[1].r);
	c.atom[1].set_bounds(1);
	bound_connect(c.atom + 1, c.atom + 0, SINGLE_BOUND);

	c.setup_hinge(0, 0, 1);
	ZMatrix(c.atom[0].r, c.atom[1].r, VECTOR3(0, 1, 0), bl2, ba, c.hinge[0].vHinge);

	c.setup_hinge(1, 0, 2);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.hinge[0].vHinge, bl2, ba, ta, c.hinge[1].vHinge);

	c.setup_hinge(2, 0, 3);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.hinge[0].vHinge, bl2, ba, -ta, c.hinge[2].vHinge);

	c.hinge[0].vHinge -= c.atom[0].r;
	c.hinge[0].vHinge.stretch(HingeLength);

	c.hinge[1].vHinge -= c.atom[0].r;
	c.hinge[1].vHinge.stretch(HingeLength);

	c.hinge[2].vHinge -= c.atom[0].r;
	c.hinge[2].vHinge.stretch(HingeLength);

	c.d = bl2 + 1.2; // overestimated ?
	c.cID = CINDX_CH;

	return nAtoms;
}

int C2H2(CLUSTER &c) { // -CH=CH-
	float bl1 = 1.356f; // C=C
	float bl2 = 1.13f; // C-H
	float ba = 120; // bound angle, H-C=C
	float ta = 180;    // torsion angle of H-C=C-H structure

	c.set_hinges(2); // 2 hinges

	int nAtoms = 4;
	c.setAtoms(nAtoms);
	c.atom[0].aindx = DB_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(3);

	c.atom[1].aindx = DB_C; //strcpy(c.atom[1].name, "C");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl1), c.atom[1].r);
	c.atom[1].set_bounds(3);
	bound_connect(c.atom + 1, c.atom + 0, DOUBLE_BOUND);

	c.atom[2].aindx = DB_H; //strcpy(c.atom[2].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, VECTOR3(0, 1, 0), bl2, ba, c.atom[2].r);
	c.atom[2].set_bounds(1);
	bound_connect(c.atom + 2, c.atom + 0);

	c.atom[3].aindx = DB_H; //strcpy(c.atom[3].name, "H");
	ZMatrix(c.atom[1].r, c.atom[0].r, c.atom[2].r, bl2, ba, ta, c.atom[3].r);
	c.atom[3].set_bounds(1);
	bound_connect(c.atom + 3, c.atom + 1);

	c.setup_hinge(0, 0, 2);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[3].r, bl2, ba, 0, c.hinge[0].vHinge);  //actually here we get the hing end point
	c.hinge[0].vHinge -= c.atom[0].r;
	c.hinge[0].vHinge.stretch(HingeLength);
	
	c.setup_hinge(1, 1, 2);
	ZMatrix(c.atom[1].r, c.atom[0].r, c.atom[2].r, bl2, ba, 0, c.hinge[1].vHinge);  //actually here we get the hing end point
	c.hinge[1].vHinge -= c.atom[1].r;
	c.hinge[1].vHinge.stretch(HingeLength);

	c.d = 4.5 + 1; // ?
	c.cID = CINDX_C2H2;

	return nAtoms;
}

int C2H(CLUSTER &c) { // -CH=C-- ; one parent, 2 children
	float bl1 = 1.356f; // C=C
	float bl2 = 1.13f; // C-H
	float ba = 120; // bound angle, H-C=C
	float ta = 180;  // torsion angle of H-C=C-H structure
	float bl = 1.52f;

	c.set_hinges(3);  // 3 hinges

	int nAtoms = 3;
	c.setAtoms(nAtoms);
	c.atom[0].aindx = DB_C; //strcpy(c.atom[0].name, "C");
	c.atom[0].r = VECTOR3(0, 0, 0);
	c.atom[0].set_bounds(3);

	c.atom[1].aindx = DB_C; //strcpy(c.atom[1].name, "C");
	ZMatrix(c.atom[0].r, VECTOR3(0, 0, bl1), c.atom[1].r);
	c.atom[1].set_bounds(3);
	bound_connect(c.atom + 1, c.atom + 0, DOUBLE_BOUND);

	// H connects to C0
	c.atom[2].aindx = DB_H; //strcpy(c.atom[2].name, "H");
	ZMatrix(c.atom[0].r, c.atom[1].r, VECTOR3(0, 1, 0), bl2, ba, c.atom[2].r);
	c.atom[2].set_bounds(1);
	bound_connect(c.atom + 2, c.atom + 0, SINGLE_BOUND);

	//parent hinge on C0
	c.setup_hinge(0, 0, 2);
	ZMatrix(c.atom[0].r, c.atom[1].r, c.atom[2].r, bl, ba, ta, c.hinge[0].vHinge);
	c.hinge[0].vHinge -= c.atom[0].r;
	c.hinge[0].vHinge.stretch(HingeLength);

	// child hinge #1
	c.setup_hinge(1, 1, 1);
	ZMatrix(c.atom[1].r, c.atom[0].r, c.atom[2].r, bl, ba, 0, c.hinge[1].vHinge);
	c.hinge[1].vHinge -= c.atom[1].r;
	c.hinge[1].vHinge.stretch(HingeLength);

	// child hinge #1
	c.setup_hinge(2, 1, 2);
	ZMatrix(c.atom[1].r, c.atom[0].r, c.atom[2].r, bl, ba, ta, c.hinge[2].vHinge);
	c.hinge[2].vHinge -= c.atom[1].r;
	c.hinge[2].vHinge.stretch(HingeLength);

	c.d = 4.5 + 1; // ?
	c.cID = CINDX_C2H;

	return nAtoms;
}

int isoprene_monomer(MONOMER &m) { //-CH2-C2H-[CH3]CH2-
	// 4 cluster
	m.reset(4, 3); // 3 relationship
	m.setHinges(2); // 2 hinges
	int na = 0;
	na = CH2(m.cluster[0]); // cluster 0 -- CH2
	na = C2H(m.cluster[1]); // cluster 1 -- C2H : -CH=C=
	na = CH3(m.cluster[2]); // cluster 2 -- CH3
	na = CH2(m.cluster[3]); // cluster 3 -- CH2

	m.cluster[0].rotate(m.cluster[0].hinge[0].vHinge, -45, false);
	m.cluster[1].rotate(m.cluster[1].hinge[0].vHinge, 30, false);
	m.cluster[3].rotate(m.cluster[3].hinge[0].vHinge, -165, false);

	// cluster 0 has a free hinge # 0
	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 1, 0, m.cluster, 1, SINGLE_BOUND);
	m.rship[0].set_relationship(1, 0, 0, 1, SINGLE_BOUND);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 1, 1, m.cluster + 2, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(1, 1, 2, 0, SINGLE_BOUND);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 1, 2, m.cluster + 3, 0, SINGLE_BOUND);
	m.rship[2].set_relationship(1, 2, 3, 0, SINGLE_BOUND);
	// cluster 0 has free hinge # 0; cluster 3 has free hinge # 1

	m.nc[0] = 0; m.iHinge[0] = 0;
	m.nc[1] = 3; m.iHinge[1] = 1;

	m.nbase = 0;

	return 4; // the base cluster
}

int isoprene1_monomer(MONOMER &m) { //CH3-C2H-[CH3]CH2-
	// 4 cluster
	m.reset(4, 3); // 3 relationship
	int na = 0;
	na = CH3(m.cluster[0]); // cluster 0 -- CH3
	na = CH3(m.cluster[1]); // cluster 1 -- CH3
	na = C2H(m.cluster[2]); // cluster 2 -- C2H : -CH=C=
	na = CH2(m.cluster[3]); // cluster 3 -- CH2

	m.cluster[0].rotate(m.cluster[0].hinge[0].vHinge, -45, false);
	m.cluster[3].rotate(m.cluster[3].hinge[0].vHinge, -165, false);

	//m.cluster[3].rotate(m.cluster[3].hinge[0].vHinge, -90, false);

	m.setHinges(1); // 1 free hinges

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 2, 0, m.cluster, 0, SINGLE_BOUND);
	m.rship[0].set_relationship(2, 0, 0, 0, SINGLE_BOUND);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 2, 1, m.cluster + 1, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(2, 1, 1, 0, SINGLE_BOUND);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 2, 2, m.cluster + 3, 0, SINGLE_BOUND);
	m.rship[2].set_relationship(2, 2, 3, 0, SINGLE_BOUND);
	// cluster 3 has a free hinge # 1

	m.nc[0] = 3;
	m.iHinge[0] = 1;

	m.nbase = 0;

	return 4; // base cluster
}

int isoprene2_monomer(MONOMER &m) { //-CH2-C2H-[CH3]CH3
	// 4 cluster
	m.reset(4, 3); // 3 relationship
	int na = 0;
	na = CH3(m.cluster[0]); // cluster 0 -- CH3
	na = CH3(m.cluster[1]); // cluster 1 -- CH3
	na = C2H(m.cluster[2]); // cluster 2 -- C2H : -CH=C=
	na = CH2(m.cluster[3]); // cluster 3 -- CH2

	m.setHinges(1); // 1 free hinge

	//m.cluster[0].rotate(m.cluster[0].hinge[0].vHinge, -45, false);
	m.cluster[2].rotate(m.cluster[2].hinge[0].vHinge, 30, false);
	m.cluster[3].rotate(m.cluster[3].hinge[1].vHinge, -60, false);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 2, 0, m.cluster + 3, 1, SINGLE_BOUND);
	m.rship[0].set_relationship(2, 0, 3, 1, SINGLE_BOUND);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 2, 1, m.cluster, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(2, 1, 0, 0, SINGLE_BOUND);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 2, 2, m.cluster + 1, 0, SINGLE_BOUND);
	m.rship[2].set_relationship(2, 2, 1, 0, SINGLE_BOUND);
	// cluster 3 has free hinge # 0

	m.nc[0] = 3; m.iHinge[0] = 0;

	m.nbase = 0;

	return 4; // the base cluster
}

int styrene_monomer(MONOMER &m) { //-CH-[C6H5]CH2-
	// 3 cluster
	m.reset(3, 2); // 2 relationship
	m.setHinges(2); // 2 free hinge
	int na = 0;
	na = CH(m.cluster[0]); // cluster 0 -- -CH--
	na = Benzene1(m.cluster[1]); // cluster 1 -- C6H5 : -C6H5
	na = CH2(m.cluster[2]); // cluster 2 -- -CH2-

	m.cluster[1].rotate(m.cluster[1].hinge[0].vHinge, -15, false);
	// cluster 0 has free hinge # 0
	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster, 1, m.cluster + 1, 0, SINGLE_BOUND);
	m.rship[0].set_relationship(0, 1, 1, 0, SINGLE_BOUND);

	m.cluster[2].rotate(m.cluster[2].hinge[0].vHinge, 60, false);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster, 2, m.cluster + 2, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(0, 2, 2, 0, SINGLE_BOUND);
	// cluster 2 has free hinge # 1

	m.nc[0] = 0; m.iHinge[0] = 0;
	m.nc[1] = 2; m.iHinge[1] = 1;

	m.nbase = 0;

	return 3; // the base cluster
}

int styrene1_monomer(MONOMER &m) { //-CH2-[C6H5]CH2-
	// 3 cluster
	m.reset(3, 2); // 2 relationship
	int na = 0;
	na = CH2(m.cluster[0]); // cluster 0 -- -CH2-
	na = Benzene1(m.cluster[1]); // cluster 1 -- C6H5 : -C6H5
	na = CH2(m.cluster[2]); // cluster 2 -- -CH2-

	m.setHinges(1); // 1 free hinge

	//m.cluster[1].rotate(m.cluster[1].hinge[0].vHinge, -30, false);
	
	m.cluster[0].rotate(m.cluster[0].hinge[0].vHinge, 120, false);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 1, 0, m.cluster, 0, SINGLE_BOUND);
	m.rship[0].set_relationship(1, 0, 0, 0, SINGLE_BOUND);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster, 1, m.cluster + 2, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(0, 1, 2, 0, SINGLE_BOUND);
	// cluster 2 has a free hinge # 1

	m.nc[0] = 2;
	m.iHinge[0] = 1;

	m.nbase = 0;

	return 3; // the base cluster
}

int styrene2_monomer(MONOMER &m) { //-CH-[C6H5]CH3
	// 3 cluster
	m.reset(3, 2); // 2 relationship
	int na = 0;
	na = Benzene1(m.cluster[0]); // cluster 0 -- C6H5 : -C6H5
	na = CH3(m.cluster[1]); // cluster 1 -- -CH3
	na = CH(m.cluster[2]); // cluster 2 -- -CH--

	m.setHinges(1); // 1 free hinge

	m.cluster[0].rotate(m.cluster[0].hinge[0].vHinge, -30, false);

	// cluster 2 has a free hinge #0
	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 2, 1, m.cluster, 0, SINGLE_BOUND);
	m.rship[0].set_relationship(2, 1, 0, 0, SINGLE_BOUND);

	// no parent-child relation is required, 2nd cluster is single
	ConnectClusters(m.cluster + 2, 2, m.cluster + 1, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(2, 2, 1, 0, SINGLE_BOUND);

	m.nc[0] = 2;
	m.iHinge[0] = 0;

	m.nbase = 0;

	return 3; // the base cluster
}

int ethylene_oxide_monomer(MONOMER &m) { //-CH2-CH2-O-
	// 3 cluster
	m.reset(3, 2); // 2 relationship
	int na = 0;
	na = CH2_template(m.cluster[0], PEO_C, PEO_H, CH1_BLENGTH, CINDX_PEO_CH2); // cluster 0 -- -CH2-
	na = CH2_template(m.cluster[1], PEO_C, PEO_H, CH1_BLENGTH, CINDX_PEO_CH2); // cluster 1 -- -CH2-
	na = DoubleBoundAtom(m.cluster[2], PEO_O, CO_BLENGTH, 112, CINDX_PEO_O); // cluster 2 -- -O-

	m.setHinges(2); // 2 free hinge

	m.cluster[1].rotate(m.cluster[1].hinge[0].vHinge, 180, false);
	// cluster 0 has free hinge # 0
	ConnectClusters(m.cluster, 1, m.cluster + 1, 0, SINGLE_BOUND);
	m.rship[0].set_relationship(0, 1, 1, 0, SINGLE_BOUND);

	m.cluster[2].rotate(m.cluster[2].hinge[0].vHinge, 60, false);

	ConnectClusters(m.cluster + 1, 1, m.cluster + 2, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(1, 1, 2, 0, SINGLE_BOUND);
	// cluster 2 has free hinge # 1

	m.nc[0] = 0; m.iHinge[0] = 0;
	m.nc[1] = 2; m.iHinge[1] = 1;

	m.nbase = 0;

	return 3; // number of cluster
}

int ethylene_oxide1_monomer(MONOMER &m) { //-CH2-CH2-OH
	// 3 cluster
	m.reset(3, 2); // 2 relationship
	int na = 0;
	na = CH2_template(m.cluster[0], PEO_C, PEO_H, CH1_BLENGTH, CINDX_PEO_CH2); // cluster 0 -- -CH2-
	na = CH2_template(m.cluster[1], PEO_C, PEO_H, CH1_BLENGTH, CINDX_PEO_CH2); // cluster 1 -- -CH2-
	na = OH_template(m.cluster[2], PEO_O, SB_H, CINDX_PEO_OH); // cluster 2 -- -OH
	// SB_H does not have any charge

	m.setHinges(1); // 1 free hinge

	m.cluster[1].rotate(m.cluster[1].hinge[0].vHinge, 180, false);
	// cluster 0 has free hinge # 0
	ConnectClusters(m.cluster, 1, m.cluster + 1, 0, SINGLE_BOUND);
	m.rship[0].set_relationship(0, 1, 1, 0, SINGLE_BOUND);

	m.cluster[2].rotate(m.cluster[2].hinge[0].vHinge, 60, false);

	ConnectClusters(m.cluster + 1, 1, m.cluster + 2, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(1, 1, 2, 0, SINGLE_BOUND);
	// cluster 2 has no hinge

	m.nc[0] = 0; m.iHinge[0] = 0;

	m.nbase = 0;

	return 3; // number of cluster
}

int propylene_oxide_monomer(MONOMER &m) { //-CH2-CH(CH3)-O-
	// 4 cluster
	m.reset(4, 3); // 3 relationship
	int na = 0;
	na = CH2(m.cluster[0]); // cluster 0 -- -CH2-
	na = CH(m.cluster[1]); // cluster 1 -- -CH--
	na = CH3(m.cluster[2]); // cluster 2 -- -CH3
	na = DoubleBoundAtom(m.cluster[3], NEUTRAL_O, CO_BLENGTH, 112, CINDX_PPO_O); // cluster 2 -- -O-
	// NEUTRAL_O has no charge

	m.setHinges(2); // 2 free hinge

	m.cluster[1].rotate(m.cluster[1].hinge[0].vHinge, 60, false);
	// cluster 0 has free hinge # 0
	ConnectClusters(m.cluster, 1, m.cluster + 1, 0, SINGLE_BOUND);
	m.rship[0].set_relationship(0, 1, 1, 0, SINGLE_BOUND);

	m.cluster[2].rotate(m.cluster[2].hinge[0].vHinge, 60, false);
	ConnectClusters(m.cluster + 1, 1, m.cluster + 2, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(1, 1, 2, 0, SINGLE_BOUND);

	m.cluster[3].rotate(m.cluster[3].hinge[0].vHinge, 60, false);
	ConnectClusters(m.cluster + 1, 2, m.cluster + 3, 0, SINGLE_BOUND);
	m.rship[2].set_relationship(1, 2, 3, 0, SINGLE_BOUND);

	// cluster 3 has free hinge # 1

	m.nc[0] = 0; m.iHinge[0] = 0;
	m.nc[1] = 3; m.iHinge[1] = 1;

	m.nbase = 0;

	return 4; // number of cluster
}

int propylene_oxide1_monomer(MONOMER &m) { //-CH2-CH(CH3)-OH
	// 4 cluster
	m.reset(4, 3); // 3 relationship
	int na = 0;
	na = CH2(m.cluster[0]); // cluster 0 -- -CH2-
	na = CH(m.cluster[1]); // cluster 1 -- -CH--
	na = CH3(m.cluster[2]); // cluster 2 -- -CH3
	na = OH_template(m.cluster[3], NEUTRAL_O, SB_H, CINDX_PPO_OH); // cluster 2 -- -OH
	// SB_H does not have any charge, NEUTRAL_O has no charge

	m.setHinges(1); // 1 free hinge

	m.cluster[1].rotate(m.cluster[1].hinge[0].vHinge, -15, false);
	// cluster 0 has free hinge # 0
	ConnectClusters(m.cluster, 1, m.cluster + 1, 0, SINGLE_BOUND);
	m.rship[0].set_relationship(0, 1, 1, 0, SINGLE_BOUND);

	m.cluster[2].rotate(m.cluster[2].hinge[0].vHinge, 60, false);
	ConnectClusters(m.cluster + 1, 1, m.cluster + 2, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(1, 1, 2, 0, SINGLE_BOUND);

	m.cluster[3].rotate(m.cluster[3].hinge[0].vHinge, 60, false);
	ConnectClusters(m.cluster + 1, 2, m.cluster + 3, 0, SINGLE_BOUND);
	m.rship[2].set_relationship(1, 2, 3, 0, SINGLE_BOUND);

	m.nc[0] = 0; m.iHinge[0] = 0;

	m.nbase = 0;

	return 4; // number of cluster
}

int propylene_oxide2_monomer(MONOMER &m) { //CH3-CH(CH3)-O-
	// 4 cluster
	m.reset(4, 3); // 3 relationship
	int na = 0;
	na = CH3(m.cluster[0]); // cluster 0 -- CH3-
	na = CH(m.cluster[1]); // cluster 1 -- -CH--
	na = CH3(m.cluster[2]); // cluster 2 -- -CH3
	na = DoubleBoundAtom(m.cluster[3], NEUTRAL_O, CO_BLENGTH, 112, CINDX_PPO_O); // cluster 2 -- -O-
	// NEUTRAL_O has no charge

	m.setHinges(1); // 1 free hinge

	m.cluster[1].rotate(m.cluster[1].hinge[0].vHinge, -15, false);
	// cluster 0 has free hinge # 0
	ConnectClusters(m.cluster, 0, m.cluster + 1, 0, SINGLE_BOUND);
	m.rship[0].set_relationship(0, 0, 1, 0, SINGLE_BOUND);

	m.cluster[2].rotate(m.cluster[2].hinge[0].vHinge, 60, false);
	ConnectClusters(m.cluster + 1, 1, m.cluster + 2, 0, SINGLE_BOUND);
	m.rship[1].set_relationship(1, 1, 2, 0, SINGLE_BOUND);

	m.cluster[3].rotate(m.cluster[3].hinge[0].vHinge, 60, false);
	ConnectClusters(m.cluster + 1, 2, m.cluster + 3, 0, SINGLE_BOUND);
	m.rship[2].set_relationship(1, 2, 3, 0, SINGLE_BOUND);

	m.nc[0] = 0; m.iHinge[0] = 1;

	m.nbase = 0;

	return 4; // number of cluster
}

bool construct_defined_simple_molecule(char *mol, BASIC_CLUSTER* p) {
	if (p == NULL) return false;
	if (strcmp(mol, "BENZENE") == 0) Benzene(*p);
	else if (strcmp(mol, "SPC_WATER") == 0) SPCwater(*p);
	else if (strcmp(mol, "SPCE_WATER") == 0) SPCEwater(*p);
	else return false;
	return true;
}

bool construct_defined_cluster(char *unit, CLUSTER* p) {
	if (p == NULL) return false;
	if (strcmp(unit, "CH3") == 0) CH3(*p);
	else if (strcmp(unit, "CH2") == 0) CH2(*p);
	else if (strcmp(unit, "C2H2") == 0) C2H2(*p);
	else if (strcmp(unit, "CH") == 0) CH(*p);
	else if (strcmp(unit, "C2H") == 0) C2H(*p);
	else if (strcmp(unit, "BENZENE1") == 0) Benzene1(*p);
	else if (strcmp(unit, "BENZENE12") == 0) Benzene12(*p);
	else if (strcmp(unit, "BENZENE13") == 0) Benzene13(*p);
	else if (strcmp(unit, "BENZENE14") == 0) Benzene14(*p);
	else if (strcmp(unit, "SPC_O") == 0) SPC_O(*p);
	else if (strcmp(unit, "SPC_OH") == 0) SPC_OH(*p);
	else return false;
	return true;
}

int construct_defined_monomer(char *monomer, MONOMER &m) {
	if (strcmp(monomer, "STYRENE") == 0) return styrene_monomer(m);
	else if (strcmp(monomer, "STYRENE1") == 0) return styrene1_monomer(m);
	else if (strcmp(monomer, "STYRENE2") == 0) return styrene2_monomer(m);
	if (strcmp(monomer, "ISOPRENE") == 0) return isoprene_monomer(m);
	else if (strcmp(monomer, "ISOPRENE1") == 0) return isoprene1_monomer(m);
	else if (strcmp(monomer, "ISOPRENE2") == 0) return isoprene2_monomer(m);
	else if (strcmp(monomer, "OXYETHYLENE") == 0 || strcmp(monomer, "EO") == 0) return ethylene_oxide_monomer(m);
	else if (strcmp(monomer, "OXYETHYLENE1") == 0 || strcmp(monomer, "EO1") == 0) return ethylene_oxide1_monomer(m);
	else if (strcmp(monomer, "OXYPROPYLENE") == 0 || strcmp(monomer, "PO") == 0) return propylene_oxide_monomer(m);
	else if (strcmp(monomer, "OXYPROPYLENE1") == 0 || strcmp(monomer, "PO1") == 0) return propylene_oxide1_monomer(m);
	else if (strcmp(monomer, "OXYPROPYLENE2") == 0 || strcmp(monomer, "PO2") == 0) return propylene_oxide2_monomer(m);
	else return 0;
}

*/
// twist the torsion angle of cluster along its parent hinge
void tweak(MMOLECULE &mm, bool bsave) {
	char buffer[250] = "\0";
	int nCluster = 0;
	float omega = 0;
	while (1) {
		sprintf(buffer, "[cluster < %d, omega] : ", mm.nCluster);
		show_msg(buffer, false);
		input(buffer, 250);
		if (sscanf(buffer, "%d %f", &nCluster, &omega) != 2) break;
		else {
			mm.twTA(nCluster, omega, false); // in degree
			if (bsave) {
				{
					char fname[200] = "\0";
					sprintf(fname, "%s.xyz", mm.mol);
					ofstream out;
					out.open(fname);
					if (!out.is_open()) {
						sprintf(errmsg, "can not open file %s to save the cordinates", fname);
						show_msg(buffer);
						return;
					}
					mm.xyz_save(out);
					out.close();
				}
			}
		}
	}
}

// copy the clusters in monomer to cluster array;
// and setup the same relationship between the clusters in monomer
bool copy_monomer(CLUSTER* cluster, MONOMER* monomer) {
	int i, nc;
	// copy the clusters
	for (nc = 0; nc < monomer->nCluster; nc++) {
		cluster[nc].setAtoms(monomer->cluster[nc].nAtoms);
		if (!cp(cluster + nc, monomer->cluster + nc)) return false;
	}
	// setup connections
	for (i = 0; i < monomer->nrship; i++) {
		if (!ConnectClusters(cluster + monomer->rship[i].n1, monomer->rship[i].iHinge1, 
			cluster + monomer->rship[i].n2, monomer->rship[i].iHinge2, monomer->rship[i].bound_type)) return false;
	}
	return true;
}

bool get_period(char **mol, char *period) {
	char *str = NULL;
	char *start = NULL, *end = NULL, *pch = NULL;
	char buff[50] = "\0";
	int i = 0;

	strcpy(period, "\0");
	if (strlen(*mol) == 0) return false;

	str = *mol;
	bool bracket = false, status = true;

	start = NULL; end = NULL;
	if (!pass(&str, ' ')) {*mol = str; return false;}
	if (!pass(&str, '-')) {*mol = str; return false;}
	bracket = (*str == '[' ? true : false);
	if (bracket) if (!pass(&str, '[')) {*mol = str; return false;}
	start = str;
	if (bracket) {
		status = NextSection(&str, ']');
		*mol = str; 
		if (!status) str--;
	}
	else {
		status = NextSection(&str, '[');
		if (status) *mol = str;
		else {
			str = start;
			status = NextSection(&str, '-');
			if (status) *mol = str;
			else {str--; *mol = NULL;}
		}
	}
	end = str;

	if (start == end || start == NULL) {*mol = NULL; return false;}
	pch = start; i = 0;
	while (pch != end) {period[i] = *pch; pch++; i++;}
	period[i] = *end; period[i + 1] = '\0';
	for (i = strlen(period) - 1; i >= 0; i--) {
		if (period[i] == ' ' || period[i] == '-' || period[i] == ']' || period[i] == '[')
			period[i] = '\0';
		else break;
	}
	return true;
}

bool get_period_number(char **str, int &n) {
	if (*str == NULL) return false;
	if (!pass(str, ']')) {n = 1; return true;}
	char buffer[250] = "\0";
	int i = 0, nunit = 0;
	char *pch = *str;
	while (strlen(pch) > 0 && *pch != '-' && *pch != '\0' && *pch != '[' && *pch != ']') {
		buffer[i] = *pch;
		pch++; i++;
	}
	buffer[i] = '\0';
	*str = pch;
	if (strlen(buffer) > 0) {
		if (sscanf(buffer, "%d", &nunit) == 1) {n = nunit; return true;}
		else {n = 1; return false;}
	}
	else {n = 1; return true;}
}

bool get_cluster_unit(char **str, char *unit) {
	char *start = NULL, *end = NULL, *pch = NULL;
	int i = 0;
	bool status;
	strcpy(unit, "\0");
	if (strlen(*str) == 0) return false;

	if (!pass(str, ' ')) return false;
	if (!pass(str, '-')) return false;
	start = *str;
	status = NextSection(str, '-');
	end = (status ? *(str--) : *str);

	if (!status) {
		if (sscanf(start, "%s", unit) != 1) return false;
	}
	else {
		pch = start; i = 0;
		while (pch != end) {unit[i] = *pch; pch++; i++;}
		unit[i] = *end; unit[i+1] = '\0';
	}
	kick_off_tail(unit, ' ');
	kick_off_tail(unit, '-');
	return true;
}

bool construct_monomer(char *monomer_name, MONOMER *m) {
	strcpy(m->name, monomer_name);

	char *name = monomer_name;
	char cunit[20] = "\0";
	bool status = false;
	int ncluster = 0, n = 0;
	char *pch = NULL;

	pch = monomer_name; ncluster = 0;
	while (get_cluster_unit(&pch, cunit) && strlen(cunit) > 0) ncluster++;

	if (ncluster == 0) return false;
	if (ncluster == 1) {
		pch = monomer_name;
		get_cluster_unit(&pch, cunit);
		if (defined_monomer(cunit)) {
			if (construct_defined_monomer(cunit, *m) >= 0) return true;
			else {m->reset(0, 0); return false;}
		}
		else {
			m->reset(1, 0);
			if (!construct_defined_cluster(cunit, m->cluster)) {
				m->reset(0, 0); return false;
			}
			else {
				m->nbase = 0;
				m->setHinges(m->cluster->nHinges); // we asumme that this is cluster-chain
				for (n = 0; n < m->cluster->nHinges; n++) {
					m->nc[n] = 0; // only one cluster in the monomer
					m->iHinge[n] = n;
				}
				m->uid = m->cluster->cID; // monomer uid is the cluster id
				return true;
			}
		}
	}
	else {
		// in this case, we have to asumme the cluster are linearized, one parent and one child only
		pch = monomer_name;
		m->reset(ncluster, ncluster - 1);
		n = 0;
		while (get_cluster_unit(&pch, cunit) && strlen(cunit) > 0) {
			if (!construct_defined_cluster(cunit, m->cluster + n)) {m->reset(0, 0); m->setHinges(0); return false;}
			n++;
		}
		for (n = 1; n < ncluster; n++) {
			ConnectClusters(m->cluster + n - 1, 1, m->cluster + n, 0, SINGLE_BOUND);
			m->rship[n-1].set_relationship(n-1, 1, n, 0, SINGLE_BOUND);
		}
		m->nc[0] = 0;
		m->iHinge[0] = 0;
		m->nc[1] = ncluster - 1;
		m->iHinge[1] = 1;
		m->nbase = 0;
		return true;
	}
}

int nMonomerFromPolymerName(char *polymer) {
	char *str = polymer;
	char buffer[256] = "\0";
	int n = 0, nblock = 0;
	while (get_period(&str, buffer) && strlen(buffer) > 0) {
		n++; if (str == NULL || strlen(str) == 0) break;
		get_period_number(&str, nblock);
		if (str == NULL || strlen(str) == 0) break;
	}
	return n;
}

bool ParseSingleChainPolymer(char *polymer, MONOMER **m, int &n_monomer, int **nperiod) {
	char buffer[256] = "\0";
	char mhead[250] = "\0", mtail[250] = "\0", mperiod[250] = "\0";
	int nblock = 0, n = 0, nbase = 0;

	char *str = polymer;
	while (get_period(&str, buffer) && strlen(buffer) > 0) {
		n++; //if (str == NULL || strlen(str) == 0) break;
		get_period_number(&str, nblock);
		if (str == NULL || strlen(str) == 0) break;
	}
	n_monomer = n; nblock = n;
	if (nblock == 0) return false;
	if (*m != NULL) delete[] *m;
	*m = new MONOMER[nblock];
	if (*nperiod != NULL) {delete[] *nperiod; *nperiod = NULL;}
	if (nblock > 2) *nperiod = new int[nblock - 2];

	str = polymer;
	get_period(&str, mhead);
	if (!construct_monomer(mhead, *m)) { // this is head
		delete[] *m; *m = NULL; 
		if (*nperiod != NULL) {delete[] *nperiod; *nperiod = NULL;}
		return false;
	}

	n = 0;
	while (n < nblock - 2) {
		if (get_period(&str, mperiod) && get_period_number(&str, (*nperiod)[n])) {
			// here true will not affect the monomer because it is not head or tail
			if (!construct_monomer(mperiod, (*m) + n + 1)) {
				delete[] *m; *m = NULL; 
				if (*nperiod != NULL) {delete[] *nperiod; *nperiod = NULL;}
				return false;
			}
			n++;
		}
	}

	if (!get_period(&str, mtail) || !construct_monomer(mtail, *m + nblock - 1)) { // this is tail
		delete[] *m; *m = NULL; 
		if (*nperiod != NULL) {delete[] *nperiod; *nperiod = NULL;}
		return false;
	}

	return true;
}

// unit is composed with cluster c in order, start with head and end with tail
// using C-C as hing bound to connect the clusters
void SingleChainPolymer(MMOLECULE &mm, MONOMER *head, MONOMER *m, int n_monomer, int *nperiod, float *ta, MONOMER *tail, float ta_tail, int &nbase, int nCG) {
	float ta_polymer = 0;
	float bl = 1.54f; // Op & Om must be C, C-C distance
	int i = 0, j = 0, k = 0, n = 1;

	// deuterize the monomer
	for (n = 0; n < n_monomer; n++) monomer_deuterization(m + n);

	int nCluster = head->nCluster + tail->nCluster;
	for (n = 0; n < n_monomer; n++) nCluster += m[n].nCluster * nperiod[n];
	CLUSTER *cluster = NULL;
	int nChild = 0, iChildHinge = 0, n2Parent = 0, iParentHinge = 0;
	int ncbase = 0, i_monomer = 0;
	int indx_monomer = 0, nCG_monomer = nCG, iCG_monomer = 0;

	int nMonomer = 2; // monomer in the macromolecule, head and tail
	int indx_Monomer = 0;
	for (n = 0; n < n_monomer; n++) nMonomer += nperiod[n];

	mm.reset();
	mm.nCluster = nCluster;
	mm.cluster = new CLUSTER[nCluster];

	mm.monomer = new _GROUP<CLUSTER>[nMonomer]; mm.nMonomer = nMonomer;

	n = 0;
	cluster = mm.cluster + n;

	copy_monomer(cluster, head); // copy the monomer in head to molecule
	mm.monomer[0].set(head->nCluster);
	for (k = 0; k < head->nCluster; k++) {
		cluster[k].cmIndx = n + k;
		cluster[k].cgIndx = 0; //labeling the cgIndx
		mm.monomer[0].u[k] = cluster + k; // setup monomer
		cluster[k].monIndx = 0;
	}
	mm.monomer[0].uid = head->uid;
	// setup parent-child relationship inside the monomer
	setup_parent_child_relationship(cluster + head->nbase); 
	// the cluster index to connect to next monomer
	// head monomer has one free child and we asumme this is single-chain polymer
	nChild = n + head->nc[0]; iChildHinge = head->iHinge[0];
	cluster = cluster + head->nCluster; // cluster for next monomer
	n += head->nCluster; // cluster index for next monomer

	i_monomer = 1;
	if (i_monomer < nbase) ncbase += head->nCluster;

	iCG_monomer = 0;
	for (i = 0; i < n_monomer; i++) {
		for (j = 0; j < nperiod[i]; j++) {
			indx_Monomer++; // a new monomer in macromolecule
			copy_monomer(cluster, m + i); // copy the monomer in head to molecule
			mm.monomer[indx_Monomer].set(m[i].nCluster);
			mm.monomer[indx_Monomer].uid = m[i].uid;
			// labeling the cluster's cgIndx
			for (k = 0; k < m[i].nCluster; k++) {
				cluster[k].cmIndx = n + k;
				cluster[k].cgIndx = indx_monomer;
				mm.monomer[indx_Monomer].u[k] = cluster + k; // setup monomer
				cluster[k].monIndx = indx_Monomer;
			}
			iCG_monomer++; if (iCG_monomer == nCG_monomer) {indx_monomer++; iCG_monomer = 0;}
			// very important
			// setup parent-child relationship from the parent of the cluster
			cluster[m[i].nc[0]].set_parent_hinge(m[i].iHinge[0]);
			setup_parent_child_relationship(cluster + m[i].nc[0]); 
			// connect the monomer to its parent
			ConnectClusters(mm.cluster + nChild, iChildHinge, cluster + m[i].nc[0], m[i].iHinge[0], SINGLE_BOUND);
			// tweak the monomer with some torsion angle
			ta_polymer = ta[i] * Angle2Radian - cluster[m[i].nc[0]].km.ta;
			mm.twTA(n + m[i].nc[0], ta_polymer, true);
			cluster[m[i].nc[0]].km.ta = ta[i] * Angle2Radian;
			PERIOD_RANGE(cluster[m[i].nc[0]].km.ta, -PI, PI, PI2)
			// the cluster index to connect to next monomer
			nChild = n + m[i].nc[1]; iChildHinge = m[i].iHinge[1];
			cluster = cluster + m[i].nCluster; // cluster for next monomer
			n += m[i].nCluster; // cluster index for next monomer

			i_monomer += 1;
			if (i_monomer < nbase) ncbase += m[i].nCluster;
		}
	}
	// tail monomer
	indx_Monomer++; // a new monomer in macromolecule
	copy_monomer(cluster, tail); // copy the monomer in head to molecule
	mm.monomer[indx_Monomer].set(tail->nCluster);
	mm.monomer[indx_Monomer].uid = tail->uid;
	if (iCG_monomer == 0) indx_monomer -= 1; // tail is attached to the last cluster
	for (k = 0; k < tail->nCluster; k++) {
		cluster[k].cmIndx = n + k;
		cluster[k].cgIndx = indx_monomer; // labeling the cgIndx 
		mm.monomer[indx_Monomer].u[k] = cluster + k; // setup monomer
		cluster[k].monIndx = indx_Monomer;
	}
	cluster[tail->nc[0]].set_parent_hinge(tail->iHinge[0]);
	// very important
	// setup parent-child relationship from the parent of the cluster
	setup_parent_child_relationship(cluster + tail->nc[0]);
	// connect the monomer to its parent
	ConnectClusters(mm.cluster + nChild, iChildHinge, cluster + tail->nc[0], tail->iHinge[0], SINGLE_BOUND);
	// tweak the monomer with some torsion angle
	ta_polymer = ta_tail * Angle2Radian - cluster[tail->nc[0]].km.ta;
	mm.twTA(n + tail->nc[0], ta_polymer, true);
	cluster[tail->nc[0]].km.ta = ta_tail * Angle2Radian;
	PERIOD_RANGE(cluster[tail->nc[0]].km.ta, -PI, PI, PI2)

	i_monomer += 1;
	if (i_monomer < nbase) ncbase += tail->nCluster;

	nbase = ncbase;
}

bool ConstructSingleChainPolymerMonomerChain(MMOLECULE &mm, char *polymer, float* ta, float ta_tail, int &nbase, int nCG) {
	mm.reset();
	MONOMER *m = NULL;
	int n_monomer = 0, nblock = 0, *nperiod = NULL;
	if (!ParseSingleChainPolymer(polymer, &m, n_monomer, &nperiod)) {
		sprintf(errmsg, "polymer %s can not be recognized");
		if (m != NULL) {delete[] m; m = NULL;}
		if (nperiod != NULL) {delete[] nperiod; nperiod = NULL;}
		return false;
	}
	nblock = n_monomer - 2; // kick off head and tail
	SingleChainPolymer(mm, m, m + 1, nblock, nperiod, ta, m + n_monomer - 1, ta_tail, nbase, nCG);
	if (m != NULL) {delete[] m; m = NULL;}
	if (nperiod != NULL) {delete[] nperiod; nperiod = NULL;}

	return true;
}

bool ConstructSingleChainPolymerFromParamFile(MMOLECULE &mm, char *fname, bool save) {
	char buffer[256] = "\0", mol[250] = "\0", msg[256] = "\0";
	int struct_model = -1;
	if (!read_msg(fname, buffer, "[MOL_STRUCT_STYLE]")) struct_model = 0;
	if (sscanf(buffer, "%d", &struct_model) != 1) {
		sprintf(msg, "Molecular structural definition style is not given, assuming the structure is defined as single-chain structure");
		show_msg(msg); struct_model = 0;
	}
	else if (struct_model == 1) {return ReadMacroMol(fname, mm, save);}
	else {
		sprintf(msg, "Unknown Molecular structural definition style: %d", struct_model); show_msg(msg);
		sprintf(msg, "  style 0 -- single-chain structure"); show_msg(msg);
		sprintf(msg, "  style 1 -- structure with monomer(group) given in file, and given links"); show_msg(msg);
		return false;
	}

	float *ta = NULL;
	int nta = 0;

	if (!read_msg(fname, mol, "[MOLECULE]")) return false;
	strcpy(mm.mol, mol);
	if (!read_msg(fname, buffer, "[POLYMER_INIT_TORSION_ANGLE]", true)) return false;
	if ((nta = read_float(buffer, &ta)) == 0) {
		sprintf(errmsg, "can not get torsion angle for monomers from %s", buffer);
		return false;
	}

	int nbase = 0;
	if (!read_msg(fname, buffer, "[BASE_CLUSTER]") || sscanf(buffer, "%d", &nbase) != 1 ||
		nbase < 0) nbase = 0;

	int nCG = 1;
	if (!read_msg(fname, buffer, "[COARSE_GRAIN]") || sscanf(buffer, "%d", &nCG) != 1 || nCG <= 0)
		nCG = 1;

	bool status = ConstructSingleChainPolymerMonomerChain(mm, mol, ta, ta[nta-1], nbase, nCG);
	if (ta != NULL) {delete[] ta; ta = NULL;}
	if (!status) return false;

	if (nbase >= mm.nCluster) nbase = mm.nCluster - 1;

	mm.base = mm.cluster + nbase;
	mm.base->disable_parent();
	setup_parent_child_relationship(mm.base, true); // report error when free hinge exist
	mm.checkMolecule(true);

	if (!check_parent_child_setup(mm.cluster, mm.nCluster)) return false;

	if (!mm.linkPars(atompar_db, errmsg)) {show_infor(errmsg); return false;}

	mm.calTorsionAngle();

	mm.setup_cluster_matom(); // require the atomic mass, and the parent-child relationship

	// set mass center to {0, 0, 0}
	VECTOR<3> cm;
	mm.SetMassCenter(cm.v, true, true);

	char template_fname[256] = "\0";
	if (save) {
		if (read_msg(fname, template_fname, "[MOL_TEMPLATE_FILENAME]")) {
			::save_xyz_struct(mm, template_fname);
		}
	}
	return true;
}

bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol) {
	sm->reset();
	if (!defined_simple_molecule(mol)) return false;
	sm->c = new BASIC_CLUSTER;
	if (!construct_defined_simple_molecule(mol, sm->c)) {sm->reset(); return false;}
	sm->c->mType = 1; // simple molecule
	sm->xaxis.v[0] = 1; sm->xaxis.v[1] = 0; sm->xaxis.v[2] = 0;
	sm->zaxis.v[0] = 0; sm->zaxis.v[1] = 0; sm->zaxis.v[2] = 1;
	if (!sm->c->linkPars(atompar_db, errmsg)) {show_infor(errmsg); return false;}
	sm->check_atom_cordinate();
	//sm->c->vp = &(sm->r); sm->c->Op = &(sm->r);  // done in check_atom_cordinate()
	strcpy(sm->mol, mol);
	sm->c->setup_cluster_matom();
	return true;
}

bool ConstructPointMolecule(PMOLECULE* pm, char *mol) {
	pm->reset();
	if (!defined_point_molecule(mol)) return false;
	pm->c = new BASIC_CLUSTER;
	if (!construct_defined_point_molecule(mol, pm->c)) {pm->reset(); return false;}
	pm->c->mType = 2; // point molecule
	if (!pm->c->linkPars(atompar_db, errmsg)) {show_infor(errmsg); return false;}
	pm->c->M = pm->c->atom->par->m;
	pm->q = pm->c->atom->c;
	pm->check_coordinate();
	strcpy(pm->mol, mol);
	pm->c->setup_cluster_matom();
	return true;
}

/*************************************************************************************/
/**           define a structure with multi-monomers in a file                       */
/*************************************************************************************/
bool ReadMonomersFromFile(char *fname, char *mname, MONOMER_GROUP &mgroup) {
	char title[256] = "\0", buffer[512] = "\0", msg[256] = "\0";
	strcpy(mgroup.name, mname);
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to construct %s from %s: failure to open the file", mname, fname); show_msg(msg);
		in.close(); return false;
	}
	
	int nm = 0, nlink = 0, nHinges = 0, i = 0, indx = 0, m1 = 0, iHinge1 = 0, m2 = 0, iHinge2 = 0;
	sprintf(title, "[%s]", mname);
	if (!search(in, title, buffer)) {
		sprintf(msg, "%s is not defined in file %s", title, fname); show_msg(msg);
		in.close(); return false;
	}
	if (!search(in, "MONOMER", buffer)) {
		sprintf(msg, "MONOMER is not defined in file %s", fname); show_msg(msg);
		in.close(); return false;
	}
	if (sscanf(buffer, "%s %d", title, &nm) != 2 || nm < 1) {
		sprintf(msg, "see %s -- number of monomers [>0] is not given in file %s", buffer, fname); show_msg(msg);
		in.close(); return false;
	}
	mgroup.monomer.SetArray(nm);
	i = 0;
	while (i < nm) {
		if (in.eof()) {
			sprintf(msg, "failure to get enough nomomers from file %s, get %d expect %d", fname, i, nm); show_msg(msg);
			mgroup.reset(); in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %s", &indx, title) != 2) {
			sprintf(msg, "UNKNOWN format of monomer in %s: %s  -- right format [indx, monomer name]", fname, buffer); show_msg(msg);
			mgroup.reset(); in.close(); return false;
		}
		if (indx < 0 || indx >= nm) continue;
		if (!construct_monomer(title, mgroup.monomer.m + indx)) {
			sprintf(msg, "unknown monomer %d %s in %s", indx, title, fname); show_msg(msg);
			mgroup.reset(); in.close(); return false;
		}
		i++;
	}

	if (!search(in, "LINK", buffer)) {
		sprintf(msg, "LINK is not defined in file %s", fname); show_msg(msg);
		mgroup.mLink.release(); in.close(); return true;
	}
	if (sscanf(buffer, "%s %d", title, &nlink) != 2) {
		sprintf(msg, "see %s -- assuming that number of link is 0", buffer, fname); show_msg(msg);
		mgroup.mLink.release(); in.close(); return true;
	}
	if (nlink < 1) {mgroup.mLink.release(); in.close(); return true;}
	else mgroup.mLink.SetArray(nlink);
	i = 0;
	while (i < nlink) {
		if (in.eof()) {
			sprintf(msg, "failure to get enough link-relationship from file %s, get %d expect %d", fname, i, nlink); show_msg(msg);
			mgroup.reset(); in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d %d %d %d", &indx, &m1, &iHinge1, &m2, &iHinge2) != 5) {
			sprintf(msg, "UNKNOWN format of link-relationship in %s: %s ", fname, buffer); show_msg(msg);
			show_msg(" -- right format: [indx] [monomer indx1] [hinge indx1] [monomer indx2] [hinge indx2]");
			mgroup.reset(); in.close(); return false;
		}
		if (indx < 0 || indx >= nlink) continue;
		mgroup.mLink.m[indx].mIndx1 = m1; mgroup.mLink.m[indx].iHinge1 = iHinge1;
		mgroup.mLink.m[indx].mIndx2 = m2; mgroup.mLink.m[indx].iHinge2 = iHinge2;
		i++;
	}

	if (!search(in, "FREE-HINGES", buffer)) {
		sprintf(msg, "FREE-HINGES is not defined in file %s", fname); show_msg(msg);
		mgroup.mHinge.release(); in.close(); return true;
	}
	if (sscanf(buffer, "%s %d", title, &nHinges) != 2) {
		sprintf(msg, "see %s -- assuming that number of free-hinge is 0", buffer, fname); show_msg(msg);
		mgroup.mHinge.release(); in.close(); return true;
	}
	if (nHinges < 1) {mgroup.mHinge.release(); in.close(); return true;}
	mgroup.mHinge.SetArray(nHinges);
	i = 0;
	while (i < nHinges) {
		if (in.eof()) {
			sprintf(msg, "failure to get enough hinges from file %s, get %d expect %d", fname, i, nHinges); show_msg(msg);
			mgroup.reset(); in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d %d ", &indx, &m1, &iHinge1) != 3) {
			sprintf(msg, "UNKNOWN format of HINGE [mIndx, hIndx] in %s: %s ", fname, buffer); show_msg(msg);
			mgroup.reset(); in.close(); return false;
		}
		if (indx < 0 || indx >= nHinges) continue;
		mgroup.mHinge.m[indx].mIndx = m1; mgroup.mHinge.m[indx].iHinge = iHinge1;
		i++;
	}
	in.close(); return true;
}

bool ReadMacroMol(char *fname, MMOLECULE &mm, bool save) {
	char buffer[256] = "\0", msg[256] = "\0", title[256] = "\0", mname[256] = "\0", sfname[256] = "\0";
	int n_monomer = 0, nlink = 0, i, indx, npars = 0, nc1 = 0, iHinge1 = 0, nc2 = 0, iHinge2 = 0;
	ifstream in; in.open(fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open file %s", fname); show_msg(msg); return false;
	}
	if (!search(in, "[COMPONENT]", buffer)) {
		sprintf(msg, "[COMPONENT] is not given in file %s", fname); show_msg(msg); in.close(); return false;
	}
	if (sscanf(buffer, "%s %d", title, &n_monomer) != 2) {
		sprintf(msg, "In file %s, number of component is not given, see %s", fname, buffer); show_msg(msg); in.close(); return false;
	}
	if (n_monomer < 1) {
		sprintf(msg, "ERROR: number of component has to be > 0, check the definition in file %s", fname); show_msg(msg); in.close(); return false;
	}
	ARRAY<MONOMER_GROUP> mgroup; mgroup.SetArray(n_monomer);
	i = 0;
	while (i < n_monomer) {
		if (in.eof()) {
			sprintf(msg, "ERROR: do not get enough component -- get %d expect %d, check the definition in file %s", i, n_monomer, fname); 
			show_msg(msg); in.close(); return false;
		}
		in.getline(buffer, 250);
		if ((npars = sscanf(buffer, "%d %s %s", &indx, mname, sfname)) < 2) {
			sprintf(msg, "ERROR: failure to get component defined in %s, see:\n%s", fname, buffer); show_msg(msg);
			in.close(); return false;
		}
		if (indx < 0 || indx >= n_monomer) continue;

		if (npars == 3 && ((sfname[0] >= 'A' && sfname[0] <= 'Z') || (sfname[0] >= 'a' && sfname[0] <= 'z') || (sfname[0] >= '0' && sfname[0] <= '9'))) {
			if (!ReadMonomersFromFile(sfname, mname, mgroup.m[indx])) {in.close(); return false;}
		}
		else {
			mgroup.m[indx].monomer.SetArray(1);
			if (!construct_monomer(mname, mgroup.m[indx].monomer.m)) {
				sprintf(msg, "UNKNOWN component %d %s defined in %s", indx, mname, fname); show_msg(msg);
				in.close(); return false;
			}
			// indx the monomer-group hinge based on the monomer's free-hinge
			mgroup.m[indx].mHinge.SetArray(mgroup.m[indx].monomer.m[0].nHinges);
			for (iHinge1 = 0; iHinge1 < mgroup.m[indx].monomer.m[0].nHinges; iHinge1++) {
				mgroup.m[indx].mHinge.m[iHinge1].mIndx = 0;
				mgroup.m[indx].mHinge.m[iHinge1].iHinge = iHinge1;
			}
		}
		i++;
	}


	ARRAY<MONOMER_LINK_RELSHIP> linkgroup;
	if (n_monomer > 1) { // no link when monomers is 1
		if (!search(in, "[LINK]", buffer)) {
			sprintf(msg, "[LINK] is not given in file %s", fname); show_msg(msg); in.close(); return false;
		}
		if (sscanf(buffer, "%s %d", title, &nlink) != 2) {
			sprintf(msg, "number of links is not given in %s, see %s", fname, buffer); show_msg(msg); in.close(); return false;
		}
		if (nlink < 1) {
			sprintf(msg, "number of link is expected > 0, see %s", fname, buffer); show_msg(msg); in.close(); return false;
		}
		linkgroup.SetArray(nlink);
		i = 0;
		while (i < nlink) {
			if (in.eof()) {
				sprintf(msg, "do not get enough links in %s, get %d expcet %d", fname, i, nlink); show_msg(msg); in.close(); return false;
			}
			in.getline(buffer, 250);
			if (sscanf(buffer, "%d %d %d %d %d", &indx, &nc1, &iHinge1, &nc2, &iHinge2) != 5) {
				sprintf(msg, "ERROR format of link [indx, nc1, iHinge1, nc2, iHinge2] in %s: %s", fname, buffer); show_msg(msg); in.close(); return false;
			}
			if (indx < 0 || indx >= nlink) continue;
			linkgroup.m[indx].mIndx1 = nc1; linkgroup.m[indx].iHinge1 = iHinge1;
			linkgroup.m[indx].mIndx2 = nc2; linkgroup.m[indx].iHinge2 = iHinge2;
			i++;
		}
	}
	in.close();

	// constructing macromolecules
	int img, im, ic, ilink, m_indx, c_indx = 0;
	// sort the monomer index
	n_monomer = 0;
	for (img = 0; img < mgroup.n; img++) n_monomer += mgroup.m[img].monomer.n;

	ARRAY<MONOMER*> pm;
	pm.SetArray(n_monomer);
	ARRAY<int> mIndx, cmIndx; 
	mIndx.SetArray(mgroup.n); mIndx.m[0] = 0;
	cmIndx.SetArray(n_monomer); cmIndx.m[0] = 0;
	m_indx = 0; 
	for (img = 0; img < mgroup.n; img++) {
		for (im = 0; im < mgroup.m[img].monomer.n; im++) {
			pm.m[m_indx] = mgroup.m[img].monomer.m + im; 
			if (m_indx > 0) cmIndx.m[m_indx] = cmIndx.m[m_indx-1] + pm.m[m_indx-1]->nCluster;
			m_indx++;
		}
		if (img != 0) mIndx.m[img] = mIndx.m[img-1] + mgroup.m[img-1].monomer.n;
	}

	mm.nCluster = 0;
	for (im = 0; im < pm.n; im++) mm.nCluster += pm.m[im]->nCluster;
	mm.cluster = new CLUSTER[mm.nCluster];
	mm.monomer = new _GROUP<CLUSTER>[n_monomer]; mm.nMonomer = n_monomer;
	CLUSTER *cluster = mm.cluster;
	c_indx = 0;
	for (im = 0; im < pm.n; im++) {
		copy_monomer(cluster, pm.m[im]);
		mm.monomer[im].set(pm.m[im]->nCluster);
		for (ic = 0; ic < pm.m[im]->nCluster; ic++) {
			cluster[ic].cmIndx = c_indx;
			cluster[ic].cgIndx = 0;
			mm.monomer[im].u[ic] = cluster + ic; // setup monomer
			cluster[ic].monIndx = im;
		}
		mm.monomer[im].uid = pm.m[im]->uid;
		cluster = cluster + pm.m[im]->nCluster;
	}

	//setup connection defined inside each monomer
	int ic1, ic2, im1, im2;
	for (img = 0; img < mgroup.n; img++) {
		for (im = 0; im < mgroup.m[img].monomer.n; im++) {
			for (ilink = 0; ilink < mgroup.m[img].monomer.m[im].nrship; ilink++) {
				im1 = mIndx.m[img] + im;
				ic1 = cmIndx.m[im1] + mgroup.m[img].monomer.m[im].rship[ilink].n1;
				ic2 = cmIndx.m[im1] + mgroup.m[img].monomer.m[im].rship[ilink].n2;
				if (!ConnectClusters(mm.cluster + ic1, mgroup.m[img].monomer.m[im].rship[ilink].iHinge1, 
					mm.cluster + ic2, mgroup.m[img].monomer.m[im].rship[ilink].iHinge2, 
					mgroup.m[img].monomer.m[im].rship[ilink].bound_type)) {
					sprintf(msg, "failure to connect cluster %d with %d, in monomer %d of monomer-group %d %s", mgroup.m[img].monomer.m[im].rship[ilink].n1,
						mgroup.m[img].monomer.m[im].rship[ilink].n2, im, img, mgroup.m[img].name);
					show_msg(msg); mm.reset(); return false;
				}
			}
		}
	}
	// setup link defined between monomers inside the each monomer-group
	for (img = 0; img < mgroup.n; img++) {
		for (ilink = 0; ilink < mgroup.m[img].mLink.n; ilink++) {
			im1 = mIndx.m[img] + mgroup.m[img].mLink.m[ilink].mIndx1;
			if (mgroup.m[img].mLink.m[ilink].iHinge1 < 0 || mgroup.m[img].mLink.m[ilink].iHinge1 >= pm.m[im1]->nHinges) {
				sprintf(msg, "link-hinge %d in monomer-group %d %s, is not avaliable", ilink, img, mgroup.m[img].name); show_msg(msg); 
				mm.reset(); return false;
			}
			ic1 = cmIndx.m[im1] + pm.m[im1]->nc[mgroup.m[img].mLink.m[ilink].iHinge1];
			iHinge1 = pm.m[im1]->iHinge[mgroup.m[img].mLink.m[ilink].iHinge1];

			im2 = mIndx.m[img] + mgroup.m[img].mLink.m[ilink].mIndx2;
			if (mgroup.m[img].mLink.m[ilink].iHinge2 < 0 || mgroup.m[img].mLink.m[ilink].iHinge2 >= pm.m[im2]->nHinges) {
				sprintf(msg, "link-hinge %d in monomer-group %d %s, is not avaliable", ilink, img, mgroup.m[img].name); show_msg(msg); 
				mm.reset(); return false;
			}
			ic2 = cmIndx.m[im2] + pm.m[im2]->nc[mgroup.m[img].mLink.m[ilink].iHinge2];
			iHinge2 = pm.m[im2]->iHinge[mgroup.m[img].mLink.m[ilink].iHinge2];

			if (!ConnectClusters(mm.cluster + ic1, iHinge1, mm.cluster + ic2, iHinge2, SINGLE_BOUND)) {
				sprintf(msg, "failure to setup the %d link defined in monomer-group %d", ilink, img); show_msg(msg);
				mm.reset(); return false;
			}
		}
	}
	
	// setup link defined between monomer-groups
	MONOMER_GROUP_HINGE *pmg1_hinge = NULL, *pmg2_hinge = NULL;
	for (ilink = 0; ilink < linkgroup.n; ilink++) {
		if (linkgroup.m[ilink].mIndx1 < 0 || linkgroup.m[ilink].mIndx1 >= mgroup.n 
			|| linkgroup.m[ilink].iHinge1 < 0 || linkgroup.m[ilink].iHinge1 >= mgroup.m[linkgroup.m[ilink].mIndx1].mHinge.n
			|| linkgroup.m[ilink].mIndx2 < 0 || linkgroup.m[ilink].mIndx2 >= mgroup.n 
			|| linkgroup.m[ilink].iHinge2 < 0 || linkgroup.m[ilink].iHinge2 >= mgroup.m[linkgroup.m[ilink].mIndx2].mHinge.n) {
				sprintf(msg, "link-hinge %d between monomer-group %d and %d, is not avaliable", ilink,
					linkgroup.m[ilink].mIndx1, linkgroup.m[ilink].mIndx2); show_msg(msg); mm.reset(); return false;
		}
		pmg1_hinge = mgroup.m[linkgroup.m[ilink].mIndx1].mHinge.m + linkgroup.m[ilink].iHinge1;
		im1 = mIndx.m[linkgroup.m[ilink].mIndx1] + pmg1_hinge->mIndx;
		ic1 = cmIndx.m[im1] + pm.m[im1]->nc[pmg1_hinge->iHinge];
		iHinge1 = pm.m[im1]->iHinge[pmg1_hinge->iHinge];

		pmg2_hinge = mgroup.m[linkgroup.m[ilink].mIndx2].mHinge.m + linkgroup.m[ilink].iHinge2;
		im2 = mIndx.m[linkgroup.m[ilink].mIndx2] + pmg2_hinge->mIndx;
		ic2 = cmIndx.m[im2] + pm.m[im2]->nc[pmg2_hinge->iHinge];
		iHinge2 = pm.m[im2]->iHinge[pmg2_hinge->iHinge];

		if (!ConnectClusters(mm.cluster + ic1, iHinge1, mm.cluster + ic2, iHinge2, SINGLE_BOUND)) {
			sprintf(msg, "failure to setup the %d link defined between monomer-groups: %d [%s] -- %d [%s] ", ilink, 
				linkgroup.m[ilink].mIndx1, mgroup.m[linkgroup.m[ilink].mIndx1].name, linkgroup.m[ilink].mIndx2, mgroup.m[linkgroup.m[ilink].mIndx2].name); show_msg(msg);
			mm.reset(); return false;
		}
	}

	int nbase = 0, nbase_monomer = 0;
	if (!read_msg(fname, buffer, "[BASE_CLUSTER]") || sscanf(buffer, "%d", &nbase_monomer) != 1 ||
		nbase_monomer < 0 || nbase_monomer >= pm.n) nbase = mm.nCluster / 2; // setup the base in default
	else {
		nbase = cmIndx.m[nbase_monomer] + pm.m[nbase_monomer]->nbase;
	}
	mm.base = mm.cluster + nbase;
	mm.base->disable_parent();
	setup_parent_child_relationship(mm.base, true); // report error when free hinge exist
	mm.checkMolecule(true);

	if (!check_parent_child_setup(mm.cluster, mm.nCluster)) return false;

	if (!mm.linkPars(atompar_db, errmsg)) {show_infor(errmsg); return false;}

	mm.calTorsionAngle();

	mm.setup_cluster_matom(); // require the atomic mass, and the parent-child relationship

	// set mass center to {0, 0, 0}
	VECTOR<3> cm;
	mm.SetMassCenter(cm.v, true, true);

	strcpy(mm.mol, fname);

	char template_fname[256] = "\0";
	if (save) {
		if (read_msg(fname, template_fname, "[MOL_TEMPLATE_FILENAME]")) {
			::save_xyz_struct(mm, template_fname);
		}
	}
	return true;
}
