#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
#include "show.h"
#include "vector.h"

// functions in Z-matrix to define a new atom cordinate

// r1 -- adjacent atom; d -- shift vector
void ZMatrix(VECTOR3 &r1, VECTOR3 &d, VECTOR3 &res) {
	res.v[0] = r1.v[0] + d.v[0];
	res.v[1] = r1.v[1] + d.v[1];
	res.v[2] = r1.v[2] + d.v[2];
}

// r1, r2 are atoms' position, bl -- bound distance 
// axis -- normal to new-1-2 plane, R_1-2 X R_new-1
// ba -- bound angle between new-1-2, in degree
void ZMatrix(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &axis, float bl, float ba, VECTOR3 &res) {
	VECTOR3 v1 = r2 - r1;
	VECTOR3 v2;
	rotate(axis, ba, v1, v2, false); //rotate v1 along axis with angle ba; ==> v2
	float len = (float)v2.abs();
	res = v2 * (bl / len);
	res += r1;
}

// bl -- bound distance; ba -- bound angle between new-1-2
// ta -- torsion angle between new-1-2 & 3-2-1 along axis 2-1
// r1, r2 & r3 -- atoms' position
void ZMatrix(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &r3, float bl, float ba, float ta, VECTOR3 &res) {
	VECTOR3 v1 = r2 - r1;
	VECTOR3 v2 = r3 - r2;
	VECTOR3 v3;
	cp(v1, v2, v3); // v3 = v1 x v2
	rotate(v3, ba, v1, v2, false); //rotate v1 along v3 with angle ba; ==> v2
	float len = (float)v2.abs();
	v2 *= (bl / len);
	rotate(v1, ta, v2, res, false); // rotate v2 along v1 with angle ta; ==> res
	res += r1;
}

// Torsion Angle -- to get r1 with r2-r3-r4 axis
// reverse procedure of ZMatrix(r1, r2, r3, bound-length, bound-angle)
// return : arc degree [-PI, PI]
double TorsionAngle(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &r3, VECTOR3 &r4) {
	VECTOR3 dr1 = r3 - r2, dr2 = r1 - r2;
	VECTOR3 n1, n2;
	cp(dr2, dr1, n1); // n1 = r21 x r23
	dr1 = r2 - r3;
	dr2 = r4 - r3;
	cp(dr1, dr2, n2); // n2 = r32 x r34
	double ta = angle(n1, n2); // angle between n1 & n2
	if (ta == 0 || ta == PI) return ta;
	VECTOR3 axis;
	cp(n1, n2, axis); // axis = n1 x n2
	dr1 = r2 - r3; // r23
	double c2 = axis * dr1;
	if (c2 >= 0) return ta; // axis is same as r23
	else return -ta; // axis is anti-parallel to r23
}

// bound angle between r21 and r23 -- 1-2-3
// return : arc degree [0, PI]
double BoundAngle(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &r3) {
	VECTOR3 v1 = r1 - r2;
	VECTOR3 v2 = r3 - r2;
	double ba = angle(v1, v2);
	return ba;
}

