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
#include "Matrix.h"
#include "bound.h"


void DIHEDRAL_INTERACT::init() {
	VECT3((*r1), (*r2), r12)
	VECT3((*r2), (*r3), r23)
	VECT3((*r3), (*r4), r34)
	V3PV3(r12, r23, rp123)
	V3PV3(r23, r34, rp234)
	V3ABS2(rp123, R2p123) Rp123 = sqrt(R2p123);
	V3ABS2(rp234, R2p234) Rp234 = sqrt(R2p234);
	inv_Rp = 1 / (Rp123 * Rp234);
	scale_uv3(rp123, rp234, B) B /= Rp123 * Rp234;
	if (B > 1) {B = 1; phi = 0;}
	else if (B < -1) {B = -1; phi = PI;}
	else phi = acos(B);

	com_12_12.v[0] = DIH_V3_COMM_0(r12.v, r12.v);
	com_12_12.v[1] = DIH_V3_COMM_1(r12.v, r12.v);
	com_12_12.v[2] = DIH_V3_COMM_2(r12.v, r12.v);

	com_12_23.v[0] = DIH_V3_COMM_0(r12.v, r23.v);
	com_12_23.v[1] = DIH_V3_COMM_1(r12.v, r23.v);
	com_12_23.v[2] = DIH_V3_COMM_2(r12.v, r23.v);

	com_12_34.v[0] = DIH_V3_COMM_0(r12.v, r34.v);
	com_12_34.v[1] = DIH_V3_COMM_1(r12.v, r34.v);
	com_12_34.v[2] = DIH_V3_COMM_2(r12.v, r34.v);

	com_23_23.v[0] = DIH_V3_COMM_0(r23.v, r23.v);
	com_23_23.v[1] = DIH_V3_COMM_1(r23.v, r23.v);
	com_23_23.v[2] = DIH_V3_COMM_2(r23.v, r23.v);

	com_23_34.v[0] = DIH_V3_COMM_0(r23.v, r34.v);
	com_23_34.v[1] = DIH_V3_COMM_1(r23.v, r34.v);
	com_23_34.v[2] = DIH_V3_COMM_2(r23.v, r34.v);

	com_34_34.v[0] = DIH_V3_COMM_0(r34.v, r34.v);
	com_34_34.v[1] = DIH_V3_COMM_1(r34.v, r34.v);
	com_34_34.v[2] = DIH_V3_COMM_2(r34.v, r34.v);
}

void DIHEDRAL_INTERACT::derive_B() {
	double dv = 0;
	double dv1, dv2, dv3;
	int i;
	// atom 1
	for (i = 0; i < 3; i++) {
		dv1 = -r23.v[i] * com_23_34.v[i] + r34.v[i] * com_23_23.v[i];
		dv2 = (-r12.v[i] * com_23_23.v[i] + r23.v[i] * com_12_23.v[i]) * 2; 
		dv2 /= R2p123;
		dv3 = 0;
		dB1.v[i] = inv_Rp * dv1 - 0.5 * B * (dv2 + dv3);
	}

	// atom 2
	for (i = 0; i < 3; i++) {
		dv1 = -r12.v[i] * com_23_34.v[i] + r23.v[i] * com_23_34.v[i] + r34.v[i] * (-com_12_23.v[i] -com_23_23.v[i]) + 2 * r23.v[i] * com_12_34.v[i];
		dv2 = (r12.v[i] * (com_23_23.v[i] + com_12_23.v[i]) + r23.v[i] * (-com_12_12.v[i] - com_12_23.v[i])) * 2; 
		dv2 /= R2p123;
		dv3 = (r34.v[i] * com_23_34.v[i] - r23.v[i] * com_34_34.v[i]) * 2;
		dv3 /= R2p234;
		dB2.v[i] = inv_Rp * dv1 - 0.5 * B * (dv2 + dv3);
	}

	// atom 3
	for (i = 0; i < 3; i++) {
		dv1 = r12.v[i] * (com_23_23.v[i] + com_23_34.v[i]) - r23.v[i] * com_12_23.v[i] + r34.v[i] * com_12_23.v[i] - 2 * r23.v[i] * com_12_34.v[i];
		dv2 = (-r12.v[i] * com_12_23.v[i] + r23.v[i] * com_12_12.v[i]) * 2; 
		dv2 /= R2p123;
		dv3 = (r34.v[i] * (-com_23_23.v[i] - com_23_34.v[i]) + r23.v[i] * (com_34_34.v[i] + com_23_34.v[i])) * 2;
		dv3 /= R2p234;
		dB3.v[i] = inv_Rp * dv1 - 0.5 * B * (dv2 + dv3);
	}

	// atom 4
	for (i = 0; i < 3; i++) {
		dv1 = -r12.v[i] * com_23_23.v[i] + r23.v[i] * com_12_23.v[i];
		dv2 = 0; 
		dv3 = (r34.v[i] * com_23_23.v[i] - r23.v[i] * com_23_34.v[i]) * 2;
		dv3 /= R2p234;
		dB4.v[i] = inv_Rp * dv1 - 0.5 * B * (dv2 + dv3);
	}
}

void DIHEDRAL_INTERACT::cal_virial() { // this virial function is not right, although it is inherited from DL_POLY manual
	double p1, p4, p23, g1, g3, h2, h4;
	int i, j;
	for (i = 0; i < 3; i++) {
		j = i;
		p1 = (r23.v[j] * com_23_34.v[j] - r34.v[j] * com_23_23.v[j]) * inv_Rp;
		p4 = (r23.v[j] * com_12_23.v[j] - r12.v[j] * com_23_23.v[j]) * inv_Rp;
		p23 = (r12.v[j] * com_23_34.v[j] + r34.v[j] * com_12_23.v[j] - 2 * r23.v[j] * com_12_34.v[j]) * inv_Rp;
		g1 = 2 * (r12.v[j] * com_23_23.v[j] - r23.v[j] * com_12_23.v[j]) / R2p123;
		g3 = 2 * (r23.v[j] * com_12_12.v[j] - r12.v[j] * com_12_23.v[j]) / R2p123;
		h2 = 2 * (r23.v[j] * com_34_34.v[j] - r34.v[j] * com_23_34.v[j]) / R2p234;
		h4 = 2 * (r34.v[j] * com_34_34.v[j] - r23.v[j] * com_23_34.v[j]) / R2p234;
		
		virial.v[i] = r12.v[i] * p1 + r23.v[i] * p23 + r34.v[i] * p4 - 0.5 * B * (r12.v[i] * g1 + r23.v[i] * (g3 + h2) + r34.v[i] * h4);
	}
}


bool VECTOR_decompose(double f, double *uf, double *u1, double *u2, double &f1, double &f2) {
	// IMPORTANT: u1, u2 and uf has to be in the same plane
	double C1 = u1[0] * uf[0] + u1[1] * uf[1] + u1[2] * uf[2]; // cos(u1^uf)
	double S1 = 1 - C1 * C1;
	if (S1 > 0) S1 = sqrt(S1);
	else {C1 = sign(C1); S1 = 0;}

	double C2 = u2[0] * uf[0] + u2[1] * uf[1] + u2[2] * uf[2]; // cos(u2^uf)
	double S2 = 1 - C2 * C2;
	if (S2 > 0) S2 = sqrt(S2);
	else {C2 = sign(C2); S2 = 0;}

	if (C1 == 0 && C2 == 0) {
		if (f == 0) {f1 = 0; f2 = 0; return true;}
		else return false;
	}

	if (S1 == 0) {
		if (C1 == 1) f1 = f;
		else f1 = -f;
		f2 = 0; return true;
	}
	else if (S2 == 0) {
		if (C2 == 1) f2 = f;
		else f2 = -f;
		f1 = 0; return true;
	}

	double axis1[3], axis2[3];
	DV3PV3(u1, uf, axis1) DV3PV3(u2, uf, axis2)
	double dir = axis1[0] * axis2[0] + axis1[1] * axis2[1] + axis1[2] * axis2[2];

	if (dir < 0) { // uf is between u1 and u2
		f1 = f / (C1 + C2 * S1 / S2);
		f2 = f1 * S1 / S2;
	}
	else {
		f1 = f / (C1 - C2 * S1 / S2);
		f2 = -f1 * S1 / S2;
	}

	// check whether the result is right
	/*
	double df[3];
	df[0] = f1 * u1[0] + f2 * u2[0] - f * uf[0];
	df[1] = f1 * u1[1] + f2 * u2[1] - f * uf[1];
	df[2] = f1 * u1[2] + f2 * u2[2] - f * uf[2];
	*/

	return true;
}

bool HINGE_CONSTRAINT_FORCE::cal_constraint_force() {
	// The additional torque between the cluster, from the constraints here, is 0 (ZERO) along the hinge.
	// This is is CONFIRMED by NEIMO method, that T = H * f_hinge 
	// see the literature, "A Fast Recursive Althorithm for Molecular Dynamics Simulation", by A. Jain et al

	// the forces from 10 to 01 and 02 do NOT contribute to the torque from cluster #1 to cluster #0, 
	// because the torque is relative to position 10
	// So, torque on hinge, from cluster #1 to #0, is only from the froces from 11 and 12 on 00
	char msg[256] = "\0";
	SVECTOR<double, 3> torq1, torq2;
	V3PV3(r00, u_11_00, torq1) V3PV3(r00, u_12_00, torq2)
	double tq1, tq2;
	SVECTOR<double, 3> u_torq1, u_torq2;
	V3ABS2(torq1, tq1) tq1 = sqrt(tq1);
	V3ABS2(torq2, tq2) tq2 = sqrt(tq2);
	if (tq1 > 0) {
		u_torq1.v[0] = torq1.v[0] / tq1;
		u_torq1.v[1] = torq1.v[1] / tq1;
		u_torq1.v[2] = torq1.v[2] / tq1;
	}
	else {
		sprintf(msg, "strange: hinge-constraint points on child cluster were chosen dependent"); show_log(msg, true); return false;
	}
	if (tq2 > 0) {
		u_torq2.v[0] = torq2.v[0] / tq2;
		u_torq2.v[1] = torq2.v[1] / tq2;
		u_torq2.v[2] = torq2.v[2] / tq2;
	}
	else {
		sprintf(msg, "strange: hinge-constraint points on parent cluster were chosen dependent"); show_log(msg, true); return false;
	}

	// since torque has direction normal to the hinge, torq1 and torq2 are also normal to the hinge, torq, torq1 and torq2 are in the same plane
	// and we know the directions of torq1 and torq2, the forces from 11 and 12 to 00 can be calculated with torque : torq = f_11_00 * torq1 +  f_12_00 * torq2
	
	bool status = VECTOR_decompose(Torq, u_torq.v, u_torq1.v, u_torq2.v, f_11_00, f_12_00);
	if (!status) return false;
	f_11_00 /= tq1; f_12_00 /= tq2;

	// check whether the result is right
	/*
	SVECTOR<double, 3> tq;
	tq.v[0] = torq1.v[0] * f_11_00 + torq2.v[0] * f_12_00 - torq.v[0];
	tq.v[1] = torq1.v[1] * f_11_00 + torq2.v[1] * f_12_00 - torq.v[1];
	tq.v[2] = torq1.v[2] * f_11_00 + torq2.v[2] * f_12_00 - torq.v[2];
	sprintf(msg, "[%f, %f, %f]", tq.v[0], tq.v[1], tq.v[2]); show_log(msg, true);
	*/

	//double f_10_00, f_10_01, f_10_02; // forces on atom 00, 01, 02 from 10
	//These three forces are to be calculated with matrix
	// variable (force) indexed 0 -- 10==>00, 1 -- 10==>01, 2 -- 10==>02
	// X axis
	m.m[0][0] = u_10_00.v[0];
	m.m[0][1] = u_10_01.v[0];
	m.m[0][2] = u_10_02.v[0];
	// Y axis
	m.m[1][0] = u_10_00.v[1];
	m.m[1][1] = u_10_01.v[1];
	m.m[1][2] = u_10_02.v[1];
	// Z axis
	m.m[2][0] = u_10_00.v[2];
	m.m[2][1] = u_10_01.v[2];
	m.m[2][2] = u_10_02.v[2];

	b.v[0] = f.v[0] - f_11_00 * u_11_00.v[0] - f_12_00 * u_12_00.v[0];
	b.v[1] = f.v[1] - f_11_00 * u_11_00.v[1] - f_12_00 * u_12_00.v[1];
	b.v[2] = f.v[2] - f_11_00 * u_11_00.v[2] - f_12_00 * u_12_00.v[2];

	status = InverseMatrix<3>(m, invm);
	if (!status) return false;
	SVECTOR<double, 3> ft;
	MpV3(invm, b, ft)
	f_10_00 = ft.v[0]; f_10_01 = ft.v[1]; f_10_02 = ft.v[2];
	return true;
}

