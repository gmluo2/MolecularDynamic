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

void rand_vect3(double *v) {
	double theta = 0, p = 0;
	while (1) {
		theta = PI * ranf();
		p = sin(theta); // theta is [0, PI], sin(theta) [0, 1]
		if (p < 0) p = -p;
		if (ranf() <= p) break; // the possibility at theta is sin(theta)
	}
	double omega = PI2 * ranf();
	v[2] = cos(theta);
	double xy = sin(theta);
	v[0] = xy * cos(omega);
	v[1] = xy * sin(omega);
}

double rand_circle_radius() {
	double r = ranf();
	while (r < ranf()) 	r = ranf();
	return r;
}

void rand_circle_vect(double *v) {
	double omega = ranf() * PI2;
	v[0] = cos(omega);
	v[1] = sin(omega);
}


double rand_sphere_shell_radius() {
	double r = ranf();
	while (r * r < ranf()) r = ranf();
	return r;
}

void project(double *u, double *axis, double *p) {
	// projection : p = (u * axis) / (|u| * |axis|) * axis.udir() = (u * axis) / (|u| * |A|^2) * axis
	double A2 = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]; 
	if (A2 == 0) {memset(p, 0, SIZE_V3); return;}
	double uabs = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
	if (uabs == 0) {memset(p, 0, SIZE_V3); return;}
	double c = (u[0] * axis[0] + u[1] * axis[1] + u[2] * axis[2]);
	c /= uabs * A2;
	p[0] = axis[0] * c; p[1] = axis[1] * c; p[2] = axis[2] * c;
}

void rotate(VECTOR3 &axis, double omega, VECTOR3 &v, VECTOR3 &result, bool arc)
{
	double theta, phi;
	theta = axis.theta();
	phi = axis.phi();
	result = R(3, -phi, true) * v;
	result = R(2, -theta, true) * result;
	result = R(3, omega, arc) * result;
	result = R(2, theta, true) * result;
	result = R(3, phi, true) * result;
}

void rotate(double theta, double phi, double omega, VECTOR3 &v, VECTOR3 &result, bool arc)
{
	result = R(3, -phi, arc) * v;
	result = R(2, -theta, arc) * result;
	result = R(3, omega, arc) * result;
	result = R(2, theta, arc) * result;
	result = R(3, phi, arc) * result;
}

//angle between u & v, in arc degree
double angle(VECTOR3 &u, VECTOR3 &v) {
	double uv = (u.v[0] * u.v[0] + u.v[1] * u.v[1] + u.v[2] * u.v[2]) * 
		(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]);
	if (uv == 0) return 0;
	uv = sqrt(uv);
	double c = (u.v[0] * v.v[0] + u.v[1] * v.v[1] + u.v[2] * v.v[2]);
	c /= uv;
	if (c >= 1.0) return 0; //sometime, this could happen accidently, slightly > 1
	else if (c <= -1.0) return PI; //sometime, this could happen accidently, slightly < -1
	else return acos(c);
}

void RotMatrix(int axis, double angle, bool arc, MATRIX<3> &r, MATRIX<3> *inv_r) {
	int i, j;
	double arc_angle = angle;
	if (!arc) arc_angle = angle * Angle2Radian;
	switch (axis) {
	case 1: // x axis
			r.m[0][0] = 1.0;
			r.m[0][1] = 0.0;
			r.m[0][2] = 0.0;
			r.m[1][0] = 0.0;
			r.m[1][1] = cos(arc_angle);
			r.m[1][2] = -sin(arc_angle);
			r.m[2][0] = 0.0;
			r.m[2][1] = -r.m[1][2];
			r.m[2][2] = r.m[1][1];
			if (inv_r != NULL) {
				inv_r->m[0][0] = 1.0;
				inv_r->m[0][1] = 0.0;
				inv_r->m[0][2] = 0.0;
				inv_r->m[1][0] = 0.0;
				inv_r->m[1][1] = r.m[1][0]; // cos(-angle) = cos(angle)
				inv_r->m[1][2] = -r.m[1][2]; // sin(-angle) = -sin(angle)
				inv_r->m[2][0] = 0.0;
				inv_r->m[2][1] = -inv_r->m[1][2];
				inv_r->m[2][2] = inv_r->m[1][1];
			}
			break;
	case 2: // y axis
			r.m[0][0] = cos(arc_angle);
			r.m[0][1] = 0.0;
			r.m[0][2] = sin(arc_angle);
			r.m[1][0] = 0.0;
			r.m[1][1] = 1.0;
			r.m[1][2] = 0.0;
			r.m[2][0] = -r.m[0][2];
			r.m[2][1] = 0.0;
			r.m[2][2] = r.m[0][0];
			if (inv_r != NULL) {
				inv_r->m[0][0] = r.m[0][0]; // cos(-angle)
				inv_r->m[0][1] = 0.0;
				inv_r->m[0][2] = -r.m[0][2]; //sin(-angle);
				inv_r->m[1][0] = 0.0;
				inv_r->m[1][1] = 1.0;
				inv_r->m[1][2] = 0.0;
				inv_r->m[2][0] = -inv_r->m[0][2];
				inv_r->m[2][1] = 0.0;
				inv_r->m[2][2] = inv_r->m[0][0];
			}
			break;
	case 3: // z axis
			r.m[0][0] = cos(arc_angle);
			r.m[0][1] = -sin(arc_angle);
			r.m[0][2] = 0.0;
			r.m[1][0] = -r.m[0][1];
			r.m[1][1] = r.m[0][0];
			r.m[1][2] = 0.0;
			r.m[2][0] = 0.0;
			r.m[2][1] = 0.0;
			r.m[2][2] = 1.0;
			if (inv_r != NULL) {
				inv_r->m[0][0] = r.m[0][0]; //cos(-arc_angle);
				inv_r->m[0][1] = -r.m[0][1]; //sin(-arc_angle);
				inv_r->m[0][2] = 0.0;
				inv_r->m[1][0] = -inv_r->m[0][1];
				inv_r->m[1][1] = inv_r->m[0][0];
				inv_r->m[1][2] = 0.0;
				inv_r->m[2][0] = 0.0;
				inv_r->m[2][1] = 0.0;
				inv_r->m[2][2] = 1.0;
			}
			break;
	default: //undefined
			for (i = 0; i < 3; i++) {
				for(j = 0; j < 3; j++) {
					r.m[i][j] = 0.0;
				}
			}
			break;
	}
}

// create the MATRIX for rotation
void RMatrix(VECTOR3 &axis, double omega, R &r, bool arc)
{
	int i, j, k;
	if (FABS(omega) < MIN_ROT) {I3Matrix(r)  return;}
	if (axis.v[0] == 0 && axis.v[1] == 0 && axis.v[2] == 0) {I3Matrix(r) return;}
	double theta, phi;
	theta = axis.theta();
	phi = axis.phi();
	/*
	result = R(3, -phi) * v;
	result = R(2, -theta) * result;
	result = R(3, omega) * result;
	result = R(2, theta) * result;
	result = R(3, phi) * result;
	*/
	/*
	r = R(3, phi, true) * R(2, theta, true) * R(3, omega, arc) * R(2, -theta, true) * R(3, -phi, true);
	*/


	MATRIX<3> R3_phi, R2_theta, R3_omega, inv_R2_theta, inv_R3_phi, m1, m2, m3;
	RotMatrix(3, phi, true, R3_phi, &inv_R3_phi);
	RotMatrix(2, theta, true, R2_theta, &inv_R2_theta);
	RotMatrix(3, omega, arc, R3_omega, NULL);
	
	MpM(R3_phi, R2_theta, m1, i, j, k, 3)
	MpM(m1, R3_omega, m2, i, j, k, 3)
	MpM(inv_R2_theta, inv_R3_phi, m3, i, j, k, 3)
	MpM(m2, m3, r, i, j, k, 3)

}
