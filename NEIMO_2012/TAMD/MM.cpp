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
#include "ZMatrix.h"
#include "MM.h"
#include "var.h"

#include "Matrix.h"

void BASIC_CLUSTER::calElectroMultipole(bool bEwaldSum) {
#if _EWALD_SUM == _MultiPole_EWALD_SUM
	if (this->eNeutral || emp == NULL) return;
	emp->d = this->d;
	int i, j, k, n;
	VECTOR3 r, rc;
	float q = 0, q0 = 0;
	emp->q = 0; V3zero(emp->mu) V3zero(emp->r0) Mzero(emp->Q, i, j)
	for (n = 0; n < nAtoms; n++) {
		q = atom[n].c; q = FABS(q);
		for (i = 0; i < 3; i++) rc.v[i] += q * atom[n].r0.v[i];
		q0 += q;
	}
	if (q0 == 0) {this->eNeutral = true; return;} //we do not need to calculate other terms, since no charge in this cluster
	// important, multipole is in the cell cordinate
	for (i = 0; i < 3; i++) {
		rc.v[i] /= q0; 
		emp->r0.v[i] = rc.v[i] + this->Op->v[i] + this->dr_cmm.v[i]; // in experiment space
	}
	for (n = 0; n < nAtoms; n++) {
		q = atom[n].c;
		emp->q += q;
		for (i = 0; i < 3; i++) r.v[i] = atom[n].r0.v[i] - rc.v[i];
		for (i = 0; i < 3; i++) {
			emp->mu.v[i] += q * r.v[i];
#if _ATOM_INDUCED_DIPOLE == 1
			emp->mu.v[i] += atom[n].ind_mu.v[i]; // induced dipole -- polarized dipole
#endif
			for (j = 0; j < 3; j++) {
				// we kick off the term in Q, -1/2 * r^2 * delta_i,j, since this term is useless in final calculation
				emp->Q.m[i][j] += 1.5 * q * r.v[i] * r.v[j];
#if _ATOM_INDUCED_DIPOLE == 1
				// quadrupole from induced dipole
				emp->Q.m[i][j] += 1.5 * (r.v[i] * atom[n].ind_mu.v[j] + r.v[j] * atom[n].ind_mu.v[i]);
#endif
			}
		}
	}

	double hR = 0, t1 = 0, t2 = 0;
	int ti, tj;
	if (bEwaldSum) {
		for (i = 0; i < NH; i++) {for (j = 0; j < NH; j++) {for (k = 0; k < NH; k++) {
			scale_uv3(emp->r0, ::h[i][j][k], hR)
			//PERIOD(hR, 0, PI2, t1) COS(emp->fcos.m[i][j][k], t1, n, t1) SIN(emp->fsin.m[i][j][k], t1, ti, t2)
			emp->fcos.m[i][j][k] = cos(hR); emp->fsin.m[i][j][k] = sin(hR);
			scale_uv3(emp->mu, ::h[i][j][k], emp->muh.m[i][j][k])
			scale_VMV(::h[i][j][k], emp->Q, ::h[i][j][k], emp->Qh.m[i][j][k], ti, tj, 3)
		}}}
	}
#endif
}

void CLUSTER::shift(VECTOR3 &d) {
	int i = 0, j;
	SVECTOR<double, 3> dr; V32V3(d, dr);
	for (i = 0; i < nAtoms; i++) {for (j = 0; j < 3; j++) atom[i].r.v[j] += dr.v[j];}
	V3plusV3(cm, dr, cm) // mass center
	if (bImplicitSolvent) {
		V3plusV3(rc_solv, dr, rc_solv) // geometrical center
	}
}

void CLUSTER::rotate(VECTOR3 &axis, double omega, bool arc) { // omega in degree ?
	if (arc) {PERIOD_RANGE(omega, -PI, PI, PI2)}
	else {PERIOD_RANGE(omega, -180, 180, 360)}
	if (FABS(omega) < MIN_ROT) return;
	R r;
	RMatrix(axis, omega, r, arc);
	rotate(r);
}

void CLUSTER::rotate(R &r) { // given rotation matrix
	int i = 0;
	SVECTOR<double, 3> v;
	VECTOR3 *pr;
	for (i = 0; i < nAtoms; i++) {
		pr = &(atom[i].r);
		V32V3((*pr), v); // atom r
		//atom[i].r = r * v;
		MpV3(r, v, (*pr))
		atom[i].rotate1(r);
	}
	//v = up;  // up
	//up = r * v;
	V32V3(up, v);
	MpV3(r, v, up)
	// hinges
	for (i = 0; i < nHinges; i++) {
		//v = hinge[i].vHinge; // hinge vector
		//hinge[i].vHinge = r * v;
		pr = &(hinge[i].vHinge);
		V32V3((*pr), v);
		MpV3(r, v, (*pr));
	}
	V32V3(cm, v)
	//cm = r * v;
	MpV3(r, v, cm)
	if (Op != NULL) cm0 = cm - (*Op);
	if (bImplicitSolvent) {
		V32V3(rc_solv, v)
		//rc_solv = r * v;
		MpV3(r, v, rc_solv);
		if (Op != NULL) {
			//rc0_solv = rc_solv - (*Op);
			V3minusV3(rc_solv, (*Op), rc_solv);
		}
	}
}

void CLUSTER::rotate(VECTOR3 &r1, VECTOR3 &r2, double omega, bool arc) { // r1 is base point, omega in degree ?
	if (FABS(omega) < MIN_ROT) return;
	int i;
	VECTOR3 r01, r02;
	for (i = 0; i < 3; i++) {r01.v[i] = r1.v[i]; r02.v[i] = r2.v[i];}
	VECTOR3 v0;
	for (i = 0; i < 3; i++) v0.v[i]= -(r01.v[i]);
	VECTOR3 axis = r02 - r01;
	shift(v0);
	rotate(axis, omega, arc);
	shift(r01);
}

void CLUSTER::shift_rotate(VECTOR3 &dv, VECTOR3 &v0, R &r) {
	int i = 0, j = 0, k = 0;
	VECTOR3 v;
	VECTOR3 ds, ds0;
	for (i = 0; i < 3; i++) {ds0.v[i] = -(v0.v[i]); ds.v[i] = dv.v[i] + v0.v[i];}

	for (i = 0; i < nAtoms; i++) {
		//for (j = 0; j < 3; j++) v.v[j] = atom[i].r.v[j] + ds0.v[j];
		V3plusV3(atom[i].r, ds0, v)
		//atom[i].r = r * v;
		MpV3(r, v, atom[i].r)
		//for (j = 0; j < 3; j++) atom[i].r.v[j] += ds.v[j];
		V3plusV3(atom[i].r, ds, atom[i].r)

		atom[i].rotate1(r);
	}
	// hinges
	for (i = 0; i < nHinges; i++) {
		//for (j = 0; j < 3; j++) v.v[j] = hinge[i].vHinge.v[j]; // hinge vector
		V32V3(hinge[i].vHinge, v)
		//hinge[i].vHinge = r * v;
		MpV3(r, v, hinge[i].vHinge)
	}
	//for (j = 0; j < 3; j++) v.v[j] = up.v[j];  // up
	V32V3(up, v)
	//up = r * v;
	MpV3(r, v, up)

	V32V3(cm, v) // mass center
	MpV3(r, v, cm)
	if (Op != NULL) {
		V3minusV3(cm, (*Op), cm0)
	}

	if (bImplicitSolvent) {
		V32V3(rc_solv, v)
		//rc_solv = r * v;
		MpV3(r, v, rc_solv)
		if (Op != NULL) {
			//rc0_solv = rc_solv - (*Op);
			V3minusV3(rc_solv, (*Op), rc0_solv)
		}
	}
}

void non_parallel(VECTOR3 &v, VECTOR3 &nv) {
	int i = 0, n = 0;
	//for (i = 0; i < 3; i++) nv.v[i] = v.v[i];
	V32V3(v, nv)
	for (i = 0; i < 3; i++) {if (FABS(v.v[i]) > 1e-7) {n = i; break;}}
	for (i = 0; i < 3; i++) {if (i != n) nv.v[i] += 1;}
}

void non_parallel(VECTOR<3> &v, VECTOR<3> &nv) {
	int i = 0, n = 0;
	//for (i = 0; i < 3; i++) nv.v[i] = v.v[i];
	V32V3(v, nv)
	for (i = 0; i < 3; i++) {if (FABS(v.v[i]) > 1e-7) {n = i; break;}}
	for (i = 0; i < 3; i++) {if (i != n) nv.v[i] += 1;}
}

void CLUSTER::align(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &r0, VECTOR3 &u) // align r1==>r0 & r12 || u
{
	VECTOR3 R0 = r0, R1 = r1;

	VECTOR3 v1 = r2 - r1;
	VECTOR3 dr = r1 * (-1);

	double cos = 0;
	int i = 0;
	scale_uv3(v1, u, cos) // c = v1 * u
	VECTOR3 nu;

	// rotate r12 ==> u
	VECTOR3 axis;
	double omega = 0;
	cp(v1, u, axis);
	if (axis.abs() > 1e-7){
		omega = angle(v1, u); // in arc degree
	}
	else { // the two axises are parallel
		if (cos < 0) { // inverse direction
			non_parallel(u, nu); // nu is an axis which is not parallel to u
			cp(u, nu, axis);
			omega = PI;
		}
		else omega = 0;
	}

	// shift r1 ==> r0
	dr = R0 - R1;
	/*
	if (FABS(omega) > MIN_ROT) tweak_cluster_and_children(this, R1, axis, omega, true); // need parent / child relationship with the clusters associated with it
	shift_cluster_and_children(this, dr);
	*/
	R rmatrix;
	RMatrix(axis, omega, rmatrix, true);
	series_shift_rotate(this, -1, dr, R1, rmatrix);
}

bool CLUSTER::xyz_save(ofstream &out, int indx) {
	int i = 0;
	char buffer[256] = "\0";
	char buff[20] = "\0";
	char atom_type[5] = "\0";
	for (i = 0; i < nAtoms; i++) {
		strcpy(atom_type, atom[i].par->atom);
		sprintf(buffer, "%s %.4f %.4f %.4f", atom_type, atom[i].r.v[0], atom[i].r.v[1], atom[i].r.v[2]);

		//sprintf(buffer, "%s %.4f %.4f %.4f %s", atom[i].atom, atom[i].r.v[0], atom[i].r.v[1], atom[i].r.v[2], atom_type);
		//if (indx >= 0) {sprintf(buff, " %d", indx); strcat(buffer, buff);}

		out<<buffer<<endl;
	}
	return true;
}

void CLUSTER::InertialCord() {
	int i = 0;
	VECTOR3 *v0 = NULL;
	v0 = Op;

	// inertial cordinates of this cluster transfer to parent's inertial cordinates
	// it is actually the inertial vector of Om [to this cluster] in parent's inertial cordinate
	VECTOR3 *v00 = NULL;
	if (parent != NULL) v00 = parent->Op;
	if (v00 == NULL) {
		//for (i = 0; i < 3; i++) l.v[i] = 0;
		V3zero(l)
	}
	else {
		//for (i = 0; i < 3; i++) l.v[i] = v0->v[i] - v00->v[i];
		V3minusV3((*v0), (*v00), l)
	}

	double r = 0;
	VECTOR3 *pr = NULL;
	if (parent != NULL) {
		//r = hinge[n2Parent].vHinge.abs();
		pr = &(hinge[n2Parent].vHinge);
		scale_uv3((*pr), (*pr), r) r = sqrt(r);
		if (r > 1e-6) {
			DIV3((*pr), r, up)
		}
		else {
			V3zero(up)
		}
	}
	else {V3zero(up)}

	for (i = 0; i < nAtoms; i++) {
		atom[i].r0.v[0] = atom[i].r.v[0] - v0->v[0];
		atom[i].r0.v[1] = atom[i].r.v[1] - v0->v[1];
		atom[i].r0.v[2] = atom[i].r.v[2] - v0->v[2];
	}
	if (bImplicitSolvent) {
		if (Op != NULL) {
			//rc0_solv = rc_solv - (*Op);
			V3minusV3(rc_solv, (*Op), rc0_solv)
		}
	}
}

void CLUSTER::MassCenter() {
	int i = 0;
	float m = 0;
	cm.v[0] = 0; cm.v[1] = 0; cm.v[2] = 0;
	M = 0;
	double *v;
	if (matom.m == NULL) {
		for (i = 0; i < nAtoms; i++) {
			m = atom[i].par->m;
			v = atom[i].r.v;
			M += m;
			cm.v[0] += m * v[0];
			cm.v[1] += m * v[1];
			cm.v[2] += m * v[2];
		}
	}
	else {
		for (i = 0; i < matom.n; i++) {
			m = matom.m[i]->par->m;
			v = matom.m[i]->r.v;
			M += m;
			cm.v[0] += m * v[0];
			cm.v[1] += m * v[1];
			cm.v[2] += m * v[2];
		}
	}
	cm.v[0] /= M;
	cm.v[1] /= M;
	cm.v[2] /= M;

	// mass center in global cordinate, local cordinate is relative to the Op
	cm0.v[0] = cm.v[0] - this->Op->v[0];
	cm0.v[1] = cm.v[1] - this->Op->v[1];
	cm0.v[2] = cm.v[2] - this->Op->v[2];
}

void CLUSTER::InertiaTensor() {
	int ni = 0, nj = 0;
	int i = 0, j = 0, k = 0;

	MATRIX<3> m1, m2;
	int na = 0;

	double r2 = 0;
	double *v, mass, rsize;
	SVECTOR<double, 3> r0;
	double *v0 = Op->v;

	M3zero(invNE->I)
	if (matom.m == NULL) { // use the atom
		for (na = 0; na < nAtoms; na++) {
			/*
			u2Matrix3(atom[na].r0, m1.m)
			MpM(m1, m1, m2, i, j, k, 3);
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) 
				invNE->I.m[i][j] -= atom[na].par->m * m2.m[i][j];
			*/
			rsize = atom[na].par->rLJ;
			v = atom[na].r.v;
			DVECT3(v0, v, r0.v)
			v = r0.v;
			r2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
			mass = atom[na].par->m;
			rsize = atom[na].par->rLJ;
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
				invNE->I.m[i][j] -= mass * v[i] * v[j];
				if (i == j) invNE->I.m[i][j] += mass * r2  + 0.4 * mass * rsize * rsize; // last term is the inertia tensor of sphere relative to mass center -- parallel-axis theorem for rotational inertia
			}
		}
	}
	else { // use matom
		for (na = 0; na < matom.n; na++) {
			/*
			u2Matrix3(matom.m[na]->r0, m1.m)
			MpM(m1, m1, m2, i, j, k, 3);
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) 
				invNE->I.m[i][j] -= matom.m[na]->par->m * m2.m[i][j];
			*/
			rsize = matom.m[na]->par->rLJ;
			v = matom.m[na]->r.v;
			DVECT3(v0, v, r0.v)
			v = r0.v;
			r2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
			mass = matom.m[na]->par->m;
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
				invNE->I.m[i][j] -= mass * v[i] * v[j];
				if (i == j) invNE->I.m[i][j] += mass * r2 + 0.4 * mass * rsize * rsize; // last term is the inertia tensor of sphere relative to mass center -- parallel-axis theorem for rotational inertia
			}
		}
	}

	for (i = 0; i < 3; i++) {
		//for (j = 0; j < 3; j++) invNE->M.m[i][j] = invNE->I.m[i][j];
		memcpy(invNE->M.m[i], invNE->I.m[i], SIZE_V3);
	}

	ni = 0; nj = 3;
	SVECTOR<double, 3> mp;
	for (i = 0; i < 3; i++) mp.v[i] = cm0.v[i] * M;
	u2Matrix(mp, invNE->M.m, ni, nj)

	ni = 3; nj = 0;
	for (i = 0; i < 3; i++) mp.v[i] = -(mp.v[i]);
	u2Matrix(mp, invNE->M.m, ni, nj)

	ni = 3; nj = 3;
	me(invNE->M.m, ni, nj, 0, 0) = M; me(invNE->M.m, ni, nj, 1, 1) = M; me(invNE->M.m, ni, nj, 2, 2) = M;

	SVECTOR<double, 3> u1;
	if (parent == NULL) Iw = 0;
	else {
		MpV3(invNE->I, up, u1)
		scale_uv3(u1, up, this->Iw)
	}

	if (mcDyn != NULL) cal_cluster_inertia_tensor_masscenter();
}

void CLUSTER::calJacobianMatrix() {
	// phai & H are Jacobian matrixs of this cluster to its parent
	// In this program, the cordinate frames in all clusters have same direction
	// The only difference is their original points are different (parent Om[])
	// So, H is a constant matrix, H[i][i] = 1. So, we do not have specific variable for H
	int i = 0, j = 0;
	M6zero(km.phai)
	//M6zero(km.phai_T)
	for (i = 0; i < 6; i++) km.phai.m[i][i] = 1;
	u2Matrix(l, km.phai.m, 0, 3) 
	
	M6T2M6(km.phai, km.phai_T, i, j)  // phai_T = phai*

	// the rest elements are 0, matrix was initilized as 0
	/*
	km.H.v[0] = up.v[0];
	km.H.v[1] = up.v[1];
	km.H.v[2] = up.v[2];
	km.H.v[3] = 0;
	km.H.v[4] = 0;
	km.H.v[5] = 0;
	*/
	memcpy(km.H.v, up.v, SIZE_V3);
}

void CLUSTER::ab() {
	MATRIX<3> omega, omega_p;
	SVECTOR<double, 3> w, v, w_p, v_p;

	SVECTOR<double, 3> a1, a2;

	V2wv(w, v, bkm->V0)
	u2Matrix3(w, omega.m)

	V6zero(invNE->a) V6zero(invNE->b)

	// a
	if (parent == NULL) {
		// H = I, H * tp' = V
		// since w(k+1) == 0, the first term can be ignored
		// for the second term,  omega(k) x omega(k) == 0, but omega(k) x v(k) is not zero
		V2wv(w, v, bkm->V0)
		V3PV3(w, v, a1) // a1 = w x v, cross product
		memcpy(invNE->a.v + 3, a1.v, SIZE_V3);
	}
	else {
		V2wv(w_p, v_p, ((BASIC_CLUSTER*)parent)->bkm->V0)
		u2Matrix3(w_p, omega_p.m)

		MpV3(omega_p, up, a1)   // a1 = omega_p * up
		invNE->a.v[0] = a1.v[0] * km.tp; 
		invNE->a.v[1] = a1.v[1] * km.tp;
		invNE->a.v[2] = a1.v[2] * km.tp;

		//for (i = 0; i < 3; i++) a2.v[i] = v.v[i] - v_p.v[i];
		V3minusV3(v, v_p, a2)
		MpV3(omega_p, a2, a1)   //a1 = omega_p * a2
		invNE->a.v[3] = a1.v[0]; 
		invNE->a.v[4] = a1.v[1]; 
		invNE->a.v[5] = a1.v[2];
	}

	// b
	SVECTOR<double, 3> zero; // u
	//V2wv(u, zero, km.H)
	//MpV(omega, u, a1, i, j, 3)
	//for (i = 0; i < 3; i++) a1.v[i] *= km.tp;
	//invNE->a.v[0] += a1.v[0]; invNE->a.v[1] += a1.v[1]; invNE->a.v[2] += a1.v[2];

	MpV3(invNE->I, w, a2)   // a2 = invNE.I * w
	MpV3(omega, a2, a1)    // a1 = omega * a2;

	MpV3(omega, cm0, zero)  // zero = omega * l_cm;
	MpV3(omega, zero, a2)  //a2 = omega * zero;
	TIME3(a2, M, a2) // a2 *= M;

	wv2V(a1, a2, invNE->b)
}

bool MMOLECULE::xyz_save(ofstream &out) {
	int i = 0;
	bool status = true;
	for (i = 0; i < nCluster; i++) status = status & cluster[i].xyz_save(out, i);
	return status;
}

void cp(MATOM *d, MATOM *s) {
	//strcpy(d->name, s->name);
	d->aindx = s->aindx;
	d->c = s->c; d->c0 = s->c0;
	d->r = s->r; d->rg = s->rg;
	d->r0 = s->r0;
	d->par = s->par;
	// copying an atom, information is not enough to form the bound connection
	d->set_bounds(s->nb);
	for (int i = 0; i < d->nb; i++) {d->ba[i] = NULL; d->bv[i] = s->bv[i]; d->btype[i] = s->btype[i];}
	d->eNeutral = s->eNeutral;
	d->mMP = s->mMP;
	d->bIMu = s->bIMu;
	d->intr_mu = s->intr_mu;
	d->ind_mu = s->ind_mu;
	d->mu = s->mu;
	// bound will be connected in cluster copying
	//d->psolv = s->psolv;
}

bool cp(BASIC_CLUSTER *d, BASIC_CLUSTER *s) {
	//if (d->nAtoms != s->nAtoms) return false;
	d->reset();
	if (s->nAtoms <= 0) return true;
	d->nAtoms = s->nAtoms;
	d->atom = new MATOM[d->nAtoms];
	int i = 0, j = 0, n = 0, nb = 0;
	for (i = 0; i < d->nAtoms; i++) cp(d->atom + i, s->atom + i);
	// forming the bound 
	for (i = 0; i < d->nAtoms; i++) {
		for (j = 0; j < d->atom[i].nb; j++) {
			if (s->atom[i].ba[j] != NULL) {
				n = s->iAtom((MATOM*)(s->atom[i].ba[j]));
				if (n < 0) { // this bound could be in other cluster
					//sprintf(errmsg, "ERROR : Program is not consistent. Can not find atom defined in cluster");
					//return false;
				}
				else {
					nb = s->atom[n].ibound(s->atom + i);
					if (nb < 0) {
						sprintf(errmsg, "ERROR: Program is not consisten to find out the bound connection in defined cluster");
						return false;
					}
					if (!bound_connect<MATOM>(d->atom + i, j, d->atom + n, nb, s->atom[i].bv[j], s->atom[i].btype[j])) return false;
				}
			}
		}
	}
	d->set_hinges(s->nHinges);

	for (i = 0; i < s->nHinges; i++) {
		d->hinge[i].vHinge = s->hinge[i].vHinge;
		if (s->hinge[i].indx >= 0) d->setup_hinge(i, s->hinge[i].indx, s->hinge[i].ibound);
		// the connection can not be copied
	}
	d->n2Parent = s->n2Parent;
	d->d = s->d;
	d->eNeutral = s->eNeutral;
	d->cID = s->cID; d->mType = s->mType;

	d->deuterized = s->deuterized;

	//short mIndx; // the molecule index in the whole system, where this cluster is
	d->cgIndx = s->cgIndx; // index of coarse-grained cluster where this cluster is in
	d->monIndx = s->monIndx; // monomer index inside the macro-molecule
	d->cmIndx = s->cmIndx; // cluster index inside macro-molecule
	d->cTypeIndx = s->cTypeIndx;

	d->psolv = s->psolv;
	V32V3(s->rc_solv, d->rc_solv)
	V32V3(s->rc0_solv, d->rc0_solv)

	return true;
}

bool cp(CLUSTER *d, CLUSTER *s) {
	if (!cp((BASIC_CLUSTER*)d, (BASIC_CLUSTER*)s)) return false; 
	d->M = s->M;
	d->cm0 = s->cm0;
	d->cm = s->cm;

	d->r_solve = s->r_solve;
	d->ita_Brownian = s->ita_Brownian;
	d->cd_rot = s->cd_rot;
	d->cd_trans = s->cd_trans;

	return true;
}

// torsion angle is viewing from c2 to c1
double torsion_angle(BASIC_CLUSTER *c1, BASIC_CLUSTER *c2) { // c1 & c2 could not has parent-child relationship, but linked
	double tangle = 0;
	if (c1->atom == NULL || c2->atom == NULL) return tangle;
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge1 = NULL, *hinge2 = NULL;
	int i = 0;
	for (i = 0; i < c2->nHinges; i++) {
		if (c2->hinge[i].plink == (_BASIC_CLUSTER<MATOM>*)c1) {hinge2 = c2->hinge + i; break;}
	}
	hinge1 = c1->hinge + hinge2->indx_link;
	if (hinge1 == NULL || hinge2 == NULL) return 0;

	VECTOR3 r1, r2, r3, r4;

	MATOM *patom2 = c1->atom + hinge1->indx; // the atom where hinge is connected in c1
	MATOM *patom3 = c2->atom + hinge2->indx; // the atom where hinge is connected in c2

	if (patom2 == NULL || patom3 == NULL) {show_log("hige-atoms are null. can not calculate torsion angle.", true); return tangle;}

	V32V3(patom2->r, r2) V32V3(patom3->r, r3)

	MATOM *patom1 = (MATOM*)patom2->ba[0]; // the first atom connecting to patom2 in the cluster

	if (c1->nAtoms == 1) {
		for (i = 0; i < c1->nHinges; i++) {
			if (i != hinge2->indx_link) break; // the first hinge which is not connected to c2
		}
		if (i == c1->nHinges) {
			show_log("connected atom in cluster has only one bound. can not calculate torsion angle", true); return tangle;
		}
		V3plusV3(patom2->r, c1->hinge[i].vHinge, r1)
	}
	else {
		if (patom1 == patom3) {
			show_log("STRANGE : defined bound atom is in another cluster. calculated torsion angle is not right.", true);
			return tangle;
		}
		else {V32V3(patom1->r, r1)} // since c1 has more than 1 atom, so there must be a defined atom in the c1 connected with patom2
	}

	MATOM *patom4 = (MATOM*)patom3->ba[0]; // the first atom connecting to patom3 in the cluster
	if (c2->nAtoms == 1) {
		for (i = 0; i < c2->nHinges; i++) {
			if (i != hinge1->indx_link) break; // the first hinge which is not connected to c2
		}
		if (i == c2->nHinges) {
			show_log("connected atom in cluster has only one bound. can not calculate torsion angle", true); return tangle;
		}
		V3plusV3(patom3->r, c2->hinge[i].vHinge, r4)
	}
	else {
		if (patom4 == patom2) {
			show_log("STRANGE : defined bound atom is in another cluster. calculated torsion angle is not right.", true);
			return tangle;
		}
		else {V32V3(patom4->r, r4)} // since c1 has more than 1 atom, so there must be a defined atom in the c1 connected with patom2
	}

	// nornally, there must be atoms connecting to both side of hinge
	// since the atom on hinge must has more than 1 bound,
	// the atom inside cluster, or the atom in connecting cluster
	/*
	if (patom1 == NULL || patom4 == NULL) {
		show_log("atoms on hinge are null. can not calculate the torsion angle.", true);
		return tangle;
	}

	for (i = 0; i < 3; i++) {
		r1.v[i] = patom1->r.v[i];// + c1->dr_cmm.v[i];
		r2.v[i] = patom2->r.v[i];// + c1->dr_cmm.v[i];
		r3.v[i] = patom3->r.v[i];// + c2->dr_cmm.v[i];
		r4.v[i] = patom4->r.v[i];// + c2->dr_cmm.v[i];
	}
	*/
	tangle = TorsionAngle(r4, r3, r2, r1);
	return tangle;
}

double IdeaTorsionAngleFromSymetry(MATOM *patom) {
	int nfold = patom->nb - 1; 
	double tangle = 0;
	switch (nfold) {
	case 0:
		tangle = PI; break;
	default:
		tangle = PI / nfold; break;
	}
	return tangle;
}

double IdeaTorsionAngleFromSymetry(BASIC_CLUSTER *c1, BASIC_CLUSTER *c2) {
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge2 = NULL;
	for (int i = 0; i < c2->nHinges; i++) {
		if (c2->hinge[i].plink == c1) {hinge2 = c2->hinge + i; break;}
	}
	if (hinge2 == NULL) return 0;
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge1 = c1->hinge + hinge2->indx_link;
	MATOM *patom2 = c2->atom + hinge2->indx;
	MATOM *patom1 = c1->atom + hinge1->indx;
	double tangle = 0;
	if (patom1->nb > patom2->nb) tangle = IdeaTorsionAngleFromSymetry(patom1);
	else tangle = IdeaTorsionAngleFromSymetry(patom2);
	return tangle;
}

// connect c2 [child] to c1 [parent] with vcHinge [in c1]
// align vpHinge in c2 parallel to vcHinge in c1
bool ConnectClusters(CLUSTER *c1, int iHinge1, CLUSTER *c2, int iHinge2, char bound_type) {
	if (iHinge1 >= c1->nHinges) {
		sprintf(errmsg, "hinge index %d does not exist [HINGES : %d]", iHinge1, c1->nHinges); show_msg(errmsg); return false;
	}
	if (iHinge2 >= c2->nHinges) {
		sprintf(errmsg, "hinge index %d does not exist [HINGES : %d]", iHinge2, c2->nHinges); show_msg(errmsg); return false;
	}

	HINGE<_BASIC_CLUSTER<MATOM> > *hinge1 = c1->hinge + iHinge1;
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge2 = c2->hinge + iHinge2;
	MATOM *patom2 = c2->atom + hinge2->indx;
	MATOM *patom1 = c1->atom + hinge1->indx;

	//float hinge_length = (c1->hinge[iHinge1].vHinge.abs() + c2->hinge[iHinge2].vHinge.abs()) / 2;
	float hinge_length = 0;
	if (!find_bound_pars(force_field_fname, patom1->par->atom, patom2->par->atom, hinge_length, ::errmsg)) {
		show_msg(::errmsg); return false;
	}

	c1->hinge[iHinge1].vHinge.stretch(hinge_length);
	c2->hinge[iHinge2].vHinge.stretch(hinge_length);

	VECTOR3 O2 = (*(c2->hinge[iHinge2].O));
	VECTOR3 Op2 = c2->hinge[iHinge2].vHinge + O2;
	VECTOR3 O1 = (*(c1->hinge[iHinge1].O));
	c2->align(Op2, O2, O1, c1->hinge[iHinge1].vHinge);

	if (!c1->setup_link_relationship(iHinge1, c2, iHinge2, bound_type)) return false;

	c2->km.ta = torsion_angle(c1, c2);

	return true;
}

void shift_cluster_and_children(CLUSTER *c, VECTOR3 &d) {
	if (d.v[0] == 0 && d.v[1] == 0 && d.v[2] == 0) return;
	VECTOR3 dr = d;

	int i = 0, nc = 0;
	LIST<CLUSTER> *plist = NULL, *clist = NULL, *child_list = NULL;
	CLUSTER *pc = NULL, *child = NULL;

	clist = new LIST<CLUSTER>; clist->p = c;
	plist = clist;
	while (plist != NULL) {
		pc = plist->p;
		for (nc = 0; nc < pc->nHinges; nc++) {
			if (nc == pc->n2Parent) continue; // this is the hinge to parent
			child = (CLUSTER*)(pc->hinge[nc].plink);
			if (child == NULL) continue;
			// add child to child_list
			if (child_list == NULL) {
				child_list = new LIST<CLUSTER>; 
				child_list->p = child;
			}
			else child_list->add(child, false);
		}
		// shift the cluster only
		pc->shift(dr);

		plist = plist->next;

		if (plist == NULL) { // this level is done, go to child level
			release_list<CLUSTER>(&clist, false);
			clist = child_list; child_list = NULL;
			plist = clist;
		}
	}
	if (clist != NULL) release_list<CLUSTER>(&clist, false);
}

// shift the cluster and its children with given shift
void MMOLECULE::shift(CLUSTER *c, VECTOR3 &d) {
	return shift_cluster_and_children(c, d);
}

void tweak_cluster_and_children(CLUSTER *c, VECTOR3 &v0, R &r) {
	VECTOR3 R0 = v0;
	VECTOR3 invR0 = R0 * (-1);

	int i = 0, nc = 0;
	LIST<CLUSTER> *plist = NULL, *clist = NULL, *child_list = NULL;
	CLUSTER *pc = NULL, *child = NULL;

	clist = new LIST<CLUSTER>; clist->p = c;
	plist = clist;
	while (plist != NULL) {
		pc = plist->p;
		for (nc = 0; nc < pc->nHinges; nc++) {
			if (nc == pc->n2Parent) continue; // this is the hinge to parent
			child = (CLUSTER*)(pc->hinge[nc].plink);
			if (child == NULL) continue;
			// add child to child_list
			if (child_list == NULL) {
				child_list = new LIST<CLUSTER>; 
				child_list->p = child;
			}
			else child_list->add(child, false);
		}
		// tweak the cluster only
		pc->shift(invR0);
		pc->rotate(r);
		pc->shift(R0);

		plist = plist->next;

		if (plist == NULL) { // this level is done, go to child level
			release_list<CLUSTER>(&clist, false);
			clist = child_list; child_list = NULL;
			plist = clist;
		}
	}
	if (clist != NULL) release_list<CLUSTER>(&clist, false);
}

void tweak_cluster_and_children(CLUSTER *c, VECTOR3 &v0, VECTOR3 &axis, double omega, bool arc) {
	if (FABS(omega) < MIN_ROT) return;
	R r;
	RMatrix(axis, omega, r, arc);
	tweak_cluster_and_children(c, v0, r);
}

void shift_tweak_cluster_and_children(CLUSTER *c, VECTOR3 &dv, VECTOR3 &v0, R &r) {
	VECTOR3 R0 = v0;
	VECTOR3 ds = dv;

	int i = 0, nc = 0;
	LIST<CLUSTER> *plist = NULL, *clist = NULL, *child_list = NULL;
	CLUSTER *pc = NULL, *child = NULL;

	clist = new LIST<CLUSTER>; clist->p = c;
	plist = clist;
	while (plist != NULL) {
		pc = plist->p;
		for (nc = 0; nc < pc->nHinges; nc++) {
			if (nc == pc->n2Parent) continue; // this is the hinge to parent
			child = (CLUSTER*)(pc->hinge[nc].plink);
			if (child == NULL) continue;
			// add child to child_list
			if (child_list == NULL) {
				child_list = new LIST<CLUSTER>; 
				child_list->p = child;
			}
			else child_list->add(child, false);
		}
		// tweak the cluster only
		pc->shift_rotate(ds, R0, r);

		plist = plist->next;

		if (plist == NULL) { // this level is done, go to child level
			release_list<CLUSTER>(&clist, false);
			clist = child_list; child_list = NULL;
			plist = clist;
		}
	}
	if (clist != NULL) release_list<CLUSTER>(&clist, false);
}

void shift_tweak_cluster_and_children(CLUSTER *c, VECTOR3 &dv, VECTOR3 &v0, VECTOR3 &axis, double omega, bool arc) {
	R r;
	RMatrix(axis, omega, r, arc);
	return shift_tweak_cluster_and_children(c, dv, v0, r);
}

// tweak cluster and its children with given shift and rotation matrix
void MMOLECULE::tw(CLUSTER *c, VECTOR3 &v0, R &r) {
	return tweak_cluster_and_children(c, v0, r);
}

// twist nth cluster along its vpHinge with domega
void MMOLECULE::twTA(int n, double domega, bool arc) {
	if (domega == 0) return;
	if (n < 0 || n >= nCluster) return;
	if (cluster[n].n2Parent < 0) return; // parent hinge is not defined even
	if (cluster[n].hinge[cluster[n].n2Parent].O == NULL) return; // hinge was not set-up

	R r;
	VECTOR3 axis = cluster[n].hinge[cluster[n].n2Parent].vHinge;
	VECTOR3 O = (*(cluster[n].hinge[cluster[n].n2Parent].O)); // same as (*(cluster[n].hinge[cluster[n].n2Parent].Olink));

	RMatrix(axis, domega, r, arc);
	tweak_cluster_and_children(cluster + n, O, r);
}

void MMOLECULE::resetDynPars() {
	int i = 0, j = 0, nc = 0;
	for (nc = 0; nc < nCluster; nc++) {
		V6zero(cluster[nc].dyn->fc)
		V6zero(cluster[nc].invNE->a)
		V6zero(cluster[nc].invNE->b)
		//V6zero(cluster[nc].invNE->f)
		M6zero(cluster[nc].invNE->P);
		cluster[nc].invNE->T = 0;
		cluster[nc].invNE->gamma = 0;
	}
}

void setup_parent_child_relationship(CLUSTER *parent, CLUSTER *child) {
	int n;
	for (n = 0; n < parent->nHinges; n++) {
		if (n == parent->n2Parent) continue; // this is hinge to parent
		if (parent->hinge[n].plink == child) break;
	}
	if (n >= parent->nHinges) return;
	child->n2Parent = parent->hinge[n].indx_link;
	// setup parent, Op position defined in cluster
	// these two parameters are redundent
	child->SetupParent();
}

void setup_parent_child_relationship(CLUSTER *base, bool linked) {
	int i = 0, nc = 0;
	double r = 0;
	LIST<CLUSTER> *plist = NULL, *clist = NULL, *child_list = NULL;
	CLUSTER *pc = NULL, *child = NULL;

	base->SetupParent();
	clist = new LIST<CLUSTER>; clist->p = base;
	plist = clist;
	while (plist != NULL) {
		pc = plist->p;
		for (nc = 0; nc < pc->nHinges; nc++) {
			if (nc == pc->n2Parent) continue; // this is the hinge to parent
			child = (CLUSTER*)(pc->hinge[nc].plink);
			if (child == NULL) {
				if (linked) {
					sprintf(errmsg, "ERROR: hinge is free in linked cluster. Cluster chain is not constructed properly!");
					show_infor(errmsg); strcpy(errmsg, "\0");
				}
				continue;
			}
			
			// indx_link is the hinge index in linked cluster to this cluster
			child->n2Parent = pc->hinge[nc].indx_link;
			// setup parent, Op position defined in cluster
			// these two parameters are redundent
			child->SetupParent();
			
			// setup up direction
			if (child->parent != NULL) {
				r = child->hinge[child->n2Parent].vHinge.abs();
				if (r > 1e-6) {for (i = 0; i < 3; i++) child->up.v[i] = child->hinge[child->n2Parent].vHinge.v[i] / r;}
			}
			else {for (i = 0; i < 3; i++) child->up.v[i] = 0;}

			// add child to child_list
			if (child_list == NULL) {
				child_list = new LIST<CLUSTER>; 
				child_list->p = child;
			}
			else child_list->add(child, false);
		}
		// do something on pc ?
		/***** here nothing has to be done *****/
		// end of the work

		plist = plist->next;

		if (plist == NULL) { // this level is done, go to child level
			release_list<CLUSTER>(&clist, false);
			clist = child_list; child_list = NULL;
			plist = clist;
		}
	}
	if (clist != NULL) release_list<CLUSTER>(&clist, false);
}

bool check_parent_child_setup(CLUSTER *cluster, int nc) {
	bool status = true;
	int i, n, nbase = 0, nsingle = 0;
	char msg[256] = "\0";
	for (n = 0; n < nc; n++) {
		if (cluster[n].n2Parent < 0) nbase++;
		if (cluster[n].nHinges == 0) nsingle++;

		if (cluster[n].n2Parent >= 0) {
			if (cluster[n].parent == NULL || cluster[n].Op == NULL) {
				sprintf(errmsg, "ERROR : parent cluster or Op was not set properly");
				show_infor(errmsg); status = false;
			}
		}

		for (i = 0; i < cluster[n].nHinges; i++) {
			if (cluster[n].hinge[i].O == NULL || cluster[n].hinge[i].plink == NULL || cluster[n].hinge[i].Olink == NULL) {
				sprintf(errmsg, "ERROR : hinge was not setup or connected properly");
				show_infor(errmsg); status = false;
			}
		}
	}
	if (nbase == 0) {
		sprintf(errmsg, "ERROR : no base cluster is found. This is a closed chain.");
		show_infor(errmsg);
		status = false;
	}
	else if (nbase > 1) {
		sprintf(errmsg, "ERROR : %d base clusters are found.\nparent-child relationship between clusters are defined properly.", nbase);
		show_infor(errmsg);
		status = false;
	}
	if (nsingle > 0) {
		sprintf(errmsg, "ERROR : %d stand-alone clusters are found. cluster-chain is not defined properly.", nsingle);
		show_infor(errmsg);
		status = false;
	}
	return status;
}

void MMOLECULE::checkMolecule(bool tree) {
	int i = 0, j = 0;
	for (i = 0; i < nCluster; i++) cluster[i].mType = 0; // macromolecule
	reset_tips();
	for (i = 0; i < nCluster; i++) {
		if (cluster[i].nHinges == 1 && cluster[i].n2Parent >= 0) ntips++;
	}
	if (ntips == 0) {
		sprintf(errmsg, "No tip is found in the macromolecule. This is not consistent.");
		//show_msg(errmsg);
		show_infor(errmsg);
		return;
	}
	tip = new CLUSTER*[ntips];
	j = 0;
	for (i = 0; i < nCluster; i++) {
		if (cluster[i].nHinges == 1 && cluster[i].n2Parent >= 0) {
			tip[j] = cluster + i;
			j++;
		}
	}
	if (tree) construct_tree();

	LIST< LIST<CLUSTER> >*ptl = this->ctree.tl;
	LIST<CLUSTER> *plist = NULL;
	int nc = 0;
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			if (plist->p == NULL) {show_infor("ERROR: in the cluster-tree of macromol, an empty link exists");}
			else nc++; 
			plist = plist->next;
		}
		ptl = ptl->next;
	}
	if (nc != nCluster) {
		sprintf(errmsg, "number of cluster in tree : %d, macromol %s has cluster %d", nc, this->mol, nCluster);
		show_infor(errmsg);
	}
}

void MMOLECULE::construct_tree() {
	ctree.release_tree(false);
	ctree.new_level(); 
	ctree.add(ctree.tl, this->base, false);

	CLUSTER *pc = NULL;
	CLUSTER* child = NULL;
	LIST<CLUSTER> *parent_list = NULL;
	LIST< LIST<CLUSTER> > *ptl = NULL;
	int plevel = 0, i = 0;

	parent_list = ctree.tl->p;
	while (parent_list != NULL) {
		if (ptl == NULL) ptl = ctree.new_level();
		pc = parent_list->p;
		for (i = 0; i < pc->nHinges; i++) {
			if (i == pc->n2Parent) continue; // this hinge connects to parent
			child = (CLUSTER*)(pc->hinge[i].plink);
			if (child != NULL) ctree.add(ptl, child, false);
		}
		parent_list = parent_list->next;

		if (parent_list == NULL) { // end of this level
			if (ptl->p == NULL) detach< LIST<CLUSTER> >(&(ctree.tl), ptl, false);
			else {
				parent_list = ptl->p;
				ptl = NULL; // add child again
			}
		}
	}

	// construct a tree about the index of the cluster in ctree
	indx_tree.SetArray(ctree.nlevels);
	int nc = 0, j;
	LIST<CLUSTER> *plist = NULL;
	ptl = ctree.tl;
	for (i = 0; i < ctree.nlevels; i++) {
		plist = ptl->p;
		nc = plist->nList();
		indx_tree.m[i].SetArray(nc);
		j = 0;
		while (plist != NULL) {
			indx_tree.m[i].m[j] = plist->p->cmIndx;
			plist = plist->next; j++;
		}
		ptl = ptl->next; // next level of the tree
	}
	/*
	char msg[256] = "\0";
	sprintf(msg, "macromolecule %s -- cluster tree :", this->mol); show_log(msg, true);
	for (i = 0; i < indx_tree.n; i++) {
		sprintf(msg, "%d : ", i); show_log(msg, false);
		for (j = 0; j < indx_tree.m[i].n; j++) {
			sprintf(msg, "%d ", indx_tree.m[i].m[j]); show_log(msg, false);
		}
		show_log("\0", true);
	}
	show_log("---------------------- over ----------------------", true);
	*/
}

void MMOLECULE::CalMolInertiaMoment() {
	int nc, i, j, k;
	CLUSTER *pc = NULL;
	VECTOR<3> dr;
	MATRIX<6> m1, m2, phai, phai_T;

	// inertial momentum relative to base cluster
	M6zero(M0)
	for (nc = 0; nc < this->nCluster; nc++) {
		pc = cluster + nc;
		// the position of the origin of the cluster's inertial cordinate
		// relative to the base 
		//for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - base->Op->v[i];
		// dr from base's origin to cluster
		//for (i = 0; i < 3; i++) dr.v[i] = dr.v[i];

		//PHAI(phai, dr, i, j)
		//M6T2M6(phai, phai_T, i, j)
		MpM(pc->invNE->M, pc->invNE->phai_T_base2me, m1, i, j, k, 6)
		MpM(pc->invNE->phai_base2me, m1, m2, i, j, k, 6) // m2 = phai(base=>cluster) * M * phai_T(base=>cluster)
		MplusM(M0, m2, M0, i, j)
		/*
		{
			VECTOR<3> r;
			Matrix2Vect(m2.m, 0, 3, r.v);
			for (i = 0; i < 3; i++) r.v[i] /= pc->M;
			VECTOR<3> r1;
			for (i = 0; i < 3; i++) r1.v[i] = dr.v[i] + pc->cm.v[i];
			i = i;
		}
		*/
	}

	/*
	{
		VECTOR<3> r;
		Matrix2Vect(M0.m, 0, 3, r.v);
		for (i = 0; i < 3; i++) r.v[i] /= mDynKin.M;
		VECTOR<3> r1;
		for (i = 0; i < 3; i++) r1.v[i] = mDynKin.cm.v[i] - base->Op->v[i];
		i = i;
	}
	*/

	InverseMatrix<6>(M0, invM0); // m1 will be changed

	// inertial momentum relative to the mass center
	M6zero(mDynKin.I)
	
	// dr from molecular mass center to base cluster
	for (i = 0; i < 3; i++) dr.v[i] = r0.v[i];

	PHAI(phai, dr, i, j)
	M6T2M6(phai, phai_T, i, j)
	MpM(M0, phai_T, m1, i, j, k, 6)
	MpM(phai, m1, mDynKin.I, i, j, k, 6)

	InverseMatrix<6>(mDynKin.I, mDynKin.invI); // m1 will be changed
}

void MMOLECULE::CalMassCenterVelocity_MassCenterMomentum() {
	// calculate the velocity of mass center
	MpV6(mDynKin.invI, mDynKin.P, mDynKin.V)
}

void MMOLECULE::SetBaseForce() {
	int i, j;
	SVECTOR<double, 3> dr;
	MATRIX<6> phai;

	// r0 is the base cluster's position relative to the mass center 
	// dr is from base cluster to molecular mass center
	for (i = 0; i < 3; i++) dr.v[i] = -r0.v[i];

	PHAI(phai, dr, i, j)

	// calculate the molecular acceleration moment with cordinate at base
	MpV6(phai, mDynKin.f, f)
	
	MpV6(phai, mDynKin.b, b)
}

void SumGyroscopicSpatialForce(MMOLECULE *mm, int n1, int n2, VECTOR<6> &res) {
	V6zero(res)
	VECTOR<6> sum, ft;
	CLUSTER *pc = NULL;
	int nc = 0;
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		// dr is from mass center to cluster
		//for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - mDynKin.cm.v[i];
		//PHAI(phai, dr, i, j)
		// transfer GyroscopicSpatialForce on cluster to the mass center's cordinate
		MpV6(pc->invNE->phai_cm2me, pc->invNE->b, ft)
		V6plusV6(sum, ft, sum)
	}
	V62V6(sum, res)
	//return res;
}

void MMOLECULE::CalGyroscopicSpatialForce() {
	int nc = 0;
	CLUSTER *pc = NULL;
	//VECTOR3 dr;
	VECTOR<6> ft;
	//MATRIX<6> phai;
	V6zero(mDynKin.b)
	if (nCluster > NCLUSTERS_PARALLEL) {
		MacroMoleculeOperate2<MMOLECULE, VECTOR<6> >((void*)(&SumGyroscopicSpatialForce), (*this), 0, nCluster, MAX_THREADS, mDynKin.b);
		return;
	}
	else {
		for (nc = 0; nc < nCluster; nc++) {
			pc = cluster + nc;
			// dr is from mass center to cluster
			//for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - mDynKin.cm.v[i];
			//PHAI(phai, dr, i, j)
			// transfer GyroscopicSpatialForce on cluster to the mass center's cordinate
			MpV6(pc->invNE->phai_cm2me, pc->invNE->b, ft)
			V6plusV6(mDynKin.b, ft, mDynKin.b)
		}
	}
}

void MMOLECULE::CalFreeBaseVelocity(VECTOR<6> &Iw0_base) {
	VECTOR<6> dIw, dv;
	//for (i = 0; i < 6; i++) dIw.v[i] = Iw.v[i] - Iw0.v[i];
	V6minusV6(Iw0_base, Iw, dIw)
	MpV6(invM0, dIw, dv)
	//for (i = 0; i < 6; i++) this->V0.v[i] += dv.v[i];
	V6plusV6(this->V0, dv, this->V0)
}

void MMOLECULE::CalFreeBaseVelocity_GivenMassCenterMomentum() {
	int i, j;
	VECTOR<6> Iw1, dIw, dv;
	VECTOR3 r_base2mc;
	r_base2mc.v[0] = -r0.v[0]; // r0 : mass center to base cluster
	r_base2mc.v[1] = -r0.v[1]; // r0 : mass center to base cluster
	r_base2mc.v[2] = -r0.v[2]; // r0 : mass center to base cluster
	MATRIX<6> phai;
	// mm.r0 is from molecular mass center to base cluster
	PHAI(phai, r_base2mc, i, j)
	MpV6(phai, mDynKin.P, Iw1)
	V6minusV6(Iw1, Iw, dIw)
	MpV6(invM0, dIw, dv)
	V6plusV6(V0, dv, V0)

	//Iw = Iw1;
	memcpy(Iw.v, Iw1.v, SIZE_V6); //update right Iw
}

void MMOLECULE::twMol(VECTOR3 &vr0, VECTOR3 &axis, double domega, bool arc) {
	if (FABS(domega) < 1e-10) return;
	int i;
	double raxis = 0;
	V3ABS2(axis, raxis)
	if (raxis < 1e-10) return; // axis is [0, 0, 0]

	VECTOR3 nr0, r0;
	for (i = 0; i < 3; i++) {nr0.v[i] = -(vr0.v[i]); r0.v[i] = vr0.v[i];}
	shiftMol(nr0);

	R rot;
	RMatrix(axis, domega, rot, arc);

	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *pc = NULL;

	ptl = this->ctree.tl;
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;

			pc->rotate(rot);

			plist = plist->next;
		}
		ptl = ptl->next;
	}
	VECTOR<3> r;
	V32V3(mDynKin.cm, r)
	MpV3(rot, r, mDynKin.cm)

	V32V3(zaxis, r)
	MpV3(rot, r, zaxis)

	V32V3(xaxis, r)
	MpV3(rot, r, xaxis)

	shiftMol(r0);
}

// twist molecule with drot
void MMOLECULE::twMol(VECTOR3 &vr0, VECTOR3 &drot, bool arc) {
	VECTOR3 axis;
	double domega = 0;
	int i;
	for (i = 0; i < 3; i++) domega += drot.v[i] * drot.v[i];
	domega = sqrt(domega);
	if (domega < MIN_ROT) return;
	for (i = 0; i < 3; i++) axis.v[i] = drot.v[i] / domega;
	twMol(vr0, axis, domega, arc);
}

// shift the whole molecule
void MMOLECULE::shiftMol(VECTOR3 &d) {
	if (d.v[0] == 0 && d.v[1] == 0 && d.v[2] == 0) return;
	VECTOR3 dr = d;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	int nc = 0, na;

	for (nc = 0; nc < nCluster; nc++) {
		pc = cluster + nc;
		//pc->shift(dr);
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			V3plusV3(patom->r, dr, patom->r)
		}
		V3plusV3(pc->cm, dr, pc->cm)
		pc->update_vp();
	}
	V3plusV3(mDynKin.cm, dr, mDynKin.cm)
}

void MMOLECULE::shift_tw_mol(VECTOR3 &dv, VECTOR3 &v0, VECTOR3 &axis, double domega, bool arc) {
	int i;
	VECTOR3 O, ds;
	for (i = 0; i < 3; i++) {O.v[i] = v0.v[i]; ds.v[i] = dv.v[i];}

	R rot;
	RMatrix(axis, domega, rot, arc);

	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *pc = NULL;

	ptl = this->ctree.tl;
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;

			pc->shift_rotate(ds, O, rot);

			plist = plist->next;
		}
		ptl = ptl->next;
	}

	VECTOR<3> r;
	V32V3(mDynKin.cm, r)
	MpV3(rot, r, mDynKin.cm)

	V32V3(zaxis, r)
	MpV3(rot, r, zaxis)

	V32V3(xaxis, r)
	MpV3(rot, r, xaxis)
}

// twist molecule with dw
void MMOLECULE::shift_tw_mol(VECTOR3 &dv, VECTOR3 &v0, VECTOR3 &drot, bool arc) {
	VECTOR3 axis;
	double domega = 0;
	int i;
	for (i = 0; i < 3; i++) domega += drot.v[i] * drot.v[i];
	domega = sqrt(domega);
	if (domega < MIN_ROT) return;
	for (i = 0; i < 3; i++) axis.v[i] = drot.v[i] / domega;
	shift_tw_mol(dv, v0, axis, domega, arc);
}

void MMOLECULE::calMassCenter(bool cal_cluster_mc) {
	int i, nc;
	VECTOR<3> cm;

	mDynKin.M = 0;
	V3zero(mDynKin.cm)
	for (nc = 0; nc < nCluster; nc++) {
		if (cal_cluster_mc) cluster[nc].MassCenter();
		mDynKin.M += cluster[nc].M;
		mDynKin.cm.v[0] += cluster[nc].M * cluster[nc].cm.v[0];
		mDynKin.cm.v[1] += cluster[nc].M * cluster[nc].cm.v[1];
		mDynKin.cm.v[2] += cluster[nc].M * cluster[nc].cm.v[2];
	}
	for (i = 0; i < 3; i++) mDynKin.cm.v[i] /= mDynKin.M;
	V32V3(mDynKin.cm, mDynKin.r0)

	// r0 is from molecular mass center to base cluster
	for (i = 0; i < 3; i++) r0.v[i] = base->Op->v[i] - mDynKin.cm.v[i];
}

void MMOLECULE::SetMassCenter(double *v, bool cal_mc, bool cal_cluster_mc) {
	int i;

	if (cal_mc) calMassCenter(cal_cluster_mc);
	VECTOR3 dr;
	for (i = 0; i < 3; i++) dr.v[i] = v[i] - mDynKin.cm.v[i];
	shiftMol(dr);
}

void MMOLECULE::SetBasePosZero() {
	int i;
	VECTOR3 dr;
	for (i = 0; i < 3; i++) dr.v[i] = -(base->Op->v[i]);
	shiftMol(dr);
}

void MMOLECULE::setMolIndx(int indx) {
	int nc = 0;
	this->mIndx = indx;
	for (nc = 0; nc < this->nCluster; nc++) cluster[nc].mIndx = indx;
}

void calClusterVelocity0(MMOLECULE &mm, bool zeroBaseVelocity) {
	int i = 0, j = 0, k = 0;
	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *p = NULL, *parent = NULL;

	ptl = mm.ctree.tl;
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			p = plist->p;
			parent = (CLUSTER*)(p->parent);
			if (parent != NULL) {MpV6(p->km.phai_T, parent->bkm->V0, p->bkm->V0)} //p->km.V0 = p->km.phai * parent->km.V;
			else {
				//for (i = 0; i < 6; i++) p->bkm->V0.v[i] = mm.V0.v[i];
				if (zeroBaseVelocity) { memset(p->bkm->V0.v, 0, SIZE_V6); }
				else {V62V6(mm.V0, p->bkm->V0)}
			}
			//V3plusV3(p->bkm->V0, p->km.V, p->bkm->V0)
			for (i = 0; i < 3; i++) p->bkm->V0.v[i] += p->up.v[i] * p->km.tp;
			plist = plist->next;
		}
		ptl = ptl->next;
	}
}

void calClusterAlpha(MMOLECULE &mm) {
	int i = 0, j = 0, k = 0;
	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *p = NULL;

	ptl = mm.ctree.tl;
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			p = plist->p;
			if (mm.free_base && p == mm.base) {
				//for (i = 0; i < 6; i++) p->bkm->alpha.v[i] = mm.alpha0.v[i];
				V62V6(mm.alpha0, p->bkm->alpha)
			}
			else {
				MpV6(p->km.phai_T, ((BASIC_CLUSTER*)p->parent)->bkm->alpha, p->bkm->alpha)
				for (i = 0; i < 6; i++) p->bkm->alpha.v[i] += p->km.tpp * p->km.H.v[i] + p->invNE->a.v[i]; // H[4,5,6] = 0
			}
			plist = plist->next;
		}
		ptl = ptl->next;
	}
}

 // return energy in kT
double calKineticEnergy(MMOLECULE &mm) {
	double Ek = 0;
	if (mm.nCluster > NCLUSTERS_PARALLEL) Ek = ParallelCalKineticEnergy(mm);
	else calClustersKineticEnergy(&mm, 0, mm.nCluster, Ek);
	mm.Ek = Ek;
	return Ek;
}

 // return kinetic energy in kT, of mass center
double calFrameKineticEnergy(MMOLECULE &mm) {
	double Ek = 0;
	VECTOR<6> Iw;
	int i = 0, j = 0;

	MpV6(mm.mDynKin.I, mm.mDynKin.V, Iw) // Iw = I * V0;
	V62V6(Iw, mm.mDynKin.P)
	//scale_uv6(mm.mDynKin.V, Iw, Ek) // Ek = V0 * I * V0;
	//Ek /= 2;
	//return Ek * ::unit_Ek_kT;

	mm.Ek_frame_rot = (Iw.v[0] * mm.mDynKin.V.v[0] + Iw.v[1] * mm.mDynKin.V.v[1] + Iw.v[2] * mm.mDynKin.V.v[2]) * 0.5 * ::unit_Ek_kT;
	mm.Ek_frame_trans = (Iw.v[3] * mm.mDynKin.V.v[3] + Iw.v[4] * mm.mDynKin.V.v[4] + Iw.v[5] * mm.mDynKin.V.v[5]) * 0.5 * ::unit_Ek_kT;
	mm.Ek_frame = mm.Ek_frame_rot + mm.Ek_frame_trans;

	return mm.Ek_frame;
}

// set the cluster, its children and grand-children's V0, by assuming its parent is null
void ClusterVelocity(CLUSTER *pc) {
	int i = 0, j = 0, k = 0;
	LIST<CLUSTER> *plist = NULL, *clist = NULL, *child_list = NULL;
	CLUSTER *p = NULL, *child = NULL;

	clist = new LIST<CLUSTER>; clist->p = pc;
	plist = clist;
	while (plist != NULL) {
		p = plist->p;
		for (i = 0; i < p->nHinges; i++) {
			if (i == p->n2Parent) continue; // this hinge connects to parent
			child = (CLUSTER*)(p->hinge[i].plink);
			if (child == NULL) continue;
			if (child_list == NULL) {
				child_list = new LIST<CLUSTER>; child_list->p = child;
			}
			else child_list->add(child, false);
		}
		if (p == pc) {}
		else {
			MpV6(p->km.phai_T, ((BASIC_CLUSTER*)p->parent)->bkm->V0, p->bkm->V0) //p->km.V0 = p->km.phai * p->parent->km.V;
			for (i = 0; i < 3; i++) p->bkm->V0.v[i] += p->km.tp * pc->km.H.v[i]; // H[4,5,6] = 0
		}
		plist = plist->next;
		
		if (plist == NULL) { // end of list, start child_list
			release_list<CLUSTER>(&clist, false);
			clist = child_list; child_list = NULL;
			plist = clist;
		}
	}
	if (clist != NULL) release_list<CLUSTER>(&clist, false);
}

// return energy in kT
// calculate kinetic energy of the cluster with its children + grand-children
double calClusterKineticEnergy(CLUSTER *pc) {
	double Ek = 0, Ekc = 0;
	VECTOR<6> Iw;
	
	int i = 0, j = 0, k = 0;
	LIST<CLUSTER> *plist = NULL, *clist = NULL, *child_list = NULL;
	CLUSTER *p = NULL, *child = NULL;

	clist = new LIST<CLUSTER>; clist->p = pc;
	plist = clist;
	while (plist != NULL) {
		p = plist->p;
		for (i = 0; i < p->nHinges; i++) {
			if (i == p->n2Parent) continue; // this hinge connects to parent
			child = (CLUSTER*)(p->hinge[i].plink);
			if (child == NULL) continue;
			if (child_list == NULL) {
				child_list = new LIST<CLUSTER>; child_list->p = child;
			}
			else child_list->add(child, false);
		}
			
		MpV6(p->invNE->M, p->bkm->V0, Iw) // Iw = M * V0;
		scale_uv6(p->bkm->V0, Iw, Ekc) // Ekc = V0 * M * V0;
		Ek += Ekc;

		plist = plist->next;
		
		if (plist == NULL) { // end of list, start child_list
			release_list<CLUSTER>(&clist, false);
			clist = child_list; child_list = NULL;
			plist = clist;
		}
	}
	if (clist != NULL) release_list<CLUSTER>(&clist, false);

	return Ek * ::unit_Ek_kT * 0.5;
}

// rotation of cluster is initilized based on the velocity of base cluster
double InitMolClusterVelocity(MMOLECULE &mm) {
	int i = 0, j = 0, k = 0;
	if (mm.ctree.tl == NULL) return 0;

	double Eth = 0.5; // Eth in unit kT
	double unit_E = kT / U_InertiaMoment;
	double tq = 0;
	VECTOR<3> torque, force;

	double E_total = Eth; // base has 0.5kT
	double tp1 = 0, tp2 = 0, tp = 0;
	double E0 = 0, E1 = 0, E2 = 0, Ek = 0, dE = 0;
	double dE1 = 0, dE2 = 0;
	int niter = 0;

	VECTOR<6> V0;
	MATRIX<6> phai, phai_T;
	VECTOR<3> dr;

	LIST< LIST<CLUSTER> > *ptl = mm.ctree.tl->get_tail();
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *pc = NULL;

	 // give a random direction to base up
	if (mm.free_base) rand_vect3(mm.base->up.v);
	//if (mm.free_base) mm.base->up.v[0] = 1; // [100]

	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;
			V2wv(torque, force, pc->dyn->fc)
			scale_uv3(torque, pc->up, tq)

			for (i = 0; i < 3; i++) {
				if (pc->parent == NULL) dr.v[i] = 0;
				else dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm.base->Op->v[i];
			}
			// dr from base to cluster
			PHAI(phai, dr, i, j)
			M6T2M6(phai, phai_T, i, j)
			// mm.V0 is the velocity of base, V0 is the velocity at cluster
			MpV6(phai_T, mm.V0, V0)

			if (NO_CHILD(pc)) {
				Ek = Eth;
				tp = 1e-6;
				pc->km.tp = (tq > 0 ? tp : -tp);
				for (i = 0; i < 6; i++) pc->bkm->V0.v[i] = V0.v[i];
				ClusterVelocity(pc);
				//Ek = calClusterKineticEnergy(pc);
			}
			else {
				pc->km.tp = 0; tp1 = 0;
				for (i = 0; i < 6; i++) pc->bkm->V0.v[i] = V0.v[i];
				ClusterVelocity(pc);
				E0 = calClusterKineticEnergy(pc);
				E1 = E0;

				tp2 = 1e-6; 
				tp2 = (tq > 0 ? tp2 : -tp2);
				pc->km.tp = tp2;
				for (i = 0; i < 6; i++) pc->bkm->V0.v[i] = V0.v[i];
				ClusterVelocity(pc);
				E2 = calClusterKineticEnergy(pc);

				dE = E2 - E1;
				tp = tp1 + (tp2 - tp1) / dE * (Eth + E0 - E1);
				
				niter = 0;
				Ek = E2;
				dE1 = E1 - E0 - Eth; dE1 = FABS(dE1);
				dE2 = E2 - E0 - Eth; dE2 = FABS(dE2);
				if (dE1 > dE2) {
					dE = E2; E2 = E1; E1 = dE;
					dE = tp2; tp2 = tp1; tp1 = dE;
				}
				while (FABS(Ek - Eth - E0) > 0.01f) {
					dE = E2 - E1;
					tp = tp1 + (tp2 - tp1) / dE * (Eth + E0 - E1);

					pc->km.tp = tp;
					for (i = 0; i < 6; i++) pc->bkm->V0.v[i] = V0.v[i];
					ClusterVelocity(pc);
					Ek = calClusterKineticEnergy(pc);
					dE = Ek - E0 - Eth; dE = FABS(dE);
					dE1 = E1 - E0 - Eth; dE1 = FABS(dE1);
					dE2 = E2 - E0 - Eth; dE2 = FABS(dE2);
					if (dE < dE1) {E2 = E1; E1 = Ek; tp2 = tp1; tp1 = tp;}
					else {E2 = Ek; tp2 = tp;}
					
					niter++;
					if (niter == 8) break;
				}
			}
			plist = plist->next;
		}
		ptl = ptl->pre;
	}
	double tp_avg = 0;
	int nc = 0;
	for (nc = 0; nc < mm.nCluster; nc++) {
		if (mm.cluster + nc == mm.base) continue;
		tp_avg = FABS(mm.cluster[nc].km.tp);
	}
	tp_avg /= (mm.nCluster - 1);
	for (nc = 0; nc < mm.nCluster; nc++) {
		if (mm.cluster + nc == mm.base) continue;
		mm.cluster[nc].km.tp = sign(mm.cluster[nc].km.tp) * tp_avg;
	}

	if (mm.free_base) {
		V62V6(mm.base->bkm->V0, mm.V0);
		mm.base->km.tp = 0;
		for (i = 0; i < 3; i++) mm.base->up.v[i] = 0; // no up
	}
	calClusterVelocity0(mm);
	calSpatialMomentum(mm); // momentum on base cluster

	// translation spatial moment = 0
	V6zero(mm.mDynKin.P);
	mm.CalMassCenterVelocity_MassCenterMomentum();
	mm.CalFreeBaseVelocity_GivenMassCenterMomentum();

	calClusterVelocity0(mm);

	// total spatial moment with cordinate at base cluster
	calSpatialMomentum0(mm, true); 

	calKineticEnergy(mm);
	calFrameKineticEnergy(mm);
	mm.Ek_internal = mm.Ek - mm.Ek_frame;

	return mm.Ek;
}

// Spatial Momentum is relative to the base cluster's cordinate
void calSpatialMomentum(MMOLECULE &mm) {
	if (mm.nCluster > NCLUSTERS_PARALLEL) return ParallelCalSpatialMomentum(mm);

	int i = 0, j = 0, nc = 0;
	CLUSTER *pc = NULL;
	VECTOR<6> iw, iw0;

	//VECTOR<3> dr;
	//MATRIX<6> phai;

	V6zero(mm.Iw);
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		MpV6(pc->invNE->M, pc->bkm->V0, iw)
		memcpy(pc->bkm->P0.v, iw.v, SIZE_V6); // update pc->bkm->P0

		//MpV6(pc->invNE->P, pc->bkm->V0, pc->bkm->P0) // update pc->bkm->P0
		
		// dr is from base to cluster
		/*
		for (i = 0; i < 3; i++) {
			dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm.base->Op->v[i];
		}
		PHAI(phai, dr, i, j)
		*/
		MpV6(pc->invNE->phai_base2me, iw, iw0)
		//for (i = 0; i < 6; i++) Iw.v[i] += iw0.v[i];
		V6plusV6(mm.Iw, iw0, mm.Iw)
	}
}

// Spatial Moment is relative to molecular Mass Center
void calSpatialMomentum0(MMOLECULE &mm, bool bCalBaseMomentum) {
	int i, j;
	if (bCalBaseMomentum) calSpatialMomentum(mm);

	MATRIX<6> phai;
	// mm.r0 is from molecular mass center to base cluster
	PHAI(phai, mm.r0, i, j)
	MpV6(phai, mm.Iw, mm.mDynKin.P)
}

// Spatial Acceleration Moment is relative to the base cluster's cordinate
void calSpatialAccel(MMOLECULE &mm, VECTOR<6> &Ialpha) {
	if (mm.nCluster > NCLUSTERS_PARALLEL) return ParallelCalSpatialAccel(mm, Ialpha);

	int i = 0, j = 0, nc = 0;
	CLUSTER *pc = NULL;
	VECTOR<6> iw, iw0;

	//VECTOR<3> dr;
	//MATRIX<6> phai;

	V6zero(Ialpha);
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		MpV6(pc->invNE->M, pc->bkm->alpha, iw)
		// dr is from base to cluster
		/*
		for (i = 0; i < 3; i++) {
			dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm.base->Op->v[i];
		}
		PHAI(phai, dr, i, j)
		*/
		MpV6(pc->invNE->phai_base2me, iw, iw0)
		//for (i = 0; i < 6; i++) Ialpha.v[i] += iw0.v[i];
		V6plusV6(Ialpha, iw0, Ialpha)
	}
}

// Spatial Acceleration Moment is relative to molecular Mass Center
void calSpatialAccel0(MMOLECULE &mm, VECTOR<6> &Ialpha) {
	int i, j;
	VECTOR<6> Iac_base;
	calSpatialAccel(mm, Iac_base);

	MATRIX<6> phai;
	// mm.r0 is from molecular mass center to base cluster
	PHAI(phai, mm.r0, i, j)
	MpV6(phai, Iac_base, Ialpha)
}

void ScaleVelocity(MMOLECULE &mm, double c) {
	int i = 0, j = 0, k = 0, nc = 0;

	mm.V0 *= c;
	for (nc = 0; nc < mm.nCluster; nc++) {
		mm.cluster[nc].km.tp *= c;
		mm.cluster[nc].bkm->V0 *= c;
		mm.cluster[nc].bkm->P0 *= c;
	}
	mm.mDynKin.V *= c;
	mm.mDynKin.P *= c;
}

double ScaleKineticEnergy(MMOLECULE &mm, double Ek0, bool reset_total_SpacialMoment) {
	int i = 0, j = 0, k = 0, nc = 0;

	if (mm.free_base && reset_total_SpacialMoment) {
		calSpatialMomentum(mm);
		V6zero(mm.mDynKin.P)
		mm.CalFreeBaseVelocity_GivenMassCenterMomentum();
	}

	double mdiff = 0, max_diff = 0.05; 
	double tp = 0, tp_max = 0, tp_min = 0, tp_av = 0;

	int nloops = 0;
	while (nloops < 1) {
		tp_max = 0; tp_min = 0; tp_av = 0;
		for (nc = 0; nc < mm.nCluster; nc++) {
			tp = FABS(mm.cluster[nc].km.tp);
			if (tp_min == 0) tp_min = tp;
			else if (tp_min > tp) tp_min = tp;
			if (tp_max == 0) tp_max = tp;
			else if (tp_max < tp) tp_max = tp;
			tp_av += tp;
		}
		tp_av /= (mm.nCluster - 1);
		if ((tp_max - tp_min) < tp_av / 5) {
			tp_av *= 0.8;
			tp_min = tp_av * 0.8;
			tp_max = tp_av * 1.05;
		}
		else {
			//tp_max = tp;
			tp_max = tp_min + (tp_max - tp_min) * 0.8;
		}
		for (nc = 0; nc < mm.nCluster; nc++) {
			if (mm.cluster[nc].parent == NULL) continue;
			if (FABS(mm.cluster[nc].km.tp) > tp_max) mm.cluster[nc].km.tp = sign(mm.cluster[nc].km.tp) * tp_max;
			/*
			for (i = 0; i < 3; i++) {
				mm.cluster[nc].km.V.v[i] = mm.cluster[nc].km.tp * mm.cluster[nc].up.v[i];
				mm.cluster[nc].km.V.v[i + 3] = 0;
			}
			*/
		}
		calClusterVelocity0(mm);

		if (mm.free_base) {
			calSpatialMomentum(mm);
			mm.CalFreeBaseVelocity_GivenMassCenterMomentum();
		}
		calClusterVelocity0(mm);
		// total spatial moment with cordinate at base cluster
		calSpatialMomentum0(mm); 
		mm.CalMassCenterVelocity_MassCenterMomentum();
		calKineticEnergy(mm);
		calFrameKineticEnergy(mm);
		mm.Ek_internal = mm.Ek - mm.Ek_internal;
		if (mm.Ek > 1) {
			ScaleVelocity(mm, sqrt(Ek0 / mm.Ek));
			calKineticEnergy(mm);
			calFrameKineticEnergy(mm);
			mm.Ek_internal = mm.Ek - mm.Ek_internal;
		}
		if (mm.Ek <= Ek0) break;
		nloops ++;
	}

	/*
	{
		char buffer[256] = "\0";
		sprintf(buffer, "kinetic energy is scaled to %f, expected to %f", Ek, Ek0);
		show_log(buffer, true);
	}
	*/

	return mm.Ek;
}

void ClusterMaxVelocity(MMOLECULE &mm, double &vmax, double &wmax) {
	double v = 0, w = 0;
	int i = 0, nc = 0;
	VECTOR<3> V, W;
	CLUSTER *p = NULL;
	for (nc = 0; nc < mm.nCluster; nc++) {
		p = mm.cluster + nc;
		V2wv(W, V, p->bkm->V0)
		w = W.v[0] * W.v[0] + W.v[1] * W.v[1] + W.v[2] * W.v[2];
		v = V.v[0] * V.v[0] + V.v[1] * V.v[1] + V.v[2] * V.v[2];
		if (nc == 0) { vmax = v; wmax = w; }
		else {
			vmax = (vmax < v ? v : vmax);
			wmax = (wmax < w ? w : wmax);
		}
	}
	vmax = sqrt(vmax);
}

bool cp(MMOLECULE *m2, MMOLECULE* m1, bool setup_tree) {
	int nc, i;
	m2->reset();
	m2->mType = m1->mType; m2->mIndx = m1->mIndx; strcpy(m2->mol, m1->mol);
	if (m1->nCluster <= 0) return true;
	m2->cluster = new CLUSTER[m1->nCluster];
	m2->nCluster = m1->nCluster;
	for (nc = 0; nc < m1->nCluster; nc++) cp(m2->cluster + nc, m1->cluster + nc);

	int nbase = 0, n = 0;
	for (n = 0; n < m1->nCluster; n++) if (m1->base == m1->cluster + n) {nbase = n; break;}
	
	// monomer's definition
	m2->monomer = new _GROUP<CLUSTER>[m1->nMonomer];
	m2->nMonomer = m1->nMonomer;
	for (n = 0; n < m2->nMonomer; n++) {
		m2->monomer[n].set(m1->monomer[n].nUnit);
		for (i = 0; i < m2->monomer[n].nUnit; i++) {
			// clusters are defined in the same way in m1 and m2
			m2->monomer[n].u[i] = m2->cluster + int(m1->monomer[n].u[i] - m1->cluster);
		}
		m2->monomer[n].uid = m1->monomer[n].uid;
	}
	
	// setup the relationship
	int iHinge1 = 0, iHinge2 = 0, nc2;
	CLUSTER *pc1 = NULL, *pc2 = NULL;
	for (nc = 0; nc < m2->nCluster; nc++) {
		pc1 = m2->cluster + nc;
		for (i = 0; i < m1->cluster[nc].nHinges; i++) {
			if (pc1->hinge[i].plink != NULL) continue; // the link was setup already
			if ((pc2 = (CLUSTER*)(m1->cluster[nc].hinge[i].plink)) == NULL) continue; // the link in original molecule is free
			for (nc2 = 0; nc2 < m1->nCluster; nc2++) {if (m1->cluster + nc2 == pc2) break;}
			pc2 = m2->cluster + nc2;
			iHinge1 = i;
			iHinge2 = m1->cluster[nc].hinge[i].indx_link;
			if (!(pc1->setup_link_relationship(iHinge1, pc2, iHinge2, SINGLE_BOUND))) {
				sprintf(errmsg, "failure to setup link between cluster %d & %d", nc, nc2);
				show_msg(errmsg); strcpy(errmsg, "\0"); 
				m2->reset();
				return false;
			}
		}
	}

	for (n = 0; n < 3; n++) {m2->zaxis.v[n] = m1->zaxis.v[n]; m2->xaxis.v[n] = m1->xaxis.v[n];}
	if (setup_tree && m1->base != NULL) {
		for (n = 0; n < m1->nCluster; n++) if (m1->base == m1->cluster + n) {nbase = n; break;}
		m2->base = m2->cluster + nbase;
		m2->base->disable_parent();
		setup_parent_child_relationship(m2->base, true); // report error when free hinge exist
		m2->checkMolecule(true);
		if (!check_parent_child_setup(m2->cluster, m2->nCluster)) {m2->reset(); return false;}
		if (!m2->linkPars(atompar_db, errmsg)) {m2->reset(); show_infor(errmsg); return false;}
		m2->calTorsionAngle();
		m2->setup_cluster_matom(); // require the atomic mass, and the parent-child relationship
	}
	return true;
}

bool atom_indx(MMOLECULE &mm, MATOM* patom, int &nc, int &na, int &indx) {
	nc = 0; na = 0; indx = 0;
	CLUSTER *pc = NULL;
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		for (na = 0; na < pc->nAtoms; na++) {
			if (pc->atom + na == patom) return true;
			indx++;
		}
	}
	return false;
}

extern char force_field_fname[256];

// FOLLOWING FUNCTIONS ARE WORKING WITH TORSION TO CALCULATE THE TORSION ON THE HINGE IN IMPLICIT WAY, WHICH IS WRONG
#if _DISABLE_ == 0
// ******************   calculating torsion, from here on  *******************
void copyHingeAtom2Cluster(CLUSTER *d, MATOM *patom1, MATOM *patom2) {
	//copy the atoms connected to the atoms on hinge
	d->setAtoms(patom1->nb);

	MATOM *patom = d->atom;
	cp(patom, patom1); // main atom
	int i = 0, n = 1; //atom0 in d is patom1
	for (i = 0; i < patom1->nb; i++) {
		if (patom1->ba[i] == patom2) continue; // ignore this atom
		cp(d->atom + n, (MATOM*)patom1->ba[i]);
		d->atom[n].set_bounds(1); //connect to atom 0 only
		bound_connect<MATOM>(d->atom + n, 0, d->atom, n - 1, patom1->bv[i], patom1->btype[i]);
		n++;
	}
	d->set_hinges(1);
	d->hinge[0].vHinge = patom2->r - patom1->r;
	d->setup_hinge(0, 0, patom->nb - 1);
}

// c1 & c2 could not has parent-child relationship, but linked, torsion angle is viewing from c1 to c2
void calTorsionEnergyTorque(BASIC_CLUSTER *c1, BASIC_CLUSTER *c2, FLEX_ASYM_MATRIX<TORSION_PARS> &tpar, double &E, double &torque) {
	E = 0, torque = 0;
	int n = 0;
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge1 = NULL;
	for (n = 0; n < c1->nHinges; n++) {
		if (c1->hinge[n].plink == c2) {hinge1 = c1->hinge + n; break;}
	}
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge2 = c2->hinge + hinge1->indx_link;
	MATOM *patom2 = c1->atom + hinge1->indx;
	MATOM *patom3 = c2->atom + hinge2->indx;
	MATOM *patom1 = NULL, *patom4 = NULL;
	double torsion_angle = 0, alpha = 0;

	TORSION_PARS *pTorPar = NULL;

	int na1, na4;
	for (na1 = 0; na1 < patom2->nb - 1; na1++) {
		patom1 = (MATOM*)patom2->ba[na1]; if (patom1 == patom3) continue;
		for (na4 = 0; na4 < patom3->nb - 1; na4++) {
			patom4 = (MATOM*)patom3->ba[na4]; if (patom4 == patom2) continue;
			torsion_angle = TorsionAngle(patom1->r, patom2->r, patom3->r, patom4->r);

			pTorPar = &(tpar.m[na1].m[na4]);
			while (pTorPar != NULL) {
				alpha = pTorPar->pn * torsion_angle - pTorPar->phase;
				E += pTorPar->v0 * (1 + cos(alpha)) * U_eps / kT; //kJ/mol ==> kT
				torque -= pTorPar->pn * pTorPar->v0 * sin(alpha) * U_LJForce / U_InertiaMoment; // kJ/mol ==> atom inertial moment
				pTorPar = pTorPar->next;
			}
		}
	}
	return;
}

// c1 & c2 hinge is on atoms #0 in c1 & c2, torsion angle is viewing from c1 to c2, with given torsion parameters
void calTorsionEnergyTorque1(BASIC_CLUSTER *c1, BASIC_CLUSTER *c2, FLEX_ASYM_MATRIX<TORSION_PARS> &tpar, double &E, double &torque) {
	//char msg[256] = "\0";
	E = 0, torque = 0;
	MATOM *patom2 = c1->atom;
	MATOM *patom3 = c2->atom;
	MATOM *patom1 = NULL, *patom4 = NULL;
	double torsion_angle = 0, alpha = 0;

	TORSION_PARS *pTorPar = NULL;

	int na1, na4, nsum = 0;
	for (na1 = 0; na1 < patom2->nb; na1++) {
		patom1 = (MATOM*)patom2->ba[na1]; if (patom1 == patom3) continue;
		for (na4 = 0; na4 < patom3->nb; na4++) {
			patom4 = (MATOM*)patom3->ba[na4]; if (patom4 == patom2) continue;
			torsion_angle = TorsionAngle(patom1->r, patom2->r, patom3->r, patom4->r);

			pTorPar = &(tpar.m[na1].m[na4]);
			while (pTorPar != NULL) {
				/*
				if (FABS(pTorPar->v0) > 0.0001) {
					sprintf(msg, "torsion angle = %f", torsion_angle);
					show_log(msg, true);
					nsum ++;
				}
				*/
				alpha = pTorPar->pn * torsion_angle - pTorPar->phase;
				E += pTorPar->v0 * (1 + cos(alpha)) * U_eps / kT; //kJ/mol ==> kT
				torque -= pTorPar->pn * pTorPar->v0 * sin(alpha) * U_LJForce / U_InertiaMoment; // kJ/mol ==> atom inertial moment
				pTorPar = pTorPar->next;
			}
		}
	}
	//sprintf(msg, "sum %d", nsum); show_log(msg, true);
	return;
}

// c1 & c2 could not has parent-child relationship, but linked, torsion angle is viewing from c2 to c1
void calTorsionProfile(BASIC_CLUSTER *pc1, BASIC_CLUSTER *pc2) {
	char msg[256] = "\0";
	CLUSTER c1, c2;
	int n1 = 0, n2 = 0, i;
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge1 = NULL;
	for (i = 0; i < pc1->nHinges; i++) {
		if (pc1->hinge[i].plink == pc2) {hinge1 = pc1->hinge + i; break;}
	}
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge2 = pc2->hinge + hinge1->indx_link;
	MATOM *patom1 = pc1->atom + hinge1->indx;
	MATOM *patom2 = pc2->atom + hinge2->indx;

	copyHingeAtom2Cluster(&c1, patom1, patom2);
	copyHingeAtom2Cluster(&c2, patom2, patom1);
	ConnectClusters(&c1, 0, &c2, 0, SINGLE_BOUND);
	setup_parent_child_relationship(&c1, true); // c1 is parent, c2 is child
	c2.km.ta = torsion_angle(&c1, &c2);

	TORSION *pTor = pc2->pTors;
	int npts = 90;
	double step = PI2 / npts;
	pTor->tProf.set(npts, step, -PI);
	//pTor->torqueProf.set(npts, step, -PI);

	sprintf(msg, "cluster %d==>%d :", pc2->cID, pc1->cID); show_log(msg, true);

	FLEX_ASYM_MATRIX<TORSION_PARS> tpar;
	TORSION_PARS *pTorPar = NULL;
	tpar.set(c1.atom[0].nb - 1, c2.atom[0].nb - 1);
	patom1 = c1.atom; patom2 = c2.atom;
	for (n1 = 0; n1 < patom1->nb - 1; n1++) {
		if (patom1->ba[n1] == NULL) {
			sprintf(errmsg, "torsion profile calculation wrong : atoms duplication is not consistent!");
			show_msg(errmsg); return;
		}
		for (n2 = 0; n2 < patom2->nb - 1; n2++) {
			if (patom2->ba[n2] == NULL) {
				sprintf(errmsg, "torsion profile calculation wrong : atoms duplication is not consistent!");
				show_msg(errmsg); return;
			}

			if (general_find_torsion_pars(force_field_fname, false, 
				((MATOM*)(patom1->ba[n1]))->par->atom, patom1->par->atom, 
				patom2->par->atom, ((MATOM*)(patom2->ba[n2]))->par->atom, 
				tpar.m[n1].m[n2], ::errmsg)) {
					pTorPar = &(tpar.m[n1].m[n2]);
					while (pTorPar != NULL) {
						sprintf(msg, "%s-%s-%s-%s : V = %f kCal/mol, phase = %f, pn = %f", ((MATOM*)(patom1->ba[n1]))->par->atom,
							patom1->par->atom, patom2->par->atom, ((MATOM*)(patom2->ba[n2]))->par->atom, pTorPar->v0, pTorPar->phase, pTorPar->pn);
						show_log(msg, true);

						pTorPar->v0 *= 4.176f; //kCal/mol ==> kJ/mol
						pTorPar->phase *= fAngle2Radian; // to arc
						pTorPar->pn = FABS(pTorPar->pn); // |pn|  

						pTorPar = pTorPar->next;
					}
			}
			else {
				//mlog.show(errmsg);
				tpar.m[n1].m[n2].v0 = 0; tpar.m[n1].m[n2].phase = 0; tpar.m[n1].m[n2].pn = 0;
			}
		}
	}

	VECTOR3 r0, r1;
	for (i = 0; i < pTor->tProf.npts; i++) {
		r0 = *(c2.hinge[c2.n2Parent].O); r1 = r0 + c2.hinge[c2.n2Parent].vHinge;
		c2.rotate(r0, r1, pTor->tProf.x[i] - c2.km.ta, true); // rotate to current x[i] torsion angle
		c2.km.ta = pTor->tProf.x[i];
		calTorsionEnergyTorque1(&c1, &c2, tpar, pTor->tProf.y[0][i], pTor->tProf.y[1][i]); // pTor->torqueProf.y[i]);
	}
	tpar.release();
	c1.reset(); c2.reset();
}

void MMOLECULE::SetupHingeTorsion() {
	int nc = 0;
	BASIC_CLUSTER *parent = NULL, *pc = NULL;
	TORSION *tor = NULL;

	char ofname[256] = "\0", buffer[256] = "\0";
	ofstream *out = NULL;
	int i;

	show_infor("calculating torsion profiles ... ");
	for (nc = 0; nc < this->nCluster; nc++) {
		pc = (BASIC_CLUSTER*)(&(cluster[nc]));
		parent = (BASIC_CLUSTER*)pc->parent;
		if (parent == NULL) continue;
		tor = torsion_db.get_torsion(parent->cID, pc->cID);
		if (tor == NULL) {
			tor = new TORSION(parent->cID, pc->cID);
			pc->pTors = tor;
			torsion_db.attach(tor);
			calTorsionProfile(parent, pc);

			sprintf(ofname, "tor-%d-%d.dat", pc->cID, parent->cID);
			out = new ofstream;
			out->open(ofname);
			if (out->is_open()) {
				for (i = 0; i < tor->tProf.npts; i++) {
					//(*out)<<tor->tProf.x[i]<<" "<<tor->tProf.y[i]<<" "<<tor->torqueProf.y[i]<<endl;
					(*out)<<tor->tProf.x[i]<<" "<<tor->tProf.y[0][i]<<" "<<tor->tProf.y[1][i]<<endl;
				}
				out->close();
			}
			delete out; out = NULL;
			sprintf(buffer, "torsion energy %d ==> %d is saved in %s", pc->cID, parent->cID, ofname);
			show_log(buffer, true);
			show_log("\n", true);
		}
		else pc->pTors = tor;
	}
	nc = number(torsion_db.ch_torsion);
	sprintf(buffer, "%d torsion profiles are calculated", nc);
	show_infor(buffer);
}

void MMOLECULE::calTorsion(double &E) {
	int nc = 0;
	double Etor = 0, torque = 0;
	CLUSTER *pc = NULL, *pc1 = NULL;
	E = 0;
	VECTOR3 t;
	for (nc = 0; nc < nCluster; nc++) {
		pc = cluster + nc;
		if (pc->pTors == NULL) continue;
		pc->pTors->get_torsion(pc->km.ta, Etor, torque);
		E += Etor; // in kT
		t.v[0] = -torque * pc->up.v[0]; // in atom inertia moment unit
		t.v[1] = -torque * pc->up.v[1];
		t.v[2] = -torque * pc->up.v[2];

		pc->dyn->fc.v[0] += t.v[0];
		pc->dyn->fc.v[1] += t.v[1];
		pc->dyn->fc.v[2] += t.v[2];

		// because torsion induced net force is 0, the torque from torsion is transferable
		// from one point to any other point
		// to parent, the torsion induced torque is negative to that of child
		if (pc->parent != NULL) {
			pc1 = (CLUSTER*)(pc->parent);
			pc1->dyn->fc.v[0] -= t.v[0];
			pc1->dyn->fc.v[1] -= t.v[1];
			pc1->dyn->fc.v[2] -= t.v[2];
		}
	}
}

void calTorsion(MMOLECULE*mm, int nc1, int nc2, double &E) { // same as MMOLECULE::calTorsion(double &E) above
	E = 0;
	int nc = 0;
	double Etor = 0, torque = 0;
	CLUSTER *pc = NULL, *pc1 = NULL;
	E = 0;
	VECTOR3 t;
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		if (pc->pTors == NULL) continue;
		pc->pTors->get_torsion(pc->km.ta, Etor, torque);
		E += Etor; // in kT
		t.v[0] = -torque * pc->up.v[0]; // in atom inertia moment unit
		t.v[1] = -torque * pc->up.v[1];
		t.v[2] = -torque * pc->up.v[2];

		pc->dyn->fc.v[0] += t.v[0];
		pc->dyn->fc.v[1] += t.v[1];
		pc->dyn->fc.v[2] += t.v[2];

		// because torsion induced net force is 0, the torque from torsion is transferable
		// from one point to any other point
		// to parent, the torsion induced torque is negative to that of child
		if (pc->parent != NULL) {
			pc1 = (CLUSTER*)(pc->parent);
			pc1->dyn->fc.v[0] -= t.v[0];
			pc1->dyn->fc.v[1] -= t.v[1];
			pc1->dyn->fc.v[2] -= t.v[2];
		}
	}
	return;
}
#endif  // _DSIABLE_ TORSION

void MMOLECULE::calTorsionAngle() { //calculate torsion angle with parent-child relationship setup
	int nc = 0;
	BASIC_CLUSTER *parent = NULL;
	CLUSTER *pc = NULL;
	for (nc = 0; nc < nCluster; nc++) {
		pc = cluster + nc;
		parent = (BASIC_CLUSTER*)pc->parent;
		if (parent != NULL) pc->km.ta = torsion_angle(parent, (BASIC_CLUSTER*)pc);
	}
}

// **********************   calculating torsion, above  ***********************

// rotate all the clusters related with pc, except the one linked with parent_hinge
// with given rotation matrix, original point -- r0, and after the rotation, shift dr
void series_shift_rotate(CLUSTER *pc, int parent_hinge, VECTOR3 &dr, VECTOR3 &r0, R &rmatrix) {
	int n;
	pc->shift_rotate(dr, r0, rmatrix);
	BASIC_CLUSTER *child = NULL;
	for (n = 0; n < pc->nHinges; n++) {
		if (n == parent_hinge) continue;
		child = (BASIC_CLUSTER*)(pc->hinge[n].plink);
		if (child != NULL) series_shift_rotate((CLUSTER*)child, pc->hinge[n].indx_link, dr, r0, rmatrix);
	}
}

// rotate pc + all clusters associated with pc except the one linked with parent_hinge, along parent_hinge with domega
void series_rotate(CLUSTER *pc, int parent_hinge, double domega, bool arc) { 
	if (parent_hinge >= pc->nHinges) return;
	R rmatrix;
	RMatrix(pc->hinge[parent_hinge].vHinge, domega, rmatrix, arc);
	VECTOR3 dr(0, 0, 0); // no shift
	VECTOR3 r0 = *(pc->hinge[parent_hinge].O);
	series_shift_rotate(pc, parent_hinge, dr, r0, rmatrix);
}


void MMOLECULE::setup_base_movement(VECTOR3 &dv, VECTOR3 &v0, VECTOR3 &drot, bool arc) {
	int i;
	V2V(dv, base->rot.dr, i, 3)
	double domega = drot.abs();
	VECTOR3 axis, v1;
	if (domega > MIN_ROT) {
		for (i = 0; i < 3; i++) axis.v[i] = drot.v[i] / domega;
		RMatrix(axis, domega, base->rot.R0, true);
	}
	else {
		I3Matrix(base->rot.R0)
	}
	MpV3(base->rot.R0, v0, v1)
	V3minusV3(v0, v1, base->rot.dr) // base->dr = v0 - base->R0 * v0
	V3plusV3(base->rot.dr, dv, base->rot.dr) // base->dr += dv
}

void setup_LocalRotMatrix(MMOLECULE *mm, int nc1, int nc2) {
	int nc = 0;
	CLUSTER *pc = NULL;
	BASIC_CLUSTER *parent = NULL;
	VECTOR3 v0, v1;
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		parent = (BASIC_CLUSTER*)pc->parent;
		// we do not calculate rotation matrix of base cluster
		if (parent == NULL) continue; 
		RMatrix(pc->up, pc->rot.dta, pc->rot.R0, true);
		//RMatrix(pc->hinge[pc->n2Parent].vHinge, pc->rot.dta, pc->rot.R0, true);
		V32V3((*(pc->Op)), v0)
		MpV3(pc->rot.R0, v0, v1)
		V3minusV3(v0, v1, pc->rot.dr) // pc->dr = -pc->R0 * v0 + v0
	}
}

void setup_all_LocalRotMatrix(MMOLECULE &mm) {
	MacroMoleculeOperate<MMOLECULE>((void*)(&setup_LocalRotMatrix), mm, 0, mm.nCluster, MAX_THREADS);
}

void setup_GlobalRotMatrix(MMOLECULE &mm) {
	int i, j, k;
	VECTOR3 v1;
	LIST< LIST<CLUSTER> > *cl = mm.ctree.tl;
	LIST<CLUSTER> *clist = NULL;
	CLUSTER *pc = NULL, *parent = NULL;
	while (cl != NULL) {
		clist = cl->p;
		while (clist != NULL) {
			pc = clist->p;
			parent = (CLUSTER*)(pc->parent);
			
			if (parent == NULL) {
				M2M(pc->rot.R0, pc->rot.Rg, i, j, 3) // Rg = R0
				V32V3(pc->rot.dr, pc->rot.drg)  // drg = dr
			}
			else {
				MpM(parent->rot.Rg, pc->rot.R0, pc->rot.Rg, i, j, k, 3) // Rg = parent->Rg * pc->R0
				MpV3(parent->rot.Rg, pc->rot.dr, v1)  // v1 = parent->Rg * pc->dr
				V3plusV3(v1, parent->rot.drg, pc->rot.drg)  // pc->drg = parent->drg + parent->Rg * pc->dr
			}
			clist = clist->next;
		}
		cl = cl->next;
	}
}

void tweak_clusters(MMOLECULE *mm, int nc1, int nc2) {
	int nc = 0, na = 0;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 v0, v1;
	DMATRIX<3> m0, m1;
	VECTOR3 *pr = NULL;
	double rabs = 0;
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			V32V3(patom->r, v0)
			MpV3(pc->rot.Rg, v0, v1)
			V3plusV3(pc->rot.drg, v1, patom->r) // r = pc->Rg * r + pc->drg
		}
		/*
		// calculate hinges, and up vector
		for (na = 0; na < pc->nHinges; na++) {
			V32V3(pc->hinge[na].vHinge, v0) // hinge vector
			MpV3(pc->rot.Rg, v0, pc->hinge[na].vHinge) // hinge[na].vHinge = pc->Rg * hinge[i].vHinge
		}
		if (pc->parent != NULL) {
			V32V3(pc->up, v0) // hinge vector
			MpV3(pc->rot.Rg, v0, pc->up) // hinge[na].vHinge = pc->Rg * hinge[i].vHinge
		}
		*/
		// mass center, global cordinate
		V32V3(pc->cm, v0)
		MpV3(pc->rot.Rg, v0, v1)
		V3plusV3(pc->rot.drg, v1, pc->cm) // r = pc->Rg * r + pc->drg

		if (bImplicitSolvent) {
			// geometrical center, global cordinate
			V32V3(pc->rc_solv, v0)
			MpV3(pc->rot.Rg, v0, v1)
			V3plusV3(pc->rot.drg, v1, pc->rc_solv) // r = pc->Rg * r + pc->drg
		}
	}
}

void check_hinges(MMOLECULE *mm, int nc1, int nc2) {
	int nc = 0, na = 0;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 *pr = NULL;
	VECTOR3 *Op, *Op_parent;
	double rabs = 0;
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		for (na = 0; na < pc->nHinges; na++) {
			VECT3((*(pc->hinge[na].O)), (*(pc->hinge[na].Olink)), pc->hinge[na].vHinge)
		}
		if (pc->parent != NULL) {
			pr = &(pc->hinge[pc->n2Parent].vHinge);
			scale_uv3((*pr), (*pr), rabs);
			if (rabs > 0.01) {
				rabs = sqrt(rabs);
				DIV3((*pr), rabs, pc->up)
			}
			Op = pc->Op; Op_parent = pc->parent->Op; pr = &(pc->l);
			VECT3((*Op_parent), (*Op), (*pr)) // l = Op_parent ==> Op
		}
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			V3minusV3(patom->r, (*(pc->Op)), patom->r0)
		}
		// mas center, local cordinate
		V3minusV3(pc->cm, (*(pc->Op)), pc->cm0)
		if (bImplicitSolvent) {
			// geometrical center, local cordinate
			V3minusV3(pc->rc_solv, (*(pc->Op)), pc->rc0_solv)
		}
	}
}

void tweak_all_clusters(MMOLECULE &mm) {
	if (mm.nCluster > NCLUSTERS_PARALLEL) {
		MacroMoleculeOperate<MMOLECULE>((void*)(&tweak_clusters), mm, 0, mm.nCluster, MAX_THREADS);
		MacroMoleculeOperate<MMOLECULE>((void*)(&check_hinges), mm, 0, mm.nCluster, MAX_THREADS); // if tweak_clusters does not configure the hinges, do it here
	}
	else {
		tweak_clusters(&mm, 0, mm.nCluster);
		check_hinges(&mm, 0, mm.nCluster);
	}
	VECTOR3 v0;
	V32V3(mm.xaxis, v0)
	MpV3(mm.base->rot.Rg, v0, mm.xaxis) // rotate xaxis
	V32V3(mm.zaxis, v0)
	MpV3(mm.base->rot.Rg, v0, mm.zaxis) // rotate xaxis
}


// return energy in kT
void calClustersKineticEnergy(MMOLECULE *mm, int n1, int n2, double& Ek) {
	Ek = 0;
	double unit_E = ::unit_Ek_kT * 0.5; // important, here unit_E is divided by 2
	double Ekc = 0;
	VECTOR<6> Iw;
	int i = 0, j = 0, n = 0;
	CLUSTER *pc = NULL;
	for (n = n1; n <= n2; n++) {
		if (n >= mm->nCluster) break;
		pc = mm->cluster + n;
		MpV6(pc->invNE->M, pc->bkm->V0, Iw) // Iw = M * V0;
		scale_uv6(pc->bkm->V0, Iw, Ekc) // Ekc = V0 * M * V0;
		Ekc *= unit_E;
		pc->Ek = (float)Ekc;
		Ek += Ekc;
	}
	//return Ek;
}

// return energy in kT
double ParallelCalKineticEnergy(MMOLECULE &mm) {
	double Ek = 0;

	int nThreads = MAX_THREADS;
	MacroMoleculeOperate2<MMOLECULE, double>((void*)(&calClustersKineticEnergy), mm, 0, mm.nCluster, nThreads, Ek);
	return Ek;
}


// Spatial Moment is relative to the base cluster's cordinate
void calClusterSpatialMomentum(MMOLECULE *mm, int n1, int n2, VECTOR<6> &res) {
	int i = 0, j = 0, nc = 0;
	CLUSTER *pc = NULL;
	VECTOR<6> iw, iw0, Iw;

	//VECTOR<3> dr;
	//MATRIX<6> phai;

	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		MpV6(pc->invNE->M, pc->bkm->V0, iw)
		// dr is from base to cluster
		/*
		for (i = 0; i < 3; i++) {
			dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm->base->Op->v[i];
		}
		PHAI(phai, dr, i, j)
		*/
		MpV6(pc->invNE->phai_base2me, iw, iw0)
		//for (i = 0; i < 6; i++) Iw.v[i] += iw0.v[i];
		V6plusV6(Iw, iw0, Iw)
	}
	V62V6(Iw, res)
}

void ParallelCalSpatialMomentum(MMOLECULE &mm) {
	V6zero(mm.Iw)
	int nThreads = MAX_THREADS;
	MacroMoleculeOperate2<MMOLECULE, VECTOR<6> >((void*)(&calClusterSpatialMomentum), mm, 0, mm.nCluster, nThreads, mm.Iw);
}

// Spatial Acceleration Moment is relative to the base cluster's cordinate
void calClusterSpatialAccel(MMOLECULE *mm, int n1, int n2, VECTOR<6> &res) {
	int i = 0, j = 0, nc = 0;
	CLUSTER *pc = NULL;
	VECTOR<6> iw, iw0, Ialpha;

	//VECTOR<3> dr;
	//MATRIX<6> phai;

	V6zero(Ialpha);
	//for (nc = 0; nc < mm.nCluster; nc++) {
	for (nc = n1; nc <= n2; nc++) {
		if (nc >=  mm->nCluster) break;
		pc = mm->cluster + nc;
		MpV6(pc->invNE->M, pc->bkm->alpha, iw)
		// dr is from base to cluster
		/*
		for (i = 0; i < 3; i++) {
			dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm->base->Op->v[i];
		}
		PHAI(phai, dr, i, j)
		*/
		MpV6(pc->invNE->phai_base2me, iw, iw0)
		//for (i = 0; i < 6; i++) Ialpha.v[i] += iw0.v[i];
		V6plusV6(Ialpha, iw0, Ialpha)
	}
	V62V6(Ialpha, res)
}

void ParallelCalSpatialAccel(MMOLECULE &mm, VECTOR<6> &Iw) {
	int nThreads = MAX_THREADS;
	MacroMoleculeOperate2<MMOLECULE, VECTOR<6> >((void*)(&calClusterSpatialAccel), mm, 0, mm.nCluster, nThreads, Iw);
}


CLUSTER_IMPSOLV_PARS* CLUSTER_IMPSOLV_PARS_DATABASE::attach(short cID) {
	CHAIN<CLUSTER_IMPSOLV_PARS> *pch = ch;
	CLUSTER_IMPSOLV_PARS *par = NULL;
	while (pch != NULL) {
		par = pch->p;
		if (par->cID == cID) return par; // this par is in the chain
		else if (pch->next == NULL) break;
		else pch = pch->next;
	};
	par = new CLUSTER_IMPSOLV_PARS; par->cID = cID;
	//if (!get_solvation_energy(par)) {delete par; return NULL;}
	if (pch == NULL) {
		ch = new CHAIN<CLUSTER_IMPSOLV_PARS>; ch->p = par;
	}
	else {
		pch->next = new CHAIN<CLUSTER_IMPSOLV_PARS>; pch->next->p = par; par->id = pch->p->id + 1;
	}
	return par;
}
/*
void CLUSTER::setup_solvation_neighbors(int N) {
	TREE<CLUSTER> tree;
	LIST<CLUSTER> *clist = NULL, *parent_list = NULL, *child_list = NULL;
	LIST< LIST<CLUSTER> > *plevel = NULL, *parent_level = NULL, *child_level = NULL;
	CLUSTER *pc = NULL, *pc1 = NULL;

	int n = 0, i;
	tree.new_level();
	tree.add(tree.tl, this, false); // this cluster is the root
	plevel = tree.tl;

	int n_neighbors = 0;
	for (n = 1; n <= N; n++) {
		child_level = tree.new_level();
		if (plevel->pre != NULL) parent_list = plevel->pre->p;
		else parent_list = NULL;
		clist = plevel->p;
		while (clist != NULL) {
			pc = clist->p;
			if (pc != NULL) {
				for (i = 0; i < pc->nHinges; i++) {
					pc1 = (CLUSTER*)(pc->hinge[i].plink);
					if (parent_list != NULL && parent_list->get_list(pc1) != NULL) continue; //its parent is in the parent level
					else {
						if (child_level->p == NULL) {
							child_level->p = new LIST<CLUSTER>;
							child_level->p->p = pc1;
						}
						else child_level->p->add(pc1, false);
						n_neighbors++;
					}
				}
			}
			clist = clist->next;
		}
		if (child_level->p == NULL) {
			parent_level = child_level->pre;
			parent_level->next = NULL;
			delete child_level; child_level = NULL;
			break;
		}
		plevel = child_level;
	}
	for (i = 0; i < neighbor_Solve.n; i++) neighbor_Solve.m[i] = NULL;
	neighbor_Solve.release();
	if (n_neighbors == 0) return;
	neighbor_Solve.SetArray(n_neighbors);
	plevel = tree.tl->next;
	n = 0;
	while (plevel != NULL) {
		clist = plevel->p;
		while (clist != NULL) {
			if (n < n_neighbors) neighbor_Solve.m[n] = clist->p;
			n++;
			clist = clist->next;
		}
		plevel = plevel->next;
	}
	tree.release_tree(false);
}

void MMOLECULE::setup_Solvation_neighbors(int N) {
	for (int n = 0; n < nCluster; n++) cluster[n].setup_solvation_neighbors(N);
}

class CLUSTER_DIS {
public:
	CLUSTER *pbase;
	float dis;
	CLUSTER_DIS() {pbase = NULL; dis = 0;};
	CLUSTER_DIS(CLUSTER *pc, float dis) {this->pbase = pc; this->dis = dis;};
	~CLUSTER_DIS() {pbase = NULL; dis = 0;};
};

// setup solvation-neighbors with cut-distance
void setup_solvation_neighbors(CLUSTER* rpc, float dis_cut) {
	TREE<CLUSTER_DIS> tree;
	LIST<CLUSTER_DIS> *clist = NULL, *parent_list = NULL, *child_list = NULL;
	LIST< LIST<CLUSTER_DIS> > *plevel = NULL, *parent_level = NULL, *child_level = NULL;
	
	CLUSTER_DIS *pdis = NULL;
	CLUSTER *pc = NULL, *pc1 = NULL;

	int n = 0, i;
	tree.new_level();
	//pdis = new CLUSTER_DIS(this, 0);
	pdis = new CLUSTER_DIS(rpc, 0);
	tree.add(tree.tl, pdis, false); // this cluster is the root
	plevel = tree.tl;

	int n_neighbors = 0;
	while (true) {
		child_level = tree.new_level();
		if (plevel->pre != NULL) parent_list = plevel->pre->p;
		else parent_list = NULL;
		clist = plevel->p;
		while (clist != NULL) {
			//if (clist->p->dis > dis_cut) {clist = clist->next; continue;} // this branch is too far
			pc = clist->p->pbase;
			if (pc != NULL) {
				for (i = 0; i < pc->nHinges; i++) {
					pc1 = (CLUSTER*)(pc->hinge[i].plink);
					if (clist->p->dis > dis_cut + pc->r_solve + pc1->r_solve) continue; // this branch is too far
					if (parent_list != NULL && get_wrap_list<CLUSTER_DIS, CLUSTER>(parent_list, pc1) != NULL) continue; //its parent is in the parent level
					else {
						pdis = new CLUSTER_DIS(pc1, clist->p->dis + pc1->r_solve + pc->r_solve);
						if (child_level->p == NULL) {
							child_level->p = new LIST<CLUSTER_DIS>;
							child_level->p->p = pdis;
						}
						else {
							child_level->p->add(pdis, false);
						}
						n_neighbors++;
					}
				}
			}
			clist = clist->next;
		}
		if (child_level->p == NULL) {
			parent_level = child_level->pre;
			parent_level->next = NULL;
			delete child_level; child_level = NULL;
			break;
		}
		plevel = child_level;
	}

	for (i = 0; i < rpc->neighbor_Solve.n; i++) rpc->neighbor_Solve.m[i] = NULL;
	rpc->neighbor_Solve.release();
	if (n_neighbors == 0) return;
	rpc->neighbor_Solve.SetArray(n_neighbors);
	plevel = tree.tl->next;
	n = 0;
	while (plevel != NULL) {
		clist = plevel->p;
		while (clist != NULL) {
			if (n < n_neighbors) rpc->neighbor_Solve.m[n] = clist->p->pbase;
			n++;
			clist = clist->next;
		}
		plevel = plevel->next;
	}
	tree.release_tree(true);
}

void MMOLECULE::setup_Solvation_neighbors(float dis_cut) {
	for (int n = 0; n < nCluster; n++) setup_solvation_neighbors(cluster + n, dis_cut);
}
*/

void MMOLECULE::calTotalForce() {
	int n, i;
	CLUSTER *pc = NULL;
	VECTOR<6> *f, *fc, *f_Brown, *f_hoover, *fc_cm = &(mDynKin.fc), fct;
	MATRIX<6> *phai;
	V6zero((*fc_cm))

	bool bForceMassCenter = (::MM_SpeedVerlet_method == 0 ? true : false); // if we check momentum conservation, force at mass center has to be calculated
	for (n = 0; n < nCluster; n++) {
		pc = cluster + n;
		fc = &(pc->dyn->fc); f_hoover = &(pc->dyn->f_Hoover); 
		f = &(pc->dyn->f); f_Brown = &(pc->dyn->f_Brown);
		for (i = 0; i < 6; i++) f->v[i] = fc->v[i] + f_Brown->v[i] + f_hoover->v[i];

		if (bForceMassCenter) {
			phai = &(pc->invNE->phai_cm2me);
			MpV6((*phai), (*fc), fct)
			V6plus(fc_cm->v, fct.v, fc_cm->v)
		}
	}
	if (bForceMassCenter) {
		V6plusV6(mDynKin.fc, mDynKin.f_Hoover, mDynKin.f)
		V6plusV6(mDynKin.f, mDynKin.f_Brown, mDynKin.f) // Brownian force
	}
}

float BASIC_CLUSTER::cal_geometry_radius() {
	int na;
	double d, dmax = 0;
	VECTOR3 r, dr;
	MATOM *patom;
	if (matom.m == NULL) {
		for (na = 0; na < nAtoms; na++) {
			patom = atom + na;
			r.v[0] += patom->r.v[0];
			r.v[1] += patom->r.v[1];
			r.v[2] += patom->r.v[2];
		}
		r.v[0] /= nAtoms; r.v[1] /= nAtoms; r.v[2] /= nAtoms;
	}
	else {
		for (na = 0; na < matom.n; na++) {
			patom = matom.m[na];
			r.v[0] += patom->r.v[0];
			r.v[1] += patom->r.v[1];
			r.v[2] += patom->r.v[2];
		}
		r.v[0] /= matom.n; r.v[1] /= matom.n; r.v[2] /= matom.n;
	}

	if (bImplicitSolvent) {
		V32V3(r, this->rc_solv)
		V3minusV3(rc_solv, (*Op), rc0_solv)
	}

	dmax = 0;
	if (matom.m == NULL) {
		for (na = 0; na < nAtoms; na++) {
			patom = atom + na;
			VECT3(r, patom->r, dr) V3ABS2(dr, d) d = sqrt(d);
			if (dmax < d) dmax = d;
		}
	}
	else {
		for (na = 0; na < matom.n; na++) {
			patom = matom.m[na];
			VECT3(r, patom->r, dr) V3ABS2(dr, d) d = sqrt(d);
			if (dmax < d) dmax = d;
		}
	}
	return (float)dmax;
}

void BASIC_CLUSTER::calSolvationRadius() {
	float d2 = 0;
	if (this->matom.n == 1) {
		r_solve = 0.6f;
		this->d = 1.2f;
	} // single atom in cluster, set size 1.2 Angs.
	else {
		r_solve = cal_geometry_radius();
		this->d = r_solve + 1;  // H has LJ distance, say 0.5
	}
	d2 = this->d * this->d;

	if (::bBrownian) {
		this->ita_Brownian = ::ita_Brownian * r_solve;
		this->cd_rot = ::cd_drag_rot * d2;
		this->cd_trans = ::cd_drag_trans * d2;
	}
}

void MMOLECULE::show_cluster_BrownianForcePars(CHAIN<short> **ch_id) {
	char msg[256] = "\0";
	bool show = true;
	CHAIN<short> *pch = NULL, *tail = NULL;
	for (int nc = 0; nc < nCluster; nc++) {
		show = true;
		if (ch_id != NULL && *ch_id != NULL) {
			if (exist_in_chain<short>(*ch_id, cluster[nc].cID)) show = false;
		}
		if (show) {
			sprintf(msg, "cluster [%d] : %f, %f, %f", cluster[nc].cID, cluster[nc].ita_Brownian, cluster[nc].cd_trans, cluster[nc].cd_rot);
			show_log(msg, true);
			if (ch_id != NULL) {
				pch = new CHAIN<short>; pch->p = new short; *(pch->p) = cluster[nc].cID;
				if (*ch_id == NULL) *ch_id = pch;
				else {
					if (tail == NULL) tail = (*ch_id)->get_tail(); 
					tail->next = pch; tail = pch;
				}
			}
		}
	}
}

void BASIC_CLUSTER::show_BrownianForcePars(CHAIN<short> **ch_id) {
	char msg[256] = "\0";
	bool show = true;
	CHAIN<short> *pch = NULL, *tail = NULL;

	show = true;
	if (ch_id != NULL && *ch_id != NULL) {
		if (exist_in_chain<short>(*ch_id, this->cID)) show = false;
	}
	if (show) {
		sprintf(msg, "cluster [%d] : %f, %f, %f", cID, ita_Brownian, cd_trans, cd_rot);
		show_log(msg, true);
		if (ch_id != NULL) {
			pch = new CHAIN<short>; pch->p = new short; *(pch->p) = cID;
			if (*ch_id == NULL) *ch_id = pch;
			else {
				if (tail == NULL) tail = (*ch_id)->get_tail(); 
				tail->next = pch; tail = pch;
			}
		}
	}
}

void CalClusterElectroMultipole(MMOLECULE *mm, int nc1, int nc2, bool bEwaldSum) {
	CLUSTER *pc = NULL;
	for (int nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		pc->calElectroMultipole(bEwaldSum);
	}
}

void MMOLECULE::Init_MolFrameRot2(VECTOR3& dr, VECTOR3& zaxis2, VECTOR3& xaxis2) {
	VECTOR3 raxis, xaxis;
	double domega = 0;
	R Rz, Rx;
	int i, j, k;

	cp(zaxis, zaxis2, raxis);
	domega = angle(zaxis, zaxis2);
	RMatrix(raxis, domega, Rz, true);

	MpV3(Rz, this->xaxis, xaxis); // xaxis of molecule changes after Rz operation
	cp(xaxis, xaxis2, raxis);
	domega = angle(xaxis, xaxis2);
	RMatrix(raxis, domega, Rx, true);

	MpM(Rx, Rz, this->base->rot.R0, i, j, k, 3)

	// important: the rotation above uses origin at base->Op, so, base->Op does not change after the rotating operation above
	VECTOR3 v0, v1;
	V32V3((*(this->base->Op)), v0)

	MpV3(base->rot.R0, v0, v1)
	V3minusV3(v0, v1, base->rot.dr) // base->dr = v0 - base->R0 * v0
	V3plusV3(base->rot.dr, dr, base->rot.dr) // base->dr += dr
}

void MMOLECULE::Init_MolFrameRot(VECTOR3& dr, VECTOR3& drot) {
	int i;

	//V2V(dv, base->rot.dr, i, 3)
	double domega = drot.abs();
	VECTOR3 axis;
	if (domega > MIN_ROT) {
		for (i = 0; i < 3; i++) axis.v[i] = drot.v[i] / domega;
		RMatrix(axis, domega, base->rot.R0, true);
	}
	else {
		I3Matrix(base->rot.R0)
	}

	VECTOR3 v0, v1;
	V32V3((*(this->base->Op)), v0)

	MpV3(base->rot.R0, v0, v1)
	V3minusV3(v0, v1, base->rot.dr) // base->dr = v0 - base->R0 * v0
	V3plusV3(base->rot.dr, dr, base->rot.dr) // base->dr += dr
}

void MacroMolMove(MMOLECULE *mm) { // using the parameter in MOL_BASE_MOVEMENT (mbm) & CLUSTER_ROT (rot)
	//suppose mdpar.nCluster == mm.nCluster
	int nc = 0;
	CLUSTER *pc = NULL;
	VECTOR3 raxis, ZeroV;
	double domega = 0;

	if (mm->nCluster < MAX_THREADS * 100) {
	/*******************************************************
				rotation all clusters in series
			NOTE: lines started with // is commented
	*******************************************************/
	
		if (mm->mbm.bRot) {
			//mm.twMol(*(mm.base->Op), drot, true);
			//mm.shiftMol(dr); // translation was done when shift molecule
			mm->shift_tw_mol(mm->mbm.dr, *(mm->base->Op), mm->mbm.drot, true);
		}
		else {
			mm->shiftMol(mm->mbm.dr);
			cp(mm->zaxis, mm->mbm.zaxis, raxis);
			domega = angle(mm->zaxis, mm->mbm.zaxis);
			mm->twMol((*(mm->base->Op)), raxis, domega, true);

			cp(mm->xaxis, mm->mbm.xaxis, raxis);
			domega = angle(mm->xaxis, mm->mbm.xaxis);
			mm->twMol((*(mm->base->Op)), raxis, domega, true);
		}

		mm->setup_base_movement(ZeroV, *(mm->base->Op), ZeroV, true);
	
		// move other clusters and guess velocity for next step
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = mm->cluster + nc;
			if (mm->free_base && pc == mm->base) continue;
			pc->km.ta += pc->rot.dta;
			mm->twTA(nc, pc->rot.dta, true); // in arc degree
		
			//delta = torsion_angle(p->parent, p); // check torsion angle
		}
	}
	else {
	/*******************************************************
				rotation all clusters parallelly
			NOTE: lines started with // is commented
	*******************************************************/
		// move all clusters 

		//mm.setup_base_movement(ZeroV, *(mm.base->Op), ZeroV, true);
		if (mm->mbm.bRot) mm->Init_MolFrameRot(mm->mbm.dr, mm->mbm.drot);
		else mm->Init_MolFrameRot2(mm->mbm.dr, mm->mbm.zaxis, mm->mbm.xaxis);

		setup_all_LocalRotMatrix(*mm);
		setup_GlobalRotMatrix(*mm);
		tweak_all_clusters(*mm);

		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = mm->cluster + nc;
			if (pc == mm->base) continue;
			pc->km.ta += pc->rot.dta;
		
			//delta = torsion_angle(p->parent, p); // check torsion angle
		}
	}
}

// calculate effective dipole moment of multi-macromolecules
void polarize_MM(MMOLECULE *mm, int n) {
	int imol, nc, iatom;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	double *mu = NULL, *E = NULL, *mut = NULL, *mu0 = NULL;
	ATOM_COMM_PARS *par = NULL;
	double fp = ::fp;
	double d = 0;
	for (imol = 0; imol < n; imol++) {
		for (nc = 0; nc < mm[imol].nCluster; nc++) {
			pc = mm[imol].cluster + nc;
			if (pc->eNeutral) continue;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom; 
				if (patom->mMP > 0) {
					mu = patom->ind_mu.v; E = patom->E.v; par = patom->par;
					d = par->alpha * E[0]; mu[0] = mu[0] + fp * (d - mu[0]);
					d = par->alpha * E[1]; mu[1] = mu[1] + fp * (d - mu[1]);
					d = par->alpha * E[2]; mu[2] = mu[2] + fp * (d - mu[2]);

					mu0 = patom->intr_mu.v; mut = patom->mu.v;
					mut[0] = mu0[0] + mu[0]; mut[1] = mu0[1] + mu[1]; mut[2] = mu0[2] + mu[2];
				}
			}
		}
	}
}

void polarize_VMM(MMOLECULE **mm, int n) {
	int imol, nc, iatom;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	double *mu = NULL, *E = NULL, *mut = NULL, *mu0 = NULL;
	ATOM_COMM_PARS *par = NULL;
	double fp = ::fp;
	double d = 0;
	for (imol = 0; imol < n; imol++) {
		for (nc = 0; nc < mm[imol]->nCluster; nc++) {
			pc = mm[imol]->cluster + nc;
			if (pc->eNeutral) continue;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom; 
				if (patom->mMP > 0) {
					mu = patom->ind_mu.v; E = patom->E.v; par = patom->par;
					d = par->alpha * E[0]; mu[0] = mu[0] + fp * (d - mu[0]);
					d = par->alpha * E[1]; mu[1] = mu[1] + fp * (d - mu[1]);
					d = par->alpha * E[2]; mu[2] = mu[2] + fp * (d - mu[2]);

					mu0 = patom->intr_mu.v; mut = patom->mu.v;
					mut[0] = mu0[0] + mu[0]; mut[1] = mu0[1] + mu[1]; mut[2] = mu0[2] + mu[2];
				}
			}
		}
	}
}

void CLUSTER::cal_cluster_inertia_tensor_masscenter() {
	int ia;
	MATOM *patom = NULL;
	MATRIX<3> *I = &(mcDyn->I);
	VECTOR3 dr;
	float m = 0, a = 0;
	SVECTOR<double, 3> r0;
	
	M3zero((*I))
	if (matom.m == NULL) {
		for (ia = 0; ia < nAtoms; ia++) {
			patom = atom + ia; m = patom->par->m;
			VECT3((*Op), patom->r, r0)
			dr.v[0] = r0.v[0] - cm0.v[0];
			dr.v[1] = r0.v[1] - cm0.v[1];
			dr.v[2] = r0.v[2] - cm0.v[2];
			I->m[0][0] += m * (dr.v[1] * dr.v[1] + dr.v[2] * dr.v[2]);
			I->m[1][1] += m * (dr.v[0] * dr.v[0] + dr.v[2] * dr.v[2]);
			I->m[2][2] += m * (dr.v[0] * dr.v[0] + dr.v[1] * dr.v[1]);
			a = m * dr.v[0] * dr.v[1];
			I->m[0][1] -= a; I->m[1][0] -= a;
			a = m * dr.v[0] * dr.v[2];
			I->m[0][2] -= a; I->m[2][0] -= a;
			a = m * dr.v[1] * dr.v[2];
			I->m[1][2] -= a; I->m[2][1] -= a;
		}
	}
	else {
		for (ia = 0; ia < matom.n; ia++) {
			patom = matom.m[ia]; m = patom->par->m;
			VECT3((*Op), patom->r, r0)
			dr.v[0] = r0.v[0] - cm0.v[0];
			dr.v[1] = r0.v[1] - cm0.v[1];
			dr.v[2] = r0.v[2] - cm0.v[2];
			I->m[0][0] += m * (dr.v[1] * dr.v[1] + dr.v[2] * dr.v[2]);
			I->m[1][1] += m * (dr.v[0] * dr.v[0] + dr.v[2] * dr.v[2]);
			I->m[2][2] += m * (dr.v[0] * dr.v[0] + dr.v[1] * dr.v[1]);
			a = m * dr.v[0] * dr.v[1];
			I->m[0][1] -= a; I->m[1][0] -= a;
			a = m * dr.v[0] * dr.v[2];
			I->m[0][2] -= a; I->m[2][0] -= a;
			a = m * dr.v[1] * dr.v[2];
			I->m[1][2] -= a; I->m[2][1] -= a;
		}
	}
}


// Dihedral torsion relating to the hinge, calculating with force on each atom of the dihedral angle
void initTorsionPars(MMOLECULE* mm) { // initilize the database for torsion pars relating to the macromolecule
	int nc = 0, na, na1;
	CLUSTER *pc = NULL, *pc1 = NULL;
	MATOM *patom1, *patom2, *patom3, *patom4;
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge, *hinge1;
	TORSION_PARS *tpar = NULL, *pTorPar = NULL;

	char msg[256] = "\0";

	for (nc = 0; nc < mm->nCluster; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		if (pc == mm->base) continue;

		hinge = pc->hinge + pc->n2Parent;
		pc1 = (CLUSTER*)(hinge->plink);
		hinge1 = pc1->hinge + hinge->indx_link;
		patom2 = pc->atom + hinge->indx;
		patom3 = pc1->atom + hinge1->indx;

		for (na = 0; na < patom2->nb; na++) {
			patom1 = (MATOM*)(patom2->ba[na]);
			if (patom1 == patom3) continue;
			for (na1 = 0; na1 < patom3->nb; na1++) {
				patom4 = (MATOM*)(patom3->ba[na1]);
				if (patom4 == patom2) continue;
				tpar = dh_db.get_torsion_ch(patom1->aindx, patom2->aindx, patom3->aindx, patom4->aindx);
				if (tpar == NULL) {
					tpar = new TORSION_PARS;
					if (general_find_torsion_pars(force_field_fname, false, 
						patom1->par->atom, patom2->par->atom, patom3->par->atom, patom4->par->atom, *tpar, msg)) {
							pTorPar = tpar;
							while (pTorPar != NULL) {
								sprintf(msg, "%s-%s-%s-%s : V = %f kCal/mol, phase = %f, pn = %f", patom1->par->atom,
									patom2->par->atom, patom3->par->atom, patom4->par->atom, pTorPar->v0, pTorPar->phase, pTorPar->pn);
								show_log(msg, true);

								pTorPar->v0 *= 4.176f; //kCal/mol ==> kJ/mol
								pTorPar->phase *= fAngle2Radian; // to arc
								pTorPar->pn = FABS(pTorPar->pn); // |pn|  
								
								pTorPar->aindx[0] = patom1->aindx; pTorPar->aindx[1] = patom2->aindx;
								pTorPar->aindx[2] = patom3->aindx; pTorPar->aindx[3] = patom4->aindx;

								pTorPar = pTorPar->next;
							}
					}
					else {
						//mlog.show(errmsg);
						tpar->v0 = 0; tpar->phase = 0; tpar->pn = 0;
						tpar->aindx[0] = patom1->aindx; tpar->aindx[1] = patom2->aindx;
						tpar->aindx[2] = patom3->aindx; tpar->aindx[3] = patom4->aindx;
					}
					dh_db.attach(tpar);
				}
			}
		}
	}
}

bool setup_hinge_constraint(CLUSTER *pc) {
	if (pc->parent == NULL) return false;
	if (pc->nAtoms < 3 || pc->parent->nAtoms < 3) return false; // each cluster has to have 2 atom in the cluster, to construct the constraint
	CLUSTER *parent = (CLUSTER*)(pc->parent);
	MATOM *patom_00 = pc->atom + pc->hinge[pc->n2Parent].indx;
	MATOM *patom_10 = parent->atom + parent->hinge[pc->hinge[pc->n2Parent].indx_link].indx;
	MATOM *patom_01 = NULL, *patom_02 = NULL;
	MATOM *patom_11 = NULL, *patom_12 = NULL;
	bool status = false, status1 = true;
	int i, j;
	for (i = 0; i < patom_00->nb; i++) {
		if (patom_00->ba[i] == patom_10) continue;
		status1 = false;
		for (j = 0; j < pc->nAtoms; j++) {
			if ((MATOM*)(patom_00->ba[i]) == pc->atom + j) {status1 = true; break;}
		}
		if (!status1) continue; // NOT in cluster pc
		if (patom_01 == NULL) {
			patom_01 = (MATOM*)(patom_00->ba[i]); continue;
		}
		else if (patom_02 == NULL) {
			patom_02 = (MATOM*)(patom_00->ba[i]); break;
		}
	}
	if (patom_01 == NULL || patom_02 == NULL) return false;

	
	for (i = 0; i < patom_10->nb; i++) {
		if (patom_10->ba[i] == patom_00) continue;
		status1 = false;
		for (j = 0; j < parent->nAtoms; j++) {
			if ((MATOM*)(patom_10->ba[i]) == parent->atom + j) {status1 = true; break;}
		}
		if (!status1) continue; // NOT in cluster parent
		if (patom_11 == NULL) {
			patom_11 = (MATOM*)(patom_10->ba[i]); continue;
		}
		else if (patom_12 == NULL) {
			patom_12 = (MATOM*)(patom_10->ba[i]); break;
		}
	}
	if (patom_11 == NULL || patom_12 == NULL) return false;

	if (pc->hinge_constraint == NULL) pc->hinge_constraint = new HINGE_CONSTRAINT;
	pc->hinge_constraint->patom_00 = patom_00;
	pc->hinge_constraint->patom_10 = patom_10;
	pc->hinge_constraint->patom_01 = patom_01;
	pc->hinge_constraint->patom_02 = patom_02;
	pc->hinge_constraint->patom_11 = patom_11;
	pc->hinge_constraint->patom_12 = patom_12;
	return true;
}

bool MMOLECULE::setup_hinge_constraint() {
	int ic;
	char msg[256] = "\0";
	for (ic = 0; ic < nCluster; ic++) {
		if (cluster + ic == base) continue;
		if (!::setup_hinge_constraint(cluster +ic)) {
			sprintf(msg, "failure to setup hinge-constraint for cluster #%d, check whether the hinge atoms has 2 conncected atoms in both clusters of the parent hinge", ic);
			show_log(msg, true); return false;
		}
	}
	return true;
}

void CLUSTER::setup_hinge_dihedral() {
	if (parent == NULL) return;

	MATOM *patom1, *patom2, *patom3, *patom4;
	HINGE< _BASIC_CLUSTER<MATOM> > *hinge;
	int ia1, ia4;
	CHAIN< DIHEDRAL<MATOM, TORSION_PARS> > *dh_ch = NULL, *tail_ch = NULL, *tch = NULL;
	DIHEDRAL<MATOM, TORSION_PARS> *dh = NULL;
	TORSION_PARS *tpar = NULL;
	int ndh;

	ndh = 0;
	hinge = this->hinge + this->n2Parent;
	patom2 = this->atom + hinge->indx;
	patom3 = this->parent->atom + hinge->indx_link;
	for (ia1 = 0; ia1 < patom2->nb; ia1++) {
		patom1 = (MATOM*)(patom2->ba[ia1]);
		if (patom1 == patom3) continue;
		for (ia4 = 0; ia4 < patom3->nb; ia4++) {
			patom4 = (MATOM*)(patom3->ba[ia4]);
			if (patom4 == patom2) continue;
			if ((tpar = dh_db.get_torsion(patom1->aindx, patom2->aindx, patom3->aindx, patom4->aindx)) == NULL || tpar->v0 == 0) continue;
			dh = new DIHEDRAL<MATOM, TORSION_PARS>;
			dh->dh[0] = patom1; dh->dh[1] = patom2;
			dh->dh[2] = patom3; dh->dh[3] = patom4;
			dh->tpar = tpar;
			tch = new CHAIN< DIHEDRAL<MATOM, TORSION_PARS> >;
			tch->p = dh;
			if (dh_ch == NULL) dh_ch = tch;
			else tail_ch->next = tch; 
			tail_ch = tch;
		}
	}
	chain_array_context_cp< DIHEDRAL<MATOM, TORSION_PARS> >(&dh_ch, this->dh);
	release_chain< DIHEDRAL<MATOM, TORSION_PARS> >(&dh_ch, true);
	tail_ch = NULL; tch = NULL; dh = NULL;
}

void MMOLECULE::setup_cluster_matom() {
	ARRAY< ARRAY<bool> > status;
	CLUSTER *pc, *pc1;
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge;
	int i, j, ic, ia, ihinge, ic1, ia1;
	MATOM *patom, *patom1;
	bool *pstatus;
	status.SetArray(nCluster);
	for (i = 0; i < nCluster; i++) {
		status.m[i].SetArray(cluster[i].nAtoms);
		pstatus = status.m[i].m;
		for (j = 0; j < cluster[i].nAtoms; j++) pstatus[j] = true;
	}

	int nm = 0; // the atom involving in each cluster -- gatom
	char msg[256] = "\0";
	for (i = 0; i < indx_tree.n; i++) {
		for (j = 0; j < indx_tree.m[i].n; j++) {
			ic = indx_tree.m[i].m[j];
			pc = cluster + ic; nm = 0;
			for (ia = 0; ia < pc->nAtoms; ia++) {
				if (status.m[ic].m[ia]) nm++;
			}
			nm += pc->nHinges;
			if (pc->parent != NULL) nm -= 1; // the atom connecting to parent cluster will be included in its parent
			if (nm < 1) {
				sprintf(msg, "STRANGE: m %d cluster %d has NO atom with mass", this->mIndx, pc->cmIndx); show_msg(msg); continue;
			}
			pc->gatom.SetArray(nm); nm = 0;
			for (ia = 0; ia < pc->nAtoms; ia++) {
				if (status.m[ic].m[ia]) {
					pc->gatom.m[nm] = pc->atom + ia; nm++;
				}
			}

			for (ihinge = 0; ihinge < pc->nHinges; ihinge++) {
				if (ihinge == pc->n2Parent) continue;
				pc1 = (CLUSTER*)(pc->hinge[ihinge].plink);
				ic1 = pc1->cmIndx;
				ia1 = pc->hinge[ihinge].indx_link;
				patom1 = pc1->atom + ia1;
				status.m[ic1].m[ia1] = false;

				pc->gatom.m[nm] = patom1; nm++;
			}

			//sprintf(msg, "cluster %d has %d matom", ic, pc->matom.n); show_log(msg, true);
		}
	}

	// matom of each cluster
	for (ic = 0; ic < nCluster; ic++) {
		pc = cluster + ic; nm = 0;
		for (ia = 0; ia < pc->gatom.n; ia++) {
			patom = pc->gatom.m[ia];
			if (patom->par->m != 0) nm++;
		}
		pc->matom.SetArray(nm);
		nm = 0;
		for (ia = 0; ia < pc->gatom.n; ia++) {
			patom = pc->gatom.m[ia];
			if (patom->par->m != 0) {pc->matom.m[nm] = patom; nm++;}
		}
	}

	//int na = 0; nm = 0;
	//for (ic = 0; ic < nCluster; ic++) {na += cluster[ic].nAtoms; nm += cluster[ic].matom.n;};
	//sprintf(msg, "total atom %d, total used atom : %d", na, nm); show_log(msg, true);
}

void BASIC_CLUSTER::setup_cluster_matom() {
	int nm = 0, ia;
	for (ia = 0; ia < nAtoms; ia++) {
		if (atom[ia].par->m != 0) nm++;
	}
	gatom.SetArray(nAtoms);
	matom.SetArray(nm);
	nm = 0; M = 0;
	for (ia = 0; ia < nAtoms; ia++) {
		gatom.m[ia] = atom + ia;
		if (atom[ia].par->m != 0) {matom.m[nm] = atom + ia; M += atom[ia].par->m; nm++;}
	}
}

void BASIC_CLUSTER::setup_cluster_fatom() { // be sure the eneutral is checked
	int nf = 0, ia;
	MATOM *patom;
	for (ia = 0; ia < gatom.n; ia++) {
		patom = gatom.m[ia];
		if (patom->par->epsLJ != 0 || !patom->eNeutral) nf++;
	}
	fatom.SetArray(nf);
	nf = 0;
	for (ia = 0; ia < gatom.n; ia++) {
		patom = gatom.m[ia];
		if (patom->par->epsLJ != 0 || !patom->eNeutral) {
			fatom.m[nf] = patom; nf++;
		}
	}
}



void assignJobIndx(JobIndx *job, int nThreads, int n1, int n2) {
	int iThread, nleft = n2 - n1 + 1, njleft = nThreads;
	int nstep = nleft / njleft;
	int indx1 = n1, indx2;
	for (iThread = 0; iThread < nThreads; iThread++) {
		nleft = n2 - indx1 + 1; njleft = nThreads - iThread;
		nstep = nleft / njleft; if (nstep * njleft < nleft) nstep += 1;
		indx2 = indx1 + nstep - 1;

		job[iThread].n1 = indx1; job[iThread].n2 = indx2;

		indx1 += nstep;
	}
}

void assignJobIndx1(JobIndx *job, int &nThreads, int n1, int n2, int nStep, int nMaxThreads) {
	int iThread, nleft = n2 - n1 + 1, njleft = nThreads;
	int nstep = nleft / njleft;

	int nthreads = nThreads;
	if (nThreads < 1) {
		nThreads = nMaxThreads;
		if (nstep < nStep) nstep = nStep;
		if (nleft < nstep * (nThreads + 1)) nthreads = nleft / nstep;
		if (nthreads < 1) nthreads = 1;
		nThreads = nthreads;
		nstep = nleft / nThreads;
	}

	int indx1 = n1, indx2;
	for (iThread = 0; iThread < nThreads; iThread++) {
		nleft = n2 - indx1 + 1; njleft = nThreads - iThread;
		nstep = nleft / njleft; if (nstep * njleft < nleft) nstep += 1;
		indx2 = indx1 + nstep - 1;

		job[iThread].n1 = indx1; job[iThread].n2 = indx2;

		indx1 += nstep;
	}
}
