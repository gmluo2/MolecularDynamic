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
#include "MM.h"

#include "var.h"

#include "Mol.h"

bool cp(SMOLECULE *d, SMOLECULE *s, bool setup) {
	d->mType = s->mType; d->mIndx = s->mIndx;  strcpy(d->mol, s->mol);
	if (d->c == NULL) d->c = new BASIC_CLUSTER;
	else d->c->reset();
	bool status = cp(d->c, s->c);
	if (!status) return false;
	d->c->vp = &(d->r);
	d->c->vcm = &(d->r);
	int i = 0;
	for (i = 0; i < 3; i++) {
		d->r.v[i] = s->r.v[i];
		d->xaxis.v[i] = s->xaxis.v[i]; 
		d->zaxis.v[i] = s->zaxis.v[i];
	}
	if (setup) {
		if (!d->c->linkPars(::atompar_db, errmsg)) {show_infor(errmsg); return false;}
		d->check_atom_cordinate();
		d->c->setup_cluster_matom();
	}
	return true;
}

// set cluster's mass center at mass center
void SMOLECULE::check_atom_cordinate() {
	if (c == NULL) return;
	int na = 0, i = 0;
	float m = 0;
	V3zero(r)
	for (na = 0; na < c->nAtoms; na++) {
		m += c->atom[na].par->m;
		for (i = 0; i < 3; i++) r.v[i] += c->atom[na].r.v[i] * c->atom[na].par->m;
	}
	for (i = 0; i < 3; i++) r.v[i] /= m;
	for (na = 0; na < c->nAtoms; na++) {
		for (i = 0; i < 3; i++) {c->atom[na].r0.v[i] = c->atom[na].r.v[i] - r.v[i];}
	}
	c->Op = &r;
	c->vp = &r;
	c->vcm = &r;
	c->M = m;
}

void rotate_SM(SMOLECULE& m, R &r) {
	if (m.c == NULL) return;
	int i = 0, n = 0;
	VECTOR3 v;
	for (n = 0; n < m.c->nAtoms; n++) {
		//for (i = 0; i < 3; i++) v.v[i] = m.c->atom[n].r0.v[i]; // atom r
		V32V3(m.c->atom[n].r0, v)
		//m.c->atom[n].r0 = r * v;
		MpV3(r, v, m.c->atom[n].r0)

		//for (i = 0; i < 3; i++) m.c->atom[n].r.v[i] = m.c->atom[n].r0.v[i] + m.r0.v[i];
		V3plusV3(m.c->atom[n].r0, m.r, m.c->atom[n].r)  // Op is at r0, mass center

		m.c->atom[n].rotate1(r);
	}
	// hinges -- is there hinges ? NO!
	for (n = 0; n < m.c->nHinges; n++) {
		//for (i = 0; i < 3; i++) v.v[i] = m.c->hinge[n].vHinge.v[i]; // hinge vector
		V32V3(m.c->hinge[n].vHinge, v)
		//m.c->hinge[n].vHinge = r * v;
		MpV3(r, v, m.c->hinge[n].vHinge)
	}

	//for (i = 0; i < 3; i++) v.v[i] = m.xaxis.v[i];
	V32V3(m.xaxis, v)
	MpV3(r, v, m.xaxis)
	//for (i = 0; i < 3; i++) v.v[i] = m.zaxis.v[i];
	V32V3(m.zaxis, v)
	MpV3(r, v, m.zaxis)
}

void rotate_SM(SMOLECULE &m, VECTOR3 &axis, double omega, bool arc) {
	if (FABS(omega) < MIN_ROT) return;
	R r;
	RMatrix(axis, omega, r, arc);
	rotate_SM(m, r);
}

void rotate_SM(SMOLECULE &m, VECTOR3 &drot, bool arc) {
	VECTOR3 axis;
	int i = 0;
	double domega = 0;
	for (i = 0; i < 3; i++) domega += drot.v[i] * drot.v[i];
	if (domega < MIN_ROT) return;
	domega = sqrt(domega);
	for (i = 0; i < 3; i++) axis.v[i] = drot.v[i] / domega;
	
	R r;
	RMatrix(axis, domega, r, arc);
	rotate_SM(m, r);
}

void SMOLECULE::calInertiaMoment() {
	if (c == NULL) return;
	int ni = 0, nj = 0;
	int i = 0, j = 0, k = 0;

	MATRIX<3> m1, m2;
	int na = 0;
	double mass, rsize;

	M3zero(this->I)
	for (na = 0; na < c->nAtoms; na++) {
		u2Matrix3(c->atom[na].r0, m1.m)
		MpM(m1, m1, m2, i, j, k, 3)
		mass = c->atom[na].par->m;
		rsize = c->atom[na].par->rLJ;
		for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
			this->I.m[i][j] -= c->atom[na].par->m * m2.m[i][j];
			if (i == j) this->I.m[i][j] += 0.4 * mass * rsize * rsize;  // this term is the inertia tensor of sphere relative to mass center -- parallel-axis theorem for rotational inertia
		}}
	}

	InverseMatrix<3>(this->I, this->invI);
}

void SMolMove(SMOLECULE &sm) { // using the parameter in MOL_BASE_MOVEMENT (mbm)
	sm.shiftMol(sm.mbm.dr);

	VECTOR3 xaxis, raxis;
	R Rz, Rx, R0;
	int i, j, k;
	double domega = 0;

	if (sm.mbm.bRot) rotate_SM(sm, sm.mbm.drot, true);
	else {
		cp(sm.zaxis, sm.mbm.zaxis, raxis);
		domega = angle(sm.zaxis, sm.mbm.zaxis);
		RMatrix(raxis, domega, Rz, true);

		MpV3(Rz, sm.xaxis, xaxis); // xaxis of molecule changes after Rz operation
		cp(xaxis, sm.mbm.xaxis, raxis);
		domega = angle(xaxis, sm.mbm.xaxis);
		RMatrix(raxis, domega, Rx, true);

		MpM(Rx, Rz, R0, i, j, k, 3)
		rotate_SM(sm, R0);
	}
}

double calKineticEnergy(SMOLECULE &m) {
	int i = 0, j = 0;
	double Ek = 0;
	double v2 = 0;
	for (i = 3; i < 6; i++) v2 += m.c->bkm->V0.v[i] * m.c->bkm->V0.v[i];
	Ek += 0.5 * m.c->M * v2;

	m.Ekt = Ek * ::unit_Ek_kT;

	VECTOR3 w, Iw;
	double Ew = 0;
	//for (i = 0; i < 3; i++) w.v[i] = m.c->bkm.V0.v[i];
	V32V3(m.c->bkm->V0, w)
	MpV3(m.I, w, Iw)
	scale_uv3(w, Iw, Ew)
	Ek += Ew * 0.5;

	m.Ek = Ek * ::unit_Ek_kT;

	return m.Ek;
}

void InitVelocity(SMOLECULE &m) {
	int i = 0, j = 0;
	float unit_E = kT / U_InertiaMoment;
	double Ethermal = 1.5 * unit_E; // translation energy 1.5 kT; rotation energy 1.5 kT
	
	double v0 = sqrt(Ethermal * 2 / m.c->M);
	VECTOR3 F, Torq;
	double f = 0, t = 0;
	if (m.c->dyn != NULL) {
		for (i = 0; i < 3; i++) {F.v[i] = m.c->dyn->fc.v[i+3]; Torq.v[i] = m.c->dyn->fc.v[i];}
		f = F.abs(); t = Torq.abs();
	}
	if (f == 0) rand_vect3(m.c->bkm->V0.v + 3); // translation
	else {
		for (i = 0; i < 3; i++) m.c->bkm->V0.v[i+3] = F.v[i] / f;
	}
	for (i = 3; i < 6; i++) {
		m.c->bkm->V0.v[i] *= v0;
		m.c->bkm->P0.v[i] = m.c->M * m.c->bkm->V0.v[i];
	}

	VECTOR3 Iw, w;
	double I = 0, w0 = 0;
	if (t == 0) rand_vect3(w.v); // rotation
	else {
		for (i = 0; i < 3; i++) w.v[i] = Torq.v[i] / t;
	}
	MpV3(m.I, w, Iw)
	scale_uv3(w, Iw, I)
	w0 = sqrt(Ethermal * 2 / I);
	for (i = 0; i < 3; i++) {
		w.v[i] *= w0;
		m.c->bkm->V0.v[i] = w.v[i];
	}

	calKineticEnergy(m);
	MpV3(m.I, w, Iw)
	memcpy(m.c->bkm->P0.v, Iw.v, SIZE_V3);
	return;
}

void ScaleVelocity(SMOLECULE &m, double c) {
	m.c->bkm->V0 *= c;
	m.c->bkm->P0 *= c;
	m.Ek *= c * c;
}

double ScaleKineticEnergy(SMOLECULE &m, double Ethermal) {
	VECTOR<6> V0;
	V62V6(m.c->bkm->V0, V0);
	VECTOR3 v, w;
	V2wv(w, v, V0)

	Ethermal /= ::unit_Ek_kT * 2; // translation energy 1.5 kT; rotation energy 1.5 kT

	VECTOR3 axis;
	double vabs = v.abs();
	if (vabs > 1e-6) {
		axis.v[0] = v.v[0] / vabs;
		axis.v[1] = v.v[1] / vabs;
		axis.v[2] = v.v[2] / vabs;
	}
	else rand_vect3(axis.v);

	double v0 = sqrt(Ethermal * 2 / m.c->M);
	v.v[0] = axis.v[0] * v0;
	v.v[1] = axis.v[1] * v0;
	v.v[2] = axis.v[2] * v0;
	
	VECTOR3 Iw;
	double I = 0, w0 = 0;
	
	vabs = v.abs();
	if (vabs > 1e-6) {
		axis.v[0] = v.v[0] / vabs;
		axis.v[1] = v.v[1] / vabs;
		axis.v[2] = v.v[2] / vabs;
	}
	else rand_vect3(axis.v);

	MpV3(m.I, axis, Iw)
	scale_uv3(axis, Iw, I)
	w0 = sqrt(Ethermal * 2 / I);
	w.v[0] = axis.v[0] * w0;
	w.v[1] = axis.v[1] * w0;
	w.v[2] = axis.v[2] * w0;

	wv2V(w, v, m.c->bkm->V0)

	memcpy(m.c->bkm->P0.v, Iw.v, SIZE_V3);
	m.c->bkm->P0.v[3] = m.c->M * v.v[0];
	m.c->bkm->P0.v[4] = m.c->M * v.v[1];
	m.c->bkm->P0.v[5] = m.c->M * v.v[2];
	calKineticEnergy(m);
	return m.Ek;
}

void calAccel(SMOLECULE &m, float ksi) { // Op of the cluster is at mass center
	int i = 0, j = 0;
	VECTOR3 a0, alpha0;
	VECTOR3 torque;

	double *fc = m.c->dyn->fc.v, *fb = m.c->dyn->f_Brown.v;
	double *f_hoover = m.c->dyn->f_Hoover.v, *f = m.c->dyn->f.v;

	for (i = 0; i < 6; i++) a0.v[i] = -ksi * m.c->bkm->P0.v[i];
	for (i = 0; i < 6; i++) f[i] = fc[i] + fb[i] + f_hoover[i];

	torque.v[0] = f[0];
	torque.v[1] = f[1];
	torque.v[2] = f[2];

	MpV3(m.invI, torque, alpha0)

	for (i = 0; i < 3; i++) a0.v[i] = f[i+3] / m.c->M;

	wv2V(alpha0, a0, m.c->bkm->alpha)
}

void calAccel_f(SMOLECULE &m) { // Op of the cluster is at mass center
	int i = 0, j = 0;
	VECTOR3 a0, alpha0;
	VECTOR3 fc, torque;
	V2wv(torque, fc, m.c->dyn->f)

	MpV3(m.invI, torque, alpha0)

	for (i = 0; i < 3; i++) a0.v[i] = fc.v[i] / m.c->M;

	wv2V(alpha0, a0, m.c->bkm->alpha)
}

bool atom_indx(SMOLECULE &sm, MATOM* patom, int &na, int &indx) {
	na = 0; indx = 0;
	BASIC_CLUSTER *pc = sm.c;
	if (pc == NULL || pc->nAtoms == 0 || na >= pc->nAtoms) return false;
	for (na = 0; na < pc->nAtoms; na++) {
		if (pc->atom + na == patom) return true;
		indx++;
	}
	return false;
}




void CalSMElectroMultipole(SMOLECULE *sm, int nSM, bool bEwaldSum) {
	BASIC_CLUSTER *pc = NULL;
	for (int iMol = 0; iMol < nSM; iMol++) {
		pc = sm[iMol].c;
		pc->calElectroMultipole(bEwaldSum);
	}
}




// FUNCTIONS FOR PMOLECULE

double calKineticEnergy(PMOLECULE &m) {
	int i = 0, j = 0;
	double Ek = 0;
	double v2 = 0;
	V3ABS2(m.v, v2)
	Ek = 0.5 * m.c->M * v2;
	m.Ek = Ek * ::unit_Ek_kT;
	m.Ekt = m.Ek; // no rotation

	return m.Ek;
}

void InitVelocity(PMOLECULE &m) {
	float unit_E = kT / U_InertiaMoment;
	double Ethermal = 1.5 * unit_E; // translation energy 1.5 kT; rotation energy 1.5 kT
	
	double v0 = sqrt(Ethermal * 2 / m.c->M);

	int i;
	VECTOR3 F;
	double f = 0;
	if (m.c->dyn != NULL) {
		for (i = 0; i < 3; i++) {F.v[i] = m.c->dyn->fc.v[i+3];}
		f = F.abs();
	}
	if (f == 0) rand_vect3(m.v.v ); // translation
	else {
		for (i = 0; i < 3; i++) m.v.v[i] = F.v[i] / f;
	}

	m.v.v[0] *= v0; m.v.v[1] *= v0; m.v.v[2] *= v0;

	m.c->bkm->P0.v[3] = m.c->M * m.v.v[0];
	m.c->bkm->P0.v[4] = m.c->M * m.v.v[1];
	m.c->bkm->P0.v[5] = m.c->M * m.v.v[2];

	calKineticEnergy(m);
	return;
}

void ScaleVelocity(PMOLECULE &m, double c) {
	m.v *= c;
	m.c->bkm->P0 *= c;
	m.Ek *= c * c;
}

double ScaleKineticEnergy(PMOLECULE &m, double Ek0) {
	float unit_E = kT / U_InertiaMoment;
	double Ethermal = 1.5 * unit_E; // translation energy 1.5 kT; rotation energy 1.5 kT
	
	double v = sqrt(Ethermal * 2 / m.c->M);
	double v0 = m.v.abs();
	if (v0 < 1e-6) {
		rand_vect3(m.v.v);
		m.v.v[0] *= v; m.v.v[1] *= v; m.v.v[2] *= v;
	}
	else {
		v /= v0;
		m.v.v[0] *= v; m.v.v[1] *= v; m.v.v[2] *= v;
	}

	m.c->bkm->P0.v[3] = m.c->M * m.v.v[0];
	m.c->bkm->P0.v[4] = m.c->M * m.v.v[1];
	m.c->bkm->P0.v[5] = m.c->M * m.v.v[2];

	double Ek = calKineticEnergy(m);
	return Ek;
}

void calAccel(PMOLECULE &m, float ksi) {
	double M = m.c->M;
	double *fc = m.c->dyn->fc.v, *fb = m.c->dyn->f_Brown.v;
	double *f_hoover = m.c->dyn->f_Hoover.v, *f = m.c->dyn->f.v;
	if (ksi != 0) {
		f_hoover[3] = -ksi * m.v.v[0] * M;
		f_hoover[4] = -ksi * m.v.v[1] * M;
		f_hoover[5] = -ksi * m.v.v[2] * M;
	}
	f[3] = fc[3] + fb[3] + f_hoover[3];
	f[4] = fc[4] + fb[4] + f_hoover[4];
	f[5] = fc[5] + fb[5] + f_hoover[5];

	m.alpha.v[0] = f[3] / M;
	m.alpha.v[1] = f[4] / M;
	m.alpha.v[2] = f[5] / M;
}

void calAccel_f(PMOLECULE &m) {
	double M = m.c->M;
	m.alpha.v[0] = m.c->dyn->f.v[3] / M;
	m.alpha.v[1] = m.c->dyn->f.v[4] / M;
	m.alpha.v[2] = m.c->dyn->f.v[5] / M;
}

void polarize_SM(SMOLECULE *sm, int n) {
	int imol, iatom;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	double *mu = NULL, *E = NULL, *mut = NULL, *mu0 = NULL;
	ATOM_COMM_PARS *par = NULL;
	double fp = ::fp;
	double d = 0, alpha = 0;
	for (imol = 0; imol < n; imol++) {
		pc = sm[imol].c;
		if (pc->eNeutral) continue;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom; 
			if (patom->mMP > 0) {
				mu = patom->ind_mu.v; E = patom->E.v; par = patom->par; alpha = par->alpha;
				d = alpha * E[0]; mu[0] = mu[0] + fp * (d - mu[0]);
				d = alpha * E[1]; mu[1] = mu[1] + fp * (d - mu[1]);
				d = alpha * E[2]; mu[2] = mu[2] + fp * (d - mu[2]);

				mu0 = patom->intr_mu.v; mut = patom->mu.v;
				mut[0] = mu0[0] + mu[0]; mut[1] = mu0[1] + mu[1]; mut[2] = mu0[2] + mu[2];
			}
		}
	}
}

void polarize_VSM(SMOLECULE **sm, int n) {
	int imol, iatom;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	double *mu = NULL, *E = NULL, *mut = NULL, *mu0 = NULL;
	ATOM_COMM_PARS *par = NULL;
	double fp = ::fp;
	double d = 0, alpha = 0;
	for (imol = 0; imol < n; imol++) {
		pc = sm[imol]->c;
		if (pc->eNeutral) continue;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom; 
			if (patom->mMP > 0) {
				mu = patom->ind_mu.v; E = patom->E.v; par = patom->par; alpha = par->alpha;
				d = alpha * E[0]; mu[0] = mu[0] + fp * (d - mu[0]);
				d = alpha * E[1]; mu[1] = mu[1] + fp * (d - mu[1]);
				d = alpha * E[2]; mu[2] = mu[2] + fp * (d - mu[2]);

				mu0 = patom->intr_mu.v; mut = patom->mu.v;
				mut[0] = mu0[0] + mu[0]; mut[1] = mu0[1] + mu[1]; mut[2] = mu0[2] + mu[2];
			}
		}
	}
}

// calculate effective dipole moment of multi-macromolecules
void polarize_PM(PMOLECULE *pm, int n) {
	int imol;
	MATOM *patom = NULL;
	double *mu = NULL, *E = NULL, *mut = NULL, *mu0 = NULL;
	ATOM_COMM_PARS *par = NULL;
	double fp = ::fp;
	double d = 0, alpha = 0;
	for (imol = 0; imol < n; imol++) {
		patom = pm[imol].c->atom;
		if (patom->mMP > 0) {
			mu = patom->ind_mu.v; E = patom->E.v; par = patom->par; alpha = par->alpha;
			d = alpha * E[0]; mu[0] = mu[0] + fp * (d - mu[0]);
			d = alpha * E[1]; mu[1] = mu[1] + fp * (d - mu[1]);
			d = alpha * E[2]; mu[2] = mu[2] + fp * (d - mu[2]);

			mu0 = patom->intr_mu.v; mut = patom->mu.v;
			mut[0] = mu0[0] + mu[0]; mut[1] = mu0[1] + mu[1]; mut[2] = mu0[2] + mu[2];
		}
	}
}

void polarize_VPM(PMOLECULE **pm, int n) {
	int imol;
	MATOM *patom = NULL;
	double *mu = NULL, *E = NULL, *mut = NULL, *mu0 = NULL;
	ATOM_COMM_PARS *par = NULL;
	double fp = ::fp;
	double d = 0, alpha = 0;
	for (imol = 0; imol < n; imol++) {
		patom = pm[imol]->c->atom;
		if (patom->mMP > 0) {
			mu = patom->ind_mu.v; E = patom->E.v; par = patom->par; alpha = par->alpha;
			d = alpha * E[0]; mu[0] = mu[0] + fp * (d - mu[0]);
			d = alpha * E[1]; mu[1] = mu[1] + fp * (d - mu[1]);
			d = alpha * E[2]; mu[2] = mu[2] + fp * (d - mu[2]);

			mu0 = patom->intr_mu.v; mut = patom->mu.v;
			mut[0] = mu0[0] + mu[0]; mut[1] = mu0[1] + mu[1]; mut[2] = mu0[2] + mu[2];
		}
	}
}

bool cp(PMOLECULE *d, PMOLECULE *s, bool setup) {
	d->mType = s->mType; d->mIndx = s->mIndx;  strcpy(d->mol, s->mol);
	if (d->c == NULL) {
		d->c = new BASIC_CLUSTER; d->c->vp = &(d->r); d->c->Op = &(d->r); d->c->vcm = &(d->r);
	}
	else d->c->reset();
	bool status = cp(d->c, s->c);
	if (!status) return false;

	int i = 0;
	for (i = 0; i < 3; i++) {
		d->r.v[i] = s->r.v[i];
	}
	if (setup) {
		if (!d->c->linkPars(::atompar_db, errmsg)) {show_infor(errmsg); return false;}
		d->c->M = d->c->atom->par->m;
		d->q = d->c->atom->c;
		d->c->setup_cluster_matom();
	}
	else {d->c->M = s->c->M; d->q = s->q;}

	d->check_coordinate();

	return true;
}
