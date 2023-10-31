#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "ranlib.h"

#include "def.h"
#include "show.h"
#include "vector.h"
#include "bound.h"
#include "ZMatrix.h"
#include "Matrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "cluster.h"

#define INTERACT_DEF 1
#define _INTERACTION 1
#include "Interaction.h"
#include "Interact2.h"

#include "var.h"

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

namespace _rand_ {
	GaussRand2 gauss_rand2;

	// GaussianRand1 is not a good Gaussian distribution !!!
double GaussianRand1() {
	// Gaussian Function : y = exp[-x^2 / 2 ]
	double x, w = 0;
	while (w == 0) w = ranf(); // w is [0, 1] but not 0
	//x = sqrt( (-2.0 * log( w ) ) / w);
	x = sqrt( (-2.0 * log( w ) ));
	return x - 1;
}

// polar form of Box¨CMuller transform to generate two Gaussian distributed non-correlated randoms [-1, 1]
double GaussianRand2() {
	gauss_rand2.irand++;
	if (gauss_rand2.irand == 2) {gauss_rand2.irand = 0; return gauss_rand2.v[1];}

	// Gaussian Function : y = exp[-(x / width)^2 / 2]
	double x1, x2, w = 2;
	while(w >= 1.0) {
		x1 = 2.0 * ranf() - 1.0;
		x2 = 2.0 * ranf() - 1.0;
		w = x1 * x1 + x2 * x2;
	}
	w = sqrt( (-2.0 * log( w ) ) / w );
	gauss_rand2.v[0] = x1 * w;
	gauss_rand2.v[1] = x2 * w;
	gauss_rand2.irand = 1;
	return gauss_rand2.v[0];
}

#include "distribution.h"
void test_GaussRand(float width, int max) {
	using namespace _distribution_;
	DISTRIBT<float, float> dt;
	dt.set_distribt_x(-4, 4, 80);
	int i;
	double r = 0;
	for (i = 0; i < max; i++) {
		r = GaussianRand2() * width;
		dt.sample(r, true);
	}
	dt.eb_raw();
	ofstream out;
	out.open("rand.dat");
	dt.save_x_y_eb(out);
	out.close();
}

void init_array(RAND_ARRAY<int> &randa) {
	for (int i = 0; i < randa.nw; i++) randa.buf.m[i] = i;
}

void init_random_array(RAND_ARRAY<int> &rand, int n) {
	rand.init_elements(n);
	init_array(rand);
	rand.random_order();
}

}

void test_GaussRand(float width, int max) {
	_rand_::test_GaussRand(width, max);
}

void test_rand_array(int n) {
	using namespace _rand_;
	RAND_ARRAY<int> randa;
	init_random_array(randa, n);
	ofstream out;
	out.open("rand_array.dat");
	for (int i = 0; i < randa.nw; i++) out<<i<<"  "<<randa.a.m[i]<<"  "<<randa.buf.m[randa.nw -1 - i]<<endl;
	out.close();
}

void gauss_random_force(float f0, float sigma, double *f, float ctheta) {
	using namespace _rand_;
	rand_vect3(f);
	double gf = GaussianRand2() * sigma;
	double force = f0 + gf;
	f[0] *= force; f[1] *= force; f[2] *= force;
}

void gen_add_random_force(MMOLECULE &mm, int ncs, float f0, float f_sigma, _rand_::RAND_ARRAY<int> &order) {
	SVECTOR<double, 3> f, t, dr;
	double d = 0;
	int i, j;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	for (i = 0; i < ncs;  i++) {
		if (i >= order.nw) return;
		if (order.a.m[i] >= mm.nCluster) continue;
		pc = mm.cluster + order.a.m[i];
		// random force and torque on each cluster
		gauss_random_force(f0, f_sigma, f.v, 1);

		VECT3((*(pc->Op)), pc->cm, dr)
		V3PV3(dr, f, t)
		
		for (j = 0; j < 3; j++) {pc->dyn->fc.v[j] += t.v[j]; pc->dyn->fc.v[3 + j] += f.v[j];}
	}
}

void gen_add_random_force(MMOLECULE &mm, int ncs, float f0, float f_sigma, float torq, float t_sigma, _rand_::RAND_ARRAY<int> &order) {
	SVECTOR<double, 3> f, t, dr;
	double d = 0;
	int i, j;
	CLUSTER *pc = NULL, *parent;
	MATOM *patom = NULL;
	for (i = 0; i < ncs;  i++) {
		if (i >= order.nw) return;
		if (order.a.m[i] >= mm.nCluster) continue;
		pc = mm.cluster + order.a.m[i];

		if (pc->parent == NULL) {
			// random force and torque on each cluster
			gauss_random_force(f0, f_sigma, f.v, 1);
			VECT3((*(pc->Op)), pc->cm, dr)
			V3PV3(dr, f, t)

			for (j = 0; j < 3; j++) {pc->dyn->fc.v[3 + j] += f.v[j]; pc->dyn->fc.v[j] += t.v[j];}
		}
		else { // this is a torque
			//gauss_random_force(torq, t_sigma, t.v, 1);
			d = torq + _rand_::GaussianRand2() * t_sigma;
			d *= (ranf() > 0.5 ? 1 : -1);
			t.v[0] = d * pc->km.H.v[0];
			t.v[1] = d * pc->km.H.v[1];
			t.v[2] = d * pc->km.H.v[2];

			V3plus(pc->dyn->fc.v, t.v, pc->dyn->fc.v)

			parent = (CLUSTER*)(pc->parent);
			INV_V3(t)
			V3plus(parent->dyn->fc.v, t.v, parent->dyn->fc.v)
		}
		
		
	}
}

extern EF_DPROF LJ12_6;
// LJ 12-6 profiles has energy unit in kT, and force is [atomic mass, Angstrom, fs]
void construct_Alternative_LJ12_6_prof() {
	/****************************************************************
		Alternative format of L-J (used in Amber): 
			U = U0 * (1/x^12 - 2 / x^6)
	*****************************************************************/
	double ELJ_unit = U_eps / kT; // unit kT
	double FLJ_unit = U_LJForce / U_InertiaMoment; // unit mass
	double dx = 0.01, x2, x6;
	int npts = 401, i; 
	LJ12_6.set(npts, dx, 0);
	for (i = 2; i < npts; i++) {
		LJ12_6.x[i] = i * dx; x2 = LJ12_6.x[i] * LJ12_6.x[i]; x6 = 1 / (x2 * x2 * x2);
		LJ12_6.y[0][i] =  x6 * x6 - 2 * x6; LJ12_6.y[0][i] *= ELJ_unit;
		LJ12_6.y[1][i] = 12 * (x6 * x6 - x6) / LJ12_6.x[i]; LJ12_6.y[1][i] *= FLJ_unit;
	}
	for (i = 0; i < 2; i++) {LJ12_6.y[0][i] = LJ12_6.y[0][2]; LJ12_6.y[1][i] = LJ12_6.y[1][2];}

	// cut-off distance
	int icut = int(rcut_LJ_norm / dx + 0.5); icut = (icut >= LJ12_6.npts ? LJ12_6.npts - 1 : icut);
	double Ucut = LJ12_6.y[0][icut], fcut = 0; //fcut = LJ12_6.y[1][icut];
	for (i = 0; i <= icut; i++) {LJ12_6.y[0][i] -= Ucut; LJ12_6.y[1][i] -= fcut;}
	for (i = icut; i < LJ12_6.npts; i++) {LJ12_6.y[0][i] = 0; LJ12_6.y[1][i] = 0;}

	//LJ12_6.y[0][LJ12_6.npts - 1] = 0; // set the last point and infinit 0
	//LJ12_6.y[1][LJ12_6.npts - 1] = 0; // set the last point and infinit 0

	/*
	ofstream out;
	out.open("LJ12_6.dat");
	double rLJ = 3, r2LJ = 9, r = 0, r2 = 0, t1, t2;
	double V = 0;
	VECTOR3 u(0, 0, 1), f, dr;
	for (r = 1.2; r < 15; r += 0.1) {
		r2 = r * r;
		LJ_V3FORCE(1, rLJ, r, u, i, t1, t2, V, f)
		out<<r<<" "<<V<<" "<<f.v[2]<<" ";
		dr.v[2] = r;
		_LJ_V3FORCE(1, r2LJ, dr, r2, t1, f, V)
		out<<V<<" "<<f.v[2]<<endl;
	}
	out.close();
	*/
}

void construct_Standard_LJ12_6_prof() {
	/****************************************************************
		L-J : U = 4 * U0 * (1/x^12 - 1 / x^6), r* = 1.122462 
	*****************************************************************/
	double ELJ_unit = 4 * U_eps / kT; // unit kT
	double FLJ_unit = 4 * U_LJForce / U_InertiaMoment; // unit mass
	double dx = 0.01, x2, x6;
	int npts = 401, i; 
	LJ12_6.set(npts, dx, 0);
	for (i = 2; i < npts; i++) {
		LJ12_6.x[i] = i * dx; x2 = LJ12_6.x[i] * LJ12_6.x[i]; x6 = 1 / (x2 * x2 * x2);
		LJ12_6.y[0][i] =  x6 * x6 - x6; LJ12_6.y[0][i] *= ELJ_unit;
		LJ12_6.y[1][i] = (12 * x6 * x6 - 6 * x6) / LJ12_6.x[i]; LJ12_6.y[1][i] *= FLJ_unit;
	}
	for (i = 0; i < 2; i++) {LJ12_6.y[0][i] = LJ12_6.y[0][2]; LJ12_6.y[1][i] = LJ12_6.y[1][2];}

	// cut-off distance
	int icut = int(rcut_LJ_norm * 1.122462 / dx + 0.5); icut = (icut >= LJ12_6.npts ? LJ12_6.npts - 1 : icut);
	double Ucut = LJ12_6.y[0][icut], fcut = 0; //fcut = LJ12_6.y[1][icut];
	for (i = 0; i <= icut; i++) {LJ12_6.y[0][i] -= Ucut; LJ12_6.y[1][i] -= fcut;}
	for (i = icut; i < LJ12_6.npts; i++) {LJ12_6.y[0][i] = 0; LJ12_6.y[1][i] = 0;}

	//LJ12_6.y[0][LJ12_6.npts - 1] = 0; // set the last point and infinit 0
	//LJ12_6.y[1][LJ12_6.npts - 1] = 0; // set the last point and infinit 0

	/*
	ofstream out;
	out.open("LJ12_6.dat");
	double rLJ = 3, r2LJ = 9, r = 0, r2 = 0, t1, t2;
	double V = 0;
	VECTOR3 u(0, 0, 1), f, dr;
	for (r = 1.2; r < 15; r += 0.1) {
		r2 = r * r;
		LJ_V3FORCE(0.1, rLJ, r2LJ, r, u, i, t1, t2, V, f)
		out<<r<<" "<<V<<" "<<f.v[2]<<" ";
		dr.v[2] = r;
		_LJ_V3FORCE(0.1, r2LJ, dr, r2, t1, f, V)
		out<<V<<" "<<f.v[2]<<endl;
	}
	out.close();
	*/
}

void construct_LJ12_6_prof(int mode) {
	switch (mode) {
	case 0: // standard L-J 12-6
		construct_Standard_LJ12_6_prof();
		break;
	case 1: // Alternative format, used in Amber force field
		construct_Alternative_LJ12_6_prof();
		break;
	}
}

// force on the molecular mass center
// the interactions between the clusters inside the molecule,
// give NO net force on the molecular mass center
void CalNetForce_Thread(MMOLECULE *mm, int nc1, int nc2, VECTOR<6>& fc) {
	V6zero(fc)
	int nc = 0;
	CLUSTER *pc = NULL;
	VECTOR<6> ft;

	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		MpV6(pc->invNE->phai_cm2me, pc->dyn->fc, ft)
		V6plusV6(fc, ft, fc)
	}
	//return fc;
}

void CalNetForce(MMOLECULE &mm, int nc1, int nc2) {
	V6zero(mm.mDynKin.fc)
	MacroMoleculeOperate2<MMOLECULE, VECTOR<6> >((void*)(&CalNetForce_Thread), mm, nc1, nc2, MAX_THREADS, mm.mDynKin.fc);
}

#if _DISABLE_ == 0
void FreeCellClusterInteraction(MMOLECULE &mm, InteractVar &var, InteractRes &res) { // V in kT
	VECTOR<3> f, t, l;
	VECTOR<3> u, r;
	double f0 = 0, dis = 0, dis2 = 0, t1, t2;
	bool bforce = true;
	double E_LJ = 0, E_estat = 0;
	double Et_LJ = 0, Et_estat = 0;
	double Uext = 0;

	CLUSTER *p1 = NULL, *p2 = NULL;
	MATOM *patom1 = NULL, *patom2 = NULL;

	LJ_PARS *pLJ = NULL;

	int i = 0, j = 0, k = 0;
	int nc1 = 0, nc2 = 0;
	int na1 = 0, na2 = 0;

	bool dis_not_cut = true;
	bool bound = false, bound14 = false;
	
	float safe_dis = 1.2f;
	
	double U = 0;

	res.STxx = 0; res.STyy = 0; res.STzz = 0; res.STtr = 0;

	for (nc1 = 0; nc1 < mm.nCluster; nc1++) {
		p1 = mm.cluster + nc1;
		for (na1 = 0; na1 < p1->nAtoms; na1++) {
			patom1 = p1->atom + na1;

			for (nc2 = 0; nc2 < mm.nCluster; nc2++) {
				if (nc1 == nc2) continue; // in same cluster
				p2 = mm.cluster + nc2;
				for (na2 = 0; na2 < p2->nAtoms; na2++) {
					patom2 = p2->atom + na2;

					V3zero(f)
					V3zero(t)
					f0 = 0;
					pLJ = NULL;
					VECT(patom2->r, patom1->r, r, i)  // VECTOR from atom2 to atom1
					V3ABS2(r, dis2) // dis2 = |r|^2

					if (dis2 < d2_14) {
						BOUND(patom1, patom2, i, bound)
						if (bound) continue; // bound atom
						BOUND_SAME_ATOM(patom1, patom2, i, j, bound)
						if (bound) continue; // atom2 - atom - atom1 is a bound angle
						BOUND_CLOSE_ATOMS(patom1, patom2, i, j, k, bound14)
					}
					else {bound = false; bound14 = false;}

					if (dis2 < _R2min_Interact) continue; // the two atoms are too close

					dis = sqrt(dis2);
					DIV(r, dis, u, i) // u = r / dis; unit vector of r; from atom2 to atom1

					bforce = false;
					pLJ = &(LJPARS.m[patom1->aindx].m[patom2->aindx]);
					//pLJ = NULL;
					//search_LJPARS(patom1->aindx, patom2->aindx, LJPARS, nLJPARS, i, pLJ);
					if (pLJ != NULL && dis2 < r2cut_LJ && pLJ->epsLJ != 0) {
						//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, r, dis2, temp, f, E_LJ)
						LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, pLJ->rLJ2, dis, u, i, t1, t2, E_LJ, f)
						if (bound14) {E_LJ *= A14_LJ; t2 *= A14_LJ; for (i = 0; i < 3; i++) f.v[i] *= A14_LJ;}

						if (var.bVirial) {
							if (var.bST_diagonal) Accumulate_Virial_Diagonal(res, f.v, r.v);
							else Accumulate_Virial_Trace(res, t2, dis);
						}

						bforce = true;
						Et_LJ += E_LJ;
					}
					if (patom1->c != 0 && patom2->c != 0) {
						E_estat = patom1->c * patom2->c / dis;
						if (bound14) E_estat *= A14_E;
						f0 = E_estat / dis2 * fUnit_estat_mass;
						for (i = 0; i < 3; i++) f.v[i] += f0 * r.v[i]; // f -- from atom 2 ==> 1

						if (var.bVirial) {
							if (var.bST_diagonal) Accumulate_Virial_Diagonal(res, f.v, r.v);
							else Accumulate_Virial_Trace(res, f0, dis2);
						}

						bforce = true;
						Et_estat += E_estat * eUnit_estat_kT;
					}

					if (bforce) {
						V3PV3(patom1->r0, f, t)   // t = r0 x f on atom1
						for (i = 0; i < 3; i++) p1->dyn->fc.v[i] += t.v[i];
						for (i = 0; i < 3; i++) p1->dyn->fc.v[3 + i] += f.v[i];
						//V3PV3(patom2->r0, f, t)   // t = r0 x f on atom2, real torque is -t because force on atom 2 = -f
						//for (i = 0; i < 3; i++) p2->dyn->fc.v[i] -= t.v[i]; 
						//for (i = 0; i < 3; i++) p2->dyn->fc.v[3 + i] -= f.v[i];
					}
				}
			}
			if (bEext && patom1->c != 0) { // external electric field
				f.v[0] = 0; f.v[1] = 0; f.v[2] = patom1->c0 * Eext * fUnit_ext_mass;
				V3PV3(patom1->r0, f, t);
				for (i = 0; i < 3; i++) p1->dyn->fc.v[i] += t.v[i];
				for (i = 0; i < 3; i++) p1->dyn->fc.v[3 + i] += f.v[i];
				Uext += patom1->c0 * patom1->r.v[2] * Eext * eUnit_ext_kT;
			}
		}
	}

	U = (Et_LJ + Et_estat) * 0.5 + Uext;

	res.U = U;
}

// calculate the interaction energy, Electric field, gradient electric field on each atom and LJ force
void EFIELD_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, InteractVar &var, InteractRes &res) {
	MATOM *patom = NULL;
	VECTOR3 E;
	DMATRIX<3> Eg;
	VECTOR3 dr, u;
	int i, j, k, l, na, na1;
	VECTOR3 r0, r1;
	double t1, t2, t3;

	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
#if _EWALD_SUM == _MultiPole_EWALD_SUM
	CHAIN<ELECTRO_MULTIPOLE> *ch_emp = NULL;
	ELECTRO_MULTIPOLE *emp = NULL;
#endif

	RELSHIP_ITERATOR<BASIC_CLUSTER, LJ_RELSHIP> LJit;

	CHAIN<BASIC_CLUSTER> *ch_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	VECTOR3 mT1;
	DMATRIX<3> mT2, mD2;
	CMATRIX<3> mT3;
	QMATRIX<3> mT4;

	double R, R2, R3, R4, R5, R7, R9;
	VECTOR3 *h = NULL;

	bool bound = false, bound14 = false;

	double U_LJ = 0;
	VECTOR3 force, torque;
	LJ_PARS *pLJ = NULL;
	
	double Utotal = 0;
	UTotal = 0;

	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		Utotal = 0;
		for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		V3zero(E)
		Mzero(Eg, i, j)

	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 1
		// ****** ignore the electric interaction between atoms *******
		goto _EOVER_EFIELD_CMM;
		// ************************************************************
	#endif

	#if _ATOM_INDUCED_DIPOLE == 0
		if (FABS(patom->c) < 0.0001) goto _EOVER_EFIELD_CMM;
	#endif

	#if _EWALD_SUM == _MultiPole_EWALD_SUM
#error the interaction results is not realized for multipole interaction 
		// using relation-ship
		if (local_interact) ch_emp = NULL; // ignore the interaction with multipoles
		else ch_emp = pc->ERship.ch_mpl;
		while (ch_emp != NULL) {
			emp = ch_emp->p;
			//for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - emp->r0.v[i];
			//ABS2(dr, i, R2)
			VECT3(emp->r0, r0, dr)
			scale_uv3(dr, dr, R2)

			R = sqrt(R2); R3 = R2 * R; R4 = R2 * R2; R5 = R2 * R3; R7 = R5 * R2;
			MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
			MTP_T3(mT3, dr, R2, R7, i, j, k, t1)
		#if _ATOM_INDUCED_DIPOLE == 1
			R9 = R7 * R2;
			MTP_T4(mT4, dr, R2, R4, R9, i, j, k, l, t1, t2)
		#endif
			// energy
			t1 = (emp->q == 0 ? 0 : patom->c * emp->q / R); // q-q
			for (i = 0; i < 3; i++) {
				t1 += -patom->c * mT1.v[i] * emp->mu.v[i]; // q-mu
			#if _ATOM_INDUCED_DIPOLE == 1
				if (emp->q != 0) t1 += patom->ind_mu.v[i] * mT1.v[i] * emp->q; // mu-q
			#endif
				for (j = 0; j < 3; j++) {
					t1 += patom->c * mT2.m[i][j] * emp->Q.m[i][j] / 3; // q-quadrupole
				#if _ATOM_INDUCED_DIPOLE == 1
					t1 -= patom->ind_mu.v[i] * mT2.m[i][j] * emp->mu.v[j]; // mu-mu
					for (k = 0; k < 3; k++) t1 += patom->ind_mu.v[k] * mT3.m[i][j][k] * emp->Q.m[i][j] / 3; // mu-quadrupole
				#endif
				}
			}
			Utotal += t1 * eUnit_estat_kT;

			// electric field
			for (i = 0; i < 3; i++) {
				if (emp->q != 0) E.v[i] -= emp->q * mT1.v[i];
				for (j = 0; j < 3; j++) {
					E.v[i] += emp->mu.v[j] * mT2.m[i][j];
					for (k = 0; k < 3; k++) E.v[i] -= emp->Q.m[j][k] * mT3.m[i][j][k] / 3;
				}
			}
		#if _ATOM_INDUCED_DIPOLE == 1
			// electric field gradient
			for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
				if (emp->q != 0) Eg.m[i][j] -= emp->q * mT2.m[i][j];
				for (k = 0; k < 3; k++) {
					Eg.m[i][j] += emp->mu.v[k] * mT3.m[i][j][k];
					for (l = 0; l < 3; l++) Eg.m[i][j] -= emp->Q.m[k][l] * mT4.m[i][j][k][l] / 3;
				}
			}}
		#endif
			ch_emp = ch_emp->next;
		}

		ch_cluster = pc->ERship.ch_cluster;
		while (ch_cluster != NULL) {
			pc1 = ch_cluster->p;
			if (pc1 == pc) {ch_cluster = ch_cluster->next; continue;} // same cluster

			if (::local_relax) { // do parent-child interaction only
				if (pc1->parent != pc && pc->parent != pc1) {ch_cluster = ch_cluster->next; continue;}
			}

			for (na1 = 0; na1 < pc1->nAtoms; na1++) {
				patom1 = pc1->atom + na1;
			#if _ATOM_INDUCED_DIPOLE == 0
				if (FABS(patom1->c) < 0.00001) continue; // this atom has no charge
			#endif
				for (i = 0; i < 3; i++) dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i];
				//ABS2(dr, i, R2)
				scale_uv3(dr, dr, R2)

				R = sqrt(R2); R3 = R * R2; 
			#if _ATOM_INDUCED_DIPOLE == 1
				R4 = R2 * R2; R5 = R3 * R2; R7 = R5 * R2;
				MTP_T1(mT1, dr, R3) MTP_T2(mT2, dr, R2, R5, i, j) MTP_D2(mD2, dr, i, j)
				MTP_T3(mT3, dr, R2, R7, i, j, k, t1)
			#endif
				if (R2 < d2_14 && pc1->mIndx == pc->mIndx) {
					BOUND(patom, patom1, i, bound) // bound atom
					if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
					if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
					else bound14 = false;
				}
				else {bound = false; bound14 = false;}

				if (bound || pc == pc1) { // bounding atoms or they are in the same cluster
				}
				else if (bound14) { // 1-4 interaction
					// electric field
					for (i = 0; i < 3; i++) {
						t1 = patom1->c * dr.v[i] / R3;
						E.v[i] += t1 * A14_E;
					#if _ATOM_INDUCED_DIPOLE == 1
						for (j = 0; j < 3; j++) {
							E.v[i] += patom1->ind_mu.v[j] * mT2.m[i][j] * A14_E;
						}
					#endif
					}
				#if _ATOM_INDUCED_DIPOLE == 1
					// electric field gradient
					for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
						Eg.m[i][j] -= patom1->c * mT2.m[i][j] * A14_E;
						for (k = 0; k < 3; k++) {
							Eg.m[i][j] += patom1->ind_mu.v[k] * mT3.m[i][j][k] * A14_E;
						}
					}}
				#endif
					// energy
					Utotal += patom->c * patom1->c / R * A14_E * eUnit_estat_kT; // 1-4 interaction is corrected by * A14_E
				#if _ATOM_INDUCED_DIPOLE == 1
					for (i = 0; i < 3; i++) {
						Utotal += (-patom->c * patom1->ind_mu.v[i] + patom1->c * patom->ind_mu.v[i]) * 
							mT1.v[i] * A14_E * eUnit_estat_kT;
						for (j = 0; j < 3; j++) Utotal -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mT2.m[i][j] * A14_E * eUnit_estat_kT;
					}
				#endif
				}
				else {
					// electric field
					for (i = 0; i < 3; i++) {
						E.v[i] += patom1->c * dr.v[i] / R3;
				#if _ATOM_INDUCED_DIPOLE == 1
						for (j = 0; j < 3; j++) E.v[i] -= patom1->ind_mu.v[j] * mT2.m[i][j];
				#endif
					}
				#if _ATOM_INDUCED_DIPOLE == 1
					// electric field gradient
					for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
						Eg.m[i][j] -= patom1->c * mT2.m[i][j];
						for (k = 0; k < 3; k++) Eg.m[i][j] += patom1->ind_mu.v[k] * mT3.m[i][j][k];
					}}
				#endif
					// energy
					Utotal += patom->c * patom1->c / R * eUnit_estat_kT;
				#if _ATOM_INDUCED_DIPOLE == 1
					for (i = 0; i < 3; i++) {
						Utotal += (patom1->c * patom->ind_mu.v[i] - patom1->ind_mu.v[i] * patom->c) * mT1.v[i] * eUnit_estat_kT;
						for (j = 0; j < 3; j++) Utotal -= patom->ind_mu.v[i] * patom1->ind_mu.v[j] * mT2.m[i][j] * eUnit_estat_kT;
					}
				#endif
				}
			}
			ch_cluster = ch_cluster->next;
		}
	#endif  //#if _EWALD_SUM == _MultiPole_EWALD_SUM

_EOVER_EFIELD_CMM:
		for (i = 0; i < 3; i++) {
			patom->E.v[i] = E.v[i];
#if _ATOM_INDUCED_DIPOLE == 1
			for (j = 0; j < 3; j++) patom->Eg.m[i][j] = Eg.m[i][j];
#endif
		}

		// LJ interaction
		if (patom->par->epsLJ > 0.0001) {
			//ch_imag_cluster = pc->LJRship.ch;
			//while (ch_imag_cluster != NULL) {
			LJit.set(&(pc->LJRship)); LJit.start_iterate();
			while (!LJit.eoi()) {
				imag_cluster = LJit.current();
				if (imag_cluster == NULL) {LJit.next_iterate(); continue;}
				pc1 = imag_cluster->p;
				if (pc1 == pc) {
					//ch_imag_cluster = ch_imag_cluster->next; 
					LJit.next_iterate();
					continue;
				}

				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					patom1 = pc1->atom + na1;
					if (patom1->par->epsLJ < 0.0001) continue;
					for (i = 0; i < 3; i++) {
						dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i];
					}
					//ABS2(dr, i, R2)
					R2 = dr.v[0] * dr.v[0]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[1] * dr.v[1]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[2] * dr.v[2]; if (R2 > r2cut_LJ) continue;

					if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}
					if (bound) continue; // ignore bound interaction

					if (R2 < _R2min_Interact) continue; // the two atoms are too close

					pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					if (pLJ != NULL && pLJ->epsLJ > 0.0001) {
						//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
						R = sqrt(R2); u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;
						LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ, pLJ->rLJ2, R, u, i, t1, t3, t2, force)
						if (bound14) { // 1-4 interaction, LJ /= 2
							t2 *= A14_LJ;
							for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
						}
						if (var.bVirial) {
							if (var.bST_diagonal) Accumulate_Virial_Diagonal(res, force.v, dr.v);
							else Accumulate_Virial_Trace(res, t3, R);
						}
						U_LJ += t2 * 0.5; // all clusters are in same system
						//U_LJ += t2; // will be corrected later
						V3PV3(patom->r0, force, torque)
						for (i = 0; i < 3; i++) {pc->dyn->fc.v[i+3] += force.v[i]; pc->dyn->fc.v[i] += torque.v[i];}
					}
				}
				//ch_imag_cluster = ch_imag_cluster->next;
				LJit.next_iterate();
			}
		}

		// energy
		UTotal += Utotal;
	}
	// total energy in kT
	UTotal = UTotal + U_LJ;

	res.U = UTotal
}

void EFIELD_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, InteractVar& var, InteractRes &res) {
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;
	InteractRes tres;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			EFIELD_CMM(cell, pc, var, tres);
			res += tres;
			it.next_iterate();
		}
	}
	Utotal *= 0.5; // electrostatic and LJ are pairwise interaction
}

void FreeCellCMMClusterInteraction(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, InteractVar &var, InteractRes &res) { // V in kT
	int nc, na, i;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	InteractRes &tres;
#error needs to be checked. Force is better to be calculated for each atom, while not accumulate with the electric field

	float coeff = sqrt(::eps_dc);

	// calculate the electric field on each atom in these cells
	EFIELD_CMM(cell, n1Cell, n2Cell, V, OVERLAP);

	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			for (na = 0; na < pc->nAtoms; na++) {
				patom = pc->atom + na;
				if (patom->c != 0) {
					for (i = 0; i < 3; i++) {
						F.v[i] = patom->c * patom->E.v[i] * fUnit_estat_mass;
						if (bEext) {if (i == 2) F.v[i] += patom->c0 * Eext * fUnit_ext_mass;} // external electric field
#if _ATOM_INDUCED_DIPOLE == 1
						for (j = 0; j < 3; j++) {
							F.v[i] += patom->ind_mu.v[j] * patom->Eg.m[i][j] * fUnit_estat_mass;
						}
#endif
						pc->dyn->fc.v[3+i] += F.v[i];
					}
					V3PV3(patom->r0, F, torque)
					for (i = 0; i < 3; i++) pc->dyn->fc.v[i] += torque.v[i];
					if (bEext) {
						V -= patom->c0 * patom->r0.v[2] * Eext * U_EextFr / kT;
					}
				}
			}
			it.next_iterate();
		}
		//if (bEext) {
		//	V -= coeff * (cell.cell[nc].emp->mu.v[2] + cell.cell[nc].emp->q * cell.cell[nc].emp->r0.v[2]) * Eext * U_EextFr / kT; // external electric field
		//}
	}
}

#endif //_DISABLE_

extern PAIRWISE_DB< PAIRWISE<float> > TholeRadius_db;
// calculate the electrostatic interaction energy, Electric field and force on each atom
void EInteract(BASIC_CLUSTER *pc, EInteractVar &evar, InteractRes &res) {
	_EwaldSum_real_::EwaldSumRealVars1 *esv = &(evar.esv1); 
	MATOM *patom = NULL;
	VECTOR3 dr, r0, r1;
	int i, j, k, na, na1;
	double R, R2;
	
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;
	bool bound = false, bound14 = false;
	VECTOR3 E, F, E_Thole, F_Thole;
	double U, U_Thole;

	char aindx[2] = {0x00, 0x00};
	unsigned int key = 0;
	PAIRWISE<float> *TholeRadius = NULL;
	double a_Thole = 0;

	RELSHIP_ITERATOR<BASIC_CLUSTER, SPME_RELSHIP> rshipIt;

	res.reset();

	//for (na = 0; na < pc->nAtoms; na++) {
		//patom = pc->atom + na;
	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];
	#if _ATOM_INDUCED_DIPOLE == 0
		if (FABS(patom->c) < 0.0001) continue;
	#endif

		V3zero(patom->E)

		V32V3(patom->rg, r0)

		rshipIt.set(&pc->spmeRship); rshipIt.start_iterate();

		while (!rshipIt.eoi()) {
			imag_cluster = rshipIt.current(); if (imag_cluster == NULL) {rshipIt.next_iterate(); continue;}
			pc1 = imag_cluster->p;
			if (pc1->eNeutral) {rshipIt.next_iterate(); continue;}
			if (pc1 == pc && imag_cluster->local) {rshipIt.next_iterate(); continue;}

			//for (na1 = 0; na1 < pc1->nAtoms; na1++) {
				//patom1 = pc1->atom + na1;
			for (na1 = 0; na1 < pc1->fatom.n; na1++) {
				patom1 = pc1->fatom.m[na1];
				V32V3(patom1->rg, r1)
				VECT3(r1, r0, dr)
				scale_uv3(dr, dr, R2)
				R = sqrt(R2); 
				if (R2 < d2_14 && pc1->mIndx == pc->mIndx) {
					BOUND(patom, patom1, i, bound) // bound atom
					if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
					if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
					else bound14 = false;
				}
				else {bound = false; bound14 = false;}
				if (bound) continue;

				if (patom->mMP == 0) {
					esv->EwaldSumRealVars0::init_vars(dr, R, R2);
					MTP_Interact(*esv, patom->c, patom1->c, E, F, U, false, 0);
				}
				else if (patom->mMP == 1) {
					esv->init_vars(dr, R, R2);
					MTP_Interact(*esv, patom->c, patom->mu.v, patom1->c, patom1->mu.v, E, F, U, false, 0);
				}
	
				if (bound14) { // 1-4 interaction
					U *= A14_E;
					// electric field
					for (i = 0; i < 3; i++) {
						E.v[i] *= A14_E; F.v[i] *= A14_E;
					}
				}

				if (evar.bTholePolarize) {
					aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
					TholeRadius = ::TholeRadius_db.search(key);
					if (TholeRadius != NULL) {
						a_Thole = TholeRadius->pv;
						if (evar.tv.init_vars(a_Thole, R, R2, dr.v)) {
							::_Thole_polarize_::TholeCorr(evar.tv, a_Thole, patom->c, patom->mu.v, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v, F_Thole.v, U_Thole);
							if (bound14) {
								E_Thole.v[0] *= A14_E; E_Thole.v[1] *= A14_E; E_Thole.v[2] *= A14_E;
								F_Thole.v[0] *= A14_E; F_Thole.v[1] *= A14_E; F_Thole.v[2] *= A14_E;
								U_Thole *= A14_E;
							}
							V3plusV3(E, E_Thole, E) V3plusV3(F, F_Thole, F) U += U_Thole;
						}
					}
				}

				F.v[0] *= ::fUnit_estat_mass;
				F.v[1] *= ::fUnit_estat_mass;
				F.v[2] *= ::fUnit_estat_mass;
				V3plusV3(patom->E, E, patom->E)
				V3plusV3(patom->F, F, patom->F)
				res.U += U * ::eUnit_estat_kT;
				if (evar.bVirial) Accumulate_Virial_Diagonal(res, F.v, dr.v);
			}
					
			rshipIt.next_iterate();
		}
	}
}


// do not call this function, as scale_force_cluster rescale force on clusters within 10-cluster residule or more for macromolecule
/*
void mm_scale_force(MMOLECULE *mm, int nc1, int nc2) {
	CLUSTER *pc = NULL;
	double f = 0, t = 0;
	VECTOR3 force, torque;
	int nc = 0, i;
	bool scale = false;
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		scale_force_cluster(pc);
	}
}
*/

void scale_force_mmolecule(MMOLECULE &mm) {
	CLUSTER *pc = NULL;
	double f = 0, t = 0;
	VECTOR3 force, torque;
	int nc = 0, i;
	bool scale = false;
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		scale_force_cluster(pc);
	}
}

void get_child(CHAIN<BASIC_CLUSTER> *ch, CHAIN<BASIC_CLUSTER> **child) {
	if (*child != NULL) release_chain<BASIC_CLUSTER>(child, false);
	if (ch == NULL) return;
	CHAIN<BASIC_CLUSTER> *tch = ch, *tail = NULL;
	int ih;
	BASIC_CLUSTER *pc;
	while (tch != NULL) {
		pc = tch->p;
		if (pc->nHinges > 0) {
			for (ih = 0; ih < pc->nHinges; ih++) {
				if (ih == pc->n2Parent) continue;
				if (*child == NULL) {
					*child = new CHAIN<BASIC_CLUSTER>();
					(*child)->p = (BASIC_CLUSTER*)(pc->hinge[ih].plink);
					tail = *child;
				}
				else {
					tail->attach2tail( (BASIC_CLUSTER*)(pc->hinge[ih].plink)) ;
					tail = tail->next;
				}
			}
		}
		tch = tch->next;
	}
}

void scale_force_cluster(BASIC_CLUSTER *pc) {
	int i;
	double f2 = 0, t2 = 0, c;
	double *f = pc->dyn->fc.v, *t, *P;
	bool scale = pc->bShocked;
	if (!scale && ::scale_force) {
		for (i = 0; i < 3; i++) if (fabs(f[i]) > ::torque_max) {scale = true; break;}
		for (i = 3; i < 6; i++) if (fabs(f[i]) > ::force_max) {scale = true; break;}
	}
	if (!scale) return;

	//show_infor("scaled force on cluster", true);
	VECTOR3 force, torque, vt;
	MATRIX<6> *phai, *phai_T;
	CLUSTER *parent, *child, *tc;
	VECTOR<6> vp, dv, vnew, alpha;
	MATRIX<3> invI;

	V2wv(torque, force, pc->dyn->fc)

	f = pc->dyn->fc.v;
	// translation
	P = pc->bkm->P0.v + 3;
	project(P, force.v, vt.v);
	c = vt.v[0] * force.v[0] + vt.v[1] * force.v[1] + vt.v[2] * force.v[2];
	if (c < 0) {// force and P-projection are in same direction
		vt.v[0] = -vt.v[0]; vt.v[1] = -vt.v[1]; vt.v[2] = -vt.v[2];
	}
	else if (c == 0) goto _torque_;
	f2 = vt.abs();
	if (f2 / pc->M < 0.01) vt.stretch(0.01); // shift 0.4 Angstrom in a timestep at least
	memcpy(vnew.v + 3, vt.v, SIZE_V3); // new translation speed of cluster
	// force
	f2 = timestep_md * 0.6;// * 0.5;
	for (i = 0; i < 3; i++) {
		f[i+3] = vt.v[i] / f2;
		alpha.v[i+3] = f[i+3] / pc->M;
	}

_torque_:
	// this force related torque
	if (pc->parent != NULL) {
		VECT3((*(pc->Op)), (*(pc->vcm)), vt)
		V3PV3(vt, force, torque)
		memcpy(f, torque.v, SIZE_V3);
		tc = (CLUSTER*)pc;
		InverseMatrix<3>(tc->invNE->I, invI);
		MpV3(invI, torque, vt)
		memcpy(alpha.v, vt.v, SIZE_V3);
	}
	else if (/*pc->nHinges == 0 && */pc->nAtoms > 1) { // isolated cluster, SMOLECULE, not PMOLECULE
		memset(f, 0, SIZE_V3);
		// rotation speed change
		P = pc->bkm->P0.v;
		f = pc->dyn->fc.v;
		project(P, torque.v, vt.v);
		c = vt.v[0] * torque.v[0] + vt.v[1] * torque.v[1] + vt.v[2] * torque.v[2];
		if (c < 0) {// torque and P-projection are in same direction
			vt.v[0] = -vt.v[0]; vt.v[1] = -vt.v[1]; vt.v[2] = -vt.v[2];
		}
		memcpy(vnew.v, vt.v, SIZE_V3);
		// related torque
		t2 = timestep_md * 0.6;// * 0.5;
		for (i = 0; i < 3; i++) {
			vt.v[i] /= t2;
			f[i] += vt.v[i]; 
		}

		if (pc->nHinges > 0) { // this is the base of a macromolecule
			tc = (CLUSTER*)pc;
			InverseMatrix<3>(tc->invNE->I, invI);
			MpV3(invI, vt, alpha)
		}
	}

/*
	f = pc->dyn->fc.v;
	for (i = 0; i < 3; i++) {
		//if (f[i] > ::torque_max) f[i] = ::torque_max;
		//else if (f[i] < -::torque_max) f[i] = -::torque_max;

		if (f[i+3] > ::force_max) f[i+3] = ::force_max;
		else if (f[i+3] < -::force_max) f[i+3] = -::force_max;
	}
	*/

	
	CHAIN< BASIC_CLUSTER > *ch_this = NULL, *ch_child = NULL, *tch, *tthis;
	BASIC_CLUSTER *pbc;
	int nchild = 0, max_child = 10;
	VECTOR<6> *alpha0;
	if (pc->nHinges > 0) {
		// to have the cluster's speed changing to another direction, the rotation speed has to be changed also
		// so, we have to figure out the new rotation speed, relative to its parent, and even child
		// to setup the torque required !
		memcpy(((CLUSTER*)pc)->invNE->dalpha_H.v, alpha.v, SIZE_V6); // use dalpha_H to save the induced acceleration
		ch_this = new CHAIN<BASIC_CLUSTER>; ch_this->p = pc;
		while (ch_this != NULL) {
			nchild++; // nchild is the distance from pc
			if (nchild >= max_child) {
				if (ch_this != NULL) release_chain<BASIC_CLUSTER>(&ch_this, false);
				if (ch_child != NULL) release_chain<BASIC_CLUSTER>(&ch_child, false);
				break;
			}
			get_child(ch_this, &ch_child); 
			if (ch_child == NULL) break;
			tch = ch_child;
			while (tch != NULL) {
				pbc = tch->p; tc = (CLUSTER*)pbc;
				parent = (CLUSTER*)(tc->hinge[tc->n2Parent].plink);
				phai_T = &(tc->km.phai_T);
				phai = &(tc->km.phai);
				alpha0 = &(parent->invNE->dalpha_H);
				MpV6((*phai_T), (*alpha0), tc->invNE->dalpha_H)

				c = 1 - double(nchild) / max_child; // we compensate the induced acceleration gradually
				for (i = 0; i < 6; i++) {
					tc->invNE->dalpha_H.v[i] *= c;
					dv.v[i] = tc->invNE->dalpha_H.v[i] * (c - 1); // compensating acceleration, inverse direction
				}
				MpV6(tc->invNE->M, dv, vp) // required force
				
				f = tc->dyn->fc.v;
				V6plus(f, vp.v, f)

				tch = tch->next;
			}

			release_chain<BASIC_CLUSTER>(&ch_this, false);
			ch_this = ch_child; ch_child = NULL;
		}
	}
/*
	if (::scale_force) {
		f = pc->dyn->fc.v;
		for (i = 0; i < 3; i++) {
			if (f[i] > ::torque_max) f[i] = ::torque_max;
			else if (f[i] < -::torque_max) f[i] = -::torque_max;
		}
		for (i = 3; i < 6; i++) {
			if (f[i] > ::force_max) f[i] = ::force_max;
			else if (f[i] < -::force_max) f[i] = -::force_max;
		}
	}
*/
	/*
	int nc = 0, i;
	double f = 0, t = 0;
	bool scale = false;

	for (i = 0; i < 3; i++) {
		if (force.v[i] > force_max) pc->dyn->fc.v[i+3] = force_max;
		else if (force.v[i] < -force_max) pc->dyn->fc.v[i+3] = -force_max;

		if (torque.v[i] > torque_max) pc->dyn->fc.v[i] = torque_max;
		else if (torque.v[i] < -torque_max) pc->dyn->fc.v[i] = -torque_max;
	}
	*/
	/* scale in absolute force
	V3ABS2(force, f)
	if (f > force2_max) {
		f = force_max / sqrt(f);
		force.v[0] *= f; force.v[1] *= f; force.v[2] *= f;
		scale = true;
	}
	V3ABS2(torque, t)
	if (t > torque2_max) {
		t = torque_max / sqrt(t);
		torque.v[0] *= t; torque.v[1] *= t; torque.v[2] *= t;
		scale = true;
	}
	if (scale) {
		wv2V(torque, force, pc->dyn->fc)
	}
	*/
}

void SCALE_FORCE_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell) {
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		if (nc >= cell.acell.n) break;
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			scale_force_cluster(pc);
			it.next_iterate();
		}
	}
}

void EInteract_SINGLE_MM(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, EInteractVar &evar, InteractRes &res) {
	evar.esv1.bEwald = false; // direct interaction

	int nc;
	BASIC_CLUSTER *pc = NULL;
	InteractRes tres;
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = (BASIC_CLUSTER*)(mm.cluster + nc);
		EInteract(pc, evar, tres);
		res += tres;
	}
}

#if _DISABLE_  == 0
void FreeCellCMMInteraction_SingleMM(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, InteractVar &var, InteractRes &res) { // V in kT
	int nc, na, i, j;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	float coeff = sqrt(::eps_dc);

#error force is better to be calculated in EFIELD, while not accumulating here

	// calculate the electric field on each atom in these cells
	EFIELD_CMM_SINGLE_MM(cell, mm, n1, n2, V, OVERLAP);

	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = (BASIC_CLUSTER*)(mm.cluster + nc);
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			if (patom->c != 0) {
				for (i = 0; i < 3; i++) {
					F.v[i] = patom->c * patom->E.v[i] * fUnit_estat_mass;
					if (bEext) {if (i == 2) F.v[i] += patom->c0 * Eext * fUnit_ext_mass;} // external electric field
#if _ATOM_INDUCED_DIPOLE == 1
					for (j = 0; j < 3; j++) {
						F.v[i] += patom->ind_mu.v[j] * patom->Eg.m[i][j] * fUnit_estat_mass;
					}
#endif
					pc->dyn->fc.v[3+i] += F.v[i];
				}
				V3PV3(patom->r0, F, torque)
				for (i = 0; i < 3; i++) pc->dyn->fc.v[i] += torque.v[i];
				if (bEext) {
					V -= patom->c0 * patom->r0.v[2] * Eext * U_EextFr / kT;
				}
			}
		}

		//if (bEext) {
		//	V -= coeff * (pc->emp->mu.v[2] + pc->emp->q * pc->emp->r0.v[2]) * Eext * U_EextFr / kT; // external electric field
		//}
	}
}
#endif // _DISABLE_

// total solvent accessible surfacial area of the cluster -- for implicit solvent model
void Cluster_SASA(BASIC_CLUSTER *pc, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *bpc1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, IMPSOLV_RELSHIP<BASIC_CLUSTER> > ImpIt;
	
	double S0 = pc->psolv->S0, S1 = 0, S = 0, b = 0, bp = 0, cs = 0, fb, fbp;
	double radius0, radius1, rw = ::rSolvent, delta_r = 0;
	double s = rw + rw; // hard sphere atom, see Proc. Nat. Acad. Sci. USA 77(1980), 1737
	VECTOR3 f, torque, dr, r0, r1;
	double d = 0, d2, d_cut = 0, r0_w = 0, r1_w = 0, t;
	int i, n_neighbors = 0;

	pc->impsolv_fpar.release();

	CLUSTER_IMPSOLV_PARS *solv_par = NULL, *solv1_par = NULL;

	solv_par = pc->psolv; radius0 = solv_par->r0; r0_w = radius0 + rw;
	double r0_w_2 = r0_w * r0_w, r1_w_2 = 0;

	pc->E_ImplicitSolvation = 0;
	V3zero(f) V3zero(torque)
	if (cell->periodic_cell) {
		V3plusV3(pc->rc_solv, pc->dr_cmm, r0) // r0 = rc_solv + dr_cmm of pc
	}
	else {
		V32V3(pc->rc_solv, r0)
	}
	//ch = pc->impsolv_rship.ch;
	n_neighbors = pc->impsolv_rship.nrship();

	cs = 1;

	if (n_neighbors == 0) {pc->SASA = S0; return;}
	else {
		pc->impsolv_fpar.SetArray(n_neighbors);
		i = 0;
		cs = 1;
		//while (ch != NULL) {
		ImpIt.set(&(pc->impsolv_rship)); ImpIt.start_iterate();
		while (!ImpIt.eoi()) {
			imag_cluster = ImpIt.current(); if (imag_cluster == NULL) {ImpIt.next_iterate(); continue;}
			//bpc1 = ch->p->p;
			bpc1 = imag_cluster->p;
			if (bpc1 == NULL) {
				pc->impsolv_fpar.m[i].fb = 0; pc->impsolv_fpar.m[i].fbp = 0;
				pc->impsolv_fpar.m[i].fb1 = 0; pc->impsolv_fpar.m[i].fbp1 = 0;
				i++; 
				//ch = ch->next; 
				ImpIt.next_iterate();
				continue;
			}
			solv1_par = bpc1->psolv; radius1 = solv1_par->r0; S1 = solv1_par->S0;
			d_cut = radius0 + radius1 + rw + rw;
			delta_r = radius1 - radius0;
			r1_w = radius1 + rw;
			r1_w_2 = r1_w * r1_w;

			if (cell->periodic_cell) {
				if (imag_cluster->local && bpc1 == pc) { // same cluster on itself, ignore it
					pc->impsolv_fpar.m[i].fb = 0; pc->impsolv_fpar.m[i].fbp = 0;
					pc->impsolv_fpar.m[i].fb1 = 0; pc->impsolv_fpar.m[i].fbp1 = 0;
					//ch = ch->next; 
					ImpIt.next_iterate();
					i++; continue;
				}

				V3plusV3(bpc1->rc_solv, bpc1->dr_cmm, r1)  // r1 = r_solv + dr_cmm of bpc1
				if (!imag_cluster->local) {
					r1.v[0] -= imag_cluster->nx * cell->xd[0];
					r1.v[1] -= imag_cluster->ny * cell->xd[1];
					r1.v[2] -= imag_cluster->nz * cell->xd[2];
				}
			}
			else {
				V32V3(bpc1->rc_solv, r1)
			}
			
			VECT3(r0, r1, dr) // dr = r1 - r0
			V3ABS2(dr, d2) d = sqrt(d2); 
			if (d < 0.1) {
				pc->impsolv_fpar.m[i].fb = 0; pc->impsolv_fpar.m[i].fbp = 0;
				pc->impsolv_fpar.m[i].fb1 = 0; pc->impsolv_fpar.m[i].fbp1 = 0;
				//ch = ch->next; 
				ImpIt.next_iterate();
				i++; continue;
			}
			pc->impsolv_fpar.m[i].u.v[0] = dr.v[0] / d; 
			pc->impsolv_fpar.m[i].u.v[1] = dr.v[1] / d; 
			pc->impsolv_fpar.m[i].u.v[2] = dr.v[2] / d;

			if (d < d_cut) {
				b = PI * r0_w * (d_cut - d) * (1 + delta_r / d);
				if (b > 0) fb = -PI * r0_w * ((d_cut - d) * delta_r / d2 + (1 + delta_r / d));
				else {b = 0; fb = 0;}

				bp = PI * r0_w * (r0_w + r1_w - s - d) * (1 + (radius1 - s - radius0) / d);
				if (bp > 0) fbp = -PI * r0_w * ((r0_w + r1_w - s - d) * (radius1 - s - radius0) / d2 + (1 + (radius1 - s - radius0) / d));
				else {bp = 0; fbp = 0;}

				cs *= (1 - (b - bp) / S0);
				if (cs < 0) cs = 0;

				fb = (fb - fbp) / (S0 - (b - bp));
				pc->impsolv_fpar.m[i].fb = fb * solv_par->f_sigma;
				pc->impsolv_fpar.m[i].fbp = fbp * solv_par->f_sigma;

				// for bpc1, covered by pc
				b = PI * r1_w * (d_cut - d) * (1 - delta_r / d);
				if (b > 0) fb = PI * r1_w * ((d_cut - d) * delta_r / d2 - (1 - delta_r / d));
				else {b = 0; fb = 0;}

				bp = PI * r1_w * (r1_w + r0_w - s - d) * (1 + (radius0 - s - radius1) / d);
				if (bp > 0) fbp = -PI * r1_w * ((r1_w + r0_w - s - d) * (radius0 - s - radius1) / d2 + (1 + (radius0 - s - radius1) / d));
				else {bp = 0; fbp = 0;}

				fb = (fb - fbp) / (S1 - (b - bp));
				pc->impsolv_fpar.m[i].fb1 = fb * solv1_par->f_sigma;
				pc->impsolv_fpar.m[i].fbp1 = fbp * solv1_par->f_sigma;
			}
			else {
				pc->impsolv_fpar.m[i].fb = 0; pc->impsolv_fpar.m[i].fbp = 0;
				pc->impsolv_fpar.m[i].fb1 = 0; pc->impsolv_fpar.m[i].fbp1 = 0;
			}
			i++;
			//ch = ch->next;
			ImpIt.next_iterate();
		}
		S = S0 * cs; pc->SASA = S;
		pc->E_ImplicitSolvation = solv_par->sigma * S;
	}
}

void MM_ImpSolvSASA(CMM_CELL3D<BASIC_CLUSTER> *cell, MMOLECULE *mm, int nc1, int nc2) {
	BASIC_CLUSTER *bpc = NULL;
	for (int nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		bpc = (BASIC_CLUSTER*)(mm->cluster + nc);
		Cluster_SASA(bpc, cell);
	}
}

void MM_ImpSolvForce(MMOLECULE *mm, int nc1, int nc2) {
	BASIC_CLUSTER *bpc = NULL;
	for (int nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		bpc = (BASIC_CLUSTER*)(mm->cluster + nc);
		ClusterImpSolvForce(bpc);
	}
}

// Implicit Solvation Energy & Force, calculated with given impsolv_par
void ClusterImpSolvForce(BASIC_CLUSTER *pc) {
	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch = NULL, *ch1 = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *bpc1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, IMPSOLV_RELSHIP<BASIC_CLUSTER> > ImpIt;
	
	VECTOR3 f, torque;
	int i;
	IMPSOLV_F *fpar = pc->impsolv_fpar.m;
	double S = pc->SASA, S1 = 0, f0 = 0;

	V3zero(f) V3zero(torque)
	if (pc->impsolv_fpar.n == 0) return;
	else {
		i = 0; 
		//ch = pc->impsolv_rship.ch;
		ImpIt.set(&(pc->impsolv_rship)); ImpIt.start_iterate();
		V3zero(f)
		// we use the fact that: pc->impsolv_fpar was created with same amount of neighbors, and it was only created when the neighbors was changed
		for (i = 0; i < pc->impsolv_fpar.n; i++) {
			//bpc1 = ch->p->p;
			imag_cluster = ImpIt.current();
			bpc1 = imag_cluster->p;
			/*
			if (bpc1->mIndx == pc->mIndx) {
				if (bpc1->monIndx == pc->monIndx || pc->parent == (_BASIC_CLUSTER<MATOM>*)bpc1 || bpc1->parent == (_BASIC_CLUSTER<MATOM>*)pc) {
					ch = ch->next; continue;
				}


				//if (bpc1->monIndx == pc->monIndx || FABS(pc->monIndx - bpc1->monIndx) == 1) {
				//	ch = ch->next; continue;
				//}

				if (pc->parent == (_BASIC_CLUSTER<MATOM>*)bpc1 || bpc1->parent == (_BASIC_CLUSTER<MATOM>*)pc) {
					ch = ch->next; continue;
				}
			}
			*/

			S1 = bpc1->SASA;
			f0 = fpar[i].fb * S + fpar[i].fbp + fpar[i].fb1 * S1 + fpar[i].fbp1;

			f.v[0] -= f0 * fpar[i].u.v[0];
			f.v[1] -= f0 * fpar[i].u.v[1];
			f.v[2] -= f0 * fpar[i].u.v[2];

			//ch = ch->next;
			ImpIt.next_iterate();
		}
	}
	V3PV3(pc->rc0_solv, f, torque)
	for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += f.v[i];}
}

// Implicit Solvation Energy & Force
void ClusterImplicitSolvationForce(BASIC_CLUSTER *pc, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *bpc1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, IMPSOLV_RELSHIP<BASIC_CLUSTER> > ImpIt;
	
	double S0 = 0, S = 0, b = 0, cs = 0, fb;
	double radius0, radius1, rw = ::rSolvent, delta_r = 0;
	VECTOR3 f, torque, dr, r0, r1;
	double d = 0, d2, d_cut = 0, r0_w = 0, t;
	int i, n_neighbors = 0;

	IMPSOLV_F *fpar = NULL;

	CLUSTER_IMPSOLV_PARS *solv_par = NULL, *solv1_par = NULL;

	solv_par = pc->psolv; radius0 = solv_par->r0; r0_w = radius0 + rw;
	double r0_w_2 = (radius0 + rw) * (radius0 + rw), r1_w_2 = 0;

	pc->E_ImplicitSolvation = 0;
	V3zero(f) V3zero(torque)
	if (cell->periodic_cell) {
		V3plusV3(pc->rc_solv, pc->dr_cmm, r0) // r0 = rc_solv + dr_cmm of pc
	}
	else {
		V32V3(pc->rc_solv, r0)
	}
	//ch = pc->impsolv_rship.ch;
	//n_neighbors = number< CMM_IMAGINE<BASIC_CLUSTER> >(ch);
	n_neighbors = pc->impsolv_rship.nrship();

	S0 = 4 * PI * r0_w * r0_w;
	cs = 1;

	if (n_neighbors == 0) return;
	else {
		fpar = new IMPSOLV_F[n_neighbors];
		i = 0;
		cs = 1;
		//while (ch != NULL) {
		ImpIt.set(&(pc->impsolv_rship)); ImpIt.start_iterate();
		while (!ImpIt.eoi()) {
			imag_cluster = ImpIt.current();
			if (imag_cluster == NULL) {ImpIt.next_iterate(); continue;}
			//bpc1 = ch->p->p;
			bpc1 = imag_cluster->p;
			if (bpc1 == NULL) {
				fpar[i].fb = 0; fpar[i].fbp = 0;
				fpar[i].fb1 = 0; fpar[i].fbp1 = 0;
				//ch = ch->next; 
				ImpIt.next_iterate();
				i++; continue;
			}
			solv1_par = bpc1->psolv; radius1 = solv1_par->r0;
			d_cut = radius0 + radius1 + rw + rw;
			delta_r = radius1 - radius0;
			r1_w_2 = (radius1 + rw) * (radius1 + rw);

			if (cell->periodic_cell) {
				if (imag_cluster->local && bpc1 == pc) { // same cluster on itself, ignore it
					fpar[i].fb = 0; fpar[i].fbp = 0;
					fpar[i].fb1 = 0; fpar[i].fbp1 = 0;
					//ch = ch->next; 
					ImpIt.next_iterate();
					i++; continue;
				}

				V3plusV3(bpc1->rc_solv, bpc1->dr_cmm, r1)  // r1 = r_solv + dr_cmm of bpc1
				if (!imag_cluster->local) {
					r1.v[0] -= imag_cluster->nx * cell->xd[0];
					r1.v[1] -= imag_cluster->ny * cell->xd[1];
					r1.v[2] -= imag_cluster->nz * cell->xd[2];
				}
			}
			else {
				V32V3(bpc1->rc_solv, r1)
			}
			
			VECT3(r0, r1, dr) // dr = r1 - r0
			V3ABS2(dr, d2) d = sqrt(d2); 
			if (d < 0.1) {
				fpar[i].fb = 0; fpar[i].fbp = 0;
				fpar[i].fb1 = 0; fpar[i].fbp1 = 0;
				//ch = ch->next; 
				ImpIt.next_iterate();
				i++; continue;
			}
			fpar[i].u.v[0] = dr.v[0] / d; fpar[i].u.v[1] = dr.v[1] / d; fpar[i].u.v[2] = dr.v[2] / d;

			if (d < d_cut) {
				b = PI * r0_w * (d_cut - d) * (1 + delta_r / d);
				if (d >= FABS(delta_r)) fb = -PI * r0_w * ((d_cut - d) * delta_r / d2 + (1 + delta_r / d));
				else fb = -PI * r0_w * ((d_cut - d) / delta_r + 1 + sign(delta_r));
				cs *= (1 - b / S0);
				if (cs < 0) cs = 0;

				fb /= (S0 - b);
				fpar[i].fb = fb;
			}
			else {
				fpar[i].fb = 0;
			}
			i++;
			//ch = ch->next;
			ImpIt.next_iterate();
		}
		S = S0 * cs;
		pc->E_ImplicitSolvation = solv_par->sigma * S;

		V3zero(f)
		for (i = 0; i < n_neighbors; i++) {
			fpar[i].fb *= solv_par->f_sigma * S;
			f.v[0] -= fpar[i].fb * fpar[i].u.v[0];
			f.v[1] -= fpar[i].fb * fpar[i].u.v[1];
			f.v[2] -= fpar[i].fb * fpar[i].u.v[2];
		}
	}
	if (fpar != NULL) {delete[] fpar; fpar = NULL;}
	V3PV3(pc->rc0_solv, f, torque)
	for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += f.v[i];}
}

double ClusterImplicitSolvationEnergy(BASIC_CLUSTER *pc, BASIC_CLUSTER *pc1, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *bpc1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, IMPSOLV_RELSHIP<BASIC_CLUSTER> > ImpIt;
	
	double S0 = 0, S = 0, b = 0, cs = 0, fbi;
	double radius0, radius1, rw = ::rSolvent, delta_r = 0;
	VECTOR3 f, torque, dr, r0, r1;
	double d = 0, d2, d_cut = 0, r0_w = 0, t;
	int i, n_neighbors = 0;

	double eImpSolv = 0;

	CLUSTER_IMPSOLV_PARS *solv_par = NULL, *solv1_par = NULL;

	solv_par = pc->psolv; radius0 = solv_par->r0; r0_w = radius0 + rw;

	double r0_w_2 = (radius0 + rw) * (radius0 + rw), r1_w_2 = 0;

	V3zero(f) V3zero(torque)
	if (cell->periodic_cell) {
		V3plusV3(pc->rc_solv, pc->dr_cmm, r0) // r0 = rc_solv + dr_cmm of pc
	}
	else {
		V32V3(pc->rc_solv, r0)
	}
	//ch = pc->impsolv_rship.ch;
	//n_neighbors = number< CMM_IMAGINE<BASIC_CLUSTER> >(ch);
	n_neighbors = pc->impsolv_rship.nrship();

	S0 = 4 * PI * r0_w * r0_w;
	cs = 1;

	//if (n_neighbors == 0 || in_imagine_rship(ch, pc1) == NULL) return eImpSolv;
	if (n_neighbors == 0 || pc->impsolv_rship.get_rship(pc1) == NULL) return eImpSolv;
	else {
		i = 0;
		cs = 1;
		//while (ch != NULL) {
		ImpIt.set(&(pc->impsolv_rship)); ImpIt.start_iterate();
		while (!ImpIt.eoi()) {
			//bpc1 = ch->p->p;
			imag_cluster = ImpIt.current(); if (imag_cluster == NULL) {ImpIt.next_iterate(); continue;}
			if (bpc1 == NULL) {
				//ch = ch->next; 
				ImpIt.next_iterate();
				i++; continue;
			}
			solv1_par = bpc1->psolv; radius1 = solv1_par->r0;
			d_cut = radius0 + radius1 + rw + rw;
			delta_r = radius1 - radius0;
			r1_w_2 = (radius1 + rw) * (radius1 + rw);

			if (cell->periodic_cell) {
				if (imag_cluster->local && bpc1 == pc) { // same cluster on itself, ignore it
					//ch = ch->next; 
					ImpIt.next_iterate();
					i++; continue;
				}

				V3plusV3(bpc1->rc_solv, bpc1->dr_cmm, r1)  // r1 = r_solv + dr_cmm of bpc1
				if (!imag_cluster->local) {
					r1.v[0] += imag_cluster->nx * cell->xd[0];
					r1.v[1] += imag_cluster->ny * cell->xd[1];
					r1.v[2] += imag_cluster->nz * cell->xd[2];
				}
			}
			else {
				V32V3(bpc1->rc_solv, r1)
			}
			
			VECT3(r0, r1, dr) // dr = r1 - r0
			V3ABS2(dr, d2) d = sqrt(d2); 
			if (d < 0.1) {
				//ch = ch->next; 
				ImpIt.next_iterate();
				i++; continue;
			}

			if (d < d_cut) {
				b = PI * r0_w * (d_cut - d) * (1 + delta_r / d);
				cs *= (1 - b / S0);
				if (cs < 0) cs = 0;
			}
			i++;
			//ch = ch->next;
			ImpIt.next_iterate();
		}
		S = S0 * cs;
		eImpSolv = solv_par->sigma * S;
		pc->E_ImplicitSolvation = solv_par->sigma * S;
	}
	return eImpSolv;
}

double ClusterImplicitSolvationEnergy_Direct(BASIC_CLUSTER *pc, BASIC_CLUSTER *pc1, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	//CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	BASIC_CLUSTER *bpc1 = NULL;
	
	double S0 = 0, S = 0, b = 0, cs = 0, fbi;
	double radius0, radius1, rw = ::rSolvent, delta_r = 0;
	VECTOR3 f, torque, dr, r0, r1;
	double d = 0, d2, d_cut = 0, r0_w = 0, t;
	int i, n_neighbors = 0;

	double eImpSolv = 0;

	CLUSTER_IMPSOLV_PARS *solv_par = NULL, *solv1_par = NULL;

	solv_par = pc->psolv; radius0 = solv_par->r0; r0_w = radius0 + rw;

	double r0_w_2 = (radius0 + rw) * (radius0 + rw), r1_w_2 = 0;

	V3zero(f) V3zero(torque)
	if (cell->periodic_cell) {
		V3plusV3(pc->rc_solv, pc->dr_cmm, r0) // r0 = rc_solv + dr_cmm of pc
	}
	else {
		V32V3(pc->rc_solv, r0)
	}
	//ch = pc->impsolv_rship.ch;
	//n_neighbors = number< CMM_IMAGINE<BASIC_CLUSTER> >(ch);
	n_neighbors = pc->impsolv_rship.nrship();

	S0 = 4 * PI * r0_w * r0_w;
	cs = 1;

	bpc1 = pc1;

	solv1_par = bpc1->psolv; radius1 = solv1_par->r0;
	d_cut = radius0 + radius1 + rw + rw;
	delta_r = radius1 - radius0;
	r1_w_2 = (radius1 + rw) * (radius1 + rw);

	V3plusV3(bpc1->rc_solv, bpc1->dr_cmm, r1)  // r1 = r_solv + dr_cmm of bpc1
	VECT3(r0, r1, dr) // dr = r1 - r0
	if (cell->periodic_cell) {
		if (dr.v[0] > cell->hw[0]) dr.v[0] -= cell->xd[0];
		else if (dr.v[0] < -(cell->hw[0])) dr.v[0] += cell->xd[0] ;
		if (dr.v[1] > cell->hw[1]) dr.v[1] -= cell->xd[1];
		else if (dr.v[1] < -(cell->hw[1])) dr.v[1] += cell->xd[1];
		if (dr.v[2] > cell->hw[2]) dr.v[2] -= cell->xd[2];
		else if (dr.v[2] < -(cell->hw[2])) dr.v[2] += cell->xd[2];
	}

	V3ABS2(dr, d2)d = sqrt(d2); 

	if (d < d_cut) {
		b = PI * r0_w * (d_cut - d) * (1 + delta_r / d);
		cs = (1 - b / S0);
		if (cs < 0) cs = 0;
	}
	else cs = 1;
	S = S0 * cs;
	eImpSolv = solv_par->sigma * S;
	pc->E_ImplicitSolvation = solv_par->sigma * S;
	return eImpSolv;
}

// solvent effect between cluster in the same macro-molecules
void CMM_ImplicitSolvationForce(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell) {
	BASIC_CLUSTER *basic_pc = NULL;
	CLUSTER *pc = NULL;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	int i, j, k;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;

	for (i = n1Cell; i <= n2Cell; i++) {
		if (i >= cell.bcell.nx) continue;
		for (j = 0; j < cell.bcell.ny; j++) {for (k = 0; k < cell.bcell.nz; k++) {
			pcell = cell.bcell.m[i].m[j].m + k; it.set(pcell); it.start_iterate();
			while (!it.eoi()) {
				basic_pc = it.current(); if (basic_pc == NULL) {it.next_iterate(); continue;}
				ClusterImplicitSolvationForce(basic_pc, &cell);
				it.next_iterate();
			}
		}}
	}
}

void MM_ImplicitSolvationForce(CMM_CELL3D<BASIC_CLUSTER> *cell, MMOLECULE *mm, int nc1, int nc2) {
	BASIC_CLUSTER *bpc = NULL;
	for (int nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		bpc = (BASIC_CLUSTER*)(mm->cluster + nc);
		ClusterImplicitSolvationForce(bpc, cell);
	}
}



void Brownian_force(MMOLECULE *mm, int nc1, int nc2, VECTOR<6> &f0) { // f0 -- total Brownian force at mass center
	V6zero(f0)
	int nc = 0, i;
	VECTOR<6> ft, t0;
	VECTOR3 f, t, dr;
	CLUSTER *pc = NULL;

	double v2 = 0, v = 0;
	double *u = NULL;
	bool bForceMassCenter = (::MM_SpeedVerlet_method == 0 ? true : false);
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) return;// f0;
		pc = mm->cluster + nc;
		//TIME(pc->bkm.V0, (-ita_Brownian), pc->dyn->f_Brown, i)
		u = pc->bkm->V0.v; v2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]; v = sqrt(v2);
		
		for (i = 0; i < 6; i++) {
			//pc->dyn->f_Brown.v[i] = -pc->ita_Brownian * pc->bkm->V0.v[i] - pc->cd_rot * pc->bkm->V0.v[i] * v;
			ft.v[i] = -pc->ita_Brownian * pc->bkm->V0.v[i] - pc->cd_rot * pc->bkm->V0.v[i] * v;
		}
		MpV6(pc->invNE->M, ft, t0)
		memcpy(pc->dyn->f_Brown.v, t0.v, SIZE_V6);
		//memcpy(pc->dyn->f_Brown.v, ft.v, SIZE_V3);
		/*
		u = pc->bkm->V0.v + 3; v2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]; v = sqrt(v2);
		for (i = 0; i < 3; i++) {
			//pc->dyn->f_Brown.v[i] = -pc->ita_Brownian * pc->bkm->V0.v[i] - pc->cd_trans * pc->bkm->V0.v[i] * v;
			f.v[i] = -pc->ita_Brownian * pc->bkm->V0.v[i + 3] - pc->cd_trans * pc->bkm->V0.v[i + 3] * v;
			f.v[i] *= pc->M;
		}
		VECT3((*(pc->Op)), pc->cm, dr)
		V3PV3(dr, f, t)
		for (i = 0; i < 3; i++) {pc->dyn->f_Brown.v[i] = t.v[i] + t.v[i]; pc->dyn->f_Brown.v[i+3] = f.v[i]; }
		*/
	
		if (bForceMassCenter) {
			MpV6(pc->invNE->phai_cm2me, pc->dyn->f_Brown, ft)
			V6plusV6(f0, ft, f0)
		}
	}
	//return f0;
}

void SM_Brownian_force(SMOLECULE *sm) {
	V6zero(sm->c->dyn->f_Brown)
	int nc = 0, i;
	VECTOR<6> ft;
	BASIC_CLUSTER *bpc = sm->c;

	double v2 = 0, v = 0;
	double *u = NULL;
	u = bpc->bkm->V0.v; v2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]; v = sqrt(v2);
	for (i = 0; i < 3; i++) {
		bpc->dyn->f_Brown.v[i] = -bpc->ita_Brownian * bpc->bkm->V0.v[i] - bpc->cd_rot * bpc->bkm->V0.v[i] * v;
	}
	u = bpc->bkm->V0.v + 3; v2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]; v = sqrt(v2);
	for (i = 3; i < 6; i++) {
		bpc->dyn->f_Brown.v[i] = -bpc->ita_Brownian * bpc->bkm->V0.v[i] - bpc->cd_trans * bpc->bkm->V0.v[i] * v;
	}
}

void PM_Brownian_force(PMOLECULE *pm) {
	V6zero(pm->c->dyn->f_Brown)
	int nc = 0, i;
	VECTOR<6> ft;
	BASIC_CLUSTER *bpc = pm->c;

	double v2 = 0, v = 0;
	double *u = pm->v.v; 
	if (bpc->cd_trans != 0) {
		v2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]; v = sqrt(v2);
	}
	for (i = 0; i < 3; i++) {
		bpc->dyn->f_Brown.v[i+3] = -bpc->ita_Brownian * u[i];
		if (bpc->cd_trans != 0) bpc->dyn->f_Brown.v[i+3]-= bpc->cd_trans * u[i] * v;
	}
}


// Dihedral torsion relating to the hinge, calculating with force on each atom of the dihedral angle
void TorsionInteract(MMOLECULE*mm, int nc1, int nc2, TorsionVar &var, InteractRes& ires) { // same as MMOLECULE::calTorsion(double &E) above
	int nc = 0, na, na1, i, ih;
	double Etor = 0, E = 0, alpha, f0, torque;
	CLUSTER *pc = NULL, *pc1 = NULL;
	MATOM *patom1, *patom2, *patom3, *patom4;
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge, *hinge1;
	TORSION_PARS *tpar = NULL, *pTorPar = NULL;

	VECTOR3 f1, f2, f3, f4;
	double virial[3] = {0, 0, 0};
	double v1[3] = {0, 0, 0};
	double sin_phi;

	VECTOR3 t;

	char msg[256] = "\0";

	bool bAtomicForce = false;

	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		if (pc == mm->base) continue;
/*
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
				tpar = dh_db.get_torsion(patom1->aindx, patom2->aindx, patom3->aindx, patom4->aindx);
				if (tpar == NULL || (tpar->v0 == 0 && tpar->next == NULL)) continue;
				*/
		for (ih = 0; ih < pc->dh.n; ih++) {
			patom1 = pc->dh.m[ih].dh[0]; patom2 = pc->dh.m[ih].dh[1];
			patom3 = pc->dh.m[ih].dh[2]; patom4 = pc->dh.m[ih].dh[3];
			tpar = pc->dh.m[ih].tpar;

			var.dh.r1 = &(patom1->r); var.dh.r2 = &(patom2->r); var.dh.r3 = &(patom3->r); var.dh.r4 = &(patom4->r);
			var.dh.init(); var.dh.derive_B();
			f0 = 0; Etor = 0; // f0 = -d_U/d_cos(phi) = 1 / sin(phi) * dU / d_phi = - 1/sin(phi) * SUM[ n * v0 * sin(n * phi - phase) ]
			pTorPar = tpar; torque = 0;
			while (pTorPar != NULL) {
				if (pTorPar->v0 == 0) {pTorPar = pTorPar->next; continue;}
				alpha = pTorPar->pn * var.dh.phi - pTorPar->phase;
				Etor += pTorPar->v0 * (1 + cos(alpha));
				f0 += tpar->pn * tpar->v0 * sin(alpha);
				pTorPar = pTorPar->next;
			}
			E += Etor * eUnit_LJ_kT; //kJ/mol ==> kT
			torque = -f0 * fUnit_LJ_mass;
			sin_phi = sin(var.dh.phi);
			if (FABS(sin_phi) < 1e-5) f0 = 0;
			else f0 *= -1.0 / sin_phi * fUnit_LJ_mass; // kJ/mol ==> atom inertial moment
			// F = -d_U / d_cos(phi) * d_cos(phi) / d_r = f0 * d_B / d_r = f0 * dh.dB
			for (i = 0; i < 3; i++) {
				f1.v[i] = f0 * var.dh.dB1.v[i];
				f2.v[i] = f0 * var.dh.dB2.v[i];
				f3.v[i] = f0 * var.dh.dB3.v[i];
				f4.v[i] = f0 * var.dh.dB4.v[i];

				if (bAtomicForce) {
					patom1->F.v[i] += f1.v[i];
					patom2->F.v[i] += f2.v[i];
					patom3->F.v[i] += f3.v[i];
					patom4->F.v[i] += f4.v[i];
				}
			}
			// virial calculation of dh.cal_virial() is not correct!!
			if (var.bVirial) {
				// virial calculation of dh.cal_virial() is not correct!!
				//var.dh.cal_virial();
				//virial[0] += f0 * var.dh.virial.v[0];
				//virial[1] += f0 * var.dh.virial.v[1];
				//virial[2] += f0 * var.dh.virial.v[2];
			}


			// calculate virial in explicit way
			if (var.bVirial) {
				v1[0] -= f1.v[0] * patom1->r.v[0] + f2.v[0] * patom2->r.v[0] + f3.v[0] * patom3->r.v[0] + f4.v[0] * patom4->r.v[0];
				v1[1] -= f1.v[1] * patom1->r.v[1] + f2.v[1] * patom2->r.v[1] + f3.v[1] * patom3->r.v[1] + f4.v[1] * patom4->r.v[1];
				v1[2] -= f1.v[2] * patom1->r.v[2] + f2.v[2] * patom2->r.v[2] + f3.v[2] * patom3->r.v[2] + f4.v[2] * patom4->r.v[2];
			}

			//t.v[0] = f1.v[0] + f2.v[0] + f3.v[0] + f4.v[0];
			//t.v[1] = f1.v[1] + f2.v[1] + f3.v[1] + f4.v[1];
			//t.v[2] = f1.v[2] + f2.v[2] + f3.v[2] + f4.v[2];

			if (!bAtomicForce) { // calculate torque directly
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
	}
	ires.U = E;
	if (var.bVirial) {
		//ires.STxx = virial[0];
		//ires.STyy = virial[1];
		//ires.STzz = virial[2];
		//ires.STtr = virial[0] + virial[1] + virial[2];

		ires.STxx = v1[0];
		ires.STyy = v1[1];
		ires.STzz = v1[2];
		ires.STtr = v1[0] + v1[1] + v1[2];
	}
}

bool HingeConstraintForce(MMOLECULE &mm, HingeConstraintVar &var, InteractRes& ires) {
	int ic, ia;
	CLUSTER *pc;
	HINGE_CONSTRAINT_FORCE *hc = &(var.hc);
	char msg[256] = "\0";

	double *fv;
	VECTOR3 *f00, *f01, *f02, *f10, *f11, *f12;

	SVECTOR<double, 3> ft, tq;
	ires.reset();

	for (ic = 0; ic < mm.nCluster; ic++) {
		pc = mm.cluster + ic;
		if (pc == mm.base) continue; // base cluster has no parent hinge
		hc->init(&(pc->hinge_constraint->patom_00->r0), &(pc->hinge_constraint->patom_01->r0), &(pc->hinge_constraint->patom_02->r0),
			&(pc->hinge_constraint->patom_10->r0), &(pc->hinge_constraint->patom_11->r0), &(pc->hinge_constraint->patom_12->r0));
		hc->hinge_force(pc->fphinge);
		if (!hc->cal_constraint_force()) {
			sprintf(msg, "failure to construct the constraint force on the parent hinge of cluster %d", ic); show_log(msg, true);
			return false;
		}

		// force on atom 00
		fv = pc->hinge_constraint->fc0[0].v; f00 = pc->hinge_constraint->fc0;
		fv[0] = hc->f_10_00 * hc->u_10_00.v[0] + hc->f_11_00 * hc->u_11_00.v[0] + hc->f_12_00 * hc->u_12_00.v[0];
		fv[1] = hc->f_10_00 * hc->u_10_00.v[1] + hc->f_11_00 * hc->u_11_00.v[1] + hc->f_12_00 * hc->u_12_00.v[1];
		fv[2] = hc->f_10_00 * hc->u_10_00.v[2] + hc->f_11_00 * hc->u_11_00.v[2] + hc->f_12_00 * hc->u_12_00.v[2];

		// force on atom 01
		fv = pc->hinge_constraint->fc0[1].v; f01 = pc->hinge_constraint->fc0 + 1;
		fv[0] = hc->f_10_01 * hc->u_10_01.v[0];
		fv[1] = hc->f_10_01 * hc->u_10_01.v[1];
		fv[2] = hc->f_10_01 * hc->u_10_01.v[2];

		// force on atom 02
		fv = pc->hinge_constraint->fc0[2].v; f02 = pc->hinge_constraint->fc0 + 2;
		fv[0] = hc->f_10_02 * hc->u_10_02.v[0];
		fv[1] = hc->f_10_02 * hc->u_10_02.v[1];
		fv[2] = hc->f_10_02 * hc->u_10_02.v[2];

		// force on atom 10
		fv = pc->hinge_constraint->fc1[0].v; f10 = pc->hinge_constraint->fc1;
		fv[0] = -hc->f_10_00 * hc->u_10_00.v[0] - hc->f_10_01 * hc->u_10_01.v[0] - hc->f_10_02 * hc->u_10_02.v[0];
		fv[1] = -hc->f_10_00 * hc->u_10_00.v[1] - hc->f_10_01 * hc->u_10_01.v[1] - hc->f_10_02 * hc->u_10_02.v[1];
		fv[2] = -hc->f_10_00 * hc->u_10_00.v[2] - hc->f_10_01 * hc->u_10_01.v[2] - hc->f_10_02 * hc->u_10_02.v[2];

		// force on atom 11
		fv = pc->hinge_constraint->fc1[1].v; f11 = pc->hinge_constraint->fc1 + 1;
		fv[0] = -hc->f_11_00 * hc->u_11_00.v[0];
		fv[1] = -hc->f_11_00 * hc->u_11_00.v[1];
		fv[2] = -hc->f_11_00 * hc->u_11_00.v[2];

		// force on atom 12
		fv = pc->hinge_constraint->fc1[2].v; f12 = pc->hinge_constraint->fc1 + 2;
		fv[0] = -hc->f_12_00 * hc->u_12_00.v[0];
		fv[1] = -hc->f_12_00 * hc->u_12_00.v[1];
		fv[2] = -hc->f_12_00 * hc->u_12_00.v[2];

		if (var.bVirial) {
			// important: r10 is zero, and the total force from cluster #0 to cluster #1 is same as that from cluster #1 to cluster @0
			// So, total constraint forces on cluster #0 and cluster #1 is zero!
			ires.STxx += f00->v[0] * hc->r00.v[0] + f01->v[0] * hc->r01.v[0] + f02->v[0] * hc->r02.v[0] + f11->v[0] * hc->r11.v[0] + f12->v[0] * hc->r12.v[0];
			ires.STyy += f00->v[1] * hc->r00.v[1] + f01->v[1] * hc->r01.v[1] + f02->v[1] * hc->r02.v[1] + f11->v[1] * hc->r11.v[1] + f12->v[1] * hc->r12.v[1];
			ires.STzz += f00->v[2] * hc->r00.v[2] + f01->v[2] * hc->r01.v[2] + f02->v[2] * hc->r02.v[2] + f11->v[2] * hc->r11.v[2] + f12->v[2] * hc->r12.v[2];
		}

		// check whether the result is right
		/*
		ft.v[0] = f00->v[0] + f01->v[0] + f02->v[0];
		ft.v[1] = f00->v[1] + f01->v[1] + f02->v[1];
		ft.v[2] = f00->v[2] + f01->v[2] + f02->v[2];
		sprintf(msg, "hinge force : [%f, %f, %f]; ", ft.v[0] - hc->f.v[0], ft.v[1] - hc->f.v[1], ft.v[2] - hc->f.v[2]); show_log(msg, false);

		ft.v[0] = f00->v[0] + f01->v[0] + f02->v[0] + f10->v[0] + f11->v[0] + f12->v[0];
		ft.v[1] = f00->v[1] + f01->v[1] + f02->v[1] + f10->v[1] + f11->v[1] + f12->v[1];
		ft.v[2] = f00->v[2] + f01->v[2] + f02->v[2] + f10->v[2] + f11->v[2] + f12->v[2];
		sprintf(msg, "net force : [%f, %f, %f]", ft.v[0], ft.v[1], ft.v[2]); show_log(msg, true);

		V3zero(tq)
		V3PV3(hc->r00, (*f00), ft) V3plusV3(ft, tq, tq);
		V3PV3(hc->r01, (*f01), ft) V3plusV3(ft, tq, tq);
		V3PV3(hc->r02, (*f02), ft) V3plusV3(ft, tq, tq);
		sprintf(msg, "torq : [%f, %f, %f]", tq.v[0] - hc->torq.v[0], tq.v[1] - hc->torq.v[1], tq.v[2] - hc->torq.v[2]); show_log(msg, true);
		sprintf(msg, "torq : [%f, %f, %f]", hc->torq.v[0], hc->torq.v[1], hc->torq.v[2]); show_log(msg, true);
		*/
	}
	//show_log("", true);
	if (var.bVirial) ires.STtr = ires.STxx + ires.STyy + ires.STzz;
	return true;
}





// External Field on the molecules

void ExternalEField(BASIC_CLUSTER *bpc, float &Eext, double &Uext) {
	MATOM *patom = NULL;
	int ia;
	Uext = 0;
	//if (::bEext) {
		for (ia = 0; ia < bpc->nAtoms; ia++) {
			patom = bpc->atom + ia;
			patom->E.v[2] += Eext * fUnit_ext_estat; // need to convert to unit e/A^2 
			//if (!var.bEfieldOnly) {
				patom->F.v[2] += patom->c0 * Eext * fUnit_ext_mass;
				Uext -= patom->c0 * patom->r0.v[2] * Eext * U_EextFr / kT;  // external electric field
			//}
		}
	//}
}

void CLUSTER_ExternalEField(CLUSTER *cluster, int ncs, float &Eext, double &Uext) {
	Uext = 0;
	double Ut = 0;
	int ic;
	for (ic = 0; ic < ncs; ic++) {
		ExternalEField(cluster + ic, Eext, Ut);
		Uext += Ut;
	}
}

void MM_ExternalEField(MMOLECULE *mm, int nM, float &Eext, double &Uext) {
	Uext = 0;
	double Ut = 0;
	int im, ic;
	for (im = 0; im < nM; im++) {
		for (ic = 0; ic < mm[im].nCluster; ic++) {
			ExternalEField(mm[im].cluster + ic, Eext, Ut);
			Uext += Ut;
		}
	}
}

void VMM_ExternalEField(MMOLECULE **mm, int nM, float &Eext, double &Uext) {
	Uext = 0;
	double Ut = 0;
	int im, ic;
	for (im = 0; im < nM; im++) {
		for (ic = 0; ic < mm[im]->nCluster; ic++) {
			ExternalEField(mm[im]->cluster + ic, Eext, Ut);
			Uext += Ut;
		}
	}
}


void SM_ExternalEField(SMOLECULE *sm, int nM, float &Eext, double &Uext) {
	Uext = 0;
	double Ut = 0;
	int im;
	for (im = 0; im < nM; im++) {
		ExternalEField(sm[im].c, Eext, Ut);
		Uext += Ut;
	}
}

void VSM_ExternalEField(SMOLECULE **sm, int nM, float &Eext, double &Uext) {
	Uext = 0;
	double Ut = 0;
	int im;
	for (im = 0; im < nM; im++) {
		ExternalEField(sm[im]->c, Eext, Ut);
		Uext += Ut;
	}
}

void PM_ExternalEField(PMOLECULE *pm, int nM, float &Eext, double &Uext) {
	Uext = 0;
	double Ut = 0;
	int im;
	for (im = 0; im < nM; im++) {
		ExternalEField(pm[im].c, Eext, Ut);
		Uext += Ut;
	}
}

void VPM_ExternalEField(PMOLECULE **pm, int nM, float &Eext, double &Uext) {
	Uext = 0;
	double Ut = 0;
	int im;
	for (im = 0; im < nM; im++) {
		ExternalEField(pm[im]->c, Eext, Ut);
		Uext += Ut;
	}
}



