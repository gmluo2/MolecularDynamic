#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "ranlib.h"

#include "def.h"
#include "show.h"
#include "ranlib.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"
#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "var.h"

#include "Interaction.h"
#include "Interact1.h"
#include "Interact2.h"
#include "Interact3.h"
#include "NEIMO.h"
#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"
#include "MD_POLYMER.h"
#include "MD_SM.h"
#include "md_cell.h"

#include "vmd_cell.h"

#include "complex.h"
#include <fftw3.h>

#include "EwaldSum.h"
#include "spme.h"
#include "spme_2d.h"

extern BSplineFunc<2> bsp; // with 2nd order derivation

namespace _spme_2d_ {
#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

void cal_Q_spme(VSPME<_VQatom>& vspme, int nThreads) { // Q value is initialized
	//vspme.Q = 0; // important
	FLEX_ASYM_CMATRIX<double> zeroM(0.0);
	RECP_SPME_RES res[MAX_THREADS];
	int iThread;
	for (iThread = 0; iThread < nThreads; iThread++) {
		res[iThread].mMP = 0; // no multipole
		res[iThread].Q = vspme.Qbuf + iThread;
		*(res[iThread].Q) = zeroM; //(FLEX_ASYM_CMATRIX<double>)0;
	}

	MacroMoleculeOperate1_mv<VSPME<_VQatom>, RECP_SPME_RES >((void*)(&VSPME_Q_accum_q), vspme, 0, vspme.va.n, res, nThreads);
	int i, j, k;
	if (nThreads == 1) {
		for (i = 0; i < vspme.Q.nx; i++) {for (j = 0; j < vspme.Q.ny; j++) {
			memcpy(vspme.Q.m[i].m[j].m, res[0].Q->m[i].m[j].m, vspme.Q.nz * SIZE_DOUBLE);
		}}
	}
	else {
		for (i = 0; i < vspme.Q.nx; i++) {for (j = 0; j < vspme.Q.ny; j++) {for (k = 0; k < vspme.Q.nz; k++) {
			vspme.Q.m[i].m[j].m[k] = res[0].Q->m[i].m[j].m[k];
			for (iThread = 1; iThread < nThreads; iThread += 1) vspme.Q.m[i].m[j].m[k] += res[iThread].Q->m[i].m[j].m[k];
		}}}
	}
}

void cal_Q_spme(VSPME<_VMUatom>& vspme, int nThreads) { // Q value is initialized
	//vspme.Q = 0; // important
	FLEX_ASYM_CMATRIX<double> zeroM(0.0);
	RECP_SPME_RES res[MAX_THREADS];
	
#if _CHECK_TIME_STAMP_
	TIME mt;
	mt.start();
	show_log("  (", false);
#endif

	int iThread;
	for (iThread = 0; iThread < nThreads; iThread+=1) {
		res[iThread].mMP = 1; // dipole
		res[iThread].Q = vspme.Qbuf + iThread; 
		*(res[iThread].Q) = zeroM; //(FLEX_ASYM_CMATRIX<double>)0;
		res[iThread].f_mu_accumulate = (void*)(&MU_accumulate);
	}
#if _CHECK_TIME_STAMP_
	mt.elapse(true, 1); show_log(", ", false);
#endif

	MacroMoleculeOperate1_mv<VSPME<_VMUatom>, RECP_SPME_RES >((void*)(&VSPME_Q_accum_mu), vspme, 0, vspme.va.n, res, nThreads);
	int i, j, k;
	if (nThreads == 1) {
		for (i = 0; i < vspme.Q.nx; i++) {for (j = 0; j < vspme.Q.ny; j++) {
			memcpy(vspme.Q.m[i].m[j].m, res[0].Q->m[i].m[j].m, vspme.Q.nz * SIZE_DOUBLE);
		}}
	}
	else {
		for (i = 0; i < vspme.Q.nx; i++) {for (j = 0; j < vspme.Q.ny; j++) {for (k = 0; k < vspme.Q.nz; k++) {
			vspme.Q.m[i].m[j].m[k] = res[0].Q->m[i].m[j].m[k];
			for (iThread = 1; iThread < nThreads; iThread += 1) vspme.Q.m[i].m[j].m[k] += res[iThread].Q->m[i].m[j].m[k];
		}}}
	}

#if _CHECK_TIME_STAMP_
	mt.elapse(true, 2); show_log("), ", false);
#endif
}





void cal_TQ_spme(VSPME<_VQatom>& vspme, int nThreads) { // Q value is initialized
	//vspme.TQ = 0; // important
	FLEX_ASYM_CMATRIX<double> zeroM(0.0);
	RECP_SPME_RES res[MAX_THREADS];
	int iThread;
	for (iThread = 0; iThread < nThreads; iThread++) {
		res[iThread].mMP = 0; // no multipole
		res[iThread].Q = vspme.Qbuf + iThread;
		*(res[iThread].Q) = zeroM; //(FLEX_ASYM_CMATRIX<double>)0;
	}

	MacroMoleculeOperate1_mv<VSPME<_VQatom>, RECP_SPME_RES >((void*)(&VSPME_TQ_accum_q), vspme, 0, vspme.va.n, res, nThreads);
	int i, j, k;
	if (nThreads == 1) {
		for (i = 0; i < vspme.TQ.nx; i++) {for (j = 0; j < vspme.TQ.ny; j++) {
			memcpy(vspme.TQ.m[i].m[j].m, res[0].Q->m[i].m[j].m, vspme.TQ.nz * SIZE_DOUBLE);
		}}
	}
	else {
		for (i = 0; i < vspme.TQ.nx; i++) {for (j = 0; j < vspme.TQ.ny; j++) {for (k = 0; k < vspme.TQ.nz; k++) {
			vspme.TQ.m[i].m[j].m[k] = res[0].Q->m[i].m[j].m[k];
			for (iThread = 1; iThread < nThreads; iThread += 1) vspme.TQ.m[i].m[j].m[k] += res[iThread].Q->m[i].m[j].m[k];
		}}}
	}
}

void cal_TQ_spme(VSPME<_VMUatom>& vspme, int nThreads) { // Q value is initialized
	//vspme.TQ = 0; // important
	FLEX_ASYM_CMATRIX<double> zeroM(0.0);
	RECP_SPME_RES res[MAX_THREADS];
	
#if _CHECK_TIME_STAMP_
	TIME mt;
	mt.start();
	show_log("  (", false);
#endif

	int iThread;
	for (iThread = 0; iThread < nThreads; iThread+=1) {
		res[iThread].mMP = 1; // dipole
		res[iThread].Q = vspme.Qbuf + iThread; 
		*(res[iThread].Q) = zeroM; //(FLEX_ASYM_CMATRIX<double>)0;
		res[iThread].f_mu_accumulate = (void*)(&TMU_accumulate);
	}
#if _CHECK_TIME_STAMP_
	mt.elapse(true, 1); show_log(", ", false);
#endif

	MacroMoleculeOperate1_mv<VSPME<_VMUatom>, RECP_SPME_RES >((void*)(&VSPME_TQ_accum_mu), vspme, 0, vspme.va.n, res, nThreads);
	int i, j, k;
	if (nThreads == 1) {
		for (i = 0; i < vspme.TQ.nx; i++) {for (j = 0; j < vspme.TQ.ny; j++) {
			memcpy(vspme.TQ.m[i].m[j].m, res[0].Q->m[i].m[j].m, vspme.Q.nz * SIZE_DOUBLE);
		}}
	}
	else {
		for (i = 0; i < vspme.TQ.nx; i++) {for (j = 0; j < vspme.TQ.ny; j++) {for (k = 0; k < vspme.TQ.nz; k++) {
			vspme.TQ.m[i].m[j].m[k] = res[0].Q->m[i].m[j].m[k];
			for (iThread = 1; iThread < nThreads; iThread += 1) vspme.TQ.m[i].m[j].m[k] += res[iThread].Q->m[i].m[j].m[k];
		}}}
	}

#if _CHECK_TIME_STAMP_
	mt.elapse(true, 2); show_log("), ", false);
#endif
}

/***************************************************************************************
	CALCULATE E_recp -- E_recp = F_recp / q
	Mesh dimension is assumed to be bigger than the order of the B-spline function
	In another word, the spline points is within the nearest neighbor cell only
***************************************************************************************/

void cal_E_recp(VSPME<_VQatom>& vspme, bool bF, int nThreads) {// E of atom is accumulated ,without initilized
	ERECP_SPME_VAR var;
	var.bF = bF; // calculate force at the same time
	var.fFunc = (void*)(&F_spme_accumulate_mu);
	var.BCQ = &(vspme.BCQ);
	void* f_Erecp = (void*)(&E_recp_q);
	MacroMoleculeOperate1<VSPME<_VQatom>, ERECP_SPME_VAR>(f_Erecp, vspme, 0, vspme.va.n, var, nThreads);
}

void cal_E_recp(VSPME<_VMUatom>& vspme, bool bF, int nThreads) {// E of atom is accumulated ,without initilized
	ERECP_SPME_VAR var;
	var.bF = bF; // calculate force at the same time
	var.fFunc = (void*)(&F_spme_accumulate_mu);
	var.BCQ = &(vspme.BCQ);
	void* f_Erecp = (void*)(&E_recp_mu);
	MacroMoleculeOperate1<VSPME<_VMUatom>, ERECP_SPME_VAR>(f_Erecp, vspme, 0, vspme.va.n, var, nThreads);
}

void cal_F_recp(VSPME<_VMUatom>& vspme, int nThreads) {// F of atom is accumulated ,without initilized
	ERECP_SPME_VAR var;
	var.bE = false;
	var.bF = true; // calculate force at the same time
	var.fFunc = (void*)(&F_spme_accumulate_mu);
	var.BCQ = &(vspme.BCQ);
	void* f_Erecp = (void*)(&E_recp_mu);
	MacroMoleculeOperate1<VSPME<_VMUatom>, ERECP_SPME_VAR>(f_Erecp, vspme, 0, vspme.va.n, var, nThreads);
}

void cal_E_recp(VSPME<_VQatom>& vspme, int na1, int na2, bool bF, int nThreads) {// E of atom is accumulated ,without initilized
	ERECP_SPME_VAR var;
	var.bF = bF; // calculate force at the same time
	var.fFunc = (void*)(&F_spme_accumulate_mu);
	var.BCQ = &(vspme.BCQ);
	void* f_Erecp = (void*)(&E_recp_q);
	MacroMoleculeOperate1<VSPME<_VQatom>, ERECP_SPME_VAR>(f_Erecp, vspme, na1, na2, var, nThreads);
}

void cal_E_recp(VSPME<_VMUatom>& vspme, int na1, int na2, bool bF, int nThreads) {// E of atom is accumulated ,without initilized
	ERECP_SPME_VAR var;
	var.bF = bF; // calculate force at the same time
	var.fFunc = (void*)(&F_spme_accumulate_mu);
	var.BCQ = &(vspme.BCQ);
	void* f_Erecp = (void*)(&E_recp_mu);
	MacroMoleculeOperate1<VSPME<_VMUatom>, ERECP_SPME_VAR>(f_Erecp, vspme, na1, na2, var, nThreads);
}

void cal_F_recp(VSPME<_VMUatom>& vspme, int na1, int na2, int nThreads) {// F of atom is accumulated ,without initilized
	ERECP_SPME_VAR var;
	var.bE = false;
	var.bF = true; // calculate force at the same time
	var.fFunc = (void*)(&F_spme_accumulate_mu);
	var.BCQ = &(vspme.BCQ);
	void* f_Erecp = (void*)(&E_recp_mu);
	MacroMoleculeOperate1<VSPME<_VMUatom>, ERECP_SPME_VAR>(f_Erecp, vspme, na1, na2, var, nThreads);
}


void setup_vqatom(_VQatom *va, MATOM *patom) {
	va->eNeutral = patom->eNeutral;
	va->r = &(patom->rg);
	va->q = &(patom->c);
	va->E = &(patom->E);
	if (::bMTForce) va->F = &(patom->f.m[_IF_ES]);
	else va->F = &(patom->F);
	va->mMP = 0;
}

void init_spme(MMOL_MD_CELL *mdcell, VSPME<_VQatom> &vspme) {
	init_spme_tmpl<_VQatom>(mdcell, vspme, (void*)(&setup_vqatom));
}

void setup_vmuatom(_VMUatom *va, MATOM *patom) {
	va->eNeutral = patom->eNeutral;
	va->r = &(patom->rg);
	va->q = &(patom->c);
	va->E = &(patom->E);
	if (::bMTForce) va->F = &(patom->f.m[_IF_ES]);
	else va->F = &(patom->F);

	va->mMP = patom->mMP;
	if (patom->mMP == 1) {
		va->mu = &(patom->mu);
	}
	else {
		va->mu = &(patom->mu);
	}
}

void init_spme(MMOL_MD_CELL *mdcell, VSPME<_VMUatom> &vspme) {
	init_spme_tmpl<_VMUatom>(mdcell, vspme, (void*)(&setup_vmuatom));
}


void init_spme(MMOL_VMD_CELL *mdcell, VSPME<_VQatom> &vspme) {
	init_spme_tmpl_vcell<_VQatom>(mdcell, vspme, (void*)(&setup_vqatom));
}
void init_spme(MMOL_VMD_CELL *mdcell, VSPME<_VMUatom> &vspme) {
	init_spme_tmpl_vcell<_VMUatom>(mdcell, vspme, (void*)(&setup_vmuatom));
}



void setup_vmuatom_induced(_VMUatom *va, MATOM *patom, float *q) {
	va->eNeutral = patom->eNeutral;
	va->r = &(patom->rg);
	//va->q = &(patom->c);
	va->q = q;
	//va->E = &(patom->E);
	if (::bMTForce) va->F = &(patom->f.m[_IF_ES]);
	else va->F = &(patom->F);

	va->mMP = patom->mMP;
	if (patom->mMP == 1) {
		va->mu = &(patom->ind_mu);
	}
	else {
		//va->mu = NULL;
		va->mu = &(patom->ind_mu);
	}
}

void init_spme_induced(MMOL_MD_CELL *mdcell, VSPME<_VMUatom> &vspme, float* q) {
	init_spme_induced_tmpl<_VMUatom>(mdcell, vspme, (void*)(&setup_vmuatom_induced), q);
}



void SressTensor1(VSPME<_VMUatom> *vspme, int na1, int na2, VECTOR<4> &st) {
	int ia;
	_VMUatom *patom;
	double STxx = 0, STyy = 0, STzz = 0, STtr = 0;
	double *mu = NULL, *E = NULL;
	for (ia = na1; ia <= na2; ia++) {
		if (ia >= vspme->va.n) break;
		patom = vspme->va.m + ia;
		mu = patom->mu->v;
		E = patom->E_recp.v;
		STxx -= mu[0] * E[0];
		STyy -= mu[1] * E[1];
		STzz -= mu[2] * E[2];
	}
	STxx *= ::fUnit_estat_mass;
	STyy *= ::fUnit_estat_mass;
	STzz *= ::fUnit_estat_mass;
	STtr = STxx + STyy + STzz;
	st.v[0] = STxx; st.v[1] = STyy; st.v[2] = STzz; st.v[3] = STtr;
}

void cal_StressTensor(VSPME<_VMUatom>& vspme, int nThreads) {  // Receiprocal Energy is calculated !
	cal_TQ_spme(vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
	vspme.cal_FTQ();
	vspme.cal_StressTensor_recp_0();

	VECTOR<4> st1;
	MacroMoleculeOperate2<VSPME<_VMUatom>, VECTOR<4> >((void*)(&SressTensor1), vspme, 0, vspme.va.n, nThreads, st1);
	if (vspme.bST_diagonal) {
		vspme.STxx += st1.v[0];
		vspme.STyy += st1.v[1];
		vspme.STzz += st1.v[2];
	}
	vspme.STtr += st1.v[3];
}




//***********************************************************
//     2-d Slab with Image Charge / Atom 
//***********************************************************


void init_spme(ARRAY<MATOM*> &atom_r, ARRAY<MATOM*> &img1, ARRAY<MATOM*> &img2, VSPME_ImageSlab<_VQatom> &vspme) {
	init_spme_tmpl_ImageSlab<_VQatom>(atom_r, img1, img2, vspme, (void*)(&setup_vqatom));
}
void init_spme(ARRAY<MATOM*> &atom_r, ARRAY<MATOM*> &img1, ARRAY<MATOM*> &img2, VSPME_ImageSlab<_VMUatom> &vspme) {
	init_spme_tmpl_ImageSlab<_VMUatom>(atom_r, img1, img2, vspme, (void*)(&setup_vmuatom));
}


// E & F calculation has to be done for both real and imagine atoms, and we assume they are done already
//f_A = f_A_ri1 + (f_Ai_ri1 - f_Ai_i1)_T + f_A_ri2 + (f_Ai_ri2 - f_Ai_i2)_T - f_A_B_r
void ImageSlab_dump_EF_Q(VSPME_ImageSlab<_VQatom> &vspme, int n1, int n2, bool& bF, bool bCleanBuff) {
	int i, i1, i2, j;
	double F, v;
	for (i = n1; i <= n2; i++) {
		if (i >= vspme.Ivspme1.n_real) break;
		i1 = i + vspme.Ivspme1.n_real; i2 = i + vspme.Ivspme2.n_real;
		vspme.Ivspme1.va.m[i].E->v[0] += vspme.Ivspme1.va.m[i].E_recp.v[0] + vspme.Ivspme2.va.m[i].E_recp.v[0] - vspme.vspme.va.m[i].E_recp.v[0];
		vspme.Ivspme1.va.m[i].E->v[0] += vspme.Ivspme1.va.m[i1].E_recp.v[0] + vspme.Ivspme2.va.m[i2].E_recp.v[0]; // on image charge
		vspme.Ivspme1.va.m[i].E->v[0] -= vspme.img1_vspme.va.m[i].E_recp.v[0] + vspme.img2_vspme.va.m[i].E_recp.v[0];
		v = 0;
		for (j = 0; j < vspme.Ivspme1.va.m[i].aE.n; j++) {
			v += vspme.Ivspme1.va.m[i].aE.m[j].v[0] + vspme.Ivspme2.va.m[i].aE.m[j].v[0] - vspme.vspme.va.m[i].aE.m[j].v[0];
			v += vspme.Ivspme1.va.m[i1].aE.m[j].v[0] + vspme.Ivspme2.va.m[i2].aE.m[j].v[0]; // on image charge
			v -= vspme.img1_vspme.va.m[i].aE.m[j].v[0] + vspme.img2_vspme.va.m[i].aE.m[j].v[0];
		}
		vspme.Ivspme1.va.m[i].E->v[0] += v;

		vspme.Ivspme1.va.m[i].E->v[1] += vspme.Ivspme1.va.m[i].E_recp.v[1] + vspme.Ivspme2.va.m[i].E_recp.v[1] - vspme.vspme.va.m[i].E_recp.v[1];
		vspme.Ivspme1.va.m[i].E->v[1] += vspme.Ivspme1.va.m[i1].E_recp.v[1] + vspme.Ivspme2.va.m[i2].E_recp.v[1]; // on image charge
		vspme.Ivspme1.va.m[i].E->v[1] -= vspme.img1_vspme.va.m[i].E_recp.v[1] + vspme.img2_vspme.va.m[i].E_recp.v[1];
		v = 0;
		for (j = 0; j < vspme.Ivspme1.va.m[i].aE.n; j++) {
			v += vspme.Ivspme1.va.m[i].aE.m[j].v[1] + vspme.Ivspme2.va.m[i].aE.m[j].v[1] - vspme.vspme.va.m[i].aE.m[j].v[1];
			v += vspme.Ivspme1.va.m[i1].aE.m[j].v[1] + vspme.Ivspme2.va.m[i2].aE.m[j].v[1]; // on image charge
			v -= vspme.img1_vspme.va.m[i].aE.m[j].v[1] + vspme.img2_vspme.va.m[i].aE.m[j].v[1];
		}
		vspme.Ivspme1.va.m[i].E->v[1] += v;

		vspme.Ivspme1.va.m[i].E->v[2] += vspme.Ivspme1.va.m[i].E_recp.v[2] + vspme.Ivspme2.va.m[i].E_recp.v[2] - vspme.vspme.va.m[i].E_recp.v[2];
		vspme.Ivspme1.va.m[i].E->v[2] -= vspme.Ivspme1.va.m[i1].E_recp.v[2] + vspme.Ivspme2.va.m[i2].E_recp.v[2]; // on image charge
		vspme.Ivspme1.va.m[i].E->v[2] += vspme.img1_vspme.va.m[i].E_recp.v[2] + vspme.img2_vspme.va.m[i].E_recp.v[2];
		v = 0;
		for (j = 0; j < vspme.Ivspme1.va.m[i].aE.n; j++) {
			v += vspme.Ivspme1.va.m[i].aE.m[j].v[2] + vspme.Ivspme2.va.m[i].aE.m[j].v[2] - vspme.vspme.va.m[i].aE.m[j].v[2];
			v -= vspme.Ivspme1.va.m[i1].aE.m[j].v[2] + vspme.Ivspme2.va.m[i2].aE.m[j].v[2]; // on image charge
			v += vspme.img1_vspme.va.m[i].aE.m[j].v[2] + vspme.img2_vspme.va.m[i].aE.m[j].v[2];
		}
		vspme.Ivspme1.va.m[i].E->v[2] += v;

		if (bF) {
			F = vspme.Ivspme1.va.m[i].F_recp.v[0] + vspme.Ivspme2.va.m[i].F_recp.v[0] - vspme.vspme.va.m[i].F_recp.v[0];
			F += vspme.Ivspme1.va.m[i1].F_recp.v[0] + vspme.Ivspme2.va.m[i2].F_recp.v[0]; // on image charge
			F -= vspme.img1_vspme.va.m[i].F_recp.v[0] + vspme.img2_vspme.va.m[i].F_recp.v[0];
			for (j = 0; j < vspme.Ivspme1.va.m[i].aF.n; j++) {
				F += vspme.Ivspme1.va.m[i].aF.m[j].v[0] + vspme.Ivspme2.va.m[i].aF.m[j].v[0] - vspme.vspme.va.m[i].aF.m[j].v[0];
				F += vspme.Ivspme1.va.m[i1].aF.m[j].v[0] + vspme.Ivspme2.va.m[i2].aF.m[j].v[0]; // on image charge
				F -= vspme.img1_vspme.va.m[i].aF.m[j].v[0] + vspme.img2_vspme.va.m[i].aF.m[j].v[0];
			}
			vspme.Ivspme1.va.m[i].F->v[0] += F * fUnit_estat_mass;

			F = vspme.Ivspme1.va.m[i].F_recp.v[1] + vspme.Ivspme2.va.m[i].F_recp.v[1] - vspme.vspme.va.m[i].F_recp.v[1];
			F += vspme.Ivspme1.va.m[i1].F_recp.v[1] + vspme.Ivspme2.va.m[i2].F_recp.v[1]; // on image charge
			F -= vspme.img1_vspme.va.m[i].F_recp.v[1] + vspme.img2_vspme.va.m[i].F_recp.v[1];
			for (j = 0; j < vspme.Ivspme1.va.m[i].aF.n; j++) {
				F += vspme.Ivspme1.va.m[i].aF.m[j].v[1] + vspme.Ivspme2.va.m[i].aF.m[j].v[1] - vspme.vspme.va.m[i].aF.m[j].v[1];
				F += vspme.Ivspme1.va.m[i1].aF.m[j].v[1] + vspme.Ivspme2.va.m[i2].aF.m[j].v[1]; // on image charge
				F -= vspme.img1_vspme.va.m[i].aF.m[j].v[1] + vspme.img2_vspme.va.m[i].aF.m[j].v[1];
			}
			vspme.Ivspme1.va.m[i].F->v[1] += F * fUnit_estat_mass;

			F = vspme.Ivspme1.va.m[i].F_recp.v[2] + vspme.Ivspme2.va.m[i].F_recp.v[2] - vspme.vspme.va.m[i].F_recp.v[2];
			F -= vspme.Ivspme1.va.m[i1].F_recp.v[2] + vspme.Ivspme2.va.m[i2].F_recp.v[2]; // on image charge
			F += vspme.img1_vspme.va.m[i].F_recp.v[2] + vspme.img2_vspme.va.m[i].F_recp.v[2];
			for (j = 0; j < vspme.Ivspme1.va.m[i].aF.n; j++) {
				F += vspme.Ivspme1.va.m[i].aF.m[j].v[2] + vspme.Ivspme2.va.m[i].aF.m[j].v[2] - vspme.vspme.va.m[i].aF.m[j].v[2];
				F -= vspme.Ivspme1.va.m[i1].aF.m[j].v[2] + vspme.Ivspme2.va.m[i2].aF.m[j].v[2]; // on image charge
				F += vspme.img1_vspme.va.m[i].aF.m[j].v[2] + vspme.img2_vspme.va.m[i].aF.m[j].v[2];
			}
			vspme.Ivspme1.va.m[i].F->v[2] += F * fUnit_estat_mass;
		}
		if (bCleanBuff) {
			vspme.Ivspme1.va.m[i].reset_local_EF(bF);
			vspme.Ivspme2.va.m[i].reset_local_EF(bF);
			vspme.img1_vspme.va.m[i].reset_local_EF(bF);
			vspme.img2_vspme.va.m[i].reset_local_EF(bF);
			vspme.vspme.va.m[i].reset_local_EF(bF);
		}
	}
}

void ImageSlab_dump_EF_MU(VSPME_ImageSlab<_VMUatom> &vspme, int n1, int n2, bool& bF, bool bCleanBuff) {
	int i, i1, i2, j;
	double F, v;
	for (i = n1; i <= n2; i++) {
		if (i >= vspme.Ivspme1.n_real) break;
		i1 = i + vspme.Ivspme1.n_real; i2 = i + vspme.Ivspme2.n_real;
		vspme.Ivspme1.va.m[i].E->v[0] += vspme.Ivspme1.va.m[i].E_recp.v[0] + vspme.Ivspme2.va.m[i].E_recp.v[0] - vspme.vspme.va.m[i].E_recp.v[0];
		vspme.Ivspme1.va.m[i].E->v[0] += vspme.Ivspme1.va.m[i1].E_recp.v[0] + vspme.Ivspme2.va.m[i2].E_recp.v[0]; // on image charge
		vspme.Ivspme1.va.m[i].E->v[0] -= vspme.img1_vspme.va.m[i].E_recp.v[0] + vspme.img2_vspme.va.m[i].E_recp.v[0];
		v = 0;
		for (j = 0; j < vspme.Ivspme1.va.m[i].aE.n; j++) {
			v += vspme.Ivspme1.va.m[i].aE.m[j].v[0] + vspme.Ivspme2.va.m[i].aE.m[j].v[0] - vspme.vspme.va.m[i].aE.m[j].v[0];
			v += vspme.Ivspme1.va.m[i1].aE.m[j].v[0] + vspme.Ivspme2.va.m[i2].aE.m[j].v[0]; // on image charge
			v -= vspme.img1_vspme.va.m[i].aE.m[j].v[0] + vspme.img2_vspme.va.m[i].aE.m[j].v[0];
		}
		vspme.Ivspme1.va.m[i].E->v[0] += v;

		vspme.Ivspme1.va.m[i].E->v[1] += vspme.Ivspme1.va.m[i].E_recp.v[1] + vspme.Ivspme2.va.m[i].E_recp.v[1] - vspme.vspme.va.m[i].E_recp.v[1];
		vspme.Ivspme1.va.m[i].E->v[1] += vspme.Ivspme1.va.m[i1].E_recp.v[1] + vspme.Ivspme2.va.m[i2].E_recp.v[1]; // on image charge
		vspme.Ivspme1.va.m[i].E->v[1] -= vspme.img1_vspme.va.m[i].E_recp.v[1] + vspme.img2_vspme.va.m[i].E_recp.v[1];
		v = 0;
		for (j = 0; j < vspme.Ivspme1.va.m[i].aE.n; j++) {
			v += vspme.Ivspme1.va.m[i].aE.m[j].v[1] + vspme.Ivspme2.va.m[i].aE.m[j].v[1] - vspme.vspme.va.m[i].aE.m[j].v[1];
			v += vspme.Ivspme1.va.m[i1].aE.m[j].v[1] + vspme.Ivspme2.va.m[i2].aE.m[j].v[1]; // on image charge
			v -= vspme.img1_vspme.va.m[i].aE.m[j].v[1] + vspme.img2_vspme.va.m[i].aE.m[j].v[1];
		}
		vspme.Ivspme1.va.m[i].E->v[1] += v;

		vspme.Ivspme1.va.m[i].E->v[2] += vspme.Ivspme1.va.m[i].E_recp.v[2] + vspme.Ivspme2.va.m[i].E_recp.v[2] - vspme.vspme.va.m[i].E_recp.v[2];
		vspme.Ivspme1.va.m[i].E->v[2] -= vspme.Ivspme1.va.m[i1].E_recp.v[2] + vspme.Ivspme2.va.m[i2].E_recp.v[2]; // on image charge
		vspme.Ivspme1.va.m[i].E->v[2] += vspme.img1_vspme.va.m[i].E_recp.v[2] + vspme.img2_vspme.va.m[i].E_recp.v[2];
		v = 0;
		for (j = 0; j < vspme.Ivspme1.va.m[i].aE.n; j++) {
			v += vspme.Ivspme1.va.m[i].aE.m[j].v[2] + vspme.Ivspme2.va.m[i].aE.m[j].v[2] - vspme.vspme.va.m[i].aE.m[j].v[2];
			v -= vspme.Ivspme1.va.m[i1].aE.m[j].v[2] + vspme.Ivspme2.va.m[i2].aE.m[j].v[2]; // on image charge
			v += vspme.img1_vspme.va.m[i].aE.m[j].v[2] + vspme.img2_vspme.va.m[i].aE.m[j].v[2];
		}
		vspme.Ivspme1.va.m[i].E->v[2] += v;

		if (bF) {
			F = vspme.Ivspme1.va.m[i].F_recp.v[0] + vspme.Ivspme2.va.m[i].F_recp.v[0] - vspme.vspme.va.m[i].F_recp.v[0];
			F += vspme.Ivspme1.va.m[i1].F_recp.v[0] + vspme.Ivspme2.va.m[i2].F_recp.v[0]; // on image charge
			F -= vspme.img1_vspme.va.m[i].F_recp.v[0] + vspme.img2_vspme.va.m[i].F_recp.v[0];
			for (j = 0; j < vspme.Ivspme1.va.m[i].aF.n; j++) {
				F += vspme.Ivspme1.va.m[i].aF.m[j].v[0] + vspme.Ivspme2.va.m[i].aF.m[j].v[0] - vspme.vspme.va.m[i].aF.m[j].v[0];
				F += vspme.Ivspme1.va.m[i1].aF.m[j].v[0] + vspme.Ivspme2.va.m[i2].aF.m[j].v[0]; // on image charge
				F -= vspme.img1_vspme.va.m[i].aF.m[j].v[0] + vspme.img2_vspme.va.m[i].aF.m[j].v[0];
			}
			vspme.Ivspme1.va.m[i].F->v[0] += F * fUnit_estat_mass;

			F = vspme.Ivspme1.va.m[i].F_recp.v[1] + vspme.Ivspme2.va.m[i].F_recp.v[1] - vspme.vspme.va.m[i].F_recp.v[1];
			F += vspme.Ivspme1.va.m[i1].F_recp.v[1] + vspme.Ivspme2.va.m[i2].F_recp.v[1]; // on image charge
			F -= vspme.img1_vspme.va.m[i].F_recp.v[1] + vspme.img2_vspme.va.m[i].F_recp.v[1];
			for (j = 0; j < vspme.Ivspme1.va.m[i].aF.n; j++) {
				F += vspme.Ivspme1.va.m[i].aF.m[j].v[1] + vspme.Ivspme2.va.m[i].aF.m[j].v[1] - vspme.vspme.va.m[i].aF.m[j].v[1];
				F += vspme.Ivspme1.va.m[i1].aF.m[j].v[1] + vspme.Ivspme2.va.m[i2].aF.m[j].v[1]; // on image charge
				F -= vspme.img1_vspme.va.m[i].aF.m[j].v[1] + vspme.img2_vspme.va.m[i].aF.m[j].v[1];
			}
			vspme.Ivspme1.va.m[i].F->v[1] += F * fUnit_estat_mass;

			F = vspme.Ivspme1.va.m[i].F_recp.v[2] + vspme.Ivspme2.va.m[i].F_recp.v[2] - vspme.vspme.va.m[i].F_recp.v[2];
			F -= vspme.Ivspme1.va.m[i1].F_recp.v[2] + vspme.Ivspme2.va.m[i2].F_recp.v[2]; // on image charge
			F += vspme.img1_vspme.va.m[i].F_recp.v[2] + vspme.img2_vspme.va.m[i].F_recp.v[2];
			for (j = 0; j < vspme.Ivspme1.va.m[i].aF.n; j++) {
				F += vspme.Ivspme1.va.m[i].aF.m[j].v[2] + vspme.Ivspme2.va.m[i].aF.m[j].v[2] - vspme.vspme.va.m[i].aF.m[j].v[2];
				F -= vspme.Ivspme1.va.m[i1].aF.m[j].v[2] + vspme.Ivspme2.va.m[i2].aF.m[j].v[2]; // on image charge
				F += vspme.img1_vspme.va.m[i].aF.m[j].v[2] + vspme.img2_vspme.va.m[i].aF.m[j].v[2];
			}
			vspme.Ivspme1.va.m[i].F->v[2] += F * fUnit_estat_mass;
		}
		if (bCleanBuff) {
			vspme.Ivspme1.va.m[i].reset_local_EF(bF);
			vspme.Ivspme2.va.m[i].reset_local_EF(bF);
			vspme.img1_vspme.va.m[i].reset_local_EF(bF);
			vspme.img2_vspme.va.m[i].reset_local_EF(bF);
			vspme.vspme.va.m[i].reset_local_EF(bF);
		}
	}
}

/*
void ImageSlab_dump_EF_MU(VSPME_ImageSlab<_VMUatom> &vspme, int n1, int n2, bool& bF) {
	int i, i1, i2;
	double F;
	for (i = n1; i <= n2; i++) {
		if (i >= vspme.Ivspme1.n_real) break;
		i1 = i + vspme.Ivspme1.n_real; i2 = i + vspme.Ivspme2.n_real;
		vspme.Ivspme1.va.m[i].E->v[0] += vspme.Ivspme1.va.m[i].E_recp.v[0] + vspme.Ivspme2.va.m[i].E_recp.v[0] - vspme.vspme.va.m[i].E_recp.v[0];
		vspme.Ivspme1.va.m[i].E->v[0] += vspme.Ivspme1.va.m[i1].E_recp.v[0] + vspme.Ivspme2.va.m[i2].E_recp.v[0]; // on image charge
		vspme.Ivspme1.va.m[i].E->v[0] -= vspme.img1_vspme.va.m[i].E_recp.v[0] + vspme.img2_vspme.va.m[i].E_recp.v[0];

		vspme.Ivspme1.va.m[i].E->v[1] += vspme.Ivspme1.va.m[i].E_recp.v[1] + vspme.Ivspme2.va.m[i].E_recp.v[1] - vspme.vspme.va.m[i].E_recp.v[1];
		vspme.Ivspme1.va.m[i].E->v[1] += vspme.Ivspme1.va.m[i1].E_recp.v[1] + vspme.Ivspme2.va.m[i2].E_recp.v[1]; // on image charge
		vspme.Ivspme1.va.m[i].E->v[1] -= vspme.img1_vspme.va.m[i].E_recp.v[1] + vspme.img2_vspme.va.m[i].E_recp.v[1];

		vspme.Ivspme1.va.m[i].E->v[2] += vspme.Ivspme1.va.m[i].E_recp.v[2] + vspme.Ivspme2.va.m[i].E_recp.v[2] - vspme.vspme.va.m[i].E_recp.v[2];
		vspme.Ivspme1.va.m[i].E->v[2] -= vspme.Ivspme1.va.m[i1].E_recp.v[2] + vspme.Ivspme2.va.m[i2].E_recp.v[2]; // on image charge
		vspme.Ivspme1.va.m[i].E->v[2] += vspme.img1_vspme.va.m[i].E_recp.v[2] + vspme.img2_vspme.va.m[i].E_recp.v[2];

		if (bF) {
			F = vspme.Ivspme1.va.m[i].F_recp.v[0] + vspme.Ivspme2.va.m[i].F_recp.v[0] - vspme.vspme.va.m[i].F_recp.v[0];
			F += vspme.Ivspme1.va.m[i1].F_recp.v[0] + vspme.Ivspme2.va.m[i2].F_recp.v[0]; // on image charge
			F -= vspme.img1_vspme.va.m[i].F_recp.v[0] + vspme.img2_vspme.va.m[i].F_recp.v[0];
			vspme.Ivspme1.va.m[i].F->v[0] += F * fUnit_estat_mass;

			F = vspme.Ivspme1.va.m[i].F_recp.v[1] + vspme.Ivspme2.va.m[i].F_recp.v[1] - vspme.vspme.va.m[i].F_recp.v[1];
			F += vspme.Ivspme1.va.m[i1].F_recp.v[1] + vspme.Ivspme2.va.m[i2].F_recp.v[1]; // on image charge
			F -= vspme.img1_vspme.va.m[i].F_recp.v[1] + vspme.img2_vspme.va.m[i].F_recp.v[1];
			vspme.Ivspme1.va.m[i].F->v[1] += F * fUnit_estat_mass;

			F = vspme.Ivspme1.va.m[i].F_recp.v[2] + vspme.Ivspme2.va.m[i].F_recp.v[2] - vspme.vspme.va.m[i].F_recp.v[2];
			F -= vspme.Ivspme1.va.m[i1].F_recp.v[2] + vspme.Ivspme2.va.m[i2].F_recp.v[2]; // on image charge
			F += vspme.img1_vspme.va.m[i].F_recp.v[2] + vspme.img2_vspme.va.m[i].F_recp.v[2];
			vspme.Ivspme1.va.m[i].F->v[2] += F * fUnit_estat_mass;
		}
	}
}
*/

void Slab_dump_EF_Q(VSPME<_VQatom> &vspme, int n1, int n2, bool& bF, bool bCleanBuff) {
	int i, j;
	double v;
	for (i = n1; i <= n2; i++) {
		if (i >= vspme.va.n) break;
		vspme.va.m[i].E->v[0] += vspme.va.m[i].E_recp.v[0];
		vspme.va.m[i].E->v[1] += vspme.va.m[i].E_recp.v[1];
		vspme.va.m[i].E->v[2] += vspme.va.m[i].E_recp.v[2];
		for (j = 0; j < vspme.va.m[i].aE.n; j++) {
			vspme.va.m[i].E->v[0] += vspme.va.m[i].aE.m[j].v[0];
			vspme.va.m[i].E->v[1] += vspme.va.m[i].aE.m[j].v[1];
			vspme.va.m[i].E->v[2] += vspme.va.m[i].aE.m[j].v[2];
		}
		if (bF) {
			v = vspme.va.m[i].F_recp.v[0];
			for (j = 0; j < vspme.va.m[i].aF.n; j++) v += vspme.va.m[i].aF.m[j].v[0];
			vspme.va.m[i].F->v[0] += v * fUnit_estat_mass;

			v = vspme.va.m[i].F_recp.v[1];
			for (j = 0; j < vspme.va.m[i].aF.n; j++) v += vspme.va.m[i].aF.m[j].v[1];
			vspme.va.m[i].F->v[1] += v * fUnit_estat_mass;

			v = vspme.va.m[i].F_recp.v[2];
			for (j = 0; j < vspme.va.m[i].aF.n; j++) v += vspme.va.m[i].aF.m[j].v[2];
			vspme.va.m[i].F->v[2] += v * fUnit_estat_mass;
		}
		if (bCleanBuff) vspme.va.m[i].reset_local_EF(bF);
	}
}

void Slab_dump_EF_MU(VSPME<_VMUatom> &vspme, int n1, int n2, bool& bF, bool bCleanBuff) {
	int i, j;
	double v;
	for (i = n1; i <= n2; i++) {
		if (i >= vspme.va.n) break;
		vspme.va.m[i].E->v[0] += vspme.va.m[i].E_recp.v[0];
		vspme.va.m[i].E->v[1] += vspme.va.m[i].E_recp.v[1];
		vspme.va.m[i].E->v[2] += vspme.va.m[i].E_recp.v[2];
		for (j = 0; j < vspme.va.m[i].aE.n; j++) {
			vspme.va.m[i].E->v[0] += vspme.va.m[i].aE.m[j].v[0];
			vspme.va.m[i].E->v[1] += vspme.va.m[i].aE.m[j].v[1];
			vspme.va.m[i].E->v[2] += vspme.va.m[i].aE.m[j].v[2];
		}
		if (bF) {
			v = vspme.va.m[i].F_recp.v[0];
			for (j = 0; j < vspme.va.m[i].aF.n; j++) v += vspme.va.m[i].aF.m[j].v[0];
			vspme.va.m[i].F->v[0] += v * fUnit_estat_mass;

			v = vspme.va.m[i].F_recp.v[1];
			for (j = 0; j < vspme.va.m[i].aF.n; j++) v += vspme.va.m[i].aF.m[j].v[1];
			vspme.va.m[i].F->v[1] += v * fUnit_estat_mass;

			v = vspme.va.m[i].F_recp.v[2];
			for (j = 0; j < vspme.va.m[i].aF.n; j++) v += vspme.va.m[i].aF.m[j].v[2];
			vspme.va.m[i].F->v[2] += v * fUnit_estat_mass;
		}
		if (bCleanBuff) vspme.va.m[i].reset_local_EF(bF);
	}
}

void Slab_dump_F_exclude_induced_mu(VSPME<_VMUatom> &vspme, int n1, int n2, bool &status, bool bCleanBuff) {
	int i, j;
	double v;
	for (i = n1; i <= n2; i++) {
		if (i >= vspme.va.n) break;
		//vspme.va.m[i].E->v[0] += vspme.va.m[i].E_recp.v[0];
		//vspme.va.m[i].E->v[1] += vspme.va.m[i].E_recp.v[1];
		//vspme.va.m[i].E->v[2] += vspme.va.m[i].E_recp.v[2];
		//if (bF) {
			//vspme.va.m[i].F->v[0] -= vspme.va.m[i].F_recp.v[0] * fUnit_estat_mass;
			//vspme.va.m[i].F->v[1] -= vspme.va.m[i].F_recp.v[1] * fUnit_estat_mass;
			//vspme.va.m[i].F->v[2] -= vspme.va.m[i].F_recp.v[2] * fUnit_estat_mass;

			v = vspme.va.m[i].F_recp.v[0];
			for (j = 0; j < vspme.va.m[i].aF.n; j++) v += vspme.va.m[i].aF.m[j].v[0];
			vspme.va.m[i].F->v[0] -= v * fUnit_estat_mass;

			v = vspme.va.m[i].F_recp.v[1];
			for (j = 0; j < vspme.va.m[i].aF.n; j++) v += vspme.va.m[i].aF.m[j].v[1];
			vspme.va.m[i].F->v[1] -= v * fUnit_estat_mass;

			v = vspme.va.m[i].F_recp.v[2];
			for (j = 0; j < vspme.va.m[i].aF.n; j++) v += vspme.va.m[i].aF.m[j].v[2];
			vspme.va.m[i].F->v[2] -= v * fUnit_estat_mass;
		//}
			if (bCleanBuff) vspme.va.m[i].reset_local_EF(true);
	}
}

void Slab_dumpJob_F_exclude_induced_mu(JobIndx *job, VSPME<_VMUatom> *vspme) {
	bool status = true;
	Slab_dump_F_exclude_induced_mu(*vspme, job->n1, job->n2, status, true);
}









void Slab_dump_EF_Q_Job(JobIndx *job, VSPME<_VQatom> *vspme, bool *bE, bool *bF) {
	int i, j;
	double v;
	for (i = job->n1; i <= job->n2; i++) {
		if (i >= vspme->va.n) break;
		if (*bE) {
			vspme->va.m[i].E->v[0] += vspme->va.m[i].E_recp.v[0];
			vspme->va.m[i].E->v[1] += vspme->va.m[i].E_recp.v[1];
			vspme->va.m[i].E->v[2] += vspme->va.m[i].E_recp.v[2];

			for (j = 0; j < vspme->va.m[i].aE.n; j++) {
				vspme->va.m[i].E->v[0] += vspme->va.m[i].aE.m[j].v[0];
				vspme->va.m[i].E->v[1] += vspme->va.m[i].aE.m[j].v[1];
				vspme->va.m[i].E->v[2] += vspme->va.m[i].aE.m[j].v[2];
			}
		}
		if (*bF) {
			v = vspme->va.m[i].F_recp.v[0];
			for (j = 0; j < vspme->va.m[i].aF.n; j++) v += vspme->va.m[i].aF.m[j].v[0];
			vspme->va.m[i].F->v[0] += v * fUnit_estat_mass;

			v = vspme->va.m[i].F_recp.v[1];
			for (j = 0; j < vspme->va.m[i].aF.n; j++) v += vspme->va.m[i].aF.m[j].v[1];
			vspme->va.m[i].F->v[1] += v * fUnit_estat_mass;

			v = vspme->va.m[i].F_recp.v[2];
			for (j = 0; j < vspme->va.m[i].aF.n; j++) v += vspme->va.m[i].aF.m[j].v[2];
			vspme->va.m[i].F->v[2] += v * fUnit_estat_mass;
		}
		vspme->va.m[i].reset_local_EF(*bF);
	}
}

void Slab_dump_EF_MU_Job(JobIndx *job, VSPME<_VMUatom> *vspme, bool *bE, bool *bF) {
	int i, j;
	double v;
	for (i = job->n1; i <= job->n2; i++) {
		if (i >= vspme->va.n) break;
		if (*bE) {
			vspme->va.m[i].E->v[0] += vspme->va.m[i].E_recp.v[0];
			vspme->va.m[i].E->v[1] += vspme->va.m[i].E_recp.v[1];
			vspme->va.m[i].E->v[2] += vspme->va.m[i].E_recp.v[2];

			for (j = 0; j < vspme->va.m[i].aE.n; j++) {
				vspme->va.m[i].E->v[0] += vspme->va.m[i].aE.m[j].v[0];
				vspme->va.m[i].E->v[1] += vspme->va.m[i].aE.m[j].v[1];
				vspme->va.m[i].E->v[2] += vspme->va.m[i].aE.m[j].v[2];
			}
		}
		if (*bF) {
			v = vspme->va.m[i].F_recp.v[0];
			for (j = 0; j < vspme->va.m[i].aF.n; j++) v += vspme->va.m[i].aF.m[j].v[0];
			vspme->va.m[i].F->v[0] += v * fUnit_estat_mass;

			v = vspme->va.m[i].F_recp.v[1];
			for (j = 0; j < vspme->va.m[i].aF.n; j++) v += vspme->va.m[i].aF.m[j].v[1];
			vspme->va.m[i].F->v[1] += v * fUnit_estat_mass;

			v = vspme->va.m[i].F_recp.v[2];
			for (j = 0; j < vspme->va.m[i].aF.n; j++) v += vspme->va.m[i].aF.m[j].v[2];
			vspme->va.m[i].F->v[2] += v * fUnit_estat_mass;
		}
		vspme->va.m[i].reset_local_EF(*bF);
	}
}


#endif

} // end of namespace _spme_

