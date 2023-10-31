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

BSplineFunc<2> bsp; // with 2nd order derivation

void BSplineFuncs::init(int indx) {
	if (indx < 2) return; // B-spline function has indx >= 2
	this->N = indx;
	this->nprof = indx - 1;
	this->xmax = N; // M2 -- [0, 2], M3 -- [0, 3], M4 -- [0, 4] ...
	this->npts = int(xmax / dx + 0.5) + 1;
	bsp.set(nprof, npts); dbsp.set(nprof, npts); x.SetArray(npts);

	int n, ipt;
	int mpt = int(1.0 / dx + 0.1); // max at this point for M2
	for (ipt = 0; ipt < npts; ipt++) x.m[ipt] = (double)ipt / mpt;
	for (ipt = 0; ipt <= mpt; ipt++) {
		bsp.m[0].m[ipt] = (double)ipt / mpt;
		dbsp.m[0].m[ipt] = 1.0;
	}
	for (ipt = mpt + 1; ipt <= mpt + mpt; ipt++) {
		bsp.m[0].m[ipt] = bsp.m[0].m[mpt + mpt - ipt];
		dbsp.m[0].m[ipt] = -1;
	}
	for (ipt = mpt + mpt + 1; ipt < npts; ipt++) {
		bsp.m[0].m[ipt] = 0;
		dbsp.m[0].m[ipt] = 0;
	}

	double u;
	for (n = 3; n <= N; n++) {
		for (ipt = 0; ipt < npts; ipt++) {
			u = (double)ipt / mpt;
			bsp.m[n-2].m[ipt] = u / (n - 1) * bsp.m[n-3].m[ipt];
			dbsp.m[n-2].m[ipt] = bsp.m[n-3].m[ipt];
			if (ipt >= mpt) {
				bsp.m[n-2].m[ipt] += (n - u) / (n - 1) * bsp.m[n-3].m[ipt - mpt];
				dbsp.m[n-2].m[ipt] -= bsp.m[n-3].m[ipt - mpt];
			}
		}
	}
}

double diff(BSplineFuncs &bsps, int indx, int nth, double x) { // indx B-Spline function in bsp, nth-th derivation
	if (indx > bsps.N) return 0;
	else if (indx < 2) return 0;

	double y = 0;
	switch (nth) {
	case 0:
		y = bsps.bsp_interpolate(indx, bsps.get_pt(x), x);
		break;
	case 1:
		y = bsps.dbsp_interpolate(indx, bsps.get_pt(x), x);
		break;
	default:
		y = diff(bsps, indx - 1, nth - 1, x) - diff(bsps, indx - 1, nth - 1, x - 1);
		break;
	}
	return y;
}

namespace _spme_ {
#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

void cal_Q_spme(VSPME<_VQatom>& vspme, int nThreads) { // Q value is initialized
	//vspme.Q = 0; // important
	FLEX_ASYM_CMATRIX<double> zeroM(0.0);
	RECP_SPME_RES res[MAX_THREADS];
	int iThread;
	for (iThread = 0; iThread < nThreads; iThread++) {
		res[iThread].mMP = 0; // no multipole
		res[iThread].Q = vspme.Qbuf + iThread;
		*(res[iThread].Q) = zeroM; //FLEX_ASYM_CMATRIX<double>(0.0);
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
	int iThread;
	for (iThread = 0; iThread < nThreads; iThread+=1) {
		res[iThread].mMP = 1; // dipole
		res[iThread].Q = vspme.Qbuf + iThread; 
		*(res[iThread].Q) = zeroM; //FLEX_ASYM_CMATRIX<double>(0.0);
		res[iThread].f_mu_accumulate = (void*)(&MU_accumulate);
	}

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

void cal_F_recp(VSPME<_VMUatom>& vspme, int nThreads) {// E of atom is accumulated ,without initilized
	ERECP_SPME_VAR var;
	var.bE = false; // disable electric calculation
	var.bF = true; // calculate force at the same time
	var.fFunc = (void*)(&F_spme_accumulate_mu);
	var.BCQ = &(vspme.BCQ);
	void* f_Erecp = (void*)(&E_recp_mu);
	MacroMoleculeOperate1<VSPME<_VMUatom>, ERECP_SPME_VAR>(f_Erecp, vspme, 0, vspme.va.n, var, nThreads);
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
		//va->mu = NULL;
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

void init_spme_induced(MMOL_MD_CELL *mdcell, VSPME<_VMUatom> &vspme, float *q) {
	init_spme_induced_tmpl<_VMUatom>(mdcell, vspme, (void*)(&setup_vmuatom_induced), q);
}


void dump_F_exclude_induced_mu(VSPME<_VMUatom> &vspme, int n1, int n2, bool &status, bool bCleanBuff) {
	int i;
	for (i = n1; i <= n2; i++) {
		if (i >= vspme.va.n) break;
		//vspme.va.m[i].E->v[0] += vspme.va.m[i].E_recp.v[0];
		//vspme.va.m[i].E->v[1] += vspme.va.m[i].E_recp.v[1];
		//vspme.va.m[i].E->v[2] += vspme.va.m[i].E_recp.v[2];
		//if (bF) {
			vspme.va.m[i].F->v[0] -= vspme.va.m[i].F_recp.v[0] * fUnit_estat_mass;
			vspme.va.m[i].F->v[1] -= vspme.va.m[i].F_recp.v[1] * fUnit_estat_mass;
			vspme.va.m[i].F->v[2] -= vspme.va.m[i].F_recp.v[2] * fUnit_estat_mass;
		//}
			if (bCleanBuff) vspme.va.m[i].reset_local_EF(true);
	}
}


void dumpJob_F_exclude_induced_mu(JobIndx *job, VSPME<_VMUatom> *vspme) {
	bool status = true;
	dump_F_exclude_induced_mu(*vspme, job->n1, job->n2, status, true);
}

void dump_EF_Q_Job(JobIndx *job, VSPME<_VQatom> *vspme, bool *bE, bool *bF) {
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

void dump_EF_MU_Job(JobIndx *job, VSPME<_VMUatom> *vspme, bool *bE, bool *bF) {
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

#endif

} // end of namespace _spme_

