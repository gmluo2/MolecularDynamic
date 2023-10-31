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
#include "Matrix.h"
#include "nhc.h"

#include "ZMatrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "CMM_2d.h"
#include "cluster.h"

#include "var.h"

#define INTERACT_DEF 1
#define _INTERACTION 0
#include "Interaction.h"
#include "Interact1.h"
#include "Interact2.h"

#include <fftw3.h>
#include "complex.h"

#include "NEIMO.h"
#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"
#include "MD_POLYMER.h"
#include "MD_SM.h"
#include "md_cell.h"

#include "vmd_cell.h"

#include "EwaldSum.h"
#include "spme.h"
#include "spme_2d.h"
#include "spme_interact.h"
#include "spme_interact_2d.h"

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
//long t_real = 0, t_recp = 0, t_emp = 0, t_cluster = 0, t_LJ = 0;
extern ITER_SAVE tsave;
#endif

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, _EwaldSum_real_::SURF_DIPOLE &mu_surf);
extern void cal_dipole(MMOL_VMD_CELL &mmcell, int nThreads, _EwaldSum_real_::SURF_DIPOLE &mu_surf);

extern PAIRWISE_DB< PAIRWISE<EF_DPROF> > sr_db; // pair-wise short-range interaction database

// pair-wise database about Thole-radius between ions in Thole polarization model, s = a * (alpha_i * alpha_j)^1/6
extern PAIRWISE_DB< PAIRWISE<float> > TholeRadius_db; //


#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

using namespace _evatom_;
using namespace _EwaldSum_real_;
namespace _EwaldSum_2d_ {

void EwaldSumRecpKzero(EwaldSumRecpK0Var &var, double q1, double q2, double &Ez, double &Fz, double &U) {
	Ez = -q2 * var.ferf * var.PI2_over_A;
	if (var.bEfieldOnly) return;

	double q12 = q1 * q2;
	U = -q12 * var.W * var.PI_over_A; // half of interaction
	Fz = -q12 * var.ferf * var.PI2_over_A;
}

void EwaldSumRecpKzero(EwaldSumRecpK0Var &var, double q1, double mu1z, double q2, double mu2z, double &Ez, double &Fz, double &U) {
	//Ez = -(q2 * var.ferf + mu2z * var.KAPPA2_over_SQRT_PI * var.fexp) * var.PI2_over_A;
	Ez = -q2 * var.ferf;
	if (var.fexp != 0) {
		Ez -= mu2z * var.KAPPA2_over_SQRT_PI * var.fexp;
	}
	Ez *= var.PI2_over_A;
	if (var.bEfieldOnly) return;

	double q12 = q1 * q2;
	if (var.iExIndDipole == 1) {
		U = -(q12 * var.W + (q1 * mu2z - mu1z * q2) * var.ferf /*- mu1z * mu2z * var.KAPPA2_over_SQRT_PI * var.fexp*/) * var.PI_over_A; // half of interaction
		Fz = -(q12 * var.ferf + (q1 * mu2z - mu1z * q2) * var.KAPPA2_over_SQRT_PI * var.fexp /*+ mu1z * mu2z * var.KAPPA2_over_SQRT_PI * 2 * var.kappa2 * var.dz * var.fexp*/) * var.PI2_over_A;
	}
	else {
		U = -(q12 * var.W + (q1 * mu2z - mu1z * q2) * var.ferf - mu1z * mu2z * var.KAPPA2_over_SQRT_PI * var.fexp) * var.PI_over_A; // half of interaction
		Fz = -(q12 * var.ferf + (q1 * mu2z - mu1z * q2) * var.KAPPA2_over_SQRT_PI * var.fexp + mu1z * mu2z * var.KAPPA2_over_SQRT_PI * 2 * var.kappa2 * var.dz * var.fexp) * var.PI2_over_A;
	}
}

// for 2d system, long-range interaction in Ewald Sum has non-zero term for k = 0
void EwaldSum_2D_K0_MU(ARRAY<_VMUatom>& sys, _VMUatom *atom, int natoms, EwaldSumRecpK0Var& var, InteractRes& res) {
	//using namespace _Thole_polarize_;
	//EwaldSumRealVars1 *esv = &(var.esv1);
	//TholeCorrTensor *tv = &(var.tv);

	_VMUatom *patom = NULL, *patom1 = NULL;
	double Ez, Fz; //E_Thole, F_Thole;
	double dz;
	int na, na1;

	double t1, t2;

	//double alpha_Thole = 0, a_Thole = ::_Thole_polarize_::a;

	double Utotal = 0, Uestat = 0, Urecp = 0;//, U_Thole = 0;

	/*
	char aindx[2] = {0x00, 0x00};
	unsigned int key;
	PAIRWISE<float> *TholeRadius = NULL;
	*/
#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	res.reset();

	//char msg[256] = "\0";

	Utotal = 0;

	for (na = 0; na < natoms; na++) {
		patom = atom + na;
		if (patom->eNeutral) continue;
		Uestat = 0;

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		V3zero(patom->aE.m[0]) 
		if (!var.bEfieldOnly) {V3zero(patom->aF.m[0])}

		for (na1 = 0; na1 < sys.n; na1++) {
			patom1 = sys.m + na1;
			if (patom1->eNeutral) continue; // ignore this atom for charge interaction
			dz = patom1->r->v[2] - patom->r->v[2];

			var.init_var(dz);
			EwaldSumRecpKzero(var, *(patom->q), patom->mu->v[2], *(patom1->q), patom1->mu->v[2], Ez, Fz, Urecp);

			/*
				if (var.bEfieldOnly) {
					EwaldReal_E(*esv, patom, patom1, dr, R, R2, E_real, bEwaldKickoff, t1);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->aindx; aindx[1] = patom1->par->aindx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr_E(*tv, a_Thole, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v)) {
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
								}
								V3plusV3(E_real, E_Thole, E_real)
							}
						}
					}
				}
				else {
					EwaldReal(*esv, patom, patom1, dr, R, R2, E_real, F_real, U_real, bEwaldKickoff, t1);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->aindx; aindx[1] = patom1->par->aindx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr(*tv, a_Thole, patom->c, patom->mu.v, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v, F_Thole.v, U_Thole)) {
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
									F_Thole.v[0] *= t2; F_Thole.v[1] *= t2; F_Thole.v[2] *= t2;
									U_Thole *= t2;
								}
								V3plusV3(E_real, E_Thole, E_real) V3plusV3(F_real, F_Thole, F_real) U_real += U_Thole;
							}
						}
					}
				}
				*/
			//patom->E_recp.v[2] += Ez;
			patom->aE.m[0].v[2] += Ez;
			if (!var.bEfieldOnly) {
				//patom->F_recp.v[2] += Fz;
				patom->aF.m[0].v[2] += Fz;
				Uestat += Urecp; // Urecp is half of interaction energy calculated above
				if (var.bVirial) {//Accumulate_Virial_Diagonal(res, F_real.v, dr.v);
					t1 = -Fz * dz * 0.5 * fUnit_estat_mass;
					res.STzz += t1;
					t2 = Urecp * fUnit_estat_mass;
					res.STxx += t2; res.STyy += t2;
					res.STtr += t1 + t2 + t2;
				}
			}
		}
		
		Utotal += Uestat;
	}

	if (var.bEfieldOnly) return;
	res.U = Utotal * eUnit_estat_kT;
}


// for 2d system, long-range interaction in Ewald Sum has non-zero term for k = 0
void EwaldSum_2D_K0_VMU(ARRAY<_VMUatom*>& sys, _VMUatom *atom, int natoms, EwaldSumRecpK0Var& var, InteractRes& res) {
	_VMUatom *patom = NULL, *patom1 = NULL;
	double Ez, Fz; //E_Thole, F_Thole;
	double dz;
	int na, na1;

	double t1, t2;

	double Utotal = 0, Uestat = 0, Urecp = 0;//, U_Thole = 0;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	res.reset();

	//char msg[256] = "\0";

	Utotal = 0;

	for (na = 0; na < natoms; na++) {
		patom = atom + na;
		if (patom->eNeutral) continue;
		Uestat = 0;
		V3zero(patom->aE.m[0])
		if (!var.bEfieldOnly) {V3zero(patom->aF.m[0])}

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		for (na1 = 0; na1 < sys.n; na1++) {
			patom1 = sys.m[na1];
			if (patom1->eNeutral) continue; // ignore this atom for charge interaction
			dz = patom1->r->v[2] - patom->r->v[2];

			var.init_var(dz);
			EwaldSumRecpKzero(var, *(patom->q), patom->mu->v[2], *(patom1->q), patom1->mu->v[2], Ez, Fz, Urecp);

			//patom->E_recp.v[2] += Ez;
			patom->aE.m[0].v[2] += Ez;
			if (!var.bEfieldOnly) {
				//patom->F_recp.v[2] += Fz;
				patom->aF.m[0].v[2] += Fz;
				Uestat += Urecp; // Urecp is half of interaction energy calculated above
				if (var.bVirial) {//Accumulate_Virial_Diagonal(res, F_real.v, dr.v);
					t1 = -Fz * dz * 0.5 * fUnit_estat_mass;
					res.STzz += t1;
					t2 = Urecp * fUnit_estat_mass;
					res.STxx += t2; res.STyy += t2;
					res.STtr += t1 + t2 + t2;
				}
			}
		}
		Utotal += Uestat;
	}

	if (var.bEfieldOnly) return;
	res.U = Utotal * eUnit_estat_kT;
}


// for 2d system, long-range interaction in Ewald Sum has non-zero term for k = 0
void EwaldSum_2D_K0_Q(ARRAY<_VQatom>& sys, _VQatom *atom, int natoms, EwaldSumRecpK0Var& var, InteractRes& res) {
	_VQatom *patom = NULL, *patom1 = NULL;
	double Ez, Fz; //E_Thole, F_Thole;
	double dz;
	int na, na1;

	double t1, t2;

	double Utotal = 0, Uestat = 0, Urecp = 0;//, U_Thole = 0;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	res.reset();

	//char msg[256] = "\0";

	Utotal = 0;

	for (na = 0; na < natoms; na++) {
		patom = atom + na;
		if (patom->eNeutral) continue;
		Uestat = 0;
		V3zero(patom->aE.m[0])
		if (!var.bEfieldOnly) { V3zero(patom->aF.m[0]) }

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		for (na1 = 0; na1 < sys.n; na1++) {
			patom1 = sys.m + na1;
			if (patom1->eNeutral) continue; // ignore this atom for charge interaction
			dz = patom1->r->v[2] - patom->r->v[2];

			var.init_var(dz);
			EwaldSumRecpKzero(var, *(patom->q), *(patom1->q), Ez, Fz, Urecp);

			//patom->E_recp.v[2] += Ez;
			patom->aE.m[0].v[2] += Ez;
			if (!var.bEfieldOnly) {
				//patom->F_recp.v[2] += Fz;
				patom->aF.m[0].v[2] += Fz;
				Uestat += Urecp; // Urecp is half of interaction energy calculated above
				if (var.bVirial) {//Accumulate_Virial_Diagonal(res, F_real.v, dr.v);
					t1 = -Fz * dz * 0.5 * fUnit_estat_mass;
					res.STzz += t1;
					t2 = Urecp * fUnit_estat_mass;
					res.STxx += t2; res.STyy += t2;
					res.STtr += t1 + t2 + t2;
				}
			}
		}
		Utotal += Uestat;
	}

	if (var.bEfieldOnly) return;
	res.U = Utotal * eUnit_estat_kT;
}

// for 2d system, long-range interaction in Ewald Sum has non-zero term for k = 0
void EwaldSum_2D_K0_VQ(ARRAY<_VQatom*>& sys, _VQatom *atom, int natoms, EwaldSumRecpK0Var& var, InteractRes& res) {
	_VQatom *patom = NULL, *patom1 = NULL;
	double Ez, Fz; //E_Thole, F_Thole;
	double dz;
	int na, na1;

	double t1, t2;

	double Utotal = 0, Uestat = 0, Urecp = 0;//, U_Thole = 0;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	res.reset();

	//char msg[256] = "\0";

	Utotal = 0;

	for (na = 0; na < natoms; na++) {
		patom = atom + na;
		if (patom->eNeutral) continue;
		Uestat = 0;
		V3zero(patom->aE.m[0])
		if (!var.bEfieldOnly) { V3zero(patom->aF.m[0])}

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		for (na1 = 0; na1 < sys.n; na1++) {
			patom1 = sys.m[na1];
			if (patom1->eNeutral) continue; // ignore this atom for charge interaction
			dz = patom1->r->v[2] - patom->r->v[2];

			var.init_var(dz);
			EwaldSumRecpKzero(var, *(patom->q), *(patom1->q), Ez, Fz, Urecp);

			//patom->E_recp.v[2] += Ez;
			patom->aE.m[0].v[2] += Ez;
			if (!var.bEfieldOnly) {
				//patom->F_recp.v[2] += Fz;
				patom->aF.m[0].v[2] += Fz;
				Uestat += Urecp; // Urecp is half of interaction energy calculated above
				if (var.bVirial) {//Accumulate_Virial_Diagonal(res, F_real.v, dr.v);
					t1 = -Fz * dz * 0.5 * fUnit_estat_mass;
					res.STzz += t1;
					t2 = Urecp * fUnit_estat_mass;
					res.STxx += t2; res.STyy += t2;
					res.STtr += t1 + t2 + t2;
				}
			}
		}
		Utotal += Uestat;
	}

	if (var.bEfieldOnly) return;
	res.U = Utotal * eUnit_estat_kT;
}


} // end of namespace _EwaldSum_2d_

namespace _spme_2d_ {
using namespace _EwaldSum_real_;
using namespace _EwaldSum_2d_;

void EwaldSum_2D_K0(_VSPME<_VQatom> &vspme, int nThreads, _EwaldSum_2d_::EwaldSumRecpK0Var &var, InteractRes &res) {
	MRunSys< ARRAY<_VQatom>, _VQatom, _EwaldSum_2d_::EwaldSumRecpK0Var, InteractRes>((void*)(&_EwaldSum_2d_::EwaldSum_2D_K0_Q), vspme.va, vspme.va.m, vspme.va.n, nThreads, var, res);
}

void EwaldSum_2D_K0(_VSPME<_VMUatom> &vspme, int nThreads, _EwaldSum_2d_::EwaldSumRecpK0Var &var, InteractRes &res) {
	MRunSys< ARRAY<_VMUatom>, _VMUatom, _EwaldSum_2d_::EwaldSumRecpK0Var, InteractRes>((void*)(&_EwaldSum_2d_::EwaldSum_2D_K0_MU), vspme.va, vspme.va.m, vspme.va.n, nThreads, var, res);
}


#if _SYS_ == _WINDOWS_SYS_
int ES_2DK0_Q(LPVOID *vp) {
#elif _SYS_ == _LINUX_SYS_
void* ES_2DK0_Q(void *vp) {
#endif
	EwaldSum_2D_K0_VARS_Q *par = (EwaldSum_2D_K0_VARS_Q*)vp;
	par->res->reset();

	EwaldSum_2D_K0(*(par->vspme), par->nThreads, *(par->var), *(par->res));
	
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->thread_var.hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void Asyn_EwaldSum_2D_K0_Q(_VSPME<_evatom_::_VQatom> &vspme, int nThreads, _EwaldSum_2d_::EwaldSumRecpK0Var &var, InteractRes &res, EwaldSum_2D_K0_VARS_Q &av) {
	av.thread_var.set_pars(0);
	av.vspme = &vspme;
	av.var = &var;
	av.res = &res;
	av.nThreads = nThreads;
#if _SYS_ == _WINDOWS_SYS_
	ResetEvent(av.thread_var.hMMThread);
	av.thread_var.thread = AfxBeginThread((AFX_THREADPROC)ES_2DK0_Q, (LPVOID)(&av), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
#elif _SYS_ == _LINUX_SYS_
	pthread_create(&(av.thread_var.thread), NULL, &ES_2DK0_Q, (void *)(&av));
#endif
}



#if _SYS_ == _WINDOWS_SYS_
int ES_2DK0_MU(LPVOID *vp) {
#elif _SYS_ == _LINUX_SYS_
void* ES_2DK0_MU(void *vp) {
#endif
	EwaldSum_2D_K0_VARS_MU *par = (EwaldSum_2D_K0_VARS_MU*)vp;
	par->res->reset();

	EwaldSum_2D_K0(*(par->vspme), par->nThreads, *(par->var), *(par->res));
	
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->thread_var.hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void Asyn_EwaldSum_2D_K0_MU(_VSPME<_evatom_::_VMUatom> &vspme, int nThreads, _EwaldSum_2d_::EwaldSumRecpK0Var &var, InteractRes &res, EwaldSum_2D_K0_VARS_MU &av) {
	av.thread_var.set_pars(0);
	av.vspme = &vspme;
	av.var = &var;
	av.res = &res;
	av.nThreads = nThreads;
#if _SYS_ == _WINDOWS_SYS_
	ResetEvent(av.thread_var.hMMThread);
	av.thread_var.thread = AfxBeginThread((AFX_THREADPROC)ES_2DK0_MU, (LPVOID)(&av), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
#elif _SYS_ == _LINUX_SYS_
	pthread_create(&(av.thread_var.thread), NULL, &ES_2DK0_MU, (void *)(&av));
#endif
}

/* multi-thread function to calculate real part of Ewald Sum for electrostatic interaction, to calculate Electric field and F on each atom */
void LocalEF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	InteractRes tres_real;

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
}

/* multi-thread function to calculate real part of Ewald Sum for electrostatic interaction, to calculate Electric field and F on each atom */
void LocalEF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	InteractRes tres_real;

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	_spme_::SPME_REAL_VARS svar;
	InteractRes tres_real;
	if (nThreads == 1) {
		if (mdcell->mm.m != NULL) {
			MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var.spme_3d_var, tres);
			res += tres;
		}
		if (mdcell->sm.m != NULL) {
			MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var.spme_3d_var, tres);
			res += tres;
		}
		if (mdcell->pm.m != NULL) {
			MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var.spme_3d_var, tres);
			res += tres;
		}

		// do we need this term if the system is not neutral ?
		//if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;

		//sprintf(errmsg, "Real : %f kT", res.U); show_log(errmsg, true);
	}
	else _spme_::Asyn_SPME_EF_Real(mdcell, cmm_cell, nThreads, var.spme_3d_var, tres_real, svar);

	//vspme.reset_local_EF(); // reset the variable for E & F

	EwaldSum_2D_K0_VARS_Q vk0;
	InteractRes tres_k0;
	if (nThreads == 1) {
		tres.reset();
		EwaldSum_2D_K0(vspme, nThreads, var.K0var, tres);
		res += tres;
		//sprintf(errmsg, "Recip [k=0]: %f kT", tres.U); show_log(errmsg, true);
	}
	else Asyn_EwaldSum_2D_K0_Q(vspme, nThreads, var.K0var, tres_k0, vk0);

	// accumulate the electric field from the image part
	init_normalized_coordinates(vspme, nThreads);
	cal_Q_spme(vspme, nThreads); // SPME Q matrix
	double U = vspme.cal_U_recp(!var.spme_3d_var.bEfieldOnly); //calculate the energy from the image part 
	res.U += U;
	cal_E_recp(vspme, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.cal_FTQ();
		vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {res.STxx += vspme.STxx; res.STyy += vspme.STyy; res.STzz += vspme.STzz;}
		res.STtr += vspme.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0]: %f kT", U); show_log(errmsg, true);

	// virial coming from the surface dipole, U(V), from receiprcal part
	if (!var.spme_3d_var.bEfieldOnly && var.bVirial) {
		EwaldSum_SurfaceDipole(mdcell->mu_surf, (EwaldSumRealVars0*)(&(var.spme_3d_var.esv1)), tres);
		res += tres;
	}
	

	//sprintf(errmsg, "After surface dipole correction, virial changes to %f", res.STtr); show_log(errmsg, true);

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VQatom> *pvcell = (_VE_CELL<_VQatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	if (nThreads > 1) {
		svar.thread_var.WaitUntilThreadOver();
		res += tres_real;
		vk0.thread_var.WaitUntilThreadOver();
		res += tres_k0;
	}

	bool bDumpF = !var.K0var.bEfieldOnly;
	if (nThreads < 2 || vspme.va.n < N4PARALLEL) vspme.dump_EF(bDumpF); // transfer the calculated E & F to the variables in real atom
	else MRunSys2< VSPME<_VQatom>, bool>((void*)(&Slab_dump_EF_Q), vspme, 0, vspme.va.n, nThreads, bDumpF);
}

void SPME_EF(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE*, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&VMM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE*, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&VSM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE*, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&VPM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}

	// do we need this term if the system is not neutral ?
	//if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;

	//sprintf(errmsg, "Real : %f, %f", res.U, res.STtr); show_log(errmsg, true);

	// accumulate the electric field from the image part
	//vspme.reset_local_EF(); // reset the variable for E & F
	init_normalized_coordinates(vspme, nThreads);
	cal_Q_spme(vspme, nThreads); // SPME Q matrix
	double U = vspme.cal_U_recp(!var.spme_3d_var.bEfieldOnly); //calculate the energy from the image part 
	res.U += U;
	cal_E_recp(vspme, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.cal_FTQ();
		vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {res.STxx += vspme.STxx; res.STyy += vspme.STyy; res.STzz += vspme.STzz;}
		res.STtr += vspme.STtr;
	}
	//sprintf(errmsg, "Recip : %f, %f", U, vspme.STtr); show_log(errmsg, true);

	tres.reset();
	EwaldSum_2D_K0(vspme, nThreads, var.K0var, tres);
	res += tres;

	// virial coming from the surface dipole, U(V), from receiprcal part
	if (!var.spme_3d_var.bEfieldOnly && var.bVirial) {
		EwaldSum_SurfaceDipole(mdcell->mu_surf, (EwaldSumRealVars0*)(&(var.spme_3d_var.esv1)), tres);
		res += tres;
	}

	//sprintf(errmsg, "After surface dipole correction, virial changes to %f", res.STtr); show_log(errmsg, true);

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VQatom> *pvcell = (_VE_CELL<_VQatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	//vspme.dump_EF(!var.K0var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
	bool bDumpF = !var.K0var.bEfieldOnly;
	if (nThreads < 2 || vspme.va.n < N4PARALLEL) vspme.dump_EF(bDumpF); // transfer the calculated E & F to the variables in real atom
	else MRunSys2< VSPME<_VQatom>, bool>((void*)(&Slab_dump_EF_Q), vspme, 0, vspme.va.n, nThreads, bDumpF);
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
void SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res);

	if (var.K0var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}

void SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res);

	if (var.K0var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE*>(VMM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE*>(VSM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE*>(VPM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}

// with dipole
/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

#if _CHECK_TIME_STAMP_
	TIME mt;
	mt.start();
	show_log("  ||  ", true);
#endif

	_spme_::SPME_REAL_VARS svar;
	InteractRes tres_real;
	if (nThreads == 1) {
		if (mdcell->mm.m != NULL) {
			MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var.spme_3d_var, tres);
			res += tres;
		}
		if (mdcell->sm.m != NULL) {
			MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var.spme_3d_var, tres);
			res += tres;
		}
		if (mdcell->pm.m != NULL) {
			MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var.spme_3d_var, tres);
			res += tres;
		}

		// do we need this term if the system is not neutral ?
		//if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;
	}
	else _spme_::Asyn_SPME_EF_Real(mdcell, cmm_cell, nThreads, var.spme_3d_var, tres_real, svar);

	//vspme.reset_local_EF(); // reset the variable for E & F

	EwaldSum_2D_K0_VARS_MU vk0;
	InteractRes tres_k0;
	if (nThreads == 1) {
		tres.reset();
		EwaldSum_2D_K0(vspme, nThreads, var.K0var, tres);
		res += tres;
	}
	else Asyn_EwaldSum_2D_K0_MU(vspme, nThreads, var.K0var, tres_k0, vk0);


#if _CHECK_TIME_STAMP_
	mt.elapse(true, 1); show_log(", ", false);
#endif
	// accumulate the electric field from the image part
#if _CHECK_TIME_STAMP_
	mt.elapse(true, 2); show_log(", ", false);
#endif

	//init_normalized_coordinates(vspme, nThreads); -- vspme is initilized outside before this funciton is called
	cal_Q_spme(vspme, nThreads); // SPME Q matrix
#if _CHECK_TIME_STAMP_
	mt.elapse(true, 3); show_log(", ", false);
#endif

	double U = vspme.cal_U_recp(!var.spme_3d_var.bEfieldOnly); //calculate the energy from the image part 
	res.U += U;
#if _CHECK_TIME_STAMP_
	mt.elapse(true, 4); show_log(", ", false);
#endif

	cal_E_recp(vspme, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

#if _CHECK_TIME_STAMP_
	mt.elapse(true, 5); show_log(", ", false);
#endif

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_StressTensor(vspme, nThreads);
		if (var.bST_diagonal) {
			res.STxx += vspme.STxx; res.STyy += vspme.STyy; res.STzz += vspme.STzz;
		}
		res.STtr += vspme.STtr;
	}

#if _CHECK_TIME_STAMP_
	mt.elapse(true, 6); show_log(", ", false);
#endif

	// virial coming from the surface dipole, U(V), from receiprical part
	if (!var.spme_3d_var.bEfieldOnly && var.bVirial) {
		EwaldSum_SurfaceDipole(mdcell->mu_surf, (EwaldSumRealVars0*)(&(var.spme_3d_var.esv1)), tres);
		res += tres;
	}

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VMUatom> *pvcell = (_VE_CELL<_VMUatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	if (nThreads > 1) {
		svar.thread_var.WaitUntilThreadOver();
		res += tres_real;
		vk0.thread_var.WaitUntilThreadOver();
		res += tres_k0;
	}

	//vspme.dump_EF(!var.K0var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
	bool bDumpF = !var.K0var.bEfieldOnly;
	if (nThreads < 2 || vspme.va.n < N4PARALLEL) vspme.dump_EF(bDumpF); // transfer the calculated E & F to the variables in real atom
	else MRunSys2< VSPME<_VMUatom>, bool>((void*)(&Slab_dump_EF_MU), vspme, 0, vspme.va.n, nThreads, bDumpF);
}

void SPME_EF(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE*, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&VMM_EFIELD_SPME_REAL), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE*, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&VSM_EFIELD_SPME_REAL), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE*, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&VPM_EFIELD_SPME_REAL), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}

	// do we need this term if the system is not neutral ?
	//if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;

	// accumulate the electric field from the image part
	//vspme.reset_local_EF(); // reset the variable for E & F
	//init_normalized_coordinates(vspme, nThreads); -- vspme is initilized outside before this funciton is called
	cal_Q_spme(vspme, nThreads); // SPME Q matrix
	double U = vspme.cal_U_recp(!var.spme_3d_var.bEfieldOnly); //calculate the energy from the image part 
	res.U += U;
	cal_E_recp(vspme, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_StressTensor(vspme, nThreads);
		if (var.bST_diagonal) {
			res.STxx += vspme.STxx; res.STyy += vspme.STyy; res.STzz += vspme.STzz;
		}
		res.STtr += vspme.STtr;
	}

	tres.reset();
	EwaldSum_2D_K0(vspme, nThreads, var.K0var, tres);
	res += tres;

	// virial coming from the surface dipole, U(V), from receiprical part
	if (!var.spme_3d_var.bEfieldOnly && var.bVirial) {
		EwaldSum_SurfaceDipole(mdcell->mu_surf, (EwaldSumRealVars0*)(&(var.spme_3d_var.esv1)), tres);
		res += tres;
	}

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VMUatom> *pvcell = (_VE_CELL<_VMUatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	//vspme.dump_EF(!var.K0var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
	bool bDumpF = !var.K0var.bEfieldOnly;
	if (nThreads < 2 || vspme.va.n < N4PARALLEL) vspme.dump_EF(bDumpF); // transfer the calculated E & F to the variables in real atom
	else MRunSys2< VSPME<_VMUatom>, bool>((void*)(&Slab_dump_EF_MU), vspme, 0, vspme.va.n, nThreads, bDumpF);
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
void SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	init_normalized_coordinates(vspme, nThreads); // vspme has to be initialized before SPME_EF for atom
	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res);

	if (var.K0var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}

void SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	init_normalized_coordinates(vspme, nThreads); // vspme has to be initialized before SPME_EF for atom
	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res);

	if (var.K0var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}

void dmax(double &x1, double &x2, double &res) {
	res = (x1 > x2 ? x1 : x2);
}
/*
void AtomicEFieldCorrect_SM(SMOLECULE *sm, int nSM) {
	EwaldSumRealVars1 esv1;
	esv1.bEwald = false; // direct interaction
	VECTOR3 E, dr;
	double R, R2;
	int im, ia, ia1;
	BASIC_CLUSTER *pc;
	MATOM *patom, *patom1;
	for (im = 0; im < nSM; im++) {
		pc = sm[im].c;
		for (ia = 0; ia < pc->nAtoms; ia++) {
			patom = pc->atom + ia;
			for (ia1 = 0; ia1 < pc->nAtoms; ia1++) {
				if (ia == ia1) continue;
				patom1 = pc->atom + ia1;
				VECT3(patom->r, patom1->r, dr)
				V3ABS2(dr, R2) R = sqrt(R2);
				esv1.init_vars(dr, R, R2);
				MTP_E(esv1, patom1->c, patom1->mu.v, E, false, 0);
				patom->E.v[0] += E.v[0];
				patom->E.v[1] += E.v[1];
				patom->E.v[2] += E.v[2];
			}
		}
	}
}
*/
/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
bool Polarize_SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field

	bool bEfieldOnly = var.K0var.bEfieldOnly;
	res.reset();
	InteractRes tres;

	init_normalized_coordinates(vspme, nThreads); // vspme has to be initialized before SPME_EF for atom with dipole
	double dmu_max = 0, vt = 0, muConvg = ::muConvg; // eA
	int iloop = 0, max_loops = ::nmax_loops_polarize;
	bool status = false;
	var.K0var.bEfieldOnly = true;
	var.spme_3d_var.bEfieldOnly = true;
	while (bPolarize && iloop < max_loops) {
		SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, tres); // calculate the dipole only

		// add the electric field from the atom inside the cluster
		//if (mdcell->sm.n > 0) MOperate<SMOLECULE>((void*)(&AtomicEFieldCorrect_SM), mdcell->sm.m, mdcell->sm.n, nThreads);

		// polarization
		if (mdcell->mm.n > 0) MOperate<MMOLECULE>((void*)(&polarize_MM), mdcell->mm.m, mdcell->mm.n, nThreads);
		if (mdcell->sm.n > 0) MOperate<SMOLECULE>((void*)(&polarize_SM), mdcell->sm.m, mdcell->sm.n, nThreads);
		if (mdcell->pm.n > 0) MOperate<PMOLECULE>((void*)(&polarize_PM), mdcell->pm.m, mdcell->pm.n, nThreads);
		cal_dipole(*mdcell, nThreads, mdcell->mu_surf); // accumulate dipole of each atom, and the whole mdcell
		memcpy(&(cmm_cell->mu_surf), &(mdcell->mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
		dmu_max = 0; vt = 0;
		if (mdcell->ind_mu_mm.n > 0) MOperate2<MOL_POL<MMOLECULE>, double>((void*)(&MM_maxdiff_backup_polarization), mdcell->ind_mu_mm.m, mdcell->ind_mu_mm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_sm.n > 0) MOperate2<MOL_POL<SMOLECULE>, double>((void*)(&SM_maxdiff_backup_polarization), mdcell->ind_mu_sm.m, mdcell->ind_mu_sm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_pm.n > 0) MOperate2<MOL_POL<PMOLECULE>, double>((void*)(&PM_maxdiff_backup_polarization), mdcell->ind_mu_pm.m, mdcell->ind_mu_pm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;

		iloop++;
		if (dmu_max <= muConvg && iloop > 0) {status = true; break;}
	}

	//char msg[256] = "\0";
	//if (status) sprintf(msg, "Convergent polarized dipole is obtainted with %d loops", iloop); 
	//else sprintf(msg, "No convergent polarized dipole is obtained with %d loops", iloop);
	//show_log(msg, true);

	var.K0var.bEfieldOnly = bEfieldOnly;
	var.spme_3d_var.bEfieldOnly = bEfieldOnly;

	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res); // calculate the dipole and force

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}

	return status;
}

bool Polarize_SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field

	bool bEfieldOnly = var.K0var.bEfieldOnly;
	res.reset();
	InteractRes tres;

	init_normalized_coordinates(vspme, nThreads); // vspme has to be initialized before SPME_EF for atom with dipole
	double dmu_max = 0, vt = 0, muConvg = ::muConvg; // eA
	int iloop = 0, max_loops = ::nmax_loops_polarize;
	bool status = false;
	var.K0var.bEfieldOnly = true;
	var.spme_3d_var.bEfieldOnly = true;
	while (bPolarize && iloop < max_loops) {
		SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, tres); // calculate the dipole only

		// add the electric field from the atom inside the cluster
		//if (mdcell->sm.n > 0) MOperate<SMOLECULE>((void*)(&AtomicEFieldCorrect_SM), mdcell->sm.m, mdcell->sm.n, nThreads);

		// polarization
		if (mdcell->mm.n > 0) MOperate<MMOLECULE*>((void*)(&polarize_VMM), mdcell->mm.m, mdcell->mm.n, nThreads);
		if (mdcell->sm.n > 0) MOperate<SMOLECULE*>((void*)(&polarize_VSM), mdcell->sm.m, mdcell->sm.n, nThreads);
		if (mdcell->pm.n > 0) MOperate<PMOLECULE*>((void*)(&polarize_VPM), mdcell->pm.m, mdcell->pm.n, nThreads);
		cal_dipole(*mdcell, nThreads, mdcell->mu_surf); // accumulate dipole of each atom, and the whole mdcell
		memcpy(&(cmm_cell->mu_surf), &(mdcell->mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
		dmu_max = 0; vt = 0;
		if (mdcell->ind_mu_mm.n > 0) MOperate2<MOL_POL<MMOLECULE>, double>((void*)(&MM_maxdiff_backup_polarization), mdcell->ind_mu_mm.m, mdcell->ind_mu_mm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_sm.n > 0) MOperate2<MOL_POL<SMOLECULE>, double>((void*)(&SM_maxdiff_backup_polarization), mdcell->ind_mu_sm.m, mdcell->ind_mu_sm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_pm.n > 0) MOperate2<MOL_POL<PMOLECULE>, double>((void*)(&PM_maxdiff_backup_polarization), mdcell->ind_mu_pm.m, mdcell->ind_mu_pm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;

		iloop++;
		if (dmu_max <= muConvg && iloop > 0) {status = true; break;}
	}

	//char msg[256] = "\0";
	//if (status) sprintf(msg, "Convergent polarized dipole is obtainted with %d loops", iloop); 
	//else sprintf(msg, "No convergent polarized dipole is obtained with %d loops", iloop);
	//show_log(msg, true);

	var.K0var.bEfieldOnly = bEfieldOnly;
	var.spme_3d_var.bEfieldOnly = bEfieldOnly;

	SPME_EF(mdcell, cmm_cell, vspme, nThreads, var, res); // calculate the dipole and force

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}

	return status;
}



// polarized SPME, with the interaction between induced dipole excluded
// in some polarized model, the interaction between induced dipole is excluded, e.g. only q-induced dipole is considered.
void VSPME_F_exclude_induced(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme_induced, int nThreads, SPME_VAR &var, InteractRes &res, bool bInitABSP_buff) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	// accumulate the electric field from the image part
	//vspme_induced.reset_local_EF(); // reset the variable for E & F

	// dipole-dipole interaction should be excluded before, so disable this calculation here
	/*
	EwaldSum_2D_K0_VARS_MU vk0;
	InteractRes tres_k0;
	if (nThreads == 1) {
		tres.reset();
		EwaldSum_2D_K0(vspme_induced, nThreads, var.K0var, tres);
		res += tres;
	}
	else Asyn_EwaldSum_2D_K0_MU(vspme_induced, nThreads, var.K0var, tres_k0, vk0);
	*/

	init_normalized_coordinates(vspme_induced, nThreads); //-- vspme is initilized outside before this funciton is called
	if (!bInitABSP_buff) vspme_induced.bReady_absp = true; // copied from vspme
	cal_Q_spme(vspme_induced, nThreads); // SPME Q matrix
	double U = vspme_induced.cal_U_recp(!var.spme_3d_var.bEfieldOnly); //calculate the energy from the image part 
	res.U = U;
	cal_E_recp(vspme_induced, true, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (var.bVirial) {
		cal_StressTensor(vspme_induced, nThreads);
		if (var.bST_diagonal) {
			res.STxx += vspme_induced.STxx; res.STyy += vspme_induced.STyy; res.STzz += vspme_induced.STzz;
		}
		res.STtr += vspme_induced.STtr;
	}

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VMUatom> *pvcell = (_VE_CELL<_VMUatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/
	/*
	if (nThreads > 1) {
		vk0.thread_var.WaitUntilThreadOver();
		res += tres_k0;
	}
	*/

	//vspme.dump_EF(!var.K0var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
	//bool bDumpF = !var.K0var.bEfieldOnly;
	bool bDumpF = true; // useless
	if (nThreads < 2 || vspme_induced.va.n < N4PARALLEL) Slab_dump_F_exclude_induced_mu(vspme_induced, 0, vspme_induced.va.n - 1, bDumpF); // transfer the calculated E & F to the variables in real atom
	else MRunSys2< VSPME<_VMUatom>, bool >((void*)(&Slab_dump_F_exclude_induced_mu), vspme_induced, 0, vspme_induced.va.n, nThreads, bDumpF);
}

bool Polarize_SPME_Interact2(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, VSPME<_VMUatom>& vspme_induced, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
	vspme_induced.release_plans();
	vspme.init_fft_plans();
#if _CHECK_TIME_STAMP_
	show_log("Polarize: ", false);
	show_vlog<int>(mt.elapse(), false);
#endif
	bool status = Polarize_SPME_Interact(mdcell, cmm_cell, vspme, nThreads, bPolarize, var, res);
	InteractRes tres;
	bool bABSP_Buff = false;
	if (bPolarize && bExcludeInducedDipole) {
		vspme.release_plans();
		vspme_induced.init_fft_plans();
		if (vspme.absp.n == vspme_induced.absp.n) {
			cp_absp_buff((__AtomBSP_BUFF*)(&vspme), (__AtomBSP_BUFF*)(&vspme_induced));
			bABSP_Buff = true;
		}
		else bABSP_Buff = false;
		VSPME_F_exclude_induced(mdcell, cmm_cell, vspme_induced, nThreads, var, tres, !bABSP_Buff);
		res -= tres;
	}
#if _CHECK_TIME_STAMP_
	show_log("", true);
#endif
	return status;
}

// for 2d system to calculate electrostatic interaction with explicitly
void ExplicitEInteract_2d(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, E2dInteractVar& evar, InteractRes &res) {
	_EwaldSum_real_::EwaldSumRealVars1 *esv = &(evar.esv1);
	esv->bEwald = false; // direct interaction
	MATOM *patom = NULL;
	VECTOR3 E, F, dr;//, E_Thole, F_Thole;
	int i, j, k, l, na, na1;
	VECTOR3 r0, r1;

	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;

	double R, R2, U;
	int npc[2], nCell, NC0 = 64, NC = 128;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, intrinsic_cell = true;

	double Uestat = 0;//, U_Thole = 0;

	/*
	char aindx[2] = {0x00, 0x00};
	unsigned int key = 0;
	PAIRWISE<float> *TholeRadius = NULL;
	double a_Thole = 0;
	*/

	res.reset();

	//for (na = 0; na < pc->nAtoms; na++) {
		//patom = pc->atom + na;
	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];

		V3zero(patom->E) V3zero(patom->F)

		V32V3(patom->rg, r0)

		// periodic cell iteration
		for (npc[0] = -NC; npc[0] <= NC; npc[0]++) {drCell.v[0] = npc[0] * cell.xd[0]; 
		for (npc[1] = -NC; npc[1] <= NC; npc[1]++) {drCell.v[1] = npc[1] * cell.xd[1]; 
			intrinsic_cell = ((npc[0] == 0 && npc[1] == 0) ? true : false);
			// image cell which is far from the cluster
			if (FABS(npc[0]) < NC0 && FABS(npc[1]) < NC0) { // calculate explicitly for each cluster
			for (nCell = 0; nCell < cell.acell.n; nCell++) {
				pcell = cell.acell.m[nCell]; it.set(pcell); it.start_iterate();
				while (!it.eoi()) {
					pc1 = it.current(); if (pc1 == NULL) {it.next_iterate(); continue;}
					if (intrinsic_cell && pc1 == pc) { // same cluster in same cell
						it.next_iterate();
						continue;
					}
					//for (na1 = 0; na1 < pc1->nAtoms; na1++) {
						//patom1 = pc1->atom + na1;
					for (na1 = 0; na1 < pc1->fatom.n; na1++) {
						patom1 = pc1->fatom.m[na1];
						V32V3(patom1->rg, r1)
						for (i = 0; i < 2; i++) dr.v[i] = r0.v[i] - r1.v[i] - drCell.v[i];
						dr.v[2] = r0.v[2] - r1.v[2];
						V3ABS2(dr, R2)
						R = sqrt(R2);

						if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
							BOUND(patom, patom1, i, bound) // bound atom
							if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
							if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
							else bound14 = false;
						}
						else {bound = false; bound14 = false;}

						if (bound || (pc == pc1 && intrinsic_cell)) { // bounding atoms, or inside the same cluster
							continue;
						}

						if (patom->mMP == 0) {
							esv->EwaldSumRealVars0::init_vars(dr, R, R2);
							MTP_Interact(*esv, patom->c, patom1->c, E, F, U, false, 0);
						}
						else if (patom->mMP == 1) {
							esv->init_vars(dr, R, R2);
							MTP_Interact(*esv, patom->c, patom->mu.v, patom1->c, patom1->mu.v, E, F, U, false, 0);
						}

						if (bound14 && pc != pc1) { // 1-4 interaction
							U *= A14_E;
							// electric field
							for (i = 0; i < 3; i++) {
								E.v[i] *= A14_E; F.v[i] *= A14_E;
							}
						}

						/*
						if (evar.bTholePolarize) {
							aindx[0] = patom->par->aindx; aindx[1] = patom1->par->aindx; key = construct_char2_key_pairwise(aindx);
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
						*/

						Uestat += U;
						for (i = 0; i < 3; i++) {
							patom->E.v[i] += E.v[i];
							patom->F.v[i] += F.v[i] * ::fUnit_estat_mass;
						}
						if (evar.bVirial) Accumulate_Virial_Diagonal(res, F.v, dr.v);
					}

					it.next_iterate();
				}
			}
			}
			else {
				r1.v[0] = 0; r1.v[1] = 0; r1.v[2] = 0; // center of whole CELL
				for (i = 0; i < 2; i++) dr.v[i] = r0.v[i] - r1.v[i] - drCell.v[i];
				dr.v[2] = r0.v[2] - r1.v[2];
				V3ABS2(dr, R2)
				R = sqrt(R2);
				esv->init_vars(dr, R, R2);
				MTP_Interact(*esv, patom->c, patom->mu.v, cell.mu_surf.q, cell.mu_surf.mu.v, E, F, U, false, 0);
				Uestat += U;
				for (i = 0; i < 3; i++) {
					patom->E.v[i] += E.v[i];
					patom->F.v[i] += F.v[i] * ::fUnit_estat_mass;
				}
				if (evar.bVirial) Accumulate_Virial_Diagonal(res, F.v, dr.v);
			}
		}} // end of periodic cell iteration

		//V3PV3(patom->r0, force, torque);
		//for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
	}

	// total energy in kT
	res.U = Uestat * ::eUnit_estat_kT;
	if (evar.bVirial) {
		if (evar.bST_diagonal) {
			res.STxx *= fUnit_estat_mass;
			res.STyy *= fUnit_estat_mass;
			res.STzz *= fUnit_estat_mass;
		}
		res.STtr *= fUnit_estat_mass;
	}
}

void ExplicitEInteract_2d_U(CMM_CELL3D<BASIC_CLUSTER> &cell, E2dInteractVar& evar, InteractRes &res) {
	_EwaldSum_real_::EwaldSumRealVars1 *esv = &(evar.esv1);
	esv->bEwald = false; // direct interaction
	VECTOR3 E, F, dr;//, E_Thole, F_Thole;
	int i, j, k, l;
	VECTOR3 r0, r1;

	double R, R2, U;
	int npc[2], nCell, NC = 360, NC0 = 64;
	VECTOR3 drCell;
	double Uestat = 0;

	// periodic cell iteration
	for (npc[0] = -NC; npc[0] <= NC; npc[0]++) {drCell.v[0] = npc[0] * cell.xd[0]; 
	for (npc[1] = -NC; npc[1] <= NC; npc[1]++) {drCell.v[1] = npc[1] * cell.xd[1]; 
		if (FABS(npc[0]) <= NC0 && FABS(npc[1]) <= NC0) { // calculate explicitly for each cluster
		}
		else {
			V32V3(drCell, dr);
			V3ABS2(dr, R2)
			R = sqrt(R2);
			esv->init_vars(dr, R, R2);
			MTP_Interact(*esv, cell.mu_surf.q, cell.mu_surf.mu.v, cell.mu_surf.q, cell.mu_surf.mu.v, E, F, U, false, 0);
			Uestat += U;
		}
	}}
	res.U = Uestat * ::eUnit_estat_kT;
}

double ShapeEnergy_2d(CMM_CELL3D<BASIC_CLUSTER> &cell) {
	double U2 = 0, R2 = 0;
	// sphere geometry
	V3ABS2(cell.mu_surf.mu, R2)
	double Uestat1 = 2 * PI / (3 * cell.xd[0] * cell.xd[1] * cell.xd[2]) * R2; 
	// infinit thin geometry
	R2 = cell.mu_surf.mu.v[2] * cell.mu_surf.mu.v[2];
	double Uestat2 = 2 * PI / (cell.xd[0] * cell.xd[1] * cell.xd[2]) * R2; 

	double Uestat = Uestat1 - Uestat2;
	return Uestat * ::eUnit_estat_kT;
}

void CMM_Explicit_EInteract_2d(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, E2dInteractVar& evar, InteractRes &res) {
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;
	InteractRes tres;
	for (nc = n1Cell; nc <= n2Cell; nc++) {
		pcell = cell.acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			ExplicitEInteract_2d(cell, pc, evar, tres);
			res += tres;
			it.next_iterate();
		}
	}
}






// for 2d system to calculate electrostatic interaction with explicitly, excluding the interaction between induced dipoles
void ExplicitEInteract_2d_ExcludeInducedDiploe(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, E2dInteractVar& evar, InteractRes &res) {
	_EwaldSum_real_::EwaldSumRealVars1 *esv = &(evar.esv1);
	esv->bEwald = false; // direct interaction
	MATOM *patom = NULL;
	VECTOR3 E, F, dr;//, E_Thole, F_Thole;
	VECTOR3 Ei, Fi;
	int i, j, k, l, na, na1;
	VECTOR3 r0, r1;

	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;

	double R, R2, U, Ui;
	int npc[2], nCell, NC0 = 64, NC = 128;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, intrinsic_cell = true;

	double Uestat = 0;//, U_Thole = 0;

	/*
	char aindx[2] = {0x00, 0x00};
	unsigned int key = 0;
	PAIRWISE<float> *TholeRadius = NULL;
	double a_Thole = 0;
	*/

	res.reset();

	//for (na = 0; na < pc->nAtoms; na++) {
		//patom = pc->atom + na;
	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];

		V3zero(patom->E) V3zero(patom->F)

		V32V3(patom->rg, r0)

		// periodic cell iteration
		for (npc[0] = -NC; npc[0] <= NC; npc[0]++) {drCell.v[0] = npc[0] * cell.xd[0]; 
		for (npc[1] = -NC; npc[1] <= NC; npc[1]++) {drCell.v[1] = npc[1] * cell.xd[1]; 
			intrinsic_cell = ((npc[0] == 0 && npc[1] == 0) ? true : false);
			// image cell which is far from the cluster
			if (FABS(npc[0]) < NC0 && FABS(npc[1]) < NC0) { // calculate explicitly for each cluster
			for (nCell = 0; nCell < cell.acell.n; nCell++) {
				pcell = cell.acell.m[nCell]; it.set(pcell); it.start_iterate();
				while (!it.eoi()) {
					pc1 = it.current(); if (pc1 == NULL) {it.next_iterate(); continue;}
					if (intrinsic_cell && pc1 == pc) { // same cluster in same cell
						it.next_iterate();
						continue;
					}
					//for (na1 = 0; na1 < pc1->nAtoms; na1++) {
						//patom1 = pc1->atom + na1;
					for (na1 = 0; na1 < pc1->fatom.n; na1++) {
						patom1 = pc1->fatom.m[na1];
						V32V3(patom1->rg, r1)
						for (i = 0; i < 2; i++) dr.v[i] = r0.v[i] - r1.v[i] - drCell.v[i];
						dr.v[2] = r0.v[2] - r1.v[2];
						V3ABS2(dr, R2)
						R = sqrt(R2);

						if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
							BOUND(patom, patom1, i, bound) // bound atom
							if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
							if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
							else bound14 = false;
						}
						else {bound = false; bound14 = false;}

						if (bound || (pc == pc1 && intrinsic_cell)) { // bounding atoms, or inside the same cluster
							continue;
						}

						if (patom->mMP == 0) {
							esv->EwaldSumRealVars0::init_vars(dr, R, R2);
							MTP_Interact(*esv, patom->c, patom1->c, E, F, U, false, 0);
						}
						else if (patom->mMP == 1) {
							esv->init_vars(dr, R, R2);
							MTP_Interact(*esv, patom->c, patom->mu.v, patom1->c, patom1->mu.v, E, F, U, false, 0);

							MTP_Interact(*esv, 0, patom->ind_mu.v, 0, patom1->ind_mu.v, Ei, Fi, Ui, false, 0);
							//V3minusV3(E, Ei, E) 
							V3minusV3(F, Fi, F)
							U -= Ui;
						}

						if (bound14 && pc != pc1) { // 1-4 interaction
							U *= A14_E;
							// electric field
							for (i = 0; i < 3; i++) {
								E.v[i] *= A14_E; F.v[i] *= A14_E;
							}
						}

						/*
						if (evar.bTholePolarize) {
							aindx[0] = patom->par->aindx; aindx[1] = patom1->par->aindx; key = construct_char2_key_pairwise(aindx);
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
						*/

						Uestat += U;
						for (i = 0; i < 3; i++) {
							patom->E.v[i] += E.v[i];
							patom->F.v[i] += F.v[i] * ::fUnit_estat_mass;
						}
						if (evar.bVirial) Accumulate_Virial_Diagonal(res, F.v, dr.v);
					}

					it.next_iterate();
				}
			}
			}
			else {
				r1.v[0] = 0; r1.v[1] = 0; r1.v[2] = 0; // center of whole CELL
				for (i = 0; i < 2; i++) dr.v[i] = r0.v[i] - r1.v[i] - drCell.v[i];
				dr.v[2] = r0.v[2] - r1.v[2];
				V3ABS2(dr, R2)
				R = sqrt(R2);
				esv->init_vars(dr, R, R2);
				MTP_Interact(*esv, patom->c, patom->mu.v, cell.mu_surf.q, cell.mu_surf.mu.v, E, F, U, false, 0);
				MTP_Interact(*esv, 0, patom->ind_mu.v, 0, cell.mu_surf.mu_ind.v, Ei, Fi, Ui, false, 0);
				//V3minusV3(E, Ei, E)
				V3minusV3(F, Fi, F)
				U -= Ui;
				Uestat += U;
				for (i = 0; i < 3; i++) {
					patom->E.v[i] += E.v[i];
					patom->F.v[i] += F.v[i] * ::fUnit_estat_mass;
				}
				if (evar.bVirial) Accumulate_Virial_Diagonal(res, F.v, dr.v);
			}
		}} // end of periodic cell iteration

		//V3PV3(patom->r0, force, torque);
		//for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
	}

	// total energy in kT
	res.U = Uestat * ::eUnit_estat_kT;
	if (evar.bVirial) {
		if (evar.bST_diagonal) {
			res.STxx *= fUnit_estat_mass;
			res.STyy *= fUnit_estat_mass;
			res.STzz *= fUnit_estat_mass;
		}
		res.STtr *= fUnit_estat_mass;
	}
}



// for 2d system to calculate electrostatic interaction with explicitly
// cell includes the image cluster inside
void ExplicitEInteract_2d_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, E2dInteractVar& evar, InteractRes &res) {
	_EwaldSum_real_::EwaldSumRealVars1 *esv = &(evar.esv1);
	esv->bEwald = false; // direct interaction
	MATOM *patom = NULL;
	VECTOR3 E, F, dr;//, E_Thole, F_Thole;
	int i, j, k, l, na, na1;
	VECTOR3 r0, r1;

	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;

	double R, R2, U;
	int npc[2], nCell, NC = 16;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, intrinsic_cell = true;

	double Uestat = 0;//, U_Thole = 0;

	/*
	char aindx[2] = {0x00, 0x00};
	unsigned int key = 0;
	PAIRWISE<float> *TholeRadius = NULL;
	double a_Thole = 0;
	*/

	res.reset();

	//for (na = 0; na < pc->nAtoms; na++) {
		//patom = pc->atom + na;
	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];

		V32V3(patom->rg, r0)

		// periodic cell iteration
		for (npc[0] = -NC; npc[0] <= NC; npc[0]++) {drCell.v[0] = npc[0] * cell.xd[0]; 
		for (npc[1] = -NC; npc[1] <= NC; npc[1]++) {drCell.v[1] = npc[1] * cell.xd[1]; 
			intrinsic_cell = ((npc[0] == 0 && npc[1] == 0) ? true : false);
			// image cell which is far from the cluster
			for (nCell = 0; nCell < cell.acell.n; nCell++) {
				pcell = cell.acell.m[nCell]; it.set(pcell); it.start_iterate();
				while (!it.eoi()) {
					pc1 = it.current(); if (pc1 == NULL) {it.next_iterate(); continue;}
					if (intrinsic_cell && pc1 == pc) { // same cluster in same cell
						it.next_iterate();
						continue;
					}
					//for (na1 = 0; na1 < pc1->nAtoms; na1++) {
						//patom1 = pc1->atom + na1;
					for (na1 = 0; na1 < pc1->fatom.n; na1++) {
						patom1 = pc1->fatom.m[na1];
						V32V3(patom1->rg, r1)
						for (i = 0; i < 2; i++) dr.v[i] = r0.v[i] - r1.v[i] - drCell.v[i];
						dr.v[2] = r0.v[2] - r1.v[2];
						V3ABS2(dr, R2)
						R = sqrt(R2);

						if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
							BOUND(patom, patom1, i, bound) // bound atom
							if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
							if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
							else bound14 = false;
						}
						else {bound = false; bound14 = false;}

						if (bound || (pc == pc1 && intrinsic_cell)) { // bounding atoms, or inside the same cluster
							continue;
						}

						if (patom->mMP == 0) {
							esv->EwaldSumRealVars0::init_vars(dr, R, R2);
							MTP_Interact(*esv, patom->c, patom1->c, E, F, U, false, 0);
						}
						else if (patom->mMP == 1) {
							esv->init_vars(dr, R, R2);
							MTP_Interact(*esv, patom->c, patom->mu.v, patom1->c, patom1->mu.v, E, F, U, false, 0);
						}

						if (bound14 && pc != pc1) { // 1-4 interaction
							U *= A14_E;
							// electric field
							for (i = 0; i < 3; i++) {
								E.v[i] *= A14_E; F.v[i] *= A14_E;
							}
						}

						/*
						if (evar.bTholePolarize) {
							aindx[0] = patom->par->aindx; aindx[1] = patom1->par->aindx; key = construct_char2_key_pairwise(aindx);
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
						*/

						// pc is assumed to be real cluster in slab
						if (pc1->cID >= 0) {// pc1 is real cluster also
							Uestat += U;
							for (i = 0; i < 3; i++) {
								patom->E.v[i] += E.v[i];
								patom->F.v[i] += F.v[i] * ::fUnit_estat_mass;
							}
							if (evar.bVirial) Accumulate_Virial_Diagonal(res, F.v, dr.v);
						}
						else { // pc1 is image cluster
							Uestat += U * 2; // MTP function return U * 0.5
							// this interaction is A with B_i, also it can be used for B with A_i with force on A
							// because A_i's coordinate is associate with A, so E & F are is on A also
							for (i = 0; i < 3; i++) {
								patom->E.v[i] += E.v[i] * 2;
								patom->F.v[i] += F.v[i] * ::fUnit_estat_mass * 2;
							}
							if (evar.bVirial) Accumulate_Virial_Diagonal(res, F.v, dr.v);
						}
					}

					it.next_iterate();
				}
			}
		}} // end of periodic cell iteration

		//V3PV3(patom->r0, force, torque);
		//for (i = 0; i < 3; i++) {pc->dyn->fc.v[i] += torque.v[i]; pc->dyn->fc.v[i+3] += force.v[i];}
	}
	// total energy in kT
	res.U = Uestat * ::eUnit_estat_kT;
	if (evar.bVirial) {
		if (evar.bST_diagonal) {
			res.STxx *= fUnit_estat_mass;
			res.STyy *= fUnit_estat_mass;
			res.STzz *= fUnit_estat_mass;
		}
		res.STtr *= fUnit_estat_mass;
	}
}

#endif // _IGNORE_ELECTROSTATIC_INTERACTION_

} // end of namespace _spme_

#define TEST 1
#if TEST == 1
using namespace _spme_2d_;

void get_atoms(MMOL_MD_CELL &mdcell, ARRAY<MATOM*>& atom) {
	int im, ic, ia, n;
	atom.set_array(nAtoms(mdcell));
	n = 0;
	for (im = 0; im < mdcell.mm.n; im++) {
		for (ic = 0; ic < mdcell.mm.m[im].nCluster; ic++) {
			for (ia = 0; ia < mdcell.mm.m[im].cluster[ic].nAtoms; ia++) {
				atom.m[n] = mdcell.mm.m[im].cluster[ic].atom + ia;
				atom.m[n]->aIndx = n;
				n++;
			}
		}
	}
	for (im = 0; im < mdcell.sm.n; im++) {
		for (ia = 0; ia < mdcell.sm.m[im].c->nAtoms; ia++) {
			atom.m[n] = mdcell.sm.m[im].c->atom + ia;
			atom.m[n]->aIndx = n;
			n++;
		}
	}
	for (im = 0; im < mdcell.pm.n; im++) {
		for (ia = 0; ia < mdcell.pm.m[im].c->nAtoms; ia++) {
			atom.m[n] = mdcell.pm.m[im].c->atom + ia;
			atom.m[n]->aIndx = n;
			n++;
		}
	}
}

void get_bclusters(MMOL_MD_CELL &mdcell, ARRAY<BASIC_CLUSTER*>& bcluster) {
	int im, ic, ia, n;
	n = 0;
	for (im = 0; im < mdcell.mm.n; im++) n += mdcell.mm.m[im].nCluster;
	n += mdcell.sm.n + mdcell.pm.n;
	bcluster.set_array(n);
	
	n = 0;
	for (im = 0; im < mdcell.mm.n; im++) {
		for (ic = 0; ic < mdcell.mm.m[im].nCluster; ic++) {
			bcluster.m[n] = (BASIC_CLUSTER*)(mdcell.mm.m[im].cluster + ic);
			n++;
		}
	}
	for (im = 0; im < mdcell.sm.n; im++) {
		bcluster.m[n] = mdcell.sm.m[im].c;
		n++;
	}
	for (im = 0; im < mdcell.pm.n; im++) {
		bcluster.m[n] = mdcell.pm.m[im].c;
		n++;
	}
}


double EwaldSum_Recip_2d_Q(MMOL_MD_CELL &mdcell) {
#define _NH 4
#define NH 9
#define AH_MIN 5e-3

	float AhMin = (float)1e-8;

	float A = mdcell.h[0] * mdcell.h[1] * 4;
	double kappa = _KAPPA / ::rcut_Ewd;

	//InitEwaldSumPars(V, inv_3V, kappa, kappa3);

	ARRAY<MATOM*> atom;
	get_atoms(mdcell, atom);

	VECTOR3 h, dr;
	double fcos, h_abs, hR, hz;
	double t1, t2, ferfc1, ferfc2;
	double U_recip, f_recip;
	int nH = int(kappa * mdcell.h[0] * 2 + 0.5); if (nH < 32) nH = 32;
	int nc[3];
	U_recip = 0;
	int n, i;
	int m = 0;
	BASIC_CLUSTER *pc = NULL, *pc1;
	MATOM *patom = NULL, *patom1;

	double EF = 0;

	for (n = 0; n < atom.n; n++) {
		patom = atom.m[n];
		for (i = 0; i < 3; i++) {
			V3zero(patom->E);
		}
	}

	//nH = 20;
	VECTOR3 r;

	double q, q1, U;

	for (n = 0; n < atom.n; n++) {
		patom = atom.m[n]; q = patom->c;
		for (m = 0; m < atom.n; m++) {
			patom1 = atom.m[m]; q1 = patom1->c;
			VECT3(patom->rg, patom1->rg, dr)

			U = 0;
			for (nc[0] = -nH; nc[0] <= nH; nc[0]++) {for (nc[1] = -nH; nc[1] <= nH; nc[1]++) {
				if (nc[0] == 0 && nc[1] == 0) continue;
				for (i = 0; i < 2; i++) h.v[i] = PI / mdcell.h[i] * nc[i];
				h_abs = sqrt(h.v[0] * h.v[0] + h.v[1] * h.v[1]);

				
				hR = dr.v[0] * h.v[0] + dr.v[1] * h.v[1];
				hz = h_abs * dr.v[2];
				fcos = cos(hR);
				t1 = kappa * dr.v[2] + h_abs / (2 * kappa); t2 = -kappa * dr.v[2] + h_abs / (2 * kappa);
				ERFC(ferfc1, t1, i) ERFC(ferfc2, t2, i)
				U += fcos / h_abs * (exp(hz) * ferfc1 + exp(-hz) * ferfc2);
			}}
			U_recip += q * q1 * U;
/*
			//if (patom != patom1) {
				t1 = kappa * dr.v[2]; ERF(t2, t1, i)
				U_recip -= q * q1 * (dr.v[2] * t2 + 1 / kappa * INV_SQRT_PI * exp(-t1 * t1)) * 2;
			//}
			*/
		}
	}
	U_recip *= PI / A * 0.5;


	return U_recip * eUnit_estat_kT;
}

void test_SPME_EwaldSum_2d() {
	int nThreads = MAX_THREADS;

	using namespace _cmm_2d_;

	bool bExplicit = true; // calculate interaction explicitly

	::eps_dc = 1;
	MMOL_MD_CELL mdcell;
	//int nmol = 2;
	//int nmol = 8;
	int nmol = 10;
	//int nmol = 52;
	//int nmol = 104;
	mdcell.SetMol(0, nmol, 0);
	CMM_CELL3D<BASIC_CLUSTER> cmm;
	::mlog.init("test.log", false);

	::rcut_Ewd = 6;
	//float h[3] = {16, 16, 16};
	//float h[3] = {8, 8, 5};
	float h[3] = {5, 5, 5};
	//float h[3] = {4, 4, 4};
	//float h[3] = {32, 32, 16};
	cmm.set_cell(5, 5, 5);
	cmm.set_cell_range(true, -h[0], h[0], -h[1], h[1], -h[2], h[2]);
	mdcell.h[0] = h[0]; mdcell.h[1] = h[1]; mdcell.h[2] = h[2];  mdcell.period = true;
	cmm.set_cell_pos();
	int ncs = 10;
	CMM_array_set_storage<BASIC_CLUSTER>(&cmm, ncs);

	int imol = 0, ia, i;
	SMOLECULE cw;
	VECTOR3 dr, axis;
	axis.v[2] = 1;
	//axis.v[0] = 1;

	extern bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol);
	extern bool cp(SMOLECULE *d, SMOLECULE *s, bool setup = true);
	extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, _EwaldSum_real_::SURF_DIPOLE &mu_surf);
	extern void init_LJ_Pars();
	extern void construct_LJ12_6_prof(int mode);
	extern double cal_charge(MMOL_MD_CELL &mmcell);

	::bPolarize = true; // polarizable water molecule

	ConstructSimpleMolecule(&cw, "SPC");
	cw.c->cal_geometry_radius();
	dr.v[0] = -cw.r.v[0]; dr.v[1] = -cw.r.v[1]; dr.v[2] = -cw.r.v[2];
	cw.shiftMol(dr);
	::iFormat_LJ = 1;  // Amber
	init_LJ_Pars();
	construct_LJ12_6_prof(::iFormat_LJ);
	for (ia = 0; ia < cw.c->nAtoms; ia++) cw.c->atom[ia].par = NULL;

	double Ep = 0;
	
	//cmm_init_subcell<BASIC_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size

	float hgap = h[2] * 0.9, dhs = h[2] - hgap;
/*
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		//dr.v[2] = (h[2] * 0.8) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (ranf() > 0.5 ? 1: -1) * (hgap + dhs * ranf());
		mdcell.sm.m[imol].shiftMol(dr);
		//rand_vect3(axis.v);
		//rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI, true);

		rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI* 0.2, true);

		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
*/	

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (hgap + dhs * ranf());
		mdcell.sm.m[imol].shiftMol(dr);
		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();

		imol++;
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[2] = -dr.v[2];
		mdcell.sm.m[imol].shiftMol(dr);
		axis.v[0] = 1; axis.v[1] = 0; axis.v[2] = 0; // along x axis
		rotate_SM(mdcell.sm.m[imol], axis, PI, true); // rotate along x axis 180 degree
		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
	get_atoms(mdcell, mdcell.atom);
	get_bclusters(mdcell, mdcell.bcluster);
	mdcell.check_eneutral();

	mdcell.set_atom_multipole((bPolarize ? 1 : 0));
	mdcell.setMolIndx();
	molecule_check_periodic_cell(mdcell);
	check_atom_rg(mdcell, nThreads);

	cmm_init_distribute_cluster(mdcell.sm.m, mdcell.sm.n, cmm, true); // add simple solid molecule into cmm

	cmm_check_cluster<BASIC_CLUSTER>(cmm);
	cmm_check<BASIC_CLUSTER>(cmm, MAX_THREADS);
	_cmm_2d_::SMCheckRelshipInCell(mdcell.sm.m, mdcell.sm.n, cmm, MAX_THREADS);

	/*
	{
		for (int imol = 0; imol < mdcell.sm.n; imol++) {
			sprintf(errmsg, "m %d [%f, %f, %f]", imol, mdcell.sm.m[imol].r0.v[0], mdcell.sm.m[imol].r0.v[1], mdcell.sm.m[imol].r0.v[2]);
			show_log(errmsg, true);
			sprintf(errmsg, "m %d has %d of neighbors", imol, number<CMM_IMAGINE<BASIC_CLUSTER> >(mdcell.sm.m[imol].c->spmeRship.ch));
			show_log(errmsg, true);
		}
	}
	*/

	// calculat the total dipole of the whole cell
	cal_dipole(mdcell, 1, mdcell.mu_surf);
	memcpy(&(cmm.mu_surf), &(mdcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));

	// SPME parameters
	extern BSplineFunc<2> bsp;

	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 32;
	//::dim_spme[0] = 20; ::dim_spme[1] = 20; ::dim_spme[2] = 20;
	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 64;
	::dim_spme[0] = 64; ::dim_spme[1] = 64; ::dim_spme[2] = 64;
	//::dim_spme[0] = 128; ::dim_spme[1] = 128; ::dim_spme[2] = 128;

	//float rspline = 4;
	//int nSpline = int(0.5 * rspline / h[0] * dim_spme[0] + 0.5);
	//init_BSplineFunc(bsp, nSpline);
	init_BSplineFunc(bsp, 6);


	VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	vspme_q.bST = true; vspme_q.bST_diagonal = true;
	VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	vspme_mu.bST = true; vspme_mu.bST_diagonal = true;
	if (bPolarize) {
		mdcell.init_polarization_buff();
		mdcell.init_dipole_hist();
		mdcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		//vspme_mu.bsp = &(bsp);
		init_spme(&mdcell, vspme_mu); // init the viral atoms
		vspme_mu.set_BSplineFunc(&bsp);
		//vspme_mu.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);
		vspme_mu.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true);
		//vspme_mu.xl[0] = cmm.xl[0]; vspme_mu.xl[1] = cmm.xl[1]; vspme_mu.xl[2] = cmm.xl[2];
		//vspme_mu.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();

		vspme_mu.init_k();
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();
	}
	else { // charge only
		mdcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		//vspme_q.bsp = &(bsp);
		vspme_q.set_BSplineFunc(&bsp);
		init_spme(&mdcell, vspme_q); // init the viral atoms
		//vspme_q.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);
		vspme_q.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true);

		//vspme_q.xl[0] = cmm.xl[0]; vspme_q.xl[1] = cmm.xl[1]; vspme_q.xl[2] = cmm.xl[2];
		//vspme_q.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();

		vspme_q.init_k();
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}

	_spme_2d_::SPME_VAR spme_var;
	spme_var.set_virial(true, true);

	spme_var.spme_3d_var.esv1.bEwald = true; 
	//spme_var.spme_3d_var.esv1.init_EwaldSum(mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8, ::rcut_Ewd);
	if (bPolarize) {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_mu.V, ::rcut_Ewd);
		spme_var.spme_3d_var.esv1.iSurface = 1; // slab
	}
	else {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_q.V, ::rcut_Ewd);
		spme_var.spme_3d_var.esv1.iSurface = 1; // slab
	}
	//spme_var.spme_3d_var.bVirial = spme_var.bVirial; spme_var.spme_3d_var.bST_diagonal = spme_var.bST_diagonal;
	spme_var.spme_3d_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.spme_3d_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	spme_var.spme_3d_var.bEfieldOnly = false;

	double kappa = spme_var.spme_3d_var.esv1.kappa;
	spme_var.K0var.bEfieldOnly = false;
	spme_var.K0var.bVirial = spme_var.bVirial;
	spme_var.K0var.set(kappa, cmm.xd[0] * cmm.xd[1]);

	spme_var.spme_3d_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	spme_var.K0var.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);

	{
		double Ushape = _spme_2d_::ShapeEnergy_2d(cmm);
		sprintf(errmsg, "shape energy: %f [kT]", Ushape); show_log(errmsg, true);
	}

	//V3zero(cmm.mu_surf.mu) cmm.mu_surf.q = 0; V3zero(cmm.mu_surf.qr) V3zero(cmm.mu_surf.mu0)
	InteractRes ires;
	if (bPolarize) {
		if (!Polarize_SPME_Interact(&mdcell, &cmm, vspme_mu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else SPME_Interact(&mdcell, &cmm, vspme_q, nThreads, spme_var, ires);
	Ep = ires.U;

	sprintf(errmsg, "total dipole: [%f, %f, %f]; from q: [%f, %f, %f], from induced dipole: [%f, %f, %f]", mdcell.mu_surf.mu.v[0], mdcell.mu_surf.mu.v[1], mdcell.mu_surf.mu.v[2], mdcell.mu_surf.qr.v[0], mdcell.mu_surf.qr.v[1], mdcell.mu_surf.qr.v[2], mdcell.mu_surf.mu_ind.v[0], mdcell.mu_surf.mu_ind.v[1], mdcell.mu_surf.mu_ind.v[2]); show_log(errmsg, true);

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);
	
	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-SPME-EwaldSum.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);


	{
	ofstream out;
	out.open("test-atomPos.dat");
	out<<"atom pos. : "<<endl;
	int ia, i;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"[ SM "<<imol<<" ] : "<<mdcell.sm.m[imol].r.v[0]<<"  "<<mdcell.sm.m[imol].r.v[1]<<"  "<<mdcell.sm.m[imol].r.v[2]<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<mdcell.sm.m[imol].c->atom[ia].r.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[2]<<endl;
		}
		out<<endl;
	}
	out.close();
	}



	{
	if (!bPolarize) {
		Ep = EwaldSum_Recip_2d_Q(mdcell);
		sprintf(errmsg, "Restrict Ewald Sum [reciprocal] : %f kT", Ep); show_log(errmsg, true); 

		init_normalized_coordinates(vspme_q, nThreads);
		cal_Q_spme(vspme_q, nThreads); // SPME Q matrix
		double U = vspme_q.cal_U_recp(); //calculate the energy from the image part 	

		sprintf(errmsg, "SPME Ewald Sum [reciprocal] : %f kT", U); show_log(errmsg, true);
	}
	}


	if (bPolarize) vspme_mu.release_plans();
	else vspme_q.release_plans();

	if (!bExplicit) return;

	// explicit calculation
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		V6zero(mdcell.sm.m[imol].c->dyn->fc)
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			V3zero(mdcell.sm.m[imol].c->atom[ia].E)
			V3zero(mdcell.sm.m[imol].c->atom[ia].F)
			for (i = 0; i < mdcell.sm.m[imol].c->atom[ia].f.n; i++) {V3zero(mdcell.sm.m[imol].c->atom[ia].f.m[i]);}
		}
	}

	E2dInteractVar evar_2d;
	evar_2d.bVirial = true; evar_2d.bST_diagonal = true; 
	evar_2d.bEfieldOnly = false; 
	evar_2d.esv1.bEwald = false; // direct interaction, not EwaldSum
	InteractRes tres;

	ires.reset();
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		tres.reset();
		_spme_2d_::ExplicitEInteract_2d(cmm, mdcell.sm.m[imol].c, evar_2d, tres);
		ires += tres;
	}
	Ep = ires.U;

	tres.reset();
	_spme_2d_::ExplicitEInteract_2d_U(cmm, evar_2d, tres);
	Ep += tres.U;

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);

	sprintf(errmsg, "explicit calculation : total energy %f [kT], with dipole-dipole interaction: %f [kT]", Ep, tres.U); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-explicit.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	
	::mlog.close();
}



void test_SPME_EwaldSum_2d_Exclude_InducedDipole() {
	int nThreads = MAX_THREADS;

	using namespace _cmm_2d_;

	bool bExplicit = true; // calculate interaction explicitly

	::eps_dc = 1;
	MMOL_MD_CELL mdcell;
	//int nmol = 2;
	//int nmol = 8;
	int nmol = 10;
	//int nmol = 52;
	//int nmol = 104;
	mdcell.SetMol(0, nmol, 0);
	CMM_CELL3D<BASIC_CLUSTER> cmm;
	::mlog.init("test.log", false);

	::rcut_Ewd = 6;
	//float h[3] = {16, 16, 16};
	//float h[3] = {8, 8, 5};
	float h[3] = {5, 5, 5};
	//float h[3] = {4, 4, 4};
	//float h[3] = {32, 32, 16};
	cmm.set_cell(5, 5, 5);
	cmm.set_cell_range(true, -h[0], h[0], -h[1], h[1], -h[2], h[2]);
	mdcell.h[0] = h[0]; mdcell.h[1] = h[1]; mdcell.h[2] = h[2];  mdcell.period = true;
	cmm.set_cell_pos();
	int ncs = 10;
	CMM_array_set_storage<BASIC_CLUSTER>(&cmm, ncs);

	int imol = 0, ia, i;
	SMOLECULE cw;
	VECTOR3 dr, axis;
	axis.v[2] = 1;
	//axis.v[0] = 1;

	extern bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol);
	extern bool cp(SMOLECULE *d, SMOLECULE *s, bool setup = true);
	extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, _EwaldSum_real_::SURF_DIPOLE &mu_surf);
	extern void init_LJ_Pars();
	extern void construct_LJ12_6_prof(int mode);
	extern double cal_charge(MMOL_MD_CELL &mmcell);

	::bPolarize = true; // polarizable water molecule

	ConstructSimpleMolecule(&cw, "SPC");
	cw.c->cal_geometry_radius();
	dr.v[0] = -cw.r.v[0]; dr.v[1] = -cw.r.v[1]; dr.v[2] = -cw.r.v[2];
	cw.shiftMol(dr);
	::iFormat_LJ = 1;  // Amber
	init_LJ_Pars();
	construct_LJ12_6_prof(::iFormat_LJ);
	for (ia = 0; ia < cw.c->nAtoms; ia++) cw.c->atom[ia].par = NULL;

	double Ep = 0;
	
	//cmm_init_subcell<BASIC_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size

	float hgap = h[2] * 0.9, dhs = h[2] - hgap;

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		//dr.v[2] = (h[2] * 0.8) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (ranf() > 0.5 ? 1: -1) * (hgap + dhs * ranf());
		mdcell.sm.m[imol].shiftMol(dr);
		//rand_vect3(axis.v);
		//rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI, true);

		rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI* 0.2, true);

		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
	get_atoms(mdcell, mdcell.atom);
	get_bclusters(mdcell, mdcell.bcluster);
	mdcell.check_eneutral();

	mdcell.set_atom_multipole((bPolarize ? 1 : 0));
	mdcell.setMolIndx();
	molecule_check_periodic_cell(mdcell);
	check_atom_rg(mdcell, nThreads);

	cmm_init_distribute_cluster(mdcell.sm.m, mdcell.sm.n, cmm, true); // add simple solid molecule into cmm

	cmm_check_cluster<BASIC_CLUSTER>(cmm);
	cmm_check<BASIC_CLUSTER>(cmm, MAX_THREADS);
	_cmm_2d_::SMCheckRelshipInCell(mdcell.sm.m, mdcell.sm.n, cmm, MAX_THREADS);

	/*
	{
		for (int imol = 0; imol < mdcell.sm.n; imol++) {
			sprintf(errmsg, "m %d [%f, %f, %f]", imol, mdcell.sm.m[imol].r0.v[0], mdcell.sm.m[imol].r0.v[1], mdcell.sm.m[imol].r0.v[2]);
			show_log(errmsg, true);
			sprintf(errmsg, "m %d has %d of neighbors", imol, number<CMM_IMAGINE<BASIC_CLUSTER> >(mdcell.sm.m[imol].c->spmeRship.ch));
			show_log(errmsg, true);
		}
	}
	*/

	// calculat the total dipole of the whole cell
	cal_dipole(mdcell, 1, mdcell.mu_surf);
	memcpy(&(cmm.mu_surf), &(mdcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));

	::bExcludeInducedDipole = true;

	// SPME parameters
	extern BSplineFunc<2> bsp;

	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 32;
	::dim_spme[0] = 20; ::dim_spme[1] = 20; ::dim_spme[2] = 20;
	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 64;
	//::dim_spme[0] = 64; ::dim_spme[1] = 64; ::dim_spme[2] = 64;
	//::dim_spme[0] = 128; ::dim_spme[1] = 128; ::dim_spme[2] = 128;

	//float rspline = 4;
	//int nSpline = int(0.5 * rspline / h[0] * dim_spme[0] + 0.5);
	//init_BSplineFunc(bsp, nSpline);
	init_BSplineFunc(bsp, 6);


	VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	vspme_q.bST = true; vspme_q.bST_diagonal = true;
	VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	vspme_mu.bST = true; vspme_mu.bST_diagonal = true;
	float q0 = 0;
	VSPME<_VMUatom> vspme_mu_induced; vspme_mu_induced.mMP = 1; // with dipole
	vspme_mu_induced.bST = true; vspme_mu_induced.bST_diagonal = true;
	if (bPolarize) {
		mdcell.init_polarization_buff();
		mdcell.init_dipole_hist();
		mdcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		//vspme_mu.bsp = &(bsp);
		init_spme(&mdcell, vspme_mu); // init the viral atoms
		vspme_mu.set_BSplineFunc(&bsp);
		//vspme_mu.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);
		vspme_mu.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true);
		//vspme_mu.xl[0] = cmm.xl[0]; vspme_mu.xl[1] = cmm.xl[1]; vspme_mu.xl[2] = cmm.xl[2];
		//vspme_mu.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();

		vspme_mu.init_k();
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();

		init_spme_induced(&mdcell, vspme_mu_induced, &(q0)); // init the viral atoms
		vspme_mu_induced.set_BSplineFunc(&bsp);
		vspme_mu_induced.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true); // important, after init();
		vspme_mu_induced.init_k(); vspme_mu_induced.init_b(); vspme_mu_induced.init_C(); vspme_mu_induced.init_vars();
	}
	else { // charge only
		mdcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		//vspme_q.bsp = &(bsp);
		vspme_q.set_BSplineFunc(&bsp);
		init_spme(&mdcell, vspme_q); // init the viral atoms
		//vspme_q.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);
		vspme_q.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true);

		//vspme_q.xl[0] = cmm.xl[0]; vspme_q.xl[1] = cmm.xl[1]; vspme_q.xl[2] = cmm.xl[2];
		//vspme_q.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();

		vspme_q.init_k();
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}

	_spme_2d_::SPME_VAR spme_var;
	spme_var.bVirial = true; spme_var.bST_diagonal = true;

	spme_var.spme_3d_var.esv1.bEwald = true; 
	//spme_var.spme_3d_var.esv1.init_EwaldSum(mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8, ::rcut_Ewd);
	if (bPolarize) {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_mu.V, ::rcut_Ewd);
		spme_var.spme_3d_var.esv1.iSurface = 1; // slab
	}
	else {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_q.V, ::rcut_Ewd);
		spme_var.spme_3d_var.esv1.iSurface = 1; // slab
	}
	spme_var.spme_3d_var.bVirial = spme_var.bVirial; spme_var.spme_3d_var.bST_diagonal = spme_var.bST_diagonal;
	spme_var.spme_3d_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.spme_3d_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	spme_var.spme_3d_var.bEfieldOnly = false;

	double kappa = spme_var.spme_3d_var.esv1.kappa;
	spme_var.K0var.bEfieldOnly = false;
	spme_var.K0var.bVirial = spme_var.bVirial;
	spme_var.K0var.set(kappa, cmm.xd[0] * cmm.xd[1]);

	spme_var.spme_3d_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	spme_var.K0var.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);

	{
		double Ushape = _spme_2d_::ShapeEnergy_2d(cmm);
		sprintf(errmsg, "shape energy: %f [kT]", Ushape); show_log(errmsg, true);
	}

	//V3zero(cmm.mu_surf.mu) cmm.mu_surf.q = 0; V3zero(cmm.mu_surf.qr) V3zero(cmm.mu_surf.mu0)
	InteractRes ires;
	if (bPolarize) {
		if (!Polarize_SPME_Interact2(&mdcell, &cmm, vspme_mu, vspme_mu_induced, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else SPME_Interact(&mdcell, &cmm, vspme_q, nThreads, spme_var, ires);
	Ep = ires.U;

	sprintf(errmsg, "total dipole: [%f, %f, %f]; from q: [%f, %f, %f], from induced dipole: [%f, %f, %f]", mdcell.mu_surf.mu.v[0], mdcell.mu_surf.mu.v[1], mdcell.mu_surf.mu.v[2], mdcell.mu_surf.qr.v[0], mdcell.mu_surf.qr.v[1], mdcell.mu_surf.qr.v[2], mdcell.mu_surf.mu_ind.v[0], mdcell.mu_surf.mu_ind.v[1], mdcell.mu_surf.mu_ind.v[2]); show_log(errmsg, true);

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);
	
	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-SPME-EwaldSum.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);


	{
	ofstream out;
	out.open("test-atomPos.dat");
	out<<"atom pos. : "<<endl;
	int ia;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"[ SM "<<imol<<" ] : "<<mdcell.sm.m[imol].r.v[0]<<"  "<<mdcell.sm.m[imol].r.v[1]<<"  "<<mdcell.sm.m[imol].r.v[2]<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<mdcell.sm.m[imol].c->atom[ia].r.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[2]<<endl;
		}
		out<<endl;
	}
	out.close();
	}



	{
	if (!bPolarize) {
		Ep = EwaldSum_Recip_2d_Q(mdcell);
		sprintf(errmsg, "Restrict Ewald Sum [reciprocal] : %f kT", Ep); show_log(errmsg, true); 

		init_normalized_coordinates(vspme_q, nThreads);
		cal_Q_spme(vspme_q, nThreads); // SPME Q matrix
		double U = vspme_q.cal_U_recp(); //calculate the energy from the image part 	

		sprintf(errmsg, "SPME Ewald Sum [reciprocal] : %f kT", U); show_log(errmsg, true);
	}
	}


	if (bPolarize) vspme_mu.release_plans();
	else vspme_q.release_plans();

	if (!bExplicit) return;

	// explicit calculation
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		V6zero(mdcell.sm.m[imol].c->dyn->fc)
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			V3zero(mdcell.sm.m[imol].c->atom[ia].E)
			V3zero(mdcell.sm.m[imol].c->atom[ia].F)
			for (i = 0; i < mdcell.sm.m[imol].c->atom[ia].f.n; i++) {V3zero(mdcell.sm.m[imol].c->atom[ia].f.m[i]);}
		}
	}

	E2dInteractVar evar_2d;
	evar_2d.bVirial = true; evar_2d.bST_diagonal = true; 
	evar_2d.bEfieldOnly = false; 
	evar_2d.esv1.bEwald = false; // direct interaction, not EwaldSum
	InteractRes tres;

	ires.reset();
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		tres.reset();
		_spme_2d_::ExplicitEInteract_2d_ExcludeInducedDiploe(cmm, mdcell.sm.m[imol].c, evar_2d, tres);
		ires += tres;
	}
	Ep = ires.U;

	tres.reset();
	_spme_2d_::ExplicitEInteract_2d_U(cmm, evar_2d, tres);
	Ep += tres.U;

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);

	sprintf(errmsg, "explicit calculation : total energy %f [kT], with dipole-dipole interaction: %f [kT]", Ep, tres.U); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-explicit.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	
	::mlog.close();
}



//***********************************************************
//     VSPME_ImageSlab for 2-d Slab with Image Charge / Atom 
//  VSPME calculate interaction between all the atoms
//  So, a slab with image charges needs to exclude the  
//  interactions between image charges(atoms), and the 
//  interactions between real atoms with their image charges(atoms)
//  respectively with same x-y, non-periodical image-cell
//***********************************************************

//********************************************************************
//                02/03/2012
//   above description is wrong! the interaction between a charge with its image charge
//   can NOT be excluded, which is the interaction of charge with the planner interface,
//   a part of the interaction energy. However, interactions between image charges needs
//   to be excluded. Anyway, see below, where A is real slab + imag1, B is real slab + img2
//********************************************************************
//********************************************************************
//   !!!! Important Modification on 02/01/2012  !!!!
//  energy of a slab with two interfaces (image charges), the electrostatic interaction
//  energy is: U = U_r + (U_A - U_r - U_img1) + (U_B - U_r - U_img2)
//               = U_A + U_B - U_r - U_img1 - U_img2
// The force calculation above is also not correct, because the coordination of image charge is associated with real charge.
// Correct way to calculate force is:
//    electrostatic interaction energy between real charges A and B with one interface is:
//      U_AB = U_A_B + U_A_Ai + U_A_Bi + U_Ai_B + U_B_Bi, where Ai and Bi are image charge of A and B respectively
//    because r_Ai = r_x_A - r_z_A, the force on A is:
//   f_A = f_A_B + f_A_Ai + f_A_Bi + f_Ai_B_T + f_Ai_A_T, 
// here f_A_ri = f_A_B + f_A_Ai + f_A_Bi is the calculated force on A within system (real + imagine)
// f_Ai_B_T is imagine of f_Ai_B, f_x_Ai_B_T = f_Ai_B and f_z_Ai_B_T = -f_z_Ai_B
// f_Ai_A_T is imagine of f_Ai_A, f_x_Ai_A_T = f_Ai_A and f_z_Ai_A_T = -f_z_Ai_A
// and, f_Ai_ri = f_Ai_A + f_Ai_B + f_Ai_Bi is the calculated force on A within system (real + imagine)
// so, f_Ai_A + f_Ai_B = f_Ai_ri - f_Ai_Bi, where f_Ai_Bi can be calculated within imagine system
// this means the final force on A should be calculated:
// f_A = f_A_ri + (f_Ai_ri - f_Ai_i)_T, here _T means imagine transform of the force on Ai, keeping xy but inverse of z

// with two interfaces, f_A = f_A_ri1 + (f_Ai_ri1 - f_Ai_i1)_T + f_A_ri2 + (f_Ai_ri2 - f_Ai_i2)_T - f_A_B_r

// same expression has to be done for E
//********************************************************************

//********************************************************************
//                02/03/2012
// Also, in case the slab has thickness thicker than 2 * r_cut, each real atom in slab, or imag1 or img2,
// will have short-range neighbors for real part EwardSum and k=0 receipt part, in slab + img1, or slab + img2, but NEVER in slab + img1 + img2
// This property can be used for real part EwardSum and k=0 receipt part. As shown above
// U = U_A + U_B - U_r - U_img1 - U_img2
// for an atom in real slab, the short-range neighbors in A can be the atoms in slab and img1
// and the short-range neighbors in B can be the atoms in slab and img2. So, for real atom in slab,
// U = U_A + U_B - U_r = U in a unit system slab + img1 + img2

// for atoms in img1, U = U_A - U_img1 = U of image charge with real atom in slab only. Similary,
// for atoms in img2, U = U_A - U_img2 = U of image charge with real atom in slab only
//
// same is for E and F calculation
//********************************************************************

namespace _EwaldSum_2d_ {
// cell has real charge, image1 and image2
// also we assume that the cluster of image charge has cluster ID < 0
void CLUSTER_EFIELD_SPME_REAL_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, _EwaldSum_real_::SPME_VAR& var, InteractRes& res) {
	using namespace _Thole_polarize_;
	EwaldSumRealVars1 *esv = &(var.esv1);
	TholeCorrTensor *tv = &(var.tv);

	MATOM *patom = NULL;
	VECTOR3 E_real, F_real, E_Thole, F_Thole;
	VECTOR3 dr, u;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	CHAIN<BASIC_CLUSTER> *ch_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, SPME_RELSHIP> rshipIt;

	double R, R2;
	double t1, t2;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, local = true, bEwaldCut = false;
	bool bEwaldKickoff = false, bSameCluster = false;

	double alpha_Thole = 0, a_Thole = ::_Thole_polarize_::a;

	double Utotal = 0, Uestat = 0, Uext = 0, U_real = 0, U_Thole = 0;

	char aindx[2] = {0x00, 0x00};
	unsigned int key;
	PAIRWISE<float> *TholeRadius = NULL;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	res.reset();

	//char msg[256] = "\0";
	//sprintf(msg, "CLUSTER %d :", pc->mIndx); show_log(msg, true);

	//for (na = 0; na < pc->nAtoms; na++) {
		//patom = pc->atom + na;
	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];
		// reset E & F of the atom
		V3zero(patom->E) 
		//if (!bEfieldOnly) {V3zero(patom->F)}
		if (patom->eNeutral) continue;
		Uestat = 0;
		//for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		V32V3(patom->rg, r0)

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		//ch_imag_cluster = pc->spmeRship.ch;
		//while (ch_imag_cluster != NULL) {
		rshipIt.set(&(pc->spmeRship)); rshipIt.start_iterate();
		while (!rshipIt.eoi()) {
			imag_cluster = rshipIt.current(); if (imag_cluster == NULL) {rshipIt.next_iterate(); continue;}
			pc1 = imag_cluster->p;
			local = imag_cluster->local;
			drCell.v[0] = imag_cluster->nx * cell.xd[0];
			drCell.v[1] = imag_cluster->ny * cell.xd[1];
			drCell.v[2] = imag_cluster->nz * cell.xd[2];
			//for (na1 = 0; na1 < pc1->nAtoms; na1++) {
				//patom1 = pc1->atom + na1;
			for (na1 = 0; na1 < pc1->fatom.n; na1++) {
				patom1 = pc1->fatom.m[na1];
				bSameCluster = (local && pc1 == pc ? true : false);
				if (bSameCluster && na1 == na) continue; // same atom
				if (patom1->eNeutral) continue; // ignore this atom for charge interaction
				for (i = 0; i < 3; i++) {
					//dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
					dr.v[i] = r0.v[i] - patom1->rg.v[i] - drCell.v[i];
				}
				if (FABS(dr.v[0]) > rcut_Ewd || FABS(dr.v[1]) > rcut_Ewd || FABS(dr.v[2]) > rcut_Ewd) {
					bEwaldCut = true;
					continue;
				}
				V3ABS2(dr, R2)
				//if (R2 < _R2min_Interact) continue;
				bEwaldCut = (R2 > r2cut_Ewd ? true : false);
				if (bEwaldCut) continue;

				R = sqrt(R2); 

				//if (R < 3) {
				//	sprintf(msg, "C %d M %d -- C %d M %d,  CELL [%d, %d, %d] : |R| = %f", pc->mIndx, na, pc1->mIndx, na1, imag_cluster->nx, imag_cluster->ny, imag_cluster->nz, R); show_log(msg, true);
				//}

				if (pc->mIndx == pc1->mIndx && R2 < d2_14) {
					BOUND(patom, patom1, i, bound) // bound atom
					if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
					if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
					else bound14 = false;
				}
				else {bound = false; bound14 = false;}

				if (bound || bSameCluster) { // bounding atoms or they are in the same cluster
					/*
					for (i = 0; i < 3; i++) {
						t1 = patom1->c * dr.v[i] / R3;
						E_real.v[i] -= t1 * (1 - f1); // bound interaction has to be kicked off in Ewald Sum
					}
					// energy, no self interaction
					t1 = patom1->c / R;
					Uestat +=  0.5 * t1 * (f0 - 1); // bound interaction has to be kicked off in Ewald Sum
					*/
					bEwaldKickoff = true; t1 = 1;
				}
				else if (bound14 && pc != pc1) { // 1-4 interaction
					/*
					for (i = 0; i < 3; i++) {
						t1 = patom1->c * dr.v[i] / R3;
						E_real.v[i] += t1 * f1;
						E_14.v[i] += (A14_E - 1) * t1; // 1-4 interaction will to be corrected later
					}
					// energy
					t1 = patom1->c / R;
					Uestat +=  0.5 * t1 * (f0 - 1 + A14_E); // 1-4 interaction is corrected by * A14_E
					*/
					bEwaldKickoff = true; t1 = 1 - A14_E;
				}
				else {
					/*
					for (i = 0; i < 3; i++) E_real.v[i] += f1 * patom1->c * dr.v[i] / R3;
					// energy
					t1 = patom1->c / R;
					Uestat += 0.5 * t1 * f0;
					*/
					bEwaldKickoff = false; t1 = 0;
				}
				if (var.bEfieldOnly) {
					EwaldReal_E(*esv, patom, patom1, dr, R, R2, E_real, bEwaldKickoff, t1);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr_E(*tv, a_Thole, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v)) {
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
								}
								V3plusV3(E_real, E_Thole, E_real)
							}
						}
					}
				}
				else {
					EwaldReal(*esv, patom, patom1, dr, R, R2, E_real, F_real, U_real, bEwaldKickoff, t1, bExcludeInducedDipole);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr(*tv, a_Thole, patom->c, patom->mu.v, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v, F_Thole.v, U_Thole)) {
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
									F_Thole.v[0] *= t2; F_Thole.v[1] *= t2; F_Thole.v[2] *= t2;
									U_Thole *= t2;
								}
								V3plusV3(E_real, E_Thole, E_real) V3plusV3(F_real, F_Thole, F_real) U_real += U_Thole;
							}
						}
					}
				}

				// important, if the interaction is with the image charge, the contribution of Ewald Sum is in half
				// and we assume the image cluster has cluster ID < 0 
				// Why is the contribution NOT half? 
				/*
				if (pc1->cID < 0) {
					U_real *= 0.5;
					E_real.v[0] *= 0.5; E_real.v[1] *= 0.5; E_real.v[2] *= 0.5;
					F_real.v[0] *= 0.5; F_real.v[1] *= 0.5; F_real.v[2] *= 0.5;
				}
				*/

				patom->E.v[0] += E_real.v[0]; patom->E.v[1] += E_real.v[1]; patom->E.v[2] += E_real.v[2];
				if (!var.bEfieldOnly) {
					// the forece on atom has unit on mass unit
					F_real.v[0] *= fUnit_estat_mass; 
					F_real.v[1] *= fUnit_estat_mass; 
					F_real.v[2] *= fUnit_estat_mass;

					patom->F.v[0] += F_real.v[0]; 
					patom->F.v[1] += F_real.v[1]; 
					patom->F.v[2] += F_real.v[2];
					Uestat += 0.5 * U_real;
					if (var.bVirial) Accumulate_Virial_Diagonal(res, F_real.v, dr.v);
				}
			}
			//ch_imag_cluster = ch_imag_cluster->next;
			rshipIt.next_iterate();
		}
		//show_log("\n", true);
		// in 2d system, no contribution from surface induced dipole
		// the electric field from the surface induced dipole
		for (i = 0; i < 3; i++) {
			t1 = esv->inv_3V * cell.mu_surf.mu.v[i];
			patom->E.v[i] -= t1;
			if (!var.bEfieldOnly) patom->F.v[i] -= patom->c * t1 * fUnit_estat_mass;
		}

		if (var.bEfieldOnly) {
			// VERY IMPORTANT :
			// electric field calculation stop here
			// if more code is required for electric field calculation, the code needs to added before this line
			continue;
		}

		Uestat -= esv->kappa * INV_SQRT_PI * patom->c * patom->c; // self energy from charge
		if (patom->mMP > 0) { // self energy from dipole
			V3ABS2(patom->mu, t1)
			Uestat -= 2 * esv->kappa3 / 3 * INV_SQRT_PI * t1;
		}

		// no surface dipole effect for 2d system
		// energy from the surface induced dipole 
		/*
		switch (patom->mMP) {
		case 0:
			Uestat += 0.5 * esv->inv_3V * patom->c * (patom->r0.v[0] * cell.mu_surf.mu.v[0] + patom->r0.v[1] * cell.mu_surf.mu.v[1] + patom->r0.v[2] * cell.mu_surf.mu.v[2]);
			break;
		case 1:
			Uestat += 0.5 * esv->inv_3V * ((patom->c * patom->r0.v[0] + patom->mu.v[0]) * cell.mu_surf.mu.v[0] 
			+ (patom->c * patom->r0.v[1] + patom->mu.v[1]) * cell.mu_surf.mu.v[1] 
			+ (patom->c * patom->r0.v[2] + patom->mu.v[2]) * cell.mu_surf.mu.v[2]);
			break;
		default:
			break;
		}
		*/

		Utotal += Uestat;

	#if _CHECK_TIME_STAMP_ == 1
		//t_real += mt.elapse();
		//time.start();
	#endif
	}
	if (var.bEfieldOnly) return;

	// virial coming from the surface dipole will be calculated later

	// total energy in kT
	res.U = Utotal * eUnit_estat_kT + Uext;
}

// *********  02/03/2012  **********
// this is to calculate the U, E & F on the image charge
// as discussed before, on image charge, only the short-range interaction with real atoms are required
// also, because they are image charge, there is no difference about whether the interacted atoms in same cluster or not
// cell has real charge, image1 and image2
// also we assume that the cluster of image charge has cluster ID < 0
void ImageCluster_EFIELD_SPME_REAL_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, _EwaldSum_real_::SPME_VAR& var, InteractRes& res) {
	if (pc->cID >= 0) {
		return; // this is not a image cluster
	}
	using namespace _Thole_polarize_;
	EwaldSumRealVars1 *esv = &(var.esv1);
	TholeCorrTensor *tv = &(var.tv);

	MATOM *patom = NULL;
	VECTOR3 E_real, F_real, E_Thole, F_Thole;
	VECTOR3 dr, u;
	int i, j, k, na, na1;
	VECTOR3 r0, r1;

	CMM_IMAGINE<BASIC_CLUSTER> *imag_cluster = NULL;
	CHAIN<BASIC_CLUSTER> *ch_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	RELSHIP_ITERATOR<BASIC_CLUSTER, SPME_RELSHIP> rshipIt;

	double R, R2;
	double t1, t2;
	VECTOR3 drCell;

	bool bound = false, bound14 = false, local = true, bEwaldCut = false;
	bool bEwaldKickoff = false, bSameCluster = false;

	double alpha_Thole = 0, a_Thole = ::_Thole_polarize_::a;

	double Utotal = 0, Uestat = 0, Uext = 0, U_real = 0, U_Thole = 0;

	char aindx[2] = {0x00, 0x00};
	unsigned int key;
	PAIRWISE<float> *TholeRadius = NULL;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	res.reset();

	//char msg[256] = "\0";
	//sprintf(msg, "CLUSTER %d :", pc->mIndx); show_log(msg, true);

	//for (na = 0; na < pc->nAtoms; na++) {
		//patom = pc->atom + na;
	for (na = 0; na < pc->fatom.n; na++) {
		patom = pc->fatom.m[na];
		// reset E & F of the atom
		V3zero(patom->E) 
		if (!var.bEfieldOnly) {V3zero(patom->F)}
		if (patom->eNeutral) continue;
		Uestat = 0;
		//for (i = 0; i < 3; i++) r0.v[i] = patom->r.v[i] + pc->dr_cmm.v[i];
		V32V3(patom->rg, r0)

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		//ch_imag_cluster = pc->spmeRship.ch;
		//while (ch_imag_cluster != NULL) {
		rshipIt.set(&(pc->spmeRship)); rshipIt.start_iterate();
		while (!rshipIt.eoi()) {
			imag_cluster = rshipIt.current(); if (imag_cluster == NULL) {rshipIt.next_iterate(); continue;}
			pc1 = imag_cluster->p;
			if (pc1->cID < 0) {// pc1 is image cluster, ignore it
				rshipIt.next_iterate(); continue;
			}
			local = imag_cluster->local;
			drCell.v[0] = imag_cluster->nx * cell.xd[0];
			drCell.v[1] = imag_cluster->ny * cell.xd[1];
			drCell.v[2] = imag_cluster->nz * cell.xd[2];
			//for (na1 = 0; na1 < pc1->nAtoms; na1++) {
				//patom1 = pc1->atom + na1;
			for (na1 = 0; na1 < pc1->fatom.n; na1++) {
				patom1 = pc1->fatom.m[na1];
//				bSameCluster = (local && pc1 == pc ? true : false);
//				if (bSameCluster && na1 == na) continue; // same atom
				if (patom1->eNeutral) continue; // ignore this atom for charge interaction
				for (i = 0; i < 3; i++) {
					//dr.v[i] = r0.v[i] - patom1->r.v[i] - pc1->dr_cmm.v[i] - drCell.v[i];
					dr.v[i] = r0.v[i] - patom1->rg.v[i] - drCell.v[i];
				}
				if (FABS(dr.v[0]) > rcut_Ewd || FABS(dr.v[1]) > rcut_Ewd || FABS(dr.v[2]) > rcut_Ewd) {
					bEwaldCut = true;
					continue;
				}
				V3ABS2(dr, R2)
				//if (R2 < _R2min_Interact) continue;
				bEwaldCut = (R2 > r2cut_Ewd ? true : false);
				if (bEwaldCut) continue;

				R = sqrt(R2); 

				//if (R < 3) {
				//	sprintf(msg, "C %d M %d -- C %d M %d,  CELL [%d, %d, %d] : |R| = %f", pc->mIndx, na, pc1->mIndx, na1, imag_cluster->nx, imag_cluster->ny, imag_cluster->nz, R); show_log(msg, true);
				//}

//				if (pc->mIndx == pc1->mIndx && R2 < d2_14) {
//					BOUND(patom, patom1, i, bound) // bound atom
//					if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
//					if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
//					else bound14 = false;
//				}
//				else {bound = false; bound14 = false;}

//				if (bound || bSameCluster) { // bounding atoms or they are in the same cluster
					/*
					for (i = 0; i < 3; i++) {
						t1 = patom1->c * dr.v[i] / R3;
						E_real.v[i] -= t1 * (1 - f1); // bound interaction has to be kicked off in Ewald Sum
					}
					// energy, no self interaction
					t1 = patom1->c / R;
					Uestat +=  0.5 * t1 * (f0 - 1); // bound interaction has to be kicked off in Ewald Sum
					*/
//					bEwaldKickoff = true; t1 = 1;
//				}
//				else if (bound14 && pc != pc1) { // 1-4 interaction
					/*
					for (i = 0; i < 3; i++) {
						t1 = patom1->c * dr.v[i] / R3;
						E_real.v[i] += t1 * f1;
						E_14.v[i] += (A14_E - 1) * t1; // 1-4 interaction will to be corrected later
					}
					// energy
					t1 = patom1->c / R;
					Uestat +=  0.5 * t1 * (f0 - 1 + A14_E); // 1-4 interaction is corrected by * A14_E
					*/
//					bEwaldKickoff = true; t1 = 1 - A14_E;
//				}
//				else {
					/*
					for (i = 0; i < 3; i++) E_real.v[i] += f1 * patom1->c * dr.v[i] / R3;
					// energy
					t1 = patom1->c / R;
					Uestat += 0.5 * t1 * f0;
					*/
					bEwaldKickoff = false; t1 = 0;
//				}
				if (var.bEfieldOnly) {
					EwaldReal_E(*esv, patom, patom1, dr, R, R2, E_real, bEwaldKickoff, t1);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr_E(*tv, a_Thole, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v)) {
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
								}
								V3plusV3(E_real, E_Thole, E_real)
							}
						}
					}
				}
				else {
					EwaldReal(*esv, patom, patom1, dr, R, R2, E_real, F_real, U_real, bEwaldKickoff, t1, bExcludeInducedDipole);
					if (var.bTholePolarize && !bSameCluster) {
						aindx[0] = patom->par->indx; aindx[1] = patom1->par->indx; key = construct_char2_key_pairwise(aindx);
						TholeRadius = ::TholeRadius_db.search(key);
						if (TholeRadius != NULL) {
							a_Thole = TholeRadius->pv;
							if (::_Thole_polarize_::TholeCorr(*tv, a_Thole, patom->c, patom->mu.v, patom1->c, patom1->mu.v, R, R2, dr.v, E_Thole.v, F_Thole.v, U_Thole)) {
								if (bEwaldKickoff) {
									t2 = 1 - t1;
									E_Thole.v[0] *= t2; E_Thole.v[1] *= t2; E_Thole.v[2] *= t2;
									F_Thole.v[0] *= t2; F_Thole.v[1] *= t2; F_Thole.v[2] *= t2;
									U_Thole *= t2;
								}
								V3plusV3(E_real, E_Thole, E_real) V3plusV3(F_real, F_Thole, F_real) U_real += U_Thole;
							}
						}
					}
				}

				// important, if the interaction is with the image charge, the contribution of Ewald Sum is in half
				// and we assume the image cluster has cluster ID < 0 
				// Why is the contribution NOT half? 
				/*
				if (pc1->cID < 0) {
					U_real *= 0.5;
					E_real.v[0] *= 0.5; E_real.v[1] *= 0.5; E_real.v[2] *= 0.5;
					F_real.v[0] *= 0.5; F_real.v[1] *= 0.5; F_real.v[2] *= 0.5;
				}
				*/

				patom->E.v[0] += E_real.v[0]; patom->E.v[1] += E_real.v[1]; patom->E.v[2] += E_real.v[2];
				if (!var.bEfieldOnly) {
					// the forece on atom has unit on mass unit
					F_real.v[0] *= fUnit_estat_mass; 
					F_real.v[1] *= fUnit_estat_mass; 
					F_real.v[2] *= fUnit_estat_mass;

					patom->F.v[0] += F_real.v[0]; 
					patom->F.v[1] += F_real.v[1]; 
					patom->F.v[2] += F_real.v[2];
					Uestat += 0.5 * U_real;
					if (var.bVirial) Accumulate_Virial_Diagonal(res, F_real.v, dr.v);
				}
			}
			//ch_imag_cluster = ch_imag_cluster->next;
			rshipIt.next_iterate();
		}
		//show_log("\n", true);
		// in 2d system, no contribution from surface induced dipole
		// the electric field from the surface induced dipole
		/*
		for (i = 0; i < 3; i++) {
			t1 = esv->inv_3V * cell.mu_surf.mu.v[i];
			patom->E.v[i] -= t1;
			if (!var.bEfieldOnly) patom->F.v[i] -= patom->c * t1 * fUnit_estat_mass;
		}
		*/

		if (var.bEfieldOnly) {
			// VERY IMPORTANT :
			// electric field calculation stop here
			// if more code is required for electric field calculation, the code needs to added before this line
			continue;
		}

		// self energy needs to be excluded
		/*
		Uestat -= esv->kappa * INV_SQRT_PI * patom->c * patom->c; // self energy from charge
		if (patom->mMP > 0) { // self energy from dipole
			V3ABS2(patom->mu, t1)
			Uestat -= 2 * esv->kappa3 / 3 * INV_SQRT_PI * t1;
		}
		*/

		// no surface dipole effect for 2d system
		// energy from the surface induced dipole 
		/*
		switch (patom->mMP) {
		case 0:
			Uestat += 0.5 * esv->inv_3V * patom->c * (patom->r0.v[0] * cell.mu_surf.mu.v[0] + patom->r0.v[1] * cell.mu_surf.mu.v[1] + patom->r0.v[2] * cell.mu_surf.mu.v[2]);
			break;
		case 1:
			Uestat += 0.5 * esv->inv_3V * ((patom->c * patom->r0.v[0] + patom->mu.v[0]) * cell.mu_surf.mu.v[0] 
			+ (patom->c * patom->r0.v[1] + patom->mu.v[1]) * cell.mu_surf.mu.v[1] 
			+ (patom->c * patom->r0.v[2] + patom->mu.v[2]) * cell.mu_surf.mu.v[2]);
			break;
		default:
			break;
		}
		*/

		Utotal += Uestat;

	#if _CHECK_TIME_STAMP_ == 1
		//t_real += time.elapse();
		//time.start();
	#endif
	}
	if (var.bEfieldOnly) return;

	// virial coming from the surface dipole will be calculated later

	// total energy in kT
	res.U = Utotal * eUnit_estat_kT + Uext;
}

void SingleMM_EFIELD_SPME_REAL_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, _EwaldSum_real_::SPME_VAR& var, InteractRes& res) {
	int nc;
	BASIC_CLUSTER *pc = NULL;

	res.reset();
	InteractRes tres; 
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm.nCluster) break;
		pc = (BASIC_CLUSTER*)(mm.cluster + nc);
		CLUSTER_EFIELD_SPME_REAL_ImageSlab(cell, pc, var, tres);
		res += tres;
	}
}

void MM_EFIELD_SPME_REAL_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE *mm, int nMM, _EwaldSum_real_::SPME_VAR& var, InteractRes& res) {
	int nm, nc;
	
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < nMM; nm++) {
		//SingleMM_EFIELD_SPME_REAL(cell, mm[nm], 0, mm[nm].nCluster, U, bEfieldOnly);
		for (nc = 0; nc < mm[nm].nCluster; nc++) {
			CLUSTER_EFIELD_SPME_REAL_ImageSlab(cell, mm[nm].cluster + nc, var, tres);
			res += tres;
		}
	}
}

void SM_EFIELD_SPME_REAL_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE *sm, int nSM, _EwaldSum_real_::SPME_VAR& var, InteractRes& res) {
	int nm;
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < nSM; nm++) {
		CLUSTER_EFIELD_SPME_REAL_ImageSlab(cell, sm[nm].c, var, tres);
		res += tres;
	}
}

void PM_EFIELD_SPME_REAL_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, PMOLECULE *pm, int nPM, _EwaldSum_real_::SPME_VAR& var, InteractRes& res) {
	int nm;
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < nPM; nm++) {
		CLUSTER_EFIELD_SPME_REAL_ImageSlab(cell, pm[nm].c, var, tres);
		res += tres;
	}
}

void ImageClusters_EFIELD_SPME_REAL_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, ImageCluster *m, int n, _EwaldSum_real_::SPME_VAR& var, InteractRes& res) {
	int nm;
	res.reset();
	InteractRes tres;
	for (nm = 0; nm < n; nm++) {
		ImageCluster_EFIELD_SPME_REAL_ImageSlab(cell, (BASIC_CLUSTER*)(m + nm), var, tres);
		res += tres;
	}
}

void DumpEF_ImageClusters_EFIELD_SPME_REAL_ImageSlab(MMOL_MD_CELL &mdcell, ARRAY<ImageCluster> &cimg1, ARRAY<ImageCluster> &cimg2, bool bF) {
	int i, ia;
	BASIC_CLUSTER *pc = NULL, *pc1 = NULL, *pc2 = NULL;
	MATOM *pa = NULL, *pa1 = NULL, *pa2 = NULL;
	for (i = 0; i < mdcell.ccluster.n; i++) {
		pc = mdcell.ccluster.m[i]; pc1 = (BASIC_CLUSTER*)(cimg1.m + i); pc2 = (BASIC_CLUSTER*)(cimg2.m + i);
		for (ia = 0; ia < pc->nAtoms; ia++) {
			pa = pc->atom + ia; pa1 = pc1->atom + ia; pa2 = pc2->atom + ia;
			pa->E.v[0] += pa1->E.v[0] + pa2->E.v[0];
			pa->E.v[1] += pa1->E.v[1] + pa2->E.v[1];
			pa->E.v[2] -= pa1->E.v[2] + pa2->E.v[2];

			if (bF) {
				pa->F.v[0] += pa1->F.v[0] + pa2->F.v[0];
				pa->F.v[1] += pa1->F.v[1] + pa2->F.v[1];
				pa->F.v[2] -= pa1->F.v[2] + pa2->F.v[2];
			}
		}
	}
}

}

namespace _spme_2d_ {
using namespace _EwaldSum_real_;
using namespace _EwaldSum_2d_;

void Kickoff_SelfImageInteraction_Qatoms(VSPME_ImageSlab<_VQatom> &vspme, int na1, int na2, _spme_2d_::SPME_VAR &var, InteractRes &res) {
	int ia, ja, natoms = vspme.vspme.va.n;
	_VQatom *patom = NULL, *patom1 = NULL;
	double Utotal = 0, U = 0, R, R2;
	VECTOR3 E, F, dr;
	bool bEwald = var.spme_3d_var.esv1.bEwald;
	var.spme_3d_var.esv1.bEwald = false;
	_EwaldSum_real_::EwaldSumRealVars0 *esv = (_EwaldSum_real_::EwaldSumRealVars0*)(&(var.spme_3d_var.esv1));

	res.reset();

	for (ia = na1; ia <= na2; ia++) {
		if (ia >= natoms) break;
		ja = natoms + ia;
		if (ja < vspme.Ivspme1.va.n) {
			patom = vspme.Ivspme1.va.m + ia;
			patom1 = vspme.Ivspme1.va.m + ja;
			dr.v[2] = patom->r->v[2] - patom1->r->v[2];
			R = FABS(dr.v[2]); R2 = R * R;
			esv->init_vars(dr, R, R2);
			MTP_Interact(*esv, *(patom->q), *(patom1->q), E, F, U, false, 0);
			Utotal += U;
			patom->E_recp.v[0] -= E.v[0];
			patom->E_recp.v[1] -= E.v[1];
			patom->E_recp.v[2] -= E.v[2];
			patom->F_recp.v[0] -= F.v[0];
			patom->F_recp.v[1] -= F.v[1];
			patom->F_recp.v[2] -= F.v[2];

			if (var.bVirial) {
				if (var.bST_diagonal) Accumulate_Virial_Diagonal(res, F.v, dr.v);
				else Accumulate_Virial_Trace(res, F.v, dr.v);
			}
		}
		if (ja < vspme.Ivspme2.va.n) {
			patom = vspme.Ivspme2.va.m + ia;
			patom1 = vspme.Ivspme2.va.m + ja;
			dr.v[2] = patom->r->v[2] - patom1->r->v[2];
			R = FABS(dr.v[2]); R2 = R * R;
			esv->init_vars(dr, R, R2);
			MTP_Interact(*esv, *(patom->q), *(patom1->q), E, F, U, false, 0);
			Utotal += U;
			patom->E_recp.v[0] -= E.v[0];
			patom->E_recp.v[1] -= E.v[1];
			patom->E_recp.v[2] -= E.v[2];
			patom->F_recp.v[0] -= F.v[0];
			patom->F_recp.v[1] -= F.v[1];
			patom->F_recp.v[2] -= F.v[2];

			if (var.bVirial) {
				if (var.bST_diagonal) Accumulate_Virial_Diagonal(res, F.v, dr.v);
				else Accumulate_Virial_Trace(res, F.v, dr.v);
			}
		}
	}
	res.U = Utotal * eUnit_estat_kT;
	if (var.bVirial) {
		res.STtr *= fUnit_estat_mass;
		if (var.bST_diagonal) {
			res.STxx *= fUnit_estat_mass;
			res.STyy *= fUnit_estat_mass;
			res.STzz *= fUnit_estat_mass;
		}
	}
	esv->bEwald = bEwald;
}

void Kickoff_SelfImageInteraction_MUatoms(VSPME_ImageSlab<_VMUatom> &vspme, int na1, int na2, _spme_2d_::SPME_VAR &var, InteractRes &res) {
	int ia, ja, natoms = vspme.vspme.va.n;
	_VMUatom *patom = NULL, *patom1 = NULL;
	double Utotal = 0, U = 0, R, R2;
	VECTOR3 E, F, dr;
	bool bEwald = var.spme_3d_var.esv1.bEwald;
	var.spme_3d_var.esv1.bEwald = false;
	_EwaldSum_real_::EwaldSumRealVars1 *esv = &(var.spme_3d_var.esv1);

	res.reset();

	for (ia = na1; ia <= na2; ia++) {
		if (ia >= natoms) break;
		ja = natoms + ia;
		if (ja < vspme.Ivspme1.va.n) {
			patom = vspme.Ivspme1.va.m + ia;
			patom1 = vspme.Ivspme1.va.m + ja;
			dr.v[2] = patom->r->v[2] - patom1->r->v[2];
			R = FABS(dr.v[2]); R2 = R * R;
			esv->init_vars(dr, R, R2);
			MTP_Interact(*esv, *(patom->q), patom->mu->v, *(patom1->q), patom1->mu->v, E, F, U, false, 0);
			Utotal += U;
			patom->E_recp.v[0] -= E.v[0];
			patom->E_recp.v[1] -= E.v[1];
			patom->E_recp.v[2] -= E.v[2];
			patom->F_recp.v[0] -= F.v[0];
			patom->F_recp.v[1] -= F.v[1];
			patom->F_recp.v[2] -= F.v[2];

			if (var.bVirial) {
				if (var.bST_diagonal) Accumulate_Virial_Diagonal(res, F.v, dr.v);
				else Accumulate_Virial_Trace(res, F.v, dr.v);
			}
		}
		if (ja < vspme.Ivspme2.va.n) {
			patom = vspme.Ivspme2.va.m + ia;
			patom1 = vspme.Ivspme2.va.m + ja;
			dr.v[2] = patom->r->v[2] - patom1->r->v[2];
			R = FABS(dr.v[2]); R2 = R * R;
			esv->init_vars(dr, R, R2);
			MTP_Interact(*esv, *(patom->q), *(patom1->q), E, F, U, false, 0);
			Utotal += U;
			patom->E_recp.v[0] -= E.v[0];
			patom->E_recp.v[1] -= E.v[1];
			patom->E_recp.v[2] -= E.v[2];
			patom->F_recp.v[0] -= F.v[0];
			patom->F_recp.v[1] -= F.v[1];
			patom->F_recp.v[2] -= F.v[2];

			if (var.bVirial) {
				if (var.bST_diagonal) Accumulate_Virial_Diagonal(res, F.v, dr.v);
				else Accumulate_Virial_Trace(res, F.v, dr.v);
			}
		}
	}
	res.U = Utotal * eUnit_estat_kT;
	if (var.bVirial) {
		res.STtr *= fUnit_estat_mass;
		if (var.bST_diagonal) {
			res.STxx *= fUnit_estat_mass;
			res.STyy *= fUnit_estat_mass;
			res.STzz *= fUnit_estat_mass;
		}
	}
	
	esv->bEwald = bEwald;
}

void Kickoff_SelfImageInteraction(VSPME_ImageSlab<_VQatom> &vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	MRunSys1< VSPME_ImageSlab<_VQatom>, SPME_VAR, InteractRes>((void*)(&Kickoff_SelfImageInteraction_Qatoms), vspme, 0, vspme.vspme.va.n - 1, nThreads, var, res);
}

void Kickoff_SelfImageInteraction(VSPME_ImageSlab<_VMUatom> &vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	MRunSys1< VSPME_ImageSlab<_VMUatom>, SPME_VAR, InteractRes>((void*)(&Kickoff_SelfImageInteraction_MUatoms), vspme, 0, vspme.vspme.va.n - 1, nThreads, var, res);
}


/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(MMOL_MD_CELL *mdcell, ARRAY<ImageCluster> &cimg1, ARRAY<ImageCluster> &cimg2, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME_ImageSlab<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field

	// Important! cmm_cell has atoms of real, image1 and image2
	res.reset();
	InteractRes tres, res_recp;

	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}

	MultiMolInteraction<ImageCluster, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&ImageClusters_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, cimg1.m, cimg1.n, nThreads, var.spme_3d_var, tres);
	res += tres;
	MultiMolInteraction<ImageCluster, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&ImageClusters_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, cimg2.m, cimg2.n, nThreads, var.spme_3d_var, tres);
	res += tres;
	DumpEF_ImageClusters_EFIELD_SPME_REAL_ImageSlab(*mdcell, cimg1, cimg2, !var.spme_3d_var.bEfieldOnly);

	// do we need this term if the system is not neutral ?
	//if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;

	//sprintf(errmsg, "Real : %f kT", res.U); show_log(errmsg, true);
	// real + image1
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	// accumulate the electric field from the image part
	//vspme.Ivspme1.reset_local_EF(); // reset the variable for E & F
	vspme.Ivspme1.init_fft_plans();
	init_normalized_coordinates(vspme.Ivspme1, nThreads);
	cal_Q_spme(vspme.Ivspme1, nThreads); // SPME Q matrix
	double U = vspme.Ivspme1.cal_U_recp(); //calculate the energy from the image part 
	res_recp.U += U;
	cal_E_recp(vspme.Ivspme1, 0, vspme.Ivspme1.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.Ivspme1, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.Ivspme1.cal_FTQ();
		vspme.Ivspme1.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx += vspme.Ivspme1.STxx; 
			res_recp.STyy += vspme.Ivspme1.STyy; 
			res_recp.STzz += vspme.Ivspme1.STzz;
		}
		res_recp.STtr += vspme.Ivspme1.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of whole sys: %f kT", U); show_log(errmsg, true);

	tres.reset();
	EwaldSum_2D_K0(vspme.Ivspme1, nThreads, var.K0var, tres);
	res_recp += tres;

	//sprintf(errmsg, "Recip [k=0] of whole sys: %f kT", tres.U); show_log(errmsg, true);


	// real + image2
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	// accumulate the electric field from the image part
	//vspme.Ivspme2.reset_local_EF(); // reset the variable for E & F
	vspme.Ivspme2.init_fft_plans();
	init_normalized_coordinates(vspme.Ivspme2, nThreads);
	cal_Q_spme(vspme.Ivspme2, nThreads); // SPME Q matrix
	U = vspme.Ivspme2.cal_U_recp(); //calculate the energy from the image part 
	res_recp.U += U;
	cal_E_recp(vspme.Ivspme2, 0, vspme.Ivspme2.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.Ivspme2, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.Ivspme2.cal_FTQ();
		vspme.Ivspme2.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx += vspme.Ivspme2.STxx; 
			res_recp.STyy += vspme.Ivspme2.STyy; 
			res_recp.STzz += vspme.Ivspme2.STzz;
		}
		res_recp.STtr += vspme.Ivspme2.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of whole sys: %f kT", U); show_log(errmsg, true);

	tres.reset();
	EwaldSum_2D_K0(vspme.Ivspme2, nThreads, var.K0var, tres);
	res_recp += tres;

	//sprintf(errmsg, "Recip [k=0] of whole sys: %f kT", tres.U); show_log(errmsg, true);

	// interaction between image1 charges
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	vspme.img1_vspme.init_fft_plans();
	init_normalized_coordinates(vspme.img1_vspme, nThreads);
	cal_Q_spme(vspme.img1_vspme, nThreads); // SPME Q matrix
	U = vspme.img1_vspme.cal_U_recp(); //calculate the energy from the image part 
	res_recp.U -= U; // energies of image1
	cal_E_recp(vspme.img1_vspme, 0, vspme.img1_vspme.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.img1_vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.img1_vspme.cal_FTQ();
		vspme.img1_vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx -= vspme.img1_vspme.STxx; 
			res_recp.STyy -= vspme.img1_vspme.STyy; 
			res_recp.STzz -= vspme.img1_vspme.STzz;
		}
		res_recp.STtr -= vspme.img1_vspme.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of slab: %f kT", U); show_log(errmsg, true);

	vspme.release_fftw_plan(); // make sure FFTW has no plan

	tres.reset();
	EwaldSum_2D_K0(vspme.img1_vspme, nThreads, var.K0var, tres);
	res_recp -= tres;


	// interaction between image2 charges
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	vspme.img2_vspme.init_fft_plans();
	init_normalized_coordinates(vspme.img2_vspme, nThreads);
	cal_Q_spme(vspme.img2_vspme, nThreads); // SPME Q matrix
	U = vspme.img2_vspme.cal_U_recp(); //calculate the energy from the image part 
	res_recp.U -= U; // energies of image1
	cal_E_recp(vspme.img2_vspme, 0, vspme.img2_vspme.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.img2_vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.img2_vspme.cal_FTQ();
		vspme.img2_vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx -= vspme.img2_vspme.STxx; 
			res_recp.STyy -= vspme.img2_vspme.STyy; 
			res_recp.STzz -= vspme.img2_vspme.STzz;
		}
		res_recp.STtr -= vspme.img2_vspme.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of slab: %f kT", U); show_log(errmsg, true);

	vspme.release_fftw_plan(); // make sure FFTW has no plan

	tres.reset();
	EwaldSum_2D_K0(vspme.img2_vspme, nThreads, var.K0var, tres);
	res_recp -= tres;

	// interaction between real charges
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	vspme.vspme.init_fft_plans();
	init_normalized_coordinates(vspme.vspme, nThreads);
	cal_Q_spme(vspme.vspme, nThreads); // SPME Q matrix
	U = vspme.vspme.cal_U_recp(); //calculate the energy from the image part 
	res_recp.U -= U; // energies of image1 and image2
	cal_E_recp(vspme.vspme, 0, vspme.vspme.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.vspme.cal_FTQ();
		vspme.vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx -= vspme.vspme.STxx; 
			res_recp.STyy -= vspme.vspme.STyy; 
			res_recp.STzz -= vspme.vspme.STzz;
		}
		res_recp.STtr -= vspme.vspme.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of slab: %f kT", U); show_log(errmsg, true);

	vspme.release_fftw_plan(); // make sure FFTW has no plan

	tres.reset();
	EwaldSum_2D_K0(vspme.vspme, nThreads, var.K0var, tres);
	res_recp -= tres;


	//sprintf(errmsg, "Recip [k=0] of real charge: %f kT", tres.U); show_log(errmsg, true);

	//vspme.dump_EF(!var.K0var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
	bool bDumpF = !var.K0var.bEfieldOnly;
	if (nThreads < 2 || vspme.Ivspme1.n_real < N4PARALLEL) vspme.dump_EF(bDumpF); // transfer the calculated E & F to the variables in real atom
	else MRunSys2< VSPME_ImageSlab<_VQatom>, bool>((void*)(&ImageSlab_dump_EF_Q), vspme, 0, vspme.Ivspme1.n_real, nThreads, bDumpF);

//	res_recp.U *= 0.5;
//	if (var.bVirial) {
//		res_recp.STtr *= 0.5;
//		if (var.bST_diagonal) {
//			res_recp.STxx *= 0.5;
//			res_recp.STyy *= 0.5;
//			res_recp.STzz *= 0.5;
//		}
//	}
	res += res_recp;


	//sprintf(errmsg, "Interaction of charges with their image charges: %f kT", tres.U); show_log(errmsg, true);

	// virial coming from the surface dipole, U(V), from receiprcal part
	/*
	double *v1, *v2;
	double STxx = 0, STyy = 0, STzz = 0;
	if (!var.K0var.bEfieldOnly && var.bVirial) {
		v1 = cmm_cell->mu_surf.qr.v; v2 = cmm_cell->mu_surf.mu0.v;
		STxx = (v1[0] * v1[0] * 0.5 + 2 * v1[0] * v2[0] + 1.5 * v2[0] * v2[0]) * vspme.inv_V / 3 * ::fUnit_estat_mass;
		STyy = (v1[1] * v1[1] * 0.5 + 2 * v1[1] * v2[1] + 1.5 * v2[1] * v2[1]) * vspme.inv_V / 3 * ::fUnit_estat_mass;
		STzz = (v1[2] * v1[2] * 0.5 + 2 * v1[2] * v2[2] + 1.5 * v2[2] * v2[2]) * vspme.inv_V / 3 * ::fUnit_estat_mass;
		if (var.bST_diagonal) {
			res.STxx += STxx; res.STyy += STyy; res.STzz += STzz;
		}
		res.STtr += STxx + STyy + STzz;
	}
	*/
	//sprintf(errmsg, "After surface dipole correction, virial changes to %f", res.STtr); show_log(errmsg, true);

	//vspme.dump_EF(!var.K0var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom


	// interaction between charge with its image charge directly
	/*
	vspme.Ivspme1.reset_local_EF();
	vspme.Ivspme2.reset_local_EF();
	vspme.vspme.reset_local_EF();

	tres.reset();
	Kickoff_SelfImageInteraction(vspme, nThreads, var, tres);
	vspme.Ivspme1.dump_EF();
	vspme.Ivspme2.dump_EF();
	res -= tres;
	*/
}


/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
void SPME_Interact(MMOL_MD_CELL *mdcell, ARRAY<ImageCluster> &cimg1, ARRAY<ImageCluster> &cimg2, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME_ImageSlab<_VQatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	SPME_EF(mdcell, cimg1, cimg2, cmm_cell, vspme, nThreads, var, res);

	if (var.K0var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}


// with dipole
/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(MMOL_MD_CELL *mdcell, ARRAY<ImageCluster> &cimg1, ARRAY<ImageCluster> &cimg2, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME_ImageSlab<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres, res_recp;

	// Important! cmm_cell has atoms of real, image1 and image2
	if (mdcell->mm.m != NULL) {
		MultiMolInteraction<MMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&MM_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, mdcell->mm.m, mdcell->mm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->sm.m != NULL) {
		MultiMolInteraction<SMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&SM_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, mdcell->sm.m, mdcell->sm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}
	if (mdcell->pm.m != NULL) {
		MultiMolInteraction<PMOLECULE, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&PM_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, mdcell->pm.m, mdcell->pm.n, nThreads, var.spme_3d_var, tres);
		res += tres;
	}

	MultiMolInteraction<ImageCluster, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&ImageClusters_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, cimg1.m, cimg1.n, nThreads, var.spme_3d_var, tres);
	res += tres;
	MultiMolInteraction<ImageCluster, BASIC_CLUSTER, _EwaldSum_real_::SPME_VAR, InteractRes>((void*)(&ImageClusters_EFIELD_SPME_REAL_ImageSlab), *cmm_cell, cimg2.m, cimg2.n, nThreads, var.spme_3d_var, tres);
	res += tres;
	DumpEF_ImageClusters_EFIELD_SPME_REAL_ImageSlab(*mdcell, cimg1, cimg2, !var.spme_3d_var.bEfieldOnly);

	// do we need this term if the system is not neutral ?
	//if (cmm_cell->mu_surf.q != 0) res.U -= 1.5 * var.esv1.inv_3V / (var.esv1.kappa * var.esv1.kappa) * cmm_cell->mu_surf.q * cmm_cell->mu_surf.q * eUnit_estat_kT;

	// real + image1
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	// accumulate the electric field from the image part
	//vspme.Ivspme1.reset_local_EF(); // reset the variable for E & F
	vspme.Ivspme1.init_fft_plans();
	//init_normalized_coordinates(vspme.Ivspme1, nThreads);
	cal_Q_spme(vspme.Ivspme1, nThreads); // SPME Q matrix
	double U = vspme.Ivspme1.cal_U_recp(); //calculate the energy from the image part 
	res_recp.U += U;
	cal_E_recp(vspme.Ivspme1, 0, vspme.Ivspme1.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.Ivspme1, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.Ivspme1.cal_FTQ();
		vspme.Ivspme1.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx += vspme.Ivspme1.STxx; 
			res_recp.STyy += vspme.Ivspme1.STyy; 
			res_recp.STzz += vspme.Ivspme1.STzz;
		}
		res_recp.STtr += vspme.Ivspme1.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of whole sys: %f kT", U); show_log(errmsg, true);

	tres.reset();
	EwaldSum_2D_K0(vspme.Ivspme1, nThreads, var.K0var, tres);
	res_recp += tres;

	//sprintf(errmsg, "Recip [k=0] of whole sys: %f kT", tres.U); show_log(errmsg, true);


	// real + image2
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	// accumulate the electric field from the image part
	//vspme.Ivspme2.reset_local_EF(); // reset the variable for E & F
	vspme.Ivspme2.init_fft_plans();
	//init_normalized_coordinates(vspme.Ivspme2, nThreads);
	cal_Q_spme(vspme.Ivspme2, nThreads); // SPME Q matrix
	U = vspme.Ivspme2.cal_U_recp(); //calculate the energy from the image part 
	res_recp.U += U;
	cal_E_recp(vspme.Ivspme2, 0, vspme.Ivspme2.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.Ivspme2, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.Ivspme2.cal_FTQ();
		vspme.Ivspme2.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx += vspme.Ivspme2.STxx; 
			res_recp.STyy += vspme.Ivspme2.STyy; 
			res_recp.STzz += vspme.Ivspme2.STzz;
		}
		res_recp.STtr += vspme.Ivspme2.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of whole sys: %f kT", U); show_log(errmsg, true);

	tres.reset();
	EwaldSum_2D_K0(vspme.Ivspme2, nThreads, var.K0var, tres);
	res_recp += tres;

	//sprintf(errmsg, "Recip [k=0] of whole sys: %f kT", tres.U); show_log(errmsg, true);


	// interaction between image1 charges
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	vspme.img1_vspme.init_fft_plans();
	//init_normalized_coordinates(vspme.img1_vspme, nThreads);
	cal_Q_spme(vspme.img1_vspme, nThreads); // SPME Q matrix
	U = vspme.img1_vspme.cal_U_recp(); //calculate the energy from the image part 
	res_recp.U -= U; // energies of image1
	cal_E_recp(vspme.img1_vspme, 0, vspme.img1_vspme.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.img1_vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.img1_vspme.cal_FTQ();
		vspme.img1_vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx -= vspme.img1_vspme.STxx; 
			res_recp.STyy -= vspme.img1_vspme.STyy; 
			res_recp.STzz -= vspme.img1_vspme.STzz;
		}
		res_recp.STtr -= vspme.img1_vspme.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of slab: %f kT", U); show_log(errmsg, true);

	vspme.release_fftw_plan(); // make sure FFTW has no plan

	tres.reset();
	EwaldSum_2D_K0(vspme.img1_vspme, nThreads, var.K0var, tres);
	res_recp -= tres;


	// interaction between image2 charges
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	vspme.img2_vspme.init_fft_plans();
	//init_normalized_coordinates(vspme.img2_vspme, nThreads);
	cal_Q_spme(vspme.img2_vspme, nThreads); // SPME Q matrix
	U = vspme.img2_vspme.cal_U_recp(); //calculate the energy from the image part 
	res_recp.U -= U; // energies of image1
	cal_E_recp(vspme.img2_vspme, 0, vspme.img2_vspme.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.img2_vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.img2_vspme.cal_FTQ();
		vspme.img2_vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx -= vspme.img2_vspme.STxx; 
			res_recp.STyy -= vspme.img2_vspme.STyy; 
			res_recp.STzz -= vspme.img2_vspme.STzz;
		}
		res_recp.STtr -= vspme.img2_vspme.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of slab: %f kT", U); show_log(errmsg, true);

	vspme.release_fftw_plan(); // make sure FFTW has no plan

	tres.reset();
	EwaldSum_2D_K0(vspme.img2_vspme, nThreads, var.K0var, tres);
	res_recp -= tres;

	// interaction between real charges, for energy correction only
	vspme.release_fftw_plan(); // make sure FFTW has no plan
	vspme.vspme.init_fft_plans();
	//init_normalized_coordinates(vspme.vspme, nThreads);
	cal_Q_spme(vspme.vspme, nThreads); // SPME Q matrix
	U = vspme.vspme.cal_U_recp(); //calculate the energy from the image part 
	//res_recp.U -= U * (vspme.f1 * vspme.f1 + vspme.f2 * vspme.f2); // energies of image1 and image2
	res_recp.U -= U;
	cal_E_recp(vspme.vspme, 0, vspme.vspme.va.n - 1, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		cal_TQ_spme(vspme.vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.vspme.cal_FTQ();
		vspme.vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {
			res_recp.STxx -= vspme.vspme.STxx; 
			res_recp.STyy -= vspme.vspme.STyy; 
			res_recp.STzz -= vspme.vspme.STzz;
		}
		res_recp.STtr -= vspme.vspme.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0] of slab: %f kT", U); show_log(errmsg, true);

	vspme.release_fftw_plan(); // make sure FFTW has no plan

	tres.reset();
	EwaldSum_2D_K0(vspme.vspme, nThreads, var.K0var, tres);
	res_recp -= tres;

	//sprintf(errmsg, "Recip [k=0] of real charge: %f kT", tres.U); show_log(errmsg, true);

//	res_recp.U *= 0.5;
//	if (var.bVirial) {
//		res_recp.STtr *= 0.5;
//		if (var.bST_diagonal) {
//			res_recp.STxx *= 0.5;
//			res_recp.STyy *= 0.5;
//			res_recp.STzz *= 0.5;
//		}
//	}
	res += res_recp;



	//sprintf(errmsg, "Interaction of charges with their image charges: %f kT", tres.U); show_log(errmsg, true);

	// virial coming from the surface dipole, U(V), from receiprical part
	/*
	double *v1, *v2;
	double STxx = 0, STyy = 0, STzz = 0;
	if (!var.K0var.bEfieldOnly && var.bVirial) {
		v1 = cmm_cell->mu_surf.qr.v; v2 = cmm_cell->mu_surf.mu0.v;
		STxx = (v1[0] * v1[0] * 0.5 + 2 * v1[0] * v2[0] + 1.5 * v2[0] * v2[0]) * vspme.inv_V / 3 * ::fUnit_estat_mass;
		STyy = (v1[1] * v1[1] * 0.5 + 2 * v1[1] * v2[1] + 1.5 * v2[1] * v2[1]) * vspme.inv_V / 3 * ::fUnit_estat_mass;
		STzz = (v1[2] * v1[2] * 0.5 + 2 * v1[2] * v2[2] + 1.5 * v2[2] * v2[2]) * vspme.inv_V / 3 * ::fUnit_estat_mass;
		if (var.bST_diagonal) {
			res.STxx += STxx; res.STyy += STyy; res.STzz += STzz;
		}
		res.STtr += STxx + STyy + STzz;
	}
	*/

	/*
	{
		int iatom = 10;
		char msg[256] = "\0";
		sprintf(msg, "Total E : %f", U); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);

		_VE_CELL<_VMUatom> *pvcell = (_VE_CELL<_VMUatom> *)(&vspme);
		double U1 = EwaldSum_Recip(*pvcell);
		sprintf(msg, "Total E : %f", U1); show_log(msg, true);
		sprintf(msg, "%f %f %f", vspme.va.m[iatom].E_recp.v[0], vspme.va.m[iatom].E_recp.v[1], vspme.va.m[iatom].E_recp.v[2]); 
		show_log(msg, true);
	}
	*/

	//vspme.dump_EF(!var.K0var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
	bool bDumpF = !var.K0var.bEfieldOnly;
	if (nThreads < 2 || vspme.Ivspme1.n_real < N4PARALLEL) vspme.dump_EF(bDumpF); // transfer the calculated E & F to the variables in real atom
	else MRunSys2< VSPME_ImageSlab<_VMUatom>, bool>((void*)(&ImageSlab_dump_EF_MU), vspme, 0, vspme.Ivspme1.n_real, nThreads, bDumpF);


	// interaction between charge with its image charge directly
	/*
	vspme.Ivspme1.reset_local_EF();
	vspme.Ivspme2.reset_local_EF();
	vspme.vspme.reset_local_EF();

	tres.reset();
	Kickoff_SelfImageInteraction(vspme, nThreads, var, tres);
	vspme.Ivspme1.dump_EF();
	vspme.Ivspme2.dump_EF();
	res -= tres;
	*/
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
void SPME_Interact(MMOL_MD_CELL *mdcell, ARRAY<ImageCluster> &img1, float f1, ARRAY<ImageCluster> &img2, float f2, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME_ImageSlab<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	init_normalized_coordinates(vspme.Ivspme1, nThreads); // vspme has to be initialized before SPME_EF for atom
	init_normalized_coordinates(vspme.Ivspme2, nThreads); // vspme has to be initialized before SPME_EF for atom
	init_normalized_coordinates(vspme.img1_vspme, nThreads); // vspme has to be initialized before SPME_EF for atom
	init_normalized_coordinates(vspme.img2_vspme, nThreads); // vspme has to be initialized before SPME_EF for atom
	init_normalized_coordinates(vspme.vspme, nThreads); // vspme has to be initialized before SPME_EF for atom

	UpdateDipole_ImageCluster(img1, f1, img2, f2, nThreads);

	SPME_EF(mdcell, img1, img2, cmm_cell, vspme, nThreads, var, res);

	if (var.K0var.bEfieldOnly) return;

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
bool Polarize_SPME_Interact(MMOL_MD_CELL *mdcell, ARRAY<ImageCluster> &img1, float f1, ARRAY<ImageCluster> &img2, float f2, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME_ImageSlab<_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field

	bool bEfieldOnly = var.K0var.bEfieldOnly;
	res.reset();
	InteractRes tres;

	init_normalized_coordinates(vspme.Ivspme1, nThreads); // vspme has to be initialized before SPME_EF for atom with dipole
	init_normalized_coordinates(vspme.Ivspme2, nThreads); // vspme has to be initialized before SPME_EF for atom with dipole
	init_normalized_coordinates(vspme.img1_vspme, nThreads); // vspme has to be initialized before SPME_EF for atom
	init_normalized_coordinates(vspme.img2_vspme, nThreads); // vspme has to be initialized before SPME_EF for atom
	init_normalized_coordinates(vspme.vspme, nThreads); // vspme has to be initialized before SPME_EF for atom

	double dmu_max = 0, vt = 0, muConvg = ::muConvg; // eA
	int iloop = 0, max_loops = ::nmax_loops_polarize;
	bool status = false;
	var.K0var.bEfieldOnly = true;
	var.spme_3d_var.bEfieldOnly = true;
	while (bPolarize && iloop < max_loops) {
		SPME_EF(mdcell, img1, img2, cmm_cell, vspme, nThreads, var, tres); // calculate the dipole only

		// add the electric field from the atom inside the cluster
		//if (mdcell->sm.n > 0) MOperate<SMOLECULE>((void*)(&AtomicEFieldCorrect_SM), mdcell->sm.m, mdcell->sm.n, nThreads);

		// polarization
		if (mdcell->mm.n > 0) MOperate<MMOLECULE>((void*)(&polarize_MM), mdcell->mm.m, mdcell->mm.n, nThreads);
		if (mdcell->sm.n > 0) MOperate<SMOLECULE>((void*)(&polarize_SM), mdcell->sm.m, mdcell->sm.n, nThreads);
		if (mdcell->pm.n > 0) MOperate<PMOLECULE>((void*)(&polarize_PM), mdcell->pm.m, mdcell->pm.n, nThreads);
		cal_dipole(*mdcell, nThreads, mdcell->mu_surf); // accumulate dipole of each atom, and the whole mdcell
		memcpy(&(cmm_cell->mu_surf), &(mdcell->mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
		UpdateDipole_ImageCluster(img1, f1, img2, f2, nThreads);
		dmu_max = 0; vt = 0;
		if (mdcell->ind_mu_mm.n > 0) MOperate2<MOL_POL<MMOLECULE>, double>((void*)(&MM_maxdiff_backup_polarization), mdcell->ind_mu_mm.m, mdcell->ind_mu_mm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_sm.n > 0) MOperate2<MOL_POL<SMOLECULE>, double>((void*)(&SM_maxdiff_backup_polarization), mdcell->ind_mu_sm.m, mdcell->ind_mu_sm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell->ind_mu_pm.n > 0) MOperate2<MOL_POL<PMOLECULE>, double>((void*)(&PM_maxdiff_backup_polarization), mdcell->ind_mu_pm.m, mdcell->ind_mu_pm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;

		iloop++;
		if (dmu_max <= muConvg && iloop > 0) {status = true; break;}
	}

	//char msg[256] = "\0";
	//if (status) sprintf(msg, "Convergent polarized dipole is obtainted with %d loops", iloop); 
	//else sprintf(msg, "No convergent polarized dipole is obtained with %d loops", iloop);
	//show_log(msg, true);

	var.K0var.bEfieldOnly = bEfieldOnly;
	var.spme_3d_var.bEfieldOnly = bEfieldOnly;

	SPME_EF(mdcell, img1, img2, cmm_cell, vspme, nThreads, var, res); // calculate the dipole and force

	// add the force and torque from atom to cluster
	//if (mdcell->mm.m != NULL) {
	//	MOperate<MMOLECULE>(MM_TOTAL_FORCE, mdcell->mm.m, mdcell->mm.n, nThreads);
	//}
	//if (mdcell->sm.m != NULL) {
	//	MOperate<SMOLECULE>(SM_TOTAL_FORCE, mdcell->sm.m, mdcell->sm.n, nThreads);
	//}
	//if (mdcell->pm.m != NULL) {
	//	MOperate<PMOLECULE>(PM_TOTAL_FORCE, mdcell->pm.m, mdcell->pm.n, nThreads);
	//}

	return status;
}


// for 2d system to calculate electrostatic interaction with explicitly
void ExplicitEInteract_ImageSlab(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOL_MD_CELL &mdcell, ARRAY<MATOM*> &img1, ARRAY<MATOM*> &img2, E2dInteractVar& evar, InteractRes &res) {
	int im, ic, ia;
	res.reset();
	InteractRes tres;
	for (im = 0; im < mdcell.mm.n; im++) {
		for (ic = 0; ic < mdcell.mm.m[im].nCluster; ic++) {
			tres.reset();
			ExplicitEInteract_2d_ImageSlab(cell, mdcell.mm.m[im].cluster + ic, evar, tres);
			res += tres;
		}
	}
	for (im = 0; im < mdcell.sm.n; im++) {
		tres.reset();
		ExplicitEInteract_2d_ImageSlab(cell, mdcell.sm.m[im].c, evar, tres);
		res += tres;
	}
	for (im = 0; im < mdcell.pm.n; im++) {
		tres.reset();
		ExplicitEInteract_2d_ImageSlab(cell, mdcell.pm.m[im].c, evar, tres);
		res += tres;
	}

	// we have to remove the interaction of atoms with their image atoms
	/*
	_EwaldSum_real_::EwaldSumRealVars1 esv;
	_EwaldSum_real_::EwaldSumRealVars0 *esv0 = (_EwaldSum_real_::EwaldSumRealVars0 *)(&esv);
	esv.bEwald = false;
	VECTOR3 dr, E, F;
	double U, R, R2;
	MATOM *patom, *patom1;
	tres.reset();
	for (ia = 0; ia < mdcell.atom.n; ia++) {
		patom = mdcell.atom.m[ia];
		patom1 = img1.m[ia];
		dr.v[0] = patom->rg.v[0] - patom1->rg.v[0];
		dr.v[1] = patom->rg.v[1] - patom1->rg.v[1];
		dr.v[2] = patom->rg.v[2] - patom1->rg.v[2];
		V3ABS2(dr, R2); R = sqrt(R2);
		if (patom->mMP == 0) {
			esv0->init_vars(dr, R, R2);
			MTP_Interact(*esv0, patom->c, patom1->c, E, F, U, false, 0);
		}
		else if (patom->mMP == 1) {
			esv.init_vars(dr, R, R2);
			MTP_Interact(esv, patom->c, patom->mu.v, patom1->c, patom1->mu.v, E, F, U, false, 0);
		}
		tres.U += U * eUnit_estat_kT;
		patom->E.v[0] -= E.v[0];
		patom->E.v[1] -= E.v[1];
		patom->E.v[2] -= E.v[2];
		F.v[0] *= fUnit_estat_mass;
		F.v[1] *= fUnit_estat_mass;
		F.v[2] *= fUnit_estat_mass;
		patom->F.v[0] -= F.v[0];
		patom->F.v[1] -= F.v[1];
		patom->F.v[2] -= F.v[2];
		if (evar.bVirial) {
			if (evar.bST_diagonal) Accumulate_Virial_Diagonal(tres, F.v, dr.v);
			else Accumulate_Virial_Trace(tres, F.v, dr.v);
		}

		patom1 = img2.m[ia];
		dr.v[0] = patom->rg.v[0] - patom1->rg.v[0];
		dr.v[1] = patom->rg.v[1] - patom1->rg.v[1];
		dr.v[2] = patom->rg.v[2] - patom1->rg.v[2];
		V3ABS2(dr, R2); R = sqrt(R2);
		if (patom->mMP == 0) {
			esv0->init_vars(dr, R, R2);
			MTP_Interact(*esv0, patom->c, patom1->c, E, F, U, false, 0);
		}
		else if (patom->mMP == 1) {
			esv.init_vars(dr, R, R2);
			MTP_Interact(esv, patom->c, patom->mu.v, patom1->c, patom1->mu.v, E, F, U, false, 0);
		}
		tres.U += U * eUnit_estat_kT;
		patom->E.v[0] -= E.v[0];
		patom->E.v[1] -= E.v[1];
		patom->E.v[2] -= E.v[2];
		F.v[0] *= fUnit_estat_mass;
		F.v[1] *= fUnit_estat_mass;
		F.v[2] *= fUnit_estat_mass;
		patom->F.v[0] -= F.v[0];
		patom->F.v[1] -= F.v[1];
		patom->F.v[2] -= F.v[2];
		if (evar.bVirial) {
			if (evar.bST_diagonal) Accumulate_Virial_Diagonal(tres, F.v, dr.v);
			else Accumulate_Virial_Trace(tres, F.v, dr.v);
		}
	}
	res -= tres;
	*/
}

void CreateImageCluster(MMOL_MD_CELL &mdcell, float z, float f, ARRAY<ImageCluster> &c_img, ARRAY<MATOM*> &a_img, int mIndx_start) {
	if (mdcell.catom.n == 0 || mdcell.ccluster.n == 0) mdcell.check_eneutral();

	int nc = 0, ic, ia;
	c_img.set_array(mdcell.ccluster.n);

	float z2 = z + z;
	
	BASIC_CLUSTER *pc;
	nc = 0;
	for (ic = 0; ic < mdcell.ccluster.n; ic++) {
		pc = mdcell.ccluster.m[ic];
		cp((BASIC_CLUSTER*)(c_img.m + nc), pc);
		c_img.m[nc].s = pc;
		c_img.m[nc].vp = &(c_img.m[nc].r); c_img.m[nc].Op = &(c_img.m[nc].r);
		c_img.m[nc].vp->v[0] = pc->vp->v[0];
		c_img.m[nc].vp->v[1] = pc->vp->v[1];
		c_img.m[nc].vp->v[2] = z2 - (pc->vp->v[2]);
		V3zero(c_img.m[nc].dr_cmm) c_img.m[nc].bCMMShift = false;
		c_img.m[nc].eNeutral = pc->eNeutral;
		c_img.m[nc].mIndx = mIndx_start + nc;
		c_img.m[nc].cID = -1;
		c_img.m[nc].mType = -1;
		c_img.m[nc].fatom.set_array(c_img.m[nc].nAtoms);
		for (ia = 0; ia < c_img.m[nc].nAtoms; ia++) {
			c_img.m[nc].fatom.m[ia] = c_img.m[nc].atom + ia;
		}
		for (ia = 0; ia < c_img.m[nc].nAtoms; ia++) {
			c_img.m[nc].atom[ia].rg.v[0] = pc->atom[ia].rg.v[0];
			c_img.m[nc].atom[ia].rg.v[1] = pc->atom[ia].rg.v[1];
			c_img.m[nc].atom[ia].rg.v[2] = z2 - pc->atom[ia].rg.v[2];

			c_img.m[nc].atom[ia].r.v[0] = c_img.m[nc].atom[ia].rg.v[0];
			c_img.m[nc].atom[ia].r.v[1] = c_img.m[nc].atom[ia].rg.v[1];
			c_img.m[nc].atom[ia].r.v[2] = c_img.m[nc].atom[ia].rg.v[2];

			c_img.m[nc].atom[ia].c = pc->atom[ia].c * f;
			c_img.m[nc].atom[ia].c0 = pc->atom[ia].c0 * f;
			c_img.m[nc].atom[ia].eNeutral = pc->atom[ia].eNeutral;
			c_img.m[nc].atom[ia].par = NULL; // NO LJ interaction
			c_img.m[nc].atom[ia].aindx = -1; // image atom
			c_img.m[nc].atom[ia].mMP = pc->atom[ia].mMP;

			if (pc->atom[ia].mMP > 0) {
				c_img.m[nc].atom[ia].mu.v[0] = pc->atom[ia].mu.v[0] * f;
				c_img.m[nc].atom[ia].mu.v[1] = pc->atom[ia].mu.v[1] * f;
				c_img.m[nc].atom[ia].mu.v[2] = -f * pc->atom[ia].mu.v[2]; // image along z axis
			}
		}
		nc++;
	}

	int n = 0;
	for (ic = 0; ic < c_img.n; ic++) n += c_img.m[ic].nAtoms;
	a_img.set_array(n);
	n = 0;
	for (ic = 0; ic < c_img.n; ic++) {
		for (ia = 0; ia < c_img.m[ic].nAtoms; ia++) {
			a_img.m[n] = c_img.m[ic].atom + ia;
			n++;
		}
	}
}

void Update_ImageCluster_Pos(ImageCluster* cimg, int n, float &z) {
	int ic, ia;

	float z2 = z + z;
	
	BASIC_CLUSTER *pc;
	for (ic = 0; ic < n; ic++) {
		pc = cimg[ic].s; if (pc == NULL) continue;
		cimg[ic].vp->v[0] = pc->vp->v[0];
		cimg[ic].vp->v[1] = pc->vp->v[1];
		cimg[ic].vp->v[2] = z2 - (pc->vp->v[2]);
		for (ia = 0; ia < cimg[ic].nAtoms; ia++) {
			cimg[ic].atom[ia].rg.v[0] = pc->atom[ia].rg.v[0];
			cimg[ic].atom[ia].rg.v[1] = pc->atom[ia].rg.v[1];
			cimg[ic].atom[ia].rg.v[2] = z2 - pc->atom[ia].rg.v[2];

			cimg[ic].atom[ia].r.v[0] = cimg[ic].atom[ia].rg.v[0];
			cimg[ic].atom[ia].r.v[1] = cimg[ic].atom[ia].rg.v[1];
			cimg[ic].atom[ia].r.v[2] = cimg[ic].atom[ia].rg.v[2];
		}
	}
}

void UpdatePos_ImageCluster(ARRAY<ImageCluster> &img1, float z1, ARRAY<ImageCluster> &img2, float z2, int nThreads) {
	MOperate1<ImageCluster, float>((void*)(&Update_ImageCluster_Pos), img1.m, img1.n, z1, nThreads);
	MOperate1<ImageCluster, float>((void*)(&Update_ImageCluster_Pos), img2.m, img2.n, z2, nThreads);
}

void Update_ImageCluster_dipole(ImageCluster* cimg, int n, float &f) {
	int ic, ia;

	BASIC_CLUSTER *pc;
	for (ic = 0; ic < n; ic++) {
		pc = cimg[ic].s; if (pc == NULL) continue;
		for (ia = 0; ia < cimg[ic].nAtoms; ia++) {
			cimg[ic].atom[ia].mMP = pc->atom[ia].mMP;
			if (pc->atom[ia].mMP > 0) {
				cimg[ic].atom[ia].mu.v[0] = pc->atom[ia].mu.v[0] * f;
				cimg[ic].atom[ia].mu.v[1] = pc->atom[ia].mu.v[1] * f;
				cimg[ic].atom[ia].mu.v[2] = -f * pc->atom[ia].mu.v[2]; // image along z axis
			}
		}
	}
}


void UpdateDipole_ImageCluster(ARRAY<ImageCluster> &img1, float f1, ARRAY<ImageCluster> &img2, float f2, int nThreads) {
	MOperate1<ImageCluster, float>((void*)(&Update_ImageCluster_dipole), img1.m, img1.n, f1, nThreads);
	MOperate1<ImageCluster, float>((void*)(&Update_ImageCluster_dipole), img2.m, img2.n, f2, nThreads);
}


void cmm_init_distribute_cluster(ImageCluster* c, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain) {
	if (reset_cell_chain) cell.reset_cell_chain();

	int nm = 0, nc = 0;
	int ni = 0, nj = 0, nk = 0;
	VECTOR<3> dr, rc;
	int vi = 0;

	BASIC_CLUSTER* pc = NULL;
	for (nm = 0; nm < n; nm++) {
		pc = (BASIC_CLUSTER*)(c + nm);
		//V3plusV3((*(pc->vp)), pc->dr_cmm, rc)  // use vp position
		V32V3((*(pc->vp)), rc)

		CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
		//cell_n(cell, ni, nj, nk, nc)
		//if (nc >= cell.nCells) nc = cell.nCells - 1;
		//cell.cell[nc].attach(pc);
		cell.bcell.m[ni].m[nj].m[nk].attach(pc);
	}
}

} // _spme_2d_


void test_SPME_EwaldSum_ImageSlab() {
	//int nThreads = MAX_THREADS;
	int nThreads = 2;

	using namespace _cmm_2d_;

	extern bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol);
	extern bool cp(SMOLECULE *d, SMOLECULE *s, bool setup = true);
	extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, _EwaldSum_real_::SURF_DIPOLE &mu_surf);
	extern void init_LJ_Pars();
	extern void construct_LJ12_6_prof(int mode);
	extern double cal_charge(MMOL_MD_CELL &mmcell);

	int imol = 0, ia;
	SMOLECULE cw;
	VECTOR3 dr, axis;
	axis.v[2] = 1;

	::eps_dc = 72;
	::bPolarize = false; // polarizable water molecule

	::rcut_Ewd = 8;
	::rcut_LJ = 25;
	ConstructSimpleMolecule(&cw, "SPC");
	cw.c->cal_geometry_radius();
	dr.v[0] = -cw.r.v[0]; dr.v[1] = -cw.r.v[1]; dr.v[2] = -cw.r.v[2];
	cw.shiftMol(dr);
	::iFormat_LJ = 1;  // Amber
	init_LJ_Pars();
	construct_LJ12_6_prof(::iFormat_LJ);
	for (ia = 0; ia < cw.c->nAtoms; ia++) {
		cw.c->atom[ia].par->bLJ = false;
		cw.c->atom[ia].par->epsLJ = 0;
	}

	MMOL_MD_CELL mdcell;
	//int nmol = 25;
	int nmol = 55;
	mdcell.SetMol(0, nmol, 0);
	float h[3] = {10, 10, 10};
	mdcell.h[0] = h[0]; mdcell.h[1] = h[1]; mdcell.h[2] = h[2];  mdcell.period = true;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (h[2] - 1) * (ranf() > 0.5 ? 1 : -1) * ranf();
		mdcell.sm.m[imol].shiftMol(dr);
		//rand_vect3(axis.v);
		rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI * 0.5, true);

		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
	get_atoms(mdcell, mdcell.atom);
	get_bclusters(mdcell, mdcell.bcluster);
	mdcell.check_eneutral();

	mdcell.set_atom_multipole((bPolarize ? 1 : 0));
	mdcell.setMolIndx();
	molecule_check_periodic_cell(mdcell);
	check_atom_rg(mdcell, nThreads);
	mdcell.init_groups();

	{
		char fname[256] = "test_EwaldSum_atom_pos.dat";
		ofstream out;
		out.open(fname);
		for (imol = 0; imol < mdcell.sm.n; imol++) {
			out<<"SM "<<imol<<" : "<<endl;
			for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
				out<<ia<<" : "<<mdcell.sm.m[imol].c->atom[ia].rg.v[0]<<", "<<mdcell.sm.m[imol].c->atom[ia].rg.v[1]<<", "<<mdcell.sm.m[imol].c->atom[ia].rg.v[2]<<endl;
			}
		}
		out.close();
	}

	double f1 = 1, z1 = -h[2], f2 = 1, z2 = h[2];
	ARRAY<ImageCluster> cimg1, cimg2;
	ARRAY<MATOM*> aimg1, aimg2;
	CreateImageCluster(mdcell, z1, f1, cimg1, aimg1, mdcell.ccluster.n);
	CreateImageCluster(mdcell, z2, f2, cimg2, aimg2, mdcell.ccluster.n * 2);
	//UpdatePos_ImageCluster(cimg1, z1, cimg2, z2, nThreads); // position was updated in CreateImageCluster already

	CMM_CELL3D<BASIC_CLUSTER> cmm;
	::mlog.init("test.log", false);

	float h_cmm[3] = {h[0], h[1], h[2]* 2};
	cmm.set_cell(10, 10, 30);
	cmm.set_cell_range(true, -h_cmm[0], h_cmm[0], -h_cmm[1], h_cmm[1], -h_cmm[2], h_cmm[2]);

	mdcell.h[2] = h_cmm[2];
	
	cmm.set_cell_pos();
	int ncs = 10;
	CMM_array_set_storage<BASIC_CLUSTER>(&cmm, ncs);

	double Ep = 0;
	
	//cmm_init_subcell<BASIC_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size

	

	cmm_init_distribute_cluster(mdcell.sm.m, mdcell.sm.n, cmm, true); // add simple solid molecule into cmm
	cmm_init_distribute_cluster(cimg1.m, cimg1.n, cmm, false);
	cmm_init_distribute_cluster(cimg2.m, cimg2.n, cmm, false);

	int ncs0 = cmm.number_atoms();

	cmm_check_cluster<BASIC_CLUSTER>(cmm);
	int ncs1 = cmm.number_atoms();
	cmm_check<BASIC_CLUSTER>(cmm, nThreads);
	int ncs2 = cmm.number_atoms();
	_cmm_2d_::SMCheckRelshipInCell(mdcell.sm.m, mdcell.sm.n, cmm, nThreads);
	_cmm_2d_::ImageClusterCheckRelshipInCell(cimg1.m, cimg1.n, cmm, nThreads);
	_cmm_2d_::ImageClusterCheckRelshipInCell(cimg2.m, cimg2.n, cmm, nThreads);


	// calculat the total dipole of the whole cell
	cal_dipole(mdcell, 1, mdcell.mu_surf);
	memcpy(&(cmm.mu_surf), &(mdcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
	// for system with image charge, surface dipole effect was not considered properly at least
	// so disable surface dipole effect

	bool bExplicit = true; // calculate interaction explicitly, after SPME calculation ?

	// SPME parameters
	extern BSplineFunc<2> bsp;
	init_BSplineFunc(bsp, 8);

	::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 128;
	//::dim_spme[0] = 64; ::dim_spme[1] = 64; ::dim_spme[2] = 512;

	VSPME_ImageSlab<_VQatom> vspme_q; 
	vspme_q.set_MP(0); // charge only
	vspme_q.vspme.bST = true; vspme_q.vspme.bST_diagonal = true;
	vspme_q.Ivspme1.bST = true; vspme_q.Ivspme1.bST_diagonal = true;
	vspme_q.Ivspme2.bST = true; vspme_q.Ivspme2.bST_diagonal = true;
	vspme_q.img1_vspme.bST = true; vspme_q.img1_vspme.bST_diagonal = true;
	vspme_q.img2_vspme.bST = true; vspme_q.img2_vspme.bST_diagonal = true;

	vspme_q.f1 = f1; vspme_q.f2 = f2;
	vspme_q.z1 = z1; vspme_q.z2 = z2;
	

	VSPME_ImageSlab<_VMUatom> vspme_mu; 
	vspme_mu.set_MP(1); // with dipole
	vspme_mu.vspme.bST = true; vspme_mu.vspme.bST_diagonal = true;
	vspme_mu.Ivspme1.bST = true; vspme_mu.Ivspme1.bST_diagonal = true;
	vspme_mu.Ivspme2.bST = true; vspme_mu.Ivspme2.bST_diagonal = true;

	vspme_mu.f1 = f1; vspme_mu.f2 = f2;
	vspme_mu.z1 = z1; vspme_mu.z2 = z2;

	if (bPolarize) {
		mdcell.init_polarization_buff();
		mdcell.init_dipole_hist();
		mdcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		//vspme_mu.bsp = &(bsp);
		vspme_mu.set_BSplineFunc(&bsp);
		init_spme(mdcell.catom, aimg1, aimg2, vspme_mu); // init the viral atoms
		//vspme_mu.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);

		//vspme_mu.set_unit_cells(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2]);
		vspme_mu.set_cells(-h[0], -h[1], -h[2], h[0] * 2, h[1] * 2, h[2] * 2, dim_spme[0], dim_spme[1], dim_spme[2]);

		//vspme_mu.xl[0] = cmm.xl[0]; vspme_mu.xl[1] = cmm.xl[1]; vspme_mu.xl[2] = cmm.xl[2];
		//vspme_mu.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();
	}
	else { // charge only
		mdcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		//vspme_q.bsp = &(bsp);
		vspme_q.set_BSplineFunc(&bsp);
		init_spme(mdcell.catom, aimg1, aimg2, vspme_q); // init the viral atoms
		//vspme_q.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);
		
		//vspme_q.set_unit_cells(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2]);
		vspme_q.set_cells(-h[0], -h[1], -h[2], h[0] * 2, h[1] * 2, h[2] * 2, dim_spme[0], dim_spme[1], dim_spme[2]);
		vspme_q.set_cells(-h[0], -h[1], -h[2] * 1.5, h[0] * 2, h[1] * 2, h[2] * 3, dim_spme[0], dim_spme[1], dim_spme[2]);

		//vspme_q.xl[0] = cmm.xl[0]; vspme_q.xl[1] = cmm.xl[1]; vspme_q.xl[2] = cmm.xl[2];
		//vspme_q.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();
	}

	_spme_2d_::SPME_VAR spme_var;
	spme_var.bVirial = true; spme_var.bST_diagonal = true;

	spme_var.spme_3d_var.esv1.bEwald = true; 
	//spme_var.spme_3d_var.esv1.init_EwaldSum(cmm.xd[0] * cmm.xd[1] * cmm.xd[2], ::rcut_Ewd);
	if (bPolarize) {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_mu.vspme.V, ::rcut_Ewd);
	}
	else {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_q.vspme.V, ::rcut_Ewd);
	}
	spme_var.spme_3d_var.bVirial = spme_var.bVirial; spme_var.spme_3d_var.bST_diagonal = spme_var.bST_diagonal;
	spme_var.spme_3d_var.bTholePolarize = false;//::_Thole_polarize_::bTholePolarize;
	spme_var.spme_3d_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	spme_var.spme_3d_var.bEfieldOnly = false;

	double kappa = spme_var.spme_3d_var.esv1.kappa;
	spme_var.K0var.bEfieldOnly = false;
	spme_var.K0var.bVirial = spme_var.bVirial;
	spme_var.K0var.set(kappa, cmm.xd[0] * cmm.xd[1]);

	//V3zero(cmm.mu_surf.mu) cmm.mu_surf.q = 0; V3zero(cmm.mu_surf.qr) V3zero(cmm.mu_surf.mu0)
	InteractRes ires;
	if (bPolarize) {
		if (!Polarize_SPME_Interact(&mdcell, cimg1, f1, cimg2, f2, &cmm, vspme_mu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else SPME_Interact(&mdcell, cimg1, cimg2, &cmm, vspme_q, nThreads, spme_var, ires);
	Ep = ires.U;

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);
	
	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-SPME-EwaldSum.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);


/*
	{
	if (!bPolarize) {
		Ep = EwaldSum_Recip_2d_Q(mdcell);
		sprintf(errmsg, "Restrict Ewald Sum [reciprocal] : %f kT", Ep); show_log(errmsg, true); 

		init_normalized_coordinates(vspme_q, nThreads);
		cal_Q_spme(vspme_q, nThreads); // SPME Q matrix
		double U = vspme_q.cal_U_recp(); //calculate the energy from the image part 	

		sprintf(errmsg, "SPME Ewald Sum [reciprocal] : %f kT", U); show_log(errmsg, true);
	}
	}
*/

	// explicit calculation
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		V6zero(mdcell.sm.m[imol].c->dyn->fc)
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			V3zero(mdcell.sm.m[imol].c->atom[ia].E)
			V3zero(mdcell.sm.m[imol].c->atom[ia].F)
		}
	}

	E2dInteractVar evar_2d;
	evar_2d.bVirial = true; evar_2d.bST_diagonal = true; 
	evar_2d.bEfieldOnly = false; 
	evar_2d.esv1.bEwald = false; // direct interaction, not EwaldSum
	InteractRes tres;

	ires.reset();
	if (bExplicit) {
		_spme_2d_::ExplicitEInteract_ImageSlab(cmm, mdcell, aimg1, aimg2, evar_2d, ires);
		Ep = ires.U;

		if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
		if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
		if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);

		sprintf(errmsg, "explicit calculation : total energy %f [kT]", Ep); show_log(errmsg, true);
		{
		ofstream out;
		out.open("test-explicit.dat");
		out<<"cluster force : "<<endl;
		for (imol = 0; imol < mdcell.sm.n; imol++) {
			out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
			out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
		}

		out<<endl<<"dipole : "<<endl;
		for (imol = 0; imol < mdcell.sm.n; imol++) {
			out<<"SM : "<<imol<<endl;
			for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
				out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
			}
		}

		out<<endl<<"electric field : "<<endl;
		for (imol = 0; imol < mdcell.sm.n; imol++) {
			out<<"SM : "<<imol<<endl;
			for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
				out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
			}
		}

		out<<endl<<"atomic force : "<<endl;
		for (imol = 0; imol < mdcell.sm.n; imol++) {
			out<<"SM : "<<imol<<endl;
			for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
				out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
			}
		}

		out<<endl<<"virial : "<<endl;
		out<<"Total virial : "<<ires.STtr<<endl;
		out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

		out.close();
		}
	}
	
	::mlog.close();
}


#endif

