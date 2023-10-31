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

using namespace _evatom_;
using namespace _EwaldSum_real_;
using namespace _EwaldSum_2d_;

#include "atom_sr.h"
#include "atom_sr_interact.h"


#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
//long t_real = 0, t_recp = 0, t_emp = 0, t_cluster = 0, t_LJ = 0;
#endif

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
extern void cal_dipole(MMOL_VMD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);

extern PAIRWISE_DB< PAIRWISE<EF_DPROF> > sr_db; // pair-wise short-range interaction database

#ifdef __GPU__
namespace _gpu_md_ {
	extern __host__ void GPU_EwaldSumReal_2d(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev, cudaStream_t &stream);
}
#endif

namespace _atomic_cmm_ {
	#define LJ_FORCE(eps, rLJ, dis, u, i, t1, t2, V, f) \
		t1 = dis / rLJ; interpolate2(V, f[0], LJ12_6, t1, i, t2) \
		V *= eps; t2 = f[0] * eps / rLJ; \
		/*if (t2 > SR_FORCE_MAX) t2 = SR_FORCE_MAX; 	else if (t2 > force_max && scale_force) t2 = force_max; */\
		f[0] = t2 * u[0]; f[1] = t2 * u[1]; f[2] = t2 * u[2];  


void GPU_LJ_INTERACT_ATOM(JobIndx *job, SR_AtCELL *cell) {
	int ia1, ia, ia2, cIndx1, cIndx2, i3a1, i3a2;
	int icx, icy, icz, ncx, ncy, ncz, dcx, dcy, dcz;
	int Ncx, Ncy, Nc;
	int i;
	int *abound;
	bool bBound;

	double R, R2;
	double t1, t2, t3;
	double xd[3] = {cell->cw[0] * cell->Nx, cell->cw[1] * cell->Ny, cell->cw[2] * cell->Nz};
	double doff[3];

	double rcut = cell->rcut, r2cut = rcut * rcut;
	dcx = int(rcut / cell->cw[0]) + 2;
	dcy = int(rcut / cell->cw[1]) + 2;
	dcz = int(rcut / cell->cw[2]) + 2;

	double cw_max = cell->cw[0];
	if (cell->cw[1] > cw_max) cw_max = cell->cw[1];
	if (cell->cw[2] > cw_max) cw_max = cell->cw[2];
	double rcut_m = rcut + cw_max * 2, r2cut_m = rcut_m * rcut_m;
	double d2x, d2y, d2z, d2xy, d2;

	double r1min = 0, r2min = 0;
	double U_LJ = 0;
	double force[3], dr[3], u[3];
	LJ_PARS *pLJ = NULL;

	PAIRWISE<EF_DPROF> *pw = NULL;
	unsigned int key;
	char aindx[2];

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
	time.start();
#endif

	//V6zero(pc->dyn->fc)

	for (ia1 = job->n1; ia1 <= job->n2; ia1++) {
		//cell->Fx_sr.m[ia1] = 0; cell->Fy_sr.m[ia1] = 0; cell->Fz_sr.m[ia1] = 0; cell->U_LJ.m[ia1] = 0;
		// should be set to 0 before this function is called
		if (cell->aID.m[ia1] < 0) continue;
		cIndx1 = cell->cIndx.m[ia1];
		r1min = cell->sr_par.m[ia1 * cell->n_srpars + 2];

		i3a1 = x_idx(ia1);

		U_LJ = 0;
		for (icx = -dcx; icx <= dcx; icx++) {
			ncx = cell->ir.m[i3a1] + icx; 
			if (ncx < 0) {ncx += cell->Nx; doff[0] = -xd[0];}
			else if (ncx >= cell->Nx) {ncx -= cell->Nx; doff[0] = xd[0];}
			else doff[0] = 0;
			Ncx = ncx * cell->Nyz;
			
			d2x = icx * cell->cw[0]; d2x *= d2x;

			for (icy = -dcy; icy <= dcy; icy++) {
				ncy = cell->ir.m[i3a1 + 1] + icy; 
				if (ncy < 0) {ncy += cell->Ny; doff[1] = -xd[1];}
				else if (ncy >= cell->Ny) {ncy -= cell->Ny; doff[1] = xd[1];}
				else doff[1] = 0;
				Ncy = ncy * cell->Nz;

				d2y = icy * cell->cw[1]; d2y *= d2y;
				d2xy = d2x + d2y;
				if (d2xy > r2cut_m) continue;

				for (icz = -dcz; icz <= dcz; icz++) {
					ncz = cell->ir.m[i3a1 + 2] + icz; 
					if (ncz < 0) {ncz += cell->Nz; doff[2] = -xd[2];}
					else if (ncz >= cell->Nz) {ncz -= cell->Nz; doff[2] = xd[2];}
					else doff[2] = 0;

					d2z = icz * cell->cw[2]; d2z *= d2z;
					d2 = d2xy + d2z;
					if (d2 > r2cut_m) continue;

					Nc = (Ncx + Ncy + ncz) * (cell->NC + 1);
					for (ia = 0; ia < cell->cmm.m[Nc]; ia++) {
						ia2 = cell->cmm.m[Nc + ia + 1];
						cIndx2 = cell->cIndx.m[ia2];
						if (cIndx1 == cIndx2) continue; // in the same cluster
						if (cell->aID.m[ia2] < 0) continue;
						i3a2 = x_idx(ia2);
						r2min = cell->sr_par.m[ia2 * cell->n_srpars + 2];
						dr[0] = cell->r.m[i3a1] - cell->r.m[i3a2] - doff[0];
						dr[1] = cell->r.m[i3a1 + 1] - cell->r.m[i3a2 + 1] - doff[1];
						dr[2] = cell->r.m[i3a1 + 2] - cell->r.m[i3a2 + 2] - doff[2];
						R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
						if (R2 > r2cut) continue;
						R = sqrt(R2);

						// check bound and bShocked -- 2014/06/07
						bBound = false;
						abound = cell->abound.m + ia * (cell->Nm_b + 1);
						for (i = 0; i < abound[0]; i++) {
							if (abound[i+1] == ia2) {// 1-2 or 1-3 bound 
								bBound = true;
							}
						}
						// ignore bound interaction for short-range interaction
						if (bBound) continue;
						if (R < (r1min + r2min)) cell->bShocked.m[ia1] = 0x01; // not bound atom
						//******** 2014/06/07 *********

						u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

						// calculation of the interaction
						memset(force, 0, SIZE_V3);
						if ((cell->sr_type.m[ia1] & 0x01) && (cell->sr_type.m[ia2] & 0x01)) { // L-J interaction
							pLJ = &(LJPARS.m[cell->aID.m[ia1]].m[cell->aID.m[ia2]]);
							if (pLJ != NULL && pLJ->epsLJ > 0.0001 && R2 < pLJ->rcut2) {
								//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
								LJ_FORCE(pLJ->epsLJ, pLJ->rLJ, R, u, i, t1, t3, t2, force)
								//if (scale_force && t3 > force_max) {
								//	t3 = force_max / t3;
								//	force[0] *= t3; force[1] *= t3; force[[2] *= t3;
								//}
								/*
								if (bound14 && pc != pc1) { //1-4 LJ interaction / 2
									t2 *= A14_LJ; t3 *= A14_LJ;
									for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
								}
								*/
								U_LJ += t2 * 0.5;
							}
						}
						if ((cell->sr_type.m[ia1] & 0x02) && (cell->sr_type.m[ia2] & 0x02)) { // given profile -- bSR
							aindx[0] = cell->aID.m[ia1]; aindx[1] = cell->aID.m[ia2];
							key = construct_char2_key_pairwise(aindx);
							pw = sr_db.search(key);
							if (pw != NULL) {
								interpolate2(t1, t2, pw->pv, R, i, t3)
								U_LJ += t1 * 0.5;
								if (t2 > force_max && scale_force) t2 = force_max; 
								force[0] += t2 * u[0]; force[1] += t2 * u[1]; force[2] += t2 * u[2];
							}
						}
						//if (bVirial) {
							//if (bST_diagonal) Accumulate_Virial_Diagonal(res, force, dr);
							//else Accumulate_Virial_Trace(res, t3, R);
						//}
						// virial
						cell->vsr.m[i3a1] += force[0] * dr[0] * 0.5;
						cell->vsr.m[i3a1 + 1] += force[1] * dr[1] * 0.5;
						cell->vsr.m[i3a1 + 2] += force[2] * dr[2] * 0.5;

						cell->Fsr.m[i3a1] += force[0]; // force has mass unit, was converted in LJ12_6 profile
						cell->Fsr.m[i3a1 + 1] += force[1];
						cell->Fsr.m[i3a1 + 2] += force[2];
					}

				}
			}

			// 1-4 interactions needs to be calculated
			//for (ia = 0; ia < cell->r
		}

		cell->U_sr.m[ia1] = U_LJ;

#if _CHECK_TIME_STAMP_ == 1
	//t_LJ += time.elapse();
#endif

	}
}

// ****** 2014/06/07 ******
// short-range interaction for 1-2 and 1-3 bound are NOT calculated, so we do not need this function anymore
// and we change the function GPU_LJ_INTERACT_BOUND to GPU_LJ_INTERACT_BOUND14, since it was not calculated before !
void GPU_LJ_INTERACT_BOUND14(JobIndx *job, SR_AtCELL *cell) {
	int ib, ia1, ia2, i3a1, i3a2, i3b;//, cIndx1, cIndx2;
	int ix1, iy1, iz1, ix2, iy2, iz2, dx, dy, dz;
	//int Nx, Ny, Nz;
	int Nm;
	int i;

	int hx = cell->Nx / 2, hy = cell->Ny, hz = cell->Nz / 2;

	double R, R2;
	double t1, t2, t3;
	double xd[3] = {cell->cw[0] * cell->Nx, cell->cw[1] * cell->Ny, cell->cw[2] * cell->Nz};
	double doff[3];

	double U_LJ = 0;
	double force[3], dr[3], u[3];
	LJ_PARS *pLJ = NULL;

	PAIRWISE<EF_DPROF> *pw = NULL;
	unsigned int key;
	char aindx[2];

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
	time.start();
#endif

	//V6zero(pc->dyn->fc)

	for (ib = job->n1; ib <= job->n2; ib++) {
		if (ib >= cell->bound_14.n) break;
		Nm = ib * 2;
		ia1 = cell->bound_14.m[Nm]; ia2 = cell->bound_14.m[Nm + 1];
		i3a1 = x_idx(ia1); i3a2 = x_idx(ia2);
		if (cell->cIndx.m[ia1] == cell->cIndx.m[ia2]) continue; // they are in the same cluster, which was ignored in short-range interaction calculation
		
		ix1 = cell->ir.m[i3a1]; iy1 = cell->ir.m[i3a1 + 1]; iz1 = cell->ir.m[i3a1 + 2];
		ix2 = cell->ir.m[i3a2]; iy2 = cell->ir.m[i3a2 + 1]; iz2 = cell->ir.m[i3a2 + 2];
		dx = ix2 - ix1; dy = iy2 - iy1; dz = iz2 - iz1;
		if (dx > hx) doff[0] = -xd[0]; else if (dx < -hx) doff[0] = xd[0]; else doff[0] = 0;
		if (dy > hy) doff[1] = -xd[1]; else if (dy < -hy) doff[1] = xd[1]; else doff[1] = 0;
		if (dz > hz) doff[2] = -xd[2]; else if (dz < -hz) doff[2] = xd[2]; else doff[2] = 0;

		dr[0] = cell->r.m[i3a2] - cell->r.m[i3a1] + doff[0];
		dr[1] = cell->r.m[i3a2 + 1] - cell->r.m[i3a1 + 1] + doff[1];
		dr[2] = cell->r.m[i3a2 + 2] - cell->r.m[i3a1 + 2] + doff[2];
		R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
		R = sqrt(R2);
		//if (R > cell->sr_cut) continue;
		u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

		U_LJ = 0;

		// calculation of the interaction
		if ((cell->sr_type.m[ia1] & 0x01) && (cell->sr_type.m[ia1] & 0x01)) { // L-J interaction
			pLJ = &(LJPARS.m[cell->aID.m[ia1]].m[cell->aID.m[ia1]]);
			if (pLJ != NULL && pLJ->epsLJ > 0.0001 && R2 < pLJ->rcut2) {
				//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
				LJ_FORCE(pLJ->epsLJ, pLJ->rLJ, R, u, i, t1, t3, t2, force)
				//if (scale_force && t3 > force_max) {
				//	t3 = force_max / t3;
				//	force[0] *= t3; force[1] *= t3; force[[2] *= t3;
				//}
				/*
				if (bound14 && pc != pc1) { //1-4 LJ interaction / 2
					t2 *= A14_LJ; t3 *= A14_LJ;
					for (i = 0; i < 3; i++) force.v[i] *= A14_LJ;
				}
				*/
				U_LJ += t2;
			}
		}
		if ((cell->sr_type.m[ia1] & 0x02) && (cell->sr_type.m[ia2] & 0x02)) { // given profile -- bSR
			aindx[0] = cell->aID.m[ia1]; aindx[1] = cell->aID.m[ia2];
			key = construct_char2_key_pairwise(aindx);
			pw = sr_db.search(key);
			if (pw != NULL) {
				interpolate2(t1, t2, pw->pv, R, i, t3)
				U_LJ += t1;
				if (t2 > force_max && scale_force) t2 = force_max; 
				force[0] += t2 * u[0]; force[1] += t2 * u[1]; force[2] += t2 * u[2];
			}
		}
		i3b = x_idx(ib);
		cell->vsr_14.m[i3b] += force[0] * dr[0];
		cell->vsr_14.m[i3b + 1] += force[1] * dr[1];
		cell->vsr_14.m[i3b + 2] += force[2] * dr[2];

		cell->U14_sr.m[ib] = U_LJ;


#if _CHECK_TIME_STAMP_ == 1
	//t_LJ += time.elapse();
#endif

	}
}




// interaction on patom1 from patom2, 
// dr -- vector from patom2 to patom1
// R -- distance |dr|
// R -- |dr|^2
// E -- electric field from patom2
// F -- force from patom2 to patom1
// U -- interaction energy
void EwaldReal_Q(EwaldSumRealVars1 &esv, double c1, double c2, double* dr, double R, double R2, /*double *E, */double *F, double& U) {
	double t, c12 = c1 * c2;
	esv.EwaldSumRealVars0::init_vars(dr, R, R2);

	// U
	t = c12 / R; // q1 - q2
	U = t * esv.f0;

	// E
	/*
	E.v[0] = -c2 * esv.mET1.v[0]; // from q2
	E.v[1] = -c2 * esv.mET1.v[1]; // from q2
	E.v[2] = -c2 * esv.mET1.v[2]; // from q2
	*/

	// F
	F[0] = -c12 * esv.mET1.v[0]; // q1 - q2
	F[1] = -c12 * esv.mET1.v[1]; // q1 - q2
	F[2] = -c12 * esv.mET1.v[2]; // q1 - q2
}

void EInteract_Q(EwaldSumRealVars1 &esv, double c1, double c2, double* dr, double R, double R2, /*double *E, */double *F, double& U, bool init_esv) {
	double t, c12 = c1 * c2;

	if (init_esv) esv.EwaldSumRealVars0::init_vars(dr, R, R2);

	// U
	U = c12 / R; // q1 - q2
	//U = t * esv.f0;

	// E
	/*
	E.v[0] = -c2 * esv.mET1.v[0]; // from q2
	E.v[1] = -c2 * esv.mET1.v[1]; // from q2
	E.v[2] = -c2 * esv.mET1.v[2]; // from q2
	*/

	// F
	F[0] = -c12 * esv.mT1.v[0]; // q1 - q2
	F[1] = -c12 * esv.mT1.v[1]; // q1 - q2
	F[2] = -c12 * esv.mT1.v[2]; // q1 - q2
}

void EwaldReal_MU(EwaldSumRealVars1 &esv, double c1, double *mu1, double c2, double *mu2, double* dr, double R, double R2, /*double *E,*/ double *F, double& U, char bExcludeDipole) {
	int i, j, k;
	double t;

	if (bExcludeDipole) esv.init_vars_E(dr, R, R2); // no induced dipole and dipole-dipole interaction will be excluded
	else esv.init_vars(dr, R, R2);

	// U
	t = c1 * c2 / R; // q1 - q2
	U = t * esv.f0;

	U += c2 * (mu1[0] * esv.mET1.v[0] + mu1[1] * esv.mET1.v[1] + mu1[2] * esv.mET1.v[2]);
	U -= c1 * (mu2[0] * esv.mET1.v[0] + mu2[1] * esv.mET1.v[1] + mu2[2] * esv.mET1.v[2]);

	if (bExcludeDipole == 0) {
		U -= mu1[0] * ( esv.mET2.m[0][0] * mu2[0] + esv.mET2.m[0][1] * mu2[1] + esv.mET2.m[0][2] * mu2[2])
			+ mu1[1] * ( esv.mET2.m[1][0] * mu2[0] + esv.mET2.m[1][1] * mu2[1] + esv.mET2.m[1][2] * mu2[2])
			+ mu1[2] * ( esv.mET2.m[2][0] * mu2[0] + esv.mET2.m[2][1] * mu2[1] + esv.mET2.m[2][2] * mu2[2]);
	}

	// E
	/*
	for (i = 0; i < 3; i++) {
		E.v[i] -= c2 * esv.mET1.v[i]; // from q2
		for (j = 0; j < 3; j++) E.v[i] += esv.mET2.m[i][j] * mu2[j];
	}
	*/

	// F
	for (i = 0; i < 3; i++) {
		F[i] = -c1 * c2 * esv.mET1.v[i]; // q1 - q2
		for (j = 0; j < 3; j++) {
			F[i] += c1 * esv.mET2.m[i][j] * mu2[j];
			F[i] -= c2 * esv.mET2.m[i][j] * mu1[j];
		}
		if (!bExcludeDipole) { // mu1 - mu2
			for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
				F[i] += mu1[j] * esv.mET3.m[i][j][k] * mu2[k];
			}}
		}
	}
}

void EInteract_MU(EwaldSumRealVars1 &esv, double c1, double *mu1, double c2, double *mu2, double* dr, double R, double R2, /*double *E,*/ double *F, double& U, char bExcludeDipole, bool init_esv) {
	int i, j, k;
	double t;

	if (init_esv) {
		if (bExcludeDipole) esv.init_vars_E(dr, R, R2); // no induced dipole and dipole-dipole interaction will be excluded
		else esv.init_vars(dr, R, R2);
	}

	// U
	U = c1 * c2 / R; // q1 - q2
	//U = t * esv.f0;

	U += c2 * (mu1[0] * esv.mT1.v[0] + mu1[1] * esv.mT1.v[1] + mu1[2] * esv.mT1.v[2]);
	U -= c1 * (mu2[0] * esv.mT1.v[0] + mu2[1] * esv.mT1.v[1] + mu2[2] * esv.mT1.v[2]);

	if (bExcludeDipole == 0) {
		U -= mu1[0] * ( esv.mT2.m[0][0] * mu2[0] + esv.mT2.m[0][1] * mu2[1] + esv.mT2.m[0][2] * mu2[2])
			+ mu1[1] * ( esv.mT2.m[1][0] * mu2[0] + esv.mT2.m[1][1] * mu2[1] + esv.mT2.m[1][2] * mu2[2])
			+ mu1[2] * ( esv.mT2.m[2][0] * mu2[0] + esv.mT2.m[2][1] * mu2[1] + esv.mT2.m[2][2] * mu2[2]);
	}

	// E
	/*
	for (i = 0; i < 3; i++) {
		E.v[i] -= c2 * esv.mET1.v[i]; // from q2
		for (j = 0; j < 3; j++) E.v[i] += esv.mET2.m[i][j] * mu2[j];
	}
	*/

	// F
	for (i = 0; i < 3; i++) {
		F[i] = -c1 * c2 * esv.mT1.v[i]; // q1 - q2
		for (j = 0; j < 3; j++) {
			F[i] += c1 * esv.mT2.m[i][j] * mu2[j];
			F[i] -= c2 * esv.mT2.m[i][j] * mu1[j];
		}
		if (!bExcludeDipole) { // mu1 - mu2
			for (j = 0; j < 3; j++) {for (k = 0; k < 3; k++) {
				F[i] += mu1[j] * esv.mT3.m[i][j][k] * mu2[k];
			}}
		}
	}
}


// calculate the electric field only about the real part of Ewald Sum
void EwaldReal_E_Q(EwaldSumRealVars1 &esv, double c2, double* dr, double R, double R2, double *E) {
	esv.EwaldSumRealVars0::init_vars(dr, R, R2);

	// E
	E[0] = -c2 * esv.mET1.v[0]; // from q2
	E[1] = -c2 * esv.mET1.v[1];
	E[2] = -c2 * esv.mET1.v[2];
}

void EInteract_E_Q(EwaldSumRealVars1 &esv, double c2, double* dr, double R, double R2, double *E, bool init_esv) {
	if (init_esv) esv.EwaldSumRealVars0::init_vars(dr, R, R2);

	// E
	E[0] = -c2 * esv.mT1.v[0]; // from q2
	E[1] = -c2 * esv.mT1.v[1];
	E[2] = -c2 * esv.mT1.v[2];
}

void EwaldReal_E_MU(EwaldSumRealVars1 &esv, double c2, double *mu2, double* dr, double R, double R2, double *E) {
	int i, j;
	esv.init_vars_E(dr, R, R2);

	// E
	for (i = 0; i < 3; i++) {
		E[i] = -c2 * esv.mET1.v[i]; // from q2
		for (j = 0; j < 3; j++) E[i] += esv.mET2.m[i][j] * mu2[j];
	}
}

void EInteract_E_MU(EwaldSumRealVars1 &esv, double c2, double *mu2, double* dr, double R, double R2, double *E, bool init_esv) {
	int i, j;
	if (init_esv) esv.init_vars_E(dr, R, R2);

	// E
	for (i = 0; i < 3; i++) {
		E[i] = -c2 * esv.mT1.v[i]; // from q2
		for (j = 0; j < 3; j++) E[i] += esv.mT2.m[i][j] * mu2[j];
	}
}

/*
void MTP_E_Q(EwaldSumRealVars0 &esv, double c2, double *dr, double R, double R2, double *E) {
	esv.init_vars(VECTOR3(dr), R, R2);
	E[0] = -c2 * esv.mT1.v[0];
	E[1] = -c2 * esv.mT1.v[1];
	E[2] = -c2 * esv.mT1.v[2];
}

void MTP_E_MU(EwaldSumRealVars1 &esv, double c2, double *mu2, double *dr, double R, double R2, double *E) {
	esv.init_vars(VECTOR3(dr), R, R2);
	int i, j;
	// E, - q2, mu2
	for (i = 0; i < 3; i++) {
		// E from q2
		E[i] -= c2 * esv.mT1.v[i];
		for (j = 0; j < 3; j++) {
			// E from mu2
			E[i] += esv.mT2.m[i][j] * mu2[j];
		}
	}
}

void MTP_Interact_Q(EwaldSumRealVars0 &esv, double c1, double c2, double *dr, double R, double R2, double *F, double& U) {
	esv.init_vars(VECTOR3(dr), R, R2);
	double c12 = c1 * c2;
	U = c12 / esv.R; // q1 - q2
	F[0] = -c12 * esv.mT1.v[0];
	F[1] = -c12 * esv.mT1.v[1];
	F[2] = -c12 * esv.mT1.v[2];
}


void MTP_Interact_MU(EwaldSumRealVars1 &esv, double c1, double *mu1, double c2, double *mu2, double *dr, double R, double R2, double *F, double& U, char bExcludeDipole) {
	esv.init_vars(VECTOR3(dr), R, R2);
	int i, j, k;
	double c12 = c1 * c2;
	U = c12 / esv.R; // q1 - q2

	U += c2 * (mu1[0] * esv.mT1.v[0] + mu1[1] * esv.mT1.v[1] + mu1[2] * esv.mT1.v[2]); // mu1 - q2
	U -= c1 * (mu2[0] * esv.mT1.v[0] + mu2[1] * esv.mT1.v[1] + mu2[2] * esv.mT1.v[2]); // q1 - mu2
	if (bExcludeDipole == 1) {
		// mu1 - mu2
		U -=  mu1[0] * (esv.mT2.m[0][0] * mu2[0] + esv.mT2.m[0][1] * mu2[1] + esv.mT2.m[0][2] * mu2[2])
			+ mu1[1] * (esv.mT2.m[1][0] * mu2[0] + esv.mT2.m[1][1] * mu2[1] + esv.mT2.m[1][2] * mu2[2])
			+ mu1[2] * (esv.mT2.m[2][0] * mu2[0] + esv.mT2.m[2][1] * mu2[1] + esv.mT2.m[2][2] * mu2[2]);
	}

	// F, q1, mu1 - q2, mu2
	for (i = 0; i < 3; i++) {
		// force
		F[i] -= c12 * esv.mT1.v[i]; // q1 - q2
		for (j = 0; j < 3; j++) {
			// force
			F[i] += c1 * esv.mT2.m[i][j] * mu2[j]; // q1 - mu2
			F[i] -= c2 * esv.mT2.m[i][j] * mu1[j]; // q2 - mu1
			if (bExcludeDipole == 1) {
				for (k = 0; k < 3; k++) { // mu1 - mu2
					F[i] += mu1[j] * esv.mT3.m[i][j][k] * mu2[k];
				}
			}
		}
	}
}
*/


// real part of Ewald-Sum
void GPU_EwdSumReal_INTERACT_ATOM(JobIndx *job, EAtCELL *cell) {
	int ia1, ia, ia2, cIndx1, cIndx2, i3a1, i3a2;
	int icx, icy, icz, ncx, ncy, ncz, dcx, dcy, dcz;
	int Ncx, Ncy, Nc;

	double R, R2;
	double t1, t2, t3;
	double xd[3] = {cell->cw[0] * cell->Nx, cell->cw[1] * cell->Ny, cell->cw[2] * cell->Nz};
	double doff[3];

	double rcut = cell->rcut, r2cut = rcut * rcut;
	dcx = int(rcut / cell->cw[0]) + 2;
	dcy = int(rcut / cell->cw[1]) + 2;
	dcz = int(rcut / cell->cw[2]) + 2;

	double cw_max = cell->cw[0];
	if (cell->cw[1] > cw_max) cw_max = cell->cw[1];
	if (cell->cw[2] > cw_max) cw_max = cell->cw[2];
	double rcut_m = rcut + cw_max * 2, r2cut_m = rcut_m * rcut_m;
	double d2x, d2y, d2z, d2xy, d2;

	double U_estat = 0;
	double dr[3];//, u[3];
	double c1, c2, mu1[3], mu2[3], E[3], Et[3], F[3], Ft[3];

	double U, Ut, Uestat = 0;

	EwaldSumRealVars1 esv;
	esv.cp_EwaldSum_vars(*(EwaldSumRealVars0*)(&cell->esv));

	short iSurfaceBoundary = esv.iSurfaceBoundary;

	short bExcludeDipole = esv.iExIndDipole;
	int cal_mode = cell->bCalE;
	if (cal_mode <= 0 || cal_mode >= 3) return; // with charge only, so far choice 1 -- electric field; 2 -- force and energy
	if (cell->bDipole) cal_mode += 10; // with dipole

	//ATOM_COMM_PARS *apar1 = NULL, *apar2 = NULL;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

	//V6zero(pc->dyn->fc)

	for (ia1 = job->n1; ia1 <= job->n2; ia1++) {
		//cell->Fx_sr.m[ia1] = 0; cell->Fy_sr.m[ia1] = 0; cell->Fz_sr.m[ia1] = 0; cell->U_LJ.m[ia1] = 0;
		// should be set to 0 before this function is called
		//if (cell->aID.m[ia1] < 0) continue;
		//if (::atompar_db.par.m[cell->aID.m[ia1]]->epsLJ < 1e-4) continue;
		i3a1 = x_idx(ia1);
		c1 = cell->c.m[ia1]; mu1[0] = cell->mu.m[i3a1]; mu1[1] = cell->mu.m[i3a1 + 1]; mu1[2] = cell->mu.m[i3a1 + 2];
		cIndx1 = cell->cIndx.m[ia1];

		Uestat = 0;
		for (icx = -dcx; icx <= dcx; icx++) {
			ncx = cell->ir.m[i3a1] + icx; 
			if (ncx < 0) {ncx += cell->Nx; doff[0] = -xd[0];}
			else if (ncx >= cell->Nx) {ncx -= cell->Nx; doff[0] = xd[0];}
			else doff[0] = 0;
			Ncx = ncx * cell->Nyz;
			
			d2x = icx * cell->cw[0]; d2x *= d2x;

			for (icy = -dcy; icy <= dcy; icy++) {
				ncy = cell->ir.m[i3a1 + 1] + icy; 
				if (ncy < 0) {ncy += cell->Ny; doff[1] = -xd[1];}
				else if (ncy >= cell->Ny) {ncy -= cell->Ny; doff[1] = xd[1];}
				else doff[1] = 0;
				Ncy = ncy * cell->Nz;

				d2y = icy * cell->cw[1]; d2y *= d2y;
				d2xy = d2x + d2y;
				if (d2xy > r2cut_m) continue;

				for (icz = -dcz; icz <= dcz; icz++) {
					ncz = cell->ir.m[i3a1 + 2] + icz; 
					if (ncz < 0) {ncz += cell->Nz; doff[2] = -xd[2];}
					else if (ncz >= cell->Nz) {ncz -= cell->Nz; doff[2] = xd[2];}
					else doff[2] = 0;

					d2z = icz * cell->cw[2]; d2z *= d2z;
					d2 = d2xy + d2z;
					if (d2 > r2cut_m) continue;

					Nc = (Ncx + Ncy + ncz) * (cell->NC + 1);
					for (ia = 0; ia < cell->cmm.m[Nc]; ia++) {
						ia2 = cell->cmm.m[Nc + ia + 1];
						if (ia2 == ia1) continue;
						i3a2 = x_idx(ia2);
						cIndx2 = cell->cIndx.m[ia2];
						//if (cIndx1 == cIndx2) continue; // in the same cluster
						//if (cell->aID.m[ia2] < 0) continue;
						c2 = cell->c.m[ia2]; mu2[0] = cell->mu.m[i3a2]; mu2[1] = cell->mu.m[i3a2 + 1]; mu2[2] = cell->mu.m[i3a2 + 2];
						dr[0] = cell->r.m[i3a1] - cell->r.m[i3a2] - doff[0];
						dr[1] = cell->r.m[i3a1 + 1] - cell->r.m[i3a2 + 1] - doff[1];
						dr[2] = cell->r.m[i3a1 + 2] - cell->r.m[i3a2 + 2] - doff[2];
						R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
						if (R2 > r2cut) continue;
						R = sqrt(R2);
						//u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

						// calculation of the interaction
						switch (cal_mode) {
						case 1: // electric field only
							EwaldReal_E_Q(esv, c2, dr, R, R2, E);
							/*
							if (cell->cIndx.m[ia1] == cell->cIndx.m[ia2]) { // excluding that from the atoms inside the cluster
								esv.bEwald = false;
								EInteract_E_Q(esv, c2, dr, R, R2, E, false);
								E[0] -= Et[0]; E[1] -= Et[1]; E[2] -= Et[2];
								esv.bEwald = true;
							}
							*/
							cell->E.m[i3a1] += E[0]; cell->E.m[i3a1 + 1] += E[1]; cell->E.m[i3a1 + 2] += E[2];
							break;
						case 2: // force and U
							EwaldReal_Q(esv, c1, c2, dr, R, R2, F, U);
							/*
							if (cell->cIndx.m[ia1] == cell->cIndx.m[ia2]) { // excluding that from the atoms inside the cluster
								esv.bEwald = false;
								EInteract_Q(esv, c1, c2, dr, R, R2, Ft, Ut, false);
								F[0] -= Ft[0]; F[1] -= Ft[1]; F[2] -= Ft[2];
								U -= Ut;
								esv.bEwald = true;
							}
							*/
							cell->Fe.m[i3a1] += F[0]; cell->Fe.m[i3a1 + 1] += F[1]; cell->Fe.m[i3a1 + 2] += F[2];
							Uestat += U * 0.5;
							// virial
							cell->ve.m[i3a1] += F[0] * dr[0] * 0.5;
							cell->ve.m[i3a1 + 1] += F[1] * dr[1] * 0.5;
							cell->ve.m[i3a1 + 2] += F[2] * dr[2] * 0.5;
							break;
						case 11: 
							EwaldReal_E_MU(esv, c2, mu2, dr, R, R2, E);
							/*
							if (cell->cIndx.m[ia1] == cell->cIndx.m[ia2]) { // excluding that from the atoms inside the cluster
								esv.bEwald = false;
								EInteract_E_MU(esv, c2, mu2, dr, R, R2, Et, false);
								E[0] -= Et[0]; E[1] -= Et[1]; E[2] -= Et[2];
								esv.bEwald = true;
							}
							*/
							cell->E.m[i3a1] += E[0]; cell->E.m[i3a1 + 1] += E[1]; cell->E.m[i3a1 + 2] += E[2];
							break;
						case 12: // force and U
							EwaldReal_MU(esv, c1, mu1, c2, mu2, dr, R, R2, F, U, (char)bExcludeDipole);
							/*
							if (cell->cIndx.m[ia1] == cell->cIndx.m[ia2]) { // excluding that from the atoms inside the cluster
								esv.bEwald = false;
								EInteract_MU(esv, c1, mu1, c2, mu2, dr, R, R2, Ft, Ut, (char)bExcludeDipole, false);
								F[0] -= Ft[0]; F[1] -= Ft[1]; F[2] -= Ft[2];
								U -= Ut;
								esv.bEwald = true;
							}
							*/
							cell->Fe.m[i3a1] += F[0]; cell->Fe.m[i3a1 + 1] += F[1]; cell->Fe.m[i3a1 + 2] += F[2];
							Uestat += U * 0.5;
							// virial
							cell->ve.m[i3a1] += F[0] * dr[0] * 0.5;
							cell->ve.m[i3a1 + 1] += F[1] * dr[1] * 0.5;
							cell->ve.m[i3a1 + 2] += F[2] * dr[2] * 0.5;
							break;
						default:
							break;
						}
					}

				}
			}
		}

		switch (cal_mode) {
		case 1: // q, electric field
			// no correction
			break;
		case 2:
			Uestat -= esv.kappa * INV_SQRT_PI * c1 * c1; // self energy from charge
			break;
		case 11: // electric field with dipole
			// self term of electric field from dipole
			t1 = 4 * esv.kappa3 / 3 * INV_SQRT_PI;
			cell->E.m[i3a1] += t1 * mu1[0];
			cell->E.m[i3a1 + 1] += t1 * mu1[1];
			cell->E.m[i3a1 + 2] += t1 * mu1[2];
			break;
		case 12:
			Uestat -= esv.kappa * INV_SQRT_PI * c1 * c1; // self energy from charge
			// self energy from dipole
			Uestat -= 2 * esv.kappa3 / 3 * INV_SQRT_PI * (mu1[0] * mu1[0] + mu1[1] * mu1[1] + mu1[2] * mu1[2]);
			break;
		default:
			break;
		}

		// the electric field from the surface induced dipole, virial of surface dipole will be calculated later
		if (esv.iSurface == 0) { // cubic
			switch (esv.iSurfaceBoundary) { // normal
			case 0:
				if (cal_mode == 1 || cal_mode == 11) { // electric field only
					E[0] = -esv.inv_3V * cell->mu_surf.mu.v[0];
					E[1] = -esv.inv_3V * cell->mu_surf.mu.v[1];
					E[2] = -esv.inv_3V * cell->mu_surf.mu.v[2];

					cell->E.m[i3a1] += E[0];
					cell->E.m[i3a1 + 1] += E[1];
					cell->E.m[i3a1 + 2] += E[2];
				}
				else if (cal_mode == 2 || cal_mode == 12) {// force only. energy and virial will be calculated somewhere else
					if (esv.iExIndDipole == 0) {
						F[0] = -c1 * esv.inv_3V * cell->mu_surf.mu.v[0];// * fUnit_estat_mass;
						F[1] = -c1 * esv.inv_3V * cell->mu_surf.mu.v[1];// * fUnit_estat_mass;
						F[2] = -c1 * esv.inv_3V * cell->mu_surf.mu.v[2];// * fUnit_estat_mass;
					}
					else if (esv.iExIndDipole == 1) {
						F[0] = -c1 * esv.inv_3V * (cell->mu_surf.mu.v[0] - cell->mu_surf.mu_ind.v[0]);// * fUnit_estat_mass;
						F[1] = -c1 * esv.inv_3V * (cell->mu_surf.mu.v[1] - cell->mu_surf.mu_ind.v[1]);// * fUnit_estat_mass;
						F[2] = -c1 * esv.inv_3V * (cell->mu_surf.mu.v[2] - cell->mu_surf.mu_ind.v[2]);// * fUnit_estat_mass;
					}

					cell->Fe.m[i3a1] += F[0];
					cell->Fe.m[i3a1 + 1] += F[1];
					cell->Fe.m[i3a1 + 2] += F[2];
				}
				break;
			case 1: // disable surface dipole in all direction
				break;
			case 2: // disable surface dipole along z axis
				if (cal_mode == 1 || cal_mode == 11) { // electric field only
					E[0] = -esv.inv_3V * cell->mu_surf.mu.v[0];
					E[1] = -esv.inv_3V * cell->mu_surf.mu.v[1];
					//E[2] = -esv->inv_3V * cell->mu_surf.mu.v[2];

					cell->E.m[i3a1] += E[0];
					cell->E.m[i3a1 + 1] += E[1];
					//cell->E.m[i3a1 + 2] += E[2];
				}
				else if (cal_mode == 2 || cal_mode == 12) {// force only. energy and virial will be calculated somewhere else
					if (esv.iExIndDipole == 0) {
						F[0] = -c1 * esv.inv_3V * cell->mu_surf.mu.v[0];// * fUnit_estat_mass;
						F[1] = -c1 * esv.inv_3V * cell->mu_surf.mu.v[1];// * fUnit_estat_mass;
						//F[2] = -c1 * esv->inv_3V * cell->mu_surf.mu.v[2];// * fUnit_estat_mass;
					}
					else if (esv.iExIndDipole == 1) {
						F[0] = -c1 * esv.inv_3V * (cell->mu_surf.mu.v[0] - cell->mu_surf.mu_ind.v[0]);// * fUnit_estat_mass;
						F[1] = -c1 * esv.inv_3V * (cell->mu_surf.mu.v[1] - cell->mu_surf.mu_ind.v[1]);// * fUnit_estat_mass;
						//F[2] = -c1 * esv.inv_3V * (cell->mu_surf.mu.v[2] - cell->mu_surf.mu_ind.v[2]);// * fUnit_estat_mass;
					}

					cell->Fe.m[i3a1] += F[0];
					cell->Fe.m[i3a1 + 1] += F[1];
					//cell->Fe.m[i3a1 + 2] += F[2];
				}
				break;
			default:
				break;
			}
		}
		//else if (esv.iSurface == 1) { // slab
			//switch (esv.iSurfaceBoundary) {
			//case 0: // normal
				// confirmed with explicit calculation, this boundary condition should be disabled
			//	break;
			//default: // do not do anything
			//	break;
			//}
		//}

		if (cal_mode == 2 || cal_mode == 12) cell->U_estat.m[ia1] += Uestat;

#if _CHECK_TIME_STAMP_ == 1
	//t_LJ += time.elapse();
#endif

	}
}


void GPU_EwdSumReal_INTERACT_InCluster(JobIndx *job, EAtCELL *cell) {
	int ic, ica1, ica2, ia1, ia2, i3a1, i3a2;//, cIndx1, cIndx2;
	int ix1, iy1, iz1, ix2, iy2, iz2, dx, dy, dz;
	//int Nx, Ny, Nz;
	int Nm;

	int hx = cell->Nx / 2, hy = cell->Ny, hz = cell->Nz / 2;

	double R, R2;
	double t1, t2, t3;
	double E[3], F[3], mu1[3], mu2[3];

	double c1, c2;
	double U = 0, Uestat = 0;
	double dr[3];

	char bExcludeDipole = (char)cell->esv.iExIndDipole;
	int cal_mode = cell->bCalE;
	if (cal_mode <= 0 || cal_mode >= 3) return; // with charge only, so far choice 1 -- electric field; 2 -- force and energy
	if (cell->bDipole) cal_mode += 10; // with dipole

	EwaldSumRealVars1 esv;
	esv.cp_EwaldSum_vars(*(EwaldSumRealVars0*)(&cell->esv));
	esv.bEwald = false;


#if _CHECK_TIME_STAMP_ == 1
	TIME time;
	time.start();
#endif

	//V6zero(pc->dyn->fc)

	for (ic = job->n1; ic <= job->n2; ic++) {
		Nm = ic * (cell->Nca + 1);
		if (cell->cluster.m[Nm] <= 1) continue; // has only one atom
		for (ica1 = 0; ica1 < cell->cluster.m[Nm]; ica1++) {
			ia1 = cell->cluster.m[Nm + ica1 + 1];
			i3a1 = x_idx(ia1);
			c1 = cell->c.m[ia1]; mu1[0] = cell->mu.m[i3a1]; mu1[1] = cell->mu.m[i3a1 + 1]; mu1[2] = cell->mu.m[i3a1 + 2];
			Uestat = 0;
			for (ica2 = 0; ica2 < cell->cluster.m[Nm]; ica2++) {
				if (ica2 == ica1) continue;
				ia2 = cell->cluster.m[Nm + ica2 + 1]; i3a2 = x_idx(ia2);
				c2 = cell->c.m[ia2]; mu2[0] = cell->mu.m[i3a2]; mu2[1] = cell->mu.m[i3a2 + 1]; mu2[2] = cell->mu.m[i3a2 + 2];
		
				dr[0] = cell->r.m[i3a1] - cell->r.m[i3a2];
				dr[1] = cell->r.m[i3a1 + 1] - cell->r.m[i3a2 + 1];
				dr[2] = cell->r.m[i3a1 + 2] - cell->r.m[i3a2 + 2];
				R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
				R = sqrt(R2);
				//if (R > cell->sr_cut) continue;
				//u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

				// calculate electrostatic interaction, excluding the interactions between atoms inside the same cluster
				switch (cal_mode) {
				case 1: // electric field only
					EInteract_E_Q(esv, c2, dr, R, R2, E, true);
					cell->E.m[i3a1] -= E[0]; cell->E.m[i3a1 + 1] -= E[1]; cell->E.m[i3a1 + 2] -= E[2]; // excluding
					break;
				case 2: // force and U
					EInteract_Q(esv, c1, c2, dr, R, R2, F, U, true);
					cell->Fe.m[i3a1] -= F[0]; cell->Fe.m[i3a1 + 1] -= F[1]; cell->Fe.m[i3a1 + 2] -= F[2]; // excluding
					Uestat -= U * 0.5;
					// virial
					cell->ve.m[i3a1] -= F[0] * dr[0] * 0.5;
					cell->ve.m[i3a1 + 1] -= F[1] * dr[1] * 0.5;
					cell->ve.m[i3a1 + 2] -= F[2] * dr[2] * 0.5;
					break;
				case 11: 
					EInteract_E_MU(esv, c2, mu2, dr, R, R2, E, true);
					cell->E.m[i3a1] -= E[0]; cell->E.m[i3a1 + 1] -= E[1]; cell->E.m[i3a1 + 2] -= E[2];
					break;
				case 12: // force and U
					EInteract_MU(esv, c1, mu1, c2, mu2, dr, R, R2, F, U, (char)bExcludeDipole, true);
					cell->Fe.m[i3a1] -= F[0]; cell->Fe.m[i3a1 + 1] -= F[1]; cell->Fe.m[i3a1 + 2] -= F[2];
					Uestat -= U * 0.5;
					// virial
					cell->ve.m[i3a1] -= F[0] * dr[0] * 0.5;
					cell->ve.m[i3a1 + 1] -= F[1] * dr[1] * 0.5;
					cell->ve.m[i3a1 + 2] -= F[2] * dr[2] * 0.5;
					break;
				default:
					break;
				}
			}
			cell->U_estat.m[ia1] += Uestat; // Uestat is excluded energy
		}
	}
#if _CHECK_TIME_STAMP_ == 1
	//t_LJ += time.elapse();
#endif

}

void GPU_EwdSumReal_INTERACT_BOUND(JobIndx *job, EAtCELL *cell) {
	int ib, ia1, ia2, i3a1, i3a2, i3b;//, cIndx1, cIndx2;
	int ix1, iy1, iz1, ix2, iy2, iz2, dx, dy, dz;
	//int Nx, Ny, Nz;
	int Nm;

	int hx = cell->Nx / 2, hy = cell->Ny, hz = cell->Nz / 2;

	double R, R2;
	double t1, t2, t3;
	double xd[3] = {cell->cw[0] * cell->Nx, cell->cw[1] * cell->Ny, cell->cw[2] * cell->Nz};
	double doff[3], E[3], F[3], mu1[3], mu2[3];

	double c1, c2;
	double U = 0, Ut = 0;
	double dr[3];

	char bExcludeDipole = (char)cell->esv.iExIndDipole;
	int cal_mode = cell->bCalE;
	if (cal_mode <= 0 || cal_mode >= 3) return; // with charge only, so far choice 1 -- electric field; 2 -- force and energy
	if (cell->bDipole) cal_mode += 10; // with dipole

	EwaldSumRealVars1 esv;
	esv.cp_EwaldSum_vars(*(EwaldSumRealVars0*)(&cell->esv));
	esv.bEwald = false;


#if _CHECK_TIME_STAMP_ == 1
	TIME time;
	time.start();
#endif

	//V6zero(pc->dyn->fc)

	for (ib = job->n1; ib <= job->n2; ib++) {
		Nm = ib * 2; i3b = x_idx(ib);
		ia1 = cell->bound.m[Nm]; ia2 = cell->bound.m[Nm + 1]; i3a1 = x_idx(ia1); i3a2 = x_idx(ia2);
		c1 = cell->c.m[ia1]; c2 = cell->c.m[ia2];
		if (cell->cIndx.m[ia1] == cell->cIndx.m[ia2]) continue; // they are in the same cluster, which was ignored in short-range interaction calculation
		
		ix1 = cell->ir.m[i3a1]; iy1 = cell->ir.m[i3a1 + 1]; iz1 = cell->ir.m[i3a1 + 2];
		ix2 = cell->ir.m[i3a2]; iy2 = cell->ir.m[i3a2 + 1]; iz2 = cell->ir.m[i3a2 + 2];
		dx = ix2 - ix1; dy = iy2 - iy1; dz = iz2 - iz1;
		if (dx > hx) doff[0] = -xd[0]; else if (dx < -hx) doff[0] = xd[0]; else doff[0] = 0;
		if (dy > hy) doff[1] = -xd[1]; else if (dy < -hy) doff[1] = xd[1]; else doff[1] = 0;
		if (dz > hz) doff[2] = -xd[2]; else if (dz < -hz) doff[2] = xd[2]; else doff[2] = 0;

		dr[0] = cell->r.m[i3a2] - cell->r.m[i3a1] + doff[0];
		dr[1] = cell->r.m[i3a2 + 1] - cell->r.m[i3a1 + 1] + doff[1];
		dr[2] = cell->r.m[i3a2 + 2] - cell->r.m[i3a1 + 2] + doff[2];
		R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
		R = sqrt(R2);
		//if (R > cell->sr_cut) continue;
		//u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

		U = 0;

		// calculate electrostatic interaction
		switch (cal_mode) {
		case 1: // electric field only
			//MTP_E_Q(esv, c2, dr, R, R2, E);
			EInteract_E_Q(esv, c2, dr, R, R2, E, true);
			cell->E_b.m[i3b] = E[0]; cell->E_b.m[i3b + 1] = E[1]; cell->E_b.m[i3b + 2] = E[2];
			break;
		case 2: // force and U
			//MTP_Interact_Q(esv, c1, c2, dr, R, R2, F, U);
			EInteract_Q(esv, c1, c2, dr, R, R2, F, U, true);
			cell->Fe_b.m[i3b] = F[0]; cell->Fe_b.m[i3b + 1] = F[1]; cell->Fe_b.m[i3b + 2] = F[2];
			cell->Ub_estat.m[ib] = U;
			// virial
			cell->ve_b.m[i3b] = F[0] * dr[0];
			cell->ve_b.m[i3b + 1] = F[1] * dr[1];
			cell->ve_b.m[i3b + 2] = F[2] * dr[2];
			break;
		case 11: 
			//MTP_E_MU(esv, c2, mu2, dr, R, R2, E);
			EInteract_E_MU(esv, c2, mu2, dr, R, R2, E, true);
			cell->E_b.m[i3b] = E[0]; cell->E_b.m[i3b + 1] = E[1]; cell->E_b.m[i3b + 2] = E[2];
			break;
		case 12: // force and U
			//MTP_Interact_MU(esv, c1, mu1, c2, mu2, dr, R, R2, F, U, bExcludeDipole);
			EInteract_MU(esv, c1, mu1, c2, mu2, dr, R, R2, F, U, bExcludeDipole, true);
			cell->Fe_b.m[i3b] = F[0]; cell->Fe_b.m[i3b + 1] = F[1]; cell->Fe_b.m[i3b + 2] = F[2];
			cell->Ub_estat.m[ib] = U;
			// virial
			cell->ve_b.m[i3b] = F[0] * dr[0];
			cell->ve_b.m[i3b + 1] = F[1] * dr[1];
			cell->ve_b.m[i3b + 2] = F[2] * dr[2];
			break;
		default:
			break;
		}
	}
#if _CHECK_TIME_STAMP_ == 1
	//t_LJ += time.elapse();
#endif

}


// 2d Ewald-Sum, local part only

void EwaldSumRecpKzero_E_Q(EwaldSumRecpK0Var &var, double q2, double &Ez) {
	Ez = -q2 * var.ferf * var.PI2_over_A;
}

void EwaldSumRecpKzero_E_MU(EwaldSumRecpK0Var &var, double q2, double mu2z, double &Ez) {
	Ez = -(q2 * var.ferf + mu2z * var.KAPPA2_over_SQRT_PI * var.fexp) * var.PI2_over_A;
}


void EwaldSumRecpKzero_Q(EwaldSumRecpK0Var &var, double q1, double q2, /*double &Ez, */double &Fz, double &U) {
	//Ez = -q2 * var.ferf * var.PI2_over_A;
	//if (var.bEfieldOnly) return;

	double q12 = q1 * q2;
	U = -q12 * var.W * var.PI_over_A; // half of interaction
	Fz = -q12 * var.ferf * var.PI2_over_A;
}

void EwaldSumRecpKzero_MU(EwaldSumRecpK0Var &var, double q1, double mu1z, double q2, double mu2z, /*double &Ez, */double &Fz, double &U) {
	//Ez = -(q2 * var.ferf + mu2z * var.KAPPA2_over_SQRT_PI * var.fexp) * var.PI2_over_A;
	//Ez = -q2 * var.ferf;
	//if (var.fexp != 0) {
	//	Ez -= mu2z * var.KAPPA2_over_SQRT_PI * var.fexp;
	//}
	//Ez *= var.PI2_over_A;
	//if (var.bEfieldOnly) return;

	double q12 = q1 * q2;
	if (var.iExIndDipole == 1) {
		U = -(q12 * var.W + (q1 * mu2z - mu1z * q2) * var.ferf/* - mu1z * mu2z * var.KAPPA2_over_SQRT_PI * var.fexp*/) * var.PI_over_A; // half of interaction
		Fz = -(q12 * var.ferf + (q1 * mu2z - mu1z * q2) * var.KAPPA2_over_SQRT_PI * var.fexp /*+ mu1z * mu2z * var.KAPPA2_over_SQRT_PI * 2 * var.kappa2 * var.dz * var.fexp*/) * var.PI2_over_A;
	}
	else {
		U = -(q12 * var.W + (q1 * mu2z - mu1z * q2) * var.ferf - mu1z * mu2z * var.KAPPA2_over_SQRT_PI * var.fexp) * var.PI_over_A; // half of interaction
		Fz = -(q12 * var.ferf + (q1 * mu2z - mu1z * q2) * var.KAPPA2_over_SQRT_PI * var.fexp + mu1z * mu2z * var.KAPPA2_over_SQRT_PI * 2 * var.kappa2 * var.dz * var.fexp) * var.PI2_over_A;
	}
}

// for 2d system, long-range interaction in Ewald Sum has non-zero term for k = 0
void GPU_EwaldSum_2D_K0(JobIndx *job, EAtCELL *cell) {
	EwaldSumRecpK0Var var = cell->k0var;

	double q1, q2, muz1, muz2;
	double Ez, Fz, Ezt, Fzt;
	double dz, z1;
	int ia1, ia2, i3az1, i3az2;

	double t1;

	double Uestat = 0, STzz, STxx, STyy;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	int cal_mode = cell->bCalE;
	var.bEfieldOnly = (cal_mode == 1 ? true : false);
	if (cal_mode <= 0 || cal_mode >= 3) return; // with charge only, so far choice 1 -- electric field; 2 -- force and energy
	if (cell->bDipole) cal_mode += 10; // with dipole

	//char msg[256] = "\0";

#if _CHECK_TIME_STAMP_ == 1
	time.start();
#endif

	for (ia1 = job->n1; ia1 <= job->n2; ia1++) {
		//if (eNeutral) continue;
		i3az1 = z_idx(ia1);
		q1 = cell->c.m[ia1]; muz1 = cell->mu.m[i3az1]; z1 = cell->r.m[i3az1];

		Uestat = 0; Ez = 0; Fz = 0; STxx = 0; STyy = 0; STzz = 0;

		for (ia2 = 0; ia2 < cell->Na; ia2++) {
			//if (patom1->eNeutral) continue; // ignore this atom for charge interaction
			i3az2 = z_idx(ia2);
			dz = cell->r.m[i3az2] - z1;
			q2 = cell->c.m[ia2]; muz2 = cell->mu.m[i3az2];

			var.init_var(dz);
			switch (cal_mode) {
			case 1: // electric field only with q
				EwaldSumRecpKzero_E_Q(var, q2, Ezt);
				Ez += Ezt;
				break;
			case 2: // force and energy with q
				EwaldSumRecpKzero_Q(var, q1, q2, Fzt, t1); //half interaction energy is calculated
				Fz += Fzt; Uestat += t1;
				STzz -= Fzt * dz * 0.5; // * fUnit_estat_mass;
				STxx += t1; STyy += t1;
				break;
			case 11: // electric field only with q and mu
				EwaldSumRecpKzero_E_MU(var, q2, muz2, Ezt);
				Ez += Ezt;
				break;
			case 12: // force and energy with q and mu
				EwaldSumRecpKzero_MU(var, q1, muz1, q2, muz2, Fzt, t1); //half interaction energy is calculated
				Fz += Fzt; Uestat += t1;
				STzz -= Fzt * dz * 0.5; // * fUnit_estat_mass;
				STxx += t1; STyy += t1;
				break;
			default:
				break;
			}
		}
		
		switch (cal_mode) {
		case 1: // electric field only with q
			cell->E.m[i3az1] += Ez;
			break;
		case 2: // force and energy with q
			cell->U_estat.m[ia1] += Uestat;
			cell->Fe.m[i3az1] += Fz;
			cell->ve.m[i3az1 - 2] += STxx;
			cell->ve.m[i3az1 - 1] += STyy;
			cell->ve.m[i3az1] += STzz;
			break;
		case 11: // electric field only with q and mu
			cell->E.m[i3az1] += Ez;
			break;
		case 12: // force and energy with q and mu
			cell->U_estat.m[ia1] += Uestat;
			cell->Fe.m[i3az1] += Fz;
			cell->ve.m[i3az1 - 2] += STxx;
			cell->ve.m[i3az1 - 1] += STyy;
			cell->ve.m[i3az1] += STzz;
			break;
		default:
			break;
		}
	}

#if _CHECK_TIME_STAMP_ == 1
	time.stop();
#endif

}




/*
void collect_efield(JobIndx* job, AtCELL *cell, ARRAY<MATOM*> &atom) {
	int ia;
	for (ia = job->n1; ia <= job->n2; ia++) {
		atom.m[ia]->E.v[0] += cell->Ex.m[ia];
		atom.m[ia]->E.v[1] += cell->Ey.m[ia];
		atom.m[ia]->E.v[2] += cell->Ez.m[ia];
	}
}
*/
void collect_sr_fresult(JobIndx* job, SR_AtCELL *cell) {
	int ia, i3a;
	MATOM *atom;
	double *F;
	for (ia = job->n1; ia <= job->n2; ia++) {
		i3a = x_idx(ia);
		atom = cell->satom.m[ia].m;
		if (::bMTForce) F = atom->f.m[_IF_SR].v;
		else F = atom->F.v;
		F[0] += cell->Fsr.m[i3a]; // has unit in mass
		F[1] += cell->Fsr.m[i3a + 1];
		F[2] += cell->Fsr.m[i3a + 2];

		atom->bShocked = (cell->bShocked.m[ia] == 0x00 ? false : true) ;
	}
}

void collect_Uvirial(JobIndx* job, EAtCELL *cell, SVECTOR<double, 4> *res) {
	int ia, iThread = job->iThread, i3a;
	memset(res[iThread].v, 0, 4 * SIZE_DOUBLE);
	/*
	if (cell->bCalSR == 1) {
		for (ia = job->n1; ia <= job->n2; ia++) {
			i3a = x_idx(ia);
			res[iThread].v[0] += cell->U_sr.m[ia]; // has unit in mass
			res[iThread].v[1] += cell->vsr.m[i3a];
			res[iThread].v[2] += cell->vsr.m[i3a + 1];
			res[iThread].v[3] += cell->vsr.m[i3a + 2];
		}
	}
	*/

	if (cell->bCalE == 2) {
		for (ia = job->n1; ia <= job->n2; ia++) {
			i3a = x_idx(ia);
			res[iThread].v[0] += cell->U_estat.m[ia] * ::eUnit_estat_kT; // has unit in mass
			res[iThread].v[1] += cell->ve.m[i3a] * ::fUnit_estat_mass;
			res[iThread].v[2] += cell->ve.m[i3a + 1] * ::fUnit_estat_mass;
			res[iThread].v[3] += cell->ve.m[i3a + 2] * ::fUnit_estat_mass;
		}
	}
}

void collect_Uvirial_sr(JobIndx* job, SR_AtCELL *cell, SVECTOR<double, 4> *res) {
	int ia, iThread = job->iThread, i3a;
	memset(res[iThread].v, 0, 4 * SIZE_DOUBLE);
	if (cell->bCalSR == 1) {
		for (ia = job->n1; ia <= job->n2; ia++) {
			i3a = x_idx(ia);
			res[iThread].v[0] += cell->U_sr.m[ia]; // has unit in mass
			res[iThread].v[1] += cell->vsr.m[i3a];
			res[iThread].v[2] += cell->vsr.m[i3a + 1];
			res[iThread].v[3] += cell->vsr.m[i3a + 2];
		}
	}
}

#ifdef __GPU__
void cudaEwaldSum2d_K0(EAtCELL *acell) { // asynchrotron thread-function
	if (cudaSetDevice(acell->igpuDev) != cudaSuccess) {
		show_infor("GPU is not enabled!", true); return;
	}
	acell->gpu_prior_cal(2);
	cudaStream_t cudaStream;
	cudaStreamCreate(&cudaStream);
	_gpu_md_::GPU_EwaldSumReal_2d(&(acell->gpu_cell), acell->gpu_cell_dev, cudaStream);
	cudaStreamSynchronize(cudaStream);
	cudaStreamDestroy(cudaStream);
}

void cpuEwaldSum2d_K0(EAtCELL *cell) {
	int nThreads = MAX_THREADS, nthreads = nThreads;
	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	
	// ***************** following should be done in GPU  ***************

	//if (cell->bCalE == 1 || cell->bCalE == 11) cell->reset_efield(); // efield calculation
	//else if (cell->bCalE == 2 || cell->bCalE == 12) cell->reset_eforce(); // force and energy calculation

	assignJobIndx(job, nThreads, 0, cell->Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, cell, true);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwaldSum_2D_K0), job, nThreads, cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell->Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, cell, true);
	
	if (cell->Nb > 0) { // bound
		if (cell->Nb > MAX_THREADS * 10) {
			nthreads = -1; assignJobIndx1(job, nthreads, 0, cell->Nb - 1, 50, nThreads);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nthreads, cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell->Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, cell);}
	}
	// need to do 1-4 neighbors here

	/*
	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
	*/
	// ***************** end of jobs done in GPU  ***************
}


extern void cudaEwaldSumReal_3d(EAtCELL *acell); // asynchrotron thread-function
extern void cpuEwaldSumReal_3d(EAtCELL *cell); // using CPU 
extern void gpu_reset_esr_neighbors(EAtCELL *acell);

#endif

void LocalF(EAtCELL &cell, _spme_2d_::VSPME<_evatom_::_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	res.reset();
	InteractRes tres;

	cell.bCalE = 2; // force only

	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	
	// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
	int nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
	if (cell.bCalE == 1 || cell.bCalE == 11) MTOperate1<JobIndx, _spme_2d_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::reset_global_efield_q_3d), job, nthreads, &vspme, true);
	else if (cell.bCalE == 2 || cell.bCalE == 12) MTOperate1<JobIndx, _spme_2d_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::reset_global_force_q_3d), job, nthreads, &vspme, true);

	// ***************** following should be done in GPU  ***************

	if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
	else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

#ifdef __GPU__
	// calculate short-range interaction with GPU
	cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.esv1));
	AsynRun<EAtCELL> asynRun;
	asynRun.set_par(&cell);
	if (cell.bCalUseGPU) asynRun.start_func((void*)&cudaEwaldSumReal_3d);
	else asynRun.start_func((void*)&cpuEwaldSumReal_3d);
	//asynRun.wait();
#else
	assignJobIndx(job, nThreads, 0, cell.Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, &cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell.Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, &cell, true);
	
	if (cell.Nb > 0) { // bound
		if (cell.Nb > MAX_THREADS * 10) {
			nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Nb - 1, 200, nThreads);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nthreads, &cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell.Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, &cell);}
	}
	// need to do 1-4 neighbors here

	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif
	// ***************** end of jobs done in GPU  ***************

#ifdef __GPU__ // we have to wait for the GPU result on real part ewald-sum
	asynRun.wait();
	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif

	// then we needs to dump the efield or force of local interaction
	SVECTOR<double, 4> sr_res[MAX_THREADS];
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_2d_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::collect_efield_sr_q_3d), job, nthreads, &cell, &vspme, true);
	else {
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_2d_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::collect_fresult_sr_q_3d), job, nthreads, &cell, &vspme, true);
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, SVECTOR<double, 4> >((void*)(&collect_Uvirial), job, nthreads, &cell, sr_res, true);
		for (i = 0; i < nthreads; i++) {
			res.U += sr_res[i].v[0];
			res.STxx += sr_res[i].v[1];
			res.STyy += sr_res[i].v[2];
			res.STzz += sr_res[i].v[3];
			res.STtr += sr_res[i].v[1] + sr_res[i].v[2] + sr_res[i].v[3];
		}
	}
}

void LocalF(EAtCELL &cell, _spme_2d_::VSPME<_evatom_::_VMUatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res) {
	res.reset();
	InteractRes tres;

	cell.bCalE = 2; // force only

	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	
	// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
	int nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
	if (cell.bCalE == 1 || cell.bCalE == 11) MTOperate1<JobIndx, _spme_2d_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::reset_global_efield_mu_3d), job, nthreads, &vspme, true);
	else if (cell.bCalE == 2 || cell.bCalE == 12) MTOperate1<JobIndx, _spme_2d_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::reset_global_force_mu_3d), job, nthreads, &vspme, true);

	// ***************** following should be done in GPU  ***************

	if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
	else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

#ifdef __GPU__
	// calculate short-range interaction with GPU
	cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.esv1));
	AsynRun<EAtCELL> asynRun;
	asynRun.set_par(&cell);
	if (cell.bCalUseGPU) asynRun.start_func((void*)&cudaEwaldSumReal_3d);
	else asynRun.start_func((void*)&cpuEwaldSumReal_3d);
	//asynRun.wait();
#else
	assignJobIndx(job, nThreads, 0, cell.Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, &cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell.Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, &cell, true);
	
	if (cell.Nb > 0) { // bound
		if (cell.Nb > MAX_THREADS * 10) {
			assignJobIndx(job, nThreads, 0, cell.Nb - 1);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nThreads, &cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell.Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, &cell);}
	}
	// need to do 1-4 neighbors here

	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif
	// ***************** end of jobs done in GPU  ***************

#ifdef __GPU__ // we have to wait for the GPU result on real part ewald-sum
	asynRun.wait();
	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif


	// then we needs to dump the efield or force of local interaction
	SVECTOR<double, 4> sr_res[MAX_THREADS];
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_2d_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::collect_efield_sr_mu_3d), job, nthreads, &cell, &vspme, true);
	else {
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_2d_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::collect_fresult_sr_mu_3d), job, nthreads, &cell, &vspme, true);
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, SVECTOR<double, 4> >((void*)(&collect_Uvirial), job, nthreads, &cell, sr_res, true);
		for (i = 0; i < nthreads; i++) {
			res.U += sr_res[i].v[0];
			res.STxx += sr_res[i].v[1];
			res.STyy += sr_res[i].v[2];
			res.STzz += sr_res[i].v[3];
			res.STtr += sr_res[i].v[1] + sr_res[i].v[2] + sr_res[i].v[3];
		}
	}
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(EAtCELL &cell, _spme_2d_::VSPME<_evatom_::_VQatom>& vspme, int nThreads, _spme_2d_::SPME_VAR& var, InteractRes& res) {
	res.reset();
	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	
	// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
	int nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
	if (cell.bCalE == 1 || cell.bCalE == 11) MTOperate1<JobIndx, _spme_2d_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::reset_global_efield_q_2d), job, nthreads, &vspme, true);
	else if (cell.bCalE == 2 || cell.bCalE == 12) MTOperate1<JobIndx, _spme_2d_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::reset_global_force_q_2d), job, nthreads, &vspme, true);

	// ***************** following should be done in GPU  ***************

	if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
	else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

#ifdef __GPU__
	// calculate short-range interaction with GPU
	cell.k0var = var.K0var;
	cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.spme_3d_var.esv1));
	AsynRun<EAtCELL> asynRun;
	asynRun.set_par(&cell);
	if (cell.bCalUseGPU) asynRun.start_func((void*)&cudaEwaldSum2d_K0);
	else asynRun.start_func((void*)&cpuEwaldSum2d_K0);
	//asynRun.wait();
#else
	assignJobIndx(job, nThreads, 0, cell.Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, &cell, true);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwaldSum_2D_K0), job, nThreads, &cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell.Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, &cell, true);
	
	if (cell.Nb > 0) { // bound
		if (cell.Nb > MAX_THREADS * 10) {
			nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Nb - 1, 50, nThreads);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nthreads, &cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell.Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, &cell);}
	}
	// need to do 1-4 neighbors here

	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif
	// ***************** end of jobs done in GPU  ***************



	// accumulate the electric field from the image part
	_spme_2d_::init_normalized_coordinates(vspme, nThreads);
	_spme_2d_::cal_Q_spme(vspme, nThreads); // SPME Q matrix
	double U = vspme.cal_U_recp(!var.spme_3d_var.bEfieldOnly); //calculate the energy from the image part 
	res.U += U;
	_spme_2d_::cal_E_recp(vspme, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		_spme_2d_::cal_TQ_spme(vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.cal_FTQ();
		vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {res.STxx += vspme.STxx; res.STyy += vspme.STyy; res.STzz += vspme.STzz;}
		res.STtr += vspme.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0]: %f kT", U); show_log(errmsg, true);

	// virial coming from the surface dipole, U(V), from receiprcal part
	InteractRes tres;
	if (!var.spme_3d_var.bEfieldOnly && var.bVirial) {
		EwaldSum_SurfaceDipole(cell.mu_surf, (EwaldSumRealVars0*)(&(var.spme_3d_var.esv1)), tres);
		res += tres;
	}


	//sprintf(errmsg, "After surface dipole correction, virial changes to %f", res.STtr); show_log(errmsg, true);


	bool bDumpF = (cell.bCalE == 2 ? true : false);
	bool bDumpE = (cell.bCalE == 1 ? true : false);
	if (nThreads < 2 || vspme.va.n < N4PARALLEL) vspme.dump_EF(bDumpF, true, bDumpE); // transfer the calculated E & F to the variables in real atom
	else {
		nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
		MTOperate3< JobIndx, _spme_2d_::VSPME<_evatom_::_VQatom>, bool, bool>((void*)(&_spme_2d_::Slab_dump_EF_Q_Job), job, nthreads, &vspme, &bDumpE, &bDumpF, true);
	}




#ifdef __GPU__ // we have to wait for the GPU result on real part ewald-sum
	asynRun.wait();
	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif


	// then we needs to dump the efield or force of local interaction
	SVECTOR<double, 4> sr_res[MAX_THREADS];
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (var.spme_3d_var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_2d_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::collect_efield_sr_q_2d), job, nthreads, &cell, &vspme, true);
	else {
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_2d_::VSPME<_evatom_::_VQatom> >((void*)(&_atomic_cmm_::collect_fresult_sr_q_2d), job, nthreads, &cell, &vspme, true);
		MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, SVECTOR<double, 4> >((void*)(&collect_Uvirial), job, nthreads, &cell, sr_res, true);
		for (i = 0; i < nthreads; i++) {
			res.U += sr_res[i].v[0];
			res.STxx += sr_res[i].v[1];
			res.STyy += sr_res[i].v[2];
			res.STzz += sr_res[i].v[3];
			res.STtr += sr_res[i].v[1] + sr_res[i].v[2] + sr_res[i].v[3];
		}
	}
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum, to calculate Electric field and F on each atom */
void SPME_EF(EAtCELL &cell, _spme_2d_::VSPME<_evatom_::_VMUatom>& vspme, int nThreads, _spme_2d_::SPME_VAR& var, InteractRes& res, bool bInitCoordination) {
	res.reset();
	int i;

	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	int nthreads;

	// reset the electric field or force on the real atom atom, for electrostatic interaction, which is recorded on vspme.va
	nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
	if (cell.bCalE == 1 || cell.bCalE == 11) MTOperate1<JobIndx, _spme_2d_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::reset_global_efield_mu_2d), job, nthreads, &vspme, true);
	else if (cell.bCalE == 2 || cell.bCalE == 12) MTOperate1<JobIndx, _spme_2d_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::reset_global_force_mu_2d), job, nthreads, &vspme, true);

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
	time.start();
#endif

	// ***************** following should be done in GPU  ***************

	if (cell.bCalE == 1 || cell.bCalE == 11) cell.reset_efield(); // efield calculation
	else if (cell.bCalE == 2 || cell.bCalE == 12) cell.reset_eforce(); // force and energy calculation

#ifdef __GPU__
	// calculate short-range interaction with GPU
	cell.k0var = var.K0var;
	cell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&var.spme_3d_var.esv1));
	AsynRun<EAtCELL> asynRun;
	asynRun.set_par(&cell);
	if (cell.bCalUseGPU) asynRun.start_func((void*)&cudaEwaldSum2d_K0);
	else asynRun.start_func((void*)&cpuEwaldSum2d_K0);
	//asynRun.wait();
#else
	assignJobIndx(job, nThreads, 0, cell.Na - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_ATOM), job, nThreads, &cell, true);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwaldSum_2D_K0), job, nThreads, &cell, true);
	// excluding the electrostatic interaction between the atoms in the same cluster
	assignJobIndx(job, nThreads, 0, cell.Nc - 1);
	MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_InCluster), job, nThreads, &cell, true);

	if (cell.Nb > 0) { // bound
		if (cell.Nb > MAX_THREADS * 10) {
			nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Nb - 1, 100, nThreads);
			MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_EwdSumReal_INTERACT_BOUND), job, nthreads, &cell, true);
		}
		else {job[0].n1 = 0; job[0].n2 = cell.Nb - 1; GPU_EwdSumReal_INTERACT_BOUND(job, &cell);}
	}
	// need to do 1-4 neighbors here

	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, EAtCELL>((void*)(&GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif
	// ***************** end of jobs done in GPU  ***************

#if _CHECK_TIME_STAMP_ == 1
	long dt = time.elapse();
#endif

	// accumulate the electric field from the image part
	if (bInitCoordination) _spme_2d_::init_normalized_coordinates(vspme, nThreads);
	_spme_2d_::cal_Q_spme(vspme, nThreads); // SPME Q matrix
	double U = vspme.cal_U_recp(!var.spme_3d_var.bEfieldOnly); //calculate the energy from the image part 
	res.U += U;
	_spme_2d_::cal_E_recp(vspme, !var.K0var.bEfieldOnly, nThreads); // calculate the electric field, and force for each atom in reciprocal space

	if (!var.K0var.bEfieldOnly && var.bVirial) {
		_spme_2d_::cal_TQ_spme(vspme, nThreads); // SPME TQ matrix for stress sensor element STzz
		vspme.cal_FTQ();
		vspme.cal_StressTensor_recp_0();
		if (var.bST_diagonal) {res.STxx += vspme.STxx; res.STyy += vspme.STyy; res.STzz += vspme.STzz;}
		res.STtr += vspme.STtr;
	}
	//sprintf(errmsg, "Recip [k!=0]: %f kT", U); show_log(errmsg, true);

#if _CHECK_TIME_STAMP_ == 1
	long dt1 = time.elapse();
	dt1 -= dt;

	char msg[256] = "\0";
	sprintf(msg, "electrostatic interaction calculation: short-range used %d ms, reciprcal used %d ms.", dt, dt1); show_log(msg, true);
#endif


	// virial coming from the surface dipole, U(V), from receiprcal part
	InteractRes tres;
	if (!var.spme_3d_var.bEfieldOnly && var.bVirial) {
		EwaldSum_SurfaceDipole(cell.mu_surf, (EwaldSumRealVars0*)(&(var.spme_3d_var.esv1)), tres);
		res += tres;
	}


	//sprintf(errmsg, "After surface dipole correction, virial changes to %f", res.STtr); show_log(errmsg, true);


	//bool bDumpF = !var.K0var.bEfieldOnly;
	bool bDumpF = (cell.bCalE == 2 ? true : false);
	bool bDumpE = (cell.bCalE == 1 ? true : false);
	if (nThreads < 2 || vspme.va.n < N4PARALLEL) vspme.dump_EF(bDumpF, true, bDumpE); // transfer the calculated E & F to the variables in real atom
	else {
		//MRunSys2< _spme_2d_::VSPME<_evatom_::_VMUatom>, bool>((void*)(&_spme_2d_::Slab_dump_EF_MU), vspme, 0, vspme.va.n, nThreads, bDumpF);
		nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme.va.n - 1, 400, nThreads);
		MTOperate3< JobIndx, _spme_2d_::VSPME<_evatom_::_VMUatom>, bool, bool>((void*)(&_spme_2d_::Slab_dump_EF_MU_Job), job, nthreads, &vspme, &bDumpE, &bDumpF, true);
	}





#ifdef __GPU__ // we have to wait for the GPU result on real part ewald-sum
	asynRun.wait();
	// accumulating the electric field or force in CELL, or GPU
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (cell.bCalE == 1) MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_efield), job, nthreads, &cell, true);
	else if (cell.bCalE == 2) { // force
		MTOperate1<JobIndx, _atomic_cmm_::EAtCELL>((void*)(&_atomic_cmm_::GPU_collecting_eforce), job, nthreads, &cell, true);
	}
#endif


	// then we needs to dump the efield or force of local interaction
	SVECTOR<double, 4> sr_res[MAX_THREADS];
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell.Na - 1, 400, nThreads);
	if (var.K0var.bEfieldOnly) MTOperate2<JobIndx, _atomic_cmm_::EAtCELL, _spme_2d_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::collect_efield_sr_mu_2d), job, nthreads, &cell, &vspme, true);
	else {
		MTOperate2<JobIndx, EAtCELL, _spme_2d_::VSPME<_evatom_::_VMUatom> >((void*)(&_atomic_cmm_::collect_fresult_sr_mu_2d), job, nthreads, &cell, &vspme, true);
		MTOperate2<JobIndx, EAtCELL, SVECTOR<double, 4> >((void*)(&collect_Uvirial), job, nthreads, &cell, sr_res, true);
		for (i = 0; i < nthreads; i++) {
			res.U += sr_res[i].v[0];
			res.STxx += sr_res[i].v[1];
			res.STyy += sr_res[i].v[2];
			res.STzz += sr_res[i].v[3];
			res.STtr += sr_res[i].v[1] + sr_res[i].v[2] + sr_res[i].v[3];
		}
	}
}

void dmax(double &x1, double &x2, double &res) {
	res = (x1 > x2 ? x1 : x2);
}

/* multi-thread function to calculate the interaction with SPME method for Ewald Sum */
bool Polarize_SPME_Interact(MMOL_MD_CELL &mdcell, EAtCELL &cell, _spme_2d_::VSPME<_evatom_::_VMUatom>& vspme, int nThreads, bool bPolarize, _spme_2d_::SPME_VAR &var, InteractRes &res) {
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
	cell.bCalE = 1; // electric field only
	cell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

#ifdef __GPU__
	gpu_reset_esr_neighbors(&cell);
#endif
	while (bPolarize && iloop < max_loops) {
		cal_dipole(mdcell, nThreads, mdcell.mu_surf); // accumulate dipole of each atom, and the whole mdcell
		memcpy(&(cell.mu_surf), &(mdcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
		_atomic_cmm_::update_atomic_charge_dipole(&cell);

		SPME_EF(cell, vspme, nThreads, var, tres, false); // calculate the dipole only

		// add the electric field from the atom inside the cluster
		//if (mdcell.sm.n > 0) MOperate<SMOLECULE>((void*)(&AtomicEFieldCorrect_SM), mdcell.sm.m, mdcell.sm.n, nThreads);

		// polarization
		if (mdcell.mm.n > 0) MOperate<MMOLECULE>((void*)(&polarize_MM), mdcell.mm.m, mdcell.mm.n, nThreads);
		if (mdcell.sm.n > 0) MOperate<SMOLECULE>((void*)(&polarize_SM), mdcell.sm.m, mdcell.sm.n, nThreads);
		if (mdcell.pm.n > 0) MOperate<PMOLECULE>((void*)(&polarize_PM), mdcell.pm.m, mdcell.pm.n, nThreads);
		dmu_max = 0; vt = 0;
		if (mdcell.ind_mu_mm.n > 0) MOperate2<MOL_POL<MMOLECULE>, double>((void*)(&MM_maxdiff_backup_polarization), mdcell.ind_mu_mm.m, mdcell.ind_mu_mm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell.ind_mu_sm.n > 0) MOperate2<MOL_POL<SMOLECULE>, double>((void*)(&SM_maxdiff_backup_polarization), mdcell.ind_mu_sm.m, mdcell.ind_mu_sm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;
		if (mdcell.ind_mu_pm.n > 0) MOperate2<MOL_POL<PMOLECULE>, double>((void*)(&PM_maxdiff_backup_polarization), mdcell.ind_mu_pm.m, mdcell.ind_mu_pm.n, nThreads, vt, false, (void*)(&dmax));
		dmu_max = (dmu_max > vt ? dmu_max : vt); vt = 0;

		iloop++;
		if (dmu_max <= muConvg && iloop > 0) {status = true; break;}
	}

	//char msg[256] = "\0";
	//if (status) sprintf(msg, "Convergent polarized dipole is obtainted with %d loops", iloop); 
	//else sprintf(msg, "No convergent polarized dipole is obtained with %d loops", iloop);
	//show_log(msg, true);

	cal_dipole(mdcell, nThreads, mdcell.mu_surf); // accumulate dipole of each atom, and the whole mdcell
	memcpy(&(cell.mu_surf), &(mdcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
	_atomic_cmm_::update_atomic_charge_dipole(&cell);
	if (bEfieldOnly) { // calculate f
		SPME_EF(cell, vspme, nThreads, var, res, false); // calculate the dipole and force
	}
	else { // calculate f
		var.K0var.bEfieldOnly = bEfieldOnly;
		var.spme_3d_var.bEfieldOnly = bEfieldOnly;
		cell.bCalE = (bEfieldOnly ? 1 : 2);

		SPME_EF(cell, vspme, nThreads, var, res, false); // calculate the dipole and force
	}

	return status;
}


// polarized SPME, with the interaction between induced dipole excluded
// in some polarized model, the interaction between induced dipole is excluded, e.g. only q-induced dipole is considered.
void VSPME_F_exclude_induced(MMOL_MD_CELL &mdcell, EAtCELL &cell, _spme_2d_::VSPME<_evatom_::_VMUatom>& vspme_induced, int nThreads, _spme_2d_::SPME_VAR &var, InteractRes &res, bool bInitABSP_buff) {
	// we have to calculate the real part of electro-static interaction, and LJ interaction
	// in which case, the electric field of each atom will be initialized/refreshed at the beginning
	// then we can calculate the image part of SPME Ewald Sum
	// and then to accumulate the force from electric field
	res.reset();
	InteractRes tres;

	// accumulate the electric field from the image part
	//vspme_induced.reset_local_EF(); // reset the variable for E & F
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
	//vspme.dump_EF(!var.K0var.bEfieldOnly); // transfer the calculated E & F to the variables in real atom
	//bool bDumpF = !var.K0var.bEfieldOnly;

	JobIndx job[MAX_THREADS];
	for (int i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	int nthreads;

	bool bDumpF = true; // useless
	bool bDumpE = false;
	//if (nThreads < 2 || vspme_induced.va.n < N4PARALLEL) Slab_dump_F_exclude_induced_mu(vspme_induced, 0, vspme_induced.va.n - 1, bDumpF); // transfer the calculated E & F to the variables in real atom
	//else MRunSys2< VSPME<_evatom_::_VMUatom>, bool >((void*)(&Slab_dump_F_exclude_induced_mu), vspme_induced, 0, vspme_induced.va.n, nThreads, bDumpF);
	if (nThreads < 2 || vspme_induced.va.n < N4PARALLEL) _spme_2d_::Slab_dump_F_exclude_induced_mu(vspme_induced, 0, vspme_induced.va.n - 1, bDumpF); // transfer the calculated E & F to the variables in real atom
	else {
		//MRunSys2< _spme_2d_::VSPME<_evatom_::_VMUatom>, bool>((void*)(&_spme_2d_::Slab_dump_EF_MU), vspme_induced, 0, vspme_induced.va.n, nThreads, bDumpF);
		nthreads = -1; assignJobIndx1(job, nthreads, 0, vspme_induced.va.n - 1, 400, nThreads);
		MTOperate1< JobIndx, _spme_2d_::VSPME<_evatom_::_VMUatom> >((void*)(&_spme_2d_::Slab_dumpJob_F_exclude_induced_mu), job, nthreads, &vspme_induced, true);
	}
}

bool Polarize_SPME_Interact2(MMOL_MD_CELL &mdcell, EAtCELL &cell, _spme_2d_::VSPME<_evatom_::_VMUatom>& vspme, _spme_2d_::VSPME<_evatom_::_VMUatom>& vspme_induced, int nThreads, bool bPolarize, _spme_2d_::SPME_VAR &var, InteractRes &res) {
	vspme_induced.release_plans();
	vspme.init_fft_plans();
#if _CHECK_TIME_STAMP_
	show_log("Polarize: ", false);
	show_vlog<int>(mt.elapse(), false);
#endif
	bool status = Polarize_SPME_Interact(mdcell, cell, vspme, nThreads, bPolarize, var, res);
	InteractRes tres;
	bool bABSP_Buff = false;
	if (bPolarize && var.spme_3d_var.esv1.iExIndDipole == 1) {
		vspme.release_plans();
		vspme_induced.init_fft_plans();
		if (vspme.absp.n == vspme_induced.absp.n) {
			cp_absp_buff((__AtomBSP_BUFF*)(&vspme), (__AtomBSP_BUFF*)(&vspme_induced));
			bABSP_Buff = true;
		}
		else bABSP_Buff = false;
		VSPME_F_exclude_induced(mdcell, cell, vspme_induced, nThreads, var, tres, !bABSP_Buff);
		res -= tres;
	}
#if _CHECK_TIME_STAMP_
	show_log("", true);
#endif
	return status;
}


#ifdef __GPU__
void GPU_LJ_Interact(SR_AtCELL *cell, int nThreads, InteractRes &res) {
	res.reset();
	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	int nthreads;

	// ********* following should run on GPU  ****
	cell->reset_sr();
	/*
	assignJobIndx(job, nThreads, 0, cell->Na - 1);
	MTOperate1<JobIndx, SR_AtCELL>((void*)(&GPU_LJ_INTERACT_ATOM), job, nThreads, cell, true);
	
	if (cell->Nb14 > 0) {
		if (cell->Nb14 < 10) {job[0].n1 = 0; job[0].n2 = cell->Nb14 - 1; GPU_LJ_INTERACT_BOUND14(job, cell);}
		else {
			nthreads = -1; assignJobIndx1(job, nthreads, 0, cell->Nb14 - 1, 100, nThreads);
			MTOperate1<JobIndx, SR_AtCELL>((void*)(&GPU_LJ_INTERACT_BOUND14), job, nthreads, cell, true);
		}
	}
	
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell->Na - 1, 400, nThreads);
	MTOperate1<JobIndx, SR_AtCELL>((void*)(&GPU_collecting_srforce), job, nthreads, cell, true);
	*/
	GPU_LJ_interact(&(cell->gpu_cell), cell->gpu_cell_dev);
	// ********* end here on GPU  ****

	SVECTOR<double, 4> sr_res[MAX_THREADS];
	// collecting force to ATOM, and sum the energy and virial coefficients
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell->Na - 1, 400, nThreads);
	MTOperate1<JobIndx, _atomic_cmm_::SR_AtCELL>((void*)(&collect_sr_fresult), job, nthreads, cell, true);
	MTOperate2<JobIndx, _atomic_cmm_::SR_AtCELL, SVECTOR<double, 4> >((void*)(&collect_Uvirial_sr), job, nthreads, cell, sr_res, true);
	for (i = 0; i < nthreads; i++) {
		res.U += sr_res[i].v[0];
		res.STxx += sr_res[i].v[1];
		res.STyy += sr_res[i].v[2];
		res.STzz += sr_res[i].v[3];
		res.STtr += sr_res[i].v[1] + sr_res[i].v[2] + sr_res[i].v[3];
	}
}





void cudaLJ(AsynLJVar *var) { // asynchrotron thread-function
	if (cudaSetDevice(var->cell->igpuDev) != cudaSuccess) {
		show_infor("GPU is not enabled!", true); return;
	}
	var->cell->reset_sr();
	var->cell->gpu_prior_cal();
	// ********* following should run on GPU  ****
	GPU_LJ_interact(&(var->cell->gpu_cell), var->cell->gpu_cell_dev);
}

#endif //#ifdef __GPU__

//collecting LJ force calculated with cudaLJ above or cpuLJ below
void collect_LJ_res(SR_AtCELL *cell, int nThreads, InteractRes &res) {
	res.reset();
	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	int nthreads;

	SVECTOR<double, 4> sr_res[MAX_THREADS];
	// collecting force to ATOM, and sum the energy and virial coefficients
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell->Na - 1, 400, nThreads);
	MTOperate1<JobIndx, _atomic_cmm_::SR_AtCELL>((void*)(&collect_sr_fresult), job, nthreads, cell, true);
	MTOperate2<JobIndx, _atomic_cmm_::SR_AtCELL, SVECTOR<double, 4> >((void*)(&collect_Uvirial_sr), job, nthreads, cell, sr_res, true);
	for (i = 0; i < nthreads; i++) {
		res.U += sr_res[i].v[0];
		res.STxx += sr_res[i].v[1];
		res.STyy += sr_res[i].v[2];
		res.STzz += sr_res[i].v[3];
		res.STtr += sr_res[i].v[1] + sr_res[i].v[2] + sr_res[i].v[3];
	}
}


//#else
void LJ_Interact(SR_AtCELL *cell, int nThreads, InteractRes &res) {
	res.reset();
	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	int nthreads;

	// ********* following should run on GPU  ****
	cell->reset_sr();
	assignJobIndx(job, nThreads, 0, cell->Na - 1);
	MTOperate1<JobIndx, SR_AtCELL>((void*)(&GPU_LJ_INTERACT_ATOM), job, nThreads, cell, true);

	if (cell->Nb14 > 0) {
		if (cell->Nb14 < 10) {job[0].n1 = 0; job[0].n2 = cell->Nb14 - 1; GPU_LJ_INTERACT_BOUND14(job, cell);}
		else {
			nthreads = -1; assignJobIndx1(job, nthreads, 0, cell->Nb14 - 1, 100, nThreads);
			MTOperate1<JobIndx, SR_AtCELL>((void*)(&GPU_LJ_INTERACT_BOUND14), job, nthreads, cell, true);
		}
	}

	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell->Na - 1, 400, nThreads);
	MTOperate1<JobIndx, SR_AtCELL>((void*)(&GPU_collecting_srforce), job, nthreads, cell, true);
	// ********* end here on GPU  ****

	SVECTOR<double, 4> sr_res[MAX_THREADS];
	// collecting force to ATOM, and sum the energy and virial coefficients
	nthreads = -1; assignJobIndx1(job, nthreads, 0, cell->Na - 1, 400, nThreads);
	MTOperate1<JobIndx, _atomic_cmm_::SR_AtCELL>((void*)(&collect_sr_fresult), job, nthreads, cell, true);
	MTOperate2<JobIndx, _atomic_cmm_::SR_AtCELL, SVECTOR<double, 4> >((void*)(&collect_Uvirial_sr), job, nthreads, cell, sr_res, true);
	for (i = 0; i < nthreads; i++) {
		res.U += sr_res[i].v[0];
		res.STxx += sr_res[i].v[1];
		res.STyy += sr_res[i].v[2];
		res.STzz += sr_res[i].v[3];
		res.STtr += sr_res[i].v[1] + sr_res[i].v[2] + sr_res[i].v[3];
	}
}

void cpuLJ(AsynLJVar *var) { // asynchrotron thread-function
	var->res->reset();
	// ********* following should run on GPU  ****
	int i;
	JobIndx job[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
	int nthreads;

	// ********* following should run on GPU  ****
	var->cell->reset_sr();
	assignJobIndx(job, var->nThreads, 0, var->cell->Na - 1);
	MTOperate1<JobIndx, SR_AtCELL>((void*)(&GPU_LJ_INTERACT_ATOM), job, var->nThreads, var->cell, true);

	if (var->cell->Nb14 > 0) {
		if (var->cell->Nb14 < 10) {job[0].n1 = 0; job[0].n2 = var->cell->Nb14 - 1; GPU_LJ_INTERACT_BOUND14(job, var->cell);}
		else {
			nthreads = -1; assignJobIndx1(job, nthreads, 0, var->cell->Nb14 - 1, 100, var->nThreads);
			MTOperate1<JobIndx, SR_AtCELL>((void*)(&GPU_LJ_INTERACT_BOUND14), job, nthreads, var->cell, true);
		}
	}

	nthreads = -1; assignJobIndx1(job, nthreads, 0, var->cell->Na - 1, 400, var->nThreads);
	MTOperate1<JobIndx, SR_AtCELL>((void*)(&GPU_collecting_srforce), job, nthreads, var->cell, true);
	// ********* end here on GPU  ****
}

//****************************   3d Ewald-Sum  *************************


} // end of namespace
