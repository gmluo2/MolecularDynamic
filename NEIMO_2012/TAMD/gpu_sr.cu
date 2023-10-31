#include "project.h"

#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "def.h"
/*
#include "show.h"
#include "ranlib.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"
*/

#ifdef __GPU__
//#include "atom_sr.h"
#include "gpu_vector.h"
#include "gpu_sr.h"

namespace _gpu_md_ {
#define RBOUND_MAX  5   // maximum distance between 1-3 bound

	__device__ void standard_LJ(DOUBLE eps, DOUBLE rLJ, DOUBLE r, DOUBLE& U, DOUBLE& f) {
		// U = 4 * eps * [(rLJ / r)^12 - (rLJ / r)^6]
		// minimum U locates at r = pow(2, 1/6) * rLJ = 1.1225 * rLJ
		eps *= 4;
		DOUBLE r0 = rLJ / r;
		r0 = (r0 > r0LJ_max ? r0LJ_max : r0);
		DOUBLE r6 = pow(r0, 6), r12 = r6 * r6;
		U = eps * (r12 - r6);
		f = eps * (12 * r12 - 6 * r6) * r0 / rLJ;
	};

	__device__ void alternative_LJ(DOUBLE eps, DOUBLE rLJ, DOUBLE r, DOUBLE &U, DOUBLE &f) {
		// U = eps * [(rLJ / r)^12 - 2 * (rLJ / r)^6]
		// comparing to standard LJ, rLJ_alternative = 1.1225 * rLJ_standard, and the minimum of U locates at r = rLJ_alternative
		DOUBLE r0 = rLJ / r;
		r0 = (r0 > r0LJ_max ? r0LJ_max : r0);
		DOUBLE r6 = pow(r0, 6), r12 = r6 * r6;
		U = eps * (r12 - 2 * r6);
		f = 12 * eps * (r12 - r6) * r0 / rLJ;
	};

	__global__ void LJ_interact(GPU_SR_AtCELL *acell) {
		// works on 1 dimensional grid and one dimensional blocks
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int Na = acell->Na;
		int nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;


		int i, ia1, ia2;
		int ia, cIndx1, cIndx2, i3a1, i3a2;
		int icx, icy, icz, ncx, ncy, ncz, dcx, dcy, dcz;
		int Ncx, Ncy, Nc;

		int isr1, isr2, nsr = acell->npars;

		DOUBLE R, R2;
		//double t1, t2, t3;
		DOUBLE xd[3] = {acell->cw[0] * acell->Nx, acell->cw[1] * acell->Ny, acell->cw[2] * acell->Nz};
		DOUBLE doff[3];

		int *abound;
		bool bBound;		
		double r1min = 0, r2min = 0;

		DOUBLE eps, rLJ, U, Ut, f;
		DOUBLE force[3], ft[3], dr[3], u[3], vsr[3];

		DOUBLE rcut = acell->rcut, r2cut = rcut * rcut;
		DOUBLE fUnit_LJ_mass = acell->mgvar.fUnit_LJ_mass;
		DOUBLE eUnit_LJ_kT = acell->mgvar.eUnit_LJ_kT;
		
		dcx = int(rcut / acell->cw[0]) + 1;
		dcy = int(rcut / acell->cw[1]) + 1;
		dcz = int(rcut / acell->cw[2]) + 1;

		int iLJ = acell->sr_mode; //== ::iFormat_LJ, 0 -- standard L-J; 1 -- alternative L-J / Amber

		for (i = 0; i < nEachThread; i++) {
			ia1 = tid * nEachThread + i; i3a1 = x_idx(ia1);
			if (ia1 >= Na) break;
			icx = acell->ir_dev[i3a1]; icy = acell->ir_dev[i3a1 + 1]; icz = acell->ir_dev[i3a1 + 2];

			Ut = 0;
			//memset(ft, 0, 3 * sizeof(double)); 
			ft[0] = 0; ft[1] = 0; ft[2] = 0;
			//memset(vsr, 0, 3 * sizeof(double));
			vsr[0] = 0; vsr[1] = 0; vsr[2] = 0;
	
			if (acell->aID_dev[ia1] < 0) continue;
			cIndx1 = acell->cIndx_dev[ia1];
			isr1 = nsr * ia1;
			
			r1min = acell->sr_par_dev[isr1 + 2];

			for (icx = -dcx; icx <= dcx; icx++) {
				ncx = acell->ir_dev[i3a1] + icx; 
				if (ncx < 0) {ncx += acell->Nx; doff[0] = -xd[0];}
				else if (ncx >= acell->Nx) {ncx -= acell->Nx; doff[0] = xd[0];}
				else doff[0] = 0;
				Ncx = ncx * acell->Nyz;
				
				for (icy = -dcy; icy <= dcy; icy++) {
					ncy = acell->ir_dev[i3a1 + 1] + icy; 
					if (ncy < 0) {ncy += acell->Ny; doff[1] = -xd[1];}
					else if (ncy >= acell->Ny) {ncy -= acell->Ny; doff[1] = xd[1];}
					else doff[1] = 0;
					Ncy = ncy * acell->Nz;

					for (icz = -dcz; icz <= dcz; icz++) {
						ncz = acell->ir_dev[i3a1 + 2] + icz; 
						if (ncz < 0) {ncz += acell->Nz; doff[2] = -xd[2];}
						else if (ncz >= acell->Nz) {ncz -= acell->Nz; doff[2] = xd[2];}
						else doff[2] = 0;
						
						Nc = (Ncx + Ncy + ncz) * (acell->NC + 1);
						for (ia = 0; ia < acell->cmm_dev[Nc]; ia++) {
							ia2 = acell->cmm_dev[Nc + ia + 1];
							if (ia2 >= Na) continue;
							cIndx2 = acell->cIndx_dev[ia2];
							if (cIndx1 == cIndx2) continue; // in the same cluster
							if (acell->aID_dev[ia2] < 0) continue;
							isr2 = nsr * ia2;
							i3a2 = x_idx(ia2);
							
							dr[0] = acell->r_dev[i3a1] - acell->r_dev[i3a2] - doff[0];
							dr[1] = acell->r_dev[i3a1 + 1] - acell->r_dev[i3a2 + 1] - doff[1];
							dr[2] = acell->r_dev[i3a1 + 2] - acell->r_dev[i3a2 + 2] - doff[2];
							R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]; 
							if (R2 < 0.1) R2 = 0.1; 
							if (R2 > r2cut) continue;
							R = sqrt(R2);
							
							// check 1-2 and/or 1-3 bound and bShocked -- 2014/06/07
							if (R < RBOUND_MAX) {
								r2min = acell->sr_par_dev[isr2 + 2];
								bBound = false;
								abound = acell->device_abound(ia1);
								for (i = 0; i < abound[0]; i++) {
									if (abound[i+1] == ia2) {// 1-2 or 1-3 bound 
										bBound = true;
									}
								}
								// ignore bound interaction for short-range interaction
								if (bBound) continue;
								if (R < (r1min + r2min)) acell->bShocked_dev[ia1] = 0x01; // not bound atom
							}
							//******** 2014/06/07 *********
							
							u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

							// calculation of the interaction
							//memset(force, 0, SIZE_V3);
							if ((acell->sr_type_dev[ia1] & 0x01) && (acell->sr_type_dev[ia2] & 0x01)) { // L-J interaction
								eps = sqrt(acell->sr_par_dev[isr1] * acell->sr_par_dev[isr2]);
								rLJ = (acell->sr_par_dev[isr1 + 1] + acell->sr_par_dev[isr2 + 1]) * 0.5;
								
								if (iLJ == 0) standard_LJ(eps, rLJ, R, U, f);
								else alternative_LJ(eps, rLJ, R, U, f);
								
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
								Ut += U * 0.5;
								force[0] = f * u[0]; force[1] = f * u[1]; force[2] = f * u[2];
							}
							/*  so far we do not enable it
							if ((acell->sr_type[ia1] & 0x02) && (acell->sr_type[ia2] & 0x02)) { // given profile -- bSR
								aindx[0] = acell->aID_dev[ia1]; aindx[1] = acell->aID_dev.m[ia2];
								key = construct_char2_key_pairwise(aindx);
								pw = sr_db.search(key);
								if (pw != NULL) {
									interpolate2(t1, t2, pw->pv, R, i, t3)
									U_LJ += t1 * 0.5;
									if (t2 > force_max && scale_force) t2 = force_max; 
									force[0] += t2 * u[0]; force[1] += t2 * u[1]; force[2] += t2 * u[2];
								}
							}
							*/
							vsr[0] += force[0] * dr[0] * 0.5;
							vsr[1] += force[1] * dr[1] * 0.5;
							vsr[2] += force[2] * dr[2] * 0.5;

							ft[0] += force[0]; ft[1] += force[1]; ft[2] += force[2];
						}
					}
				}
			}

			acell->Fsr_dev[i3a1] = ft[0] * fUnit_LJ_mass; // force has mass unit, was converted in LJ12_6 profile
			acell->Fsr_dev[i3a1 + 1] = ft[1] * fUnit_LJ_mass;
			acell->Fsr_dev[i3a1 + 2] = ft[2] * fUnit_LJ_mass;

			acell->vsr_dev[i3a1] = vsr[0] * fUnit_LJ_mass;
			acell->vsr_dev[i3a1 + 1] = vsr[1] * fUnit_LJ_mass;
			acell->vsr_dev[i3a1 + 2] = vsr[2] * fUnit_LJ_mass;

			acell->U_sr_dev[ia1] = Ut * eUnit_LJ_kT;
		}
		
		// should we calculate SR of bound and bound14 ? for those of macromolecule

		//tid = threadIdx.x + blockIdx.x * blockDim.x;
		//Na = acell->Na;
		DOUBLE rmax_bound = 3, Rmax_boundLength;
		int Nb = acell->Nb, ib, i3b;
		Rmax_boundLength = 3 * 2; // bound include 1-2 and 1-3 bound, latter could be twice of closest bound

	// 1-2 and 1-3 bound are avoided in the code above		
#ifdef _OLD_VERSION_
		if (Nb > 0) {
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;
			for (i = 0; i < nEachThread; i++) {
				ib = tid * nEachThread + i; 
				if (ib > acell->Nb) break;
				i3b = x_idx(ib);
				//memset(acell->Fsr_b_dev + i3b, 0, SIZE_V3); 
				acell->Fsr_b_dev[i3b] = 0; acell->Fsr_b_dev[i3b + 1] = 0; acell->Fsr_b_dev[i3b + 2] = 0;
				//memset(acell->vsr_b_dev + i3b, 0, SIZE_V3); 
				acell->vsr_b_dev[i3b] = 0; acell->vsr_b_dev[i3b + 1] = 0; acell->vsr_b_dev[i3b + 2] = 0;
				acell->Ub_sr_dev[ib] = 0;
				ia1 = acell->bound_dev[ib * 2]; ia2 = acell->bound_dev[ib * 2 + 1];
				i3a1 = x_idx(ia1); i3a2 = x_idx(ia2);
				dr[0] = acell->r_dev[i3a1] - acell->r_dev[i3a2];
				dr[1] = acell->r_dev[i3a1 + 1] - acell->r_dev[i3a2 + 1];
				dr[2] = acell->r_dev[i3a1 + 2] - acell->r_dev[i3a2 + 2];
				if (dr[0] > Rmax_boundLength) dr[0] -= xd[0];
				else if (dr[0] < -Rmax_boundLength) dr[0] += xd[0];
				if (dr[1] > Rmax_boundLength) dr[1] -= xd[1];
				else if (dr[1] < -Rmax_boundLength) dr[1] += xd[1];
				if (dr[2] > Rmax_boundLength) dr[2] -= xd[2];
				else if (dr[2] < -Rmax_boundLength) dr[2] += xd[2];

				R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]; if (R2 < 0.2) R2 = 0.2;
				R = sqrt(R2);
			
				if (R > rcut) continue;
				u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

				// calculation of the interaction
				//memset(force, 0, SIZE_V3);
				if ((acell->sr_type_dev[ia1] & 0x01) && (acell->sr_type_dev[ia2] & 0x01)) { // L-J interaction
					eps = sqrt(acell->sr_par_dev[isr1] * acell->sr_par_dev[isr2]);
					rLJ = (acell->sr_par_dev[isr1 + 1] + acell->sr_par_dev[isr2 + 1]) * 0.5;
								
					if (iLJ == 0) standard_LJ(eps, rLJ, R, U, f);
					else alternative_LJ(eps, rLJ, R, U, f);

					Ut += U * 0.5;
					force[0] = f * u[0]; force[1] = f * u[1]; force[2] = f * u[2];
					/*  so far we do not enable it
					if ((acell->sr_type[ia1] & 0x02) && (acell->sr_type[ia2] & 0x02)) { // given profile -- bSR
						aindx[0] = acell->aID_dev[ia1]; aindx[1] = acell->aID_dev.m[ia2];
						key = construct_char2_key_pairwise(aindx);
						pw = sr_db.search(key);
						if (pw != NULL) {
							interpolate2(t1, t2, pw->pv, R, i, t3)
							U_LJ += t1 * 0.5;
							if (t2 > force_max && scale_force) t2 = force_max; 
							force[0] += t2 * u[0]; force[1] += t2 * u[1]; force[2] += t2 * u[2];
						}
					}
					*/
					acell->Fsr_b_dev[i3b] = force[0];
					acell->Fsr_b_dev[i3b + 1] = force[1];
					acell->Fsr_b_dev[i3b + 2] = force[2];
					//memcpy(acell->Fsr_b_dev + i3b, force, SIZE_V3);
					acell->vsr_b_dev[i3b] = force[0] * dr[0];
					acell->vsr_b_dev[i3b + 1] = force[1] * dr[1];
					acell->vsr_b_dev[i3b + 2] = force[2] * dr[2];
					acell->Ub_sr_dev[ib] = U;
				}
			}
		}
#endif // _OLD_VERSION

		// bound14
		Nb = acell->Nb14;
		Rmax_boundLength = 3 * 3; // bound 1-4
		if (Nb > 0) {
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;
			for (i = 0; i < nEachThread; i++) {
				ib = tid * nEachThread + i; 
				if (ib > acell->Nb14) break;
				i3b = x_idx(ib);
				//memset(acell->Fsr_14_dev + i3b, 0, SIZE_V3); 
				acell->Fsr_14_dev[i3b] = 0; acell->Fsr_14_dev[i3b + 1] = 0; acell->Fsr_14_dev[i3b + 2] = 0;
				//memset(acell->vsr_14_dev + i3b, 0, SIZE_V3); 
				acell->vsr_14_dev[i3b] = 0; acell->vsr_14_dev[i3b + 1] = 0; acell->vsr_14_dev[i3b + 2] = 0;
				acell->U14_sr_dev[ib] = 0;
				ia1 = acell->bound_14_dev[ib * 2]; ia2 = acell->bound_14_dev[ib * 2 + 1];
				i3a1 = x_idx(ia1); i3a2 = x_idx(ia2);
				dr[0] = acell->r_dev[i3a1] - acell->r_dev[i3a2];
				dr[1] = acell->r_dev[i3a1 + 1] - acell->r_dev[i3a2 + 1];
				dr[2] = acell->r_dev[i3a1 + 2] - acell->r_dev[i3a2 + 2];
				if (dr[0] > Rmax_boundLength) dr[0] -= xd[0];
				else if (dr[0] < -Rmax_boundLength) dr[0] += xd[0];
				if (dr[1] > Rmax_boundLength) dr[1] -= xd[1];
				else if (dr[1] < -Rmax_boundLength) dr[1] += xd[1];
				if (dr[2] > Rmax_boundLength) dr[2] -= xd[2];
				else if (dr[2] < -Rmax_boundLength) dr[2] += xd[2];

				R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]; if (R2 < 0.2) R2 = 0.2;
				R = sqrt(R2);
			
				if (R > rcut) continue;
				u[0] = dr[0] / R; u[1] = dr[1] / R; u[2] = dr[2] / R;

				// calculation of the interaction
				//memset(force, 0, SIZE_V3);
				if ((acell->sr_type_dev[ia1] & 0x01) && (acell->sr_type_dev[ia2] & 0x01)) { // L-J interaction
					eps = sqrt(acell->sr_par_dev[isr1] * acell->sr_par_dev[isr2]);
					rLJ = (acell->sr_par_dev[isr1 + 1] + acell->sr_par_dev[isr2 + 1]) * 0.5;
								
					if (iLJ == 0) standard_LJ(eps, rLJ, R, U, f);
					else alternative_LJ(eps, rLJ, R, U, f);

					Ut += U * 0.5;
					force[0] = f * u[0]; force[1] = f * u[1]; force[2] = f * u[2];
					/*  so far we do not enable it
					if ((acell->sr_type[ia1] & 0x02) && (acell->sr_type[ia2] & 0x02)) { // given profile -- bSR
						aindx[0] = acell->aID_dev[ia1]; aindx[1] = acell->aID_dev.m[ia2];
						key = construct_char2_key_pairwise(aindx);
						pw = sr_db.search(key);
						if (pw != NULL) {
							interpolate2(t1, t2, pw->pv, R, i, t3)
							U_LJ += t1 * 0.5;
							if (t2 > force_max && scale_force) t2 = force_max; 
							force[0] += t2 * u[0]; force[1] += t2 * u[1]; force[2] += t2 * u[2];
						}
					}
					*/
					acell->Fsr_14_dev[i3b] = force[0];
					acell->Fsr_14_dev[i3b + 1] = force[1];
					acell->Fsr_14_dev[i3b + 2] = force[2];
					//memcpy(acell->Fsr_14_dev + i3b, force, SIZE_V3);
					acell->vsr_14_dev[i3b] = force[0] * dr[0];
					acell->vsr_14_dev[i3b + 1] = force[1] * dr[1];
					acell->vsr_14_dev[i3b + 2] = force[2] * dr[2];
					acell->U14_sr_dev[ib] = U;
				}
			}
		}
	}

	__global__ void collect_LJ_res(GPU_SR_AtCELL *acell) {
		// works on 1 dimensional grid and one dimensional blocks
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int Na = acell->Na;
		int nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;

		int i, j, ia1, ia2, ia, ib, iba1, iba2;
		int i3a1, i3a2;

		int Nb, Nm;

		if (acell->Na_b > 0) {
			Nb = acell->Na_b;
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;

			for (i = 0; i < nEachThread; i++) {
				ia = tid * nEachThread + i;
				if (ia >= Nb) break;
				Nm = ia * (acell->Nm_b + 2);
				ia1 = acell->b_relship_dev[Nm];
				i3a1 = x_idx(ia1);
				for (j = 0; j < acell->b_relship_dev[Nm + 1]; j++) {
					ib = acell->b_relship_dev[Nm + 2 + j];
					i3a2 = x_idx(ib);
					iba1 = acell->bound_dev[2 * ib];
					iba2 = acell->bound_dev[2 * ib + 1];
					if (iba1 == ia1) { // int this bound, ia1 is the first atom
						acell->Fsr_dev[i3a1] -= acell->Fsr_b_dev[i3a2];
						acell->Fsr_dev[i3a1 + 1] -= acell->Fsr_b_dev[i3a2 + 1];
						acell->Fsr_dev[i3a1 + 2] -= acell->Fsr_b_dev[i3a2 + 2];

						// collecting energy and virial term at this case only
						acell->U_sr_dev[ia1] -= acell->Ub_sr_dev[ib];
						acell->vsr_dev[i3a1] -= acell->vsr_b_dev[i3a2];
						acell->vsr_dev[i3a1 + 1] -= acell->vsr_b_dev[i3a2 + 1];
						acell->vsr_dev[i3a1 + 2] -= acell->vsr_b_dev[i3a2 + 2];
					}
					else if (iba2 == ia1) { // ia1 is the second atom of the bound
						acell->Fsr_dev[i3a1] += acell->Fsr_b_dev[i3a2];
						acell->Fsr_dev[i3a1 + 1] += acell->Fsr_b_dev[i3a2 + 1];
						acell->Fsr_dev[i3a1 + 2] += acell->Fsr_b_dev[i3a2 + 2];
					}
				}
			}
		}
		//__syncthreads();

		if (acell->Na_b14 > 0) {
			Nb = acell->Na_b14;
			nEachThread = Nb / (gridDim.x * blockDim.x);
			if (nEachThread == 0) nEachThread = 1;
			if ((nEachThread * gridDim.x * blockDim.x) < Nb) nEachThread += 1;

			for (i = 0; i < nEachThread; i++) {
				ia = tid * nEachThread + i;
				if (ia >= Nb) break;
				Nm = ia * (acell->Nm_b14 + 2);
				ia1 = acell->b14_relship_dev[Nm];
				i3a1 = x_idx(ia1);
				for (j = 0; j < acell->b14_relship_dev[Nm + 1]; j++) {
					ib = acell->b14_relship_dev[Nm + 2 + j];
					i3a2 = x_idx(ib);
					iba1 = acell->bound_14_dev[2 * ib];
					iba2 = acell->bound_14_dev[2 * ib + 1];
					if (iba1 == ia1) { // int this bound, ia1 is the first atom
						acell->Fsr_dev[i3a1] -= acell->Fsr_14_dev[i3a2] * acell->c14;
						acell->Fsr_dev[i3a1 + 1] -= acell->Fsr_14_dev[i3a2 + 1] * acell->c14;
						acell->Fsr_dev[i3a1 + 2] -= acell->Fsr_14_dev[i3a2 + 2] * acell->c14;

						// collecting energy and virial term at this case only
						acell->U_sr_dev[ia1] -= acell->U14_sr_dev[ib] * acell->c14;
						acell->vsr_dev[i3a1] -= acell->vsr_14_dev[i3a2] * acell->c14;
						acell->vsr_dev[i3a1 + 1] -= acell->vsr_14_dev[i3a2 + 1] * acell->c14;
						acell->vsr_dev[i3a1 + 2] -= acell->vsr_14_dev[i3a2 + 2] * acell->c14;
					}
					else if (iba2 == ia1) { // ia1 is the second atom of the bound
						acell->Fsr_dev[i3a1] += acell->Fsr_14_dev[i3a2] * acell->c14;
						acell->Fsr_dev[i3a1 + 1] += acell->Fsr_14_dev[i3a2 + 1] * acell->c14;
						acell->Fsr_dev[i3a1 + 2] += acell->Fsr_14_dev[i3a2 + 2] * acell->c14;
					}
				}
			}
		}
		//__syncthreads();

		Na = acell->Na;
		nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;

		for (i = 0; i < nEachThread; i++) {
			ia1 = tid * nEachThread + i; i3a1 = x_idx(ia1);
			if (ia1 >= Na) break;
			acell->bShocked_map[ia1] = acell->bShocked_dev[ia1];
			//memcpy(acell->Fsr_map + i3a1, acell->Fsr_dev + i3a1, SIZE_V3);
			acell->Fsr_map[i3a1] = acell->Fsr_dev[i3a1];
			acell->Fsr_map[i3a1 + 1] = acell->Fsr_dev[i3a1 + 1];
			acell->Fsr_map[i3a1 + 2] = acell->Fsr_dev[i3a1 + 2];
			//memcpy(acell->vsr_map + i3a1, acell->vsr_dev + i3a1, SIZE_V3);
			acell->vsr_map[i3a1] = acell->vsr_dev[i3a1];
			acell->vsr_map[i3a1 + 1] = acell->vsr_dev[i3a1 + 1];
			acell->vsr_map[i3a1 + 2] = acell->vsr_dev[i3a1 + 2];
			acell->U_sr_map[ia1] = acell->U_sr_dev[ia1];
		}
	}
	
	__host__ void GPU_LJ_interact(GPU_SR_AtCELL *acell_host, GPU_SR_AtCELL *acell_dev) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		
		int dimBlock = 8;
		int dimGrid = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(acell_host->Na, dimBlock, dimGrid);
		LJ_interact<<< dimGrid, dimBlock, 0, cuStream >>>(acell_dev);
		cudaStreamSynchronize(cuStream);
		
		collect_LJ_res<<< dimGrid, dimBlock, 0, cuStream >>>(acell_dev);
		cudaStreamSynchronize(cuStream);
		
		cudaStreamDestroy(cuStream);
	}

} // end of namespace _gpu_md_



namespace _gpu_md_ {
	__global__ void prior_relocate_atom_kernel(GPU_AtCELL *acell) {
		// work on one cubic of the meshed cell, gridDim is [Nx, Ny] while blockDim is Nz
		int icx = blockIdx.x;
		int icy = blockIdx.y;
		int icz = threadIdx.x;

		int *cubic = acell->device_cubic(icx, icy, icz);
		int i;
		for (i = 0; i <= acell->NC; i++) cubic[i] = 0;
	};
	__global__ void update_location_at_host_kernel(GPU_AtCELL *acell) {
		// work on one cubic of the meshed cell, gridDim is [Nx, Ny] while blockDim is Nz
		int icx = blockIdx.x;
		int icy = blockIdx.y;
		int icz = threadIdx.x;

		int *cubic = acell->device_cubic(icx, icy, icz);
		int *cubic_host = acell->host_cubic(icx, icy, icz);
		int i;
		for (i = 0; i <= acell->NC; i++) cubic_host[i] = cubic[i];
	};
	__global__ void relocate_atom_kernel(GPU_AtCELL *acell) {
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(acell->Na, gridDim.x, blockDim.x);
		
		double *r_map = acell->r_map, *r_dev = acell->r_dev;
		int *ir_map = acell->ir_map, *ir_dev = acell->ir_dev;
		int *cubic;
		double du;
		
		double xl[3] = {acell->cx[0], acell->cx[1], acell->cx[2]};
		double md[3] = {acell->cw[0], acell->cw[1], acell->cw[2]};
		int Nx = acell->Nx, Ny = acell->Ny, Nz = acell->Nz, Na = acell->Na;
		//double xd[3] = {md[0] * Nx, md[1] * Ny, md[2] * Nz};

		int i, iatom, i3;
		int ix, iy, iz, nc, NC = acell->NC;
		for (i = 0; i < nEachThread; i++) {
			iatom = tid * nEachThread + i;
			if (iatom >= Na) break;
			i3 = iatom * 3;
			r_dev[i3] = r_map[i3]; r_dev[i3 + 1] = r_map[i3 + 1]; r_dev[i3 + 2] = r_map[i3 + 2];
			du = (r_dev[i3] - xl[0]) / md[0]; ix = int(du);
			du = (r_dev[i3 + 1] - xl[1]) / md[1]; iy = int(du);
			du = (r_dev[i3 + 2] - xl[2]) / md[2]; iz = int(du);
			
			// if atom is out of the whole cell, we do not want to locate it to another side of the box
			// because short-range interaction will need the atoms in one cluster has same coordination !!
			if (ix < 0) ix = 0; else if (ix >= Nx) ix = Nx - 1;
			if (iy < 0) iy = 0; else if (iy >= Ny) iy = Ny - 1;
			if (iz < 0) iz = 0; else if (iz >= Nz) iz = Nz - 1;
			
			ir_dev[i3] = ix; ir_dev[i3 + 1] = iy; ir_dev[i3 + 2] = iz;
			ir_map[i3] = ix; ir_map[i3 + 1] = iy; ir_map[i3 + 2] = iz;

			cubic = acell->device_cubic(ix, iy, iz);
			nc = atomicAdd(cubic, 1); // do this first, else we can not get unit index
			if (nc < (NC-1)) {
				atomicExch(&(cubic[nc + 1]), iatom);
			}
			else {
				atomicSub(cubic, 1); //-- reached limit, ignore this atom
				atomicAdd(acell->ivar_map, 1);
			}
		}
		__syncthreads(); // not really necessary
	};

	__host__ bool relocate_atom(GPU_AtCELL *acell_host, GPU_AtCELL *acell_dev, bool bUpdateHost) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		
		acell_host->ivar_map[0] = 0;
		
		dim3 dim3Grid(acell_host->Nx, acell_host->Ny);
		int dimBlock = acell_host->Nz;
		prior_relocate_atom_kernel<<< dim3Grid, dimBlock, 0, cuStream >>>(acell_dev);
		
		dimBlock = 128;
		int dimGrid = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(acell_host->Na, dimBlock, dimGrid);
		relocate_atom_kernel<<< dimGrid, dimBlock, 0, cuStream >>>(acell_dev);
		
		dimBlock = acell_host->Nz;
		bUpdateHost = true;
		if (bUpdateHost) update_location_at_host_kernel<<< dim3Grid, dimBlock, 0, cuStream >>>(acell_dev);
		
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		
		if (acell_host->ivar_map[0] > 0) return false;
		else return true;
	};
	
	__global__ void check_neighbors_kernel(GPU_EAtCELL *acell) {
		// works on 1 dimensional grid and one dimensional blocks
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int Na = acell->Na;
		int nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;

		if (acell->srnb_dev == NULL) return;

		int i, ia1, ia, ia2, i3a1, i3a2;
		int icx, icy, icz, ncx, ncy, ncz, dcx, dcy, dcz;
		int Nc;

		DOUBLE R, R2;
		DOUBLE xd[3] = {acell->cw[0] * acell->Nx, acell->cw[1] * acell->Ny, acell->cw[2] * acell->Nz};
		DOUBLE doff[3];

		DOUBLE dr[3];//, u[3];

		int hx = acell->Nx / 2, hy = acell->Ny, hz = acell->Nz / 2;
		int *cubic;

		DOUBLE rcut = acell->rcut, r2cut = rcut * rcut;
		DOUBLE cw_max = acell->cw[0];
		if (acell->cw[1] > cw_max) cw_max = acell->cw[1];
		if (acell->cw[2] > cw_max) cw_max = acell->cw[2];
		DOUBLE rcut_m = rcut + cw_max * 2, r2cut_m = rcut_m * rcut_m;
		DOUBLE d2x, d2y, d2z, d2xy, d2;

		dcx = int(rcut / acell->cw[0]) + 1;//2;
		dcy = int(rcut / acell->cw[1]) + 1;
		dcz = int(rcut / acell->cw[2]) + 1;

		int drIdx, srIdx, nsr, n_srnb = acell->n_srnb;
		int *isr = NULL;
		DOUBLE *drsr = NULL;

		for (i = 0; i < nEachThread; i++) {
			ia1 = tid * nEachThread + i; i3a1 = x_idx(ia1);
			if (ia1 >= Na) break;

			srIdx = acell->srIndx(ia1);
			drIdx = acell->drIndx(ia1);
			isr = acell->srnb_dev + srIdx;
			drsr = acell->dr_dev + drIdx;

			for (icx = -dcx; icx <= dcx; icx++) {
				ncx = acell->ir_dev[i3a1] + icx;
				if (ncx < 0) {ncx += acell->Nx; doff[0] = -xd[0];}
				else if (ncx >= acell->Nx) {ncx -= acell->Nx; doff[0] = xd[0];}
				else doff[0] = 0;

				d2x = icx * acell->cw[0]; d2x *= d2x;

				for (icy = -dcy; icy <= dcy; icy++) {
					ncy = acell->ir_dev[i3a1 + 1] + icy;
					if (ncy < 0) {ncy += acell->Ny; doff[1] = -xd[1];}
					else if (ncy >= acell->Ny) {ncy -= acell->Ny; doff[1] = xd[1];}
					else doff[1] = 0;

					d2y = icy * acell->cw[1]; d2y *= d2y;
					d2xy = d2x + d2y;
					if (d2xy > r2cut_m) continue;

					for (icz = -dcz; icz <= dcz; icz++) {
						ncz = acell->ir_dev[i3a1 + 2] + icz;
						if (ncz < 0) {ncz += acell->Nz; doff[2] = -xd[2];}
						else if (ncz >= acell->Nz) {ncz -= acell->Nz; doff[2] = xd[2];}
						else doff[2] = 0;

						d2z = icz * acell->cw[2]; d2z *= d2z;
						d2 = d2xy + d2z;
						if (d2 > r2cut_m) continue;

						cubic = acell->device_cubic(ncx, ncy, ncz);
						Nc = cubic[0];
						for (ia = 0; ia < Nc; ia++) {
							ia2 = cubic[ia + 1];
							if (ia2 >= Na) continue;
							if (ia2 == ia1) continue;

							i3a2 = x_idx(ia2);
							dr[0] = acell->r_dev[i3a1] - acell->r_dev[i3a2] - doff[0];
							dr[1] = acell->r_dev[i3a1 + 1] - acell->r_dev[i3a2 + 1] - doff[1];
							dr[2] = acell->r_dev[i3a1 + 2] - acell->r_dev[i3a2 + 2] - doff[2];
							R2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
							if (R2 > r2cut) continue;
							R = GPUSQRT(R2);

							nsr = atomicAdd(isr, 1);
							if (nsr < n_srnb) {
								isr[nsr + 1] = ia2;
								drsr[nsr * 4] = -dr[0]; drsr[nsr * 4 + 1] = -dr[1]; drsr[nsr * 4 + 2] = -dr[2]; drsr[nsr * 4 + 3] = R;
								//atomicAdd(acell->ivar_map + 1, 1);
							}
							else {
								atomicSub(isr, 1); // overflowed
								atomicAdd(acell->ivar_map, 1);
							}

						}
						/*
						atomicAdd(acell->ivar_map + 2, Nc);
						if (ia1 == 1) {
							acell->ivar_map[5] = acell->ir_dev[i3a1];
							acell->ivar_map[6] = acell->ir_dev[i3a1 + 1];
							acell->ivar_map[7] = acell->ir_dev[i3a1 + 2];
						}
						*/
					}
				}
			}
		}
	}
	
	__global__ void update_neighbors_status_kernel(GPU_EAtCELL *acell_dev, bool status) {
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		if (tid == 0) {acell_dev->bSRInitialized = status;}
	}

	__host__ void check_neighbors(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);

		int dimBlock = 128;
		int dimGrid = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(acell_host->Na, dimBlock, dimGrid);
		check_neighbors_kernel<<< dimGrid, dimBlock, 0, cuStream >>>(acell_dev);
		update_neighbors_status_kernel<<< 1, 32, 0, cuStream >>>(acell_dev, true);

		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	}

	__host__ void update_neighbors_status(GPU_EAtCELL *acell_dev, bool status) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);

		update_neighbors_status_kernel<<< 1, 32, 0, cuStream >>>(acell_dev, status);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	}

	__global__ void update_charge_dipole_kernel(GPU_EAtCELL *acell) {
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(acell->Na, gridDim.x, blockDim.x);
		
		double *mu_map = acell->mu_map, *mu_dev = acell->mu_dev;
		double *c_map = acell->c_map, *c_dev = acell->c_dev;
		
		int i, iatom, i3;
		int Na = acell->Na;
		for (i = 0; i < nEachThread; i++) {
			iatom = tid * nEachThread + i;
			if (iatom >= Na) break;
			i3 = iatom * 3;
			mu_dev[i3] = mu_map[i3]; mu_dev[i3 + 1] = mu_map[i3 + 1]; mu_dev[i3 + 2] = mu_map[i3 + 2];
			c_dev[iatom] = c_map[iatom];
		}
	};
	
	__host__ void update_charge_dipole(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		
		int dimBlock = 128;
		int dimGrid = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(acell_host->Na, dimBlock, dimGrid);
		update_charge_dipole_kernel<<< dimGrid, dimBlock, 0, cuStream >>>(acell_dev);

		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	};
	
	__global__ void reset_efield_kernel(GPU_EAtCELL *acell) {
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(acell->Na, gridDim.x, blockDim.x);
		
		double *E_map = acell->E_map, *E_dev = acell->E_dev;
		
		int i, iatom, i3;
		int Na = acell->Na;
		for (i = 0; i < nEachThread; i++) {
			iatom = tid * nEachThread + i;
			if (iatom >= Na) break;
			i3 = iatom * 3;
			E_dev[i3] = 0; E_dev[i3 + 1] = 0; E_dev[i3 + 2] = 0;
			E_map[i3] = 0; E_map[i3 + 1] = 0; E_map[i3 + 2] = 0;
		}
	};
	__host__ void reset_efield(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		
		int dimBlock = 128;
		int dimGrid = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(acell_host->Na, dimBlock, dimGrid);
		reset_efield_kernel<<< dimGrid, dimBlock, 0, cuStream >>>(acell_dev);

		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	};
}

#endif // __GPU__
