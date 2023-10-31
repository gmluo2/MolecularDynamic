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
#include "gpu_esr.h"

#include "gpu_sr.h"


namespace _gpu_md_ {
	namespace _gpu_ewald_ {
		// following is copied from atom_sr_interact.cpp
		// 2d Ewald-Sum, local part only

		__device__ void EwaldSumRecpKzero_E_Q(EwaldSumRecpK0Var &var, DOUBLE q2, DOUBLE &Ez) {
			Ez = -q2 * var.ferf * var.PI2_over_A;
		}

		__device__ void EwaldSumRecpKzero_E_MU(EwaldSumRecpK0Var &var, DOUBLE q2, DOUBLE mu2z, DOUBLE &Ez) {
			Ez = -(q2 * var.ferf + mu2z * var.KAPPA2_over_SQRT_PI * var.fexp) * var.PI2_over_A;
		}


		__device__ void EwaldSumRecpKzero_Q(EwaldSumRecpK0Var &var, DOUBLE q1, DOUBLE q2, /*DOUBLE &Ez, */DOUBLE &Fz, DOUBLE &U) {
			//Ez = -q2 * var.ferf * var.PI2_over_A;
			//if (var.bEfieldOnly) return;

			DOUBLE q12 = q1 * q2;
			U = -q12 * var.W * var.PI_over_A; // half of interaction
			Fz = -q12 * var.ferf * var.PI2_over_A;
		}

		__device__ void EwaldSumRecpKzero_MU(EwaldSumRecpK0Var &var, DOUBLE q1, DOUBLE mu1z, DOUBLE q2, DOUBLE mu2z, /*DOUBLE &Ez, */DOUBLE &Fz, DOUBLE &U) {
			//Ez = -(q2 * var.ferf + mu2z * var.KAPPA2_over_SQRT_PI * var.fexp) * var.PI2_over_A;
			//Ez = -q2 * var.ferf;
			//if (var.fexp != 0) {
			//	Ez -= mu2z * var.KAPPA2_over_SQRT_PI * var.fexp;
			//}
			//Ez *= var.PI2_over_A;
			//if (var.bEfieldOnly) return;

			DOUBLE q12 = q1 * q2;
			if (var.iExIndDipole == 1) {
				U = -(q12 * var.W + (q1 * mu2z - mu1z * q2) * var.ferf/* - mu1z * mu2z * var.KAPPA2_over_SQRT_PI * var.fexp*/) * var.PI_over_A; // half of interaction
				Fz = -(q12 * var.ferf + (q1 * mu2z - mu1z * q2) * var.KAPPA2_over_SQRT_PI * var.fexp /*+ mu1z * mu2z * var.KAPPA2_over_SQRT_PI * 2 * var.kappa2 * var.dz * var.fexp*/) * var.PI2_over_A;
			}
			else {
				U = -(q12 * var.W + (q1 * mu2z - mu1z * q2) * var.ferf - mu1z * mu2z * var.KAPPA2_over_SQRT_PI * var.fexp) * var.PI_over_A; // half of interaction
				Fz = -(q12 * var.ferf + (q1 * mu2z - mu1z * q2) * var.KAPPA2_over_SQRT_PI * var.fexp + mu1z * mu2z * var.KAPPA2_over_SQRT_PI * 2 * var.kappa2 * var.dz * var.fexp) * var.PI2_over_A;
			}
		}

	} // namespace _gpu_ewald_

	// for 2d system, long-range interaction in Ewald Sum has non-zero term for k = 0
	// the force, efield and other results are accumulated to the device memory
	__global__ void EwaldSum_2D_K0(GPU_EAtCELL *cell) {
		// works on 1 dimensional grid and one dimensional blocks
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int Na = cell->Na;
		int nEachThread = Na / (gridDim.x * blockDim.x);
		if (nEachThread == 0) nEachThread = 1;
		if ((nEachThread * gridDim.x * blockDim.x) < Na) nEachThread += 1;


		_gpu_ewald_::EwaldSumRecpK0Var var;
		var.init(((cell->bCalE == 1 || cell->bCalE == 11) ? true: false), cell->esv.iExIndDipole, cell->esv.kappa, cell->esv.A);

		DOUBLE q1, q2, muz1, muz2;
		DOUBLE Ez, Fz, Ezt, Fzt;
		DOUBLE dz, z1;
		int i, ia1, ia2, i3az1, i3az2, i3a1;

		DOUBLE t1;

		DOUBLE Uestat = 0, STzz, STxx, STyy;

		int cal_mode = cell->bCalE;
		if (cal_mode <= 0 || cal_mode >= 3) return; // with charge only, so far choice 1 -- electric field; 2 -- force and energy
		if (cell->bDipole) cal_mode += 10; // with dipole


		for (i = 0; i < nEachThread; i++) {
			ia1 = tid * nEachThread + i; i3a1 = x_idx(ia1);
			if (ia1 >= Na) break;
			//if (eNeutral) continue;
			i3az1 = z_idx(ia1);
			q1 = cell->c_dev[ia1]; muz1 = cell->mu_dev[i3az1]; z1 = cell->r_dev[i3az1];

			Uestat = 0; Ez = 0; Fz = 0; STxx = 0; STyy = 0; STzz = 0;

			for (ia2 = 0; ia2 < cell->Na; ia2++) {
				//if (patom1->eNeutral) continue; // ignore this atom for charge interaction
				i3az2 = z_idx(ia2);
				dz = cell->r_dev[i3az2] - z1;
				q2 = cell->c_dev[ia2]; muz2 = cell->mu_dev[i3az2];

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
				cell->E_dev[i3az1] += Ez;
				break;
			case 2: // force and energy with q
				cell->U_estat_dev[ia1] += Uestat;
				cell->Fe_dev[i3az1] += Fz;
				cell->ve_dev[i3a1] += STxx;
				cell->ve_dev[i3a1 + 1] += STyy;
				cell->ve_dev[i3az1] += STzz;
				break;
			case 11: // electric field only with q and mu
				cell->E_dev[i3az1] += Ez;
				break;
			case 12: // force and energy with q and mu
				cell->U_estat_dev[ia1] += Uestat;
				cell->Fe_dev[i3az1] += Fz;
				cell->ve_dev[i3a1] += STxx;
				cell->ve_dev[i3a1 + 1] += STyy;
				cell->ve_dev[i3az1] += STzz;
				break;
			default:
				break;
			}
		}
	}

} // end of namespace _gpu_md_

#endif // __GPU__
