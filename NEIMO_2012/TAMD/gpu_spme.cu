#include "project.h"

//#include <iostream>
//#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

//#include "ranlib.h"

#include "def.h"
//#include "vector.h"
#include "gpu_vector.h"

//#include "var.h"

//#include "complex.h"
//#include <fftw3.h>


#include "gpu_spme.h"

namespace _gpu_md_ {
namespace _spme_ {
#ifdef __GPU__

	__global__ void FT_prior(VSPME *vspme) { // 1d block, Q ==> ftd
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		if (tid == 0) {
			memcpy(vspme->ftd.m, vspme->Q.mb.m, vspme->ftd.n * sizeof(DOUBLE));
		}
	};
	__global__ void FT_prior_2d(VSPME *vspme) { // Q ==> ftd
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(vspme->Nv, gridDim.x, blockDim.x);
		int ith, i, j, k, n, tn;
		DOUBLE vb, vc, vbc;
		//COMPLEX *fq = NULL;

		for (ith = 0; ith < nEachThread; ith++) {
			n = nEachThread * tid + ith;
			if (n >= vspme->Nv) break;
			i = n / vspme->Nyz; tn = n - i * vspme->Nyz;
			j = tn / vspme->Nz; k = tn - j * vspme->Nz;
		
			if (vspme->bUseHostFFT == 0 || vspme->ftd_host == NULL) vspme->ftd.m[n] = vspme->Q.mb.m[n];
			else vspme->ftd_host[n] = vspme->Q.mb.m[n];
		}
	};
	
	__global__ void FT_post1_2d(VSPME *vspme) { // 2d -- block dimension, ftc ==> FQ
		// FT have half dimension only
		int Nz = vspme->nhz + 1, Ny = vspme->Ny, Nx = vspme->Nx;
		int Nyz = Ny * Nz, Nv = Nx * Nyz;

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(Nv, gridDim.x, blockDim.x);
		int ith, i, j, k, n, m, tn;
		
		DOUBLE vb, vc, vbc;
		COMPLEX *fq = NULL, *ftc;
		cuDoubleComplex *dfq = NULL, *dftc;
		
		int inx, iny, inz;
		for (ith = 0; ith < nEachThread; ith++) {
			n = nEachThread * tid + ith;
			if (n >= Nv) break;
			i = n / Nyz; tn = n - i * Nyz;
			j = tn / Nz; k = tn - j * Nz;

			//inx = (i == 0 ? 0 : vspme->Nx - i); // the inverse index
			//iny = (j == 0 ? 0 : (vspme->Ny - j));
			inx = (i == 0 ? 0 : (Nx - i)); // the inverse index
			iny = (j == 0 ? 0 : (Ny - j));
			inz = (k == 0 ? 0 : (vspme->Nz - k));
			
				/*
				for (k = 0; k <= _VSPME<atom>::nhz; k++) {
					FQ.m[i].m[j].m[k].x = _VSPME<atom>::ftc.m[n][0];
					FQ.m[i].m[j].m[k].y = _VSPME<atom>::ftc.m[n][1];
					n++;
				}
				*/
				
			//n = (i * vspme->Ny + j) * (vspme->nhz + 1);
			//m = i * vspme->Nyz + j * vspme->Nz;
			n = (i * Ny + j) * Nz + k;
			m = (i * Ny + j) * vspme->Nz + k;
			if (vspme->bUseHostFFT == 0 || vspme->ftc_host == NULL) {
				fq = vspme->FQ.mb.m + m; ftc = vspme->ftc.m + n;
				fq->x = ftc->x; fq->y = ftc->y;
				if (k != 0) {
					vspme->FQ.e(inx, iny, inz)->x = ftc->x;
					vspme->FQ.e(inx, iny, inz)->y = -ftc->y;
				}
			}
			else {
				fq = vspme->FQ.mb.m + m; dftc = vspme->ftc_host + n;
				fq->x = DOUBLE(dftc->x); fq->y = DOUBLE(dftc->y);
				if (k != 0) {
					vspme->FQ.e(inx, iny, inz)->x = DOUBLE(dftc->x);
					vspme->FQ.e(inx, iny, inz)->y = -DOUBLE(dftc->y);
				}
			}
			
			// FQ ==> FBCQ
			vc = *(vspme->C.e(i, j, k)); //vc = vspme->C.m[i].m[j].m[k];
			vb = vspme->b[0].m[i] * vspme->b[1].m[j] * vspme->b[2].m[k]; vbc = vc * vb;
			fq = vspme->FQ.e(i, j, k);
			vspme->FBCQ.e(i, j, k)->x = vbc * fq->x; vspme->FBCQ.e(i, j, k)->y = vbc * fq->y;
				
			if (k != 0) {
				vc = *(vspme->C.e(inx, iny, inz)); //vc = vspme->C.m[i].m[j].m[k]; 
				vb = vspme->b[0].m[inx] * vspme->b[1].m[iny] * vspme->b[2].m[inz]; vbc = vc * vb;
				fq = vspme->FQ.e(inx, iny, inz);
				vspme->FBCQ.e(inx, iny, inz)->x = vbc * fq->x; vspme->FBCQ.e(inx, iny, inz)->y = vbc * fq->y;
			}
			
			// FBCQ ==>ftc
			if (vspme->bUseHostFFT == 0 || vspme->ftc_host == NULL) {
				ftc = vspme->ftc.m + n; fq = vspme->FBCQ.mb.m + m;
				ftc->x = fq->x; ftc->y = fq->y;
			}
			else {
				dftc = vspme->ftc_host + n; fq = vspme->FBCQ.mb.m + m;
				dftc->x = double(fq->x); dftc->y = double(fq->y);
			}
			//n += _VSPME<atom>::nhz + 1;
		}
	};

	__global__ void IFT_post_2d(VSPME *vspme) { // ftd ==> BCQ
		int Nv = vspme->Nv;
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(vspme->Nv, gridDim.x, blockDim.x);
		int ith, n;

		for (ith = 0; ith < nEachThread; ith++) {
			n = nEachThread * tid + ith;
			if (n >= Nv) break;

			if (vspme->bUseHostFFT == 0 || vspme->ftd_host == NULL) {
				//memcpy(vspme->BCQ.mb.m + n, vspme->ftd.m + n, vspme->Nz * sizeof(double));
				vspme->BCQ.mb.m[n] = vspme->ftd.m[n];
			}
			else {
				//memcpy(vspme->BCQ.mb.m + n, vspme->ftd_host + n, vspme->Nz * sizeof(double));
				vspme->BCQ.mb.m[n] = DOUBLE(vspme->ftd_host[n]);
			}
		}
	};

/*
	__global__ void IFT_post(VSPME *vspme) { // ftd ==> BCQ
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		if (tid == 0) {
			memcpy(vspme->BCQ.mb.m, vspme->ftd.m, vspme->ftd.n * sizeof(DOUBLE));
		}
	};
*/
	__host__ void vspme_FT_IFT(VSPME *vspme_host, VSPME *vspme_dev) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		
		//cufftPlan3d(&vspme_host->ft_plan, vspme_host->Nx, vspme_host->Ny, vspme_host->Nz, CUFFT_R2C);
		//cufftPlan3d(&vspme_host->ift_plan, vspme_host->Nx, vspme_host->Ny, vspme_host->Nz, CUFFT_C2R);
		
		int blockDim = 256;
		int gridDim_hv = 0, gridDim = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(vspme_host->Nv, blockDim, gridDim);
		_gpu_dstruct_::_gpu_util_::gpu_grid(vspme_host->Nx * vspme_host->Ny * (vspme_host->nhz + 1), blockDim, gridDim_hv);

		cufftSetStream(vspme_host->ft_plan, cuStream);
		cufftSetStream(vspme_host->ift_plan, cuStream);

		//FT_prior<<< 1, 1, 0, cuStream >>>(vspme_dev);
		FT_prior_2d<<< gridDim, blockDim, 0, cuStream >>>(vspme_dev);
		//cudaStreamSynchronize(cuStream);
		// cuFFT is host-API library
		CUFFT_ExecD2Z(vspme_host->ft_plan, vspme_host->ftd.m, vspme_host->ftc.m); // ftd ==FFT==> ftc, double ==> cuDoubleComplex
		//cudaStreamSynchronize(cuStream);

		// the job assigned with half of Nz, or full dimension : Nx * Ny * (nhz + 1)
		FT_post1_2d<<< gridDim_hv, blockDim, 0, cuStream >>>(vspme_dev);
		//cudaStreamSynchronize(cuStream);
		
		// cuFFT is host-API library
		CUFFT_ExecZ2D(vspme_host->ift_plan, vspme_host->ftc.m, vspme_host->ftd.m);
		//cudaStreamSynchronize(cuStream);
		//IFT_post<<< 1, 1, 0, cuStream >>>(vspme_dev);
		IFT_post_2d<<< gridDim, blockDim, 0, cuStream >>>(vspme_dev);
		
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	};
	
	
	
	// use CPU-FFTW
	__host__ void vspme_FT_prior(VSPME *vspme_host, VSPME *vspme_dev) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		
		int blockDim = 256;
		int gridDim = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(vspme_host->Nv, blockDim, gridDim);

		//FT_prior<<< 1, 1, 0, cuStream >>>(vspme_dev);
		FT_prior_2d<<< gridDim, blockDim, 0, cuStream >>>(vspme_dev);

		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	};
	
	__host__ void vspme_FT_post(VSPME *vspme_host, VSPME *vspme_dev) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		
		
		int blockDim = 256;
		int gridDim = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(vspme_host->Nx * vspme_host->Ny * (vspme_host->nhz + 1), blockDim, gridDim);

		// the job assigned with half of Nz, or full dimension : Nx * Ny * (nhz + 1)
		FT_post1_2d<<< gridDim, blockDim, 0, cuStream >>>(vspme_dev);
		
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	};
	
	__host__ void vspme_IFT_post(VSPME *vspme_host, VSPME *vspme_dev) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		
		int blockDim = 256;
		int gridDim = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(vspme_host->Nv, blockDim, gridDim);

		//IFT_post<<< 1, 1, 0, cuStream >>>(vspme_dev);
		IFT_post_2d<<< gridDim, blockDim, 0, cuStream >>>(vspme_dev);
		
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	};
	
	

	__host__ double cal_U_recp(VSPME *vspme_host, VSPME *vspme_dev) {
		//double E = _gpu_dstruct_::_gpu_util_::sum_array_cuml<double>(vspme_host->Q.mb.m, vspme_host->ftd.m, vspme_host->ftd.n);
		double E = 0;
		if (vspme_host->bUseHostFFT == 0 || vspme_host->ftd_host == NULL) {
			E = (double)_gpu_dstruct_::_gpu_util_::SUM_ARRAY2DD_CUML(vspme_host->Q.mb.m, vspme_host->ftd.m, vspme_host->ftd.n, vspme_host->dbuff_dev, vspme_host->dbuff_host, vspme_host->ndbuff);
		}
		else {
			E = (double)_gpu_dstruct_::_gpu_util_::SUM_ARRAY2Dd_CUML(vspme_host->Q.mb.m, vspme_host->ftd_host, vspme_host->ftd.n, vspme_host->dbuff_dev, vspme_host->dbuff_host, vspme_host->ndbuff);
		}
		E *= PI2 / vspme_host->V;
		return E;
	};


	__global__ void cal_ST_recp_0_kernel(VSPME* vspme, int iST, DOUBLE *s_block) {
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(vspme->Nv, gridDim.x, blockDim.x);

		__shared__ DOUBLE sum[512];
		int ith, i, j, k, n, tn;

		_gpu_dstruct_::FLEX_CMATRIX<DOUBLE> *c;
		switch(iST) {
		case 0:
			c = &(vspme->Ctr); break;
		case 1:
			c = &(vspme->Cxx); break;
		case 2:
			c = &(vspme->Cyy); break;
		case 3:
			c = &(vspme->Czz); break;
		default:
			c = &(vspme->Ctr); break;
		}

		DOUBLE vc, vb, vbc, sv = 0;
		COMPLEX *fq = NULL;
		//double ST = 0;
		for (ith = 0; ith < nEachThread; ith++) {
			n = nEachThread * tid + ith;
			if (n >= vspme->Nv) break;
			i = n / vspme->Nyz; tn = n - i * vspme->Nyz;
			j = tn / vspme->Nz; k = tn - j * vspme->Nz;

			vc = *(c->e(i, j, k)); //vc = c->m[i].m[j].m[k];
			if (vc == 0) continue;
			vb = vspme->b[0].m[i] * vspme->b[1].m[j] * vspme->b[2].m[k]; vbc = vc * vb;
			fq = vspme->FQ.e(i, j, k);
			//FBTQ.m[i].m[j].m[k].x = vbc * fq->x; FBTQ.m[i].m[j].m[k].y = vbc * fq->y;
			sv += vbc * (fq->x * fq->x + fq->y * fq->y);
		}
		sum[threadIdx.x] = sv;

		// sum up the value in the same block
		__syncthreads();
		if (threadIdx.x == 0) {// each block
			sv = 0;
			for (i = 0; i < blockDim.x; i++) sv += sum[i];
			s_block[blockIdx.x] = sv;
		}
		__syncthreads();
	};

	__host__ double cal_ST_recp_0(VSPME *vspme_host, VSPME* vspme_dev, int iST) {
		int blockDim = 256;
		int gridDim = vspme_host->Nv / blockDim;
		if (gridDim == 0) gridDim = 1;
		else if (gridDim > 1024) gridDim = 1024; // do not have too many blocks, but also can not be just few blocks
		if (gridDim > vspme_host->ndbuff) gridDim = vspme_host->ndbuff;

		cudaStream_t cuStream = NULL;
		cudaStreamCreate(&cuStream);
		cal_ST_recp_0_kernel<<< gridDim, blockDim, 512, cuStream >>>(vspme_dev, vspme_host->Nv, vspme_host->dbuff_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		DOUBLE *buff_host = vspme_host->dbuff_host;
		DOUBLE sum = 0;
		for (int i = 0; i < gridDim; i++) sum += buff_host[i];

		double ST = sum * PI2 / vspme_host->V;
		return ST;
	};

	__host__ void cal_ST_recp(VSPME *vspme_host, VSPME *vspme_dev, double& STxx, double& STyy, double& STzz) {
		STxx = cal_ST_recp_0(vspme_host, vspme_dev, 1); // xx
		STyy = cal_ST_recp_0(vspme_host, vspme_dev, 2); // yy
		STzz = cal_ST_recp_0(vspme_host, vspme_dev, 3); // zz
	}

	__host__ void cal_STtr_recp(VSPME *vspme_host, VSPME *vspme_dev, double& STtr) {
		STtr = cal_ST_recp_0(vspme_host, vspme_dev, 0); // trace
	}

	__global__ void prior_init_normalized_coordinates_kernel(_VSPME *vspme) {
		// work on one cubic of the meshed cell, gridDim is [Nx, Ny] while blockDim is Nz
		int icx = blockIdx.x;
		int icy = blockIdx.y;
		int icz = threadIdx.x;

		int *cubic = vspme->cmm_cubic(icx, icy, icz);
		int i;
		for (i = 0; i <= vspme->Nca; i++) cubic[i] = 0;
	};
	__global__ void init_normalized_coordinates_kernel(VSPME *vspme) {
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(vspme->Na, gridDim.x, blockDim.x);
		
		double *u = vspme->u.m, *r = vspme->r;
		int *cubic, *ic = vspme->ic.m;
		double *absp_M, *absp_dM, *absp_d2M, du;
		
		double xl[3] = {vspme->xl[0], vspme->xl[1], vspme->xl[2]};
		double md[3] = {vspme->md[0], vspme->md[1], vspme->md[2]};
		int Nx = vspme->Nx, Ny = vspme->Ny, Nz = vspme->Nz;

		int i, iatom, i3, j, ipt;
		int ix, iy, iz, nc, Nca = vspme->Nca;
		for (i = 0; i < nEachThread; i++) {
			iatom = tid * nEachThread + i;
			if (iatom >= vspme->Na) break;
			i3 = iatom * 3;
			u[i3] = (r[i3] - xl[0]) / md[0];
			u[i3 + 1] = (r[i3 + 1] - xl[1]) / md[1];
			u[i3 + 2] = (r[i3 + 2] - xl[2]) / md[2];

			ix = int(u[i3]); iy = int(u[i3+1]); iz = int(u[i3+2]);

			if (ix < 0) {ix += Nx; u[i3] += Nx;}
			else if (ix >= Nx) {ix -= Nx; u[i3] -= Nx;}

			if (iy < 0) {iy += Ny; u[i3 + 1] += Ny;}
			else if (iy >= Ny) {iy -= Ny; u[i3 + 1] -= Ny;}

			if (iz < 0) {iz += Nz; u[i3 + 2] += Nz;}
			else if (iz >= Nz) {iz -= Nz; u[i3 + 2] -= Nz;}

			ic[i3] = ix; ic[i3 + 1] = iy; ic[i3 + 2] = iz;
			cubic = vspme->cmm_cubic(ix, iy, iz);
			nc = atomicAdd(cubic, 1); // do this first, else we can not get unit index
			if (nc < (Nca-1)) {
				atomicExch(&(cubic[nc + 1]), iatom);
			}
			else {
				atomicSub(cubic, 1); //-- reached limit, ignore this atom
				atomicAdd(vspme->status_host, 1);
			}

			du = u[i3] - ix; //int(u[i3]);
			absp_M = vspme->absp_M(iatom, 0, 0); 
			absp_dM = vspme->absp_M(iatom, 1, 0);
			absp_d2M = vspme->absp_M(iatom, 2, 0);
			for (j = 0; j < vspme->N; j++) {
				ipt = vspme->bsp.get_pt(du);
				absp_M[j] = vspme->bsp.bsp_interpolate(ipt, du);
				absp_dM[j] = vspme->bsp.dbsp_interpolate(ipt, du, 1);
				absp_d2M[j] = vspme->bsp.dbsp_interpolate(ipt, du, 2);
				du += 1;
			}

			du = u[i3 + 1] - iy; //int(u[i3 + 1]);
			absp_M = vspme->absp_M(iatom, 0, 1); 
			absp_dM = vspme->absp_M(iatom, 1, 1);
			absp_d2M = vspme->absp_M(iatom, 2, 1);
			for (j = 0; j < vspme->N; j++) {
				ipt = vspme->bsp.get_pt(du);
				absp_M[j] = vspme->bsp.bsp_interpolate(ipt, du);
				absp_dM[j] = vspme->bsp.dbsp_interpolate(ipt, du, 1);
				absp_d2M[j] = vspme->bsp.dbsp_interpolate(ipt, du, 2);
				du += 1;
			}

			du = u[i3 + 2] - iz; //int(u[i3 + 2]);
			absp_M = vspme->absp_M(iatom, 0, 2); 
			absp_dM = vspme->absp_M(iatom, 1, 2);
			absp_d2M = vspme->absp_M(iatom, 2, 2);
			for (j = 0; j < vspme->N; j++) {
				ipt = vspme->bsp.get_pt(du);
				absp_M[j] = vspme->bsp.bsp_interpolate(ipt, du);
				absp_dM[j] = vspme->bsp.dbsp_interpolate(ipt, du, 1);
				absp_d2M[j] = vspme->bsp.dbsp_interpolate(ipt, du, 2);
				du += 1;
			}
		}
		__syncthreads(); // not really necessary
	};

	__host__ void init_normalized_coordinates(VSPME *vspme_host, VSPME *vspme_dev) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		
		vspme_host->status_host[0] = 0;
		
		dim3 dim3Grid(vspme_host->Nx, vspme_host->Ny);
		int dimBlock = vspme_host->Nz;
		prior_init_normalized_coordinates_kernel<<< dim3Grid, dimBlock, 0, cuStream >>>(vspme_dev);
		cudaStreamSynchronize(cuStream);
		
		dimBlock = 128;
		int dimGrid = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(vspme_host->Na, dimBlock, dimGrid);
		init_normalized_coordinates_kernel<<< dimGrid, dimBlock, 0, cuStream >>>(vspme_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	};

#ifdef __DOUBLE_PRECISION__
	__global__ void VSPME_Q_accum(VSPME *vspme, bool bQ, bool bMU) {
		// work on one cubic of the meshed cell, gridDim is [Nx, Ny] while blockDim is Nz
		int icx = blockIdx.x;
		int icy = blockIdx.y;
		int icz = threadIdx.x;

		int *cubic, nc;

		int mMP = vspme->mMP;
		double q;
		int i, j, k, kx, ky, kz;//, nx, ny, nz;
		int iatom, ipt, i3a, ia;
		double *u = vspme->u.m, *mu = NULL;
		int ib;
		double Mn[3], dMn[3], ux, uy, uz, dux, duy, duz;
		BSplineFunc *bsp = &(vspme->bsp);
		int bsp_indx = vspme->N;
		int nb = bsp_indx;
		int d[3];
		double md[3] = {vspme->md[0], vspme->md[1], vspme->md[2]};
		int Nx = vspme->Nx, Ny = vspme->Ny, Nz = vspme->Nz;

		DOUBLE *Q = vspme->Q.e(icx, icy, icz);
		Q[0] = 0;

		//for (i = -nb; i <= nb; i++) {
		//for (i = -1; i <= nb; i++) {
		for (i = 0; i < nb; i++) {
			kx = icx + i; d[0] = 0;
			if (kx < 0) {kx += Nx; d[0] = -Nx;} else if (kx >= Nx) {kx -= Nx; d[0] = Nx;}
			//for (j = -nb; j <= nb; j++) {
			//for (j = -1; j <= nb; j++) {
			for (j = 0; j < nb; j++) {
				ky = icy + j; d[1] = 0;
				if (ky < 0) {ky += Ny; d[1] = -Ny;} else if (ky >= Ny) {ky -= Ny; d[1] = Ny;}
				//for (k = -nb; k <= nb; k++) {
				//for (k = -1; k <= nb; k++) {
				for (k = 0; k < nb; k++) {
					kz = icz + k; d[2] = 0;
					if (kz < 0) {kz += Nz; d[2] = -Nz;} else if (kz >= Nz) {kz -= Nz; d[2] = Nz;}
					
					cubic = vspme->cmm_cubic(kx, ky, kz);
					nc = cubic[0];
					for (iatom = 1; iatom <= nc; iatom++) {
						ia = cubic[iatom];
						if (ia >= vspme->Na) continue;
						i3a = ia * 3;
						q = vspme->q.m[ia]; mu = vspme->mu.m + i3a;
						ux = u[i3a];// nx = int(ux);
						uy = u[i3a + 1];// ny = int(uy);
						uz = u[i3a + 2];// nz = int(uz);

						mu = vspme->mu.m + i3a;
		
						dux = ux + d[0] - icx; 
						ib = int(dux); if (ib < 0 || ib >= bsp_indx) continue;
						
						//ipt = bsp->get_pt(dux); Mn[0] = bsp->bsp_interpolate(ipt, dux); if (Mn[0] == 0) continue;
						//dMn[0] = bsp->dbsp_interpolate(ipt, dux);
						Mn[0] = vspme->absp_Mn(ia, 0, 0, ib); dMn[0] = vspme->absp_Mn(ia, 1, 0, ib);
						
						duy = uy + d[1] - icy;
						ib = int(duy); if (ib < 0 || ib >= bsp_indx) continue;
						//ipt = bsp->get_pt(duy); Mn[1] = bsp->bsp_interpolate(ipt, duy); if (Mn[1] == 0) continue;
						//dMn[1] = bsp->dbsp_interpolate(ipt, duy);
						Mn[1] = vspme->absp_Mn(ia, 0, 1, ib); dMn[1] = vspme->absp_Mn(ia, 1, 1, ib);

						duz = uz + d[2] - icz;
						ib = int(duz); if (ib < 0 || ib >= bsp_indx) continue;
						//ipt = bsp->get_pt(duz); Mn[2] = bsp->bsp_interpolate(ipt, duz); if (Mn[2] == 0) continue;
						//dMn[2] = bsp->dbsp_interpolate(ipt, duz);
						Mn[2] = vspme->absp_Mn(ia, 0, 2, ib); dMn[2] = vspme->absp_Mn(ia, 1, 2, ib);
						
						if (bQ) {
							//vspme->Q->m[kx].m[ky].m[kz] += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
							Q[0] += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
						}
						if (bMU && mMP > 0) {
							//((FUNC)res.f_mu_accumulate)(patom, vspme, Mn, dMn, res, kx, ky, kz);
							Q[0] += mu[0] * dMn[0] * Mn[1] * Mn[2] / md[0] + mu[1] * Mn[0] * dMn[1] * Mn[2] / md[1] + mu[2] * Mn[0] * Mn[1] * dMn[2] / md[2];
						}
					}
				}
			}
		}

		__syncthreads();
	};

	__host__ void cal_Q_spme(VSPME *vspme_host, VSPME *vspme_dev, bool bQ, bool bMU) {
		dim3 dimGrid(vspme_host->Nx, vspme_host->Ny);
		int dimBlock = vspme_host->Nz;
		
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		VSPME_Q_accum<<<dimGrid, dimBlock, 0, cuStream>>>(vspme_dev, bQ, bMU);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	};

#else
	__global__ void prior_Q_accum_kernel(VSPME *vspme) {
		// work on one cubic of the meshed cell, gridDim is [Nx, Ny] while blockDim is Nz
		int icx = blockIdx.x;
		int icy = blockIdx.y;
		int icz = threadIdx.x;
		int i;

		DOUBLE *Q;// = vspme->Q.e(icx, icy, icz); Q[0] = 0;
		for (i = 0; i < Nflows; i++) {
			Q = vspme->Qbuf[i].e(icx, icy, icz); Q[0] = 0;
		}
	};

	__global__ void sum_Q_kernel(VSPME *vspme) {
		// work on one cubic of the meshed cell, gridDim is [Nx, Ny] while blockDim is Nz
		int icx = blockIdx.x;
		int icy = blockIdx.y;
		int icz = threadIdx.x;
		int i;

		DOUBLE *Q = vspme->Q.e(icx, icy, icz); Q[0] = 0;
		DOUBLE *Qbuf = NULL;
		for (i = 0; i < Nflows; i++) {
			Qbuf = vspme->Qbuf[i].e(icx, icy, icz); Q[0] += Qbuf[0];
		}
	};

	__global__ void VSPME_Q_accum(VSPME *vspme, bool bQ, bool bMU) {
		// work on each atom, send its contribution to close cell
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(vspme->Na, gridDim.x, blockDim.x);

		int icx, icy, icz;

		int *ic = vspme->ic.m;

		int mMP = vspme->mMP;
		double q;
		int i, j, k, kx, ky, kz;//, nx, ny, nz;
		int iatom, i3a, ia;
		double *u = vspme->u.m, *mu = NULL;
		int ib;
		double Mn[3], dMn[3], ux, uy, uz, dux, duy, duz;
		BSplineFunc *bsp = &(vspme->bsp);
		int bsp_indx = vspme->N;
		int nb = bsp_indx;
		double md[3] = {vspme->md[0], vspme->md[1], vspme->md[2]};
		int Nx = vspme->Nx, Ny = vspme->Ny, Nz = vspme->Nz;

		DOUBLE *Q = NULL;
		double Qa = 0;

		double *absp_M, *absp_dM, *absp_d2M;

		double xl[3] = {vspme->xl[0], vspme->xl[1], vspme->xl[2]};

		int iflow = threadIdx.x % Nflows;

		for (iatom = 0; iatom < nEachThread; iatom++) {
			ia = tid * nEachThread + iatom;
			if (ia >= vspme->Na) break;
			i3a = ia * 3;
			icx = ic[i3a]; icy = ic[i3a + 1]; icz = ic[i3a + 2];
			q = vspme->q.m[ia]; mu = vspme->mu.m + i3a;

			for (i = 0; i < nb; i++) {
				kx = icx - i; //dux = ux + i;
				if (kx < 0) {kx += Nx;} //else if (kx >= Nx) {kx -= Nx;}

				ib = i;//int(dux); if (ib < 0 || ib >= bsp_indx) continue;
				Mn[0] = vspme->absp_Mn(ia, 0, 0, ib); dMn[0] = vspme->absp_Mn(ia, 1, 0, ib);

				for (j = 0; j < nb; j++) {
					ky = icy - j; //duy = uy + j;
					if (ky < 0) {ky += Ny;} //else if (ky >= Ny) {ky -= Ny;}

					ib = j;//int(duy); if (ib < 0 || ib >= bsp_indx) continue;
					Mn[1] = vspme->absp_Mn(ia, 0, 1, ib); dMn[1] = vspme->absp_Mn(ia, 1, 1, ib);

					for (k = 0; k < nb; k++) {
						kz = icz - k; //duz = uz + k;
						if (kz < 0) {kz += Nz;} //else if (kz >= Nz) {kz -= Nz;}

						ib = k;//int(duz); if (ib < 0 || ib >= bsp_indx) continue;
						Mn[2] = vspme->absp_Mn(ia, 0, 2, ib); dMn[2] = vspme->absp_Mn(ia, 1, 2, ib);

						//Q = vspme->Q.e(kx, ky, kz);
						Q = vspme->Qbuf[iflow].e(kx, ky, kz);

						Qa = 0;
						if (bQ) {
							//vspme->Q->m[kx].m[ky].m[kz] += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
							Qa += q * Mn[0] * Mn[1] * Mn[2]; // accumulating
						}
						if (bMU && mMP > 0) {
							//((FUNC)res.f_mu_accumulate)(patom, vspme, Mn, dMn, res, kx, ky, kz);
							Qa += mu[0] * dMn[0] * Mn[1] * Mn[2] / md[0] + mu[1] * Mn[0] * dMn[1] * Mn[2] / md[1] + mu[2] * Mn[0] * Mn[1] * dMn[2] / md[2];
						}

						atomicAdd(Q, DOUBLE(Qa));
					}
				}
			}
		}
	};

	__host__ void cal_Q_spme(VSPME *vspme_host, VSPME *vspme_dev, bool bQ, bool bMU) {
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);

		dim3 dim3Grid(vspme_host->Nx, vspme_host->Ny);
		int dimBlock = vspme_host->Nz;
		prior_Q_accum_kernel<<< dim3Grid, dimBlock, 0, cuStream >>>(vspme_dev);
		//cudaStreamSynchronize(cuStream);

		dimBlock = 128;
		int dimGrid = 0;
		_gpu_dstruct_::_gpu_util_::gpu_grid(vspme_host->Na, dimBlock, dimGrid);
		VSPME_Q_accum<<< dimGrid, dimBlock, 0, cuStream >>>(vspme_dev, bQ, bMU);
		//cudaStreamSynchronize(cuStream);

		//dim3Grid(vspme_host->Nx, vspme_host->Ny);
		dimBlock = vspme_host->Nz;
		sum_Q_kernel<<< dim3Grid, dimBlock, 0, cuStream >>>(vspme_dev);
		cudaStreamSynchronize(cuStream);

		cudaStreamDestroy(cuStream);
	};

#endif


/***************************************************************************************
	CALCULATE E_recp -- E_recp = F_recp / q
	Mesh dimension is assumed to be bigger than the order of the B-spline function
	In another word, the spline points is within the nearest neighbor cell only
***************************************************************************************/

	__global__ void EF_recp(VSPME *vspme, bool bE, bool bF, bool bQ, bool bMU) { // E of atom is accumulated ,without initilized
		// work on each atom, with 1d dimension
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = _gpu_dstruct_::_gpu_util_::gpuJob_nEachThread(vspme->Na, gridDim.x, blockDim.x);

		int i, j, k, nx, ny, nz, kx, ky, kz;
		int ia, i3a, iatom, ipt;
		int ib;
		double Mn[3], dMn[3], d2Mn[3];
		double ux, uy, uz, dux, duy, duz;
		BSplineFunc *bsp = &(vspme->bsp);
		int bsp_indx = vspme->N;
		int nb = bsp_indx - 1;

		double md[3] = {vspme->md[0], vspme->md[1], vspme->md[2]};
		double *E, *F, *u;
		double tv_q[3] = {0, 0, 0}, tv_mu[3] = {0, 0, 0};
		double q, *mu, bcq;
		int Nx = vspme->Nx, Ny = vspme->Ny, Nz = vspme->Nz;

		double inv_V = vspme->inv_V;
		_gpu_dstruct_::FLEX_CMATRIX<DOUBLE> *BCQ = &(vspme->BCQ);
		int mMP = vspme->mMP;

		for (iatom = 0; iatom < nEachThread; iatom++) {
			ia = tid * nEachThread + iatom;
			if (ia >= vspme->Na) break;
			i3a = ia * 3;
			if (bE) {E = vspme->E.m + i3a; E[0] = 0; E[1] = 0; E[2] = 0;}
			if (bF) {F = vspme->F.m + i3a; F[0] = 0; F[1] = 0; F[2] = 0;}
			u = vspme->u.m + i3a;
			
			ux = u[0]; nx = int(ux);
			uy = u[1]; ny = int(uy);
			uz = u[2]; nz = int(uz);
			
			//for (i = -nb; i <= nb; i++) {
			//for (i = -nb; i <= 1; i++) {
			for (i = -nb; i <= 0; i++) {
				kx = nx + i; dux = ux - kx; //if (dux < 0 || dux >= bsp_indx) continue;
				//ipt = bsp->get_pt(dux); Mn[0] = bsp->bsp_interpolate(ipt, dux);
				//dMn[0] = bsp->dbsp_interpolate(ipt, dux);
				//if (bF && mMP > 0) d2Mn[0] = bsp->dbsp_interpolate(ipt, dux, 2);
				ipt = int(dux); if (ipt < 0 || ipt >= bsp_indx) continue;
				Mn[0] = vspme->absp_Mn(ia, 0, 0, ipt); dMn[0] = vspme->absp_Mn(ia, 1, 0, ipt); d2Mn[0] = vspme->absp_Mn(ia, 2, 0, ipt);
				if (kx < 0) kx += Nx; else if (kx >= Nx) kx -= Nx;

				//for (j = -nb; j <= nb; j++) {
				//for (j = -nb; j <= 1; j++) {
				for (j = -nb; j <= 0; j++) {
					ky = ny + j; duy = uy - ky; //if (duy < 0 || duy >= bsp_indx) continue;
					//ipt = bsp->get_pt(duy); Mn[1] = bsp->bsp_interpolate(ipt, duy);
					//dMn[1] = bsp->dbsp_interpolate(ipt, duy);
					//if (bF && mMP > 0) d2Mn[1] = bsp->dbsp_interpolate(ipt, duy, 2);
					ipt = int(duy); if (ipt < 0 || ipt >= bsp_indx) continue;
					Mn[1] = vspme->absp_Mn(ia, 0, 1, ipt); dMn[1] = vspme->absp_Mn(ia, 1, 1, ipt); d2Mn[1] = vspme->absp_Mn(ia, 2, 1, ipt);
					if (ky < 0) ky += Ny; else if (ky >= Ny) ky -= Ny;

					//for (k = -nb; k <= nb; k++) {
					//for (k = -nb; k <= 1; k++) {
					for (k = -nb; k <= 0; k++) {
						kz = nz + k; duz = uz - kz; //if (duz < 0 || duz >= bsp_indx) continue;
						//ipt = bsp->get_pt(duz); Mn[2] = bsp->bsp_interpolate(ipt, duz);
						//dMn[2] = bsp->dbsp_interpolate(ipt, duz);
						//if (bF && mMP > 0) d2Mn[2] = bsp->dbsp_interpolate(ipt, duz, 2);
						ipt = int(duz); if (ipt < 0 || ipt >= bsp_indx) continue;
						Mn[2] = vspme->absp_Mn(ia, 0, 2, ipt); dMn[2] = vspme->absp_Mn(ia, 1, 2, ipt); d2Mn[2] = vspme->absp_Mn(ia, 2, 2, ipt);
						if (kz < 0) kz += Nz; else if (kz >= Nz) kz -= Nz;

						bcq = (*(BCQ->e(kx, ky, kz)));
						if (bE) {// accumulating E
							E[0] -= dMn[0] / md[0] * Mn[1] * Mn[2] * bcq; //BCQ->m[kx].m[ky].m[kz]; // accumulating
							E[1] -= Mn[0] * dMn[1] / md[1] * Mn[2] * bcq; //BCQ->m[kx].m[ky].m[kz]; // accumulating
							E[2] -= Mn[0] * Mn[1] * dMn[2] / md[2] * bcq; //BCQ->m[kx].m[ky].m[kz]; // accumulating
						}

						if (bF) {// accumulating F
							//((FUNC)var.fFunc)(patom, vspme, Mn, dMn, d2Mn, var, kx, ky, kz);
							//memset(tv_q, 0, SIZE_V3); memset(tv_mu, 0, SIZE_V3);
							tv_q[0] = 0; tv_q[1] = 0; tv_q[2] = 0; tv_mu[0] = 0; tv_mu[1] = 0; tv_mu[2] = 0;
							q = vspme->q.m[ia];
							mu = vspme->mu.m + i3a;
							if (bQ && q != 0) {
								tv_q[0] = q * dMn[0] / md[0] * Mn[1] * Mn[2];
								tv_q[1] = q * Mn[0] * dMn[1] / md[1] * Mn[2];
								tv_q[2] = q * Mn[0] * Mn[1] * dMn[2] / md[2];
							}
							if (bMU && mMP > 0) {
								tv_mu[0] = mu[0] * d2Mn[0] * Mn[1] * Mn[2] / (md[0] * md[0]) + mu[1] * dMn[0] * dMn[1] * Mn[2] / (md[0] * md[1]) + mu[2] * dMn[0] * Mn[1] * dMn[2] / (md[0] * md[2]);
								tv_mu[1] = mu[0] * dMn[0] * dMn[1] * Mn[2] / (md[0] * md[1]) + mu[1] * Mn[0] * d2Mn[1] * Mn[2] / (md[1] * md[1]) + mu[2] * Mn[0] * dMn[1] * dMn[2] / (md[1] * md[2]);
								tv_mu[2] = mu[0] * dMn[0] * Mn[1] * dMn[2] / (md[0] * md[2]) + mu[1] * Mn[0] * dMn[1] * dMn[2] / (md[1] * md[2]) + mu[2] * Mn[0] * Mn[1] * d2Mn[2] / (md[2] * md[2]);
							}
							F[0] -= (tv_q[0] + tv_mu[0]) * bcq; //BCQ->m[kx].m[ky].m[kz]; 
							F[1] -= (tv_q[1] + tv_mu[1]) * bcq; //BCQ->m[kx].m[ky].m[kz];
							F[2] -= (tv_q[2] + tv_mu[2]) * bcq; //BCQ->m[kx].m[ky].m[kz];
						}
					}
				}
			}
			if (bE) {
				//memcpy(vspme->E_host + i3a, E, SIZE_V3);
				vspme->E_host[i3a] = inv_V * E[0];
				vspme->E_host[i3a + 1] = inv_V * E[1];
				vspme->E_host[i3a + 2] = inv_V * E[2];
			}
			if (bF) {
				//memcpy(vspme->F_host + i3a, F, SIZE_V3);
				vspme->F_host[i3a] = inv_V * F[0];
				vspme->F_host[i3a + 1] = inv_V * F[1];
				vspme->F_host[i3a + 2] = inv_V * F[2];
			}
		}
	};

	__host__ void cal_EF_spme(VSPME *vspme_host, VSPME *vspme_dev, bool bQ, bool bMU, bool bE, bool bF) {
		int dimBlock = 128, dimGrid;
		_gpu_dstruct_::_gpu_util_::gpu_grid(vspme_host->Na, dimBlock, dimGrid);
		

		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);
		EF_recp<<<dimGrid, dimBlock, 0, cuStream>>>(vspme_dev, bE, bF, bQ, bMU);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
	};
	
	
	__global__ void _init_device_spme(VSPME *vspme) {
		int tid = threadIdx.x + blockDim.x * blockIdx.x;
		if (tid == 0) {
			vspme->_init();
		}
		__syncthreads();
	};
	
	__global__ void init_device_spme(VSPME* vspme) {
		int tid = threadIdx.x + blockDim.x * blockIdx.x;
		if (tid == 0) {
			vspme->init_k(); vspme->init_b(); vspme->init_C(); 
		}
		__syncthreads();
	};

	// bsp_diff2_mapped_host is the 2nd differential profile of BSplineFunc
	// because GPU can not run self-iterative function, the 2nd differential profile can not be obtained
	// it will have to be copied from the host function
	// r_hostmap is the GPU mapped memory, on host
	__host__ bool setup_vspme_host(VSPME *vspme_host, VSPME** vspme_dev, bool b_init, float *status_host, int bsp_indx, double *bsp_x, double *bsp_y0, double *bsp_y1, double *bsp_y2, int bsp_npts, float *xl, float *xd, double kappa, int Nx, int Ny, int Nz, int nBSP, short mMP, double *r_hostmap, double *q_hostmap, double *mu_hostmap, double *E_hostmap, double *F_hostmap, int Na) {
		if (b_init) vspme_host->_init(); // important
		
		vspme_host->status_host = status_host;
		vspme_host->mMP = mMP;
		vspme_host->r = r_hostmap;
		vspme_host->q_host = q_hostmap;
		vspme_host->mu_host = mu_hostmap;
		vspme_host->E_host = E_hostmap;
		vspme_host->F_host = F_hostmap;
		vspme_host->Nx = Nx;
		vspme_host->Ny = Ny;
		vspme_host->Nz = Nz;
		vspme_host->Nyz = Ny * Nz;
		vspme_host->Nv = Nx * Ny * Nz;
		vspme_host->xl[0] = xl[0];
		vspme_host->xl[1] = xl[1];
		vspme_host->xl[2] = xl[2];
		vspme_host->xd[0] = xd[0];
		vspme_host->xd[1] = xd[1];
		vspme_host->xd[2] = xd[2];
		vspme_host->V = xd[0] * xd[1] * xd[2];
		vspme_host->kappa = kappa;
		
		vspme_host->cell_dim(xd[0], xd[1], xd[2]);

		vspme_host->N = nBSP;
		vspme_host->Na = Na;
		
		// BSplineFunc
		vspme_host->bsp._init();
		vspme_host->bsp.set_range(bsp_indx, bsp_npts);
		vspme_host->bsp.set_profiles_host(bsp_npts);
		cudaMemcpy(vspme_host->bsp.x, bsp_x, bsp_npts * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vspme_host->bsp.y0, bsp_y0, bsp_npts * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vspme_host->bsp.y1, bsp_y1, bsp_npts * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vspme_host->bsp.y2, bsp_y2, bsp_npts * sizeof(double), cudaMemcpyHostToDevice);

		
		vspme_host->Nca = int(xd[0] * xd[1] * xd[2] / (Nx * Ny * Nz) + 0.1) + 5;
		if (vspme_host->Nca < 5) vspme_host->Nca = 5;
		
		size_t size_free = 0, size_total = 0;
		//cudaMemGetInfo(&size_free, &size_total);
		
		if (!vspme_host->bUseHostFFT) vspme_host->init_fft_plans();
		if (!vspme_host->init_buffer_host()) return false;
		
		//cudaMemGetInfo(&size_free, &size_total);
		
		cudaStream_t cuStream;
		cudaStreamCreate(&cuStream);

		cudaMalloc(vspme_dev, sizeof(VSPME));
		_init_device_spme<<<1, 1, 0, cuStream>>>(*vspme_dev);
		cudaStreamSynchronize(cuStream);
		
		cudaMemcpy(*vspme_dev, vspme_host, sizeof(VSPME), cudaMemcpyHostToDevice);
		
		//cudaMemGetInfo(&size_free, &size_total);
		
		init_device_spme<<<1, 1, 0, cuStream>>>(*vspme_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		
		//cudaMemGetInfo(&size_free, &size_total);

		// we need to copy the device vspme back to host, to make sure the parameters are consistent, eg. the device address of arrays, and their dimensions
		// so that we can update the device memory
		size_total = sizeof(VSPME);
		cudaMemcpy(vspme_host, *vspme_dev, sizeof(VSPME), cudaMemcpyDeviceToHost);

		cudaMemcpy(vspme_host->q.m, q_hostmap, Na * sizeof(double), cudaMemcpyHostToDevice);
		
		cudaMemGetInfo(&size_free, &size_total);
		if (size_free == 0) {
			// some problem to allocate memory on GPU, during setting up the vspme on GPU
			return false;
		}
		else return true;
	};

	__host__ void update_dipole_device_spme(VSPME *vspme_host) {
		cudaMemcpy(vspme_host->mu.m, vspme_host->mu_host, vspme_host->mu.n * sizeof(double), cudaMemcpyDeviceToDevice);
	};

	__host__ void update_q_device_spme(VSPME *vspme_host) {
		cudaMemcpy(vspme_host->q.m, vspme_host->q_host, vspme_host->q.n * sizeof(double), cudaMemcpyHostToDevice);
	};

	__host__ void release_device_spme(VSPME *vspme_host, VSPME *vspme_dev) {
		vspme_host->release_plans();
		vspme_host->release_host();
		
		//cudaMemcpy(vspme_host, vspme_dev, sizeof(VSPME), cudaMemcpyDeviceToHost);
		cudaFree(vspme_dev);
	};

#endif // __GPU__
} // end of namespace _spme_

} // _gpu_md_

