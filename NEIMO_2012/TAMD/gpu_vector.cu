#include "project.h"

//#include <iostream>
//#include <fstream>
#include <cstdio>
//#include <cmath>

using namespace std;

#ifdef __GPU__

#include "def.h"
#include "gpu_vector.h"
namespace _gpu_dstruct_ {

/*
void rand_vect3(double *v) {
	double theta = 0, p = 0;
	while (1) {
		theta = PI * ranf();
		p = sin(theta); // theta is [0, PI], sin(theta) [0, 1]
		if (p < 0) p = -p;
		if (ranf() <= p) break; // the possibility at theta is sin(theta)
	}
	double omega = PI2 * ranf();
	v[2] = (float)cos(theta);
	double xy = sin(theta);
	v[0] = float(xy * cos(omega));
	v[1] = float(xy * sin(omega));
}
*/




namespace _gpu_util_ {
	__global__ void sum_array_kernel(double *a_dev, int n, double *s_block) {
		__shared__ double sum[512]; // threads in block is less than 512
		sum[threadIdx.x] = 0;

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = gpuJob_nEachThread(n, gridDim.x, blockDim.x);
		int i, ia;
		double sv = 0;
		
		for (i = 0; i < nEachThread; i++) { // each thread
			ia = tid * nEachThread + i;
			if (ia >= n) break;
			sum[threadIdx.x] += a_dev[ia];
		}
		__syncthreads();
		if (threadIdx.x == 0) {// each block
			for (i = 0; i < blockDim.x; i++) sv += sum[i];
			s_block[blockIdx.x] = sv;
		}
		__syncthreads();
	};

	__global__ void sum_array_cuml_kernel(double *a1_dev, double *a2_dev, int n, double *s_block) {
		__shared__ double sum[512];  // threads in block is less than 512
		sum[threadIdx.x] = 0;

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = gpuJob_nEachThread(n, gridDim.x, blockDim.x);
		int i, ia;
		double sv = 0;
		
		for (i = 0; i < nEachThread; i++) { // each thread
			ia = tid * nEachThread + i;
			if (ia >= n) break;
			sum[threadIdx.x] += a1_dev[ia] * a2_dev[ia];
		}
		__syncthreads();
		if (threadIdx.x == 0) {// each block
			for (i = 0; i < blockDim.x; i++) sv += sum[i];
			s_block[blockIdx.x] = sv;
		}
		__syncthreads();
	};

	__host__ double sum_array(double *a_dev, int n, double *buf_dev, double *buf_host, int dim_buf) {
		int blockDim = 256;
		int gridDim = n / blockDim;
		if (gridDim == 0) gridDim = 1;
		else if (gridDim > 4096) gridDim = 4096; // do not have too many blocks, but also can not be just few blocks
		if (gridDim > dim_buf) gridDim = dim_buf;

		cudaStream_t cuStream = NULL;
		cudaStreamCreate(&cuStream);
		sum_array_kernel<<< gridDim, blockDim, 512, cuStream >>>(a_dev, n, buf_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		double sum = 0;
		for (int i = 0; i < gridDim; i++) sum += buf_host[i];
		return sum;
	};

	__host__ double sum_array_cuml(double *a1_dev, double *a2_dev, int n, double *buf_dev, double *buf_host, int dim_buf) {
		int blockDim = 256;
		int gridDim = n / blockDim;
		if (gridDim == 0) gridDim = 1;
		else if (gridDim > 4096) gridDim = 4096; // do not have too many blocks, but also can not be just few blocks
		if (gridDim > dim_buf) gridDim = dim_buf;

		cudaStream_t cuStream = NULL;
		cudaStreamCreate(&cuStream);
		sum_array_cuml_kernel<<< gridDim, blockDim, 512, cuStream >>>(a1_dev, a2_dev, n, buf_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		double sum = 0;
		for (int i = 0; i < gridDim; i++) sum += buf_host[i];
		return sum;
	};







	__global__ void sum_float_array_kernel(float *a_dev, int n, float *s_block) {
		__shared__ double sum[512]; // threads in block is less than 512
		sum[threadIdx.x] = 0;

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = gpuJob_nEachThread(n, gridDim.x, blockDim.x);
		int i, ia;
		float sv = 0;

		for (i = 0; i < nEachThread; i++) { // each thread
			ia = tid * nEachThread + i;
			if (ia >= n) break;
			sum[threadIdx.x] += a_dev[ia];
		}
		__syncthreads();
		if (threadIdx.x == 0) {// each block
			for (i = 0; i < blockDim.x; i++) sv += sum[i];
			s_block[threadIdx.x] = sv;
		}
		__syncthreads();
	};

	__global__ void sum_float_array_cuml_kernel(float *a1_dev, float *a2_dev, int n, float *s_block) {
		__shared__ float sum[512];  // threads in block is less than 512
		sum[threadIdx.x] = 0;

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = gpuJob_nEachThread(n, gridDim.x, blockDim.x);
		int i, ia;
		float sv = 0;

		for (i = 0; i < nEachThread; i++) { // each thread
			ia = tid * nEachThread + i;
			if (ia >= n) break;
			sum[threadIdx.x] += a1_dev[ia] * a2_dev[ia];
		}
		__syncthreads();
		if (threadIdx.x == 0) {// each block
			for (i = 0; i < blockDim.x; i++) sv += sum[i];
			s_block[blockIdx.x] = sv;
		}
		__syncthreads();
	};
	__global__ void sum_float_double_array_cuml_kernel(float *a1_dev, double *a2_dev, int n, float *s_block) {
		__shared__ float sum[512];  // threads in block is less than 512
		sum[threadIdx.x] = 0;

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = gpuJob_nEachThread(n, gridDim.x, blockDim.x);
		int i, ia;
		float sv = 0;

		for (i = 0; i < nEachThread; i++) { // each thread
			ia = tid * nEachThread + i;
			if (ia >= n) break;
			sum[threadIdx.x] += a1_dev[ia] * a2_dev[ia];
		}
		__syncthreads();
		if (threadIdx.x == 0) {// each block
			for (i = 0; i < blockDim.x; i++) sv += sum[i];
			s_block[blockIdx.x] = sv;
		}
		__syncthreads();
	};
	__host__ float sum_float_array(float *a_dev, int n, float *buf_dev, float *buf_host, int dim_buf) {
		int blockDim = 256;
		int gridDim = n / blockDim;
		if (gridDim == 0) gridDim = 1;
		else if (gridDim > 512) gridDim = 512; // do not have too many blocks, but also can not be just few blocks
		if (gridDim > dim_buf) gridDim = dim_buf;

		cudaStream_t cuStream = NULL;
		cudaStreamCreate(&cuStream);
		sum_float_array_kernel<<< gridDim, blockDim, 512, cuStream >>>(a_dev, n, buf_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		float sum = 0;
		for (int i = 0; i < gridDim; i++) sum += buf_host[i];
		return sum;
	};
	__host__ float sum_float_array_cuml(float *a1_dev, float *a2_dev, int n, float *buf_dev, float *buf_host, int dim_buf) {
		int blockDim = 256;
		int gridDim = n / blockDim;
		if (gridDim == 0) gridDim = 1;
		else if (gridDim > 4096) gridDim = 4096; // do not have too many blocks, but also can not be just few blocks
		if (gridDim > dim_buf) gridDim = dim_buf;

		cudaStream_t cuStream = NULL;
		cudaStreamCreate(&cuStream);
		sum_float_array_cuml_kernel<<< gridDim, blockDim, 512, cuStream >>>(a1_dev, a2_dev, n, buf_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		double sum = 0;
		for (int i = 0; i < gridDim; i++) sum += buf_host[i];
		return sum;
	};

	__host__ float sum_float_double_array_cuml(float *a1_dev, double *a2_dev, int n, float *buf_dev, float *buf_host, int dim_buf) {
		int blockDim = 256;
		int gridDim = n / blockDim;
		if (gridDim == 0) gridDim = 1;
		else if (gridDim > 4096) gridDim = 4096; // do not have too many blocks, but also can not be just few blocks
		if (gridDim > dim_buf) gridDim = dim_buf;

		cudaStream_t cuStream = NULL;
		cudaStreamCreate(&cuStream);
		sum_float_double_array_cuml_kernel<<< gridDim, blockDim, 512, cuStream >>>(a1_dev, a2_dev, n, buf_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		double sum = 0;
		for (int i = 0; i < gridDim; i++) sum += buf_host[i];
		return sum;
	};



} // _gpu_util_

__device__ double atomicAdd(double *address, double val) {
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	while(assumed != old) {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
	}
	return __longlong_as_double(old);
};

} // end of namespace _gpu_dstruct_

#endif // __GPU__
