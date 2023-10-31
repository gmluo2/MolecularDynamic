/*
 * gpu_def.h
 *
 *  Created on: May 20, 2013
 *      Author: gmluo
 */

#ifdef __GPU__

// for the cases like Q calculation of meshed grid, requiring lots of atomic operation
// normally, it can be warpSize, or smaller
#define Nflows  32

// using double precision for GPU calculation ?
//#define __DOUBLE_PRECISION__
#ifdef __DOUBLE_PRECISION__
#define GPU_DOUBLE double
#define DOUBLE double
#define COMPLEX cuDoubleComplex
#define CUFFT_PLAN_D2Z  CUFFT_D2Z
#define CUFFT_PLAN_Z2D  CUFFT_Z2D
#define CUFFT_ExecD2Z  cufftExecD2Z
#define CUFFT_ExecZ2D  cufftExecZ2D
#define SUM_ARRAY2DD_CUML   sum_array_cuml
#define SUM_ARRAY2Dd_CUML   sum_array_cuml

#define GPUSQRT  sqrt
#define GPUERFC  erfc
#define GPUEXP   exp

#else

#define GPU_DOUBLE float
#define DOUBLE float
#define COMPLEX cuComplex
#define CUFFT_PLAN_D2Z  CUFFT_R2C
#define CUFFT_PLAN_Z2D  CUFFT_C2R
#define CUFFT_ExecD2Z  cufftExecR2C
#define CUFFT_ExecZ2D  cufftExecC2R
#define SUM_ARRAY2DD_CUML   sum_float_array_cuml
#define SUM_ARRAY2Dd_CUML   sum_float_double_array_cuml

#define GPUSQRT  sqrtf
#define GPUERFC  erfcf
#define GPUEXP   expf

#endif  //__DOUBLE_PRECISION__

#endif  // __GPU__
