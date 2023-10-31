#define MAX_THREADS   4

#define _DOS_SYS_       0
#define _WINDOWS_SYS_   1
#define _LINUX_SYS_     2

//#define _USE_GPU_

#define _SYS_   1

// if disable GPU, remove the definition of __GPU__
//#define __GPU__
#ifdef __GPU__
// use GPU for real EwaldSum calculation ? true : false
#define bGPURealEwaldSum true

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
//#include <cuFFT.h>
#include <cufft.h>

#define _CudaDevice_ 0  // the device

// use CPU for SPME calculation ? enable -- define __CPU_FFT__, disable -- disable the definition of __CPU_FFT__
//#define __CPU_FFT__
#endif //__GPU__

#include "gpu_def.h"

#if _SYS_ == _WINDOWS_SYS_
//#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
//#define _CRT_SECURE_CPP_OVERLOAD_SECURE_NAMES 1
#include "..\stdafx.h"
#include "..\NEIMOMD_DLL.h"
#include "..\MDInforDlg.h"

//#define _CRT_SECURE_DEPRECATE_MEMORY
//#include <memory.h>
//#define memcpy CopyMemory
#include <string.h>
//#define _CRT_SECURE_NO_DEPRECATE 1
//#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1

#elif _SYS_ == _LINUX_SYS_
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>

#define LPVOID void

#endif

#include <iostream>
#include <sstream>
#include <ostream>
#include <streambuf>
// use MPI calculation
#define _MPI_    0

// ignore the electrostatic interaction between atoms ? [0/1]
#define _IGNORE_ELECTROSTATIC_INTERACTION_   0
// distinguish intra-molecule and inter-molecule interaction for Coarse-Grained MD [0/1]
#define _INTRA_INTER_INTERACTS      1

// disable some old functions, which could always conflict with the useful code

// always keep it ON
#define _DISABLE_  1
