//#if _MPI_ == 0 // NON-MPI

#if _SYS_ == _WINDOWS_SYS_
#include <stdio.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <time.h>

#define MYTIME struct _timeb
#elif _SYS_ == _LINUX_SYS_
#include <stdio.h>
//#include <time.h>
#include <sys/timeb.h>
#define MYTIME struct timeb
#endif
/*
#elif _MPI_ == 1 //MPI

#ifdef SEEK_SET
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#endif
#include "mpi.h"

#define MYTIME  double

#endif //_MPI_
*/
#ifndef _DEF_TIME_STAMP_
#define _DEF_TIME_STAMP_
class TIME {
	MYTIME t0, t1;
public:
	void start();
	void stop();
	long glance(bool bLOG = false, int iSTAMP = 0);
	long elapse(bool bLOG = false, int iSTAMP = 0);
};



#ifndef _CHECK_TIME_STAMPE_
#define _CHECK_TIME_STAMPE_  1
#endif

#ifdef _MAIN_TIME_
	#define _DEF_TIME_STAMP_
	TIME mt;
#elif _CHECK_TIME_STAMP_
	#ifndef _DEF_GTIME_STAMP_
	#define _DEF_GTIME_STAMP_
	extern TIME mt;
	#endif
#endif


#endif
