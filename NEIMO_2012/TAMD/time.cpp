#include "project.h"

//#if _MPI_ == 0  // NON-MPI

#if _SYS_ == _WINDOWS_SYS_

#elif _SYS_ == _LINUX_SYS_

#endif

#define _MAIN_TIME_ 1
#include "time.h"

extern void show_log(char *msg, bool endline);

#if _SYS_ == _WINDOWS_SYS_
void TIME::start() {
	_ftime_s(&t0);
}

void TIME::stop() {
	_ftime_s(&t1);
}

long TIME::glance(bool bLOG, int iSTAMP) {
	MYTIME now;
	_ftime_s(&now);
	long dt = (now.time - t0.time) * 1000 + (now.millitm - t0.millitm);
	char msg[256] = "\0";
	if (bLOG) {
		sprintf(msg, "%d : %d", iSTAMP, dt); show_log(msg, false);
	}
	return dt;
}

long TIME::elapse(bool bLOG, int iSTAMP) {
	_ftime_s(&t1);
	long dt = (t1.time - t0.time) * 1000 + (t1.millitm - t0.millitm);
	char msg[256] = "\0";
	if (bLOG) {
		sprintf(msg, "%d : %d", iSTAMP, dt); show_log(msg, false);
	}
	return dt;
}

#elif _SYS_ == _LINUX_SYS_
void TIME::start() {
  ftime(&t0);
}

void TIME::stop() {
  ftime(&t1);
}

long TIME::glance(bool bLOG, int iSTAMP) {
  MYTIME now;
  ftime(&now);
 // return (now.time - t0.time) * 1000 + (now.millitm - t0.millitm);
  long dt = (now.time - t0.time) * 1000 + (now.millitm - t0.millitm);
	char msg[256] = "\0";
	if (bLOG) {
		sprintf(msg, "%d : %d", iSTAMP, dt); show_log(msg, false);
	}
	return dt;
}

long TIME::elapse(bool bLOG, int iSTAMP) {
  ftime(&t1);
  //return (t1.time - t0.time) * 1000 + (t1.millitm - t0.millitm);
  long dt = (t1.time - t0.time) * 1000 + (t1.millitm - t0.millitm);
	char msg[256] = "\0";
	if (bLOG) {
		sprintf(msg, "%d : %d", iSTAMP, dt); show_log(msg, false);
	}
	return dt;
}

#endif
/*
#elif _MPI_ == 1
#include "time.h"

void TIME::start() {
	t0 = MPI_Wtime() * 1000; // to ms
}

void TIME::stop() {
  t1 = MPI_Wtime() * 1000; // to ms
}

long TIME::glance() {
  MYTIME now;
  now = MPI_Wtime() * 1000; // to ms
  return long(now - t0);
}

long TIME::elapse() {
  t1 = MPI_Wtime() * 1000; // to ms
  return long(t1 - t0);
}
#endif //_MPI
*/


