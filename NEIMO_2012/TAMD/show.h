#include "log.h"

extern void show_msg(char *msg, bool endline = true); // show the message in MessageBox and log file
void input(char *msg, int length = 250);

void show_loop_infor(int nloop, int nsaved, int nLoops, float Ek, float Ep);
void show_infor(char *msg, bool endline = true); // showing to MD dialog and log file
void show_log(char *msg, bool endline); // show the message in log file
void show_log1(char *msg, bool endline, bool bReplaceLog = false);

template<class T> void show_vlog(T v, bool endline = true) {
	char msg[256] = "\0";
	extern LOG mlog;
	if (mlog.out != NULL) {
		(*mlog.out)<<v;
		if (endline) (*mlog.out)<<endl;
	}
};

