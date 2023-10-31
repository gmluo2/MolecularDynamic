#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
using namespace std;

#include "show.h"
extern LOG mlog;

#if _SYS_ == _DOS_SYS_ || _SYS_ == _LINUX_SYS_

void show_msg(char *msg, bool endline) {
	mlog.show(msg, endline);
	//if (endline) mlog.show("\n", false);
}

void show_log(char *msg, bool endline) {
	mlog.show(msg, endline);
	//if (endline) mlog.show("\n", false);
}

void show_log1(char *msg, bool endline, bool bReplaceLog) {
	mlog.show(msg, endline, false, bReplaceLog);
}

void input(char *msg, int length) {
	cin.getline(msg, length);
}
#elif _SYS_ == _WINDOWS_SYS_
void show_msg(char *msg, bool endline) {
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	AfxMessageBox(LPCTSTR(CString(msg)));
	mlog.show(msg, endline, false, endline);
	//if (endline) mlog.show("\n", false);
};

void show_log(char *msg, bool endline) {
	mlog.show(msg, endline, false, endline);
	//if (endline) mlog.show("\n", false);
}

void show_log1(char *msg, bool endline, bool bReplaceLog) {
	mlog.show(msg, endline, false, bReplaceLog);
}

void input(char *msg, int length) {
}
#endif
