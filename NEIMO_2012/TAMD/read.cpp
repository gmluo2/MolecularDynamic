#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
#include "show.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"

#include "ZMatrix.h"

#include "MM.h"
#include "Mol.h"
#include "MD.h"

#include "read.h"
#include "var.h"

extern bool search(ifstream &in, char *title, char *buffer);

void rewind(ifstream &in) {
	in.clear();
#if _SYS_ == _WINDOWS_SYS_
	in.seekg(0, ios_base::beg);
#elif _SYS_ == _LINUX_SYS_
	in.seekg(0, ios::beg);
#endif
}

file_pos current_fpos(ifstream &in) {
	return in.tellg();
}

void set_fpos(ifstream &in, file_pos &pos) {
	rewind(in);
	in.seekg(pos);
}

void Upper(char* str)
{
	int len = int(strlen(str));
	for (int i = 0; i < len; i++) {
		if (*(str + i) <= 'z' && *(str + i) >= 'a') *(str + i) =  *(str + i) + 'A' - 'a';
	}
}

bool NextSection(char** original, char ch)  //default ch = ' '
{
	//go to next section seperated by " " in str
	int i = 0;
	char* str = *original;
	while(strlen(str) > 0 && *str == ch) str++;
	if (strlen(str) == 0) {*original = str; return false;}
	while(*str != ch && strlen(str) > 0) str++;
	*original = str;
	if (strlen(str) == 0) return false;
	else return true;
}

bool pass(char** original, char ch) {
	int i = 0;
	char *str = *original;
	while (strlen(str) > 0 && *str == ch) str++;
	*original = str;
	if (strlen(str) == 0) return false; // reach end
	else return true;
}

void kick_off_tail(char *str, char ch) {
	int n = 0;
	while ((n = strlen(str) - 1) >= 0 && *(str + n) == ch) {
		*(str + n) = '\0';
	}
}

int read_float(char *buffer, float *v, int max_v, char* format) {
	int n = 0;
	char *buff = NULL;
	buff = buffer;
	while (n < max_v && buff != NULL && sscanf(buff, format, v + n) == 1) {
		NextSection(&buff, ' '); n++;
	}
	return n;
}

int read_float(char *buffer, float **ta) {
	if (*ta != NULL) {delete[] *ta; *ta = NULL;}
	float v = 0;
	int n = 0;
	char *pch = buffer;
	while (pch != NULL && sscanf(pch, "%f", &v) == 1) {NextSection(&pch); n++;}
	if (n == 0) return 0;
	*ta = new float[n];
	pch = buffer; n = 0;
	while (pch != NULL && sscanf(pch, "%f", &v) == 1) {(*ta)[n] = v; NextSection(&pch); n++;}
	return n;
}
/*
void Upper(char* str)
{
	int len = int(strlen(str));
	for (int i = 0; i < len; i++) {
		if (*(str + i) <= 'z' && *(str + i) >= 'a') *(str + i) =  *(str + i) + 'A' - 'a';
	}
}

bool NextSection(char** original, char ch)  //default ch = ' '
{
	//go to next section seperated by " " in str
	int i = 0;
	char* str = *original;
	while(strlen(str) > 0 && *str == ch) str++;
	if (strlen(str) == 0) {*original = str; return false;}

	while(*str != ch && strlen(str) > 0) str++;

	while(strlen(str) > 0 && *str == ch) str++;
	if (strlen(str) == 0) {*original = str; return false;}

	*original = str;
	if (strlen(str) == 0) return false;
	else return true;
}
*/
bool InString(char ch, char *str) {
	int nstr = strlen(str);
	for (int i = 0; i < nstr; i++) {
		if (ch == str[i]) return true;
	}
	return false;
}

bool NextSection(char** original, char* ch)  
{
	//go to next section seperated by " " in str
	int i = 0;
	char* str = *original;
	while(strlen(str) > 0 && InString(*str, ch)) str++;
	if (strlen(str) == 0) {*original = NULL; return false;}

	while(!InString(*str, ch) && strlen(str) > 0) str++;

	while(strlen(str) > 0 && InString(*str, ch)) str++;
	if (strlen(str) == 0) {*original = NULL; return false;}

	*original = str;
	if (strlen(str) == 0) return false;
	else return true;
}
/*
bool pass(char** original, char ch) {
	int i = 0;
	char *str = *original;
	while (strlen(str) > 0 && *str == ch) str++;
	*original = str;
	if (strlen(str) == 0) return false; // reach end
	else return true;
}
*/
bool pass(char** original, char* ch) {
	int i = 0;
	char *str = *original;
	while (strlen(str) > 0 && InString(*str, ch)) str++;
	*original = str;
	if (strlen(str) == 0) return false; // reach end
	else return true;
}
/*
void kick_off_tail(char *str, char ch) {
	int n = 0;
	while ((n = strlen(str) - 1) >= 0 && *(str + n) == ch) {
		*(str + n) = '\0';
	}
}
*/
void kick_off_tail(char *str, char* ch) {
	int n = 0;
	while ((n = strlen(str) - 1) >= 0 && InString(str[n], ch)) {
		*(str + n) = '\0';
	}
}

void pickup_str(char *source, int nstart, int len, char *dest) {
	int n_s = strlen(source);
	if (nstart >= n_s) {strcpy(dest, "\0"); return;}
	int nend = 0;
	if (len < 0) nend = n_s - nstart;
	else nend = nstart + len;
	if (nend >= n_s) nend = n_s;

	//if (nend >= nstart) {strcpy(dest, "\0"); return;}
	for (int i = nstart; i < nend; i++) dest[i-nstart] = source[i];
	dest[nend-nstart] = '\0';
}

bool next_section(char** original, char* ch) {
	if (original == NULL || strlen(*original) == 0) return false;
	int nlen = strlen(ch);
	char item[2560] = "\0";
	char *bf = NULL;
	for (int i = 0; i <= strlen(*original)- nlen; i++) {
		bf = *original + i;
		pickup_str(bf, 0, nlen, item);
		if (strcmp(item, ch) == 0) {
			*original = *original + i;
			return true;
		}
	}
	*original = NULL;
	return false;
}

bool pass_next_section(char** original, char* ch) {
	if (!next_section(original, ch)) return false;
	*original += strlen(ch);
	return true;
}

int read_float(char *buffer, float *v, int max_v, char* format, char seperator) {
	int n = 0;
	char *buff = NULL;
	buff = buffer;
	while (n < max_v && buff != NULL && sscanf(buff, format, v + n) == 1) {
		n++;
		if (!NextSection(&buff, seperator)) break; 
	}
	return n;
}

int read_float(char *buffer, float *v, int max_v, char* format, char* seperator) {
	int n = 0, nsep = strlen(seperator);
	char *buff = NULL;
	buff = buffer;
	while (n < max_v && buff != NULL && sscanf(buff, format, v + n) == 1) {
		n++;
		if (!NextSection(&buff, seperator)) break; 
	}
	return n;
}

int read_float(char *buffer, int ncol, float *v, char *format, char* seperator) {
	int n = 0, nsep = strlen(seperator);
	char *buff = NULL;
	buff = buffer;
	while (n < ncol && buff != NULL) {
		n++;
		if (!NextSection(&buff, seperator)) break; 
	}
	if (buff != NULL && sscanf(buff, format, v) == 1) return 1;
	else return 0;
}

bool read_profile(char *fname, int nIgnoreLines, float **x, int ix, float **y, int iy, int& npts) {
	char buffer[256] = "\0", msg[256] = "\0", seperator[6] = " ,;\t|";
	float vx = 0, vy = 0;
	int i = 0;
	ifstream in;
	in.open(fname); if (!in.is_open()) {sprintf(msg, "failure to open file %s", fname); show_log(msg, true); return false;}
	for (i = 0; i < nIgnoreLines; i++) in.getline(buffer, 250); // ignore some lines
	npts = 0;
	while (!in.eof()) {
		in.getline(buffer, 250); if (strlen(buffer) == 0) continue;
		if (read_float(buffer, ix, &vx, "%f", seperator) == 1 && read_float(buffer, iy, &vy, "%f", seperator) == 1) npts++;
	}
	if (npts == 0) {in.close(); return false;}
	if (*x != NULL) delete[] *x; *x = new float[npts];
	if (*y != NULL) delete[] *y; *y = new float[npts];

	rewind(in);
	for (i = 0; i < nIgnoreLines; i++) in.getline(buffer, 250); // ignore some lines
	i = 0;
	while (!in.eof()) {
		in.getline(buffer, 250); if (strlen(buffer) == 0) continue;
		if (read_float(buffer, ix, &vx, "%f", seperator) == 1 && read_float(buffer, iy, &vy, "%f", seperator) == 1) {
			(*x)[i] = vx; (*y)[i] = vy; i++;
		}
	}
	in.close(); return true;
}

bool read_profile(char *fname, int nIgnoreLines, ARRAY<float*> &va, ARRAY<int> &icol, int& npts) {
	char buffer[256] = "\0", msg[256] = "\0", seperator[6] = " ,;\t|";
	ARRAY<float> v; v.SetArray(va.n);
	int i = 0, ipt;
	bool status;
	ifstream in;
	in.open(fname); if (!in.is_open()) {sprintf(msg, "failure to open file %s", fname); show_log(msg, true); return false;}
	for (i = 0; i < nIgnoreLines; i++) in.getline(buffer, 250); // ignore some lines
	npts = 0;
	while (!in.eof()) {
		in.getline(buffer, 250); if (strlen(buffer) == 0) continue;
		status = true;
		for (i = 0; i < va.n; i++) {
			if (read_float(buffer, icol.m[i], v.m + i, "%f", seperator) != 1) status = false;
		}
		if (status) npts++;
	}
	if (npts == 0) {in.close(); return false;}
	for (i = 0; i < va.n; i++) {
		if (va.m[i] != NULL) delete[] va.m[i]; va.m[i] = new float[npts];
	}

	rewind(in);
	for (i = 0; i < nIgnoreLines; i++) in.getline(buffer, 250); // ignore some lines
	ipt = 0; 
	while (!in.eof()) {
		in.getline(buffer, 250); if (strlen(buffer) == 0) continue;
		status = true;
		for (i = 0; i < va.n; i++) {
			if (read_float(buffer, icol.m[i], v.m + i, "%f", seperator) != 1) status = false;
		}
		if (status) {
			for (i = 0; i < va.n; i++) (va.m[i])[ipt] = v.m[i];
			ipt++;
		}
	}
	in.close(); return true;
}

/*
bool search(ifstream &in, char *title, char *buffer) {
	char buff[2560] = "\0", *bf = "\0";
	strcpy(buffer, "\0"); 
	while (!in.eof()) {
		in.getline(buffer, 2560);
		if (strlen(buffer) == 0) continue;
		bf = buffer; if (bf[0] == ' ') pass_next_section(&bf, " ");
		pickup_str(bf, 0, strlen(title), buff);
		if (strcmp(buff, title) == 0) return true;
	}
	//if (in.eof()) return false;
	return false;
}

bool search(ifstream &in, char *title, int indx, char *buffer) {
	char buff[256] = "\0";
	int n = indx - 1;
	strcpy(buffer, "\0"); 
	while (!in.eof()) {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s %d", buff, &n) == 2 &&
			(strcmp(buff, title) == 0 && n == indx)) return true;
	}
	//if (in.eof()) return false;
	return false;
}
*/

int nItem(char* str, char ch) {
	char *ps = str;
	char buff[256] = "\0";
	int n = 0;
	pass(&ps, ch);
	while (sscanf(ps, "%s", buff) == 1 && strlen(ps) > 0) {
		n++; 
		if (!NextSection(&ps, ch) || ps == NULL) break;
	}
	return n;
}

int nItem(char* str, char* ch) {
	char *ps = str;
	char buff[256] = "\0";
	int n = 0;
	pass(&ps, ch);
	while (sscanf(ps, "%s", buff) == 1 && strlen(ps) > 0) {
		n++; 
		if (!NextSection(&ps, ch) || ps == NULL) break;
	}
	return n;
}


// assuming the macromolecule structure is given
// and reading the torsion angle to reconstruct the new configuration
bool read_mol_MD_TA_struct(ifstream &in, MMOLECULE &mol, int njump, char *errmsg) {
	char buffer[256] = "\0", buff[250] = "\0";
	char TA_title[5] = "[MD]";

	int i = 0;
	float ta = 0;

	i = 0;
	while (!in.eof() && i < njump) {
		strcpy(buffer, "\0"); strcpy(buff, "\0");
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s" , buff) == 1 && strcmp(buff, TA_title) == 0)
			i++;
	}
	if (in.eof()) {sprintf(errmsg, "can not get the line %s", TA_title); return false;}

	i = 0;
	while (!in.eof()) {
		strcpy(buffer, "\0");
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %f", &i, &ta) == 2) {
			mol.twTA(i, ta - mol.cluster[i].km.ta, true);
			mol.cluster[i].km.ta = ta;
		}
		else break;
	}
	return true;
}

bool read_msg(char *fname, char *msg, char *item, bool full_line) {
	char buffer[256] = "\0", buff[250] = "\0";
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(errmsg, "Can not open file %s", fname);
		return false;
	}
	while (strcmp(buff, item) != 0) {
		if (in.eof()) {
			sprintf(errmsg, "Can not find item %s in %s", item, fname);
			in.close(); return false;
		}
		in.getline(buffer, 250);
		sscanf(buffer, "%s", buff);
	}
	strcpy(buffer, "\0");
	in.getline(buffer, 250);
	in.close();
	if (full_line) {strcpy(msg, buffer); return (strlen(msg) > 0 ? true: false);}
	else if (sscanf(buffer, "%s", msg) == 1) return true;
	else {
		sprintf(errmsg, "No message is given after %s in %s", item, fname);
		return false;
	}
}

bool read_MD_TIMESTEP(char *fname, float &dt) {
	char buffer[256] = "\0";
	char item[50] = "\0";
	float v = 0;

	strcpy(item, "[MD_TIME_STEP]");
	if (!read_msg(fname, buffer, item)) return false;
	if (sscanf(buffer, "%f", &v) == 1) {dt = v; return true;}
	else return false;
}

bool read_MD_LOOPS(char *fname, long &nloops) {
	char buffer[256] = "\0";
	char item[50] = "\0";
	long int v = 0;

	strcpy(item, "[TOTAL_LOOPS]");
	if (!read_msg(fname, buffer, item)) return false;
	if (sscanf(buffer, "%ld", &v) == 1) {nloops = (long)v; return true;}
	else return false;
}

bool read_MD_BASE_Pars(char *fname, MD_BASE* mdpar) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(errmsg, "Can not open md parameter file %f", fname);
		return false;
	}
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	float v = 0;
	int n = 0;

	strcpy(item, "[MD_TIME_STEP]");
	while (!in.eof() && strcmp(title, item) != 0) {
		in.getline(buffer, 250);
		sscanf(buffer, "%s", title);
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f", &v) != 1) {
		sprintf(errmsg, "Can not get time step for MD after %s", item);
		in.close(); return false;
	}
	mdpar->set_dt(v);

	strcpy(item, "[TOTAL_LOOPS]");
	while (!in.eof() && strcmp(title, item) != 0) {
		in.getline(buffer, 250);
		sscanf(buffer, "%s", title);
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &n) != 1) {
		sprintf(errmsg, "Can not get total MD loops after %s", item);
		in.close(); return false;
	}
	mdpar->LOOPS = n;

	in.close();
	return true;
}

bool read_MD_SAVE_Pars(char *fname, MD_SAVE &msave) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(errmsg, "Can not open md parameter file %f", fname);
		return false;
	}
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int n = 0, maxEachFile = 100000;

	strcpy(item, "[MD_TITLE]");
	while (!in.eof() && strcmp(title, item) != 0) {
		in.getline(buffer, 250);
		sscanf(buffer, "%s", title);
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%s", title) != 1) {
		sprintf(errmsg, "Can not get title after %s", item);
		in.close();
		return false;
	}
	msave.set_title(title);

	strcpy(item, "[SAVE]");
	while (!in.eof() && strcmp(title, item) != 0) {
		in.getline(buffer, 250);
		sscanf(buffer, "%s", title);
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d", &n, &maxEachFile) != 2) {
		sprintf(errmsg, "Can not get [save every loop] [max. save in each file] after %s", item);
		in.close(); return false;
	}
	msave.set_save(n, maxEachFile);
	
	in.close();
	return true;
}

bool read_HOOVER_Pars(char *fname, float &tao, float &max_ksi) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(errmsg, "Can not open md parameter file %f", fname);
		return false;
	}
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	float v = 0;

	strcpy(item, "[HOOVER_CONST]");
	while (!in.eof() && strcmp(title, item) != 0) {
		in.getline(buffer, 250);
		sscanf(buffer, "%s", title);
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f", &v) != 1) {
		sprintf(errmsg, "Can not get Hoover const. after %s", item);
		in.close(); return false;
	}
	tao = v;

	strcpy(item, "[MAX_KSI_HOOVER]");
	while (!in.eof() && strcmp(title, item) != 0) {
		in.getline(buffer, 250);
		sscanf(buffer, "%s", title);
	}
	if (in.eof()) {
		sprintf(errmsg, "[MAX_KSI_HOOVER] is not defined in %s.\n default maximum value is used : %f", fname, ::max_ksi_Hoover);
		show_msg(errmsg); strcpy(errmsg, "\0");
		in.close(); return true;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f", &v) != 1) {
		sprintf(errmsg, "Can not get maximum ksi value for Hoover Dynamic after %s", item);
		in.close(); return false;
	}
	max_ksi = FABS(v);

	in.close();
	return true;
}

bool readExternalElectricField(char *fname, bool &bEext, float &Eext) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		//sprintf(errmsg, "Can not open md parameter file %f", fname);
		bEext = false; Eext = 0;
		return true;
	}
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	float v = 0;
	int nflag = 0;

	strcpy(item, "[ExternalElectricField]");
	if (!search(in, item, buffer)) {in.close(); bEext = false; Eext = 0; return true;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %f", &nflag, &v) != 2) {
		bEext = false; Eext = 0;
		in.close(); return true;
	}
	else {
		bEext = (nflag == 1 ? true : false);
		if (bEext) Eext = v;
	}

	in.close();
	return true;
}

bool readForceScale(char *fname, bool &scale_force, float &force_max, float &torque_max) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		//sprintf(errmsg, "Can not open md parameter file %f", fname);
		scale_force = false;
		return true;
	}
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[SCALE_FORCE]");
	if (!search(in, item, buffer)) {in.close(); scale_force = false; return true;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %f %f", &nflag, &force_max, &torque_max) != 3) {
		scale_force = false; force_max = 0; torque_max = 0;
		in.close(); return true;
	}
	else {
		scale_force = (nflag == 1 ? true : false);
	}

	in.close();
	return true;
}

bool readDielectricConst(char *fname, float &eps) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		//sprintf(errmsg, "Can not open md parameter file %f", fname);
		eps = 1;
		return true;
	}
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[DIELECTRIC_CONST]");
	if (!search(in, item, buffer)) {
		sprintf(errmsg, "dielectric const is set to 1");
		show_infor(errmsg); strcpy(errmsg, "\0");
		in.close(); eps = 1; return true;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f", &eps) != 1) {
		sprintf(errmsg, "dielectric const is set to 1");
		show_infor(errmsg); strcpy(errmsg, "\0");
		eps = 1; in.close(); return true;
	}
	else {
		if (eps <= 1) {
			eps = 1;
			sprintf(errmsg, "dielectric const is set to 1");
			show_infor(errmsg); strcpy(errmsg, "\0");
		}
		else {
			sprintf(errmsg, "dielectric const is %f", eps);
			show_infor(errmsg); strcpy(errmsg, "\0");
		}
	}

	in.close();
	return true;
}

bool readLocalRelax(char *fname, bool &local_relax, float &dta_max_relax) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		//sprintf(errmsg, "Can not open md parameter file %f", fname);
		local_relax = false;
		return true;
	}
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[LOCAL-RELAX]");
	if (!search(in, item, buffer)) {in.close(); local_relax = false; return true;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %f", &nflag, &dta_max_relax) != 2) {
		local_relax = false;
		in.close(); return true;
	}
	else {
		local_relax = (nflag == 1 ? true : false);
	}

	in.close();
	return true;
}

bool readLocalInteract(char *fname, bool &local_interact) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		//sprintf(errmsg, "Can not open md parameter file %f", fname);
		local_interact = false;
		return true;
	}
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[LOCAL-INTERACT]");
	if (!search(in, item, buffer)) {in.close(); local_interact = false; return true;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nflag) != 1) {
		local_interact = false;
		in.close(); return true;
	}
	else {
		local_interact = (nflag == 1 ? true : false);
	}

	in.close();
	return true;
}

bool readImplicitSolvent(char *fname, bool &bImplicitSolvent) {
	int iflag = 0;
	char buffer[256] = "\0", msg[256] = "\0";
	if (!read_msg(fname, buffer, "[IMPLICIT_SOLVENT]", false) || sscanf(buffer, "%d", &iflag) != 1) {
		sprintf(msg, "[IMPLICIT_SOLVENT] is not defined in %s. Implicit-solvent force is disabled.", fname); show_log(msg, true);
		bImplicitSolvent = false;
	}
	else {
		bImplicitSolvent = (iflag > 0 ? true : false);
		if (bImplicitSolvent) {
			show_log("Implicit-Solvent force will be applied", true);
		}
	}
	return true;
}

bool readInternalKinEnergy(char *fname, float &Ek0_Internal) {
	Ek0_Internal = 0.5f; //in theory is 0.5kT
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return true;

	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	float v = 0;

	strcpy(item, "[INTERNAL_KINETIC_ENERGY]");
	if (!search(in, item, buffer)) {in.close(); return true;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%f", &v) != 1) {in.close(); return true;}
	else {
		if (v > 0) Ek0_Internal = v;
	}

	in.close();
	return true;
}

bool readMacroMolFreeBase(char *fname, bool &free_base) {
	free_base = true;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return true;

	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[FREE_MACROMOLECULE]");
	if (!search(in, item, buffer)) {in.close(); return true;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nflag) != 1) {in.close(); return true;}
	else {
		free_base = (nflag == 1 ? true : false);
	}

	in.close();
	return true;
}

bool readPolarization(char *fname, bool &bPolarize, int &iPolarize) {
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0, imodel = 0;

	bPolarize = false;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(buffer, "failure to open %s. ", fname); show_log(buffer, false);
		show_log("polarization is disabled in default.", true);
		return true;
	}

	strcpy(item, "[POLARIZATION]");
	if (!search(in, item, buffer)) {
		sprintf(buffer, "[POLARIZATION] is not defined in %s. ", fname); show_log(buffer, false);
		show_log("polarization is disabled in default.", true);
		in.close(); return true;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d", &nflag, &imodel) < 1) {
		bPolarize = false; 
		sprintf(title, "[POLARIZATION] is not defined properly in %s, see : ", fname); show_log(title, false); show_log(buffer, true);
		show_log("polarization is disabled in default.", true);
		in.close(); return true;
	}
	else {
		bPolarize = (nflag > 0 ? true : false);
		iPolarize = imodel;
		if (!bPolarize) {
			iPolarize = 0; ::_Thole_polarize_::bTholePolarize = false;
			show_log("Non-polarization!", true);
		}
		else if (iPolarize <= 0) {
			iPolarize = 0; ::_Thole_polarize_::bTholePolarize = false;
			show_log("Polarizable atom as point charge/dipole", true);
		}
		else if (iPolarize == 1 || iPolarize == 2) {
			iPolarize = 1;
			::_Thole_polarize_::bTholePolarize = true;
			::_Thole_polarize_::iTholeModel = iPolarize - 1;
			in.getline(buffer, 250);
			if (sscanf(buffer, "%f", &(::_Thole_polarize_::a)) != 1 || ::_Thole_polarize_::a < 0) {
				if (::_Thole_polarize_::iTholeModel == 0) ::_Thole_polarize_::a = 2.1304f;
				else if (::_Thole_polarize_::iTholeModel == 1) ::_Thole_polarize_::a = 0.572f;
			}
			sprintf(title, "Use Thole-exponential polarization model, with a = %f", ::_Thole_polarize_::a); show_log(title, true);
			show_log("Check spme_interact.cpp and spme_interact_2d.cpp! Right now Thole polarization is used when another atom is trivalent ion!", true);
		}
		else {
			show_log("So far only Thole-exponential polarization model is realized in this program.", true);
			show_log("Since the model is not realized, use point charge/dipole", true);
			iPolarize = 0; ::_Thole_polarize_::bTholePolarize = false;
		}
	}

	in.close();
	return true;
}

bool readPolarization1(char *fname, bool &bExcludeInducedDipole) {
	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0, imodel = 0;

	bExcludeInducedDipole = false;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(buffer, "failure to open %s. ", fname); show_log(buffer, false);
		show_log("interaction between induced dipole is to be included in default.", true);
		return true;
	}

	strcpy(item, "[INDUCED_DIPOLE_INTERACTION]");
	if (!search(in, item, buffer)) {
		sprintf(buffer, "[INDUCED_DIPOLE_INTERACTION] is not defined in %s. ", fname); show_log(buffer, false);
		show_log("interaction between induced dipole is to be included in default.", true);
		in.close(); return true;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d", &nflag, &imodel) < 1) {
		bExcludeInducedDipole = false; 
		sprintf(title, "[INDUCED_DIPOLE_INTERACTION] is not defined properly in %s, see : ", fname); show_log(title, false); show_log(buffer, true);
		show_log("interaction between induced dipole is to be included in default.", true);
		in.close(); return true;
	}
	else {
		bExcludeInducedDipole = (nflag > 0 ? false : true);
		show_log("interaction between induced dipole is to be ", false);
		if (bExcludeInducedDipole) show_log("excluded!", true);
		else show_log("included!", true);
	}

	in.close();
	return true;
}

bool readCoarseGrainMD(char *fname, bool &coarse_grain, char *cgpar_fname, int &nNormalMD, int &nCG) {
	coarse_grain = false;
	nNormalMD = 1; nCG = 1; strcpy(cgpar_fname, "\0");

	ifstream in;
	in.open(fname);
	if (!in.is_open()) return true;

	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[COARSE_GRAIN]");
	if (!search(in, item, buffer)) {in.close(); return true;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %s", &nflag, cgpar_fname) != 2) {in.close(); return true;}
	else {
		coarse_grain = (nflag == 1 ? true : false);
	}
	in.getline(buffer,250);
	if (sscanf(buffer, "%d", &nNormalMD) != 1) {
		sprintf(errmsg, "can not get normal MD step before coarse-grain MD from : %s", buffer);
		show_msg(errmsg); in.close(); return false;
	}
	in.getline(buffer,250);
	if (sscanf(buffer, "%d", &nCG) != 1) {
		sprintf(errmsg, "can not get coarse-grain MD steps from : %s", buffer);
		show_msg(errmsg);in.close(); return false;
	}

	in.close();
	return true;
}

bool read_coarse_grain_md_var(char *fname, float &dt, float &tao, float &max_ksi, int &nCMMRefresh) {
	char title[250] = "\0", msg[256] = "\0";

	strcpy(title, "[CG_MD_TIME_STEP]");
	if (!read_msg(fname, msg, title, true) || sscanf(msg, "%f", &dt) != 1) {
		sprintf(errmsg, "can not get coarse-grain-md-step time in %s from : %s", fname, msg);
		show_msg(errmsg); return false;
	}
	strcpy(title, "[HOOVER_CONST]");
	if (!read_msg(fname, msg, title, true) || sscanf(msg, "%f", &tao) != 1) {
		sprintf(errmsg, "can not get coarse-grain-md tao in %s from : %s", fname, msg);
		show_msg(errmsg); return false;
	}
	strcpy(title, "[MAX_KSI_HOOVER]");
	if (!read_msg(fname, msg, title, true) || sscanf(msg, "%f", &max_ksi) != 1) {
		sprintf(errmsg, "can not get coarse-grain-md max-ksi in %s from : %s", fname, msg);
		show_msg(errmsg); return false;
	}
	strcpy(title, "[CMM_REFRESH]");
	if (!read_msg(fname, msg, title, true) || sscanf(msg, "%d", &nCMMRefresh) != 1) {
		sprintf(errmsg, "can not get cmm-checking md steps in %s from : %s", fname, msg);
		show_msg(errmsg); return false;
	}
	return true;
}

bool readBrownianPars(char *fname, bool &Brownian, float &ita_Bronian, float &cd_trans, float &cd_rot) {
	Brownian = false;

	ifstream in;
	in.open(fname);
	if (!in.is_open()) return true;

	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[BROWNIAN_FORCE]");
	if (!search(in, item, buffer)) {in.close(); return true;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %f %f %f", &nflag, &ita_Bronian, &cd_trans, &cd_rot) != 4) {in.close(); return true;}
	else {
		Brownian = (nflag == 1 ? true : false);
	}

	in.close();
	return true;
}

bool readSPMEPars(char *fname, int &indx_bspline, int &dim_spme_x, int &dim_spme_y, int &dim_spme_z, bool &bExpandSPMECellZ) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return true;

	char buffer[256] = "\0", title[250] = "\0";
	char *bf = buffer;
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[SPME]");
	if (!search(in, item, buffer)) {in.close(); return true;}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d %d %d", &indx_bspline, &dim_spme_x, &dim_spme_y, &dim_spme_z) != 4) {
		sprintf(errmsg, "SPME parameters are not given properly [B-Spline func. indx], mesh dimension [x, y, z]"); show_msg(errmsg);
		in.close(); return false;
	}
	else {
		bf = buffer;
		NextSection(&bf); NextSection(&bf); NextSection(&bf); NextSection(&bf);
		if (bf != NULL && sscanf(bf, "%d", &nflag) == 1 && nflag == 1) bExpandSPMECellZ = true;
		else bExpandSPMECellZ = false;
		if (indx_bspline / 2 * 2 != indx_bspline) {
			sprintf(errmsg, "B-Spline Func. indx has to be even"); show_msg(errmsg);
			in.close(); return false;
		}
	}

	in.close();
	return true;
}

bool readNPTPars(char *fname, float& W_npt, float &tao_npt) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return true;

	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[NPT]");
	if (!search(in, item, buffer)) {
		sprintf(errmsg, "[NPT] is not defined in %s", fname); show_msg(errmsg);
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if ((nflag = sscanf(buffer, "%f %f", &W_npt, &tao_npt)) == 1) {
		tao_npt = W_npt;
		sprintf(errmsg, "tao of barostat is not given, use the same value for the weight, %f fs", tao_npt); show_log(errmsg, true);
	}
	else if (nflag == 0) {
		sprintf(errmsg, "Weight of barostat is not given in %s following [NPT]", buffer); show_msg(errmsg);
		return false;
	}

	in.close();
	return true;
}

bool readNPT_RF_Exclude(char *fname, bool& bNPT_RF_Exclude, int &nmax, float &Fmin) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return true;

	char buffer[256] = "\0", title[250] = "\0";
	char item[50] = "\0";
	int nflag = 0;

	strcpy(item, "[NPT_RF_EXCLUDE]");
	if (!search(in, item, buffer)) {
		sprintf(errmsg, "[NPT_RF_EXCLUDE] is not defined in %s. Disable it!", fname); show_log(errmsg, true);
		bNPT_RF_Exclude = false;
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d %d %f", &nflag, &nmax, &Fmin) == 3 && nflag > 0) {
		bNPT_RF_Exclude = true;//indepedent Barostat 
	}
	else {
		bNPT_RF_Exclude = false;
		if (nflag > 0) {
			sprintf(errmsg, "NPT_RF_EXCLUDE is not given properly [flag, nmax_exclude, Fmin], see %s", buffer); show_log(errmsg, true);
			show_log("NPT_RF_EXCLUDE is disabled", true);
		}
	}

	in.close();
	return true;
}


bool readSurfaceBoundary(char *fname) {
	char buffer[256] = "\0", msg[256] = "\0";
	if (!read_item(fname, "[SURFACE_DIPOLE]", buffer)) {
		show_log("[SURFACE_DIPOLE] boundary condition is not set in ", false); show_log(fname, true);
		show_log("enable surface term in Ewarld-Sum !", true); iSurfaceBoundary = 0; return false;
	}
	int flag = 0;
	if (sscanf(buffer, "%d", &flag) != 1 || flag == 0) {
		show_log("enable surface term in Ewarld-Sum !", true); iSurfaceBoundary = 0; return true;
	}
	else if (flag == 1) {
		show_log("using tinfoil boundary condition, disable surface term in Ewarld-Sum !", true); iSurfaceBoundary = flag; return true;
	}
	else if (flag == 2) {
		show_log("using tinfoil boundary condition along z direction, disable surface term in Ewarld-Sum along z axis!", true); iSurfaceBoundary = flag; return true;
	}
	else {
		show_log("unknown surface model. enable surface term in Ewarld-Sum !", true); iSurfaceBoundary = 0; return true;
	}
}

bool readNetForceCorrectionMode(char *fname) {
	char buffer[256] = "\0", msg[256] = "\0";
	if (!read_item(fname, "[NETFORCE_CORRECTION_MODE]", buffer)) {
		::mode_NetForce = 0; show_log("[NETFORCE_CORRECTION_MODE] is not given, disable net-force-correction.", true); return true;
	}
	int flag = 0;
	if (sscanf(buffer, "%d", &flag) == 1) {
		::mode_NetForce = flag; 
		switch (::mode_NetForce) {
		case 1:
			show_log("correct net-force along x/y/z axises.", true); return true;
			break;
		case 2:
			show_log("correct net-force along z-axis.", true); return true;
			break;
		default:
			::mode_NetForce = 0;
			show_log("[NETFORCE_CORRECTION_MODE] is not given, disable net-force-correction.", true); return true;
		}
		return true;
	}
	else {
		::mode_NetForce = 0;
		show_log("[NETFORCE_CORRECTION_MODE] value is not expected, disable net-force-correction.", true); return true;
	}
}
