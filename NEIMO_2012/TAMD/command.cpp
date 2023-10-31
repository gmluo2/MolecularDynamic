#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
#include "show.h"
#include "ranlib.h"
#include "vector.h"
#include "bound.h"
#include "MM.h"
#include "Mol.h"
#include "EwaldSumVar.h"
#include "CMM.h"
#include "var.h"

#include "distribution.h"

#include "command.h"

extern LIST<MMOLECULE> *mlist;

COMMAND* pCmd = NULL;
COMMAND* pCurrentCmd = NULL;
int nCurMol = -1;
char cmd_msg[256] = "\0";

extern void show_all_molecules();
extern bool search(ifstream &in, char *title, char *buffer);

#if _SYS_ == _WINDOWS_SYS_
extern CWnd *pMainFrame;
#endif

using namespace _distribution_;
using namespace _cluster_distribution_;

void UpperString(char* str) {
	int len = int(strlen(str));
	for (int i = 0; i < len; i++) {
		if (*(str + i) <= 'z' && *(str + i) >= 'a') *(str + i) =  *(str + i) + 'A' - 'a';
	}
}

COMMAND* AddCommand(char* c1, char* c2, bool (*func)())
{
	COMMAND *p = new COMMAND;
	strcpy(p->cmd1, c1);
	strcpy(p->cmd2, c2);
	p->func = func;

	COMMAND *tmp = pCmd;
	if (pCmd == NULL) pCmd = p;
	else {
		while(tmp->next != NULL) tmp = tmp->next;
		tmp->next = p;
	}
	return p;
}

void DefCommands()
{
	pCmd = NULL;
	COMMAND* cmd = NULL;
	char help_string[1024] = "\0";
	char format[256] = "\0";

	cmd = AddCommand("HELP", "", ShowCommand);
	cmd->Help_String("HELP", "");
	
	cmd = AddCommand("H", "", ShowCommand);
	cmd->Help_String("H", "");

	AddCommand("?", "", ShowCommand);
	cmd->Help_String("?", "");

	cmd = AddCommand("SHOW", "", show_cur_mol);
	cmd->Help_String("show", "show current molecule");

	cmd = AddCommand("TWEAK", "MOL", tweak_mol);
	cmd->Help_String("tweak mol [wx] [wy] [wz]", "tweak current molecule at its base with [wx] [wy] [wz]");

	cmd = AddCommand("TWEAK", "", tweak);
	cmd->Help_String("tweak", "tweak cluster [n] of current molecule along its hinge");

	cmd = AddCommand("SHIFT", "", shift);
	cmd->Help_String("shift [dx] [dy] [dz]", "shift current molecule with [dx] [dy] [dz]");

	cmd = AddCommand("DELETE", "", delete_cur_mol);
	cmd->Help_String("delete", "delete current molecule");

	cmd = AddCommand("COPY", "", copy_cur_mol);
	cmd->Help_String("copy", "copy current molecule to a new molecule");

	cmd = AddCommand("SAVE", "", save_cur_mol);
	cmd->Help_String("save [ofname]", "save current molecule to a file");

	cmd = AddCommand("USE", "", set_cur_mol);
	cmd->Help_String("USE [n]", "set molecule #n as current molecule");

	cmd = AddCommand("HIGHLIGHT", "", highlight);
	cmd->Help_String("HIGHLIGHT [n]", "highlight cluster #n in current molecule");

	cmd = AddCommand("READ", "CHARGE", read_charge_mol2);
	cmd->Help_String("read charge [filename]", "read atom charges from a file calculated with antechamber, in mol2 format");
	
	cmd = AddCommand("SHOW", "CHARGE", show_charge_cluster);
	cmd->Help_String("show charge [n]", "show atomic charges in cluster n ");
	
	cmd = AddCommand("AVERAGE", "CHARGE", average_charge_cluster);
	cmd->Help_String("average charge [n]", "average atomic charges in same kind of clusters as cluster #n ");

	cmd = AddCommand("MOL", "DIPOLE", total_dipole);
	cmd->Help_String("mol dipole", "calculate the total net charge & dipole of the macro-molecule");

	cmd = AddCommand("MONOMER", "DIPOLE", monomer_dipole);
	cmd->Help_String("monomer dipole [n]", "calculate the net charge & dipole of the monomer in macro-molecule");

	cmd = AddCommand("DIPOLE", "DISTRIBUTION", monomer_dipole_distribution);
	cmd->Help_String("dipole distribution [mode] [monomer] [from] [to] [npts] [ofname]", "ananyze the monomer's dipole moment distribution");

	cmd = AddCommand("CLUSTER", "RADIUS", cal_cluster_radius);
	cmd->Help_String("cluster radius [n]", "calculate the radius of cluster (atom center)");
}

void DelCommands() {
	COMMAND *p = ::pCmd, *next = NULL;
	while (p != NULL) {
		next = p->next;
		p->next = NULL;
		delete p;
		p = next;
	}
	::pCmd = NULL;
}

bool ShowCommand()
{
	char msg[4800] = "\0";
	char buff[256] = "\0";
	sprintf(buff, "? / h /help -- show commands\n"); strcat(msg, buff);
	sprintf(buff, "show -- show current molecule #\n"); strcat (msg, buff);
	sprintf(buff, "tweak mol [dx] [dy] [dz]\n"); strcat(msg, buff);
	sprintf(buff, "tweak [nCluster] [domega] -- tweak cluster and its children along its parent hinge\n"); strcat(msg, buff);
	sprintf(buff, "shift [dx] [dy] [dz] -- shift molecule \n"); strcat(msg, buff);
	sprintf(buff, "use [n] -- set molecule #n as the one working at\n"); strcat(msg, buff);
	sprintf(buff, "save [output filename] -- save current molecule struct. to a file: torsion-angles\n"); strcat(msg, buff);
	sprintf(buff, "copy -- make a copy of current molecule and set it as current one\n"); strcat(msg, buff);
	sprintf(buff, "delete -- delete current molecule\n"); strcat(msg, buff);
	sprintf(buff, "highlight [n] -- highlight the cluster in current molecule\n"); strcat(msg, buff);
	sprintf(buff, "read charge [ifname] -- read atom charges from a file calculated with antechamber, in mol2 format\n"); strcat(msg, buff);
	sprintf(buff, "show charge [n] -- show atom charges in cluster #n\n"); strcat(msg, buff);
	sprintf(buff, "average charge [n] -- average atomic charges in same kind of clusters as cluster #n \n"); strcat(msg, buff);
	sprintf(buff, "mol dipole -- calculate the total net charge & dipole of the macro-molecule \n"); strcat(msg, buff);
	sprintf(buff, "monomer dipole [n] -- calculate the net charge & dipole of the monomer in macro-molecule \n"); strcat(msg, buff);
	sprintf(buff, "dipole distribution [mode] [from] [to] [npts] [ofname] \n       -- ananyze the monomer's dipole moment distribution \n"); strcat(msg, buff);
	sprintf(buff, "       [mode] -- 0 : |mu|; 1 : mu_x; 2 : mu_y; 3 : mu_z"); strcat(msg, buff);
	sprintf(buff, "cluster radius [n] -- calculate radius of cluster n \n"); strcat(msg, buff);

	show_msg(msg);
	return true;
}

COMMAND* GetCommand(char* c1, char* c2)
{
	COMMAND *p = pCmd;
	bool getcommand = false;
	char cmd1[256] = "\0", cmd2[256] = "\0";
	strcpy(cmd1, c1); UpperString(cmd1);
	strcpy(cmd2, c2); UpperString(cmd2);
	while (p != NULL) {
		if (strcmp(p->cmd1, cmd1) == 0) {
			if (strcmp(p->cmd2, cmd2) == 0) {getcommand = true; break;}
			else if (strcmp(p->cmd2, "\0") == 0 || strlen(p->cmd2) == 0) {
				getcommand= true;
				// do not break, there could be other command match it
				p = p->next;
			}
			else p = p->next;
		}
		else p = p->next;
	}
	if (getcommand) {::pCurrentCmd = p; return p;}
	else return NULL;
}

bool get_cur_mol(MMOLECULE **mm) {
	if (mlist == NULL || nCurMol < 0) return false;
	LIST<MMOLECULE> *ml = mlist->get_list(nCurMol);
	if (ml == NULL || ml->p == NULL) {
		sprintf(errmsg, "molecule # %d is not defined", nCurMol);
		show_msg(errmsg); strcpy(errmsg, "\0"); return false;
	}
	MMOLECULE *m = ml->p;
	if (m->nCluster <= 0) {
		sprintf(errmsg, "no cluster is defined in molecule # %d", nCurMol);
		show_msg(errmsg); strcpy(errmsg, "\0"); return false;
	}
	if (m->base == NULL) {
		sprintf(errmsg, "base cluster is not assigned in molecule # %d", nCurMol);
		show_msg(errmsg); strcpy(errmsg, "\0"); return false;
	}
	*mm = m;
	return true;
}

bool tweak_mol() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0", cmd2[20] = "\0";
	float dx = 0, dy = 0, dz = 0;
	if (sscanf(cmd_msg, "%s %s %f %f %f", cmd1, cmd2, &dx, &dy, &dz) != 5) return false;
	VECTOR3 rot;
	rot.v[0] = dx * Angle2Radian; rot.v[1] = dy * Angle2Radian; rot.v[2] = dz * Angle2Radian;

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;
	m->twMol(*(m->base->Op), rot, true);

	show_all_molecules();
	return true;
}

bool tweak() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0", cmd2[20] = "\0";
	float domega = 0;
	int nCluster = 0;
	if (sscanf(cmd_msg, "%s %d %f", cmd1, &nCluster, &domega) != 3) return false;
	domega *= fAngle2Radian;

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;
	m->twTA(nCluster, domega, true); // in degree
	m->cluster[nCluster].km.ta += domega;
	PERIOD_RANGE(m->cluster[nCluster].km.ta, -PI, PI, PI2)

	show_all_molecules();
	return true;
}

bool shift() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0", cmd2[20] = "\0";
	float dx[3] = {0, 0, 0};
	if (sscanf(cmd_msg, "%s %f %f %f", cmd1, dx, dx + 1, dx + 2) != 4) return false;
	VECTOR3 ds(dx[0], dx[1], dx[2]);

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;
	m->shiftMol(ds);

	show_all_molecules();
	return true;
}

bool show_cur_mol() {
	char buffer[256] = "\0";
	int nTotal = 0;
	if (mlist != NULL) nTotal = mlist->nList(true);
	sprintf(buffer, "%d molecules; current molecule # %d", nTotal, nCurMol);
	show_msg(buffer);
	return true;
}

bool copy_cur_mol() {
	if (strlen(cmd_msg) == 0) return true;
	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;
	MMOLECULE *dm = new MMOLECULE;
	if (!cp(dm, m, true)) { // molecule will be copied, the tree will be setup
		sprintf(errmsg, "failure to copy molecule"); show_msg(errmsg); strcpy(errmsg, "\0"); return false;
	}
	
	LIST<MMOLECULE> *plist = new LIST<MMOLECULE>;
	plist->p = dm;
	int indx = mlist->indx(plist);
#if _SYS_ == _WINDOWS_SYS
	if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_ADD_MOLECULE, 0, (LPARAM)indx);
#endif
	return true;
}

bool delete_cur_mol() {
	if (strlen(cmd_msg) == 0) return true;
	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;

	int indx = mlist->indx(m);
	detach(&mlist, m, true);
	if (mlist != NULL) nCurMol = 0;
	else nCurMol = -1;

#if _SYS_ == _WINDOWS_SYS
	if (pMainFrame != NULL) ::SendMessage(pMainFrame->m_hWnd, WM_DELETE_MOLECULE, 0, (LPARAM)indx);	
#endif
	return true;
}

bool set_cur_mol() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0";
	int n = 0;
	if (sscanf(cmd_msg, "%s %d", cmd1, &n) != 2) return false;
	int nTotal = 0;
	if (mlist != NULL) nTotal = mlist->nList(false);
	if (n < nTotal) nCurMol = n;
	else if (nTotal <= 0) nCurMol = -1;
	else nCurMol = nTotal - 1;
	return true;
}

bool save_cur_mol() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd[20] = "\0", ofname[256] = "\0";
	if (sscanf(cmd_msg, "%s %s", cmd, ofname) != 2) return false;

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;
	ofstream out;
	out.open(ofname);
	if (!out.is_open()) {
		sprintf(errmsg, "Can not open file %s to save molecule structure", ofname);
		show_msg(errmsg);
		return true;
	}
	CLUSTER *pc = m->cluster;
	out<<m->base->Op->v[0]<<"  "<<m->base->Op->v[1]<<"  "<<m->base->Op->v[2];
	out<<"  "<<m->zaxis.v[0]<<"  "<<m->zaxis.v[1]<<"  "<<m->zaxis.v[2];
	out<<"  "<<m->xaxis.v[0]<<"  "<<m->xaxis.v[1]<<"  "<<m->xaxis.v[2]<<endl;
	for (int n = 0; n < m->nCluster; n++) {
		pc = m->cluster + n;
		out<<n<<"  "<<pc->km.ta<<endl;
	}
	out.close();
	return true;
}

bool highlight() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd[20] = "\0";
	int n = 0;
	if (sscanf(cmd_msg, "%s %d", cmd, &n) != 2) return false;
	
	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;

	if (n >= m->nCluster) {
		sprintf(errmsg, "current molecule has %d cluster only", m->nCluster); show_msg(errmsg);
		return false;
	}
	m->highlight(false);
	m->cluster[n].highlight = true;

	show_all_molecules();
	return true;
}

bool read_antech_atom(char *buffer, char *amber_atom, float &charge) {
	int n = 0;
	float c= 0;
	char amber_alias[50] = "\0";
	if (strlen(buffer) < 75 || sscanf(buffer + 47, "%s", amber_alias) != 1 || sscanf(buffer + 67, "%f", &c) != 1) {
		sprintf(errmsg, "this is not generated by antechamber : %s", buffer);
		show_msg(errmsg);
		return false;
	}
	else {
		if (strcmp(amber_alias, amber_atom) != 0) {
			sprintf(errmsg, "atom given from antechamber : %s, while in molecule it is %s\n %s", amber_alias, amber_atom, buffer);
			show_msg(errmsg);
		}
		charge = c;
		return true;
	}
}

bool read_antech_atoms(char *fname, MMOLECULE &m) {
	ifstream in;
	char buffer[256] = "\0";
	in.open(fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", fname); show_msg(errmsg); return false;}
	if (!search(in, "@<TRIPOS>ATOM", buffer)) {
		sprintf(errmsg, "can not find defined atoms in %s", fname); show_msg(errmsg); in.close(); return false;
	}
	int nc = 0, na;
	float c = 0;
	MATOM *patom = NULL;
	CLUSTER *pc = NULL;
	for (nc = 0; nc < m.nCluster; nc++) {
		pc = m.cluster + nc;
		for (na = 0; na < pc->nAtoms; na++) {
			if (in.eof()) {
				sprintf(errmsg, "%s does not match to this macro-molecule", fname);
				show_msg(errmsg); in.close(); return false;
			}
			patom = pc->atom + na;
			in.getline(buffer, 250);
			if (!read_antech_atom(buffer, patom->par->atom, c)) {
				sprintf(errmsg, "%s does not match to this macro-molecule", fname);
				show_msg(errmsg); in.close(); return false;
			}
			else {patom->c = c; patom->c0 = c;}
		}
	}
	in.close();
	return true;
}

bool read_charge_mol2() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0", cmd2[20] = "\0", ifname[256] = "\0";
	if (sscanf(cmd_msg, "%s %s %s", cmd1, cmd2, ifname) != 3) return false;

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;

	bool status = read_antech_atoms(ifname, *m);

	return status;
}

bool show_charge_cluster() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0", cmd2[20] = "\0";
	int nc = 0;
	if (sscanf(cmd_msg, "%s %s %d", cmd1, cmd2, &nc) != 3) return false;

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;
	if (nc < 0) nc = m->nCluster + nc;

	if (nc < 0 || nc >= m->nCluster) {
		sprintf(errmsg, "macro-molecule has only %d clusters", m->nCluster);
		show_msg(errmsg); return false;
	}

	// highlight the cluster
	m->highlight(false);
	m->cluster[nc].highlight = true;
	show_all_molecules();

	char buffer[1024] = "\0", buff[250] = "\0";
	CLUSTER *pc = m->cluster + nc;
	MATOM *patom = NULL;
	int na = 0;
	float c_total = 0;
	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		sprintf(buff, "%6s : %8.4f\n", patom->par->atom, patom->c0);
		strcat(buffer, buff);
		c_total += patom->c0;
	}
	sprintf(buff, "net charge : %8.4f\n", c_total); strcat(buffer, buff);
	show_msg(buffer);
	return true;
}

bool average_charge_cluster() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0", cmd2[20] = "\0";
	int nc = 0;
	if (sscanf(cmd_msg, "%s %s %d", cmd1, cmd2, &nc) != 3) return false;

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;
	if (nc < 0) nc = m->nCluster + nc;

	if (nc < 0 || nc >= m->nCluster) {
		sprintf(errmsg, "macro-molecule has only %d clusters", m->nCluster);
		show_msg(errmsg); return false;
	}

	int cindx = m->cluster[nc].cID;

	m->highlight(false); // disable highlight 
	
	float *c = new float[m->cluster[nc].nAtoms];

	char buffer[1024] = "\0", buff[250] = "\0";
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	int na = 0, nc1 = 0, nSameCluster = 0;
	pc = m->cluster + nc;
	for (na = 0; na < pc->nAtoms; na++) c[na] = 0;
	for (nc1 = 0; nc1 < m->nCluster; nc1++) {
		pc = m->cluster + nc1;
		if (pc->cID != cindx) continue;
		nSameCluster++; // one more cluster
		pc->highlight = true; //highlight this cluster
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			c[na] += patom->c0;
		}
	}

	float c_total = 0;
	pc = m->cluster + nc;
	for (na = 0; na < pc->nAtoms; na++) {
		patom = pc->atom + na;
		c[na] /= nSameCluster;
		sprintf(buff, "%6s : %8.4f\n", patom->par->atom, c[na]);
		strcat(buffer, buff);
		c_total += c[na];
	}
	sprintf(buff, "net charge : %8.4f\n", c_total); strcat(buffer, buff);

	delete[] c; c = NULL;
	show_all_molecules(); // refresh view

	show_msg(buffer);
	return true;
}

bool total_dipole() {
	if (strlen(cmd_msg) == 0) return true;
	//char cmd1[20] = "\0", cmd2[20] = "\0";
	int nc = 0;
	//if (sscanf(cmd_msg, "%s %s %d", cmd1, cmd2, &nc) != 3) return false;

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;
/*
	if (nc < 0) nc = m->nCluster + nc;
	if (nc < 0 || nc >= m->nCluster) {
		sprintf(errmsg, "macro-molecule has only %d clusters", m->nCluster);
		show_msg(errmsg); return false;
	}
*/
	int cindx = m->cluster[nc].cID;

	m->highlight(false); // disable highlight 
	
	char buffer[1024] = "\0", buff[250] = "\0";
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	int na = 0, nc1 = 0, nSameCluster = 0, i;

	VECTOR3 mu;
	float c_total = 0;

	for (nc = 0; nc < m->nCluster; nc++) {
		pc = m->cluster + nc;
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			for (i = 0; i < 3; i++) mu.v[i] += patom->c0 * patom->r.v[i];
			c_total += patom->c0;
		}
	}

	sprintf(buff, "net charge : %8.4f; \ndipole : [%8.4f, %8.4f, %8.4f]", c_total, mu.v[0], mu.v[1], mu.v[2]); strcat(buffer, buff);

	show_all_molecules(); // refresh view

	show_msg(buffer);
	return true;
}

bool monomer_dipole() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0", cmd2[20] = "\0";
	int n_monomer = 0;
	if (sscanf(cmd_msg, "%s %s %d", cmd1, cmd2, &n_monomer) != 3) return false;

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;

	if (n_monomer < 0) n_monomer = m->nMonomer + n_monomer;
	if (n_monomer < 0 || n_monomer >= m->nMonomer) {
		sprintf(errmsg, "macro-molecule has only %d clusters", m->nCluster);
		show_msg(errmsg); return false;
	}

	m->highlight(false); // disable highlight 
	
	char buffer[1024] = "\0", buff[250] = "\0";
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	int na = 0, i, nc;

	VECTOR3 mu;
	float c_total = 0;

	_GROUP<CLUSTER> *monomer = m->monomer + n_monomer;

	for (nc = 0; nc < monomer->nUnit; nc++) {
		pc = monomer->u[nc];
		pc->highlight = true;
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			for (i = 0; i < 3; i++) mu.v[i] += patom->c0 * patom->r.v[i];
			c_total += patom->c0;
		}
	}

	float mu_abs = mu.abs();

	sprintf(buff, "net charge : %8.4f; \ndipole : [%8.4f, %8.4f, %8.4f]; \n|mu| = %8.4f", c_total, mu.v[0], mu.v[1], mu.v[2], mu_abs); strcat(buffer, buff);

	show_all_molecules(); // refresh view

	show_msg(buffer);
	return true;
}

extern bool get_monomer_uid(char *monomer, int &uid);

bool monomer_dipole_distribution() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0", cmd2[20] = "\0", title[256] = "\0", ofname[256] = "\0", monomer[256] = "\0";
	float xf = 0, xt = 0;
	int npts = 0, mode = 0, monomer_uid = 0;
	if (sscanf(cmd_msg, "%s %s %s %f %f %d %s", cmd1, cmd2, monomer, &xf, &xt, &npts, title) != 8) return false;
	if (npts < 0 || xf >= xt) {
		sprintf(errmsg, "from < to, npts > 0"); show_msg(errmsg);
		return false;
	}
	if (!get_monomer_uid(monomer, monomer_uid)) {
		sprintf(errmsg, "unknown monomer %s", monomer); show_msg(errmsg); return false;
	}

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;

	char buffer[1024] = "\0", buff[250] = "\0";
	int n = 0;

	int i = 0;
	DISTRIBT<float, long> dt[N_MU_DT];
	dt[MU_ABS].set_distribt_x(xf, xt, npts);
	dt[MU_X].set_distribt_x(-xt, xt, npts + npts);
	dt[MU_Y].set_distribt_x(-xt, xt, npts + npts);
	dt[MU_Z].set_distribt_x(-xt, xt, npts + npts);
	dt[MU_XY].set_distribt_x(0, xt, npts + npts);
	dt[MU_THETA].set_distribt_x(0, fPI, 90);
	dt[MU_OMEGA].set_distribt_x(0, fPI + fPI, 180);
	monomer_dipole_distribution(*m, false, monomer_uid, true, xf, xt, dt);

	dt[MU_ABS].eb_raw();
	ofstream *out = NULL;
	sprintf(ofname, "%s-mu_abs.dat", title);
	out->open(ofname);
	dt[MU_ABS].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;

	dt[MU_X].eb_raw();
	sprintf(ofname, "%s-mu_x.dat", title);
	out->open(ofname);
	dt[MU_X].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;

	dt[MU_Y].eb_raw();
	sprintf(ofname, "%s-mu_y.dat", title);
	out->open(ofname);
	dt[MU_Y].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;

	dt[MU_Z].eb_raw();
	sprintf(ofname, "%s-mu_z.dat", title);
	out->open(ofname);
	dt[MU_Z].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;

	dt[MU_XY].eb_raw();
	sprintf(ofname, "%s-mu_xy.dat", title);
	out->open(ofname);
	dt[MU_XY].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;

	dt[MU_THETA].eb_raw();
	sprintf(ofname, "%s-mu_theta.dat", title);
	out->open(ofname);
	dt[MU_THETA].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;

	dt[MU_OMEGA].eb_raw();
	sprintf(ofname, "%s-mu_omega.dat", title);
	out->open(ofname);
	dt[MU_OMEGA].save_x_y_eb(*out);
	out->close(); delete out; out = NULL;

	sprintf(errmsg, "distribution is saved in files titled with %s", title);
	show_msg(errmsg);
	return true;
}

bool cal_cluster_radius() {
	if (strlen(cmd_msg) == 0) return true;
	char cmd1[20] = "\0", cmd2[20] = "\0";
	int nc = 0;
	if (sscanf(cmd_msg, "%s %s %d", cmd1, cmd2, &nc) != 3) return false;

	MMOLECULE *m = NULL;
	if (!get_cur_mol(&m)) return true;
	if (nc < 0) nc = m->nCluster + nc;

	if (nc < 0 || nc >= m->nCluster) {
		sprintf(errmsg, "macro-molecule has only %d clusters", m->nCluster);
		show_msg(errmsg); return false;
	}

	// highlight the cluster
	m->highlight(false);
	m->cluster[nc].highlight = true;
	show_all_molecules();

	char buffer[1024] = "\0", buff[250] = "\0";
	CLUSTER *pc = m->cluster + nc;

	float radius = pc->cal_geometry_radius();
	
	sprintf(buff, "radius : %8.4f\n", radius); strcat(buffer, buff);
	show_msg(buffer);
	return true;
}
