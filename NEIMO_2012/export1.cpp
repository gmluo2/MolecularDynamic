#include "TAMD\project.h"
//#include "project.h"

#if _SYS_ == _WINDOWS_SYS_
#include "stdafx.h"
#endif

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#define _READ_MM_MD_ 1

using namespace std;
#if _SYS_ == _WINDOWS_SYS_
#include "TAMD\def.h"
#include "TAMD\show.h"
#include "TAMD\ranlib.h"
#include "TAMD\vector.h"
#include "TAMD\bound.h"
#include "TAMD\Matrix.h"
#include "TAMD\nhc.h"

#include "TAMD\EwaldSumVar.h"

#include "TAMD\MM.h"
#include "TAMD\Mol.h"
#include "TAMD\CMM.h"
#include "TAMD\cluster.h"
#include "TAMD\MD.h"
#include "TAMD\cg-mm.h"
#include "TAMD\cg-md.h"
#include "TAMD\cg-cmm.h"
#include "TAMD\cg-cluster.h"
#include "TAMD\Interaction.h"
#include "TAMD\Interact1.h"
#include "TAMD\md_cell.h"
#include "TAMD\cg-md-cell.h"
#include "TAMD\cg-mpi-md-cell.h"
#include "TAMD\test.h"
#include "TAMD\read.h"
#include "TAMD\var.h"
#include "TAMD\ReadMDProc.h"
#include "TAMD\mol-cell.h"
#include "TAMD\show_md.h"
#include "export1.h"

#include "TAMD\vmd_cell.h"

#elif _SYS_ == _LINUX_SYS_
#include "def.h"
#include "show.h"
#include "ranlib.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "cluster.h"
#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"
#include "cg-cmm.h"
#include "cg-cluster.h"
#include "Interaction.h"
#include "Interact1.h"
#include "md_cell.h"
#include "cg-md-cell.h"
#include "cg-mpi-md-cell.h"
#include "test.h"
#include "read.h"
#include "var.h"
#include "ReadMDProc.h"
#include "mol-cell.h"
#include "show_md.h"
#include "export1.h"

#include "vmd_cell.h"

#endif  // _SYS_

using namespace _coarse_grain_;

extern MMOLECULE mm;
extern CG_MMOLECULE cgmm;
extern MMOL_MD_CELL mm_md_cell;
extern FLEX_MMOL_CELL flex_mmcell;
extern MMOL_VMD_CELL vmdcell;
extern CGMM_MD_CELL cgmm_md_cell;
extern LIST<MMOLECULE> *mlist;
extern LIST<CG_MMOLECULE> *cgmlist;

extern int JOB;
extern int MD_MODE;

extern bool md_run;

#if _SYS_ == _WINDOWS_SYS_
extern CWnd *pMainFrame;
extern HANDLE hWaitGetMolStructEvent;
extern HANDLE hWaitEvent;
#endif

#define _X_RAY  0
#define _NEUTRON 1
int LIGHT_SOURCE = _X_RAY;
PICK_SCATT_STRUCT pick_struct;

#define UID_H   1

extern bool read_item(char *fname, char *title, char *item);

extern bool get_molecule_upto_job(int nMol, GENERAL_MOLECULE& gm); // defined in export.cpp
extern int get_macromol_atom_number(MMOLECULE *mm); // defined in export.cpp
extern int get_macromol_atom_number(CG_MMOLECULE *cgmm);
extern int get_total_atom_number(); // defined in export.cpp

extern bool ReadCGMMStruct(ifstream &in, CG_MMOLECULE &mm);
extern bool ReadCGSMStruct(ifstream &in, CGSM &sm);

extern "C" void set_scatt(int light_source) {
	::LIGHT_SOURCE = light_source;
}

int get_macromol_scatt_atom_number(MMOLECULE *mm, bool all_atom) {
	if (all_atom) return get_macromol_atom_number(mm);

	int n = 0, nc = 0, na = 0;
	CLUSTER *pc = NULL;
	for (nc = 0; nc < mm->nCluster; nc++) {
		pc = mm->cluster + nc;
		for (na = 0; na < pc->nAtoms; na++) {
			if (all_atom) n++;
			else if (pc->atom[na].par->uid != UID_H) n++; // 1 is H
		}
	}
	return n;
}

int get_simple_molecule_scatt_atom_number(SMOLECULE *sm, bool all_atom) {
	if (all_atom) return sm->c->nAtoms;

	int n = 0, na = 0;
	BASIC_CLUSTER *pc = sm->c;
	for (na = 0; na < pc->nAtoms; na++) {
		if (all_atom) n++;
		else if (pc->atom[na].par->uid != UID_H) n++; // 1 is H
	}
	return n;
}

int get_point_molecule_scatt_atom_number(PMOLECULE *pm, bool all_atom) {
	if (all_atom) return pm->c->nAtoms;

	int n = 0, na = 0;
	BASIC_CLUSTER *pc = pm->c;
	for (na = 0; na < pc->nAtoms; na++) {
		if (all_atom) n++;
		else if (pc->atom[na].par->uid != UID_H) n++; // 1 is H
	}
	return n;
}

// get the atoms in current job
extern "C" int get_total_scatt_atom_number(bool all_atom) {
	// specific type of monomer ?
	if (pick_struct.monomer_uid.n > 0) return get_monomer_total_scatt_atom_number(pick_struct.monomer_uid, all_atom);

	int n = 0, nMol = 0;
	LIST<MMOLECULE> *plist = NULL;

	switch (JOB) {
	case GET_SINGLE_MM_MD_PROC:
		return get_macromol_scatt_atom_number(&mm, all_atom);
		break;
	case SHOW_SINGLE_MM_MD_PROC:
		return get_macromol_scatt_atom_number(&mm, all_atom);
		break;
	case GET_MULTI_MM_MD_PROC:
		for (nMol = 0; nMol < mm_md_cell.mm.n; nMol++) {
			n += get_macromol_scatt_atom_number(mm_md_cell.mm.m + nMol, all_atom);
		}
		for (nMol = 0; nMol < mm_md_cell.sm.n; nMol++) {
			n += get_simple_molecule_scatt_atom_number(mm_md_cell.sm.m + nMol, all_atom);
		}
		for (nMol = 0; nMol < mm_md_cell.pm.n; nMol++) {
			n += get_point_molecule_scatt_atom_number(mm_md_cell.pm.m + nMol, all_atom);
		}
		return n;
		break;
	case SHOW_MULTI_MM_MD_PROC:
		for (nMol = 0; nMol < mm_md_cell.mm.n; nMol++) {
			n += get_macromol_scatt_atom_number(mm_md_cell.mm.m + nMol, all_atom);
		}
		for (nMol = 0; nMol < mm_md_cell.sm.n; nMol++) {
			n += get_simple_molecule_scatt_atom_number(mm_md_cell.sm.m + nMol, all_atom);
		}
		for (nMol = 0; nMol < mm_md_cell.pm.n; nMol++) {
			n += get_point_molecule_scatt_atom_number(mm_md_cell.pm.m + nMol, all_atom);
		}
		return n;
		break;
	case SINGLE_MM_FREE_CELL_MD:
		return get_macromol_scatt_atom_number(&mm, all_atom);
		break;
	case MULTI_MM_FREE_CELL_MD:
		for (nMol = 0; nMol < mm_md_cell.mm.n; nMol++) {
			n += get_macromol_scatt_atom_number(mm_md_cell.mm.m + nMol, all_atom);
		}
		for (nMol = 0; nMol < mm_md_cell.sm.n; nMol++) {
			n += get_simple_molecule_scatt_atom_number(mm_md_cell.sm.m + nMol, all_atom);
		}
		for (nMol = 0; nMol < mm_md_cell.pm.n; nMol++) {
			n += get_point_molecule_scatt_atom_number(mm_md_cell.pm.m + nMol, all_atom);
		}
		return n;
		break;
	case MULTI_MM_PERIDIC_CELL_MD:
		for (nMol = 0; nMol < mm_md_cell.mm.n; nMol++) {
			n += get_macromol_scatt_atom_number(mm_md_cell.mm.m + nMol, all_atom);
		}
		for (nMol = 0; nMol < mm_md_cell.sm.n; nMol++) {
			n += get_simple_molecule_scatt_atom_number(mm_md_cell.sm.m + nMol, all_atom);
		}
		for (nMol = 0; nMol < mm_md_cell.pm.n; nMol++) {
			n += get_point_molecule_scatt_atom_number(mm_md_cell.pm.m + nMol, all_atom);
		}
		return n;
		break;
	// coarse-grained
	case GET_SINGLE_CGMM_MD_PROC:
		return cgmm.nCluster;
		break;
	case GET_MULTI_CGMM_MD_PROC:
		for (nMol = 0; nMol < cgmm_md_cell.mm.n; nMol++) {
			n += cgmm_md_cell.mm.m->nCluster;
		}
		n += cgmm_md_cell.sm.n;
		return n;
		break;
	case SHOW_MULTI_CGMM_MD_PROC:
		for (nMol = 0; nMol < cgmm_md_cell.mm.n; nMol++) {
			n += cgmm_md_cell.mm.m->nCluster;
		}
		n += cgmm_md_cell.sm.n;
		return n;
		break;
	case MOLECULE_OPERATION:
		return 0;
	default:
		return 0;
	}
}

float x_ray_scatt(char atype) {
	switch (atype) {
	case 1: // H 
		return 1;
	case 3: // Li 
		return 3;
	case 6: // C
		return 6;
	case 7: // N
		return 7;
	case 8: // O
		return 8;
	case 9: // F
		return 9;
	case 11: // Na
		return 11;
	case 12: // Mg
		return 12;
	case 15: // P
		return 15;
	case 16: // S
		return 16;
	case 17: // Cl
		return 17;
	case 19: // K
		return 19;
	case 20: // Ca
		return 20;
	case 35: // Br
		return 35;
	case 53: // I
		return 53;
	default:
		return 0;
	}
}

float neutron_scatt(char atype, bool deuterized) {
	switch (atype) {
	case 1: // H 
		return (deuterized ? 2 : 1);
	case 3: // Li 
		return 3;
	case 6: // C
		return 6;
	case 7: // N
		return 7;
	case 8: // O
		return 8;
	case 9: // F
		return 9;
	case 11: // Na
		return 11;
	case 12: // Mg
		return 12;
	case 15: // P
		return 15;
	case 16: // S
		return 16;
	case 17: // Cl
		return 17;
	case 19: // K
		return 19;
	case 20: // Ca
		return 20;
	case 35: // Br
		return 35;
	case 53: // I
		return 53;
	default:
		return 0;
	}
}

float get_atom_scatt(int SOURCE, MATOM *atom, bool deuterized) {
	switch(SOURCE) {
	case _X_RAY: 
		return x_ray_scatt(atom->par->uid);
	case _NEUTRON:
		return neutron_scatt(atom->par->uid, deuterized);
	default:
		return 0;
	}
}

bool CGATOM_SCATT_PARS::get_scatt(char *fname, int SOURCE) {
	char msg[256] = "\0", atom[250] = "\0", *bf = NULL;
	int i = 0;
	float scf = 0;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(msg, "can not open %s to get coarse-grained atomic scatt. factor", fname);
		show_log(msg, true); return false;
	}
	while (!in.eof()) {
		in.getline(msg, 250);
		if (sscanf(msg, "%s", atom) == 1 && strcmp(atom, this->atom) == 0) {
			bf = msg; NextSection(&bf);
			for (i = 0; i < SOURCE; i++) {
				if (bf == NULL || !NextSection(&bf)) {this->scf = 0; in.close(); return false;}
			}
			if (sscanf(bf, "%f", &scf) == 1) {
				this->scf = scf; in.close(); return true;
			}
			else {
				this->scf = 0; sprintf(msg, "failure to get scatt. factor of %s with source %d in %s", this->atom, SOURCE, fname);
				show_log(msg, true); in.close(); return false;
			}
		}
	}
	sprintf(msg, "failure to get scatt. factor of %s with source %d in %s", this->atom, SOURCE, fname);
	show_log(msg, true); in.close(); return false;
}

ARRAY<CGATOM_SCATT_PARS> cgatom_scf_db;
char cgatom_scatt_fname[256] = "\0";

extern "C" bool init_cgatom_scatt() {
	int natoms = number<CGATOM_COMM_PARS>(cgatompar_db.ch);
	cgatom_scf_db.SetArray(natoms);
	if (natoms == 0) return true;
	bool status = true;
	CHAIN<CGATOM_COMM_PARS> *ch = cgatompar_db.ch;
	for (int iatom = 0; iatom < natoms; iatom++) {
		cgatom_scf_db.m[iatom].indx = ch->p->indx;
		cgatom_scf_db.m[iatom].uid = ch->p->uid;
		strcpy(cgatom_scf_db.m[iatom].atom, ch->p->atom);

		status = status & cgatom_scf_db.m[iatom].get_scatt(cgatom_scatt_fname, ::LIGHT_SOURCE);
		if (cgatom_scf_db.m[iatom].indx != iatom) {
			sprintf(::errmsg, "strange : when constructing cg-atomic scatt. factor database, atomic index in NEIMO is not consistent");
			show_log(::errmsg, true);
		}
		ch = ch->next;
	}
	return status;
}

float get_cgatom_scatt(int SOURCE, CGATOM *pcgatom) {
	if (cgatom_scf_db.n <= pcgatom->aindx) return 0;
	else return cgatom_scf_db.m[pcgatom->aindx].scf;
}

void get_implicit_atom_scatt(int LIGHT_SOURCE, MATOM *atom, bool deuterize, float &fsc, VECTOR3 &vsc) {
	fsc = 0; V3zero(vsc)
	VECTOR3 r;
	float f = 0;

	if (atom->par->uid == UID_H) return; // H
	f = get_atom_scatt(LIGHT_SOURCE, atom, false); // this atom is not H, so can not be deuterized
	fsc += f;
	V3PC(atom->r, f, r) V3plusV3(vsc, r, vsc)
	
	for (int n = 0; n < atom->nb; n++) {
		if (((MATOM*)atom->ba[n])->par->uid == UID_H) { // H
			f = get_atom_scatt(LIGHT_SOURCE, (MATOM*)atom->ba[n], deuterize);
			fsc += f;
			V3PC(atom->ba[n]->r, f, r) V3plusV3(vsc, r, vsc)
		}
	}
	vsc.v[0] /= fsc; vsc.v[1] /= fsc; vsc.v[2] /= fsc;
}

extern "C" int get_molecule_scatt_struct(bool all_atom, float *fscatt, float *x, float *y, float *z, int MAX, int nMol, int nstart) {
	// specific type of monomer ?
	if (pick_struct.monomer_uid.n > 0) return get_monomer_scatt_struct(all_atom, fscatt, x, y, z, MAX, nMol, nstart);

	int i = 0, j = 0;
	int iatom = 0;
	int nc = 0, natom = 0;
	int n = 0;
	VECTOR3 rsc;
	MATOM *patom = NULL;
	CGATOM *pcgatom = NULL;
	MMOLECULE *m = NULL;
	CG_MMOLECULE *cgmm = NULL;
	SMOLECULE *sm = NULL;
	CGSM *cgsm = NULL;
	GENERAL_MOLECULE gm;
	bool status = get_molecule_upto_job(nMol, gm);
	if (!status) {pick_struct.reset(); return -1;}

	if (gm.mol_type == 1) { // macro-molecule
		m = gm.mm;
		if (nstart == 0) {nc = 0; natom = 0;}
		else if (nMol == pick_struct.nMol) {
			nc = pick_struct.nCluster;
			natom = pick_struct.natom; if (natom >= m->cluster[nc].nAtoms) {nc++; natom = 0;}
			if (nc >= m->nCluster) {pick_struct.reset(); return -1;}
		}
		else {
			if (all_atom) {
				natom = 0; nc = 0; 
				while (nc < m->nCluster && natom < nstart) {
					if (natom + m->cluster[nc].nAtoms < nstart) {
						natom += m->cluster[nc].nAtoms; nc++;
						if (nc == m->nCluster) {pick_struct.reset(); return -1;}
					}
					else {natom = nstart - natom; break;}
				}
			}
			else {
				natom = 0; nc = 0;
				while (nc < m->nCluster) {
					for (i = 0; i < m->cluster[nc].nAtoms; i++) {
						if (m->cluster[nc].atom[i].par->uid != UID_H) natom++; // not H
						if (natom == nstart + 1) break;
					}
					if (natom == nstart + 1) {natom = i; break;}
					else nc++;
				}
			}
		}

		if (nc == m->nCluster) {pick_struct.reset(); n = 0; return n;} // molecule has no more atoms
		n = 0;
		while (n < MAX) {
			for (iatom = natom; iatom < m->cluster[nc].nAtoms; iatom++) {
				patom = m->cluster[nc].atom + iatom;
				if (all_atom) {
					fscatt[n] = get_atom_scatt(LIGHT_SOURCE, patom, m->cluster[nc].deuterized);
					x[n] = float(patom->r.v[0] + m->cluster[nc].dr_cmm.v[0]);
					y[n] = float(patom->r.v[1] + m->cluster[nc].dr_cmm.v[1]);
					z[n] = float(patom->r.v[2] + m->cluster[nc].dr_cmm.v[2]);
					n++;
				}
				else {
					if (patom->par->uid == UID_H) continue; // H
					get_implicit_atom_scatt(LIGHT_SOURCE, patom, m->cluster[nc].deuterized, fscatt[n], rsc);
					x[n] = float(rsc.v[0] + m->cluster[nc].dr_cmm.v[0]);
					y[n] = float(rsc.v[1] + m->cluster[nc].dr_cmm.v[1]);
					z[n] = float(rsc.v[2] + m->cluster[nc].dr_cmm.v[2]);
					n++;
				}
				if (n == MAX) {pick_struct.set(nMol, *m, nc, iatom + 1); return n;} // no rest atom in this macro-molecule is not picked up
			}
			nc++;
			if (nc == m->nCluster) break; //over
			natom = 0; // in new cluster, natom starts from 0
		}
		pick_struct.set(nMol, *m, nc, iatom + 1);
		return n;
	}
	else if (gm.mol_type == 0) { // simple-molecule
		sm = gm.sm;
		pick_struct.reset(); // not macro-molecule anymore
		if (sm->c == NULL || nstart >= sm->c->nAtoms) return -1;
		natom = nstart; n = 0;
		while (n < MAX) {
			for (iatom = natom; iatom < sm->c->nAtoms; iatom++) {
				patom = sm->c->atom + iatom;
				if (all_atom) {
					//atom[i] = patom->par->uid;
					fscatt[n] = get_atom_scatt(LIGHT_SOURCE, patom, sm->c->deuterized);
					x[n] = float(patom->r.v[0] + sm->c->dr_cmm.v[0]); 
					y[n] = float(patom->r.v[1] + sm->c->dr_cmm.v[1]);
					z[n] = float(patom->r.v[2] + sm->c->dr_cmm.v[2]);
					n++;
				}
				else {
					if (patom->par->uid == UID_H) continue; // H
					get_implicit_atom_scatt(LIGHT_SOURCE, patom, sm->c->deuterized, fscatt[n], rsc);
					x[n] = float(rsc.v[0] + sm->c->dr_cmm.v[0]);
					y[n] = float(rsc.v[1] + sm->c->dr_cmm.v[1]);
					z[n] = float(rsc.v[2] + sm->c->dr_cmm.v[2]);
					n++;
				}
				if (n == MAX) return n;
			}
		}
		return n;
	}
	else if (gm.mol_type == 2) { // coarse-grained macro-molecule
		cgmm = gm.cgmm; natom = nstart;
		if (nMol == pick_struct.nMol) {
			if (natom >= cgmm->nCluster) {pick_struct.reset(); n = 0; return n;} // molecule has no more atoms
		}

		n = 0;
		for (iatom = natom; iatom < cgmm->nCluster; iatom++) {
			pcgatom = cgmm->cluster[iatom].atom;
			fscatt[n] = get_cgatom_scatt(LIGHT_SOURCE, pcgatom);
			x[n] = float(pcgatom->r.v[0] + cgmm->cluster[iatom].dr_cmm.v[0]);
			y[n] = float(pcgatom->r.v[1] + cgmm->cluster[iatom].dr_cmm.v[1]);
			z[n] = float(pcgatom->r.v[2] + cgmm->cluster[iatom].dr_cmm.v[2]);
			n++;
			if (n == MAX) {pick_struct.set(nMol, *cgmm, iatom + 1); return n;} 
		}
		pick_struct.reset();
		return n;
	}
	else if (gm.mol_type == 3) { // coarse-grained simple-molecule
		cgsm = gm.cgsm; 
		if (nstart != 0) {pick_struct.reset(); n = 0; return 0;} // cgmm has only one atom
		pcgatom = cgsm->atom;
		fscatt[0] = get_cgatom_scatt(LIGHT_SOURCE, pcgatom);
		x[0] = float(pcgatom->r.v[0] + cgsm->dr_cmm.v[0]);
		y[0] = float(pcgatom->r.v[1] + cgsm->dr_cmm.v[1]);
		z[0] = float(pcgatom->r.v[2] + cgsm->dr_cmm.v[2]);
		pick_struct.reset();
		return 1;
	}
	else return -1;
}

bool READ_MD_PROC_FILE::init(char *mdpar_fname) {
	extern int MD_ENSEMBLE, MD_MODE;

	char buffer[256] = "\0", msg[256] = "\0";
	if (!read_item(mdpar_fname, "[MD_CELL]", buffer)) return false;
	if (sscanf(buffer, "%d", &::MD_MODE) != 1) {
		sprintf(msg, "failure to obtain the type of CELL, see %s", buffer); show_msg(msg, true); return false;
	}
	if (!read_item(mdpar_fname, "[MD_MODE]", buffer)) MD_ENSEMBLE = MD_NVT;
	else if (sscanf(buffer, "%d", &MD_ENSEMBLE) != 1) {
		MD_ENSEMBLE = MD_NVT;
		sprintf(msg, "it is assumed to be NVT essemble %s", buffer); show_msg(msg, true); 
		//return false;
	}

	if (MD_ENSEMBLE == MD_NPT) {
		this->bNPT = true;
		if (MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
			this->xd[0] = ::mm_md_cell.h[0] * 2;
			this->xd[1] = ::mm_md_cell.h[1] * 2;
			this->xd[2] = ::mm_md_cell.h[2] * 2;
		}
		else if (MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			this->xd[0] = ::cgmm_md_cell.h[0] * 2;
			this->xd[1] = ::cgmm_md_cell.h[1] * 2;
			this->xd[2] = ::cgmm_md_cell.h[2] * 2;
		}
		else {
			show_msg("ERROR: cell is NOT periodical, while simulation is NPT", true);
			this->bNPT = false;
			return false;
		}
	}
	else if (MD_ENSEMBLE == MD_NVT) this->bNPT = false;

	return true;
}

bool READ_MD_PROC_FILE::search_proc(int isaved, char* buffer) {
	int ifile = isaved / msave.max_save;
	//isave -= ifile * msave.max_save;
	char fname[256] = "\0";

	if (in == NULL || nfile != ifile || nsaved >= isaved) {
		reset_opened_file();
		sprintf(fname, "%s-%d.kin", msave.title, ifile);
		in = new ifstream;
		in->open(fname);
		if (!in->is_open()) {
			sprintf(errmsg, "Can not open file %s", fname);
			reset_opened_file();
			return false;
		}

		if (bNPT) {
			sprintf(fname, "%s-pv-%d.kin", msave.title, ifile);
			PVin = new ifstream;
			PVin->open(fname);
			if (!PVin->is_open()) {
				sprintf(errmsg, "Can not open file %s", fname);
				reset_opened_file();
				return false;
			}
		}

		this->nfile = ifile;
	}

	bool status = true;
	status = search(*in, "[MD]", isaved, buffer);
	if (bNPT && PVin != NULL) {
		status = status & search(*PVin, "[MD]", isaved, buffer);
	}
	if (status) {
		this->nsaved = isaved; return true;
	}
	else {
		reset_opened_file();
		sprintf(errmsg, "Can not find MD Proc. # %d from %s", isaved, fname);
		return false;
	}
}

READ_MD_PROC_FILE rmdproc;

extern "C" bool init_read_mdproc(char *mdpar_fname) {
	rmdproc.reset_opened_file();
	return rmdproc.init(mdpar_fname) && read_MD_SAVE_Pars(mdpar_fname, rmdproc.msave);
}

extern "C" void clear_read_mdproc() {
	rmdproc.reset_opened_file();
	return;
}

bool serial_single_mol_md_proc(long isave) {
	char buffer[256] = "\0";
	if (!rmdproc.search_proc(isave, buffer)) return false;
	ReadMoleculeStruct(*(rmdproc.in), mm);
#if _SYS_ == _WINDOWS_SYS_
	if (pMainFrame != NULL) ::PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
#endif
	return true;
}

bool serial_single_cgmm_md_proc(long isave) {
	char buffer[256] = "\0";
	if (!rmdproc.search_proc(isave, buffer)) return false;
	ReadCGMMStruct(*(rmdproc.in), cgmm);
#if _SYS_ == _WINDOWS_SYS_
	if (pMainFrame != NULL) ::PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
#endif
	return true;
}

bool read_NPT_PV(ifstream *in, float &eps) {
	char buffer[256] = "\0";
	if (in == NULL) return false;
	if (in->eof()) return false;
	in->getline(buffer, 250);  // V and P
	in->getline(buffer, 250);  // eps, veps, alpha and E, see function, ISOTROPIC_P_MD_SAVE::save_kinetic_pars, in NPT_MD.cpp
	if (sscanf(buffer, "%f", &eps) == 1) return true;
	else return false;
}

bool serial_multi_mol_md_proc(long isave) {
	if (mm_md_cell.mm.n <= 0 && mm_md_cell.sm.n <= 0 && mm_md_cell.pm.n <= 0) return false;

	char buffer[256] = "\0", buff[256] = "\0";
	long nloop = -1;
	if (!rmdproc.search_proc(isave, buffer)) return false;
	sscanf(buffer, "%s %d %ld", buff, &isave, &nloop);

	float eps = 1;
	if (rmdproc.bNPT) {
		if (!read_NPT_PV(rmdproc.PVin, eps)) {
			sprintf(buffer, "failure to get PV of config. %d", isave); show_log(buffer, true); return false;
		}
		else {
			mm_md_cell.h[0] = rmdproc.xd[0] * exp(eps) * 0.5;
			mm_md_cell.h[1] = rmdproc.xd[1] * exp(eps) * 0.5;
			mm_md_cell.h[2] = rmdproc.xd[2] * exp(eps) * 0.5;
		}
	}

	int nm = 0;
	for (nm = 0; nm < mm_md_cell.pm.n; nm++) {
		if (!search(*(rmdproc.in), "PM", nm, buffer)) {
			sprintf(errmsg, "can not find PM %d", nm);
			return false;
		}
		ReadPointMoleculeStruct(*(rmdproc.in), mm_md_cell.pm.m[nm]);
	}
	for (nm = 0; nm < mm_md_cell.sm.n; nm++) {
		if (!search(*(rmdproc.in), "SM", nm, buffer)) {
			sprintf(errmsg, "can not find SM %d", nm);
			return false;
		}
		ReadSolidMoleculeStruct(*(rmdproc.in), mm_md_cell.sm.m[nm]);
	}
	for (nm = 0; nm < mm_md_cell.mm.n; nm++) {
		if (!search(*(rmdproc.in), "MM", nm, buffer)) {
			sprintf(errmsg, "can not find MM %d", nm);
			return false;
		}
		ReadMoleculeStruct(*(rmdproc.in), mm_md_cell.mm.m[nm]);
	}

	molecule_check_periodic_cell(mm_md_cell);

#if _SYS_ == _WINDOWS_SYS_
	if (pMainFrame != NULL) {
		::PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
	}
#endif
	return true;
}


bool serial_multi_cgmm_md_proc(long isave) {
	if (cgmm_md_cell.mm.n <= 0 && cgmm_md_cell.sm.n <= 0 && cgmm_md_cell.pm.n <= 0) return false;
	char buffer[256] = "\0", buff[256] = "\0";
	int nloop = -1;
	if (!rmdproc.search_proc(isave, buffer)) return false;
	sscanf(buffer, "%s %d %ld", buff, &isave, &nloop);

	float eps = 1;
	if (rmdproc.bNPT) {
		if (!read_NPT_PV(rmdproc.PVin, eps)) {
			sprintf(buffer, "failure to get PV of config. %d", isave); show_log(buffer, true); return false;
		}
		else {
			mm_md_cell.h[0] = rmdproc.xd[0] * eps * 0.5;
			mm_md_cell.h[1] = rmdproc.xd[1] * eps * 0.5;
			mm_md_cell.h[2] = rmdproc.xd[2] * eps * 0.5;
		}
	}

	int nm = 0;
	for (nm = 0; nm < cgmm_md_cell.pm.n; nm++) {
		if (!search(*(rmdproc.in), "PM", nm, buffer)) {
			sprintf(errmsg, "can not find PM %d", nm);
			return false;
		}
		ReadCGSMStruct(*(rmdproc.in), cgmm_md_cell.pm.m[nm]);
	}
	for (nm = 0; nm < cgmm_md_cell.sm.n; nm++) {
		if (!search(*(rmdproc.in), "SM", nm, buffer)) {
			sprintf(errmsg, "can not find SM %d", nm);
			return false;
		}
		ReadCGSMStruct(*(rmdproc.in), cgmm_md_cell.sm.m[nm]);
	}
	for (nm = 0; nm < cgmm_md_cell.mm.n; nm++) {
		if (!search(*(rmdproc.in), "MM", nm, buffer)) {
			sprintf(errmsg, "can not find MM %d", nm);
			return false;
		}
		ReadCGMMStruct(*(rmdproc.in), cgmm_md_cell.mm.m[nm]);
	}

	_coarse_grain_::molecule_check_periodic_cell(cgmm_md_cell);

#if _SYS_ == _WINDOWS_SYS_
	if (pMainFrame != NULL) {
		::PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
	}
#endif
	return true;
}

extern "C" bool serial_md_proc(long isave) {
	switch (JOB) {
	case GET_SINGLE_MM_MD_PROC:
		return serial_single_mol_md_proc(isave);
		break;
	case GET_MULTI_MM_MD_PROC:
		return serial_multi_mol_md_proc(isave);
		break;
	case GET_SINGLE_CGMM_MD_PROC:
		return serial_single_cgmm_md_proc(isave);
		break;
	case GET_MULTI_CGMM_MD_PROC:
		return serial_multi_cgmm_md_proc(isave);
		break;
	default:
		return false;
	}
}


/******************************************************************
******               SPECIFIC TYPE OF MONOMER                ******
******************************************************************/

extern char struct_fname[256];

bool get_defined_monomer_uid(char *monomer, int& monomer_uid) {
	ifstream in;
	char title[256] = "\0", buffer[256] = "\0";
	in.open(struct_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open monomer structure database"); show_log(errmsg, true); in.close(); return false;}
	sprintf(title, "[%s]", monomer);
	if (!search(in, title, buffer)) {sprintf(errmsg, "monomer %s is not defined", monomer); show_log(errmsg, true); in.close(); return false;}
	sprintf(title, "UID");
	if (!search(in, title, buffer)) {
		sprintf(errmsg, "UID of monomer %s is not given", monomer); show_log(errmsg, true); in.close(); return false;
	}
	if (sscanf(buffer, "%s %d", title, &monomer_uid) != 2) {
		sprintf(errmsg, "UID of monomer %s is not given", monomer); show_log(errmsg, true); in.close(); return false;
	}
	in.close();
	return true;
}

extern "C" void set_pickup_monomer(char *monomer) {
	pick_struct.monomer_uid.SetArray(0);
	int monomer_uid = -1;
	if (monomer == NULL || strlen(monomer) == 0) return;

	if (JOB == SINGLE_MM_FREE_CELL_MD || JOB == MULTI_MM_FREE_CELL_MD || JOB == MULTI_MM_PERIDIC_CELL_MD ||
		JOB == GET_SINGLE_MM_MD_PROC || JOB == SHOW_SINGLE_MM_MD_PROC || 
		JOB == GET_MULTI_MM_MD_PROC || JOB == SHOW_MULTI_MM_MD_PROC) {
			if (!get_defined_monomer_uid(monomer, monomer_uid)) return;
	}
	else if (JOB == SINGLE_CGMM_FREE_CELL_MD || JOB == MULTI_CGMM_FREE_CELL_MD ||
		JOB == MULTI_CGMM_PERIDIC_CELL_MD || JOB == GET_SINGLE_CGMM_MD_PROC ||
		JOB == SHOW_SINGLE_CGMM_MD_PROC || JOB == GET_MULTI_CGMM_MD_PROC ||
		JOB == SHOW_MULTI_CGMM_MD_PROC) {
			if (!get_cgatom_uid(monomer, monomer_uid)) return;
	}

	pick_struct.monomer_uid.SetArray(1);
	pick_struct.monomer_uid.m[0] = monomer_uid;
}

extern "C" void add_pickup_monomer(char *monomer) {
	int monomer_uid = -1;
	if (monomer == NULL || strlen(monomer) == 0) return;

	if (JOB == SINGLE_MM_FREE_CELL_MD || JOB == MULTI_MM_FREE_CELL_MD || JOB == MULTI_MM_PERIDIC_CELL_MD ||
		JOB == GET_SINGLE_MM_MD_PROC || JOB == SHOW_SINGLE_MM_MD_PROC || 
		JOB == GET_MULTI_MM_MD_PROC || JOB == SHOW_MULTI_MM_MD_PROC) {
			if (!get_defined_monomer_uid(monomer, monomer_uid)) return;
	}
	else if (JOB == SINGLE_CGMM_FREE_CELL_MD || JOB == MULTI_CGMM_FREE_CELL_MD ||
		JOB == MULTI_CGMM_PERIDIC_CELL_MD || JOB == GET_SINGLE_CGMM_MD_PROC ||
		JOB == SHOW_SINGLE_CGMM_MD_PROC || JOB == GET_MULTI_CGMM_MD_PROC ||
		JOB == SHOW_MULTI_CGMM_MD_PROC) {
			if (!get_cgatom_uid(monomer, monomer_uid)) return;
	}

	int n = pick_struct.monomer_uid.n;
	int *uid = pick_struct.monomer_uid.m;
	pick_struct.monomer_uid.m = NULL; pick_struct.monomer_uid.n = 0;
	pick_struct.monomer_uid.SetArray(n+1);
	for (int i = 0; i < pick_struct.monomer_uid.n; i++) {
		if (i < n) pick_struct.monomer_uid.m[i] = uid[i];
		else pick_struct.monomer_uid.m[i] = monomer_uid;
	}
	delete[] uid; uid = NULL;
}

int get_macromol_monomer_scatt_atom_number(MMOLECULE *mm, ARRAY<int>& monomer_uid, bool all_atom) {
	if (all_atom) return get_macromol_atom_number(mm);

	int n = 0, nc = 0, na = 0, nm = 0;
	_GROUP<CLUSTER> *pm = NULL;
	CLUSTER *pc = NULL;
	for (nm = 0; nm < mm->nMonomer; nm++) {
		pm = mm->monomer + nm;
		if (monomer_uid.m != NULL && (!monomer_uid.exist(pm->uid))) continue;
		for (nc = 0; nc < pm->nUnit; nc++) {
			pc = pm->u[nc];
			for (na = 0; na < pc->nAtoms; na++) {
				if (all_atom) n++;
				else if (pc->atom[na].par->uid != UID_H) n++; // 1 is H
			}
		}
	}
	return n;
}

int get_CGMM_monomer_scatt_atom_number(CG_MMOLECULE *mm, ARRAY<int>& uid) {
	int n = 0, nc = 0, na = 0;
	CG_CLUSTER *pc = NULL;
	for (nc = 0; nc < mm->nCluster; nc++) {
		pc = mm->cluster + nc;
		if (uid.m != NULL && (!uid.exist(pc->cID))) continue;
		n += pc->nAtoms;
	}
	return n;
}

int get_CGSM_monomer_scatt_atom_number(CGSM *sm, ARRAY<int>& uid) {
	int n = 0;
	if (uid.m != NULL && !uid.exist(sm->cID)) return 0;
	else return sm->nAtoms;
}


// get the atoms in current job
int get_monomer_total_scatt_atom_number(ARRAY<int>& monomer_uid, bool all_atom) {
	int n = 0, nMol = 0;
	LIST<MMOLECULE> *plist = NULL;

	switch (JOB) {
	case GET_SINGLE_MM_MD_PROC:
		return get_macromol_monomer_scatt_atom_number(&mm, monomer_uid, all_atom);
		break;
	case SHOW_SINGLE_MM_MD_PROC:
		return get_macromol_monomer_scatt_atom_number(&mm, monomer_uid, all_atom);
		break;
	case GET_MULTI_MM_MD_PROC:
		for (nMol = 0; nMol < mm_md_cell.mm.n; nMol++) {
			n += get_macromol_monomer_scatt_atom_number(mm_md_cell.mm.m + nMol, monomer_uid, all_atom);
		}
		return n;
		break;
	case SHOW_MULTI_MM_MD_PROC:
		for (nMol = 0; nMol < mm_md_cell.mm.n; nMol++) {
			n += get_macromol_monomer_scatt_atom_number(mm_md_cell.mm.m + nMol, monomer_uid, all_atom);
		}
		return n;
		break;
	case SINGLE_MM_FREE_CELL_MD:
		return get_macromol_monomer_scatt_atom_number(&mm, monomer_uid, all_atom);
		break;
	case MULTI_MM_FREE_CELL_MD:
		for (nMol = 0; nMol < mm_md_cell.mm.n; nMol++) {
			n += get_macromol_monomer_scatt_atom_number(mm_md_cell.mm.m + nMol, monomer_uid, all_atom);
		}
		return n;
		break;
	case MULTI_MM_PERIDIC_CELL_MD:
		for (nMol = 0; nMol < mm_md_cell.mm.n; nMol++) {
			n += get_macromol_monomer_scatt_atom_number(mm_md_cell.mm.m + nMol, monomer_uid, all_atom);
		}
		return n;
		break;
	// coarse-grained
	case GET_SINGLE_CGMM_MD_PROC:
		return get_CGMM_monomer_scatt_atom_number(&cgmm, monomer_uid);
		break;
	case GET_MULTI_CGMM_MD_PROC:
		for (nMol = 0; nMol < cgmm_md_cell.mm.n; nMol++) {
			n += get_CGMM_monomer_scatt_atom_number(cgmm_md_cell.mm.m + nMol, monomer_uid);
		}
		for (nMol = 0; nMol < cgmm_md_cell.sm.n; nMol++) {
			n += get_CGSM_monomer_scatt_atom_number(cgmm_md_cell.sm.m + nMol, monomer_uid);
		}
		return n;
		break;
	case MOLECULE_OPERATION:
		return 0;
	default:
		return 0;
	}
}

int get_monomer_scatt_struct(bool all_atom, float *fscatt, float *x, float *y, float *z, int MAX, int nMol, int nstart) {
	int i = 0, j = 0;
	int iatom = 0;
	int nm = 0, nc = 0, natom = 0;
	int n = 0;
	VECTOR3 rsc;
	MATOM *patom = NULL;
	CGATOM *pcgatom = NULL;
	MMOLECULE *m = NULL;
	CG_MMOLECULE *cgmm = NULL;
	CGSM *cgsm = NULL;
	_GROUP<CLUSTER> *pm = NULL;
	GENERAL_MOLECULE gm;
	bool status = get_molecule_upto_job(nMol, gm);
	if (!status) {pick_struct.reset(); return -1;}

	if (gm.mol_type == 1) { // macro-molecule
		m = gm.mm;
		if (nstart == 0) {
			nm = 0; nc = 0; natom = 0; pick_struct.reset();
			while (nm < m->nMonomer && !pick_struct.exist_uid(m->monomer[nm].uid)) nm++;
			if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
		}
		else if (nMol == pick_struct.nMol) {
			nm = pick_struct.nMonomer; nc = pick_struct.nCluster; natom = pick_struct.natom; 
			if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
			if (nc >= m->monomer[nm].nUnit) {nm++; nc = 0; natom = 0;}
			if (natom >= m->monomer[nm].u[nc]->nAtoms) {nc++; natom = 0;}
			
			if (nm < m->nMonomer) {
				if (!pick_struct.exist_uid(m->monomer[nm].uid)) {
					while (nm < m->nMonomer && m->monomer[nm].uid) nm++;
					if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
				}
			}
			else {pick_struct.reset(); return -1;}
		}
		else {
			if (all_atom) {
				nm = 0; nc = 0; natom = 0; 
				while (natom < nstart) {
					if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
					else if (!pick_struct.exist_uid(m->monomer[nm].uid)) {
						nm++; nc = 0;
					}
					else if (nc >= m->monomer[nm].nUnit) {nm++; nc = 0;}
					else if (natom + m->monomer[nm].u[nc]->nAtoms < nstart) {
						natom += m->monomer[nm].u[nc]->nAtoms; nc++;
						if (nc == m->monomer[nm].nUnit) {nm++; nc = 0;}
					}
					else {natom = nstart - natom; break;}
				}
			}
			else {
				nm = 0; nc = 0; natom = 0; 
				while (natom < nstart) {
					if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
					else if (!pick_struct.exist_uid(m->monomer[nm].uid)) {
						nm++; nc = 0; natom = 0;
					}
					else if (nc >= m->monomer[nm].nUnit) {nm++; nc = 0; natom = 0;}
					else {
						for (i = 0; i < m->monomer[nm].u[nc]->nAtoms; i++) {
							if (m->monomer[nm].u[nc]->atom[i].par->uid != UID_H) natom++; // not H
							if (natom == nstart + 1) break;
						}
						if (natom == nstart + 1) {natom = i; break;}
						else {
							nc++; if (nc == m->monomer[nm].nUnit) {nm++; nc = 0; natom = 0;}
						}
					}
				}
			}
		}

		if (nm >= m->nMonomer) {pick_struct.reset(); n = 0; return n;} // molecule has no more atoms
		n = 0;
		while (n < MAX) {
			if (nm >= m->nMonomer) break;
			if (!pick_struct.exist_uid(m->monomer[nm].uid)) {
				nm++; nc = 0; natom = 0; continue;
			}
			for (iatom = natom; iatom < m->monomer[nm].u[nc]->nAtoms; iatom++) {
				patom = m->monomer[nm].u[nc]->atom + iatom;
				if (all_atom) {
					fscatt[n] = get_atom_scatt(LIGHT_SOURCE, patom, m->monomer[nm].u[nc]->deuterized);
					x[n] = float(patom->r.v[0] + m->monomer[nm].u[nc]->dr_cmm.v[0]);
					y[n] = float(patom->r.v[1] + m->monomer[nm].u[nc]->dr_cmm.v[1]);
					z[n] = float(patom->r.v[2] + m->monomer[nm].u[nc]->dr_cmm.v[2]);
					n++;
				}
				else {
					if (patom->par->uid == UID_H) continue; // H
					get_implicit_atom_scatt(LIGHT_SOURCE, patom, m->monomer[nm].u[nc]->deuterized, fscatt[n], rsc);
					x[n] = float(rsc.v[0] + m->monomer[nm].u[nc]->dr_cmm.v[0]);
					y[n] = float(rsc.v[1] + m->monomer[nm].u[nc]->dr_cmm.v[1]);
					z[n] = float(rsc.v[2] + m->monomer[nm].u[nc]->dr_cmm.v[2]);
					n++;
				}
				if (n == MAX) {pick_struct.set(nMol, *m, nm, nc, iatom + 1); return n;} // no rest atom in this macro-molecule is not picked up
			}
			nc++;
			if (nc >= m->monomer[nm].nUnit) {nm++; nc = 0;}
			if (nm >= m->nMonomer) break; //over
			natom = 0; // in new cluster, natom starts from 0
		}
		pick_struct.set(nMol, *m, nm, nc, iatom + 1);
		return n;
	}
	else if (gm.mol_type == 2) { // coarse-grained macro-molecule
		cgmm = gm.cgmm;
		if (nMol == pick_struct.nMol) {
			if (natom >= cgmm->nCluster) {pick_struct.reset(); n = 0; return n;} // molecule has no more atoms
			else natom = pick_struct.natom;
		}
		else natom = 0;

		n = 0;
		for (iatom = natom; iatom < cgmm->nCluster; iatom++) {
			if (pick_struct.monomer_uid.m != NULL && !pick_struct.exist_uid(cgmm->cluster[iatom].cID)) continue;
			pcgatom = cgmm->cluster[iatom].atom;
			fscatt[n] = get_cgatom_scatt(LIGHT_SOURCE, pcgatom);
			x[n] = float(pcgatom->r.v[0] + cgmm->cluster[iatom].dr_cmm.v[0]);
			y[n] = float(pcgatom->r.v[1] + cgmm->cluster[iatom].dr_cmm.v[1]);
			z[n] = float(pcgatom->r.v[2] + cgmm->cluster[iatom].dr_cmm.v[2]);
			n++;
			if (n == MAX) {pick_struct.set(nMol, *cgmm, iatom + 1); return n;} 
		}
		pick_struct.reset();
		return n;
	}
	else if (gm.mol_type == 3) { // coarse-grained simple-molecule
		cgsm = gm.cgsm; 
		if (pick_struct.monomer_uid.m != NULL && !pick_struct.exist_uid(cgsm->cID)) {
			pick_struct.reset(); return 0;
		}
		if (nstart != 0) {pick_struct.reset(); n = 0; return 0;} // cgmm has only one atom
		pcgatom = cgsm->atom;
		fscatt[0] = get_cgatom_scatt(LIGHT_SOURCE, pcgatom);
		x[0] = float(pcgatom->r.v[0] + cgsm->dr_cmm.v[0]);
		y[0] = float(pcgatom->r.v[1] + cgsm->dr_cmm.v[1]);
		z[0] = float(pcgatom->r.v[2] + cgsm->dr_cmm.v[2]);
		pick_struct.reset();
		return 1;
	}
	else return -1;
}


extern "C" int get_macromol_scatt_atom_number(int iMol, bool all_atom) {
	GENERAL_MOLECULE gm;
	bool status = get_molecule_upto_job(iMol, gm);
	if (!status) return -1;
	else if (gm.mol_type == 1 && gm.mm != NULL) 
		return get_macromol_monomer_scatt_atom_number(gm.mm, pick_struct.monomer_uid, all_atom);
	else if (gm.mol_type == 2 && gm.cgmm != NULL)
		return get_CGMM_monomer_scatt_atom_number(gm.cgmm, pick_struct.monomer_uid);
	else if (gm.mol_type == 3 && gm.cgsm != NULL)
		return get_CGSM_monomer_scatt_atom_number(gm.cgsm, pick_struct.monomer_uid);
	else return 0;
}

// get macro-molecule's scatters
extern "C" int get_macromol_scatt_struct(int iMol, bool bInCMM, bool all_atom, float *fscatt, float *x, float *y, float *z, int MAX, int nstart) {
	int i = 0, j = 0;
	int iatom = 0;
	int nm = 0, nc = 0, natom = 0;
	int n = 0;
	VECTOR3 rsc;
	MATOM *patom = NULL;
	CGATOM *pcgatom = NULL;
	MMOLECULE *m = NULL;
	CG_MMOLECULE *cgmm = NULL;
	CGSM *cgsm = NULL;
	_GROUP<CLUSTER> *pm = NULL;
	GENERAL_MOLECULE gm;
	bool status = get_molecule_upto_job(iMol, gm);
	if (!status) {pick_struct.reset(); return -1;}

	if (gm.mol_type == 1) { // macro-molecule
		m = gm.mm;
		if (nstart == 0) {
			nm = 0; nc = 0; natom = 0; pick_struct.reset(); pick_struct.nMol = iMol;
			while (nm < m->nMonomer && !pick_struct.exist_uid(m->monomer[nm].uid)) nm++;
			if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
		}
		else if (iMol == pick_struct.nMol) {
			nm = pick_struct.nMonomer; nc = pick_struct.nCluster; natom = pick_struct.natom; 
			if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
			if (nc >= m->monomer[nm].nUnit) {nm++; nc = 0; natom = 0;}
			if (natom >= m->monomer[nm].u[nc]->nAtoms) {nc++; natom = 0;}
			
			if (nm < m->nMonomer) {
				if (!pick_struct.exist_uid(m->monomer[nm].uid)) {
					while (nm < m->nMonomer && m->monomer[nm].uid) nm++;
					if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
				}
			}
			else {pick_struct.reset(); return -1;}
		}
		else {
			if (all_atom) {
				nm = 0; nc = 0; natom = 0; 
				while (natom < nstart) {
					if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
					else if (!pick_struct.exist_uid(m->monomer[nm].uid)) {
						nm++; nc = 0;
					}
					else if (nc >= m->monomer[nm].nUnit) {nm++; nc = 0;}
					else if (natom + m->monomer[nm].u[nc]->nAtoms < nstart) {
						natom += m->monomer[nm].u[nc]->nAtoms; nc++;
						if (nc == m->monomer[nm].nUnit) {nm++; nc = 0;}
					}
					else {natom = nstart - natom; break;}
				}
			}
			else {
				nm = 0; nc = 0; natom = 0; 
				while (natom < nstart) {
					if (nm >= m->nMonomer) {pick_struct.reset(); return -1;}
					else if (!pick_struct.exist_uid(m->monomer[nm].uid)) {
						nm++; nc = 0; natom = 0;
					}
					else if (nc >= m->monomer[nm].nUnit) {nm++; nc = 0; natom = 0;}
					else {
						for (i = 0; i < m->monomer[nm].u[nc]->nAtoms; i++) {
							if (m->monomer[nm].u[nc]->atom[i].par->uid != UID_H) natom++; // not H
							if (natom == nstart + 1) break;
						}
						if (natom == nstart + 1) {natom = i; break;}
						else {
							nc++; if (nc == m->monomer[nm].nUnit) {nm++; nc = 0; natom = 0;}
						}
					}
				}
			}
		}

		if (nm >= m->nMonomer) {pick_struct.reset(); n = 0; return n;} // molecule has no more atoms
		n = 0;
		while (n < MAX) {
			if (nm >= m->nMonomer) break;
			if (!pick_struct.exist_uid(m->monomer[nm].uid)) {
				nm++; nc = 0; natom = 0; continue;
			}
			for (iatom = natom; iatom < m->monomer[nm].u[nc]->nAtoms; iatom++) {
				patom = m->monomer[nm].u[nc]->atom + iatom;
				if (all_atom) {
					fscatt[n] = get_atom_scatt(LIGHT_SOURCE, patom, m->monomer[nm].u[nc]->deuterized);
					x[n] = float(patom->r.v[0]); if (bInCMM) x[n] += float(m->monomer[nm].u[nc]->dr_cmm.v[0]);
					y[n] = float(patom->r.v[1]); if (bInCMM) y[n] += float(m->monomer[nm].u[nc]->dr_cmm.v[1]);
					z[n] = float(patom->r.v[2]); if (bInCMM) z[n] += float(m->monomer[nm].u[nc]->dr_cmm.v[2]);
					n++;
				}
				else {
					if (patom->par->uid == UID_H) continue; // H
					get_implicit_atom_scatt(LIGHT_SOURCE, patom, m->monomer[nm].u[nc]->deuterized, fscatt[n], rsc);
					x[n] = float(rsc.v[0]); if (bInCMM) x[n] += float(m->monomer[nm].u[nc]->dr_cmm.v[0]);
					y[n] = float(rsc.v[1]); if (bInCMM) y[n] += float(m->monomer[nm].u[nc]->dr_cmm.v[1]);
					z[n] = float(rsc.v[2]); if (bInCMM) z[n] += float(m->monomer[nm].u[nc]->dr_cmm.v[2]);
					n++;
				}
				if (n == MAX) {pick_struct.set(iMol, *m, nm, nc, iatom + 1); return n;} // no rest atom in this macro-molecule is not picked up
			}
			nc++;
			if (nc >= m->monomer[nm].nUnit) {nm++; nc = 0;}
			if (nm >= m->nMonomer) break; //over
			natom = 0; // in new cluster, natom starts from 0
		}
		pick_struct.set(iMol, *m, nm, nc, iatom + 1);
		return n;
	}
	else if (gm.mol_type == 2) { // coarse-grained macro-molecule
		cgmm = gm.cgmm;
		if (iMol == pick_struct.nMol) {
			if (natom >= cgmm->nCluster) {pick_struct.reset(); n = 0; return n;} // molecule has no more atoms
			else natom = pick_struct.natom;
		}
		else natom = 0;

		n = 0;
		for (iatom = natom; iatom < cgmm->nCluster; iatom++) {
			if (pick_struct.monomer_uid.m != NULL && !pick_struct.exist_uid(cgmm->cluster[iatom].cID)) continue;
			pcgatom = cgmm->cluster[iatom].atom;
			fscatt[n] = get_cgatom_scatt(LIGHT_SOURCE, pcgatom);
			x[n] = float(pcgatom->r.v[0]); if (bInCMM) x[n] += float(cgmm->cluster[iatom].dr_cmm.v[0]);
			y[n] = float(pcgatom->r.v[1]); if (bInCMM) y[n] += float(cgmm->cluster[iatom].dr_cmm.v[1]);
			z[n] = float(pcgatom->r.v[2]); if (bInCMM) z[n] += float(cgmm->cluster[iatom].dr_cmm.v[2]);
			n++;
			if (n == MAX) {pick_struct.set(iMol, *cgmm, iatom + 1); return n;} 
		}
		pick_struct.reset();
		return n;
	}
	else if (gm.mol_type == 3) { // coarse-grained simple-molecule
		cgsm = gm.cgsm; 
		if (pick_struct.monomer_uid.m != NULL && !pick_struct.exist_uid(cgsm->cID)) {
			pick_struct.reset(); return 0;
		}
		if (nstart != 0) {pick_struct.reset(); n = 0; return 0;} // cgmm has only one atom
		pcgatom = cgsm->atom;
		fscatt[0] = get_cgatom_scatt(LIGHT_SOURCE, pcgatom);
		x[0] = float(pcgatom->r.v[0]); if (bInCMM) x[0] += float(cgsm->dr_cmm.v[0]);
		y[0] = float(pcgatom->r.v[1]); if (bInCMM) y[0] += float(cgsm->dr_cmm.v[1]);
		z[0] = float(pcgatom->r.v[2]); if (bInCMM) z[0] += float(cgsm->dr_cmm.v[2]);
		pick_struct.reset();
		return 1;
	}
	else return -1;
}

