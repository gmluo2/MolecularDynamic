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

#include "Interaction.h"
#include "Interact1.h"

#include "md_cell.h"
#include "cg-md-cell.h"
#include "test.h"
#define _READ_MM_MD_ 1
#include "read.h"
#include "var.h"
#include "ReadMDProc.h"
#include "mol-cell.h"

using namespace _coarse_grain_;

extern MMOLECULE mm;
extern CG_MMOLECULE cgmm;
extern MMOL_MD_CELL mm_md_cell;
extern CGMM_MD_CELL cgmm_md_cell;
extern LIST<MMOLECULE> *mlist;

extern bool ReadCGMMStruct(ifstream& in, CG_MMOLECULE& mm);
extern bool ReadCGSMStruct(ifstream& in, CGSM& sm);

extern int JOB;
extern int MD_MODE;

extern bool md_run;

#if _SYS_ == _WINDOWS_SYS_
extern CWnd *pMainFrame;
extern HANDLE hWaitGetMolStructEvent;
extern HANDLE hWaitEvent;
#endif

extern int delay_time; // in milisecond
extern int max_time_pickup_molstruct;  // in milisecond
extern int nshow; // the molecule is shown every nshow loops
extern int ishow;

extern int MD_MODE, MD_ENSEMBLE;

template <class MM> bool get_single_mol_md_proc(char *mdpar_fname, int isave, void *read_molstruct, MM *m) {
	typedef bool (*READ_MOL_STRUCT)(ifstream&, MM&);

	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {
		show_msg(::errmsg);
		return false;
	}
	int ifile = isave / msave.max_save;
	//isave -= ifile * msave.max_save;
	char fname[256] = "\0";
	char buffer[256] = "\0";
	ifstream in;
	sprintf(fname, "%s-%d.kin", msave.title, ifile);
	in.open(fname);
	if (!in.is_open()) {
		sprintf(errmsg, "Can not open file %s", fname);
		show_msg(::errmsg);
		return false;
	}
	if (!search(in, "[MD]", isave, buffer)) {
		in.close();
		sprintf(errmsg, "Can not find MD Proc. # %d from %s", isave, fname);
		show_msg(::errmsg);
		return false;
	}
	//ReadMoleculeStruct(in, mm);
	((READ_MOL_STRUCT)read_molstruct)(in, *m);
	in.close();
#if _SYS_ == _WINDOWS_SYS_
	if (pMainFrame != NULL) ::PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
#endif
	return true;
}

bool get_single_mm_md_proc(char *mdpar_fname, int isave) {
	return get_single_mol_md_proc<MMOLECULE>(mdpar_fname, isave, (void*)(&ReadMoleculeStruct), &mm);
}

bool get_single_cgmm_md_proc(char *mdpar_fname, int isave) {
	return get_single_mol_md_proc<CG_MMOLECULE>(mdpar_fname, isave, (void*)(&ReadCGMMStruct), &cgmm);
}

template <class MM> bool show_single_mol_md_proc(char *mdpar_fname, void *read_molstruct, MM *m, int nf) {
	typedef bool (*READ_MOL_STRUCT)(ifstream&, MM&);
	if (nf < 0) nf = 0;

	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {
		show_msg(errmsg);
		return false;
	}
	int ifile = 0, nfile = 0;
	int isave = nf, nloop = 0;
	char fname[256] = "\0";
	char buffer[256] = "\0", buff[20] = "\0";
	ifstream *in = NULL;
	ifstream *Ein = NULL;

	bool proc = true;
	int n = 0;
	long LOOPS = 0;

	ifstream par_in;
	par_in.open(mdpar_fname);
	if (search(par_in, "[TOTAL_LOOPS]", buffer)) {
		par_in.getline(buffer, 250);
		sscanf(buffer, "%ld", &LOOPS);
		LOOPS *= msave.nsave;
		par_in.close();
	}
	else {
		par_in.close();
		sprintf(errmsg, "Can not get [TOTAL_LOOPS] in %s", mdpar_fname);
		show_msg(errmsg);
		return false;
	}

	float Ek = 0, Ep = 0;

	COMMAND_MD_STOP = false;

	while (proc) {
		if (COMMAND_MD_STOP) break; // stop with given command
		ifile = isave / msave.max_save;
		if (in == NULL || ifile != nfile) {
			if (in != NULL) {if (in->is_open()) in->close(); delete in; in = NULL;}
			in = new ifstream;
			sprintf(fname, "%s-%d.kin", msave.title, ifile);
			in->open(fname);
			if (!in->is_open()) {
				proc = false;
				delete in; in = NULL;
				if (n == 0) {
					sprintf(errmsg, "Can not open file %s", fname);
					show_msg(::errmsg);
					return false;
				}
			}
		}
		if (Ein == NULL || ifile != nfile) {
			if (Ein != NULL) {if (Ein->is_open()) Ein->close(); delete Ein; Ein = NULL;}
			Ein = new ifstream;
			sprintf(fname, "%s-%d.E", msave.title, ifile);
			Ein->open(fname);
			if (!Ein->is_open()) {
				proc = false;
				delete Ein; Ein = NULL;
				if (n == 0) {
					sprintf(errmsg, "Can not open file %s", fname);
					show_msg(errmsg);
					return false;
				}
			}
			nfile = ifile;
		}
		if (search((*Ein), "[MD]", isave, buffer)) {
			Ein->getline(buffer, 250);
			sscanf(buffer, "%f", &Ek);
			Ein->getline(buffer, 250);
			sscanf(buffer, "%f", &Ep);
		}
		if (!search((*in), "[MD]", isave, buffer)) proc = false;
		else {
			sscanf(buffer, "%s %d %d", buff, &isave, &nloop);
			show_loop_infor(nloop, isave, LOOPS - nloop, Ek, Ep);
			//ReadMoleculeStruct((*in), mm);
			((READ_MOL_STRUCT)read_molstruct)(*in, *m);


		#if _SYS_ == _WINDOWS_SYS_
			if (pMainFrame != NULL) {
				if (delay_time > 0) WaitForSingleObject(hWaitEvent, delay_time);
				ResetEvent(hWaitGetMolStructEvent);
				PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
				WaitForSingleObject(hWaitGetMolStructEvent, max_time_pickup_molstruct);
			}
		#endif

			isave += nshow;
		}
	}
	if (in != NULL) {if (in->is_open()) in->close(); delete in; in = NULL;}
	if (Ein != NULL) {if (Ein->is_open()) Ein->close(); delete Ein; Ein = NULL;}
	return true;
}

bool show_single_mm_md_proc(char *mdpar_fname, int nf) {
	return show_single_mol_md_proc<MMOLECULE>(mdpar_fname, (void*)(&ReadMoleculeStruct), &mm, nf);
}

bool show_single_cgmm_md_proc(char *mdpar_fname, int nf) {
	return show_single_mol_md_proc<CG_MMOLECULE>(mdpar_fname, (void*)(&ReadCGMMStruct), &cgmm, nf);
}

ifstream* open_multi_mol_md_file(ifstream *in, int ncur_file, char *title, char *ext, int ifile) {
	char fname[256] = "\0";
	sprintf(fname, "%s-%d.%s", title, ifile, ext);
	if (ncur_file == ifile) return in;
	else if (in != NULL) {if (in->is_open()) in->close(); delete in; in = NULL;}

	if (in == NULL) in = new ifstream;
	in->open(fname);
	if (!in->is_open()) {sprintf(errmsg, "can not open file %s", fname); delete in; in = NULL; return NULL;}
	return in;
}

template <class PM, class SM, class MM> bool get_multi_mol_md_confg(ifstream** in, char *title, int ifile, int ncur_file, int isave, int& nloop, ARRAY<PM> &pm, ARRAY<SM> &sm, ARRAY<MM> &mm, void *read_pmstruct, void *read_smstruct, void *read_mmstruct) {
	typedef bool (*READ_MM_STRUCT)(ifstream&, MM&);
	typedef bool (*READ_SM_STRUCT)(ifstream&, SM&);
	typedef bool (*READ_PM_STRUCT)(ifstream&, PM&);

	char buffer[256] = "\0", buff[20] = "\0";

	*in = open_multi_mol_md_file(*in, ncur_file, title, "kin", ifile);
	if (*in == NULL) {show_msg(errmsg); return false;}

	if (!search(**in, "[MD]", isave, buffer)) {
		sprintf(errmsg, "Can not find MD Proc. # %d", isave);
		show_msg(errmsg); delete *in; *in = NULL; return false;
	}
	sscanf(buffer, "%s %d %d", buff, &isave, &nloop);
	int nm = 0;
	for (nm = 0; nm < pm.n; nm++) {
		if (!search(**in, "PM", nm, buffer)) {
			(*in)->close(); delete (*in); (*in) = NULL;
			sprintf(errmsg, "can not find PM %d", nm); show_msg(errmsg);
			return false;
		}
		((READ_PM_STRUCT)read_pmstruct)(**in, pm.m[nm]);
	}
	for (nm = 0; nm < sm.n; nm++) {
		if (!search(**in, "SM", nm, buffer)) {
			(*in)->close(); delete (*in); (*in) = NULL;
			sprintf(errmsg, "can not find SM %d", nm); show_msg(errmsg);
			return false;
		}
		((READ_SM_STRUCT)read_smstruct)(**in, sm.m[nm]);
	}
	for (nm = 0; nm < mm.n; nm++) {
		if (!search(**in, "MM", nm, buffer)) {
			(*in)->close(); delete (*in); (*in) = NULL;
			sprintf(errmsg, "can not find MM %d", nm); show_msg(errmsg);
			return false;
		}
		//ReadMoleculeStruct(**in, mm.m[nm]);
		((READ_MM_STRUCT)read_mmstruct)(**in, mm.m[nm]);
	}
	return true;
}

bool get_multi_mol_md_Energy(ifstream** Ein, char *title, int ifile, int ncur_file, int isave, int &nloop, ARRAY<float> &Ek, float& Ep) {
	if ((*Ein = open_multi_mol_md_file(*Ein, ncur_file, title, "E", ifile)) == NULL) {show_msg(errmsg); return false;}

	char buffer[256] = "\0", buff[20] = "\0";
	float ek = 0, ep = 0;
	if (!search(**Ein, "[MD]", isave, buffer)) {
		sprintf(errmsg, "Can not find MD Proc. # %d", isave);
		show_msg(errmsg); delete *Ein; *Ein = NULL; return false;
	}
	sscanf(buffer, "%s %d %d", buff, &isave, &nloop);

	int i = 0;
	for (i = 0; i < Ek.n; i++) {
		(**Ein).getline(buffer, 250);
		if (sscanf(buffer, "%f", &ek) != 1) {
			sprintf(errmsg, "Can not get energy from %s", buffer);
			show_msg(errmsg); delete *Ein; *Ein = NULL; return false;
		}
		else Ek.m[i] = ek;
	}
	(**Ein).getline(buffer, 250);
	if (sscanf(buffer, "%f", &ep) != 1) {
		sprintf(errmsg, "Can not get total potential energy from %s", buffer);
		show_msg(errmsg); delete *Ein; *Ein = NULL; return false;
	}
	Ep = ep;
	return true;
}

bool get_multi_mol_md_NPT_PV(ifstream** pv_in, char *title, int ifile, int ncur_file, int isave, int &nloop, float &P, float& V, float &eps) {
	char pv_title[256] = "\0";
	sprintf(pv_title, "%s-pv", title);
	if ((*pv_in = open_multi_mol_md_file(*pv_in, ncur_file, pv_title, "", ifile)) == NULL) {show_msg(errmsg); return false;}

	char buffer[256] = "\0", buff[20] = "\0";
	if (!search(**pv_in, "[MD]", isave, buffer)) {
		sprintf(errmsg, "Can not find MD Proc. # %d", isave);
		show_msg(errmsg); delete *pv_in; *pv_in = NULL; return false;
	}
	sscanf(buffer, "%s %d %d", buff, &isave, &nloop);

	(**pv_in).getline(buffer, 250);
	if (sscanf(buffer, "%f %f", &V, &P) != 2) {
		sprintf(errmsg, "Can not get V & P from %s", buffer);
		show_msg(errmsg); delete *pv_in; *pv_in = NULL; return false;
	}
	(**pv_in).getline(buffer, 250);
	if (sscanf(buffer, "%f", &eps) != 1) {
		sprintf(errmsg, "Can not get eps from %s", buffer);
		show_msg(errmsg); delete *pv_in; *pv_in = NULL; return false;
	}
	return true;
}

template <class PM, class SM, class MM> bool get_multi_mol_md_proc(char *mdpar_fname, long isave, ARRAY<PM> &pm, ARRAY<SM> &sm, ARRAY<MM> &mm, void *read_pmstruct, void *read_smstruct, void *read_mmstruct) {
	if (mm_md_cell.mm.n <= 0 && mm_md_cell.sm.n <= 0 && mm_md_cell.pm.n <= 0) return true;

	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {
		show_msg(::errmsg);
		return false;
	}
	int ifile = isave / msave.max_save;
	int nloop = 0;

	ifstream *in = NULL;
	//bool status = get_multi_mol_md_confg<SM, MM>(&in, msave.title, ifile, -1, isave, nloop, mm_md_cell.sm, mm_md_cell.mm);
	bool status = get_multi_mol_md_confg<PM, SM, MM>(&in, msave.title, ifile, -1, isave, nloop, pm, sm, mm, read_pmstruct, read_smstruct, read_mmstruct);
	if (in != NULL) {if (in->is_open()) in->close(); delete in; in = NULL;}
	return status;
}

bool get_NPT_PV_proc(char *mdpar_fname, long isave, float &P, float &V) {
	char buffer[256] = "\0", pv_title[256] = "\0";
	if (read_msg(mdpar_fname, buffer, "[MD_MODE]") && sscanf(buffer, "%d", &MD_ENSEMBLE) == 1) {
		if (MD_ENSEMBLE != MD_NPT) MD_ENSEMBLE = MD_NVT;
	}
	else MD_ENSEMBLE = MD_NVT;

	if (MD_ENSEMBLE != MD_NPT) return false;

	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {
		show_msg(::errmsg);
		return false;
	}
	int ifile = isave / msave.max_save;
	int nloop = 0;

	float veps;
	ifstream *pv_in = NULL;
	if (!get_multi_mol_md_NPT_PV(&pv_in, msave.title, ifile, -1, isave, nloop, P, V, veps)) return false;
	else return true;
}

void mdcell_size_update(float &x, float &y, float &z, float V) {
	float V0 = x * y * z;
	float veps = powf(V / V0, 1.0/3); 
	x *= veps; y *= veps; z *= veps;
}

bool get_multi_mm_md_proc(char *mdpar_fname, int isave) {
	bool status = get_multi_mol_md_proc<PMOLECULE, SMOLECULE, MMOLECULE>(mdpar_fname, isave, mm_md_cell.pm, mm_md_cell.sm, mm_md_cell.mm, (void*)(&ReadPointMoleculeStruct), (void*)(&ReadSolidMoleculeStruct), (void*)(&ReadMoleculeStruct));
	if (!status) return false;

	float P = 0, V = 0;
	if (MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
		if (MD_ENSEMBLE == MD_NPT) {
			if (!get_NPT_PV_proc(mdpar_fname, isave, P, V)) return false;
			mdcell_size_update(mm_md_cell.h[0], mm_md_cell.h[1], mm_md_cell.h[2], V / 8);
		}
		molecule_check_periodic_cell(mm_md_cell);
	}
#if _SYS_ == _WINDOWS_SYS_
	if (pMainFrame != NULL) {
		if (status) ::PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
	}
#endif
	return true;
}

bool get_multi_cgmm_md_proc(char *mdpar_fname, int isave) {
	bool status = get_multi_mol_md_proc<CGSM, CGSM, CG_MMOLECULE>(mdpar_fname, isave, cgmm_md_cell.pm, cgmm_md_cell.sm, cgmm_md_cell.mm, (void*)(&ReadCGSMStruct), (void*)(&ReadCGSMStruct), (void*)(&ReadCGMMStruct));
	if (!status) return false;
	if (MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
		molecule_check_periodic_cell(mm_md_cell);
	}
#if _SYS_ == _WINDOWS_SYS_
	if (pMainFrame != NULL) {
		if (status) ::PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
	}
#endif
	return true;
}

extern void molecule_check_periodic_cell(MMOL_MD_CELL &mmcell);


bool show_multi_mm_md_proc(char *mdpar_fname, int nf) {
	if (mm_md_cell.mm.n <= 0 && mm_md_cell.sm.n <= 0 && mm_md_cell.pm.n <= 0) return true;
	if (nf < 0) nf = 0;

	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {
		show_msg(::errmsg); return false;
	}

	ARRAY<float> Ek;
	float Ep = 0;
	Ek.SetArray(mm_md_cell.mm.n + mm_md_cell.sm.n + mm_md_cell.pm.n);

	int ifile = 0, nfile = 0;
	int isave = nf, nloop = 0;
	char buffer[256] = "\0", buff[20] = "\0";
	
	ifstream *in = NULL;
	ifstream *Ein = NULL;
	ifstream *pv_in = NULL;

	bool proc = true;
	int n = 0, nm = 0;
	long LOOPS = 0;

	ifstream par_in;
	par_in.open(mdpar_fname);
	if (search(par_in, "[MD_CELL]", buffer)) {
		par_in.getline(buffer, 250);
		sscanf(buffer, "%d", &(MD_MODE));
	}
	if (search(par_in, "[MD_MODE]", buffer)) {
		if (sscanf(buffer, "%d", &MD_ENSEMBLE) != 1) MD_ENSEMBLE = MD_NVT;
		else if (MD_ENSEMBLE != MD_NVT && MD_ENSEMBLE != MD_NPT) MD_ENSEMBLE = MD_NVT;
	}
	if (search(par_in, "[TOTAL_LOOPS]", buffer)) {
		par_in.getline(buffer, 250);
		sscanf(buffer, "%ld", &LOOPS);
		LOOPS *= msave.nsave;
		par_in.close();
	}
	else {
		par_in.close();
		sprintf(errmsg, "Can not get [TOTAL_LOOPS] in %s", mdpar_fname);
		show_msg(errmsg); return false;
	}

	float Ek_total = 0, Ep_total = 0;
	float P = 0, V = 0, veps = 0;
	char pv_title[256] = "\0";
	sprintf(pv_title, "%s-pv", msave.title);
	bool success = true;

	COMMAND_MD_STOP = false;

	ifile = 0; nfile = -1;
	while (proc) {
		if (COMMAND_MD_STOP) break; // stop with given command
		ifile = isave / msave.max_save;

		if (MD_ENSEMBLE == MD_NPT) {
			if (!get_multi_mol_md_NPT_PV(&pv_in, pv_title, ifile, nfile, isave, nloop, P, V, veps)) {
				proc = false; success = false; break;
			}
			else {
				mdcell_size_update(mm_md_cell.h[0], mm_md_cell.h[1], mm_md_cell.h[2], V / 8);
			}
		}

		if (!get_multi_mol_md_confg<PMOLECULE, SMOLECULE, MMOLECULE>(&in, msave.title, ifile, nfile, isave, nloop, mm_md_cell.pm, mm_md_cell.sm, mm_md_cell.mm, (void*)(&ReadPointMoleculeStruct), (void*)(&ReadSolidMoleculeStruct), (void*)(&ReadMoleculeStruct)) ||
			!get_multi_mol_md_Energy(&Ein, msave.title, ifile, nfile, isave, nloop, Ek, Ep)) {
			proc = false; success = false; break;
		}
		nfile = ifile;

		if (MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
			molecule_check_periodic_cell(mm_md_cell);
		}

		Ek_total = 0; Ep_total = 0;
		for (nm = 0; nm < mm_md_cell.mm.n + mm_md_cell.sm.n + mm_md_cell.pm.n; nm++) {
			Ek_total += Ek.m[nm]; Ep_total = Ep;
		}

		show_loop_infor(nloop, isave, LOOPS - nloop, Ek_total, Ep_total);

	#if _SYS_ == _WINDOWS_SYS_
		if (pMainFrame != NULL) {
			if (delay_time > 0) WaitForSingleObject(hWaitEvent, delay_time);
			ResetEvent(hWaitGetMolStructEvent);
			PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
			WaitForSingleObject(hWaitGetMolStructEvent, max_time_pickup_molstruct);
		}
	#endif

		isave += nshow;
	}

	if (in != NULL) {if (in->is_open()) in->close(); delete in; in = NULL;}
	if (Ein != NULL) {if (Ein->is_open()) Ein->close(); delete Ein; Ein = NULL;}
	Ek.release();
	return success;
}

extern void molecule_check_periodic_cell(CGMM_MD_CELL &mmcell);

bool show_multi_cgmm_md_proc(char *mdpar_fname, int nf) {
	if (cgmm_md_cell.mm.n <= 0 && cgmm_md_cell.sm.n <= 0 && cgmm_md_cell.pm.n <= 0) return true;
	if (nf < 0) nf = 0;

	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {
		show_msg(::errmsg); return false;
	}

	ARRAY<float> Ek;
	float Ep = 0;
	Ek.SetArray(cgmm_md_cell.mm.n + cgmm_md_cell.sm.n + cgmm_md_cell.pm.n);

	int ifile = 0, nfile = 0;
	int isave = nf, nloop = 0;
	char buffer[256] = "\0", buff[20] = "\0";
	
	ifstream *in = NULL;
	ifstream *Ein = NULL;

	bool proc = true;
	int n = 0, nm = 0;
	long LOOPS = 0;
	int MD_MODE = 0;

	ifstream par_in;
	par_in.open(mdpar_fname);
	if (search(par_in, "[MD_CELL]", buffer)) {
		par_in.getline(buffer, 250);
		sscanf(buffer, "%d", &(MD_MODE));
	}
	if (search(par_in, "[MD_MODE]", buffer)) {
		if (sscanf(buffer, "%d", &MD_ENSEMBLE) != 1) MD_ENSEMBLE = MD_NVT;
		else if (MD_ENSEMBLE != MD_NVT && MD_ENSEMBLE != MD_NPT) MD_ENSEMBLE = MD_NVT;
	}
	if (search(par_in, "[TOTAL_LOOPS]", buffer)) {
		par_in.getline(buffer, 250);
		sscanf(buffer, "%ld", &LOOPS);
		LOOPS *= msave.nsave;
		par_in.close();
	}
	else {
		par_in.close();
		sprintf(errmsg, "Can not get [TOTAL_LOOPS] in %s", mdpar_fname);
		show_msg(errmsg); return false;
	}

	float Ek_total = 0, Ep_total = 0;
	bool success = true;

	COMMAND_MD_STOP = false;

	ifile = 0; nfile = -1;
	while (proc) {
		if (COMMAND_MD_STOP) break; // stop with given command
		ifile = isave / msave.max_save;

		if (!get_multi_mol_md_confg<CGSM, CGSM, CG_MMOLECULE>(&in, msave.title, ifile, nfile, isave, nloop, cgmm_md_cell.pm, cgmm_md_cell.sm, cgmm_md_cell.mm, (void*)(&ReadCGSMStruct), (void*)(&ReadCGSMStruct), (void*)(&ReadCGMMStruct)) ||
			!get_multi_mol_md_Energy(&Ein, msave.title, ifile, nfile, isave, nloop, Ek, Ep)) {
			proc = false; success = false; break;
		}
		nfile = ifile;

		if (MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			_coarse_grain_::molecule_check_periodic_cell(cgmm_md_cell);
		}

		Ek_total = 0; Ep_total = 0;
		for (nm = 0; nm < cgmm_md_cell.mm.n + cgmm_md_cell.sm.n; nm++) {
			Ek_total += Ek.m[nm]; Ep_total = Ep;
		}

		show_loop_infor(nloop, isave, LOOPS - nloop, Ek_total, Ep_total);

	#if _SYS_ == _WINDOWS_SYS_
		if (pMainFrame != NULL) {
			if (delay_time > 0) WaitForSingleObject(hWaitEvent, delay_time);
			ResetEvent(hWaitGetMolStructEvent);
			PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
			WaitForSingleObject(hWaitGetMolStructEvent, max_time_pickup_molstruct);
		}
	#endif

		isave += nshow;
	}

	if (in != NULL) {if (in->is_open()) in->close(); delete in; in = NULL;}
	if (Ein != NULL) {if (Ein->is_open()) Ein->close(); delete Ein; Ein = NULL;}
	Ek.release();
	return success;
}

void ClusterSaveXYZ(ofstream *out, BASIC_CLUSTER *pc, bool bInCMM, double *dv) {
	if (out == NULL) return;
	if (!out->is_open()) return;
	bool bSave = false;
	VECTOR3 r;
	int i;
	for (i = 0; i < pc->nAtoms; i++) {
		if (bInCMM) {V3plusV3(pc->atom[i].r, pc->dr_cmm, r)}
		else {V32V3(pc->atom[i].r, r)}
		if (dv != NULL) {r.v[0] += dv[0]; r.v[1] += dv[1]; r.v[2] += dv[2];}
		(*out)<<pc->atom[i].par->atom<<" "<<r.v[0]<<" "<<r.v[1]<<" "<<r.v[2]<<endl;
	}
	return;
}

void MMSaveXYZ(ofstream *out, MMOLECULE *mm, ARRAY<short> &id_monomer, ARRAY<short> &id_cluster, bool bInCMM, double *dv) {
	if (out == NULL) return;
	if (!out->is_open()) return;
	int imonomer = 0, icluster = 0, i;
	_GROUP<CLUSTER> *pmonomer = NULL;
	CLUSTER *pc = NULL;
	bool bSave = false;
	VECTOR3 r;
	for (imonomer = 0; imonomer < mm->nMonomer; imonomer++) {
		pmonomer = mm->monomer + imonomer;
		bSave = (id_monomer.m == NULL || id_monomer.exist(pmonomer->uid) ? true : false);
		if (!bSave) continue;
		for (icluster = 0; icluster < pmonomer->nUnit; icluster++) {
			pc = pmonomer->u[icluster];
			bSave = (id_cluster.m == NULL || id_cluster.exist(pc->cID) ? true : false);
			if (!bSave) continue;
			for (i = 0; i < pc->nAtoms; i++) {
				if (bInCMM) {
					V3plusV3(pc->atom[i].r, pc->dr_cmm, r)
					if (dv != NULL) {r.v[0] += dv[0]; r.v[1] += dv[1]; r.v[2] += dv[2];}
				}
				else {
					V32V3(pc->atom[i].r, r)
				}
				(*out)<<pc->atom[i].par->atom<<" "<<r.v[0]<<" "<<r.v[1]<<" "<<r.v[2]<<endl;
			}
		}
	}
	return;
}

namespace _coarse_grain_ {

void CGClusterSaveXYZ(ofstream *out, CG_CLUSTER *pc, bool bInCMM, double *dv) {
	if (out == NULL) return;
	if (!out->is_open()) return;
	bool bSave = false;
	VECTOR3 r;
	int i;
	for (i = 0; i < pc->nAtoms; i++) {
		if (bInCMM) {V3plusV3(pc->atom[i].r, pc->dr_cmm, r)}
		else {V32V3(pc->atom[i].r, r)}
		if (dv != NULL) {r.v[0] += dv[0]; r.v[1] += dv[1]; r.v[2] += dv[2];}
		(*out)<<pc->atom[i].par->atom<<" "<<r.v[0]<<" "<<r.v[1]<<" "<<r.v[2]<<endl;
	}
	return;
}

void CGMMSaveXYZ(ofstream *out, CG_MMOLECULE *mm, ARRAY<short> &id_cgatom, bool bInCMM, double *dv) {
	if (out == NULL) return;
	if (!out->is_open()) return;
	int icluster = 0, i = 0;
	CG_CLUSTER *pc = NULL;
	bool bSave = false;
	VECTOR3 r;
	for (icluster = 0; icluster < mm->nCluster; icluster++) {
		pc = mm->cluster + icluster;
		bSave = (id_cgatom.m == NULL || id_cgatom.exist(pc->cID) ? true : false);
		if (!bSave) continue;
		for (i = 0; i < pc->nAtoms; i++) {
			if (bInCMM) {
				V3plusV3(pc->atom[i].r, pc->dr_cmm, r)
				if (dv != NULL) {r.v[0] += dv[0]; r.v[1] += dv[1]; r.v[2] += dv[2];}
			}
			else {
				V32V3(pc->atom[i].r, r)
			}
			(*out)<<pc->atom[i].par->pdb_alias<<" "<<r.v[0]<<" "<<r.v[1]<<" "<<r.v[2]<<endl;
		}
	}
	return;
}

} // end of namespace _coarse_grain_

void SMSaveXYZ(ofstream *out, SMOLECULE *sm, ARRAY<short> &id_cluster) {
	if (out == NULL) return;
	if (!out->is_open()) return;
	BASIC_CLUSTER *pc = sm->c;
	bool bSave = false;
	bSave = (id_cluster.m == NULL || id_cluster.exist(pc->cID) ? true : false);
	if (!bSave) return;
	for (int i = 0; i < pc->nAtoms; i++) {
		(*out)<<pc->atom[i].par->atom<<" "<<pc->atom[i].r.v[0]<<" "<<pc->atom[i].r.v[1]<<" "<<pc->atom[i].r.v[2]<<endl;
	}
	return;
}

void PMSaveXYZ(ofstream *out, PMOLECULE *pm, ARRAY<short> &id_cluster) {
	if (out == NULL) return;
	if (!out->is_open()) return;
	BASIC_CLUSTER *pc = pm->c;
	bool bSave = false;
	bSave = (id_cluster.m == NULL || id_cluster.exist(pc->cID) ? true : false);
	if (!bSave) return;
	for (int i = 0; i < pc->nAtoms; i++) {
		(*out)<<pc->atom[i].par->atom<<" "<<pc->atom[i].r.v[0]<<" "<<pc->atom[i].r.v[1]<<" "<<pc->atom[i].r.v[2]<<endl;
	}
	return;
}

namespace _coarse_grain_ {
void CGSMSaveXYZ(ofstream *out, CGSM *sm, ARRAY<short> &id_cgatom) {
	if (out == NULL) return;
	if (!out->is_open()) return;
	bool bSave = false;
	bSave = (id_cgatom.m == NULL || id_cgatom.exist(sm->cID) ? true : false);
	if (!bSave) return;
	for (int i = 0; i < sm->nAtoms; i++) {
		(*out)<<sm->atom[i].par->pdb_alias<<" "<<sm->atom[i].r.v[0]<<" "<<sm->atom[i].r.v[1]<<" "<<sm->atom[i].r.v[2]<<endl;
	}
	return;
}

} // end of namespace _coarse_grain_
