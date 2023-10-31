#include "TAMD\project.h"
//#include "project.h"

#if _SYS_ == _WINDOWS_SYS_
#include "stdafx.h"
#endif

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#include <fftw3.h>

using namespace std;

#define _READ_MM_MD_  1

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
#include "TAMD\fmd_cell.h"
#include "TAMD\cg-md-cell.h"
#include "TAMD\cg-mpi-md-cell.h"
#include "TAMD\test.h"
#include "TAMD\read.h"
#include "TAMD\var.h"
#include "TAMD\ReadMDProc.h"
#include "TAMD\mol-cell.h"
#include "TAMD\show_md.h"
#include "TAMD\mdcell_nhc.h"
#include "TAMD\NPT_MD.h"
#include "TAMD\md_module.h"
#include "TAMD\sr.h"
#include "export1.h"

#include "TAMD\vmd_cell.h"
#include "TAMD\interface.h"

#include "TAMD\complex.h"

#define _READ_MM_MD_ 1

#if _MPI_ == 1
#include "TAMD\mpi_md.h"
#include "TAMD\mpi_md_cell.h"
#endif

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
#include "fmd_cell.h"
#include "cg-md-cell.h"
#include "cg-mpi-md-cell.h"
#include "test.h"
#include "read.h"
#include "var.h"
#include "ReadMDProc.h"
#include "mol-cell.h"
#include "show_md.h"
#include "mdcell_nhc.h"
#include "NPT_MD.h"
#include "md_module.h"
#include "sr.h"
#include "export1.h"

#include "vmd_cell.h"
#include "interface.h"

#include "complex.h"

#if _MPI_ == 1
#include "mpi_md.h"
#include "mpi_md_cell.h"
#endif

#endif

#if _MPI_ == 1
#ifdef SEEK_SET
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include "mpi.h"
#endif
#endif

using namespace _coarse_grain_;

MMOLECULE mm;
CG_MMOLECULE cgmm;
MMOL_MD_CELL mm_md_cell;
CGMM_MD_CELL cgmm_md_cell;
FLEX_MMOL_CELL flex_mmcell; 
MMOL_VMD_CELL vmdcell;
LIST<MMOLECULE> *mlist = NULL;
LIST<CG_MMOLECULE> *cgmlist = NULL;

int JOB = UNKNOWN_JOB;
int MD_MODE = SINGLE_MM_FREE_CELL_MD;
int MD_ENSEMBLE = MD_NVT;

bool md_run = false;

#if _SYS_ == _WINDOWS_SYS_
CWnd *pMainFrame = NULL;
HANDLE hWaitGetMolStructEvent = NULL;
HANDLE hWaitEvent = NULL;
#endif

int delay_time = 2; // in milisecond
int max_time_pickup_molstruct = 100;  // in milisecond
int nshow = 1; // the molecule is shown every nshow loops
int ishow = 0;

extern void molecule_check_periodic_cell(MMOL_MD_CELL &mmcell);

extern PARS_DATABASE< COMM_PARS<50> > cluster_db;

#if _SYS_ == _WINDOWS_SYS_
void set_main_frame(CWnd* pMainFrame) {
	::pMainFrame = pMainFrame;
}
#endif

#if _MPI_ == 1
void reset_mpi_pars() {
	if (::n1MM != NULL) {delete[] ::n1MM; ::n1MM = NULL;}
	if (::n2MM != NULL) {delete[] ::n2MM; ::n2MM = NULL;}
	if (::n1SM != NULL) {delete[] ::n1SM; ::n1SM = NULL;}
	if (::n2SM != NULL) {delete[] ::n2SM; ::n2SM = NULL;}
	//if (::n1CMM != NULL) {delete[] ::n1CMM; ::n1CMM = NULL;}
	//if (::n2CMM != NULL) {delete[] ::n2CMM; ::n2CMM = NULL;}
}

void set_mpi(bool bMM, bool bSM, bool bCMM) {
	reset_mpi_pars(); if (::nMPI <= 0) return; // this can not happen
	if (bMM) {::n1MM = new int[::nMPI]; ::n2MM = new int[::nMPI];}
	if (bSM) {::n1SM = new int[::nMPI]; ::n2SM = new int[::nMPI];}
	//if (bCMM) {::n1CMM = new int[::nMPI]; ::n2CMM = new int[::nMPI];}
}

void MPI_SingleMM_JobAssign(int nCluster) {
	char msg[256] = "\0";
	sprintf(msg, "Macromolecule's clusters are distributed to :"); show_log(msg, true); 
	set_mpi(true, false, false);
	int dn = int(float(nCluster) / ::nMPI + 1.f / ::nMPI);
	for (int i = 0; i < ::nMPI; i++) {
		::n1MM[i] = i * dn;
		if (i == ::nMPI - 1) ::n2MM[i] = nCluster - 1;
		else ::n2MM[i] = ::n1MM[i] + dn - 1;
		sprintf(msg, "proc %d: [%d, %d]", i, n1MM[i], n2MM[i]); show_log(msg, true);
	}
}

bool MPI_MultiMM_JobAssign_Auto(int nMM, int nSM) {
	char msg[256] = "\0";
	sprintf(msg, "Automatically distribute molecules to different processors"); show_log(msg, true); 
	sprintf(msg, "molecules are distributed to :"); show_log(msg, true); 
	bool bMM = (nMM > 0 ? true : false);
	bool bSM = (nSM > 0 ? true : false);
	set_mpi(bMM, bSM, false);
	int dnMM = int(float(nMM) / ::nMPI + 1.0f / ::nMPI), dnSM = int(float(nSM) / ::nMPI + 1.0f / ::nMPI);
	for (int i = 0; i < ::nMPI; i++) {
		if (bMM) ::n1MM[i] = i * dnMM;
		if (bSM) ::n1SM[i] = i * dnSM;
		if (i == ::nMPI - 1) {
			if (bMM) ::n2MM[i] = nMM - 1;
			if (bSM) ::n2SM[i] = nSM - 1;
		}
		else {
			if (bMM) ::n2MM[i] = ::n1MM[i] + dnMM - 1;
			if (bSM) ::n2SM[i] = ::n1SM[i] + dnSM - 1;
		}
		if (bMM) {
			sprintf(msg, "%s proc %d: macro-molecules [%d, %d]", ::proc_name[i], i, n1MM[i], n2MM[i]); show_log(msg, false);
		}
		if (bSM) {
			sprintf(msg, ", small-molecules [%d, %d]", i, n1SM[i], n2SM[i]); show_log(msg, true);
		}
		else show_log("", true);
	}
	return true;
}

bool MPI_MultiMM_JobAssign_FromFile(int nMM, int nSM) {
	char msg[256] = "\0", buffer[256] = "\0";
	char job_fname[256] = "job-assign";
	ifstream in;
	in.open(job_fname);
	if (!in.is_open()) {
		sprintf(msg, "file %s does not exist, to assigne job on each processor", job_fname); show_log(msg, true); return false;
	}
	int nlen = 0;
	sprintf(msg, "molecules are distributed to :"); show_log(msg, true); 
	bool bMM = (nMM > 0 ? true : false);
	bool bSM = (nSM > 0 ? true : false);
	set_mpi(bMM, bSM, false);
	int n1MM = 0, n2MM = 0, n1SM = 0, n2SM = 0;
	for (int i = 0; i < ::nMPI; i++) {
		rewind(in); 
		if (!search(in, ::proc_name[i], buffer)) {
			sprintf(msg, "job for processor %s is not given in %s", ::proc_name[i], job_fname); show_log(msg, true); 
			in.close(); return false;
		}
		if (sscanf(buffer, "%s %d %d %d %d", msg, &n1MM, &n2MM, &n1SM, &n2SM) != 5) {
			sprintf(msg, "job assignment is not correct [process name] [n1MM] [n2MM] [n1SM] [n2SM] : %s", buffer); show_log(msg, true);
			in.close(); return false;
		}
		if (bMM) {::n1MM[i] = n1MM; ::n2MM[i] = n2MM;}
		if (bSM) {::n1SM[i] = n1SM; ::n2SM[i] = n2SM;}
		if (bMM) {
			sprintf(msg, "%s proc %d: macro-molecules [%d, %d]", ::proc_name[i], i, ::n1MM[i], ::n2MM[i]); show_log(msg, false);
		}
		if (bSM) {
			sprintf(msg, ", small-molecules [%d, %d]", i, ::n1SM[i], ::n2SM[i]); show_log(msg, true);
		}
		else show_log("", true);
	}
	in.close(); return true;
}

bool MPI_MultiMM_JobAssign(int nMM, int nSM) {
	if (!MPI_MultiMM_JobAssign_FromFile(nMM, nSM)) return MPI_MultiMM_JobAssign_Auto(nMM, nSM);
	else return true;
}

#endif

extern "C" void set_md_mode(int mode) {MD_MODE = mode; JOB = MD_MODE;}
extern "C" void set_job(int job) {JOB = job;}
// get the molecule in current job
bool get_molecule_upto_job(int nMol, GENERAL_MOLECULE& gm) {
	int n = 0;
	LIST<MMOLECULE> *plist = NULL;
	LIST<CG_MMOLECULE> *pcglist = NULL;

	switch (JOB) {
	case GET_SINGLE_MM_MD_PROC:
		if (nMol == 0) {gm.set_mmol(&mm); return true;}
		else return false;
		break;
	case SHOW_SINGLE_MM_MD_PROC:
		if (nMol == 0) {gm.set_mmol(&(mm)); return true;}
		else return false;
		break;
	case GET_MULTI_MM_MD_PROC:
		if (nMol < mm_md_cell.mm.n && nMol >= 0) {
			gm.set_mmol(mm_md_cell.mm.m + nMol);
			return true;
		}
		else if (nMol >= 0 && nMol < mm_md_cell.mm.n + mm_md_cell.sm.n) {
			gm.set_smol(mm_md_cell.sm.m + nMol - mm_md_cell.mm.n);
			return true;
		}
		else if (nMol >= 0 && nMol < mm_md_cell.mm.n + mm_md_cell.sm.n + mm_md_cell.pm.n) {
			gm.set_pmol(mm_md_cell.pm.m + nMol - mm_md_cell.mm.n - mm_md_cell.sm.n);
			return true;
		} 
		else return false;
		break;
	case SHOW_MULTI_MM_MD_PROC:
		if (nMol < mm_md_cell.mm.n && nMol >= 0) {
			gm.set_mmol(mm_md_cell.mm.m + nMol);
			return true;
		}
		else if (nMol >= 0 && nMol < mm_md_cell.mm.n + mm_md_cell.sm.n) {
			gm.set_smol(mm_md_cell.sm.m + nMol - mm_md_cell.mm.n);
			return true;
		}
		else if (nMol >= 0 && nMol < mm_md_cell.mm.n + mm_md_cell.sm.n + mm_md_cell.pm.n) {
			gm.set_pmol(mm_md_cell.pm.m + nMol - mm_md_cell.mm.n - mm_md_cell.sm.n);
			return true;
		} 
		else return false;
		break;
	case SINGLE_MM_FREE_CELL_MD:
		if (nMol == 0) {gm.set_mmol(&mm); return true;}
		else return false;
		break;
	case MULTI_MM_FREE_CELL_MD:
		if (nMol < mm_md_cell.mm.n && nMol >= 0) {
			gm.set_mmol(mm_md_cell.mm.m + nMol);
			return true;
		}
		else if (nMol >= 0 && nMol < mm_md_cell.mm.n + mm_md_cell.sm.n) {
			gm.set_smol(mm_md_cell.sm.m + nMol - mm_md_cell.mm.n);
			return true;
		}
		else if (nMol >= 0 && nMol < mm_md_cell.mm.n + mm_md_cell.sm.n + mm_md_cell.pm.n) {
			gm.set_pmol(mm_md_cell.pm.m + nMol - mm_md_cell.mm.n - mm_md_cell.sm.n);
			return true;
		} 
		else return false;
		break;
	case MULTI_MM_PERIDIC_CELL_MD:
		if (nMol < mm_md_cell.mm.n && nMol >= 0) {
			gm.set_mmol(mm_md_cell.mm.m + nMol);
			return true;
		}
		else if (nMol >= 0 && nMol < mm_md_cell.mm.n + mm_md_cell.sm.n) {
			gm.set_smol(mm_md_cell.sm.m + nMol - mm_md_cell.mm.n);
			return true;
		}
		else if (nMol >= 0 && nMol < mm_md_cell.mm.n + mm_md_cell.sm.n + mm_md_cell.pm.n) {
			gm.set_pmol(mm_md_cell.pm.m + nMol - mm_md_cell.mm.n - mm_md_cell.sm.n);
			return true;
		} 
		else return false;
		break;

		/*
		#define SINGLE_CGMM_FREE_CELL_MD   3
		#define MULTI_CGMM_FREE_CELL_MD    4
		#define MULTI_CGMM_PERIDIC_CELL_MD 5
		*/
	case SINGLE_CGMM_FREE_CELL_MD:
		if (nMol == 0) {gm.set_cgmm(&cgmm); return true;}
		else return false;
		break;
	case GET_SINGLE_CGMM_MD_PROC:
		if (nMol == 0) {gm.set_cgmm(&cgmm); return true;}
		else return false;
		break;
	case SHOW_SINGLE_CGMM_MD_PROC:
		if (nMol == 0) {gm.set_cgmm(&(cgmm)); return true;}
		else return false;
		break;
	case MULTI_CGMM_FREE_CELL_MD:
		if (nMol < cgmm_md_cell.mm.n && nMol >= 0) {
			gm.set_cgmm(cgmm_md_cell.mm.m + nMol);
			return true;
		}
		else if (nMol >= 0 && nMol < cgmm_md_cell.mm.n + cgmm_md_cell.sm.n) {
			gm.set_cgsm(cgmm_md_cell.sm.m + nMol - cgmm_md_cell.mm.n);
			return true;
		}
		else return false;
		break;
	case MULTI_CGMM_PERIDIC_CELL_MD:
		if (nMol < cgmm_md_cell.mm.n && nMol >= 0) {
			gm.set_cgmm(cgmm_md_cell.mm.m + nMol);
			return true;
		}
		else if (nMol >= 0 && nMol < cgmm_md_cell.mm.n + cgmm_md_cell.sm.n) {
			gm.set_cgsm(cgmm_md_cell.sm.m + nMol - cgmm_md_cell.mm.n);
			return true;
		}
		else return false;
		break;
	case GET_MULTI_CGMM_MD_PROC:
		if (nMol < cgmm_md_cell.mm.n && nMol >= 0) {
			gm.set_cgmm(cgmm_md_cell.mm.m + nMol);
			return true;
		}
		else if (nMol >= 0 && nMol < cgmm_md_cell.mm.n + cgmm_md_cell.sm.n) {
			gm.set_cgsm(cgmm_md_cell.sm.m + nMol - cgmm_md_cell.mm.n);
			return true;
		}
		else return false;
		break;
	case SHOW_MULTI_CGMM_MD_PROC:
		if (nMol < cgmm_md_cell.mm.n && nMol >= 0) {
			gm.set_cgmm(cgmm_md_cell.mm.m + nMol);
			return true;
		}
		else if (nMol >= 0 && nMol < cgmm_md_cell.mm.n + cgmm_md_cell.sm.n) {
			gm.set_cgsm(cgmm_md_cell.sm.m + nMol - cgmm_md_cell.mm.n);
			return true;
		}
		else return false;
		break;
	case MOLECULE_OPERATION:
		plist = mlist; n = 0;
		if (plist != NULL) {
			while (n < nMol) {
				if (plist == NULL) break;
				plist = plist->next; n++;
			}
			if (n == nMol && plist != NULL) {
				gm.set_mmol(plist->p); return true;
			}
		}

		pcglist = cgmlist; 
		while (n < nMol) {
			if (pcglist == NULL) break;
			pcglist = pcglist->next; n++;
		}
		if (pcglist == NULL) return false;
		else {gm.set_cgmm(pcglist->p); return true;}

		break;

	case GCMD:
		if (nMol < vmdcell.mm.n && nMol >= 0) {
			gm.set_mmol(vmdcell.mm.m[nMol]);
			return true;
		}
		else if (nMol >= 0 && nMol < vmdcell.mm.n + vmdcell.sm.n) {
			gm.set_smol(vmdcell.sm.m[nMol - vmdcell.mm.n]);
			return true;
		}
		else if (nMol >= 0 && nMol < vmdcell.mm.n + vmdcell.sm.n + vmdcell.pm.n) {
			gm.set_pmol(vmdcell.pm.m[nMol - vmdcell.mm.n - vmdcell.sm.n]);
			return true;
		} 
		else return false;
		break;
	case GET_GCMD_PROC:
		if (nMol < vmdcell.mm.n && nMol >= 0) {
			gm.set_mmol(vmdcell.mm.m[nMol]);
			return true;
		}
		else if (nMol >= 0 && nMol < vmdcell.mm.n + vmdcell.sm.n) {
			gm.set_smol(vmdcell.sm.m[nMol - vmdcell.mm.n]);
			return true;
		}
		else if (nMol >= 0 && nMol < vmdcell.mm.n + vmdcell.sm.n + vmdcell.pm.n) {
			gm.set_pmol(vmdcell.pm.m[nMol - vmdcell.mm.n - vmdcell.sm.n]);
			return true;
		} 
		else return false;
		break;
	case SHOW_GCMD_PROC:
		if (nMol < vmdcell.mm.n && nMol >= 0) {
			gm.set_mmol(vmdcell.mm.m[nMol]);
			return true;
		}
		else if (nMol >= 0 && nMol < vmdcell.mm.n + vmdcell.sm.n) {
			gm.set_smol(vmdcell.sm.m[nMol - vmdcell.mm.n]);
			return true;
		}
		else if (nMol >= 0 && nMol < vmdcell.mm.n + vmdcell.sm.n + vmdcell.pm.n) {
			gm.set_pmol(vmdcell.pm.m[nMol - vmdcell.mm.n - vmdcell.sm.n]);
			return true;
		} 
		else return false;
		break;
	default:
		return false;
	}
}

int get_macromol_atom_number(MMOLECULE *mm) {
	int n = 0, nc = 0;
	for (nc = 0; nc < mm->nCluster; nc++) {
		n += mm->cluster[nc].nAtoms;
	}
	return n;
}

int get_macromol_atom_number(CG_MMOLECULE *cgmm) {
	int n = 0, nc = 0;
	for (nc = 0; nc < cgmm->nCluster; nc++) {
		n += cgmm->cluster[nc].nAtoms;
	}
	return n;
}

// get the atoms in current job
extern "C" int get_total_atom_number() {
	int n = 0, nMol = 0;
	LIST<MMOLECULE> *plist = NULL;
	LIST<MMOLECULE> *pcglist = NULL;

	switch (JOB) {
	case GET_SINGLE_MM_MD_PROC:
		return get_macromol_atom_number(&mm);
		break;
	case SHOW_SINGLE_MM_MD_PROC:
		return get_macromol_atom_number(&mm);
		break;
	case GET_SINGLE_CGMM_MD_PROC:
		return get_macromol_atom_number(&cgmm);
		break;
	case SHOW_SINGLE_CGMM_MD_PROC:
		return get_macromol_atom_number(&cgmm);
		break;
	case GET_MULTI_MM_MD_PROC:
		return nAtoms(mm_md_cell);
		break;
	case SHOW_MULTI_MM_MD_PROC:
		return nAtoms(mm_md_cell);
		break;
	case GET_MULTI_CGMM_MD_PROC:
		for (nMol = 0; nMol < cgmm_md_cell.mm.n; nMol++) {
			n += get_macromol_atom_number(cgmm_md_cell.mm.m + nMol);
		}
		for (nMol = 0; nMol < cgmm_md_cell.sm.n; nMol++) {
			n += cgmm_md_cell.sm.m[nMol].nAtoms;
		}
		return n;
		break;
	case SHOW_MULTI_CGMM_MD_PROC:
		for (nMol = 0; nMol < cgmm_md_cell.mm.n; nMol++) {
			n += get_macromol_atom_number(cgmm_md_cell.mm.m + nMol);
		}
		for (nMol = 0; nMol < cgmm_md_cell.sm.n; nMol++) {
			n += cgmm_md_cell.sm.m[nMol].nAtoms;
		}
		return n;
		break;
	case SINGLE_MM_FREE_CELL_MD:
		return get_macromol_atom_number(&mm);
		break;
	case SINGLE_CGMM_FREE_CELL_MD:
		return get_macromol_atom_number(&cgmm);
		break;
	case MULTI_MM_FREE_CELL_MD:
		return nAtoms(mm_md_cell);
		break;
	case MULTI_MM_PERIDIC_CELL_MD:
		return nAtoms(mm_md_cell);
		break;
	case MULTI_CGMM_FREE_CELL_MD:
		for (nMol = 0; nMol < cgmm_md_cell.mm.n; nMol++) {
			n += get_macromol_atom_number(cgmm_md_cell.mm.m + nMol);
		}
		for (nMol = 0; nMol < cgmm_md_cell.sm.n; nMol++) {
			n += cgmm_md_cell.sm.m[nMol].nAtoms;
		}
		return n;
		break;
	case MULTI_CGMM_PERIDIC_CELL_MD:
		for (nMol = 0; nMol < cgmm_md_cell.mm.n; nMol++) {
			n += get_macromol_atom_number(cgmm_md_cell.mm.m + nMol);
		}
		for (nMol = 0; nMol < cgmm_md_cell.sm.n; nMol++) {
			n += cgmm_md_cell.sm.m[nMol].nAtoms;
		}
		return n;
		break;
	case MOLECULE_OPERATION:
		return 0;
	case GCMD:
		return nAtoms(vmdcell);
	case GET_GCMD_PROC:
		return nAtoms(vmdcell);
		break;
	case SHOW_GCMD_PROC:
		return nAtoms(vmdcell);
		break;
	default:
		return 0;
	}
}

extern bool ConstructCGMMStruct(CG_MMOLECULE *m, char *fname);
//extern bool ConstructCGSM(CGSM* cgsm, char *mol);
extern bool ReadCGMMStruct(ifstream &in, CG_MMOLECULE &mm);
extern bool ReadCGSMStruct(ifstream &in, CGSM &sm);

extern void molecule_check_periodic_cell(CGMM_MD_CELL &mmcell);


extern void test_SPME_EwaldSum();
extern void test_SPME_EwaldSum_2d();
extern void test_SPME_EwaldSum_2d_Exclude_InducedDipole();
extern void test_SPME_EwaldSum_ImageSlab();
extern void construct_profiles(); // some functions, exp2, erfc, cos, sin ...

extern void test1_SPME_EwaldSum_2d();
extern void test1_SPME_EwaldSum_2d_Exclude_InducedDipole();
extern void test1_SPME_EwaldSum_3d();
extern void test1_SPME_EwaldSum_3d_Exclude_InducedDipole();

extern "C" bool construct_macromolecule(char *fname) {
	if (strlen(::mlog.fname) == 0) ::mlog.init("mol_construct.log");

	construct_profiles(); // some functions, exp2, erfc, cos, sin ...

	//test_SPME_EwaldSum(); ::mlog.close();
	//test_SPME_EwaldSum_2d(); ::mlog.close();
	//test_SPME_EwaldSum_2d_Exclude_InducedDipole(); ::mlog.close();
	//test_SPME_EwaldSum_ImageSlab(); ::mlog.close();
	//test1_SPME_EwaldSum_2d(); ::mlog.close();
	//test1_SPME_EwaldSum_2d_Exclude_InducedDipole(); ::mlog.close();
	//test1_SPME_EwaldSum_3d(); ::mlog.close();
	//test1_SPME_EwaldSum_3d_Exclude_InducedDipole(); mlog.close();
	//return false;

	bool status = 0;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(errmsg, "Can not open file %s", fname);
		show_msg(errmsg);
		return false;
	}
	char item[100] = "\0", buffer[256] = "\0";
	char mol_template_fname[250] = "\0";
	int n = 0;

	readDielectricConst(fname, eps_dc);

	mm.reset();
	cgmm.reset();
	mm_md_cell.release();
	atompar_db.release();
	LJPARS.release();
	//torsion_db.release();
	ESolv_db.release(); // solvation energy
	cluster_db.release();
	_interface_constraint_::clusterIC_db.release();

	if (MD_MODE == SINGLE_MM_FREE_CELL_MD) {
		strcpy(item, "[MOLECULE_TEMPLATE]");
		if (search(in, item, buffer)) {
			in.getline(buffer, 250);
			if (sscanf(buffer, "%d %s", &n, mol_template_fname) == 2 && n == 1) {
				status = ConstructMoleculeStruct(&mm, mol_template_fname);
			}
			else status = ConstructSingleChainPolymerFromParamFile(mm, fname);
		}
		else {
			status = ConstructSingleChainPolymerFromParamFile(mm, fname);
		}
		if (!status) {mm.reset(); show_msg(errmsg);}

		mm.calMassCenter(true);

		// atompar_db was setup during constructing the molecules;
		// setup LJ parameters of all atoms
		init_LJ_Pars();
	}
	else if (MD_MODE == SINGLE_CGMM_FREE_CELL_MD) {
		strcpy(item, "[MOLECULE_TEMPLATE]");
		if (search(in, item, buffer)) {
			in.getline(buffer, 250);
			if (sscanf(buffer, "%d %s", &n, mol_template_fname) == 2 && n == 1) {
				status = ConstructCGMMStruct(&cgmm, mol_template_fname);
			}
			else status = ConstructSingleChainPolymerFromParamFile(cgmm, fname, false);
		}
		else {
			status = ConstructSingleChainPolymerFromParamFile(cgmm, fname, false);
		}
		if (!status) {cgmm.reset(); show_msg(errmsg);}

		// atompar_db was setup during constructing the molecules;
		// setup LJ parameters of all atoms
		init_LJ_Pars_CG();
	}
	else if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
		strcpy(item, "[MD_CELL_TEMPLATE]");
		if (!search(in, item, buffer)) {
			sprintf(errmsg, "Can not find %s in %s", item, fname);
			show_msg(errmsg); in.close(); 
			atompar_db.release(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", &mol_template_fname) != 1 || strlen(mol_template_fname) < 1) {
			sprintf(errmsg, "number of molecules can not be find from %s", buffer);
			show_msg(errmsg); in.close(); 
			atompar_db.release(); return false;
		}

		mm_md_cell.SetMol(0, 0); // reset mm_md_cell

		if (!(status = ConstructDefinedMDCell<MMOLECULE, SMOLECULE, PMOLECULE, MMOL_MD_CELL>(mol_template_fname, mm_md_cell, (void*)(&ConstructMoleculeStruct),
			(void*)(&ConstructSimpleMolecule), (void*)(&ConstructPointMolecule), (void*)(&ReadMoleculeStruct), (void*)(&ReadSolidMoleculeStruct), (void*)(&ReadPointMoleculeStruct)))) {
			mm_md_cell.SetMol(0, 0, 0);
			show_msg(errmsg); in.close(); 
			atompar_db.release(); return false;
		}
		mm_md_cell.setMolIndx(); // label the molecule index for all clusters

		for (n = 0; n < mm_md_cell.mm.n; n++) mm_md_cell.mm.m[n].calMassCenter(true);
		for (n = 0; n < mm_md_cell.sm.n; n++) mm_md_cell.sm.m[n].calInertiaMoment();

		if (MD_MODE == MULTI_MM_PERIDIC_CELL_MD) molecule_check_periodic_cell(mm_md_cell);

		// atompar_db was setup during constructing the molecules;
		// setup LJ parameters of all atoms
		init_LJ_Pars();
	}
	else if (MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
		strcpy(item, "[MD_CELL_TEMPLATE]");
		if (!search(in, item, buffer)) {
			sprintf(errmsg, "Can not find %s in %s", item, fname);
			show_msg(errmsg); in.close(); 
			cgatompar_db.release(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%s", &mol_template_fname) != 1 || strlen(mol_template_fname) < 1) {
			sprintf(errmsg, "number of molecules can not be find from %s", buffer);
			show_msg(errmsg); in.close(); 
			cgatompar_db.release(); return false;
		}

		cgmm_md_cell.SetMol(0, 0); // reset mm_md_cell

		if (!(status = ConstructDefinedMDCell<CG_MMOLECULE, CGSM, CGSM, CGMM_MD_CELL>(mol_template_fname, cgmm_md_cell, (void*)(&ConstructCGMMStruct),
			(void*)(&ConstructCGSM), (void*)(&ConstructCGSM), (void*)(&ReadCGMMStruct), (void*)(&ReadCGSMStruct), (void*)(&ReadCGSMStruct)))) {
			cgmm_md_cell.SetMol(0, 0);
			show_msg(errmsg); in.close(); 
			cgatompar_db.release(); return false;
		}
		cgmm_md_cell.setMolIndx(); // label the molecule index for all clusters

		if (MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) _coarse_grain_::molecule_check_periodic_cell(cgmm_md_cell);

		// atompar_db was setup during constructing the molecules;
		// setup LJ parameters of all atoms
		init_LJ_Pars_CG();
	}

	in.close();
	show_log("molecules are constructed!", true);

	return status;
}

extern EF_DPROF LJ12_6;
PAIRWISE_DB< PAIRWISE<EF_DPROF> > sr_db; // pair-wise short-range interaction database
extern char sr_fname[256];

PAIRWISE_DB< PAIRWISE<float> > TholeRadius_db; // pair-wise database about Thole-radius between ions in Thole polarization model, s = a * (alpha_i * alpha_j)^1/6

void MD_SIM_SINGLE_MOL(char *mdpar_fname) {
#ifdef __GPU__
	cudaDeviceReset(); 
#endif

	if (mm.nCluster <= 0) {show_msg("molecule is not constructed");return;}

	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {show_msg(errmsg); return;}

	LFVERLET_MM_MD mdpar;
	if (!read_MD_BASE_Pars(mdpar_fname, (MD_BASE*)(&mdpar))) {show_msg(errmsg); return;}
	::timestep_md = mdpar.dt;

	float tao = 0, max_ksi = 0;
	HOOVER_DYN ihDyn, ft_hDyn, fr_hDyn; 
	ihDyn.nhc = new NHC<MNHC>; ft_hDyn.nhc = new NHC<MNHC>; fr_hDyn.nhc = new NHC<MNHC>;
	if (!read_HOOVER_Pars(mdpar_fname, tao, max_ksi)) {
		ihDyn.release(true); ft_hDyn.release(true); fr_hDyn.release(true);
		show_msg(errmsg); return;
	}
	::max_ksi_Hoover = 0.5f / mdpar.dt;
	if (max_ksi > ::max_ksi_Hoover) max_ksi = ::max_ksi_Hoover;
	ihDyn.set_vars(tao, mdpar.dt, max_ksi);
	ft_hDyn.set_vars(tao, mdpar.dt, max_ksi);
	fr_hDyn.set_vars(tao, mdpar.dt, max_ksi);
	ihDyn.nhc->set_freedom(mm.nCluster - 1); ft_hDyn.nhc->set_freedom(3); fr_hDyn.nhc->set_freedom(3);
	mdpar.set_HDyn(&ft_hDyn, &fr_hDyn, &ihDyn);

	mdpar.mmsave.mdsave = &msave;
	mdpar.LOOPS *= msave.nsave;

	CMM_CELL3D<BASIC_CLUSTER> cmm;
	int interact_cal_method = 0;
	int n = 0, nc = 0, nCheckCell = 5;
	int nx, ny, nz;
	//float size_cluster = 0;//, rLJ_cut = 0, rshell = 0;
	float dx, dy, dz;
	char buffer[256] = "\0", md_cell_fname[256] = "\0";
	char* pch = NULL;

	if (!read_msg(mdpar_fname, buffer, "[CLUSTER_SIZE]", true) || sscanf(buffer, "%f", &size_cluster) != 1) size_cluster = 2;

	if (!read_msg(mdpar_fname, buffer, "[INTERACTION_CALCULATION_METHOD]", true)) interact_cal_method = _PAIRWISE;
	else if (sscanf(buffer, "%d", &interact_cal_method) != 1) interact_cal_method = _PAIRWISE;

	//if (interact_cal_method == _CMM) {
		pch = buffer;
		NextSection(&pch);
		if (sscanf(pch, "%d", &nCheckCell) != 1) {
			sprintf(errmsg, "can not get the loop number to re-check CMM during MD from %s", buffer);
			show_msg(errmsg);
			md_run = false;
			ihDyn.release(true); ft_hDyn.release(true); fr_hDyn.release(true);
			return;
		}

		//rLJ_cut = rLJ_cut_max();
		//rshell = (size_cluster > rLJ_cut ? size_cluster : rLJ_cut);

		//cmm.nshell = int(rshell / (hw + hw) + 0.5f);
		//if (cmm.nshell * (hw + hw) < rshell) cmm.nshell++;

		::min_cmm_size = size_cluster * 2;

		cmm.nCheckCluster = nCheckCell;

		if (!read_msg(mdpar_fname, buffer, "[MD_CELL_TEMPLATE]", true) || sscanf(buffer, "%s", md_cell_fname) != 1) {
			ihDyn.release(true); ft_hDyn.release(true); fr_hDyn.release(true);
			sprintf(errmsg, "can not get [MD_CELL_TEMPLATE] for cell size in %s", mdpar_fname); return;
		}
		if (!ReadDefinedCellSize(md_cell_fname, dx, dy, dz)) {
			ihDyn.release(true); ft_hDyn.release(true); fr_hDyn.release(true); return;
		}
		dx /= 2; dy /= 2; dz /= 2;
		nx = int(2 * dx / min_cmm_size); nx = (nx < 1 ? 1 : nx);
		ny = int(2 * dy / min_cmm_size); ny = (ny < 1 ? 1 : ny);
		nz = int(2 * dz / min_cmm_size); nz = (nz < 1 ? 1 : nz);
		cmm.set_cell(nx, ny, nz);
		cmm.set_cell_range(false, -dx, dx, -dy, dy, -dz, dz);
		cmm.set_cell_pos();
	//}
#if _MPI_ == 1
	sprintf(buffer, "%s%d.log", mdpar.mmsave.mdsave->title, ::mpi_id);
#else
	sprintf(buffer, "%s.log", mdpar.mmsave.mdsave->title);
#endif
	if (!mlog.init(buffer)) {
		sprintf(errmsg, "can not open file %s for log\nshall we continue with no log record ?", buffer);
#if _SYS_ == _WINDOWS_SYS_
		if (AfxMessageBox(LPCTSTR(CString(errmsg)), MB_YESNO) == IDNO) {
			ihDyn.release(true); ft_hDyn.release(true); fr_hDyn.release(true);
			md_run = false; return;
		}
#else
		ihDyn.release(true); ft_hDyn.release(true); fr_hDyn.release(true);
		show_msg(errmsg); md_run = false; return;
#endif
	}

	if (!read_msg(mdpar_fname, buffer, "[RANDOM_FORCE]", true) ||
		sscanf(buffer, "%d %f", &n, &randf) != 2 || n <= 0) bRandomForce = false;
	else {bRandomForce = true; randf *= kT / U_MassAccel / mdpar.dt; fwhm_randf = randf / 3; } // randf has unit kT/Angs.

	readExternalElectricField(mdpar_fname, bEext, Eext); 
	readForceScale(mdpar_fname, scale_force, force_max, torque_max);
	if (scale_force) {
		sprintf(buffer, "scale force / torque when they are more than %7.4f / %7.4f", force_max, torque_max);
		show_infor(buffer);
		force2_max = force_max * force_max;
		torque2_max = torque_max * torque_max;
	}
	readLocalRelax(mdpar_fname, local_relax, dta_max_relax);
	readLocalInteract(mdpar_fname, local_interact);
	readInternalKinEnergy(mdpar_fname, Ek0Internal);
	readBrownianPars(mdpar_fname, ::bBrownian, ::ita_Brownian, ::cd_drag_trans, ::cd_drag_rot);
	if (::bBrownian) {
		//::ita_Brownian *= kT / U_MassAccel / mdpar.dt;
		//::cd_drag_trans *= kT / U_MassAccel / mdpar.dt;
		//::cd_drag_rot *= kT / U_MassAccel / mdpar.dt;
		::ita_Brownian /= mdpar.dt;
		::cd_drag_rot /= mdpar.dt;
		::cd_drag_trans /= mdpar.dt;
	}

	readImplicitSolvent(mdpar_fname, bImplicitSolvent);
	if (bImplicitSolvent) init_ImplicitSolvation_Pars();

	// shall we free the macro-molecule base ?
	readMacroMolFreeBase(mdpar_fname, mm.free_base);

	// setup torsion of the molecule
	dh_db.release();
	//torsion_db.release();
	if (!local_relax) {
		//mm.SetupHingeTorsion();
		initTorsionPars(&mm);
		dh_db.convert_array();
	}
	for (nc = 0; nc < mm.nCluster; nc++) mm.cluster[nc].check_eneutral();
#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	for (nc = 0; nc < mm.nCluster; nc++) {
		mm.cluster[nc].calElectroMultipole(false);
	}
#endif

	/*
	bool coarse_grain = false;
	char cgpar_fname[256] = "\0";
	int nNormalMD = 1, nCG = 1;
	CG_MD_VARS cgvar;
	if (!readCoarseGrainMD(mdpar_fname, coarse_grain, cgpar_fname, nNormalMD, nCG)) return;
	if (coarse_grain) {
		if (!read_coarse_grain_md_var(cgpar_fname, cgvar.dt, cgvar.tao, cgvar.max_ksi, cgvar.nCMMRefresh)) return;
	}
	*/

	//if (interact_cal_method == _CMM) {
		::interact_cal_method = _CMM;
		cmm_init_distribute_cluster(mm, cmm, true); // set base cluster at the center cell
		cmm_init_subcell<BASIC_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
	//}

	COMMAND_MD_STOP = false;
	md_run = true;
	construct_LJ12_6_prof(::iFormat_LJ);
	sprintf(errmsg, "normalized cutoff distance for 12-6 L-J interaction is set to %f", ::rcut_LJ_norm); show_log(errmsg, true);
	float rcut_SR = 0, rmin_SR = 0;
	if (strlen(sr_fname) > 0 && initPairwiseShortRangeInteract(sr_fname, ::atompar_db, ::sr_db, rcut_SR, rmin_SR)) {
		if (rcut_SR > ::rcut_LJ) {
			::rcut_LJ = rcut_SR; 
			::r2cut_LJ = ::rcut_LJ * ::rcut_LJ;
		}
	}
	sprintf(errmsg, "maximum short-range-interaction cut-off distance is %f Angs.", ::rcut_LJ); show_log(errmsg, true);
	if (strlen(sr_fname) > 0) update_atomic_pars(sr_fname, ::atompar_db); // disable Lennard-Jones ? and update atomic polarizability

	sprintf(::errmsg, "contraint the translating and rotating kinetic energies of mass center to the thermal energy");
	show_log(::errmsg, true);
#if _MPI_ == 1
	MPI_SingleMM_JobAssign(mm.nCluster);
	mdpar.resetClusters(mm.nCluster);
	if (interact_cal_method == _CMM) {
		MPI_POLYMER(mm, mdpar, &cmm);
	}
	else {
		sprintf(::errmsg, "MPI has to use CMM method to calculate interaction, method choosed : %d", interact_cal_method); show_msg(::errmsg);
	}
#else
	int nNormalMD = 1000;
	MD(mm, mdpar, &cmm, interact_cal_method, nNormalMD);
#endif

	md_run = false;
	mlog.close();
	cmm.set_cell(0, 0, 0);
	cgmm.reset();
	// release database of atom par, LJ pars, torsion
	atompar_db.release();
	LJPARS.release();
	ESolv_db.release();
	//ESolvPars.release(); // solvation energy
	//torsion_db.release(); dh_db.release();
	ihDyn.release(true); ft_hDyn.release(true); fr_hDyn.release(true);
	cluster_db.release();
	_interface_constraint_::clusterIC_db.release();

#ifdef __GPU__
	cudaDeviceReset(); 
#endif

	return;
}

namespace _coarse_grain_ {
	CGBOUND_DB cgbound_db;
	CGBOUND_DB cg_non_neighbor_db;
	CGBOUND_DB cg_intramol_non_neighbor_db;
	CGBOUND_DIPOLE_DISTRBT_DB cg_dipole_distrbt_db;

	extern void MD_POLYMER(CG_MMOLECULE &cgmm, LFVERLET_CGMM_MD& mdpar, CMM_CELL3D<CG_CLUSTER> *cmm);
	extern bool get_neighbor_level(char *fname, int &n_neighbors);
}

void MD_SIM_SINGLE_CGMM(char *mdpar_fname) {
#ifdef __GPU__
	cudaDeviceReset(); 
#endif

	if (cgmm.nCluster <= 0) {show_msg("molecule is not constructed");return;}

	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {show_msg(errmsg); return;}

	// Leapfrog Verlet algorithm
	LFVERLET_CGMM_MD mdpar;
	if (!read_MD_BASE_Pars(mdpar_fname, (MD_BASE*)(&mdpar))) {show_msg(errmsg); return;}
	::timestep_md = mdpar.dt;

	float tao = 0, max_ksi = 0;
	HOOVER_DYN hDyn; hDyn.nhc = new NHC<MNHC>;
	if (!read_HOOVER_Pars(mdpar_fname, tao, max_ksi)) {show_msg(errmsg); hDyn.release(true); return;}
	::max_ksi_Hoover = 0.5f / mdpar.dt;
	if (max_ksi > ::max_ksi_Hoover) max_ksi = ::max_ksi_Hoover;
	hDyn.set_vars(tao, mdpar.dt, max_ksi);
	mdpar.set_HDyn(&hDyn, NULL);

	mdpar.mmsave.mdsave = &msave;
	mdpar.LOOPS *= msave.nsave;

	CMM_CELL3D<CG_CLUSTER> cmm;
	int interact_cal_method = 0;
	int n = 0, nc = 0, nCheckCell = 5;
	int nx, ny, nz;
	//float size_cluster = 0;//, rLJ_cut = 0, rshell = 0;
	float dx, dy, dz;
	char buffer[256] = "\0", md_cell_fname[256] = "\0";
	char* pch = NULL;

	if (!read_msg(mdpar_fname, buffer, "[CLUSTER_SIZE]", true) || sscanf(buffer, "%f", &size_cluster) != 1) size_cluster = 2;
	
	if (!read_msg(mdpar_fname, buffer, "[INTERACTION_CALCULATION_METHOD]", true)) interact_cal_method = _CMM;
	else if (sscanf(buffer, "%d", &interact_cal_method) != 1) {
		sprintf(errmsg, "INTERACTION METHOD has to be CMM"); show_msg(errmsg); md_run = false; hDyn.release(true); return;
	}
	::interact_cal_method = interact_cal_method;

	if (interact_cal_method == _CMM) {
		pch = buffer;
		NextSection(&pch);
		if (sscanf(pch, "%d", &nCheckCell) != 1) {
			sprintf(errmsg, "can not get the loop number to re-check CMM during MD from %s", buffer);
			show_msg(errmsg);
			md_run = false;
			hDyn.release(true); 
			return;
		}

		::min_cmm_size = size_cluster * 2;

		cmm.nCheckCluster = nCheckCell;

		if (!read_msg(mdpar_fname, buffer, "[MD_CELL_TEMPLATE]", true) || sscanf(buffer, "%s", md_cell_fname) != 1) {
			sprintf(errmsg, "can not get [MD_CELL_TEMPLATE] for cell size in %s", mdpar_fname); hDyn.release(true); return;
		}
		if (!ReadDefinedCellSize(md_cell_fname, dx, dy, dz)) {hDyn.release(true); return;}
		dx /= 2; dy /= 2; dz /= 2;
		nx = int(2 * dx / min_cmm_size); nx = (nx < 1 ? 1 : nx);
		ny = int(2 * dy / min_cmm_size); ny = (ny < 1 ? 1 : ny);
		nz = int(2 * dz / min_cmm_size); nz = (nz < 1 ? 1 : nz);
		cmm.set_cell(nx, ny, nz);
		cmm.set_cell_range(false, -dx, dx, -dy, dy, -dz, dz);
		cmm.set_cell_pos();
	}
	else {
		sprintf(errmsg, "interaction calculation has to use CMM"); show_msg(errmsg);
		md_run = false; hDyn.release(true); return;
	}
#if _MPI_ == 1
	sprintf(buffer, "%s%d.log", mdpar.mmsave.mdsave->title, ::mpi_id);
#else
	sprintf(buffer, "%s.log", mdpar.mmsave.mdsave->title);
#endif
	if (!mlog.init(buffer)) {
		sprintf(errmsg, "can not open file %s for log\nshall we continue with no log record ?", buffer);
#if _SYS_ == _WINDOWS_SYS_
		if (AfxMessageBox(LPCTSTR(CString(errmsg)), MB_YESNO) == IDNO) {
			md_run = false; hDyn.release(true); return;
		}
#else
		show_msg(errmsg); md_run = false; return;
#endif
	}

	if (!read_msg(mdpar_fname, buffer, "[RANDOM_FORCE]", true) ||
		sscanf(buffer, "%d %f", &n, &randf) != 2 || n <= 0) bRandomForce = false;
	else {bRandomForce = true; randf *= kT / U_MassAccel / mdpar.dt; fwhm_randf = randf / 3;} // randf has unit kT/Angs.

	readExternalElectricField(mdpar_fname, bEext, Eext); 
	if (!bEext || Eext == 0) ::bCGdipole = false;
	readForceScale(mdpar_fname, scale_force, force_max, torque_max);
	if (scale_force) {
		sprintf(buffer, "scale force / torque when they are more than %7.4f / %7.4f", force_max, torque_max);
		show_infor(buffer);
		force2_max = force_max * force_max;
		torque2_max = torque_max * torque_max;
	}
	//readLocalRelax(mdpar_fname, local_relax, dta_max_relax);
	//readLocalInteract(mdpar_fname, local_interact);
	readInternalKinEnergy(mdpar_fname, Ek0Internal);
	readBrownianPars(mdpar_fname, ::bBrownian, ::ita_Brownian, ::cd_drag_trans, ::cd_drag_rot);
	if (::bBrownian) {
		//::ita_Brownian *= kT / U_MassAccel / mdpar.dt;
		//::cd_drag_trans *= kT / U_MassAccel / mdpar.dt;
		//::cd_drag_rot *= kT / U_MassAccel / mdpar.dt;
		::ita_Brownian /= mdpar.dt;
		::cd_drag_rot /= mdpar.dt;
		::cd_drag_trans /= mdpar.dt;
	}

	readImplicitSolvent(mdpar_fname, bImplicitSolvent);
	if (bImplicitSolvent) init_ImplicitSolvation_Pars();

	// setup torsion of the molecule
	//if (!local_relax) cgmm.SetupHingeTorsion();
	//for (nc = 0; nc < cgmm.nCluster; nc++) {
	//	cgmm.cluster[nc].calElectroMultipole(false);
	//}

	// bound database -- neighbors were setup when cgmm is constructed
	// we have to setup neighbors of the cgcluster first
	cgbound_db.releaseChain(true); cgbound_db.releaseArray(true);
#if _CG_NEIGHBOR_INTERACT_ == 0
	// if we use FENE bound interaction, each cluster has nearest neighbor only
	cgmm.setup_cluster_neighbors(1); // 
	attach_bound_db(cgmm, cgbound_db, 0); // 0 -- FENE, 1 -- FROM FILE
	// using given potential profile, unit kT
#elif _CG_NEIGHBOR_INTERACT_ == 1
	int n_neighbors = 1;
	if (!get_neighbor_level(cgatom_par_fname, n_neighbors)) {md_run = false; hDyn.release(true); return;}
	cgmm.setup_cluster_neighbors(n_neighbors); // 
	attach_bound_db(cgmm, cgbound_db, 1); // 0 -- FENE, 1 -- FROM FILE
	#if _INTRA_INTER_INTERACTS == 1
	attach_non_neighbor_db(::cgatompar_db, cg_intramol_non_neighbor_db);
	#else
	attach_non_neighbor_db(::cgatompar_db, cg_non_neighbor_db);
	#endif
#endif

	cgbound_db.Chain2Array();
	cg_intramol_non_neighbor_db.Chain2Array();
	cg_non_neighbor_db.Chain2Array();

	// basic cell of CMM
	cmm_init_distribute_cluster(cgmm, cmm, true); // set base cluster at the center cell
	cmm_init_subcell<CG_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size

	COMMAND_MD_STOP = false;
	md_run = true;
	construct_LJ12_6_prof(::iFormat_LJ);
	sprintf(errmsg, "normalized cutoff distance for 12-6 L-J interaction is set to %f", ::rcut_LJ_norm); show_log(errmsg, true);

	sprintf(::errmsg, "contraint the translating and rotating kinetic energies of mass center to the thermal energy");
	show_log(::errmsg, true);
//#if _MPI_ == 1
	//MPI_SingleMM_JobAssign(mm.nCluster);
	//mdpar.resetClusters(mm.nCluster);
	//MPI_POLYMER(mm, mdpar, &cmm);
//#else
	MD_POLYMER(cgmm, mdpar, &cmm);
//#endif

	md_run = false;
	mlog.close();
	cmm.set_cell(0, 0, 0);
	mm.reset();
	// release database of atom par, LJ pars, torsion
	atompar_db.release();
	LJPARS.release();
	ESolv_db.release();
	//ESolvPars.release(); // solvation energy
	//torsion_db.release();
	cgbound_db.releaseArray(true); cgbound_db.releaseChain(true);
	cg_dipole_distrbt_db.releaseArray(true); cg_dipole_distrbt_db.releaseChain(true);
	hDyn.release(true); 
	cluster_db.release();
	_interface_constraint_::clusterIC_db.release();

#ifdef __GPU__
	cudaDeviceReset(); 
#endif

	return;
}

namespace _Thole_polarize_ {
	extern void initTholeRadius(ATOM_COMM_PARS_DATABASE &atompar_db, PAIRWISE_DB< PAIRWISE<float> > &TholeRadius_db);
}

extern void GPU_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, int md_mode);
extern void SLAB_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm);
extern void ImageSlab_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm);
char parfname_IConstraint[256];

extern "C" bool cluster_spatial_distribution(char* par_fname);

void MD_SIM_MULTI_MOL_CELL(char *mdpar_fname) {
#ifdef __GPU__
	cudaDeviceReset(); 
#endif

	// MD ESSEMBLE
/*
	char pf[256] = "distribt_Br.ini";
	cluster_spatial_distribution(pf);
	return;
*/
	show_log("initialize md ... ", true);
	char buffer[256] = "\0";
	if (read_msg(mdpar_fname, buffer, "[MD_MODE]")) {
		if (sscanf(buffer, "%d", &MD_ENSEMBLE) != 1) MD_ENSEMBLE = MD_NVT;
		else if (MD_ENSEMBLE != MD_NVT && MD_ENSEMBLE != MD_NPT) MD_ENSEMBLE = MD_NVT;
	}

	// constraint
	int mconstraint = 0;
	char constraint_parfname[256] = "\0";
	if (read_msg(mdpar_fname, buffer, "[MD_CONSTRAINT]", true)) {
		if (sscanf(buffer, "%d", &mconstraint) == 1) {
			if (mconstraint > 0) {
				if (mconstraint == MD_PMF_CONSTRAINT) {
					if (sscanf(buffer, "%d %s", &mconstraint, constraint_parfname) != 2) {
						sprintf(errmsg, "constraint parameter file is not given, see : %s", buffer);
						show_log(errmsg, true); return;
					}
					else strcpy(parfname_IConstraint, constraint_parfname);
				}
				else if (mconstraint == MD_SLAB || mconstraint == MD_IMAGE_SLAB) {
					if (sscanf(buffer, "%d %s", &mconstraint, constraint_parfname) != 2) {
						sprintf(errmsg, "constraint parameter file is not given, see : %s", buffer);
						show_log(errmsg, true); return;
					}
					else strcpy(parfname_IConstraint, constraint_parfname);
				}
				else if (mconstraint == MD_NO_CONSTRAINT) {
					if (sscanf(buffer, "%d %s", &mconstraint, constraint_parfname) == 2)
						strcpy(parfname_IConstraint, constraint_parfname);
				}
				else {
					sprintf(errmsg, "ERROR: Unknown type of constraint. see : %s", buffer); 
					show_log(errmsg, true); return;
				}
			}
			else {
				show_log("no constraint on MD", true); mconstraint = 0;
			}
		}
		else {
			sprintf(errmsg, "no constraint model is given, see : %s", buffer); show_log(errmsg, true);
			show_log("disable constraint", true); mconstraint = 0;
		}
	}
	else mconstraint = 0;

	// check local relax or local interact here, since they can be used for mass center constraint
	readLocalRelax(mdpar_fname, local_relax, dta_max_relax);
	if (::local_relax) {
		sprintf(buffer, "calculate the local relationship within neighbor CMM"); show_infor(buffer);
	}
	readLocalInteract(mdpar_fname, local_interact);

	// MD SAVE
	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {show_msg(errmsg); mlog.close(); return;}

	// Leapfrog Verlet algorithm
	//LFVERLET_MM_MD mdpar;
	//if (!read_LFVerlet_MD_Pars(mdpar_fname, mdpar)) {show_msg(errmsg); mlog.close(); return;}

	int i;

	// Nose-Hoover Chains
	float dt = 0, tao = 0, max_ksi = 0;
	if (!read_MD_TIMESTEP(mdpar_fname, dt)) {show_msg(errmsg); mlog.close(); return;}
	::timestep_md = dt;
	if (!read_HOOVER_Pars(mdpar_fname, tao, max_ksi)) {show_msg(errmsg); mlog.close(); return;}
	//::max_ksi_Hoover = 0.5f / dt;
	::max_ksi_Hoover = 0.8f / dt;
	if (max_ksi > ::max_ksi_Hoover) max_ksi = ::max_ksi_Hoover;
	mm_md_cell.set_md_timestep(dt);

	if (!ReadMDCellNHC<MMOLECULE, SMOLECULE, PMOLECULE, MMOL_MD_CELL>("NHC.INI", mm_md_cell, dt)) {
		mlog.close(); return;
	}
	mm_md_cell.hdyn.m[mm_md_cell.hdyn.n - 1].set_vars(tao, dt, max_ksi); // #the last one -- default from MD.INI
	mm_md_cell.init_nhc_mols();
	sprintf(errmsg, "%d Nose-Hoover chain are used!", mm_md_cell.hdyn.n); show_log(errmsg, true); 
	// end of Nose-Hoover Chain

	mm_md_cell.TotalMass();

	// Nose-Hoover Chain for the constraint on mass-center velocity of the whole cell
	if (mconstraint == MD_SLAB || mconstraint == MD_IMAGE_SLAB) {
		if (::local_interact) {
			mm_md_cell.bHDynMC_T = false; // z axis
		}
		else {
			mm_md_cell.bHDynMC_T = false; //true; // z axis
		}
		{
			int v = 0, nv = 0;
			float ita = 0;
			if (read_item(mdpar_fname, "[SLAB_Z]", buffer) && (nv = sscanf(buffer, "%d %f", &v, &ita) >= 1)) {
				if (nv != 2) ita = 0.1;
			}
			else v = 0;
			if (v > 0) {
				sprintf(errmsg, "each MD time step decrease speed %f %%", ita * 100); show_log(errmsg, true);
				::ita_mc = ita / ::timestep_md;
			}
			if (v == 1) {// thermalstat on translation along z
				mm_md_cell.bHDynMC_T = false;
				mm_md_cell.ix_mct = 0; mm_md_cell.iy_mct = 0; mm_md_cell.iz_mct = 1; // check movement along x, y, z axis
				show_log("cell mass center : disable translation thermostat along z axis.", true);
			}
			else if (v == 2) {
				mm_md_cell.bHDynMC_T = false;
				mm_md_cell.ix_mct = 0; mm_md_cell.iy_mct = 0; mm_md_cell.iz_mct = 1; // check movement along x, y, z axis
				show_log("cell mass center : disable translation thermostat along x, y and z axis.", true);
			}
			else {
				mm_md_cell.bHDynMC_T = false;
				show_log("cell mass center : disable thermostat along x, y and z axis.", true);
			}
		}
		mm_md_cell.check_freedom();
		mm_md_cell.setup_mcHDyn(tao, dt, max_ksi);

		show_log("translation movement of the whole cell is constrainted!", true);
	}
	else {
		if (::local_interact) {
			mm_md_cell.bHDynMC_T = false; // z axis
		}
		else {
			// translation constraint is enabled for 3d cell, when interface constraint is applied
			// which means the whole system will be confined in a local area along z axis
			if (strlen(::parfname_IConstraint) > 0) mm_md_cell.bHDynMC_T = false;
			else mm_md_cell.bHDynMC_T = false;
		}
		mm_md_cell.check_freedom();
		mm_md_cell.setup_mcHDyn(tao, dt, max_ksi);
		mm_md_cell.ix_mct = 0; mm_md_cell.iy_mct = 0; mm_md_cell.iz_mct = 1; // check movement along x, y, z axis
		show_log("translation movement of the whole cell is constraint!", true);
		{
			int v = 0, nv = 0;
			float ita = 0;
			if (read_item(mdpar_fname, "[SLAB_Z]", buffer) && (nv = sscanf(buffer, "%d %f", &v, &ita)) >= 1) {
				if (nv == 1) ita = 0.1;
			}
			else v = 0;
			if (v > 0) {
				sprintf(errmsg, "each MD time step decrease speed %f %%", ita * 100); show_log(errmsg, true);
				::ita_mc = ita / ::timestep_md;
			}
			if (v == 1) {// thermalstat on translation along z
				mm_md_cell.bHDynMC_T = false;
				mm_md_cell.ix_mct = 0; mm_md_cell.iy_mct = 0; mm_md_cell.iz_mct = 1; // check movement along x, y, z axis
				show_log("cell mass center : translation thermostat along z axis.", true);
			}
			else if (v == 2) {
				mm_md_cell.bHDynMC_T = false;
				mm_md_cell.ix_mct = 0; mm_md_cell.iy_mct = 0; mm_md_cell.iz_mct = 1; // check movement along x, y, z axis
				show_log("cell mass center : translation thermostat along x, y and z axis.", true);
			}
			else {
				mm_md_cell.bHDynMC_T = false;
				show_log("cell mass center : disable thermostat along x, y and z axis.", true);
			}
		}
	}
	mm_md_cell.Ethermal_mct = 0;
	if (mm_md_cell.ix_mct) mm_md_cell.Ethermal_mct += 0.5;
	if (mm_md_cell.iy_mct) mm_md_cell.Ethermal_mct += 0.5;
	if (mm_md_cell.iz_mct) mm_md_cell.Ethermal_mct += 0.5;
	mm_md_cell.Mthermal_mc = mm_md_cell.Ethermal_mct;

	{
			int n = 0, nv = 0;
			float v0 = 0, v1 = 0;
			if (read_item(mdpar_fname, "[SLAB_Z_SPRING]", buffer) && (nv = sscanf(buffer, "%d %f %f", &n, &v0, &v1)) >= 2) {
				if (n <= 0) {::bZSpring = false; ::k_zspring[0] = 0; ::k_zspring[1] = 0;}
				else {
					::bZSpring = true;
					::k_zspring[0] = v0 * kT / U_MassAccel;
					if (nv == 3) ::k_zspring[1] = v1 * kT / U_MassAccel;
					else ::k_zspring[1] = 0;
				}
			}
			else ::bZSpring = false;

			if (!::bZSpring) show_log("NO spring is applied to the displacement of mass center!", true);
			else {
				show_log("A spring is applied to the displacement of mass center along z-axis!", true);
				sprintf(::errmsg, "    with spring-constant: k[0] = %f, k[1] = %f kT/Angs.", v0, v1);
				show_log(::errmsg, true); strcpy(::errmsg, "\0");
			}
	}

	long nloops = 0;
	if (!read_MD_LOOPS(mdpar_fname, nloops)) {show_msg(errmsg); mlog.close(); return;}
	nloops = nloops * (long)msave.nsave;
	mm_md_cell.setMD(nloops);
	
	mm_md_cell.SetSave(&msave);

	CMM_CELL3D<BASIC_CLUSTER> cmm;
	int interact_cal_method = 0;
	int n = 0, nc = 0, nCheckCell = 5;
	int nx, ny, nz;
	//float size_cluster = 0;//, rLJ_cut = 0, rshell = 0;;
	char* pch = NULL;

	if (!read_msg(mdpar_fname, buffer, "[CLUSTER_SIZE]", true) || sscanf(buffer, "%f", &size_cluster) != 1) size_cluster = 2;

	if (!read_msg(mdpar_fname, buffer, "[INTERACTION_CALCULATION_METHOD]", true) || 
		sscanf(buffer, "%d", &interact_cal_method) != 1 || interact_cal_method != 1) {
		sprintf(errmsg, "Interaction Calculation method for multi-molecule in cell has to be CMM"); show_msg(errmsg); 
		sprintf(errmsg, "CMM is used"); show_msg(errmsg); interact_cal_method = 1;
		//mlog.close(); md_run = false; return;
	}

#if _MPI_ == 1
	if (::MD_MODE != MULTI_MM_PERIDIC_CELL_MD) {
		sprintf(::errmsg, "MPI works on periodical cell only"); show_msg(::errmsg);
		mlog.close(); md_run = false; return;
	}
#endif

	if (interact_cal_method == _CMM) {
		pch = buffer;
		NextSection(&pch);
		if (sscanf(pch, "%d", &nCheckCell) != 1) {
			sprintf(errmsg, "can not get the loop number to re-check CMM during MD from %s", buffer);
			show_msg(errmsg); md_run = false; mlog.close(); return;
		}

		::min_cmm_size = size_cluster * 2;

		nx = int(2 * mm_md_cell.h[0] / min_cmm_size); nx = (nx < 1 ? 1 : nx);
		ny = int(2 * mm_md_cell.h[1] / min_cmm_size); ny = (ny < 1 ? 1 : ny);
		nz = int(2 * mm_md_cell.h[2] / min_cmm_size); nz = (nz < 1 ? 1 : nz);
		cmm.set_cell(nx, ny, nz);
		if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
			cmm.set_cell_range(true, -mm_md_cell.h[0], mm_md_cell.h[0], -mm_md_cell.h[1], mm_md_cell.h[1], -mm_md_cell.h[2], mm_md_cell.h[2]);
			cmm.set_cell_pos();
		}
		else if (::MD_MODE == MULTI_MM_FREE_CELL_MD) {
			cmm.set_cell_range(false, -mm_md_cell.h[0], mm_md_cell.h[0], -mm_md_cell.h[1], mm_md_cell.h[1], -mm_md_cell.h[2], mm_md_cell.h[2]);
			cmm.set_cell_pos();
		}

		cmm.nCheckCluster = nCheckCell;
	}

	readSurfaceBoundary(mdpar_fname);
	readNetForceCorrectionMode(mdpar_fname);
	if (::bZSpring && ::mode_NetForce == 0) {
		show_log("Enable net-force compensation along z-axis, because z-spring is enabled!", true); ::mode_NetForce = 2;
	}

	ISOTROPIC_NPT_CELL npt_cell;
	float W_npt = 0, tao_npt = 0;
	if (MD_ENSEMBLE == MD_NPT) {
		npt_cell.mdcell = &mm_md_cell; npt_cell.cmm = &cmm;
		npt_cell.NF_pressure = check_NF_pressure(npt_cell.mdcell);
		if (!readNPTPars(mdpar_fname, W_npt, tao_npt)) {md_run = false; mlog.close(); return;}
		npt_cell.pdyn.W = W_npt * W_npt * (npt_cell.NF_pressure + 3);
		npt_cell.setup_pressure_NHC();
		npt_cell.phdyn.nhc->N = 9; // d^2
		npt_cell.pmdpar.LOOPS = nloops; 
		npt_cell.phdyn.set_vars(tao_npt, dt, ::max_ksi_Hoover);
		npt_cell.pmdpar.set_dt(dt);

		//npt_cell.phdyn.set_vars(tao_npt, 0, ::max_ksi_Hoover);
		//npt_cell.pmdpar.set_dt(0);

		npt_cell.Pext = 1; //200; // external pressure : 1 atmosphere

		if (cmm.nCheckCluster != 1) {
			sprintf(errmsg, "In NPT simulation, CMM will be updated for each step!"); show_log(errmsg, true);
			cmm.nCheckCluster = 1;
		}
	}

#if _MPI_ == 1
	sprintf(buffer, "%s%d.log", msave.title, ::mpi_id);
#else
	sprintf(buffer, "%s.log", msave.title);
#endif

	if (!mlog.init(buffer)) {
		sprintf(errmsg, "can not open file %s for log\nshall we continue with no log record ?", buffer);
#if _SYS_ == _WINDOWS_SYS_
		if (AfxMessageBox(LPCTSTR(CString(errmsg)), MB_YESNO) == IDNO) {
			md_run = false; mlog.close(); return;
		}
#else
		show_msg(errmsg); md_run = false; mlog.close(); return;
#endif
	}

	MD_SAVE pv_msave;
	if (MD_ENSEMBLE == MD_NPT) {
		sprintf(buffer, "%s-pv", msave.title);
		pv_msave.set_title(buffer);
		pv_msave.set_save(msave.nsave, msave.max_save);
		//pv_msave.set_save(1, msave.max_save * msave.nsave); // save each step for pressure
		npt_cell.pmdpar.psave.mdsave = &pv_msave;
	}

	if (!read_msg(mdpar_fname, buffer, "[RANDOM_FORCE]", true) ||
		sscanf(buffer, "%d %f", &n, &randf) != 2 || n <= 0) bRandomForce = false;
	else {bRandomForce = true; randf *= kT / U_MassAccel / dt; fwhm_randf = randf / 3;} // randf has unit kT/Angs.

	readExternalElectricField(mdpar_fname, bEext, Eext);
	readForceScale(mdpar_fname, scale_force, force_max, torque_max);
	if (scale_force) {
		sprintf(buffer, "scale force / torque when they are more than %7.4f / %7.4f", force_max, torque_max);
		show_infor(buffer);
		force2_max = force_max * force_max;
		torque2_max = torque_max * torque_max;
	}
	readLocalRelax(mdpar_fname, local_relax, dta_max_relax);
	if (::local_relax) {
		sprintf(buffer, "calculate the local relationship within neighbor CMM"); show_infor(buffer);
	}
	readLocalInteract(mdpar_fname, local_interact);
	::bRealEwaldSumOnly = (local_interact == 1 ? true : false);
	if (::bRealEwaldSumOnly) show_infor("real part Ewald-Sum is used for electrostatic interaction");
	readInternalKinEnergy(mdpar_fname, Ek0Internal);
	readBrownianPars(mdpar_fname, ::bBrownian, ::ita_Brownian, ::cd_drag_trans, ::cd_drag_rot);
	if (::bBrownian) {
		//::ita_Brownian *= kT / U_MassAccel / dt;
		//::cd_drag_trans *= kT / U_MassAccel / dt;
		//::cd_drag_rot *= kT / U_MassAccel / dt;
		::ita_Brownian /= dt;
		::cd_drag_rot /= dt;
		::cd_drag_trans /= dt;
	}

	readImplicitSolvent(mdpar_fname, ::bImplicitSolvent);
	if (::bImplicitSolvent) init_ImplicitSolvation_Pars();

	readPolarization(mdpar_fname, ::bPolarize, ::iPolarize);
	if (::bPolarize) readPolarization1(mdpar_fname, ::bExcludeInducedDipole);
	if (::bPolarize) mm_md_cell.eNeutral = false;
	mm_md_cell.check_eneutral();

	if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
		molecule_check_periodic_cell(mm_md_cell); // check the molecule mass center is inside the cell
	}

	// setup torsion of the macro-molecules
	dh_db.release();
	for (n = 0; n < mm_md_cell.mm.n; n++) {
		//mm_md_cell.mm.m[n].SetupHingeTorsion();
		initTorsionPars(mm_md_cell.mm.m + n);
#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0 && _EWALD_SUM == _MultiPole_EWALD_SUM
		for (nc = 0; nc < mm_md_cell.mm.m[n].nCluster; nc++) {
			mm_md_cell.mm.m[n].cluster[nc].calElectroMultipole(false);
		}
#endif
	}
	dh_db.convert_array();
	
	COMMAND_MD_STOP = false;
	md_run = true;
	construct_LJ12_6_prof(::iFormat_LJ);
	sprintf(errmsg, "normalized cutoff distance for 12-6 L-J interaction is set to %f", ::rcut_LJ_norm); show_log(errmsg, true);

	float rcut_SR = 0, rmin_SR = 0;
	if (strlen(sr_fname) > 0 && initPairwiseShortRangeInteract(sr_fname, ::atompar_db, ::sr_db, rcut_SR, rmin_SR)) {
		if (rcut_SR > ::rcut_LJ) {
			::rcut_LJ = rcut_SR; 
			::r2cut_LJ = ::rcut_LJ * ::rcut_LJ;
		}
	}
	sprintf(errmsg, "short-range-interaction cut-off distance is %f Angs.", ::rcut_LJ); show_log(errmsg, true);
	if (strlen(sr_fname) > 0) update_atomic_pars(sr_fname, ::atompar_db); // disable Lennard-Jones ? and update atomic polarizability

	if (::bPolarize) {
		switch (iPolarize) {
		case 1: // Thole polarization with exponential screening function
			_Thole_polarize_::initTholeRadius(::atompar_db, ::TholeRadius_db);
			break;
		}
	}

	mm_md_cell.SetMDMemory();
#if _MPI_ == 1
	MPI_MultiMM_JobAssign(mm_md_cell.mm.n, mm_md_cell.sm.n);
	mpi_multi_mol_md_proc(mm_md_cell, &cmm, ::MD_MODE);
#endif

	if (!readSPMEPars(mdpar_fname, ::indx_bspline, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], bExpandSPMECellZ)) return;

	extern void NPT_MD_PROC(ISOTROPIC_NPT_CELL &npt_cell);
	extern void PMF_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, char *pmf_parfname);

	extern void GPU_MD_MCONFS(MMOL_MD_CELL &mmcell, MC_SPECIES &mcsp, MD_PRO *mdpro, int npro);
	extern void SLAB_MD_PROCS(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, MD_PRO *mdpro, int npro);
	extern void SLAB_MD_MCONFS(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, MC_SPECIES &mcsp, MD_PRO *mdpro, int npro);
	extern bool read_mdprocs(char *fname, ARRAY<MD_PRO> &mdpro);
#ifdef __GPU__
	extern void GPU_SLAB_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm);
#endif


	ARRAY<MD_PRO> mdpro;
	read_mdprocs(mdpar_fname, mdpro);

	MC_SPECIES mcsp;
	bool bMCMD = read_mc_species(mdpar_fname, mcsp);
	if (bMCMD) {
		mcsp.xf[0] = -mm_md_cell.h[0]; mcsp.xt[0] = mm_md_cell.h[0];
		mcsp.xf[1] = -mm_md_cell.h[1]; mcsp.xt[1] = mm_md_cell.h[1];
	}

	if (mconstraint == 0) {
#ifdef __GPU__
		if (MD_ENSEMBLE == MD_NVT) {
			if (bMCMD && mcsp.sp.n > 0 && mdpro.n > 0) GPU_MD_MCONFS(mm_md_cell, mcsp, mdpro.m, mdpro.n);
			else GPU_MD_PROC(mm_md_cell, &cmm, ::MD_MODE); //MD_PROC(mm_md_cell, &cmm, ::MD_MODE);
		}
#else
#ifdef _USE_GPU_
		if (MD_ENSEMBLE == MD_NVT) {
			if (bMCMD && mcsp.sp.n > 0 && mdpro.n > 0) GPU_MD_MCONFS(mm_md_cell, mcsp, mdpro.m, mdpro.n);
			else GPU_MD_PROC(mm_md_cell, &cmm, ::MD_MODE);
		}
#else
		if (MD_ENSEMBLE == MD_NVT) {
			if (bMCMD && mcsp.sp.n > 0 && mdpro.n > 0) {
				show_msg("to use multi-configuration MD simulation, program has to enable _USE_GPU_ in project.h");
				//GPU_MD_PROCS(mm_md_cell, mdpro.m, mdpro.n);
			}
			else MD_PROC(mm_md_cell, &cmm, ::MD_MODE);
		}
#endif //_USE_GPU_
#endif //__GPU__
		else if (MD_ENSEMBLE == MD_NPT) NPT_MD_PROC(npt_cell);
	}
	else if (mconstraint == MD_PMF_CONSTRAINT) {
#ifdef __GPU__
		if (MD_ENSEMBLE == MD_NVT) {
			if (bMCMD && mcsp.sp.n > 0 && mdpro.n > 0) GPU_MD_MCONFS(mm_md_cell, mcsp, mdpro.m, mdpro.n);
			else GPU_MD_PROC(mm_md_cell, &cmm, ::MD_MODE);
		}
#else
#ifdef _USE_GPU_
		if (MD_ENSEMBLE == MD_NVT) {
			if (bMCMD && mcsp.sp.n > 0 && mdpro.n > 0) GPU_MD_MCONFS(mm_md_cell, mcsp, mdpro.m, mdpro.n);
			else GPU_MD_PROC(mm_md_cell, &cmm, ::MD_MODE); //PMF_MD_PROC(mm_md_cell, &cmm, constraint_parfname);
		}
#else
		if (MD_ENSEMBLE == MD_NVT) {
			if (bMCMD && mcsp.sp.n > 0 && mdpro.n > 0) {
				show_msg("to use multi-configuration MD simulation, program has to enable _USE_GPU_ in project.h");
				//GPU_MD_MCONFS(mm_md_cell, mcsp, mdpro.m, mdpro.n);
			}
			else MD_PROC(mm_md_cell, &cmm, ::MD_MODE); //PMF_MD_PROC(mm_md_cell, &cmm, constraint_parfname);
		}
#endif //_USE_GPU_
#endif //__GPU__
		else {
			show_log("For PMF-simulation, NPT essemble is not avaliable yet!", true);
		}
	}
	else if (mconstraint == MD_SLAB) {
		if (MD_ENSEMBLE == MD_NVT) {
			if (bMCMD && mcsp.sp.n > 0 && mdpro.n > 0) {
				SLAB_MD_MCONFS(mm_md_cell, &cmm, mcsp, mdpro.m, mdpro.n);
			}
#ifdef __GPU__
			else if (mdpro.n > 0) SLAB_MD_PROCS(mm_md_cell, &cmm, mdpro.m, mdpro.n);
			else GPU_SLAB_MD_PROC(mm_md_cell, &cmm);
#else
			else if (mdpro.n > 0) SLAB_MD_PROCS(mm_md_cell, &cmm, mdpro.m, mdpro.n);
			else SLAB_MD_PROC(mm_md_cell, &cmm);
#endif
		}
		else {
			show_log("For simulation of 2d slab, NPT essemble is not avaliable yet!", true);
		}
	}
	else if (mconstraint == MD_IMAGE_SLAB) {
		if (MD_ENSEMBLE == MD_NVT) ImageSlab_MD_PROC(mm_md_cell, &cmm);
		else {
			show_log("For simulation of 2d slab, NPT essemble is not avaliable yet!", true);
		}
	}
	/*
	if (mm_md_cell.sm.n > 0 || mm_md_cell.pm.n > 0) {
		FMD_PROC(mm_md_cell, &cmm);
	}
	else MD_PROC(mm_md_cell, &cmm, ::MD_MODE);
	*/

_MULTI_MD_OVER:
	md_run = false;
	cmm.set_cell(0, 0, 0);
	mm_md_cell.release();
	atompar_db.release(); // since all atoms are released, the database has to be released also
	LJPARS.release(); // release LJ pars 
	ESolv_db.release();
	//ESolvPars.release(); // solvation energy
	//torsion_db.release(); // release torsion database
	dh_db.release();
	mlog.close();
	cmm.reset_basic_cell_matrix();
	cluster_db.release();
	_interface_constraint_::clusterIC_db.release();

#ifdef __GPU__
	cudaDeviceReset(); 
#endif

	return;
}

void MD_SIM_MULTI_CGMM_CELL(char *mdpar_fname) {
#ifdef __GPU__
	cudaDeviceReset(); 
#endif

	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {show_msg(errmsg); mlog.close(); return;}

	// Leapfrog Verlet algorithm
	/*
	LFVERLET_CGMM_MD mdpar;
	if (!read_MD_BASE_Pars(mdpar_fname, (MD_BASE*)(&mdpar))) {show_msg(errmsg); return;}
	if (!read_HOOVER_Pars(mdpar_fname, (HOOVER_DYN*)(&mdpar))) {show_msg(errmsg); return;}
	mdpar.mmsave.mdsave = &msave;
	mdpar.LOOPS *= msave.nsave;
	::max_ksi_Hoover = 0.5f / mdpar.dt;
	if (mdpar.max_ksi_Hoover > ::max_ksi_Hoover) mdpar.max_ksi_Hoover = ::max_ksi_Hoover;

	cgmm_md_cell.SetMD(float(mdpar.dt), float(mdpar.tao_s), float(mdpar.max_ksi_Hoover), mdpar.LOOPS, ::bMMFrameConstraint);
	*/

	cgmm_md_cell.SetSave(&msave);

	// Nose-Hoover Chains
	float dt = 0, tao = 0, max_ksi = 0;
	if (!read_MD_TIMESTEP(mdpar_fname, dt)) {show_msg(errmsg); mlog.close(); return;}
	::timestep_md = dt;
	if (!read_HOOVER_Pars(mdpar_fname, tao, max_ksi)) {show_msg(errmsg); mlog.close(); return;}
	::max_ksi_Hoover = 0.6f / dt;
	if (max_ksi > ::max_ksi_Hoover) max_ksi = ::max_ksi_Hoover;
	cgmm_md_cell.set_md_timestep(dt);

	if (!ReadMDCellNHC<MMOLECULE, SMOLECULE, PMOLECULE, MMOL_MD_CELL>("NHC.INI", mm_md_cell, dt)) {
		mlog.close(); return;
	}
	cgmm_md_cell.hdyn.m[0].set_vars(tao, dt, max_ksi); // #0 -- default from MD.INI
	cgmm_md_cell.init_nhc_mols();
	cgmm_md_cell.check_freedom();
	// end of Nose-Hoover Chain

	long nloops = 0;
	if (!read_MD_LOOPS(mdpar_fname, nloops)) {show_msg(errmsg); mlog.close(); return;}
	nloops *= msave.nsave;
	cgmm_md_cell.SetMD(nloops);

	CMM_CELL3D<CG_CLUSTER> cmm;
	int interact_cal_method = 0;
	int n = 0, nc = 0, nCheckCell = 5;
	int nx, ny, nz;
	//float size_cluster = 0;//, rLJ_cut = 0, rshell = 0;;
	char buffer[256] = "\0";
	char* pch = NULL;

	if (!read_msg(mdpar_fname, buffer, "[CLUSTER_SIZE]", true) || sscanf(buffer, "%f", &size_cluster) != 1) size_cluster = 2;

	if (!read_msg(mdpar_fname, buffer, "[INTERACTION_CALCULATION_METHOD]", true) || 
		sscanf(buffer, "%d", &interact_cal_method) != 1 || interact_cal_method != 1) {
		sprintf(errmsg, "Interaction Calculation method for multi-molecule in cell has to be CMM");
		show_msg(errmsg); mlog.close(); md_run = false; return;
	}

#if _MPI_ == 1
	if (::MD_MODE != MULTI_CGMM_PERIDIC_CELL_MD) {
		sprintf(::errmsg, "MPI works on periodical cell only"); show_msg(::errmsg);
		mlog.close(); md_run = false; return;
	}
#endif

	if (interact_cal_method == _CMM) {
		pch = buffer;
		NextSection(&pch);
		if (sscanf(pch, "%d", &nCheckCell) != 1) {
			sprintf(errmsg, "can not get the loop number to re-check CMM during MD from %s", buffer);
			show_msg(errmsg); md_run = false; mlog.close(); return;
		}

		::min_cmm_size = size_cluster * 2;

		nx = int(2 * mm_md_cell.h[0] / min_cmm_size); nx = (nx < 1 ? 1 : nx);
		ny = int(2 * mm_md_cell.h[1] / min_cmm_size); ny = (ny < 1 ? 1 : ny);
		nz = int(2 * mm_md_cell.h[2] / min_cmm_size); nz = (nz < 1 ? 1 : nz);
		cmm.set_cell(nx, ny, nz);
		if (::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
			cmm.set_cell_range(true, -cgmm_md_cell.h[0], cgmm_md_cell.h[0], -cgmm_md_cell.h[1], cgmm_md_cell.h[1], -cgmm_md_cell.h[2], cgmm_md_cell.h[2]);
			cmm.set_cell_pos();
		}
		else if (::MD_MODE == MULTI_CGMM_FREE_CELL_MD) {
			cmm.set_cell_range(false, -cgmm_md_cell.h[0], cgmm_md_cell.h[0], -cgmm_md_cell.h[1], cgmm_md_cell.h[1], -cgmm_md_cell.h[2], cgmm_md_cell.h[2]);
			cmm.set_cell_pos();
		}

		cmm.nCheckCluster = nCheckCell;
	}

	readSurfaceBoundary(mdpar_fname);

#if _MPI_ == 1
	sprintf(buffer, "%s%d.log", msave.title, ::mpi_id);
#else
	sprintf(buffer, "%s.log", msave.title);
#endif

	if (!mlog.init(buffer)) {
		sprintf(errmsg, "can not open file %s for log\nshall we continue with no log record ?", buffer);
#if _SYS_ == _WINDOWS_SYS_
		if (AfxMessageBox(LPCTSTR(CString(errmsg)), MB_YESNO) == IDNO) {
			md_run = false; mlog.close(); return;
		}
#else
		show_msg(errmsg); md_run = false; mlog.close(); return;
#endif
	}

	if (!read_msg(mdpar_fname, buffer, "[RANDOM_FORCE]", true) ||
		sscanf(buffer, "%d %f", &n, &randf) != 2 || n <= 0) bRandomForce = false;
	else {bRandomForce = true; randf *= kT / U_MassAccel / dt; fwhm_randf = randf / 3;} // randf has unit kT/Angs.

	readExternalElectricField(mdpar_fname, bEext, Eext);
	if (!bEext || Eext == 0) ::bCGdipole = false;
	readForceScale(mdpar_fname, scale_force, force_max, torque_max);
	if (scale_force) {
		sprintf(buffer, "scale force / torque when they are more than %7.4f / %7.4f", force_max, torque_max);
		show_infor(buffer);
		force2_max = force_max * force_max;
		torque2_max = torque_max * torque_max;
	}
	//readLocalRelax(mdpar_fname, local_relax, dta_max_relax);
	readLocalInteract(mdpar_fname, local_interact);
	readInternalKinEnergy(mdpar_fname, Ek0Internal);
	readBrownianPars(mdpar_fname, ::bBrownian, ::ita_Brownian, ::cd_drag_trans, ::cd_drag_rot);
	if (::bBrownian) {
		//::ita_Brownian *= kT / U_MassAccel / dt;
		//::cd_drag_trans *= kT / U_MassAccel / dt;
		//::cd_drag_rot *= kT / U_MassAccel / dt;
		::ita_Brownian /= dt;
		::cd_drag_rot /= dt;
		::cd_drag_trans /= dt;
	}

	readImplicitSolvent(mdpar_fname, bImplicitSolvent);
	if (bImplicitSolvent) init_ImplicitSolvation_Pars();

	if (::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) {
		_coarse_grain_::molecule_check_periodic_cell(cgmm_md_cell); // check the molecule mass center is inside the cell
	}

	// setup torsion of the macro-molecules
	/*
	for (n = 0; n < cgmm_md_cell.mm.n; n++) {
		cgmm_md_cell.mm.m[n].SetupHingeTorsion();
		for (nc = 0; nc < cgmm_md_cell.mm.m[n].nCluster; nc++) {
			cgmm_md_cell.mm.m[n].cluster[nc].calElectroMultipole(false);
		}
	}
	*/

	// bound database -- neighbors were setup when cgmm is constructed
	// we have to setup neighbors of the cgcluster first
	cgbound_db.releaseChain(true); cgbound_db.releaseArray(true);
#if _CG_NEIGHBOR_INTERACT_ == 1
	int n_neighbors = 1;
	if (!get_neighbor_level(cgatom_par_fname, n_neighbors)) {md_run = false; return;}
#endif
	for (n = 0; n < cgmm_md_cell.mm.n; n++) {
		// if we use FENE bound interaction, each cluster has nearest neighbor only
	#if _CG_NEIGHBOR_INTERACT_ == 0
		cgmm_md_cell.mm.m[n].setup_cluster_neighbors(1); // 
		attach_bound_db(cgmm_md_cell.mm.m[n], cgbound_db, 0); // 0 -- FENE, 1 -- FROM FILE
	#elif _CG_NEIGHBOR_INTERACT_ == 1
		cgmm_md_cell.mm.m[n].setup_cluster_neighbors(n_neighbors); // 
		attach_bound_db(cgmm_md_cell.mm.m[n], cgbound_db, 1); // 0 -- FENE, 1 -- FROM FILE
		#if _INTRA_INTER_INTERACTS == 1
		attach_intramol_non_neighbor_db(::cgatompar_db, cg_intramol_non_neighbor_db);
		#endif
		attach_non_neighbor_db(::cgatompar_db, cg_non_neighbor_db);
	#endif
	}
	cgbound_db.Chain2Array();
	cg_intramol_non_neighbor_db.Chain2Array();
	cg_non_neighbor_db.Chain2Array();

	/*
	bool reset_cell_chain = true;
	for (n = 0; n < cgmm_md_cell.mm.n; n++) {
		reset_cell_chain = (n == 0 ? true : false);
		cmm_init_distribute_cluster(cgmm_md_cell.mm.m[n], cmm, reset_cell_chain); // set base cluster at the center cell
	}
	reset_cell_chain = (cgmm_md_cell.mm.n == 0 ? true : false);
	cmm_init_distribute_cluster(cgmm_md_cell.sm.m, cgmm_md_cell.sm.n, cmm, reset_cell_chain); // add simple solid molecule into cmm
	cmm_init_subcell<CG_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
	*/
	COMMAND_MD_STOP = false;
	md_run = true;
	construct_LJ12_6_prof(::iFormat_LJ);
	sprintf(errmsg, "normalized cutoff distance for 12-6 L-J interaction is set to %f", ::rcut_LJ_norm); show_log(errmsg, true);

	cgmm_md_cell.SetMDMemory();
#if _MPI_ == 1
	MPI_MultiMM_JobAssign(cgmm_md_cell.mm.n, cgmm_md_cell.sm.n);
	mpi_multi_mol_md_proc(cgmm_md_cell, &cmm, ::MD_MODE);
#else
	MD_PROC(cgmm_md_cell, &cmm, ::MD_MODE);
#endif

_MULTI_MD_OVER:
	md_run = false;
	cmm.set_cell(0, 0, 0);
	cgmm_md_cell.release();
	cgatompar_db.release(); // since all atoms are released, the database has to be released also
	LJPARS.release(); // release LJ pars 
	cgbound_db.releaseArray(true); cgbound_db.releaseChain(true);
	cg_dipole_distrbt_db.releaseArray(true); cg_dipole_distrbt_db.releaseChain(true);
	ESolv_db.release();
	//ESolvPars.release(); // solvation energy
	//torsion_db.release(); // release torsion database
	cluster_db.release();
	_interface_constraint_::clusterIC_db.release();
	mlog.close();
	return;
}
#if _DISABLE_
void MD_SIM_MULTI_VMOL_CELL(char *mdpar_fname) {
	// MD ESSEMBLE
	char buffer[256] = "\0";
	if (read_msg(mdpar_fname, buffer, "[MD_MODE]")) {
		if (sscanf(buffer, "%d", &MD_ENSEMBLE) != 1) MD_ENSEMBLE = MD_NVT;
		else if (MD_ENSEMBLE != MD_NVT && MD_ENSEMBLE != MD_NPT) MD_ENSEMBLE = MD_NVT;
	}

	// constraint
	int mconstraint = 0;
	char constraint_parfname[256] = "\0";
	if (read_msg(mdpar_fname, buffer, "[MD_CONSTRAINT]", true)) {
		if (sscanf(buffer, "%d", &mconstraint) == 1) {
			if (mconstraint > 0) {
				if (mconstraint == MD_PMF_CONSTRAINT) {
					if (sscanf(buffer, "%d %s", &mconstraint, constraint_parfname) != 2) {
						sprintf(errmsg, "constraint parameter file is not given, see : %s", buffer);
						show_log(errmsg, true); return;
					}
				}
				else if (mconstraint == MD_SLAB || mconstraint == MD_IMAGE_SLAB) {
					if (sscanf(buffer, "%d %s", &mconstraint, constraint_parfname) != 2) {
						sprintf(errmsg, "constraint parameter file is not given, see : %s", buffer);
						show_log(errmsg, true); return;
					}
					else strcpy(parfname_IConstraint, constraint_parfname);
				}
				else if (mconstraint == MD_NO_CONSTRAINT) {
					if (sscanf(buffer, "%d %s", &mconstraint, constraint_parfname) == 2)
						strcpy(parfname_IConstraint, constraint_parfname);
				}
				else {
					sprintf(errmsg, "ERROR: Unknown type of constraint. see : %s", buffer); 
					show_log(errmsg, true); return;
				}
			}
			else {
				show_log("no constraint on MD", true); mconstraint = 0;
			}
		}
		else {
			sprintf(errmsg, "no constraint model is given, see : %s", buffer); show_log(errmsg, true);
			show_log("disable constraint", true); mconstraint = 0;
		}
	}
	else mconstraint = 0;

	// MD SAVE
	MD_SAVE msave;
	if (!read_MD_SAVE_Pars(mdpar_fname, msave)) {show_msg(errmsg); mlog.close(); return;}

	// Leapfrog Verlet algorithm
	//LFVERLET_MM_MD mdpar;
	//if (!read_LFVerlet_MD_Pars(mdpar_fname, mdpar)) {show_msg(errmsg); mlog.close(); return;}

	int i;

	// Nose-Hoover Chains
	float dt = 0, tao = 0, max_ksi = 0;
	if (!read_MD_TIMESTEP(mdpar_fname, dt)) {show_msg(errmsg); mlog.close(); return;}
	::timestep_md = dt;
	if (!read_HOOVER_Pars(mdpar_fname, tao, max_ksi)) {show_msg(errmsg); mlog.close(); return;}
	//::max_ksi_Hoover = 0.5f / dt;
	::max_ksi_Hoover = 0.8f / dt;
	if (max_ksi > ::max_ksi_Hoover) max_ksi = ::max_ksi_Hoover;
	mm_md_cell.set_md_timestep(dt);

	if (!ReadMDCellNHC<MMOLECULE, SMOLECULE, PMOLECULE, MMOL_MD_CELL>("NHC.INI", mm_md_cell, dt)) {
		mlog.close(); return;
	}
	mm_md_cell.hdyn.m[mm_md_cell.hdyn.n - 1].set_vars(tao, dt, max_ksi); // #the last one -- default from MD.INI
	mm_md_cell.init_nhc_mols();
	mm_md_cell.check_freedom();
	// end of Nose-Hoover Chain

	// Nose-Hoover Chain for the constraint on mass-center velocity of the whole cell
	if (mconstraint != MD_SLAB && mconstraint != MD_IMAGE_SLAB) {
		mm_md_cell.setup_mcHDyn(tao, dt, max_ksi);
		mm_md_cell.mcHDyn.nhc->set_freedom(3);
	}
	mm_md_cell.TotalMass();

	long nloops = 0;
	if (!read_MD_LOOPS(mdpar_fname, nloops)) {show_msg(errmsg); mlog.close(); return;}
	nloops *= msave.nsave;
	mm_md_cell.setMD(nloops);
	
	mm_md_cell.SetSave(&msave);

	CMM_CELL3D<BASIC_CLUSTER> cmm;
	int interact_cal_method = 0;
	int n = 0, nc = 0, nCheckCell = 5;
	int nx, ny, nz;
	//float size_cluster = 0;//, rLJ_cut = 0, rshell = 0;;
	char* pch = NULL;

	if (!read_msg(mdpar_fname, buffer, "[CLUSTER_SIZE]", true) || sscanf(buffer, "%f", &size_cluster) != 1) size_cluster = 2;

	if (!read_msg(mdpar_fname, buffer, "[INTERACTION_CALCULATION_METHOD]", true) || 
		sscanf(buffer, "%d", &interact_cal_method) != 1 || interact_cal_method != 1) {
		sprintf(errmsg, "Interaction Calculation method for multi-molecule in cell has to be CMM"); show_msg(errmsg); 
		sprintf(errmsg, "CMM is used"); show_msg(errmsg); interact_cal_method = 1;
		//mlog.close(); md_run = false; return;
	}

#if _MPI_ == 1
	if (::MD_MODE != MULTI_MM_PERIDIC_CELL_MD) {
		sprintf(::errmsg, "MPI works on periodical cell only"); show_msg(::errmsg);
		mlog.close(); md_run = false; return;
	}
#endif

	if (interact_cal_method == _CMM) {
		pch = buffer;
		NextSection(&pch);
		if (sscanf(pch, "%d", &nCheckCell) != 1) {
			sprintf(errmsg, "can not get the loop number to re-check CMM during MD from %s", buffer);
			show_msg(errmsg); md_run = false; mlog.close(); return;
		}

		::min_cmm_size = size_cluster * 2;

		nx = int(2 * mm_md_cell.h[0] / min_cmm_size); nx = (nx < 1 ? 1 : nx);
		ny = int(2 * mm_md_cell.h[1] / min_cmm_size); ny = (ny < 1 ? 1 : ny);
		nz = int(2 * mm_md_cell.h[2] / min_cmm_size); nz = (nz < 1 ? 1 : nz);
		cmm.set_cell(nx, ny, nz);
		if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
			cmm.set_cell_range(true, -mm_md_cell.h[0], mm_md_cell.h[0], -mm_md_cell.h[1], mm_md_cell.h[1], -mm_md_cell.h[2], mm_md_cell.h[2]);
			cmm.set_cell_pos();
		}
		else if (::MD_MODE == MULTI_MM_FREE_CELL_MD) {
			cmm.set_cell_range(false, -mm_md_cell.h[0], mm_md_cell.h[0], -mm_md_cell.h[1], mm_md_cell.h[1], -mm_md_cell.h[2], mm_md_cell.h[2]);
			cmm.set_cell_pos();
		}

		cmm.nCheckCluster = nCheckCell;
	}

	readSurfaceBoundary(mdpar_fname);

	ISOTROPIC_NPT_CELL npt_cell;
	float W_npt = 0, tao_npt = 0;
	if (MD_ENSEMBLE == MD_NPT) {
		npt_cell.mdcell = &mm_md_cell; npt_cell.cmm = &cmm;
		npt_cell.NF_pressure = check_NF_pressure(npt_cell.mdcell);
		if (!readNPTPars(mdpar_fname, W_npt, tao_npt)) {md_run = false; mlog.close(); return;}
		npt_cell.pdyn.W = W_npt * W_npt * (npt_cell.NF_pressure + 3);
		npt_cell.setup_pressure_NHC();
		npt_cell.phdyn.nhc->N = 9; // d^2
		npt_cell.pmdpar.LOOPS = nloops; 
		npt_cell.phdyn.set_vars(tao_npt, dt, ::max_ksi_Hoover);
		npt_cell.pmdpar.set_dt(dt);

		//npt_cell.phdyn.set_vars(tao_npt, 0, ::max_ksi_Hoover);
		//npt_cell.pmdpar.set_dt(0);

		npt_cell.Pext = 1; //200; // external pressure : 1 atmosphere

		if (cmm.nCheckCluster != 1) {
			sprintf(errmsg, "In NPT simulation, CMM will be updated for each step!"); show_log(errmsg, true);
			cmm.nCheckCluster = 1;
		}
	}

#if _MPI_ == 1
	sprintf(buffer, "%s%d.log", msave.title, ::mpi_id);
#else
	sprintf(buffer, "%s.log", msave.title);
#endif

	if (!mlog.init(buffer)) {
		sprintf(errmsg, "can not open file %s for log\nshall we continue with no log record ?", buffer);
#if _SYS_ == _WINDOWS_SYS_
		if (AfxMessageBox(LPCTSTR(CString(errmsg)), MB_YESNO) == IDNO) {
			md_run = false; mlog.close(); return;
		}
#else
		show_msg(errmsg); md_run = false; mlog.close(); return;
#endif
	}

	MD_SAVE pv_msave;
	if (MD_ENSEMBLE == MD_NPT) {
		sprintf(buffer, "%s-pv", msave.title);
		pv_msave.set_title(buffer);
		pv_msave.set_save(msave.nsave, msave.max_save);
		//pv_msave.set_save(1, msave.max_save * msave.nsave); // save each step for pressure
		npt_cell.pmdpar.psave.mdsave = &pv_msave;
	}

	if (!read_msg(mdpar_fname, buffer, "[RANDOM_FORCE]", true) ||
		sscanf(buffer, "%d %f", &n, &randf) != 2 || n <= 0) bRandomForce = false;
	else {bRandomForce = true; randf *= kT / U_MassAccel / dt; fwhm_randf = randf / 3;} // randf has unit kT/Angs.

	readExternalElectricField(mdpar_fname, bEext, Eext);
	readForceScale(mdpar_fname, scale_force, force_max, torque_max);
	if (scale_force) {
		sprintf(buffer, "scale force / torque when they are more than %7.4f / %7.4f", force_max, torque_max);
		show_infor(buffer);
		force2_max = force_max * force_max;
		torque2_max = torque_max * torque_max;
	}
	readLocalRelax(mdpar_fname, local_relax, dta_max_relax);
	if (::local_relax) {
		sprintf(buffer, "calculate the local relationship within neighbor CMM"); show_infor(buffer);
	}
	readLocalInteract(mdpar_fname, local_interact);
	readInternalKinEnergy(mdpar_fname, Ek0Internal);
	readBrownianPars(mdpar_fname, ::bBrownian, ::ita_Brownian, ::cd_drag_trans, ::cd_drag_rot);
	if (::bBrownian) {
		//::ita_Brownian *= kT / U_MassAccel / dt;
		//::cd_drag_trans *= kT / U_MassAccel / dt;
		//::cd_drag_rot *= kT / U_MassAccel / dt;
		::ita_Brownian /= dt;
		::cd_drag_rot /= dt;
		::cd_drag_trans /= dt;
	}

	readImplicitSolvent(mdpar_fname, ::bImplicitSolvent);
	if (::bImplicitSolvent) init_ImplicitSolvation_Pars();

	readPolarization(mdpar_fname, ::bPolarize, ::iPolarize);
	if (::bPolarize) readPolarization1(mdpar_fname, ::bExcludeInducedDipole);
	if (::bPolarize) mm_md_cell.eNeutral = false;
	mm_md_cell.check_eneutral();

	if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD) {
		molecule_check_periodic_cell(mm_md_cell); // check the molecule mass center is inside the cell
	}

	// setup torsion of the macro-molecules
	dh_db.release();
	for (n = 0; n < mm_md_cell.mm.n; n++) {
		//mm_md_cell.mm.m[n].SetupHingeTorsion();
		initTorsionPars(mm_md_cell.mm.m + n);
#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0 && _EWALD_SUM == _MultiPole_EWALD_SUM
		for (nc = 0; nc < mm_md_cell.mm.m[n].nCluster; nc++) {
			mm_md_cell.mm.m[n].cluster[nc].calElectroMultipole(false);
		}
#endif
	}
	dh_db.convert_array();
	
	COMMAND_MD_STOP = false;
	md_run = true;
	construct_LJ12_6_prof(::iFormat_LJ);
	sprintf(errmsg, "normalized cutoff distance for 12-6 L-J interaction is set to %f", ::rcut_LJ_norm); show_log(errmsg, true);

	float rcut_SR = 0, rmin_SR = 0;
	if (strlen(sr_fname) > 0 && initPairwiseShortRangeInteract(sr_fname, ::atompar_db, ::sr_db, rcut_SR, rmin_SR)) {
		if (rcut_SR > ::rcut_LJ) {
			::rcut_LJ = rcut_SR; 
			::r2cut_LJ = ::rcut_LJ * ::rcut_LJ;
		}
	}
	sprintf(errmsg, "short-range-interaction cut-off distance is %f Angs.", ::rcut_LJ); show_log(errmsg, true);
	if (strlen(sr_fname) > 0) update_atomic_pars(sr_fname, ::atompar_db); // disable Lennard-Jones ? and update atomic polarizability

	if (::bPolarize) {
		switch (iPolarize) {
		case 1: // Thole polarization with exponential screening function
			_Thole_polarize_::initTholeRadius(::atompar_db, ::TholeRadius_db);
			break;
		}
	}

	mm_md_cell.SetMDMemory();
#if _MPI_ == 1
	MPI_MultiMM_JobAssign(mm_md_cell.mm.n, mm_md_cell.sm.n);
	mpi_multi_mol_md_proc(mm_md_cell, &cmm, ::MD_MODE);
#else

	if (!readSPMEPars(mdpar_fname, ::indx_bspline, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], bExpandSPMECellZ)) return;

	extern void NPT_MD_PROC(ISOTROPIC_NPT_CELL &npt_cell);
	extern void PMF_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, char* constraint_parfname);

	if (mconstraint == 0) {
		if (MD_ENSEMBLE == MD_NVT) MD_PROC(mm_md_cell, &cmm, ::MD_MODE);
		else if (MD_ENSEMBLE == MD_NPT) NPT_MD_PROC(npt_cell);
	}
	else if (mconstraint == MD_PMF_CONSTRAINT) {
		if (MD_ENSEMBLE == MD_NVT) {
			show_log("For PMF-simulation, NVT essemble is not avaliable yet!", true);
			//PMF_MD_PROC(mm_md_cell, &cmm, constraint_parfname);
		}
		else {
			show_log("For PMF-simulation, NPT essemble is not avaliable yet!", true);
		}
	}
	/*
	if (mm_md_cell.sm.n > 0 || mm_md_cell.pm.n > 0) {
		FMD_PROC(mm_md_cell, &cmm);
	}
	else MD_PROC(mm_md_cell, &cmm, ::MD_MODE);
	*/
#endif

_MULTI_MD_OVER:
	md_run = false;
	cmm.set_cell(0, 0, 0);
	mm_md_cell.release();
	atompar_db.release(); // since all atoms are released, the database has to be released also
	LJPARS.release(); // release LJ pars 
	ESolv_db.release();
	//ESolvPars.release(); // solvation energy
	//torsion_db.release(); // release torsion database
	dh_db.release();
	mlog.close();
	cmm.reset_basic_cell_matrix();
	cluster_db.release();
	_interface_constraint_::clusterIC_db.release();

#ifdef __GPU__
	cudaDeviceReset(); 
#endif

	return;
}
#endif

extern "C" void MD_SIM(char *mdpar_fname) {
	switch (MD_MODE) {
	case SINGLE_MM_FREE_CELL_MD:
		return MD_SIM_SINGLE_MOL(mdpar_fname);
		break;
	case SINGLE_CGMM_FREE_CELL_MD:
		return MD_SIM_SINGLE_CGMM(mdpar_fname);
		break;
	case MULTI_MM_FREE_CELL_MD:
		return MD_SIM_MULTI_MOL_CELL(mdpar_fname);
		break;
	case MULTI_MM_PERIDIC_CELL_MD:
		return MD_SIM_MULTI_MOL_CELL(mdpar_fname);
		break;
	case MULTI_CGMM_FREE_CELL_MD:
		return MD_SIM_MULTI_CGMM_CELL(mdpar_fname);
		break;
	case MULTI_CGMM_PERIDIC_CELL_MD:
		return MD_SIM_MULTI_CGMM_CELL(mdpar_fname);
		break;
	}
}

extern "C" bool running_md(void) {
	return md_run;
}

extern "C" void stop_md(bool &status) {
	COMMAND_MD_STOP = true;
}

extern "C" int get_molecule_struct(float *x, float *y, float *z, bool *highlight, int MAX, int nMol, int nstart) {
	int i = 0, j = 0;
	int iatom = 0;
	int nc = 0;
	int n = 0;
	MATOM *patom = NULL;
	MMOLECULE *m = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	CGATOM *pcgatom = NULL;
	CG_MMOLECULE *cgmm = NULL;
	CGSM *cgsm = NULL;
	GENERAL_MOLECULE gm;
	bool status = get_molecule_upto_job(nMol, gm);
	if (!status) return -1;

	if (gm.mol_type == 1) { // macro-molecule
		m = gm.mm;
		while (nc < m->nCluster && j + m->cluster[nc].nAtoms <= nstart) {j += m->cluster[nc].nAtoms; nc++;}
		if (nc == m->nCluster) {n = 0; return n;} // molecule has no more atoms
		j = nstart - j; i = 0;
		while (i < MAX) {
			for (iatom = j; iatom < m->cluster[nc].nAtoms; iatom++) {
				patom = m->cluster[nc].atom + iatom;
				x[i] = float(patom->r.v[0] + m->cluster[nc].dr_cmm.v[0]); 
				y[i] = float(patom->r.v[1] + m->cluster[nc].dr_cmm.v[1]);
				z[i] = float(patom->r.v[2] + m->cluster[nc].dr_cmm.v[2]);
				highlight[i] = m->cluster[nc].highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			nc++;
			if (nc == m->nCluster) break; //over
			j = 0; // for new cluster, atom start from 0
		}
		return n;
	}
	else if (gm.mol_type == 0) { // simple molecule
		sm = gm.sm;
		if (sm->c == NULL || nstart >= sm->c->nAtoms) return -1;
		j = nstart; i = 0;
		while (i < MAX) {
			for (iatom = j; iatom < sm->c->nAtoms; iatom++) {
				patom = sm->c->atom + iatom;
				x[i] = float(patom->r.v[0] + sm->c->dr_cmm.v[0]); 
				y[i] = float(patom->r.v[1] + sm->c->dr_cmm.v[1]); 
				z[i] = float(patom->r.v[2] + sm->c->dr_cmm.v[2]);
				highlight[i] = sm->c->highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			break;
		}
		return n;
	}
	else if (gm.mol_type == 4) { // point molecule
		pm = gm.pm;
		if (pm->c == NULL || nstart >= pm->c->nAtoms) return -1;
		j = nstart; i = 0;
		while (i < MAX) {
			for (iatom = j; iatom < pm->c->nAtoms; iatom++) {
				patom = pm->c->atom + iatom;
				x[i] = float(patom->r.v[0] + pm->c->dr_cmm.v[0]); 
				y[i] = float(patom->r.v[1] + pm->c->dr_cmm.v[1]); 
				z[i] = float(patom->r.v[2] + pm->c->dr_cmm.v[2]);
				highlight[i] = pm->c->highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			break;
		}
		return n;
	}
	else if (gm.mol_type == 2) { // coarse-grained macromolecule
		cgmm = gm.cgmm;
		while (nc < cgmm->nCluster && j + cgmm->cluster[nc].nAtoms <= nstart) {j += cgmm->cluster[nc].nAtoms; nc++;}
		if (nc == cgmm->nCluster) {n = 0; return n;} // molecule has no more atoms
		j = nstart - j; i = 0;
		while (i < MAX) {
			for (iatom = j; iatom < cgmm->cluster[nc].nAtoms; iatom++) {
				pcgatom = cgmm->cluster[nc].atom + iatom;
				x[i] = float(pcgatom->r.v[0] + cgmm->cluster[nc].dr_cmm.v[0]); 
				y[i] = float(pcgatom->r.v[1] + cgmm->cluster[nc].dr_cmm.v[1]);
				z[i] = float(pcgatom->r.v[2] + cgmm->cluster[nc].dr_cmm.v[2]);
				highlight[i] = cgmm->cluster[nc].highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			nc++;
			if (nc == cgmm->nCluster) break; //over
			j = 0; // for new cluster, atom start from 0
		}
		return n;
	}
	else if (gm.mol_type == 3) { // coarse-grained simple molecule
		cgsm = gm.cgsm;
		j = nstart; i = 0;
		while (i < MAX) {
			for (iatom = j; iatom < sm->c->nAtoms; iatom++) {
				pcgatom = cgsm->atom + iatom;
				x[i] = float(pcgatom->r.v[0] + cgsm->dr_cmm.v[0]); 
				y[i] = float(pcgatom->r.v[1] + cgsm->dr_cmm.v[1]); 
				z[i] = float(pcgatom->r.v[2] + cgsm->dr_cmm.v[2]);
				highlight[i] = cgsm->highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			break;
		}
		return n;
	}
	else return -1;
}

extern "C" int get_molecule(int *atom, float *x, float *y, float *z, bool *highlight, int MAX, int nMol, int nstart) {
	int i = 0, j = 0;
	int iatom = 0;
	int nc = 0;
	int n = 0;
	MATOM *patom = NULL;
	MMOLECULE *m = NULL;
	CGATOM *pcgatom = NULL;
	CG_MMOLECULE *cgmm = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	CGSM *cgsm = NULL;
	GENERAL_MOLECULE gm;
	bool status = get_molecule_upto_job(nMol, gm);
	if (!status) return -1;

	if (gm.mol_type == 1) { // macro-molecule
		m = gm.mm;
		while (nc < m->nCluster && j + m->cluster[nc].nAtoms <= nstart) {j += m->cluster[nc].nAtoms; nc++;}
		if (nc == m->nCluster) {n = 0; return n;} // molecule has no more atoms
		j = nstart - j;
		while (i < MAX) {
			for (iatom = j; iatom < m->cluster[nc].nAtoms; iatom++) {
				patom = m->cluster[nc].atom + iatom;
				//atom[i] = patom->aindx;
				atom[i] = patom->par->uid;
				x[i] = float(patom->r.v[0] + m->cluster[nc].dr_cmm.v[0]);
				y[i] = float(patom->r.v[1] + m->cluster[nc].dr_cmm.v[1]);
				z[i] = float(patom->r.v[2] + m->cluster[nc].dr_cmm.v[2]);
				highlight[i] = m->cluster[nc].highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			nc++;
			if (nc == m->nCluster) break; //over
			j = 0; // for new cluster, atom start from 0
		}
		return n;
	}
	else if (gm.mol_type == 0) { // simple-molecule
		sm = gm.sm;
		if (sm->c == NULL || nstart >= sm->c->nAtoms) return -1;
		j = nstart;
		while (i < MAX) {
			for (iatom = j; iatom < sm->c->nAtoms; iatom++) {
				patom = sm->c->atom + iatom;
				//atom[i] = patom->aindx;
				atom[i] = patom->par->uid;
				x[i] = float(patom->r.v[0] + sm->c->dr_cmm.v[0]); 
				y[i] = float(patom->r.v[1] + sm->c->dr_cmm.v[1]);
				z[i] = float(patom->r.v[2] + sm->c->dr_cmm.v[2]);
				highlight[i] = sm->c->highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			break;
		}
		return n;
	}
	else if (gm.mol_type == 4) { // point-molecule
		pm = gm.pm;
		if (pm->c == NULL || nstart >= pm->c->nAtoms) return -1;
		j = nstart;
		while (i < MAX) {
			for (iatom = j; iatom < pm->c->nAtoms; iatom++) {
				patom = pm->c->atom + iatom;
				//atom[i] = patom->aindx;
				atom[i] = patom->par->uid;
				x[i] = float(patom->r.v[0] + pm->c->dr_cmm.v[0]); 
				y[i] = float(patom->r.v[1] + pm->c->dr_cmm.v[1]);
				z[i] = float(patom->r.v[2] + pm->c->dr_cmm.v[2]);
				highlight[i] = pm->c->highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			break;
		}
		return n;
	}
	else if (gm.mol_type == 2) { // coarse-grained macro-molecule
		cgmm = gm.cgmm;
		while (nc < cgmm->nCluster && j + cgmm->cluster[nc].nAtoms <= nstart) {j += cgmm->cluster[nc].nAtoms; nc++;}
		if (nc == cgmm->nCluster) {n = 0; return n;} // molecule has no more atoms
		j = nstart - j;
		while (i < MAX) {
			for (iatom = j; iatom < cgmm->cluster[nc].nAtoms; iatom++) {
				pcgatom = cgmm->cluster[nc].atom + iatom;
				//atom[i] = patom->aindx;
				atom[i] = pcgatom->par->uid;
				x[i] = float(pcgatom->r.v[0] + cgmm->cluster[nc].dr_cmm.v[0]);
				y[i] = float(pcgatom->r.v[1] + cgmm->cluster[nc].dr_cmm.v[1]);
				z[i] = float(pcgatom->r.v[2] + cgmm->cluster[nc].dr_cmm.v[2]);
				highlight[i] = cgmm->cluster[nc].highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			nc++;
			if (nc == cgmm->nCluster) break; //over
			j = 0; // for new cluster, atom start from 0
		}
		return n;
	}
	else if (gm.mol_type == 3) { // coarse-grained simple-molecule
		cgsm = gm.cgsm;
		j = nstart;
		while (i < MAX) {
			for (iatom = j; iatom < cgsm->nAtoms; iatom++) {
				pcgatom = cgsm->atom + iatom;
				//atom[i] = patom->aindx;
				atom[i] = pcgatom->par->uid;
				x[i] = float(pcgatom->r.v[0] + cgsm->dr_cmm.v[0]); 
				y[i] = float(pcgatom->r.v[1] + cgsm->dr_cmm.v[1]);
				z[i] = float(pcgatom->r.v[2] + cgsm->dr_cmm.v[2]);
				highlight[i] = cgsm->highlight;
				i++; n++;
				if (i == MAX) return n;
			}
			break;
		}
		return n;
	}
	else return -1;
}

extern "C" int get_atom_bound(int *nb, char *bv, char *cvb_type, int MAX, int nMol, int natom) {
	int i = 0, j = 0;
	int iatom = 0;
	int nc = 0;
	int n = 0;
	MATOM *patom = NULL;
	MMOLECULE *m = NULL;
	CGATOM *pcgatom = NULL;
	CG_MMOLECULE *cgmm = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	CGSM *cgsm = NULL;
	GENERAL_MOLECULE gm;
	bool status = get_molecule_upto_job(nMol, gm);
	if (!status) return -1;

	if (gm.mol_type == 1) { // macro-molecule
		m = gm.mm;
		while (nc < m->nCluster && j + m->cluster[nc].nAtoms <= natom) {j += m->cluster[nc].nAtoms; nc++;}
		if (nc == m->nCluster) {n = 0; return n;} // molecule has no more atoms
		j = natom - j;
		patom = m->cluster[nc].atom + j;
		for (i = 0; i < patom->nb; i++) {
			if (i == MAX) return i;
			if (!atom_indx(*m, (MATOM*)patom->ba[i], nc, iatom, nb[i])) {
				sprintf(errmsg, "ERROR : can not find atom defined in cluster");
				show_msg(errmsg);
			}
			bv[i] = patom->bv[i]; cvb_type[i] = patom->btype[i];
		}
		return patom->nb;
	}
	else if (gm.mol_type == 0) { // simple-molecule
		sm = gm.sm;
		if (sm->c == NULL || natom >= sm->c->nAtoms) return -1;

		j = natom;
		patom = sm->c->atom + j;
		for (i = 0; i < patom->nb; i++) {
			if (i == MAX) return i;
			if (!atom_indx(*sm, (MATOM*)patom->ba[i], iatom, nb[i])) {
				sprintf(errmsg, "ERROR : can not find atom defined in cluster");
				show_msg(errmsg);
			}
			bv[i] = patom->bv[i]; cvb_type[i] = patom->btype[i];
		}
		return patom->nb;
	}
	else if (gm.mol_type == 4) { // point-molecule
		pm = gm.pm;
		if (pm->c == NULL || natom >= pm->c->nAtoms) return -1;

		j = natom;
		patom = pm->c->atom + j;
		for (i = 0; i < patom->nb; i++) {
			if (i == MAX) return i;
			if (!atom_indx(*sm, (MATOM*)patom->ba[i], iatom, nb[i])) {
				sprintf(errmsg, "ERROR : can not find atom defined in cluster");
				show_msg(errmsg);
			}
			bv[i] = patom->bv[i]; cvb_type[i] = patom->btype[i];
		}
		return patom->nb;
	}
	else if (gm.mol_type == 2) { // coarse-grained macro-molecule
		cgmm = gm.cgmm;
		while (nc < cgmm->nCluster && j + cgmm->cluster[nc].nAtoms <= natom) {j += cgmm->cluster[nc].nAtoms; nc++;}
		if (nc == cgmm->nCluster) {n = 0; return n;} // molecule has no more atoms
		j = natom - j;
		pcgatom = cgmm->cluster[nc].atom + j;
		for (i = 0; i < pcgatom->nb; i++) {
			if (i == MAX) return i;
			if (!atom_indx(*cgmm, (CGATOM*)pcgatom->ba[i], nc, iatom, nb[i])) {
				sprintf(errmsg, "ERROR : can not find atom defined in cluster");
				show_msg(errmsg);
			}
			bv[i] = pcgatom->bv[i]; cvb_type[i] = pcgatom->btype[i];
		}
		return pcgatom->nb;
	}
	else if (gm.mol_type == 3) { // coarse-grained simple-molecule
		/*
		cgsm = gm.cgsm;
		j = natom;
		pcgatom = cgsm->atom + j;
		for (i = 0; i < pcgatom->nb; i++) {
			if (i == MAX) return i;
			if (!atom_indx(*sm, (MATOM*)patom->ba[i], iatom, nb[i])) {
				sprintf(errmsg, "ERROR : can not find atom defined in cluster");
				show_msg(errmsg);
			}
			bv[i] = patom->bv[i]; cvb_type[i] = patom->btype[i];
		}
		return patom->nb;
		*/
		return 0;
	}
	else return -1;
}

extern "C" void get_errmsg(char *errmsg) {
	strcpy(errmsg, ::errmsg);
}

extern "C" void get_mol_struct_over() {
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(hWaitGetMolStructEvent);
#endif
}

extern "C" void set_show_infor(int min_time, int max_time, int nshow) {
	::delay_time = min_time;
	::max_time_pickup_molstruct = max_time - min_time;
	if (::max_time_pickup_molstruct < 100) ::max_time_pickup_molstruct = 100;
	::nshow = nshow;
}

extern "C" void get_show_infor(int* min_time, int *max_time, int *nshow) {
	*min_time = ::delay_time;
	*max_time = ::max_time_pickup_molstruct;
	*nshow = ::nshow;
}

extern "C" bool get_single_md_proc(char *mdpar_fname, int isave) {
	switch (JOB) {
	case GET_SINGLE_MM_MD_PROC:
		return get_single_mm_md_proc(mdpar_fname, isave);
		break;
	case GET_MULTI_MM_MD_PROC:
		return get_multi_mm_md_proc(mdpar_fname, isave);
		break;
	case GET_SINGLE_CGMM_MD_PROC:
		return get_single_cgmm_md_proc(mdpar_fname, isave);
		break;
	case GET_MULTI_CGMM_MD_PROC:
		return get_multi_cgmm_md_proc(mdpar_fname, isave);
		break;
	default:
		return false;
	}
}

extern "C" bool show_single_md_proc(char *mdpar_fname, int nf) {
	switch (JOB) {
	case SHOW_SINGLE_MM_MD_PROC:
		return show_single_mm_md_proc(mdpar_fname, nf);
		break;
	case SHOW_MULTI_MM_MD_PROC:
		return show_multi_mm_md_proc(mdpar_fname, nf);
		break;
	case SHOW_SINGLE_CGMM_MD_PROC:
		return show_single_cgmm_md_proc(mdpar_fname, nf);
		break;
	case SHOW_MULTI_CGMM_MD_PROC:
		return show_multi_cgmm_md_proc(mdpar_fname, nf);
		break;
	default:
		return false;
	}
}

extern READ_MD_PROC_FILE rmdproc;
extern PICK_SCATT_STRUCT pick_struct;

void clear_memory() {
	mm.reset();
	mm_md_cell.release();
	release_list<MMOLECULE>(&mlist, true);

	rmdproc.reset_opened_file();
	pick_struct.reset();

#if _MPI_ == 1
	reset_mpi_pars();
#endif
}

extern "C" void periodic_cell_check() {
	if (::MD_MODE == MULTI_MM_PERIDIC_CELL_MD) molecule_check_periodic_cell(mm_md_cell);
	else if (::MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) _coarse_grain_::molecule_check_periodic_cell(cgmm_md_cell);
}
