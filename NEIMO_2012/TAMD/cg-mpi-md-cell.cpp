#include "project.h"
#if _MPI_ == 1
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
#include "cg-mm.h"
#include "CMM.h"
#include "cg-cmm.h"
#include "var.h"

#include "Interaction.h"
#include "Interact1.h"
#include "Interact2.h"
#include "Interact3.h"
#include "Interact4.h"
#include "MD.h"
#include "cg-md.h"
//#include "MD_POLYMER.h"
#include "cgmd_polymer.h"
#include "md_cell.h"
#include "cg-md-cell.h"

#ifdef SEEK_SET
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include "mpi.h"
#endif

#include "mpi_md_cell.h"

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

namespace _coarse_grain_ {

void MPI_GatherCGClusterForce(CGMM_MD_CELL &mmcell) {
	if (mmcell.mm.m == NULL) return;
#define NC_MAX   200
#define V_MAX    2000
	CG_MMOLECULE *mm = NULL;

	int iMol = 0, n1Mol = 0, n2Mol = 0;
	int nc = 0, nv = 0, ncs = 0;
	double dbf[V_MAX];
	CG_CLUSTER *pc = NULL;
	double *fc = NULL;

	int isource = 0, ipkg = 0, npkgs = 0;
	int nSend[3] = {0, 0, 0};
	int total_clusters = 0;
	MPI_Status status;
	int nget = 0;
	char msg[256] = "\0";

	if (::mpi_id != 0) { //this is not the master
		// macromolecules in this process
		nSend[0] = ::n1MM[::mpi_id]; nSend[1] = ::n2MM[::mpi_id];
		MPI_Send(nSend, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);

		for (iMol = ::n1MM[::mpi_id]; iMol <= ::n2MM[::mpi_id]; iMol++) {
			// macromolecule indx, number of clusters, number of sending packages
			nSend[0] = iMol; nSend[1] = mmcell.mm.m[iMol].nCluster;
			nSend[2] = mmcell.mm.m[iMol].nCluster / NC_MAX;
			if (nSend[2] * NC_MAX < total_clusters) nSend[2] ++;
			MPI_Send(nSend, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);

			mm = mmcell.mm.m + iMol;
			nc = 0;
			while (nc < mm->nCluster) {
				ncs = 0; nv = 0;
				while (ncs < NC_MAX && nc < mm->nCluster) {
					if (nc >= mm->nCluster) break;
					pc = mm->cluster + nc; fc = pc->dyn->fc.v;
					dbf[nv+0] = fc[0]; dbf[nv+1] = fc[1]; dbf[nv+2] = fc[2]; 
					nv += 3; ncs++; nc++;
				}
				MPI_Send(dbf, nv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // send to the master
			}
		}
	}
	else {
		for (isource = 1; isource < ::nMPI; isource++) {
			MPI_Recv(nSend, 2, MPI_INT, isource, 0, MPI_COMM_WORLD, &status);
			n1Mol = nSend[0]; n2Mol = nSend[1];
			
			for (iMol = n1Mol; iMol <= n2Mol; iMol++) {
				mm = mmcell.mm.m + iMol;
				MPI_Recv(nSend, 3, MPI_INT, isource, 0, MPI_COMM_WORLD, &status);
				if (iMol != nSend[0] || mm->nCluster != nSend[1]) {
					sprintf(msg, "ERROR : the molecule information is not consistent, contacting with processor #%d, at macromolecule #%d", isource, iMol);
					show_msg(msg);
				}
				npkgs = nSend[2];
				nc = 0;
				for (ipkg = 0; ipkg < npkgs; ipkg++) {
					MPI_Recv(dbf, V_MAX, MPI_DOUBLE, isource, 0, MPI_COMM_WORLD, &status);
					ncs = 0; MPI_Get_count(&status, MPI_DOUBLE, &nget);
					while (ncs < nget) {
						if (nc >= mm->nCluster) {
							sprintf(msg, "ERROR : MPI get force of cluster %d", nc); show_log(msg, true);
						}
						else {
							pc = mm->cluster + nc; fc = pc->dyn->fc.v;
							fc[0] = dbf[ncs]; fc[1] = dbf[ncs+1]; fc[2] = dbf[ncs+2];
							nc++;
						}
						ncs += 3;
					}
				}
			}
		}
	}

	
#undef NC_MAX
#undef V_MAX
}


bool MPI_ExchangeStructInfor(CGMM_MD_CELL &mmcell) {
	CG_MMOLECULE *mm = NULL;
	CG_CLUSTER *pc = NULL;
	CGATOM *patom = NULL;
	double *dp = NULL;
#define NC_MAX   2000
	int iMol = 0, nc = 0, proc_id = 0;
	int i = 0, iv = 0;
	double dbf[NC_MAX];

	MPI_Status status;
	char msg[256] = "\0";

	for (iMol = 0; iMol < mmcell.mm.n; iMol++) {
		proc_id = MM_mpi_process(iMol);
		mm = mmcell.mm.m + iMol;
		//sprintf(msg, "mm %d with cluster %d : broadcast from %d", iMol, mm->nCluster, proc_id); show_log(msg, true);
		// cluster's parameters
		for (nc = 0; nc < mm->nCluster; nc++) {
			//MPI_Barrier(MPI_COMM_WORLD);
			pc = mm->cluster + nc;
			iv = 0;
			if (::mpi_id == proc_id) {
				for (i = 0; i < pc->nAtoms; i++) {
					patom = pc->atom + i;
					dp = dbf + iv; v3cp(dp, patom->r.v) iv += 3;
				}
				/*
				for (i = 0; i < pc->nHinges; i++) {
					dp = dbf + iv; v3cp(dp, pc->hinge[i].vHinge.v) iv += 3;
				}
				*/
				//dp = dbf + iv; v3cp(dp, pc->dr_cmm.v) iv += 3; // dr_cmm
				//dp = dbf + iv; v3cp(dp, pc->bkm->v.v) iv += 3; // v
				if (::bImplicitSolvent) {
					dp = dbf + iv; v3cp(dp, pc->rc_solv.v) iv += 3; // rc_solv
				}
			}
			else {
				// following iv calculation has to be consistent as above in ::mpi_id == proc_id
				iv += pc->nAtoms * 3; // atoms r
				//iv += pc->nHinges * 3; // hinges, vHinge
				//iv += 3; // dr_cmm
				//iv += 3; // v
				if (bImplicitSolvent) iv += 3; // rc_solv
			}
			if (iv >= NC_MAX) {
				show_log("struct. parameters is more than 2000", true);
				return false;
			}
			MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
			if (::mpi_id != proc_id) {
				iv = 0;
				// following iv calculation has to be consistent as above in ::mpi_id == proc_id
				for (i = 0; i < pc->nAtoms; i++) {
					patom = pc->atom + i;
					dp = dbf + iv; v3cp(patom->r.v, dp) iv += 3;
				}
				/*
				for (i = 0; i < pc->nHinges; i++) {
					dp = dbf + iv; v3cp(pc->hinge[i].vHinge.v, dp) iv += 3;
				}
				*/
				//dp = dbf + iv; v3cp(pc->dr_cmm.v, dp) iv += 3; // dr_cmm
				//dp = dbf + iv; v3cp(pc->bkm->v.v, dp) iv += 3; // v
				if (bImplicitSolvent) {
					dp = dbf + iv; v3cp(pc->rc_solv.v, dp) iv += 3; // rc_solv
				}
			}
		}
		//sprintf(msg, "clusters of mol #%d are exchanged", iMol); show_log(msg, true);
		if (proc_id != 0) { // send Ek to master
			iv = 0;
			if (::mpi_id == proc_id) {
				dp = dbf + iv; dp[0] = mmcell.mm.m[iMol].Ek; iv += 1;  // Ek
				MPI_Send(dbf, iv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
			else if (::mpi_id == 0) { // master
				iv += 1;  // Ek
				MPI_Recv(dbf, iv, MPI_DOUBLE, proc_id, 0, MPI_COMM_WORLD, &status);
				mmcell.mm.m[iMol].Ek = (float)dbf[0];
			}
		}
	}

	if (mmcell.sm.m == NULL) return true;
	CGSM *sm = NULL;
	for (iMol = 0; iMol < mmcell.sm.n; iMol++) {
		proc_id = SM_mpi_process(iMol);
		sm = mmcell.sm.m + iMol;
		// cluster's parameters
		iv = 0;
		if (::mpi_id == proc_id) {
			for (i = 0; i < sm->nAtoms; i++) {
				patom = sm->atom + i;
				dp = dbf + iv; v3cp(dp, patom->r.v) iv += 3;
			}
			//dp = dbf + iv; v3cp(dp, sm->dr_cmm.v) iv += 3; // dr_cmm
			//dp = dbf + iv; v3cp(dp, sm->bkm->v.v) iv += 3; // v
			if (bImplicitSolvent) {
				dp = dbf + iv; v3cp(dp, sm->rc_solv.v) iv += 3; // rc_solv
			}
		}
		else {
			// following iv calculation has to be consistent as above in ::mpi_id == proc_id
			iv += sm->nAtoms * 3; // atoms, r
			//iv += 3; // dr_cmm
			//iv += 3; // v
			iv += 3; // rc_solv
		}
		if (iv >= NC_MAX) return false;
		MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
		if (::mpi_id != proc_id) {
			iv = 0;
			// following iv calculation has to be consistent as above in ::mpi_id == proc_id
			for (i = 0; i < pc->nAtoms; i++) {
				patom = pc->atom + i;
				dp = dbf + iv; v3cp(patom->r.v, dp) iv += 3;
			}
			//dp = dbf + iv; v3cp(sm->dr_cmm.v, dp) iv += 3; // dr_cmm
			//dp = dbf + iv; v3cp(sm->bkm->v.v, dp) iv += 3; // v
			if (bImplicitSolvent) {
				dp = dbf + iv; v3cp(sm->rc_solv.v, dp) iv += 3; // rc_solve
			}
		}
	}

#undef NC_MAX
	return true;
}


bool MPI_GatherMDInfor(LFVERLET_CGMM_MD *lfmdpar, int nMM, LFVERLET_CGSM_MD *lfsmdpar, int nSM) {
	double *dp = NULL;
#define NC_MAX   100
	int iMol = 0, iv = 0, proc_id = -1;
	double dbf[NC_MAX];
	MPI_Status status;

	for (iMol = 0; iMol < nMM; iMol++) {
		proc_id = MM_mpi_process(iMol);
		if (proc_id == 0) continue;
		iv = 2;
		if (::mpi_id == proc_id) { // this process is working on this macro-molecule
			dbf[0] = lfmdpar[iMol].ksi_Hoover;
			dbf[1] = lfmdpar[iMol].vksi_Hoover;
			MPI_Send(dbf, iv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		else if (::mpi_id == 0) {// master
			MPI_Recv(dbf, iv, MPI_DOUBLE, proc_id, 0, MPI_COMM_WORLD, &status);
			lfmdpar[iMol].ksi_Hoover = dbf[0];
			lfmdpar[iMol].vksi_Hoover = dbf[1];
		}
	}
	if (nSM <= 0) return true;
	for (iMol = 0; iMol < nSM; iMol++) {
		proc_id = SM_mpi_process(iMol);
		iv = 2;
		if (proc_id != 0 && ::mpi_id == proc_id) { // this process is working on this macro-molecule
			dbf[0] = lfsmdpar[iMol].ksi_Hoover;
			dbf[1] = lfsmdpar[iMol].vksi_Hoover;
			MPI_Send(dbf, iv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		else if (proc_id != 0 && ::mpi_id == 0) {// master
			MPI_Recv(dbf, iv, MPI_DOUBLE, proc_id, 0, MPI_COMM_WORLD, &status);
			lfsmdpar[iMol].ksi_Hoover = dbf[0];
			lfsmdpar[iMol].vksi_Hoover = dbf[1];
		}
	}
#undef NC_MAX
	return true;
}

void mpi_cluster_check_periodic_cell(CGMM_MD_CELL &mmcell) {
	CG_MMOLECULE *m = NULL;
	CG_CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
		if (nm >= mmcell.mm.n) break;
		m = mmcell.mm.m + nm;
		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			for (i = 0; i < 3; i++) {
				xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
				x = (float)(pc->r->v[i]);
				PERIOD_CHECK(x, xl, xr, period, pc->dr_cmm.v[i])
			}
		}
	}
	/* not necessary, because sm has only one cluster
	SMOLECULE *sm = NULL;
	for (nm = ::n1SM[::mpi_id]; nm <= ::n2SM[::mpi_id]; nm++) {
		if (nm >= mmcell.sm.n) break;
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(sm->r->v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
	}
	*/
}

void mpi_molecule_check_periodic_cell(CGMM_MD_CELL &mmcell) {
	CG_MMOLECULE *m = NULL;
	int nm = 0, i = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
		if (nm >= mmcell.mm.n) break;
		m = mmcell.mm.m + nm;
		m->calMassCenter();
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(m->cm.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) m->shiftMol(ds);
	}
	CGSM *sm = NULL;
	if (mmcell.sm.m != NULL) {
		for (nm = ::n1SM[::mpi_id]; nm <= ::n2SM[::mpi_id]; nm++) {
			if (nm >= mmcell.sm.n) break;
			sm = mmcell.sm.m + nm;
			for (i = 0; i < 3; i++) {
				xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
				x = (float)(sm->r->v[i]);
				PERIOD_CHECK(x, xl, xr, period, ds.v[i])
			}
			if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
		}
	}

	mpi_cluster_check_periodic_cell(mmcell);
}


void multi_mol_distribute_force(CG_MMOLECULE *mm, int nMM, CGSM *sm, int nSM, CMM_CELL3D<CG_CLUSTER> *cmm, CGMM_MD_CELL* mmcell, int nThreads, MM_EF_VAR& var, MM_EF_RES &res) {
	int i = 0, nc = 0, nm = 0;
	double Ep = 0, Ep_t = 0, Etorsion = 0;
	bool Overlap = false;

	if (mm != NULL && nMM > 0) {
		for (nm = 0; nm < nMM; nm++) {
			for (nc = 0; nc < mm[nm].nCluster; nc++) {
				V6zero(mm[nm].cluster[nc].dyn->fc)
				mm[nm].cluster[nc].dyn->overlap = false;
			}
		}
	}
	if (sm != NULL && nSM > 0) {
		for (nm = 0; nm < nSM; nm++) {
			V6zero(sm[nm].dyn->fc)
		}
	}

#if _INTRA_INTER_INTERACTS == 1
	if (var.bCheckClusterInCMM) {
		if (mm != NULL && nMM > 0) {
			for (nm = 0; nm < nMM; nm++) {
				mm[nm].check_non_neighbor_sm(::r2cut_intramol);
			}
		}
	}
#endif
#if _CHECK_TIME_STAMP_ == 1
	sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
	show_infor(errmsg);
#endif

	Ep = 0;
	if (cmm->periodic_cell) {
	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		#error need to realize the multipole calculation of each cg-macro-molecule
		//calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, true, MAX_THREADS);
		// calculate the multipole of each cell
		//calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, false, true, MAX_THREADS);
	#endif

		// check relationship of each cluster, electrostatic and LJ
		if (var.checkCMM) {
		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
			if (mm != NULL && nMM > 0) {
				for (nm = 0; nm < nMM; nm++) {
					// imethod == 0 : LJ interaction for all clusters, ignoring the neighbors
					// imethod == 1 : LJ interaction for the clusters in different neighbors
					// here we use imethod = 0
					CGMMCheckLocalRelshipInCell(*cmm, mm[nm], 0, mm[nm].nCluster, 0, nThreads); // setup relationship
				}
			}
			if (sm != NULL) {
				CGSMCheckLocalRelshipInCell(*cmm, sm, nSM, nThreads); // setup relationship
			}
			#if _CHECK_TIME_STAMP_ == 1
				sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
				show_infor(errmsg);
			#endif
		}

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		Ep = 0;
		if (mm != NULL && nMM > 0) {
			for (nm = 0; nm < nMM; nm++) {
				MMClusterInteraction<CG_MMOLECULE, CG_CLUSTER>((void*)(&FreeCellCMMInteraction_SingleCGMM), *cmm, mm[nm], 0, mm[nm].nCluster, Ep_t, Overlap, false, nThreads);
				Ep += Ep_t;
			}
		}
		if (sm!= NULL && nSM > 0) {
			MultiMolInteraction<CGSM, CG_CLUSTER>((void*)(&FreeCellCMMInteraction_CGSM), *cmm, sm, nSM, Ep_t, Overlap, nThreads);
			Ep += Ep_t;
		}
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "Interaction takes %d ms", time.glance());
		show_infor(errmsg);
	#endif
	}
	else { //case MULTI_CGMM_FREE_CELL_MD:
	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, false, MAX_THREADS);
	#endif
		// check relationship of each cluster, electrostatic and LJ
		if (var.checkCMM) {
		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
			if (mm != NULL && nMM > 0) {
				for (nm = 0; nm < nMM; nm++) {
					// imethod == 0 : LJ interaction for all clusters, ignoring the neighbors
					// imethod == 1 : LJ interaction for the clusters in different neighbors
					// here we use imethod = 0
					CGMMCheckLocalRelshipInCell(*cmm, mm[nm], 0, mm[nm].nCluster, 0, nThreads); // setup relationship
				}
			}
			if (sm!= NULL && nSM > 0) {
				CGSMCheckLocalRelshipInCell(*cmm, sm, nSM, nThreads); // setup relationship
			}
		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
			show_infor(errmsg);
		#endif
		}
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		Ep = 0;
		if (mm != NULL && nMM > 0) {
			for (nm = 0; nm < nMM; nm++) {
				MMClusterInteraction<CG_MMOLECULE, CG_CLUSTER>((void*)(&FreeCellCMMInteraction_SingleCGMM), *cmm, mm[nm], 0, mm[nm].nCluster, Ep_t, Overlap, false, nThreads);
				Ep += Ep_t;
				if (Overlap) res.OVERLAP = true;
			}
		}
		if (sm != NULL && nSM > 0) {
			MultiMolInteraction<CGSM, CG_CLUSTER>((void*)(&FreeCellCMMInteraction_CGSM), *cmm, sm, nSM, Ep_t, Overlap, nThreads);
			Ep += Ep_t;
		}
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "Interaction takes %d ms", time.glance());
		show_infor(errmsg);
	#endif
	}

	#if _CHECK_TIME_STAMP_ == 1
	time.start();
	#endif
	// intra-molecule non-bound interaction
	#if _INTRA_INTER_INTERACTS == 1
	if (mm != NULL && nMM > 0) {
		for (nm = 0; nm < nMM; nm++) {
			MacroMoleculeOperate2<CG_MMOLECULE, double>((void*)(&cgmm_IntraMolInteract), mm[nm], 0, mm[nm].nCluster, nThreads, Ep_t);
			Ep += Ep_t;
		}
	}
	#endif
	// bound-interaction
	if (mm != NULL && nMM > 0) {
		for (nm = 0; nm < nMM; nm++) {
			MacroMoleculeOperate2<CG_MMOLECULE, double>((void*)(&cgbound_interact), mm[nm], 0, mm[nm].nCluster, nThreads, Ep_t);
			Ep += Ep_t;
			if (::bCGdipole && ::Eext != 0) {
				mm[nm].bound_dipole();
				MacroMoleculeOperate2<CG_MMOLECULE, double>((void*)(&calDipoleEffect), mm[nm], 0, mm[nm].nBound, nThreads, Ep_t);
				Ep += Ep_t;
			}
		}
	}
	#if _CHECK_TIME_STAMP_ == 1
	sprintf(::errmsg, "bound-interaction takes %d ms", time.glance());
	show_infor(errmsg);
	#endif

	if (bImplicitSolvent) {
		// neighbors of cluster for implicit solvation 
		if (mm != NULL && nMM > 0) {
			for (nm = 0; nm < nMM; nm++) {
				CGMMClusterImpSolvNeighborsInCell(cmm, mm + nm, 0, mm[nm].nCluster - 1);
				// solvation energy
				MacroMolCellOperate((void*)(&CGMM_ImplicitSolvationForce), *cmm, mm[nm], 0, mm[nm].nCluster - 1, nThreads);
				for (nc = 0; nc < mm[nm].nCluster; nc++) Ep += mm[nm].cluster[nc].E_ImplicitSolvation;
			}
		}
	}

	if (scale_force) {
		if (mm != NULL && nMM > 0) {
			for (nm = 0; nm < nMM; nm++) {
				MacroMoleculeOperate<CG_MMOLECULE>((void*)(&scale_force_cgcluster), mm[nm], 0, mm[nm].nCluster, nThreads);
			}
		}
		if (sm!= NULL && nSM > 0) {
			MOperate<CGSM>((void*)(&scale_force_cgsm), sm, nSM, nThreads);
		}
	}

	res.Ep = Ep;
}

void mpi_multi_mol_md_proc(CGMM_MD_CELL &mmcell, CMM_CELL3D<CG_CLUSTER> *cmm, int md_mode) {
	if (::bCGdipole) {
		sprintf(errmsg, "external E-field : %8.2f kV/mm", ::Eext); show_msg(errmsg);
	}

	VECTOR3 mass_center;
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}
	if (md_mode == MULTI_CGMM_FREE_CELL_MD) {
		set_free_cell_mass_center(mmcell, mass_center);
	}
	else if (md_mode == MULTI_CGMM_PERIDIC_CELL_MD) {
		molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell
		// using Ewald_Sum to calculate interaction
	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		InitEwaldSumPars(cmm->xd[0]);
	#endif
	}

	int i = 0, j = 0, nloop = 0, nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0;

	for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
		if (nm >= mmcell.mm.n) break;
		mmcell.mm.m[nm].setup_kin_dyn();
		mmcell.mm.m[nm].calBrownianPars();
		// bound-dipole
		if (::bCGdipole) mmcell.mm.m[nm].setup_bounds();
	}
	if (mmcell.sm.m != NULL) {
		for (nm = ::n1SM[::mpi_id]; nm <= ::n2SM[::mpi_id]; nm++) {
			if (nm >= mmcell.sm.n) break;
			mmcell.sm.m[nm].setup_kin_dyn();
			mmcell.sm.m[nm].calBrownianPars();
		}
	}

	// here we have to check the periodity of the cluster position
	if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);

	// since molecule mass center position could be shifted, distribute cluster in CMM again 
	bool reset_cell_chain = true;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		reset_cell_chain = (nm == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.mm.m[nm], *cmm, reset_cell_chain); // set base cluster at the center cell
	}
	reset_cell_chain = (mmcell.mm.n == 0 ? true : false);
	cmm_init_distribute_cluster(mmcell.sm.m, mmcell.sm.n, *cmm, reset_cell_chain); // add simple solid molecule into cmm
	cmm_init_subcell<CG_CLUSTER>(*cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
	if (!cmm->setup_basic_cell_matrix()) {
		sprintf(::errmsg, "failure to construct basic-cell matrix for CMM"); show_msg(::errmsg);
		return;
	}
	// *****************  over  ******************
	show_all_molecules();

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;

	bool Overlap = false;
	int nstart = 0, nEachThread = 0, nmol = 0;

	CG_MMOLECULE *mm = NULL;
	CG_CLUSTER *pc = NULL;
	LFVERLET_CGMM_MD *lfmdpar = NULL;
	int nMaxCluster = mmcell.max_ncluster();

	int nMMThreads = 0, nMMEachThread = 0;
	CGMM_MD_THREAD_PAR mm_mdpar[MAX_THREADS];

	if (mmcell.mm.n > 0) {
		nMMEachThread = (::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1) / MAX_THREADS;
		if (nMMEachThread > 0) nMMThreads = MAX_THREADS;
		else {nMMThreads = ::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1; nMMEachThread = 1;}

		nstart = ::n1MM[::mpi_id];
		for (i = 0; i < nMMThreads; i++) {
			if (i == nMMThreads - 1) nmol = ::n2MM[::mpi_id] + 1 - nstart;
			else nmol = nMMEachThread;
			mm_mdpar[i].set(mmcell.mm.m + nstart, mmcell.lfmd.m + nstart, nmol, i);
			mm_mdpar[i].set_buffer_memory(nMaxCluster);
			nstart += nmol;
		}
	}

	CGSM *sm = NULL;
	LFVERLET_CGSM_MD *lfsmdpar = NULL;

	int nSMThreads = 0, nSMEachThread = 0;
	CGSM_MD_THREAD_PAR sm_mdpar[MAX_THREADS];

	if (mmcell.sm.n > 0) {
		nSMEachThread = (::n2SM[::mpi_id] - ::n1SM[::mpi_id] + 1) / MAX_THREADS;
		if (nSMEachThread > 0) nSMThreads = MAX_THREADS;
		else {nSMThreads = ::n2SM[::mpi_id] - ::n1SM[::mpi_id] + 1; nSMEachThread = 1;}

		nstart = ::n1SM[::mpi_id];
		for (i = 0; i < nSMThreads; i++) {
			if (i == nSMThreads - 1) nmol = ::n2SM[::mpi_id] + 1 - nstart;
			else nmol = nSMEachThread;
			sm_mdpar[i].set(mmcell.sm.m + nstart, mmcell.lfsmd.m + nstart, nmol, i);
			nstart += nmol;
		}
	}

	for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
		mm = mmcell.mm.m + nm;
		InitVelocity(*mm, 3 * ::Ek0Internal);
		for (nc = 0; nc < mm->nCluster; nc++) {
			V3zero(mm->cluster[nc].dyn->fc)
			mm->cluster[nc].dyn->overlap = false;
		}
	}

	if (mmcell.sm.n > 0) {
		Ek =  3 * ::Ek0Internal;
		MOperate1<CGSM, double>((void*)(&InitCGSMVelocity), mmcell.sm.m + ::n1SM[::mpi_id], ::n2SM[::mpi_id] - ::n1SM[::mpi_id] + 1, Ek, MAX_THREADS);
	}

	if (!MPI_ExchangeStructInfor(mmcell)) {
		sprintf(errmsg, "failure to exchange the structure"); show_msg(errmsg);
		return;
	}

	mmcell.Init(); //copy the initial torsion angle, speed to the mdpar

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	long LOOPS = 0;
	int save_indx = 0;
	MD_SAVE *mdsave = NULL;
	if (mmcell.mm.n > 0) {LOOPS = mmcell.lfmd.m[0].LOOPS; mdsave = mmcell.lfmd.m[0].mmsave.mdsave;}
	else if (mmcell.sm.n > 0) {LOOPS = mmcell.lfsmd.m[0].LOOPS; mdsave = mmcell.lfsmd.m[0].smsave.mdsave;}

	bool bCheckClusterInCMM = true;

	MM_EF_VAR mef_var;
	MM_EF_RES mef_res;
	MOperateVAR<CG_CLUSTER, CG_MMOLECULE, CGSM, CGSM, CMM_CELL3D<CG_CLUSTER>, CGMM_MD_CELL> mvar;
	mvar.cell = cmm; mvar.mmcell = &mmcell;
	mvar.mm = mmcell.mm.m + ::n1MM[::mpi_id]; mvar.nMM = ::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1;
	if (mmcell.sm.m != NULL) {
		mvar.sm = mmcell.sm.m + ::n1SM[::mpi_id]; mvar.nSM = ::n2SM[::mpi_id] - ::n1SM[::mpi_id] + 1;
	}
	else {mvar.sm = NULL; mvar.nSM = 0;}
	mvar.pm = NULL; mvar.nPM = 0;
	mvar.nMainThreads = MAX_THREADS;
	mvar.nMinorThreads = 1;
	if (mvar.nMinorThreads < 2) mvar.nMinorThreads = 2;

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command
		if (iCheckCMM == 0) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) {
			switch(md_mode) {
			case MULTI_CGMM_PERIDIC_CELL_MD:
				// here we have to check the periodity of the molecule position
				MPI_Barrier(MPI_COMM_WORLD);
				// get the structure calculated in other processors
				if (!MPI_ExchangeStructInfor(mmcell)) {
					sprintf(errmsg, "failure to exchange the structure"); show_msg(errmsg);
					break;
				}
				molecule_check_periodic_cell(mmcell);
				break;
			case MULTI_CGMM_FREE_CELL_MD:
				// get the structure calculated in other processors
				MPI_Barrier(MPI_COMM_WORLD);
				if (!MPI_ExchangeStructInfor(mmcell)) {
					sprintf(errmsg, "failure to exchange the structure"); show_msg(errmsg);
					break;
				}
				set_free_cell_mass_center(mmcell, mass_center);
				break;
			}
		}
		if (::mpi_id == 0) {
			mmcell.MDSave(nloop, (float)Ep);
			show_all_molecules();
		}

		for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
			mm = mmcell.mm.m + nm;
			calInertiaMoment(mm);
		}

		//CellOperate<CG_CLUSTER>((void*)(&InitCGClusterForceInCMM), *cmm, 0, cmm->nCells - 1, MAX_THREADS);

	// re-distribute / check clusters in CMM, use vp of each cluster only
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		if (nloop == 0) cmm_check_cluster<CG_CLUSTER>(*cmm);
		else if (bCheckClusterInCMM) cmm_check<CG_CLUSTER>(*cmm, MAX_THREADS);

		// find out relationship
		// calculate the interaction, and the implicit solvation
		mef_var.checkCMM = (iCheckCMM == 0 || nloop == 0 ? true : false);
		mef_var.bCheckClusterInCMM = bCheckClusterInCMM;
		MOperateInCell<CG_CLUSTER, CG_MMOLECULE, CGSM, CMM_CELL3D<CG_CLUSTER>, CGMM_MD_CELL, MM_EF_VAR, MM_EF_RES>(mvar, 
			(void*)(&multi_mol_distribute_force), mef_var, mef_res);
		Ep = mef_res.Ep;
		Overlap = mef_res.OVERLAP;

		MPI_Barrier(MPI_COMM_WORLD);
		// sum Ep @ master
		Ep_t = Ep;
		MPI_Reduce(&Ep_t, &Ep, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		for (i = 0; i < nMMThreads; i++) mm_mdpar[i].set_loop(nloop);
		for (i = 0; i < nSMThreads; i++) sm_mdpar[i].set_loop(nloop);

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		CGMM_MULTI_THREAD_MD(mm_mdpar, nMMThreads);
		if (mmcell.sm.m != NULL) CGSM_MULTI_THREAD_MD(sm_mdpar, nSMThreads);
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "md takes %d ms", time.glance());
		show_infor(errmsg);
	#endif

		for (i = 0; i < nMMThreads; i++) {
			if (strlen(mm_mdpar[i].msg_show) > 1) show_infor(mm_mdpar[i].msg_show, true);
			if (strlen(mm_mdpar[i].msg_log) > 1) mlog.show(mm_mdpar[i].msg_log, true);
		}
		for (i = 0; i < nSMThreads; i++) {
			if (strlen(sm_mdpar[i].msg_show) > 1) show_infor(sm_mdpar[i].msg_show, true);
			if (strlen(sm_mdpar[i].msg_log) > 1) mlog.show(sm_mdpar[i].msg_log, true);
		}

	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif

		MPI_Barrier(MPI_COMM_WORLD);
		// gather mdpar into master
		MPI_GatherMDInfor(mmcell.lfmd.m, mmcell.lfmd.n, mmcell.lfsmd.m, mmcell.lfsmd.n); 

		//MPI_Barrier(MPI_COMM_WORLD);
		// get the structure calculated in other processors
		//if (!MPI_ExchangeStructInfor(mmcell)) {
		//	sprintf(errmsg, "failure to exchange the structure"); show_msg(errmsg);
		//	break;
		//}
	
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "exchanging structure takes %d ms", time.glance());
		show_infor(errmsg);
	#endif

		Ek = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) Ek += mmcell.mm.m[nm].Ek;
		for (nm = 0; nm < mmcell.sm.n; nm++) Ek += mmcell.sm.m[nm].Ek;

		if (::mpi_id == 0) {
			save_indx = mdsave->mESave.save_indx();
			show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek, (float)Ep);
		}

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		nloop++;
	}
	for (i = 0; i < nMMThreads; i++) mm_mdpar[i].reset();
	for (i = 0; i < nSMThreads; i++) sm_mdpar[i].reset();
}

} // end of namespace _coarse_grain_

#endif // _MPI_

