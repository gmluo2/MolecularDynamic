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
#include "CMM.h"
#include "var.h"

#include "Interaction.h"
#include "Interact1.h"
#include "Interact2.h"
#include "Interact3.h"
#include "NEIMO.h"
#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"
#include "MD_POLYMER.h"
#include "MD_SM.h"
#include "md_cell.h"
#include "mpi_md_cell.h"

#ifdef SEEK_SET
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include "mpi.h"
#endif

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif


void MPI_GatherClusterForce(MMOL_MD_CELL &mmcell) {
	if (mmcell.mm.m == NULL) return;
#define NC_MAX   200
#define V_MAX    2000
	MMOLECULE *mm = NULL;

	int iMol = 0, n1Mol = 0, n2Mol = 0;
	int nc = 0, nv = 0, ncs = 0;
	double dbf[V_MAX];
	CLUSTER *pc = NULL;
	double *fc = NULL;

	int isource = 0, ipkg = 0, npkgs = 0;
	int nSend[3] = {0, 0, 0};
	int total_clusters = 0;
	MPI_Status status;
	int nget = 0;
	char msg[256] = "\0";

	VECTOR<6> fc_net; // net force

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
					dbf[nv+3] = fc[3]; dbf[nv+4] = fc[4]; dbf[nv+5] = fc[5];
					nv += 6; ncs++; nc++;
				}
				MPI_Send(dbf, nv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // send to the master
			}
			MPI_Reduce(mm->mDynKin.fc.v, fc_net.v, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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
					show_infor(msg);
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
							fc[3] = dbf[ncs+3]; fc[4] = dbf[ncs+4]; fc[5] = dbf[ncs+5];
							nc++;
						}
						ncs += 6;
					}
				}
				MPI_Reduce(mm->mDynKin.fc.v, fc_net.v, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				V62V6(fc_net, mm->mDynKin.fc)
			}
		}
	}

	
#undef NC_MAX
#undef V_MAX
}

int MM_mpi_process(int imol) {
	int id = 0;
	for (id = 0; id < ::nMPI; id++) {
		if (imol >= ::n1MM[id] && imol <= ::n2MM[id]) return id;
	}
	char msg[256] = "\0";
	sprintf(msg, "ERROR: macromolecule #%d is not operated by any MPI process", imol);
	show_infor(msg);
	return -1;
}

int SM_mpi_process(int imol) {
	int id = 0;
	for (id = 0; id < ::nMPI; id++) {
		if (imol >= ::n1SM[id] && imol <= ::n2SM[id]) return id;
	}
	char msg[256] = "\0";
	sprintf(msg, "ERROR: macromolecule #%d is not operated by any MPI process", imol);
	show_infor(msg);
	return -1;
}

bool MPI_ExchangeStructInfor(MMOL_MD_CELL &mmcell) {
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	double *dp = NULL;
#define NC_MAX   2000
	int iMol = 0, nc = 0, proc_id = 0;
	int i = 0, iv = 0;
	double dbf[NC_MAX];

	//MPI_Status status;
	char msg[256] = "\0";

	for (iMol = 0; iMol < mmcell.mm.n; iMol++) {
		proc_id = MM_mpi_process(iMol);
		mm = mmcell.mm.m + iMol;
		//sprintf(msg, "mm %d : broadcast from %d", iMol, proc_id); show_log(msg, true);
		// cluster's parameters
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = mm->cluster + nc;
			iv = 0;
			if (::mpi_id == proc_id) {
				for (i = 0; i < pc->nAtoms; i++) {
					patom = pc->atom + i;
					dp = dbf + iv; v3cp(dp, patom->r0.v) iv += 3;
					dp = dbf + iv; v3cp(dp, patom->r.v) iv += 3;
				}
				for (i = 0; i < pc->nHinges; i++) {
					dp = dbf + iv; v3cp(dp, pc->hinge[i].vHinge.v) iv += 3;
				}
				dp = dbf + iv; v3cp(dp, pc->l.v) iv += 3; // l
				dp = dbf + iv; v3cp(dp, pc->up.v) iv += 3; // up
				dp = dbf + iv; v3cp(dp, pc->cm0.v) iv += 3; // cm0
				dp = dbf + iv; v3cp(dp, pc->cm.v) iv += 3; // cm
				dp = dbf + iv; dp[0] = pc->km.ta; dp[1] = pc->km.tp; dp[2] = pc->km.tpp; iv += 3; // ta, tp, tpp
				dp = dbf + iv; v3cp(dp, pc->dr_cmm.v) iv += 3; // dr_cmm
				//dp = dbf + iv; v6cp(dp, pc->bkm->V0.v) iv += 6; // V0
				//dp = dbf + iv; v6cp(dp, pc->bkm->alpha.v) iv += 6; // alpha
				dp = dbf + iv; v3cp(dp, pc->rc_solv.v) iv += 3; // rc_solv
				dp = dbf + iv; v3cp(dp, pc->rc0_solv.v) iv += 3; // rc0_solv
			}
			else {
				// following iv calculation has to be consistent as above in ::mpi_id == proc_id
				iv += pc->nAtoms * 6; // atoms, r, r0
				iv += pc->nHinges * 3; // hinges, vHinge
				iv += 3; // l
				iv += 3; // up
				iv += 3; // cm0
				iv += 3; // cm
				iv += 3; // ta, tp, tpp
				iv += 3; // dr_cmm
				//iv += 6; // V0
				//iv += 6; // alpha
				iv += 3; // rc_solv
				iv += 3; // rc0_solv
			}
			if (iv >= NC_MAX) return false;
			MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
			if (::mpi_id != proc_id) {
				iv = 0;
				// following iv calculation has to be consistent as above in ::mpi_id == proc_id
				for (i = 0; i < pc->nAtoms; i++) {
					patom = pc->atom + i;
					dp = dbf + iv; v3cp(patom->r0.v, dp) iv += 3;
					dp = dbf + iv; v3cp(patom->r.v, dp) iv += 3;
				}
				for (i = 0; i < pc->nHinges; i++) {
					dp = dbf + iv; v3cp(pc->hinge[i].vHinge.v, dp) iv += 3;
				}
				dp = dbf + iv; v3cp(pc->l.v, dp) iv += 3; // l
				dp = dbf + iv; v3cp(pc->up.v, dp) iv += 3; // up
				dp = dbf + iv; v3cp(pc->cm0.v, dp) iv += 3; // cm0
				dp = dbf + iv; v3cp(pc->cm.v, dp) iv += 3; // cm
				dp = dbf + iv; pc->km.ta = dp[0]; pc->km.tp = dp[1]; pc->km.tpp = dp[2]; iv += 3; // ta, tp, tpp
				dp = dbf + iv; v3cp(pc->dr_cmm.v, dp) iv += 3; // dr_cmm
				//dp = dbf + iv; v6cp(pc->bkm->V0.v, dp) iv += 6; // V0
				//dp = dbf + iv; v6cp(pc->bkm->alpha.v, dp) iv += 6; // alpha
				dp = dbf + iv; v3cp(pc->rc_solv.v, dp) iv += 3; // rc_solv
				dp = dbf + iv; v3cp(pc->rc0_solv.v, dp) iv += 3; // rc0_solv
			}
		}
		//sprintf(msg, "clusters are exchanged"); show_log(msg, true);

		// macromolecule parameters
		iv = 0; 
		if (::mpi_id == proc_id) {
			dp = dbf + iv; v3cp(dp, mm->xaxis.v) iv += 3; // xaxis
			dp = dbf + iv; v3cp(dp, mm->zaxis.v) iv += 3; // zaxis
			dp = dbf + iv; v3cp(dp, mm->mDynKin.cm.v) iv += 3; // mass center
			dp = dbf + iv; v3cp(dp, mm->mDynKin.r0.v) iv += 3; // r0
			dp = dbf + iv; v6cp(dp, mm->mDynKin.V.v) iv += 6; // mass center velocity
			dp = dbf + iv; v6cp(dp, mm->mDynKin.P.v) iv += 6; // moment at mass center
		
			dp = dbf + iv; v3cp(dp, mm->r0.v) iv += 3; // base position
			dp = dbf + iv; v6cp(dp, mm->V0.v) iv += 6; // base velocity
			dp = dbf + iv; v6cp(dp, mm->alpha0.v) iv += 6; // base acceleration
			dp = dbf + iv; dp[0] = mm->Ek; iv += 1; //Ek
		}
		else {
			// following iv calculation has to be consistent as above in ::mpi_id == proc_id
			iv += 3; // xaxis 
			iv += 3; // zaxis
			iv += 3; // mass center
			iv += 3; // r0
			iv += 6; // mass center velocity
			iv += 6; // moment at mass center

			iv += 3; // base position
			iv += 6; // base velocity
			iv += 6; // base acceleration
			iv += 1; // Ek
		}
		if (iv >= NC_MAX) return false;
		MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
		iv = 0;
		if (::mpi_id != proc_id) {
			// following iv calculation has to be consistent as above in ::mpi_id == proc_id
			dp = dbf + iv; v3cp(mm->xaxis.v, dp) iv += 3; // xaxis
			dp = dbf + iv; v3cp(mm->zaxis.v, dp) iv += 3; // zaxis
			dp = dbf + iv; v3cp(mm->mDynKin.cm.v, dp) iv += 3; // mass center
			dp = dbf + iv; v3cp(mm->mDynKin.r0.v, dp) iv += 3; // r0
			dp = dbf + iv; v6cp(mm->mDynKin.V.v, dp) iv += 6; // mass center velocity
			dp = dbf + iv; v6cp(mm->mDynKin.P.v, dp) iv += 6; // moment at mass center
		
			dp = dbf + iv; v3cp(mm->r0.v, dp) iv += 3; // base position
			dp = dbf + iv; v6cp(mm->V0.v, dp) iv += 6; // base velocity
			dp = dbf + iv; v6cp(mm->alpha0.v, dp) iv += 6; // base acceleration
			dp = dbf + iv; mm->Ek = (float)(dp[0]); iv += 1; // Ek
		}
	}

	if (mmcell.sm.m == NULL) return true;
	SMOLECULE *sm = NULL;
	for (iMol = 0; iMol < mmcell.sm.n; iMol++) {
		proc_id = SM_mpi_process(iMol);
		sm = mmcell.sm.m + iMol;
		// cluster's parameters
		iv = 0;
		if (::mpi_id == proc_id) {
			for (i = 0; i < sm->c->nAtoms; i++) {
				patom = sm->c->atom + i;
				dp = dbf + iv; v3cp(dp, patom->r0.v) iv += 3;
				dp = dbf + iv; v3cp(dp, patom->r.v) iv += 3;
			}
			dp = dbf + iv; v3cp(dp, sm->r0.v) iv += 3; // mass center
			dp = dbf + iv; v3cp(dp, sm->c->dr_cmm.v) iv += 3; // dr_cmm
			//dp = dbf + iv; v6cp(dp, sm->c->bkm.V0.v) iv += 6; // V0
			//dp = dbf + iv; v6cp(dp, sm->c->bkm.alpha.v) iv += 6; // alpha
			dp = dbf + iv; v3cp(dp, sm->xaxis.v) iv += 3; // xaxis
			dp = dbf + iv; v3cp(dp, sm->zaxis.v) iv += 3; // zaxis
			dp = dbf + iv; v3cp(dp, sm->c->rc0_solv.v) iv += 3; // rc0_solv
			dp = dbf + iv; v3cp(dp, sm->c->rc_solv.v) iv += 3; // rc_solv
			dp = dbf + iv; dp[0] = sm->Ek; dp[1] = sm->Ep; iv += 2; // Ek, Ep
		}
		else {
			// following iv calculation has to be consistent as above in ::mpi_id == proc_id
			iv += sm->c->nAtoms * 6; // atoms, r, r0
			iv += 3; // mass center
			iv += 3; // dr_cmm
			//iv += 6; // V0
			//iv += 6; // alpha
			iv += 3; // xaxis
			iv += 3; // zaxis
			iv += 3; // rc0_solv
			iv += 3; // rc_solv
			iv += 2; // Ek, Ep
		}
		if (iv >= NC_MAX) return false;
		MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
		if (::mpi_id != proc_id) {
			iv = 0;
			// following iv calculation has to be consistent as above in ::mpi_id == proc_id
			for (i = 0; i < pc->nAtoms; i++) {
				patom = pc->atom + i;
				dp = dbf + iv; v3cp(patom->r0.v, dp) iv += 3;
				dp = dbf + iv; v3cp(patom->r.v, dp) iv += 3;
			}
			dp = dbf + iv; v3cp(sm->r0.v, dp) iv += 3; // mass center
			dp = dbf + iv; v3cp(sm->c->dr_cmm.v, dp) iv += 3; // dr_cmm
			//dp = dbf + iv; v6cp(sm->c->bkm.V0.v, dp) iv += 6; // V0
			//dp = dbf + iv; v6cp(sm->c->bkm.alpha.v, dp) iv += 6; // alpha
			dp = dbf + iv; v3cp(sm->xaxis.v, dp) iv += 3; // xaxis
			dp = dbf + iv; v3cp(sm->zaxis.v, dp) iv += 3; // zaxis
			dp = dbf + iv; v3cp(sm->c->rc0_solv.v, dp) iv += 3; // rc0_solve
			dp = dbf + iv; v3cp(sm->c->rc_solv.v, dp) iv += 3; // rc_solve
			dp = dbf + iv; sm->Ek = (float)dp[0]; sm->Ep = (float)dp[1]; iv += 2; // Ek, Ep
		}
	}

#undef NC_MAX
	return true;
}

bool MPI_ExchangeStruct(MMOL_MD_CELL &mmcell) {
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	double *dp = NULL;
#define NC_MAX   2000
	int iMol = 0, nc = 0, proc_id = 0;
	int i = 0, iv = 0;
	double dbf[NC_MAX];
	VECTOR3 r2, zaxis2, xaxis2;

	//MPI_Status status;
	char msg[256] = "\0";

	ARRAY< MACROMOL_RUN_BKGD<MMOLECULE> > mRunBkgd;
	mRunBkgd.SetArray(mmcell.mm.n);
	for (i = 0; i < mRunBkgd.n; i++) {mRunBkgd.m[i].func = (void*)(&MacroMolMove); mRunBkgd.m[i].mm = mmcell.mm.m + i;}

	for (iMol = 0; iMol < mmcell.mm.n; iMol++) {
		proc_id = MM_mpi_process(iMol);
		mm = mmcell.mm.m + iMol;
		//sprintf(msg, "mm %d : broadcast from %d", iMol, proc_id); show_log(msg, true);
		// cluster's parameters
		dp = dbf; iv = 0; v3cp(dp, mm->base->Op->v) iv += 3;
		dp = dbf + iv; v3cp(dp, mm->zaxis.v) iv += 3;
		dp = dbf + iv; v3cp(dp, mm->xaxis.v) iv += 3;
		MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
		if (::mpi_id != proc_id) {
			// distribute the parameter to 
			dp = dbf; iv = 0; v3cp(r2.v, dp) iv += 3; 
			dp = dbf + iv; v3cp(zaxis2.v, dp) iv += 3;
			dp = dbf + iv; v3cp(xaxis2.v, dp) iv += 3;
			mm->set_base_move2(r2, zaxis2, xaxis2);
		}
		nc = 0;
		while (1) {
			dp = dbf; iv = 0;
			while (iv < NC_MAX) {
				if (::mpi_id == proc_id) {
					dp = dbf + iv; dp[0] = mm->cluster[nc].km.ta; iv += 1;
				}
				else iv += 1;
				nc += 1;
				if (nc >= mm->nCluster) break;
			}
			//if (::mpi_id == proc_id) {sprintf(msg, "distributing %d cluster of MM#%d from %s", iv, iMol, ::proc_name[proc_id]); show_log(msg, true);}
			MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
			//if (::mpi_id != proc_id) {sprintf(msg, "get %d cluster of MM#%d from %s", iv, iMol, ::proc_name[proc_id]); show_log(msg, true);}
			if (::mpi_id != proc_id) {
				// following iv calculation has to be consistent as above in ::mpi_id == proc_id
				for (i = 0; i < iv; i++) {
					mm->cluster[nc - iv + i].rot.dta = dbf[i] - mm->cluster[nc - iv + i].km.ta;
				}
			}
			if (nc >= mm->nCluster) break;
		}
		if (::mpi_id != proc_id) MacroMol_Run_bkgd<MMOLECULE>(mRunBkgd.m[iMol]);
	}
	for (i = 0; i < mRunBkgd.n; i++) mRunBkgd.m[i].WaitUntilThreadOver();
	mRunBkgd.release();

	if (mmcell.sm.m == NULL) return true;

	ARRAY< MACROMOL_RUN_BKGD<SMOLECULE> > smRunBkgd;
	smRunBkgd.SetArray(mmcell.sm.n);
	for (i = 0; i < smRunBkgd.n; i++) {smRunBkgd.m[i].func = (void*)(&SMolMove); smRunBkgd.m[i].mm = mmcell.sm.m + i;}

	SMOLECULE *sm = NULL;
	for (iMol = 0; iMol < mmcell.sm.n; iMol++) {
		proc_id = SM_mpi_process(iMol);
		sm = mmcell.sm.m + iMol;
		// cluster's parameters
		dp = dbf; iv = 0; v3cp(dp, sm->r0.v) iv += 3;
		dp = dbf + iv; v3cp(dp, sm->zaxis.v) iv += 3;
		dp = dbf + iv; v3cp(dp, sm->xaxis.v) iv += 3;
		MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
		if (::mpi_id != proc_id) {
			// distribute the parameter to 
			dp = dbf; iv = 0; v3cp(r2.v, dp) iv += 3; 
			dp = dbf + iv; v3cp(zaxis2.v, dp) iv += 3;
			dp = dbf + iv; v3cp(xaxis2.v, dp) iv += 3;
			sm->set_base_move2(r2, zaxis2, xaxis2);
			MacroMol_Run_bkgd<SMOLECULE>(smRunBkgd.m[iMol]);
		}
	}
	for (i = 0; i < smRunBkgd.n; i++) smRunBkgd.m[i].WaitUntilThreadOver();
	smRunBkgd.release();

#undef NC_MAX
	return true;
}

void MPI_ExchangeClusterImpSolvRshipPars(MMOL_MD_CELL &mmcell) {
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	double *dp = NULL;
#define NC_MAX   2000
	int iMol = 0, nc = 0, proc_id = 0, iproc = 0;
	int i = 0, iv = 0;
	double dbf[NC_MAX];

	char msg[256] = "\0";

	if (mmcell.sm.n > 0) {
		for (iproc = 0; iproc < ::nMPI; iproc++) {
			iMol = ::n1SM[iproc];
			while (1) {
				iv = 0;
				while (iv < NC_MAX) {
					if (::mpi_id == iproc) dbf[iv] = mmcell.sm.m[iMol].c->SASA;
					iv++; iMol++;
					if (iMol > ::n2SM[iproc]) break;
				}
				MPI_Bcast(dbf, iv, MPI_DOUBLE, iproc, MPI_COMM_WORLD);
				if (::mpi_id != iproc) {
					for (i = 0; i < iv; i++) mmcell.sm.m[iMol - iv + i].c->SASA = dbf[i];
				}
				if (iMol > ::n2SM[iproc]) break;
			}
		}
	}

	//show_log("exchange implicit solvation parameters of MM: ", false);
	for (iMol = 0; iMol < mmcell.mm.n; iMol++) {
		proc_id = MM_mpi_process(iMol);
		//sprintf(msg, "%d ", iMol); show_log(msg, false);
		mm = mmcell.mm.m + iMol;
		// cluster's parameters
		nc = 0;
		while (1) {
			dp = dbf; iv = 0;
			while (iv < NC_MAX) {
				if (::mpi_id == proc_id) {
					dp = dbf + iv; dp[0] = mm->cluster[nc].SASA; iv += 1;
				}
				else iv += 1;
				nc += 1;
				if (nc >= mm->nCluster) break;
			}
			MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
			if (::mpi_id != proc_id) {
				// following iv calculation has to be consistent as above in ::mpi_id == proc_id
				for (i = 0; i < iv; i++) {
					mm->cluster[nc - iv + i].SASA = dbf[i];
				}
			}
			if (nc >= mm->nCluster) break;
		}
	}
	//show_log("", true);
#undef NC_MAX
}

void MPI_ExchangeKineticEnergy(MMOL_MD_CELL &mmcell) {
	MMOLECULE *mm = NULL;
	double *dp = NULL;
#define NC_MAX   2000
	int iproc = 0, iMol = 0;
	int i = 0, iv = 0;
	double dbf[NC_MAX];

	char msg[256] = "\0";

	if (mmcell.sm.n > 0) {
		for (iproc = 0; iproc < ::nMPI; iproc++) {
			iMol = ::n1SM[iproc];
			while (1) {
				iv = 0;
				while (iv < NC_MAX) {
					if (::mpi_id == iproc) dbf[iv] = mmcell.sm.m[iMol].Ek;
					iv++; iMol++;
					if (iMol > ::n2SM[iproc]) break;
				}
				MPI_Bcast(dbf, iv, MPI_DOUBLE, iproc, MPI_COMM_WORLD);
				if (::mpi_id != iproc) {
					for (i = 0; i < iv; i++) mmcell.sm.m[iMol - iv + i].Ek = dbf[i];
				}
				if (iMol > ::n2SM[iproc]) break;
			}
		}
	}
	for (iproc = 0; iproc < ::nMPI; iproc++) {
		iMol = ::n1MM[iproc];
		while (1) {
			iv = 0;
			while (iv < NC_MAX) {
				if (::mpi_id == iproc) dbf[iv] = mmcell.mm.m[iMol].Ek;
				iv++; iMol++;
				if (iMol > ::n2MM[iproc]) break;
			}
			MPI_Bcast(dbf, iv, MPI_DOUBLE, iproc, MPI_COMM_WORLD);
			if (::mpi_id != iproc) {
				for (i = 0; i < iv; i++) mmcell.mm.m[iMol - iv + i].Ek = dbf[i];
			}
			if (iMol > ::n2MM[iproc]) break;
		}
	}
#undef NC_MAX
}

bool MPI_ExchangeClusterMultipoleInfor(int proc_id, ELECTRO_MULTIPOLE *emp, bool bEwaldSum) {
	double *dp = NULL;
#define NC_MAX   2000
	int i, j, k, n, iv = 0;
	double dbf[NC_MAX];

	iv = 0;
	if (::mpi_id == proc_id) {
		dp = dbf + iv; dp[0] = emp->q; dp[1] = emp->d; iv += 2; // q & d
		dp = dbf + iv; v3cp(dp, emp->r0.v) iv += 3; // charge center
		dp = dbf + iv; v3cp(dp, emp->mu.v) iv += 3; // mu
		dp = dbf + iv; dmatrix2array(dp, emp->Q, i, j, n) iv += n; // Q
	}
	else {
		// following iv calculation has to be consistent as above in ::mpi_id == proc_id
		iv += 2; // q & d
		iv += 3; // charge center
		iv += 3; // mu
		iv += 3 * 3; // Q -- 3x3 matrix
	}
	MPI_Bcast(dbf, iv, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
	if (::mpi_id != proc_id) {
		iv = 0;
		// following iv calculation has to be consistent as above in ::mpi_id == proc_id
		dp = dbf + iv; emp->q = (float)dp[0]; emp->d = (float)dp[1]; iv += 2; // q & d
		dp = dbf + iv; v3cp(emp->r0.v, dp) iv += 3; // charge center
		dp = dbf + iv; v3cp(emp->mu.v, dp) iv += 3; // mu
		dp = dbf + iv; array2dmatrix(emp->Q, dp, i, j, n) iv += n; // Q
	}
	if (bEwaldSum) {
		// fsin
		if (::mpi_id == proc_id) {
			cmatrix2array(dbf, emp->fsin, i, j, k, n)
		}
		else n = NH * NH * NH; 
		if (n >= NC_MAX) return false;
		MPI_Bcast(dbf, n, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
		if (::mpi_id != proc_id) {
			array2cmatrix(emp->fsin, dbf, i, j, k, n)
		}

		// fcos
		if (::mpi_id == proc_id) {
			cmatrix2array(dbf, emp->fcos, i, j, k, n)
		}
		else n = NH * NH * NH; 
		if (n >= NC_MAX) return false;
		MPI_Bcast(dbf, n, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
		if (::mpi_id != proc_id) {
			array2cmatrix(emp->fcos, dbf, i, j, k, n)
		}

		// muh
		if (::mpi_id == proc_id) {
			cmatrix2array(dbf, emp->muh, i, j, k, n)
		}
		else n = NH * NH * NH; 
		MPI_Bcast(dbf, n, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
		if (::mpi_id != proc_id) {
			array2cmatrix(emp->muh, dbf, i, j, k, n)
		}

		// Qh
		if (::mpi_id == proc_id) {
			cmatrix2array(dbf, emp->Qh, i, j, k, n)
		}
		else n = NH * NH * NH; 
		MPI_Bcast(dbf, n, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
		if (::mpi_id != proc_id) {
			array2cmatrix(emp->Qh, dbf, i, j, k, n)
		}
	}
#undef NC_MAX
	return true;
}

bool MPI_ExchangeClusterMultipoleInfor(MMOL_MD_CELL &mmcell, bool bEwaldSum) {
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	ELECTRO_MULTIPOLE *emp = NULL;

	int iMol, nc;
	int proc_id = -1;

	//MPI_Status status;
	char msg[256] = "\0";

	for (iMol = 0; iMol < mmcell.mm.n; iMol++) {
		proc_id = MM_mpi_process(iMol);
		mm = mmcell.mm.m + iMol;
		// cluster's parameters
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = mm->cluster + nc;
			if (pc->emp == NULL) continue;
			emp = pc->emp;
			if (!MPI_ExchangeClusterMultipoleInfor(proc_id, emp, bEwaldSum)) return false;
		}
	}

	if (mmcell.sm.n == 0) return true;
	SMOLECULE *sm = NULL;
	for (iMol = 0; iMol < mmcell.sm.n; iMol++) {
		proc_id = SM_mpi_process(iMol);
		sm = mmcell.sm.m + iMol;
		if (sm->c->emp == NULL) continue;
		emp = sm->c->emp;
		// cluster's parameters
		if (!MPI_ExchangeClusterMultipoleInfor(proc_id, emp, bEwaldSum)) return false;
	}
	return true;
}

bool MPI_GatherMDInfor(LFVERLET_MM_MD *lfmdpar, int nMM, LFVERLET_SM_MD *lfsmdpar, int nSM) {
	double *dp = NULL;
#define NC_MAX   100
	int iMol = 0, iv = 0, proc_id = -1;
	double dbf[NC_MAX];
	MPI_Status status;

	for (iMol = 0; iMol < nMM; iMol++) {
		proc_id = MM_mpi_process(iMol);
		iv = 8;
		if (proc_id != 0 && ::mpi_id == proc_id) { // this process is working on this macro-molecule
			dbf[0] = lfmdpar[iMol].ksi_frame_rot;
			dbf[1] = lfmdpar[iMol].ksi_frame_trans;
			dbf[2] = lfmdpar[iMol].ksi_internal;
			dbf[3] = lfmdpar[iMol].ksi_Hoover;
			dbf[4] = lfmdpar[iMol].vksi_frame_rot;
			dbf[5] = lfmdpar[iMol].vksi_frame_trans;
			dbf[6] = lfmdpar[iMol].vksi_internal;
			dbf[7] = lfmdpar[iMol].vksi_Hoover;
			MPI_Send(dbf, iv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		else if (proc_id != 0 && ::mpi_id == 0) {// master
			MPI_Recv(dbf, iv, MPI_DOUBLE, proc_id, 0, MPI_COMM_WORLD, &status);
			lfmdpar[iMol].ksi_frame_rot = dbf[0];
			lfmdpar[iMol].ksi_frame_trans = dbf[1];
			lfmdpar[iMol].ksi_internal = dbf[2];
			lfmdpar[iMol].ksi_Hoover = dbf[3];
			lfmdpar[iMol].vksi_frame_rot = dbf[4];
			lfmdpar[iMol].vksi_frame_trans = dbf[5];
			lfmdpar[iMol].vksi_internal = dbf[6];
			lfmdpar[iMol].vksi_Hoover = dbf[7];
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

void mpi_cluster_check_periodic_cell(MMOL_MD_CELL &mmcell) {
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
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
				x = (float)(pc->cm.v[i]);
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
			x = (float)(sm->r0.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
	}
	*/
}

void mpi_molecule_check_periodic_cell(MMOL_MD_CELL &mmcell) {
	MMOLECULE *m = NULL;
	int nm = 0, i = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
		if (nm >= mmcell.mm.n) break;
		m = mmcell.mm.m + nm;
		m->calMassCenter(false);
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(m->mDynKin.cm.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) m->shiftMol(ds);
	}
	SMOLECULE *sm = NULL;
	if (mmcell.sm.m != NULL) {
		for (nm = ::n1SM[::mpi_id]; nm <= ::n2SM[::mpi_id]; nm++) {
			if (nm >= mmcell.sm.n) break;
			sm = mmcell.sm.m + nm;
			for (i = 0; i < 3; i++) {
				xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
				x = (float)(sm->r0.v[i]);
				PERIOD_CHECK(x, xl, xr, period, ds.v[i])
			}
			if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
		}
	}

	mpi_cluster_check_periodic_cell(mmcell);
}

void mpi_calClusterMultipole(MMOL_MD_CELL &mmcell, bool bEwaldSum) {
	int iMol = 0;
	MMOLECULE *mm = NULL;
	for (iMol = ::n1MM[::mpi_id]; iMol <= ::n1MM[::mpi_id]; iMol++) {
		if (iMol >= mmcell.mm.n) break;
		mm = mmcell.mm.m + iMol;
		MacroMoleculeOperate1<MMOLECULE, bool>((void*)(&CalClusterElectroMultipole), *mm, 0, mm->nCluster - 1, bEwaldSum, MAX_THREADS);
	}
	if (mmcell.sm.n > 0) MOperate1<SMOLECULE, bool>((void*)(&CalSMElectroMultipole), mmcell.sm.m + ::n1SM[::mpi_id], ::n2SM[::mpi_id] - ::n1SM[::mpi_id] + 1, bEwaldSum, MAX_THREADS);
}

/*
#if _SYS_ == _WINDOWS_SYS_
int SM_LFVERLET_MD_thread_func(LPVOID *par) {
#elif _SYS_ == _LINUX_SYS_
void* SM_LFVERLET_MD_thread_func(void *par) {
#endif
	SM_MD_THREAD_PAR *sm_mdpar = (SM_MD_THREAD_PAR*)par;
	LFVERLET_SM_MD *lfsmdpar = NULL;
	SMOLECULE *sm = NULL;
	int nm = 0;

	strcpy(sm_mdpar->msg_show, "\0"); strcpy(sm_mdpar->msg_log, "\0");

	for (nm = 0; nm < sm_mdpar->n; nm++) {
		sm = sm_mdpar->sm + nm;
		sm->calInertiaMoment();

		lfsmdpar = sm_mdpar->lfsmdpar + nm;
		sm_mdpar->res[nm] = SM_LeapFrog_velocity_acceleration(sm_mdpar->bHooverDyn[nm], 
			*sm, *lfsmdpar, sm_mdpar->md_loop, sm_mdpar->nconvg_loops[nm], sm_mdpar->msg_show, sm_mdpar->msg_log, 5120);
	}
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(sm_mdpar->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void SM_MULTI_THREAD_MD(SM_MD_THREAD_PAR* sm_mdpar, int nThreads) { // we assume nThread > 1
	int i = 0;
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nThreads; i++) {
		ResetEvent(sm_mdpar[i].hMMThread);
		sm_mdpar[i].thread = AfxBeginThread((AFX_THREADPROC)SM_LFVERLET_MD_thread_func, (LPVOID)(sm_mdpar + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nThreads; i++) {
		WaitForSingleObject(sm_mdpar[i].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nThreads; i++) {
		pthread_create(&(sm_mdpar[i].thread), NULL, &SM_LFVERLET_MD_thread_func, (void *)(sm_mdpar + i));
	}
	for (i = 0; i < nThreads; i++) {
		if (sm_mdpar[i].thread != 0) pthread_join(sm_mdpar[i].thread, NULL); // wait the thread over
	}
#endif
}

#if _SYS_ == _WINDOWS_SYS_
int MM_LFVERLET_MD_thread_func(LPVOID *par) {
#elif _SYS_ == _LINUX_SYS_
void* MM_LFVERLET_MD_thread_func(void *par) {
#endif
	MM_MD_THREAD_PAR *mm_mdpar = (MM_MD_THREAD_PAR*)par;
	LFVERLET_MM_MD *lfmdpar = NULL;
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, nc = 0, j = 0;

	VECTOR<6> Iw0;
	
	strcpy(mm_mdpar->msg_show, "\0"); strcpy(mm_mdpar->msg_log, "\0");
	
	for (nm = 0; nm < mm_mdpar->n; nm++) {
		mm = mm_mdpar->mm + nm;
		if (bRandomForce) gen_add_random_force1(*mm, 0, mm->nCluster, randf, 0.35);
		CalNetForce(*mm, 0, mm->nCluster); // the net force on molecular mass center

		lfmdpar = mm_mdpar->lfmdpar + nm;
		mm_mdpar->res[nm] = LeapFrog_SelfConsistent_velocity_acceleration(mm_mdpar->bHooverDyn[nm], 
			*mm, *lfmdpar, mm_mdpar->md_loop, mm_mdpar->nconvg_loops[nm],
			mm_mdpar->vect1, mm_mdpar->vect2, mm_mdpar->dv1, mm_mdpar->dv2, mm_mdpar->local_relax, mm_mdpar->msg_show, mm_mdpar->msg_log, 5120); // translation and rotation of whole molecule are 1.5kT
	}
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(mm_mdpar->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void MM_MULTI_THREAD_MD(MM_MD_THREAD_PAR* mm_mdpar, int nThreads) { // we assume nThread > 1
	int i = 0;
#if _SYS_ == _WINDOWS_SYS_ 
	for (i = 0; i < nThreads; i++) {
		ResetEvent(mm_mdpar[i].hMMThread);
		mm_mdpar[i].thread = AfxBeginThread((AFX_THREADPROC)MM_LFVERLET_MD_thread_func, (LPVOID)(mm_mdpar + i), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (i = 0; i < nThreads; i++) {
		WaitForSingleObject(mm_mdpar[i].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (i = 0; i < nThreads; i++) {
		pthread_create(&(mm_mdpar[i].thread), NULL, &MM_LFVERLET_MD_thread_func, (void *)(mm_mdpar + i));
	}
	for (i = 0; i < nThreads; i++) {
		if (mm_mdpar[i].thread != 0) pthread_join(mm_mdpar[i].thread, NULL); // wait the thread over
	}
#endif
}
*/

extern void molecule_check_periodic_cell(MMOL_MD_CELL &mmcell);
extern void SM_MULTI_THREAD_MD(SM_MD_THREAD_PAR* sm_mdpar, int nThreads);
extern void MM_MULTI_THREAD_MD(MM_MD_THREAD_PAR* mm_mdpar, int nThreads);

void multi_mol_distribute_force(MMOLECULE *mm, int nMM, SMOLECULE *sm, int nSM, CMM_CELL3D<BASIC_CLUSTER> *cmm, MMOL_MD_CELL *mmcell, int nThreads, MM_EF_VAR& var, MM_EF_RES &res) {
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
			V6zero(sm[nm].c->dyn->fc)
		}
	}

#define _interact 0
#if _interact == 0 // CMM + EwaldSum
	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		// calculate the electro-multipole of each cluster in the cell, EwaldSum
		//CellOperate1<BASIC_CLUSTER, bool>((void*)(&calClusterElectroMultipole), *cmm, 0, cmm->nCells, true, MAX_THREADS);
		
		// calculate the cluster in this processor
		// and exchange the multipole information of clusters in other processors
		mpi_calClusterMultipole(*mmcell, true);
		MPI_Barrier(MPI_COMM_WORLD);
		if (!MPI_ExchangeClusterMultipoleInfor(*mmcell, true)) {
			sprintf(errmsg, "failure to exchange the cluster's multipole"); show_infor(errmsg);
			break;
		}

		// calculate the multipole of each cell
		calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, false, true, MAX_THREADS);
	#endif

	// cluster's relationship checking, electrostatic and LJ
	if (var.checkCMM) {
		// setup relationship for each cluster in this process
		if (mm != NULL && nMM > 0) {
			for (nm = 0; nm < nMM; nm++) {
				MMCheckRelshipInCell(mm[nm], 0, mm[nm].nCluster - 1, *cmm, nThreads);
			}
		}
		if (sm != NULL && nSM > 0) {
		//SMRelshipInCell(sm, nSM, cmm);
			SMCellOperate<SMOLECULE, BASIC_CLUSTER>((void*)(&SMRelshipInCell), *cmm, sm, nSM, nThreads);
		}
	}

#if	_ATOM_INDUCED_DIPOLE == 1
	Efield_Polarization(*cmm, MAX_THREADS, (nloop > 1 ? 3 : 1));
#endif

#if _ATOM_INDUCED_DIPOLE == 0
	if (mm != NULL && nMM > 0) {
		for (nm = 0; nm < nMM; nm++) {
			MMClusterInteraction<MMOLECULE, BASIC_CLUSTER>((void*)(&EwaldSumCMMInteraction_SingleMM), *cmm, mm[nm], 0, mm[nm].nCluster, mm[nm].Ep, Overlap, false, nThreads); // calculate the interaction with CMM + Ewald_Sum 
			Ep += mm[nm].Ep;
		}
	}
	if (sm != NULL && nSM > 0) {
		Ep_t = 0;
		MultiMolInteraction<SMOLECULE, BASIC_CLUSTER>((void*)(&EwaldSumCMMInteraction_SM), *cmm, sm, nSM, Ep_t, Overlap, nThreads);
		if (Overlap) res.OVERLAP = true;
		Ep += Ep_t;
	}

	#elif _ATOM_INDUCED_DIPOLE == 1
		#error interaction @ ATOM INDUCED DIPOLE, is not realized yet
		ClusterInteraction<BASIC_CLUSTER>((void*)(&Polar_CMMEwaldSumInteraction), *cmm, 0, cmm->nCells, Ep, Overlap, false, MAX_THREADS); // calculate the interaction with CMM + Ewald_Sum 
	#endif
#elif _interact == 1 // CMM, direct electrostatic, non-EwaldSum, check the function in Interact1.cpp
	#error direct interaction @ CMM, is not realized yet
	/**** direct calculation, non-Ewald Sum ****/
	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	// calculate the electro-multipole of each cluster in the cell, EwaldSum
	CellOperate1<BASIC_CLUSTER, bool>((void*)(&calClusterElectroMultipole), *cmm, 0, cmm->nCells, false, MAX_THREADS);
	#endif

	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		// calculate the multipole of each cell
		calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, false, false, MAX_THREADS);
	#endif

	if (iCheckCMM == 0) {
		CheckRelship_ImagCell(*cmm, 0, cmm->nCells, false, MAX_THREADS); // setup relationship for each atom with CMM
	}

	ClusterInteraction<BASIC_CLUSTER>((void*)(&CMMInteraction), *cmm, 0, cmm->nCells, Ep, Overlap, false, MAX_THREADS); // calculate the interaction with CMM
	/*************OVER***************/
#endif
#undef _interact

	if (mm != NULL && nMM > 0) {
		for (nm = 0; nm < nMM; nm++) {
			//mm[nm].calTorsion(Etorsion);
			MacroMoleculeOperate2<MMOLECULE, double>((void*)(&calTorsion), mm[nm], 0, mm[nm].nCluster, nThreads, Etorsion);
			Ep += Etorsion;
		}
	}

	/*
	if (scale_force) {
		if (mm != NULL && nMM > 0) {
			for (nm = 0; nm < nMM; nm++) {
				MacroMoleculeOperate<MMOLECULE>((void*)(&mm_scale_force), mm[nm], 0, mm[nm].nCluster, nThreads);
			}
		}
		// we do not care the small molecule here
	}
	*/

	if (bImplicitSolvent) {
		// implicit solvation energy / force, of the molecule in this process
		if (mm != NULL && nMM > 0) {
			// implicit solvation energy / force -- dependent on the cluster itself
			/*
			for (nm = 0; nm < nMM; nm++) {
				// implicit solvation relationship, neighbors of cluster for implicit solvation 
				if (var.bCheckClusterInCMM) MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MMClusterImpSolvNeighborsInCell), *cmm, mm[nm], 0, mm[nm].nCluster - 1, nThreads);
				MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MM_ImplicitSolvationForce), *cmm, mm[nm], 0, mm[nm].nCluster, nThreads);
				for (nc = 0; nc < mm[nm].nCluster; nc++) Ep += mm[nm].cluster[nc].E_ImplicitSolvation;
			}
			*/

			// implicit solvation energy / force -- force = d(U1+U2)/dr, dependent on the related two clusters
			/*
			for (nm = 0; nm < nMM; nm++) {
				// implicit solvation relationship, neighbors of cluster for implicit solvation 
				if (var.bCheckClusterInCMM) MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MMClusterImpSolvNeighborsInCell), *cmm, mm[nm], 0, mm[nm].nCluster - 1, nThreads);
				MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MM_ImpSolvSASA), *cmm, mm[nm], 0, mm[nm].nCluster, nThreads);
				for (nc = 0; nc < mm[nm].nCluster; nc++) Ep += mm[nm].cluster[nc].E_ImplicitSolvation;
			}
			//in this case, exchange of solvation parameters has to be done out of this parallel calculation
			//MPI_ExchangeClusterImpSolvRshipPars(*mmcell); // have to exchange it so that we can get the parameters for related cluster
			//for (nm = 0; nm < nMM; nm++) {
			//	MacroMoleculeOperate<MMOLECULE>((void*)(&MM_ImpSolvForce), mm[nm], 0, mm[nm].nCluster, nThreads);
			//}
			*/
		}
	}

	res.Ep = Ep;
}

void multi_mol_impsolv_initialize(MMOLECULE *mm, int nMM, SMOLECULE *sm, int nSM, CMM_CELL3D<BASIC_CLUSTER> *cmm, MMOL_MD_CELL *mmcell, int nThreads, MM_IMPSOLV_VAR& var, MM_IMPSOLV_RES &res) {
	int nm = 0;
	// check the implicit solvation relationship and their SASA calculations
	// implicit solvation energy / force -- force = d(U1+U2)/dr, dependent on the related two clusters
	// release the implicit-solv relationships and the related parameters, for all the molecules
	if (mm != NULL && nMM > 0) {
		for (nm = 0; nm < nMM; nm++) {
			// implicit solvation relationship, neighbors of cluster for implicit solvation 
			if (var.bCheckClusterInCMM) MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MMClusterImpSolvNeighborsInCell), *cmm, mm[nm], 0, mm[nm].nCluster - 1, nThreads);
			MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MM_ImpSolvSASA), *cmm, mm[nm], 0, mm[nm].nCluster, nThreads);
		}
	}
	res.res = 1;
}

void impsolv_initilize(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, bool bCheckClusterInCMM) {
	int i;

	// check the implicit-solv relationships for the molecules in this process, and cal SASA
	MOperateVAR<CLUSTER, MMOLECULE, SMOLECULE, PMOLECULE, CMM_CELL3D<BASIC_CLUSTER>, MMOL_MD_CELL> mvar;
	mvar.cell = cmm; mvar.mmcell = &mmcell;
	mvar.mm = mmcell.mm.m + ::n1MM[::mpi_id]; mvar.nMM = ::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1;
	if (mvar.sm == NULL) {mvar.sm = NULL; mvar.nSM = 0;}
	else {
		mvar.sm = mmcell.sm.m + ::n1SM[::mpi_id]; mvar.nSM = ::n2SM[::mpi_id] - ::n1SM[::mpi_id] + 1;
	}
	mvar.pm = NULL; mvar.nPM = 0;
	mvar.nMainThreads = MAX_THREADS;
	mvar.nMinorThreads = MAX_THREADS * 2 / MAX_THREADS;
	if (mvar.nMinorThreads < 2) mvar.nMinorThreads = 2;

	MM_IMPSOLV_VAR var; var.bCheckClusterInCMM = bCheckClusterInCMM;
	MM_IMPSOLV_RES res;
	MOperateInCell<CLUSTER, MMOLECULE, SMOLECULE, CMM_CELL3D<BASIC_CLUSTER>, MMOL_MD_CELL, MM_IMPSOLV_VAR, MM_IMPSOLV_RES>(mvar, 
		(void*)(&multi_mol_impsolv_initialize), var, res);

	MPI_ExchangeClusterImpSolvRshipPars(mmcell); // exchange the SASA of all clusters
}

void mpi_multi_mol_md_proc(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, int md_mode) {
	int i = 0, j = 0, nloop = 0, nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0;

	for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
		if (nm >= mmcell.mm.n) break;
		mmcell.mm.m[nm].setup_kin_dyn();
	}
	if (mmcell.sm.m != NULL) {
		for (nm = ::n1SM[::mpi_id]; nm <= ::n2SM[::mpi_id]; nm++) {
			if (nm >= mmcell.sm.n) break;
			mmcell.sm.m[nm].setup_kin_dyn();
		}
	}
	if (mmcell.pm.m != NULL) {
		for (nm = 0; nm <= mmcell.pm.n; nm++) {
			if (nm >= mmcell.pm.n) break;
			mmcell.pm.m[nm].setup_dyn();
		}
	}

	if (mmcell.mm.n == 0 && mmcell.sm.n == 0) {
		strcpy(errmsg, "no molecule is defined in the cell"); show_infor(errmsg); return;
	}
	if (md_mode == MULTI_MM_FREE_CELL_MD) {
		sprintf(errmsg, "this function works for multi-mol + periodic cell"); show_infor(errmsg); return;
	}
	else if (md_mode == MULTI_MM_PERIDIC_CELL_MD) {
		molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell
		// using Ewald_Sum to calculate interaction
	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		_EwaldSum_recp_::InitEwaldSumPars(cmm->xd[0]);
	#endif
	}

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	// since molecule mass center position could be shifted, distribute cluster in CMM again 
	bool reset_cell_chain = true;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		reset_cell_chain = (nm == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.mm.m[nm], *cmm, reset_cell_chain); // set base cluster at the center cell
	}
	reset_cell_chain = (mmcell.mm.n == 0 ? true : false);
	cmm_init_distribute_cluster(mmcell.sm.m, mmcell.sm.n, *cmm, reset_cell_chain); // add simple solid molecule into cmm
	cmm_init_subcell<BASIC_CLUSTER>(*cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
	if (!cmm->setup_basic_cell_matrix()) {
		sprintf(::errmsg, "failure to construct basic-cell matrix for CMM"); show_msg(::errmsg);
		return;
	}
	// *****************  over  ******************

	//show_all_molecules();

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;

	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;

	bool Overlap = false;
	int nstart = 0, nEachThread = 0, nmol = 0;

	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	LFVERLET_MM_MD *lfmdpar = NULL;
	int nMaxCluster = mmcell.max_ncluster();

	int nMMThreads = 0, nMMEachThread = 0;
	MM_MD_THREAD_PAR mm_mdpar[MAX_THREADS];

	if (mmcell.mm.n > 0) {
		nMMEachThread = (::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1) / MAX_THREADS;
		if (nMMEachThread > 0) nMMThreads = MAX_THREADS;
		else {nMMThreads = ::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1; nMMEachThread = 1;}

		nstart = ::n1MM[::mpi_id];
		for (i = 0; i < nMMThreads; i++) {
			if (i == nMMThreads - 1) nmol = ::n2MM[::mpi_id] + 1 - nstart;
			else nmol = nMMEachThread;
			mm_mdpar[i].set(mmcell.mm.m + nstart, mmcell.lfmd.m + nstart, false, nmol, i);
			mm_mdpar[i].set_buffer_memory(nMaxCluster);
			nstart += nmol;
		}
	}

	SMOLECULE *sm = NULL;
	LFVERLET_SM_MD *lfsmdpar = NULL;

	int nSMThreads = 0, nSMEachThread = 0;
	SM_MD_THREAD_PAR sm_mdpar[MAX_THREADS];

	if (mmcell.sm.n > 0) {
		nSMEachThread = (::n2SM[::mpi_id] - ::n1SM[::mpi_id] + 1) / MAX_THREADS;
		if (nSMEachThread > 0) nSMThreads = MAX_THREADS;
		else {nSMThreads = (::n2SM[::mpi_id] - ::n1SM[::mpi_id] + 1); nSMEachThread = 1;}

		nstart = ::n1SM[::mpi_id];
		for (i = 0; i < nSMThreads; i++) {
			if (i == nSMThreads - 1) nmol = ::n2SM[::mpi_id] + 1 - nstart;
			else nmol = nSMEachThread;
			sm_mdpar[i].set(mmcell.sm.m + nstart, mmcell.lfsmd.m + nstart, nmol, i);
			nstart += nmol;
		}
	}

	bool bHooverDyn = true;

	VECTOR<6> Iw0;

	int nConvg = 0;

	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mm = mmcell.mm.m + nm;
		mm->calClusterSolvationRadius(); // it uses the cluster's mass center position
		// the macro-molecule not running in this processor, does not have bkm, dyn and invNE for the clusters
		if (nm >= ::n1MM[::mpi_id] && nm <= ::n2MM[::mpi_id]) {
			InitVelocity(*mm);
			mm->resetDynPars(); // set force & torque to 0
		}
	}

	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		if (nm >= ::n1SM[::mpi_id] && nm <= ::n2SM[::mpi_id]) {
			InitVelocity(*sm);
			// set force & torque to 0
			V6zero(sm->c->dyn->fc)
		}
	}

	CHAIN<short> *ch_cluster_id = NULL;
	strcpy(::errmsg, "\n"); show_log(::errmsg, true);
	sprintf(::errmsg, "--------------------------------------------------------"); show_log(::errmsg, true);
	sprintf(::errmsg, "  CLUSTER UID: ita, cd_translation, cd_rotation"); show_log(::errmsg, true);
	sprintf(::errmsg, "  Brownian force: f = -ita * V - cd_trans/rot * |V| * V"); show_log(::errmsg, true);
	sprintf(::errmsg, "  Brownian force unit: mass * Angs/fs^2"); show_log(::errmsg, true);
	sprintf(::errmsg, "--------------------------------------------------------"); show_log(::errmsg, true);
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mm = mmcell.mm.m + nm;
		mm->show_cluster_BrownianForcePars(&ch_cluster_id);
	}
	release_chain<short>(&ch_cluster_id, true);
	sprintf(::errmsg, "------------------------- OVER -------------------------"); show_log(::errmsg, true);
	strcpy(::errmsg, "\n"); show_log(::errmsg, true);


	mmcell.Init(); //copy the initial torsion angle, speed to the mdpar

	long LOOPS = 0;
	int save_indx = 0;
	MD_SAVE *mdsave = NULL;
	if (mmcell.mm.n > 0) {LOOPS = mmcell.lfmd.m[0].LOOPS; mdsave = mmcell.lfmd.m[0].mmsave.mdsave;}
	else if (mmcell.sm.n > 0) {LOOPS = mmcell.lfsmd.m[0].LOOPS; mdsave = mmcell.lfsmd.m[0].smsave.mdsave;}

	NEIMO_THREAD_VARS *neimo_matrix_cal_var = new NEIMO_THREAD_VARS[::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1];

	bool bCheckClusterInCMM = true;

	MM_EF_VAR mef_var;
	MM_EF_RES mef_res;
	MOperateVAR<CLUSTER, MMOLECULE, SMOLECULE, PMOLECULE, CMM_CELL3D<BASIC_CLUSTER>, MMOL_MD_CELL> mvar;
	mvar.cell = cmm; mvar.mmcell = &mmcell;
	mvar.mm = mmcell.mm.m + ::n1MM[::mpi_id]; mvar.nMM = ::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1;
	if (mmcell.sm.m != NULL) {
		mvar.sm = mmcell.sm.m + ::n1SM[::mpi_id]; mvar.nSM = ::n2SM[::mpi_id] - ::n1SM[::mpi_id] + 1;
	}
	else {mvar.sm = NULL; mvar.nSM = 0;}
	mvar.pm = NULL; mvar.nPM = 0;
	mvar.nMainThreads = MAX_THREADS;
	mvar.nMinorThreads = 2;
	if (mvar.nMinorThreads < 2) mvar.nMinorThreads = 2;

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command
		if (iCheckCMM == 0) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		MPI_Barrier(MPI_COMM_WORLD);
		// get the structure calculated in other processors
		//if (!MPI_ExchangeStructInfor(mmcell)) {
		if (!MPI_ExchangeStruct(mmcell)) {
			sprintf(errmsg, "failure to exchange the structure"); show_infor(errmsg);
			break;
		}
		molecule_check_periodic_cell(mmcell);

		show_all_molecules();

	#define _interact 0
	#if _interact == 0 // CMM + EwaldSum
		#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		// calculate the electro-multipole of each cluster in the cell, EwaldSum
		//CellOperate1<BASIC_CLUSTER, bool>((void*)(&calClusterElectroMultipole), *cmm, 0, cmm->nCells, true, MAX_THREADS);
		
		// calculate the cluster in this processor
		// and exchange the multipole information of clusters in other processors
		mpi_calClusterMultipole(mmcell, true);
		MPI_Barrier(MPI_COMM_WORLD);
		if (!MPI_ExchangeClusterMultipoleInfor(mmcell, true)) {
			sprintf(errmsg, "failure to exchange the cluster's multipole"); show_infor(errmsg);
			break;
		}

		// calculate the multipole of each cell
		calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, false, true, MAX_THREADS);
		#endif
	#endif

		// calculate the inertia moment of each cluster and the whole molecule
		for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
			if (nm >= mmcell.mm.n) break;
			mm = mmcell.mm.m + nm;
			CalMolInertiaMoment(*mm);
			neimo_matrix_cal_var[nm - ::n1MM[::mpi_id]].set_pars(mm, 0);
			Start_NEIMO_MATRIX_CALC(*mm, neimo_matrix_cal_var[nm - ::n1MM[::mpi_id]]);
		}
		// sm inertia moment will be calculated in SM_LFVERLET_MD_thread_func before Newton law calculation
		//show_log("CMM checking ...", true);
		// check cluster in cell, using vp of each cluster only
		if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
		else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
		
		//show_log("cal interaction ...", true);
		// find out relationship
		// calculate the interaction, and the implicit solvation
		mef_var.checkCMM = (iCheckCMM == 0 ? true : false);
		mef_var.bCheckClusterInCMM = bCheckClusterInCMM;
		MOperateInCell<CLUSTER, MMOLECULE, SMOLECULE, CMM_CELL3D<BASIC_CLUSTER>, MMOL_MD_CELL, MM_EF_VAR, MM_EF_RES>(mvar, 
			(void*)(&multi_mol_distribute_force), mef_var, mef_res);
		Ep = mef_res.Ep;
		Overlap = mef_res.OVERLAP;
		if (bImplicitSolvent) {
			// implicit solvation energy / force -- force = d(U1+U2)/dr, dependent on the related two clusters
			/*
			MPI_ExchangeClusterImpSolvRshipPars(mmcell); // have to exchange it so that we can get the parameters for related cluster
			for (nm = ::n1MM[::mpi_id]; nm <= ::n1MM[::mpi_id]; nm++) {
				MacroMoleculeOperate<MMOLECULE>((void*)(&MM_ImpSolvForce), mmcell.mm.m[nm], 0, mmcell.mm.m[nm].nCluster, MAX_THREADS);
			}
			*/
			//show_log("cal implicit solvation ...", true);
			// Another way to realize force = d(U1 + dU2) / dr
			//if (::local_relax) {
			impsolv_initilize(mmcell, cmm, bCheckClusterInCMM);
			for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
				MacroMoleculeOperate<MMOLECULE>((void*)(&MM_ImpSolvForce), mmcell.mm.m[nm], 0, mmcell.mm.m[nm].nCluster - 1, MAX_THREADS);
			}
			//}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		// sum Ep @ master
		Ep_t = Ep;
		MPI_Reduce(&Ep_t, &Ep, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		//bHooverDyn = (Overlap ? false : true);
		bHooverDyn = true;

		for (i = 0; i < nMMThreads; i++) mm_mdpar[i].set_loop(nloop, bHooverDyn, Ep);
		for (i = 0; i < nSMThreads; i++) sm_mdpar[i].set_loop(nloop, bHooverDyn, Ep);

		// wait until the neimo-matrix calculation is done
		for (nm = 0; nm < ::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1; nm++) {
			if (nm >= mmcell.mm.n) break;
#if _SYS_ == _WINDOWS_SYS_
			WaitForSingleObject(neimo_matrix_cal_var[nm].hMMThread, INFINITE);
#elif _SYS_ == _LINUX_SYS_
			pthread_join(neimo_matrix_cal_var[nm].thread, NULL);
#endif
		}
		//show_log("moving mols ...", true);
		MM_MULTI_THREAD_MD(mm_mdpar, nMMThreads);
		if (mmcell.sm.m != NULL) SM_MULTI_THREAD_MD(sm_mdpar, nSMThreads);

		for (i = 0; i < nMMThreads; i++) {
			if (strlen(mm_mdpar[i].msg_show) > 1) show_infor(mm_mdpar[i].msg_show, true);
			if (strlen(mm_mdpar[i].msg_log) > 1) mlog.show(mm_mdpar[i].msg_log, true);
		}
		for (i = 0; i < nSMThreads; i++) {
			if (strlen(sm_mdpar[i].msg_show) > 1) show_infor(sm_mdpar[i].msg_show, true);
			if (strlen(sm_mdpar[i].msg_log) > 1) mlog.show(sm_mdpar[i].msg_log, true);
		}

		//show_log("exchange kinetic energy ...", true);
		MPI_ExchangeKineticEnergy(mmcell); // the kinetic energy of each molecule

		if (::mpi_id == 0) mmcell.MDSave(nloop, (float)Ep); // master

		Ek_t = 0;
		for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
			Ek_t += mmcell.mm.m[nm].Ek;
			Verlet(mmcell.lfmd.m[nm], mmcell.mm.m[nm]); // Verlet algorithm to calculate the ta and tp for next timestep
		}
		if (mmcell.sm.m != NULL) {
			for (nm = ::n1SM[::mpi_id]; nm <= ::n2SM[::mpi_id]; nm++) {
				Ek_t += mmcell.sm.m[nm].Ek;
				Verlet(mmcell.lfsmd.m[nm], mmcell.sm.m[nm]); // Verlet algorithm to calculate the ta and tp for next timestep
			}
		}

		//show_log("check torsion angle ...", true);
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			//for (nm = ::n1MM[::mpi_id]; nm <= ::n2MM[::mpi_id]; nm++) {
			for (nm = 0; nm < mmcell.mm.n; nm++) {
				mm = mmcell.mm.m + nm;
				lfmdpar = mmcell.lfmd.m + nm;
				check_angle(*mm, *lfmdpar);
			}
		}
		iCheckAngle++; if (iCheckAngle >= nCheckAngle) iCheckAngle = 0;

		MPI_Barrier(MPI_COMM_WORLD);
		// gather mdpar into master
		MPI_GatherMDInfor(mmcell.lfmd.m, mmcell.lfmd.n, mmcell.lfsmd.m, mmcell.lfsmd.n); 

		//MPI_Barrier(MPI_COMM_WORLD);
		// get the structure calculated in other processors
		//if (!MPI_ExchangeStructInfor(mmcell)) {
		//	sprintf(errmsg, "failure to exchange the structure"); show_infor(errmsg);
		//	break;
		//}

		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mmcell.CalibrateTorsionAngle();

		nloop++;
	}
	for (i = 0; i < nMMThreads; i++) mm_mdpar[i].reset();
	for (i = 0; i < nSMThreads; i++) sm_mdpar[i].reset();
	if (neimo_matrix_cal_var != NULL) {delete[] neimo_matrix_cal_var; neimo_matrix_cal_var = NULL;}
}


#endif // _MPI_
