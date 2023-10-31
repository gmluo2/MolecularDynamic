/********************************************************************
******          SINGLE POLYMER SIMULATION WITH MPICH2          ******
********************************************************************/

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

#include "ZMatrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "cluster.h"
#include "Interaction.h"
#include "Interaction1.h"
#include "Interact2.h"
#include "NEIMO.h"
#include "MD.h"
#include "MD_POLYMER.h"

#include "var.h"
#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

#ifdef SEEK_SET
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include "mpi.h"
#endif

extern bool trace;

extern void get_Frame_Hoover_force(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, double Ethermal_trans, double Ethermal_rot, double &ksi_frame_trans, double &vksi_frame_trans, double &ksi_frame_rot, double &vksi_frame_rot, VECTOR<6> &f);

// set cluster's force & torque due to the Hoover correction on acceleration of the mass center
extern void setClusterHooverForce(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, double Ethermal_trans, double Ethermal_rot, double &ksi_frame_trans, double &vksi_frame_trans, double &ksi_frame_rot, double &vksi_frame_rot, double ksi_internal);

extern double max_diff(VECTOR<6> *v1, VECTOR<6> *v2, int n);

// initilize the local cordinates and calculate the inertia moment of each cluster and the whole molecule

extern void CalClusterInertiaMoment(MMOLECULE *mm, int n1, int n2);

extern void CalMolInertiaMoment(MMOLECULE &mm);

// calculate Coriolis acceleration term and gyroscopic spatial force
extern void CalMol_ab(MMOLECULE &mm);

extern void CalInternalKinHooverForce(MMOLECULE &mm);

extern double PredictKsi(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, double dEk);

extern bool LeapFrog_SelfConsistent_velocity_acceleration(bool bHooverDyn, MMOLECULE &mm, LFVERLET_MM_MD& mdpar, int md_loop, int &nconvg_loops, double Ethermal_trans, double Ethermal_rot, VECTOR<6> *vect1, VECTOR<6> *vect2, double *dv1, double *dv2, bool local_relax, char *show_msg, char *log_msg, int n_msglen, bool init_md);

extern void InitVelocity(MMOLECULE &mm);

extern void ZeroVelocity(MMOLECULE &mm);

//adjusting free base velocity so that the total spacial moment is consistent
extern void MakeFreeBaseVelocityConsistent(MMOLECULE &mm);

void MPI_ExchangeClusterImpSolvRshipPars(MMOLECULE& mm) {
	CLUSTER *pc = NULL;
	double *dp = NULL;
#define NC_MAX   2000
	int nc = 0;
	int i = 0, iv = 0, iproc = 0;
	double dbf[NC_MAX];

	char msg[256] = "\0";

	for (iproc = 0; iproc < ::nMPI; iproc++) {
		// cluster's parameters
		nc = ::n1MM[iproc];
		while (nc <= ::n2MM[iproc]) {
			dp = dbf; iv = 0;
			while (iv < NC_MAX) {
				if (::mpi_id == iproc) {
					dp = dbf + iv; dp[0] = mm.cluster[nc].SASA; iv += 1;
				}
				else iv += 1;
				nc += 1;
				if (nc > ::n2MM[iproc]) break;
			}
			MPI_Bcast(dbf, iv, MPI_DOUBLE, iproc, MPI_COMM_WORLD);
			if (::mpi_id != iproc) {
				// following iv calculation has to be consistent as above in ::mpi_id == proc_id
				for (i = 0; i < iv; i++) {
					mm.cluster[nc - iv + i].SASA = dbf[i];
				}
			}
			if (nc > ::n2MM[iproc]) break;
		}
	}
#undef NC_MAX
}

void MPI_GatherClusterForce(MMOLECULE *mm) {
#define NC_MAX   200
#define V_MAX    1200
	int nc = 0, nv = 0, ncs = 0;
	double dbf[V_MAX];
	CLUSTER *pc = NULL;
	double *fc = NULL;

	int nSend[3] = {0, 0, 0};
	int total_clusters = 0;
	MPI_Status status;
	int nget = 0;
	char msg[256] = "\0";

	if (::mpi_id != 0) { //this is not the master
		nSend[0] = ::n1MM[::mpi_id]; nSend[1] = ::n2MM[::mpi_id];
		total_clusters = ::n2MM[::mpi_id] - ::n1MM[::mpi_id] + 1;
		nSend[2] = total_clusters / NC_MAX;
		if (nSend[2] * NC_MAX < total_clusters) nSend[2] ++;
		MPI_Send(nSend, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);

		nc = ::n1MM[::mpi_id];
		while (nc <= ::n2MM[::mpi_id]) {
			ncs = 0; nv = 0;
			while (ncs < NC_MAX && nc <= ::n2MM[::mpi_id]) {
				if (nc >= mm->nCluster) break;
				pc = mm->cluster + nc; fc = pc->dyn->fc.v;
				dbf[nv+0] = fc[0]; dbf[nv+1] = fc[1]; dbf[nv+2] = fc[2]; 
				dbf[nv+3] = fc[3]; dbf[nv+4] = fc[4]; dbf[nv+5] = fc[5];
				nv += 6; ncs++; nc++;
			}
			MPI_Send(dbf, nv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // send to the master
		}
	}
	else {
		for (int isource = 1; isource < ::nMPI; isource++) {
			MPI_Recv(nSend, 3, MPI_INT, isource, 0, MPI_COMM_WORLD, &status);
			nc = nSend[0]; 
			for (int ipkg = 0; ipkg < nSend[2]; ipkg++) {
				MPI_Recv(dbf, V_MAX, MPI_DOUBLE, isource, 0, MPI_COMM_WORLD, &status);
				ncs = 0; MPI_Get_count(&status, MPI_DOUBLE, &nget);
				while (ncs < nget) {
					if (nc >= mm->nCluster) {
						sprintf(msg, "ERROR : MPI get force of cluster %d", nc); show_log(msg, true);
					}
					else if (nc > nSend[1]) {
						sprintf(msg, "ERROR : MPI get force of cluster %d, out of range [%d, %d]", nc, nSend[0], nSend[1]); show_log(msg, true);
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
		}
	}

	// net force
	VECTOR<6> fc_net;
	MPI_Reduce(mm->mDynKin.fc.v, fc_net.v, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (::mpi_id == 0) {
		V62V6(fc_net, mm->mDynKin.fc)
	}

#undef NC_MAX
#undef V_MAX
}

void MPI_DistributeMovePars(VECTOR<6> &dv_base, double *dta, int NC) {
#define NC_MAX   2000
	int nc = 0, ncs = 0;
	double dbf[NC_MAX];

	int nSend = NC / NC_MAX;
	if (nSend * NC_MAX < NC) nSend++;
	MPI_Status status;
	int nget = 0;
	char msg[256] = "\0";

	if (::mpi_id != 0) { //this is slave
		MPI_Recv(dv_base.v, 6, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		nc = 0;
		for (int ipkg = 0; ipkg < nSend; ipkg++) {
			MPI_Recv(dbf, NC_MAX, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_DOUBLE, &nget);
			for (int iget = 0; iget < nget; iget++) {
				if (nc >= NC) {
					sprintf(msg, "ERROR: we get dta more than the size"); show_log(msg, true); 
				}
				else dta[nc] = dbf[iget];
				nc++;
			}
		}
	}
	else { // master
		for (int idest = 1; idest < ::nMPI; idest++) {
			MPI_Send(dv_base.v, 6, MPI_DOUBLE, idest, 0, MPI_COMM_WORLD);
			nc = 0;
			for (int ipkg = 0; ipkg < nSend; ipkg++) {
				ncs = 0;
				for (int i = 0; i < NC_MAX; i++) {
					dbf[i] = dta[nc];
					nc++; ncs++;
					if (nc >= NC) break;
				}
				MPI_Send(dbf, ncs, MPI_DOUBLE, idest, 0, MPI_COMM_WORLD);
			}
		}
	}

#undef NC_MAX
}

void ClusterInertialCord(MMOLECULE *mm, int n1, int n2) {
	int nc, i, j;
	VECTOR3 dr;
	CLUSTER *pc = NULL;
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		pc->InertialCord(); // calcualte the inertial cordinate of each atoms in each cluster
		V3minusV3(pc->cm, (*(pc->Op)), pc->cm0)
		//pc->MassCenter();
		pc->calJacobianMatrix();

		for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm->mDynKin.cm.v[i];
		PHAI(pc->invNE->phai_cm2me, dr, i, j)
		M6T2M6(pc->invNE->phai_cm2me, pc->invNE->phai_T_cm2me, i, j)
		for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - CLUSTER_ORIGIN(mm->base, i);
		PHAI(pc->invNE->phai_base2me, dr, i, j)
		M6T2M6(pc->invNE->phai_base2me, pc->invNE->phai_T_base2me, i, j)

		//pc->InertiaTensor();
	}
}

void MolInertialCord(MMOLECULE &mm) {
	int nThreads = 1;
	if (mm.nCluster > NCLUSTERS_PARALLEL) nThreads = MAX_THREADS;
	MacroMoleculeOperate<MMOLECULE>((void*)(&ClusterInertialCord), mm, 0, mm.nCluster, nThreads);
}

void MPI_POLYMER(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, CMM_CELL3D<BASIC_CLUSTER> *cmm) {
	mm.setup_kin_dyn();
	char msg_show[5120] = "\0", msg_log[5120] = "\0";
	int n_msglen = 5120;

	long nloop = 0;
	int nc = 0, ncell = 0;
	double Ek = 0, Ep = 0, Etorsion = 0, myEp = 0;

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;

	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;

	bool Overlap = false;

	VECTOR<6> *vect1 = new VECTOR<6>[mm.nCluster];
	VECTOR<6> *vect2 = new VECTOR<6>[mm.nCluster];

	double *dv1 = new double[mm.nCluster];
	double *dv2 = new double[mm.nCluster];

	int nConvg = 0;

	bool bHooverDyn = false;

	CLUSTER *pc = NULL;
	MD_SAVE *mdsave = mdpar.mmsave.mdsave;
	if (mdsave == NULL) {sprintf(errmsg, "no MD saving object is defined"); show_msg(errmsg); return;}
	sprintf(refresh_mol_xyz_struct_fname, "%s.xyz", mdpar.mmsave.mdsave->title);

	{
		// set mass center to {0, 0, 0}
		VECTOR<3> cm;
		mm.SetMassCenter(cm.v, true, true);
	}

	InitVelocity(mm);

	mm.calClusterSolvationRadius(); // it uses the cluster's mass center position
	//mm.setup_Solvation_neighbors(::n_neighbor_Solvation);
	//mm.setup_Solvation_neighbors(::solv_cut);
	mm.calClusterSolvationRadius(); // it uses the cluster's mass center position
	CHAIN<short> *ch_cluster_id = NULL;
	strcpy(msg_show, "\n"); show_log(msg_show, true);
	sprintf(msg_show, "--------------------------------------------------------"); show_log(msg_show, true);
	sprintf(msg_show, "  CLUSTER UID: ita, cd_translation, cd_rotation"); show_log(msg_show, true);
	sprintf(msg_show, "  Brownian force: f = -ita * V - cd_trans/rot * |V| * V"); show_log(msg_show, true);
	sprintf(msg_show, "  Brownian force unit: mass * Angs/fs^2"); show_log(msg_show, true);
	sprintf(msg_show, "--------------------------------------------------------"); show_log(msg_show, true);
	mm.show_cluster_BrownianForcePars(&ch_cluster_id);
	release_chain<short>(&ch_cluster_id, true);
	sprintf(msg_show, "------------------------ OVER --------------------------"); show_log(msg_show, true);
	strcpy(msg_show, "\n"); show_log(msg_show, true);

	mdpar.Init(mm); // copy initial torsion angle and speed to mdpar
	
	mm.resetDynPars(); // set force & torque to 0

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	bool init_md = false;
	int nCMMCheck = cmm->nCheckCluster;

	VECTOR3 mIw, mP;
	//float rand_torque = 0;

	NEIMO_THREAD_VARS neimo_matrix_cal_var;
	VECTOR<6> dv_base;
	double *dta = new double[mm.nCluster];
	for (nc = 0; nc < mm.nCluster; nc++) dta[nc] = 0;

	bool bCheckClusterInCMM = true;

	while (nloop < mdpar.LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command
		init_md = (nloop == 0 ? true : false); 
		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM) {
		{
			// set mass center to {0, 0, 0}
			VECTOR<3> cm;
			mm.SetMassCenter(cm.v, true, true);
		}
		}
		
		if (::mpi_id == 0) { // do this only on master
			// calculate the inertia moment of each cluster and the whole molecule
			CalMolInertiaMoment(mm);
		}
		else MolInertialCord(mm);

#if _CHECK_TIME_STAMP_ == 1
		time.start();
#endif
		if (::mpi_id == 0) { // do this only on master
			// we have to do it always before we start new thread, so that the thread-related variable will be refreshed
			// calculate the neimo-matrix, related to the innertial matrix only
			neimo_matrix_cal_var.set_pars(&mm, 0);
			Start_NEIMO_MATRIX_CALC(mm, neimo_matrix_cal_var);
		}

		for (nc = 0; nc < mm.nCluster; nc++) {
			V6zero(mm.cluster[nc].dyn->fc)
			mm.cluster[nc].dyn->overlap = false;
		}
		Ep = 0;

	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		// calculate the electro-multipole of each cluster in the cell, not EwaldSum
		CellOperate1<BASIC_CLUSTER, bool>((void*)(&calClusterElectroMultipole), *cmm, 0, cmm->nCells, false, MAX_THREADS);
	#endif
		// check cluster in cell and calculate the multipole of each cell
		if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
		else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);

	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, false, false, MAX_THREADS);
	#endif

		if (iCheckCMM == 0) {
			// setup relationship of the cluster used in this MPI process
			MMCheckRelshipInCell(mm, ::n1MM[::mpi_id], ::n2MM[::mpi_id], *cmm, MAX_THREADS); 
			//sprintf(::errmsg, "check relationship takes %d", time.elapse());
			//show_infor(errmsg);
		}
		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "multipole calculation takes %d ms", time.glance());
		show_infor(errmsg);
#endif
		// parallel calculation of interaction with CMM, and multipole-interaction
		MMClusterInteraction<MMOLECULE, BASIC_CLUSTER>((void*)(&FreeCellCMMInteraction1_SingleMM), *cmm, mm, ::n1MM[::mpi_id], ::n2MM[::mpi_id], Ep, Overlap, false, MAX_THREADS);

#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "cluster interact.  calculation takes %d ms", time.glance());
		show_infor(errmsg);
#endif

		// implicit solvation energy / force, for the clusters used in this MPI process
		if (bImplicitSolvent) {
			// neighbors of cluster used in this MPI process for implicit solvation

			// implicit solvation energy / force -- dependent on the cluster itself
			/*
			if (bCheckClusterInCMM) MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MMClusterImpSolvNeighborsInCell), *cmm, mm, ::n1MM[::mpi_id], ::n2MM[::mpi_id], MAX_THREADS);
			MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MM_ImplicitSolvationForce), *cmm, mm, ::n1MM[::mpi_id], ::n2MM[::mpi_id], MAX_THREADS);
			*/
			if (bCheckClusterInCMM) MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MMClusterImpSolvNeighborsInCell), *cmm, mm, ::n1MM[::mpi_id], ::n2MM[::mpi_id], MAX_THREADS);
			MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MM_ImpSolvSASA), *cmm, mm, ::n1MM[::mpi_id], ::n2MM[::mpi_id], MAX_THREADS);
			for (nc = ::n1MM[::mpi_id]; nc <= ::n2MM[::mpi_id]; nc++) {
				if (nc >= mm.nCluster) break;
				Ep += mm.cluster[nc].E_ImplicitSolvation;
			}
			MPI_ExchangeClusterImpSolvRshipPars(mm); // have to exchange it so that we can get the parameters for related cluster
			MacroMoleculeOperate<MMOLECULE>((void*)(&MM_ImpSolvForce), mm, ::n1MM[::mpi_id], ::n2MM[::mpi_id], MAX_THREADS);
		}

		if (bRandomForce) gen_add_random_force1(mm, ::n1MM[::mpi_id], ::n2MM[::mpi_id], randf, 0.35);
		CalNetForce(mm, ::n1MM[::mpi_id], ::n2MM[::mpi_id]); // the net force on molecular mass center

		MPI_GatherClusterForce(&mm); // get the force on clusters from each MPI process, and the netforce
		myEp = Ep;
		MPI_Reduce(&myEp, &Ep, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (::mpi_id == 0) {
			// Because torsion is a torque along hinge, on both clusters on the hinge,
			// the net force/torque of the whole macro-molecule will be 0.
			// So we can calculate the torsion here
			//mm.calTorsion(Etorsion);
			MacroMoleculeOperate2<MMOLECULE, double>((void*)(&calTorsion), mm, 0, mm.nCluster, MAX_THREADS, Etorsion);
			Ep += Etorsion;
		}

		if (scale_force && ::mpi_id == 0) {
			switch (interact_cal_method) {
			case _CMM:
				CellOperate<BASIC_CLUSTER>((void*)(&SCALE_FORCE_CMM), *cmm, 0, cmm->acell.n, MAX_THREADS);
				break;
			default:
				scale_force_mmolecule(mm);
				break;
			}
		}

#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "interaction calculation takes %d ms", time.glance());
		show_infor(errmsg);
#endif

		if (::mpi_id == 0) { // do this only on master
			mdpar.mmsave.save_kinetic_pars(mm, nloop, true, NULL, true);
			mdsave->mESave.one_more();
			mdsave->save_energy(mm.Ek, nloop, false, NULL, true);
			mdsave->save_energy((float)Ep, nloop, false, NULL, false);// save Ep
			if (mdsave->mESave.bsave()) mdsave->mESave.save_over();
			mdpar.save_mol_dyn_pars(mm, nloop, true, NULL, true);

			show_molecule_struct(mm, refresh_mol_xyz_struct_fname);

			//sprintf(errmsg, "[%f, %f, %f]", float(mm.mDynKin.cm.v[0]), float(mm.mDynKin.cm.v[1]), float(mm.mDynKin.cm.v[2]));
			//sprintf(errmsg, "[%f, %f, %f]", float(mm.base->atom[0].r.v[0]), float(mm.base->atom[0].r.v[1]), float(mm.base->atom[0].r.v[2]));
			//show_infor(errmsg);
		}

#if _CHECK_TIME_STAMP_ == 1
		//time.start();
#endif

		if (::mpi_id == 0) { // do this only on master
			// wait until the neimo-matrix calculation is done
#if _SYS_ == _WINDOWS_SYS_
			WaitForSingleObject(neimo_matrix_cal_var.hMMThread, INFINITE);
#elif _SYS_ == _LINUX_SYS_
			pthread_join(neimo_matrix_cal_var.thread, NULL);
#endif

#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "matrix calculation takes %d ms", time.elapse());
			show_infor(errmsg);
#endif

#if _CHECK_TIME_STAMP_ == 1
			time.start();
#endif

			//bHooverDyn = (Overlap ? false : true);
			bHooverDyn = true;
			strcpy(msg_show, "\0"); strcpy(msg_log, "\0");
			// to get self-consistent velcoty and acceleration at this step based on leap-frog verlert algorithm
			// since we are simulating single macro-molecule, set whole molecule translating thermal energy 0, rotation thermal energy 1.5 kT
			LeapFrog_SelfConsistent_velocity_acceleration(bHooverDyn, mm, mdpar, nloop, nConvg, vect1, vect2, dv1, dv2, ::local_relax, msg_show, msg_log, n_msglen, init_md);
			if (strlen(msg_show) > 1) show_infor(msg_show, true);
			if (strlen(msg_log) > 1) mlog.show(msg_log, true);
			Ek = mm.Ek;
		
			//Verlet(mdpar, mm); // Verlet algorithm to calculate the ta and tp for next timestep
			VerletStep(mdpar, mm, dv_base, dta);
		}
		MPI_DistributeMovePars(dv_base, dta, mm.nCluster);
		VerletMove(mdpar, dv_base, dta, mm);

#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "NEIMO + Verlet take %d", time.elapse());
		show_infor(errmsg);
#endif

		// make sure torsion angle in range [-PI, PI]
		iCheckAngle++;
		if (iCheckAngle >= nCheckAngle) {
			iCheckAngle = 0;
			check_angle(mm, mdpar);
		}

		if (::mpi_id == 0) {
			//show_molecule_struct(mm, refresh_mol_xyz_struct_fname);
			show_loop_infor(nloop, mdpar.mmsave.mdsave->mESave.save_indx(), mdpar.LOOPS - nloop, (float)Ek, (float)Ep);
		}

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mdpar.CalibrateTorsionAngle(mm);

		nloop++;
	}
	if (vect1 != NULL) {delete[] vect1; vect1 = NULL;}
	if (vect2 != NULL) {delete[] vect2; vect2 = NULL;}
	if (dv1 != NULL) {delete[] dv1; dv1 = NULL;}
	if (dv2 != NULL) {delete[] dv2; dv2 = NULL;}
	if (dta != NULL) {delete[] dta; dta = NULL;}
}

#endif
