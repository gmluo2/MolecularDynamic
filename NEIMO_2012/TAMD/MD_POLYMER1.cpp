
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

#include "ZMatrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
using namespace _cmm_3d_;
#include "cluster.h"
#include "Interaction.h"
#include "Interact2.h"
#include "NEIMO.h"
#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"

#include "MD_POLYMER.h"

#include "var.h"
#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

extern bool trace;


bool FixedBase_LeapFrog_SelfConsistent_velocity_acceleration(bool bHooverDyn, MMOLECULE &mm, LFVERLET_MM_MD& mdpar, int md_loop, int &nconvg_loops, VECTOR<6> *vect1, VECTOR<6> *vect2, double *dv1, double *dv2, bool local_relax, char *show_msg, char *log_msg, int n_msglen) {
	char msg[512] = "\0";

	int nc = 0, i = 0, j;
	CLUSTER *pc = NULL;
	double T = 0, Ek0 = 0, Ek = 0, dEk = 0;
	double Ethermal = (double)(mm.nCluster - 1) * ::Ek0Internal; // / 2; // in unit kT

	double dEmax = Ethermal; // maximum kinetic energy change in one step. Each freedom has kinetic energy 0.5kT only

	double vmax_diff = 0, wmax_diff = 0;
	double vConvg = 0.005 / mdpar.dt; //0.001 / mdpar.dt;
	if (vConvg < 2e-4) vConvg = 2e-4;
	double vConvg0 = 0;
	
	bool convergent = false;
	int nmax_convg = 8;
	int MD_LOOP = md_loop;

	VECTOR<6> *tp_last = vect1, *tp_new = vect2;
	double *mtp = dv1, *mtpp = dv2;
	VECTOR<6> mv0, malpha0;

	bool tp_scaled = false;

	double ksi_Hoover = mdpar.iHDyn->nhc->ksi.v[0];
	double vksi_Hoover = mdpar.iHDyn->nhc->vksi.v[0];

	int nconvg = 0;
	bool new_config = true;

	while (!convergent) {
		// speed is varying to get consisten speed and accleration based on leap-frog velert MD 
		// calculate Coriolis acceleration term and gyroscopic spatial force
		CalMol_ab(mm);

		//if (nloop == 0 && convergent == 0) InitKinetic(mdpar, mm);

		for (nc = 0; nc < mm.nCluster; nc++) {
			if (nconvg == 0) {for (i = 0; i < 6; i++) tp_last[nc].v[i] = mm.cluster[nc].bkm->V0.v[i];}
			else {for (i = 0; i < 6; i++) tp_last[nc].v[i] = tp_new[nc].v[i];}
		}

		if (local_relax) {
			Ek0 = 0;
			ksi_Hoover = 0; vksi_Hoover = 0;
		}
		else {
			Ek0 = calKineticEnergy(mm); // kinetic energy of while macro-molecule

			if (bHooverDyn) {
				T = Ek0 / Ethermal;
				mdpar.iHDyn->nhc->set_dT(T - 1);
				mdpar.iHDyn->nhc->verlet_NHC();
				ksi_Hoover = mdpar.iHDyn->nhc->ksi.v[0];
			}
			else {ksi_Hoover = 0;}

			if (bHooverDyn && max_ksi_Hoover > 1e-8 && T > 9) { //ksi_Hoover / max_ksi_Hoover >= 0.98
					/*
					sprintf(msg, " Ek is more thane 9 times of thermal energy. scale to 80%% at loop %d", MD_LOOP);
					if (n_msglen - strlen(show_msg) > 150) {if (strlen(show_msg) > 10) strcat(show_msg, " \r\n"); strcat(show_msg, msg);}
					sprintf(msg, " ksi : [%5.3f]", ksi_Hoover);
					if (n_msglen - strlen(show_msg) > 150) {if (strlen(show_msg) > 10) strcat(show_msg, " \r\n"); strcat(show_msg, msg);}
					*/
					
					// decrease the speed is not efficient enough
					// so scale the kinetic energy back to thermal energy
					// and calculate the a, b since speed changed
			
					Ek0 = ScaleKineticEnergy(mm, Ethermal, false);
					mdpar.iHDyn->nhc->reset_NHC();
					md_loop = 0; // to trigger the case that do not get consistent velocity and accleration

					continue;
			}
		}

		NEIMO_CalHooverForce(mm, ksi_Hoover, mm.mDynKin.f_Hoover);

		NEIMO_FixedBase(mm, (float)ksi_Hoover, new_config);

		if (new_config) new_config = false;

		if (local_relax) {
			Speed_Kinetic(mdpar, mm, true);
			convergent = true;
			return 1;
		}

		if (nconvg == 0) {// backup tpp
			for (nc = 0; nc < mm.nCluster; nc++) {
				mtp[nc] = mm.cluster[nc].km.tp;
				mtpp[nc] = mm.cluster[nc].km.tpp;
			}
			for (i = 0; i < 6; i++) {
				mv0.v[i] = mm.V0.v[i]; malpha0.v[i] = mm.alpha0.v[i];
			}
		}

		if (md_loop == 0) Speed_Kinetic(mdpar, mm, true);
		else {
			// get new speed at current position and +1/2 time step
			Speed_Verlet(mdpar, mm); // reset the speed in mm to the new speed
		}

		for (nc = 0; nc < mm.nCluster; nc++) {
			for (i = 0; i < 6; i++) tp_new[nc].v[i] = mm.cluster[nc].bkm->V0.v[i];
		}
		max_diff(tp_new, tp_last, mm.nCluster, vmax_diff, wmax_diff); 
		if (nconvg == 0) vConvg0 = vmax_diff;
		Ek = calKineticEnergy(mm); // kinetic energy 

		if (md_loop == 0) {convergent = true; break;}
		else if (vmax_diff < vConvg && wmax_diff < vConvg) convergent = true;
		else {
			//dEk = Ek - Ethermal;
			dEk = Ek - Ek0;
			if (dEk > dEmax) {
				//sprintf(msg, " Kinetic energy jumps from %f to %f kT with iteration %d at loop %d", Ek0, Ek, nconvg, MD_LOOP);
				//if (n_msglen - strlen(show_msg) > 150) {if (strlen(show_msg) > 10) strcat(show_msg, "\n"); strcat(show_msg, msg);}
				//if (nconvg > 0) {
				if (nconvg < nmax_convg) {
					//sprintf(msg, " scale to 0.8 of thermal energy at loop %d ", MD_LOOP);
					//if (n_msglen - strlen(show_msg) > 150) {if (strlen(show_msg) > 10) strcat(show_msg, "\n"); strcat(show_msg, msg);}
					Ek0 = ScaleKineticEnergy(mm, Ethermal, false); // do not reset the spacial moment
					Speed_Kinetic(mdpar, mm, false); // do not check free base velocity
					mdpar.iHDyn->nhc->reset_NHC();
					md_loop = 0; // to trigger the case that do not get consistent velocity and accleration
					continue;
					// using the first iteration acceleration, see below
					//nconvg = nmax_convg + 1;
				}
				else {
					sprintf(msg, " ERROR : Unable to handle this problem!");
					if (n_msglen - strlen(show_msg) > 100) {if (strlen(show_msg) > 10) strcat(show_msg, "\n"); strcat(show_msg, msg);}
					nconvg = nmax_convg + 1; 
				}
				break;
			}
		}

		nconvg++;
		if (nconvg > nmax_convg) {
		  if (vmax_diff < vConvg0) convergent = true;
		  else convergent = false;
		  break;
		}
	}

	if (!convergent && Ek > 4 * Ethermal) {
		sprintf(msg, " failure to get consistent speed & acceleration at loop %d ", MD_LOOP);
		if (n_msglen - strlen(show_msg) > 100) {if (strlen(show_msg) > 10) strcat(show_msg, "\n"); strcat(show_msg, msg);}
		// use the acceleration speed based on initial speed
		// torsion angle is not changed in the iteration

		for (nc = 0; nc < mm.nCluster; nc++) {
			mm.cluster[nc].km.tp = mtp[nc];
			mm.cluster[nc].km.tpp = mtpp[nc];		
			//mm.cluster[nc].km.tpp = 0;
		}
		mm.base->km.tp = 0; mm.base->km.tpp = 0;
		for (i = 0; i < 6; i++) {
			mm.V0.v[i] = mv0.v[i];
			mm.alpha0.v[i] = malpha0.v[i];
			//mm.alpha0.v[i] = 0;
		}
		/*
		for (nc = 0; nc < mm.nCluster; nc++) {
			pc = mm.cluster + nc;
			for (j = 0; j < 3; j++) {
				pc->km.V.v[j] = pc->km.tp * pc->up.v[j]; pc->km.V.v[j + 3] = 0;
			}
		}
		*/
		calClusterVelocity0(mm);

		Ek = ScaleKineticEnergy(mm, Ethermal, false);
		//CalMol_ab(mm);
		for (nc = 0; nc < mm.nCluster; nc++) {
			pc = mm.cluster + nc; 
			V6zero(pc->invNE->a) V6zero(pc->invNE->b)
		}
		Speed_Kinetic(mdpar, mm, true);
		mdpar.iHDyn->nhc->reset_NHC();
	}

	mm.Ek = (float)Ek;

	if (bHooverDyn) {
		mdpar.iHDyn->nhc->accept_nextstep();
	}

	return convergent;
}

void InitVelocity_FixedBase(MMOLECULE &mm, double &Ek) {
	int i = 0, j = 0;
	//float unit_E = kT / U_InertiaMoment;

	// calculate the inertia moment of each cluster and the whole molecule
	CalMolInertiaMoment(mm, MAX_THREADS);
	bool Overlap = false;

	V6zero(mm.V0); // set V0 of base 0
	for (i = 0; i < 0; i++) {
		mm.cluster[i].km.tp = 0;
		V6zero(mm.cluster[i].bkm->V0)
	}
	Ek = 0;
}

void ZeroVelocity_FixedBase(MMOLECULE &mm) {
	int nc;
	CLUSTER *pc = NULL;
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		pc->km.tp = 0; V6zero(pc->bkm->V0)
		V6zero(pc->invNE->a) V6zero(pc->invNE->b)
	}
	V6zero(mm.mDynKin.P) V6zero(mm.mDynKin.V)
	V6zero(mm.V0)
}

extern void MM_TorsionInteract(MACROMOL_THREAD_VARS<MMOLECULE> *mjob, TorsionVar *tvar, InteractRes *ires);

void FIXEDBASE_MD_LPOLYMER(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, CMM_CELL3D<BASIC_CLUSTER> *cmm, int interact_cal_method) {
	mm.setup_kin_dyn();
	::interact_cal_method = interact_cal_method; 
	char msg_show[5120] = "\0", msg_log[5120] = "\0";
	int n_msglen = 5120;

	long nloop = 0;
	int nc = 0, ncell = 0, i;
	double Ek = 0, Ep = 0, Etorsion = 0;
	//mm.setup_cluster_matom();
	mm.setup_kin_dyn();

	bool bCMM_StorageChain = false, bClusterRelshipChain = false;
	float vcluster = 1.5f * 1.5f * 1.5f;
	int nCluster_CMM = int(cmm->fw[0] * cmm->fw[1] * cmm->fw[2] / (cmm->nx * cmm->ny * cmm->nz) / vcluster * 1.2 + 0.5) * 27;
	if (nCluster_CMM < 2) nCluster_CMM = 2;
	// for single molecule structure, non-periodical cell, the total related clusters can not be more than the clusters of the macromolecule
	if (nCluster_CMM > mm.nCluster) nCluster_CMM = mm.nCluster + 2;
	if (!bCMM_StorageChain) {
		CMM_array_set_storage<BASIC_CLUSTER>(cmm, nCluster_CMM);
	}
	float rcut = ::rcut_LJ + ::size_cluster * 2;
	int nrship_LJ = int(rcut * rcut * rcut / vcluster * 1.2);
	// for single molecule structure, non-periodical cell, the total related clusters can not be more than the clusters of the macromolecule
	if (nrship_LJ > mm.nCluster) nrship_LJ = mm.nCluster + 2;
	rcut = ::rcut_Ewd + ::size_cluster * 2;
	int nrship_spme = int(rcut * rcut * rcut / vcluster * 1.2);
	// for single molecule structure, non-periodical cell, the total related clusters can not be more than the clusters of the macromolecule
	if (nrship_spme > mm.nCluster) nrship_spme = mm.nCluster + 2;
	int nrship_ImpSolv = int(::dSolv_max * ::dSolv_max * ::dSolv_max / vcluster * 1.2);
	// for single molecule structure, non-periodical cell, the total related clusters can not be more than the clusters of the macromolecule
	if (nrship_ImpSolv > mm.nCluster) nrship_ImpSolv = mm.nCluster + 2;
	if (!bClusterRelshipChain) {
		for (nc = 0; nc < mm.nCluster; nc++) {
			mm.cluster[nc].LJRship.use_array(nrship_LJ);
			mm.cluster[nc].spmeRship.use_array(nrship_spme);
			if (::bImplicitSolvent) mm.cluster[nc].impsolv_rship.use_array(nrship_ImpSolv);
		}
	}

	MACROMOL_THREAD_VARS<MMOLECULE> mmJob[MAX_THREADS];
	assignMMJob<MMOLECULE, MACROMOL_THREAD_VARS<MMOLECULE> >(&mm, mmJob, MAX_THREADS);

	CMM_CELL3D<BASIC_CLUSTER> job_cmm[MAX_THREADS];
	MM_CMM_JOB<MACROMOL_THREAD_VARS<MMOLECULE>, CMM_CELL3D<BASIC_CLUSTER> > mmCMMAllocator[MAX_THREADS];
	MM_INTERACT_JOB job_interact[MAX_THREADS];
	for (nc = 0; nc < MAX_THREADS; nc++) {
		CMM_array_construct_for_storage<BASIC_CLUSTER>(job_cmm + nc, cmm, nCluster_CMM);
		mmCMMAllocator[nc].set_pars(nc, mmJob + nc);
		mmCMMAllocator[nc].cmm = job_cmm + nc;
		
		job_interact[nc].set_pars(nc, mmJob + nc);
		job_interact[nc].cmm = cmm;
		job_interact[nc].var.bVirial = false;
		job_interact[nc].evar.esv1.bEwald = false; // not EwaldSum !
		job_interact[nc].evar.bVirial = false; 
		job_interact[nc].evar.bST_diagonal = false; 
		job_interact[nc].evar.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
		job_interact[nc].evar.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	}

	mm.setup_cluster_fatom();
	mm.setup_hinge_dihedral();
	TorsionVar torsion_var[MAX_THREADS];
	HingeConstraintVar hinge_constraint_var[MAX_THREADS];
	InteractRes mires[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		torsion_var[i].bVirial = false; torsion_var[i].bST_diagonal = true;
		hinge_constraint_var[i].bVirial = false; hinge_constraint_var[i].bST_diagonal = true;
	}

	double Uext = 0; // external electric field

	int nThreads = MAX_THREADS;

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;

	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;

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

	mm.SetBasePosZero();

	if (local_relax) {CalMolInertiaMoment(mm, MAX_THREADS); ZeroVelocity_FixedBase(mm);}
	else InitVelocity_FixedBase(mm, Ek);
	mm.calClusterSolvationRadius(); // it uses the cluster's mass center position
	//mm.setup_Solvation_neighbors(::n_neighbor_Solvation);
	//mm.setup_Solvation_neighbors(::solv_cut);

	mdpar.Init(mm); // copy initial torsion angle and speed to mdpar
	
	mm.resetDynPars(); // set force & torque to 0

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	int nCMMCheck = cmm->nCheckCluster;
	bool bCheckClusterInCMM = true;

	_rand_::RAND_ARRAY<int> randa;

	while (nloop < mdpar.LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command

		if (iCheckCMM == 0) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		// calculate the inertia moment of each cluster and the whole molecule
		CalMolInertiaMoment(mm, MAX_THREADS);

		if (local_relax) ZeroVelocity_FixedBase(mm);

#if _CHECK_TIME_STAMP_ == 1
		time.start();
#endif
		
		extern void InitClusterForce(MMOLECULE *mm);
		InitClusterForce(&mm);
		Ep = 0;
		
		switch (interact_cal_method) {
		case _CMM:
		#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		#if _EWALD_SUM == _MultiPole_EWALD_SUM
			// calculate the electro-multipole of each cluster in the cell, not EwaldSum
			CellOperate1<BASIC_CLUSTER, bool>((void*)(&calClusterElectroMultipole), *cmm, 0, cmm->nCells, false, MAX_THREADS);
		#endif
		#endif
			// check cluster in cell and calculate the multipole of each cell
			if (bCMM_StorageChain) {
				if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
				else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
			}
			else {
				MTOperate<MM_CMM_JOB<MACROMOL_THREAD_VARS<MMOLECULE>, CMM_CELL3D<BASIC_CLUSTER> > >((void*)(&MM_ClusterAllocate2CMM), mmCMMAllocator, nThreads);
				for (nc = 0; nc < nThreads; nc++) {
					if (!mmCMMAllocator[nc].status) {
						sprintf(errmsg, "loop %d thread %d: buffer of CMM cell leaks when allocating cluster", nloop, nc);
						show_log(errmsg, true);
					}
				}
				if (!combine_AtomAllocatedCMM_array(job_cmm, nThreads, cmm, true, true)) {
					sprintf(errmsg, "loop %d : buffer of CMM cell leaks when combining allocated cluster", nloop);
					show_log(errmsg, true);
				}
			}

		#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		#if _EWALD_SUM == _MultiPole_EWALD_SUM
			calColumbMultipoles<BASIC_CLUSTER>(*cmm, 0, cmm->nCells, false, false, MAX_THREADS);
		#endif
		#endif

			if (iCheckCMM == 0) {
				//time.start();
				//ParallelCheckRelshipInCells(*cmm, 0, cmm->nCells, false, MAX_THREADS); // setup relationship for each atom with CMM
				MMCheckRelshipInCell(mm, 0, mm.nCluster, *cmm, MAX_THREADS);
				//sprintf(::errmsg, "check relationship takes %d", time.elapse());
				//show_infor(errmsg);
			}
			iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;
		
#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "multipole calculation takes %d ms", time.glance());
			show_infor(errmsg);
#endif
			
			MTOperate< MM_INTERACT_JOB >((void*)(&MM_FreeCellInteract), job_interact, nThreads);
			break;
		case 0:
			MTOperate< MM_INTERACT_JOB >((void*)(&MM_FreeCellInteract), job_interact, nThreads); // calculate the interactions between the clusters
			break;
		default: // do not calculate interaction
			break; 
		}
		for (nc = 0; nc < nThreads; nc++) {
			Ep += job_interact[nc].LJ_res.U + job_interact[nc].eres.U;
		}

		MTOperate2< MACROMOL_THREAD_VARS<MMOLECULE>, TorsionVar, InteractRes >((void*)(&MM_TorsionInteract), mmJob, nThreads, torsion_var, mires);
		Etorsion = 0; 
		for (i = 0; i < nThreads; i++) Etorsion += mires[i].U;
		Ep += Etorsion;

		if (::bEext) {
			MOperate3<CLUSTER, float, double>((void*)(&CLUSTER_ExternalEField), mm.cluster, mm.nCluster, ::Eext, MAX_THREADS, Uext);
			Ep += Uext;
		}

		MTOperate< MACROMOL_THREAD_VARS<MMOLECULE> >((void*)(&MM_ClusterTotalForce), mmJob, nThreads);

		//mm.calTorsion(Etorsion);
		//Ep += Etorsion;

		//bHooverDyn = (Overlap ? false : true);
		bHooverDyn = true;

		switch (interact_cal_method) {
		case _CMM:
			CellOperate<BASIC_CLUSTER>((void*)(&SCALE_FORCE_CMM), *cmm, 0, cmm->acell.n, MAX_THREADS);
			break;
		default:
			scale_force_mmolecule(mm);
			break;
		}

		if (bRandomForce) {
			_rand_::init_random_array(randa, mm.nCluster);
			gen_add_random_force(mm, mm.nCluster / 5, randf, fwhm_randf, randa);
		}
		CalNetForce(mm, 0, mm.nCluster); // the net force on molecular mass center

#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "interaction calculation takes %d", time.elapse());
		show_infor(errmsg);
#endif

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

#if _CHECK_TIME_STAMP_ == 1
		time.start();
#endif

		strcpy(msg_show, "\0"); strcpy(msg_log, "\0");
		// to get self-consistent velcoty and acceleration at this step based on leap-frog verlert algorithm
		FixedBase_LeapFrog_SelfConsistent_velocity_acceleration(bHooverDyn, mm, mdpar, nloop, nConvg, vect1, vect2, dv1, dv2, ::local_relax, msg_show, msg_log, n_msglen);
		if (strlen(msg_show) > 1) show_infor(msg_show, true);
		if (strlen(msg_log) > 1) mlog.show(msg_log, true);
		Ek = mm.Ek;
		
		if (local_relax) Verlet_Relax(mdpar, mm, (double)dta_max_relax);
		else Verlet(mdpar, mm); // Verlet algorithm to calculate the ta and tp for next timestep

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

		//show_molecule_struct(mm, refresh_mol_xyz_struct_fname);

		show_loop_infor(nloop, mdpar.mmsave.mdsave->mESave.save_indx(), mdpar.LOOPS - nloop, (float)Ek, (float)Ep);

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mdpar.CalibrateTorsionAngle(mm);

		nloop++;
	}
	if (vect1 != NULL) {delete[] vect1; vect1 = NULL;}
	if (vect2 != NULL) {delete[] vect2; vect2 = NULL;}
	if (dv1 != NULL) {delete[] dv1; dv1 = NULL;}
	if (dv2 != NULL) {delete[] dv2; dv2 = NULL;}
}
