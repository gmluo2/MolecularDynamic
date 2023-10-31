#include "project.h"

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
using namespace _cmm_3d_;

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
#include "complex.h"
#include "fftw3.h"
#include "EwaldSum.h"
#include "spme.h"
#include "spme_2d.h"
#include "spme_interact.h"
#include "spme_interact_2d.h"
extern BSplineFunc<2> bsp;

using namespace _cmm_3d_;

using namespace _EwaldSum_real_;

#include "interface.h"
using namespace _interface_constraint_;

#include "mdvar.h"

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif


using namespace _evatom_;
using namespace _EwaldSum_real_;
using namespace _EwaldSum_2d_;

#include "atom_sr.h"
#include "atom_sr_interact.h"

#include "gpu_vector.h"
#include "gpu_spme.h"
#include "gpu_spme_2d.h"
#include "gpu_interact.h"
#include "gpu_interact_2d.h"

// end of interface constraint
extern char parfname_IConstraint[256];

extern void CompensateNetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 3> *sv3);
extern void CompensateNetForceZ(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 3> *sv3);
extern void NetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob);
extern void ZeroSpatialMomentum(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob);
extern void MDCELL_CalMolTransKineticEnergy(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4> *Ekt);
extern void InterfaceConstraints_CMM(CMM_CELL3D<BASIC_CLUSTER> *cmm, int nThreads, double &Ut);
extern void ApplyZSpring(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob, float *k);
extern void MDCellMassCenter(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob);

extern _bias_force_::BIAS_POINT mdcell_bias;

#ifndef __GPU__
void GPU_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, int md_mode) {
	using namespace _spme_;

	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}
	// we need to make sure the mass center of the cell is at zero
	VECTOR3 mass_center;
	set_free_cell_mass_center(mmcell, mass_center);

	if (md_mode == MULTI_MM_PERIDIC_CELL_MD) {
		molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell
	}

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************
	show_all_molecules();

	if (strlen(parfname_IConstraint) > 0) {
		if (init_clusterIC_DB(2, cmm->xl[2], cmm->xr[2])) {
			sprintf(errmsg, "Interface defined in %s is enabled!", parfname_IConstraint); show_log(errmsg, true);
		}
	}
	if (clusterIC_db.n == 0) { // no constraint is applied
		if (mmcell.bHDynMC_T) {
			sprintf(errmsg, "Interface constraint is disabled! So, constraint on translation along z axis is disabled also!", parfname_IConstraint); show_log(errmsg, true);
			mmcell.bHDynMC_T = false;
		}
	}

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	if (::rcut_Ewd > cmm->xd[0] / 3) ::rcut_Ewd = cmm->xd[0] / 3;
	if (::rcut_Ewd > cmm->xd[1] / 3) ::rcut_Ewd = cmm->xd[1] / 3;
	if (::rcut_Ewd > cmm->xd[2] / 3) ::rcut_Ewd = cmm->xd[2] / 3;

	bool bStorageChain = false;
	InitMD(mmcell, cmm, bStorageChain);


	NVT_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	mdvar.q0 = 0; // a charge for the interactions between induced dipoles. The charge has to be set to 0, important!!

	mdvar.set_virial(::bVirial, true);
	mdvar.spme_var.bEfieldOnly = false; 
	mdvar.spme_var.esv1.bEwald = true; 
	mdvar.spme_var.esv1.init_EwaldSum(cmm->xd[0] * cmm->xd[1] * cmm->xd[2], ::rcut_Ewd);
	mdvar.spme_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	mdvar.spme_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;

	init_BSplineFunc<2>(::bsp, ::indx_bspline);
	mdvar.vspme_q.mMP = 0; // charge only
	mdvar.vspme_mu.mMP = 1; // with dipole
	mdvar.vspme_mu_induced.mMP = 1; // with dipole of course
	if (bPolarize) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		mdvar.vspme_mu.bsp = &(::bsp);
		init_spme(&mmcell, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_mu.xl[0] = cmm->xl[0]; mdvar.vspme_mu.xl[1] = cmm->xl[1]; mdvar.vspme_mu.xl[2] = cmm->xl[2];
		mdvar.vspme_mu.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		mdvar.vspme_mu.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();

		mdvar.vspme_mu_induced.bsp = &(::bsp);
		init_spme_induced(&mmcell, mdvar.vspme_mu_induced, &(mdvar.q0)); // init the viral atoms
		mdvar.vspme_mu_induced.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_mu_induced.xl[0] = cmm->xl[0]; mdvar.vspme_mu_induced.xl[1] = cmm->xl[1]; mdvar.vspme_mu_induced.xl[2] = cmm->xl[2];
		mdvar.vspme_mu_induced.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		mdvar.vspme_mu_induced.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		mdvar.vspme_mu_induced.init_b(); mdvar.vspme_mu_induced.init_C(); mdvar.vspme_mu_induced.init_vars();
	}
	else if (::bDipole) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		mdvar.vspme_mu.bsp = &(::bsp);
		init_spme(&mmcell, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_mu.xl[0] = cmm->xl[0]; mdvar.vspme_mu.xl[1] = cmm->xl[1]; mdvar.vspme_mu.xl[2] = cmm->xl[2];
		mdvar.vspme_mu.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		mdvar.vspme_mu.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();
	}
	else { // charge only
		mmcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		mdvar.vspme_q.bsp = &(::bsp);
		init_spme(&mmcell, mdvar.vspme_q); // init the viral atoms
		mdvar.vspme_q.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_q.xl[0] = cmm->xl[0]; mdvar.vspme_q.xl[1] = cmm->xl[1]; mdvar.vspme_q.xl[2] = cmm->xl[2];
		mdvar.vspme_q.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		mdvar.vspme_q.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		mdvar.vspme_q.init_b(); mdvar.vspme_q.init_C(); mdvar.vspme_q.init_vars();
	}

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
#endif

	int i = 0;
	int nloop = 0, nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0, U_LJ = 0, U_estat = 0;


	double Uext = 0; // external electric field

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	int nMDThreads = MAX_THREADS;

	MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> mdcellHDynThreadVar[MAX_THREADS];
	AssignJob_MDCELL_HDYN<MMOL_MD_CELL>(&mmcell, mdcellHDynThreadVar, nMDThreads, false);

	mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::init_job(mmcell, cmm);
	sprintf(errmsg, "Job # : %d; working with molecules:", MAX_THREADS); show_log(errmsg, true);
	for (i = 0; i < MAX_THREADS; i++) {
		sprintf(errmsg, "%d -- MM %d, SM %d, PM %d", i, mdvar.mdcellJob[i].mmIndx.n, mdvar.mdcellJob[i].smIndx.n, mdvar.mdcellJob[i].pmIndx.n);
		show_log(errmsg, true);
		//for (nm = 0; nm < mdvar.mdcellJob[i].mmIndx.n; nm++) {
		//	sprintf(errmsg, "%d ", mdvar.mdcellJob[i].mmIndx.m[nm]); show_log(errmsg, false);
		//}
		//show_log("", true);
	}

	double djob[MAX_THREADS];
	//MATRIX<3> mjob[MAX_THREADS];
	bool status_job[MAX_THREADS];
	//INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];


	// make sure the spatial momentum of the whole cell is zero
	MDCellMassCenter(mdvar.job_mdcell, sv3job, nMDThreads);
	ZeroSpatialMomentum(mdvar.job_mdcell, sv3job, MAX_THREADS);

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;
	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;
	bool bCheckClusterInCMM = true;
	bool Overlap = false;
	int nConvg = 0;
	bool bSpeedVerlet = true;

	mmcell.Init_MD(); //copy the initial torsion angle, speed to the mdpar
	if (bPolarize) mmcell.init_dipole_hist(); // has to be done after Init_MD();

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	int LOOPS = mmcell.LOOPS, save_indx = 0;
	MD_SAVE *mdsave = mmcell.msave;
	AsyncMDSave_VARS asyncMDSave;

	bool bSaveVirial = false;
	ITER_SAVE virial_save;
	if (::bVirial) {
		virial_save.set_file(mdsave->title, "virial");
		virial_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	ITER_SAVE mu_save;
	bool bSave_mu = (bPolarize || ::bDipole ? true : false);
	if (bSave_mu) {
		mu_save.set_file(mdsave->title, "mu");
		mu_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	//LJ_AsynVARS asynVar_LJ;

	int ncs_mmcell = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) ncs_mmcell += mmcell.mm.m[nm].nCluster;
	ncs_mmcell += mmcell.sm.n + mmcell.pm.n;

	mdvar.spme_var.esv1.iSurface = 0; // 3d
	mdvar.spme_var.esv1.iSurfaceBoundary = 0;//::iSurfaceBoundary; // normal
	mdvar.spme_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);

	_atomic_cmm_::EAtCELL ecell;
	_atomic_cmm_::SR_AtCELL srcell;
	ecell.rmax_c = ::size_cluster / 2; // maximum radius of the cluster
	srcell.rmax_c = ::size_cluster / 2; // maximum radius of the cluster
	ecell.rcut = ::rcut_Ewd;
	srcell.rcut = ::rcut_LJ;

	float cell_corner[3] = {cmm->xl[0], cmm->xl[1], cmm->xl[2]};
	//int ngrid[3] = {cmm->nx, cmm->ny, cmm->nz};
	int ngrid[3] = {int(cmm->xd[0]), int(cmm->xd[1]), int(cmm->xd[2])};
	if (ngrid[0] < 10) ngrid[0] = 10;
	if (ngrid[1] < 10) ngrid[1] = 10;
	if (ngrid[2] < 10) ngrid[2] = 10;
	float cell_cw[3] = {cmm->xd[0] / ngrid[0], cmm->xd[1] / ngrid[1], cmm->xd[2] / ngrid[2]}; // is this right ?
	_atomic_cmm_::init_EAtCell(mmcell, ecell, ngrid[0], ngrid[1], ngrid[2], cell_corner, cell_cw);
	//distribute_atoms(&ecell);
	_atomic_cmm_::update_atomic_charge_dipole(&ecell);
	ecell.esv.cp_EwaldSum_vars(*((_EwaldSum_real_::EwaldSumRealVars0*)&mdvar.spme_var.esv1));
	ecell.esv.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	ecell.bCalE = 1; // calculate force
	ecell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

	_atomic_cmm_::init_SRAtCell(mmcell, srcell, ngrid[0], ngrid[1], ngrid[2], cell_corner, cell_cw);
	//distribute_atoms(&srcell);
	srcell.bCalSR = 1; // calculating short-range L-J or other interaction

	AsynRun<_atomic_cmm_::AsynLJVar> asynRun_LJ;
	_atomic_cmm_::AsynLJVar asynVar_LJ;
	asynVar_LJ.cell = &srcell; asynVar_LJ.nThreads = MAX_THREADS;
	asynVar_LJ.res = &(mdvar.MDVar_LJ::LJ_res);
	asynRun_LJ.set_par(&asynVar_LJ);

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command

		iCheckCMM = 0;

		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			switch(md_mode) {
			case MULTI_MM_PERIDIC_CELL_MD:
				// here we have to check the periodity of the molecule position
				//molecule_check_periodic_cell(mmcell);
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_molecule_PeriodicCell_check), 
					mdvar.job_mdcell, nMDThreads);
				break;
			case MULTI_MM_FREE_CELL_MD:
				set_free_cell_mass_center(mmcell, mass_center);
				break;
			}
		}

		check_atom_rg(mmcell, nMDThreads);

	// re-distribute / check clusters in CMM, use vp of each cluster only
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		
		if (bCheckClusterInCMM) { // update distribution, including coordination
			_atomic_cmm_::distribute_atoms(&ecell);
			_atomic_cmm_::distribute_atoms(&srcell);
			/*
			{
				_atomic_cmm_::check_distribution((_AtCELL*)&ecell);
				_atomic_cmm_::check_distribution((_AtCELL*)&srcell);
			}
			*/
		}
		else { // update coordination only
			_atomic_cmm_::update_atoms(&ecell);
			_atomic_cmm_::update_atoms(&srcell);
		}

	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		show_log(errmsg, true);
	#endif

		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_InitClusterForce), mdvar.job_mdcell, nMDThreads);

		Ep = 0;
		mdvar.MDVar_LJ::LJ_res.reset();
		if (!::bMTForce) {
			// important : calculate LJ interaction first, here cluster force & torque are reset
			_atomic_cmm_::LJ_Interact(&srcell, MAX_THREADS, mdvar.MDVar_LJ::LJ_res);
			U_LJ = mdvar.MDVar_LJ::LJ_res.U;
		}
		else {
			// use asynchrotron function for L-J interaction
		#ifdef __GPU__
			if (srcell.bUseGPU) asynRun_LJ.start_func((void*)(&_atomic_cmm_::cudaLJ));
			else
		#endif
				asynRun_LJ.start_func((void*)(&_atomic_cmm_::cpuLJ));
		}

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalInertiaMoment), mdvar.job_mdcell, nMDThreads);

		MDCellMassCenter(mdvar.job_mdcell, sv3job, nMDThreads);

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!mmcell.eNeutral && ::bRealEwaldSumOnly) {
			if (::bPolarize || ::bDipole) _atomic_cmm_::LocalF(ecell, mdvar.vspme_mu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			else _atomic_cmm_::LocalF(ecell, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
		}
		else if (!::local_interact && !mmcell.eNeutral) {
			ecell.bCalE = 2; // force

			// important : calculate Coulomb secondly, here cluster force & torque are accumulated
			if (::bPolarize) {
				if (nloop > 3) Guess_induced_dipole(mmcell, MAX_THREADS);
				if (::bExcludeInducedDipole) {
					if (!_atomic_cmm_::Polarize_SPME_Interact2(mmcell, ecell, mdvar.vspme_mu, mdvar.vspme_mu_induced, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
						sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
					}
				}
				else if (!_atomic_cmm_::Polarize_SPME_Interact(mmcell, ecell, mdvar.vspme_mu, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
				}
				Backup_induced_dipole_hist(mmcell, MAX_THREADS);
			}
			else if (::bDipole) {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); 
				memcpy(&(ecell.mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				_atomic_cmm_::update_atomic_charge_dipole(&ecell);
				// we have to transfer the surface dipole to GPU

				_atomic_cmm_::SPME_EF(ecell, mdvar.vspme_mu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res, true);
			}
			else {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); 
				memcpy(&(ecell.mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				//_atomic_cmm_::update_atomic_charge_dipole(&cell); // no atomic dipole
				// we have to transfer the surface dipole to GPU

				_atomic_cmm_::SPME_EF(ecell, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			}

			U_estat = mdvar.spme_res.U;
		}
#endif
		
		if (::bMTForce) {
			asynRun_LJ.wait();
			collect_LJ_res(&srcell, MAX_THREADS, mdvar.MDVar_LJ::LJ_res);
			U_LJ = mdvar.MDVar_LJ::LJ_res.U;
		}

		Ep = U_LJ + U_estat;
		if (::bVirial) {
			mdvar.virial.v[0] = mdvar.LJ_res.STxx + mdvar.spme_res.STxx;
			mdvar.virial.v[1] = mdvar.LJ_res.STyy + mdvar.spme_res.STyy;
			mdvar.virial.v[2] = mdvar.LJ_res.STzz + mdvar.spme_res.STzz;
			mdvar.virial.v[3] = mdvar.LJ_res.STtr + mdvar.spme_res.STtr;

			if (!virial_save.one_more()) break;
			bSaveVirial = virial_save.bsave(); 
			//if (virial_save.bsave()) {
			//	*(virial_save.out)<<mdvar.LJ_res.STxx<<"  "<<mdvar.LJ_res.STyy<<"  "<<mdvar.LJ_res.STzz<<"  "<<mdvar.LJ_res.STtr<<endl;
			//	*(virial_save.out)<<mdvar.spme_res.STxx<<"  "<<mdvar.spme_res.STyy<<"  "<<mdvar.spme_res.STzz<<"  "<<mdvar.spme_res.STtr<<endl;
			//}
		}
		
		//virial.v[0] = 0; virial.v[1] = 0; virial.v[2] = 0; virial.v[3] = 0; Ep = 0;

		if (bSave_mu) {
			mu_save.one_more();
			if (mu_save.bsave()) {
				*(mu_save.out) << mmcell.mu_surf.mu.v[0]<<"  "<<mmcell.mu_surf.mu.v[1]<<"  "<<mmcell.mu_surf.mu.v[2]<<endl;
			}
		}

		if (mmcell.mm.n > 0) {
			MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, TorsionVar, InteractRes>(
				(void*)(&MMCELL_MM_Torsion), mdvar.job_mdcell, nMDThreads, mdvar.torsion_var, mdvar.mires);
			Etorsion = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ep += mdvar.mires[i].U; Etorsion += mdvar.mires[i].U;
				if (mdvar.torsion_var[i].bVirial) {
					mdvar.virial.v[0] += mdvar.mires[i].STxx;
					mdvar.virial.v[1] += mdvar.mires[i].STyy;
					mdvar.virial.v[2] += mdvar.mires[i].STzz;
					mdvar.virial.v[3] += mdvar.mires[i].STtr;
				}
			}
			//sprintf(errmsg, "%d: Torsion %f kT", nloop, Etorsion); show_log(errmsg, true);
		}

		if (::bEext) {
			if (mmcell.mm.m != NULL) {
				MOperate3<MMOLECULE, float, double>((void*)(&MM_ExternalEField), mmcell.mm.m, mmcell.mm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.sm.m != NULL) {
				MOperate3<SMOLECULE, float, double>((void*)(&SM_ExternalEField), mmcell.sm.m, mmcell.sm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.pm.m != NULL) {
				MOperate3<PMOLECULE, float, double>((void*)(&PM_ExternalEField), mmcell.pm.m, mmcell.pm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
		}


		//if (mmcell.mm.m != NULL) MOperate<MMOLECULE>(MM_TOTAL_FORCE, mmcell.mm.m, mmcell.mm.n, MAX_THREADS);
		//if (mmcell.sm.m != NULL) MOperate<SMOLECULE>(SM_TOTAL_FORCE, mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
		//if (mmcell.pm.m != NULL) MOperate<PMOLECULE>(PM_TOTAL_FORCE, mmcell.pm.m, mmcell.pm.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ClusterForce), mdvar.job_mdcell, nMDThreads);

		if (::bVirial && bSaveVirial && (mmcell.mm.n > 0 || mmcell.sm.n > 0)) {
			for (i = 0; i < nMDThreads; i++) {
				sv4job[i].v[0] = 0; sv4job[i].v[1] = 0; sv4job[i].v[2] = 0; sv4job[i].v[3] = 0;
			}

			switch (::method_Virial) {
			case CLUSTER_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >(
					(void*)(&MDCELL_VirialCorrect_RigidCluster), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			case MOLECULE_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >(
					(void*)(&MDCELL_VirialCorrect_Molecule), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			}
			
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] -= sv4job[i].v[0];
				mdvar.virial.v[1] -= sv4job[i].v[1];
				mdvar.virial.v[2] -= sv4job[i].v[2];
				mdvar.virial.v[3] -= sv4job[i].v[0] + sv4job[i].v[1] + sv4job[i].v[2];
			}
		}
		
		if (::bVirial) {
			if (virial_save.bsave()) {
				(*virial_save.out)<<nloop<<"  "<<mdvar.virial.v[0] * ::unit_Ek_kT<<"  "<<mdvar.virial.v[1] * ::unit_Ek_kT
					<<"  "<<mdvar.virial.v[2] * ::unit_Ek_kT<<"  "<<mdvar.virial.v[3] * ::unit_Ek_kT;
			}
		}

		// IMPORTANT: the forces from here on needs to calculate the induced torque explicitly
		// randome force and ImplicitSolvent force, and the induced torque, will be calculated explicitly
		// so we do not need to call MDCELL_ClusterForce.

		if (bImplicitSolvent) {
			// neighbors of cluster for implicit solvation 
			//if (bCheckClusterInCMM) CellOperate<BASIC_CLUSTER>((void*)(&CMMCheckImpSolvNeighbors_Cell), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			// solvation energy
			//CellOperate<BASIC_CLUSTER>((void*)(&CMM_ImplicitSolvationForce), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			//for (nm = 0; nm < mmcell.mm.n; nm++) {
			//	mm = mmcell.mm.m + nm;
			//	for (nc = 0; nc < mm->nCluster; nc++) Ep += mm->cluster[nc].E_ImplicitSolvation;
			//}
			MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >(
				(void*)(&MMCELL_CMM_CheckImpSolvRelship), mdvar.job_mdcellCMM, nMDThreads);
			MTOperate1< MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, double>(
				(void*)(&MMCELL_cluster_impsolv_interact), mdvar.job_mdcellCMM, nMDThreads, djob);
			for (i = 0; i < nMDThreads; i++) Ep += djob[i];
		}

		//CellOperate<BASIC_CLUSTER>((void*)(&SCALE_FORCE_CMM), *cmm, 0, cmm->acell.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ScaleClusterForce), mdvar.job_mdcell, nMDThreads);

		if (bRandomForce) gen_random_force(mmcell, mdvar.rand[0], mdvar.rand[1]);

		if (_bias_force_::biasOn) apply_point_bias_force(mmcell, mdcell_bias);

		if (!::local_interact && ::mode_NetForce > 0) {
			NetForce(mdvar.job_mdcell, nMDThreads); // calculate net force and acceleration of the whole cell
			if (::bZSpring) {
				mmcell.fc.v[2] += ::k_zspring[0] * mmcell.rc.v[2] + (mmcell.rc.v[2] > 0 ? 1 : -1) * ::k_zspring[1] * mmcell.rc.v[2] * mmcell.rc.v[2];
				extern void MDCELL_ZVelocity(MDCELL_ThreadJob<MMOL_MD_CELL> *job, int nJob);
				MDCELL_ZVelocity(mdvar.job_mdcell, nMDThreads); // calculate z-velocity of cell
				mmcell.fc.v[2] += ::ita_mc * mmcell.Vmc.v[2] * mmcell.M;
			}
			if (::mode_NetForce == 1) {
				mmcell.ac.v[0] = mmcell.fc.v[0] / mmcell.M;
				mmcell.ac.v[1] = mmcell.fc.v[1] / mmcell.M;
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}
			else if (::mode_NetForce == 2) {
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}
			if (::mode_NetForce == 1) CompensateNetForce(mdvar.job_mdcell, nMDThreads, sv3job); // compensate the net force on all specises in the cell
			else if (::mode_NetForce == 2) CompensateNetForceZ(mdvar.job_mdcell, nMDThreads, sv3job); // compensate the net force on all specises in the cell
			/*
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] += sv6job[i].v[0];
				mdvar.virial.v[1] += sv6job[i].v[1];
				mdvar.virial.v[2] += sv6job[i].v[2];
				mdvar.virial.v[3] += (sv6job[i].v[0] + sv6job[i].v[1] + sv6job[i].v[2]);
			}
			*/
		}

		if (clusterIC_db.n > 0) {
			// the force from interface, and the torque
			// however, the force will not contribute to the virial
			InterfaceConstraints_CMM(cmm, nMDThreads, djob[0]); Ep += djob[0];
		}

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		// self-consistent velocity and acceleration, Nose-Hoover Chain
		bSpeedVerlet = (nloop == 0 ? false : true);
		if (!LeapFrogVerlet_SelfConsistent_SpeedAccel(&mmcell, mdvar.job_mdcell, mdcellHDynThreadVar, bSpeedVerlet, nMDThreads)) {
			sprintf(errmsg, "failure to get consistent velocity and acceleration at loop %d", nloop);
			show_log(errmsg, true);
		}

		Ep_t = Ep; Ek_t = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) Ek_t += mmcell.mm.m[nm].Ek;
		for (nm = 0; nm < mmcell.sm.n; nm++) Ek_t += mmcell.sm.m[nm].Ek;
		for (nm = 0; nm < mmcell.pm.n; nm++) Ek_t += mmcell.pm.m[nm].Ek;

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);
		
		// save md asynchrotronly
		//if (mmcell.msave->bBuffered) asyncMDSave.WaitUntilThreadOver();
		//mmcell.MDSave(nloop, (float)Ep_t, false, true);
		//AsyncMDSave(mmcell, asyncMDSave);

		mmcell.MDSave(nloop, (float)Ep_t, false);

		//if (::bPolarize && mdsave->mKinSave.bsave()) log_dipole(&mmcell); // write the dipole of cluster into log file

		if (::bVirial && bSaveVirial && ::method_Virial == CLUSTER_VIRIAL) {
			if (mmcell.mm.n > 0 && mmcell.bHingeConstraint) {
				MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, HingeConstraintVar, InteractRes>(
					(void*)(&MDCELL_ClusterVirialCorrect_HingeConstraint), mdvar.job_mdcell, nMDThreads, 
					mdvar.hinge_constraint_var, mdvar.mires);
				for (i = 0; i < nMDThreads; i++) {
					mdvar.virial.v[0] += mdvar.mires[i].STxx;
					mdvar.virial.v[1] += mdvar.mires[i].STyy;
					mdvar.virial.v[2] += mdvar.mires[i].STzz;
					mdvar.virial.v[3] += mdvar.mires[i].STtr;
				}
			}

			MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_CalMolTransKineticEnergy), mdvar.job_mdcell, nMDThreads, sv4job);
			Ek_t = 0;
			for (i = 0; i < nMDThreads; i++) Ek_t += sv4job[i].v[3];
			mdvar.Pint = (mdvar.virial.v[3] + Ek_t * 2 / unit_Ek_kT) / (::fUnit_atm_u * mmcell.h[0] * mmcell.h[1] * mmcell.h[2] * 8 * 3);
			if (bSaveVirial) {
				(*virial_save.out)<<"          "<<mdvar.virial.v[3] * ::unit_Ek_kT<<"  "<<Ek_t<<"  "<<mdvar.Pint<<endl;
			}
		}

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_Verlet), mdvar.job_mdcell, nMDThreads);
		// accept current NHC vars and prepair for next step
		for (i = 0; i < mmcell.hdyn.n; i++) mmcell.hdyn.m[i].nhc->accept_nextstep();
		if (mmcell.bHDynMC_T) mmcell.mcHDyn.nhc->accept_nextstep();

		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			if (mmcell.mm.n > 2) {
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CheckRotAngle), mdvar.job_mdcell, nMDThreads);
			}
			else {
				for (nm = 0; nm < mmcell.mm.n; nm++) {
					mm = mmcell.mm.m + nm;
					check_angle(*mm, mmcell.lfmd.m[nm]);
				}
			}
		}
		iCheckAngle++; if (iCheckAngle >= nCheckAngle) iCheckAngle = 0;

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mmcell.CalibrateTorsionAngle();

		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "checking procedure takes %d ms", time.glance());
			show_infor(errmsg);
		#endif

		if (::bVirial) virial_save.save_over();

		nloop++;
	}
	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
	mdsave->flush_buf();
	if (::bVirial) virial_save.reset();
	if (bSave_mu) mu_save.save_over();
}

#else // define __GPU__


void GPU_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, int md_mode) {
	using namespace _spme_;

	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}
	// we need to make sure the mass center of the cell is at zero
	VECTOR3 mass_center;
	set_free_cell_mass_center(mmcell, mass_center);

	if (md_mode == MULTI_MM_PERIDIC_CELL_MD) {
		molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell
	}

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************
	show_all_molecules();

	if (strlen(parfname_IConstraint) > 0) {
		if (init_clusterIC_DB(2, cmm->xl[2], cmm->xr[2])) {
			sprintf(errmsg, "Interface defined in %s is enabled!", parfname_IConstraint); show_log(errmsg, true);
			apply_slab_barrier_constraint(clusterIC_db);
		}
	}
	if (clusterIC_db.n == 0) { // no constraint is applied
		if (mmcell.bHDynMC_T) {
			sprintf(errmsg, "Interface constraint is disabled! So, constraint on translation along z axis is disabled also!", parfname_IConstraint); show_log(errmsg, true);
			mmcell.bHDynMC_T = false;
		}
	}
	if (!::bZSpring && ::mode_NetForce == 0 && mmcell.bHDynMC_T) {
		sprintf(errmsg, "although constraint on translation along z axis is enabled in MD.INI, disabled it!"); show_log(errmsg, true);
		mmcell.bHDynMC_T = false;
	}

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	if (::rcut_Ewd > cmm->xd[0] / 3) ::rcut_Ewd = cmm->xd[0] / 3;
	if (::rcut_Ewd > cmm->xd[1] / 3) ::rcut_Ewd = cmm->xd[1] / 3;
	if (::rcut_Ewd > cmm->xd[2] / 3) ::rcut_Ewd = cmm->xd[2] / 3;

	bool bStorageChain = false;
	InitMD(mmcell, cmm, bStorageChain);


	NVT_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	mdvar.q0 = 0; // a charge for the interactions between induced dipoles. The charge has to be set to 0, important!!

	mdvar.set_virial(::bVirial, true);
	mdvar.spme_var.bEfieldOnly = false; 
	mdvar.spme_var.esv1.bEwald = true; 
	mdvar.spme_var.esv1.init_EwaldSum(cmm->xd[0] * cmm->xd[1] * cmm->xd[2], ::rcut_Ewd);
	mdvar.spme_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	mdvar.spme_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;

	init_BSplineFunc<2>(::bsp, ::indx_bspline);
	mdvar.vspme_q.mMP = 0; // charge only
	mdvar.vspme_mu.mMP = 1; // with dipole
	mdvar.vspme_mu_induced.mMP = 1; // with dipole of course
	if (bPolarize) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		mdvar.vspme_mu.bsp = &(::bsp);
		init_spme(&mmcell, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_mu.xl[0] = cmm->xl[0]; mdvar.vspme_mu.xl[1] = cmm->xl[1]; mdvar.vspme_mu.xl[2] = cmm->xl[2];
		mdvar.vspme_mu.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		mdvar.vspme_mu.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();

		mdvar.vspme_mu_induced.bsp = &(::bsp);
		init_spme_induced(&mmcell, mdvar.vspme_mu_induced, &(mdvar.q0)); // init the viral atoms
		mdvar.vspme_mu_induced.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_mu_induced.xl[0] = cmm->xl[0]; mdvar.vspme_mu_induced.xl[1] = cmm->xl[1]; mdvar.vspme_mu_induced.xl[2] = cmm->xl[2];
		mdvar.vspme_mu_induced.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		mdvar.vspme_mu_induced.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		mdvar.vspme_mu_induced.init_b(); mdvar.vspme_mu_induced.init_C(); mdvar.vspme_mu_induced.init_vars();
	}
	else if (::bDipole) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		mdvar.vspme_mu.bsp = &(::bsp);
		init_spme(&mmcell, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_mu.xl[0] = cmm->xl[0]; mdvar.vspme_mu.xl[1] = cmm->xl[1]; mdvar.vspme_mu.xl[2] = cmm->xl[2];
		mdvar.vspme_mu.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		mdvar.vspme_mu.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();
	}
	else { // charge only
		mmcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		mdvar.vspme_q.bsp = &(::bsp);
		init_spme(&mmcell, mdvar.vspme_q); // init the viral atoms
		mdvar.vspme_q.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_q.xl[0] = cmm->xl[0]; mdvar.vspme_q.xl[1] = cmm->xl[1]; mdvar.vspme_q.xl[2] = cmm->xl[2];
		mdvar.vspme_q.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		mdvar.vspme_q.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		mdvar.vspme_q.init_b(); mdvar.vspme_q.init_C(); mdvar.vspme_q.init_vars();
	}

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
#endif

	int i = 0;
	long nloop = 0;
	int nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0, U_LJ = 0, U_estat = 0;


	double Uext = 0; // external electric field

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	int nMDThreads = MAX_THREADS;

	MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> mdcellHDynThreadVar[MAX_THREADS];
	AssignJob_MDCELL_HDYN<MMOL_MD_CELL>(&mmcell, mdcellHDynThreadVar, nMDThreads, false);

	mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::init_job(mmcell, cmm);
	sprintf(errmsg, "Job # : %d; working with molecules:", MAX_THREADS); show_log(errmsg, true);
	for (i = 0; i < MAX_THREADS; i++) {
		sprintf(errmsg, "%d -- MM %d, SM %d, PM %d", i, mdvar.mdcellJob[i].mmIndx.n, mdvar.mdcellJob[i].smIndx.n, mdvar.mdcellJob[i].pmIndx.n);
		show_log(errmsg, true);
		//for (nm = 0; nm < mdvar.mdcellJob[i].mmIndx.n; nm++) {
		//	sprintf(errmsg, "%d ", mdvar.mdcellJob[i].mmIndx.m[nm]); show_log(errmsg, false);
		//}
		//show_log("", true);
	}

	double djob[MAX_THREADS];
	MATRIX<3> mjob[MAX_THREADS];
	bool status_job[MAX_THREADS];
	//INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];


	// make sure the spatial momentum of the whole cell is zero
	if (mmcell.bHDynMC_T) 	MDCELL_mass_center(mmcell);
	ZeroSpatialMomentum(mdvar.job_mdcell, sv3job, MAX_THREADS);

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;
	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;
	bool bCheckClusterInCMM = true;
	bool Overlap = false;
	int nConvg = 0;
	bool bSpeedVerlet = true;

	mmcell.Init_MD(); //copy the initial torsion angle, speed to the mdpar
	if (bPolarize) mmcell.init_dipole_hist(); // has to be done after Init_MD();

#if _CHECK_TIME_STAMP_ == 1
	TIME time; int dt;
#endif

	long LOOPS = mmcell.LOOPS;
	int save_indx = 0;
	MD_SAVE *mdsave = mmcell.msave;
	AsyncMDSave_VARS asyncMDSave;

	bool bSaveVirial = false;
	ITER_SAVE virial_save;
	if (::bVirial) {
		virial_save.set_file(mdsave->title, "virial");
		virial_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	ITER_SAVE mu_save;
	bool bSave_mu = (bPolarize || ::bDipole ? true : false);
	if (bSave_mu) {
		mu_save.set_file(mdsave->title, "mu");
		mu_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	int ncs_mmcell = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) ncs_mmcell += mmcell.mm.m[nm].nCluster;
	ncs_mmcell += mmcell.sm.n + mmcell.pm.n;

	mdvar.spme_var.esv1.iSurface = 0; // 3d
	mdvar.spme_var.esv1.iSurfaceBoundary = 0;//::iSurfaceBoundary; // normal
	mdvar.spme_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);

	
	
	// ecell, srcell, spme_gpu ...
	if (cudaSetDevice(_CudaDevice_) != cudaSuccess) {
		show_infor("GPU is not enabled!", true); return;
	}
	else show_log("GPU is used!", true);
	cudaDeviceReset();

	_atomic_cmm_::EAtCELL ecell;
	_atomic_cmm_::SR_AtCELL srcell;
	ecell.rmax_c = (::size_cluster / 2 > 1 ? 1 : ::size_cluster / 2) ; // maximum radius of the cluster
	srcell.rmax_c = (::size_cluster / 2 > 1 ? 1 : ::size_cluster / 2); // maximum radius of the cluster
	ecell.rcut = ::rcut_Ewd;
	srcell.rcut = ::rcut_LJ;
	srcell.bCalUseGPU = false; // use GPU for L-J interaction calculation
	ecell.bCalUseGPU = bGPURealEwaldSum; // use CPU for real part Ewald-Sum interaction calculation

	if (srcell.bCalUseGPU) show_log("using GPU for L-J interaction.", true);
	else show_log("using CPU for L-J interaction.", true);

	if (ecell.bCalUseGPU) show_log("using GPU for real part Ewald-Sum.", true);
	else show_log("using CPU for real part Ewald-Sum.", true);

	_atomic_cmm_::AsynLJVar cudaLJVar;
	cudaLJVar.cell = &srcell; cudaLJVar.nThreads = MAX_THREADS;
	cudaLJVar.res = &(mdvar.MDVar_LJ::LJ_res);
	AsynRun<_atomic_cmm_::AsynLJVar> asynRun_LJ;
	asynRun_LJ.set_par(&cudaLJVar);

	float cell_corner[3] = {cmm->xl[0], cmm->xl[1], cmm->xl[2]};
	int ngrid[3] = {int(cmm->xd[0] + 0.5), int(cmm->xd[1] + 0.5), int(cmm->xd[2] + 0.5)};
	if (ngrid[0] < 10) ngrid[0] = 10;
	if (ngrid[1] < 10) ngrid[1] = 10;
	if (ngrid[2] < 10) ngrid[2] = 10;
	float cell_cw[3] = {cmm->xd[0] / ngrid[0], cmm->xd[1] / ngrid[1], cmm->xd[2] / ngrid[2]}; // is this right ?
	_atomic_cmm_::init_EAtCell(mmcell, ecell, ngrid[0], ngrid[1], ngrid[2], cell_corner, cell_cw);
	//distribute_atoms(&ecell);
	//_atomic_cmm_::update_atomic_charge_dipole(&ecell);
	ecell.esv.cp_EwaldSum_vars(*((_EwaldSum_real_::EwaldSumRealVars0*)&mdvar.spme_var.esv1));
	ecell.esv.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	ecell.bCalE = 1; // calculate force
	ecell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

	_atomic_cmm_::init_SRAtCell(mmcell, srcell, ngrid[0], ngrid[1], ngrid[2], cell_corner, cell_cw);
	//distribute_atoms(&srcell);
	srcell.bCalSR = 1; // calculating short-range L-J or other interaction 

	_gpu_md_::GPU_SPME_3d spme_gpu;
	spme_gpu.atom.set_array(mmcell.atom.n); _gpu_md_::init_atom(mmcell.atom, spme_gpu.atom);
	spme_gpu.init_memory(mmcell.atom.n); 
	if (!_gpu_md_::setup_gpu_spme(&spme_gpu, &::bsp, cmm->xl, cmm->xd, ::dim_spme, mdvar.spme_var.esv1.kappa, mmcell.atom.n, (bPolarize ? 1: 0))) {
		show_msg("GPU failure to initialize the memory for SPME", true); return;
	}

	show_log("using GPU for SPME, image part of Ewald-Sum.", true);


	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command

		iCheckCMM = 0;

		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			switch(md_mode) {
			case MULTI_MM_PERIDIC_CELL_MD:
				// here we have to check the periodity of the molecule position
				//molecule_check_periodic_cell(mmcell);
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_molecule_PeriodicCell_check), 
					mdvar.job_mdcell, nMDThreads);
				break;
			case MULTI_MM_FREE_CELL_MD:
				set_free_cell_mass_center(mmcell, mass_center);
				break;
			}
		}

		check_atom_rg(mmcell, nMDThreads);
		distribute_atoms(&srcell);
		distribute_atoms(&ecell);
		update_atomic_charge_dipole(&ecell);
		memcpy(spme_gpu.r.m, ecell.r.m, ecell.r.n * sizeof(double));
		//_gpu_md_::update_atomic_coordinate(&spme_gpu, (MAX_THREADS < 4 ? MAX_THREADS : 4));

		if (ecell.gpu_cell.srnb_dev == NULL && ecell.bCalUseGPU) {
			if (!ecell.check_neighbor_memory(true)) return;
			gpu_reset_esr_neighbors(&ecell);
			check_neighbors(&(ecell.gpu_cell), ecell.gpu_cell_dev);
		}


		if (nloop == 0) {
			memcpy(spme_gpu.q.m, ecell.c.m, ecell.c.n * sizeof(double));
			_gpu_md_::_spme_::update_q_device_spme(&(spme_gpu.vspme_host));
		}

		memcpy(spme_gpu.mu.m, ecell.mu.m, ecell.mu.n * sizeof(double));
		_gpu_md_::_spme_::update_dipole_device_spme(&(spme_gpu.vspme_host));

		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_InitClusterForce), mdvar.job_mdcell, nMDThreads);


		Ep = 0;
		// important : calculate LJ interaction first, here cluster force & torque are reset
		//_atomic_cmm_::LJ_Interact(&srcell, MAX_THREADS, mdvar.MDVar_LJ::LJ_res);
		//U_LJ = mdvar.MDVar_LJ::LJ_res.U;

		// use asynchrotron function for L-J interaction
		mdvar.MDVar_LJ::LJ_res.reset();
		if (srcell.bCalUseGPU) asynRun_LJ.start_func((void*)(&_atomic_cmm_::cudaLJ));
		else asynRun_LJ.start_func((void*)(&_atomic_cmm_::cpuLJ));

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalInertiaMoment), mdvar.job_mdcell, nMDThreads);
		
		MDCellMassCenter(mdvar.job_mdcell, sv3job, nMDThreads);

#if _CHECK_TIME_STAMP_ == 1
		time.start();
#endif

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		// calculat the total dipole of the whole cell
		cal_dipole(mmcell, 1, mmcell.mu_surf);
		ecell.mu_surf = mmcell.mu_surf; cmm->mu_surf = mmcell.mu_surf;

		ecell.bCalE = 2; // 1 -- Efield, 2 -- force & U
		ecell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

		//bool bstatus = _gpu_dstruct_::_gpu_util_::check_gpu();

		if (!mmcell.eNeutral && ::bRealEwaldSumOnly) {
			if (::bPolarize || ::bDipole) _atomic_cmm_::LocalF(ecell, mdvar.vspme_mu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			else _atomic_cmm_::LocalF(ecell, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
		}
		else if (!::local_interact && !mmcell.eNeutral) {
			if (bPolarize) {
				if (::bExcludeInducedDipole) {
					if (!_atomic_cmm_::Polarize_SPME_Interact2(mmcell, ecell, spme_gpu, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
						sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
					}
				}
				else if (!_atomic_cmm_::Polarize_SPME_Interact(mmcell, ecell, spme_gpu, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
				}
			}
			else _atomic_cmm_::SPME_EF(ecell, spme_gpu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);


			U_estat = mdvar.spme_res.U;
		}

#endif

#if _CHECK_TIME_STAMP_ == 1
		dt = time.elapse(); sprintf(errmsg, "interaction calculation takes: %d ms", dt); show_log(errmsg, true); time.start();
#endif
		// wait for GPU-LJ finishing
		asynRun_LJ.wait();
		collect_LJ_res(&srcell, MAX_THREADS, mdvar.MDVar_LJ::LJ_res);
		U_LJ = mdvar.MDVar_LJ::LJ_res.U;

#if _CHECK_TIME_STAMP_ == 1
		dt = time.elapse(); sprintf(errmsg, "waiting for LJ interaction takes: %d ms", dt); show_log(errmsg, true); time.start();
#endif

		Ep = U_LJ + U_estat;
		if (::bVirial) {
			mdvar.virial.v[0] = mdvar.LJ_res.STxx + mdvar.spme_res.STxx;
			mdvar.virial.v[1] = mdvar.LJ_res.STyy + mdvar.spme_res.STyy;
			mdvar.virial.v[2] = mdvar.LJ_res.STzz + mdvar.spme_res.STzz;
			mdvar.virial.v[3] = mdvar.LJ_res.STtr + mdvar.spme_res.STtr;

			if (!virial_save.one_more()) break;
			bSaveVirial = virial_save.bsave(); 
			//if (virial_save.bsave()) {
			//	*(virial_save.out)<<mdvar.LJ_res.STxx<<"  "<<mdvar.LJ_res.STyy<<"  "<<mdvar.LJ_res.STzz<<"  "<<mdvar.LJ_res.STtr<<endl;
			//	*(virial_save.out)<<mdvar.spme_res.STxx<<"  "<<mdvar.spme_res.STyy<<"  "<<mdvar.spme_res.STzz<<"  "<<mdvar.spme_res.STtr<<endl;
			//}
		}
		
		//virial.v[0] = 0; virial.v[1] = 0; virial.v[2] = 0; virial.v[3] = 0; Ep = 0;

		if (bSave_mu) {
			mu_save.one_more();
			if (mu_save.bsave()) {
				*(mu_save.out) << mmcell.mu_surf.mu.v[0]<<"  "<<mmcell.mu_surf.mu.v[1]<<"  "<<mmcell.mu_surf.mu.v[2]<<endl;
			}
		}

		if (mmcell.mm.n > 0) {
			MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, TorsionVar, InteractRes>(
				(void*)(&MMCELL_MM_Torsion), mdvar.job_mdcell, nMDThreads, mdvar.torsion_var, mdvar.mires);
			Etorsion = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ep += mdvar.mires[i].U; Etorsion += mdvar.mires[i].U;
				if (mdvar.torsion_var[i].bVirial) {
					mdvar.virial.v[0] += mdvar.mires[i].STxx;
					mdvar.virial.v[1] += mdvar.mires[i].STyy;
					mdvar.virial.v[2] += mdvar.mires[i].STzz;
					mdvar.virial.v[3] += mdvar.mires[i].STtr;
				}
			}
			//sprintf(errmsg, "%d: Torsion %f kT", nloop, Etorsion); show_log(errmsg, true);
		}

		if (::bEext) {
			if (mmcell.mm.m != NULL) {
				MOperate3<MMOLECULE, float, double>((void*)(&MM_ExternalEField), mmcell.mm.m, mmcell.mm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.sm.m != NULL) {
				MOperate3<SMOLECULE, float, double>((void*)(&SM_ExternalEField), mmcell.sm.m, mmcell.sm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
			if (mmcell.pm.m != NULL) {
				MOperate3<PMOLECULE, float, double>((void*)(&PM_ExternalEField), mmcell.pm.m, mmcell.pm.n, ::Eext, MAX_THREADS, Uext);
				Ep += Uext;
			}
		}

		if (_interface_constraint_::clusterIC_db.n > 0) {
			// the force from interface, and the torque
			// however, the force will not contribute to the virial
			//InterfaceConstraints_CMM(cmm, nMDThreads, djob[0]); Ep += djob[0];

			// constraint force will be added to atom with this function
			if (srcell.bUseGPU) cudaMemcpy(srcell.cmm.m, srcell.gpu_cell.cmm_dev, srcell.cmm.n * sizeof(int), cudaMemcpyDeviceToHost);
			_atomic_cmm_::InterfaceConstraints_AtCell((_atomic_cmm_::_AtCELL*)(&srcell), nMDThreads, djob[0]); //Ep += djob[0];
		}

		//if (mmcell.mm.m != NULL) MOperate<MMOLECULE>(MM_TOTAL_FORCE, mmcell.mm.m, mmcell.mm.n, MAX_THREADS);
		//if (mmcell.sm.m != NULL) MOperate<SMOLECULE>(SM_TOTAL_FORCE, mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
		//if (mmcell.pm.m != NULL) MOperate<PMOLECULE>(PM_TOTAL_FORCE, mmcell.pm.m, mmcell.pm.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ClusterForce), mdvar.job_mdcell, nMDThreads);

		if (::bVirial && bSaveVirial && (mmcell.mm.n > 0 || mmcell.sm.n > 0)) {
			for (i = 0; i < nMDThreads; i++) {
				sv4job[i].v[0] = 0; sv4job[i].v[1] = 0; sv4job[i].v[2] = 0; sv4job[i].v[3] = 0;
			}
			switch (::method_Virial) {
			case CLUSTER_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >(
					(void*)(&MDCELL_VirialCorrect_RigidCluster), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			case MOLECULE_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >(
					(void*)(&MDCELL_VirialCorrect_Molecule), mdvar.job_mdcell, nMDThreads, sv4job);
				break;
			}
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] -= sv4job[i].v[0];
				mdvar.virial.v[1] -= sv4job[i].v[1];
				mdvar.virial.v[2] -= sv4job[i].v[2];
				mdvar.virial.v[3] -= sv4job[i].v[0] + sv4job[i].v[1] + sv4job[i].v[2];
			}
		}
		
		if (::bVirial) {
			if (virial_save.bsave()) {
				(*virial_save.out)<<nloop<<"  "<<mdvar.virial.v[0] * ::unit_Ek_kT<<"  "<<mdvar.virial.v[1] * ::unit_Ek_kT
					<<"  "<<mdvar.virial.v[2] * ::unit_Ek_kT<<"  "<<mdvar.virial.v[3] * ::unit_Ek_kT;
			}
		}

		// IMPORTANT: the forces from here on needs to calculate the induced torque explicitly
		// randome force and ImplicitSolvent force, and the induced torque, will be calculated explicitly
		// so we do not need to call MDCELL_ClusterForce.

		if (bImplicitSolvent) {
			// neighbors of cluster for implicit solvation 
			//if (bCheckClusterInCMM) CellOperate<BASIC_CLUSTER>((void*)(&CMMCheckImpSolvNeighbors_Cell), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			// solvation energy
			//CellOperate<BASIC_CLUSTER>((void*)(&CMM_ImplicitSolvationForce), *cmm, 0, cmm->bcell.nx, MAX_THREADS);
			//for (nm = 0; nm < mmcell.mm.n; nm++) {
			//	mm = mmcell.mm.m + nm;
			//	for (nc = 0; nc < mm->nCluster; nc++) Ep += mm->cluster[nc].E_ImplicitSolvation;
			//}
			MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >(
				(void*)(&MMCELL_CMM_CheckImpSolvRelship), mdvar.job_mdcellCMM, nMDThreads);
			MTOperate1< MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, double>(
				(void*)(&MMCELL_cluster_impsolv_interact), mdvar.job_mdcellCMM, nMDThreads, djob);
			for (i = 0; i < nMDThreads; i++) Ep += djob[i];
		}

		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ScaleClusterForce), mdvar.job_mdcell, nMDThreads);

		if (bRandomForce) gen_random_force(mmcell, mdvar.rand[0], mdvar.rand[1]);

		if (_bias_force_::biasOn) apply_point_bias_force(mmcell, mdcell_bias);

		if (!::local_interact && ::mode_NetForce > 0) {
			NetForce(mdvar.job_mdcell, nMDThreads); // calculate net force and acceleration of the whole cell
			if (::bZSpring) {
				mmcell.fc.v[2] += ::k_zspring[0] * mmcell.rc.v[2] + (mmcell.rc.v[2] > 0 ? 1 : -1) * ::k_zspring[1] * mmcell.rc.v[2] * mmcell.rc.v[2];
				extern void MDCELL_ZVelocity(MDCELL_ThreadJob<MMOL_MD_CELL> *job, int nJob);
				MDCELL_ZVelocity(mdvar.job_mdcell, nMDThreads); // calculate z-velocity of cell
				mmcell.fc.v[2] += ::ita_mc * mmcell.Vmc.v[2] * mmcell.M;
			}
			if (::mode_NetForce == 1) {
				mmcell.ac.v[0] = mmcell.fc.v[0] / mmcell.M;
				mmcell.ac.v[1] = mmcell.fc.v[1] / mmcell.M;
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}
			else if (::mode_NetForce == 2) {
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}

			if (::mode_NetForce == 1) CompensateNetForce(mdvar.job_mdcell, nMDThreads, sv3job); // compensate the net force on all specises in the cell
			else if (::mode_NetForce == 2) CompensateNetForceZ(mdvar.job_mdcell, nMDThreads, sv3job); // compensate the net force on all specises in the cell
			/*
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] += sv3job[i].v[0];
				mdvar.virial.v[1] += sv3job[i].v[1];
				mdvar.virial.v[2] += sv3job[i].v[2];
				mdvar.virial.v[3] += (sv3job[i].v[0] + sv3job[i].v[1] + sv3job[i].v[2]);
			}
			*/
		}

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		// self-consistent velocity and acceleration, Nose-Hoover Chain
		bSpeedVerlet = (nloop == 0 ? false : true);
		if (!LeapFrogVerlet_SelfConsistent_SpeedAccel(&mmcell, mdvar.job_mdcell, mdcellHDynThreadVar, bSpeedVerlet, nMDThreads)) {
			sprintf(errmsg, "failure to get consistent velocity and acceleration at loop %d", nloop);
			show_log(errmsg, true);
		}

		Ep_t = Ep; Ek_t = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) Ek_t += mmcell.mm.m[nm].Ek;
		for (nm = 0; nm < mmcell.sm.n; nm++) Ek_t += mmcell.sm.m[nm].Ek;
		for (nm = 0; nm < mmcell.pm.n; nm++) Ek_t += mmcell.pm.m[nm].Ek;

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);
		
		// save md asynchrotronly
		//if (mmcell.msave->bBuffered) asyncMDSave.WaitUntilThreadOver();
		//mmcell.MDSave(nloop, (float)Ep_t, false, true);
		//AsyncMDSave(mmcell, asyncMDSave);

		mmcell.MDSave(nloop, (float)Ep_t, false);

		//if (::bPolarize && mdsave->mKinSave.bsave()) log_dipole(&mmcell); // write the dipole of cluster into log file

		if (::bVirial && bSaveVirial && ::method_Virial == CLUSTER_VIRIAL) {
			if (mmcell.mm.n > 0 && mmcell.bHingeConstraint) {
				MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, HingeConstraintVar, InteractRes>(
					(void*)(&MDCELL_ClusterVirialCorrect_HingeConstraint), mdvar.job_mdcell, nMDThreads, 
					mdvar.hinge_constraint_var, mdvar.mires);
				for (i = 0; i < nMDThreads; i++) {
					mdvar.virial.v[0] += mdvar.mires[i].STxx;
					mdvar.virial.v[1] += mdvar.mires[i].STyy;
					mdvar.virial.v[2] += mdvar.mires[i].STzz;
					mdvar.virial.v[3] += mdvar.mires[i].STtr;
				}
			}

			MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_CalMolTransKineticEnergy), mdvar.job_mdcell, nMDThreads, sv4job);
			Ek_t = 0;
			for (i = 0; i < nMDThreads; i++) Ek_t += sv4job[i].v[3];
			mdvar.Pint = (mdvar.virial.v[3] + Ek_t * 2 / unit_Ek_kT) / (::fUnit_atm_u * mmcell.h[0] * mmcell.h[1] * mmcell.h[2] * 8 * 3);
			if (bSaveVirial) {
				(*virial_save.out)<<"          "<<mdvar.virial.v[3] * ::unit_Ek_kT<<"  "<<Ek_t<<"  "<<mdvar.Pint<<"    "<<mmcell.rc.v[2]<<endl;
			}
		}

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_Verlet), mdvar.job_mdcell, nMDThreads);
		// accept current NHC vars and prepair for next step
		for (i = 0; i < mmcell.hdyn.n; i++) mmcell.hdyn.m[i].nhc->accept_nextstep();
		if (mmcell.bHDynMC_T) mmcell.mcHDyn.nhc->accept_nextstep();

		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			if (mmcell.mm.n > 2) {
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CheckRotAngle), mdvar.job_mdcell, nMDThreads);
			}
			else {
				for (nm = 0; nm < mmcell.mm.n; nm++) {
					mm = mmcell.mm.m + nm;
					check_angle(*mm, mmcell.lfmd.m[nm]);
				}
			}
		}
		iCheckAngle++; if (iCheckAngle >= nCheckAngle) iCheckAngle = 0;

		iCalibrateTorsionAngle++; if (iCalibrateTorsionAngle == nCalibrateTorsionAngle) iCalibrateTorsionAngle = 0;
		if (iCalibrateTorsionAngle == 0) mmcell.CalibrateTorsionAngle();

		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "checking procedure takes %d ms", time.glance());
			show_infor(errmsg);
		#endif

		if (::bVirial) virial_save.save_over();

		nloop++;
	}
	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
	mdsave->flush_buf();
	if (::bVirial) virial_save.reset();
	if (bSave_mu) mu_save.save_over();

	_gpu_md_::_spme_::release_device_spme(&(spme_gpu.vspme_host), spme_gpu.vspme_dev);
	spme_gpu.vspme_dev = NULL;
	spme_gpu.status.release(); spme_gpu.E.release(); spme_gpu.F.release(); spme_gpu.mu.release(); spme_gpu.q.release();
	spme_gpu.release_plans();

	ecell.gpu_cell.release(); ecell.release_device(); 
	srcell.gpu_cell.release(); srcell.release_device();
	cudaDeviceReset();
}


#endif  //__GPU__

