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
#include "CMM_2d.h"
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
#include "spme_interact.h"
#include "spme_2d.h"
#include "spme_interact_2d.h"

using namespace _evatom_;
using namespace _EwaldSum_real_;
using namespace _EwaldSum_2d_;

#include "mdvar.h"
/*
#include "slab_md.h"
*/
#include "md_module.h"

#include "atom_sr.h"
#include "atom_sr_interact.h"

#include "gpu_vector.h"
#include "gpu_spme.h"
#include "gpu_spme_2d.h"
#include "gpu_interact.h"
#include "gpu_interact_2d.h"

extern BSplineFunc<2> bsp;

#include "interface.h"
using namespace _interface_constraint_;

extern _bias_force_::BIAS_POINT mdcell_bias;

//************************************************************************************
//    3d-periodical cell
//************************************************************************************

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
extern ITER_SAVE tsave;
#endif

void init_spme_3d(MMOL_MD_CELL &mmcell, NVT_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> &mdvar) {
	mdvar.q0 = 0; // a charge for the interactions between induced dipoles. The charge has to be set to 0, important!!

	mdvar.set_virial(::bVirial, true);
	mdvar.spme_var.bEfieldOnly = false;
	mdvar.spme_var.esv1.bEwald = true;
	mdvar.spme_var.esv1.init_EwaldSum(mmcell.h[0] * mmcell.h[1] * mmcell.h[2] * 8, ::rcut_Ewd);
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
		mdvar.vspme_mu.xl[0] = -mmcell.h[0]; mdvar.vspme_mu.xl[1] = -mmcell.h[1]; mdvar.vspme_mu.xl[2] = -mmcell.h[2];
		mdvar.vspme_mu.cell_dim(mmcell.h[0] * 2, mmcell.h[1] * 2, mmcell.h[2] * 2); // important, after init();
		mdvar.vspme_mu.init_k();
		mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();

		mdvar.vspme_mu_induced.bsp = &(::bsp);
		init_spme_induced(&mmcell, mdvar.vspme_mu_induced, &(mdvar.q0)); // init the viral atoms
		mdvar.vspme_mu_induced.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_mu_induced.xl[0] = -mmcell.h[0]; mdvar.vspme_mu_induced.xl[1] = -mmcell.h[1]; mdvar.vspme_mu_induced.xl[2] = -mmcell.h[2];
		mdvar.vspme_mu_induced.cell_dim(mmcell.h[0] * 2, mmcell.h[1] * 2, mmcell.h[2] * 2); // important, after init();
		mdvar.vspme_mu_induced.init_k();
		mdvar.vspme_mu_induced.init_b(); mdvar.vspme_mu_induced.init_C(); mdvar.vspme_mu_induced.init_vars();
	}
	else if (::bDipole) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		mdvar.vspme_mu.bsp = &(::bsp);
		init_spme(&mmcell, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_mu.xl[0] = -mmcell.h[0]; mdvar.vspme_mu.xl[1] = -mmcell.h[1]; mdvar.vspme_mu.xl[2] = -mmcell.h[2];
		mdvar.vspme_mu.cell_dim(mmcell.h[0] * 2, mmcell.h[1] * 2, mmcell.h[2] * 2); // important, after init();
		mdvar.vspme_mu.init_k();
		mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();
	}
	else { // charge only
		mmcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		mdvar.vspme_q.bsp = &(::bsp);
		init_spme(&mmcell, mdvar.vspme_q); // init the viral atoms
		mdvar.vspme_q.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		mdvar.vspme_q.xl[0] = -mmcell.h[0]; mdvar.vspme_q.xl[1] = -mmcell.h[1]; mdvar.vspme_q.xl[2] = -mmcell.h[2];
		mdvar.vspme_q.cell_dim(mmcell.h[0] * 2, mmcell.h[1] * 2, mmcell.h[2] * 2); // important, after init();
		mdvar.vspme_q.init_k();
		mdvar.vspme_q.init_b(); mdvar.vspme_q.init_C(); mdvar.vspme_q.init_vars();
	}

	mdvar.spme_var.esv1.iSurface = 0; // 3d
	mdvar.spme_var.esv1.iSurfaceBoundary = ::iSurfaceBoundary; // normal
	mdvar.spme_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
}

void init_AtCell(MMOL_MD_CELL &mmcell, NVT_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> &mdvar, _atomic_cmm_::EAtCELL &ecell, _atomic_cmm_::SR_AtCELL &srcell) {
	ecell.rmax_c = ::size_cluster / 2; // maximum radius of the cluster
	srcell.rmax_c = ::size_cluster / 2; // maximum radius of the cluster
	ecell.rcut = ::rcut_Ewd;
	srcell.rcut = ::rcut_LJ;

	float cell_corner[3] = {-mmcell.h[0], -mmcell.h[1], -mmcell.h[2]};
	int ngrid[3] = {int(mmcell.h[0]), int(mmcell.h[1]), int(mmcell.h[2])};
	if (ngrid[0] < 10) ngrid[0] = 10;
	if (ngrid[1] < 10) ngrid[1] = 10;
	if (ngrid[2] < 10) ngrid[2] = 10;
	float cell_cw[3] = {mmcell.h[0] * 2 / ngrid[0], mmcell.h[1] * 2 / ngrid[1], mmcell.h[2] * 2 / ngrid[2]}; // is this right ?
	_atomic_cmm_::init_EAtCell(mmcell, ecell, ngrid[0], ngrid[1], ngrid[2], cell_corner, cell_cw);
	//distribute_atoms(&ecell);
	_atomic_cmm_::update_atomic_charge_dipole(&ecell);
	ecell.esv.cp_EwaldSum_vars(*((_EwaldSum_real_::EwaldSumRealVars0*)&mdvar.spme_var.esv1));
	ecell.esv.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	ecell.bCalE = 1; // calculate force
	ecell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

	_atomic_cmm_::init_SRAtCell(mmcell, srcell, ngrid[0], ngrid[1], ngrid[2], cell_corner, cell_cw);
	//distribute_atoms(&srcell);
	srcell.bCalSR = 1; // calculating short-range L-J or other interaction )
}

extern void init_velocity(MMOL_MD_CELL &mmcell);

extern bool init_mc_species(MMOL_MD_CELL &mmcell, MC_SPECIES &sp, _bias_force_::BIAS_POINT &mdcell_bias);

extern void relocate_mc_species(MMOL_MD_CELL &mmcell, MC_SPECIES &sp);

extern void molecule_check_periodic_cell(MMOL_MD_CELL &mmcell);
extern void MDCELL_molecule_PeriodicCell_check(MDCELL_ThreadJob<MMOL_MD_CELL> *job);

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
extern void NetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob);
extern void CompensateNetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 3> *sv3);
extern void CompensateNetForceZ(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 3> *sv3);
extern void ApplyZSpring(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob, float *k);
extern void InterfaceConstraints_CMM(CMM_CELL3D<BASIC_CLUSTER> *cmm, int nThreads, double &Ut);
extern void MDCELL_CalMolTransKineticEnergy(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4> *Ekt);
extern void MDCellMassCenter(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob);

extern void MDCELL_SpatialMomentumZ(MDCELL_ThreadJob<MMOL_MD_CELL> *job, double *mp);

extern void ZeroSpatialMomentum(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob);

extern char parfname_IConstraint[256];

// electrostatic force calculation is controlled by ::bPolarize, ::bDipole, ::bRealEwaldSumOnly, ::local_interact

void GPU_MDRun(MMOL_MD_CELL &mmcell, _atomic_cmm_::EAtCELL &ecell, _atomic_cmm_::SR_AtCELL &srcell, NVT_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> &mdvar, MDControlVars &cvar, long LOOPS, bool bZeroSpatialMomentum) {
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}

	MDControlVars cv0;
	cv0.backup_global_vars();
	cvar.update_global_vars();
	bool bSaveVirial = cvar.bSaveVirial;
	bool bSavePz = cvar.bSavePz;
	bool bSave_mu = cvar.bSave_mu;
	cvar.show_status();

	// we need to make sure the mass center of the cell is at zero
	molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell
	VECTOR3 mass_center;
	set_free_cell_mass_center(mmcell, mass_center);

	show_log("MD of 3d cell ...", true);

	molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell

	// *****************  over  ******************
	show_all_molecules();
/*
	if (strlen(parfname_IConstraint) > 0) {
		if (init_clusterIC_DB(2)) {
			sprintf(errmsg, "Interface defined in %s is enabled!", parfname_IConstraint); show_log(errmsg, true);
			apply_slab_barrier_constraint(clusterIC_db);
		}
		else if (clusterIC_db.n == 0 || clusterIC_db.n != 2) {
			sprintf(errmsg, "2 interfaces are required to define a slab. check %s !", parfname_IConstraint);
			show_log(errmsg, true); return;
		}
	}
	else {
		sprintf(errmsg, "To ensure the species are not evaporated out of FFT space, external interfacial constraint is required!");
		show_log(errmsg, true); //return;
	}
*/
#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

	double kappa = mdvar.spme_var.esv1.kappa;
	sprintf(errmsg, "Ewald-Sum kappa = %f, rcut = %f", kappa, ::rcut_Ewd); show_log(errmsg, true);
#endif

	bool bStorageChain = false;

	int i = 0;
	long nloop = 0;
	int nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0, U_LJ = 0, U_estat = 0;
	double Ektxx, Ektyy, Ektzz;

	double Uext = 0; // external electric field

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	int nMDThreads = MAX_THREADS;

	MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> mdcellHDynThreadVar[MAX_THREADS];
	AssignJob_MDCELL_HDYN<MMOL_MD_CELL>(&mmcell, mdcellHDynThreadVar, nMDThreads, false);

	double djob[MAX_THREADS];
	MATRIX<3> mjob[MAX_THREADS];
	bool status_job[MAX_THREADS];
	//INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];

	// make sure the spatial momentum of the whole cell is zero
	MDCellMassCenter(mdvar.job_mdcell, sv3job, nMDThreads);
	
	if (bZeroSpatialMomentum) {
		if (mmcell.bHDynMC_T) {
			ZeroSpatialMomentum(mdvar.job_mdcell, sv3job, MAX_THREADS);
		}
	}

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

#define CPU_SPME  0
#define GPU_SPME  1
#ifdef __GPU__
	int spme_method = GPU_SPME;

	_gpu_md_::GPU_SPME_3d spme_gpu;
#else
	int spme_method = CPU_SPME;
#endif

	// polarization could be enabled / disabled
	if (::bPolarize || ::bDipole) {
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral
		if (mdvar.vspme_mu.va.n != mmcell.catom.n) init_spme_3d(mmcell, mdvar);
	}
	else {
		mmcell.set_atom_multipole(0); // atoms are NOT polarizable
		if (mdvar.vspme_q.va.n != mmcell.catom.n) init_spme_3d(mmcell, mdvar);
	}

	if (spme_method == CPU_SPME) {
		show_log("using CPU for SPME, image part of Ewald-Sum.", true);
	}
	else if (spme_method == GPU_SPME) {
#ifdef __GPU__
		float xl[3] = {-mmcell.h[0], -mmcell.h[1], -mmcell.h[2]};
		float xd[3] = {mmcell.h[0] * 2, mmcell.h[1] * 2, mmcell.h[2] * 2};
		spme_gpu.atom.set_array(mmcell.atom.n); _gpu_md_::init_atom(mmcell.atom, spme_gpu.atom);
		spme_gpu.init_memory(mmcell.atom.n);
		if (!_gpu_md_::setup_gpu_spme(&spme_gpu, &::bsp, xl, xd, ::dim_spme, mdvar.spme_var.esv1.kappa, mmcell.atom.n, (bPolarize ? 1: 0))) {
			show_msg("GPU failure to initialize the memory for SPME", true); return;
		}
		show_log("using GPU for SPME, image part of Ewald-Sum.", true);
#endif
	}

	ecell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

	//int LOOPS = mmcell.LOOPS;
	int save_indx = 0;
	MD_SAVE *mdsave = mmcell.msave;
	AsyncMDSave_VARS asyncMDSave;

	ITER_SAVE virial_save;
	if (bSaveVirial) {
		virial_save.set_file(mdsave->title, "virial");
		virial_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	ITER_SAVE Pz_save;
	if (bSavePz) {
		Pz_save.set_file(mdsave->title, "Pz");
		Pz_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	ITER_SAVE mu_save;
	//bool bSave_mu = (bPolarize || ::bDipole ? true : false);
	if (bSave_mu) {
		mu_save.set_file(mdsave->title, "mu");
		mu_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	AsynRun<_atomic_cmm_::AsynLJVar> asynRun_LJ;
	_atomic_cmm_::AsynLJVar asynVar_LJ;
	asynVar_LJ.cell = &srcell; asynVar_LJ.nThreads = MAX_THREADS;
	asynVar_LJ.res = &(mdvar.MDVar_LJ::LJ_res);
	asynRun_LJ.set_par(&asynVar_LJ);

	int ncs_mmcell = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) ncs_mmcell += mmcell.mm.m[nm].nCluster;
	ncs_mmcell += mmcell.sm.n + mmcell.pm.n;

#if _CHECK_TIME_STAMP_
	tsave.set_file(mdsave->title, "t");
	//tsave.set_iter(mdsave->nsave, mdsave->max_save);
	tsave.set_iter(1, mdsave->max_save);
	long dt = 0;
#endif


	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command

		iCheckCMM = 0;

		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			//case MULTI_MM_PERIDIC_CELL_MD:
				// here we have to check the periodity of the molecule position
				//molecule_check_periodic_cell(mmcell);
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_molecule_PeriodicCell_check),
					mdvar.job_mdcell, nMDThreads);
			//	break;
		}

#if _CHECK_TIME_STAMP_
tsave.one_more();
mt.start();
if (tsave.bsave()) (*tsave.out)<<nloop<<"  ";
#endif

		check_atom_rg(mmcell, nMDThreads);

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalInertiaMoment), mdvar.job_mdcell, nMDThreads);

		MDCellMassCenter(mdvar.job_mdcell, sv3job, nMDThreads);

		if (bCheckClusterInCMM) { // update distribution, including coordination
			_atomic_cmm_::distribute_atoms(&ecell);
			_atomic_cmm_::distribute_atoms(&srcell);
		}
		else { // update coordination only
			_atomic_cmm_::update_atoms(&ecell);
			_atomic_cmm_::update_atoms(&srcell);
		}

#ifdef __GPU__

		if (ecell.gpu_cell.srnb_dev == NULL && ecell.bCalUseGPU) {
			if (!ecell.check_neighbor_memory(true)) return;
			gpu_reset_esr_neighbors(&ecell);
			check_neighbors(&(ecell.gpu_cell), ecell.gpu_cell_dev);
		}

		if (spme_method == GPU_SPME) {
			update_atomic_charge_dipole(&ecell);
			memcpy(spme_gpu.r.m, ecell.r.m, ecell.r.n * sizeof(double));
			//_gpu_md_::update_atomic_coordinate(&spme_gpu, (MAX_THREADS < 4 ? MAX_THREADS : 4));

			if (nloop == 0) {
				memcpy(spme_gpu.q.m, ecell.c.m, ecell.c.n * sizeof(double));
				_gpu_md_::_spme_::update_q_device_spme(&(spme_gpu.vspme_host));
			}

			memcpy(spme_gpu.mu.m, ecell.mu.m, ecell.mu.n * sizeof(double));
			_gpu_md_::_spme_::update_dipole_device_spme(&(spme_gpu.vspme_host));
		}
#endif

		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_InitClusterForce), mdvar.job_mdcell, nMDThreads);

	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		dt = mt.elapse();
		(*tsave.out)<<dt<<"  ";
	#endif

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
			if (srcell.bCalUseGPU) asynRun_LJ.start_func((void*)(&_atomic_cmm_::cudaLJ));
			else
#endif
				asynRun_LJ.start_func((void*)(&_atomic_cmm_::cpuLJ));
		}

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (spme_method == CPU_SPME) {
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
			}
			U_estat = mdvar.spme_res.U;
		}
#ifdef __GPU__
		else if (spme_method == GPU_SPME) {
			// calculat the total dipole of the whole cell
			cal_dipole(mmcell, 1, mmcell.mu_surf);
			ecell.mu_surf = mmcell.mu_surf;

			ecell.bCalE = 2; // 1 -- Efield, 2 -- force & U
			ecell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

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
			}

			U_estat = mdvar.spme_res.U;
		}
#endif // __GPU__

#endif //_IGNORE_ELECTROSTATIC_INTERACTION_


#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
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
		}

		//virial.v[0] = 0; virial.v[1] = 0; virial.v[2] = 0; virial.v[3] = 0; Ep = 0;

		if (bSavePz) {
			if (!Pz_save.one_more()) break;
		}

		if (bSave_mu) {
			mu_save.one_more();
			if (mu_save.bsave()) {
				*(mu_save.out)<<nloop<<"  "<<mmcell.mu_surf.mu.v[0]<<"  "<<mmcell.mu_surf.mu.v[1]<<endl;
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
			for (i = 0; i < nMDThreads; i++) memset(sv4job[i].v, 0, 4 * SIZE_DOUBLE);
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

		if (bSaveVirial) {
			if (virial_save.bsave()) {
				(*virial_save.out)<<nloop<<"  "<<mdvar.virial.v[0]<<"  "<<mdvar.virial.v[1]
					<<"  "<<mdvar.virial.v[2]<<"  "<<mdvar.virial.v[3];
			}
		}

		// NO Implicit-Solvent

		//MOperate<BASIC_CLUSTER*>((void*)(&SCALE_FORCE_CLUSTER), mmcell.bcluster.m, mmcell.bcluster.n, MAX_THREADS);
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

	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
		}
	#endif

		if (clusterIC_db.n > 0) {
			// the force from interface, and the torque
			// however, the force will not contribute to the virial
			InterfaceConstraints_AtCell((_atomic_cmm_::_AtCELL*)&srcell, nMDThreads, djob[0]); Ep += djob[0];
		}

	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
		}
	#endif

		iCheckCMM++; iCheckCMM = 0;

		// self-consistent velocity and acceleration, Nose-Hoover Chain
		bSpeedVerlet = (nloop == 0 ? false : true);
		if (!LeapFrogVerlet_SelfConsistent_SpeedAccel(&mmcell, mdvar.job_mdcell, mdcellHDynThreadVar, bSpeedVerlet, nMDThreads)) {
			sprintf(errmsg, "failure to get consistent velocity and acceleration at loop %d", nloop);
			show_log(errmsg, true);
		}

	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
		}
	#endif

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
/*
			calHooverVirial(mdvar.job_mdcell, nMDThreads, sv3job);
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] += sv3job[i].v[0];
				mdvar.virial.v[1] += sv3job[i].v[1];
				mdvar.virial.v[2] += sv3job[i].v[2];
				mdvar.virial.v[3] += sv3job[i].v[0] + sv3job[i].v[1] + sv3job[i].v[2];
			}
*/
			MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_CalMolTransKineticEnergy), mdvar.job_mdcell, nMDThreads, sv4job);
			Ek_t = 0; Ektxx = 0; Ektyy = 0; Ektzz = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ektxx += sv4job[i].v[0];
				Ektyy += sv4job[i].v[1];
				Ektzz += sv4job[i].v[2];
				Ek_t += sv4job[i].v[3];
			}
			mdvar.Pint = (mdvar.virial.v[3] + Ek_t * 2 / unit_Ek_kT) / (::fUnit_atm_u * mmcell.h[0] * mmcell.h[1] * mmcell.h[2] * 8 * 3);
			if (bSaveVirial) {
				(*virial_save.out)<<"          "<<(mdvar.virial.v[0] + Ektxx * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pxx -- mN/meter
				(*virial_save.out)<<"  "<<(mdvar.virial.v[1] + Ektyy * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pyy -- mN/meter
				(*virial_save.out)<<"  "<<(mdvar.virial.v[2] + Ektzz * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"  "<<0.5 * (mdvar.virial.v[2] - (mdvar.virial.v[0] + mdvar.virial.v[1]) * 0.5 + (Ektzz - 0.5 * (Ektxx + Ektyy)) * 2 / ::unit_Ek_kT) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"  "<<0.5 * (mdvar.virial.v[2] - (mdvar.virial.v[0] + mdvar.virial.v[1]) * 0.5) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"          "<<mdvar.Pint<<"      "<<mmcell.rc.v[2]<<endl;
/*
				(*virial_save.out)<<"          "<<Ektxx;  // Pxx -- mN/meter
				(*virial_save.out)<<"  "<<Ektyy;  // Pyy -- mN/meter
				(*virial_save.out)<<"  "<<Ektzz;  // Pzz -- mN/meter
				(*virial_save.out)<<"  "<<Ektzz - 0.5 * (Ektxx + Ektyy);  // Pzz -- mN/meter
				(*virial_save.out)<<"  "<<0.5 * (mdvar.virial.v[2] - (mdvar.virial.v[0] + mdvar.virial.v[1]) * 0.5) / (mmcell.h[0] * mmcell.h[1] * 4) * 1.66e6;  // Pzz -- mN/meter
				(*virial_save.out)<<"          "<<mdvar.Pint<<endl;
				*/
			}
		}

		if (bSavePz && Pz_save.bsave()) {
			{
				MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, double >((void*)(&MDCELL_SpatialMomentumZ), mdvar.job_mdcell, MAX_THREADS, djob);
				double Pz = 0, vz = 0;
				float M = mmcell.M;
				int i;
				for (i = 0; i < MAX_THREADS; i++) Pz += djob[i];
				vz = Pz / M;
				*(Pz_save.out)<<nloop<<"  ";
				*(Pz_save.out)<<mmcell.fc.v[0]<<"  "<<mmcell.fc.v[1]<<"  "<<mmcell.fc.v[2]<<"  "<<mmcell.Ekt_mc<<"  "<<Pz * Pz / M * 0.5 * unit_Ek_kT<<"  "<<Pz<<"  "<<vz;
				if (::bZSpring) *(Pz_save.out)<<"  "<<mmcell.rc.v[2];
				*(Pz_save.out)<<endl;
			}
		}

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_Verlet), mdvar.job_mdcell, nMDThreads);
		// accept current NHC vars and prepair for next step
		for (i = 0; i < mmcell.hdyn.n; i++) mmcell.hdyn.m[i].nhc->accept_nextstep();
		if (mmcell.bHDynMC_T) mmcell.mcHDyn.nhc->accept_nextstep();

		#if _CHECK_TIME_STAMP_ == 1
			//mt.start();
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
			//sprintf(::errmsg, "checking procedure takes %d ms", time.glance());
			//show_infor(errmsg);
		#endif

		if (bSaveVirial) virial_save.save_over();
		if (bSavePz) Pz_save.save_over();
		if (bSave_mu) mu_save.save_over();
#if _CHECK_TIME_STAMP_
		if (tsave.bsave()) (*tsave.out)<<endl;
		tsave.save_over();
#endif

		nloop++;
	}
	/*
	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
	*/
	mdsave->flush_buf();
	if (bSaveVirial) virial_save.reset();
	if (bSavePz) Pz_save.reset();
	if (bSave_mu) mu_save.reset();
#if _CHECK_TIME_STAMP_
	tsave.reset();
#endif

	cv0.update_global_vars();
}




void GPU_MD_PROCS(MMOL_MD_CELL &mmcell, MD_PRO *mdpro, int npro) {
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}

	// we need to make sure the mass center of the cell is at zero
	//molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell
	//VECTOR3 mass_center;
	//set_free_cell_mass_center(mmcell, mass_center);

	show_log("running molecular dynamic on NVT essemble ...", true);

	//molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************
	show_all_molecules();

	if (strlen(parfname_IConstraint) > 0) {
		if (init_clusterIC_DB(2, -mmcell.h[2], mmcell.h[2])) {
			sprintf(errmsg, "Interface defined in %s is enabled!", parfname_IConstraint); show_log(errmsg, true);
			apply_slab_barrier_constraint(clusterIC_db);
		}
		else if (clusterIC_db.n == 0 || clusterIC_db.n != 2) {
			sprintf(errmsg, "2 interfaces are required to define a slab. check %s !", parfname_IConstraint);
			show_log(errmsg, true); return;
		}
	}
	else {
		//sprintf(errmsg, "To ensure the species are not evaporated out of FFT space, external interfacial constraint is required!");
		//show_log(errmsg, true); //return;
	}

	int i;
	bool bStorageChain = false;
	extern void InitMD(MMOL_MD_CELL &mmcell, bool bStorageChain);
	InitMD(mmcell, bStorageChain);

	NVT_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	init_spme_3d(mmcell, mdvar);

	mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::init_job(mmcell, NULL);
	sprintf(errmsg, "Job # : %d; working with molecules:", MAX_THREADS); show_log(errmsg, true);
	for (i = 0; i < MAX_THREADS; i++) {
		sprintf(errmsg, "%d -- MM %d, SM %d, PM %d", i, mdvar.mdcellJob[i].mmIndx.n, mdvar.mdcellJob[i].smIndx.n, mdvar.mdcellJob[i].pmIndx.n);
		show_log(errmsg, true);
	}

	// ecell, srcell, spme_gpu ...
	
#ifdef __GPU__
	if (cudaSetDevice(_CudaDevice_) != cudaSuccess) {
		show_infor("GPU is not enabled!", true); return;
	}
	else show_log("GPU is used!", true);
	cudaDeviceReset();
#endif

	_atomic_cmm_::EAtCELL ecell;
	_atomic_cmm_::SR_AtCELL srcell;
#ifdef __GPU__
	ecell.bCalUseGPU = bGPURealEwaldSum;
	srcell.bCalUseGPU = false;

	if (srcell.bCalUseGPU) show_log("using GPU for L-J interaction", true);
	else show_log("using CPU for L-J interaction", true);

	if (ecell.bCalUseGPU) show_log("using GPU for real part Ewald-Sum.", true);
	else show_log("using CPU for real part Ewald-Sum.", true);
#else
	show_log("using CPU for L-J interaction, and real part Ewald-Sum.", true);
#endif
	init_AtCell(mmcell, mdvar, ecell, srcell);

	double kappa = mdvar.spme_var.esv1.kappa;
	sprintf(errmsg, "Ewald-Sum kappa = %f, rcut = %f", kappa, ::rcut_Ewd); show_log(errmsg, true);
//#endif

	MD_PRO mdpro0; mdpro0.cvar.backup_global_vars(); mdpro0.nloops = mmcell.LOOPS;
	char oftitle[256] = "\0", toftitle[256] = "\0"; strcpy(oftitle, mmcell.msave->title);
	char msg[2560] = "\0", oldparIConstraint[256] = "\0";
	float eThermal_old = ::Ek0Internal;
	strcpy(oldparIConstraint, ::parfname_IConstraint);
	bool bConstraintChanged = false;
	int ipro = 0;

	init_velocity(mmcell); // new velocity

	for (ipro = 0; ipro < npro; ipro++) {
		mdpro[ipro].cvar.bSave_mu = false; mdpro[ipro].cvar.bSavePz = false; mdpro[ipro].cvar.bSaveVirial = false;
		sprintf(toftitle, "%s.%dpro", oftitle, ipro);
		mmcell.msave->set_title(toftitle);

		if (strlen(mdpro[ipro].parfname_IConstraint) > 0) {
			bConstraintChanged = true;
			strcpy(::parfname_IConstraint, mdpro[ipro].parfname_IConstraint);
			if (init_clusterIC_DB(2, -mmcell.h[2], mmcell.h[2])) {
				sprintf(errmsg, "use constraints defined in %s. Interface constraint is enabled!", ::parfname_IConstraint); show_log(errmsg, true);
				apply_slab_barrier_constraint(clusterIC_db);
			}
			else if (clusterIC_db.n == 0 || clusterIC_db.n != 2) {
				sprintf(errmsg, "failure to obtain constraint from file %s. check the file !", ::parfname_IConstraint); show_log(errmsg, true); return;
			}
		}
		else {
			sprintf(errmsg, "use original constraints !");
			show_log(errmsg, true); //return;
		}

		if (mdpro[ipro].eThermal != 0) ::Ek0Internal = mdpro[ipro].eThermal;
		else ::Ek0Internal = eThermal_old;
		mmcell.check_freedom();
		sprintf(msg, "thermal energy of each freedom : %5.2f kT", mmcell.Ethermal); show_log(msg, true);

		//init_velocity(mmcell); // new velocity
		GPU_MDRun(mmcell, ecell, srcell, mdvar, mdpro[ipro].cvar, long(mdpro[ipro].nloops) * long(mmcell.msave->nsave), ipro == 0 ? true : false);
		mmcell.msave->flush_buf(); mmcell.msave->close();

		if (COMMAND_MD_STOP) break; // stop with given command
	}

	::Ek0Internal = eThermal_old; mmcell.check_freedom();
	sprintf(msg, "thermal energy of each freedom : %5.2f kT", mmcell.Ethermal); show_log(msg, true);
	if (bConstraintChanged) {
		strcpy(::parfname_IConstraint, oldparIConstraint);
		if (init_clusterIC_DB(2,  -mmcell.h[2], mmcell.h[2])) {
			sprintf(errmsg, "use constraints defined in %s. Interface constraint is enabled!", ::parfname_IConstraint); show_log(errmsg, true);
			apply_slab_barrier_constraint(clusterIC_db);
		}
	}
	else {
		sprintf(errmsg, "use original constraints !");
		show_log(errmsg, true); //return;
	}

	mdpro0.cvar.update_global_vars(); mdpro0.cvar.bSave_mu = true; mdpro0.cvar.bSavePz = true; mdpro0.cvar.bSaveVirial = true;
	mmcell.msave->set_title(oftitle);
	if (npro == 0) init_velocity(mmcell); // new velocity
	GPU_MDRun(mmcell, ecell, srcell, mdvar, mdpro0.cvar, mdpro0.nloops, npro == 0 ? true : false);
	mmcell.msave->flush_buf(); mmcell.msave->close();

	if (MAX_THREADS > 1) {
		for (int i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
#if _CHECK_TIME_STAMP_
	tsave.reset();
#endif
}

extern bool read_mdprocs(char *fname, ARRAY<MD_PRO> &mdpro);
extern bool read_mc_species(char *fname, MC_SPECIES &sp);
extern void reset_random_seeds();

void GPU_MD_MCONFS(MMOL_MD_CELL &mmcell, MC_SPECIES &mcsp, MD_PRO *mdpro, int npro) {
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}
	int i;

	reset_random_seeds();

	// we need to make sure the mass center of the cell is at zero
	//molecule_check_2d_periodic_cell(mmcell); // check the molecule mass center is inside the cell
	//VECTOR3 mass_center;
	//set_free_cell_mass_center(mmcell, mass_center);

	show_log("running molecular dynamic on slab, with molecules randomly relocated periodically ...", true);

	//molecule_check_2d_periodic_cell(mmcell); // check the molecule mass center is inside the cell

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************
	show_all_molecules();

	if (strlen(parfname_IConstraint) > 0) {
		if (init_clusterIC_DB(2, -mmcell.h[2], mmcell.h[2])) {
			sprintf(errmsg, "Interface defined in %s is enabled!", parfname_IConstraint); show_log(errmsg, true);
			apply_slab_barrier_constraint(clusterIC_db);
		}
		else if (clusterIC_db.n == 0 || clusterIC_db.n != 2) {
			sprintf(errmsg, "2 interfaces are required to define a slab. check %s !", parfname_IConstraint);
			show_log(errmsg, true); return;
		}
	}
	else {
		sprintf(errmsg, "To ensure the species are not evaporated out of FFT space, external interfacial constraint is required!");
		show_log(errmsg, true); //return;
	}

	// ecell, srcell, spme_gpu ...
	/*
	if (cudaSetDevice(_CudaDevice_) != cudaSuccess) {
		show_infor("GPU is not enabled!", true); return;
	}
	else show_log("GPU #0 is used!", true);
	cudaDeviceReset();
	*/

//#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	//if (::rcut_Ewd > cmm->xd[0] / 3) ::rcut_Ewd = cmm->xd[0] / 3;
	//if (::rcut_Ewd > cmm->xd[1] / 3) ::rcut_Ewd = cmm->xd[1] / 3;
	//if (::rcut_Ewd > cmm->xd[2] / 3) ::rcut_Ewd = cmm->xd[2] / 3;

	bool bStorageChain = false;
	extern void InitMD(MMOL_MD_CELL &mmcell, bool bStorageChain);
	InitMD(mmcell, bStorageChain);
	NVT_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	init_spme_3d(mmcell, mdvar);

	mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::init_job(mmcell, NULL);
	sprintf(errmsg, "Job # : %d; working with molecules:", MAX_THREADS); show_log(errmsg, true);
	for (i = 0; i < MAX_THREADS; i++) {
		sprintf(errmsg, "%d -- MM %d, SM %d, PM %d", i, mdvar.mdcellJob[i].mmIndx.n, mdvar.mdcellJob[i].smIndx.n, mdvar.mdcellJob[i].pmIndx.n);
		show_log(errmsg, true);
	}

#ifdef __GPU__

	if (cudaSetDevice(_CudaDevice_) != cudaSuccess) {
		show_infor("GPU is not enabled!", true); return;
	}
	else show_log("GPU #0 is used!", true);
	//cudaDeviceReset();

#endif
	_atomic_cmm_::EAtCELL ecell;
	_atomic_cmm_::SR_AtCELL srcell;
#ifdef __GPU__
	ecell.bCalUseGPU = bGPURealEwaldSum; srcell.bCalUseGPU = false;

	if (srcell.bCalUseGPU) show_log("using GPU for L-J interaction", true);
	else show_log("using CPU for L-J interaction", true);

	if (ecell.bCalUseGPU) show_log("using GPU for real part Ewald-Sum.", true);
	else show_log("using CPU for real part Ewald-Sum.", true);
#else
	show_log("using CPU for L-J interaction, and real part Ewald-Sum.", true);
#endif
	init_AtCell(mmcell, mdvar, ecell, srcell);

	double kappa = mdvar.spme_var.esv1.kappa;
	sprintf(errmsg, "Ewald-Sum kappa = %f, rcut = %f", kappa, ::rcut_Ewd); show_log(errmsg, true);
//#endif

	init_mc_species(mmcell, mcsp, mdcell_bias);
	show_log("\nMOLECULES to be relocated : ", false);
	for (i = 0; i < mcsp.sp.n; i++) {
		if (mcsp.sp.m[i].uid >= 0) {show_log(mcsp.sp.m[i].name, false); show_log(", ", false);}
	}
	show_log("", true);

	MD_PRO mdpro0; mdpro0.cvar.backup_global_vars(); mdpro0.nloops = mmcell.LOOPS;
	char oftitle[256] = "\0", md_title[256] = "\0", toftitle[256] = "\0";
	strcpy(md_title, mmcell.msave->title);
	char msg[2560] = "\0", oldparIConstraint[256] = "\0";
	float eThermal_old = ::Ek0Internal;
	strcpy(oldparIConstraint, ::parfname_IConstraint);
	bool bConstraintChanged = false;
	int ipro = 0, iconf = 0;
	while (1) {
		show_log("\nrelocating species ... ", true);
		//relocate_mc_species(mmcell, mcsp);
		//molecule_check_periodic_cell(mmcell);

		relocate_mc_species(mmcell, mcsp, ::mdcell_bias);
		_bias_force_::biasOn = true; // turn on bias

		sprintf(oftitle, "%s.conf%d", md_title, iconf);
		sprintf(msg, "\nCONFIG #%d, output-file-title : %s", iconf, oftitle); show_log(msg, true);

		init_velocity(mmcell); // new velocity

		for (ipro = 0; ipro < npro; ipro++) {
			mdpro[ipro].cvar.bSave_mu = false; mdpro[ipro].cvar.bSavePz = false; mdpro[ipro].cvar.bSaveVirial = false;
			sprintf(toftitle, "%s.pro%d", oftitle, ipro);
			mmcell.msave->set_title(toftitle);
			sprintf(msg, "PROC #%d, output-file-title : %s", ipro, toftitle); show_log(msg, true);

			if (strlen(mdpro[ipro].parfname_IConstraint) > 0) {
				bConstraintChanged = true;
				strcpy(::parfname_IConstraint, mdpro[ipro].parfname_IConstraint);
				if (init_clusterIC_DB(2,  -mmcell.h[2], mmcell.h[2])) {
					sprintf(errmsg, "use constraints defined in %s. Interface constraint is enabled!", ::parfname_IConstraint); show_log(errmsg, true);
					apply_slab_barrier_constraint(clusterIC_db);
				}
				else if (clusterIC_db.n == 0 || clusterIC_db.n != 2) {
					sprintf(errmsg, "failure to obtain constraint from file %s. check the file !", ::parfname_IConstraint); show_log(errmsg, true); return;
				}
			}
			else {
				sprintf(errmsg, "use original constraints !");
				show_log(errmsg, true); //return;
			}

			if (mdpro[ipro].eThermal != 0) ::Ek0Internal = mdpro[ipro].eThermal;
			else ::Ek0Internal = eThermal_old;
			mmcell.check_freedom();
			sprintf(msg, "thermal energy of each freedom : %5.2f kT", mmcell.Ethermal); show_log(msg, true);

			//init_velocity(mmcell); // new velocity
			GPU_MDRun(mmcell, ecell, srcell, mdvar, mdpro[ipro].cvar, long(mdpro[ipro].nloops) * long(mmcell.msave->nsave), ipro == 0 ? true : false);
			mmcell.msave->flush_buf(); mmcell.msave->close();

			if (COMMAND_MD_STOP) break; // stop with given command
		}

		::Ek0Internal = eThermal_old; mmcell.check_freedom();
		sprintf(msg, "thermal energy of each freedom : %5.2f kT", mmcell.Ethermal); show_log(msg, true);
		if (bConstraintChanged) {
			strcpy(::parfname_IConstraint, oldparIConstraint);
			if (init_clusterIC_DB(2,  -mmcell.h[2], mmcell.h[2])) {
				sprintf(errmsg, "use constraints defined in %s. Interface constraint is enabled!", ::parfname_IConstraint); show_log(errmsg, true);
				apply_slab_barrier_constraint(clusterIC_db);
			}
		}
		else {
			sprintf(errmsg, "use original constraints !");
			show_log(errmsg, true); //return;
		}

		_bias_force_::biasOn = false; // turn off bias
		mdpro0.cvar.update_global_vars(); mdpro0.cvar.bSave_mu = true; mdpro0.cvar.bSavePz = true; mdpro0.cvar.bSaveVirial = true;
		mmcell.msave->set_title(oftitle);
		sprintf(msg, "FINAL PROC, output-file-title : %s", oftitle); show_log(msg, true);
		if (npro == 0) init_velocity(mmcell); // new velocity
		GPU_MDRun(mmcell, ecell, srcell, mdvar, mdpro0.cvar, mdpro0.nloops, npro == 0 ? true : false);
		mmcell.msave->flush_buf(); mmcell.msave->close();

		if (COMMAND_MD_STOP) break; // stop with given command

		iconf += 1;
	}

	if (MAX_THREADS > 1) {
		for (int i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
#if _CHECK_TIME_STAMP_
	tsave.reset();
#endif
}


