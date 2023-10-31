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

#include "mdvar.h"
#include "slab_md.h"

#include "md_module.h"

extern BSplineFunc<2> bsp;

#include "interface.h"
using namespace _interface_constraint_;
using namespace _EwaldSum_real_;

//************************************************************************************
//    2d-periodical slab is periodical in x-y, while not periodical in z axis
//************************************************************************************

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
extern ITER_SAVE tsave;
#endif

void init_spme_2d(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, Slab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> &mdvar) {
	//Slab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	mdvar.q0 = 0; // a charge for the interactions between induced dipoles. The charge has to be set to 0, important!!

	mdvar.set_virial(::bVirial, true);
	mdvar.spme_var.spme_3d_var.bEfieldOnly = false;
	mdvar.spme_var.spme_3d_var.esv1.bEwald = true; 
	//mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(cmm->xd[0] * cmm->xd[1] * cmm->xd[2], ::rcut_Ewd);
	mdvar.spme_var.spme_3d_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	mdvar.spme_var.spme_3d_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;

	init_BSplineFunc<2>(::bsp, ::indx_bspline);
	mdvar.vspme_q.mMP = 0; // charge only
	mdvar.vspme_mu.mMP = 1; // with dipole
	mdvar.vspme_mu_induced.mMP = 1; // with dipole of course
	if (bPolarize) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		init_spme(&mmcell, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.set_BSplineFunc(&(::bsp));
		mdvar.vspme_mu.set_cell(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], true, ::bExpandSPMECellZ); // important, after init();
		mdvar.vspme_mu.init_k(); mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();

		init_spme_induced(&mmcell, mdvar.vspme_mu_induced, &(mdvar.q0)); // init the viral atoms
		mdvar.vspme_mu_induced.set_BSplineFunc(&(::bsp));
		mdvar.vspme_mu_induced.set_cell(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], true, ::bExpandSPMECellZ); // important, after init();
		mdvar.vspme_mu_induced.init_k(); mdvar.vspme_mu_induced.init_b(); mdvar.vspme_mu_induced.init_C(); mdvar.vspme_mu_induced.init_vars();

		mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(mdvar.vspme_mu.V, ::rcut_Ewd);

		if (::bExpandSPMECellZ) {
			sprintf(errmsg, "SPME: expand cell z-axis from %f to %f", cmm->xd[2], mdvar.vspme_mu.xd[2]); show_log(errmsg, true);
		}
		else {
			sprintf(errmsg, "SPME: keep cell z-axis : %f / %f", cmm->xd[2], mdvar.vspme_mu.xd[2]); show_log(errmsg, true);
		}
	}
	else if (::bDipole) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		init_spme(&mmcell, mdvar.vspme_mu); // init the viral atoms
		mdvar.vspme_mu.set_BSplineFunc(&(::bsp));
		mdvar.vspme_mu.set_cell(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], true, ::bExpandSPMECellZ); // important, after init();
		mdvar.vspme_mu.init_k(); mdvar.vspme_mu.init_b(); mdvar.vspme_mu.init_C(); mdvar.vspme_mu.init_vars();

		mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(mdvar.vspme_mu.V, ::rcut_Ewd);

		if (::bExpandSPMECellZ) {
			sprintf(errmsg, "SPME: expand cell z-axis from %f to %f", cmm->xd[2], mdvar.vspme_mu.xd[2]); show_log(errmsg, true);
		}
		else {
			sprintf(errmsg, "SPME: keep cell z-axis : %f / %f", cmm->xd[2], mdvar.vspme_mu.xd[2]); show_log(errmsg, true);
		}
	}
	else { // charge only
		mmcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		init_spme(&mmcell, mdvar.vspme_q); // init the viral atoms
		mdvar.vspme_q.set_BSplineFunc(&(::bsp));
		mdvar.vspme_q.set_cell(cmm->xl[0], cmm->xl[1], cmm->xl[2], cmm->xd[0], cmm->xd[1], cmm->xd[2], ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], true, ::bExpandSPMECellZ); // important, after init();
		mdvar.vspme_q.init_k(); mdvar.vspme_q.init_b(); mdvar.vspme_q.init_C(); mdvar.vspme_q.init_vars();

		mdvar.spme_var.spme_3d_var.esv1.init_EwaldSum(mdvar.vspme_q.V, ::rcut_Ewd);

		if (::bExpandSPMECellZ) {
			sprintf(errmsg, "SPME: expand cell z-axis from %f to %f", cmm->xd[2], mdvar.vspme_q.xd[2]); show_log(errmsg, true);
		}
		else {
			sprintf(errmsg, "SPME: keep cell z-axis : %f / %f", cmm->xd[2], mdvar.vspme_q.xd[2]); show_log(errmsg, true);
		}
	}

	double kappa = mdvar.spme_var.spme_3d_var.esv1.kappa;
	mdvar.spme_var.K0var.bEfieldOnly = false;
	mdvar.spme_var.K0var.set(kappa, cmm->xd[0] * cmm->xd[1]);

	mdvar.spme_var.spme_3d_var.esv1.iSurface = 1; // 2d
	mdvar.spme_var.spme_3d_var.esv1.iSurfaceBoundary = ::iSurfaceBoundary; // normal
	mdvar.spme_var.spme_3d_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	mdvar.spme_var.K0var.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
}

void init_velocity(MMOL_MD_CELL &mmcell) {
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;

	int nm;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mm = mmcell.mm.m + nm;
		InitVelocity(*mm);
	}

	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		InitVelocity(*sm);
	}

	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m + nm;
		InitVelocity(*pm);
	}
}

bool init_mc_species(MMOL_MD_CELL &mmcell, MC_SPECIES &sp) {
	if (sp.xf[0] < -mmcell.h[0]) sp.xf[0] = -mmcell.h[0];
	if (sp.xt[0] > mmcell.h[0]) sp.xt[0] = mmcell.h[0];
	if (sp.xf[1] < -mmcell.h[1]) sp.xf[1] = -mmcell.h[1];
	if (sp.xt[1] > mmcell.h[1]) sp.xt[1] = mmcell.h[1];
	if (sp.xf[2] < -mmcell.h[2]) sp.xf[2] = -mmcell.h[2];
	if (sp.xt[2] > mmcell.h[2]) sp.xt[2] = mmcell.h[2];

	sp.mindx.release();
	bool status = true;
	int i, j, ns = 0, np = 0, n = 0;
	COMM_PARS<50> *mpar = NULL;
	for (i = 0; i < sp.sp.n; i++) {
		if ((mpar = mmcell.mTypeIndx_db.search_par(sp.sp.m[i].name)) == NULL) {
			show_log("ERROR for randomly relocating molecules : ", false); show_log(sp.sp.m[i].name, false); show_log(" is not defined in the cell", true); status = false;
			sp.sp.m[i].uid = -1;
		}
		else sp.sp.m[i].uid = mpar->uid;
	}
	n = 0;
	for (i = 0; i < mmcell.mm.n; i++) {
		for (j = 0; j < sp.sp.n; j++) {
			if (mmcell.mm.m[i].mTypeIndx == sp.sp.m[j].uid) {n++; break;}
		}
	}
	for (i = 0; i < mmcell.sm.n; i++) {
		for (j = 0; j < sp.sp.n; j++) {
			if (mmcell.sm.m[i].mTypeIndx == sp.sp.m[j].uid) {n++; break;}
		}
	}
	for (i = 0; i < mmcell.pm.n; i++) {
		for (j = 0; j < sp.sp.n; j++) {
			if (mmcell.pm.m[i].mTypeIndx == sp.sp.m[j].uid) {n++; break;}
		}
	}
	sp.mindx.set_array(n);
	if (n == 0) return status;
	n = 0;
	for (i = 0; i < mmcell.mm.n; i++) {
		for (j = 0; j < sp.sp.n; j++) {
			if (mmcell.mm.m[i].mTypeIndx == sp.sp.m[j].uid) {sp.mindx.m[n] = mmcell.mm.m[i].mIndx; n++; break;}
		}
	}
	for (i = 0; i < mmcell.sm.n; i++) {
		for (j = 0; j < sp.sp.n; j++) {
			if (mmcell.sm.m[i].mTypeIndx == sp.sp.m[j].uid) {sp.mindx.m[n] = mmcell.sm.m[i].mIndx; n++; break;}
		}
	}
	for (i = 0; i < mmcell.pm.n; i++) {
		for (j = 0; j < sp.sp.n; j++) {
			if (mmcell.pm.m[i].mTypeIndx == sp.sp.m[j].uid) {sp.mindx.m[n] = mmcell.pm.m[i].mIndx; n++; break;}
		}
	}
	return status;
}

extern int connect_ranf_server(const char *host);
extern float get_ranf();
extern void disconnect_ranf_server();

void relocate_mc_species(MMOL_MD_CELL &mmcell, MC_SPECIES &sp) {
	int imol, i, n;
	MMOLECULE *mm = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	float w[3] = {(sp.xt[0] - sp.xf[0]), (sp.xt[1] - sp.xf[1]), (sp.xt[2] - sp.xf[2])};
	float r[3];
	VECTOR3 dr;
	// randomly reorder the molecule list
	int m1, m2;

	  float franf = 0;
	  bool bUniverseRandom = false;
	  if (connect_ranf_server("127.0.0.1") > 0) bUniverseRandom = true;
	  if (bUniverseRandom) {show_log("Use unified random number from server", true);}
	  else {show_log("Use locally generated random number.", true);}

	  for (i = 0; i < sp.mindx.n; i++) {
	    franf = (bUniverseRandom ? get_ranf() : ranf());
	    m1 = int(sp.mindx.n * franf); 
	    franf = (bUniverseRandom ? get_ranf() : ranf());
	    m2 = int(sp.mindx.n * franf);
	    if (m1 != m2) { // exchange their positions
	      n = sp.mindx.m[m1]; sp.mindx.m[m1] = sp.mindx.m[m2]; sp.mindx.m[m2] = n;
	    }
	}
	for (i = 0; i < sp.mindx.n; i++) {
		n = 0; imol = sp.mindx.m[i];
		if (imol >= mmcell.mm.n) {imol -= mmcell.mm.n; n += 1;}
		if (imol >= mmcell.sm.n) {imol -= mmcell.sm.n; n += 1;}
		if (imol >= mmcell.pm.n) {imol -= mmcell.pm.n; n += 1;}

		franf = (bUniverseRandom ? get_ranf() : ranf());
		r[0] = sp.xf[0] + w[0] * franf;
		franf = (bUniverseRandom ? get_ranf() : ranf()); 
		r[1] = sp.xf[1] + w[1] * franf;
		franf = (bUniverseRandom ? get_ranf() : ranf()); 
		r[2] = sp.xf[2] + w[2] * franf;

		switch (n) {
		case 0: // MMOLECULE
			dr.v[0] = r[0] - mmcell.mm.m[imol].r0.v[0];
			dr.v[1] = r[1] - mmcell.mm.m[imol].r0.v[1];
			dr.v[2] = r[2] - mmcell.mm.m[imol].r0.v[2];
			mmcell.mm.m[imol].shiftMol(dr);
			break;
		case 1: // SMOLECULE
			dr.v[0] = r[0] - mmcell.sm.m[imol].r.v[0];
			dr.v[1] = r[1] - mmcell.sm.m[imol].r.v[1];
			dr.v[2] = r[2] - mmcell.sm.m[imol].r.v[2];
			mmcell.sm.m[imol].shiftMol(dr);
			break;
		case 2: // PMOLECULE
			dr.v[0] = r[0] - mmcell.pm.m[imol].r.v[0];
			dr.v[1] = r[1] - mmcell.pm.m[imol].r.v[1];
			dr.v[2] = r[2] - mmcell.pm.m[imol].r.v[2];
			mmcell.pm.m[imol].shiftMol(dr);
			break;
		default:
			break;
		}
	}
	disconnect_ranf_server();
}

extern void GPU_MDRun_SLAB(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, Slab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> &mdvar, MDControlVars &cvar, long LOOPS, bool bZeroSpatialMomentum);

extern void molecule_check_2d_periodic_cell(MMOL_MD_CELL &mmcell);
extern void ZeroSpatialMomentumZ(MDCELL_ThreadJob< MMOL_MD_CELL> *job, double *mp, int nJob);
extern void MMCELL_CMM_CheckRelship_2d(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator);
extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
extern void MMCELL_CMM_CheckImpSolvRelship_2d(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator);
extern void NetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob);
extern void CompensateNetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 6> *sv6);
extern void CompensateNetForceZ(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 6> *sv6);
extern void ApplyZSpring(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 6> *mp, int nJob, float *k);
extern void InterfaceConstraints_CMM(CMM_CELL3D<BASIC_CLUSTER> *cmm, int nThreads, double &Ut);
extern void MDCELL_CalMolTransKineticEnergy(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4> *Ekt);
extern void MDCELL_SpatialMomentumZ(MDCELL_ThreadJob<MMOL_MD_CELL> *job, double *mp);
extern void MDCELL_molecule_2d_PeriodicCell_check(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
extern void MDCellMassCenter(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 6> *mp, int nJob);

extern char parfname_IConstraint[256];

// following initializations were done !
	//bool bStorageChain = false;
	//InitMD(mmcell, cmm, bStorageChain);

	//Slab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	//init_spme_2d(mmcell, mdvar);

// electrostatic force calculation is controlled by ::bPolarize, ::bDipole, ::bRealEwaldSumOnly, ::local_interact

void MDRun_SLAB(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, Slab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> &mdvar, MDControlVars &cvar, long LOOPS, bool bZeroSpatialMomentum = false) {
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
	molecule_check_2d_periodic_cell(mmcell); // check the molecule mass center is inside the cell
	VECTOR3 mass_center;
	set_free_cell_mass_center(mmcell, mass_center);

	show_log("MD of slab ...", true);

	molecule_check_2d_periodic_cell(mmcell); // check the molecule mass center is inside the cell

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
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
	//if (::rcut_Ewd > cmm->xd[0] / 3) ::rcut_Ewd = cmm->xd[0] / 3;
	//if (::rcut_Ewd > cmm->xd[1] / 3) ::rcut_Ewd = cmm->xd[1] / 3;
	//if (::rcut_Ewd > cmm->xd[2] / 3) ::rcut_Ewd = cmm->xd[2] / 3;

	//bool bStorageChain = false;
	//InitMD(mmcell, cmm, bStorageChain);


	//Slab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	init_spme_2d(mmcell, cmm, mdvar);

	double kappa = mdvar.spme_var.spme_3d_var.esv1.kappa;
	sprintf(errmsg, "Ewald-Sum kappa = %f, rcut = %f", kappa, ::rcut_Ewd); show_log(errmsg, true);
#endif

	bool bStorageChain = false;

	int i = 0;
	int nloop = 0, nc = 0, nm = 0, ncell = 0;
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

/*
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
*/

	double djob[MAX_THREADS];
	MATRIX<3> mjob[MAX_THREADS];
	bool status_job[MAX_THREADS];
	//INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 6> sv6job[MAX_THREADS];

	if (bZeroSpatialMomentum) {
		// make sure the spatial momentum of the whole cell along z is zero
		if (mmcell.bHDynMC_T || mmcell.bHDynMC_R) {
			MDCELL_mass_center(mmcell);
			if (mmcell.bHDynMC_R) calCellInertiaTensor(mdvar.job_mdcell, mjob, nMDThreads);
			ZeroSpatialMomentumZ(mdvar.job_mdcell, djob, MAX_THREADS);
		}
	}

	int nCheckAngle = 5, iCheckAngle = 0;
	int nCheckSpatialMoment = 1, iCheckSpatialMoment = 0;
	int iCheckCMM = 0;
	int nCalibrateTorsionAngle = 10000, iCalibrateTorsionAngle = 0;
	bool bCheckClusterInCMM = false;
	bool Overlap = false;
	int nConvg = 0;
	bool bSpeedVerlet = true;

	mmcell.Init_MD(); //copy the initial torsion angle, speed to the mdpar
	if (::bPolarize) mmcell.init_dipole_hist(); // has to be done after Init_MD();

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

	LJ_AsynVARS asynVar_LJ;

	int ncs_mmcell = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) ncs_mmcell += mmcell.mm.m[nm].nCluster;
	ncs_mmcell += mmcell.sm.n + mmcell.pm.n;

#if _CHECK_TIME_STAMP_
	tsave.set_file(mdsave->title, "t");
	//tsave.set_iter(mdsave->nsave, mdsave->max_save);
	tsave.set_iter(1, mdsave->max_save);
	long dt = 0;
#endif
	
	show_log("using CPU for interactions calculation.", true);

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command

		iCheckCMM = 0;

		if (iCheckCMM == 0 || Ek > mmcell.Ethermal * 2) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			//case MULTI_MM_PERIDIC_CELL_MD:
				// here we have to check the periodity of the molecule position
				//molecule_check_periodic_cell(mmcell);
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_molecule_2d_PeriodicCell_check), 
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
		
		MDCellMassCenter(mdvar.job_mdcell, sv6job, nMDThreads);
		if (mmcell.bHDynMC_T || mmcell.bHDynMC_R) {
			//MDCELL_mass_center(mmcell);
			if (mmcell.bHDynMC_R) calCellInertiaTensor(mdvar.job_mdcell, mjob, nMDThreads);
		}

		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_InitClusterForce), mdvar.job_mdcell, nMDThreads);

	// re-distribute / check clusters in CMM, use vp of each cluster only
		if (bStorageChain) {
			if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
			else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
		}
		else {
			if (nloop == 0 || bCheckClusterInCMM) {
				if (MAX_THREADS == 1) {
					mdvar.job_mdcellCMM[0].cmm = cmm; // directly allocate to global cmm
					MMCELL_CMM_cluster_allocate(mdvar.job_mdcellCMM, status_job);
				}
				else {
					mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::job_mdcell_CMM_use_internal_cmm(); // set the JOB-related CMM is the job-related CMM
					MTOperate1<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, bool>(
						(void*)(&MMCELL_CMM_cluster_allocate), mdvar.job_mdcellCMM, nMDThreads, status_job);
					for (i = 0; i < nMDThreads; i++) {
						if (!status_job[i]) {
							sprintf(errmsg, "loop %d thread %d: buffer of CMM cell leaks when allocating cluster", nloop, i);
							show_log(errmsg, true);
						}
					}
					if (!combine_AtomAllocatedCMM_array<BASIC_CLUSTER>(mdvar.job_cmm, nMDThreads, cmm, true, true)) {
						sprintf(errmsg, "loop %d : buffer of CMM cell leaks when combining allocated cluster", nloop);
						show_log(errmsg, true);
					}
				}
				mdvar.MDVar_Job<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >::job_mdcell_CMM_use_external_cmm(cmm);
				//for (i = 0; i < MAX_THREADS; i++) job_mdcellCMM[i].cmm = cmm; // reset the JOB-related CMM is the main CMM
			}
/*
			if (bCheckClusterInCMM) {
			{
				nc = cmm->number_atoms();
				if (nc != ncs_mmcell) {
					sprintf(errmsg, "%d : cluster in CMM -- %d, real -- %d", nloop, nc, ncs_mmcell); show_log(errmsg, true);
				}
			}
			}
*/
		}
	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		dt = mt.elapse();
		(*tsave.out)<<dt<<"  ";
	#endif

		Ep = 0;
		// check relationship of each cluster, electrostatic and LJ
		if (iCheckCMM == 0) {
		//#if _CHECK_TIME_STAMP_ == 1
		//	time.start();
		//#endif
			if (bStorageChain) {
				//ParallelCheckRelshipInCells(*cmm, 0, cmm->acell.n - 1, MAX_THREADS); // setup relationship for each atom with CMM
				if (mmcell.mm.n > 0) _cmm_2d_::MMCheckRelshipInCell(mmcell.mm.m, mmcell.mm.n, *cmm, nMDThreads);
				if (mmcell.sm.n > 0) _cmm_2d_::SMCheckRelshipInCell(mmcell.sm.m, mmcell.sm.n, *cmm, nMDThreads);
				if (mmcell.pm.n > 0) _cmm_2d_::PMCheckRelshipInCell(mmcell.pm.m, mmcell.pm.n, *cmm, nMDThreads);
			}
			else {
				MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >(
					(void*)(&MMCELL_CMM_CheckRelship_2d), mdvar.job_mdcellCMM, nMDThreads);
			}
		//#if _CHECK_TIME_STAMP_ == 1
		//	sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
		//	show_infor(errmsg);
		//#endif
		}

		// important : calculate LJ interaction first, here cluster force & torque are reset
		if (::bMTForce) {
			Asyn_LJ_Interact(&mmcell, cmm, MAX_THREADS, mdvar.MDVar_LJ::LJ_var, mdvar.MDVar_LJ::LJ_res, asynVar_LJ);
		}
		else {
			LJ_Interact(&mmcell, cmm, MAX_THREADS, mdvar.MDVar_LJ::LJ_var, mdvar.MDVar_LJ::LJ_res);
			U_LJ = mdvar.MDVar_LJ::LJ_res.U;
		}

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!mmcell.eNeutral && ::bRealEwaldSumOnly) {
			if (::bPolarize || ::bDipole) _spme_2d_::LocalEF(&mmcell, cmm, mdvar.vspme_mu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			else _spme_2d_::LocalEF(&mmcell, cmm, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
		}
		else if (!::local_interact && !mmcell.eNeutral) {
			// For 2d x-y periodical EwalsSum, we do not need the total net dipole for EwaldSum!!
			// And, we have to set net dipole to zero because real part 2d EwaldSum used 3d EwaldSum program
			// where net cell dipole is used for energy, E and F calculation!!
			//V3zero(mmcell.mu_surf.mu) V3zero(mmcell.mu_surf.mu0) V3zero(mmcell.mu_surf.qr) mmcell.mu_surf.q = 0;
			// calculat the total dipole of the whole cell
			//cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));

			mdvar.spme_var.spme_3d_var.bEfieldOnly = false;
			mdvar.spme_var.K0var.bEfieldOnly = false;
			// important : calculate Coulomb secondly, here cluster force & torque are accumulated
			if (::bPolarize) {
				if (nloop >= 3) Guess_induced_dipole(mmcell, MAX_THREADS);
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				if (::bExcludeInducedDipole) {
					if (!_spme_2d_::Polarize_SPME_Interact2(&mmcell, cmm, mdvar.vspme_mu, mdvar.vspme_mu_induced, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
						sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
					}
				}
				else if (!_spme_2d_::Polarize_SPME_Interact(&mmcell, cmm, mdvar.vspme_mu, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
				}
				Backup_induced_dipole_hist(mmcell, MAX_THREADS);
			}
			else if (::bDipole) {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				_spme_2d_::SPME_Interact(&mmcell, cmm, mdvar.vspme_mu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			}
			else {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				_spme_2d_::SPME_Interact(&mmcell, cmm, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			}

			U_estat = mdvar.spme_res.U;
		}
#endif
		
#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
		}
#endif

		if (::bMTForce) {
			asynVar_LJ.thread_var.WaitUntilThreadOver();
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
				*(mu_save.out)<<nloop<<"  "<<cmm->mu_surf.mu.v[0]<<"  "<<cmm->mu_surf.mu.v[1]<<endl;
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
				(void*)(&MMCELL_CMM_CheckImpSolvRelship_2d), mdvar.job_mdcellCMM, nMDThreads);
			MTOperate1< MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, double>(
				(void*)(&MMCELL_cluster_impsolv_interact), mdvar.job_mdcellCMM, nMDThreads, djob);
			for (i = 0; i < nMDThreads; i++) Ep += djob[i];
		}

		/*
		if (scale_force) { // scaled in LJ interaction calculation
			CellOperate<BASIC_CLUSTER>((void*)(&SCALE_FORCE_CMM), *cmm, 0, cmm->acell.n, MAX_THREADS);
		}
		*/

		if (bRandomForce) gen_random_force(mmcell, mdvar.rand[0], mdvar.rand[1]);

		if (!::local_interact && ::mode_NetForce > 0) {
			NetForce(mdvar.job_mdcell, nMDThreads); // calculate net force and acceleration of the whole cell
			if (::bZSpring) {
				//ApplyZSpring(mdvar.job_mdcell, sv6job, nMDThreads, ::k_zspring);
				/*
				if (::mode_NetForce == 1) {
					mmcell.fc.v[0] += ::k_zspring[0] * mmcell.rc.v[0] + (mmcell.rc.v[0] > 0 ? 1 : -1) * ::k_zspring[1] * mmcell.rc.v[0] * mmcell.rc.v[0];
					mmcell.fc.v[1] += ::k_zspring[0] * mmcell.rc.v[1] + (mmcell.rc.v[1] > 0 ? 1 : -1) * ::k_zspring[1] * mmcell.rc.v[1] * mmcell.rc.v[1];
					mmcell.fc.v[2] += ::k_zspring[0] * mmcell.rc.v[2] + (mmcell.rc.v[2] > 0 ? 1 : -1) * ::k_zspring[1] * mmcell.rc.v[2] * mmcell.rc.v[2];
				}
				else if (::mode_NetForce == 2) {
				*/
					mmcell.fc.v[2] += ::k_zspring[0] * mmcell.rc.v[2] + (mmcell.rc.v[2] > 0 ? 1 : -1) * ::k_zspring[1] * mmcell.rc.v[2] * mmcell.rc.v[2];
				//}
			}
			if (::mode_NetForce == 1) {
				mmcell.ac.v[0] = mmcell.fc.v[0] / mmcell.M;
				mmcell.ac.v[1] = mmcell.fc.v[1] / mmcell.M;
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}
			else if (::mode_NetForce == 2) {
				mmcell.ac.v[2] = mmcell.fc.v[2] / mmcell.M;
			}

			if (::mode_NetForce == 1) CompensateNetForce(mdvar.job_mdcell, nMDThreads, sv6job); // compensate the net force on all specises in the cell
			else if (::mode_NetForce == 2) CompensateNetForceZ(mdvar.job_mdcell, nMDThreads, sv6job); // compensate the net force on all specises in the cell
			/*
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] += sv6job[i].v[0];
				mdvar.virial.v[1] += sv6job[i].v[1];
				mdvar.virial.v[2] += sv6job[i].v[2];
				mdvar.virial.v[3] += (sv6job[i].v[0] + sv6job[i].v[1] + sv6job[i].v[2]);
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
			InterfaceConstraints_CMM(cmm, nMDThreads, djob[0]); Ep += djob[0];
		}

	#if _CHECK_TIME_STAMP_ == 1
		//sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		//show_infor(errmsg);
		if (tsave.bsave()) {
			dt = mt.elapse();
			(*tsave.out)<<dt<<"  ";
		}
	#endif

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

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
			calHooverVirial(mdvar.job_mdcell, nMDThreads, sv6job);
			for (i = 0; i < nMDThreads; i++) {
				mdvar.virial.v[0] += sv6job[i].v[0];
				mdvar.virial.v[1] += sv6job[i].v[1];
				mdvar.virial.v[2] += sv6job[i].v[2];
				mdvar.virial.v[3] += sv6job[i].v[0] + sv6job[i].v[1] + sv6job[i].v[2];
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
				*(Pz_save.out)<<mmcell.fc.v[0]<<"  "<<mmcell.fc.v[1]<<"  "<<mmcell.fc.v[2]<<"  "<<mmcell.Ekr_mc<<"  "<<mmcell.Ekt_mc<<"  "<<Pz * Pz / M * 0.5 * unit_Ek_kT<<"  "<<Pz<<"  "<<vz;
				if (::bZSpring) *(Pz_save.out)<<"  "<<mmcell.rc.v[2];
				*(Pz_save.out)<<endl;
			}
		}

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_Verlet), mdvar.job_mdcell, nMDThreads);
		// accept current NHC vars and prepair for next step
		if (mmcell.bIndependentFrameHDyn) {
			MDCELL_HDyn_MultiThreadOperate<MMOL_MD_CELL>((void*)(&MDCELL_AcceptCurrentDynForNextStep),
				mdcellHDynThreadVar, nMDThreads);
		}
		else {
			for (i = 0; i < mmcell.hdyn.n; i++) mmcell.hdyn.m[i].nhc->accept_nextstep();
		}
		if (mmcell.bHDynMC_T) mmcell.mcHDyn.nhc->accept_nextstep();
		if (mmcell.bHDynMC_R) mmcell.mcrHDyn.nhc->accept_nextstep();

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




void SLAB_MD_PROCS(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, MD_PRO *mdpro, int npro) {
	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}

	// we need to make sure the mass center of the cell is at zero
	//molecule_check_2d_periodic_cell(mmcell); // check the molecule mass center is inside the cell
	//VECTOR3 mass_center;
	//set_free_cell_mass_center(mmcell, mass_center);

	show_log("running molecular dynamic on slab ...", true);

	//molecule_check_2d_periodic_cell(mmcell); // check the molecule mass center is inside the cell

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************
	show_all_molecules();

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

//#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	//if (::rcut_Ewd > cmm->xd[0] / 3) ::rcut_Ewd = cmm->xd[0] / 3;
	//if (::rcut_Ewd > cmm->xd[1] / 3) ::rcut_Ewd = cmm->xd[1] / 3;
	//if (::rcut_Ewd > cmm->xd[2] / 3) ::rcut_Ewd = cmm->xd[2] / 3;

	bool bStorageChain = false;
	InitMD(mmcell, cmm, bStorageChain);

	Slab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	
	int i;
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

	double kappa = mdvar.spme_var.spme_3d_var.esv1.kappa;
	sprintf(errmsg, "Ewald-Sum kappa = %f, rcut = %f", kappa, ::rcut_Ewd); show_log(errmsg, true);
//#endif

	MD_PRO mdpro0; mdpro0.cvar.backup_global_vars(); mdpro0.nloops = mmcell.LOOPS;
	char oftitle[256] = "\0", toftitle[256] = "\0"; strcpy(oftitle, mmcell.msave->title);
	char msg[2560] = "\0", oldparIConstraint[256] = "\0";
	float eThermal_old = ::Ek0Internal;
	strcpy(oldparIConstraint, ::parfname_IConstraint);
	bool bConstraintChanged = false;
	int ipro = 0;
	for (ipro = 0; ipro < npro; ipro++) {
		mdpro[ipro].cvar.bSave_mu = false; mdpro[ipro].cvar.bSavePz = false; mdpro[ipro].cvar.bSaveVirial = false; 
		sprintf(toftitle, "%s.%dpro", oftitle, ipro); 
		mmcell.msave->set_title(toftitle);

		if (strlen(mdpro[ipro].parfname_IConstraint) > 0) {
			bConstraintChanged = true;
			strcpy(::parfname_IConstraint, mdpro[ipro].parfname_IConstraint);
			if (init_clusterIC_DB(2)) {
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

		init_velocity(mmcell); // new velocity
#ifdef _USE_GPU_
		GPU_MDRun_SLAB(mmcell, cmm, mdvar, mdpro[ipro].cvar, long(mdpro[ipro].nloops) * long(mmcell.msave->nsave), (ipro == 0 ? true : false));
#else
		MDRun_SLAB(mmcell, cmm, mdvar, mdpro[ipro].cvar, long(mdpro[ipro].nloops) * long(mmcell.msave->nsave), (ipro == 0 ? true : false));
#endif
		mmcell.msave->flush_buf(); mmcell.msave->close();

		if (COMMAND_MD_STOP) break; // stop with given command
	}

	::Ek0Internal = eThermal_old; mmcell.check_freedom();
	sprintf(msg, "thermal energy of each freedom : %5.2f kT", mmcell.Ethermal); show_log(msg, true);
	if (bConstraintChanged) {
		strcpy(::parfname_IConstraint, oldparIConstraint);
		if (init_clusterIC_DB(2)) {
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
#ifdef _USE_GPU_
	GPU_MDRun_SLAB(mmcell, cmm, mdvar, mdpro[ipro].cvar, long(mdpro[ipro].nloops) * long(mmcell.msave->nsave), (npro == 0 ? true : false));
#else
	MDRun_SLAB(mmcell, cmm, mdvar, mdpro0.cvar, mdpro0.nloops, , (npro == 0 ? true : false));
#endif
	mmcell.msave->flush_buf(); mmcell.msave->close();

	if (MAX_THREADS > 1) {
		for (int i = 0; i < MAX_THREADS; i++) mdvar.job_cmm[i].reset_basic_cell_matrix();
	}
#if _CHECK_TIME_STAMP_
	tsave.reset();
#endif
}

#include "read.h"
bool read_mdprocs(char *fname, ARRAY<MD_PRO> &mdpro) {
	mdpro.release();
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		show_log("failure to open file : ", false); show_log(fname, true); return false;
	}
	char buffer[256] = "\0", *bf = NULL, msg[2560] = "\0";
	float eThermal = 0;
	char str[256] = "\0";
	int nprocs = 0, iproc = 0;
	int istat = 0;
	if (!search(in, "[MD_PROC]", buffer)) {
		show_log("[MD_PROC] is not defined in : ", false); show_log(fname, true); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nprocs) != 1) {
		show_log("number of procedures is not defined given in ", false); show_log(fname, false); 
		show_log(", see: ", false); show_log(buffer, true); in.close(); return false;
	}
	if (nprocs <= 0) {
		show_log("no additional procedures is defined", true); in.close(); return true;
	}
	mdpro.set_array(nprocs);
	for (iproc = 0; iproc < nprocs; iproc++) {
		in.getline(buffer, 250);
		bf = buffer;
		if (bf != NULL && sscanf(bf, "%d", &istat) == 1) mdpro.m[iproc].cvar.bPolarize = (istat > 0 ? true : false);
		else {
			sprintf(msg, "failure to get the status of POLARIZE (1st) from %s", buffer); mdpro.release(); in.close(); return false;
		}
		NextSection(&bf, ::seperator);

		if (bf != NULL && sscanf(bf, "%d", &istat) == 1) mdpro.m[iproc].cvar.bRealEwaldSumOnly = (istat > 0 ? true : false);
		else {
			sprintf(msg, "failure to get the status of REAL_EwaldSum_ONLY (2nd) from %s", buffer); mdpro.release(); in.close(); return false;
		}
		NextSection(&bf, ::seperator);

		if (bf != NULL && sscanf(bf, "%d", &istat) == 1) mdpro.m[iproc].cvar.bLocalInteract = (istat > 0 ? true : false);
		else {
			sprintf(msg, "failure to get the status of LOCAL_INTERACTION (3rd) from %s", buffer); mdpro.release(); in.close(); return false;
		}
		NextSection(&bf, ::seperator);

		if (bf != NULL && sscanf(bf, "%d", &istat) == 1) mdpro.m[iproc].cvar.bScaleForce = (istat > 0 ? true : false);
		else {
			sprintf(msg, "failure to get the status of SCALE_FORCE (4th) from %s", buffer); mdpro.release(); in.close(); return false;
		}
		NextSection(&bf, ::seperator);

		if (bf != NULL && sscanf(bf, "%d", &istat) == 1) mdpro.m[iproc].nloops = istat;
		else {
			sprintf(msg, "failure to get the number of md_loops (5th) from %s", buffer); mdpro.release(); in.close(); return false;
		}
		NextSection(&bf, ::seperator);

		if (bf != NULL && sscanf(bf, "%f", &eThermal) == 1) mdpro.m[iproc].eThermal = eThermal;
		else {
			mdpro.m[iproc].eThermal = 0;
			//sprintf(msg, "thermal energy is not defined.", buffer); mdpro.release(); in.close(); return false;
		}

		if (bf != NULL) {
			NextSection(&bf, ::seperator);

			if (bf != NULL && sscanf(bf, "%s", str) == 1 && strlen(str) > 0) strcpy(mdpro.m[iproc].parfname_IConstraint, str);
			else {
				strcpy(mdpro.m[iproc].parfname_IConstraint, "\0");
				//sprintf(msg, "thermal energy is not defined.", buffer); mdpro.release(); in.close(); return false;
			}
		}
	}
	in.close(); return true;
}

bool read_mc_species(char *fname, MC_SPECIES &sp) {
	sp.sp.release(); sp.mindx.release();
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		show_log("failure to open file : ", false); show_log(fname, true); return false;
	}
	char buffer[256] = "\0", *bf = NULL, msg[2560] = "\0", mol[256] = "\0";
	int nm = 0, im = 0;
	int istat = 0;
	if (!search(in, "[MC]", buffer)) {
		show_log("[MC] is not defined in : ", false); show_log(fname, true); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nm) != 1) {
		show_log("number of species to be randomly relocated is not defined given in ", false); show_log(fname, false); 
		show_log(", see: ", false); show_log(buffer, true); in.close(); return false;
	}
	if (nm <= 0) {
		show_log("no species to be relocated randomly", true); in.close(); return true;
	}
	sp.sp.set_array(nm);
	bf = buffer; NextSection(&bf, seperator);
	for (im = 0; im < nm; im++) {
		if (bf != NULL && sscanf(bf, "%s", mol) == 1 && strlen(mol) > 0) strcpy(sp.sp.m[im].name, mol);
		else {
			sprintf(msg, "failure to get the molecule %d from %s", im, buffer); sp.sp.release(); in.close(); return false;
		}
		NextSection(&bf, ::seperator);
	}

	in.getline(buffer, 250);
	if (sscanf(buffer, "%f %f", sp.xf + 2, sp.xt + 2) != 2) {
		show_log("failure to obtain z-range for molecular relocation [format : zf  zt], see: ", false); show_log(buffer, true);
		//sp.sp.release(); 
		in.close(); return false;
	}
	in.close(); return true;
}

extern void reset_random_seeds();
void SLAB_MD_MCONFS(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, MC_SPECIES &mcsp, MD_PRO *mdpro, int npro) {
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

//#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
	//if (::rcut_Ewd > cmm->xd[0] / 3) ::rcut_Ewd = cmm->xd[0] / 3;
	//if (::rcut_Ewd > cmm->xd[1] / 3) ::rcut_Ewd = cmm->xd[1] / 3;
	//if (::rcut_Ewd > cmm->xd[2] / 3) ::rcut_Ewd = cmm->xd[2] / 3;

	bool bStorageChain = false;
	InitMD(mmcell, cmm, bStorageChain);

	Slab_MDVar<MMOL_MD_CELL, BASIC_CLUSTER> mdvar;
	
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

	double kappa = mdvar.spme_var.spme_3d_var.esv1.kappa;
	sprintf(errmsg, "Ewald-Sum kappa = %f, rcut = %f", kappa, ::rcut_Ewd); show_log(errmsg, true);
//#endif

	init_mc_species(mmcell, mcsp);
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
		relocate_mc_species(mmcell, mcsp);
		molecule_check_2d_periodic_cell(mmcell);


                init_velocity(mmcell); // new velocity

		sprintf(oftitle, "%s.conf%d", md_title, iconf);
		sprintf(msg, "\nCONFIG #%d, output-file-title : %s", iconf, oftitle); show_log(msg, true);
		for (ipro = 0; ipro < npro; ipro++) {
			mdpro[ipro].cvar.bSave_mu = false; mdpro[ipro].cvar.bSavePz = false; mdpro[ipro].cvar.bSaveVirial = false; 
			sprintf(toftitle, "%s.pro%d", oftitle, ipro); 
			mmcell.msave->set_title(toftitle);
			sprintf(msg, "PROC #%d, output-file-title : %s", ipro, toftitle); show_log(msg, true);

			if (strlen(mdpro[ipro].parfname_IConstraint) > 0) {
				bConstraintChanged = true;
				strcpy(::parfname_IConstraint, mdpro[ipro].parfname_IConstraint);
				if (init_clusterIC_DB(2)) {
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
#ifdef _USE_GPU_
			GPU_MDRun_SLAB(mmcell, cmm, mdvar, mdpro[ipro].cvar, long(mdpro[ipro].nloops) * long(mmcell.msave->nsave), (ipro == 0 ? true : false));
#else
			MDRun_SLAB(mmcell, cmm, mdvar, mdpro[ipro].cvar, long(mdpro[ipro].nloops) * long(mmcell.msave->nsave), (ipro == 0 ? true : false));
#endif
			mmcell.msave->flush_buf(); mmcell.msave->close();

			if (COMMAND_MD_STOP) break; // stop with given command
		}

		::Ek0Internal = eThermal_old; mmcell.check_freedom();
		sprintf(msg, "thermal energy of each freedom : %5.2f kT", mmcell.Ethermal); show_log(msg, true);
		if (bConstraintChanged) {
			strcpy(::parfname_IConstraint, oldparIConstraint);
			if (init_clusterIC_DB(2)) {
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
		sprintf(msg, "FINAL PROC, output-file-title : %s", oftitle); show_log(msg, true);

		if (npro == 0) init_velocity(mmcell); // new velocity
#ifdef _USE_GPU_
		GPU_MDRun_SLAB(mmcell, cmm, mdvar, mdpro[ipro].cvar, long(mdpro[ipro].nloops) * long(mmcell.msave->nsave), (npro == 0 ? true : false));
#else
		MDRun_SLAB(mmcell, cmm, mdvar, mdpro0.cvar, mdpro0.nloops, (npro == 0 ? true : false));
#endif
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


