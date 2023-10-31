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
#include "spme_interact.h"
extern BSplineFunc<2> bsp;

#include "interface.h"
using namespace _interface_constraint_;
extern void InterfaceConstraints_CMM(CMM_CELL3D<BASIC_CLUSTER> *cmm, int nThreads, double &Ut);

extern void ExplicitEInteract(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, EInteractVar& evar, InteractRes &res);
extern void ZeroSpatialMomentum(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 6> *mp, int nJob);
extern void MDCELL_CalMolTransKineticEnergy(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4> *Ekt);

// MD simulation for PMF 
// the procedure is slightly changed from the MD_PROC in md_cell.cpp
// with additional change relating to the biased potential, and the constraints on the related atoms

#include "pmf.h"


bool PMF_constraint_consistent(MMOL_MD_CELL &mmcell, _pmf_constraint_::PMF_CONSTRAINT &cpmf) {
	using namespace _pmf_constraint_;
	char msg[256] = "\0";
	MATOM *patom1 = NULL, *patom2 = NULL;
	double Rmax = cpmf.r0 + cpmf.dr_max;
	double dcell[3];
	dcell[0] = mmcell.h[0] * 2; dcell[1] = mmcell.h[1] * 2; dcell[2] = mmcell.h[2] * 2; 

	bool status = true;

	if (cpmf.cModel == PMF_CONSTRAINT_0 || cpmf.cModel == PMF_CONSTRAINT_1 || cpmf.cModel == PMF_CONSTRAINT_2 || cpmf.cModel == PMF_CONSTRAINT_3) {
		if (cpmf.v_int.m[I0] < 0 || cpmf.v_int.m[I0] >= mmcell.pm.n) {
			sprintf(msg, "The cell has PM molecule [0 - %d], PM molecule [%d] relating to PMF-contraint is out of range.", mmcell.pm.n - 1, cpmf.v_int.m[I0]);
			show_log(msg, true); status = false;
		}
		if (cpmf.v_int.m[I1] < 0 || cpmf.v_int.m[I1] >= mmcell.pm.n) {
			sprintf(msg, "The cell has PM molecule [0 - %d], PM molecule [%d] relating to PMF-contraint is out of range.", mmcell.pm.n - 1, cpmf.v_int.m[I1]);
			show_log(msg, true); status = false;
		}

		if (!status) return false;

		patom1 = mmcell.pm.m[cpmf.v_int.m[I0]].c->atom;
		patom2 = mmcell.pm.m[cpmf.v_int.m[I1]].c->atom;
		
		if (Rmax >= mmcell.h[0] || Rmax >= mmcell.h[1] || Rmax >= mmcell.h[2]) {
			sprintf(msg, "The distance between contrainted-molecules [%f Angs.] is more than half of the cell dimension [%f, %f, %f].", Rmax, dcell[0], dcell[1], dcell[2]);
			show_log(msg, true); status = false;
		}
	}
	return status;
}

void PMF_constraint_force(MMOL_MD_CELL &mmcell, _pmf_constraint_::PMF_CONSTRAINT &cpmf, double &U, ofstream &out, bool bsave) {
	using namespace _pmf_constraint_;
	MATOM *patom1 = NULL, *patom2 = NULL;
	double R, R2;
	double Rmax = cpmf.r0 + cpmf.dr_max;
	VECTOR3 dr0, dr, u;
	double dcell[3];
	dcell[0] = mmcell.h[0] * 2; dcell[1] = mmcell.h[1] * 2; dcell[2] = mmcell.h[2] * 2; 

	double ft, f[3];
	U = 0;

	if (cpmf.cModel == PMF_CONSTRAINT_0 || cpmf.cModel == PMF_CONSTRAINT_1 || cpmf.cModel == PMF_CONSTRAINT_2 || cpmf.cModel == PMF_CONSTRAINT_3) {
		patom1 = mmcell.pm.m[cpmf.v_int.m[I0]].c->atom;
		patom2 = mmcell.pm.m[cpmf.v_int.m[I1]].c->atom;
		VECT3(patom1->rg, patom2->rg, dr0)

		// in this models, we believe that the cell dimensions has its half dimension bigger than Rmax
		if (dr0.v[0] > mmcell.h[0]) dr.v[0] = dr0.v[0] - dcell[0];
		else if (dr0.v[0] < -mmcell.h[0]) dr.v[0] = dr0.v[0] + dcell[0];
		else dr.v[0] = dr0.v[0];

		if (dr0.v[1] > mmcell.h[1]) dr.v[1] = dr0.v[1] - dcell[1];
		else if (dr0.v[1] < -mmcell.h[1]) dr.v[1] = dr0.v[1] + dcell[1];
		else dr.v[1] = dr0.v[1];

		if (dr0.v[2] > mmcell.h[2]) dr.v[2] = dr0.v[2] - dcell[2];
		else if (dr0.v[2] < -mmcell.h[2]) dr.v[2] = dr0.v[2] + dcell[2];
		else dr.v[2] = dr0.v[2];

		V3ABS2(dr, R2) if (R2 == 0) return; // ignore it
		R = sqrt(R2); u.v[0] = dr.v[0] / R; u.v[1] = dr.v[1] / R; u.v[2] = dr.v[2] / R;
		cpmf.ef(R, U, ft);
		f[0] = ft * u.v[0]; f[1] = ft * u.v[1]; f[2] = ft * u.v[2];
		patom2->F.v[0] += f[0]; patom2->F.v[1] += f[1]; patom2->F.v[2] += f[2];
		patom1->F.v[0] -= f[0]; patom1->F.v[1] -= f[1]; patom1->F.v[2] -= f[2];

		if (bsave) {
			out<<R<<"  "<<U;
		}
	}
}

void PMF_constraint_force2(MMOL_MD_CELL &mmcell, _pmf_constraint_::PMF_CONSTRAINT &cpmf) {
	using namespace _pmf_constraint_;
	PMOLECULE *pm1, *pm2;
	MATOM *patom1 = NULL, *patom2 = NULL;

	if (cpmf.cModel == PMF_CONSTRAINT_2 || cpmf.cModel == PMF_CONSTRAINT_3) { // with one PM (1st) fixed at a position
		pm1 = mmcell.pm.m + cpmf.v_int.m[I0];
		pm2 = mmcell.pm.m + cpmf.v_int.m[I1];
		patom1 = pm1->c->atom;
		patom2 = pm2->c->atom;

		V3zero(patom1->F) V3zero(pm1->c->dyn->fc) V3zero(pm1->c->dyn->f)
	}
}

void PMF_constraint_kinetic(MMOL_MD_CELL &mmcell, _pmf_constraint_::PMF_CONSTRAINT &cpmf) {
	using namespace _pmf_constraint_;
	PMOLECULE *pm1, *pm2;
	MATOM *patom1 = NULL, *patom2 = NULL;

	if (cpmf.cModel == PMF_CONSTRAINT_2 || cpmf.cModel == PMF_CONSTRAINT_3) { // with one PM (1st) fixed at a position
		pm1 = mmcell.pm.m + cpmf.v_int.m[I0];
		pm2 = mmcell.pm.m + cpmf.v_int.m[I1];
		patom1 = pm1->c->atom;
		patom2 = pm2->c->atom;

		V3zero(pm1->v) V3zero(pm1->alpha)
	}
}

#ifndef _DISABLE_
extern char parfname_IConstraint[256];
void PMF_MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, char *pmf_parfname) {
	// *********************************************************
	// working on periodical cell !!
	// *********************************************************

	if (strlen(parfname_IConstraint) > 0) {
		if (init_clusterIC_DB(2, cmm->xl[2], cmm->xr[2])) {
			sprintf(errmsg, "Interface defined in %s is enabled!", parfname_IConstraint); show_log(errmsg, true);
		}
	}

	using namespace _pmf_constraint_;
	// PMF constraints
	PMF_CONSTRAINT pmf_constraint;
	if (!read_pmf_constraint(pmf_parfname, pmf_constraint)) return;
	if (!PMF_constraint_consistent(mmcell, pmf_constraint)) return; // check the constraint is reasonable or not
	pmf_constraint.construct_constraint_efprofile();

	double U_pmf = 0;
	
	ITER_SAVE pmf_save;

	// end of PMF constraint

	using namespace _spme_;

	if (mmcell.mm.n == 0 && mmcell.sm.n == 0 && mmcell.pm.n == 0) {
		sprintf(errmsg, "no molecule is defined in the cell"); show_msg(errmsg); return;
	}
	// we need to make sure the mass center of the cell is at zero
	VECTOR3 mass_center;
	set_free_cell_mass_center(mmcell, mass_center);

	molecule_check_periodic_cell(mmcell); // check the molecule mass center is inside the cell

	// here we have to check the periodity of the cluster position
	//if (md_mode == MULTI_MM_PERIDIC_CELL_MD) cluster_check_periodic_cell(mmcell);
	// *****************  over  ******************
	show_all_molecules();

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

	init_BSplineFunc<2>(::bsp, ::indx_bspline);
	VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	vspme_q.bST = ::bVirial; vspme_q.bST_diagonal = true;
	VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	vspme_mu.bST = ::bVirial; vspme_mu.bST_diagonal = true;
	if (bPolarize) {
		mmcell.init_polarization_buff();
		mmcell.init_dipole_hist();
		mmcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		vspme_mu.bsp = &(::bsp);
		init_spme(&mmcell, vspme_mu); // init the viral atoms
		vspme_mu.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_mu.xl[0] = cmm->xl[0]; vspme_mu.xl[1] = cmm->xl[1]; vspme_mu.xl[2] = cmm->xl[2];
		vspme_mu.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		vspme_mu.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();
	}
	else { // charge only
		mmcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		vspme_q.bsp = &(::bsp);
		init_spme(&mmcell, vspme_q); // init the viral atoms
		vspme_q.init(::bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_q.xl[0] = cmm->xl[0]; vspme_q.xl[1] = cmm->xl[1]; vspme_q.xl[2] = cmm->xl[2];
		vspme_q.cell_dim(cmm->xd[0], cmm->xd[1], cmm->xd[2]); // important, after init();
		vspme_q.init_k(/*cmm->xd[0], cmm->xd[1], cmm->xd[2]*/);
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}

extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
#endif

	int i = 0;

	bool bStorageChain = false;
	InitMD(mmcell, cmm, bStorageChain);
	// freedom of cell changing from the constraint
	if (pmf_constraint.cModel == PMF_CONSTRAINT_2 || pmf_constraint.cModel == PMF_CONSTRAINT_3) {
		mmcell.NF -= 3; // one PM has coordinates fixed
		mmcell.Ethermal -= 3 * ::Ek0Internal;
	}

	long nloop = 0;
	int nc = 0, nm = 0, ncell = 0;
	double Ek = 0, Ep = 0, Ek_t = 0, Ep_t = 0, Etorsion = 0, U_LJ = 0, U_estat = 0;

	InteractVar LJ_var; 
	LJ_var.bVirial = ::bVirial; LJ_var.bST_diagonal = true;

	SPME_VAR spme_var;
	spme_var.bVirial = ::bVirial; spme_var.bST_diagonal = true; spme_var.bEfieldOnly = false; 
	spme_var.esv1.bEwald = true; spme_var.esv1.init_EwaldSum(cmm->xd[0] * cmm->xd[1] * cmm->xd[2], ::rcut_Ewd);
	spme_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;

	InteractRes LJ_res, spme_res;

	TorsionVar torsion_var[MAX_THREADS];
	HingeConstraintVar hinge_constraint_var[MAX_THREADS];
	InteractRes mires[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		torsion_var[i].bVirial = ::bVirial; torsion_var[i].bST_diagonal = true;
		hinge_constraint_var[i].bVirial = ::bVirial; hinge_constraint_var[i].bST_diagonal = true;
	}

	double Uext = 0; // external electric field

	double Pint = 0;
	SVECTOR<double, 4> virial;

	MMOLECULE *mm = NULL;
	//CLUSTER *pc = NULL;
	//SMOLECULE *sm = NULL;
	//PMOLECULE *pm = NULL;

	int nMDThreads = MAX_THREADS;
	MDCELL_Job<MMOL_MD_CELL> mdcellJob[MAX_THREADS];
	AssignJob_MDCELL<MMOL_MD_CELL>(&mmcell, mdcellJob, nMDThreads);
	sprintf(errmsg, "Job # : %d; working with molecules:", MAX_THREADS); show_log(errmsg, true);
	for (i = 0; i < MAX_THREADS; i++) {
		sprintf(errmsg, "%d -- MM %d, SM %d, PM %d", i, mdcellJob[i].mmIndx.n, mdcellJob[i].smIndx.n, mdcellJob[i].pmIndx.n);
		show_log(errmsg, true);
		//for (nm = 0; nm < mdcellJob[i].mmIndx.n; nm++) {
		//	sprintf(errmsg, "%d ", mdcellJob[i].mmIndx.m[nm]); show_log(errmsg, false);
		//}
		//show_log("", true);
	}

	MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> mdcellHDynThreadVar[MAX_THREADS];
	AssignJob_MDCELL_HDYN<MMOL_MD_CELL>(&mmcell, mdcellHDynThreadVar, nMDThreads, false);

	// assuming each cluster has minimum dimension 2.5^3, and we expand the array dimension for each subcell 50% bigger than estimation
	float vcluster = 1.5 * 1.5 * 1.5;
	int ncs = int(cmm->fw[0] * cmm->fw[1] * cmm->fw[2] / vcluster * 1.2 + 0.5);
	if (ncs < 2) ncs = 2;
	CMM_array_set_storage<BASIC_CLUSTER>(cmm, ncs);
	CMM_CELL3D<BASIC_CLUSTER> job_cmm[MAX_THREADS];

	MDCELL_ThreadJob<MMOL_MD_CELL> job_mdcell[MAX_THREADS];
	MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >job_mdcellCMM[MAX_THREADS];

	double djob[MAX_THREADS];
	MATRIX<3> mjob[MAX_THREADS];
	bool status_job[MAX_THREADS];
	//INTERACT interact_job[MAX_THREADS]; // using for explicit interaction calculation
	SVECTOR<double, 4> sv4job[MAX_THREADS];
	SVECTOR<double, 3> sv3job[MAX_THREADS];

	for (i = 0; i < MAX_THREADS; i++) {
		if (MAX_THREADS > 1) {
			CMM_array_construct_for_storage<BASIC_CLUSTER>(job_cmm + i, cmm, ncs);
		}

		job_mdcell[i].set_pars(i, mdcellJob + i);
		job_mdcellCMM[i].set_pars(i, mdcellJob + i);
		job_mdcellCMM[i].cmm = cmm;
	}

	// make sure the spatial momentum of the whole cell is zero
	if (mmcell.bHDynMC_T) {
		MDCELL_mass_center(mmcell);
		ZeroSpatialMomentum(job_mdcell, sv3job, MAX_THREADS);
	}

	PMF_constraint_kinetic(mmcell, pmf_constraint);

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

	long LOOPS = 0;
	int save_indx = 0;
	MD_SAVE *mdsave = NULL;
	if (mmcell.mm.n > 0) {LOOPS = mmcell.lfmd.m[0].LOOPS; mdsave = mmcell.lfmd.m[0].mmsave.mdsave;}
	else if (mmcell.sm.n > 0) {LOOPS = mmcell.lfsmd.m[0].LOOPS; mdsave = mmcell.lfsmd.m[0].smsave.mdsave;}
	else if (mmcell.pm.m > 0) {LOOPS = mmcell.lfpmd.m[0].LOOPS; mdsave = mmcell.lfpmd.m[0].pmsave.mdsave;}

	bool bSaveVirial = false;
	ITER_SAVE virial_save;
	if (::bVirial) {
		virial_save.set_file(mdsave->title, "virial");
		virial_save.set_iter(mdsave->nsave, mdsave->max_save);
	}

	pmf_save.set_file(mdsave->title, "pmf");
	pmf_save.set_iter(mdsave->nsave, mdsave->max_save);

	int ncs_mmcell = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) ncs_mmcell += mmcell.mm.m[nm].nCluster;
	ncs_mmcell += mmcell.sm.n + mmcell.pm.n;

	_rand_::RAND_ARRAY<int> randa_m, randa_c;

	while (nloop < LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command

		iCheckCMM = 0;

		if (iCheckCMM == 0) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM && iCheckCMM == 0) { // we do re-position of mmcell coupling with CMM position check, relationship check
			// here we have to check the periodity of the molecule position
			//molecule_check_periodic_cell(mmcell);
			MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_molecule_PeriodicCell_check), job_mdcell, nMDThreads);
		}

		check_atom_rg(mmcell, nMDThreads);

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalInertiaMoment), job_mdcell, nMDThreads);
		
		//MDCellMassCenter(mdvar.job_mdcell, sv6job, nMDThreads);
		if (mmcell.bHDynMC_T) {
			MDCELL_mass_center(mmcell);
		}

		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_InitClusterForce), job_mdcell, nMDThreads);

	// re-distribute / check clusters in CMM, use vp of each cluster only
	#if _CHECK_TIME_STAMP_ == 1
		time.start();
	#endif
		if (bStorageChain) {
			if (nloop == 0) cmm_check_cluster<BASIC_CLUSTER>(*cmm);
			else if (bCheckClusterInCMM) cmm_check<BASIC_CLUSTER>(*cmm, MAX_THREADS);
		}
		else {
			if (nloop == 0 || bCheckClusterInCMM) {
				if (MAX_THREADS == 1) {
					job_mdcellCMM[0].cmm = cmm;
					cmm->reset_cell_chain();
					MMCELL_CMM_cluster_allocate(job_mdcellCMM, status_job);
				}
				else {
					for (i = 0; i < MAX_THREADS; i++) job_mdcellCMM[i].cmm = job_cmm + i; // set the JOB-related CMM is the job-related CMM
					MTOperate1<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, bool>((void*)(&MMCELL_CMM_cluster_allocate), job_mdcellCMM, nMDThreads, status_job);
					for (i = 0; i < nMDThreads; i++) {
						if (!status_job[i]) {
							sprintf(errmsg, "loop %d thread %d: buffer of CMM cell leaks when allocating cluster", nloop, i);
							show_log(errmsg, true);
						}
					}
					if (!combine_AtomAllocatedCMM_array(job_cmm, nMDThreads, cmm, true, true)) {
						sprintf(errmsg, "loop %d : buffer of CMM cell leaks when combining allocated cluster", nloop);
						show_log(errmsg, true);
					}
				}
				for (i = 0; i < MAX_THREADS; i++) job_mdcellCMM[i].cmm = cmm; // reset the JOB-related CMM is the main CMM
			}

			if (bCheckClusterInCMM) {
			{
				nc = 0;
				for (i = 0; i < cmm->acell.n; i++) {
					nc += cmm->acell.m[i]->length();
					//sprintf(errmsg, "%d %d", i, cmm->acell.m[i]->length()); show_log(errmsg, true);
				}
				if (nc != ncs_mmcell) {
					sprintf(errmsg, "%d : cluster in CMM -- %d, real -- %d", nloop, nc, ncs); show_log(errmsg, true);
				}
			}
			}
			
		}
	#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "cmm_check takes %d ms", time.glance());
		show_infor(errmsg);
	#endif

		Ep = 0;
		// check relationship of each cluster, electrostatic and LJ
		if (iCheckCMM == 0) {
		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
			if (bStorageChain) {
				ParallelCheckRelshipInCells(*cmm, 0, cmm->acell.n - 1, MAX_THREADS); // setup relationship for each atom with CMM
			}
			else {
				MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >((void*)(&MMCELL_CMM_CheckRelship), job_mdcellCMM, nMDThreads);
			}
		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
			show_infor(errmsg);
		#endif
		}

		// important : calculate LJ interaction first, here cluster force & torque are reset
		
		LJ_Interact(&mmcell, cmm, MAX_THREADS, LJ_var, LJ_res);
		U_LJ = LJ_res.U;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!::local_interact && !mmcell.eNeutral) {
			spme_var.bEfieldOnly = false;
			// important : calculate Coulomb secondly, here cluster force & torque are accumulated
			if (::bPolarize) {
				if (nloop > 3) Guess_induced_dipole(mmcell, MAX_THREADS);
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				if (!Polarize_SPME_Interact(&mmcell, cmm, vspme_mu, MAX_THREADS, true, spme_var, spme_res)) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
				}
				Backup_induced_dipole_hist(mmcell, MAX_THREADS);
			}
			else {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				SPME_Interact(&mmcell, cmm, vspme_q, MAX_THREADS, spme_var, spme_res);
			}

			U_estat = spme_res.U;
		}
#endif
		
		Ep = U_LJ + U_estat;
		if (::bVirial) {
			virial.v[0] = LJ_res.STxx + spme_res.STxx;
			virial.v[1] = LJ_res.STyy + spme_res.STyy;
			virial.v[2] = LJ_res.STzz + spme_res.STzz;
			virial.v[3] = LJ_res.STtr + spme_res.STtr;

			if (!virial_save.one_more()) break;
			else bSaveVirial = virial_save.bsave();
		}
		
		//virial.v[0] = 0; virial.v[1] = 0; virial.v[2] = 0; virial.v[3] = 0; Ep = 0;

		if (mmcell.mm.n > 0) {
			MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, TorsionVar, InteractRes>((void*)(&MMCELL_MM_Torsion), job_mdcell, nMDThreads, torsion_var, mires);
			Etorsion = 0;
			for (i = 0; i < nMDThreads; i++) {
				Ep += mires[i].U; Etorsion += mires[i].U;
				if (torsion_var[i].bVirial) {
					virial.v[0] += mires[i].STxx;
					virial.v[1] += mires[i].STyy;
					virial.v[2] += mires[i].STzz;
					virial.v[3] += mires[i].STtr;
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

		// PMF constraints
		pmf_save.one_more();
		if (pmf_save.bsave()) (*pmf_save.out)<<nloop<<"  ";
		PMF_constraint_force(mmcell, pmf_constraint, U_pmf, *pmf_save.out, pmf_save.bsave());
		if (pmf_save.bsave()) (*pmf_save.out)<<endl;
		Ep += U_pmf;

		//if (mmcell.mm.m != NULL) MOperate<MMOLECULE>(MM_TOTAL_FORCE, mmcell.mm.m, mmcell.mm.n, MAX_THREADS);
		//if (mmcell.sm.m != NULL) MOperate<SMOLECULE>(SM_TOTAL_FORCE, mmcell.sm.m, mmcell.sm.n, MAX_THREADS);
		//if (mmcell.pm.m != NULL) MOperate<PMOLECULE>(PM_TOTAL_FORCE, mmcell.pm.m, mmcell.pm.n, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ClusterForce), job_mdcell, nMDThreads);

		if (::bVirial && bSaveVirial && (mmcell.mm.n > 0 || mmcell.sm.n > 0)) {
			for (i = 0; i < nMDThreads; i++) {
				sv4job[i].v[0] = 0; sv4job[i].v[1] = 0; sv4job[i].v[2] = 0; sv4job[i].v[3] = 0;
			}
			switch (::method_Virial) {
			case CLUSTER_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_VirialCorrect_RigidCluster), job_mdcell, nMDThreads, sv4job);
				break;
			case MOLECULE_VIRIAL:
				MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_VirialCorrect_Molecule), job_mdcell, nMDThreads, sv4job);
				break;
			}
			for (i = 0; i < nMDThreads; i++) {
				virial.v[0] -= sv4job[i].v[0];
				virial.v[1] -= sv4job[i].v[1];
				virial.v[2] -= sv4job[i].v[2];
				virial.v[3] -= sv4job[i].v[0] + sv4job[i].v[1] + sv4job[i].v[2];
			}
		}
		
		if (::bVirial && bSaveVirial) {
			(*virial_save.out)<<nloop<<"  "<<virial.v[0] * ::unit_Ek_kT<<"  "<<virial.v[1] * ::unit_Ek_kT<<"  "<<virial.v[2] * ::unit_Ek_kT<<"  "<<virial.v[3] * ::unit_Ek_kT;
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
			MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >((void*)(&MMCELL_CMM_CheckImpSolvRelship), job_mdcellCMM, nMDThreads);
			MTOperate1< MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> >, double>((void*)(&MMCELL_cluster_impsolv_interact), job_mdcellCMM, nMDThreads, djob);
			for (i = 0; i < nMDThreads; i++) Ep += djob[i];
		}

		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ScaleClusterForce), mdvar.job_mdcell, nMDThreads);

		if (bRandomForce) gen_random_force(mmcell, randa_m, randa_c);

		if (clusterIC_db.n > 0) {
			// the force from interface, and the torque
			// however, the force will not contribute to the virial
			InterfaceConstraints_CMM(cmm, nMDThreads, djob[0]); Ep += djob[0];
		}

		iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

		PMF_constraint_force2(mmcell, pmf_constraint); // final constraint at the force, before the dynamic calculation
		V6zero(mmcell.pm.m[2].c->dyn->fc) V6zero(mmcell.pm.m[2].c->dyn->f)

		// self-consistent velocity and acceleration, Nose-Hoover Chain
		bSpeedVerlet = (nloop == 0 ? false : true);
		if (!LeapFrogVerlet_SelfConsistent_SpeedAccel(&mmcell, job_mdcell, mdcellHDynThreadVar, bSpeedVerlet, nMDThreads)) {
			sprintf(errmsg, "failure to get consistent velocity and acceleration at loop %d", nloop);
			show_log(errmsg, true);
		}

		V3zero(mmcell.lfpmd.m[2].v0) V3zero(mmcell.lfpmd.m[2].v0_mh)
		V3zero(mmcell.pm.m[2].v) V3zero(mmcell.pm.m[2].alpha)

		Ep_t = Ep; Ek_t = 0;
		for (nm = 0; nm < mmcell.mm.n; nm++) Ek_t += mmcell.mm.m[nm].Ek;
		for (nm = 0; nm < mmcell.sm.n; nm++) Ek_t += mmcell.sm.m[nm].Ek;
		for (nm = 0; nm < mmcell.pm.n; nm++) Ek_t += mmcell.pm.m[nm].Ek;

		show_all_molecules();
		save_indx = mdsave->mESave.save_indx();
		show_loop_infor(nloop, save_indx, LOOPS - nloop, (float)Ek_t, (float)Ep_t);
		mmcell.MDSave(nloop, (float)Ep_t, false);

		//if (::bPolarize && mdsave->mKinSave.bsave()) log_dipole(&mmcell); // write the dipole of cluster into log file

		if (::bVirial && bSaveVirial && ::method_Virial == CLUSTER_VIRIAL) {
			if (mmcell.mm.n > 0 && mmcell.bHingeConstraint) {
				MTOperate2< MDCELL_ThreadJob<MMOL_MD_CELL>, HingeConstraintVar, InteractRes>((void*)(&MDCELL_ClusterVirialCorrect_HingeConstraint), job_mdcell, nMDThreads, hinge_constraint_var, mires);
				for (i = 0; i < nMDThreads; i++) {
					virial.v[0] += mires[i].STxx;
					virial.v[1] += mires[i].STyy;
					virial.v[2] += mires[i].STzz;
					virial.v[3] += mires[i].STtr;
				}
			}

			MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 4> >((void*)(&MDCELL_CalMolTransKineticEnergy), job_mdcell, nMDThreads, sv4job);
			Ek_t = 0;
			for (i = 0; i < nMDThreads; i++) Ek_t += sv4job[i].v[3];
			Pint = (virial.v[3] + Ek_t * 2 / unit_Ek_kT) / (::fUnit_atm_u * mmcell.h[0] * mmcell.h[1] * mmcell.h[2] * 8 * 3);
			if (bSaveVirial) {
				(*virial_save.out)<<"          "<<virial.v[3] * ::unit_Ek_kT<<"  "<<Ek_t<<"  "<<Pint<<endl;
			}
		}

		// before verlet movement, we could need to add more constraint for PMF
		PMF_constraint_kinetic(mmcell, pmf_constraint);

		// Verlet move, and prepare the velocity for next step
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_Verlet), job_mdcell, nMDThreads);
		// accept current NHC vars and prepair for next step
		if (mmcell.bIndependentFrameHDyn) {
			MDCELL_HDyn_MultiThreadOperate<MMOL_MD_CELL>((void*)(&MDCELL_AcceptCurrentDynForNextStep), mdcellHDynThreadVar, nMDThreads);
		}
		else {
			for (i = 0; i < mmcell.hdyn.n; i++) mmcell.hdyn.m[i].nhc->accept_nextstep();
		}
		if (mmcell.bHDynMC_T) mmcell.mcHDyn.nhc->accept_nextstep();

		#if _CHECK_TIME_STAMP_ == 1
			time.start();
		#endif
		// make sure torsion angle in range [-PI, PI]
		if (iCheckAngle == 0) {
			if (mmcell.mm.n > 2) {
				MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CheckRotAngle), job_mdcell, nMDThreads);
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

		nloop++;
	}
	if (MAX_THREADS > 1) {
		for (i = 0; i < MAX_THREADS; i++) job_cmm[i].reset_basic_cell_matrix();
	}
	if (::bVirial) virial_save.reset();
}

#endif // _DISALBE_
