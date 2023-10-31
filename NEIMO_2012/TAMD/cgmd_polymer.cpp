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
#include "cg-mm.h"
#include "CMM.h"
using namespace _cmm_3d_;
#include "cg-cmm.h"
#include "cluster.h"
#include "Interaction.h"
#include "Interact2.h"
#include "NEIMO.h"
#include "MD.h"
#include "cg-md.h"
//#include "MD_POLYMER.h"
#include "cgmd_polymer.h"
#include "Interact4.h"

#include "var.h"
#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

extern bool trace;

extern int delay_time; // in milisecond
extern int max_time_pickup_molstruct;  // in milisecond
extern int nshow;
extern int ishow;

#if _SYS_ == _WINDOWS_SYS_
extern CWnd* pMainFrame;
extern HANDLE hWaitGetMolStructEvent;
extern HANDLE hWaitEvent;
extern CMDInforDlg *mdInforDlg;
#endif


namespace _coarse_grain_ {
/*
void get_Frame_Hoover(CG_MMOLECULE &mm, LFVERLET_CGMM_MD& mdpar, double Et_thermal, double &tksi, double Er_thermal, double &rksi) {
	double T = mm.mKD->Ekt / Et_thermal;
	double vksi = 0;
	if (T > 3) {
		tksi = ::max_ksi_Hoover; // this will slow the speed to half
	}
	else {
		vksi = mdpar.inv_tao_ss * (T - 1);
		if (T > 1.5) vksi *= 1.1;
		tksi = mdpar.ksi_frame_trans + vksi * mdpar.dt;
		if (tksi > mdpar.max_ksi_Hoover) tksi = mdpar.max_ksi_Hoover;
	}

	T = mm.mKD->Ekr / Er_thermal;
	if (T > 3) rksi = ::max_ksi_Hoover; // this will slow the speed to half
	else {
		vksi = mdpar.inv_tao_ss * (T - 1);
		if (T > 1.5) vksi *= 1.1;
		rksi = mdpar.ksi_frame_rot + vksi * mdpar.dt;
		if (rksi > mdpar.max_ksi_Hoover) rksi = mdpar.max_ksi_Hoover;
	}
	
	return;
}
*/
//void get_Hoover(CG_MMOLECULE &mm, LFVERLET_CGMM_MD& mdpar, double Ethermal, double &ksi) {
void get_Hoover(CG_MMOLECULE &mm, LFVERLET_CGMM_MD& mdpar, double Ethermal) {
	double T = mm.Ek / Ethermal;
	/*
	double vksi = mdpar.inv_tao_ss * (T - 1);
	if (T > 1.5) vksi *= 1.1;
	ksi = mdpar.ksi_Hoover + vksi * mdpar.dt;
	if (FABS(ksi) > mdpar.max_ksi_Hoover) ksi = sign(ksi) * mdpar.max_ksi_Hoover;
	*/
	mdpar.hdyn->nhc->set_dT(T - 1);
	mdpar.hdyn->nhc->verlet_NHC();
	
	return;
}

//void get_Hoover(CGSM &sm, LFVERLET_CGSM_MD& mdpar, double Ethermal, double &ksi) {
void get_Hoover(CGSM &sm, LFVERLET_CGSM_MD& mdpar, double Ethermal) {
	double T = sm.Ek / Ethermal;
	/*
	double vksi = mdpar.inv_tao_ss * (T - 1);
	if (T > 1.5) vksi *= 1.1;
	ksi = mdpar.ksi_Hoover + vksi * mdpar.dt;
	if (FABS(ksi) > mdpar.max_ksi_Hoover) ksi = sign(ksi) * mdpar.max_ksi_Hoover;
	*/
	mdpar.hdyn->nhc->set_dT(T - 1);
	mdpar.hdyn->nhc->verlet_NHC();

	return;
}

void HooverCorrectionOfAccerlation(CG_MMOLECULE *mm, int nc1, int nc2, HOOVER_DYN& hDyn) {
	int nc = 0;
	CG_CLUSTER *pc = NULL;
	double mksi = hDyn.nhc->ksi.v[0];
	VECTOR3 vw, r;
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		pc->bkm->alpha.v[0] -= mksi * pc->bkm->v.v[0];
		pc->bkm->alpha.v[1] -= mksi * pc->bkm->v.v[1];
		pc->bkm->alpha.v[2] -= mksi * pc->bkm->v.v[2];

		/*
		VECT3(mm->cm, (*(pc->r)), r) // r = mass_center ==> cluster
		V3PV3(mm->mKD->w, r, vw)
		*/
	}
}

bool CGSMLeapFrog_Verlet(CGSM &sm, LFVERLET_CGSM_MD& mdpar, double Ethermal, bool init_md, char *show_msg, char *log_msg, int n_msglen) {
	char msg[512] = "\0";

	int i = 0, iCycle = 0;

	double vConvg = 0.005 / mdpar.dt; //0.001 / mdpar.dt;
	if (vConvg < 2e-4) vConvg = 2e-4;
	VECTOR3 vbf, dv;
	double dvmax = 0;

	double ksi_Hoover = 0;

	while (iCycle < 10) {
		V32V3(sm.bkm->v, vbf)
		// calculate the kinetic energy, which will be used for Hoover acceleration estimation
		calCGSMKineticEnergy(sm); 

		get_Hoover(sm, mdpar, Ethermal); 
		ksi_Hoover = mdpar.hdyn->nhc->ksi.v[0];

		// calculate Brownian force here

		// calculate the acceleration
		calCGSMDynamic(sm);
		sm.bkm->alpha.v[0] -= ksi_Hoover * sm.bkm->v.v[0];
		sm.bkm->alpha.v[1] -= ksi_Hoover * sm.bkm->v.v[1];
		sm.bkm->alpha.v[2] -= ksi_Hoover * sm.bkm->v.v[2];

		// get new speed at current position and +1/2 time step
		if (init_md) Speed_Kinetic(mdpar, sm);
		else Speed_Verlet(mdpar, sm, true); // reset the speed in mm to the new speed

		VECT3(vbf, sm.bkm->v, dv)
		V3ABS2(dv, dvmax)
		dvmax = sqrt(dvmax);
		if (dvmax < vConvg) break;

		iCycle++;
	}

	Verlet(mdpar, sm); // go to next step

	mdpar.hdyn->nhc->accept_nextstep();

	return true;
}

bool CGLeapFrog_Verlet(CG_MMOLECULE &mm, LFVERLET_CGMM_MD& mdpar, VECTOR3 *vbf, double Ethermal, bool init_md, char *show_msg, char *log_msg, int n_msglen) {
	char msg[512] = "\0";

	int nc = 0, i = 0, iCycle = 0;
	CG_CLUSTER *pc = NULL;

	double vConvg = 0.005 / mdpar.dt; //0.001 / mdpar.dt;
	if (vConvg < 2e-4) vConvg = 0.05 / mdpar.dt;
	VECTOR3 dv;
	double dvmax = 0, dvabs = 0;

	//double ksi_Hoover = 0, ksi_frame_trans = 0, ksi_frame_rot = 0;
	NHC<MNHC> *nhc = mdpar.hdyn->nhc;

	while (iCycle < 10) {
		for (nc = 0; nc < mm.nCluster; nc++) {V32V3(mm.cluster[nc].bkm->v, vbf[nc])}

		// calculate the kinetic energy, which will be used for Hoover acceleration estimation
		ParallelCalKineticEnergy(mm); 
		get_Hoover(mm, mdpar, Ethermal * mm.nCluster);

		// calculate Brownian force here
		if (bBrownian) MacroMoleculeOperate<CG_MMOLECULE>((void*)(&CG_Brownian_force), mm, 0, mm.nCluster, MAX_THREADS);

		// calculate the acceleration
		ParallelCalDynamic(mm);
		MacroMoleculeOperate1<CG_MMOLECULE, HOOVER_DYN>((void*)(&HooverCorrectionOfAccerlation), mm, 0, mm.nCluster, *(mdpar.hdyn), MAX_THREADS);

		// get new speed at current position and +1/2 time step
		if (init_md) Speed_Kinetic(mdpar, mm);
		else Speed_Verlet(mdpar, mm, true); // reset the speed in mm to the new speed

		dvmax = 0;
		for (nc = 0; nc < mm.nCluster; nc++) {
			VECT3(vbf[nc], mm.cluster[nc].bkm->v, dv)
			V3ABS2(dv, dvabs)
			dvmax = (dvmax < dvabs ? dvabs : dvmax);
		}
		dvmax = sqrt(dvmax);
		if (dvmax < vConvg) break;

		iCycle++;
	}

	Verlet(mdpar, mm); // go to next step

	mdpar.hdyn->nhc->accept_nextstep();

	return true;
}

void ZeroVelocity(CG_MMOLECULE &mm) {
	int nc;
	CG_CLUSTER *pc = NULL;
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		V3zero(pc->bkm->v)
	}
}

void InitCGClusterVelocity(CG_MMOLECULE *mm, int nc1, int nc2, double &Ek) {
	double E_InertiaMoment = Ek * kT / U_InertiaMoment;
	int nc;
	double v = 0;
	CG_CLUSTER *pc = NULL;
	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		v = E_InertiaMoment / (*(pc->M)) * 2; // Ek = 1/2 * m * v^2
		v = sqrt(v);
		rand_vect3(pc->bkm->v.v);
		pc->bkm->v.v[0] *= v; pc->bkm->v.v[1] *= v; pc->bkm->v.v[2] *= v;
	}
}

void InitVelocity(CG_MMOLECULE &mm, double Ek) {
	MacroMoleculeOperate1<CG_MMOLECULE, double>((void*)(&InitCGClusterVelocity), mm, 0, mm.nCluster, Ek, MAX_THREADS);
}

void InitCGSMVelocity(CGSM *sm, int nSM, double &Ek) {
	double E_InertiaMoment = Ek * kT / U_InertiaMoment;
	int n;
	double v = 0;
	for (n = 0; n < nSM; n++) {
		v = E_InertiaMoment / (*(sm[n].M)) * 2; // Ek = 1/2 * m * v^2
		v = sqrt(v);
		rand_vect3(sm[n].bkm->v.v);
		sm[n].bkm->v.v[0] *= v; sm[n].bkm->v.v[1] *= v; sm[n].bkm->v.v[2] *= v;
	}
}

extern CGBOUND_DB cgbound_db;
extern CGBOUND_DIPOLE_DISTRBT_DB cg_dipole_distrbt_db;

void MD_POLYMER(CG_MMOLECULE &cgmm, LFVERLET_CGMM_MD& mdpar, CMM_CELL3D<CG_CLUSTER> *cmm) {
	char msg_show[5120] = "\0", msg_log[5120] = "\0";
	int n_msglen = 5120;

	long nloop = 0;
	int nc = 0, ncell = 0;
	double Ek = 0, Ep = 0;

	bool Overlap = false;

	CG_CLUSTER *pc = NULL;
	MD_SAVE *mdsave = mdpar.mmsave.mdsave;
	if (mdsave == NULL) {sprintf(errmsg, "no MD saving object is defined"); show_msg(errmsg); return;}
	sprintf(refresh_mol_xyz_struct_fname, "%s.xyz", mdpar.mmsave.mdsave->title);


	{
		// set mass center to {0, 0, 0}
		VECTOR<3> cm;
		cgmm.SetMassCenter(cm.v, true);
	}

	cgmm.calBrownianPars();

	// bound-dipole
	if (::bCGdipole) cgmm.setup_bounds();

	// kinetic & dynamic parameters
	cgmm.setup_kin_dyn();
	InitVelocity(cgmm, 3 * ::Ek0Internal);

	mdpar.resetClusters(cgmm.nCluster);
	mdpar.Init(cgmm); // copy initial torsion angle and speed to mdpar
	
#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	bool init_md = true;
	int iCheckCMM = 0;

	VECTOR3 *vbf = new VECTOR3[cgmm.nCluster];

	_rand_::RAND_ARRAY<int> randa;

	while (nloop < mdpar.LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command
		init_md = (nloop == 0 ? true : false);

		{
			// set mass center to {0, 0, 0}
			VECTOR<3> cm;
			cgmm.SetMassCenter(cm.v, true);
		}
		calInertiaMoment(&cgmm);
#if _CHECK_TIME_STAMP_ == 1
		time.start();
#endif
		// very important to reset the force on mm cluster, because the force on mm cluster is calculate at first
		// then transfer to cg_cluster
		for (nc = 0; nc < cgmm.nCluster; nc++) {
			V3zero(cgmm.cluster[nc].dyn->fc)
			cgmm.cluster[nc].dyn->overlap = false;
		}

		Ep = 0;
		switch (interact_cal_method) {
		case _CMM:
			if (iCheckCMM == 0 || nloop == 0) {
				//if (nloop == 0) calColumbMultipoles<CG_CLUSTER>(*cmm, 0, cmm->nCells, false, MAX_THREADS);
				cmm_check_cluster<CG_CLUSTER>(*cmm);
				// imethod == 0 : LJ interaction for all clusters, ignoring the neighbors
				// imethod == 1 : LJ interaction for the clusters in different neighbors
				// here we use imethod = 0
				CGMMCheckLocalRelshipInCell(*cmm, cgmm, 0, cgmm.nCluster, 0, MAX_THREADS); // setup relationship for each atom with CMM
			}
			//calColumbMultipoles<CG_CLUSTER>(*cmm, 0, cmm->nCells, false, MAX_THREADS);
			iCheckCMM++; if (iCheckCMM >= cmm->nCheckCluster) iCheckCMM = 0;

			if (bImplicitSolvent) {
				// neighbors of cluster for implicit solvation 
				MacroMolCellOperate<CG_MMOLECULE, CG_CLUSTER>((void*)(&CGMMClusterImpSolvNeighborsInCell), *cmm, cgmm, 0, cgmm.nCluster, MAX_THREADS);
			}
			
			// parallel calculation of interaction with CMM, assigning job on clusters of macro-molecule
#if _DISABLE_ == 0
		#error needs to be realized for interaction calculation
			MMClusterInteraction<CG_MMOLECULE, CG_CLUSTER>((void*)(&FreeCellCMMInteraction_SingleCGMM), *cmm, cgmm, 0, cgmm.nCluster, Ep, Overlap, false, MAX_THREADS);
#endif
			break;
		case 0:
			//FreeCellClusterInteraction(mm, 0, mm.nClusters, Ep, Overlap, false); // calculate the interactions between the clusters
			sprintf(errmsg, "interaction in coarse-grained macromolecule has to be calculated with CMM");
			show_msg(errmsg, false);
			break;
		default: // do not calculate interaction
			break; 
		}
		// bound-interaction
		MacroMoleculeOperate2<CG_MMOLECULE, double>((void*)(&cgbound_interact), cgmm, 0, cgmm.nCluster, MAX_THREADS, Ek);
		Ep += Ek;

		if (::bCGdipole && ::Eext != 0) {
			cgmm.bound_dipole();
			MacroMoleculeOperate2<CG_MMOLECULE, double>((void*)(&calDipoleEffect), cgmm, 0, cgmm.nBound, MAX_THREADS, Ek);
			Ep += Ek;
		}

		if (bImplicitSolvent) {
			MacroMolCellOperate<CG_MMOLECULE, CG_CLUSTER>((void*)(&CGMM_ImplicitSolvationForce), *cmm, cgmm, 0, cgmm.nCluster - 1, MAX_THREADS); // solvation effect
		}

		if (scale_force) MacroMoleculeOperate<CG_MMOLECULE>((void*)(&scale_force_cgcluster), cgmm, 0, cgmm.nCluster, MAX_THREADS);

		if (::bRandomForce) {
			_rand_::init_random_array(randa, cgmm.nCluster);
			gen_add_random_force_cgmm(cgmm, randa);
		}

#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "interaction calculation takes %d", time.elapse());
		show_infor(errmsg);
#endif

		mdpar.mmsave.save_kinetic_pars(cgmm, nloop, true, NULL, true);
		mdsave->mESave.one_more();
		mdsave->save_energy(cgmm.Ek, nloop, false, NULL, true);
		mdsave->save_energy((float)Ep, nloop, false, NULL, false);// save Ep
		if (mdsave->mESave.bsave()) mdsave->mESave.save_over();
		mdpar.save_mol_dyn_pars(cgmm, nloop, true, NULL, true);

		show_molecule_struct(cgmm, refresh_mol_xyz_struct_fname);

		//sprintf(errmsg, "[%f, %f, %f]", float(mm.mDynKin.cm.v[0]), float(mm.mDynKin.cm.v[1]), float(mm.mDynKin.cm.v[2]));
		//sprintf(errmsg, "[%f, %f, %f]", float(mm.base->atom[0].r.v[0]), float(mm.base->atom[0].r.v[1]), float(mm.base->atom[0].r.v[2]));
		//show_infor(errmsg);

#if _CHECK_TIME_STAMP_ == 1
		time.start();
#endif

		strcpy(msg_show, "\0"); strcpy(msg_log, "\0");
		// to get self-consistent velcoty and acceleration at this step based on leap-frog verlert algorithm
		// since we are simulating single macro-molecule, set whole molecule translating thermal energy 0, rotation thermal energy 1.5 kT
		CGLeapFrog_Verlet(cgmm, mdpar, vbf, 3 * ::Ek0Internal, init_md, msg_show, msg_log, n_msglen);
		if (strlen(msg_show) > 1) show_infor(msg_show, true);
		if (strlen(msg_log) > 1) mlog.show(msg_log, true);
		Ek = cgmm.Ek;
		
#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "NEIMO + Verlet take %d", time.elapse());
		show_infor(errmsg);
#endif

		//show_molecule_struct(mm, refresh_mol_xyz_struct_fname);

		show_loop_infor(nloop, mdpar.mmsave.mdsave->mESave.save_indx(), mdpar.LOOPS - nloop, (float)Ek, (float)Ep);

		nloop++;
	}
}

#if _SYS_ == _DOS_SYS_
void show_molecule_struct(CG_MMOLECULE &mm, char *fname) {
	ishow++;
	if (nshow > 0 && ishow < nshow) return;
	ishow = 0;
	if (fname == NULL || strlen(fname) == 0) return;
	save_xyz_struct(mm, fname);
}

#elif _SYS_ == _WINDOWS_SYS_

void show_molecule_struct(CG_MMOLECULE &mm, char *fname) {
	ishow++;
	if (nshow > 0 && ishow < nshow) return;
	ishow = 0;

	//if (fname != NULL && strlen(fname) > 0) save_xyz_struct(mm, fname);

	if (pMainFrame == NULL) return;

	if (delay_time > 0) WaitForSingleObject(hWaitEvent, delay_time);

	ResetEvent(hWaitGetMolStructEvent);
	PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
	WaitForSingleObject(hWaitGetMolStructEvent, max_time_pickup_molstruct);
	//WaitForSingleObject(hWaitGetMolStructEvent, INFINITE);
}

#elif _SYS_ == _LINUX_SYS_

void show_molecule_struct(CG_MMOLECULE &mm, char *fname) {
	ishow++;
	if (nshow > 0 && ishow < nshow) return;
	ishow = 0;

	if (fname != NULL && strlen(fname) > 0) save_xyz_struct(mm, fname);
}

#endif

} // end of namespace _coarse_grain_ 
