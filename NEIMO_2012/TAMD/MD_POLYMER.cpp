
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
#include "Interaction1.h"
#include "Interact2.h"
#include "NEIMO.h"
#include "MD.h"
#include "MD_POLYMER.h"

#include "var.h"
#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

bool trace = false;
/*
void get_Frame_Hoover_force(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, double Ethermal_trans, double Ethermal_rot, double &ksi_frame_trans, double &vksi_frame_trans, double &ksi_frame_rot, double &vksi_frame_rot, VECTOR<6> &f) {
	double force, torque, Ek;
	int i = 0, j = 0;
	double unit_E = kT / U_InertiaMoment;
	double dE_max = 0;
	double Ethermal = 0;
	double Emode = 0;
	double T = 0;
	double t;

	// rotation
	Ethermal = 2 * Ethermal_rot * unit_E; // here we calculated energy without dividing 2

	MATRIX<3> I;
	VECTOR3 u, w, dir;
	VECTOR3 Iw;
	double Ir = 0, dIw_r = 0, dEk = 0;
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) I.m[i][j] = mm.mDynKin.I.m[i][j];
	V2wv(w, u, mm.mDynKin.V)
	MpV3(I, w, Iw) // Iw = I * w;
	Ek = Iw.v[0] * w.v[0] + Iw.v[1] * w.v[1] + Iw.v[2] * w.v[2];
	if (Ethermal > 0.01 * unit_E) T = Ek / Ethermal;
	else T = 1.6;
	vksi_frame_rot = mdpar.inv_tao_ss * (T - 1);
	if (T > 1.2) vksi_frame_trans *= 1.1;
	//ksi_frame_rot = mdpar.ksi_frame_rot + (vksi_frame_rot + mdpar.vksi_frame_rot) * mdpar.dt * 0.5;
	ksi_frame_rot = mdpar.ksi_frame_rot + vksi_frame_rot * mdpar.dt;
	if (FABS(ksi_frame_rot) > mdpar.max_ksi_Hoover) ksi_frame_rot = sign(ksi_frame_rot) * mdpar.max_ksi_Hoover;
	
	for (i = 0; i < 3; i++) f.v[i] = -ksi_frame_rot * Iw.v[i];

	// translation
	// rotation
	Ethermal = 2 * Ethermal_trans * unit_E;  // here we calculated energy without dividing 2

	double dIw_t = 0;
	ABS2(u, i, t)
	Ek = mm.mDynKin.M * t;
	if (Ethermal > 0.01 * unit_E) T = Ek / Ethermal;
	else T = 1.6;
	vksi_frame_trans = mdpar.inv_tao_ss * (T - 1);
	if (T > 1.2) vksi_frame_trans *= 1.1;
	//ksi_frame_trans = mdpar.ksi_frame_trans + (vksi_frame_trans + mdpar.vksi_frame_trans) * mdpar.dt * 0.5;
	ksi_frame_trans = mdpar.ksi_frame_trans + vksi_frame_trans * mdpar.dt;
	if (FABS(ksi_frame_trans) > mdpar.max_ksi_Hoover) ksi_frame_trans = sign(ksi_frame_trans) * mdpar.max_ksi_Hoover;
	
	for (i = 0; i < 3; i++) f.v[i+3] = -ksi_frame_trans * mm.mDynKin.M * u.v[i];

	return;
}

// set cluster's force & torque due to the Hoover correction on acceleration of the mass center
void setClusterHooverForce(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, double Ethermal_trans, double Ethermal_rot, double &ksi_frame_trans, double &vksi_frame_trans, double &ksi_frame_rot, double &vksi_frame_rot, double ksi_internal) {
	int i;
	VECTOR<6> alpha, f;
	VECTOR3 dr;
	MATRIX<6> phai, phai_T;
	CLUSTER *pc = NULL;
	double t = 0, force = 0, torque = 0;

	V6zero(mm.mDynKin.f_FrameHoover)
	
	V6zero(f)
	get_Frame_Hoover_force(mm, mdpar, Ethermal_trans, Ethermal_rot, ksi_frame_trans, vksi_frame_trans, ksi_frame_rot, vksi_frame_rot, f);
	
	V62V6(f, mm.mDynKin.f_FrameHoover)
}
*/

void max_diff(VECTOR<6> *v1, VECTOR<6> *v2, int n, double &wmax, double &umax) {
	int nc = 0, i = 0;
	double udiff = 0, wdiff = 0;
	VECTOR<6> d;
	VECTOR3 u, w;
	for (nc = 0; nc < n; nc++) {
		V6minusV6(v1[nc], v2[nc], d)
		V2wv(w, u, d)
		MAX(u, i, udiff) MAX(w, i, wdiff)
		udiff = FABS(udiff); wdiff = FABS(wdiff);
		if (nc == 0) {wmax = wdiff; umax = udiff;}
		else {
			wmax = (wmax > wdiff ? wmax : wdiff);
			umax = (umax > udiff ? umax : udiff);
		}
	}
	return;
}

// initilize the local cordinates and calculate the inertia moment of each cluster and the whole molecule

void CalClusterInertiaMoment(MMOLECULE *mm, int n1, int n2) {
	int nc, i, j;
	VECTOR3 dr;
	CLUSTER *pc = NULL;
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		pc->InertialCord(); // calcualte the inertial cordinate of each atoms in each cluster
		pc->MassCenter();
		pc->calJacobianMatrix();

		for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm->mDynKin.cm.v[i];
		PHAI(pc->invNE->phai_cm2me, dr, i, j)
		M6T2M6(pc->invNE->phai_cm2me, pc->invNE->phai_T_cm2me, i, j)
		for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - CLUSTER_ORIGIN(mm->base, i);
		PHAI(pc->invNE->phai_base2me, dr, i, j)
		M6T2M6(pc->invNE->phai_base2me, pc->invNE->phai_T_base2me, i, j)

		pc->InertiaTensor();
	}
}

void CalMolInertiaMoment(MMOLECULE &mm, int nThreads) {
	mm.calMassCenter(true); // here to calculate molecule mass center and set mm.r0 -- from mass center to base
	/*
	int nc;
	for (nc = 0; nc < mm.nCluster; nc++) {
		mm.cluster[nc].InertialCord(); // calcualte the inertial cordinate of each atoms in each cluster
		mm.cluster[nc].MassCenter();
		mm.cluster[nc].calJacobianMatrix();
		mm.cluster[nc].InertiaTensor();
	}
	*/

	if (mm.nCluster > NCLUSTERS_PARALLEL && nThreads > 1) {
		MacroMoleculeOperate<MMOLECULE>((void*)(&CalClusterInertiaMoment), mm, 0, mm.nCluster, nThreads);
	}
	else CalClusterInertiaMoment(&mm, 0, mm.nCluster - 1);

	//if (::MM_SpeedVerlet_method == 0) mm.CalMolInertiaMoment(); // calculate whole molecule's inertia moment relative to mass center
	mm.CalMolInertiaMoment(); // calculate whole molecule's inertia moment relative to mass center
}

void CalCluster_ab(MMOLECULE *mm, int n1, int n2) {
	int nc;
	for (nc = n1; nc <= n2; nc++) {
		if (nc >= mm->nCluster) break;
		mm->cluster[nc].ab();
	}
}

// calculate Coriolis acceleration term and gyroscopic spatial force
void CalMol_ab(MMOLECULE &mm) {
	int nThreads = 1;
	if (mm.nCluster > NCLUSTERS_PARALLEL) nThreads = MAX_THREADS;
	MacroMoleculeOperate<MMOLECULE>((void*)(&CalCluster_ab), mm, 0, mm.nCluster, nThreads);

	if (::MM_SpeedVerlet_method == 0) mm.CalGyroscopicSpatialForce(); // at mass center
}

void SingleMM_NHC(LFVERLET_MM_MD& mdpar, double Ek) {
	double Ti = 0;
	Ti = Ek / (mdpar.iHDyn->nhc->N * ::Ek0Internal);
	mdpar.iHDyn->nhc->set_dT(Ti - 1);
	mdpar.iHDyn->nhc->verlet_NHC();

	mdpar.fetch_ksi();

	double ksi_c = 0, Tc = 1.8;
	if (Ti > Tc) {
		ksi_c = Ti / Tc; ksi_c *= ksi_c * 0.6;
		ksi_c = ::max_ksi_Hoover * (1 - exp(-ksi_c));
		if (mdpar.ksi < ksi_c) mdpar.ksi = ksi_c;
	}
}

bool LeapFrog_SelfConsistent_velocity_acceleration(bool bHooverDyn, MMOLECULE &mm, LFVERLET_MM_MD& mdpar, int md_loop, int &nconvg_loops, VECTOR<6> *vect1, VECTOR<6> *vect2, double *dv1, double *dv2, bool local_relax, char *show_msg, char *log_msg, int n_msglen, bool init_md) {
	char msg[512] = "\0";

	extern void ScaleVelocity(MMOLECULE &mm, double c);

	int nc = 0, i = 0, j;
	CLUSTER *pc = NULL;
	double T = 0, Ek0 = 0, dEk = 0;
	double maxEk = 9;
	if (_bias_force_::biasOn) maxEk = 64;
	double Ethermal = (mm.nCluster + 5) * ::Ek0Internal, Ekmax = Ethermal * maxEk; // 4 times of thermal kinetic energy. 
	// When Macromol has kinetic energy higher than Ekmax, we think it is too fast and the kinetic energy of the macromol will be scaled

	double wmax_diff = 0, umax_diff = 0;
	double uConvg = 0.005 / mdpar.dt, wConvg = 0.005 / mdpar.dt; //0.017 is 1 arc degree;
	if (uConvg < 2e-4) uConvg = 2e-4;
	if (wConvg < 2e-4) wConvg = 2e-4;
	
	bool convergent = false;
	int nmax_convg = 8;
	int MD_LOOP = md_loop;

	VECTOR<6> *tp_last = vect1, *tp_new = vect2;
	double *mtp = dv1, *mtpp = dv2;
	VECTOR<6> mv0, malpha0;

	bool tp_scaled = false;

	// 0 -- use the spatial moment at mass center, 1 -- use the velocity of base cluster
	if (::MM_SpeedVerlet_method == 0) MakeFreeBaseVelocityConsistent(mm); // adjusting free base velocity
	else if (::MM_SpeedVerlet_method == 1) {
		calClusterVelocity0(mm); calSpatialMomentum0(mm, true); // required for Nose-Hover force calculation of the frame
	}

	double ksi = mdpar.ksi;

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
			ksi = 0;
		}
		else {
			Ek0 = calKineticEnergy(mm); // kinetic energy of while macro-molecule

			if (bHooverDyn) {
				SingleMM_NHC(mdpar, mm.Ek);
				ksi = mdpar.ksi;
			}
			else {ksi = 0;}

			if (mm.Ek > Ekmax) { // too fast, the kinetic energy will be scaled
				sprintf(msg, "rescaled kinetic energy at loop %d", md_loop); show_log(msg, true);
				//Ek0 = ScaleKineticEnergy(mm, Ethermal);
				ScaleVelocity(mm, sqrt(Ethermal / mm.Ek));
				calClusterVelocity0(mm); calSpatialMomentum0(mm, true); 
				//calFrameKineticEnergy(mm);
				//SingleMM_NHC(mdpar, mm.Ek);
				mdpar.iHDyn->nhc->reset_NHC();
				mdpar.fetch_ksi();
				ksi = mdpar.ksi;
				md_loop = 0; // to trigger the case that do not get consistent velocity and accleration
				continue;
			}
		}

		// Brownian force
		if (bBrownian) {
			for (nc = 0; nc < mm.nCluster; nc++) {V6zero(mm.cluster[nc].dyn->f_Brown)}
			V6zero(mm.mDynKin.f_Brown)
			MacroMoleculeOperate2<MMOLECULE, VECTOR<6> >((void*)(&Brownian_force), mm, 0, mm.nCluster, MAX_THREADS, mm.mDynKin.f_Brown);
		}

		V6zero(mm.mDynKin.f_Hoover)

		NEIMO_CalHooverForce(mm, ksi, mm.mDynKin.f_Hoover);

		//mm.calTotalForce();
		//NEIMO(mm, new_config);
		NEIMO(mm, true);  // it is always true even if velocity changed,, total force is to be calculated in NEIMO(...)

		if (new_config) new_config = false;

		if (local_relax) {
			Speed_Kinetic(mdpar, mm, true);
			convergent = true;
			if (mdpar.iHDyn != NULL) mdpar.iHDyn->nhc->accept_nextstep();
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

		if (md_loop == 0 || init_md) Speed_Kinetic(mdpar, mm, true);
		else {
			// get new speed at current position and +1/2 time step
			Speed_Verlet(mdpar, mm); // reset the speed in mm to the new speed
		}

		for (nc = 0; nc < mm.nCluster; nc++) {
			for (i = 0; i < 6; i++) tp_new[nc].v[i] = mm.cluster[nc].bkm->V0.v[i];
		}
		max_diff(tp_new, tp_last, mm.nCluster, wmax_diff, umax_diff); 
		calKineticEnergy(mm); // kinetic energy 

		if (md_loop == 0 || init_md) {convergent = true; break;}
		else if (umax_diff < uConvg && wmax_diff < wConvg) convergent = true;

		nconvg++;
		if (nconvg > nmax_convg) {convergent = true; break;}
	}

	if (!convergent && mm.Ek > Ekmax) {
		sprintf(msg, " failure to get consistent speed & acceleration at loop %d ", MD_LOOP);
		if (n_msglen - strlen(show_msg) > 100) {if (strlen(show_msg) > 10) strcat(show_msg, "\n"); strcat(show_msg, msg);}
		// use the acceleration speed based on initial speed
		// torsion angle is not changed in the iteration

		ScaleVelocity(mm, sqrt(Ethermal / mm.Ek));
		//CalMol_ab(mm);
		for (nc = 0; nc < mm.nCluster; nc++) {
			pc = mm.cluster + nc; 
			V6zero(pc->invNE->a) V6zero(pc->invNE->b)
		}
		Speed_Kinetic(mdpar, mm, true);

		if (bHooverDyn) {
			//if (mdpar.ftHDyn != NULL) mdpar.ftHDyn->nhc->reset_NHC();
			//if (mdpar.frHDyn != NULL) mdpar.frHDyn->nhc->reset_NHC();
			//if (mdpar.iHDyn != NULL) mdpar.iHDyn->nhc->reset_NHC();
			//mdpar.fetch_ksi();
		}
	}

	if (bHooverDyn) {
		if (mdpar.iHDyn != NULL) mdpar.iHDyn->nhc->accept_nextstep();
	}

	return convergent;
}

void InitVelocity(MMOLECULE &mm) {
	int i = 0, j = 0;
	//float unit_E = kT / U_InertiaMoment;

	// calculate the inertia moment of each cluster and the whole molecule
	CalMolInertiaMoment(mm, MAX_THREADS);
	double Ep = 0;
	bool Overlap = false;

	V6zero(mm.V0); // set V0 of base 0
	//Ek = InitMolClusterVelocity(mm); // constructing the rest clusters
	for (int nc = 0; nc < mm.nCluster; nc++) {
		if (mm.cluster[nc].parent == NULL) continue;
		mm.cluster[nc].km.tp = (ranf() > 0.5 ? 1 : -1) * ranf() * 2e-3;
	}
	double Ek0 = (mm.nCluster - 1) * 0.5; // thermal energy
	ScaleKineticEnergy(mm, Ek0, true);
}

void ZeroVelocity(MMOLECULE &mm) {
	int nc;
	CLUSTER *pc = NULL;
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		pc->km.tp = 0; pc->km.tpp = 0; V6zero(pc->bkm->V0)
		V6zero(pc->invNE->a) V6zero(pc->invNE->b)
	}
	V6zero(mm.mDynKin.P) V6zero(mm.mDynKin.V) V6zero(mm.mDynKin.alpha)
	V6zero(mm.V0) V6zero(mm.alpha0) V6zero(mm.Iw)
}

//adjusting free base velocity so that the total spacial moment is consistent
void MakeFreeBaseVelocityConsistent(MMOLECULE &mm) {
	VECTOR<6> Iw0;
	CLUSTER *pc = NULL;

	calClusterVelocity0(mm);
	calSpatialMomentum(mm);
	mm.CalFreeBaseVelocity_GivenMassCenterMomentum();
	calClusterVelocity0(mm);
	calSpatialMomentum(mm);
	calKineticEnergy(mm);
	calFrameKineticEnergy(mm);
	mm.Ek_internal = mm.Ek - mm.Ek_frame;
}

void MM_ClusterAllocate2CMM(MM_CMM_JOB<MACROMOL_THREAD_VARS<MMOLECULE>, CMM_CELL3D<BASIC_CLUSTER> > *job) {
	MM_ClusterAllocateCMM_array<BASIC_CLUSTER, MMOLECULE>(job->cmm, job->job->mm, job->job->n1, job->job->n2);
}

extern void EInteract(BASIC_CLUSTER *pc, EInteractVar &evar, InteractRes &res);
extern void CLUSTER_LJ_INTERACT(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, InteractVar &var, InteractRes &res);
void MM_FreeCellInteract(MM_INTERACT_JOB *job) {
	int i, imol, ic;
	MMOLECULE *mm = job->job->mm;
	CMM_CELL3D<BASIC_CLUSTER> *cmm = job->cmm;
	BASIC_CLUSTER * bpc = NULL;
	InteractRes tres;
	job->LJ_res.reset(); job->eres.reset();
	for (ic = job->job->n1; ic <= job->job->n2; ic++) {
		if (ic >= mm->nCluster) break;
		bpc = (BASIC_CLUSTER*)(mm->cluster + ic);
		CLUSTER_LJ_INTERACT(*cmm, bpc, job->var, tres);
		job->LJ_res += tres;
		if (!bpc->eNeutral) {
			EInteract(bpc, job->evar, tres);
			job->eres += tres;
		}
	}
}

extern void CLUSTER_TOTAL_FORCE(BASIC_CLUSTER *pc);
void MM_ClusterTotalForce(MACROMOL_THREAD_VARS<MMOLECULE> *job) {
	MMOLECULE *mm = job->mm;
	for (int ic = job->n1; ic <= job->n2; ic++) {
		if (ic >= mm->nCluster) break;
		CLUSTER_TOTAL_FORCE((BASIC_CLUSTER*)(mm->cluster + ic));
	}
}

void MM_Cluster_atomic_rg(MMOLECULE *mm, int nc1, int nc2) {
	int ic, iatom;
	CLUSTER *pc;
	MATOM *patom;
	double *r1, *r2, *r3;
	for (ic = nc1; ic <= nc2; ic++) {
		if (ic >= mm->nCluster) break;
		pc = mm->cluster + ic;
		if (pc->bCMMShift) {
			r2 = pc->dr_cmm.v;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom;
				r1 = patom->r.v; r3 = patom->rg.v;
				//V3plusV3(patom->r, pc->dr_cmm, patom->rg)
				VplusV3(r1, r2, r3)
			}
		}
		else {
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom;
				V32V3(patom->r, patom->rg)
			}
		}
	}
}

// STRANGE : FUNCTION defined with variable (&, &) is same as (*, *)
//void MM_TorsionInteract(MACROMOL_THREAD_VARS<MMOLECULE> *mjob, TorsionVar &tvar, InteractRes &ires) {
void MM_TorsionInteract(MACROMOL_THREAD_VARS<MMOLECULE> *mjob, TorsionVar *tvar, InteractRes *ires) {
	ires->reset();
	int ic;
	CLUSTER *pc;
	MMOLECULE *mm = mjob->mm;
	TorsionInteract(mm, mjob->n1, mjob->n2, *tvar, *ires);
}


void MM_MassCenter(MACROMOL_THREAD_VARS<MMOLECULE> *mjob, VECTOR<6> *mp) {
	MMOLECULE *mm = mjob->mm;
	CLUSTER *pc = NULL;
	int ic;
	VECTOR3 Rcm;

	for (ic = mjob->n1; ic <= mjob->n2; ic++) {
		if (ic >= mm->nCluster) break;
		pc = mm->cluster + ic;
		Rcm.v[0] += pc->M * pc->cm.v[0];
		Rcm.v[1] += pc->M * pc->cm.v[1];
		Rcm.v[2] += pc->M * pc->cm.v[2];
	}

	memcpy(mp->v, Rcm.v, SIZE_V3);
}

void MM_ApplySpring(MACROMOL_THREAD_VARS<MMOLECULE> *mjob, VECTOR<6> *mp) {
	MMOLECULE *mm = mjob->mm;
	CLUSTER *pc = NULL;
	int ic;

	VECTOR3 f, t, dr;

	for (ic = mjob->n1; ic <= mjob->n2; ic++) {
		if (ic >= mm->nCluster) break;
		pc = mm->cluster + ic;
		f.v[0] = pc->M * mp->v[0];
		f.v[1] = pc->M * mp->v[1];
		f.v[2] = pc->M * mp->v[2];
		VECT3((*(pc->Op)), pc->cm, dr)
		V3PV3(dr, f, t)
		pc->dyn->fc.v[0] += t.v[0];
		pc->dyn->fc.v[1] += t.v[1];
		pc->dyn->fc.v[2] += t.v[2];
		pc->dyn->fc.v[3] += f.v[0];
		pc->dyn->fc.v[4] += f.v[1];
		pc->dyn->fc.v[5] += f.v[2];
	}
}

void ApplySpring(MACROMOL_THREAD_VARS<MMOLECULE> *mjob, VECTOR<6> *mp, int nJob, float *k) {
	int i, j;
	/*
	for (i = 0; i < nJob; i++) memset(mp[i].v, 0, SIZE_V6);
	MTOperate1<MACROMOL_THREAD_VARS<MMOLECULE>, VECTOR<6> >((void*)(&MM_MassCenter), mjob, nJob, mp);
	float M = mjob->mm->mDynKin.M;
	double Rz = 0;	
	for (j = 0; j < 3; j++) {
		Rz = 0;
		for (i = 0; i < nJob; i++) Rz += mp[i].v[j];
		Rz /= M;
		Rz *= -k;
		for (i = 0; i < nJob; i++) mp[i].v[j] = Rz; 
	}
	*/
	MMOLECULE *mm = mjob->mm;
	double d = 0;
	for (i = 0; i < nJob; i++) {
		for (j = 0; j < 3; j++) {
			d = mm->mDynKin.cm.v[j];
			mp[i].v[j] = (k[1] * (d > 0 ? -1 : -1) * d * d - k[0] * d) / mm->mDynKin.M;
		}
	}
	
	MTOperate1<MACROMOL_THREAD_VARS<MMOLECULE>, VECTOR<6> >((void*)(&MM_ApplySpring), mjob, nJob, mp);
}



void MD_POLYMER(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, CMM_CELL3D<BASIC_CLUSTER> *cmm, int interact_cal_method, int nNormalMD) {
	::interact_cal_method = interact_cal_method; 
	char msg_show[5120] = "\0", msg_log[5120] = "\0";
	int n_msglen = 5120;

	// apply spring at the mass center of mm
	::bZSpring = true;
	::k_zspring[0] = 10.0f * ::kT / ::U_MassAccel;
	::k_zspring[1] = 10.0f * ::kT / ::U_MassAccel;

	long nloop = 0;
	int nc = 0, ncell = 0;
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

	int i;

	MACROMOL_THREAD_VARS<MMOLECULE> mmJob[MAX_THREADS];
	assignMMJob<MMOLECULE, MACROMOL_THREAD_VARS<MMOLECULE> >(&mm, mmJob, MAX_THREADS);

	CMM_CELL3D<BASIC_CLUSTER> job_cmm[MAX_THREADS];
	MM_CMM_JOB<MACROMOL_THREAD_VARS<MMOLECULE>, CMM_CELL3D<BASIC_CLUSTER> > mmCMMAllocator[MAX_THREADS];
	MM_INTERACT_JOB job_interact[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		if (MAX_THREADS > 1) {
			CMM_array_construct_for_storage<BASIC_CLUSTER>(job_cmm + i, cmm, nCluster_CMM);
			mmCMMAllocator[i].cmm = job_cmm + i;
		}
		else {
			mmCMMAllocator[i].cmm = cmm;
		}
		mmCMMAllocator[i].set_pars(i, mmJob + i);

		mmJob[i]._THREAD_VARS::set_pars(i);

		job_interact[i].set_pars(i, mmJob + i);
		job_interact[i].cmm = cmm;
		job_interact[i].var.bVirial = false;
		job_interact[i].evar.esv1.bEwald = false; // not EwaldSum !
		job_interact[i].evar.bVirial = false; 
		job_interact[i].evar.bST_diagonal = false; 
		job_interact[i].evar.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
		job_interact[i].evar.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	}

	mm.setup_cluster_fatom();
	mm.setup_hinge_dihedral();
	TorsionVar torsion_var[MAX_THREADS];
	HingeConstraintVar hinge_constraint_var[MAX_THREADS];
	InteractRes mires[MAX_THREADS];
	for (i = 0; i < MAX_THREADS; i++) {
		torsion_var[i].bVirial = ::bVirial; torsion_var[i].bST_diagonal = true;
		hinge_constraint_var[i].bVirial = ::bVirial; hinge_constraint_var[i].bST_diagonal = true;
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

	{
		// set mass center to {0, 0, 0}
		VECTOR<3> cm;
		mm.SetMassCenter(cm.v, true, true);
	}
	if (local_relax) {
		CalMolInertiaMoment(mm, MAX_THREADS); 
		ZeroVelocity(mm);
	}
	else InitVelocity(mm);

	CHAIN<short> *ch_cluster_id = NULL;
	if (::bImplicitSolvent) mm.calClusterSolvationRadius(); // it uses the cluster's mass center position
	if (::bBrownian) {
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
	}

	mdpar.Init(mm); // copy initial torsion angle and speed to mdpar
	
	mm.resetDynPars(); // set force & torque to 0

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
#endif

	int iNormalMDLoop = 0;
	nNormalMD *= mdsave->nsave;
	bool init_md = false;
	int nCMMCheck = cmm->nCheckCluster;

	VECTOR3 mIw, mP;
	float rand_torque = 0;

	NEIMO_THREAD_VARS neimo_matrix_cal_var;

	bool bCheckClusterInCMM = true;

	_rand_::RAND_ARRAY<int> randa;

	for (i = 0; i < mm.nCluster; i++) mm.cluster[i].eNeutral = true;

	while (nloop < mdpar.LOOPS) {
		if (COMMAND_MD_STOP) break; // stop with given command
		init_md = (nloop == 0 ? true : false); 
		if (iCheckCMM == 0) bCheckClusterInCMM = true;
		else bCheckClusterInCMM = (nloop / N_MDSTEP_CHECK_CMM * N_MDSTEP_CHECK_CMM == nloop ? true : false);

		if (bCheckClusterInCMM) {
		{
			// set mass center to {0, 0, 0}
			VECTOR<3> cm;
			mm.SetMassCenter(cm.v, true, true);
		}
		}
		MacroMoleculeOperate((void*)(&MM_Cluster_atomic_rg), mm, 0, mm.nCluster, MAX_THREADS);
		
		// calculate the inertia moment of each cluster and the whole molecule
		CalMolInertiaMoment(mm, MAX_THREADS);

		if (local_relax) ZeroVelocity(mm);

#if _CHECK_TIME_STAMP_ == 1
		time.start();
#endif
		// we have to do it always before we start new thread, so that the thread-related variable will be refreshed
		// calculate the neimo-matrix, related to the innertial matrix only
		//neimo_matrix_cal_var.set_pars(&mm, 0);
		//Start_NEIMO_MATRIX_CALC(mm, neimo_matrix_cal_var);

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
				if (MAX_THREADS == 1) {
					cmm->reset_cell_chain();
					MM_ClusterAllocate2CMM(mmCMMAllocator);
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
		for (i = 0; i < nThreads; i++) {
			Ep += job_interact[i].LJ_res.U + job_interact[i].eres.U;
		}

//#error this funciton is strange
		MTOperate2< MACROMOL_THREAD_VARS<MMOLECULE>, TorsionVar, InteractRes >((void*)(&MM_TorsionInteract), mmJob, nThreads, torsion_var, mires);
		Etorsion = 0; 
		for (i = 0; i < nThreads; i++) Etorsion += mires[i].U;
		Ep += Etorsion;

		if (::bEext) {
			MOperate3<CLUSTER, float, double>((void*)(&CLUSTER_ExternalEField), mm.cluster, mm.nCluster, ::Eext, MAX_THREADS, Uext);
			Ep += Uext;
		}

		// before scale cluster, momentum of cluster has to be updated, because when atom are shocked, momentum is to be used to scale the interaction !
		// 0 -- use the spatial moment at mass center, 1 -- use the velocity of base cluster
		if (::MM_SpeedVerlet_method == 0) MakeFreeBaseVelocityConsistent(mm); // adjusting free base velocity
		else if (::MM_SpeedVerlet_method == 1) {
			calClusterVelocity0(mm); calSpatialMomentum0(mm, true); // required for Nose-Hover force calculation of the frame
		}
		MTOperate< MACROMOL_THREAD_VARS<MMOLECULE> >((void*)(&MM_ClusterTotalForce), mmJob, nThreads);

#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "cluster interact.  calculation takes %d ms", time.glance());
		show_infor(errmsg);
#endif

		// implicit solvation energy / force
		if (bImplicitSolvent) {
			// neighbors of cluster for implicit solvation 
			if (bCheckClusterInCMM) MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MMClusterImpSolvNeighborsInCell), *cmm, mm, 0, mm.nCluster - 1, MAX_THREADS);
			//MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MM_ImplicitSolvationForce), *cmm, mm, 0, mm.nCluster, MAX_THREADS);
			MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MM_ImpSolvSASA), *cmm, mm, 0, mm.nCluster, MAX_THREADS);
			MacroMoleculeOperate<MMOLECULE>((void*)(&MM_ImpSolvForce), mm, 0, mm.nCluster, MAX_THREADS);
			for (nc = 0; nc < mm.nCluster; nc++) Ep += mm.cluster[nc].E_ImplicitSolvation;
		}

		if (bRandomForce && nloop / 20 * 20 == nloop) {
			_rand_::init_random_array(randa, mm.nCluster);
			gen_add_random_force(mm, mm.nCluster / 5, randf, fwhm_randf, randf, fwhm_randf, randa);
		}

		switch (interact_cal_method) {
		case _CMM:
			CellOperate<BASIC_CLUSTER>((void*)(&SCALE_FORCE_CMM), *cmm, 0, cmm->acell.n, MAX_THREADS);
			break;
		default:
			scale_force_mmolecule(mm);
			break;
		}

		bHooverDyn = true;

		CalNetForce(mm, 0, mm.nCluster); // the net force on molecular mass center

#if _CHECK_TIME_STAMP_ == 1
		sprintf(::errmsg, "interaction calculation takes %d ms", time.glance());
		show_infor(errmsg);
#endif

		mdpar.mmsave.save_kinetic_pars(mm, nloop, true, NULL, true);
		mdsave->mESave.one_more();
		//mdsave->save_energy(mm.Ek, nloop, false, NULL, true);
		mdpar.mmsave.save_energy(mm, nloop, false, NULL, true);
		mdsave->save_energy((float)Ep, nloop, false, NULL, false);// save Ep
		if (mdsave->mESave.bsave()) mdsave->mESave.save_over();
		mdpar.save_mol_dyn_pars(mm, nloop, true, NULL, true);

		show_molecule_struct(mm, refresh_mol_xyz_struct_fname);

		//sprintf(errmsg, "[%f, %f, %f]", float(mm.mDynKin.cm.v[0]), float(mm.mDynKin.cm.v[1]), float(mm.mDynKin.cm.v[2]));
		//sprintf(errmsg, "[%f, %f, %f]", float(mm.base->atom[0].r.v[0]), float(mm.base->atom[0].r.v[1]), float(mm.base->atom[0].r.v[2]));
		//show_infor(errmsg);

#if _CHECK_TIME_STAMP_ == 1
		//time.start();
#endif

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

		strcpy(msg_show, "\0"); strcpy(msg_log, "\0");
		// to get self-consistent velcoty and acceleration at this step based on leap-frog verlert algorithm
		// since we are simulating single macro-molecule, set whole molecule translating thermal energy 0, rotation thermal energy 1.5 kT
		LeapFrog_SelfConsistent_velocity_acceleration(bHooverDyn, mm, mdpar, nloop, nConvg, vect1, vect2, dv1, dv2, ::local_relax, msg_show, msg_log, n_msglen, init_md);
		if (strlen(msg_show) > 1) show_infor(msg_show, true);
		if (strlen(msg_log) > 1) mlog.show(msg_log, true);
		Ek = mm.Ek;

		//Ep = mm.Ek_frame;
		
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
		iNormalMDLoop++;
		if (iNormalMDLoop >= nNormalMD) iNormalMDLoop = 0;
	}
	if (vect1 != NULL) {delete[] vect1; vect1 = NULL;}
	if (vect2 != NULL) {delete[] vect2; vect2 = NULL;}
	if (dv1 != NULL) {delete[] dv1; dv1 = NULL;}
	if (dv2 != NULL) {delete[] dv2; dv2 = NULL;}
}

extern int delay_time; // in milisecond
extern int max_time_pickup_molstruct;  // in milisecond
extern int nshow;
extern int ishow;

#if _SYS_ == _DOS_SYS_
void show_loop_infor(int nloop, int nsaved, int nLoops, float Ek, float Ep) {
	char msg[256] = "\0";
	sprintf(msg, "loop # : %d; saved : %d; total : %d", nloop, nsaved, nLoops);
	cout<<msg<<endl;
	sprintf(msg, "Ek [kT] : %f; Ep [kT] : %f; total : %f", Ek, Ep, Ek + Ep);
	cout<<msg<<endl;
}

void show_infor(char *msg, bool endline) {
	mlog.show(msg, endline);
}

void show_molecule_struct(MMOLECULE &mm, char *fname) {
	ishow++;
	if (nshow > 0 && ishow < nshow) return;
	ishow = 0;
	if (fname == NULL || strlen(fname) == 0) return;
	save_xyz_struct(mm, fname);
}

#elif _SYS_ == _WINDOWS_SYS_

extern CWnd* pMainFrame;
extern HANDLE hWaitGetMolStructEvent;
extern HANDLE hWaitEvent;
extern CMDInforDlg *mdInforDlg;

void show_loop_infor(int nloop, int nsaved, int nLoops, float Ek, float Ep) {
	if (mdInforDlg == NULL) return;
	//if (ishow != 0) return;
	mdInforDlg->md_loop = nloop;
	mdInforDlg->md_nsaved = nsaved;
	mdInforDlg->md_torun = nLoops;
	mdInforDlg->md_Ek = Ek;
	mdInforDlg->md_Ep = Ep;
	mdInforDlg->md_Etotal = Ek + Ep;
	//::SendMessage(mdInforDlg->m_hWnd, WM_DIALOG_SHOW_LOOP, 0, 0); // show information
	::PostMessage(mdInforDlg->m_hWnd, WM_DIALOG_SHOW_LOOP, 0, 0); // show information
}

void show_infor(char *msg, bool endline) {
	mlog.show(msg, endline);
	if (mdInforDlg == NULL) return;
	mdInforDlg->append_msg(msg);
	::SendMessage(mdInforDlg->m_hWnd, WM_DIALOG_SHOW_INFOR, 0, 0); // show information
}

void show_molecule_struct(MMOLECULE &mm, char *fname) {
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

void show_all_molecules() {
	ishow++;
	if (nshow > 0 && ishow < nshow) return;
	ishow = 0;

	if (pMainFrame == NULL) return;
	if (delay_time > 0) WaitForSingleObject(hWaitEvent, delay_time);

	ResetEvent(hWaitGetMolStructEvent);
	PostMessage(pMainFrame->m_hWnd, WM_PICKUP_MOL_STRUCT, 0, 0);
	WaitForSingleObject(hWaitGetMolStructEvent, max_time_pickup_molstruct);
	//WaitForSingleObject(hWaitGetMolStructEvent, INFINITE);
}

#elif _SYS_ == _LINUX_SYS_
void show_loop_infor(int nloop, int nsaved, int nLoops, float Ek, float Ep) {
}

void show_infor(char *msg, bool endline) {
	mlog.show(msg, endline);
}

void show_molecule_struct(MMOLECULE &mm, char *fname) {
	ishow++;
	if (nshow > 0 && ishow < nshow) return;
	ishow = 0;

	if (fname != NULL && strlen(fname) > 0) save_xyz_struct(mm, fname);
}

void show_all_molecules() {
}
#endif

void cal_Hoover_force(MMOLECULE &mm, LFVERLET_MM_MD& mdpar) {
	double ksi = mdpar.ksi;
	NEIMO_CalHooverForce(mm, (float)ksi, mm.mDynKin.f_Hoover);
}

bool MM_SpeedVerlet(bool local_relax, bool speed_verlet, MMOLECULE &mm, LFVERLET_MM_MD& mdpar, VECTOR<6> *vect1, VECTOR<6> *vect2, double &vmax_diff, double &wmax_diff, char *show_msg, char *log_msg, int n_msglen) {
	char msg[512] = "\0";
	int nc = 0;
	CLUSTER *pc = NULL;
	double Ti = 0, Tft = 0, Tfr = 0, T = 0;
	double Eti = (mm.nCluster - 1) * ::Ek0Internal, Etfr = 3 * ::Ek0Internal, Etft = 3 * ::Ek0Internal;
	double Ethermal = Eti + Etfr + Etft, Ekmax = Ethermal * 4;

	VECTOR<6> *v_last = vect1, *v_new = vect2;

_MM_SpeedVerlet_Start:

	for (nc = 0; nc < mm.nCluster; nc++) {
		V62V6(mm.cluster[nc].bkm->V0, v_last[nc])
	}

	if (local_relax) {
		for (nc = 0; nc < mm.nCluster; nc++) {mm.cluster[nc].km.tp = 0; mm.cluster[nc].km.tpp = 0;}
		V6zero(mm.V0) V6zero(mm.Iw) V6zero(mm.alpha0)
	}
	else {
		cal_Hoover_force(mm, mdpar);
	}

	// speed is varying to get consisten speed and accleration based on leap-frog velert MD 
	// calculate Coriolis acceleration term and gyroscopic spatial force
	CalMol_ab(mm);
	
	if (bBrownian) {
		for (nc = 0; nc < mm.nCluster; nc++) {
			V6zero(mm.cluster[nc].dyn->f_Brown)
		}
		V6zero(mm.mDynKin.f_Brown)
		//MacroMoleculeOperate2<MMOLECULE, VECTOR<6> >((void*)(&Brownian_force), mm, 0, mm.nCluster, MAX_THREADS, mm.mDynKin.f_Brown);
		Brownian_force(&mm, 0, mm.nCluster - 1, mm.mDynKin.f_Brown);
	}

	//mm.calTotalForce();
	//NEIMO(mm, new_config);
	NEIMO(mm, true);  // it is always true even if velocity changed, total force is to be calculated in NEIMO(...)

	if (local_relax || !speed_verlet) {
		Speed_Kinetic(mdpar, mm, true);
	}
	else {
		// get new speed at current position and +1/2 time step
		Speed_Verlet(mdpar, mm); // reset the speed in mm to the new speed
	}

	for (nc = 0; nc < mm.nCluster; nc++) {
		V62V6(mm.cluster[nc].bkm->V0, v_new[nc])
	}
	max_diff(v_new, v_last, mm.nCluster, vmax_diff, wmax_diff);
	calKineticEnergy(mm); // kinetic energy 
	calFrameKineticEnergy(mm);
	mm.Ek_internal = mm.Ek - mm.Ek_frame;

	bool status = true;
	return status;
}

bool MM_SpeedVerlet_Simplified(bool speed_verlet, MMOLECULE &mm, LFVERLET_MM_MD& mdpar, VECTOR<6> *vect1, VECTOR<6> *vect2, double &vmax_diff, double &wmax_diff) {
	int nc = 0;
	CLUSTER *pc = NULL;

	VECTOR<6> *v_last = vect1, *v_new = vect2;

	// speed is varying to get consisten speed and accleration based on leap-frog velert MD 
	// calculate Coriolis acceleration term and gyroscopic spatial force
	CalMol_ab(mm);

	for (nc = 0; nc < mm.nCluster; nc++) {V62V6(mm.cluster[nc].bkm->V0, v_last[nc])}
	
	// the total force has to be calculated outside
	NEIMO(mm, true);  // it is always true even if velocity changed, total force is to be calculated in NEIMO(...)

	if (!speed_verlet) Speed_Kinetic(mdpar, mm, true);
	else Speed_Verlet(mdpar, mm); // get new speed at current position and +1/2 time step, and reset the speed in mm to the new speed
	
	for (nc = 0; nc < mm.nCluster; nc++) {V62V6(mm.cluster[nc].bkm->V0, v_new[nc])}
	max_diff(v_new, v_last, mm.nCluster, vmax_diff, wmax_diff);
	calKineticEnergy(mm); // kinetic energy 
	calFrameKineticEnergy(mm);
	mm.Ek_internal = mm.Ek - mm.Ek_frame;

	return true;
}

