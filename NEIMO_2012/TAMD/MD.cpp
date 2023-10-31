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

#include "ZMatrix.h"

#include "MM.h"
#include "Mol.h"
#include "MD.h"
#include "var.h"

using namespace std;

bool save_xyz_struct(MMOLECULE &mm, char *fname) {
	ofstream out;
	out.open(fname);
	if (!out.is_open()) {sprintf(errmsg, "can not open file %s to save macromolecular structure"); show_msg(errmsg); strcpy(errmsg, "\0"); return false;}
	mm.xyz_save(out);
	out.close();
	return true;
}

bool MD_SAVE::save_energy(float Ek, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	if (check_status) {
		if (!mESave.one_more()) return false;
	}
	if (!mESave.bsave()) return true; // this loop will be not saved

	if (bAsync) {
		mESave.init_buf();
		if (save_md_title) mESave.sbuf<<"[MD] "<<mESave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) mESave.sbuf<<additional_title<<endl;
		mESave.sbuf<<Ek<<endl;
		mESave.flush_buf();
	}
	else {
		if (save_md_title) *(mESave.out)<<"[MD] "<<mESave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) *(mESave.out)<<additional_title<<endl;
		*(mESave.out)<<Ek<<endl;
	}

	if (check_status) mESave.save_over();
	return true;
}

bool MM_MD_SAVE::save_energy(MMOLECULE &mm, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	if (check_status) {
		if (!mdsave->mESave.one_more()) return false;
	}
	if (!mdsave->mESave.bsave()) return true; // this loop will be not saved

	char msg[256] = "\0";

	if (bAsync) {
		mdsave->init_buf();
		if (save_md_title) mdsave->mESave.sbuf<<"[MD] "<<mdsave->mESave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) mdsave->mESave.sbuf<<additional_title<<endl;

		//*(mdsave->mESave.out)<<mm.Ek<<"  "<<mm.Ek_frame_trans<<"  "mm.Ek_frame_rot<<"  "<<mm.Ek_internal<<endl;
		sprintf(msg, "%12.4f    %8.4f    %8.4f    %12.4f", mm.Ek, mm.Ek_frame_trans, mm.Ek_frame_rot, mm.Ek_internal);
		mdsave->mESave.sbuf<<msg<<endl;
		mdsave->mESave.flush_buf();
	}
	else {
		if (save_md_title) *(mdsave->mESave.out)<<"[MD] "<<mdsave->mESave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) *(mdsave->mESave.out)<<additional_title<<endl;

		//*(mdsave->mESave.out)<<mm.Ek<<"  "<<mm.Ek_frame_trans<<"  "mm.Ek_frame_rot<<"  "<<mm.Ek_internal<<endl;
		sprintf(msg, "%12.4f    %8.4f    %8.4f    %12.4f", mm.Ek, mm.Ek_frame_trans, mm.Ek_frame_rot, mm.Ek_internal);
		*(mdsave->mESave.out)<<msg<<endl;
	}

	if (check_status) mdsave->mESave.save_over();
	return true;
}

bool MM_MD_SAVE::save_kinetic_pars(MMOLECULE &mm, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	if (mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!mdsave->mKinSave.one_more()) return false;
	}
	if (!mdsave->mKinSave.bsave()) return true; // this loop will be not saved

	int i = 0;
	double ta = 0;

	if (bAsync) {
		mdsave->init_buf();
		if (save_md_title) mdsave->mKinSave.sbuf<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0)
			mdsave->mKinSave.sbuf<<additional_title<<endl;
		mdsave->mKinSave.sbuf<<mm.base->Op->v[0]<<" "<<mm.base->Op->v[1]<<" "<<mm.base->Op->v[2]<<" ";
		mdsave->mKinSave.sbuf<<mm.zaxis.v[0]<<" "<<mm.zaxis.v[1]<<" "<<mm.zaxis.v[2]<<" ";
		mdsave->mKinSave.sbuf<<mm.xaxis.v[0]<<" "<<mm.xaxis.v[1]<<" "<<mm.xaxis.v[2]<<endl;

		mdsave->mKinSave.sbuf<<mm.V0.v[0]<<" "<<mm.V0.v[1]<<" "<<mm.V0.v[2]<<" ";
		mdsave->mKinSave.sbuf<<mm.V0.v[3]<<" "<<mm.V0.v[4]<<" "<<mm.V0.v[5]<<endl;

		mdsave->mKinSave.sbuf<<mm.alpha0.v[0]<<" "<<mm.alpha0.v[1]<<" "<<mm.alpha0.v[2]<<" ";
		mdsave->mKinSave.sbuf<<mm.alpha0.v[3]<<" "<<mm.alpha0.v[4]<<" "<<mm.alpha0.v[5]<<endl;

		for (i = 0; i < mm.nCluster; i++) {
			ta = mm.cluster[i].km.ta;
			//*(mdsave->mKinSave.out)<<i<<" "<<ta<<" "<<mm.cluster[i].km.tp<<" "<<mm.cluster[i].km.tpp<<endl;
			mdsave->mKinSave.sbuf<<i<<" "<<ta<<endl;
		}
		mdsave->mKinSave.sbuf<<endl;

		mdsave->mKinSave.flush_buf();
	}
	else {
		if (save_md_title) *(mdsave->mKinSave.out)<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0)
			*(mdsave->mKinSave.out)<<additional_title<<endl;
		*(mdsave->mKinSave.out)<<mm.base->Op->v[0]<<" "<<mm.base->Op->v[1]<<" "<<mm.base->Op->v[2]<<" ";
		*(mdsave->mKinSave.out)<<mm.zaxis.v[0]<<" "<<mm.zaxis.v[1]<<" "<<mm.zaxis.v[2]<<" ";
		*(mdsave->mKinSave.out)<<mm.xaxis.v[0]<<" "<<mm.xaxis.v[1]<<" "<<mm.xaxis.v[2]<<endl;

		*(mdsave->mKinSave.out)<<mm.V0.v[0]<<" "<<mm.V0.v[1]<<" "<<mm.V0.v[2]<<" ";
		*(mdsave->mKinSave.out)<<mm.V0.v[3]<<" "<<mm.V0.v[4]<<" "<<mm.V0.v[5]<<endl;

		*(mdsave->mKinSave.out)<<mm.alpha0.v[0]<<" "<<mm.alpha0.v[1]<<" "<<mm.alpha0.v[2]<<" ";
		*(mdsave->mKinSave.out)<<mm.alpha0.v[3]<<" "<<mm.alpha0.v[4]<<" "<<mm.alpha0.v[5]<<endl;

		for (i = 0; i < mm.nCluster; i++) {
			ta = mm.cluster[i].km.ta;
			//*(mdsave->mKinSave.out)<<i<<" "<<ta<<" "<<mm.cluster[i].km.tp<<" "<<mm.cluster[i].km.tpp<<endl;
			*(mdsave->mKinSave.out)<<i<<" "<<ta<<endl;
		}
		*(mdsave->mKinSave.out)<<endl;
	}

	if (check_status) mdsave->mKinSave.save_over();
	return true;
}

bool MM_MD_SAVE::save_atomic_dipole(MMOLECULE &mm, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	return true;
	/*
	if (!::bPolarize && !::bDipole) return true;

	if (mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!mdsave->muSave.one_more()) return false;
	}
	if (!mdsave->muSave.bsave()) return true; // this loop will be not saved

	int i = 0;
	double ta = 0;
	int ic, ia;
	CLUSTER *c;

	if (bAsync) {
		mdsave->init_buf();
		if (save_md_title) mdsave->muSave.sbuf<<"[MD] "<<mdsave->muSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0)
			mdsave->muSave.sbuf<<additional_title<<endl;
		for (ic = 0; ic < mm.nCluster; ic++) {
			c = mm.cluster + ic;
			for (ia = 0; ia < c->nAtoms; ia++) {
				mdsave->muSave.sbuf<<c->atom[ia].mu.v[0]<<" "<<c->atom[ia].mu.v[1]<<" "<<c->atom[ia].mu.v[2]<<" ";
				if (::bDipole) {
					mdsave->muSave.sbuf<<c->atom[ia].ind_mu.v[0]<<" "<<c->atom[ia].ind_mu.v[1]<<" "<<c->atom[ia].ind_mu.v[2]<<" ";
					mdsave->muSave.sbuf<<c->atom[ia].intr_mu.v[0]<<" "<<c->atom[ia].intr_mu.v[1]<<" "<<c->atom[ia].intr_mu.v[2];
				}
				mdsave->muSave.sbuf<<endl;
			}
		}
	
		mdsave->muSave.sbuf<<endl;

		mdsave->muSave.flush_buf();
	}
	else {
		if (save_md_title) *(mdsave->muSave.out)<<"[MD] "<<mdsave->muSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0)
			*(mdsave->muSave.out)<<additional_title<<endl;
		for (ic = 0; ic < mm.nCluster; ic++) {
			c = mm.cluster + ic;
			for (ia = 0; ia < c->nAtoms; ia++) {
				*(mdsave->muSave.out)<<c->atom[ia].mu.v[0]<<" "<<c->atom[ia].mu.v[1]<<" "<<c->atom[ia].mu.v[2]<<" ";
				if (::bDipole) {
					*(mdsave->muSave.out)<<c->atom[ia].ind_mu.v[0]<<" "<<c->atom[ia].ind_mu.v[1]<<" "<<c->atom[ia].ind_mu.v[2]<<" ";
					*(mdsave->muSave.out)<<c->atom[ia].intr_mu.v[0]<<" "<<c->atom[ia].intr_mu.v[1]<<" "<<c->atom[ia].intr_mu.v[2];
				}
				*(mdsave->muSave.out)<<endl;
			}
		}
	
		*(mdsave->muSave.out)<<endl;
	}

	if (check_status) mdsave->muSave.save_over();
	return true;
	*/
}

#ifndef _DISABLE_
double KineticEnergy_Foward(MMOLECULE &mm, LFVERLET_MM_MD &mdpar) {
	double Ek = 0;
	double tmp = 0;
	VECTOR<6> v0;
	int i, n;

#define SWAP_TP(mtp, tp_ph, tmp) tmp = mtp; mtp = tp_ph; tp_ph = tmp;

	V2V(mm.V0, v0, i, 6) V2V(mdpar.v0_ph, mm.V0, i, 6)
	for (n = 0; n < mm.nCluster; n++) {
		SWAP_TP(mm.cluster[n].km.tp, mdpar.tp_ph[n], tmp)
	}
	calClusterVelocity0(mm);
	Ek = calKineticEnergy(mm);
	V2V(v0, mm.V0, i, 6)
	for (n = 0; n < mm.nCluster; n++) {
		SWAP_TP(mm.cluster[n].km.tp, mdpar.tp_ph[n], tmp)
	}
	mm.base->km.tp = 0;
	calClusterVelocity0(mm);
	return Ek;

#undef SWAP_TP
}
#endif // _DISABLE_

void check_angle(MMOLECULE &mm, LFVERLET_MM_MD &mdpar) {
	int i = 0;
	double ta = 0;
	double dt = 0;
	for (i = 0; i < mm.nCluster; i++) {
		ta = mm.cluster[i].km.ta;
		dt = 0;
		while (ta > PI) {ta -= PI2; dt -= PI2;}
		while (ta < -PI) {ta += PI2; dt += PI2;}
		mm.cluster[i].km.ta = ta;
		if (mdpar.t != NULL) mdpar.t[i] += dt;
	}
}

void Speed_Verlet(LFVERLET_MM_MD &mdpar, MMOLECULE &mm) {
	int method = ::MM_SpeedVerlet_method; // 0 -- use the spatial moment at mass center, 1 -- use the spatial momentum of base cluster
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0, j;
	CLUSTER *p = NULL;

	if (method == 1 && mm.free_base) { // use the momentum of the base cluster, I am not sure this method is implemented properly !!, somthing is wrong !
		for (i = 0; i < 6; i++) {
			mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + mm.alpha0.v[i] * mdpar.dt; // ?
			mdpar.v0.v[i] = (mdpar.v0_mh.v[i] + mdpar.v0_ph.v[i]) * 0.5;
			mm.V0.v[i] = mdpar.v0.v[i];
		}
	}

	// check velocity of clusters except base 
	for (i = 0; i < mm.nCluster; i++) {
		// Leapfrog Verlet algorithm
		p = mm.cluster + i;
		if (mm.free_base && mm.base == p) {
			mdpar.tp_mh[i] = 0; mdpar.tp_ph[i] = 0; p->km.tpp = 0; p->km.tp = 0;
			continue;
		}
		mdpar.tp_ph[i] = mdpar.tp_mh[i] + p->km.tpp * mdpar.dt;
		mdpar.tp[i] = 0.5 * (mdpar.tp_mh[i] + mdpar.tp_ph[i]);
		p->km.tp = mdpar.tp[i];
	}

	// check the spatial moment with cordinate at molecular mass center
	VECTOR<6> Iw, Iac0;
	if (method == 0) { // using mass center
		// not necessary for acceleration !
		//calClusterAlpha(mm);
		//calSpatialAccel0(mm, Iac0);

		for (i = 0; i < 6; i++) {
			mdpar.P0_ph.v[i] = mdpar.P0_mh.v[i] + (mm.mDynKin.f.v[i] - mm.mDynKin.b.v[i]) * mdpar.dt; // M * alpha + b = f
			mdpar.P0.v[i] = (mdpar.P0_mh.v[i] + mdpar.P0_ph.v[i]) * 0.5;
			mm.mDynKin.P.v[i] = mdpar.P0.v[i];
		}
		MpV6(mm.mDynKin.invI, mdpar.P0, mdpar.V0)
		mm.mDynKin.V = mdpar.V0;

		calClusterVelocity0(mm); // calculate velocity of each cluster, with current base velocity, although it is not right !
		calSpatialMomentum(mm); // calculate the total spatial momentum based on current base velocity, tp[], and save momentum to mm.Iw

		if (mm.free_base) {
			// check velocity of base cluster
			// ## we might not have to double-check the base cluster velocity from the spatial moment ##
			// ## some time even get the velocity and spatial moment not so consistent ##
			// ## but if we turn feed back on the kinetic energy, constant thermal static energy
			// ## this checking will help the molecule's kinetic energy does not increase suddently
			mm.CalFreeBaseVelocity_GivenMassCenterMomentum(); // use mm.mDynkin.P to calculate right momentum relative to base cluster, use mm.Iw (calculated above) to compensate base velocity
			//mm.CalMassCenterVelocity_MassCenterMomentum();   //update already above
			calClusterVelocity0(mm); // update velocity of each cluster
			//calSpatialMomentum(mm); // update momentum again, mm.Iw, relative to base cluster, actually mm.Iw was update in CalFreeBaseVelocity_GivenMassCenterMomentum()

			for (i = 0; i < 6; i++) {
				mdpar.v0.v[i] = mm.V0.v[i];
				mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i]  + mm.alpha0.v[i] * mdpar.dt * 0.5; // using for verlet
				//mdpar.v0_mh.v[i] = mdpar.v0.v[i]  - mm.alpha0.v[i] * mdpar.dt * 0.5;
			}
		}
	}
	else if (method == 1) {
		calClusterVelocity0(mm);
		calSpatialMomentum0(mm, true);
		mm.CalMassCenterVelocity_MassCenterMomentum();
	}

	// for debug purpose
	/*
	calClusterVelocity0(mm);
	calSpatialMoment0(mm, Iw);
	i = i;
	*/
}

void Speed_Kinetic(LFVERLET_MM_MD &mdpar, MMOLECULE &mm, bool check_freebase) {
	int method = ::MM_SpeedVerlet_method; // 0 -- use the spatial moment at mass center, 1 -- use the velocity of base cluster
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0, j;
	CLUSTER *p = NULL;

	if (method == 1 && mm.free_base) { // Iw of base cluster, not implemented properly yet !
		for (i = 0; i < 6; i++) {
			mdpar.v0.v[i] = mm.V0.v[i];
			mdpar.v0_mh.v[i] = mm.V0.v[i] - mm.alpha0.v[i] * mdpar.dt * 0.5;
			mdpar.v0_ph.v[i] = mm.V0.v[i] + mm.alpha0.v[i] * mdpar.dt * 0.5;
		}
	}

	for (i = 0; i < mm.nCluster; i++) {
		// Leapfrog Verlet algorithm
		p = mm.cluster + i;
		mdpar.tp[i] = p->km.tp;
		mdpar.tp_mh[i] = mdpar.tp[i] - p->km.tpp * mdpar.dt * 0.5;
		mdpar.tp_ph[i] = mdpar.tp[i] + p->km.tpp * mdpar.dt * 0.5;
	}

	VECTOR<6> Iw, Iac0;
	if (method == 0 && mm.free_base) { // mass center
		for (i = 0; i < 6; i++) {
			mdpar.P0_mh.v[i] = mm.mDynKin.P.v[i] - (mm.mDynKin.f.v[i] - mm.mDynKin.b.v[i])* mdpar.dt * 0.5;
			mdpar.P0_ph.v[i] = mm.mDynKin.P.v[i] +  (mm.mDynKin.f.v[i] - mm.mDynKin.b.v[i])  * mdpar.dt * 0.5;
			mdpar.P0.v[i] = mm.mDynKin.P.v[i];
			mdpar.V0.v[i] = mm.mDynKin.V.v[i]; // not necessary
		}
	}
	//char msg[256] = "\0";
	//sprintf(msg, "P [%f, %f, %f, %f, %f, %f]", mm.mDynKin.P.v[0], mm.mDynKin.P.v[1], mm.mDynKin.P.v[2], mm.mDynKin.P.v[3], mm.mDynKin.P.v[4], mm.mDynKin.P.v[5]);
	//show_infor(msg);
	//sprintf(msg, "F [%f, %f, %f, %f, %f, %f]", Iac0.v[0], Iac0.v[1], Iac0.v[2], Iac0.v[3], Iac0.v[4], Iac0.v[5]);
	//show_infor(msg);

	// check velocity of base cluster
	if (mm.free_base && check_freebase && method == 0) {
		calClusterVelocity0(mm);
		calSpatialMomentum(mm);
		mm.CalFreeBaseVelocity_GivenMassCenterMomentum();
		calClusterVelocity0(mm);

		for (i = 0; i < 6; i++) {
			mdpar.v0.v[i] = mm.V0.v[i];
			mdpar.v0_mh.v[i] = mm.V0.v[i] - mm.alpha0.v[i] * mdpar.dt * 0.5;
			mdpar.v0_ph.v[i] = mm.V0.v[i] + mm.alpha0.v[i] * mdpar.dt * 0.5;
		}
	}
	// very important to construct the velocity history again
	// because free base V0 is changed
	/*
	if (mm.free_base) {
		for (i = 0; i < 6; i++) {
			mdpar.v0.v[i] = mm.V0.v[i];
			mdpar.v0_mh.v[i] = mm.V0.v[i] - mm.alpha0.v[i] * mdpar.dt * 0.5;
			mdpar.v0_ph.v[i] = mm.V0.v[i] + mm.alpha0.v[i] * mdpar.dt * 0.5;
		}
	}
	*/
}
/*
void Verlet(LFVERLET_MM_MD &mdpar, MMOLECULE &mm) {
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *p = NULL;
	double delta = 0;
	VECTOR<6> d;
	VECTOR3 dr, drot;

	// guess the spatial moment with cordinate at molecule mass center for next step
	// do spatial moment first. speed and conformation will change
	VECTOR<6> Iw0, Iac0;
	if (mm.free_base) {
		calClusterAlpha(mm);
		calSpatialAccel0(mm, Iac0);
		//calSpatialMoment0(mm, Iw0);
		for (i = 0; i < 6; i++) {
			mdpar.Iw0_ph.v[i] = mdpar.Iw0_mh.v[i] + Iac0.v[i] * mdpar.dt;
			mm.mDynKin.P.v[i] = 1.5 * mdpar.Iw0_ph.v[i] - 0.5 * mdpar.Iw0_mh.v[i]; 
			mdpar.Iw0_mh.v[i] = mdpar.Iw0_ph.v[i];
		}
	}

	// move base cluster and guess velocity for next step
	// because the base velocity is re-configured with the spacial moment at mass center,
	// the base velocity is hard to follow the leaflog algorithem

	if (mm.free_base) {
		for (i = 0; i < 6; i++) {
			// guess base cluster velocity
			//mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + mm.alpha0.v[i] * mdpar.dt;
			mdpar.v0_ph.v[i] = mm.V0.v[i] + mm.alpha0.v[i] * mdpar.dt * 0.5;

			//mdpar.v0.v[i] = 1.5 * mdpar.v0_ph.v[i] - 0.5 * mdpar.v0_mh.v[i]; 
			//mm.V0.v[i] = mdpar.v0.v[i];
			d.v[i] = mdpar.v0_ph.v[i] * mdpar.dt;
			mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; //for next step
			//mdpar.v0_mh.v[i] = mm.V0.v[i] + mm.alpha0.v[i] * mdpar.dt * 0.5; //for next step
		}
		V2wv(drot, dr, d)
	}

	if (mm.nCluster < MAX_THREADS) {
	//*******************************************************
	//			rotation all clusters in series
	//		NOTE: lines started with // is commented
	//*******************************************************
	
		//mm.twMol(*(mm.base->Op), drot, true);
		//mm.shiftMol(dr); // translation was done when shift molecule
		mm.shift_tw_mol(dr, *(mm.base->Op), drot, true);
	
		for (i = 0; i < 6; i++) mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; 
	
		//for (i = 0; i < 3; i++) mdpar.r0.v[i] = mm.mDynKin.r0.v[i];

		// move other clusters and guess velocity for next step
		for (i = 0; i < mm.nCluster; i++) {
			p = mm.cluster + i;
			if (mm.free_base && p == mm.base) continue;
			//delta = torsion_angle(p->parent, (BASIC_CLUSTER*)p); // check torsion angle
			// Leapfrog Verlet algorithm
			delta = mdpar.tp_ph[i] * mdpar.dt;
			mdpar.t[i] += delta;
			p->km.ta = mdpar.t[i];
			mm.twTA(i, delta, true); // in arc degree
		
			//delta = torsion_angle(p->parent, p); // check torsion angle

			// guess new position velocity
			mdpar.tp[i] = 1.5 * mdpar.tp_ph[i] - 0.5 * mdpar.tp_mh[i];
			p->km.tp = mdpar.tp[i];
			mdpar.tp_mh[i] = mdpar.tp_ph[i];
		}
	}
	else {
	//*******************************************************
	//			rotation all clusters parallelly
	//		NOTE: lines started with // is commented
	//*******************************************************
		// move all clusters 

		mm.setup_base_movement(dr, *(mm.base->Op), drot, true);

		for (i = 0; i < mm.nCluster; i++) {
			p = mm.cluster + i;
			if (mm.free_base && p == mm.base) continue;
			//delta = torsion_angle(p->parent, (BASIC_CLUSTER*)p); // check torsion angle
			// Leapfrog Verlet algorithm
			delta = mdpar.tp_ph[i] * mdpar.dt;
			p->rot.dta = delta;
		}
		setup_all_LocalRotMatrix(mm);
		setup_GlobalRotMatrix(mm);
		tweak_all_clusters(mm);

		//for (i = 0; i < 3; i++) mdpar.r0.v[i] = mm.mDynKin.r0.v[i];

		for (i = 0; i < mm.nCluster; i++) {
			p = mm.cluster + i;
			if (p == mm.base) continue;
			mdpar.t[i] += p->rot.dta;
			p->km.ta = mdpar.t[i];
		
			//delta = torsion_angle(p->parent, p); // check torsion angle

			// guess new position velocity
			mdpar.tp[i] = 1.5 * mdpar.tp_ph[i] - 0.5 * mdpar.tp_mh[i];
			p->km.tp = mdpar.tp[i];
			mdpar.tp_mh[i] = mdpar.tp_ph[i];
		}
	}
}
*/

void Verlet(LFVERLET_MM_MD &mdpar, MMOLECULE &mm) {
	int method = ::MM_SpeedVerlet_method; // 0 -- use the spatial moment at mass center, 1 -- use the velocity of base cluster
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *p = NULL;
	double delta = 0;
	VECTOR<6> d;
	VECTOR3 dr, drot;

	// guess the spatial moment with cordinate at molecule mass center for next step
	// do spatial moment first. speed and conformation will change
	VECTOR<6> Iw0, Iac0;
	if (method == 0 && mm.free_base) { // mass center
		//calClusterAlpha(mm);
		//calSpatialAccel0(mm, Iac0);
		for (i = 0; i < 6; i++) {
			mdpar.P0_ph.v[i] = mdpar.P0_mh.v[i] + (mm.mDynKin.f.v[i] - mm.mDynKin.b.v[i]) * mdpar.dt;
			mm.mDynKin.P.v[i] = 1.5 * mdpar.P0_ph.v[i] - 0.5 * mdpar.P0_mh.v[i]; 
			mdpar.P0_mh.v[i] = mdpar.P0_ph.v[i];
		}
	}

	// move base cluster and guess velocity for next step
	// because the base velocity is re-configured with the spacial moment at mass center,
	// the base velocity is hard to follow the leaflog algorithem

	if (mm.free_base) {
		for (i = 0; i < 6; i++) {
			// guess base cluster velocity
			if (method == 1) mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + mm.alpha0.v[i] * mdpar.dt;
			else if (method == 0) mdpar.v0_ph.v[i] = mm.V0.v[i] + mm.alpha0.v[i] * mdpar.dt * 0.5;

			d.v[i] = mdpar.v0_ph.v[i] * mdpar.dt;

			if (method == 1) {
				mdpar.v0.v[i] = 1.5 * mdpar.v0_ph.v[i] - 0.5 * mdpar.v0_mh.v[i]; 
				mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; //for next step
			}
			else if (method == 0) {
				mdpar.v0.v[i] = 1.5 * mdpar.v0_ph.v[i] - 0.5 * mdpar.v0_mh.v[i]; 
				mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; //for next step
			}
			mm.V0.v[i] = mdpar.v0.v[i];
		}
		V62V6(mm.V0, mm.base->bkm->V0)
		V2wv(drot, dr, d)
	}
	mm.set_base_move(dr, drot);

	for (i = 0; i < mm.nCluster; i++) {
		p = mm.cluster + i;
		if (mm.free_base && p == mm.base) continue;
		//delta = torsion_angle(p->parent, (BASIC_CLUSTER*)p); // check torsion angle
		// Leapfrog Verlet algorithm
		mdpar.tp_ph[i] = mdpar.tp_mh[i] + p->km.tpp * mdpar.dt;
		delta = mdpar.tp_ph[i] * mdpar.dt;
		p->rot.dta = delta;
	}
	
	MacroMolMove(&mm);

	for (i = 0; i < mm.nCluster; i++) {
		p = mm.cluster + i;
		if (p == mm.base) continue;
		mdpar.t[i] += p->rot.dta;
		p->km.ta = mdpar.t[i];
		
		//delta = torsion_angle((BASIC_CLUSTER*)p->parent, p); // check torsion angle

		// guess new position velocity
		mdpar.tp[i] = 1.5 * mdpar.tp_ph[i] - 0.5 * mdpar.tp_mh[i];
		p->km.tp = mdpar.tp[i];
		mdpar.tp_mh[i] = mdpar.tp_ph[i];
	}
}

void Verlet_Relax(LFVERLET_MM_MD &mdpar, MMOLECULE &mm, double dta_max) {
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *p = NULL;
	double delta = 0;
	VECTOR<6> d;
	VECTOR3 dr, drot;

	double tp_max = 0, t, dt = 0;
	for (i = 0; i < mm.nCluster; i++) {
		if ((t = FABS(mdpar.tp_ph[i])) > tp_max) {tp_max = t; dt = dta_max / tp_max;}
	}

	// guess the spatial moment with cordinate at molecule mass center for next step
	// do spatial moment first. speed and conformation will change
	VECTOR<6> Iw0, Iac0;
	if (mm.free_base) {
		calClusterAlpha(mm);
		calSpatialAccel0(mm, Iac0);
		calSpatialMomentum0(mm, true);
		for (i = 0; i < 6; i++) {
			mdpar.Iw_ph.v[i] = mdpar.Iw_mh.v[i] + Iac0.v[i] * mdpar.dt;
			mm.mDynKin.P.v[i] = 1.5 * mdpar.Iw_ph.v[i] - 0.5 * mdpar.Iw_mh.v[i]; 
			mdpar.Iw_mh.v[i] = mdpar.Iw_ph.v[i];
		}
	}

	// move base cluster and guess velocity for next step
	if (mm.free_base) {
		for (i = 0; i < 6; i++) {
			// guess base cluster velocity
			mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + mm.alpha0.v[i] * mdpar.dt;
			mdpar.v0.v[i] = 1.5 * mdpar.v0_ph.v[i] - 0.5 * mdpar.v0_mh.v[i]; 
			mm.V0.v[i] = mdpar.v0.v[i];
			//d.v[i] = mdpar.v0_ph.v[i] * mdpar.dt;
			d.v[i] = mdpar.v0_ph.v[i] * dt;
		}
		V2wv(drot, dr, d)
	}

	if (mm.nCluster < MAX_THREADS) {
	/*******************************************************
				rotation all clusters in series
			NOTE: lines started with // is commented
	*******************************************************/
	
		//mm.twMol(*(mm.base->Op), drot, true);
		//mm.shiftMol(dr); // translation was done when shift molecule
		mm.shift_tw_mol(dr, *(mm.base->Op), drot, true);
	
		for (i = 0; i < 6; i++) mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; 
	
		//for (i = 0; i < 3; i++) mdpar.r0.v[i] = mm.mDynKin.r0.v[i];

		// move other clusters and guess velocity for next step
		for (i = 0; i < mm.nCluster; i++) {
			p = mm.cluster + i;
			if (mm.free_base && p == mm.base) continue;
			//delta = torsion_angle(p->parent, (BASIC_CLUSTER*)p); // check torsion angle
			// Leapfrog Verlet algorithm
			//delta = mdpar.tp_ph[i] * mdpar.dt;
			delta = mdpar.tp_ph[i] * dt;
			mdpar.t[i] += delta;
			p->km.ta = mdpar.t[i];
			mm.twTA(i, delta, true); // in arc degree
		
			//delta = torsion_angle(p->parent, p); // check torsion angle

			// guess new position velocity
			mdpar.tp[i] = 1.5 * mdpar.tp_ph[i] - 0.5 * mdpar.tp_mh[i];
			p->km.tp = mdpar.tp[i];
			mdpar.tp_mh[i] = mdpar.tp_ph[i];
		}
	}
	else {
	/*******************************************************
				rotation all clusters parallelly
			NOTE: lines started with // is commented
	*******************************************************/
		// move all clusters 

		mm.setup_base_movement(dr, *(mm.base->Op), drot, true);

		for (i = 0; i < mm.nCluster; i++) {
			p = mm.cluster + i;
			if (mm.free_base && p == mm.base) continue;
			//delta = torsion_angle(p->parent, (BASIC_CLUSTER*)p); // check torsion angle
			// Leapfrog Verlet algorithm
			
			//delta = mdpar.tp_ph[i] * mdpar.dt;
			delta = mdpar.tp_ph[i] * dt;
			
			p->rot.dta = delta;
		}
		setup_all_LocalRotMatrix(mm);
		setup_GlobalRotMatrix(mm);
		tweak_all_clusters(mm);

		// guess velocity for next step
		for (i = 0; i < 6; i++) mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; 
	
		//for (i = 0; i < 3; i++) mdpar.r0.v[i] = mm.mDynKin.r0.v[i];

		for (i = 0; i < mm.nCluster; i++) {
			p = mm.cluster + i;
			if (mm.free_base && p == mm.base) continue;
			mdpar.t[i] += p->rot.dta;
			p->km.ta = mdpar.t[i];
		
			//delta = torsion_angle(p->parent, p); // check torsion angle

			// guess new position velocity
			mdpar.tp[i] = 1.5 * mdpar.tp_ph[i] - 0.5 * mdpar.tp_mh[i];
			p->km.tp = mdpar.tp[i];
			mdpar.tp_mh[i] = mdpar.tp_ph[i];
		}
	}
}

void InitKinetic(LFVERLET_MM_MD &mdpar, MMOLECULE &mm) {
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *p = NULL;
	double dv = 0;

	// base cluster velocity
	//for (i = 0; i < 3; i++) mdpar.r0.v[i] = mm.mDynKin.r0.v[i];
	for (i = 0; i < 6; i++) {
		mdpar.v0.v[i] = mm.V0.v[i];
		dv = mm.alpha0.v[i] * mdpar.dt * 0.5;
		mdpar.v0_mh.v[i] = mm.V0.v[i] - dv;
	}
	// cluster's local parameters
	for (i = 0; i < mm.nCluster; i++) {
		p = mm.cluster + i;
		// Leapfrog Verlet algorithm
		// guess tp_mh
		mdpar.tp[i] = p->km.tp;
		mdpar.tp_mh[i] = mdpar.tp[i] - p->km.tpp * mdpar.dt * 0.5;
	}

	int method = ::MM_SpeedVerlet_method; // 0 -- use the spatial moment at mass center, 1 -- use the velocity of base cluster

	VECTOR<6> Iw0, Iac0;
	if (method == 0) { // mass center
		//calClusterAlpha(mm);
		//calSpatialAccel0(mm, Iac0);
		for (i = 0; i < 6; i++) {
			mdpar.P0.v[i] = mm.mDynKin.P.v[i];
			mdpar.P0_mh.v[i] = mm.mDynKin.P.v[i] - (mm.mDynKin.f.v[i] - mm.mDynKin.b.v[i]) * mdpar.dt * 0.5;
			mdpar.P0_ph.v[i] = mm.mDynKin.P.v[i] + (mm.mDynKin.f.v[i] - mm.mDynKin.b.v[i]) * mdpar.dt * 0.5;
		}
	}
	else if (method == 1) { // base cluster, not implemented properly !!
		//calClusterAlpha(mm);
		//calSpatialAccel0(mm, Iac0);
		for (i = 0; i < 6; i++) {
			mdpar.Iw.v[i] = mm.Iw.v[i];
			mdpar.Iw_mh.v[i] = mm.Iw.v[i] - mm.f.v[i] * mdpar.dt * 0.5;
			mdpar.Iw_ph.v[i] = mm.Iw.v[i] + mm.f.v[i] * mdpar.dt * 0.5;
		}
	}
}

bool LFVERLET_MM_MD::save_mol_dyn_pars(MMOLECULE &mm, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	if (mmsave.mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}
	if (check_status) {
		if (!mmsave.mdsave->mDynSave.one_more()) return false;
	}
	if (!mmsave.mdsave->mDynSave.bsave()) return true; // this loop will be not saved

	int i = 0;
	if (bAsync) {
		mmsave.mdsave->init_buf();
		if (save_md_title) mmsave.mdsave->mDynSave.sbuf<<"[MD] "<<mmsave.mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) mmsave.mdsave->mDynSave.sbuf<<additional_title<<endl;

	
	/*
	VECTOR<6> Iw;
	calSpatialMoment0(mm, Iw);
	*(mmsave.mdsave->mDynSave.out)<<Iw.v[0];
	for (i = 1; i < 6; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<Iw.v[i];
	*(mmsave.mdsave->mDynSave.out)<<endl;

	VECTOR<6> Iac;
	calSpatialAccel0(mm, Iac);
	*(mmsave.mdsave->mdsave->mDynSave.out)<<Iac.v[0];
	for (i = 1; i < 6; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<Iac.v[i];
	*(mmsave.mdsave->mDynSave.out)<<endl;
	*/

	/*
	int nc = 0;
	for (nc = 0; nc < mm.nCluster; nc++) {
		*(mmsave.mdsave->mDynSave.out)<<nc;
		for (i = 0; i < 6; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<mm.cluster[nc].dyn->fc.v[i];
		*(mmsave.mdsave->mDynSave.out)<<endl;
	}
	*/
	/*
	*(mmsave.mdsave->mDynSave.out)<<mm.mDynKin.cm.v[0];
	for (i = 1; i < 3; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<mm.mDynKin.cm.v[i];
	*(mmsave.mdsave->mDynSave.out)<<endl;
	*/

		mmsave.mdsave->mDynSave.sbuf<<this->ksi<<"  "<<vksi<<endl;
		mmsave.mdsave->mDynSave.sbuf<<endl;

		mmsave.mdsave->mDynSave.flush_buf();
	}
	else {
		if (save_md_title) *(mmsave.mdsave->mDynSave.out)<<"[MD] "<<mmsave.mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) *(mmsave.mdsave->mDynSave.out)<<additional_title<<endl;

	
	/*
	VECTOR<6> Iw;
	calSpatialMoment0(mm, Iw);
	*(mmsave.mdsave->mDynSave.out)<<Iw.v[0];
	for (i = 1; i < 6; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<Iw.v[i];
	*(mmsave.mdsave->mDynSave.out)<<endl;

	VECTOR<6> Iac;
	calSpatialAccel0(mm, Iac);
	*(mmsave.mdsave->mdsave->mDynSave.out)<<Iac.v[0];
	for (i = 1; i < 6; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<Iac.v[i];
	*(mmsave.mdsave->mDynSave.out)<<endl;
	*/

	/*
	int nc = 0;
	for (nc = 0; nc < mm.nCluster; nc++) {
		*(mmsave.mdsave->mDynSave.out)<<nc;
		for (i = 0; i < 6; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<mm.cluster[nc].dyn->fc.v[i];
		*(mmsave.mdsave->mDynSave.out)<<endl;
	}
	*/
	/*
	*(mmsave.mdsave->mDynSave.out)<<mm.mDynKin.cm.v[0];
	for (i = 1; i < 3; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<mm.mDynKin.cm.v[i];
	*(mmsave.mdsave->mDynSave.out)<<endl;
	*/

		*(mmsave.mdsave->mDynSave.out)<<this->ksi<<"  "<<vksi<<endl;
		*(mmsave.mdsave->mDynSave.out)<<endl;
	}

	if (check_status) mmsave.mdsave->mDynSave.save_over();
	return true;
}


// this function gives the movement parameters, and change mdpar only, but does not move the macro-molecule
void VerletStep(LFVERLET_MM_MD &mdpar, MMOLECULE &mm, VECTOR<6> &dv_base, double *dta) {
	int method = ::MM_SpeedVerlet_method; // 0 -- use the spatial moment at mass center, 1 -- use the velocity of base cluster
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *p = NULL;
	double delta = 0;
	VECTOR<6> d;
	VECTOR3 dr, drot;

	// guess the spatial moment with cordinate at molecule mass center for next step
	// do spatial moment first. speed and conformation will change
	VECTOR<6> Iw0, Iac0;
	if (method == 0 && mm.free_base) {
		calClusterAlpha(mm);
		calSpatialAccel0(mm, Iac0);
		//calSpatialMoment0(mm, Iw0);
		for (i = 0; i < 6; i++) {
			mdpar.Iw_ph.v[i] = mdpar.Iw_mh.v[i] + Iac0.v[i] * mdpar.dt;
			mm.mDynKin.P.v[i] = 1.5 * mdpar.Iw_ph.v[i] - 0.5 * mdpar.Iw_mh.v[i]; 
			mdpar.Iw_mh.v[i] = mdpar.Iw_ph.v[i];
		}
	}

	// move base cluster and guess velocity for next step
	// because the base velocity is re-configured with the spacial moment at mass center,
	// the base velocity is hard to follow the leaflog algorithem

	if (mm.free_base) {
		for (i = 0; i < 6; i++) {
			// guess base cluster velocity
			if (method == 1) mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + mm.alpha0.v[i] * mdpar.dt;
			else if (method == 0) mdpar.v0_ph.v[i] = mm.V0.v[i] + mm.alpha0.v[i] * mdpar.dt * 0.5;

			mdpar.v0.v[i] = 1.5 * mdpar.v0_ph.v[i] - 0.5 * mdpar.v0_mh.v[i]; 
			mm.V0.v[i] = mdpar.v0.v[i];
			d.v[i] = mdpar.v0_ph.v[i] * mdpar.dt;
			mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; //for next step
			//mdpar.v0_mh.v[i] = mm.V0.v[i] + mm.alpha0.v[i] * mdpar.dt * 0.5; //for next step
		}
		V2wv(drot, dr, d)
	}
	V62V6(d, dv_base)

	for (i = 0; i < 6; i++) mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; 
	
	// move other clusters and guess velocity for next step
	for (i = 0; i < mm.nCluster; i++) {
		p = mm.cluster + i;
		if (mm.free_base && p == mm.base) {dta[i] = 0; continue;}
		// Leapfrog Verlet algorithm
		delta = mdpar.tp_ph[i] * mdpar.dt;
		mdpar.t[i] += delta;
		//p->km.ta = mdpar.t[i];
		dta[i] = delta;

		// guess new position velocity
		mdpar.tp[i] = 1.5 * mdpar.tp_ph[i] - 0.5 * mdpar.tp_mh[i];
		p->km.tp = mdpar.tp[i];
		mdpar.tp_mh[i] = mdpar.tp_ph[i];
	}
}

// this function move the macro-molecule with given movement parameters, and mdpar
void VerletMove(LFVERLET_MM_MD &mdpar, VECTOR<6> &dv_base, double *dta, MMOLECULE &mm) {
	//suppose mdpar.nCluster == mm.nCluster
	int i = 0;
	CLUSTER *p = NULL;
	double delta = 0;
	VECTOR3 dr, drot;

	V2wv(drot, dr, dv_base)

	if (mm.nCluster < MAX_THREADS) {
	/*******************************************************
				rotation all clusters in series
			NOTE: lines started with // is commented
	*******************************************************/
	
		mm.shift_tw_mol(dr, *(mm.base->Op), drot, true);
	
		// move other clusters and guess velocity for next step
		for (i = 0; i < mm.nCluster; i++) {
			p = mm.cluster + i;
			if (mm.free_base && p == mm.base) continue;
			mm.twTA(i, dta[i], true); // in arc degree
			p->km.ta += dta[i];
			//p->km.tp = mdpar.tp[i];
		}
	}
	else {
	/*******************************************************
				rotation all clusters parallelly
			NOTE: lines started with // is commented
	*******************************************************/
		// move all clusters 

		mm.setup_base_movement(dr, *(mm.base->Op), drot, true);

		for (i = 0; i < mm.nCluster; i++) {
			p = mm.cluster + i;
			if (mm.free_base && p == mm.base) {p->rot.dta = 0; continue;}
			p->rot.dta = dta[i];
		}
		setup_all_LocalRotMatrix(mm);
		setup_GlobalRotMatrix(mm);
		tweak_all_clusters(mm);

		for (i = 0; i < mm.nCluster; i++) {
			p = mm.cluster + i;
			if (p == mm.base) continue;
			p->km.ta += dta[i];
		}
	}
}


// MD FUNCTIONS FOR SMOLECULE

bool SM_MD_SAVE::save_kinetic_pars(SMOLECULE &m, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	if (mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!mdsave->mKinSave.one_more()) return false;
	}
	if (!mdsave->mKinSave.bsave()) return true; // this loop will be not saved

	if (bAsync) {
		mdsave->init_buf();
		if (save_md_title) mdsave->mKinSave.sbuf<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			mdsave->mKinSave.sbuf<<additional_title<<endl;

		mdsave->mKinSave.sbuf<<m.r.v[0]<<" "<<m.r.v[1]<<" "<<m.r.v[2]<<" ";
		mdsave->mKinSave.sbuf<<m.zaxis.v[0]<<" "<<m.zaxis.v[1]<<" "<<m.zaxis.v[2]<<" ";
		mdsave->mKinSave.sbuf<<m.xaxis.v[0]<<" "<<m.xaxis.v[1]<<" "<<m.xaxis.v[2]<<endl;

		mdsave->mKinSave.sbuf<<m.c->bkm->V0.v[0]<<" "<<m.c->bkm->V0.v[1]<<" "<<m.c->bkm->V0.v[2]<<" ";
		mdsave->mKinSave.sbuf<<m.c->bkm->V0.v[3]<<" "<<m.c->bkm->V0.v[4]<<" "<<m.c->bkm->V0.v[5]<<endl;

		mdsave->mKinSave.sbuf<<m.c->bkm->alpha.v[0]<<" "<<m.c->bkm->alpha.v[1]<<" "<<m.c->bkm->alpha.v[2]<<" ";
		mdsave->mKinSave.sbuf<<m.c->bkm->alpha.v[3]<<" "<<m.c->bkm->alpha.v[4]<<" "<<m.c->bkm->alpha.v[5]<<endl;

		mdsave->mKinSave.sbuf<<endl;

		mdsave->mKinSave.flush_buf();
	}
	else {
		if (save_md_title) *(mdsave->mKinSave.out)<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			*(mdsave->mKinSave.out)<<additional_title<<endl;

		*(mdsave->mKinSave.out)<<m.r.v[0]<<" "<<m.r.v[1]<<" "<<m.r.v[2]<<" ";
		*(mdsave->mKinSave.out)<<m.zaxis.v[0]<<" "<<m.zaxis.v[1]<<" "<<m.zaxis.v[2]<<" ";
		*(mdsave->mKinSave.out)<<m.xaxis.v[0]<<" "<<m.xaxis.v[1]<<" "<<m.xaxis.v[2]<<endl;

		*(mdsave->mKinSave.out)<<m.c->bkm->V0.v[0]<<" "<<m.c->bkm->V0.v[1]<<" "<<m.c->bkm->V0.v[2]<<" ";
		*(mdsave->mKinSave.out)<<m.c->bkm->V0.v[3]<<" "<<m.c->bkm->V0.v[4]<<" "<<m.c->bkm->V0.v[5]<<endl;

		*(mdsave->mKinSave.out)<<m.c->bkm->alpha.v[0]<<" "<<m.c->bkm->alpha.v[1]<<" "<<m.c->bkm->alpha.v[2]<<" ";
		*(mdsave->mKinSave.out)<<m.c->bkm->alpha.v[3]<<" "<<m.c->bkm->alpha.v[4]<<" "<<m.c->bkm->alpha.v[5]<<endl;

		*(mdsave->mKinSave.out)<<endl;
	}

	if (check_status) mdsave->mKinSave.save_over();
	return true;
}


bool SM_MD_SAVE::save_atomic_dipole(SMOLECULE &sm, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	return true;
	/*
	if (!::bPolarize && !::bDipole) return true;

	if (mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!mdsave->muSave.one_more()) return false;
	}
	if (!mdsave->muSave.bsave()) return true; // this loop will be not saved

	int i = 0;
	double ta = 0;
	int ia;
	BASIC_CLUSTER *c = sm.c;


	if (bAsync) {
		mdsave->init_buf();
		if (save_md_title) mdsave->muSave.sbuf<<"[MD] "<<mdsave->muSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0)
			mdsave->muSave.sbuf<<additional_title<<endl;

		for (ia = 0; ia < c->nAtoms; ia++) {
			mdsave->muSave.sbuf<<c->atom[ia].mu.v[0]<<" "<<c->atom[ia].mu.v[1]<<" "<<c->atom[ia].mu.v[2]<<" ";
			//mdsave->muSave.sbuf<<c->atom[ia].F.v[0]<<" "<<c->atom[ia].F.v[1]<<" "<<c->atom[ia].F.v[2]<<" ";
			if (::bDipole) {
				mdsave->muSave.sbuf<<c->atom[ia].ind_mu.v[0]<<" "<<c->atom[ia].ind_mu.v[1]<<" "<<c->atom[ia].ind_mu.v[2]<<" ";
				mdsave->muSave.sbuf<<c->atom[ia].intr_mu.v[0]<<" "<<c->atom[ia].intr_mu.v[1]<<" "<<c->atom[ia].intr_mu.v[2];
			}
			mdsave->muSave.sbuf<<endl;
		}
		//*(mdsave->muSave.out)<<c->dyn->fc.v[0]<<" "<<c->dyn->fc.v[1]<<" "<<c->dyn->fc.v[2]<<" ";
		//*(mdsave->muSave.out)<<c->dyn->fc.v[3]<<" "<<c->dyn->fc.v[4]<<" "<<c->dyn->fc.v[5]<<" "<<endl;
	
		mdsave->muSave.sbuf<<endl;

		mdsave->muSave.flush_buf();
	}
	else {
		if (save_md_title) *(mdsave->muSave.out)<<"[MD] "<<mdsave->muSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0)
			*(mdsave->muSave.out)<<additional_title<<endl;

		for (ia = 0; ia < c->nAtoms; ia++) {
			*(mdsave->muSave.out)<<c->atom[ia].mu.v[0]<<" "<<c->atom[ia].mu.v[1]<<" "<<c->atom[ia].mu.v[2]<<" ";
			//*(mdsave->muSave.out)<<c->atom[ia].F.v[0]<<" "<<c->atom[ia].F.v[1]<<" "<<c->atom[ia].F.v[2]<<" ";
			if (::bDipole) {
				*(mdsave->muSave.out)<<c->atom[ia].ind_mu.v[0]<<" "<<c->atom[ia].ind_mu.v[1]<<" "<<c->atom[ia].ind_mu.v[2]<<" ";
				*(mdsave->muSave.out)<<c->atom[ia].intr_mu.v[0]<<" "<<c->atom[ia].intr_mu.v[1]<<" "<<c->atom[ia].intr_mu.v[2];
			}
			*(mdsave->muSave.out)<<endl;
		}
		//*(mdsave->muSave.out)<<c->dyn->fc.v[0]<<" "<<c->dyn->fc.v[1]<<" "<<c->dyn->fc.v[2]<<" ";
		//*(mdsave->muSave.out)<<c->dyn->fc.v[3]<<" "<<c->dyn->fc.v[4]<<" "<<c->dyn->fc.v[5]<<" "<<endl;
	
		*(mdsave->muSave.out)<<endl;
	}

	if (check_status) mdsave->muSave.save_over();
	return true;
	*/
}

void Speed_Verlet(LFVERLET_SM_MD &mdpar, SMOLECULE &m) {
	int i = 0;

	// we suppose force cluster calculate accerlation based on total force, sm.c->dyn->f !!
	// and momentum is well maitained in bkm->P0

	// translation --
	for (i = 0; i < 6; i++) {
		//mdpar.v0.v[i] = m.c->bkm->V0.v[i];
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + m.c->bkm->alpha.v[i] * mdpar.dt;
		mdpar.v0.v[i] = 0.5 * (mdpar.v0_mh.v[i] + mdpar.v0_ph.v[i]);
		m.c->bkm->V0.v[i] = mdpar.v0.v[i];
	}
	// angular rotation --
	for (i = 0; i < 3; i++) {
		//mdpar.p0.v[i] = m.c->bkm->P0.v[i];
		mdpar.p0_ph.v[i] = mdpar.p0_mh.v[i] + m.c->dyn->f.v[i] * mdpar.dt;
		mdpar.p0.v[i] = 0.5 * (mdpar.p0_mh.v[i] + mdpar.p0_ph.v[i]);
		m.c->bkm->P0.v[i] = mdpar.p0.v[i];
	}
	// since P0 is changed, angular momentum of cluster has to be recalculated, to make sure consistent between angular momentum and its speed !
	VECTOR3 w;
	MpV3(m.invI, mdpar.p0, w)
	memcpy(mdpar.v0.v, w.v, SIZE_V3);
	memcpy(m.c->bkm->V0.v, w.v, SIZE_V3);
}

void Speed_Kinetic(LFVERLET_SM_MD &mdpar, SMOLECULE &m) {
	int i = 0;
	double dv = 0;

	double hdt = mdpar.dt * 0.5;
	for (i = 0; i < 6; i++) {
		mdpar.v0.v[i] = m.c->bkm->V0.v[i];
		dv = m.c->bkm->alpha.v[i] * hdt;
		mdpar.v0_mh.v[i] = mdpar.v0.v[i] - dv;
		mdpar.v0_ph.v[i] = mdpar.v0.v[i] + dv;
	}
	for (i = 0; i < 3; i++) {
		mdpar.p0.v[i] = m.c->bkm->P0.v[i];
		mdpar.p0_mh.v[i] = mdpar.p0.v[i] - m.c->dyn->f.v[i] * hdt;
		mdpar.p0_ph.v[i] = mdpar.p0.v[i] + m.c->dyn->f.v[i] * hdt;
	}
}

void Verlet(LFVERLET_SM_MD &mdpar, SMOLECULE &m) {
	if (m.bFixed) {
		memset(mdpar.v0_mh.v, 0, SIZE_V6);
		memset(mdpar.v0_ph.v, 0, SIZE_V6);
		memset(mdpar.v0.v, 0, SIZE_V6);
		memset(mdpar.p0.v, 0, SIZE_V3);
		memset(mdpar.p0_ph.v, 0, SIZE_V3);
		memset(mdpar.p0_mh.v, 0, SIZE_V3);
		memset(m.c->bkm->V0.v, 0, SIZE_V6);
		return;
	}

	int i = 0;
	double delta = 0;
	VECTOR<6> d;
	VECTOR3 dr, drot;

	for (i = 0; i < 6; i++) {
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + m.c->bkm->alpha.v[i] * mdpar.dt;
		mdpar.v0.v[i] = 1.5 * mdpar.v0_ph.v[i] - 0.5 * mdpar.v0_mh.v[i]; 
		m.c->bkm->V0.v[i] = mdpar.v0.v[i];
		d.v[i] = mdpar.v0_ph.v[i] * mdpar.dt;
	}
	for (i = 0; i < 3; i++) {
		mdpar.p0_ph.v[i] = mdpar.p0_mh.v[i] + m.c->dyn->f.v[i] * mdpar.dt;
		mdpar.p0.v[i] = 1.5 * mdpar.p0_ph.v[i] - 0.5 * mdpar.p0_mh.v[i];
		m.c->bkm->P0.v[i] = mdpar.p0.v[i];
	}
	V2wv(drot, dr, d)
	m.shiftMol(dr);
	rotate_SM(m, drot, true);
	
	for (i = 0; i < 6; i++) mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; 
	for (i = 0; i < 3; i++) mdpar.p0_mh.v[i] = mdpar.p0_ph.v[i];
}

bool LFVERLET_SM_MD::save_mol_dyn_pars(SMOLECULE &m, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	if (smsave.mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!smsave.mdsave->mDynSave.one_more()) return false;
	}
	if (!smsave.mdsave->mDynSave.bsave()) return true; // this loop will be not saved

	if (bAsync) {
		smsave.mdsave->init_buf();
		if (save_md_title) smsave.mdsave->mDynSave.sbuf<<"[MD] "<<smsave.mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			smsave.mdsave->mDynSave.sbuf<<additional_title<<endl;

		//smsave.mdsave->mDynSave.sbuf<<this->ksi_Hoover<<" "<<vksi_Hoover<<endl;
		if (this->hdyn != NULL) {
			smsave.mdsave->mDynSave.sbuf<<this->hdyn->nhc->ksi.v[0]<<" "<<this->hdyn->nhc->vksi.v[0]<<endl;
		}
		smsave.mdsave->mDynSave.sbuf<<endl;

		smsave.mdsave->mDynSave.flush_buf();
	}
	else {
		if (save_md_title) *(smsave.mdsave->mDynSave.out)<<"[MD] "<<smsave.mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			*(smsave.mdsave->mDynSave.out)<<additional_title<<endl;

		//*(smsave.mdsave->mDynSave.out)<<this->ksi_Hoover<<" "<<vksi_Hoover<<endl;
		if (this->hdyn != NULL) {
			*(smsave.mdsave->mDynSave.out)<<this->hdyn->nhc->ksi.v[0]<<" "<<this->hdyn->nhc->vksi.v[0]<<endl;
		}
		*(smsave.mdsave->mDynSave.out)<<endl;
	}

	if (check_status) smsave.mdsave->mDynSave.save_over();
	return true;
}



// MD functions for PMOLECULE

bool PM_MD_SAVE::save_kinetic_pars(PMOLECULE &m, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	if (mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!mdsave->mKinSave.one_more()) return false;
	}
	if (!mdsave->mKinSave.bsave()) return true; // this loop will be not saved

	if (bAsync) {
		mdsave->init_buf();
		if (save_md_title) mdsave->mKinSave.sbuf<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			mdsave->mKinSave.sbuf<<additional_title<<endl;
		mdsave->mKinSave.sbuf<<m.r.v[0]<<" "<<m.r.v[1]<<" "<<m.r.v[2]<<endl;

		mdsave->mKinSave.sbuf<<m.v.v[0]<<" "<<m.v.v[1]<<" "<<m.v.v[2]<<endl;

		mdsave->mKinSave.sbuf<<m.alpha.v[0]<<" "<<m.alpha.v[1]<<" "<<m.alpha.v[2]<<endl;

		mdsave->mKinSave.sbuf<<endl;

		mdsave->mKinSave.flush_buf();
	}
	else {
		if (save_md_title) *(mdsave->mKinSave.out)<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			*(mdsave->mKinSave.out)<<additional_title<<endl;
		*(mdsave->mKinSave.out)<<m.r.v[0]<<" "<<m.r.v[1]<<" "<<m.r.v[2]<<endl;

		*(mdsave->mKinSave.out)<<m.v.v[0]<<" "<<m.v.v[1]<<" "<<m.v.v[2]<<endl;

		*(mdsave->mKinSave.out)<<m.alpha.v[0]<<" "<<m.alpha.v[1]<<" "<<m.alpha.v[2]<<endl;

		*(mdsave->mKinSave.out)<<endl;
	}

	if (check_status) mdsave->mKinSave.save_over();
	return true;
}


bool PM_MD_SAVE::save_atomic_dipole(PMOLECULE &pm, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	return true;
	/*
	if (!::bPolarize && !::bDipole) return true;
	
	if (mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!mdsave->muSave.one_more()) return false;
	}
	if (!mdsave->muSave.bsave()) return true; // this loop will be not saved

	int i = 0;
	int ia;
	BASIC_CLUSTER *c = pm.c;

	if (bAsync) {
		mdsave->init_buf();
		if (save_md_title) mdsave->muSave.sbuf<<"[MD] "<<mdsave->muSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0)
			mdsave->muSave.sbuf<<additional_title<<endl;
		
		for (ia = 0; ia < c->nAtoms; ia++) {
			mdsave->muSave.sbuf<<c->atom[ia].mu.v[0]<<" "<<c->atom[ia].mu.v[1]<<" "<<c->atom[ia].mu.v[2]<<" ";
			//mdsave->muSave.sbuf<<c->atom[ia].F.v[0]<<" "<<c->atom[ia].F.v[1]<<" "<<c->atom[ia].F.v[2]<<" ";
			if (::bDipole) {
				mdsave->muSave.sbuf<<c->atom[ia].ind_mu.v[0]<<" "<<c->atom[ia].ind_mu.v[1]<<" "<<c->atom[ia].ind_mu.v[2]<<" ";
				mdsave->muSave.sbuf<<c->atom[ia].intr_mu.v[0]<<" "<<c->atom[ia].intr_mu.v[1]<<" "<<c->atom[ia].intr_mu.v[2];
			}
			mdsave->muSave.sbuf<<endl;
		}
		//mdsave->muSave.sbuf<<c->dyn->fc.v[0]<<" "<<c->dyn->fc.v[1]<<" "<<c->dyn->fc.v[2]<<" ";
		//mdsave->muSave.sbuf<<c->dyn->fc.v[3]<<" "<<c->dyn->fc.v[4]<<" "<<c->dyn->fc.v[5]<<" "<<endl;
	
		mdsave->muSave.sbuf<<endl;

		mdsave->muSave.flush_buf();
	}
	else {
		if (save_md_title) *(mdsave->muSave.out)<<"[MD] "<<mdsave->muSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0)
			*(mdsave->muSave.out)<<additional_title<<endl;
		
		for (ia = 0; ia < c->nAtoms; ia++) {
			*(mdsave->muSave.out)<<c->atom[ia].mu.v[0]<<" "<<c->atom[ia].mu.v[1]<<" "<<c->atom[ia].mu.v[2]<<" ";
			//*(mdsave->muSave.out)<<c->atom[ia].F.v[0]<<" "<<c->atom[ia].F.v[1]<<" "<<c->atom[ia].F.v[2]<<" ";
			if (::bDipole) {
				*(mdsave->muSave.out)<<c->atom[ia].ind_mu.v[0]<<" "<<c->atom[ia].ind_mu.v[1]<<" "<<c->atom[ia].ind_mu.v[2]<<" ";
				*(mdsave->muSave.out)<<c->atom[ia].intr_mu.v[0]<<" "<<c->atom[ia].intr_mu.v[1]<<" "<<c->atom[ia].intr_mu.v[2];
			}
			*(mdsave->muSave.out)<<endl;
		}
		//*(mdsave->muSave.out)<<c->dyn->fc.v[0]<<" "<<c->dyn->fc.v[1]<<" "<<c->dyn->fc.v[2]<<" ";
		//*(mdsave->muSave.out)<<c->dyn->fc.v[3]<<" "<<c->dyn->fc.v[4]<<" "<<c->dyn->fc.v[5]<<" "<<endl;
	
		*(mdsave->muSave.out)<<endl;
	}

	if (check_status) mdsave->muSave.save_over();
	return true;
	*/
}

void Speed_Verlet(LFVERLET_PM_MD &mdpar, PMOLECULE &m, bool reset_mol_speed) {
	int i = 0;

	for (i = 0; i < 3; i++) {
		mdpar.v0.v[i] = m.v.v[i];
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + m.alpha.v[i] * mdpar.dt;
		mdpar.v0.v[i] = 0.5 * (mdpar.v0_mh.v[i] + mdpar.v0_ph.v[i]);
		if (reset_mol_speed) m.v.v[i] = mdpar.v0.v[i];
	}
	m.updateSpatialMomentum();
}

void Speed_Kinetic(LFVERLET_PM_MD &mdpar, PMOLECULE &m) {
	int i = 0;
	double dv = 0;

	for (i = 0; i < 3; i++) {
		mdpar.v0.v[i] = m.v.v[i];
		dv = m.alpha.v[i] * mdpar.dt * 0.5;
		mdpar.v0_mh.v[i] = mdpar.v0.v[i] - dv;
		mdpar.v0_ph.v[i] = mdpar.v0.v[i] + dv;
	}
}

void Verlet(LFVERLET_PM_MD &mdpar, PMOLECULE &m) {
	if (m.bFixed) {
		memset(mdpar.v0_ph.v, 0, SIZE_V3);
		memset(mdpar.v0.v, 0, SIZE_V3);
		memset(m.v.v, 0, SIZE_V3);
		return;
	}

	int i = 0;
	double delta = 0;
	VECTOR3 dr;

	for (i = 0; i < 3; i++) {
		mdpar.v0_ph.v[i] = mdpar.v0_mh.v[i] + m.alpha.v[i] * mdpar.dt;
		mdpar.v0.v[i] = 1.5 * mdpar.v0_ph.v[i] - 0.5 * mdpar.v0_mh.v[i]; 
		m.v.v[i] = mdpar.v0.v[i];
		dr.v[i] = mdpar.v0_ph.v[i] * mdpar.dt;
	}
	m.shiftMol(dr);
	m.updateSpatialMomentum();
	
	for (i = 0; i < 3; i++) mdpar.v0_mh.v[i] = mdpar.v0_ph.v[i]; 
}

bool LFVERLET_PM_MD::save_mol_dyn_pars(PMOLECULE &m, long nloop, bool check_status, char* additional_title, bool save_md_title, bool bAsync) {
	if (pmsave.mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!pmsave.mdsave->mDynSave.one_more()) return false;
	}
	if (!pmsave.mdsave->mDynSave.bsave()) return true; // this loop will be not saved

	if (bAsync) {
		pmsave.mdsave->init_buf();
		if (save_md_title) pmsave.mdsave->mDynSave.sbuf<<"[MD] "<<pmsave.mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			pmsave.mdsave->mDynSave.sbuf<<additional_title<<endl;

		//pmsave.mdsave->mDynSave.sbuf<<this->ksi_Hoover<<" "<<vksi_Hoover<<endl;
		if (hdyn != NULL) {
			pmsave.mdsave->mDynSave.sbuf<<this->hdyn->nhc->ksi.v[0]<<" "<<this->hdyn->nhc->vksi.v[0]<<endl;
		}
		pmsave.mdsave->mDynSave.sbuf<<endl;

		pmsave.mdsave->mDynSave.flush_buf();
	}
	else {
		if (save_md_title) *(pmsave.mdsave->mDynSave.out)<<"[MD] "<<pmsave.mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
		if (additional_title != NULL && strlen(additional_title) > 0) 
			*(pmsave.mdsave->mDynSave.out)<<additional_title<<endl;

		//*(pmsave.mdsave->mDynSave.out)<<this->ksi_Hoover<<" "<<vksi_Hoover<<endl;
		if (hdyn != NULL) {
			*(pmsave.mdsave->mDynSave.out)<<this->hdyn->nhc->ksi.v[0]<<" "<<this->hdyn->nhc->vksi.v[0]<<endl;
		}
		*(pmsave.mdsave->mDynSave.out)<<endl;
	}

	if (check_status) pmsave.mdsave->mDynSave.save_over();
	return true;
}
