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

#include "MM.h"
#include "Mol.h"
#include "MD.h"
#include "var.h"

#include "cg-mm.h"
#include "cg-md.h"

namespace _coarse_grain_ {

bool save_xyz_struct(CG_MMOLECULE &mm, char *fname) {
	ofstream out;
	out.open(fname);
	if (!out.is_open()) {sprintf(errmsg, "can not open file %s to save macromolecular structure"); show_msg(errmsg); strcpy(errmsg, "\0"); return false;}
	mm.xyz_save(out);
	out.close();
	return true;
}

bool CGMM_MD_SAVE::save_kinetic_pars(CG_MMOLECULE &cgmm, long nloop, bool check_status, char* additional_title, bool save_md_title) {
	if (mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!mdsave->mKinSave.one_more()) return false;
	}
	if (!mdsave->mKinSave.bsave()) return true; // this loop will be not saved

	if (save_md_title) *(mdsave->mKinSave.out)<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
	if (additional_title != NULL && strlen(additional_title) > 0)
		*(mdsave->mKinSave.out)<<additional_title<<endl;
	int nc = 0;
	for (nc = 0; nc < cgmm.nCluster; nc++) {
		*(mdsave->mKinSave.out)<<nc<<" "<<cgmm.cluster[nc].r->v[0]<<" "<<cgmm.cluster[nc].r->v[1]<<" "<<cgmm.cluster[nc].r->v[2];
		//*(mdsave->mKinSave.out)<<" "<<cgmm.cluster[nc].bkm->v.v[0]<<" "<<cgmm.cluster[nc].bkm->v.v[1]<<" "<<cgmm.cluster[nc].bkm->v.v[2];
		*(mdsave->mKinSave.out)<<endl;
	}
	*(mdsave->mKinSave.out)<<endl;

	if (check_status) mdsave->mKinSave.save_over();
	return true;
}

void Speed_Verlet(LFVERLET_CGMM_MD &mdpar, CG_MMOLECULE &mm, bool reset_mol_speed) {
	//suppose mdpar.nCluster == mm.nCluster
	int nc = 0, i;
	CG_CLUSTER *pc = NULL;

	// check velocity of clusters except base 
	for (nc = 0; nc < mm.nCluster; nc++) {
		// Leapfrog Verlet algorithm
		pc = mm.cluster + nc;
		for (i = 0; i < 3; i++) {
			mdpar.v_ph[nc].v[i] = mdpar.v_mh[nc].v[i] + pc->bkm->alpha.v[i] * mdpar.dt;
			mdpar.v[nc].v[i] = 0.5 * (mdpar.v_mh[nc].v[i] + mdpar.v_ph[nc].v[i]);
			if (reset_mol_speed) pc->bkm->v.v[i] = mdpar.v[nc].v[i];
		}
	}
}

void Speed_Kinetic(LFVERLET_CGMM_MD &mdpar, CG_MMOLECULE &mm) {
	//suppose mdpar.nCluster == mm.nCluster
	int nc = 0, i;
	CG_CLUSTER *pc = NULL;

	for (nc = 0; nc < mm.nCluster; nc++) {
		// Leapfrog Verlet algorithm
		pc = mm.cluster + nc;
		for (i = 0; i < 3; i++) {
			mdpar.v[nc].v[i] = pc->bkm->v.v[i];
			mdpar.v_mh[nc].v[i] = mdpar.v[nc].v[i] - pc->bkm->alpha.v[i] * mdpar.dt * 0.5;
			mdpar.v_ph[nc].v[i] = mdpar.v[nc].v[i] + pc->bkm->alpha.v[i] * mdpar.dt * 0.5;
		}
	}
}

void Verlet(LFVERLET_CGMM_MD &mdpar, CG_MMOLECULE &mm) {
	//suppose mdpar.nCluster == mm.nCluster
	int nc = 0, i = 0;
	CG_CLUSTER *pc = NULL;
	VECTOR3 dr;

	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		// Leapfrog Verlet algorithm
		for (i = 0; i < 3; i++) {
			dr.v[i] = mdpar.v_ph[nc].v[i] * mdpar.dt;
			pc->r->v[i] += dr.v[i];
			// guess new position's velocity
			pc->bkm->v.v[i] = 1.5 * mdpar.v_ph[nc].v[i] - 0.5 * mdpar.v_mh[nc].v[i];
			mdpar.v_mh[nc].v[i] = mdpar.v_ph[nc].v[i];
			mdpar.v[nc].v[i] = pc->bkm->v.v[i];
		}
		if (bImplicitSolvent) {
			V3plusV3(pc->rc_solv, dr, pc->rc_solv) // geometrical center of cluster, in real cordinate
		}
	}
}

bool LFVERLET_CGMM_MD::save_mol_dyn_pars(CG_MMOLECULE &mm, long nloop, bool check_status, char* additional_title, bool save_md_title) {
	if (mmsave.mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}
	if (check_status) {
		if (!mmsave.mdsave->mDynSave.one_more()) return false;
	}
	if (!mmsave.mdsave->mDynSave.bsave()) return true; // this loop will be not saved

	if (save_md_title) *(mmsave.mdsave->mDynSave.out)<<"[MD] "<<mmsave.mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
	if (additional_title != NULL && strlen(additional_title) > 0) *(mmsave.mdsave->mDynSave.out)<<additional_title<<endl;

	int i = 0;
	/*
	for (int nc = 0; nc < mm.nCluster; nc++) {
		*(mmsave.mdsave->mDynSave.out)<<nc;
		//for (i = 0; i < 3; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<mm.cluster[nc].invNE->f.v[i];
		for (i = 0; i < 3; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<mm.cluster[nc].dyn->fc.v[i];
		*(mmsave.mdsave->mDynSave.out)<<endl;
	}
	*/

	//*(mmsave.mdsave->mDynSave.out)<<this->ksi_Hoover<<" "<<this->ksi_frame_trans<<" "<<this->ksi_frame_rot<<endl;
	if (hdyn != NULL) {
		*(mmsave.mdsave->mDynSave.out)<<hdyn->nhc->ksi.v[0]<<" "<<hdyn->nhc->vksi.v[0]<<endl;
	}
	//*(mmsave.mdsave->mDynSave.out)<<mm.mKD->Ekt<<"  "<<mm.mKD->Ekr<<endl;
	*(mmsave.mdsave->mDynSave.out)<<endl;

	if (check_status) mmsave.mdsave->mDynSave.save_over();
	return true;
}



/**********************************************************************************
						CGSM MD leapfrog Verlet algorithm
**********************************************************************************/

bool save_xyz_struct(CGSM &sm, char *fname) {
	ofstream out;
	out.open(fname);
	if (!out.is_open()) {sprintf(errmsg, "can not open file %s to save macromolecular structure"); show_msg(errmsg); strcpy(errmsg, "\0"); return false;}
	sm.xyz_save(out);
	out.close();
	return true;
}

bool CGSM_MD_SAVE::save_kinetic_pars(CGSM &cgsm, long nloop, bool check_status, char* additional_title, bool save_md_title) {
	if (mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}

	if (check_status) {
		if (!mdsave->mKinSave.one_more()) return false;
	}
	if (!mdsave->mKinSave.bsave()) return true; // this loop will be not saved

	if (save_md_title) *(mdsave->mKinSave.out)<<"[MD] "<<mdsave->mKinSave.save_indx()<<" "<<nloop<<endl;
	if (additional_title != NULL && strlen(additional_title) > 0)
		*(mdsave->mKinSave.out)<<additional_title<<endl;

	*(mdsave->mKinSave.out)<<cgsm.r->v[0]<<" "<<cgsm.r->v[1]<<" "<<cgsm.r->v[2];
	//*(mdsave->mKinSave.out)<<" "<<cgsm.bkm->v.v[0]<<" "<<cgsm.bkm->v.v[1]<<" "<<cgsm.bkm->v.v[2];
	*(mdsave->mKinSave.out)<<endl;

	*(mdsave->mKinSave.out)<<endl;

	if (check_status) mdsave->mKinSave.save_over();
	return true;
}

void Speed_Verlet(LFVERLET_CGSM_MD &mdpar, CGSM &sm, bool reset_mol_speed) {
	for (int i = 0; i < 3; i++) {
		mdpar.v_ph.v[i] = mdpar.v_mh.v[i] + sm.bkm->alpha.v[i] * mdpar.dt;
		mdpar.v.v[i] = 0.5 * (mdpar.v_mh.v[i] + mdpar.v_ph.v[i]);
		if (reset_mol_speed) sm.bkm->v.v[i] = mdpar.v.v[i];
	}
}

void Speed_Kinetic(LFVERLET_CGSM_MD &mdpar, CGSM &sm) {
	for (int i = 0; i < 3; i++) {
		mdpar.v.v[i] = sm.bkm->v.v[i];
		mdpar.v_mh.v[i] = mdpar.v.v[i] - sm.bkm->alpha.v[i] * mdpar.dt * 0.5;
		mdpar.v_ph.v[i] = mdpar.v.v[i] + sm.bkm->alpha.v[i] * mdpar.dt * 0.5;
	}
}

void Verlet(LFVERLET_CGSM_MD &mdpar, CGSM &sm) {
	int i;
	VECTOR3 dr;

	// Leapfrog Verlet algorithm
	for (i = 0; i < 3; i++) {
		dr.v[i] = mdpar.v_ph.v[i] * mdpar.dt;
		sm.r->v[i] += dr.v[i];
		// guess new position's velocity
		sm.bkm->v.v[i] = 1.5 * mdpar.v_ph.v[i] - 0.5 * mdpar.v_mh.v[i];
		mdpar.v_mh.v[i] = mdpar.v_ph.v[i];
		mdpar.v.v[i] = sm.bkm->v.v[i];
	}
	if (bImplicitSolvent) {
		V3plusV3(sm.rc_solv, dr, sm.rc_solv) // geometrical center of cluster, in real cordinate
	}
}

bool LFVERLET_CGSM_MD::save_mol_dyn_pars(CGSM &sm, long nloop, bool check_status, char* additional_title, bool save_md_title) {
	if (smsave.mdsave == NULL) {sprintf(errmsg, "no saving object"); return false;}
	if (check_status) {
		if (!smsave.mdsave->mDynSave.one_more()) return false;
	}
	if (!smsave.mdsave->mDynSave.bsave()) return true; // this loop will be not saved

	if (save_md_title) *(smsave.mdsave->mDynSave.out)<<"[MD] "<<smsave.mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
	if (additional_title != NULL && strlen(additional_title) > 0) *(smsave.mdsave->mDynSave.out)<<additional_title<<endl;

	int i = 0;
	/*
	//for (i = 0; i < 3; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<sm.invNE->f.v[i];
	for (i = 0; i < 3; i++) *(mmsave.mdsave->mDynSave.out)<<"  "<<sm.dyn->fc.v[i];
	*(mmsave.mdsave->mDynSave.out)<<endl;
	*/

	//*(smsave.mdsave->mDynSave.out)<<this->ksi_Hoover<<"  "<<vksi_Hoover<<endl;
	if (hdyn != NULL) {
		*(smsave.mdsave->mDynSave.out)<<hdyn->nhc->ksi.v[0]<<"  "<<hdyn->nhc->vksi.v[0]<<endl;
	}
	*(smsave.mdsave->mDynSave.out)<<endl;

	if (check_status) smsave.mdsave->mDynSave.save_over();
	return true;
}

} //end of namespace _coarse_grain_ 
