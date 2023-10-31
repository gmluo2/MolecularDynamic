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
#include "spme_interact.h"
extern BSplineFunc<2> bsp;

using namespace _cmm_3d_;

using namespace _EwaldSum_real_;

#include "interface.h"
using namespace _interface_constraint_;

#include "mdvar.h"

extern void ExplicitEInteract(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, EInteractVar& evar, InteractRes &res);

#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

void setup_mol_kin(MMOL_MD_CELL &mdcell) {
	int im, ic;
	MMOLECULE *mm;
	SMOLECULE *sm;
	for (im = 0; im < mdcell.mm.n; im++) {
		mm = mdcell.mm.m + im;
		for (ic = 0; ic < mm->nCluster; ic++) {
			if (mm->cluster[ic].bkm == NULL) mm->cluster[ic].bkm = new BASIC_KINEMATIC;
		}
	}
	for (im = 0; im < mdcell.sm.n; im++) {
		sm = mdcell.sm.m + im;
		if (sm->c->bkm == NULL) sm->c->bkm = new BASIC_KINEMATIC;
	}
	// pm has the velocity in PM, while not in pm->c
}

void MMOL_MD_CELL::check_eneutral() {
	eNeutral = true;
	int i, ic, nc = 0, na = 0;
	for (i = 0; i < mm.n; i++) { 
		for (ic = 0; ic < mm.m[i].nCluster; ic++) {
			mm.m[i].cluster[ic].check_eneutral();
			if (!mm.m[i].cluster[ic].eNeutral) {nc++; na += mm.m[i].cluster[ic].nAtoms;}
		}
	}
	for (i = 0; i < sm.n; i++) { 
		sm.m[i].c->check_eneutral();
		if (!sm.m[i].c->eNeutral) {nc++; na += sm.m[i].c->nAtoms;}
	}
	for (i = 0; i < pm.n; i++) { 
		pm.m[i].c->check_eneutral();
		if (!pm.m[i].c->eNeutral) {nc++; na += pm.m[i].c->nAtoms;}
	}
	ccluster.set_array(nc); catom.set_array(na);
	if (nc > 0) eNeutral = false;
	
	nc = 0;
	for (i = 0; i < mm.n; i++) { 
		for (ic = 0; ic < mm.m[i].nCluster; ic++) {
			if (!mm.m[i].cluster[ic].eNeutral) {ccluster.m[nc] = mm.m[i].cluster + ic; nc++;}
		}
	}
	for (i = 0; i < sm.n; i++) { 
		if (!sm.m[i].c->eNeutral) {ccluster.m[nc] = sm.m[i].c; nc++;}
	}
	for (i = 0; i < pm.n; i++) { 
		if (!pm.m[i].c->eNeutral) {ccluster.m[nc] = pm.m[i].c; nc++;}
	}

	na = 0;
	for (ic = 0; ic < ccluster.n; ic++) {
		for (i = 0; i < ccluster.m[ic]->nAtoms; i++) {
			catom.m[na] = ccluster.m[ic]->atom + i;
			na++;
		}
	}
}

void MMOL_MD_CELL::init_nhc_mols() { // group molecules for each HOOVER_DYN
	int nMM = 0, nSM = 0, nPM = 0;
	int imol = 0, ihdyn = 0, indx = 0;
	for (ihdyn = 0; ihdyn < hdyn.n; ihdyn++) {
		nMM = 0; nSM = 0; nPM = 0;
		for (imol = 0; imol < lfmd.n; imol++) {
			if (lfmd.m[imol].iHDyn->nhc == hdyn.m[ihdyn].nhc) nMM++;
		}
		for (imol = 0; imol < lfsmd.n; imol++) {
			if (lfsmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) nSM++;
		}
		for (imol = 0; imol < lfpmd.n; imol++) {
			if (lfpmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) nPM++;
		}
		hdyn.m[ihdyn].mmIndx.SetArray(nMM);
		if (nMM > 0) {
			indx = 0;
			for (imol = 0; imol < lfmd.n; imol++) {
				if (lfmd.m[imol].iHDyn->nhc == hdyn.m[ihdyn].nhc) {
					hdyn.m[ihdyn].mmIndx.m[indx] = imol; indx++;
				}
			}
		}
		hdyn.m[ihdyn].smIndx.SetArray(nSM);
		if (nSM > 0) {
			indx = 0;
			for (imol = 0; imol < lfsmd.n; imol++) {
				if (lfsmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) {
					hdyn.m[ihdyn].smIndx.m[indx] = imol; indx++;
				}
			}
		}
		hdyn.m[ihdyn].pmIndx.SetArray(nPM);
		if (nPM > 0) {
			indx = 0;
			for (imol = 0; imol < lfpmd.n; imol++) {
				if (lfpmd.m[imol].hdyn->nhc == hdyn.m[ihdyn].nhc) {
					hdyn.m[ihdyn].pmIndx.m[indx] = imol; indx++;
				}
			}
		}
	}
}

void MMOL_MD_CELL::check_freedom() {
	int NF = 0;
	int imol = 0, ihdyn = 0, i;
	MDCELL_HDYN *pdyn = NULL;
	this->NF = 0;
	for (ihdyn = 0; ihdyn < hdyn.n; ihdyn++) {
		pdyn = hdyn.m + ihdyn;
		NF = 0; 
		for (i = 0; i < pdyn->mmIndx.n; i++) {
			imol = pdyn->mmIndx.m[i];
			NF += (mm.m[imol].nCluster - 1);
			NF += 6;
		}
		// each sm has 6 freedomes, translation + rotation
		NF += pdyn->smIndx.n * 6;
		//for (i = 0; i < pdyn->smIndx.n; i++) {
		//	imol = pdyn->smIndx.m[i];
		//	NF += 6;
		//}
		// each sm has 3 freedomes, translation only
		NF += pdyn->pmIndx.n * 3;
		//for (i = 0; i < pm.n; i++) NF += 3; // translation only

		pdyn->nhc->set_freedom(NF);
		this->NF += NF;
	}

	this->Ethermal = this->NF * ::Ek0Internal;

	if (this->bHDynMC_T) {
		this->Ethermal_mct = 0;
		if (this->ix_mct) {this->Ethermal_mct += ::Ek0Internal; NF -= 1;}
		if (this->iy_mct) {this->Ethermal_mct += ::Ek0Internal; NF -= 1;}
		if (this->iz_mct) {this->Ethermal_mct += ::Ek0Internal; NF -= 1;}
	}

	char msg[256 ] = "\0";
	sprintf(msg, "Total freedom of the cell : %d, thermal energy : %f kT", this->NF, this->Ethermal); show_log(msg, true);
}

bool MMOL_MD_CELL::MDSave(long nloop, float Ep, bool bEachMolDyn, bool bAsync) {
	int nsave = 0, n = 0;
	MD_SAVE *mdsave = NULL;
	if (mm.n > 0) mdsave = lfmd.m[0].mmsave.mdsave;
	else if (sm.n > 0) mdsave = lfsmd.m[0].smsave.mdsave;
	else if (pm.n > 0) mdsave = lfpmd.m[0].pmsave.mdsave;
	else return false;
	nsave = mdsave->nsave;

	char title[20] = "\0";
	bool save_md_title = true;
	bool save_status = false;

	if (!mdsave->mKinSave.one_more()) return false;
	if (!mdsave->muSave.one_more()) return false;
	if (!mdsave->mESave.one_more()) return false;
	if (!mdsave->mDynSave.one_more()) return false;

	if (bAsync) {
		msave->bBuffered = true;
		mdsave->init_buf();
	}

	if (mdsave->mKinSave.bsave()) {
		save_md_title = true;
		for (n = 0; n < pm.n; n++) {
			sprintf(title, "PM %d", n);
			lfpmd.m[n].pmsave.save_kinetic_pars(pm.m[n], nloop, false, title, save_md_title, bAsync);
			lfpmd.m[n].pmsave.save_atomic_dipole(pm.m[n], nloop, false, title, save_md_title, bAsync);
			lfpmd.m[n].pmsave.mdsave->save_energy(pm.m[n].Ek, nloop, false, NULL, save_md_title, bAsync);
			if (bEachMolDyn) lfpmd.m[n].save_mol_dyn_pars(pm.m[n], nloop, false, title, save_md_title, bAsync);
			save_md_title = false;
		}
		for (n = 0; n < sm.n; n++) {
			sprintf(title, "SM %d", n);
			lfsmd.m[n].smsave.save_kinetic_pars(sm.m[n], nloop, false, title, save_md_title, bAsync);
			lfsmd.m[n].smsave.save_atomic_dipole(sm.m[n], nloop, false, title, save_md_title, bAsync);
			lfsmd.m[n].smsave.mdsave->save_energy(sm.m[n].Ek, nloop, false, NULL, save_md_title, bAsync);
			if (bEachMolDyn) lfsmd.m[n].save_mol_dyn_pars(sm.m[n], nloop, false, title, save_md_title, bAsync);
			save_md_title = false;
		}
		for (n = 0; n < mm.n; n++) {
			sprintf(title, "MM %d", n);
			lfmd.m[n].mmsave.save_kinetic_pars(mm.m[n], nloop, false, title, save_md_title, bAsync);
			lfmd.m[n].mmsave.save_atomic_dipole(mm.m[n], nloop, false, title, save_md_title, bAsync);
			//lfmd.m[n].mmsave.mdsave->save_energy(mm.m[n].Ek, nloop, false, NULL, save_md_title);
			lfmd.m[n].mmsave.save_energy(mm.m[n], nloop, false, NULL, save_md_title, bAsync);
			if (bEachMolDyn) lfmd.m[n].save_mol_dyn_pars(mm.m[n], nloop, false, title, save_md_title, bAsync);
			save_md_title = false;
		}
		mdsave->save_energy(Ep, nloop, false, NULL, false, bAsync);

		if (!bEachMolDyn && mdsave->mDynSave.bsave()) {
			if (bAsync) {
				mdsave->init_buf();
				mdsave->mDynSave.sbuf<<"[MD] "<<mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
				//mdsave->mDynSave.sbuf<<hdyn.ksi_Hoover<<" "<<hdyn.vksi_Hoover<<endl<<endl;
				for (n = 0; n < hdyn.n; n++) {
					if (hdyn.m[n].nhc != NULL && hdyn.m[n].mol_coupled()) {
						mdsave->mDynSave.sbuf<<n<<":  "<<hdyn.m[n].nhc->ksi.v[0]<<" "<<hdyn.m[n].nhc->vksi.v[0]<<endl;
					}
				}
				if (bHDynMC_T) {
					mdsave->mDynSave.sbuf<<"mc trans:  "<<mcHDyn.nhc->ksi.v[0]<<" "<<mcHDyn.nhc->vksi.v[0]<<"  "<<mcHDyn.nhc->dT<<endl;
				}
				mdsave->mDynSave.sbuf<<endl;

				mdsave->mDynSave.flush_buf();
			}
			else {
				*(mdsave->mDynSave.out)<<"[MD] "<<mdsave->mDynSave.save_indx()<<" "<<nloop<<endl;
				//*(mdsave->mDynSave.out)<<hdyn.ksi_Hoover<<" "<<hdyn.vksi_Hoover<<endl<<endl;
				for (n = 0; n < hdyn.n; n++) {
					if (hdyn.m[n].nhc != NULL && hdyn.m[n].mol_coupled()) {
						*(mdsave->mDynSave.out)<<n<<":  "<<hdyn.m[n].nhc->ksi.v[0]<<" "<<hdyn.m[n].nhc->vksi.v[0]<<endl;
					}
				}
				if (bHDynMC_T) {
					*(mdsave->mDynSave.out)<<"mc trans:  "<<mcHDyn.nhc->ksi.v[0]<<" "<<mcHDyn.nhc->vksi.v[0]<<"  "<<mcHDyn.nhc->dT<<endl;
				}
				*(mdsave->mDynSave.out)<<endl;
			}
		}
		mdsave->mKinSave.save_over();
		mdsave->muSave.save_over();
		mdsave->mESave.save_over();
		mdsave->mDynSave.save_over();
	}
	return true;
}


#if _SYS_ == _WINDOWS_SYS_
int AsyncMDSave_thread(LPVOID *vp) {
#elif _SYS_ == _LINUX_SYS_
void* AsyncMDSave_thread(void *vp) {
#endif
	AsyncMDSave_VARS *par = (AsyncMDSave_VARS*)vp;
	MD_SAVE *mdsave = par->mdcell->msave;
	MMOL_MD_CELL *mdcell = par->mdcell;	

	mdsave->flush_buf();
	
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void AsyncMDSave(MMOL_MD_CELL &md_cell, AsyncMDSave_VARS &av) {
	if (!md_cell.msave->bBuffered) return;
	av.set_pars(0);
	av.mdcell = &md_cell;
#if _SYS_ == _WINDOWS_SYS_
	ResetEvent(av.hMMThread);
	av.thread = AfxBeginThread((AFX_THREADPROC)AsyncMDSave_thread, (LPVOID)(&av), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
#elif _SYS_ == _LINUX_SYS_
	pthread_create(&(av.thread), NULL, &AsyncMDSave_thread, (void *)(&av));
#endif
}


void MMOL_MD_CELL::init_mdcell_velocity_memory() {
	int imol;
	mvb.mmv1.SetArray(mm.n); mvb.mmv2.SetArray(mm.n); 
	mvb.mm_dvmax.SetArray(mm.n); mvb.mm_dwmax.SetArray(mm.n); 
	mvb.mm_status.SetArray(mm.n); mvb.mm_sc_status.SetArray(mm.n);

	for (imol = 0; imol < mm.n; imol++) {
		mvb.mmv1.m[imol].cv.SetArray(mm.m[imol].nCluster);
		mvb.mmv2.m[imol].cv.SetArray(mm.m[imol].nCluster);
	}

	mvb.sm_dvmax.SetArray(sm.n); mvb.sm_dwmax.SetArray(sm.n);
	mvb.sm_status.SetArray(sm.n); mvb.sm_sc_status.SetArray(sm.n);

	mvb.pm_dvmax.SetArray(pm.n); 
	mvb.pm_status.SetArray(pm.n); mvb.pm_sc_status.SetArray(pm.n);
}

void MMOL_MD_CELL::init_polarization_buff() {
	int imol, iatom, natom = 0, ic;
	MMOLECULE *pmm = NULL;
	ind_mu_mm.SetArray(mm.n);
	for (imol = 0; imol < mm.n; imol++) {
		pmm = mm.m + imol;
		ind_mu_mm.m[imol].m = pmm;
		natom = 0;
		for (ic = 0; ic < pmm->nCluster; ic++) natom += pmm->cluster[ic].nAtoms;
		ind_mu_mm.m[imol].mu.SetArray(natom);
	}
	ind_mu_sm.SetArray(sm.n);
	for (imol = 0; imol < sm.n; imol++) {
		ind_mu_sm.m[imol].m = sm.m + imol;
		ind_mu_sm.m[imol].mu.SetArray(sm.m[imol].c->nAtoms);
	}
	ind_mu_pm.SetArray(pm.n);
	for (imol = 0; imol < pm.n; imol++) {
		ind_mu_pm.m[imol].m = pm.m + imol;
		ind_mu_pm.m[imol].mu.SetArray(pm.m[imol].c->nAtoms);
	}
}

void MM_maxdiff_backup_polarization(MOL_POL<MMOLECULE> *mp, int nmol, double &dmu_max) {
	int imol, iatom, ic, imu;
	MMOLECULE *pmm = NULL;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 *mu_buf = NULL;
	VECTOR3 *mu1 = NULL, *mu2 = NULL;
	VECTOR3 dmu;
	double vt;
	dmu_max = 0;
	for (imol = 0; imol < nmol; imol++) {
		pmm = mp[imol].m; mu_buf = mp[imol].mu.m; imu = 0;
		for (ic = 0; ic < pmm->nCluster; ic++) {
			pc = pmm->cluster + ic;
			if (pc->eNeutral) continue;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom; 
				if (patom->mu_hist == NULL) continue;
				mu1 = &(patom->ind_mu); mu2 = mu_buf + imu;
				V3minusV3((*mu1), (*mu2), dmu)
				vt = FABS(dmu.v[0]); dmu_max = (dmu_max > vt ? dmu_max : vt);
				vt = FABS(dmu.v[1]); dmu_max = (dmu_max > vt ? dmu_max : vt);
				vt = FABS(dmu.v[2]); dmu_max = (dmu_max > vt ? dmu_max : vt);
				V32V3((*mu1), (*mu2))
				imu++;
			}
		}
	}
}

void SM_maxdiff_backup_polarization(MOL_POL<SMOLECULE> *mp, int nmol, double &dmu_max) {
	int imol, iatom;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 *mu_buf = NULL;
	VECTOR3 *mu1 = NULL, *mu2 = NULL;
	VECTOR3 dmu;
	double vt;
	dmu_max = 0;
	for (imol = 0; imol < nmol; imol++) {
		pc = mp[imol].m->c; mu_buf = mp[imol].mu.m;
		if (pc->eNeutral) continue;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom; 
			if (patom->mu_hist == NULL) continue;
			mu1 = &(patom->ind_mu); mu2 = mu_buf + iatom;
			V3minusV3((*mu1), (*mu2), dmu)
			vt = FABS(dmu.v[0]); dmu_max = (dmu_max > vt ? dmu_max : vt);
			vt = FABS(dmu.v[1]); dmu_max = (dmu_max > vt ? dmu_max : vt);
			vt = FABS(dmu.v[2]); dmu_max = (dmu_max > vt ? dmu_max : vt);
			V32V3((*mu1), (*mu2))
		}
	}
}

void PM_maxdiff_backup_polarization(MOL_POL<PMOLECULE> *mp, int nmol, double &dmu_max) {
	int imol, iatom;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 *mu_buf = NULL;
	VECTOR3 *mu1 = NULL, *mu2 = NULL;
	VECTOR3 dmu;
	double vt;
	dmu_max = 0;
	for (imol = 0; imol < nmol; imol++) {
		pc = mp[imol].m->c; mu_buf = mp[imol].mu.m;
		if (pc->eNeutral) continue;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom; 
			if (patom->mu_hist == NULL) continue;
			mu1 = &(patom->ind_mu); mu2 = mu_buf + iatom;
			V3minusV3((*mu1), (*mu2), dmu)
			vt = FABS(dmu.v[0]); dmu_max = (dmu_max > vt ? dmu_max : vt);
			vt = FABS(dmu.v[1]); dmu_max = (dmu_max > vt ? dmu_max : vt);
			vt = FABS(dmu.v[2]); dmu_max = (dmu_max > vt ? dmu_max : vt);
			V32V3((*mu1), (*mu2))
		}
	}
}

void MMOL_MD_CELL::set_atom_multipole(int mMP) {
	int imol, iatom, ic, imu;
	MMOLECULE *pmm = NULL;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	for (imol = 0; imol < mm.n; imol++) {
		pmm = mm.m + imol;
		for (ic = 0; ic < pmm->nCluster; ic++) {
			pc = pmm->cluster + ic;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom; 
				if (patom->par->alpha < 1e-4) {
					if (patom->bIMu) patom->mMP = mMP;
					else patom->mMP = 0;
					if (patom->eNeutral) patom->eNeutral = false;
				}
				else if (patom->c != 0) {patom->mMP = mMP; patom->eNeutral = false;}
				else {patom->mMP = 0; patom->eNeutral = true;}
			}
			pc->check_eneutral();
		}
	}
	for (imol = 0; imol < sm.n; imol++) {
		pc = sm.m[imol].c;
		for (iatom = 0; iatom < sm.m[imol].c->nAtoms; iatom++) {
			patom = pc->atom + iatom; 
			if (patom->par->alpha < 1e-4) {
				if (patom->bIMu) patom->mMP = mMP;
				else patom->mMP = 0;
				if (patom->eNeutral) patom->eNeutral = false;
			}
			else if (patom->c != 0) {patom->mMP = mMP; patom->eNeutral = false;}
			else {patom->mMP = 0; patom->eNeutral = true;}
		}
		pc->check_eneutral();
	}
	for (imol = 0; imol < pm.n; imol++) {
		pc = pm.m[imol].c;
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom; 
			if (patom->par->alpha < 1e-4) {
				if (patom->bIMu) patom->mMP = mMP;
				else patom->mMP = 0;
				if (patom->eNeutral) patom->eNeutral = false;
			}
			else if (patom->c != 0) {patom->mMP = mMP; patom->eNeutral = false;}
			else {patom->mMP = 0; patom->eNeutral = true;}
		}
		pc->check_eneutral();
	}
}

void MDCELL_mass_center(MMOL_MD_CELL &mmcell) {
	VECTOR3 mc;
	float M = 0;
	MMOLECULE *m = NULL;
	int nm = 0, i = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		m->calMassCenter(true);
		for (i = 0; i < 3; i++) mc.v[i] += m->mDynKin.M * m->mDynKin.cm.v[i];
		M += m->mDynKin.M;
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 3; i++) mc.v[i] += sm->c->M * sm->r.v[i];
		M += sm->c->M;
	}
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m + nm;
		for (i = 0; i < 3; i++) mc.v[i] += pm->c->M * pm->r.v[i];
		M += pm->c->M;
	}
	for (i = 0; i < 3; i++) mc.v[i] /= M;
	V32V3(mc, mmcell.rc)
}

void set_free_cell_mass_center(MMOL_MD_CELL &mmcell, VECTOR3 &O) {
	VECTOR3 mc;
	float M = 0;
	MMOLECULE *m = NULL;
	int nm = 0, i = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		m->calMassCenter(true);
		for (i = 0; i < 3; i++) {
			mc.v[i] += m->mDynKin.M * m->mDynKin.cm.v[i];
			M += m->mDynKin.M;
		}
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 3; i++) {
			mc.v[i] += sm->c->M * sm->r.v[i];
			M += sm->c->M;
		}
	}
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m + nm;
		for (i = 0; i < 3; i++) {
			mc.v[i] += pm->c->M * pm->r.v[i];
			M += pm->c->M;
		}
	}
	for (i = 0; i < 3; i++) mc.v[i] /= M;
	
	VECTOR3 ds;
	for (i = 0; i < 3; i++) ds.v[i] = O.v[i] - mc.v[i];
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		m->shiftMol(ds);
	}
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		sm->shiftMol(ds);
	}
	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m + nm;
		pm->shiftMol(ds);
	}
}

void cluster_check_periodic_cell(MMOL_MD_CELL &mmcell) {
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 3; i++) {
				xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
				x = (float)(pc->cm.v[i]);
				PERIOD_CHECK(x, xl, xr, period, pc->dr_cmm.v[i])
				if (pc->dr_cmm.v[i] != 0) pc->bCMMShift = true;
			}
			pc->update_vp();
		}
	}
	/* not necessary, because sm has only one cluster
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(sm->r0.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
	}
	*/
}

void molecule_check_periodic_cell(MMOL_MD_CELL &mmcell) {
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		m = mmcell.mm.m + nm;
		m->calMassCenter(false);
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(m->mDynKin.cm.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) m->shiftMol(ds);

		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 3; i++) {
				xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
				x = (float)(pc->cm.v[i]);
				PERIOD_CHECK(x, xl, xr, period, pc->dr_cmm.v[i])
				if (pc->dr_cmm.v[i] != 0) pc->bCMMShift = true;
			}
			pc->update_vp();
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(sm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m + nm;
		for (i = 0; i < 3; i++) {
			xl = -mmcell.h[i]; xr = mmcell.h[i]; period = xr - xl;
			x = (float)(pm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) pm->shiftMol(ds);
	}

	//cluster_check_periodic_cell(mmcell);
}

void MDCELL_molecule_PeriodicCell_check(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	MMOLECULE *m = NULL;
	CLUSTER *pc = NULL;
	int nm = 0, i = 0, nc = 0, imol = 0;
	VECTOR3 ds;
	float xl = 0, xr = 0, period = 0, delta = 0, x = 0;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		m = mdcell->mm.m + imol;
		m->calMassCenter(false);
		for (i = 0; i < 3; i++) {
			xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
			x = (float)(m->mDynKin.cm.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) m->shiftMol(ds);

		for (nc = 0; nc < m->nCluster; nc++) {
			pc = m->cluster + nc;
			pc->bCMMShift = false;
			for (i = 0; i < 3; i++) {
				xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
				//x = (float)(pc->cm.v[i]);
				x = (float)(pc->vp->v[i]);
				PERIOD_CHECK(x, xl, xr, period, pc->dr_cmm.v[i])
				if (pc->dr_cmm.v[i] != 0) pc->bCMMShift = true;
			}
			pc->update_vp();
		}
	}

	char msg[256] = "\0";
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; 
		for (i = 0; i < 3; i++) {
			xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
			x = (float)(sm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) sm->shiftMol(ds);
		/*
		if (sm->r.v[0] > mdcell->h[0] || sm->r.v[0] < -mdcell->h[0] ||
			sm->r.v[1] > mdcell->h[1] || sm->r.v[1] < -mdcell->h[1] ||
			sm->r.v[2] > mdcell->h[2] || sm->r.v[2] < -mdcell->h[2]){
				sprintf(msg, "strange: sm %d has coordinate out of cell [%f, %f, %f]", sm->r.v[0], sm->r.v[1], sm->r.v[2]);
				show_log(msg, true);
		}
		*/
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		for (i = 0; i < 3; i++) {
			xl = -mdcell->h[i]; xr = mdcell->h[i]; period = xr - xl;
			x = (float)(pm->r.v[i]);
			PERIOD_CHECK(x, xl, xr, period, ds.v[i])
		}
		if (ds.v[0] != 0 || ds.v[1] != 0 || ds.v[2] != 0) pm->shiftMol(ds);
	}
}

// **********************************************************************************************************************
// calculating the Hoover Force on each cluster
// **********************************************************************************************************************
void cal_Hoover_force(MMOLECULE *mm, double ksi, double ksi_mc, VECTOR3 &vmc) {
	int ic;
	CLUSTER *pc;
	MATRIX<6> *M, *phai;
	MATRIX<3> *I;
	SVECTOR<double, 6> f, vh;
	SVECTOR<double, 3> fmc, torq;
	SVECTOR<double, 3> dr;
	double *v1, *v2;
	float mass;
	double *fh;

	int method = ::MM_SpeedVerlet_method;

	NEIMO_CalHooverForce(*mm, ksi, mm->mDynKin.f_Hoover);

	for (ic = 0; ic < mm->nCluster; ic++) {
		pc = mm->cluster + ic;
		mass = pc->M; 
		M = &(pc->invNE->M);
	
		V6zero(f)

	if (ksi_mc != 0) {
		fmc.v[0] = -mass * ksi_mc * vmc.v[0];
		fmc.v[1] = -mass * ksi_mc * vmc.v[1];
		fmc.v[2] = -mass * ksi_mc * vmc.v[2];
		VECT3((*(pc->Op)), pc->cm, dr)
		V3PV3(dr, fmc, torq)

		// accumulate Hoover force
		//fh = pc->dyn->f_Hoover.v;
		f.v[0] += torq.v[0];
		f.v[1] += torq.v[1];
		f.v[2] += torq.v[2];
		f.v[3] += fmc.v[0];
		f.v[4] += fmc.v[1];
		f.v[5] += fmc.v[2];

		V6plusV6(pc->dyn->f_Hoover, f, pc->dyn->f_Hoover)

		phai = &(pc->invNE->phai_cm2me);
		MpV6((*phai), f, vh) // vh = phai(cm=>pc) * f
		v1 = mm->mDynKin.f_Hoover.v; v2 = vh.v; 
		V6plus(v1, v2, v1) // mDynKin.f_Hoover += vh
	}
	}
}

void cal_Hoover_force(SMOLECULE *sm, double ksi, double ksi_mc, VECTOR3 &vmc) {
	VECTOR<6> *V = &(sm->c->bkm->V0);
	VECTOR3 w;
	VECTOR3 force, torque;
	MATRIX<3> *I = &(sm->I);
	float mass = sm->c->M;
	double *v1, *v2;
	double *fh = sm->c->dyn->f_Hoover.v;

	memset(fh, 0, SIZE_V6);

	if (ksi != 0) {
		// we want to compensate rotation from the rotation of the cell, but keep the translation part
		w.v[0] = -ksi * V->v[0];
		w.v[1] = -ksi * V->v[1];
		w.v[2] = -ksi * V->v[2];
		MpV3((*I), w, torque)
	}

	force.v[0] = -mass * (ksi_mc * vmc.v[0] + ksi * V->v[3]);
	force.v[1] = -mass * (ksi_mc * vmc.v[1] + ksi * V->v[4]);
	force.v[2] = -mass * (ksi_mc * vmc.v[2] + ksi * V->v[5]);
	
	// total virtual force
	fh[0] = torque.v[0];
	fh[1] = torque.v[1];
	fh[2] = torque.v[2];
	fh[3] = force.v[0];
	fh[4] = force.v[1];
	fh[5] = force.v[2];
}

void cal_Hoover_force(PMOLECULE *pm, double ksi, double ksi_mc, VECTOR3 &vmc) {
	double *v = pm->v.v;
	double *v1, *v2;
	double *fh = pm->c->dyn->f_Hoover.v;
	float mass = pm->c->M;

	memset(fh, 0, SIZE_V6);

	// total virtual force
	fh[3] = -mass *(ksi_mc * vmc.v[0] + ksi * v[0]);
	fh[4] = -mass *(ksi_mc * vmc.v[1] + ksi * v[1]);
	fh[5] = -mass *(ksi_mc * vmc.v[2] + ksi * v[2]);
}

extern void ScaleVelocity(MMOLECULE &mm, double c);
extern void InitVelocity(MMOLECULE &mm);
extern void ScaleVelocity(SMOLECULE &mm, double c);
extern void ScaleVelocity(PMOLECULE &mm, double c);

void MDCELL_SpeedVerlet(MDCELL_ThreadJob<MMOL_MD_CELL> *job, bool *status) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	VECTOR3 f, torq;
	double ksi_mc = 0;
	if (mdcell->bHDynMC_T) ksi_mc = mdcell->mcHDyn.nhc->ksi.v[0];

	VECTOR3 vmc, wmc, Vw;
	V32V3(mdcell->Vmc, vmc);
	if (mdcell->ix_mct == 0) vmc.v[0] = 0;
	if (mdcell->iy_mct == 0) vmc.v[1] = 0;
	if (mdcell->iz_mct == 0) vmc.v[2] = 0;
	int nm = 0, imol = 0, ihdyn = 0, ic;

	double dt =  mdcell->hdyn.m[0].nhc->dt, ksi;
	double uConvg = 0.005 / dt, wConvg = 0.005 / dt; //0.017 is 1 arc degree;
	if (uConvg < 2e-4) uConvg = 2e-4;
	if (wConvg < 2e-4) wConvg = 2e-4;

	char msg_show[1024] = "\0", msg_log[1024] = "\0", msg[1024] = "\0";

	bool bSpeedVerlet = *status;
	bool reset;

	double maxEk = 9;
	if (_bias_force_::biasOn) maxEk = 64;

	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		if (mdcell->lfmd.m[imol].iHDyn->hdyn_status) continue; // convergent already

		mdcell->lfmd.m[imol].fetch_ksi();

		if (mm->Ek > (mm->nCluster + 5) * ::Ek0Internal * maxEk) {
			ScaleVelocity(*mm, sqrt((mm->nCluster + 5) * ::Ek0Internal / mm->Ek));
			ksi = 0; reset = true;
		}
		else {ksi = mdcell->lfmd.m[imol].iHDyn->nhc->ksi.v[0]; reset = false;}
		if (reset) {
			mdcell->lfmd.m[imol].iHDyn->nhc->reset_NHC();
			mdcell->lfmd.m[imol].Init(*mm);
		}

		// 0 -- use the spatial moment at mass center, 1 -- use the velocity of base cluster
		if (::MM_SpeedVerlet_method == 0) MakeFreeBaseVelocityConsistent(*mm); // adjusting free base velocity
		else if (::MM_SpeedVerlet_method == 1) {
			calClusterVelocity0(*mm); calSpatialMomentum0(*mm, true); // required for Nose-Hover force calculation of the frame
		}

		//extern void cal_Hoover_force(MMOLECULE &mm, LFVERLET_MM_MD& mdpar);
		//cal_Hoover_force(*mm, mdcell->lfmd.m[imol]);
		cal_Hoover_force(mm, ksi, ksi_mc, vmc);

		if (::bBrownian) {V6zero(mm->mDynKin.f_Brown) Brownian_force(mm, 0, mm->nCluster - 1, mm->mDynKin.f_Brown);}

		mdcell->mvb.mm_status.m[imol] = MM_SpeedVerlet_Simplified((reset ? false : bSpeedVerlet), mdcell->mm.m[imol], mdcell->lfmd.m[imol], 
			mdcell->mvb.mmv1.m[imol].cv.m, mdcell->mvb.mmv2.m[imol].cv.m, 
			mdcell->mvb.mm_dvmax.m[imol], mdcell->mvb.mm_dwmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.mm_status.m[imol]) {
			sprintf(msg, "Macromol # %d [%s] : Speed-Verlet fails", imol, mdcell->mm.m[imol].mol);
			show_log(msg, true); //show_infor(msg_show);
		}
		if (reset) mdcell->mvb.mm_sc_status.m[imol] = true;
		else if (mdcell->mvb.mm_dvmax.m[imol] <= uConvg && mdcell->mvb.mm_dwmax.m[imol] <= wConvg) mdcell->mvb.mm_sc_status.m[imol] = true;
		else mdcell->mvb.mm_sc_status.m[imol] = false;
	}

	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		if (mdcell->lfsmd.m[imol].hdyn->hdyn_status) continue;  // convergent already

		if (sm->Ek > 6 * ::Ek0Internal * maxEk) {
			ScaleVelocity(*sm, sqrt(6 * ::Ek0Internal / sm->Ek));
			ksi = 0; reset = true;
		}
		else {ksi = mdcell->lfsmd.m[imol].hdyn->nhc->ksi.v[0]; reset = false;}
		if (reset && mdcell->lfsmd.m[imol].hdyn->nhc->N < 24) {
			mdcell->lfsmd.m[imol].hdyn->nhc->reset_NHC();
			mdcell->lfsmd.m[imol].Init(*sm);
		}

		// adjusting the force on the cluster, because of the constraint on the velocity of the whole cell
		//V3PV3(wmc, sm->r, Vw);
		//V3PV3(wmc, sm->Rg, Vw);
		cal_Hoover_force(sm, ksi, ksi_mc, vmc);

		if (::bBrownian) SM_Brownian_force(sm);

		mdcell->mvb.sm_status.m[imol] = SM_SpeedVerlet_Simplified((reset ? false : bSpeedVerlet), mdcell->sm.m[imol], mdcell->lfsmd.m[imol], 
			mdcell->mvb.sm_dvmax.m[imol], mdcell->mvb.sm_dwmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.sm_status.m[imol]) {
			sprintf(msg, "mol # %d [%s] : Speed-Verlet fails", imol, mdcell->sm.m[imol].mol);
			show_log(msg_show, true); //show_log(msg_log, true);
		}
		if (mdcell->mvb.sm_dvmax.m[imol] <= uConvg && mdcell->mvb.sm_dwmax.m[imol] < wConvg) mdcell->mvb.sm_sc_status.m[imol] = true;
		else mdcell->mvb.sm_sc_status.m[imol] = false;
	}

	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		if (mdcell->lfpmd.m[imol].hdyn->hdyn_status) continue;  // convergent already

		if (pm->Ek > 3 * ::Ek0Internal * maxEk) {
			ScaleVelocity(*pm, sqrt(3 * ::Ek0Internal / pm->Ek));
			ksi = 0; reset = true;
		}
		else {ksi = mdcell->lfpmd.m[imol].hdyn->nhc->ksi.v[0]; reset = false;}
		if (reset && mdcell->lfsmd.m[imol].hdyn->nhc->N < 12) {
			mdcell->lfpmd.m[imol].hdyn->nhc->reset_NHC();
			mdcell->lfpmd.m[imol].Init(*pm);
		}

		// adjusting the force on the cluster, because of the constraint on the velocity of the whole cell
		cal_Hoover_force(pm, ksi, ksi_mc, vmc);

		if (::bBrownian) PM_Brownian_force(pm);

		mdcell->mvb.pm_status.m[imol] = PM_SpeedVerlet_Simplified((reset ? false : bSpeedVerlet), mdcell->pm.m[imol], mdcell->lfpmd.m[imol], 
			mdcell->mvb.pm_dvmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.pm_status.m[imol]) {
			sprintf(msg, "mol # %d [%s] : Speed-Verlet fails", imol, mdcell->pm.m[imol].mol);
			show_log(msg_show, true); //show_log(msg_log, true);
		}
		if (mdcell->mvb.pm_dvmax.m[imol] <= uConvg) mdcell->mvb.pm_sc_status.m[imol] = true;
		else mdcell->mvb.pm_sc_status.m[imol] = false;
	}
}

// **********************************************************************************************************************
// calculating the Hoover Force on each cluster
// **********************************************************************************************************************
extern void cal_Hoover_force(MMOLECULE *mm, double ksi);
extern void cal_Hoover_force(SMOLECULE *sm, double ksi);
extern void cal_Hoover_force(PMOLECULE *pm, double ksi);

void MDCELL_SpeedVerlet0(MDCELL_ThreadJob<MMOL_MD_CELL> *job, bool *status) { // without confinement on total spatial momentum 
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	VECTOR3 f, torq;
	int nm = 0, imol = 0, ihdyn = 0, ic;

	double dt =  mdcell->hdyn.m[0].nhc->dt;
	double uConvg = 0.005 / dt, wConvg = 0.005 / dt; //0.017 is 1 arc degree;
	if (uConvg < 2e-4) uConvg = 2e-4;
	if (wConvg < 2e-4) wConvg = 2e-4;

	char msg_show[1024] = "\0", msg_log[1024] = "\0", msg[1024] = "\0";

	bool bSpeedVerlet = *status;

	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		if (mdcell->lfmd.m[imol].iHDyn->hdyn_status) continue; // convergent already

		mdcell->lfmd.m[imol].fetch_ksi();

		cal_Hoover_force(mm, mdcell->lfmd.m[imol].ksi);

		if (::bBrownian) Brownian_force(mm, 0, mm->nCluster - 1, mm->mDynKin.f_Brown);

		mdcell->mvb.mm_status.m[imol] = MM_SpeedVerlet_Simplified(bSpeedVerlet, mdcell->mm.m[imol], mdcell->lfmd.m[imol], 
			mdcell->mvb.mmv1.m[imol].cv.m, mdcell->mvb.mmv2.m[imol].cv.m, 
			mdcell->mvb.mm_dvmax.m[imol], mdcell->mvb.mm_dwmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.mm_status.m[imol]) {
			sprintf(msg, "Macromol # %d [%s] : Speed-Verlet fails", imol, mdcell->mm.m[imol].mol);
			show_log(msg, true); //show_infor(msg_show);
		}
		if (mdcell->mvb.mm_dvmax.m[imol] <= uConvg && mdcell->mvb.mm_dwmax.m[imol] <= wConvg) mdcell->mvb.mm_sc_status.m[imol] = true;
		else mdcell->mvb.mm_sc_status.m[imol] = false;
	}

	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		if (mdcell->lfsmd.m[imol].hdyn->hdyn_status) continue;  // convergent already

		cal_Hoover_force(sm, mdcell->lfsmd.m[imol].hdyn->nhc->ksi.v[0]);

		if (::bBrownian) SM_Brownian_force(sm);

		mdcell->mvb.sm_status.m[imol] = SM_SpeedVerlet_Simplified(bSpeedVerlet, mdcell->sm.m[imol], mdcell->lfsmd.m[imol], 
			mdcell->mvb.sm_dvmax.m[imol], mdcell->mvb.sm_dwmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.sm_status.m[imol]) {
			sprintf(msg, "mol # %d [%s] : Speed-Verlet fails", imol, mdcell->sm.m[imol].mol);
			show_log(msg_show, true); //show_log(msg_log, true);
		}
		if (mdcell->mvb.sm_dvmax.m[imol] <= uConvg && mdcell->mvb.sm_dwmax.m[imol] < wConvg) mdcell->mvb.sm_sc_status.m[imol] = true;
		else mdcell->mvb.sm_sc_status.m[imol] = false;
	}

	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		if (mdcell->lfpmd.m[imol].hdyn->hdyn_status) continue;  // convergent already

		cal_Hoover_force(pm, mdcell->lfpmd.m[imol].hdyn->nhc->ksi.v[0]);

		if (::bBrownian) PM_Brownian_force(pm);

		mdcell->mvb.pm_status.m[imol] = PM_SpeedVerlet_Simplified(bSpeedVerlet, mdcell->pm.m[imol], mdcell->lfpmd.m[imol], 
			mdcell->mvb.pm_dvmax.m[imol]); // translation and rotation of whole molecule are 1.5kT
		if (!mdcell->mvb.pm_status.m[imol]) {
			sprintf(msg, "mol # %d [%s] : Speed-Verlet fails", imol, mdcell->pm.m[imol].mol);
			show_log(msg_show, true); //show_log(msg_log, true);
		}
		if (mdcell->mvb.pm_dvmax.m[imol] <= uConvg) mdcell->mvb.pm_sc_status.m[imol] = true;
		else mdcell->mvb.pm_sc_status.m[imol] = false;
	}
}


void MDCELL_Verlet(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0, ihdyn = 0;

	//char msg_show[1024] = "\0", msg_log[1024] = "\0";

	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		Verlet(mdcell->lfmd.m[imol], mdcell->mm.m[imol]); // Verlet algorithm to calculate the ta and tp for next timestep
	}

	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		Verlet(mdcell->lfsmd.m[imol], mdcell->sm.m[imol]); // Verlet algorithm to calculate the ta and tp for next timestep
	}

	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		Verlet(mdcell->lfpmd.m[imol], mdcell->pm.m[imol]); // Verlet algorithm to calculate the ta and tp for next timestep
	}
}

void MDCELL_CalMolKineticEnergy(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0;
	int method = ::MM_SpeedVerlet_method;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		if (method == 0) MakeFreeBaseVelocityConsistent(*mm);
		else if (method == 1) {calClusterVelocity0(*mm); calSpatialMomentum0(*mm, true);}
		calKineticEnergy(*mm);
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		calKineticEnergy(*sm);
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		calKineticEnergy(*pm);
	}
}

void MDCELL_CalMolNetForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	int nm = 0, imol = 0;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= job->job->mdcell->mm.n) continue;
		mm = job->job->mdcell->mm.m + imol;
		CalNetForce(*mm, 0, mm->nCluster); // calculate the net force also here !
	}
	/*
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < smIndx.n; nm++) {
		imol = smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		sm->Ek = calKineticEnergy(*sm);
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < pmIndx.n; nm++) {
		imol = pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		pm->Ek = calKineticEnergy(*pm);
	}
	*/
}
/*
void MDCELL_ClusterRealForceBackup(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0, ic;
	BASIC_CLUSTER *pc = NULL;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = (BASIC_CLUSTER*)(mm->cluster + ic);
			V62V6(pc->dyn->fc, pc->dyn->fc0)
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; pc = sm->c;
		V62V6(pc->dyn->fc, pc->dyn->fc0)
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol; pc = pm->c;
		V62V6(pc->dyn->fc, pc->dyn->fc0)
	}
}
*/
void MDCELL_CalInertiaMoment(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		CalMolInertiaMoment(*mm);
		//NEIMO_MOMENT_MATRIX(*mm);
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		sm->calInertiaMoment();
	}
	/*
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < pmIndx.n; nm++) {
		imol = pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		pm->Ek = calKineticEnergy(*pm);
	}
	*/
}

void MDCELL_CalCellInertiaTensor(MDCELL_ThreadJob<MMOL_MD_CELL> *job, MATRIX<3> *mI) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	MATRIX<3> Iw;
	VECTOR3 Rcm, rc;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			Rcm.v[0] = pc->cm.v[0] + pc->dr_cmm.v[0];
			Rcm.v[1] = pc->cm.v[1] + pc->dr_cmm.v[1];
			Rcm.v[2] = pc->cm.v[2] + pc->dr_cmm.v[2];

			VECT3(mdcell->rc, Rcm, rc)

			Iw.m[0][0] += pc->M * (rc.v[1] * rc.v[1] + rc.v[2] * rc.v[2]);
			Iw.m[1][1] += pc->M * (rc.v[0] * rc.v[0] + rc.v[2] * rc.v[2]);
			Iw.m[2][2] += pc->M * (rc.v[0] * rc.v[0] + rc.v[1] * rc.v[1]);
			Iw.m[0][1] -= pc->M * rc.v[0] * rc.v[1];
			Iw.m[0][2] -= pc->M * rc.v[0] * rc.v[2];
			Iw.m[1][2] -= pc->M * rc.v[1] * rc.v[2];
		}
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		M = sm->c->M;

		VECT3(mdcell->rc, sm->r, rc)

		Iw.m[0][0] += M * (rc.v[1] * rc.v[1] + rc.v[2] * rc.v[2]);
		Iw.m[1][1] += M * (rc.v[0] * rc.v[0] + rc.v[2] * rc.v[2]);
		Iw.m[2][2] += M * (rc.v[0] * rc.v[0] + rc.v[1] * rc.v[1]);
		Iw.m[0][1] -= M * rc.v[0] * rc.v[1];
		Iw.m[0][2] -= M * rc.v[0] * rc.v[2];
		Iw.m[1][2] -= M * rc.v[1] * rc.v[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		M = pm->c->M;

		VECT3(mdcell->rc, pm->r, rc)

		Iw.m[0][0] += M * (rc.v[1] * rc.v[1] + rc.v[2] * rc.v[2]);
		Iw.m[1][1] += M * (rc.v[0] * rc.v[0] + rc.v[2] * rc.v[2]);
		Iw.m[2][2] += M * (rc.v[0] * rc.v[0] + rc.v[1] * rc.v[1]);
		Iw.m[0][1] -= M * rc.v[0] * rc.v[1];
		Iw.m[0][2] -= M * rc.v[0] * rc.v[2];
		Iw.m[1][2] -= M * rc.v[1] * rc.v[2];
	}
	mI->m[0][0] = Iw.m[0][0];
	mI->m[1][1] = Iw.m[1][1];
	mI->m[2][2] = Iw.m[2][2];
	mI->m[0][1] = Iw.m[0][1];
	mI->m[0][2] = Iw.m[0][2];
	mI->m[1][2] = Iw.m[1][2];
	mI->m[1][0] = Iw.m[0][1];
	mI->m[2][0] = Iw.m[0][2];
	mI->m[2][1] = Iw.m[1][2];
}

void MDCELL_HDyn(MMOL_MD_CELL* mdcell, ARRAY<int> &hdynIndx) {
	int i = 0, ihdyn = 0, im, imol = 0;
	double Ek1, Ek2, Ek3;
	double T = 0;

	char msg[256] = "\0";
	bool status = true;

	for (i = 0; i < hdynIndx.n; i++) {
		ihdyn = hdynIndx.m[i]; if (ihdyn >= mdcell->hdyn.n) continue;
		if (mdcell->hdyn.m[ihdyn].hdyn_status) continue;
		if (!mdcell->hdyn.m[ihdyn].mol_coupled()) {
			mdcell->hdyn.m[ihdyn].hdyn_status = true; continue;
		}

		status = true;
		Ek1 = 0; Ek2 = 0; Ek3 = 0;
		for (im = 0; im < mdcell->hdyn.m[ihdyn].mmIndx.n; im++) {
			imol = mdcell->hdyn.m[ihdyn].mmIndx.m[im];
			 Ek1 += mdcell->mm.m[imol].Ek;

			if (!mdcell->mvb.mm_status.m[imol]) status = false;
		}
		for (im = 0; im < mdcell->hdyn.m[ihdyn].smIndx.n; im++) {
			imol = mdcell->hdyn.m[ihdyn].smIndx.m[im];
			Ek1 += mdcell->sm.m[imol].Ek;

			if (!mdcell->mvb.sm_status.m[imol]) status = false;
		}
		for (im = 0; im < mdcell->hdyn.m[ihdyn].pmIndx.n; im++) {
			imol = mdcell->hdyn.m[ihdyn].pmIndx.m[im];
			Ek1 += mdcell->pm.m[imol].Ek;

			if (!mdcell->mvb.pm_status.m[imol]) status = false;
		}

		//if (!status) mdcell->hdyn.m[ihdyn].nhc->reset_NHC();
		T = Ek1 / (mdcell->hdyn.m[ihdyn].nhc->N * ::Ek0Internal);
		mdcell->hdyn.m[ihdyn].nhc->set_dT(T - 1);
		mdcell->hdyn.m[ihdyn].nhc->verlet_NHC();

		//sprintf(msg, "Hdyn #%d: T[%f], ksi[%f], vksi[%f]", ihdyn, T, mdcell->hdyn.m[ihdyn].nhc->ksi.v[0], mdcell->hdyn.m[ihdyn].nhc->vksi.v[0]);
		//show_log(msg, true);

		//if (!status) {
		//	sprintf(msg, "NHC #%d is reset", ihdyn); show_log(msg, true);
		//}
	}
}

void MDCELL_CheckHDynCovergent(MMOL_MD_CELL* mdcell, ARRAY<int> &hdynIndx) {
	int i = 0, ihdyn = 0, im, imol = 0;

	for (i = 0; i < hdynIndx.n; i++) {
		ihdyn = hdynIndx.m[i]; if (ihdyn >= mdcell->hdyn.n) continue;

		if (mdcell->hdyn.m[ihdyn].hdyn_status) continue;

		mdcell->hdyn.m[ihdyn].hdyn_status = true;
		for (im = 0; im < mdcell->hdyn.m[ihdyn].mmIndx.n; im++) {
			imol = mdcell->hdyn.m[ihdyn].mmIndx.m[im];
			if (!mdcell->mvb.mm_status.m[imol] || mdcell->mvb.mm_sc_status.m[imol]) continue; // LFMD speed is rescaled, or the speed/accel is covergent
			else {mdcell->hdyn.m[ihdyn].hdyn_status = false; break;}
		}
		if (!mdcell->hdyn.m[ihdyn].hdyn_status) break;

		for (im = 0; im < mdcell->hdyn.m[i].smIndx.n; im++) {
			imol = mdcell->hdyn.m[i].smIndx.m[im];
			if (!mdcell->mvb.sm_status.m[imol] || mdcell->mvb.sm_sc_status.m[imol]) continue; // LFMD speed is rescaled, or the speed/accel is covergent
			else {mdcell->hdyn.m[ihdyn].hdyn_status = false; break;}
		}
		if (!mdcell->hdyn.m[ihdyn].hdyn_status) break;

		for (im = 0; im < mdcell->hdyn.m[i].pmIndx.n; im++) {
			imol = mdcell->hdyn.m[i].pmIndx.m[im];
			if (!mdcell->mvb.pm_status.m[imol] || mdcell->mvb.pm_sc_status.m[imol]) continue; // LFMD speed is rescaled, or the speed/accel is covergent
			else {mdcell->hdyn.m[ihdyn].hdyn_status = false; break;}
		}
	}
}

void MDCELL_AcceptCurrentDynForNextStep(MMOL_MD_CELL* mdcell, ARRAY<int> &hdynIndx) {
	int i = 0, ihdyn = 0;

	for (i = 0; i < hdynIndx.n; i++) {
		ihdyn = hdynIndx.m[i]; if (ihdyn >= mdcell->hdyn.n) continue;
		mdcell->hdyn.m[ihdyn].nhc->accept_nextstep();
	}
}

bool LeapFrogVerlet_SelfConsistent_SpeedAccel(MMOL_MD_CELL *mdcell, MDCELL_ThreadJob<MMOL_MD_CELL> *job_mdcell, MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> *hdyn_thread_var, bool bSpeedVerlet, int nThreads) {
	extern void MDCELL_SpatialMomentum(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 3> *P);

	ARRAY< SVECTOR<double, 3> > aP; aP.SetArray(nThreads);
	ARRAY<bool> status; status.SetArray(nThreads);

	// before velocity verlet, backup the real force on the cluster, because the force on cluster (fc) will be adjusted by mass center velocity constraint
	//MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ClusterRealForceBackup), job_mdcell, nThreads);

	extern void MDCELL_MMClusterVelocity(MDCELL_ThreadJob<MMOL_MD_CELL> *job);

	bool bDynConvergent = false;
	int iloop = 0, max_loops = 8, i;

	VECTOR3 Pmc;

	float Ek = 0, dT;

	for (i = 0; i < mdcell->hdyn.n; i++) mdcell->hdyn.m[i].hdyn_status = false;
	while (!bDynConvergent && iloop < max_loops) {
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_MMClusterVelocity), job_mdcell, nThreads);
		// calculate the kinetic energy of each molecule, 
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalMolKineticEnergy), job_mdcell, nThreads); 

		// run NHC to get a new ksi
		MDCELL_HDyn_MultiThreadOperate<MMOL_MD_CELL>((void*)(&MDCELL_HDyn), hdyn_thread_var, nThreads);

		// run NHC for the mass-center velocity of the whole cell
		if (mdcell->bHDynMC_T) {
			for (i = 0; i < nThreads; i++) memset(aP.m[i].v, 0, SIZE_V3);
			MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 3> >((void*)(&MDCELL_SpatialMomentum), job_mdcell, nThreads, aP.m);
			V3zero(mdcell->Vmc) V3zero(Pmc)
			for (i = 0; i < nThreads; i++) {
				Pmc.v[0] += aP.m[i].v[0];
				Pmc.v[1] += aP.m[i].v[1];
				Pmc.v[2] += aP.m[i].v[2];
			}
			mdcell->Vmc.v[0] = Pmc.v[0] / mdcell->M; mdcell->Vmc.v[1] = Pmc.v[1] / mdcell->M; mdcell->Vmc.v[2] = Pmc.v[2] / mdcell->M;


			mdcell->Ekt_mc = 0;
			if (mdcell->ix_mct) mdcell->Ekt_mc += mdcell->Vmc.v[0] * mdcell->Vmc.v[0];
			if (mdcell->iy_mct) mdcell->Ekt_mc += mdcell->Vmc.v[1] * mdcell->Vmc.v[1];
			if (mdcell->iz_mct) mdcell->Ekt_mc += mdcell->Vmc.v[2] * mdcell->Vmc.v[2];
			//V3ABS2(mdcell->Vmc, Ek) 
			mdcell->Ekt_mc *= mdcell->M * 0.5 * unit_Ek_kT;
				
			if (mdcell->Ethermal_mct < 0.1) {
				if (mdcell->iz_mct == 1) mdcell->mcHDyn.nhc->set_dT(mdcell->Ekt_mc / mdcell->Mthermal_mc);
				else if (mdcell->iy_mct == 1) mdcell->mcHDyn.nhc->set_dT(mdcell->Ekt_mc / mdcell->Mthermal_mc);
				else if (mdcell->ix_mct == 1) mdcell->mcHDyn.nhc->set_dT(mdcell->Ekt_mc / mdcell->Mthermal_mc);
			}
			else mdcell->mcHDyn.nhc->set_dT((mdcell->Ekt_mc - mdcell->Ethermal_mct) / mdcell->Mthermal_mc); // supposedly the spatial momentum should be 0
				
			mdcell->mcHDyn.nhc->verlet_NHC(1e-6);
		}
		
		// Do the Speed-Verlet
		for (i = 0; i < nThreads; i++) status.m[i] = bSpeedVerlet;
		MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, bool>((void*)(&MDCELL_SpeedVerlet), job_mdcell, nThreads, status.m); 

		// check whether the speed and accel are convergent
		MDCELL_HDyn_MultiThreadOperate<MMOL_MD_CELL>((void*)(&MDCELL_CheckHDynCovergent), hdyn_thread_var, nThreads);
		bDynConvergent = true;
		for (i = 0; i < mdcell->hdyn.n; i++) {
			if (!mdcell->hdyn.m[i].hdyn_status) {bDynConvergent = false; break;}
		}
		iloop += 1;
	}
	return bDynConvergent;
}

// without confinement on spatial momentum
bool LeapFrogVerlet0_SelfConsistent_SpeedAccel(MMOL_MD_CELL *mdcell, MDCELL_ThreadJob<MMOL_MD_CELL> *job_mdcell, MDCELL_HDYN_THREAD_VARS<MMOL_MD_CELL> *hdyn_thread_var, bool bSpeedVerlet, int nThreads) {
	ARRAY<bool> status; status.SetArray(nThreads);

	// before velocity verlet, backup the real force on the cluster, because the force on cluster (fc) will be adjusted by mass center velocity constraint
	//MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_ClusterRealForceBackup), job_mdcell, nThreads);

	extern void MDCELL_MMClusterVelocity(MDCELL_ThreadJob<MMOL_MD_CELL> *job);
	MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_MMClusterVelocity), job_mdcell, nThreads);
	// calculate the kinetic energy of each molecule, 
	// initializing the velocities of each cluster and the mass-center speed based on the Spatial Moment of frame
	MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalMolKineticEnergy), job_mdcell, nThreads); 

	bool bDynConvergent = false;
	int iloop = 0, max_loops = 8, i;

	float Ek = 0;

	for (i = 0; i < mdcell->hdyn.n; i++) mdcell->hdyn.m[i].hdyn_status = false;
	while (!bDynConvergent && iloop < max_loops) {
		// run NHC to get a new ksi
		MDCELL_HDyn_MultiThreadOperate<MMOL_MD_CELL>((void*)(&MDCELL_HDyn), hdyn_thread_var, nThreads);

		// Do the Speed-Verlet
		for (i = 0; i < nThreads; i++) status.m[i] = bSpeedVerlet;
		MTOperate1< MDCELL_ThreadJob<MMOL_MD_CELL>, bool>((void*)(&MDCELL_SpeedVerlet0), job_mdcell, nThreads, status.m); 

		// check whether the speed and accel are convergent
		MDCELL_HDyn_MultiThreadOperate<MMOL_MD_CELL>((void*)(&MDCELL_CheckHDynCovergent), hdyn_thread_var, nThreads);
		bDynConvergent = true;
		for (i = 0; i < mdcell->hdyn.n; i++) {
			if (!mdcell->hdyn.m[i].hdyn_status) {bDynConvergent = false; break;}
		}
		iloop += 1;
	}
	return bDynConvergent;
}

void InitClusterForce(MMOLECULE *mm) {
	int ic, iatom, i;
	CLUSTER *pc = NULL;
	MATOM * patom = NULL;

	for (ic = 0; ic < mm->nCluster; ic++) {
		pc = mm->cluster + ic;
		V6zero(pc->dyn->fc) V6zero(pc->dyn->f_Hoover)
		V6zero(pc->dyn->f_Brown)
		V6zero(pc->dyn->f)
		if (pc->fphinge != NULL) {V6zero((*(pc->fphinge)))}
		for (iatom = 0; iatom < pc->nAtoms; iatom++) {
			patom = pc->atom + iatom; 
			V3zero(patom->F) V3zero(patom->E)
			for (i = 0; i < patom->f.n; i++) {V3zero(patom->f.m[i])}
			patom->bShocked = false;
		}
		pc->bShocked = false;
	}
}

void MDCELL_InitClusterForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0, nc = 0, iatom, i;
	CLUSTER *pc = NULL;

	MMOLECULE *mm = NULL;
	MATOM *patom = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = mm->cluster + nc;
			V6zero(pc->dyn->fc) V6zero(pc->dyn->f_Hoover)
			if (pc->fphinge != NULL) {V6zero((*(pc->fphinge)))}
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom; 
				V3zero(patom->F) V3zero(patom->E)
				for (i = 0; i < patom->f.n; i++) {V3zero(patom->f.m[i])}
				patom->bShocked = false;
			}
			pc->bShocked = false;
		}
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		V6zero(sm->c->dyn->fc) V6zero(sm->c->dyn->f_Hoover) 
		for (iatom = 0; iatom < sm->c->nAtoms; iatom++) {
			patom = sm->c->atom + iatom; 
			V3zero(patom->F) V3zero(patom->E)
			for (i = 0; i < patom->f.n; i++) {V3zero(patom->f.m[i])}
			patom->bShocked = false;
		}
		sm->c->bShocked = false;
	}
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		V6zero(pm->c->dyn->fc) V6zero(pm->c->dyn->f_Hoover) 
		for (iatom = 0; iatom < pm->c->nAtoms; iatom++) {
			patom = pm->c->atom + iatom; 
			V3zero(patom->F) V3zero(patom->E)
			for (i = 0; i < patom->f.n; i++) {V3zero(patom->f.m[i])}
			patom->bShocked = false;
		}
		pm->c->bShocked = false;
	}
}

void MDCELL_ScaleClusterForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0, nc = 0, iatom;
	CLUSTER *pc = NULL;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = mm->cluster + nc;
			scale_force_cluster(pc);
		}
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		scale_force_cluster(sm->c);
	}
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		scale_force_cluster(pm->c);
	}
}

void MDCELL_ClusterForce(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0, nc = 0, iatom;
	CLUSTER *pc = NULL;

	MMOLECULE *mm = NULL;
	MATOM *patom = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = mm->cluster + nc;
			V6zero(pc->dyn->fc)
			CLUSTER_TOTAL_FORCE(pc);
		}
	}
	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		V6zero(sm->c->dyn->fc)
		CLUSTER_TOTAL_FORCE(sm->c);
	}
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		if (::bMTForce) {
			V6zero(pm->c->dyn->fc) 
			CLUSTER_TOTAL_FORCE(pm->c);
		}
		else {
			pm->c->dyn->fc.v[3] = pm->c->atom[0].F.v[0];
			pm->c->dyn->fc.v[4] = pm->c->atom[0].F.v[1];
			pm->c->dyn->fc.v[5] = pm->c->atom[0].F.v[2];
		}
	}
}

void InitClusterForceInCMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell) {
	CMM_CELL<BASIC_CLUSTER>* pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	int iatom = 0;
	MATOM *patom = NULL;

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > it;

	for (int nCell = n1Cell; nCell <= n2Cell; nCell++) {
		if (nCell >= cell.acell.n) break;
		pcell = cell.acell.m[nCell]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current();
			V6zero(pc->dyn->fc) V6zero(pc->dyn->f_Hoover)
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = pc->atom + iatom; V3zero(patom->F)
			}
			pc->dyn->overlap = false;
			it.next_iterate();
		}
	}
}

void MM_atomic_rg(MMOLECULE *mm, int nMM) {
	int imol, ic, iatom;
	CLUSTER *pc;
	MATOM *patom;
	for (imol = 0; imol < nMM; imol++) {
		for (ic = 0; ic < mm[imol].nCluster; ic++) {
			pc = mm[imol].cluster + ic;
			for (iatom = 0; iatom < pc->nAtoms; iatom++) {
				patom = mm[imol].cluster[ic].atom + iatom;
				V3plusV3(patom->r, pc->dr_cmm, patom->rg)
			}
		}
	}
}

void SM_atomic_rg(SMOLECULE *sm, int nSM) {
	int imol, iatom;
	MATOM *patom;
	for (imol = 0; imol < nSM; imol++) {
		for (iatom = 0; iatom < sm[imol].c->nAtoms; iatom++) {
			patom = sm[imol].c->atom + iatom;
			V3plusV3(patom->r, sm[imol].c->dr_cmm, patom->rg)
		}
	}
}

void PM_atomic_rg(PMOLECULE *pm, int nPM) {
	int imol;
	MATOM *patom;
	for (imol = 0; imol < nPM; imol++) {
		patom = pm[imol].c->atom;
		V32V3(patom->r, patom->rg)
	}
}
/*
void check_atom_rg(MMOL_MD_CELL &mmcell, int nThreads) {
	if (mmcell.mm.n > 0) {
		MOperate<MMOLECULE>((void*)(&MM_atomic_rg), mmcell.mm.m, mmcell.mm.n, (mmcell.mm.n > nThreads ? nThreads : mmcell.mm.n));
	}
	if (mmcell.sm.n > 0) {
		MOperate<SMOLECULE>((void*)(&SM_atomic_rg), mmcell.sm.m, mmcell.sm.n, (mmcell.sm.n > nThreads ? nThreads : mmcell.sm.n));
	}
	if (mmcell.pm.n > 0) {
		MOperate<PMOLECULE>((void*)(&PM_atomic_rg), mmcell.pm.m, mmcell.pm.n, (mmcell.pm.n > nThreads ? nThreads : mmcell.pm.n));
	}
}
*/
void update_atom_rg(MMOL_MD_CELL &mmcell) {
	if (mmcell.mm.n > 0) MM_atomic_rg(mmcell.mm.m, mmcell.mm.n);
	if (mmcell.sm.n > 0) SM_atomic_rg(mmcell.sm.m, mmcell.sm.n);
	if (mmcell.pm.n > 0) PM_atomic_rg(mmcell.pm.m, mmcell.pm.n);
}

void cluster_atomic_rg(BASIC_CLUSTER** bc, int ncs) {
	int ic, iatom;
	MATOM *patom;
	for (ic = 0; ic < ncs; ic++) {
		for (iatom = 0; iatom < bc[ic]->nAtoms; iatom++) {
			patom = bc[ic]->atom + iatom;
			V3plusV3(patom->r, bc[ic]->dr_cmm, patom->rg)
		}
	}
}

void check_atom_rg(MMOL_MD_CELL &mmcell, int nThreads) {
	if (mmcell.bcluster.n < 100) cluster_atomic_rg(mmcell.bcluster.m, mmcell.bcluster.n);
	else MOperate<BASIC_CLUSTER*>((void*)(&cluster_atomic_rg), mmcell.bcluster.m, mmcell.bcluster.n, nThreads);
}

void MultiMM_NEIMO_MOMENT(MMOLECULE *mm, int nmm) {
	for (int imol = 0; imol < nmm; imol++) {
		CalMolInertiaMoment(mm[imol]);
		//NEIMO_MOMENT_MATRIX(mm[imol]);
	}
}

void MultiMMP_NEIMO_MOMENT(MMOLECULE** mm, int nmm) {
	for (int imol = 0; imol < nmm; imol++) {
		CalMolInertiaMoment(*(mm[imol]));
		//NEIMO_MOMENT_MATRIX(*(mm[imol]));
	}
}

/*
void MM_Torsion(MMOLECULE *mm, int nMM, double& Etorsion) {
	int nm = 0;
	double Et = 0;
	for (nm = 0; nm < nMM; nm++) {
		mm->calTorsion(Et);
		Etorsion += Et;
	}
}
*/
void MDCELL_CheckRotAngle(MDCELL_ThreadJob<MMOL_MD_CELL> *job) {
	int nm = 0, imol = 0, ihdyn = 0;

	char msg_show[1024] = "\0", msg_log[1024] = "\0";

	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= job->job->mdcell->mm.n) continue;
		check_angle(job->job->mdcell->mm.m[imol], job->job->mdcell->lfmd.m[imol]);
	}
	/*
	for (nm = 0; nm < smIndx.n; nm++) {
		imol = smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		// work to be done
	}

	for (nm = 0; nm < pmIndx.n; nm++) {
		imol = pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		// work to be done
	}
	*/
}

void MMCELL_CMM_cluster_allocate(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator, bool *status) {
	int i, imol;
	bool tstatus;
	*status = true;
	CMM_array_reset_storage<BASIC_CLUSTER>(allocator->cmm);
	for (i = 0; i < allocator->job->mmIndx.n; i++) {
		imol = allocator->job->mmIndx.m[i];
		if (imol >= allocator->job->mdcell->mm.n) continue;
		tstatus = MM_AllocateCMM_array<BASIC_CLUSTER, MMOLECULE>(allocator->job->mdcell->mm.m + imol, allocator->cmm); 
		if (!tstatus) *status = false;
	}
	for (i = 0; i < allocator->job->smIndx.n; i++) {
		imol = allocator->job->smIndx.m[i];
		if (imol >= allocator->job->mdcell->sm.n) continue;
		tstatus = SM_AllocateCluster_CMM_array(allocator->job->mdcell->sm.m + imol, 1, allocator->cmm); 
		if (!tstatus) *status = false;
	}
	for (i = 0; i < allocator->job->pmIndx.n; i++) {
		imol = allocator->job->pmIndx.m[i];
		if (imol >= allocator->job->mdcell->pm.n) continue;
		tstatus = PM_AllocateCluster_CMM_array(allocator->job->mdcell->pm.m + imol, 1, allocator->cmm); 
		if (!tstatus) *status = false;
	}
}

void MMCELL_CMM_CheckRelship(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator) {
	int i, imol, nc;
	MMOLECULE *mm = NULL;
	BASIC_CLUSTER *pc = NULL;
	for (i = 0; i < allocator->job->mmIndx.n; i++) {
		imol = allocator->job->mmIndx.m[i];
		if (imol >= allocator->job->mdcell->mm.n) continue;
		mm = allocator->job->mdcell->mm.m + imol;
		for (nc = 0; nc < mm->nCluster; nc++) {
			CheckClusterRelship((BASIC_CLUSTER*)(mm->cluster + nc), *(allocator->cmm)); 
		}
	}
	for (i = 0; i < allocator->job->smIndx.n; i++) {
		imol = allocator->job->smIndx.m[i];
		if (imol >= allocator->job->mdcell->sm.n) continue;
		CheckClusterRelship(allocator->job->mdcell->sm.m[imol].c, *(allocator->cmm)); 
	}
	for (i = 0; i < allocator->job->pmIndx.n; i++) {
		imol = allocator->job->pmIndx.m[i];
		if (imol >= allocator->job->mdcell->pm.n) continue;
		CheckClusterRelship(allocator->job->mdcell->pm.m[imol].c, *(allocator->cmm)); 
	}
}

void MMCELL_CMM_CheckImpSolvRelship(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *allocator) {
	extern float rSolvent, dSolv_max;
	//float r0_solvent = rSolvent;
	CMM_CELL3D<BASIC_CLUSTER> *cmm = allocator->cmm;
	float xd[3] = {cmm->xd[0] / cmm->bcell.nx, cmm->xd[1] / cmm->bcell.ny, cmm->xd[2] / cmm->bcell.nz};
	int dnc[3];
	float rcut = dSolv_max;
	
	if (::local_relax) {dnc[0] = 1; dnc[1] = 1; dnc[2] = 1;}
	else {
		dnc[0] = int(rcut / xd[0] + 0.5f); if (dnc[0] == 0) dnc[0] = 1; if (dnc[0] * xd[0] < rcut) dnc[0] += 1;
		dnc[1] = int(rcut / xd[1] + 0.5f); if (dnc[1] == 0) dnc[1] = 1; if (dnc[1] * xd[1] < rcut) dnc[1] += 1;
		dnc[2] = int(rcut / xd[2] + 0.5f); if (dnc[2] == 0) dnc[2] = 1; if (dnc[2] * xd[2] < rcut) dnc[2] += 1;
		dnc[0] += 1; dnc[1] += 1; dnc[2] += 1;
	}

	int i, imol, nc;
	MMOLECULE *mm = NULL;
	BASIC_CLUSTER *pc = NULL;
	for (i = 0; i < allocator->job->mmIndx.n; i++) {
		imol = allocator->job->mmIndx.m[i];
		if (imol >= allocator->job->mdcell->mm.n) continue;
		mm = allocator->job->mdcell->mm.m + imol;
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = (BASIC_CLUSTER*)(mm->cluster + nc);
			ClusterImpSolvNeighborsInCell<BASIC_CLUSTER>(cmm, pc, dnc); 
		}
	}
	for (i = 0; i < allocator->job->smIndx.n; i++) {
		imol = allocator->job->smIndx.m[i];
		if (imol >= allocator->job->mdcell->sm.n) continue;
		pc = allocator->job->mdcell->sm.m[imol].c;
		ClusterImpSolvNeighborsInCell<BASIC_CLUSTER>(cmm, pc, dnc); 
	}
	for (i = 0; i < allocator->job->pmIndx.n; i++) {
		imol = allocator->job->pmIndx.m[i];
		if (imol >= allocator->job->mdcell->pm.n) continue;
		pc = allocator->job->mdcell->pm.m[imol].c;
		ClusterImpSolvNeighborsInCell<BASIC_CLUSTER>(cmm, pc, dnc); 
	}
}

void MMCELL_cluster_impsolv_interact(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *mdcell_impsolv, double *U) {
	int i, imol, nc;
	MMOLECULE *mm = NULL;
	BASIC_CLUSTER *pc = NULL;
	CMM_CELL3D<BASIC_CLUSTER> *cmm = mdcell_impsolv->cmm;
	MMOL_MD_CELL *mdcell = mdcell_impsolv->job->mdcell;
	// calculate SASA of each cluster
	for (i = 0; i < mdcell_impsolv->job->mmIndx.n; i++) {
		imol = mdcell_impsolv->job->mmIndx.m[i];
		if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = (BASIC_CLUSTER*)(mm->cluster + nc);
			Cluster_SASA(pc, cmm);
		}
	}
	for (i = 0; i < mdcell_impsolv->job->smIndx.n; i++) {
		imol = mdcell_impsolv->job->smIndx.m[i];
		if (imol >= mdcell->sm.n) continue;
		pc = mdcell->sm.m[imol].c;
		Cluster_SASA(pc, cmm);
	}
	for (i = 0; i < mdcell_impsolv->job->pmIndx.n; i++) {
		imol = mdcell_impsolv->job->pmIndx.m[i];
		if (imol >= mdcell->pm.n) continue;
		pc = mdcell->pm.m[imol].c;
		Cluster_SASA(pc, cmm);
	}

	// calculate the force & energy
	double E = 0;
	for (i = 0; i < mdcell_impsolv->job->mmIndx.n; i++) {
		imol = mdcell_impsolv->job->mmIndx.m[i];
		if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (nc = 0; nc < mm->nCluster; nc++) {
			pc = (BASIC_CLUSTER*)(mm->cluster + nc);
			ClusterImpSolvForce(pc);
			E += pc->E_ImplicitSolvation;
		}
	}
	for (i = 0; i < mdcell_impsolv->job->smIndx.n; i++) {
		imol = mdcell_impsolv->job->smIndx.m[i];
		if (imol >= mdcell->sm.n) continue;
		pc = mdcell->sm.m[imol].c;
		ClusterImpSolvForce(pc);
		E += pc->E_ImplicitSolvation;
	}
	for (i = 0; i < mdcell_impsolv->job->pmIndx.n; i++) {
		imol = mdcell_impsolv->job->pmIndx.m[i];
		if (imol >= mdcell->pm.n) continue;
		pc = mdcell->pm.m[imol].c;
		ClusterImpSolvForce(pc);
		E += pc->E_ImplicitSolvation;
	}
	*U = E;
}

//#error the function can work with reference as variables, while not pointer for variables, how can it work with MTOperate2 ?
//void MMCELL_MM_Torsion(MDCELL_ThreadJob<MMOL_MD_CELL> *job, TorsionVar &var, InteractRes &ires) {
void MMCELL_MM_Torsion(MDCELL_ThreadJob<MMOL_MD_CELL> *job, TorsionVar *var, InteractRes *ires) {
	int i, imol;
	MMOLECULE *mm = NULL;
	ires->reset();
	InteractRes tres;
	for (i = 0; i < job->job->mmIndx.n; i++) {
		imol = job->job->mmIndx.m[i];
		if (imol >= job->job->mdcell->mm.n) continue;
		mm = job->job->mdcell->mm.m + imol;
		TorsionInteract(mm, 0, mm->nCluster, *var, tres);
		(*ires) += tres;
	}
}

void MMCELL_ExplicitInteract(MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > *job, INTERACT *interact) {
	int i, imol, ic;
	MMOLECULE *mm = NULL;
	CMM_CELL3D<BASIC_CLUSTER> *cmm = job->cmm;
	BASIC_CLUSTER * bpc = NULL;
	InteractRes tres;
	interact->LJ_res.reset(); interact->eres.reset();
	interact->evar.esv1.bEwald = false; // NON-EwaldSum, but direct interaction
	for (i = 0; i < job->job->mmIndx.n; i++) {
		imol = job->job->mmIndx.m[i];
		if (imol >= job->job->mdcell->mm.n) continue;
		mm = job->job->mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			bpc = (BASIC_CLUSTER*)(mm->cluster + ic);
			CLUSTER_LJ_INTERACT(*cmm, bpc, interact->var, tres);
			interact->LJ_res += tres;
			if (!bpc->eNeutral) {
				ExplicitEInteract(*cmm, bpc, interact->evar, tres);
				interact->eres += tres;
			}
		}
	}
	for (i = 0; i < job->job->smIndx.n; i++) {
		imol = job->job->smIndx.m[i];
		if (imol >= job->job->mdcell->sm.n) continue;
		bpc = job->job->mdcell->sm.m[imol].c;
		CLUSTER_LJ_INTERACT(*cmm, bpc, interact->var, tres);
		interact->LJ_res += tres;
		if (!bpc->eNeutral) {
			ExplicitEInteract(*cmm, bpc, interact->evar, tres);
			interact->eres += tres;
		}
	}
	for (i = 0; i < job->job->pmIndx.n; i++) {
		imol = job->job->pmIndx.m[i];
		if (imol >= job->job->mdcell->pm.n) continue;
		bpc = job->job->mdcell->pm.m[imol].c;
		CLUSTER_LJ_INTERACT(*cmm, bpc, interact->var, tres);
		interact->LJ_res += tres;
		if (!bpc->eNeutral) {
			ExplicitEInteract(*cmm, bpc, interact->evar, tres);
			interact->eres += tres;
		}
	}
}

void MDCELL_VirialCorrect_RigidCluster(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4>* virial) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	int nm = 0, imol = 0, ic, ia;
	double STxx = 0, STyy = 0, STzz = 0;
	SVECTOR<double, 3> dr;
	SVECTOR<double, 3> F;
	double *f, *r;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			if (pc->fatom.m == NULL) {
				for (ia = 0; ia < pc->nAtoms; ia++) {
					patom = pc->atom + ia;
					VECT3(pc->cm0, patom->r0, dr) r = dr.v;
					if (::bEext) { // subtract the external force
						F.v[0] = patom->F.v[0];
						F.v[1] = patom->F.v[1];
						F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
						f = F.v;
					}
					else f = patom->F.v;
					STxx += r[0] * f[0];
					STyy += r[1] * f[1];
					STzz += r[2] * f[2];
				}
			}
			else {
				for (ia = 0; ia < pc->fatom.n; ia++) {
					patom = pc->fatom.m[ia];
					VECT3(pc->cm, patom->r, dr) r = dr.v;
					if (::bEext) { // subtract the external force
						F.v[0] = patom->F.v[0];
						F.v[1] = patom->F.v[1];
						F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
						f = F.v;
					}
					else f = patom->F.v;
					STxx += r[0] * f[0];
					STyy += r[1] * f[1];
					STzz += r[2] * f[2];
				}
			}
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		if (sm->c->fatom.m == NULL) {
			for (ia = 0; ia < sm->c->nAtoms; ia++) {
				patom = sm->c->atom + ia;
				// SM has mass center at r of molecule, so r0 is the displacement from mass center
				r = patom->r0.v;
				if (::bEext) { // subtract the external force
					F.v[0] = patom->F.v[0];
					F.v[1] = patom->F.v[1];
					F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
					f = F.v;
				}
				else f = patom->F.v;
				STxx += r[0] * f[0];
				STyy += r[1] * f[1];
				STzz += r[2] * f[2];
			}
		}
		else {
			for (ia = 0; ia < sm->c->fatom.n; ia++) {
				patom = sm->c->fatom.m[ia];
				// SM has mass center at r of molecule, so r0 is the displacement from mass center
				r = patom->r0.v;
				if (::bEext) { // subtract the external force
					F.v[0] = patom->F.v[0];
					F.v[1] = patom->F.v[1];
					F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
					f = F.v;
				}
				else f = patom->F.v;
				STxx += r[0] * f[0];
				STyy += r[1] * f[1];
				STzz += r[2] * f[2];
			}
		}
	}

	// point molecule has only one atom, no correction is required for virial
	/*
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
	}
	*/
	//STxx *= 0.5; STyy *= 0.5; STzz *= 0.5;
	virial->v[0] = STxx; virial->v[1] = STyy; virial->v[2] = STzz;
	virial->v[3] = STxx + STyy + STzz;
}

void MDCELL_ClusterVirialCorrect_HingeConstraint(MDCELL_ThreadJob<MMOL_MD_CELL> *job, HingeConstraintVar &var, InteractRes &ires) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL, *parent;
	HINGE_CONSTRAINT *hc;
	int nm = 0, imol = 0, ic, ia;
	double STxx = 0, STyy = 0, STzz = 0;
	SVECTOR<double, 3> dr;
	double *f, *r;
	double *v1, *v2, *v3;

	ires.reset();
	InteractRes tres;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		NEIMO_HingeForce(*mm);
		if (!HingeConstraintForce(*mm, var, tres)) continue;
		ires += tres;
	}
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			if (pc->hinge_constraint == NULL) continue;
			hc = pc->hinge_constraint;

			v1 = pc->cm.v; v3 = dr.v; r = dr.v;
			// 00
			v2 = hc->patom_00->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc0[0].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];
			
			// 01
			v2 = hc->patom_01->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc0[1].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];

			// 02
			v2 = hc->patom_02->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc0[2].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];

			// parent cluster
			parent = (CLUSTER*)(pc->parent);
			v1 = parent->cm.v; v3 = dr.v; r = dr.v;
			// 10
			v2 = hc->patom_10->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc1[0].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];
			
			// 11
			v2 = hc->patom_11->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc1[1].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];

			// 12
			v2 = hc->patom_12->r.v; 
			DVECT3(v1, v2, v3) 
			f = hc->fc1[2].v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];
		}
	}
	ires.STxx -= STxx; 
	ires.STyy -= STyy;
	ires.STzz -= STzz;
	ires.STtr -= STxx + STyy + STzz;
}

void MDCELL_VirialCorrect_Molecule(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4>* virial) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	int nm = 0, imol = 0, ic, ia;
	double STxx = 0, STyy = 0, STzz = 0;
	SVECTOR<double, 3> dr;
	SVECTOR<double, 3> F;
	double *f, *r;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			for (ia = 0; ia < pc->nAtoms; ia++) {
				patom = pc->atom + ia;
				VECT3(mm->mDynKin.cm, patom->r, dr) r = dr.v;
				if (::bEext) { // subtract the external force
					F.v[0] = patom->F.v[0];
					F.v[1] = patom->F.v[1];
					F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
					f = F.v;
				}
				else f = patom->F.v;
				STxx += r[0] * f[0];
				STyy += r[1] * f[1];
				STzz += r[2] * f[2];
			}
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		for (ia = 0; ia < sm->c->nAtoms; ia++) {
			patom = sm->c->atom + ia;
			// SM has mass center at r of molecule, so r0 is the displacement from mass center
			r = patom->r0.v;
			if (::bEext) { // subtract the external force
				F.v[0] = patom->F.v[0];
				F.v[1] = patom->F.v[1];
				F.v[2] = patom->F.v[2] - patom->c0 * ::Eext * ::fUnit_ext_mass; 
				f = F.v;
			}
			else f = patom->F.v;
			STxx += r[0] * f[0];
			STyy += r[1] * f[1];
			STzz += r[2] * f[2];
		}
	}

	// point molecule has only one atom, no correction is required for virial
	/*
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
	}
	*/
	//STxx *= 0.5; STyy *= 0.5; STzz *= 0.5;
	virial->v[0] = STxx; virial->v[1] = STyy; virial->v[2] = STzz;
	virial->v[3] = STxx + STyy + STzz;
}

void gen_random_force(MMOL_MD_CELL &mdcell, _rand_::RAND_ARRAY<int> &randa_m, _rand_::RAND_ARRAY<int> &randa_c) {
	int imol, nc, i, i_m, i_c;
	SVECTOR<double, 3> f;
	BASIC_CLUSTER *pc = NULL;
	if (mdcell.mm.n > 0) {
		_rand_::init_random_array(randa_m, mdcell.mm.n);
		for (i_m = 0; i_m < randa_m.nw; i_m++) {
			imol = randa_m.a.m[i_m]; if (imol >= mdcell.mm.n) continue;
			_rand_::init_random_array(randa_c, mdcell.mm.m[imol].nCluster);
			for (i_c = 0; i_c < mdcell.mm.m[imol].nCluster; i_c++) {
				nc = randa_c.a.m[i_c]; if (nc >= mdcell.mm.m[imol].nCluster) continue;
				pc = mdcell.mm.m[imol].cluster + nc;
				gauss_random_force(randf, fwhm_randf, f.v, 1);
				for (i = 0; i < 3; i++) {
					pc->dyn->f.v[i+3] += f.v[i];
					pc->dyn->fc.v[i+3] += f.v[i];
				}
			}
		}
	}

	if (mdcell.sm.n > 0) {
		_rand_::init_random_array(randa_m, mdcell.sm.n);
		for (i_m = 0; i_m < randa_m.nw; i_m++) {
			imol = randa_m.a.m[i_m]; if (imol >= mdcell.sm.n) continue;
			pc = mdcell.sm.m[imol].c;
			gauss_random_force(randf, fwhm_randf, f.v, 1);
			for (i = 0; i < 3; i++) {
				pc->dyn->f.v[i+3] += f.v[i];
				pc->dyn->fc.v[i+3] += f.v[i];
			}
		}
	}

	if (mdcell.pm.n > 0) {
		_rand_::init_random_array(randa_m, mdcell.pm.n);
		for (i_m = 0; i_m < randa_m.nw; i_m++) {
			imol = randa_m.a.m[i_m]; if (imol >= mdcell.pm.n) continue;
			pc = mdcell.pm.m[imol].c;
			gauss_random_force(randf, fwhm_randf, f.v, 1);
			for (i = 0; i < 3; i++) {
				pc->dyn->fc.v[i+3] += f.v[i];
				pc->dyn->f.v[i+3] += f.v[i];
			}
		}
	}
}

extern void calTransSpatialMomentum(MMOLECULE &mm, VECTOR3 &P);
void MDCELL_SpatialMomentum(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0;
// do we need to check the rotation momentum ?
	VECTOR3 P, p;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		calTransSpatialMomentum(*mm, p); // velocity of each cluster at its mass-center
		V3plusV3(p, P, P)
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; M = sm->c->M;
		P.v[0] += sm->c->bkm->P0.v[3];
		P.v[1] += sm->c->bkm->P0.v[4];
		P.v[2] += sm->c->bkm->P0.v[5];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol; M = pm->c->M;
		P.v[0] += M * pm->v.v[0];
		P.v[1] += M * pm->v.v[1];  
		P.v[2] += M * pm->v.v[2];
	}
	mp->v[0] = P.v[0]; mp->v[1] = P.v[1]; mp->v[2] = P.v[2];
}

void MDCELL_ZeroSpatialMomentum(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *vp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	VECTOR3 *v = &(mdcell->Vmc);

	if (FABS(v->v[0]) < 1e-8 && FABS(v->v[1]) < 1e-8 && FABS(v->v[2]) < 1e-8) return;
	else show_log("zero the total spatial momentum of the cell", true);

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			mm->cluster[ic].bkm->V0.v[3] -= vp->v[0];
			mm->cluster[ic].bkm->V0.v[4] -= vp->v[1];
			mm->cluster[ic].bkm->V0.v[5] -= vp->v[2];

			mm->cluster[ic].bkm->P0.v[3] -= mm->cluster[ic].M * vp->v[0];
			mm->cluster[ic].bkm->P0.v[4] -= mm->cluster[ic].M * vp->v[1];
			mm->cluster[ic].bkm->P0.v[5] -= mm->cluster[ic].M * vp->v[2];

			mm->cluster[ic].mcDyn->v.v[0] -= vp->v[0];
			mm->cluster[ic].mcDyn->v.v[1] -= vp->v[1];
			mm->cluster[ic].mcDyn->v.v[2] -= vp->v[2];
		}
		mm->V0.v[3] -= vp->v[0];
		mm->V0.v[4] -= vp->v[1];
		mm->V0.v[5] -= vp->v[2];
		mm->mDynKin.V.v[3] -= vp->v[0];
		mm->mDynKin.V.v[4] -= vp->v[1];
		mm->mDynKin.V.v[5] -= vp->v[2];

		mm->mDynKin.P.v[3] -= mm->mDynKin.M * vp->v[0];
		mm->mDynKin.P.v[4] -= mm->mDynKin.M * vp->v[1];
		mm->mDynKin.P.v[5] -= mm->mDynKin.M * vp->v[2];
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		sm->c->bkm->V0.v[3] -= vp->v[0];
		sm->c->bkm->V0.v[4] -= vp->v[1];
		sm->c->bkm->V0.v[5] -= vp->v[2];
		
		sm->c->bkm->P0.v[3] -= sm->c->M * vp->v[0];
		sm->c->bkm->P0.v[4] -= sm->c->M * vp->v[1];
		sm->c->bkm->P0.v[5] -= sm->c->M * vp->v[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		pm->v.v[0] -= -vp->v[0];
		pm->v.v[1] -= vp->v[1];
		pm->v.v[2] -= vp->v[2];

		pm->c->bkm->P0.v[3] -= pm->c->M * vp->v[0];
		pm->c->bkm->P0.v[4] -= pm->c->M * vp->v[1];
		pm->c->bkm->P0.v[5] -= pm->c->M * vp->v[2];
	}
}

void ZeroSpatialMomentum(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob) {
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>,SVECTOR<double, 3> >((void*)(&MDCELL_SpatialMomentum), job, nJob, mp);
	VECTOR3 P, v, Pw, w;
	float M = job->job->mdcell->M;
	bool bMCt = job->job->mdcell->bHDynMC_T;
	int i;
	for (i = 0; i < nJob; i++) {
		P.v[0] += mp[i].v[0]; P.v[1] += mp[i].v[1]; P.v[2] += mp[i].v[2];
	}
	v.v[0] = P.v[0] / M; v.v[1] = P.v[1] / M; v.v[2] = P.v[2] / M;
	
	for (i = 0; i < nJob; i++) memcpy(mp[i].v, v.v, SIZE_V3);
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 3> >((void*)(&MDCELL_ZeroSpatialMomentum), job, nJob, mp);
}


void MDCELL_MassCenter(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;
	VECTOR3 Rcm;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
			pc = mm->cluster + ic;
			Rcm.v[0] += pc->M * pc->vp->v[0]; //(pc->cm.v[0] + pc->dr_cmm.v[0]);
			Rcm.v[1] += pc->M * pc->vp->v[1]; //(pc->cm.v[1] + pc->dr_cmm.v[1]);
			Rcm.v[2] += pc->M * pc->vp->v[2]; //(pc->cm.v[2] + pc->dr_cmm.v[2]);
		}
	}

	SMOLECULE *sm = NULL;
	float M = 0;
	double *r = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol; M = sm->c->M; r = sm->r.v;
		Rcm.v[0] += M * r[0];
		Rcm.v[1] += M * r[1];
		Rcm.v[2] += M * r[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol; M = pm->c->M; r = pm->r.v;
		Rcm.v[0] += M * r[0];
		Rcm.v[1] += M * r[1];
		Rcm.v[2] += M * r[2];
	}
	memcpy(mp->v, Rcm.v, SIZE_V3);
}

void MDCELL_ApplySpring(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp) {
	MMOL_MD_CELL *mdcell = job->job->mdcell;
	CLUSTER *pc = NULL;
	int nm = 0, imol = 0, ic;

	VECTOR3 f, t, dr;
	VECTOR<6> ft, fmc;
	MATRIX<6> *phai;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		for (ic = 0; ic < mm->nCluster; ic++) {
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

			wv2V(t, f, fmc)
			phai = &(pc->invNE->phai_cm2me);
			MpV6((*phai), fmc, ft)
			V6plusV6(mm->mDynKin.fc, ft, mm->mDynKin.fc)
		}
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		sm->c->dyn->fc.v[3] += sm->c->M * mp->v[0];
		sm->c->dyn->fc.v[4] += sm->c->M * mp->v[1];
		sm->c->dyn->fc.v[5] += sm->c->M * mp->v[2];
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		pm->c->dyn->fc.v[3] += pm->c->M * mp->v[0];
		pm->c->dyn->fc.v[4] += pm->c->M * mp->v[1];
		pm->c->dyn->fc.v[5] += pm->c->M * mp->v[2];
	}
}

void ApplySpring(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob, float *k) {
	int i, j;
	for (i = 0; i < nJob; i++) memset(mp[i].v, 0, SIZE_V3);
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 3> >((void*)(&MDCELL_MassCenter), job, nJob, mp);
	float M = job->job->mdcell->M;
	double Rz = 0, fz;	
	for (j = 0; j < 3; j++) {
		Rz = 0;
		for (i = 0; i < nJob; i++) Rz += mp[i].v[j];
		Rz /= M;
		job->job->mdcell->rc.v[j] = Rz; // update mass center value of cell
		fz = -k[0] * Rz + k[1] * (Rz > 0 ? -1 : 1) * Rz * Rz;
		for (i = 0; i < nJob; i++) mp[i].v[j] = fz; 
	}
	
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 3> >((void*)(&MDCELL_ApplySpring), job, nJob, mp);
}

void MDCellMassCenter(MDCELL_ThreadJob< MMOL_MD_CELL> *job, SVECTOR<double, 3> *mp, int nJob) {
	int i, j;
	for (i = 0; i < nJob; i++) memset(mp[i].v, 0, SIZE_V3);
	MTOperate1<MDCELL_ThreadJob<MMOL_MD_CELL>, SVECTOR<double, 3> >((void*)(&MDCELL_MassCenter), job, nJob, mp);
	float M = job->job->mdcell->M;
	double Rz = 0, fz;	
	for (j = 0; j < 3; j++) {
		Rz = 0;
		for (i = 0; i < nJob; i++) Rz += mp[i].v[j];
		Rz /= M;
		job->job->mdcell->rc.v[j] = Rz; // update mass center value of cell
	}
}

void InitMD(MMOL_MD_CELL &mmcell, bool bStorageChain) {
	mmcell.init_big_small_mm();
	mmcell.init_groups();

	::bVirial = true;

	int nc = 0, nm = 0;
	double Ek = 0;

	float vcluster = 1.5 * 1.5 * 1.5, rcut = 0;
	bool bRelshipStorageChain = false;
	rcut = ::rcut_LJ + ::size_cluster * 2;
	int ncluster_LJ = int(rcut * rcut * rcut / vcluster * 1.2 + 0.5);
	rcut = ::rcut_Ewd + ::size_cluster * 2;
	int ncluster_SPME = int(rcut * rcut * rcut / vcluster * 1.2 + 0.5);
	int ncluster_ImpSolv = int(::dSolv_max * ::dSolv_max * ::dSolv_max / vcluster * 1.2 + 0.5);

	mmcell.bHingeConstraint = true;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mmcell.mm.m[nm].setup_cluster_fatom();// be sure the eneutral is checked
		mmcell.mm.m[nm].setup_mcDyn();
		mmcell.mm.m[nm].setup_kin_dyn();
		mmcell.mm.m[nm].init_hinge_force();
		if (!mmcell.mm.m[nm].setup_hinge_constraint()) {
			mmcell.bHingeConstraint = false; show_log("disable hinge constraint", true);
		}
		mmcell.mm.m[nm].setup_hinge_dihedral();
		if (!bRelshipStorageChain) {
			for (nc = 0; nc < mmcell.mm.m[nm].nCluster; nc++) {
				mmcell.mm.m[nm].cluster[nc].LJRship.use_array(ncluster_LJ);
				mmcell.mm.m[nm].cluster[nc].spmeRship.use_array(ncluster_SPME);
				if (::bImplicitSolvent) mmcell.mm.m[nm].cluster[nc].impsolv_rship.use_array(ncluster_ImpSolv);
			}
		}
	}
	if (mmcell.sm.m != NULL) {
		for (nm = 0; nm < mmcell.sm.n; nm++) {
			mmcell.sm.m[nm].c->setup_cluster_fatom();// be sure the eneutral is checked
			mmcell.sm.m[nm].setup_kin_dyn();
			if (!bRelshipStorageChain) {
				mmcell.sm.m[nm].c->LJRship.use_array(ncluster_LJ);
				mmcell.sm.m[nm].c->spmeRship.use_array(ncluster_SPME);
				if (::bImplicitSolvent) mmcell.sm.m[nm].c->impsolv_rship.use_array(ncluster_ImpSolv);
			}
		}
	}
	if (mmcell.pm.m != NULL) {
		for (nm = 0; nm < mmcell.pm.n; nm++) {
			mmcell.pm.m[nm].c->setup_cluster_fatom();// be sure the eneutral is checked
			mmcell.pm.m[nm].setup_dyn();
			if (!bRelshipStorageChain) {
				mmcell.pm.m[nm].c->LJRship.use_array(ncluster_LJ);
				mmcell.pm.m[nm].c->spmeRship.use_array(ncluster_SPME);
				if (::bImplicitSolvent) mmcell.pm.m[nm].c->impsolv_rship.use_array(ncluster_ImpSolv);
			}
		}
	}

	check_atom_rg(mmcell, MAX_THREADS);
	// since molecule mass center position could be shifted, distribute cluster in CMM again 
	//for (nm = 0; nm < mmcell.mm.n; nm++) CalMolInertiaMoment(mmcell.mm.m[nm]); // here we will calculate mass center
	//for (nm = 0; nm < mmcell.sm.n; nm++) mmcell.sm.m[nm].calInertiaMoment();

	mmcell.init_mdcell_velocity_memory();

	// ************ SETUP the KinDyn variables *******************
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mm = mmcell.mm.m + nm;
		InitVelocity(*mm);
		mm->calClusterSolvationRadius(); // it uses the cluster's mass center position
		//mm->resetDynPars(); // set force & torque to 0
	}

	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		InitVelocity(*sm);
		sm->c->calSolvationRadius();
	}

	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m + nm;
		InitVelocity(*pm);
		pm->c->calSolvationRadius();
	}

	CHAIN<short> *ch_cluster_id = NULL;
	if (::bBrownian) {
		strcpy(::errmsg, "\n"); show_log(::errmsg, true);
		sprintf(::errmsg, "----------------------  Brownian Force  ---------------------"); show_log(::errmsg, true);
		sprintf(::errmsg, "  general parameters :"); show_log(::errmsg, true);
		sprintf(::errmsg, "  ita / radiu = %f", ::ita_Brownian); show_log(::errmsg, true);
		sprintf(::errmsg, "  cd_trans / d^2 = %f", ::cd_drag_trans); show_log(::errmsg, true);
		sprintf(::errmsg, "  cd_rot / d^2 = %f", ::cd_drag_rot); show_log(::errmsg, true);
		show_log("", true);
		sprintf(::errmsg, "Eistein-type diffusion coefficient D = kT / ita of cluster"); show_log(::errmsg, true);
		sprintf(::errmsg, "  For each cluster, ita = ita * radiu"); show_log(::errmsg, true);
		sprintf(::errmsg, "                    cd_trans = cd_trans * d^2"); show_log(::errmsg, true);
		sprintf(::errmsg, "                    cd_rot = cd_rot * d^2"); show_log(::errmsg, true);
		sprintf(::errmsg, ""); show_log(::errmsg, true);

		sprintf(::errmsg, "  CLUSTER UID: ita, cd_translation, cd_rotation"); show_log(::errmsg, true);
		sprintf(::errmsg, "  Brownian force: f = -ita * V - cd_trans/rot * |V| * V"); show_log(::errmsg, true);
		sprintf(::errmsg, "  Brownian force unit: mass * Angs/fs^2"); show_log(::errmsg, true);
		sprintf(::errmsg, "--------------------------------------------------------"); show_log(::errmsg, true);
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			mm = mmcell.mm.m + nm;
			mm->show_cluster_BrownianForcePars(&ch_cluster_id);
		}
		for (nm = 0; nm < mmcell.sm.n; nm++) {
			sm = mmcell.sm.m + nm;
			sm->c->show_BrownianForcePars(&ch_cluster_id);
		}
		for (nm = 0; nm < mmcell.pm.n; nm++) {
			pm = mmcell.pm.m + nm;
			pm->c->show_BrownianForcePars(&ch_cluster_id);
		}
		release_chain<short>(&ch_cluster_id, true);
		sprintf(::errmsg, "------------------------- OVER -------------------------"); show_log(::errmsg, true);
		show_log("", true);
	}

	/*
	show_log("-------  NHC ---------", true);
	sprintf(errmsg, "CELL translation : dt -- %f, tao -- %f", mmcell.mcHDyn.nhc->dt, 1 / sqrt(mmcell.mcHDyn.nhc->inv_taos)); 
	show_log(errmsg, true);
	sprintf(errmsg, "CELL rotation : dt -- %f, tao -- %f", mmcell.mcrHDyn.nhc->dt, 1 / sqrt(mmcell.mcrHDyn.nhc->inv_taos)); 
	show_log(errmsg, true);
	for (nc = 0; nc < mmcell.hdyn.n; nc++) {
		sprintf(errmsg, "# %d : dt -- %f, tao -- %f", nc, mmcell.hdyn.m[nc].nhc->dt, 1 / sqrt(mmcell.hdyn.m[nc].nhc->inv_taos));
		show_log(errmsg, true);
	}
	for (nm = 0; nm < mmcell.lfpmd.n; nm++) {
		sprintf(errmsg, "LFPMD # %d : dt = %f", nm, mmcell.lfpmd.m[nm].dt); show_log(errmsg, true);
	}
	show_log("-------- OVER --------", true);
	show_log("", true);
	*/
}

void InitMD(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, bool bStorageChain) {
	mmcell.init_big_small_mm();
	mmcell.init_groups();

	::bVirial = true;

	int nc = 0, nm = 0;
	double Ek = 0;

	float vcluster = 1.5 * 1.5 * 1.5, rcut = 0;
	bool bRelshipStorageChain = false;
	rcut = ::rcut_LJ + ::size_cluster * 2;
	int ncluster_LJ = int(rcut * rcut * rcut / vcluster * 1.2 + 0.5);
	if (ncluster_LJ > 20) ncluster_LJ = 20;
	rcut = ::rcut_Ewd + ::size_cluster * 2;
	int ncluster_SPME = int(rcut * rcut * rcut / vcluster * 1.2 + 0.5);
	if (ncluster_SPME > 20) ncluster_SPME = 20;
	int ncluster_ImpSolv = int(::dSolv_max * ::dSolv_max * ::dSolv_max / vcluster * 1.2 + 0.5);
	if (ncluster_ImpSolv > 20) ncluster_ImpSolv = 20;

	mmcell.bHingeConstraint = true;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mmcell.mm.m[nm].setup_cluster_fatom();// be sure the eneutral is checked
		mmcell.mm.m[nm].setup_mcDyn();
		mmcell.mm.m[nm].setup_kin_dyn();
		mmcell.mm.m[nm].init_hinge_force();
		if (!mmcell.mm.m[nm].setup_hinge_constraint()) {
			mmcell.bHingeConstraint = false; show_log("disable hinge constraint", true);
		}
		mmcell.mm.m[nm].setup_hinge_dihedral();
		if (!bRelshipStorageChain) {
			for (nc = 0; nc < mmcell.mm.m[nm].nCluster; nc++) {
				mmcell.mm.m[nm].cluster[nc].LJRship.use_array(ncluster_LJ);
				mmcell.mm.m[nm].cluster[nc].spmeRship.use_array(ncluster_SPME);
				if (::bImplicitSolvent) mmcell.mm.m[nm].cluster[nc].impsolv_rship.use_array(ncluster_ImpSolv);
			}
		}
	}
	if (mmcell.sm.m != NULL) {
		for (nm = 0; nm < mmcell.sm.n; nm++) {
			mmcell.sm.m[nm].c->setup_cluster_fatom();// be sure the eneutral is checked
			mmcell.sm.m[nm].setup_kin_dyn();
			if (!bRelshipStorageChain) {
				mmcell.sm.m[nm].c->LJRship.use_array(ncluster_LJ);
				mmcell.sm.m[nm].c->spmeRship.use_array(ncluster_SPME);
				if (::bImplicitSolvent) mmcell.sm.m[nm].c->impsolv_rship.use_array(ncluster_ImpSolv);
			}
		}
	}
	if (mmcell.pm.m != NULL) {
		for (nm = 0; nm < mmcell.pm.n; nm++) {
			mmcell.pm.m[nm].c->setup_cluster_fatom();// be sure the eneutral is checked
			mmcell.pm.m[nm].setup_dyn();
			if (!bRelshipStorageChain) {
				mmcell.pm.m[nm].c->LJRship.use_array(ncluster_LJ);
				mmcell.pm.m[nm].c->spmeRship.use_array(ncluster_SPME);
				if (::bImplicitSolvent) mmcell.pm.m[nm].c->impsolv_rship.use_array(ncluster_ImpSolv);
			}
		}
	}

	check_atom_rg(mmcell, MAX_THREADS);
	// since molecule mass center position could be shifted, distribute cluster in CMM again 
	//for (nm = 0; nm < mmcell.mm.n; nm++) CalMolInertiaMoment(mmcell.mm.m[nm]); // here we will calculate mass center
	//for (nm = 0; nm < mmcell.sm.n; nm++) mmcell.sm.m[nm].calInertiaMoment();

	bool bCMM_Storage_Chain = bStorageChain;
	bool reset_cell_chain = true;
	if (bCMM_Storage_Chain) {
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			reset_cell_chain = (nm == 0 ? true : false);
			cmm_init_distribute_cluster(mmcell.mm.m[nm], *cmm, reset_cell_chain); // set base cluster at the center cell
		}
		reset_cell_chain = (mmcell.mm.n == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.sm.m, mmcell.sm.n, *cmm, reset_cell_chain); // add simple solid molecule into cmm
		reset_cell_chain = (mmcell.mm.n + mmcell.sm.n == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.pm.m, mmcell.pm.n, *cmm, reset_cell_chain); // add point molecule into cmm

		//cmm_init_subcell<BASIC_CLUSTER>(*cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size
	}
	else {
		{
		float vcluster = 1.5 * 1.5 * 1.5;
		int ncs = int(cmm->fw[0] * cmm->fw[1] * cmm->fw[2] / vcluster * 1.2 + 0.5);
		if (ncs < 5) ncs = 5;
		else if (ncs > 20) ncs = 20;
		CMM_array_set_storage<BASIC_CLUSTER>(cmm, ncs);
		}

		for (nm = 0; nm < mmcell.mm.n; nm++) {
			reset_cell_chain = (nm == 0 ? true : false);
			cmm_init_distribute_cluster(mmcell.mm.m[nm], *cmm, reset_cell_chain); // set base cluster at the center cell
		}
		reset_cell_chain = (mmcell.mm.n == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.sm.m, mmcell.sm.n, *cmm, reset_cell_chain); // add simple solid molecule into cmm
		reset_cell_chain = (mmcell.mm.n + mmcell.sm.n == 0 ? true : false);
		cmm_init_distribute_cluster(mmcell.pm.m, mmcell.pm.n, *cmm, reset_cell_chain); // add point molecule into cmm
	}

	mmcell.init_mdcell_velocity_memory();

	// ************ SETUP the KinDyn variables *******************
	MMOLECULE *mm = NULL;
	CLUSTER *pc = NULL;
	SMOLECULE *sm = NULL;
	PMOLECULE *pm = NULL;
	for (nm = 0; nm < mmcell.mm.n; nm++) {
		mm = mmcell.mm.m + nm;
		InitVelocity(*mm);
		mm->calClusterSolvationRadius(); // it uses the cluster's mass center position
		//mm->resetDynPars(); // set force & torque to 0
	}

	for (nm = 0; nm < mmcell.sm.n; nm++) {
		sm = mmcell.sm.m + nm;
		InitVelocity(*sm);
		sm->c->calSolvationRadius();
	}

	for (nm = 0; nm < mmcell.pm.n; nm++) {
		pm = mmcell.pm.m + nm;
		InitVelocity(*pm);
		pm->c->calSolvationRadius();
	}

	CHAIN<short> *ch_cluster_id = NULL;
	if (::bBrownian) {
		strcpy(::errmsg, "\n"); show_log(::errmsg, true);
		sprintf(::errmsg, "----------------------  Brownian Force  ---------------------"); show_log(::errmsg, true);
		sprintf(::errmsg, "  general parameters :"); show_log(::errmsg, true);
		sprintf(::errmsg, "  ita / radiu = %f", ::ita_Brownian); show_log(::errmsg, true);
		sprintf(::errmsg, "  cd_trans / d^2 = %f", ::cd_drag_trans); show_log(::errmsg, true);
		sprintf(::errmsg, "  cd_rot / d^2 = %f", ::cd_drag_rot); show_log(::errmsg, true);
		show_log("", true);
		sprintf(::errmsg, "Eistein-type diffusion coefficient D = kT / ita of cluster"); show_log(::errmsg, true);
		sprintf(::errmsg, "  For each cluster, ita = ita * radiu"); show_log(::errmsg, true);
		sprintf(::errmsg, "                    cd_trans = cd_trans * d^2"); show_log(::errmsg, true);
		sprintf(::errmsg, "                    cd_rot = cd_rot * d^2"); show_log(::errmsg, true);
		sprintf(::errmsg, ""); show_log(::errmsg, true);

		sprintf(::errmsg, "  CLUSTER UID: ita, cd_translation, cd_rotation"); show_log(::errmsg, true);
		sprintf(::errmsg, "  Brownian force: f = -ita * V - cd_trans/rot * |V| * V"); show_log(::errmsg, true);
		sprintf(::errmsg, "  Brownian force unit: mass * Angs/fs^2"); show_log(::errmsg, true);
		sprintf(::errmsg, "--------------------------------------------------------"); show_log(::errmsg, true);
		for (nm = 0; nm < mmcell.mm.n; nm++) {
			mm = mmcell.mm.m + nm;
			mm->show_cluster_BrownianForcePars(&ch_cluster_id);
		}
		for (nm = 0; nm < mmcell.sm.n; nm++) {
			sm = mmcell.sm.m + nm;
			sm->c->show_BrownianForcePars(&ch_cluster_id);
		}
		for (nm = 0; nm < mmcell.pm.n; nm++) {
			pm = mmcell.pm.m + nm;
			pm->c->show_BrownianForcePars(&ch_cluster_id);
		}
		release_chain<short>(&ch_cluster_id, true);
		sprintf(::errmsg, "------------------------- OVER -------------------------"); show_log(::errmsg, true);
		show_log("", true);
	}

	/*
	show_log("-------  NHC ---------", true);
	sprintf(errmsg, "CELL translation : dt -- %f, tao -- %f", mmcell.mcHDyn.nhc->dt, 1 / sqrt(mmcell.mcHDyn.nhc->inv_taos)); 
	show_log(errmsg, true);
	sprintf(errmsg, "CELL rotation : dt -- %f, tao -- %f", mmcell.mcrHDyn.nhc->dt, 1 / sqrt(mmcell.mcrHDyn.nhc->inv_taos)); 
	show_log(errmsg, true);
	for (nc = 0; nc < mmcell.hdyn.n; nc++) {
		sprintf(errmsg, "# %d : dt -- %f, tao -- %f", nc, mmcell.hdyn.m[nc].nhc->dt, 1 / sqrt(mmcell.hdyn.m[nc].nhc->inv_taos));
		show_log(errmsg, true);
	}
	for (nm = 0; nm < mmcell.lfpmd.n; nm++) {
		sprintf(errmsg, "LFPMD # %d : dt = %f", nm, mmcell.lfpmd.m[nm].dt); show_log(errmsg, true);
	}
	show_log("-------- OVER --------", true);
	show_log("", true);
	*/
}

void MMOL_MD_CELL::init_dipole_hist() {
	int ia, na = 0;
	if (atom.m == NULL) {
		show_log("error: initialize atomic induced dipole history after InitMD function is called!", true);
		return;
	}
	pol_atom.release();
	for (ia = 0; ia < atom.n; ia++) {
		if (atom.m[ia]->par->alpha > 1e-4) na++;
	}
	if (na == 0) return;
	pol_atom.set_array(na);
	mu_hist.set_array(na);
	na = 0;
	for (ia = 0; ia < atom.n; ia++) {
		if (atom.m[ia]->par->alpha > 1e-4) {
			pol_atom.m[na] = atom.m[ia];
			atom.m[ia]->mu_hist = mu_hist.m + na;
			na++;
		}
		else atom.m[ia]->mu_hist = NULL;
	}
};

void Backup_atomic_induced_dipole_hist(MATOM** atom, int n) {
	int ia;
	DIPOLE_HIST *mu_hist = NULL;
	double *mu_ind = NULL;
	for (ia = 0; ia < n; ia++) {
		mu_hist = atom[ia]->mu_hist;
		mu_ind = atom[ia]->ind_mu.v;
		memcpy(mu_hist->mu[2].v, mu_hist->mu[1].v, SIZE_V3);
		memcpy(mu_hist->mu[1].v, mu_hist->mu[0].v, SIZE_V3);
		memcpy(mu_hist->mu[0].v, mu_ind, SIZE_V3);
	}
}

void Backup_induced_dipole_hist(MMOL_MD_CELL &mdcell, int nMDThreads) {
	if (nMDThreads <= 1 || mdcell.pol_atom.n < 100) Backup_atomic_induced_dipole_hist(mdcell.pol_atom.m, mdcell.pol_atom.n);
	else MOperate<MATOM*>((void*)(&Backup_atomic_induced_dipole_hist), mdcell.pol_atom.m, mdcell.pol_atom.n, nMDThreads);
};

void Guess_atomic_induced_dipole(MATOM** atom, int n) {
	// guessing induced dipole of a new step, from the last 3 history
	int ia;
	DIPOLE_HIST *mu_hist = NULL;
	double *mu_ind = NULL;
	for (ia = 0; ia < n; ia++) {
		mu_hist = atom[ia]->mu_hist;
		mu_ind = atom[ia]->ind_mu.v;
		mu_ind[0] = (mu_hist->mu[0].v[0] - mu_hist->mu[1].v[0]) * 3 + mu_hist->mu[2].v[0];
		mu_ind[1] = (mu_hist->mu[0].v[1] - mu_hist->mu[1].v[1]) * 3 + mu_hist->mu[2].v[1];
		mu_ind[2] = (mu_hist->mu[0].v[2] - mu_hist->mu[1].v[2]) * 3 + mu_hist->mu[2].v[2];

		//mu_ind[0] = 2 * mu_hist->mu[0].v[0] - mu_hist->mu[1].v[0];
		//mu_ind[1] = 2 * mu_hist->mu[0].v[1] - mu_hist->mu[1].v[1];
		//mu_ind[2] = 2 * mu_hist->mu[0].v[2] - mu_hist->mu[1].v[2];

		if (atom[ia]->bIMu) {
			atom[ia]->mu.v[0] = atom[ia]->intr_mu.v[0] + mu_ind[0];
			atom[ia]->mu.v[1] = atom[ia]->intr_mu.v[1] + mu_ind[1];
			atom[ia]->mu.v[2] = atom[ia]->intr_mu.v[2] + mu_ind[2];
		}
		else memcpy(atom[ia]->mu.v, mu_ind, SIZE_V3);
	}
}

void Guess_induced_dipole(MMOL_MD_CELL &mdcell, int nMDThreads) {
	if (nMDThreads <= 1 || mdcell.pol_atom.n < 100) Guess_atomic_induced_dipole(mdcell.pol_atom.m, mdcell.pol_atom.n);
	else MOperate<MATOM*>((void*)(&Guess_atomic_induced_dipole), mdcell.pol_atom.m, mdcell.pol_atom.n, nMDThreads);
};

void MDCELL_CalMolTransKineticEnergy(MDCELL_ThreadJob<MMOL_MD_CELL> *job, SVECTOR<double, 4> *Ekt) {
	extern void NPT_calKineticEnergy_Cluster(MMOLECULE& mm, bool bTransOnly = false);
	extern void NPT_calKineticEnergy_Mol(MMOLECULE& mm, bool bTransOnly = false);
	extern void NPT_calKineticEnergy(SMOLECULE& sm, bool bTransOnly = false);
	extern void NPT_calKineticEnergy(PMOLECULE& pm);

	MMOL_MD_CELL *mdcell = job->job->mdcell;
	int nm = 0, imol = 0;
	double Ek = 0, Ekx = 0, Eky = 0, Ekz = 0;
	VECTOR3 P;

	MMOLECULE *mm = NULL;
	for (nm = 0; nm < job->job->mmIndx.n; nm++) {
		imol = job->job->mmIndx.m[nm]; if (imol >= mdcell->mm.n) continue;
		mm = mdcell->mm.m + imol;
		calTransSpatialMomentum(*mm, P); // velocity of each cluster at its mass-center
		if (::method_Virial == CLUSTER_VIRIAL) {
			NPT_calKineticEnergy_Cluster(*mm, true); // translation kinetic energy of all clusters
		}
		else if (::method_Virial == MOLECULE_VIRIAL) {
			NPT_calKineticEnergy_Mol(*mm, true); // translation kinetic energy of the mass center of the molecule
		}
		Ek += mm->Ekt;
		Ekx += mm->Ektx;
		Eky += mm->Ekty;
		Ekz += mm->Ektz;
	}

	SMOLECULE *sm = NULL;
	for (nm = 0; nm < job->job->smIndx.n; nm++) {
		imol = job->job->smIndx.m[nm]; if (imol >= mdcell->sm.n) continue;
		sm = mdcell->sm.m + imol;
		NPT_calKineticEnergy(*sm, true);
		Ek += sm->Ekt;
		Ekx += sm->Ektx;
		Eky += sm->Ekty;
		Ekz += sm->Ektz;
	}

	PMOLECULE *pm = NULL;
	for (nm = 0; nm < job->job->pmIndx.n; nm++) {
		imol = job->job->pmIndx.m[nm]; if (imol >= mdcell->pm.n) continue;
		pm = mdcell->pm.m + imol;
		NPT_calKineticEnergy(*pm);
		Ek += pm->Ekt;
		Ekx += pm->Ektx;
		Eky += pm->Ekty;
		Ekz += pm->Ektz;
	}
	Ekt->v[0] = Ekx;
	Ekt->v[1] = Eky;
	Ekt->v[2] = Ekz;
	Ekt->v[3] = Ek;
}

void log_dipole(MMOL_MD_CELL *mdcell) {
	char msg[256] = "\0", buffer[256] = "\0";
	VECTOR3 mu, mut;
	BASIC_CLUSTER *pc;
	int im, ia;
	for (im = 0; im < mdcell->sm.n; im++) {
		pc = mdcell->sm.m[im].c; V3zero(mu) V3zero(mut)
		sprintf(msg, "sm %d : ", im);
		for (ia = 0; ia < pc->nAtoms; ia++) {
			mu += pc->atom[ia].mu;
			mut.v[0] += pc->atom[ia].c * pc->atom[ia].r0.v[0];
			mut.v[1] += pc->atom[ia].c * pc->atom[ia].r0.v[1];
			mut.v[2] += pc->atom[ia].c * pc->atom[ia].r0.v[2];
			sprintf(buffer, "%f, ", pc->atom[ia].mu.abs());
			strcat(msg, buffer);
		}
		mut += mu;
		sprintf(buffer, "induced dipole %f eA, total dipole %f eA", mu.abs(), mut.abs());
		strcat(msg, buffer);
		show_log(msg, true);
	}
	for (im = 0; im < mdcell->pm.n; im++) {
		pc = mdcell->pm.m[im].c; V3zero(mu) V3zero(mut)
		sprintf(msg, "pm %d : ", im);
		for (ia = 0; ia < pc->nAtoms; ia++) {
			mu += pc->atom[ia].mu;
			mut.v[0] += pc->atom[ia].c * pc->atom[ia].r0.v[0];
			mut.v[1] += pc->atom[ia].c * pc->atom[ia].r0.v[1];
			mut.v[2] += pc->atom[ia].c * pc->atom[ia].r0.v[2];
			sprintf(buffer, "%f, ", pc->atom[ia].mu.abs());
			strcat(msg, buffer);
		}
		mut += mu;
		sprintf(buffer, "induced dipole %f eA, total dipole %f eA", mu.abs(), mut.abs());
		strcat(msg, buffer);
		show_log(msg, true);
	}
	show_log("", true);
}

// Interface constraint
void InterfaceConstraints_cmm_cell(CMM_CELL3D<BASIC_CLUSTER> &cmm, int ix1, int ix2, double &Ut) {
	int ix, iy, iz, ic;
	BASIC_CLUSTER *pc;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > iterator;

	IConstraint *p_ic = NULL;

	Ut = 0;
	double U, f0, z;
	VECTOR3 dr, f, torq;

	CLUSTER *pcluster = NULL;

	for (ix = ix1; ix <= ix2; ix++) {
		if (ix < 0) continue;
		else if (ix >= cmm.nx) break;
	for (iy = 0; iy < cmm.bcell.ny; iy++) { for (iz = 0; iz < cmm.bcell.nz; iz++) {
		pcell = &(cmm.bcell.m[ix].m[iy].m[iz]); iterator.set(pcell); iterator.start_iterate();
		while (!iterator.eoi()) {
			pc = iterator.current();
			if (pc == NULL || pc->dyn == NULL) {iterator.next_iterate(); continue;}
			z = pc->vp->v[2];
			f.v[2] = 0;
			for (ic = 0; ic < _interface_constraint_::clusterIC_db.m[pc->cTypeIndx].n; ic++) {
				p_ic = _interface_constraint_::clusterIC_db.m[pc->cTypeIndx].m + ic;
				if (p_ic->zt < 0 && z > p_ic->zt) continue;
				else if (p_ic->zf > 0 && z < p_ic->zf) continue;
				if (p_ic->ef_prof.n <= 0) continue;
				p_ic->ef(z, U, f0);
				f.v[2] += f0;
				Ut += U;
			}
			pc->dyn->fc.v[5] += f.v[2]; 
			if (pc->parent != NULL) {
				pcluster = (CLUSTER*)pc;
				dr.v[0] = pcluster->cm0.v[0]; // cm0 is center of cluster relative to its Op
				dr.v[1] = pcluster->cm0.v[0];
				dr.v[2] = pcluster->cm0.v[0];
				V3PV3(dr, f, torq)
				pc->dyn->fc.v[0] += torq.v[0];
				pc->dyn->fc.v[1] += torq.v[1];
				pc->dyn->fc.v[2] += torq.v[2];
			}
			iterator.next_iterate();
		}
	}}
	}
}

void InterfaceConstraints_CMM(CMM_CELL3D<BASIC_CLUSTER> *cmm, int nThreads, double &Ut) {
	Ut = 0;
	if (nThreads > cmm->nx) nThreads = cmm->nx;
	CellOperate2<BASIC_CLUSTER, double>((void*)(&InterfaceConstraints_cmm_cell), *cmm, 0, cmm->nx - 1, Ut, nThreads);
}

// end of interface constraint
extern char parfname_IConstraint[256];

extern void CompensateNetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 3> *sv3);
extern void CompensateNetForceZ(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob, SVECTOR<double, 3> *sv3);
void NetForce(MDCELL_ThreadJob< MMOL_MD_CELL> *job, int nJob);

void MD_PROC(MMOL_MD_CELL &mmcell, CMM_CELL3D<BASIC_CLUSTER> *cmm, int md_mode) {
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
	/*
	if (mmcell.bHDynMC_T) {
		sprintf(errmsg, "although constraint on translation along z axis is enabled in MD.INI, disabled it!"); show_log(errmsg, true);
		mmcell.bHDynMC_T = false;
	}
	*/

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
	if (mmcell.bHDynMC_T) {
		MDCELL_mass_center(mmcell);
		ZeroSpatialMomentum(mdvar.job_mdcell, sv3job, MAX_THREADS);
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

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
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

		// calculate the inertia moment of each cluster and the whole molecule
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_CalInertiaMoment), mdvar.job_mdcell, nMDThreads);
		
		MDCellMassCenter(mdvar.job_mdcell, sv3job, nMDThreads);

		//CellOperate<BASIC_CLUSTER>((void*)(&InitClusterForceInCMM), *cmm, 0, cmm->acell.n - 1, MAX_THREADS);
		MTOperate< MDCELL_ThreadJob<MMOL_MD_CELL> >((void*)(&MDCELL_InitClusterForce), mdvar.job_mdcell, nMDThreads);

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
					mdvar.job_mdcellCMM[0].cmm = cmm; // directly allocate to global cmm
					cmm->reset_cell_chain();
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

			if (bCheckClusterInCMM) {
			{
				nc = cmm->number_atoms();
				/*
				for (i = 0; i < cmm->acell.n; i++) {
					nc += cmm->acell.m[i]->length();
					//sprintf(errmsg, "%d %d", i, cmm->acell.m[i]->length()); show_log(errmsg, true);
				}
				*/
				if (nc != ncs_mmcell) {
					sprintf(errmsg, "%d : cluster in CMM -- %d, real -- %d", nloop, nc, ncs_mmcell); show_log(errmsg, true);
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
				MTOperate<MDCELL_CMM_ThreadJob<MMOL_MD_CELL, CMM_CELL3D<BASIC_CLUSTER> > >(
					(void*)(&MMCELL_CMM_CheckRelship), mdvar.job_mdcellCMM, nMDThreads);
			}
		#if _CHECK_TIME_STAMP_ == 1
			sprintf(::errmsg, "Relationship_check takes %d ms", time.glance());
			show_infor(errmsg);
		#endif
		}

		// important : calculate LJ interaction first, here cluster force & torque are reset
		
		LJ_Interact(&mmcell, cmm, MAX_THREADS, mdvar.MDVar_LJ::LJ_var, mdvar.MDVar_LJ::LJ_res);
		U_LJ = mdvar.MDVar_LJ::LJ_res.U;

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		if (!mmcell.eNeutral && ::bRealEwaldSumOnly) {
			if (::bPolarize || ::bDipole) _spme_::LocalEF(&mmcell, cmm, mdvar.vspme_mu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			else _spme_::LocalEF(&mmcell, cmm, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
		}
		else if (!::local_interact && !mmcell.eNeutral) {
			mdvar.spme_var.bEfieldOnly = false;
			// important : calculate Coulomb secondly, here cluster force & torque are accumulated
			if (::bPolarize) {
				if (nloop > 3) Guess_induced_dipole(mmcell, MAX_THREADS);
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				if (::bExcludeInducedDipole) {
					if (!Polarize_SPME_Interact2(&mmcell, cmm, mdvar.vspme_mu, mdvar.vspme_mu_induced, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
						sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
					}
				}
				else if (!Polarize_SPME_Interact(&mmcell, cmm, mdvar.vspme_mu, MAX_THREADS, true, mdvar.spme_var, mdvar.spme_res)) {
					sprintf(errmsg, "polarized dipole is not convergent to 0.01D @ loop %d", nloop); show_log(errmsg, true);
				}
				Backup_induced_dipole_hist(mmcell, MAX_THREADS);
			}
			else if (::bDipole) {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				SPME_Interact(&mmcell, cmm, mdvar.vspme_mu, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			}
			else {
				// calculat the total dipole of the whole cell
				cal_dipole(mmcell, MAX_THREADS, mmcell.mu_surf); memcpy(&(cmm->mu_surf), &(mmcell.mu_surf), sizeof(_EwaldSum_real_::SURF_DIPOLE));
				SPME_Interact(&mmcell, cmm, mdvar.vspme_q, MAX_THREADS, mdvar.spme_var, mdvar.spme_res);
			}

			U_estat = mdvar.spme_res.U;
		}
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

		// before scale cluster, momentum of cluster has to be updated, because when atom are shocked, momentum is to be used to scale the interaction !
		// 0 -- use the spatial moment at mass center, 1 -- use the velocity of base cluster
		{
			for (int im = 0; im < mmcell.mm.n; im++) {
				if (::MM_SpeedVerlet_method == 0) MakeFreeBaseVelocityConsistent(mmcell.mm.m[im]); // adjusting free base velocity
				else if (::MM_SpeedVerlet_method == 1) {
					calClusterVelocity0(mmcell.mm.m[im]); calSpatialMomentum0(mmcell.mm.m[im], true); // required for Nose-Hover force calculation of the frame
				}
			}
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
				mdvar.virial.v[0] += sv3job[i].v[0];
				mdvar.virial.v[1] += sv3job[i].v[1];
				mdvar.virial.v[2] += sv3job[i].v[2];
				mdvar.virial.v[3] += (sv3job[i].v[0] + sv3job[i].v[1] + sv3job[i].v[2]);
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
				(*virial_save.out)<<"          "<<mdvar.virial.v[3] * ::unit_Ek_kT<<"  "<<Ek_t<<"  "<<mdvar.Pint<<"      "<<mmcell.rc.v[2]<<endl;
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

int nAtoms(MMOL_MD_CELL &mdcell) {
	int na = 0;
	int im, ic;
	for (im = 0; im < mdcell.mm.n; im++) {
		for (ic = 0; ic < mdcell.mm.m[im].nCluster; ic++) na += mdcell.mm.m[im].cluster[ic].nAtoms;
	}
	for (im = 0; im < mdcell.sm.n; im++) na += mdcell.sm.m[im].c->nAtoms;
	for (im = 0; im < mdcell.pm.n; im++) na += 1;
	return na;
}


