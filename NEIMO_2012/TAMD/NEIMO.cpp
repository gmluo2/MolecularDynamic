#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
#include "show.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "ZMatrix.h"

#include "MM.h"
#include "Mol.h"
#include "cluster.h"
#include "var.h"
#include "NEIMO.h"

extern bool trace;

#define _LINEAR_POLYMER  0
#if _LINEAR_POLYMER == 1
// given inertial cordinate of clusters
// Jacabino matrix, inertia mass matrix, inertia spatial matrix
// torque, a, b for each cluster
// to calculate the acceleration of torque angle 
void NEIMO_LINEAR_POLYMER(MMOLECULE &mm) {
	int i = 0, j = 0, k = 0;

	CLUSTER *p = mm.base;
	CLUSTER *child = NULL;

	MATRIX<6> Pp, tao, m1, m2, phi, phi_T;
	VECTOR<6> v1, v2;

	p = mm.tip[0]; // from tip to base
	while (p != NULL) {
		CHILD_CLUSTER(child, p, 0, i, j) // this is linear polymer
		if (child != NULL) {
			// if P = phi * P * phi_T, disable following two lines 
			MpM(child->km.phai, Pp, m1, i, j, k, 6) // m1 = phai(k, k - 1) * P+(k - 1);
			MpM(m1, child->km.phai_T, p->invNE->P, i, j, k, 6) // m2 = m1 * phai*(k, k - 1);
		}
		//p->invNE->P = m2 + p->invNE->M;
		for (i = 0; i < 6; i++) for (j = 0; j < 6; j++) 
			p->invNE->P.m[i][j] += p->invNE->M.m[i][j];

		MpV6(p->invNE->P, p->km.H, p->invNE->G) // v1 = P * H
		scale_uv6(p->km.H, p->invNE->G, p->invNE->D) // D = H * v1
		DIV(p->invNE->G, p->invNE->D, p->invNE->G, i)

		//MpV(p->km.phai, p->invNE->G, p->invNE->K, i, j, 6) // K(k + 1, k) = phai(k+1, k) * G
		uv2M(p->invNE->G, p->km.H, tao, i, j, 6) // tao = G H
		for (i = 0; i < 6; i++) for (j = 0; j < 6; j++) {
			if (i == j) tao.m[i][j] = 1 - tao.m[i][j];
			else tao.m[i][j] *= -1;
		} // tao = I - tao
		//tao.Show("tao");

		// P = phi * P * phi*
		/*
		if (p->parent != NULL) {
			MpM(p->km.phai, tao, phi, i, j, k, 6)
			M6T2M6(phi, phi_T, i, j)  // phi_T = phi*
			//phi_T.Show("phi_T");
			MpM(phi, p->invNE.P, m1, i, j, k, 6)
			MpM(m1, phi_T, p->parent->invNE.P, i, j, k, 6)
		}
		*/

		// Pp = tao * P
		MpM(tao, p->invNE->P, Pp, i, j, k, 6)

		// ignore the calculation of KAI(k + 1, k) = phai(k + 1, k) * tao(k)

		p = (CLUSTER*)(p->parent);
	}

	p = mm.tip[0]; //from tip
	VECTOR<6> z, zp;
	double eps = 0;
	for (i = 0; i < 6; i++) {v1.v[i] = 0; v2.v[i] = 0;} // v1 = 0; v2 = 0;
	while (p != NULL) {
		CHILD_CLUSTER(child, p, 0, i, j)
		if (child != NULL) {
			// v1 = pahi(k, k-1) * zp
			MpV6(child->km.phai, zp, v1)
		}
		// v2 = P(k) * a(k)
		MpV6(p->invNE->P, p->invNE->a, v2)
		for (i = 0; i < 6; i++) {
			z.v[i] = v1.v[i] + v2.v[i] + p->invNE->b.v[i] - p->dyn->fc.v[i];
			//z.v[i] = v1.v[i] + v2.v[i] + p->invNE->b.v[i] - p->dyn->fc.v[i];
		}

		scale_uv6(p->km.H, z, eps)
		eps = p->invNE->T - eps;
		p->invNE->gamma = eps / p->invNE->D;
		for (i = 0; i < 6; i++) zp.v[i] = z.v[i] + p->invNE->G.v[i] * eps;
		p = (CLUSTER*)(p->parent);
	}

	p = mm.base;
	VECTOR<6> alpha, alphap;
	eps = 0;
	while (p != NULL) {
		if (p->parent == NULL) { // this is base
			for (i = 0; i < 6; i++) alphap.v[i] = mm.alpha0.v[i];
			if (mm.free_base) {
				p->km.tpp = 0; // base has no parent hinge
				for (i = 0; i < 6; i++) p->bkm->alpha.v[i] = mm.alpha0.v[i] + p->invNE->a.v[i];

				CHILD_CLUSTER(p, p, 0, i, j) // linear polymer
				continue;
			}
		}
		else {
			// alphap(k) = phai_T(k+1, k) * alpha(k+1)
			MpV6(p->km.phai_T, p->bkm->alpha, alphap)
		}
		// eps = G * alphap
		scale_uv6(p->invNE->G, alphap, eps)

		p->km.tpp = p->invNE->gamma - eps;
		//p->km.tpp -=  ksi * p->km.tp; // NEIMO-Hoover dynamic, see J. Phys. Chem. 100, 10508, 1996
		for (i = 0; i < 6; i++) p->bkm->alpha.v[i] = alphap.v[i] + p->invNE->a.v[i] + p->km.H.v[i] * p->km.tpp;

		CHILD_CLUSTER(p, p, 0, i, j) // linear polymer
	}

	if (trace) {
		//calClusterAlpha(mm);
		p = mm.tip[0];
		while (p != NULL) {
			CHILD_CLUSTER(child, p, 0, i, j) // linear polymer

			if (child != NULL) MpV(child->km.phai, child->invNE->f, p->invNE->f, i, j, 6)
			MpV6(p->invNE->M, p->bkm->alpha, v1)
			for (i = 0; i < 6; i++) p->invNE->f.v[i] += v1.v[i] + p->invNE->b.v[i] + p->dyn->fc.v[i];
			scale_uv6(p->km.H, p->invNE->f, p->invNE->T)
			p = (CLUSTER*)(p->parent);
		}
	}
}
#endif

// Hoover motion will induce correction on the acceleration of child cluster
// to make sure the inertia moment of whole molecule conservation,
// the base cluster will have correction on its acceleration
/*
void NEIMO_CalFreeBaseHooverAccel(MMOLECULE &mm, float ksi, VECTOR<6> &alpha_base) {
	int i = 0, j = 0, k = 0;

	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *pc = NULL;
	CLUSTER *parent = NULL;

	VECTOR<6> vt, f0;
	VECTOR<3> dr;
	MATRIX<6> phai, phai_T;
	MATRIX<6> m_cluster, Imc, m1;

	ptl = mm.ctree.tl; // from base
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;
			parent = (CLUSTER*)(pc->parent);
			
			// dr is from mass center to cluster Op
			for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm.mDynKin.cm.v[i];
			PHAI(phai, dr, i, j)
			MpM(phai, pc->invNE->M, Imc, i, j, k, 6)  // Imc = phai(mc==>cluster) * I_cluster
			// dr is from base to cluster
			for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - CLUSTER_ORIGIN(mm.base, i);
			PHAI(phai, dr, i, j)
			M6T2M6(phai, phai_T, i, j)
			MpM(Imc, phai_T, m1, i, j, k, 6)  //m1 = phai(mc==>cluster) * I_cluster * phai_T(base==>cluster)
			MplusM(m_cluster, m1, m_cluster, i, j)  // m_cluster += m1

			if (parent == NULL) { // this is base
				Vzero(pc->invNE->dalpha_H, i)
				plist = plist->next; continue;
			}
			MpV6(pc->km.phai_T, parent->invNE->dalpha_H, vt)
			for (i = 0; i < 6; i++) pc->invNE->dalpha_H.v[i] = vt.v[i] - pc->km.H.v[i] * ksi * pc->km.tp;
			MpV6(Imc, pc->invNE->dalpha_H, vt) // vt = Imc * dalpha_H
			for (i = 0; i < 6; i++) f0.v[i] += vt.v[i]; // f0 += Imc * dalpha_H

			plist = plist->next;
		}
		ptl = ptl->next;
	}

	// alpha_base = m_cluster)^-1 * (-f0)
	for (i = 0; i < 6; i++) f0.v[i] = -f0.v[i];
	InverseMatrix(m_cluster, m1);
	MpV6(m1, f0, alpha_base)
}
*/
// Hoover motion will induce correction on the acceleration of child cluster
// to make sure the inertia moment of whole molecule conservation,
// the base cluster will have correction on its acceleration
void NEIMO_CalHooverForce_GlobalSpeed_Thread(MMOLECULE *mm, int nc1, int nc2, double& ksi, VECTOR<6> &force) {
	int i = 0, nc = 0;
	VECTOR<6> force_mc;

	CLUSTER *pc = NULL;
	CLUSTER *parent = NULL;

	VECTOR<6> vt, f0;

	double dksi = 0;
	double Ekmax = ::Ek0Internal * 4;

	bool reset = false;;
	double ksi_t, ksi_r;
	VECTOR3 dr;

	bool bForceMassCenter = (::MM_SpeedVerlet_method == 0 ? true : false);

	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		parent = (CLUSTER*)(pc->parent);
		V6zero(pc->dyn->f_Hoover)
		dksi = ksi;
		//if (ksi < 0 && pc->Ek > 1) continue; // cluster is too fast already
		/*
		if (::bCoarseGrain_MM) {
			if (pc->Ek > 0.8 && ksi < 0.2) {reset = true; ksi_t = 0.2; ksi_r = 0.2;}
			else if (pc->Ek < 0.1) {reset = true; ksi_r = -0.2; ksi_t = dksi;}
			else reset = false;
		}

		if (reset) {
			f0.v[0] = -ksi_t * pc->bkm->V0.v[3] * pc->M;
			f0.v[1] = -ksi_t * pc->bkm->V0.v[4] * pc->M;
			f0.v[2] = -ksi_t * pc->bkm->V0.v[5] * pc->M;
			VECT3((*(pc->Op)), pc->cm, dr)
			V3PV3(dr, f0, vt)

			V6zero(pc->invNE->dalpha_H)
			for (i = 0; i < 3; i++) {
				//pc->invNE->dalpha_H.v[i] = -dksi * pc->bkm->V0.v[i];
				pc->invNE->dalpha_H.v[i] = -ksi_r * pc->km.H.v[i] * pc->km.tp;
			}
			MpV6(pc->invNE->M, pc->invNE->dalpha_H, pc->dyn->f_Hoover) // f_Hoover = M(cluster) * dalpha_H
			for (i = 0; i  < 3; i++) {pc->dyn->f_Hoover.v[i] += vt.v[i]; pc->dyn->f_Hoover.v[i + 3] += f0.v[i];}
		}
		else {
		*/
			for (i = 0; i < 6; i++) {
				pc->invNE->dalpha_H.v[i] = -dksi * pc->bkm->V0.v[i];
				//pc->invNE->dalpha_H.v[i] = -dksi * pc->km.H.v[i] * pc->km.tp;
			}
			MpV6(pc->invNE->M, pc->invNE->dalpha_H, pc->dyn->f_Hoover) // f_Hoover = M(cluster) * dalpha_H
		//}

		if (bForceMassCenter) {
			MpV6(pc->invNE->phai_cm2me, pc->dyn->f_Hoover, f0) // f0 = phai(cm=>pc) * vt
			V6plusV6(force_mc, f0, force_mc)
		}
	}
	if (bForceMassCenter) {V62V6(force_mc, force)}
	//return force_mc;
}

void NEIMO_CalHooverForce_GlobalSpeed(MMOLECULE &mm, float ksi, VECTOR<6> &force_mc) {
	double dksi = (double)ksi;
	V6zero(force_mc)
	if (mm.nCluster > NCLUSTERS_PARALLEL) {
		MacroMoleculeOperate3<MMOLECULE, double, VECTOR<6> >((void*)(&NEIMO_CalHooverForce_GlobalSpeed_Thread), mm, 0, mm.nCluster, dksi, MAX_THREADS, force_mc);
		return;
	}
	else NEIMO_CalHooverForce_GlobalSpeed_Thread(&mm, 0, mm.nCluster - 1, dksi, force_mc);
}


void NEIMO_CalHooverForce(MMOLECULE &mm, float ksi, VECTOR<6> &force_mc) {
	NEIMO_CalHooverForce_GlobalSpeed(mm, ksi, force_mc);
	/*
	if (!mm.free_base) return;
	return;

	int i, j;
	VECTOR<6> f;
	VECTOR3 dr = mm.r0; // cm ==> base
	INV_V3(dr) // dr : base ==> cm
	MATRIX<6> phai;
	PHAI(phai, dr, i, j)
	//M6T2M6(phai, phai_T, i, j)

	MpV6(phai, force_mc, f) //  f_base = phai (base ==> cm) * f_cm
	INV_V6(f)
	
	mm.base->dyn->f_Hoover = f;
	//V6plusV6(mm.base->dyn->f_Hoover, f, mm.base->dyn->f_Hoover)
	V6zero(force_mc);
	*/
}

#ifndef _DISABLE_
// Hoover motion will induce correction on the acceleration of child cluster
// to make sure the inertia moment of whole molecule conservation,
// the base cluster will have correction on its acceleration
// same as NEIMO_CalInternalHooverForce
void NEIMO_CalInternalHooverForce_1(MMOLECULE &mm, VECTOR<6> &force_mc) { // ksi = 1
	int i = 0, j = 0, k = 0;
	V6zero(force_mc)

	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *pc = NULL;
	CLUSTER *parent = NULL;

	VECTOR<6> vt, f0;
	//VECTOR<3> dr;
	//MATRIX<6> phai;

	ptl = mm.ctree.tl; // from base
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;
			parent = (CLUSTER*)(pc->parent);
			
			// dr is from mass center to cluster Op
			//for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm.mDynKin.cm.v[i];
			//PHAI(phai, dr, i, j)

			if (parent == NULL) { // this is base
				V6zero(pc->invNE->dalpha_H)
				plist = plist->next; continue;
			}
			MpV6(pc->km.phai_T, parent->invNE->dalpha_H, vt)
			for (i = 0; i < 6; i++) pc->invNE->dalpha_H.v[i] = vt.v[i] - pc->km.H.v[i] * pc->km.tp; // ksi = 1
			MpV6(pc->invNE->M, pc->invNE->dalpha_H, vt) // vt = M(cluster) * dalpha_H

			V62V6(vt, pc->dyn->f0_Hoover_InternalKinEnergy) // put Hoover force on the cluster

			MpV6(pc->invNE->phai_cm2me, vt, f0) // f0 = phai(cm=>pc) * vt
			//for (i = 0; i < 6; i++) force_mc.v[i] += f0.v[i];
			V6plusV6(force_mc, f0, force_mc)

			plist = plist->next;
		}
		ptl = ptl->next;
	}
}

void distribute_MassCenter_Hoover_force_Thread(MMOLECULE *mm, int nc1, int nc2, VECTOR<6> &alpha0) {
	VECTOR<6> alpha, f;

	CLUSTER *pc = NULL;
	int nc = 0;

	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		//for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm.mDynKin.cm.v[i]; // from mass center to cluster Op
		//PHAI(phai, dr, i, j)
		//M6T2M6(phai, phai_T, i, j)
		//MpV(phai_T, alpha0, alpha, i, j, 6)
		MpV6(pc->invNE->phai_T_cm2me, alpha0, alpha)
		MpV6(pc->invNE->M, alpha, f)
		//for (i = 0; i < 6; i++) pc->dyn->f_Hoover.v[i] += f.v[i];
		V6plusV6(pc->dyn->f_Hoover, f, pc->dyn->f_Hoover)
	}
}

void distribute_MassCenter_Hoover_force(MMOLECULE &mm, int nThreads) {
	VECTOR<6> alpha0, alpha, f;
	//VECTOR3 dr;
	//MATRIX<6> phai, phai_T;

	CLUSTER *pc = NULL;
	int nc = 0;

	VECTOR<6> fc;
	V6minusV6(mm.mDynKin.f_FrameHoover, mm.mDynKin.f_Hoover, fc)
	MpV6(mm.mDynKin.invI, fc, alpha0) // mass center acceleration induced from f_Hoover
	if (nThreads > 2 && mm.nCluster > NCLUSTERS_PARALLEL) MacroMoleculeOperate1<MMOLECULE, VECTOR<6> >((void*)(&distribute_MassCenter_Hoover_force_Thread), mm, 0, mm.nCluster, alpha0, MAX_THREADS);
	else distribute_MassCenter_Hoover_force_Thread(&mm, 0, mm.nCluster, alpha0);
	V62V6(fc, mm.mDynKin.f_Hoover)
}
#endif // _DISABLE_

/*
void Hoover_Dyn(MMOLECULE &mm, float ksi) {
	CLUSTER *pc = NULL;
	int i = 0;
	double dtpp = 0;

	for (i = 0; i < mm.nCluster; i++) {
		pc = mm.cluster + i;
		//if (FABS(p->km.tpp) < 2e-6) p->km.tpp = 0;
		dtpp = ksi * pc->km.tp;
		pc->km.tpp -=  dtpp; // NEIMO-Hoover dynamic, see J. Phys. Chem. 100, 10508, 1996
	}
	if (mm.free_base) mm.base->km.tpp = 0;

	//VECTOR<6> alpha_base_Hoover;
	//NEIMO_CalFreeBaseHooverAccel(mm, ksi, alpha_base_Hoover);
	//for (i = 0; i < 6; i++) mm.alpha0.v[i] += alpha_base_Hoover.v[i];
}
*/


// modification in the midnight of 10/30/2009: add a matrix variable Pp in INVERSE_NE, 
// and the related code on calculation of NEIMO MATRIX relating to the algorithm defined in Page 264 of paper
// "A Fast Recursive Althorithm for Molecular Dynamics Simulation", by A. Jain et al,
// Now, we follow exactly the althorithm defined in the paper!!!
// Using the code above, it is found that the torque between the cluster is non-zero along the hinge when all the forces
// is put in fc. This is confliction with Eq. 4.7 in the paper, without explicit torque on the hinge.
void NEIMO_MATRIX(MMOLECULE &mm, bool new_config, bool bAcceleration) {
	int i = 0, j = 0, k = 0, nc = 0;

	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *pc = NULL;
	CLUSTER *child = NULL;

	MATRIX<6> Pp, tao, m1, m2;//, phi, phi_T;
	VECTOR<6> v1, v2;

	double d;

	VECTOR<6> *pv1 = NULL, *pv2 = NULL;
	MATRIX<6> *pm1 = NULL, *pm2 = NULL, *pm3 = NULL;

	if (!new_config) goto _force_related;

	for (nc = 0; nc < mm.nCluster; nc++) {
		pm1 = &(mm.cluster[nc].invNE->P); M6zero((*pm1))
		pm1 = &(mm.cluster[nc].invNE->Pp); M6zero((*pm1))
	}

	ptl = mm.ctree.get_tree_level(-1); // from tip to base
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;
			M6zero(Pp) 
			if (pc == mm.base) {
				//plist = plist->next;
				if (mm.free_base) {
					memset(pc->invNE->G.v, 0, SIZE_V6);
					memset(pc->invNE->beta.v, 0, SIZE_V6);
					pc->invNE->D = 0;
					pc->invNE->gamma = 0;
				}
				else {
					/*
					for (i = 0; i < 6; i++) {
						pc->invNE->G.v[i] = 0; 
						pc->invNE->beta.v[i] = 0;
					}
					*/
					memset(pc->invNE->G.v, 0, SIZE_V6);
					memset(pc->invNE->beta.v, 0, SIZE_V6);

					pc->invNE->D = 1e20;
					pc->invNE->gamma = 0;
				}
				//continue;
			}

			for (nc = 0; nc < pc->nHinges; nc++) { // run over all the children
				if (nc == pc->n2Parent) continue; // this hinge connects to parent
				// here we assumed that all the linked clusters (except the parent) are the children of this cluster
				// This is NOT true for the case with round-ring linked clusters inside the macromolecule
				child = (CLUSTER*)(pc->hinge[nc].plink);

				pm1 = &(child->invNE->Pp); pm2 = &(child->km.phai);
				MpM((*pm2), (*pm1), m1, i, j, k, 6)  // m1 = phai * Pp -- child
				pm1 = &(child->km.phai_T);
				MpM(m1, (*pm1), m2, i, j, k, 6) // m2 = m1 * phai_T = phai * Pp * phai_T -- child
				M6plusM6(m2, Pp, Pp)  // Pp += phai * Pp * phai_T -- child
			}
			pm1 = &(pc->invNE->M); pm2 = &(pc->invNE->P);
			M6plusM6(Pp, (*pm1), (*pm2))  // pc->invNE->P = SUM[child->km.phai * child->invNE->Pp * child->km.phai_T] + pc->invNE->M

			if (pc->parent == NULL) { // base of molecule
				// INVERSE_NE for base cluster
				// because base cluster is free, H  = I (unit vector to any direction) while not Hinge !
				// So, D = P; G = P * H * D^-1 = I; tao = I - G * H = 0; Pp = 0;
				// gamma is Matrix, gamma = P^-1 * (P * a + b +/- fc)
				InverseMatrix<6>(pc->invNE->P, mm.binvNE.invD);
			}
			else {
				pm1 = &(pc->invNE->P); pv1 = &(pc->km.H);
				MpV6((*pm1), (*pv1), v1) // v1 = pc->invNE->P * H*
				scale_uv6((*pv1), v1, pc->invNE->D)  // pc->invNE->D = pc->km.H * pc->invNE->P * pc->km.H*
				pv1 = &(pc->invNE->G); d = 1.0 / pc->invNE->D;
				TIME6(v1, d, (*pv1)) // pc->invNE->G = pc->invNE->P * pc->km.H / pc->invNE->D

				pv1 = &(pc->invNE->G); pv2 = &(pc->km.H); pm1 = &(pc->invNE->tao);
				uv2TAO6((*pv1), (*pv2), (*pm1)) // tao = I - tao

				// Pp = tao * P
				pm1 = &(pc->invNE->tao); pm2 = &(pc->invNE->P); pm3 = &(pc->invNE->Pp);
				MpM((*pm1), (*pm2), (*pm3), i, j, k, 6)  // pc->invNE->Pp = pc->invNE->tao * pc->invNE->P

				//MpV(pc->km.phai, pc->invNE->G, child->invNE->K, i, j, 6) // K(k + 1, k) = phai(k+1, k) * G

				// ignore the calculation of KAI(k + 1, k) = phai(k + 1, k) * tao(k)
			}

			plist = plist->next;
		}
		ptl = ptl->pre;
	}

_force_related:
	if (!bAcceleration) return; // we do not calculate the acceleration

	VECTOR<6> z, zp;
	double eps = 0;
	ptl = mm.ctree.get_tree_level(-1); // from tip to base
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;
			V6zero(z)
			for (nc = 0; nc < pc->nHinges; nc++) { // run over child
				if (nc == pc->n2Parent) continue; // this hinge connects to parent
				// here we assume all the clusters (except the parent) connecting to this cluster are the children of this cluster
				// This assumption is NOT true when the macromolecule has a ring-like structure inside
				child = (CLUSTER*)(pc->hinge[nc].plink);

				if (child == NULL) continue;
				pm1 = &(child->km.phai); pv1 = &(child->invNE->zp);
				MpV6((*pm1), (*pv1), v1)  // v1 = phai(k, k-1) * zp
				//for (i = 0; i < z.n; i++) z.v[i] += v1.v[i]; // z += v1 (for all children)
				V6plusV6(z, v1, z)
			}
			pm1 = &(pc->invNE->P); pv1 = &(pc->invNE->a);
			MpV6((*pm1), (*pv1), v2)  // v2 = P(k) * a(k)
			pv1 = &(pc->invNE->b); pv2 = &(pc->dyn->f);
			for (i = 0; i < 6; i++) {
				z.v[i] += v2.v[i] + pv1->v[i] - pv2->v[i]; 
				//z.v[i] += v2.v[i] + pc->invNE->b.v[i] - pc->dyn->fc.v[i];
			}

			if (pc == mm.base) {
				if (mm.free_base) {
					// INVERSE_NE for base cluster
					// because base cluster is free, H  = I (unit vector to any direction) while not Hinge !
					// So, D = P; G = P * H * D^-1 = I; tao = I - G * H = 0; Pp = 0;
					// gamma is Matrix, gamma = -P^-1 * (P * a + b +/- fc)
					pc->invNE->gamma = 0;
					pm1 = &(mm.binvNE.invD); 
					MpV6((*pm1), z, mm.alpha0)
					INV_V6(mm.alpha0) // mm.alpha0 = -P^-1 * (P * a + b +/- fc)
				}
				else pc->invNE->gamma = 0; //?
			}
			else {
				scale_uv6(pc->km.H, z, eps)
				eps = pc->invNE->T - eps;
				pc->invNE->gamma = eps / pc->invNE->D;
				//pc->invNE->gamma = eps / pc->invNE->D - ksi * pc->km.tpp;
			}

			for (i = 0; i < 6; i++) pc->invNE->zp.v[i] = z.v[i] + pc->invNE->G.v[i] * eps;

			plist = plist->next;
		}
		ptl = ptl->pre;
	}
}

// calculate the last step in NEIMO procedure -- cluster acceleration
void NEIMO_CalClusterAcceleration(MMOLECULE &mm) {
	int i = 0, j = 0, k = 0;

	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *pc = NULL;
	CLUSTER *parent = NULL;

	VECTOR<6> alpha, alphap;
	double eps = 0;

	MATRIX<6> *pm = NULL;
	VECTOR<6> *pv = NULL;

	ptl = mm.ctree.tl; // from base
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;
			V6zero(alphap)
			parent = (CLUSTER*)(pc->parent);
			eps = 0;
			if (pc->parent == NULL) {
				if (mm.free_base) {
					// INVERSE_NE for base cluster
					// because base cluster is free, H  = I (unit vector to any direction) while not Hinge !
					// So, D = P; G = P * H * D^-1 = I; tao = I - G * H = 0; Pp = 0;
					// gamma is Matrix, gamma = -P^-1 * (P * a + b +/- fc)
					// alpha = gamma + a

					V6plusV6(mm.alpha0, pc->invNE->a, mm.alpha0)
					V62V6(mm.alpha0, pc->bkm->alpha)
					pc->km.tpp = 0; // base has no parent hinge
				}
				else {
					pc->km.tpp = 0; // base has no parent hinge
					V6zero(pc->bkm->alpha);
				}
			}
			else {
				pm = &(pc->km.phai_T); pv = &(parent->bkm->alpha);
				MpV6((*pm), (*pv), alphap)  // alphap(k) = phai_T(k+1, k) * alpha(k+1)
				// eps = G * alphap
				scale_uv6(pc->invNE->G, alphap, eps)

				pc->km.tpp = pc->invNE->gamma - eps;
				//pc->km.tpp -=  ksi * pc->km.tp; // NEIMO-Hoover dynamic, see J. Phys. Chem. 100, 10508, 1996
				for (i = 0; i < 6; i++) pc->bkm->alpha.v[i] = alphap.v[i] + pc->invNE->a.v[i] + pc->km.H.v[i] * pc->km.tpp;
			}
			plist = plist->next;
			continue;
		}
		ptl = ptl->next;
	}
}

#ifndef _DISABLE_
void NEIMO_CalFreeBaseAccel(MMOLECULE &mm, bool bCalMatrix) {
	int i = 0, j = 0, k = 0;

	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *pc = NULL;
	CLUSTER *parent = NULL;

	VECTOR<6> vt, alpha, alpha_total;
	//VECTOR<3> dr;
	MATRIX<6> m1, m2, mA; //phai, phai_T;

	MATRIX<6> *pm1 = NULL, *pm2 = NULL, *pm3 = NULL;
	VECTOR<6> *pv1 = NULL, *pv2 = NULL;

	if (bCalMatrix) {M6zero(mm.mDynKin.M0)}
	V6zero(mm.b) V6zero(mm.f)

	ptl = mm.ctree.tl; // from base
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;
			parent = (CLUSTER*)(pc->parent);

			// without Hoover correction
			for (i = 0; i < 6; i++) pc->invNE->beta.v[i] = pc->km.H.v[i] * pc->invNE->gamma + pc->invNE->a.v[i];
			// with Hoover correction
			//for (i = 0; i < 6; i++) pc->invNE->beta.v[i] = pc->km.H.v[i] * (pc->invNE->gamma - ksi * pc->km.tp) + pc->invNE->a.v[i];

			if (bCalMatrix) {
				if (parent != NULL) {
					// mA = I - H* x G*
					/*
					for (i = 0; i < 6; i++) for (j = 0; j < 6; j++) {
						mA.m[i][j] = -(pc->km.H.v[i] * pc->invNE->G.v[j]);
						if (i == j) mA.m[i][j] += 1;
					}
					*/
					pv1 = &(pc->km.H); pv2 = &(pc->invNE->G);
					uv2TAO6((*pv1), (*pv2), mA)

					pm1 = &(pc->km.phai_T); pm2 = &(pc->invNE->tao);
					MpM(mA, (*pm1), (*pm2), i, j, k, 6)  // tao = (I - H* G*) * phai_T
				}
			}

			if (parent != NULL) {
				// beta(k) = beta(k) + tao(k) * beta(k+1)
				//         = beta(k) + tao(k) * beta(k+1) + tao(k) * tao(k+1) * beta(k+2) + ... to base
				pm1 = &(pc->invNE->tao); pv1 = &(parent->invNE->beta);
				MpV6((*pm1), (*pv1), vt) // vt = tao(k) * parent->beta [or, beta(k+1)]

				pv1 = &(pc->invNE->beta);
				V6plusV6((*pv1), vt, (*pv1)) // pc->invNE->beta = beta(k) + tao(k) * beta(k+1)
			}

			if (bCalMatrix) {
				if (parent == NULL) {
					pm1 = &(pc->invNE->M); pm2 = &(pc->invNE->m0);
					M2M((*pm1), (*pm2), i, j, 6) // m0 = M, this is base
					pm2 = &(mm.mDynKin.M0); // total inertia moment relative to the base
					M6plusM6((*pm1), (*pm2), (*pm2)) // mm.M0 += pc->invNE->M
				}
				else {
					// dr is from base to cluster's cordinate
					//for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - CLUSTER_ORIGIN(mm.base, i);
					//PHAI(phai, dr, i, j)
					pm1 = &(pc->invNE->phai_base2me); pm2 = &(pc->invNE->M); pm3 = &(pc->invNE->m0);
					MpM((*pm1), (*pm2), (*pm3), i, j, k, 6) // pc->invNE->m0 = phai(base=>cluster) * M

					//M6T2M6(phai, phai_T, i, j)
					pm1 = &(pc->invNE->phai_T_base2me);
					MpM(mA, (*pm1), m1, i, j, k, 6)  // m1 = mA * phai_T(base=>cluster)
					pm1 = &(pc->invNE->m0);
					MpM((*pm1), m1, m2, i, j, k, 6) //m2 = phai(base=>cluster) * M * Ak * phai_T(base=>cluster)
					pm2 = &(mm.mDynKin.M0);
					M6plusM6((*pm2), m2, (*pm2)) // mm.mDynKin.M0 += m2, total inertia moment relative to the base
				}
			}

			pm1 = &(pc->invNE->m0); pv1 = &(pc->invNE->beta);
			MpV6((*pm1), (*pv1), alpha)  // alpha = phai(base==>cluster) * M(k) * beta(k) -- beta(k) is accumuated from base
			V6plusV6(alpha_total, alpha, alpha_total)

			//total b & f relative to base
			if (parent == NULL) { // base, we do not need to do transform
				pv1 = &(pc->invNE->b);
				V6plusV6(mm.b, (*pv1), mm.b)
				pv1 = &(pc->dyn->f);
				V6plusV6(mm.f, (*pv1), mm.f)
			}
			else {
				pm1 = &(pc->invNE->phai_base2me); pv1 = &(pc->invNE->b);
				MpV6((*pm1), (*pv1), vt) // vt = phai(base==>cluster) * b(k)
				V6plusV6(mm.b, vt, mm.b) // accumulate to mm.b, total b relative to the base
				pv1 = &(pc->dyn->f);
				MpV6((*pm1), (*pv1), vt) // vt = phai(base==>cluster) * f(k)
				V6plusV6(mm.f, vt, mm.f) // accumulate to mm.f, total f relative to the base
			}

			plist = plist->next;
		}
		ptl = ptl->next;
	}

	if (bCalMatrix) InverseMatrix<6>(mm.mDynKin.M0, mm.mDynKin.invM0);
	for (i = 0; i < 6; i++) vt.v[i] = mm.f.v[i] - mm.b.v[i] - alpha_total.v[i];
	MpV6(mm.mDynKin.invM0, vt, mm.alpha0)
}
#endif // _DISABLE_

void check_force_consistent_mass_center(MMOLECULE &mm) {
	VECTOR<6> alpha, f, ftotal;
	//VECTOR3 dr;
	//MATRIX<6> phai, phai_T;

	CLUSTER *pc = NULL;
	int nc = 0, i;
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		//for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - mm.mDynKin.cm.v[i]; // from mass center to cluster Op
		//PHAI(phai, dr, i, j)
		MpV6(pc->invNE->phai_cm2me, pc->dyn->fc, f)
		//for (i = 0; i < 6; i++) ftotal.v[i] += f.v[i];
		V6plusV6(ftotal, f, ftotal)
	}
	//for (i = 0; i < 6; i++) f.v[i] = ftotal.v[i] - mm.mDynKin.fc.v[i];
	V6minusV6(ftotal, mm.mDynKin.fc, f)
	// f is supposed to be 0

	V6zero(ftotal)
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		//for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - CLUSTER_ORIGIN(mm.base, i); // from base cluster to cluster Op
		//PHAI(phai, dr, i, j)
		MpV6(pc->invNE->phai_base2me, pc->dyn->fc, f)
		//for (i = 0; i < 6; i++) ftotal.v[i] += f.v[i];
		V6plusV6(ftotal, f, ftotal)
	}
	for (i = 0; i < 6; i++) f.v[i] = ftotal.v[i] - mm.mDynKin.f.v[i];
	// f is supposed to be 0
	i = i;
}

void NEIMO(MMOLECULE &mm, bool bCalMatrix) {
	// do this out of this function
	mm.calTotalForce();

	//mm.CalMolDyn(); // calculate whole molecule's dynamic, acceleration
	// set the molecular spatial moment and acceleration moment at base cluster
	//mm.SetBaseForce();
	//check_force_consistent_mass_center(mm);

	NEIMO_MATRIX(mm, bCalMatrix, true);
	//if (mm.free_base) NEIMO_CalFreeBaseAccel(mm, bCalMatrix);
	NEIMO_CalClusterAcceleration(mm);
	//Hoover_Dyn(mm, ksi);

	// for debug 
	/*
	VECTOR<6> Iac;
	calClusterAlpha(mm);
	calSpatialAccel0(mm, Iac);
	ksi = ksi;
	*/
}

void NEIMO_GivenHoover(MMOLECULE &mm, float ksi, bool new_config) {
	int i, n;
	mm.calTotalForce();
	//mm.CalMolDyn(); // calculate whole molecule's dynamic, acceleration
	// set the molecular spatial moment and acceleration moment at base cluster
	//mm.SetBaseForce();
	//check_force_consistent_mass_center(mm);

	NEIMO_MATRIX(mm, new_config, true);
	//if (mm.free_base) NEIMO_CalFreeBaseAccel(mm, new_config);
	NEIMO_CalClusterAcceleration(mm);

	// for debug 
	/*
	VECTOR<6> Iac;
	calClusterAlpha(mm);
	calSpatialAccel0(mm, Iac);
	ksi = ksi;
	*/
}

void NEIMO_FixedBase(MMOLECULE &mm, float ksi, bool new_config) {
	mm.calTotalForce();
	NEIMO_MATRIX(mm, new_config, true);
	NEIMO_CalClusterAcceleration(mm);
}

#ifndef _DISABLE_
// this is a thread to calculate the matrix in NEIMO_CalFreeBaseAccel, NON-related to the force and/or acceleration
void NEIMO_MOMENT_MATRIX1_THREAD(MMOLECULE *mm, int nc1, int nc2, MATRIX<6> &M0) {
	int i = 0, j = 0, k = 0, nc = 0;

	CLUSTER *pc = NULL, *parent = NULL;

	MATRIX<6> m1, m2, mA;

	MATRIX<6> *pm1 = NULL, *pm2 = NULL, *pm3 = NULL;
	VECTOR<6> *pv1 = NULL, *pv2 = NULL;

	M6zero(M0)

	for (nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		parent = (CLUSTER*)(pc->parent);
		if (parent != NULL) {
			// mA = I - H* x G*
			/*
			for (i = 0; i < 6; i++) for (j = 0; j < 6; j++) {
				mA.m[i][j] = -(pc->km.H.v[i] * pc->invNE->G.v[j]);
				if (i == j) mA.m[i][j] += 1;
			}
			*/
			pv1 = &(pc->km.H); pv2 = &(pc->invNE->G);
			uv2TAO6((*pv1), (*pv2), mA)
			pm1 = &(pc->km.phai_T); pm2 = &(pc->invNE->tao);
			MpM(mA, (*pm1), (*pm2), i, j, k, 6)  // tao = (I - H* G*) * phai_T
		}

		if (parent == NULL) {
			pm1 = &(pc->invNE->M); pm2 = &(pc->invNE->m0);
			M2M((*pm1), (*pm2), i, j, 6) // m0 = M, this is base
			pm2 = &(M0); // total inertia moment relative to the base
			M6plusM6((*pm1), (*pm2), (*pm2)) // mm.mDynKin.M0 += pc->invNE->M
		}
		else {
			// dr is from base to cluster's cordinate
			//for (i = 0; i < 3; i++) dr.v[i] = CLUSTER_ORIGIN(pc, i) - CLUSTER_ORIGIN(mm.base, i);
			//PHAI(phai, dr, i, j)
			pm1 = &(pc->invNE->phai_base2me); pm2 = &(pc->invNE->M); pm3 = &(pc->invNE->m0);
			MpM((*pm1), (*pm2), (*pm3), i, j, k, 6) // pc->invNE->m0 = phai(base=>cluster) * M

			//M6T2M6(phai, phai_T, i, j)
			pm1 = &(pc->invNE->phai_T_base2me);
			MpM(mA, (*pm1), m1, i, j, k, 6)  // m1 = mA * phai_T(base=>cluster)
			pm1 = &(pc->invNE->m0);
			MpM((*pm1), m1, m2, i, j, k, 6) //m2 = phai(base=>cluster) * M * Ak * phai_T(base=>cluster)
			M6plusM6(M0, m2, M0) // mm.mDynKin.M0 += m2, total inertia moment relative to the base
		}
	}
	//return M0;
}
#endif // _DISABLE_

// this is a combination of functions NEIMO_MATRIX, and NEIMO_CalFreeBaseAccel
// about the part for calculation NON-related to the force and/or acceleration
void NEIMO_MOMENT_MATRIX(MMOLECULE &mm) {
	// comparing to the same function above, the explicit codes were replaced by calling NEIMO_MATRIX function directly
	NEIMO_MATRIX(mm, true, false); // calculating the matrix, but not the acceleration
/*
	// also calculating the matrix required for the acceleration speed of base cluster 
	M6zero(mm.mDynKin.M0)
	if (mm.nCluster > NCLUSTERS_BIG_MM) {
		MacroMoleculeOperate2<MMOLECULE, MATRIX<6> >((void*)(&NEIMO_MOMENT_MATRIX1_THREAD), mm, 0, mm.nCluster, MAX_THREADS, mm.mDynKin.M0);
	}
	else {
		NEIMO_MOMENT_MATRIX1_THREAD(&mm, 0, mm.nCluster - 1, mm.mDynKin.M0);
	}

	InverseMatrix<6>(mm.mDynKin.M0, mm.mDynKin.invM0);
	*/
}

#if _SYS_ == _WINDOWS_SYS_
int NEIMO_MATRIX_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
void* NEIMO_MATRIX_Thread(void *void_par) { // thread-function
#endif
	NEIMO_THREAD_VARS *par = (NEIMO_THREAD_VARS*)void_par;
	NEIMO_MOMENT_MATRIX(*(par->mm)); // calculate the matrix only, disable acceleration calculation
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

void Start_NEIMO_MATRIX_CALC(MMOLECULE &mm, NEIMO_THREAD_VARS &var) {
	var.set_pars(&mm, 0);
#if _SYS_ == _WINDOWS_SYS_ 
	ResetEvent(var.hMMThread);
	var.thread = AfxBeginThread((AFX_THREADPROC)NEIMO_MATRIX_Thread, (LPVOID)(&var), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
#elif _SYS_ == _LINUX_SYS_
	pthread_create(&(var.thread), NULL, &NEIMO_MATRIX_Thread, (void *)(&var));
#endif
}

void NEIMO_HingeForce(MMOLECULE &mm) {
	int i = 0;

	LIST< LIST<CLUSTER> > *ptl = NULL;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER *pc = NULL;
	CLUSTER *child = NULL;

	VECTOR<6> *alpha;
	MATRIX<6> *m, *phai;//, *phai_T;

	VECTOR<6> *pv = NULL;

	VECTOR<6> f_child, f;

	//int l = 0;
	ptl = mm.ctree.get_tree_level(-1); // from tip
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;
			if (pc == mm.base) {V6zero((*(pc->fphinge))); plist = plist->next; continue;}

			V6zero(f_child)
			for (i = 0; i < pc->nHinges; i++) {
				if (i == pc->n2Parent) continue;
				child = (CLUSTER*)(pc->hinge[i].plink);
				phai = &(child->km.phai); pv = child->fphinge;
				MpV6((*phai), (*pv), f)
				V6plusV6(f_child, f, f_child)
			}
			m = &(pc->invNE->M); alpha = &(pc->bkm->alpha);
			MpV6((*m), (*alpha), f);
			for (i = 0; i < 6; i++) pc->fphinge->v[i] = f_child.v[i] + f.v[i] - pc->dyn->f.v[i] + pc->invNE->b.v[i];

			// check the torque and force between the clusters on hinge
			// The torque along the hinge is supposed to be same as that in pc->invNE->T
			// in this program, dihedral torque is calculated as the force on each atom explicitly, so pc->invNE->T = 0
			/*
			{
				char msg[256] = "\0";
				SVECTOR<double, 3> f, torq;
				SVECTOR<double, 3> vt;
				torq.v[0] = pc->fphinge->v[0]; torq.v[1] = pc->fphinge->v[1]; torq.v[2] = pc->fphinge->v[2];
				f.v[0] = pc->fphinge->v[3]; f.v[1] = pc->fphinge->v[4]; f.v[2] = pc->fphinge->v[5];
				double d1, d2;
				scale_uv3(torq, pc->up, d1) 
				V3PV3(torq, pc->up, vt)
				d2 = sqrt(vt.v[0] * vt.v[0] + vt.v[1] * vt.v[1] + vt.v[2] * vt.v[2]);
				sprintf(msg, "[P = %f, N = %f], ", d1, d2); show_log(msg, false);
				scale_uv3(f, pc->up, d1)
				V3PV3(f, pc->up, vt)
				d2 = sqrt(vt.v[0] * vt.v[0] + vt.v[1] * vt.v[1] + vt.v[2] * vt.v[2]);
				sprintf(msg, "[P = %f, N = %f]", d1, d2); show_log(msg, true);
			}
			*/

			plist = plist->next;
		}
		/*
		{
			char msg[256] = "\0";
			sprintf(msg, "level %d", l); show_log(msg, true);
			l++;
		}
		*/
		ptl = ptl->pre;
	}
	//show_log("mm over\n", true);
}


// ***********************************************************************************************
//* momentum of macromolecule is the sum of the momentum of each cluster. In the frame of internal coordination,
//* the speed of each cluster is relative to the Op of its hinge, V_i. Total momentum of macromolecule = 
//*  Sum_i {Phai[i==>base] * Pp_i * V_i}, where Pp_i = (I - G * H) * P_i, P_i = M_i + Phai_i * P_(i-1) * Phai_T_i