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
#include "ZMatrix.h"
#include "MM.h"
#include "var.h"
#include "read.h"
#include "Matrix.h"

#include "cg-mm.h"

extern bool search(ifstream &in, char *title, char *buffer);
extern bool is_bound(char *buffer, char *atom1, char *atom2, char **pars);

namespace _coarse_grain_ {
extern CGBOUND_DIPOLE_DISTRBT_DB cg_dipole_distrbt_db;

void CG_CLUSTER::setup_neighbors(int n_neighbors) {
	neighbor.release(); cgb.release();
	if (n_neighbors <= 0) return;
	int i, j;
	neighbor.SetArray(n_neighbors);
	cgb.SetArray(n_neighbors);

	CHAIN<CG_CLUSTER> *ch_neighbor = NULL;

	for (i = 0; i < n_neighbors; i++) {
		grandchildren_chain<CG_CLUSTER>(this, i + 1, &ch_neighbor);
		chain_array<CG_CLUSTER>(&ch_neighbor, neighbor.m[i]);
		release_chain<CG_CLUSTER>(&ch_neighbor, false);
		if (neighbor.m[i].n == 0) return;
		// the bound profiles
		cgb.m[i].SetArray(neighbor.m[i].n); 
		for (j = 0; j < cgb.m[i].n; j++) cgb.m[i].m[j] = NULL;
	}
}

void CG_CLUSTER::calBrownianPars()  {
	float d2 = this->d * this->d;
	this->ita_Brownian = ::ita_Brownian * d2;
	this->cd_trans = ::cd_drag_trans * d2;
	return;
}

bool CG_MMOLECULE::xyz_save(ofstream &out) {
	int nc;
	VECTOR3 *r = NULL;
	for (nc = 0; nc < nCluster; nc++) {
		r = cluster[nc].r;
		out<<cluster[nc].atom[0].aindx<<" "<<r->v[0]<<" "<<r->v[1]<<" "<<r->v[2]<<endl;
	}
	return true;
}

void CG_MMOLECULE::shiftMol(VECTOR3 &dr) {
	int i = 0, nc = 0, na = 0;
	CG_CLUSTER *pc = NULL;
	CGATOM *patom = NULL;
	for (nc = 0; nc < this->nCluster; nc++) {
		pc = cluster + nc;
		for (na = 0; na < pc->nAtoms; na++) {
			patom = pc->atom + na;
			V3plusV3(patom->r, dr, patom->r)
		}
		if (bImplicitSolvent) {
			V3plusV3(pc->rc_solv, dr, pc->rc_solv); // geometrical center of cluster, in real cordinate
		}
		pc->update_vp();
	}
	V3plusV3(cm, dr, cm)
}

void CG_MMOLECULE::calMassCenter() {
	int nc, na;
	float m;

	M = 0;
	V3zero(cm)
	for (nc = 0; nc < nCluster; nc++) {
		for (na = 0; na < cluster[nc].nAtoms; na++) {
			m = cluster[nc].atom[na].par->m;
			cm.v[0] += m * cluster[nc].atom[na].r.v[0];
			cm.v[1] += m * cluster[nc].atom[na].r.v[1];
			cm.v[2] += m * cluster[nc].atom[na].r.v[2];
			M += m;
		}
	}
	cm.v[0] /= M; cm.v[1] /= M; cm.v[2] /= M;
}

void CG_MMOLECULE::SetMassCenter(double *v, bool cal_mass_center) {
	int i;

	if (cal_mass_center) calMassCenter();
	VECTOR3 dr;
	for (i = 0; i < 3; i++) dr.v[i] = v[i] - cm.v[i];
	shiftMol(dr);
}

// connect c2 to c1 with vcHinge [in c1]
// align vpHinge in c2 parallel to vcHinge in c1
bool ConnectClusters(CG_CLUSTER *c1, int iHinge1, CG_CLUSTER *c2, int iHinge2, char bound_type) {
	if (iHinge1 >= c1->nHinges) {
		sprintf(errmsg, "hinge index %d does not exist [HINGES : %d]", iHinge1, c1->nHinges); show_msg(errmsg); return false;
	}
	if (iHinge2 >= c2->nHinges) {
		sprintf(errmsg, "hinge index %d does not exist [HINGES : %d]", iHinge2, c2->nHinges); show_msg(errmsg); return false;
	}

	HINGE<_BASIC_CLUSTER<CGATOM> > *hinge1 = c1->hinge + iHinge1;
	HINGE<_BASIC_CLUSTER<CGATOM> > *hinge2 = c2->hinge + iHinge2;
	CGATOM *patom2 = c2->atom + hinge2->indx;
	CGATOM *patom1 = c1->atom + hinge1->indx;

	float hinge_length = 4;
	extern char cgatom_par_fname[256];
	extern bool find_bound_length(char *fname, char *atom1, char *atom2, float &blen);
	if (!find_bound_length(cgatom_par_fname, patom1->par->atom, patom2->par->atom, hinge_length)) 
		hinge_length = 4;

	c1->hinge[iHinge1].vHinge.stretch(hinge_length);
	c2->hinge[iHinge2].vHinge.stretch(hinge_length);

	VECTOR3 O2 = (*(c2->hinge[iHinge2].O));
	VECTOR3 Op2 = c2->hinge[iHinge2].vHinge + O2;
	VECTOR3 O1 = (*(c1->hinge[iHinge1].O));
	c2->align(Op2, O2, O1, c1->hinge[iHinge1].vHinge);
	//c2->hinge[iHinge2].vHinge = c1->hinge[iHinge1].vHinge * (-1);

	if (!c1->setup_link_relationship(iHinge1, c2, iHinge2, bound_type)) return false;

	return true;
}

// return energy in kT
void calCGClustersKineticEnergy(CG_MMOLECULE *mm, int n1, int n2, double& Ek) {
	Ek = 0;
	double unit_E = unit_Ek_kT * 0.5; // important, unit_E is divid by 2
	double Ekc = 0, dEkc = 0;
	double EkThermal = 3 * ::Ek0Internal, Ekcut = EkThermal * 3;
	int n = 0;
	CG_CLUSTER *pc = NULL;
	for (n = n1; n <= n2; n++) {
		if (n >= mm->nCluster) break;
		pc = mm->cluster + n;
		scale_uv3(pc->bkm->v, pc->bkm->v, Ekc) 
		Ekc *= pc->atom[0].par->m; // Ekc = M * v^2; assuming each cluster has only one atom inside
		Ekc *= unit_E;
		pc->bkm->Ek = Ekc;
		Ek += Ekc;

		//dEkc = Ekc - EkThermal;
		//pc->bkm->hoover = (dEkc > EkThermal ? true : false);
		pc->bkm->hoover = (Ekc > Ekcut ? true : false);
	}
	//return Ek;
}

// return energy in kT
void ParallelCalKineticEnergy(CG_MMOLECULE &mm) {
	double Ek = 0;
	int nThreads = MAX_THREADS;
	MacroMoleculeOperate2<CG_MMOLECULE, double>((void*)(&calCGClustersKineticEnergy), mm, 0, mm.nCluster, nThreads, Ek);
	mm.Ek = (float)Ek;
}


// calculate acceleration
void calClusterDynamic(CG_MMOLECULE *mm, int n1, int n2) {
	int n = 0;
	CG_CLUSTER *pc = NULL;
	double M = 0;
	for (n = n1; n <= n2; n++) {
		if (n >= mm->nCluster) break;
		pc = mm->cluster + n;
		V3plusV3(pc->dyn->fc, pc->dyn->f_Brown, pc->dyn->f)
		M = pc->atom[0].par->m;
		pc->bkm->alpha.v[0] = pc->dyn->f.v[0] / M;
		pc->bkm->alpha.v[1] = pc->dyn->f.v[1] / M;
		pc->bkm->alpha.v[2] = pc->dyn->f.v[2] / M;
	}
}

// calculate the acceleration of each cluster in mm
void ParallelCalDynamic(CG_MMOLECULE &mm) {
	int nThreads = MAX_THREADS;
	MacroMoleculeOperate<CG_MMOLECULE>((void*)(&calClusterDynamic), mm, 0, mm.nCluster, nThreads);
}

// Note : the angular momentum in this function is not real angular momentum
// here we calculate only r0 * m * v; -- r0 is from mass center to cluster
void calFrameSpatialMomentumThread(CG_MMOLECULE *mm, int n1, int n2, VECTOR<6> &P) {
	int n = 0;
	CG_CLUSTER *pc = NULL;
	VECTOR3 r0, w; // translation of mass center, total angular moment relative to mass center
	float m = 0;
	V6zero(P)
	for (n = n1; n <= n2; n++) {
		if (n >= mm->nCluster) break;
		pc = mm->cluster + n;
		m = *(pc->M);
		P.v[0] += m * pc->bkm->v.v[0];
		P.v[1] += m * pc->bkm->v.v[1];
		P.v[2] += m * pc->bkm->v.v[2];
		// Iw = r * m * v = (r0 + cm) * m * v = r0 * m * v + cm * m * v  --  here we calculate the first term only
		VECT3(mm->cm, (*(pc->r)), r0)
		V3PV3(r0, pc->bkm->v, w)
		P.v[3] += m * w.v[0];
		P.v[4] += m * w.v[1];
		P.v[5] += m * w.v[2];
	}
}

void calFrameSpatialMomentum(CG_MMOLECULE *mm) {
	V6zero(mm->mKD->P)
	MacroMoleculeOperate2< CG_MMOLECULE, VECTOR<6> >((void*)(&calFrameSpatialMomentumThread), *mm, 0, mm->nCluster, MAX_THREADS, mm->mKD->P);
	// another term of Iw -- from calFrameSpatialMomentumThread : cm * m * v
	VECTOR3 P0, vt;
	P0.v[0] = mm->mKD->P.v[0]; P0.v[1] = mm->mKD->P.v[1]; P0.v[2] = mm->mKD->P.v[2];
	V3PV3(mm->cm, P0, vt)
	mm->mKD->P.v[3] += vt.v[0]; mm->mKD->P.v[4] += vt.v[1]; mm->mKD->P.v[5] += vt.v[2];
	mm->mKD->v.v[0] = mm->mKD->P.v[0] / mm->M;
	mm->mKD->v.v[1] = mm->mKD->P.v[1] / mm->M;
	mm->mKD->v.v[2] = mm->mKD->P.v[2] / mm->M;

	P0.v[0] = mm->mKD->P.v[3];
	P0.v[1] = mm->mKD->P.v[4];
	P0.v[2] = mm->mKD->P.v[5];
	MpV3(mm->mKD->invM, P0, mm->mKD->w)
}

// Note : the angular momentum in this function is not real angular momentum
// here we calculate only r0 * m * v; -- r0 is from mass center to cluster
void calInertiaMomentThread(CG_MMOLECULE *mm, int n1, int n2, MATRIX<3> &M) {
	int n = 0, i, j;
	CG_CLUSTER *pc = NULL;
	VECTOR3 r0; // translation of mass center, total angular moment relative to mass center
	double R2 = 0;
	float m = 0;
	for (n = n1; n <= n2; n++) {
		if (n >= mm->nCluster) break;
		pc = mm->cluster + n;
		m = *(pc->M);
		// Iw = r * m * v = (r0 + cm) * m * v = r0 * m * v + cm * m * v  --  here we calculate the first term only
		VECT3(mm->cm, (*(pc->r)), r0) scale_uv3(r0, r0, R2)
		for (i = 0; i < 3; i++) {for (j = 0; j < 3; j++) {
			M.m[i][j] -= m * r0.v[i] * r0.v[j];
			if (i == j) M.m[i][j] += m * R2;
		}}
	}
}

void calInertiaMoment(CG_MMOLECULE *mm) {
	MacroMoleculeOperate2< CG_MMOLECULE, MATRIX<3> >((void*)(&calInertiaMomentThread), *mm, 0, mm->nCluster, MAX_THREADS, mm->mKD->M);
	InverseMatrix<3>(mm->mKD->M, mm->mKD->invM);
}

void calFrameEnergy(CG_MMOLECULE *mm) {
	// assuming that spatial momentum and velocity were calculated
	double v2 = 0;
	scale_uv3(mm->mKD->v, mm->mKD->v, v2) mm->mKD->Ekt = 0.5 * mm->M * v2 * ::unit_Ek_kT; // kT
	mm->mKD->Ekr = mm->mKD->w.v[0] * mm->mKD->P.v[3] + mm->mKD->w.v[1] * mm->mKD->P.v[4] + mm->mKD->w.v[2] * mm->mKD->P.v[5];
	mm->mKD->Ekr *= 0.5 * ::unit_Ek_kT; // kT
}


/*
bool construct_coarse_grained_mm(MMOLECULE &mm, CG_MMOLECULE &cgmm) {
	CG_CLUSTER *pcg = NULL;
	CLUSTER *pc = NULL;
	if (cgmm.cluster != NULL) cgmm.reset();

	char msg[256] = "\0";

	int ncluster = 0, nc = 0, ncg = 0, n, nhinge = 0;
	int nCG = 0;
	for (nc = 0; nc < mm.nCluster; nc++) {
		if (nCG < mm.cluster[nc].cgIndx) nCG = mm.cluster[nc].cgIndx;
	}
	nCG++; // cgIndx starts from 0
	cgmm.cluster = new CG_CLUSTER[nCG]; cgmm.nCluster = nCG;

	CHAIN<CLUSTER> *ch = NULL, *tch = NULL;
	HINGE<_BASIC_CLUSTER<MATOM> > *hinge = NULL;
	int nAtoms = 0, na = 0, na_start = 0;
	for (ncg = 0; ncg < nCG; ncg++) {
		ncluster = 0; nhinge = 0;
		release_chain<CLUSTER>(&ch, false);
		for (nc = 0; nc < mm.nCluster; nc++) {
			pc = mm.cluster + nc;
			if (pc->cgIndx == ncg) {
				if (ch == NULL) {ch = new CHAIN<CLUSTER>; ch->p = pc;}
				else ch->attach2tail(pc);
				ncluster++;
				for (n = 0; n < pc->nHinges; n++) {
					if (((BASIC_CLUSTER*)(pc->hinge[n].plink))->cgIndx != pc->cgIndx) nhinge++;
				}
			}
		}
		if (ncluster <= 0) {
			sprintf(msg, "no real cluster is labeled for coarse-grained cluster #%d", ncg);
			show_msg(msg);
			cgmm.reset(); return false;
		}
		pcg = cgmm.cluster + ncg;
		pcg->cluster.SetArray(ncluster);
		tch = ch;
		for (n = 0; n < ncluster; n++) {
			if (tch == NULL) {
				sprintf(msg, "STRANGE : number of cluster is not consistent when constructing coarse-grained cluster #%d", ncg);
				show_msg(msg);
				cgmm.reset(); return false;
			}
			pcg->cluster.m[n] = tch->p;
			tch = tch->next;
		}
		nAtoms = 0;
		for (n = 0; n < ncluster; n++) nAtoms += pcg->cluster.m[n]->nAtoms;
		pcg->atom.SetArray(nAtoms);
		na_start = 0;
		for (n = 0; n < ncluster; n++) {
			for (na = 0; na < pcg->cluster.m[n]->nAtoms; na++) {
				pcg->atom.m[na_start + na] = pcg->cluster.m[n]->atom + na;
			}
			na_start += pcg->cluster.m[n]->nAtoms;
		}
		if (nhinge <= 0) {
			sprintf(msg, "STRANGE : coarse-grained cluster #%d has no hinge", ncg);
			show_msg(msg);
			cgmm.reset(); return false;
		}
		else {
			pcg->set_hinge(nhinge);
			nhinge = 0;
			for (nc = 0; nc < pcg->cluster.n; nc++) {
				pc = pcg->cluster.m[nc];
				for (n = 0; n < pc->nHinges; n++) {
					if (((BASIC_CLUSTER*)pc->hinge[n].plink)->cgIndx != pc->cgIndx) {
						pcg->hinge[nhinge].O = pc->hinge[n].O;
						pcg->hinge[nhinge].Olink = pc->hinge[n].Olink;
						pcg->hinge[nhinge].plink = (cgmm.cluster + ((BASIC_CLUSTER*)pc->hinge[n].plink)->cgIndx);
						V32V3(pc->hinge[n].vHinge, pcg->hinge[nhinge].vHinge)
						// indx_link will be setup later
						// we ignore the indx, ibound 
						// in cg_cluster, the indx is different from that in cluster
						// because the atom index is different

						nhinge++;
					}
				}
			}
		}
	}

	// setup indx_link in the hinges
	HINGE<VIRTUAL_BASIC_CLUSTER> *cghinge = NULL;
	bool status = false;
	for (ncg = 0; ncg < cgmm.nCluster; ncg++) {
		pcg = cgmm.cluster + ncg;
		for (nhinge = 0; nhinge < pcg->nHinge; nhinge++) {
			cghinge = pcg->hinge + nhinge;
			status = false;
			for (n = 0; n < cghinge->plink->nHinge; n++) {
				if (cghinge->plink->hinge[n].plink == (VIRTUAL_BASIC_CLUSTER*)pcg) {
					cghinge->indx_link = n;
					cghinge->plink->hinge[n].indx_link = nhinge;
					status = true;
					break;
				}
			}
			if (!status) {
				sprintf(msg, "failure to setup link for hinge #%d of cg_cluster #%d", nhinge, ncg);
				show_msg(msg); cgmm.reset(); return false;
			}
		}
	}

	// base
	cgmm.base = cgmm.cluster + mm.base->cgIndx;
	pcg = cgmm.base;

	cgmm.base->Op = mm.base->Op; // very important
	cgmm.xaxis = mm.xaxis;
	cgmm.zaxis = mm.zaxis;
	return true;
}
*/

void cp_cgatom(CGATOM *d, CGATOM *s) {
	//strcpy(d->name, s->name);
	d->aindx = s->aindx;
	d->c = s->c; d->c0 = s->c0;
	d->r = s->r;
	d->r0 = s->r0;
	d->par = s->par;
	// copying an atom, information is not enough to form the bound connection
	d->set_bounds(s->nb);
	for (int i = 0; i < d->nb; i++) {d->ba[i] = NULL; d->bv[i] = s->bv[i]; d->btype[i] = s->btype[i];}
	// bound will be connected in cluster copying
	if (bImplicitSolvent) d->psolv  = s->psolv;
}

bool cp(CG_CLUSTER *d, CG_CLUSTER *s) {
	d->reset();
	if (!cp_basic_cluster<CGATOM>((_BASIC_CLUSTER<CGATOM>*)d, (_BASIC_CLUSTER<CGATOM>*)s, (void*)(&cp_cgatom))) return false;
	d->cID = s->cID; d->cTypeIndx = s->cTypeIndx;
	d->eNeutral = s->eNeutral;
	d->mIndx = s->mIndx;
	d->setup_cordinate();
	d->d = s->d;

	d->rc_solv = s->rc_solv;
	d->psolv = s->psolv;
	return true;
}


bool cp(CG_MMOLECULE *m2, CG_MMOLECULE* m1) {
	int nc, i;
	m2->reset();
	m2->mType = m1->mType; m2->mIndx = m1->mIndx; strcpy(m2->mol, m1->mol);
	if (m1->nCluster <= 0) return true;
	m2->cluster = new CG_CLUSTER[m1->nCluster];
	m2->nCluster = m1->nCluster;
	for (nc = 0; nc < m1->nCluster; nc++) cp(m2->cluster + nc, m1->cluster + nc);

	// setup the relationship
	int iHinge1 = 0, iHinge2 = 0, nc2;
	CG_CLUSTER *pc1 = NULL, *pc2 = NULL;
	for (nc = 0; nc < m2->nCluster; nc++) {
		pc1 = m2->cluster + nc;
		for (i = 0; i < m1->cluster[nc].nHinges; i++) {
			if (pc1->hinge[i].plink != NULL) continue; // the link was setup already
			if ((pc2 = (CG_CLUSTER*)(m1->cluster[nc].hinge[i].plink)) == NULL) continue; // the link in original molecule is free
			for (nc2 = 0; nc2 < m1->nCluster; nc2++) {if (m1->cluster + nc2 == pc2) break;}
			pc2 = m2->cluster + nc2;
			iHinge1 = i;
			iHinge2 = m1->cluster[nc].hinge[i].indx_link;
			if (!(pc1->setup_link_relationship(iHinge1, pc2, iHinge2, SINGLE_BOUND))) {
				sprintf(errmsg, "failure to setup link between cluster %d & %d", nc, nc2);
				show_msg(errmsg); strcpy(errmsg, "\0"); 
				m2->reset();
				return false;
			}
		}
	}

	return true;
}

bool atom_indx(CG_MMOLECULE &cgmm, CGATOM* pcgatom, int &nc, int &na, int &indx) {
	nc = 0; na = 0; indx = 0;
	CG_CLUSTER *pc = NULL;
	for (nc = 0; nc < cgmm.nCluster; nc++) {
		pc = cgmm.cluster + nc;
		for (na = 0; na < pc->nAtoms; na++) {
			if (pc->atom + na == pcgatom) return true;
			indx++;
		}
	}
	return false;
}

bool find_bound_length(char *fname, char *atom1, char *atom2, float &blen) {
	char msg[256] = "\0", buffer[256] = "\0", batom1[250] = "\0", batom2[250] = "\0";
	char *pars = NULL;
	float vb = 0, bl = 0;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open %s to search bount-length", fname);
		show_msg(msg); return false;
	}
	if (!search(in, "BOUND", buffer)) {
		sprintf(msg, "BOUND is not defined in %s", fname);
		show_msg(msg); in.close(); return false;
	}
	while (!in.eof()) {
		in.getline(buffer, 250);
		if (strlen(buffer) == 0) break;
		if (is_bound(buffer, atom1, atom2, &pars) && strlen(pars) > 0
			&& sscanf(pars, "%f %f", &vb, &bl) == 2) {
				blen = bl; in.close(); return true;
		}
	}
	sprintf(msg, "bound %s-%s is not given in %s", atom1, atom2, fname); show_log(msg, true);
	in.close(); return false;
}

bool find_FENE_bound_par(char *fname, char *atom1, char *atom2, FENE_BOUND &bpar) {
	char msg[256] = "\0", buffer[256] = "\0", batom1[250] = "\0", batom2[250] = "\0";
	char *pars = NULL;
	float vb = 0, bl = 0, d = 0;
	int n = 0;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open %s to search bound parameters", fname);
		show_msg(msg); return false;
	}
	if (!search(in, "BOUND", buffer)) {
		sprintf(msg, "BOUND is not defined in %s", fname);
		show_msg(msg); in.close(); return false;
	}
	while (!in.eof()) {
		in.getline(buffer, 250); if (strlen(buffer) == 0) break;
		if (is_bound(buffer, atom1, atom2, &pars) && strlen(pars) > 0
			&& sscanf(pars, "%f %f", &vb, &bl) == 2) {
				bpar.vb = vb * 4.176; // kCal/mol ==>kJ/mol
				bpar.blen = bl;
				in.close(); return true;
		}
	}
	sprintf(msg, "bound %s-%s is not given in %s", atom1, atom2, fname); show_log(msg, true);
	in.close(); return false;
}

void FENE_BOUND::construct_prof(EF_DPROF &efprof) {
	double unit = U_LJForce / U_MassAccel;
	double dstep = this->blen / 100;
	efprof.set(199, dstep); // 0 ~ 0.99 * blen
	double v = -0.5 * vb * blen * blen;
	double x = 0;
	for (int i = 0; i < 199; i++) {
		efprof.x[i] = i * dstep;
		x = efprof.x[i] / blen - 1; x = FABS(x);
		if (x > 0.999) x = 0.999;
		efprof.y[0][i] = v * log(1 - x * x) * U_eps / kT;  // kJ/mol ==> kT
		efprof.y[1][i] = vb * (blen - efprof.x[i]) / (1 - x * x) * unit;  // kJ/mol/A ==> Atomic Mass unit
	}
}

extern char cgatom_par_fname[256];
void attach_FENE_bound_db(CG_MMOLECULE &cgmm, CGBOUND_DB &cgbound_db) {
	unsigned int key = 0;
	int nc, nn;
	CG_CLUSTER *pc = NULL, *plink = NULL;
	char aindx[2];
	FENE_BOUND bpar;
	BOUND_PROF *bprof = NULL;
	bool status = true;

	for (nc = 0; nc < cgmm.nCluster; nc++) {
		pc = cgmm.cluster + nc; aindx[0] = pc->atom[0].aindx;
		for (nn = 0; nn < pc->neighbor.m[0].n; nn++) {
			if (pc->neighbor.m[0].m[nn] == NULL) continue;
			plink = pc->neighbor.m[0].m[nn];
			aindx[1] = plink->atom[0].aindx;
			key = construct_char2_key_pairwise(aindx);
			bprof = cgbound_db.search_chain(key);
			if (bprof != NULL) {pc->cgb.m[0].m[nn] = bprof; continue;}
			else {
				status = find_FENE_bound_par(cgatom_par_fname, pc->atom[0].par->atom, plink->atom[0].par->atom, bpar);
				if (status) {
					if (::d_14 < bpar.blen * 2 * 3) {::d_14 = bpar.blen * 2 * 3; ::d2_14 = ::d_14 * ::d_14;}
				}
				else continue;
				bprof = new BOUND_PROF;
				bprof->aindx[0] = pc->atom[0].aindx; bprof->aindx[1] = plink->atom[0].aindx;
				bprof->key = key;
				bpar.construct_prof(bprof->efprof);
				cgbound_db.attach(bprof);
				pc->cgb.m[0].m[nn] = bprof;
				{
					char fname[256] = "\0";
					sprintf(fname, "%s-%s-bound.dat", pc->atom[0].par->atom, plink->atom[0].par->atom);
					ofstream out;
					int i = 0;
					out.open(fname);
					for (i = 0; i < bprof->efprof.npts; i++) {
						out<<bprof->efprof.x[i]<<" "<<bprof->efprof.y[0][i]<<" "<<bprof->efprof.y[1][i]<<endl;
					}
					out.close();

					char msg[256] = "\0";
					sprintf(msg, "bound energy & force between %s & %s are saved in %s", pc->atom[0].par->atom, plink->atom[0].par->atom, fname);
					show_log(msg, true);
				}
			}
		}
	}
}

bool get_neighbor_level(char *fname, int &n_neighbors) {
	char buffer[256] = "\0", msg[256] = "\0";
	if (!read_item(fname, "NEIGHBOR-INTERACTION", buffer)) {
		sprintf(msg, "NEIGHBOR-INTERACTION is not defined in %s\nfailure to get neighbor level", fname);
		show_msg(msg, true); return false;
	}
	if (sscanf(buffer, "%d", &n_neighbors) != 1) {
		sprintf(msg, "failure to get neighbor level from %s in %s", buffer, fname);
		show_msg(msg, true); return false;
	}
	else return true;
}

bool find_neighbor_bound_file(char *db_fname, char *atom1, char *atom2, int nth_neighbor, char* bound_fname) {
	// important : neighbor is pairwise
	char msg[256] = "\0", buffer[256] = "\0", title[256] = "\0", batom1[250] = "\0", batom2[250] = "\0";
	char *bf = NULL;
	bool status = false;
	sprintf(title, "NEIGHBOR-1%d", nth_neighbor);
	ifstream in;
	in.open(db_fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open %s to search neighbor-interaction", db_fname);
		show_msg(msg); return false;
	}
	if (!search(in, title, buffer)) {
		sprintf(msg, "%s is not given in %s", title, db_fname);
		show_msg(msg); in.close(); return false;
	}
	while (!in.eof()) {
		in.getline(buffer, 250); if (strlen(buffer) == 0) break;
		if (sscanf(buffer, "%s", title) != 1) break;
		strcpy(batom1, title); bf = batom1; next_section(&bf, "-"); bf[0] = '\0'; if (strlen(batom1) == 0) break; kick_off_tail(batom1, " ");
		bf = title; if (!NextSection(&bf, "-")) break;
		if (sscanf(bf, "%s", batom2) != 1 || strlen(batom2) == 0) break;
		if ((strcmp(batom1, atom1) != 0 || strcmp(batom2, atom2) != 0) && (strcmp(batom2, atom1) != 0 || strcmp(batom1, atom2) != 0)) continue;
		if (sscanf(buffer, "%s %s", title, bound_fname) == 2 && strlen(bound_fname) > 0) {
			in.close(); return true;
		}
	}
	//sprintf(msg, "bound %s-%s is not given in %s", atom1, atom2, db_fname); show_log(msg, true);
	in.close(); return false;
}

bool read_efprof(char *fname, EF_DPROF &efprof) {
	efprof.release();
	float *x = NULL, *y = NULL;
	int npts = 0;
	if (!read_profile(fname, 0, &x, 0, &y, 1, npts) || npts < 3) return false;
	float dx = x[1] - x[0]; if (dx > 0.2f) dx = 0.2f;
	int efp_npts = int((x[npts-1] - x[0]) / dx + 0.5f);
	efprof.set(efp_npts, dx, x[0]);
	int i, j;
	for (i = 0; i < efprof.npts; i++) {
		if (efprof.x[i] <= x[0]) efprof.y[0][i] = y[0];
		else if (efprof.x[i] >= x[npts-1]) efprof.y[0][i] = y[npts-1];
		else {
			for (j = 0; j < npts-1; j++) {if (efprof.x[i] < x[j+1]) break;}
			efprof.y[0][i] = y[j] + (y[j+1] - y[j]) / (x[j+1] - x[j]) * (efprof.x[i] - x[j]);
		}
	}
	// force
	float fUnit = kT / U_InertiaMoment;
	for (i = 1; i < efprof.npts - 1; i++) {
		efprof.y[1][i] = -(efprof.y[0][i+1] - efprof.y[0][i-1]) / (dx + dx);
		efprof.y[1][i] *= fUnit;
	}
	efprof.y[1][0] = efprof.y[1][1];
	efprof.y[1][efprof.npts-1] = efprof.y[1][efprof.npts-2];
	return true;
}

void attach_neighbor_bound_db_from_file(CG_MMOLECULE &cgmm, CGBOUND_DB &cgbound_db) {
	// important : neighbor is pairwise
	unsigned int key = 0;
	int nc, n_neighbor = 0, nn = 0;
	CG_CLUSTER *pc = NULL, *plink = NULL;
	char aindx[3];
	char neighbor_bound_fname[256] = "\0";
	BOUND_PROF *bprof = NULL;
	bool status = true;
	char msg[256] = "\0";

	for (nc = 0; nc < cgmm.nCluster; nc++) {
		pc = cgmm.cluster + nc; aindx[0] = pc->atom[0].aindx;
		for (n_neighbor = 0; n_neighbor < pc->neighbor.n; n_neighbor++) {
			for (nn = 0; nn < pc->neighbor.m[n_neighbor].n; nn++) {
				if (pc->neighbor.m[n_neighbor].m[nn] == NULL) continue;
				plink = pc->neighbor.m[n_neighbor].m[nn];
				aindx[1] = plink->atom[0].aindx; aindx[2] = char(n_neighbor);
				key = construct_char3_key_12_pairwise(aindx); // important : neighbor is pairwise
				bprof = cgbound_db.search_chain(key);
				if (bprof != NULL) {pc->cgb.m[n_neighbor].m[nn] = bprof; continue;}
				else {
					status = find_neighbor_bound_file(cgatom_par_fname, pc->atom[0].par->atom, plink->atom[0].par->atom, n_neighbor + 2, neighbor_bound_fname);
					if (!status) {
						//sprintf(msg, "bound energy & force between %s & %s, neighbor1-%d is not given in %s", pc->atom[0].par->atom, plink->atom[0].par->atom, n_neighbor+2, cgatom_par_fname);
						//show_log(msg, true); 
						continue;
					}
					bprof = new BOUND_PROF;
					bprof->aindx[0] = pc->atom[0].aindx; bprof->aindx[1] = plink->atom[0].aindx;
					bprof->key = key;
					status = read_efprof(neighbor_bound_fname, bprof->efprof);
					if (status) {
						cgbound_db.attach(bprof);
						pc->cgb.m[n_neighbor].m[nn] = bprof; 
					
					{
						char fname[256] = "\0";
						sprintf(fname, "%s-%s-neighbor1%d.dat", pc->atom[0].par->atom, plink->atom[0].par->atom, n_neighbor+2);
						ofstream out;
						int i = 0;
						out.open(fname);
						for (i = 0; i < bprof->efprof.npts; i++) {
							out<<bprof->efprof.x[i]<<" "<<bprof->efprof.y[0][i]<<" "<<bprof->efprof.y[1][i]<<endl;
						}
						out.close();

						sprintf(msg, "bound energy & force between %s & %s, neighbor1-%d are saved in %s", pc->atom[0].par->atom, plink->atom[0].par->atom, n_neighbor+2, fname);
						show_log(msg, true);
					}
					}
					else {
						sprintf(msg, "bound energy & force between %s & %s, neighbor1-%d : file %s does not exist", pc->atom[0].par->atom, plink->atom[0].par->atom, n_neighbor+2, neighbor_bound_fname);
						show_log(msg, true);

						pc->cgb.m[n_neighbor].m[nn] = NULL; 
						delete bprof; bprof = NULL;
					}
				}
			}
		}
	}
}

void attach_bound_db(CG_MMOLECULE &cgmm, CGBOUND_DB &cgbound_db, int cgbound_type) {
	switch (cgbound_type) {
	case 0: // FENE
		attach_FENE_bound_db(cgmm, cgbound_db);
		break;
	case 1: // from file
		attach_neighbor_bound_db_from_file(cgmm, cgbound_db);
		break;
	default:
		return;
	}
}

void cgbound_interact(CG_MMOLECULE *cgmm, int nc1, int nc2, double& V) {
	V = 0;
	CG_CLUSTER *pc = NULL, *plink = NULL;
	int n_neighbor, nn, i;
	BOUND_PROF *bprof = NULL;
	double v = 0, f = 0, d = 0, t = 0;
	VECTOR3 dr;
	for (int nc = nc1; nc <= nc2; nc++) {
		if (nc >= cgmm->nCluster) break;
		pc = cgmm->cluster + nc;
		for (n_neighbor = 0; n_neighbor < pc->neighbor.n; n_neighbor++) {
			for (nn = 0; nn < pc->neighbor.m[n_neighbor].n; nn++) {
				if (pc->neighbor.m[n_neighbor].m[nn] == NULL) continue;
				plink = pc->neighbor.m[n_neighbor].m[nn];
				bprof = pc->cgb.m[n_neighbor].m[nn];
				if (bprof == NULL) continue;
				// important, use the cordinate in molecule, not in the real space
				// we are calculating the distance inside the molecule
				VECT3((*(plink->r)), (*(pc->r)), dr)
				d = sqrt(dr.v[0] * dr.v[0] + dr.v[1] * dr.v[1] + dr.v[2] * dr.v[2]);
				if (d > bprof->efprof.xt && n_neighbor > 0) continue;
				interpolate2(v, f, bprof->efprof, d, i, t)
				V += v * 0.5; // this a bound-interact
				if (d > 0.1) {
					pc->dyn->fc.v[0] += f * dr.v[0] / d;
					pc->dyn->fc.v[1] += f * dr.v[1] / d;
					pc->dyn->fc.v[2] += f * dr.v[2] / d;
				}
				else {
					rand_vect3(dr.v);
					pc->dyn->fc.v[0] += f * dr.v[0];
					pc->dyn->fc.v[1] += f * dr.v[1];
					pc->dyn->fc.v[2] += f * dr.v[2];
				}
			}
		}
	}
}

/******************************************************
			NON-NEIGHBOR INTERACTION
******************************************************/
bool find_non_neighbor_interact_file(char *db_fname, char *atom1, char *atom2, char* non_neighbor_fname) {
	// important : non-neighbor is pairwise
	char msg[256] = "\0", buffer[256] = "\0", title[256] = "\0", batom1[250] = "\0", batom2[250] = "\0";
	char *bf = NULL;
	bool status = false;
	strcpy(title, "NON-NEIGHBOR-INTERACTION");
	ifstream in;
	in.open(db_fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open %s to search neighbor-interaction", db_fname);
		show_msg(msg); return false;
	}
	if (!search(in, title, buffer)) {
		sprintf(msg, "%s is not given in %s", title, db_fname); show_msg(msg); 
		in.close(); return false;
	}
	while (!in.eof()) {
		in.getline(buffer, 250); if (strlen(buffer) == 0) break;
		if (sscanf(buffer, "%s", title) != 1) break;
		strcpy(batom1, title); bf = batom1; next_section(&bf, "-"); bf[0] = '\0'; if (strlen(batom1) == 0) break; kick_off_tail(batom1, " ");
		bf = title; if (!NextSection(&bf, "-")) break;
		if (sscanf(bf, "%s", batom2) != 1 || strlen(batom2) == 0) break;
		if ((strcmp(batom1, atom1) != 0 || strcmp(batom2, atom2) != 0) && (strcmp(batom2, atom1) != 0 || strcmp(batom1, atom2) != 0)) continue;
		if (sscanf(buffer, "%s %s", title, non_neighbor_fname) == 2 && strlen(non_neighbor_fname) > 0) {
			in.close(); return true;
		}
	}
	//sprintf(msg, "non-neighbor interaction of %s-%s is not given in %s", atom1, atom2, db_fname); show_log(msg, true);
	in.close(); return false;
}

void attach_non_neighbor_db(CGATOM_COMM_PARS_DATABASE &cgatompar_db, CGBOUND_DB &cg_non_neighbor_db) {
	// important : non-neighbor is pairwise
	unsigned int key = 0;
	int natoms = number(cgatompar_db.ch); if (natoms == 0) return;
	CHAIN<CGATOM_COMM_PARS> *cg_ch1 = NULL, *cg_ch2 = NULL;
	CGATOM_COMM_PARS *cgpar1 = NULL, *cgpar2 = NULL;
	char aindx[2];
	char non_neighbor_fname[256] = "\0";
	BOUND_PROF *bprof = NULL;
	bool status = true;
	char msg[256] = "\0";
	bool get_new_prof = false;

	cg_ch1 = cgatompar_db.ch;
	while (cg_ch1 != NULL) {
		cgpar1 = cg_ch1->p; aindx[0] = cgpar1->indx;
		cg_ch2 = cgatompar_db.ch;
		while (cg_ch2 != NULL) {
			cgpar2 = cg_ch2->p; aindx[1] = cgpar2->indx;
			key = construct_char2_key_pairwise(aindx); // important : non-neighbor is pairwise
			bprof = cg_non_neighbor_db.search_chain(key);
			if (bprof != NULL) {cg_ch2 = cg_ch2->next; continue;}
			else {
				status = find_non_neighbor_interact_file(cgatom_par_fname, cgpar1->atom, cgpar2->atom, non_neighbor_fname);
				if (!status) {
					//sprintf(msg, "non-neighbor energy & force between %s & %s, is not given in %s", cgpar1->atom, cgpar2->atom, cgatom_par_fname);
					//show_log(msg, true); 
					cg_ch2 = cg_ch2->next; continue;
				}
				bprof = new BOUND_PROF;
				bprof->aindx[0] = aindx[0]; bprof->aindx[1] = aindx[1];
				bprof->key = key;
				status = read_efprof(non_neighbor_fname, bprof->efprof);
				if (status) {
					cg_non_neighbor_db.attach(bprof);
					if (cgpar1->rcut < bprof->efprof.xt) cgpar1->rcut = bprof->efprof.xt;
					if (cgpar2->rcut < bprof->efprof.xt) cgpar2->rcut = bprof->efprof.xt;
					if (::rcut_LJ < bprof->efprof.xt) {
						::rcut_LJ = bprof->efprof.xt;
						::r2cut_LJ = ::rcut_LJ * ::rcut_LJ;
						get_new_prof = true;
					}
				{
					char fname[256] = "\0";
					sprintf(fname, "%s-%s-non-neighbor.dat", cgpar1->atom, cgpar2->atom);
					ofstream out;
					int i = 0;
					out.open(fname);
					for (i = 0; i < bprof->efprof.npts; i++) {
						out<<bprof->efprof.x[i]<<" "<<bprof->efprof.y[0][i]<<" "<<bprof->efprof.y[1][i]<<endl;
					}
					out.close();

					sprintf(msg, "non-neighbor energy & force between %s & %s are saved in %s", cgpar1->atom, cgpar2->atom, fname);
					show_log(msg, true);
				}
				}
				else {
					sprintf(msg, "non-neighbor energy & force between %s & %s : file %s does not exist", cgpar1->atom, cgpar2->atom, non_neighbor_fname);
					show_log(msg, true);

					delete bprof; bprof = NULL;
				}
			}
			cg_ch2 = cg_ch2->next;
		}
		cg_ch1 = cg_ch1->next;
	}
	if (get_new_prof) {
		sprintf(msg, "inter-mol interaction cut at %8.2f", ::rcut_LJ);
		show_log(msg, true);
	}
}

/******************************************************
		INTRA-MOLECULE NON-NEIGHBOR INTERACTION
******************************************************/
bool find_intramol_non_neighbor_interact_file(char *db_fname, char *atom1, char *atom2, char* non_neighbor_fname) {
	// important : non-neighbor is pairwise
	char msg[256] = "\0", buffer[256] = "\0", title[256] = "\0", batom1[250] = "\0", batom2[250] = "\0";
	char *bf = NULL;
	bool status = false;
	strcpy(title, "INTRA-MOLECULE-NON-NEIGHBOR-INTERACTION");
	ifstream in;
	in.open(db_fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open %s to search neighbor-interaction", db_fname);
		show_msg(msg); return false;
	}
	if (!search(in, title, buffer)) {
		sprintf(msg, "%s is not given in %s", title, db_fname); show_msg(msg); 
		in.close(); return false;
	}
	while (!in.eof()) {
		in.getline(buffer, 250); if (strlen(buffer) == 0) break;
		if (sscanf(buffer, "%s", title) != 1) break;
		strcpy(batom1, title); bf = batom1; next_section(&bf, "-"); bf[0] = '\0'; if (strlen(batom1) == 0) break; kick_off_tail(batom1, " ");
		bf = title; if (!NextSection(&bf, "-")) break;
		if (sscanf(bf, "%s", batom2) != 1 || strlen(batom2) == 0) break;
		if ((strcmp(batom1, atom1) != 0 || strcmp(batom2, atom2) != 0) && (strcmp(batom2, atom1) != 0 || strcmp(batom1, atom2) != 0)) continue;
		if (sscanf(buffer, "%s %s", title, non_neighbor_fname) == 2 && strlen(non_neighbor_fname) > 0) {
			in.close(); return true;
		}
	}
	//sprintf(msg, "non-neighbor interaction of %s-%s is not given in %s", atom1, atom2, db_fname); show_log(msg, true);
	in.close(); return false;
}

void attach_intramol_non_neighbor_db(CGATOM_COMM_PARS_DATABASE &cgatompar_db, CGBOUND_DB &cg_non_neighbor_db) {
	// important : non-neighbor is pairwise
	unsigned int key = 0;
	int natoms = number(cgatompar_db.ch); if (natoms == 0) return;
	CHAIN<CGATOM_COMM_PARS> *cg_ch1 = NULL, *cg_ch2 = NULL;
	CGATOM_COMM_PARS *cgpar1 = NULL, *cgpar2 = NULL;
	char aindx[2];
	char non_neighbor_fname[256] = "\0";
	BOUND_PROF *bprof = NULL;
	bool status = true;
	char msg[256] = "\0";
	bool get_new_prof = false;

	cg_ch1 = cgatompar_db.ch;
	while (cg_ch1 != NULL) {
		cgpar1 = cg_ch1->p; aindx[0] = cgpar1->indx;
		cg_ch2 = cgatompar_db.ch;
		while (cg_ch2 != NULL) {
			cgpar2 = cg_ch2->p; aindx[1] = cgpar2->indx;
			key = construct_char2_key_pairwise(aindx); // important : non-neighbor is pairwise
			bprof = cg_non_neighbor_db.search_chain(key);
			if (bprof != NULL) {cg_ch2 = cg_ch2->next; continue;}
			else {
				status = find_intramol_non_neighbor_interact_file(cgatom_par_fname, cgpar1->atom, cgpar2->atom, non_neighbor_fname);
				if (!status) {
					//sprintf(msg, "non-neighbor energy & force between %s & %s, is not given in %s", cgpar1->atom, cgpar2->atom, cgatom_par_fname);
					//show_log(msg, true); 
					cg_ch2 = cg_ch2->next; continue;
				}
				bprof = new BOUND_PROF;
				bprof->aindx[0] = aindx[0]; bprof->aindx[1] = aindx[1];
				bprof->key = key;
				status = read_efprof(non_neighbor_fname, bprof->efprof);
				if (status) {
					cg_non_neighbor_db.attach(bprof);
					if (cgpar1->rcut < bprof->efprof.xt) cgpar1->rcut = bprof->efprof.xt;
					if (cgpar2->rcut < bprof->efprof.xt) cgpar2->rcut = bprof->efprof.xt;
					if (::rcut_intramol < bprof->efprof.xt) {
						::rcut_intramol = bprof->efprof.xt;
						::r2cut_intramol = ::rcut_intramol * ::rcut_intramol;
						get_new_prof = true;
					}
				{
					char fname[256] = "\0";
					sprintf(fname, "%s-%s-intramolecule-non-neighbor.dat", cgpar1->atom, cgpar2->atom);
					ofstream out;
					int i = 0;
					out.open(fname);
					for (i = 0; i < bprof->efprof.npts; i++) {
						out<<bprof->efprof.x[i]<<" "<<bprof->efprof.y[0][i]<<" "<<bprof->efprof.y[1][i]<<endl;
					}
					out.close();

					sprintf(msg, "intra-molecule non-neighbor energy & force between %s & %s are saved in %s", cgpar1->atom, cgpar2->atom, fname);
					show_log(msg, true);
				}
				}
				else {
					sprintf(msg, "intra-molecule non-neighbor energy & force between %s & %s : file %s does not exist", cgpar1->atom, cgpar2->atom, non_neighbor_fname);
					show_log(msg, true);

					delete bprof; bprof = NULL;
				}
			}
			cg_ch2 = cg_ch2->next;
		}
		cg_ch1 = cg_ch1->next;
	}
	if (get_new_prof) {
		sprintf(msg, "intra-mol interaction cut at %8.2f", ::rcut_intramol);
		show_log(msg, true);
	}
}

// return energy in kT
double calCGSMKineticEnergy(CGSM &sm) {
	double Ek = 0;
	double unit_E = ::unit_Ek_kT * 0.5; // important, unit_E is divid by 2

	scale_uv3(sm.bkm->v, sm.bkm->v, Ek) 
	Ek *= *(sm.M); // Ek = M * v^2; assuming each cluster has only one atom inside
	Ek *= unit_E;
	sm.Ek = Ek;
	return Ek;
}

// calculate acceleration
void calCGSMDynamic(CGSM &sm) {
	double M = 0;
	V3plusV3(sm.dyn->fc, sm.dyn->f_Brown, sm.dyn->f)
	M = *(sm.M);
	sm.bkm->alpha.v[0] = sm.dyn->f.v[0] / M;
	sm.bkm->alpha.v[1] = sm.dyn->f.v[1] / M;
	sm.bkm->alpha.v[2] = sm.dyn->f.v[2] / M;
}

bool cp(CGSM *d, CGSM *s) {
	d->mType = s->mType; d->mIndx = s->mIndx; strcpy(d->mol, s->mol);
	CG_CLUSTER *dest = (CG_CLUSTER*)d;
	CG_CLUSTER *source = (CG_CLUSTER*)s;
	cp(dest, source);
	return true;
}

bool get_cgatom_pdb_alias(CGATOM_COMM_PARS *atompar, char *fname, char *title) {
	char msg[256] = "\0", atom[250] = "\0", alias[250] = "\0";
	if (fname == NULL || strlen(fname) == 0) {
		sprintf(msg, "input PDB-alias of %s : ", atompar->atom);
		show_msg(msg, false);
		input(msg, 250);
		if (sscanf(msg, "%s", alias) == 1) {strcpy(atompar->pdb_alias, alias); return true;}
		else return false;
	}
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open file %s", fname); show_log(msg, true); return false;
	}
	if (!search(in, title, msg)) {
		sprintf(msg, "%s is not defined in %s", title, fname); 
		show_log(msg, true); in.close(); return false;
	}
	while (!in.eof()) {
		in.getline(msg, 250);
		if (sscanf(msg, "%s %s", atom, alias) != 2) break;
		if (strcmp(atom, atompar->atom) == 0) {
			strcpy(atompar->pdb_alias, alias); in.close(); return true;
		}
	}
	in.close(); return false;
}

bool find_bound_dipole_distrbt_file(char *db_fname, char *atom1, char *atom2, char* mu_distrbt_fname, char* theta_distrbt_fname) {
	// important : dipole distribt is non-pairwise
	char msg[256] = "\0", buffer[256] = "\0", title[256] = "\0", batom1[250] = "\0", batom2[250] = "\0";
	char *bf = NULL;
	bool status = false;
	strcpy(title, "BOUND-DIPOLE-DISTRBT");
	ifstream in;
	in.open(db_fname);
	if (!in.is_open()) {
		sprintf(msg, "failure to open %s to search dipole-distrbt", db_fname);
		show_msg(msg); return false;
	}
	if (!search(in, title, buffer)) {
		sprintf(msg, "%s is not given in %s", title, db_fname); show_msg(msg); 
		in.close(); return false;
	}
	while (!in.eof()) {
		in.getline(buffer, 250); if (strlen(buffer) == 0) break;
		if (sscanf(buffer, "%s", title) != 1) break;
		strcpy(batom1, title); bf = batom1; next_section(&bf, "-"); bf[0] = '\0'; if (strlen(batom1) == 0) break; kick_off_tail(batom1, " ");
		bf = title; if (!NextSection(&bf, "-")) break;
		if (sscanf(bf, "%s", batom2) != 1 || strlen(batom2) == 0) break;
		//if ((strcmp(batom1, atom1) != 0 || strcmp(batom2, atom2) != 0) && (strcmp(batom2, atom1) != 0 || strcmp(batom1, atom2) != 0)) continue;
		if ((strcmp(batom1, atom1) != 0 || strcmp(batom2, atom2) != 0)) continue; // non-pairwise
		if (sscanf(buffer, "%s %s %s", title, mu_distrbt_fname, theta_distrbt_fname) == 3 && strlen(mu_distrbt_fname) > 0 && strlen(theta_distrbt_fname) > 0) {
			in.close(); return true;
		}
	}
	//sprintf(msg, "dipole-distribution of %s-%s is not given in %s", atom1, atom2, db_fname); show_log(msg, true);
	in.close(); return false;
}

// for single-chain polymer, this function will generate the bound
// direction 1==>2, 2==>, 3==>4 ...
// for non-single-chain macromolecule, new function is needed
void CG_MMOLECULE::setup_bounds() {
	if (bound != NULL) {delete[] bound; bound = NULL;}
	nBound = 0;

	int ic, ib;
	CG_CLUSTER *pc = NULL;

	if (nCluster <= 1) return;
	int max_bounds = nCluster * (nCluster - 1);
	_BOUND_<CG_CLUSTER> *b = new _BOUND_<CG_CLUSTER>[max_bounds];
	_BOUND_<CG_CLUSTER> *tb = NULL;
	int nb = 0, ntb = 0;

	char msg[256] = "\0";

	for (ic = 0; ic < nCluster; ic++) {
		pc = cluster + ic;
		for (ib = 0; ib < pc->nHinges; ib++) {
			if (pc->hinge[ib].plink == NULL) continue;
			if (nb < 20) {tb = b; ntb = nb;}
			else {tb = b + nb - 20; ntb = 20;}
			if (exist_bound< CG_CLUSTER, _BOUND_<CG_CLUSTER> >(pc, (CG_CLUSTER*)(pc->hinge[ib].plink), true, tb, ntb)) continue;
			if (nb >= max_bounds) {
				sprintf(msg, "ERROR: bounds is more than estimated number [%d]", max_bounds); show_msg(msg, true);
				delete[] b; b = NULL; return;
			}
			b[nb].a[0] = pc; b[nb].a[1] = (CG_CLUSTER*)(pc->hinge[ib].plink); nb++;
		}
	}
	nBound = nb;
	bound = new CGBOUND[nBound];
	char key = 0, aindx[2];
	CGBOUND_DIPOLE_DISTRBT *pdt = NULL;
	bool status = false;
	char mu_fname[256] = "\0", theta_fname[256] = "\0";
	float vmax = 0;
	for (ib = 0; ib < nb; ib++) {
		bound[ib].a[0] = b[ib].a[0];
		bound[ib].a[1] = b[ib].a[1];
		aindx[0] = bound[ib].a[0]->atom[0].aindx;
		aindx[1] = bound[ib].a[1]->atom[0].aindx;
		key = construct_char2_key_non_pairwise(aindx); // important : bound is non-pairwise
		pdt = cg_dipole_distrbt_db.search_chain(key);
		if (pdt != NULL) {bound[ib].dt = pdt; continue;}
		else {
			if (::bCGdipole) {
				status = find_bound_dipole_distrbt_file(cgatom_par_fname, bound[ib].a[0]->atom[0].par->atom, bound[ib].a[1]->atom[0].par->atom, mu_fname, theta_fname);
				if (!status) {
					//sprintf(msg, "dipole distribution %s-%s, is not given in %s", bound[ib].a[0]->par->atom, bound[ib].a[1]->par->atom, cgatom_par_fname);
					//show_log(msg, true); 
					bound[ib].dt = NULL; continue;
				}
				pdt = new CGBOUND_DIPOLE_DISTRBT;
				pdt->aindx[0] = aindx[0]; pdt->aindx[1] = aindx[1];
				pdt->key = key;
				status = read_profile(mu_fname, 0, &(pdt->mu.x), 0, &(pdt->mu.y), 1, pdt->mu.N);
				status = status & read_profile(theta_fname, 0, &(pdt->theta.x), 0, &(pdt->theta.y), 1, pdt->theta.N);
				if (status) {
					pdt->normalize_distrbt();
					cg_dipole_distrbt_db.attach(pdt);
					bound[ib].dt = pdt;

				{
					char fname[256] = "\0";
					sprintf(fname, "%s-%s-mu-distrbt.dat", bound[ib].a[0]->atom[0].par->atom, bound[ib].a[1]->atom[0].par->atom);
					ofstream mu_out, theta_out;
					int i = 0;
					mu_out.open(fname);
					for (i = 0; i < pdt->mu.N; i++) {
						mu_out<<pdt->mu.x[i]<<" "<<pdt->mu.y[i]<<endl;
					}
					mu_out.close();
					sprintf(msg, "|mu| distrbt %s-%s is saved in %s", bound[ib].a[0]->atom[0].par->atom, bound[ib].a[1]->atom[0].par->atom, fname);
					show_msg(msg, true);

					sprintf(fname, "%s-%s-theta-distrbt.dat", bound[ib].a[0]->atom[0].par->atom, bound[ib].a[1]->atom[0].par->atom);
					theta_out.open(fname);
					for (i = 0; i < pdt->theta.N; i++) {
						theta_out<<pdt->theta.x[i]<<" "<<pdt->theta.y[i]<<endl;
					}
					theta_out.close();
					sprintf(msg, "dipole-theta distrbt %s-%s is saved in %s", bound[ib].a[0]->atom[0].par->atom, bound[ib].a[1]->atom[0].par->atom, fname);
					show_log(msg, true);
				}
				}
				else {
					sprintf(msg, "%s==>%s dipole-distribt : mu -- %s / theta -- %s does not exist", bound[ib].a[0]->atom[0].par->atom, bound[ib].a[1]->atom[0].par->atom, mu_fname, theta_fname);
					show_log(msg, true);
					delete pdt; pdt = NULL;
					bound[ib].dt = NULL;
				}
			}
		}
	}
	delete[] b; b = NULL;
}

// create dipole of each bound
void CG_MMOLECULE::bound_dipole() {
	CGBOUND *pb = NULL;
	for (int ib = 0; ib < nBound; ib++) {
		pb = bound + ib;
		if (pb->dt == NULL) continue;
		pb->vmu = xrand<float, float>(pb->dt->mu);
		pb->theta = xrand<float, float>(pb->dt->theta) * fAngle2Radian;
		rand_vect3(pb->randv.v);
	}
}

void calDipoleEffect(CG_MMOLECULE *mm, int nb1, int nb2, double &V) {
	double fUnit_ext_mass = U_EextForce / U_MassAccel;

	CGBOUND *pb = NULL;
	VECTOR3 vb, axis, vdipole, torque, Ev(0, 0, ::Eext), force;
	R rmatrix;
	double v = 0, U = 0;
	for (int ib = nb1; ib <= nb2; ib++) {
		if (ib >= mm->nBound) break;
		pb = mm->bound + ib;
		if (pb->dt == NULL) continue;
		VECT3((*(pb->a[0]->r)), (*(pb->a[1]->r)), vb)
		cp(vb, pb->randv, axis); // axis = vb x pb->randv
		RMatrix(axis, pb->theta, rmatrix, true);
		MpV3(rmatrix, vb, vdipole)
		v = vdipole.abs(); 
		if (v == 0) {V3zero(pb->mu)}
		else {
			v = pb->vmu / v; 
			pb->mu.v[0] = vdipole.v[0] * v; 
			pb->mu.v[1] = vdipole.v[1] * v; 
			pb->mu.v[2] = vdipole.v[2] * v;
		}
		V3PV3(pb->mu, Ev, torque) // now we get torque = dipole x E
		// we assum that a virtual force on each end, in a direction -- torque x r
		V3PV3(torque, vb, force)
		v = force.abs();
		if (v < 1e-7) continue; // torque is along r, virtual force is 0
		force.v[0] /= v; force.v[1] /= v; force.v[2] /= v; // force direction
		v = torque.abs() / vb.abs();  // now we need to get |f|. torque = r(0-->1) x force
//NOTE: force is normal to vb, smallest force, assuming that electric field will not change the bound length
		v *= fUnit_ext_mass;
		force.v[0] *= v; force.v[1] *= v; force.v[2] *= v;
		
		pb->a[1]->dyn->fc.v[0] += force.v[0];
		pb->a[1]->dyn->fc.v[1] += force.v[1];
		pb->a[1]->dyn->fc.v[2] += force.v[2];

		pb->a[0]->dyn->fc.v[0] -= force.v[0];
		pb->a[0]->dyn->fc.v[1] -= force.v[1];
		pb->a[0]->dyn->fc.v[2] -= force.v[2];

		U -= pb->mu.v[2] * ::Eext;
	}
	V = U * U_EextFr / kT;
}

void cgmm_check_non_neighbor_sm(CG_MMOLECULE *cgmm, int ncf, int nct, float r2cut) {
	int nc, nc1;
	CG_CLUSTER *pc = NULL, *pc1 = NULL;
	VECTOR3 dr;
	double r2;
	CHAIN<CG_CLUSTER> *cgch = NULL, *tail = NULL, *tch = NULL;
	for (nc = ncf; nc <= nct; nc++) {
		if (nc >= cgmm->nCluster) break;
		pc = cgmm->cluster + nc;
		pc->non_neighbor_sm.release();
		for (nc1 = 0; nc1 < cgmm->nCluster; nc1++) {
			if (nc1 == nc) continue;
			pc1 = cgmm->cluster + nc1;
			if (pc->cBound(pc1) > 0) continue; // neighbor
			VECT3(pc->atom[0].r, pc1->atom[0].r, dr);
			r2 = dr.v[0] * dr.v[0]; if (r2 > r2cut) continue;
			r2 += dr.v[1] * dr.v[1]; if (r2 > r2cut) continue;
			r2 += dr.v[2] * dr.v[2]; if (r2 > r2cut) continue;
			tch = new CHAIN<CG_CLUSTER>; tch->p = pc1;
			if (cgch == NULL) cgch = tch;
			else tail->next = tch; 
			tail = tch;
		}
		chain_array<CG_CLUSTER>(&cgch, pc->non_neighbor_sm);
		release_chain<CG_CLUSTER>(&cgch, false); tail = NULL;
	}
}

void CG_MMOLECULE::check_non_neighbor_sm(float r2cut) {
	MacroMoleculeOperate1<CG_MMOLECULE, float>((void*)(&cgmm_check_non_neighbor_sm), *this, 0, nCluster, r2cut, MAX_THREADS);
}

} // end of namespace _coarse_grain_ 
