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
#include "EwaldSumVar.h"
#include "MM.h"
#include "cg-mm.h"
#include "Mol.h"
#include "CMM.h"
using namespace _cmm_3d_;
#include "cg-cmm.h"

#include "var.h"

namespace _coarse_grain_ {

/*
// for cubic cell only
void InitEwaldSumPars(float L) {
	double V = L * L * L;
	kappa = _KAPPA / L; kappa3 = kappa * kappa * kappa;
	r_Ewd = float(_ALPHA_EWALD) * L; rcut_Ewd = r_Ewd; r2cut_Ewd = rcut_Ewd * rcut_Ewd;
	inv_V = float((PI2 + PI2) / (V * 3));
	int i, j, k;
	double dh = PI2 / L, fexp = 0, kappa2 = kappa * kappa, c = 1;
	for (i = 0; i < NH; i++) {for (j = 0; j < NH; j++) {for (k = 0; k < NH; k++) {
		h[i][j][k].v[0] = (i - _NH) * dh; 
		h[i][j][k].v[1] = (j - _NH) * dh; 
		h[i][j][k].v[2] = (k - _NH) * dh;
		V3ABS2(h[i][j][k], H2[i][j][k]) //ABS2(h[i][j][k], ti, H2[i][j][k])
		c = H2[i][j][k] / kappa2 * 0.25;
		fexp = exp(-c);
		Ah[i][j][k] = fexp / H2[i][j][k] * inv_V * 3;
	}}}
	AhMin = float(AH_MIN / (L * PI));
}

float cluster_max_size(CMM_CELL3D<BASIC_CLUSTER> &cell) {
	float dmax = 0;
	int n = 0;
	CHAIN<BASIC_CLUSTER> *pch = NULL;
	for (n = 0; n < cell.nCells; n++) {
		pch = cell.cell[n].ch;
		while (pch != NULL) {
			dmax = (pch->p->d > dmax ? pch->p->d : dmax);
			pch = pch->next;
		}
	}
	return dmax;
}
*/
void cmm_init_distribute_cluster(CGSM *sm, int nSM, CMM_CELL3D<CG_CLUSTER> &cell, bool reset_cell_chain) {
	if (reset_cell_chain) cell.reset_cell_chain();

	int nc = 0, ncell = 0;
	int ni = 0, nj = 0, nk = 0;
	VECTOR<3> dr, rc;
	int vi = 0;

	CG_CLUSTER* pc = NULL;
	for (nc = 0; nc < nSM; nc++) {
		pc = (CG_CLUSTER*)(sm + nc);
		V32V3((*(pc->vp)), rc)  // use vp position

		CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
		//cell_n(cell, ni, nj, nk, ncell)
		//if (ncell >= cell.acell.n) ncell = cell.acell.n - 1;
		//cell.acell.m[ncell]->attach((CG_CLUSTER*)pc);
		cell.bcell.m[ni].m[nj].m[nk].attach((CG_CLUSTER*)pc);
	}
}

void cmm_init_distribute_cluster(CG_MMOLECULE &mm, CMM_CELL3D<CG_CLUSTER> &cell, bool reset_cell_chain) {
	if (reset_cell_chain) cell.reset_cell_chain();

	int nc = 0, ncell = 0;
	int ni = 0, nj = 0, nk = 0;
	VECTOR<3> dr, rc;
	int vi = 0;

	CG_CLUSTER* pc = NULL;
	for (nc = 0; nc < mm.nCluster; nc++) {
		pc = mm.cluster + nc;
		V32V3((*(pc->vp)), rc)  // use vp position

		CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
		//cell_n(cell, ni, nj, nk, ncell)
		//if (ncell >= cell.acell.n) ncell = cell.acell.n - 1;
		//cell.acell.m[ncell]->attach((CG_CLUSTER*)pc);
		cell.bcell.m[ni].m[nj].m[nk].attach((CG_CLUSTER*)pc);
	}
}
/*
void cmm_init_distribute_cluster(SMOLECULE* sm, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain) {
	if (reset_cell_chain) cell.reset_cell_chain();

	int nm = 0, nc = 0;
	int ni = 0, nj = 0, nk = 0;
	VECTOR<3> dr, rc;
	int vi = 0;

	BASIC_CLUSTER* pc = NULL;
	for (nm = 0; nm < n; nm++) {
		pc = sm[nm].c;
		V32V3((*(pc->vp)), rc)  // use vp position

		CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
		cell_n(cell, ni, nj, nk, nc)
		if (nc >= cell.acell.n) nc = cell.acell.n - 1;
		cell.acell.m[nc]->attach(pc);
	}
}

void calClusterElectroMultipole(CMM_CELL3D<CG_CLUSTER> &cell, bool bEwaldSum, int n1Cell, int n2Cell) {
	CHAIN<CG_CLUSTER>* pch = NULL;
	CG_CLUSTER *pc = NULL;

	for (int nCell = n1Cell; nCell <= n2Cell; nCell++) {
		if (nCell >= cell.acell.n) break;
		pch = cell.acell.m[nCell]->ch;
		while (pch != NULL) {
			pc = pch->p;
			if (pc->eNeutral) {pch = pch->next; continue;}
			pc->calElectroMultipole(bEwaldSum); // calculate multipole of the cluster
			pch = pch->next;
		}
	}
}

void CMM_CELL<CG_CLUSTER>::calElectroMultipole(bool bEachCluster, bool bEwaldSum) {
	int i = 0, j = 0, nc[3];
	CHAIN<CG_CLUSTER>* pch = NULL;
	CG_CLUSTER *pc = NULL;

	VECTOR<3> dr;

	ELECTRO_MULTIPOLE *mp = this->emp, *cmp = NULL;

	mp->q = 0;
	V3zero(mp->mu)
	Mzero(mp->Q, i, j)
	//CMzero(mp->O, i, j, k)

	for (i = 0; i < 3; i++) mp->r0.v[i] = this->r0.v[i];
	
	if (this->subCell == NULL) {
		pch = this->ch;
		while (pch != NULL) {
			pc = pch->p;
			if (pc->eNeutral) {pch = pch->next; continue;}

			if (bEachCluster) pc->calElectroMultipole(bEwaldSum); // calculate multipole of the cluster
			cmp = pc->emp;
			mp->q += pc->emp->q;
			for (i = 0; i < 3; i++) {
				dr.v[i] = cmp->r0.v[i] - mp->r0.v[i];
				mp->mu.v[i] += cmp->q * dr.v[i] + cmp->mu.v[i];
			}
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
				// we kicked off the term in Q, - 1/2 * r^2 * delta_i,j, since it is useless in electric field and force calculation
				mp->Q.m[i][j] += cmp->Q.m[i][j] + 1.5 * (cmp->mu.v[i] * dr.v[j] 
					+ cmp->mu.v[j] * dr.v[i] + cmp->q * dr.v[i] * dr.v[j]);
			}
			pch = pch->next;
		}
	}
	else {
		for (nc[0] = 0; nc[0] < this->subCell->n; nc[0]++) {for (nc[1] = 0; nc[1] < this->subCell->n; nc[1]++) {for (nc[2] = 0; nc[2] < this->subCell->n; nc[2]++) {
			// calculate multipole of subcell first
			this->subCell->cell[nc[0]][nc[1]][nc[2]].calElectroMultipole(bEachCluster, bEwaldSum);

			cmp = this->subCell->cell[nc[0]][nc[1]][nc[2]].emp;
			mp->q += cmp->q;
			for (i = 0; i < 3; i++) {
				dr.v[i] = cmp->r0.v[i] - mp->r0.v[i];
				mp->mu.v[i] += cmp->q * dr.v[i] + cmp->mu.v[i];
			}
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
				// we kicked off the term in Q, - 1/2 * r^2 * delta_i,j, since it is useless in electric field and force calculation
				mp->Q.m[i][j] += cmp->Q.m[i][j] + 1.5 * (cmp->mu.v[i] * dr.v[j] 
					+ cmp->mu.v[j] * dr.v[i] + cmp->q * dr.v[i] * dr.v[j]);
			}
		}}}
	}

	this->eNeutral = mp->neutral();
	if (this->eNeutral) return; // no charge or multipole 

	// the terms for Ewald Sum
	double hR = 0, t1 = 0, t2 = 0;
	if (bEwaldSum) {
		for (nc[0] = 0; nc[0] < NH; nc[0]++) {for (nc[1] = 0; nc[1] < NH; nc[1]++) {for (nc[2] = 0; nc[2] < NH; nc[2]++) {
			scale_uv3(mp->r0, ::h[nc[0]][nc[1]][nc[2]], hR)
			PERIOD(hR, 0, PI2, t1) COS(mp->fcos.m[nc[0]][nc[1]][nc[2]], t1, i, t1) SIN(mp->fsin.m[nc[0]][nc[1]][nc[2]], t1, i, t2)
			scale_uv3(mp->mu, ::h[nc[0]][nc[1]][nc[2]], mp->muh.m[nc[0]][nc[1]][nc[2]])
			scale_VMV(::h[nc[0]][nc[1]][nc[2]], mp->Q, ::h[nc[0]][nc[1]][nc[2]], mp->Qh.m[nc[0]][nc[1]][nc[2]], i, j, 3)
		}}}
	}
}
*/
/************************************************************************************
                        ELECTROSTATIC & LJ relationships
************************************************************************************/
/*
#define _N_MULTP  16 // 36 // 16
void CheckCellRelship(CMM_CELL3D<CG_CLUSTER> &cell, VECTOR3 &r0, float radius, CMM_CELL<CG_CLUSTER> *pcell, float *hw, bool imagine_cell, CMM_RELSHIP &rship, LJ_RELSHIP &LJrship, bool coarse_grain, int mIndx, int cgIndx, bool mandatory_divid_cell = false) {
	if (pcell->ch == NULL) return; // this cell is empty

	int i, j, k, nx, ny, nz;
	VECTOR3 dr, u, dR, r1;
	float d2, R2, rcut_LJ, r2cut_LJ;
	float t1, t2, t3;
	CMM_CELL<BASIC_CLUSTER> *pcl = NULL;
	d2 = hw[0] * hw[0] + hw[1] * hw[1] + hw[2] * hw[2];
	dr.v[0] = pcell->r0.v[0] - r0.v[0]; dr.v[1] = pcell->r0.v[1] - r0.v[1]; dr.v[2] = pcell->r0.v[2] - r0.v[2];
	if (imagine_cell) {PeriodCellCloestDist2(dr, cell, R2, u, nx, ny, nz)}
	else {V3ABS2(dr, R2) nx = 0; ny = 0; nz = 0;}

	CHAIN<BASIC_CLUSTER> *pch1 = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	
	if (R2 < d2 * _N_MULTP || mandatory_divid_cell) { // this cell is too close to the atom
		if (pcell->subCell == NULL) {
			// this the deepest level
			pch1 = pcell->ch;
			while (pch1 != NULL) {
				pc1 = pch1->p;
				if (!imagine_cell && coarse_grain && pc1->mIndx == mIndx && pc1->cgIndx == cgIndx) {pch1 = pch1->next; continue;} // same molecule and same cgIndx

				//d2 = pc1->d * pc1->d / 4;
				d2 = pc1->d + radius + radius;
				d2 = d2 * d2;
				if (pc1->eNeutral) {
					//V3plusV3(pc1->atom[0].r, pc1->dr_cmm, r1)
					V32V3((*(pc1->vp)), r1)
				}
				else {V32V3(pc1->emp->r0, r1)} 
				VECT3(r0, r1, dr) // very important for LJ to check the image cell
				if (imagine_cell) {
					PeriodCellCloestDist2(dr, cell, R2, u, nx, ny, nz)
				}
				else {
					//ABS2(dr, i, R2) 
					scale_uv3(dr, dr, R2)
					nx = 0; ny = 0; nz = 0;
				}
				if (!pc1->eNeutral) {
					if (R2 < d2 * _N_MULTP) rship.attach(pc1); // this cluster is too close to the atom, and has electrostatic interaction
					else rship.attach(pc1->emp); // use the multipole method on this cluster
				}
				// this cluster is in LJ interaction region ?
				rcut_LJ = ::rcut_LJ + pc1->d * 0.5 + radius; r2cut_LJ = rcut_LJ * rcut_LJ;
				if (imagine_cell) {
					if (R2 < r2cut_LJ && R2 > _R2min_Interact) LJrship.attach(pc1, nx, ny, nz); // this cluster in [i, j, k] image cell is in LJ interaction region
				}
				else {if (R2 < r2cut_LJ && R2 > _R2min_Interact) LJrship.attach(pc1, 0, 0, 0);}

				pch1 = pch1->next;
			}
		}
		else {
			for (i = 0; i < pcell->subCell->n; i++) for (j = 0; j < pcell->subCell->n; j++) for (k = 0; k < pcell->subCell->n; k++) {
				pcl = &(pcell->subCell->cell[i][j][k]);
				if (pcl->ch != NULL) CheckCellRelship(cell, r0, radius, pcl, pcell->subCell->hw, imagine_cell, rship, LJrship, coarse_grain, mIndx, cgIndx, false);
			}
		}
	}
	else { // cell is far away from r0
		if (!imagine_cell && coarse_grain && pcell->ch != NULL) {
			pch1 = pcell->ch;
			while (pch1 != NULL) {
				pc1 = pch1->p;
				if (pc1->mIndx == mIndx && pc1->cgIndx == cgIndx) {// same molecule and same cgIndx
					//one cluster in this cell is in the same coarse-grain
					if (pcell->ch != NULL) CheckCellRelship(cell, r0, radius, pcell, hw, imagine_cell, rship, LJrship, coarse_grain, mIndx, cgIndx, true); // mandatory recheck the cluster in the cell 
					return; // over
				}
				pch1 = pch1->next;
			}
		}
		if (pcell->ch != NULL && !pcell->eNeutral) rship.attach(pcell->emp);
		// cluster inside the cell can not has LJ interaction with r0, because it is too far away
	}
}

// check the relationship of cluster at r0 
void CheckCellRelship(CMM_CELL3D<BASIC_CLUSTER> &cell, VECTOR3 &r0, float radius, bool imagine_cell, CMM_RELSHIP &rship, LJ_RELSHIP &LJrship, bool coarse_grain, int mIndx, int cgIndx) {
	rship.release_relationship();
	LJrship.release_relationship();
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcl = NULL;
	for (nc = 0; nc < cell.nCells; nc++) {
		pcl = cell.cell + nc;
		if (pcl->ch != NULL) CheckCellRelship(cell, r0, radius, pcl, cell.hw, imagine_cell, rship, LJrship, coarse_grain, mIndx, cgIndx, false);
	}
}
*/
/************************************************************************************
                        Local relationship only
	// imethod == 0 : LJ interaction for all clusters, ignoring the neighbors
	// imethod == 1 : LJ interaction for the clusters in different neighbors
************************************************************************************/
void CheckCell_LocalRelship(CMM_CELL3D<CG_CLUSTER> &cell, CG_CLUSTER *pc, LOCAL_RELSHIP *pship, int nc1, int nc2, int nc3, int imethod) {
	if (bImplicitSolvent) {
		show_msg("neighbors / clusters in same molecule, are ignored for local interaction. Implicit solvation effect is not considered properly.");
	}
	int i, j, k, nx, ny, nz, d1cmm, d2cmm, d3cmm;
	VECTOR3 dr, r0, r1;
	float d2, r2cut;

	V32V3((*(pc->vp)), r0)

	CMM_CELL<CG_CLUSTER> *pcell = NULL;
	CG_CLUSTER *pc2 = NULL;
	int dnx, dny, dnz;
	float xd[3] = {cell.xd[0] / cell.bcell.nx, cell.xd[1] / cell.bcell.ny, cell.xd[2] / cell.bcell.nz};
	dnx = int(::rcut_LJ / xd[0] + 0.5f) + 1; 
	dny = int(::rcut_LJ / xd[1] + 0.5f) + 1; 
	dnz = int(::rcut_LJ / xd[2] + 0.5f) + 1; 

	ITERATOR<CG_CLUSTER, CMM_CELL<CG_CLUSTER> > it;

	for (i = -dnx; i <= dnx; i++) {
		nx = nc1 + i; d1cmm = 0;
		if (nx < 0) {
			if (cell.periodic_cell) {while (nx < 0) {nx += cell.bcell.nx; d1cmm -= 1;}}
			else continue;
		}
		else if (nx >= cell.bcell.nx) {
			if (cell.periodic_cell) {while (nx >= cell.bcell.nx) {nx -= cell.bcell.nx; d1cmm += 1;}}
			else continue;
		}

		for (j = -dny; j <= dny; j++) {
			ny = nc2 + j; d2cmm = 0;
			if (ny < 0) {
				if (cell.periodic_cell) {while (ny < 0) {ny += cell.bcell.ny; d2cmm -= 1;}}
				else continue;
			}
			else if (ny >= cell.bcell.ny) {
				if (cell.periodic_cell) {while (ny >= cell.bcell.ny) {ny -= cell.bcell.ny; d2cmm += 1;}}
				else continue;
			}

			for (k = -dnz; k <= dnz; k++) {
				nz = nc3 + k; d3cmm = 0;
				if (nz < 0) {
					if (cell.periodic_cell) {while (nz < 0) {nz += cell.bcell.nz; d3cmm -= 1;}}
					else continue;
				}
				else if (nz >= cell.bcell.nz) {
					if (cell.periodic_cell) {while (nz >= cell.bcell.nz) {nz -= cell.bcell.nz; d3cmm += 1;}}
					else continue;
				}

				// check the relationship within the cell
				pcell = &(cell.bcell.m[nx].m[ny].m[nz]); it.set(pcell); it.start_iterate();
				while (!it.eoi()) {
					pc2 = it.current(); if (pc2 == NULL) {it.next_iterate(); continue;}

					if (imethod == 1) { // ignore the cluster in same macromolecule
						if (pc2->mIndx == pc->mIndx) {
							it.next_iterate();
							continue;
						} 
					}

					if (pc2 == pc && d1cmm == 0 && d2cmm == 0 && d3cmm == 0) {
						it.next_iterate();
						continue;
					}
					if (cell.periodic_cell) {
						r1.v[0] = pc2->vp->v[0] + d1cmm * cell.xd[0];
						r1.v[1] = pc2->vp->v[1] + d2cmm * cell.xd[1];
						r1.v[2] = pc2->vp->v[2] + d3cmm * cell.xd[2];
					}
					VECT3(r0, r1, dr) scale_uv3(dr, dr, d2) // dr = r1 - r0, d2 = |dr|^2
					switch (imethod) {
					case 0: // LJ
						r2cut = (pc->atom[0].par->rcut > pc2->atom[0].par->rcut ? pc->atom[0].par->rcut : pc2->atom[0].par->rcut);
						break;
					case 1: // LJ
						r2cut = (pc->atom[0].par->rcut > pc2->atom[0].par->rcut ? pc->atom[0].par->rcut : pc2->atom[0].par->rcut);
						break;
					default: // or case 0 -- LJ
						r2cut = (pc->atom[0].par->rcut > pc2->atom[0].par->rcut ? pc->atom[0].par->rcut : pc2->atom[0].par->rcut);
						break;
					}
					r2cut = r2cut * r2cut;
					if (d2 <= r2cut) {
						// important : for coarse-grained cluster, the neighbor interaction
						// is treated specially. The local-relationship is specially for
						// the non-neighbor interaction
						if (imethod == 0) {if (pc->cBound(pc2) < 0) pship->attach(pc2, d1cmm, d2cmm, d3cmm);}
						else if (imethod == 1) pship->attach(pc2, d1cmm, d2cmm, d3cmm);
					}
					it.next_iterate();
				}
	}}}
}

/*****************************************************************************
               FUNCTIONS TO CHECK RELATIONSHIP
******************************************************************************/
// check the relationship of all clusters in CELL [nCell1 ~ nCell2], with imagine cells
void CheckLocalRelshipCells(CMM_CELL3D<CG_CLUSTER> &cell, int imethod_rcut, int nCell1, int nCell2) {
	int nc, i, j;
	CMM_CELL<CG_CLUSTER> *pcell = NULL;
	CG_CLUSTER *pc = NULL;
	ITERATOR<CG_CLUSTER, CMM_CELL<CG_CLUSTER> > it;
	for (nc = nCell1; nc <= nCell2; nc++) {
		if (nc >= cell.bcell.nx) break;
		for (i = 0; i < cell.bcell.ny; i++) {for (j = 0; j < cell.bcell.nz; j++) {
			pcell = &(cell.bcell.m[nc].m[i].m[j]); it.set(pcell); it.start_iterate();
			while (!it.eoi()) {
				pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
				pc->LRship.release_relationship();
				CheckCell_LocalRelship(cell, pc, &(pc->LRship), nc, i, j, imethod_rcut);
				it.next_iterate();
			}
		}}
	}
}

// check the relationship of all clusters in WHOLE CELL
void CheckLocalRelship(CMM_CELL3D<CG_CLUSTER> &cell, int n1Cell, int n2Cell, int imethod_rcut, int nThreads) {
	return CellOperate1<CG_CLUSTER, int>((void*)(&CheckLocalRelshipCells), cell, n1Cell, n2Cell, imethod_rcut, nThreads);
}

void CGMMClusterImpSolvNeighborsInCell(CMM_CELL3D<CG_CLUSTER> *cell, CG_MMOLECULE *cgmm, int nc1, int nc2) {
	return _cmm_3d_::MacroMolClusterImpSolvNeighborsInCell<CG_MMOLECULE, CG_CLUSTER>(cell, cgmm, nc1, nc2);
};

/***************************************************************************
******      SPECIALLY USING FOR SINGLE MACRO-MOLECULE SIMULATION      ******
****************************************************************************/

void CheckCellLocalRelship(CMM_CELL3D<CG_CLUSTER> &cell, CG_CLUSTER *pc, int imethod) {
	VECTOR3 r0;
	V32V3((*(pc->vp)), r0)
	int n1 = int((r0.v[0] - cell.xl[0]) / (cell.xd[0] / cell.bcell.nx) + 0.5f);
	int n2 = int((r0.v[1] - cell.xl[1]) / (cell.xd[1] / cell.bcell.ny) + 0.5f);
	int n3 = int((r0.v[2] - cell.xl[2]) / (cell.xd[2] / cell.bcell.nz) + 0.5f);
	if (n1 < 0) n1 = 0;
	else if (n1 >= cell.bcell.nx) n1 = cell.bcell.nx - 1;
	if (n2 < 0) n2 = 0;
	else if (n2 >= cell.bcell.ny) n2 = cell.bcell.ny - 1;
	if (n3 < 0) n3 = 0;
	else if (n3 >= cell.bcell.nz) n3 = cell.bcell.nz - 1;

	pc->LRship.release_relationship();

	CheckCell_LocalRelship(cell, pc, &(pc->LRship), n1, n2, n3, imethod);
}

void CGClusterLocalRelshipInCell(CMM_CELL3D<CG_CLUSTER> *cell, CG_MMOLECULE *mm, int nc1, int nc2, int& imethod) {
	CG_CLUSTER *pc = NULL;
	for (int nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = mm->cluster + nc;
		CheckCellLocalRelship(*cell, pc, imethod);
	}
}

// check the relationship of all clusters in WHOLE CELL
void CGMMCheckLocalRelshipInCell(CMM_CELL3D<CG_CLUSTER> &cell, CG_MMOLECULE &mm, int nc1, int nc2, int imethod, int nThreads) {
	return MacroMolCellOperate1<CG_MMOLECULE, CG_CLUSTER, int>((void*)(&CGClusterLocalRelshipInCell), cell, mm, nc1, nc2, imethod, nThreads);
}

void CGSMLocalRelshipInCell(CMM_CELL3D<CG_CLUSTER> *cell, CGSM *sm, int nSM) {
	for (int nc = 0; nc < nSM; nc++) CheckCellLocalRelship(*cell, (CG_CLUSTER*)(sm + nc), 0);
}

// check the relationship of all clusters in WHOLE CELL
void CGSMCheckLocalRelshipInCell(CMM_CELL3D<CG_CLUSTER> &cell, CGSM *sm, int nSM, int nThreads) {
	return SMCellOperate<CGSM, CG_CLUSTER>((void*)(&CGSMLocalRelshipInCell), cell, sm, nSM, nThreads);
}

} // end of namespace _coarse_grain_ 


