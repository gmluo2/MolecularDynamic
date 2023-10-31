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
#include "MM.h"
#include "Mol.h"
#include "EwaldSumVar.h"
#include "CMM.h"

#include "var.h"

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif

namespace _EwaldSum_recp_ {
#if _EWALD_SUM == _MultiPole_EWALD_SUM
// for cubic cell only
void InitEwaldSumPars(float Lx, float Ly, float Lz) {
	double V = Lx * Ly * Lz;
	r_Ewd = float(_ALPHA_EWALD) * L; rcut_Ewd = r_Ewd; r2cut_Ewd = rcut_Ewd * rcut_Ewd;
	kappa = _KAPPA / r_Ewd;
	inv_V = float((PI2 + PI2) / (V * 3));
	kappa3 = kappa * kappa * kappa;
	int i, j, k;
	double dh = PI2 / L, fexp = 0, kappa2 = kappa * kappa, c = 1;
	for (i = 0; i < NHX; i++) {for (j = 0; j < NHY; j++) {for (k = 0; k < NHZ; k++) {
		h[i][j][k].v[0] = (i - _NHX) * dh; 
		h[i][j][k].v[1] = (j - _NHY) * dh; 
		h[i][j][k].v[2] = (k - _NHZ) * dh;
		V3ABS2(h[i][j][k], H2[i][j][k]) //ABS2(h[i][j][k], ti, H2[i][j][k])
		c = H2[i][j][k] / kappa2 * 0.25;
		fexp = exp(-c);
		Ah[i][j][k] = fexp / H2[i][j][k] * inv_V * 3;
	}}}
	AhMin = float(AH_MIN / (L * PI));
}
#elif _EWALD_SUM == _SPME_EWALD_SUM
void InitEwaldSumPars(float V, double &inv_3V, double &kappa, double &kappa3) {
	inv_3V = float((PI2 + PI2) / (V * 3));
	kappa = _KAPPA / ::rcut_Ewd;
	kappa3 = kappa * kappa * kappa;
	r2cut_Ewd = rcut_Ewd * rcut_Ewd;
}

#endif
} // end of namespace _EwaldSum_recp_ 

float cluster_max_size(CMM_CELL3D<BASIC_CLUSTER> &cell) {
	float dmax = 0;
	int n = 0;
	BASIC_CLUSTER *pc = NULL;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > iterator;
	for (n = 0; n < cell.acell.n; n++) {
		iterator.set(cell.acell.m[n]); iterator.start_iterate();
		while (!iterator.eoi()) {
			pc = iterator.current();
			dmax = (pc->d > dmax ? pc->d : dmax);
			iterator.next_iterate();
		}
	}
	return dmax;
}

void cmm_init_distribute_cluster(MMOLECULE &mm, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain) {
	if (reset_cell_chain) cell.reset_cell_chain();

	int nc = 0;
	int ni = 0, nj = 0, nk = 0;
	VECTOR<3> dr, rc;
	int vi = 0;

	LIST< LIST<CLUSTER> > *ptl = mm.ctree.tl;
	LIST<CLUSTER> *plist = NULL;
	CLUSTER* pc = NULL;
	while (ptl != NULL) {
		plist = ptl->p;
		while (plist != NULL) {
			pc = plist->p;
			V32V3((*(pc->vp)), rc)  // use vp position

			CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
			//cell_n(cell, ni, nj, nk, nc)
			//if (nc >= cell.nCells) nc = cell.nCells - 1;
			//cell.cell[nc].attach((BASIC_CLUSTER*)pc);
			cell.bcell.m[ni].m[nj].m[nk].attach((BASIC_CLUSTER*)pc);

			plist = plist->next;
		}
		ptl = ptl->next;
	}
}

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
		//cell_n(cell, ni, nj, nk, nc)
		//if (nc >= cell.nCells) nc = cell.nCells - 1;
		//cell.cell[nc].attach(pc);
		cell.bcell.m[ni].m[nj].m[nk].attach(pc);
	}
}

void cmm_init_distribute_cluster(PMOLECULE* pm, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain) {
	if (reset_cell_chain) cell.reset_cell_chain();

	int nm = 0, nc = 0;
	int ni = 0, nj = 0, nk = 0;
	VECTOR<3> dr, rc;
	int vi = 0;

	BASIC_CLUSTER* pc = NULL;
	for (nm = 0; nm < n; nm++) {
		pc = pm[nm].c;
		V32V3((*(pc->vp)), rc)  // use vp position

		CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
		//cell_n(cell, ni, nj, nk, nc)
		//if (nc >= cell.nCells) nc = cell.nCells - 1;
		//cell.cell[nc].attach(pc);
		cell.bcell.m[ni].m[nj].m[nk].attach(pc);
	}
}

#if _EWALD_SUM == _MultiPole_EWALD_SUM
void calClusterElectroMultipole(CMM_CELL3D<BASIC_CLUSTER> &cell, bool bEwaldSum, int n1Cell, int n2Cell) {
	CHAIN<BASIC_CLUSTER>* pch = NULL;
	BASIC_CLUSTER *pc = NULL;

	for (int nCell = n1Cell; nCell <= n2Cell; nCell++) {
		if (nCell >= cell.acell.n) break;
		pch = cell.acell.m[nCell].ch;
		while (pch != NULL) {
			pc = pch->p;
			if (pc->eNeutral) {pch = pch->next; continue;}
			pc->calElectroMultipole(bEwaldSum); // calculate multipole of the cluster
			pch = pch->next;
		}
	}
}
#endif

#define _N_MULTP  16 // 36 // 16

#if _EWALD_SUM == _MultiPole_EWALD_SUM
/************************************************************************************
                        ELECTROSTATIC & LJ relationships
************************************************************************************/

void CheckClusterRelship(CMM_CELL3D<BASIC_CLUSTER> &cell, VECTOR3 &r0, float radius, CMM_CELL<BASIC_CLUSTER> *pcell, float *hw, CMM_RELSHIP &rship, LJ_RELSHIP &LJrship, bool mandatory_divid_cell = false) {
	if (pcell->ch == NULL) return; // this cell is empty

	bool imagine_cell = cell.periodic_cell;

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
		//if (pcell->subCell == NULL) {
		// this the deepest level
		pch1 = pcell->ch;
		while (pch1 != NULL) {
			pc1 = pch1->p;

			//d2 = pc1->d * pc1->d / 4;
			d2 = pc1->d + radius + radius;
			d2 = d2 * d2;
			if (pc1->eNeutral) {
				//V3plusV3(pc1->atom[0].r, pc1->dr_cmm, r1)
				V32V3((*(pc1->vp)), r1)
			}
			else {V32V3(pc1->emp->r0, r1)} 
			VECT3(r0, r1, dr) // very important for LJ to check the image cell

			nx = 0; ny = 0; nz = 0;
			if (imagine_cell) {
				//PeriodCellCloestDist2(dr, cell, R2, u, nx, ny, nz)
				u.v[0] = dr.v[0]; t1 = cell.xd[0] * 0.5;
				while (u.v[0] > t1) {
					if (u.v[0] > 0) {u.v[0] -= cell.xd[0]; nx -= 1;}
					else {u.v[0] += cell.xd[0]; nx += 1;}
				};
				u.v[1] = dr.v[1]; t1 = cell.xd[1] * 0.5;
				while (u.v[1] > t1) {
					if (u.v[1] > 0) {u.v[1] -= cell.xd[1]; ny -= 1;}
					else {u.v[1] += cell.xd[1]; ny += 1;}
				};
				u.v[2] = dr.v[2]; t1 = cell.xd[2] * 0.5;
				while (u.v[2] > t1) {
					if (u.v[2] > 0) {u.v[2] -= cell.xd[2]; nz -= 1;}
					else {u.v[2] += cell.xd[2]; nz += 1;}
				};
				scale_uv3(u, u, R2)
			}
			else {
				//ABS2(dr, i, R2) 
				scale_uv3(dr, dr, R2)
			}

			if (!pc1->eNeutral) {
				if (R2 < d2 * _N_MULTP) rship.attach(pc1); // this cluster is too close to the atom, and has electrostatic interaction
				else rship.attach(pc1->emp); // use the multipole method on this cluster
			}
			// this cluster is in LJ interaction region ?
			rcut_LJ = ::rcut_LJ + pc1->d * 0.5 + radius; r2cut_LJ = rcut_LJ * rcut_LJ;
			if (R2 < r2cut_LJ && R2 > _R2min_Interact) LJrship.attach(pc1, nx, ny, nz); // this cluster in [i, j, k] image cell is in LJ interaction region

			pch1 = pch1->next;
		}
		//}
		/*
		else {
			for (i = 0; i < pcell->subCell->n; i++) for (j = 0; j < pcell->subCell->n; j++) for (k = 0; k < pcell->subCell->n; k++) {
				pcl = &(pcell->subCell->cell[i][j][k]);
				if (pcl->ch != NULL) CheckClusterRelship(cell, r0, radius, pcl, pcell->subCell->hw, rship, LJrship);
			}
		}
		*/
	}
	else { // cell is far away from r0
		if (pcell->ch != NULL && !pcell->eNeutral) rship.attach(pcell->emp);
		// cluster inside the cell can not has LJ interaction with r0, because it is too far away
	}
}

// check the relationship of cluster at r0 
void CheckClusterRelship(CMM_CELL3D<BASIC_CLUSTER> &cell, VECTOR3 &r0, float radius, CMM_RELSHIP &rship, LJ_RELSHIP &LJrship) {
	rship.release_relationship();
	LJrship.release_relationship();
	int nc;
	CMM_CELL<BASIC_CLUSTER> *pcl = NULL;
	for (nc = 0; nc < cell.nCells; nc++) {
		pcl = cell.cell + nc;
		if (pcl->ch != NULL) CheckClusterRelship(cell, r0, radius, pcl, cell.hw, rship, LJrship);
	}
}

#endif //_EWALD_SUM == _MultiPole_EWALD_SUM

namespace _cmm_3d_ {

/************************************************************************************
                        LJ relationship only
************************************************************************************/
void CheckCluster_LJRelship(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, LJ_RELSHIP &LJrship, SPME_RELSHIP &spmeRship, int nc1, int nc2, int nc3) {
	int i, j, k, nx, ny, nz, d1cmm, d2cmm, d3cmm;
	VECTOR3 dr, r0, r1;
	float d2, rcut = 0, r2cut = 0;
	float rmax = (::rcut_LJ > ::rcut_Ewd ? ::rcut_LJ : ::rcut_Ewd) + cell.fw_max; //::min_cmm_size;

	LJrship.release_relationship(); spmeRship.release_relationship();
	bool bAttach_LJ = true, bAttach_spme = true;

	V32V3((*(pc->vp)), r0)

	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc2 = NULL;
	int dnx, dny, dnz;
	//float xd[3] = {cell.xd[0] / cell.bcell.nx, cell.xd[1] / cell.bcell.ny, cell.xd[2] / cell.bcell.nz};
	float xd[3] = {cell.hw[0] * 2, cell.hw[1] * 2, cell.hw[2] * 2};

	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > iterator;

	if (::local_relax) {dnx = 1; dny = 1; dnz = 1;}
	else {
		dnx = int(rmax / xd[0] + 0.5) + 1; 
		dny = int(rmax / xd[1] + 0.5) + 1; 
		dnz = int(rmax / xd[2] + 0.5) + 1; 
	}

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
				pcell = &(cell.bcell.m[nx].m[ny].m[nz]); iterator.set(pcell); iterator.start_iterate();
				while (!iterator.eoi()) {
					pc2 = iterator.current();
					if (pc2 == NULL) {iterator.next_iterate(); continue;}
					if (pc2 == pc && d1cmm == 0 && d2cmm == 0 && d3cmm == 0) {
						if (bAttach_spme) bAttach_spme = spmeRship.attach(pc2, d1cmm, d2cmm, d3cmm); // it has to be considered for the electrostatic interaction
						iterator.next_iterate();
						continue;
					}
					if (cell.periodic_cell) {
						r1.v[0] = pc2->vp->v[0] + d1cmm * cell.xd[0];
						r1.v[1] = pc2->vp->v[1] + d2cmm * cell.xd[1];
						r1.v[2] = pc2->vp->v[2] + d3cmm * cell.xd[2];
					}
					else {V32V3((*(pc2->vp)), r1)}
					VECT3(r0, r1, dr) scale_uv3(dr, dr, d2) // dr = r1 - r0, d2 = |dr|^2
					rcut = ::rcut_LJ + (pc->d + pc2->d); r2cut = rcut * rcut;
					if (d2 <= r2cut) {
						if (bAttach_LJ) bAttach_LJ = LJrship.attach(pc2, d1cmm, d2cmm, d3cmm);
					}
				#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
					if (pc->eNeutral || pc2->eNeutral) {}
					else {
						if (!local_interact) {
							rcut = ::rcut_Ewd + (pc->d + pc2->d); r2cut = rcut * rcut;
							if (d2 <= r2cut) {
								if (bAttach_spme) bAttach_spme = spmeRship.attach(pc2, d1cmm, d2cmm, d3cmm);
							}
						}
					}
				#endif
					iterator.next_iterate();
				}
	}}}
}

// check the relationship of a cluster in imagine cells, non-corase-grain
void CheckClusterRelship(BASIC_CLUSTER *pc, CMM_CELL3D<BASIC_CLUSTER> &cell) {
	VECTOR3 r0;
	int n1, n2, n3;
	//if (pc->eNeutral) {
		V32V3((*(pc->vp)), r0)
		n1 = int((r0.v[0] - cell.xl[0]) / (cell.xd[0] / cell.bcell.nx) + 0.5f);
		n2 = int((r0.v[1] - cell.xl[1]) / (cell.xd[1] / cell.bcell.ny) + 0.5f);
		n3 = int((r0.v[2] - cell.xl[2]) / (cell.xd[2] / cell.bcell.nz) + 0.5f);
		if (n1 < 0) n1 = 0;
		else if (n1 >= cell.bcell.nx) n1 = cell.bcell.nx - 1;
		if (n2 < 0) n2 = 0;
		else if (n2 >= cell.bcell.ny) n2 = cell.bcell.ny - 1;
		if (n3 < 0) n3 = 0;
		else if (n3 >= cell.bcell.nz) n3 = cell.bcell.nz - 1;
		CheckCluster_LJRelship(cell, pc, pc->LJRship, pc->spmeRship, n1, n2, n3);
	//}
	//else {
	//	V32V3(pc->emp->r0, r0)
	//	CheckClusterRelship(cell, r0, pc->d * 0.5f, pc->ERship, pc->LJRship);
	//}
}

/*****************************************************************************
               FUNCTIONS TO CHECK RELATIONSHIP
******************************************************************************/
// check the relationship of all clusters in CELL [nCell1 ~ nCell2]
void CheckRelshipInCells(CMM_CELL3D<BASIC_CLUSTER> &cell, int nCell1, int nCell2) {
	int nc;// = number(cell.cell[0].ch);
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	BASIC_CLUSTER *pc = NULL;
	VECTOR3 r0;
	ITERATOR<BASIC_CLUSTER, CMM_CELL<BASIC_CLUSTER> > iterator;
	for (nc = nCell1; nc <= nCell2; nc++) {
		pcell = cell.acell.m[nc]; iterator.set(pcell); iterator.start_iterate();
		while (!iterator.eoi()) {
			pc = iterator.current();
			CheckClusterRelship(pc, cell);
			iterator.next_iterate();
		}
	}
}

// check the relationship of all clusters in WHOLE CELL
void ParallelCheckRelshipInCells(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, int nThreads) {
	return CellOperate<BASIC_CLUSTER>((void*)(&CheckRelshipInCells), cell, n1Cell, n2Cell, nThreads);
}


/***************************************************************************
******      SPECIALLY USING FOR SINGLE MACRO-MOLECULE SIMULATION      ******
****************************************************************************/

void MMClusterRelshipInCell(CMM_CELL3D<BASIC_CLUSTER> *cell, MMOLECULE *mm, int nc1, int nc2) {
	BASIC_CLUSTER *pc = NULL;
	for (int nc = nc1; nc <= nc2; nc++) {
		if (nc >= mm->nCluster) break;
		pc = (BASIC_CLUSTER*)(mm->cluster + nc);
		CheckClusterRelship(pc, *cell);
	}
}

void SMRelshipInCell(SMOLECULE *sm, int nSM, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	for (int nc = 0; nc < nSM; nc++) CheckClusterRelship(sm[nc].c, *cell);
}

void PMRelshipInCell(PMOLECULE *pm, int nPM, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	for (int nc = 0; nc < nPM; nc++) CheckClusterRelship(pm[nc].c, *cell);
}

// check the relationship of all clusters in WHOLE CELL
void MMCheckRelshipInCell(MMOLECULE &mm, int nc1, int nc2, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads) {
	return MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MMClusterRelshipInCell), cell, mm, nc1, nc2, nThreads);
}

// check the relationship of all clusters in WHOLE CELL
void SMCheckRelshipInCell(SMOLECULE *sm, int nSM, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads) {
	return SMCellOperate<SMOLECULE, BASIC_CLUSTER>((void*)(&SMRelshipInCell), cell, sm, nSM, nThreads);
}

// check the relationship of all clusters in WHOLE CELL
void PMCheckRelshipInCell(PMOLECULE *pm, int nPM, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads) {
	return SMCellOperate<PMOLECULE, BASIC_CLUSTER>((void*)(&PMRelshipInCell), cell, pm, nPM, nThreads);
}


void CMMCheckImpSolvNeighbors_Cell(CMM_CELL3D<BASIC_CLUSTER> &cell, int nc1, int nc2) {
	CheckImpSolvNeighbors_Cell<BASIC_CLUSTER>(cell, nc1, nc2);
}

void MMClusterImpSolvNeighborsInCell(CMM_CELL3D<BASIC_CLUSTER> *cell, MMOLECULE *mm, int nc1, int nc2) {
	return MacroMolClusterImpSolvNeighborsInCell<MMOLECULE, BASIC_CLUSTER>(cell, mm, nc1, nc2);
}

} // end of namespace _cmm_3d_

#undef _N_MULTP

/**************************************************************************************************************************
	following functions are specially designed for the CMM using ARRAY to store the cluster inside each CMM_CELL
**************************************************************************************************************************/

bool MM_AllocateCluster_CMM_array(MMOLECULE *mm, int nmol, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	bool status = true;
	for (int imol = 0; imol < nmol; imol++) {
		if (!MM_AllocateCMM_array<BASIC_CLUSTER, MMOLECULE>(mm + imol, cell)) status = false;
	}
	return status;
};

bool SM_AllocateCluster_CMM_array(SMOLECULE *sm, int nmol, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	char msg[256] = "\0";
	bool status = true;
	int ix = 0, iy = 0, iz = 0;
	VECTOR3 rc, dr;
	BASIC_CLUSTER *pc = NULL;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	for (int imol = 0; imol < nmol; imol++) {
		pc = sm[imol].c;
		V32V3(sm[imol].r, rc)  // SMOLECULE has only one cluster, no dr_cmm has to be used
		CMMCELL_r2INDX(rc, (*cell), ix, iy, iz, dr)
		pcell = &(cell->bcell.m[ix].m[iy].m[iz]);
		if (!pcell->attach(pc)) {
			sprintf(msg, "subcell [%d, %d, %d] leaks", ix, iy, iz); show_log(msg, true); status = false;
		}
	}
	return status;
};

bool PM_AllocateCluster_CMM_array(PMOLECULE *pm, int nmol, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	char msg[256] = "\0";
	bool status = true;
	int ix = 0, iy = 0, iz = 0;
	VECTOR3 rc, dr;
	BASIC_CLUSTER *pc = NULL;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	for (int imol = 0; imol < nmol; imol++) {
		pc = pm[imol].c;
		V32V3(pm[imol].r, rc)  // SMOLECULE has only one cluster, no dr_cmm has to be used
		CMMCELL_r2INDX(rc, (*cell), ix, iy, iz, dr)
		pcell = &(cell->bcell.m[ix].m[iy].m[iz]);
		if (!pcell->attach(pc)) {
			sprintf(msg, "subcell [%d, %d, %d] leaks", ix, iy, iz); show_log(msg, true); status = false;
		}
	}
	return status;
};


bool AllocateCluster_CMM_array(ARRAY<BASIC_CLUSTER*> &bc, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	char msg[256] = "\0";
	bool status = true;
	int ix = 0, iy = 0, iz = 0;
	VECTOR3 rc, dr;
	BASIC_CLUSTER *pc = NULL;
	CMM_CELL<BASIC_CLUSTER> *pcell = NULL;
	for (int ic = 0; ic < bc.n; ic++) {
		pc = bc.m[ic];
		V32V3((*(pc->vp)), rc)
		CMMCELL_r2INDX(rc, (*cell), ix, iy, iz, dr)
		pcell = &(cell->bcell.m[ix].m[iy].m[iz]);
		if (!pcell->attach(pc)) {
			sprintf(msg, "subcell [%d, %d, %d] leaks", ix, iy, iz); show_log(msg, true); status = false;
		}
	}
	return status;
};

void AllocateCluster_CMM(JOB_THREAD_VARS<CMM_AtomJob<BASIC_CLUSTER> > *job) {
	if (job->job->atoms.m != NULL && job->job->cmm != NULL) AllocateCluster_CMM_array(job->job->atoms, job->job->cmm);
}
