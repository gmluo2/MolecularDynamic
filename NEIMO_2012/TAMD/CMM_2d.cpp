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
#include "CMM_2d.h"

#include "var.h"

#if _SYS_ == _WINDOWS_SYS_
extern HANDLE hWaitEvent;
#endif


namespace _cmm_2d_ {
	/*
void CMMCheckImpSolvNeighbors_Cell(CMM_CELL3D<BASIC_CLUSTER> &cell, int nc1, int nc2) {
	CheckImpSolvNeighbors_Cell<BASIC_CLUSTER>(cell, nc1, nc2);
};
*/


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
				if (nz < 0 || nz >= cell.bcell.nz) continue; // z direction is not periodical

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
						r1.v[2] = pc2->vp->v[2];// + pc2->dr_cmm.v[2] + d3cmm * cell.xd[2];
					}
					else {V32V3((*(pc2->vp)), r1)}
					VECT3(r0, r1, dr) scale_uv3(dr, dr, d2) // dr = r1 - r0, d2 = |dr|^2
					rcut = ::rcut_LJ + (pc->d + pc2->d); r2cut = rcut * rcut;
					//if (pc2->cID < 0) {
					//	pc2 = pc2;
					//}
					if (d2 <= r2cut && pc2->cID >= 0) {
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


void MMRelshipInCell(MMOLECULE *mm, int nMM, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	BASIC_CLUSTER *pc = NULL;
	int nm, nc;
	for (nm = 0; nm < nMM; nm++) {
		for (nc = 0; nc < mm[nm].nCluster; nc++) {
			pc = (BASIC_CLUSTER*)(mm[nm].cluster + nc);
			CheckClusterRelship(pc, *cell);
		}
	}
}

void SMRelshipInCell(SMOLECULE *sm, int nSM, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	for (int nc = 0; nc < nSM; nc++) CheckClusterRelship(sm[nc].c, *cell);
}

void PMRelshipInCell(PMOLECULE *pm, int nPM, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	for (int nc = 0; nc < nPM; nc++) CheckClusterRelship(pm[nc].c, *cell);
}


// check the relationship of all clusters in WHOLE CELL
void SingleMMCheckRelshipInCell(MMOLECULE &mm, int nc1, int nc2, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads) {
	return MacroMolCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MMClusterRelshipInCell), cell, mm, nc1, nc2, nThreads);
}

void MMCheckRelshipInCell(MMOLECULE *mm, int nM, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads) {
	return SMCellOperate<MMOLECULE, BASIC_CLUSTER>((void*)(&MMRelshipInCell), cell, mm, nM, nThreads);
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







// modified on 02/06/2012
// for image clusters, checking only the local relation with real cluster for electrostatic interaction
// image cluster has cID < 0, real cluster has cID >= 0
void CheckImageCluster_Relship(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, SPME_RELSHIP &spmeRship, int nc1, int nc2, int nc3) {
	int i, j, k, nx, ny, nz, d1cmm, d2cmm, d3cmm;
	VECTOR3 dr, r0, r1;
	float d2, rcut = 0, r2cut = 0;
	float rmax = ::rcut_Ewd + cell.fw_max; //::min_cmm_size;

	spmeRship.release_relationship();
	bool bAttach_spme = true;

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
				if (nz < 0 || nz >= cell.bcell.nz) continue; // z direction is not periodical

				// check the relationship within the cell
				pcell = &(cell.bcell.m[nx].m[ny].m[nz]); iterator.set(pcell); iterator.start_iterate();
				while (!iterator.eoi()) {
					pc2 = iterator.current();
					if (pc2 == NULL || pc2->cID < 0) {iterator.next_iterate(); continue;} // ignore image cluster also, since only the relationship with real cluster is considered
					if (pc2 == pc && d1cmm == 0 && d2cmm == 0 && d3cmm == 0) {
						if (bAttach_spme) bAttach_spme = spmeRship.attach(pc2, d1cmm, d2cmm, d3cmm); // it has to be considered for the electrostatic interaction
						iterator.next_iterate();
						continue;
					}
					if (cell.periodic_cell) {
						r1.v[0] = pc2->vp->v[0] + d1cmm * cell.xd[0];
						r1.v[1] = pc2->vp->v[1] + d2cmm * cell.xd[1];
						r1.v[2] = pc2->vp->v[2];// + pc2->dr_cmm.v[2] + d3cmm * cell.xd[2];
					}
					else {V32V3((*(pc2->vp)), r1)}
					VECT3(r0, r1, dr) scale_uv3(dr, dr, d2) // dr = r1 - r0, d2 = |dr|^2
					//rcut = ::rcut_LJ + (pc->d + pc2->d); r2cut = rcut * rcut;
					//if (pc2->cID < 0) {
					//	pc2 = pc2;
					//}
					//if (d2 <= r2cut && pc2->cID >= 0) {
					//	if (bAttach_LJ) bAttach_LJ = LJrship.attach(pc2, d1cmm, d2cmm, d3cmm);
					//}
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
void CheckImageClusterRelship(BASIC_CLUSTER *pc, CMM_CELL3D<BASIC_CLUSTER> &cell) {
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
		CheckImageCluster_Relship(cell, pc, pc->spmeRship, n1, n2, n3);
	//}
	//else {
	//	V32V3(pc->emp->r0, r0)
	//	CheckClusterRelship(cell, r0, pc->d * 0.5f, pc->ERship, pc->LJRship);
	//}
}

void ImageClusterRelshipInCell(ImageCluster *ic, int n, CMM_CELL3D<BASIC_CLUSTER> *cell) {
	for (int nc = 0; nc < n; nc++) CheckImageClusterRelship((BASIC_CLUSTER*)(ic + nc), *cell);
}

// check the relationship of all clusters in WHOLE CELL
void ImageClusterCheckRelshipInCell(ImageCluster *ic, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads) {
	return SMCellOperate<ImageCluster, BASIC_CLUSTER>((void*)(&ImageClusterRelshipInCell), cell, ic, n, nThreads);
}

} // end of namespace _cmm_2d_ {
