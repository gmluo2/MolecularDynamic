#ifdef PeriodCellCloestDist2
#undef PeriodCellCloestDist2
#endif
#define PeriodCellCloestDist2_2d(dr, cell, R2, u, nx, ny, nz) u.v[0] = FABS(dr.v[0]); u.v[1] = FABS(dr.v[1]); u.v[2] = FABS(dr.v[2]); \
	if (u.v[0] > cell.xd[0] * 0.5) {u.v[0] = cell.xd[0] - u.v[0]; nx = (dr.v[0] > 0 ? -1 : 1);} else nx = 0;\
	if (u.v[1] > cell.xd[1] * 0.5) {u.v[1] = cell.xd[1] - u.v[1]; ny = (dr.v[1] > 0 ? -1 : 1);} else ny = 0;\
	nz = 0;\
	R2 = u.v[0] * u.v[0] + u.v[1] * u.v[1] + u.v[2] * u.v[2];

namespace _cmm_2d_ {
/************************************************************************************
                        implicit solvation for clusters in x-y 2d CMM
************************************************************************************/
// pc0 is the center cluster, nc0[3] is the cell's (where pc0 is) index in basic_cell_matrix
// dnc[3] is the offset of cell index, where the cluster inside could be neighbors
// xd[3] is the dimension of basic_cell
template <class atom> void CheckImpSolvNeighbors_Cluster(CMM_CELL3D<atom> &cell, atom *pc0, int *nc0, int *dnc, float r0_solvent, IMAG_RELSHIP<atom> &impsolv_rship) {
	bool imagine_cell = cell.periodic_cell;
	int i, j, k, ni, nj, nk, nx, ny, nz;
	SVECTOR<double, 3> r0, dr, r1;
	if (pc0->bCMMShift) {V3plusV3(pc0->rc_solv, pc0->dr_cmm, r0)} // r0 is pc0's rc
	else {V32V3(pc0->rc_solv, r0)}

	float d2, rcut = 0, r2cut;
	CMM_CELL<atom> *pcell = NULL;
	atom *pc = NULL;
	int n_relationships = 0;
	ITERATOR<atom, CMM_CELL<atom> > iterator;

	for (i = nc0[0] - dnc[0]; i <= nc0[0] + dnc[0]; i++) {
		ni = i; nx = 0;
		if (ni < 0) {
			if (imagine_cell) {while (ni < 0) {ni += cell.nx; nx -= 1;}}
			else continue;
		}
		else if (ni >= cell.bcell.nx) {
			if (imagine_cell) {while (ni >= cell.bcell.nx) {ni -= cell.bcell.nx; nx += 1;}}
			else continue;
		}
		for (j = nc0[1] - dnc[1]; j <= nc0[1] + dnc[1]; j++) {
			nj = j; ny = 0;
			if (nj < 0) {
				if (imagine_cell) {while (nj < 0) {nj += cell.bcell.ny; ny -= 1;}}
				else continue;
			}
			else if (nj >= cell.bcell.ny) {
				if (imagine_cell) {while (nj >= cell.bcell.ny) {nj -= cell.bcell.ny; ny += 1;}}
				else continue;
			}
			for (k = nc0[2] - dnc[2]; k <= nc0[2] + dnc[2]; k++) {
				nk = k; nz = 0; // z is not periodical
				if (nk < 0) continue;
				else if (nk >= cell.bcell.nz) continue;

				pcell = &(cell.bcell.m[ni].m[nj].m[nk]);
				// this the deepest level
				iterator.set(pcell); iterator.start_iterate();
				while (!iterator.eoi()) {
					pc = iterator.current(); if (pc == NULL) {iterator.next_iterate(); continue;}
					if (pc == pc0) {
						if (!imagine_cell) {
							iterator.next_iterate();
							continue;
						}
						else if (nx == 0 && ny == 0 && nz == 0) {
							iterator.next_iterate();
							continue;
						}
					}
					if (imagine_cell) {
						if (pc->bCMMShift) {V3plusV3(pc->rc_solv, pc->dr_cmm, r1)}
						else {V32V3(pc->rc_solv, r1)}
						if (nx != 0) r1.v[0] += nx * cell.xd[0];
						if (ny != 0) r1.v[1] += ny * cell.xd[1];
					}
					else {V32V3(pc->rc_solv, r1)}
					VECT3(r0, r1, dr) // dr = r1 - r0
					V3ABS2(dr, d2)
					rcut = pc0->psolv->r0 + pc->psolv->r0 + r0_solvent + r0_solvent;
					r2cut = rcut * rcut;
					if (d2 <= r2cut) {
						if (::local_relax) {
							if (d2 + d2 < r2cut) {
								impsolv_rship.attach(pc, nx, ny, nz);
								n_relationships++;
								if (n_relationships > 20) return; // local_relax case: do not consider too much
							}
						}
						else impsolv_rship.attach(pc, nx, ny, nz);
					}
					iterator.next_iterate();
				}
			}
		}
	}
};

template <class atom> void ClusterImpSolvNeighborsInCell(CMM_CELL3D<atom> *cell, atom *pc, int *dnc) {
	// nc is the subcell index of the cluster
	// dnc is the neighbor subcell to search around
	int nc[3];
	if (pc->psolv == NULL) return;
	SVECTOR<double, 3> r0, dr;
	// position of cluster
	V32V3((*(pc->vp)), r0)
	// cell index where cluster in the basic_cell_matrix
	CMMCELL_r2INDX(r0, (*cell), nc[0], nc[1], nc[2], dr)

	pc->impsolv_rship.release_relationship();
	CheckImpSolvNeighbors_Cluster<atom>(*cell, pc, nc, dnc, ::rSolvent, (*((IMAG_RELSHIP<atom>*)&(pc->impsolv_rship))));
};

template <class atom> void CheckImpSolvNeighbors_Cell(CMM_CELL3D<atom> &cell, int nc1, int nc2) {
	float r0_solvent = rSolvent;
	float xd[3] = {cell.xd[0] / cell.bcell.nx, cell.xd[1] / cell.bcell.ny, cell.xd[2] / cell.bcell.nz};
	int nc[3], dnc[3];
	float rcut = dSolv_max;

	if (::local_relax) {dnc[0] = 1; dnc[1] = 1; dnc[2] = 1;}
	else {
		dnc[0] = int(rcut / xd[0] + 0.5f); if (dnc[0] == 0) dnc[0] = 1; if (dnc[0] * xd[0] < rcut) dnc[0] += 1;
		dnc[1] = int(rcut / xd[1] + 0.5f); if (dnc[1] == 0) dnc[1] = 1; if (dnc[1] * xd[1] < rcut) dnc[1] += 1;
		dnc[0] += 1; dnc[1] += 1;
		dnc[2] = 0;
	}

	CMM_CELL<atom> *pcell = NULL;
	atom *pc = NULL;
	int i, j, k;
	ITERATOR<atom, CMM_CELL<atom> > iterator;
	for (i = nc1; i <= nc2; i++) {
		if (i >= cell.bcell.nx) break;
		for (j = 0; j < cell.bcell.ny; j++) {for (k = 0; k < cell.bcell.nz; k++) {
			pcell = &(cell.bcell.m[i].m[j].m[k]); iterator.set(pcell);
			nc[0] = i; nc[1] = j; nc[2] = k;
			iterator.start_iterate();
			while (!iterator.eoi()) {
				pc = iterator.current();
				if (pc == NULL) {iterator.next_iterate(); continue;}
				if (pc->psolv == NULL) continue;
				pc->impsolv_rship.release_relationship();
				CheckImpSolvNeighbors_Cluster<atom>(cell, pc, nc, dnc, r0_solvent, pc->impsolv_rship);
				iterator.next_iterate();
			}
		}}
	}
	return;
};

// check the implicit-solvation neighbors of each cluster in cell, in parallel

void CMMCheckImpSolvNeighbors_Cell(CMM_CELL3D<BASIC_CLUSTER> &cell, int nc1, int nc2);


template <class MACROMOL, class atom> void MacroMolClusterImpSolvNeighborsInCell(CMM_CELL3D<atom> *cell, MACROMOL *mm, int nc1, int nc2) {
	atom *pc = NULL;
	//float r0_solvent = rSolvent;
	float xd[3] = {cell->xd[0] / cell->bcell.nx, cell->xd[1] / cell->bcell.ny, cell->xd[2] / cell->bcell.nz};
	int dnc[3];
	float rcut = dSolv_max;

	SVECTOR<double, 3> r0;
	
	if (::local_relax) {dnc[0] = 1; dnc[1] = 1; dnc[2] = 1;}
	else {
		dnc[0] = int(rcut / xd[0] + 0.5f); if (dnc[0] == 0) dnc[0] = 1; if (dnc[0] * xd[0] < rcut) dnc[0] += 1;
		dnc[1] = int(rcut / xd[1] + 0.5f); if (dnc[1] == 0) dnc[1] = 1; if (dnc[1] * xd[1] < rcut) dnc[1] += 1;
		dnc[0] += 1; dnc[1] += 1;
		dnc[2] = 0;
	}

	for (int mnc = nc1; mnc <= nc2; mnc++) {
		if (mnc >= mm->nCluster) break;
		pc = (atom*)(mm->cluster + mnc);
		if (pc->psolv == NULL) continue;
		// position of cluster
		//V32V3((*(pc->vp)), r0)
		
		ClusterImpSolvNeighborsInCell<atom>(cell, pc, dnc);
	}
};

void MMClusterImpSolvNeighborsInCell(CMM_CELL3D<BASIC_CLUSTER> *cell, MMOLECULE *mm, int nc1, int nc2);



void CheckClusterRelship(BASIC_CLUSTER *pc, CMM_CELL3D<BASIC_CLUSTER> &cell);
void SingleMMCheckRelshipInCell(MMOLECULE &mm, int nc1, int nc2, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads);
void MMCheckRelshipInCell(MMOLECULE *mm, int nM, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads);
void SMCheckRelshipInCell(SMOLECULE *sm, int nSM, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads);
void PMCheckRelshipInCell(PMOLECULE *pm, int nPM, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads);

void ImageClusterCheckRelshipInCell(ImageCluster *ic, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads);
}; // end of namespace _cmm_2d_ {
