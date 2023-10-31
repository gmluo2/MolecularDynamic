namespace _coarse_grain_ {

void cmm_init_distribute_cluster(CG_MMOLECULE &mm, CMM_CELL3D<CG_CLUSTER> &cell, bool reset_cell_chain);
void cmm_init_distribute_cluster(CGSM *sm, int nSM, CMM_CELL3D<CG_CLUSTER> &cell, bool reset_cell_chain);

// imethod == 0 : LJ interaction for all clusters, ignoring the neighbors
// imethod == 1 : LJ interaction for the clusters in different neighbors
void CGMMCheckLocalRelshipInCell(CMM_CELL3D<CG_CLUSTER> &cell, CG_MMOLECULE &mm, int nc1, int nc2, int imethod, int nThreads);
void CGSMCheckLocalRelshipInCell(CMM_CELL3D<CG_CLUSTER> &cell, CGSM *sm, int nSM, int nThreads);
void CheckLocalRelship(CMM_CELL3D<CG_CLUSTER> &cell, int n1Cell, int n2Cell, int imethod_rcut, int nThreads);

void CGMMClusterImpSolvNeighborsInCell(CMM_CELL3D<CG_CLUSTER> *cell, CG_MMOLECULE *cgmm, int nc1, int nc2);

} // end of namespace _coarse_grain_ 
