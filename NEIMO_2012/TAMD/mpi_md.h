void MPI_GatherClusterForce(MMOLECULE *mm);
void MPI_DistributeMovePars(VECTOR<6> &dv_base, double *dta, int NC);

void MPI_POLYMER(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, CMM_CELL3D<BASIC_CLUSTER> *cmm);
