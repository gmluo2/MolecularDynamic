void init_LJ_Pars();
// maximum distance for L-J interaction cut-off
float rLJ_cut_max();
void init_LJ_Pars_CG();

void init_ImplicitSolvation_Pars();

void MD(MMOLECULE &mm, char *title, float dt, int nsave, long nloops);

void MD(MMOLECULE &mm, LFVERLET_MM_MD &mdpar, CMM_CELL3D<BASIC_CLUSTER> *cmm, int interact_method, int nNormalMD);

