namespace _coarse_grain_ {
void get_Hoover(CG_MMOLECULE &mm, LFVERLET_CGMM_MD& mdpar, double Ethermal);
void HooverCorrectionOfAccerlation(CG_MMOLECULE *mm, int nc1, int nc2, HOOVER_DYN& hDyn);
bool CGLeapFrog_Verlet(CG_MMOLECULE &mm, LFVERLET_CGMM_MD& mdpar, VECTOR3 *vbf, double Ethermal, bool init_md, char *show_msg, char *log_msg, int n_msglen);

bool CGSMLeapFrog_Verlet(CGSM &sm, LFVERLET_CGSM_MD& mdpar, double Ethermal, bool init_md, char *show_msg, char *log_msg, int n_msglen);

void ZeroVelocity(CG_MMOLECULE &mm);
// Initilize the clusters' velocity
void InitVelocity(CG_MMOLECULE &mm, double Ek);
void InitCGSMVelocity(CGSM *sm, int nSM, double &Ek);

void MD_POLYMER(CG_MMOLECULE &cgmm, LFVERLET_CGMM_MD& mdpar, CMM_CELL3D<CG_CLUSTER> *cmm, int interact_cal_method);

void show_molecule_struct(CG_MMOLECULE &mm, char *fname);

} // end of namespace _coarse_grain_ 
