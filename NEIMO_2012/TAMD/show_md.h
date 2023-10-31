bool get_single_mm_md_proc(char *mdpar_fname, int isave);
bool get_single_cgmm_md_proc(char *mdpar_fname, int isave);
bool show_single_mm_md_proc(char *mdpar_fname, int nf);
bool show_single_cgmm_md_proc(char *mdpar_fname, int nf);

bool get_multi_mm_md_proc(char *mdpar_fname, int isave);
bool get_multi_cgmm_md_proc(char *mdpar_fname, int isave);
bool show_multi_mm_md_proc(char *mdpar_fname, int nf);
bool show_multi_cgmm_md_proc(char *mdpar_fname, int nf);

void ClusterSaveXYZ(ofstream *out, BASIC_CLUSTER *pc, bool bInCMM = true, double *dv = NULL);

void MMSaveXYZ(ofstream *out, MMOLECULE *mm, ARRAY<short> &id_monomer, ARRAY<short> &id_cluster, bool bInCMM = true, double *dv = NULL);
void SMSaveXYZ(ofstream *out, SMOLECULE *sm, ARRAY<short> &id_cluster);
void PMSaveXYZ(ofstream *out, PMOLECULE *pm, ARRAY<short> &id_cluster);

namespace _coarse_grain_ {

void CGClusterSaveXYZ(ofstream *out, CG_CLUSTER *pc, bool bInCMM = true, double *dv = NULL);
void CGMMSaveXYZ(ofstream *out, CG_MMOLECULE *mm, ARRAY<short> &id_cgatom, bool bInCMM, double *dv = NULL);
void CGSMSaveXYZ(ofstream *out, CGSM *sm, ARRAY<short> &id_cgatom);

} // end of namespace _coarse_grain
