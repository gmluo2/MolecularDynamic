void rewind(ifstream &in);
#ifndef file_pos
#define file_pos streampos
#endif
file_pos current_fpos(ifstream &in);
void set_fpos(ifstream &in, file_pos &pos);

void Upper(char* str);
bool NextSection(char** original, char ch = ' ');  //default ch = ' '
bool pass(char** original, char ch = ' '); // go to next non-ch character
void kick_off_tail(char *str, char ch = ' ');
int read_float(char *buffer, float **ta);
/*
void Upper(char* str);
bool NextSection(char** original, char ch = ' ');  
*/
bool NextSection(char** original, char* ch);
int nItem(char* str, char ch = ' ');
int nItem(char* str, char* ch);
//bool pass(char** original, char ch = ' '); // go to next non-ch character
bool pass(char** original, char* ch);
//void kick_off_tail(char *str, char ch);
void kick_off_tail(char *str, char* ch);
//bool search(ifstream &in, char *title, char *buffer);
//bool search(ifstream &in, char *title, int indx, char *buffer);
int read_float(char *buffer, float *v, int max_v, char* format, char seperator = ' ');
int read_float(char *buffer, float *v, int max_v, char* format, char* seperator);

void pickup_str(char *source, int nstart, int len, char *dest);
bool next_section(char** original, char* ch);
bool pass_next_section(char** original, char* ch);

int read_float(char *buffer, int ncol, float *v, char *format, char* seperator);
bool read_profile(char *fname, int nIgnoreLines, float **x, int ix, float **y, int iy, int& npts);
bool read_profile(char *fname, int nIgnoreLines, ARRAY<float*> &va, ARRAY<int> &icol, int& npts);

bool read_msg(char *fname, char *msg, char *item, bool full_line = false);
extern bool read_item(char *fname, char *title, char *item);
extern bool search(ifstream &in, char *title, char *buffer);

#if _READ_MM_MD_ == 1
// assuming the macromolecule structure is given
// and reading the torsion angle to reconstruct the new configuration
bool read_mol_MD_TA_struct(ifstream &in, MMOLECULE &mol, int njump, char *errmsg);

bool read_MD_SAVE_Pars(char *fname, MD_SAVE &msave);
bool read_MD_BASE_Pars(char *fname, MD_BASE *mdpar);
bool read_MD_TIMESTEP(char *fname, float &dt);
bool read_HOOVER_Pars(char *fname, float &tao, float &max_ksi);
bool read_MD_LOOPS(char *fname, long &nloops);

bool readExternalElectricField(char *fname, bool &bEext, float &Eext);

bool readForceScale(char *fname, bool &scale_force, float &force_max, float &torque_max);

bool readDielectricConst(char *fname, float &eps);

bool readLocalRelax(char *fname, bool &local_relax, float &dta_max_relax);
bool readLocalInteract(char *fname, bool &local_interact);
bool readInternalKinEnergy(char *fname, float &Ek0_Internal);
bool readImplicitSolvent(char *fname, bool &bImplicitSolvent);

bool readMacroMolFreeBase(char *fname, bool &free_base);

bool readPolarization(char *fname, bool &bPolarize, int &iPolarize);
bool readPolarization1(char *fname, bool &bExcludeInducedDipole);

bool readCoarseGrainMD(char *fname, bool &coarse_grain, char *cgpar_fname, int &nNormalMD, int &nCG);
bool read_coarse_grain_md_var(char *fname, float &dt, float &tao, float &max_ksi, int &nCMMRefresh);

bool readBrownianPars(char *fname, bool &Brownian, float &ita_Bronian, float &cd_trans, float &cd_rot);

bool readSPMEPars(char *fname, int &indx_bspline, int &dim_spme_x, int &dim_spme_y, int &dim_spme_z, bool &bExpandSPMECellZ);

bool readNPTPars(char *fname, float& W_npt, float &tao_npt);

bool readNPT_RF_Exclude(char *fname, bool& bNPT_RF_Exclude, int &nmax, float &Fmin);

bool readSurfaceBoundary(char *fname);
bool readNetForceCorrectionMode(char *fname);

#endif

