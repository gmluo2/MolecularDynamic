
extern "C" void set_md_mode(int mode);
extern "C" void set_job(int job);
// get the molecule in current job
//bool get_molecule_upto_job;
extern "C" int get_total_atom_number();
extern "C" bool construct_macromolecule(char *fname);

void MD_SIM_SINGLE_MOL(char *mdpar_fname);
void MD_SIM_MULTI_MOL_CELL(char *mdpar_fname);

extern "C" void MD_SIM(char *mdpar_fname) ;

extern "C" bool running_md(void);
extern "C" void stop_md(bool &status);

extern "C" int get_molecule_struct(float *x, float *y, float *z, bool *highlight, int MAX, int nMol, int nstart);
extern "C" int get_molecule(int *atom, float *x, float *y, float *z, bool *highlight, int MAX, int nMol, int nstart);

extern "C" int get_atom_bound(int *nb, char *cvb_type, int MAX, int nMol, int natom);

extern "C" void get_errmsg(char *errmsg) ;

extern "C" void get_mol_struct_over();

extern "C" void set_show_infor(int min_time, int max_time, int nshow);

extern "C" void get_show_infor(int* min_time, int *max_time, int *nshow);

extern "C" bool get_single_md_proc(char *mdpar_fname, int isave);
extern "C" bool show_single_md_proc(char *mdpar_fname, int nf);

void clear_memory();

extern "C" bool init_env();
extern "C" void clear_env();
