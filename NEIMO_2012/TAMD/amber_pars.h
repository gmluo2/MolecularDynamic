class TORSION_PARS {
public:
	short aindx[4]; // atom indexes of the 4 atoms relating to the dihedral angle
	float v0, phase, pn;
	TORSION_PARS *next;
	TORSION_PARS() {
		v0 = 0; phase = 0; pn = 0; next = NULL;
	};
	void release_tail() {
		if (next == NULL) return;
		delete next; next = NULL;
	};
	~TORSION_PARS() {
		if (next != NULL) {delete next; next = NULL;}
	};
};

bool find_bound_pars(char *fname, char *atom1, char *atom2, float &blen, char *errmsg);
bool find_bound_angle_pars(char *fname, char *atom1, char *atom2, char *atom3, float &angle, char *errmsg);
bool general_find_torsion_pars(char *fname, bool improper, char *atom1, char *atom2, char *atom3, char *atom4, TORSION_PARS &tpar, char *errmsg);

extern char struct_fname[256], force_field_fname[256], atom_index_fname[256];
