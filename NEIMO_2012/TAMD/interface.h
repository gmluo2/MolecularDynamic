namespace _interface_constraint_ {

class INTERFACE_CONSTRAINT {
public:
	// 2d-interface is defined within a region [z1, z2], with given energy profile and force with unit in kT and kT/Angs
	float z0, zf, zt;
	EF_DPROF ef_prof;

	void construct_efprofile_erf(float U0, float z0, float sigma);
	void ef(float z, double &U, double &f);
};

class IConstraint : public COMM_PARS<50>, public INTERFACE_CONSTRAINT {
public:
	char title[50];
	int I_indx; // index of the interface
	IConstraint() {I_indx = 0; z0 = 0; zf = 0; zt = 0; strcpy(title, "\0");};
};

class Barrier {
public:
	float zf, zt, z0;
};


#ifdef _INTERFACE_INIT_
	FLEX_ASYM_MATRIX<IConstraint> clusterIC_db;
	Barrier barrier[2];
#else
	extern FLEX_ASYM_MATRIX<IConstraint> clusterIC_db;
	extern Barrier barrier[2];
#endif

bool init_clusterIC_DB(int nInterface, float z1, float z2);

bool read_interface_constraints(IConstraint *ic);

bool read_interface_image_constant(char *par_fname, float &f1_img, float &z1_img, float &f2_img, float &z2_img);

void apply_slab_barrier_constraint(FLEX_ASYM_MATRIX<IConstraint> &clusterIC_db);
} // end of namespace _pmf_constraint_
