
class MDControlVars {
public:
	bool bPolarize, bDipole, bRealEwaldSumOnly, bLocalInteract;
	bool bSaveVirial, bSavePz, bSave_mu;
	bool bScaleForce;

	MDControlVars() {bPolarize = false; bDipole = false; bRealEwaldSumOnly = false; bLocalInteract = false; bScaleForce = false;}
	void show_status() {
		char msg[2560] = "\0", tmsg[256] = "\0";
		strcat(msg, "ATOMIC_DIPOLE: "); if (bDipole) strcat(msg, "Enabled; "); else strcat(msg, "Disabled; ");
		strcat(msg, "ATOMIC_POLARIZATION: "); if (bPolarize) strcat(msg, "Enabled; "); else strcat(msg, "Disabled; ");
		strcat(msg, "EwaldSum_RealPartOnly: "); if (bRealEwaldSumOnly) strcat(msg, "true; "); else strcat(msg, "false; ");
		strcat(msg, "LOCAL_INTERACTION: "); if (bLocalInteract) strcat(msg, "Enabled; "); else strcat(msg, "Disabled; ");
		strcat(msg, "SCALE_FORCE: "); if (bScaleForce) strcat(msg, "Enabled; "); else strcat(msg, "Disabled; ");
		strcat(msg, "SAVE_VIRIAL: "); if (bSaveVirial) strcat(msg, "Enabled; "); else strcat(msg, "Disabled; ");
		strcat(msg, "SAVE_PZ: "); if (bSavePz) strcat(msg, "Enabled; "); else strcat(msg, "Disabled; ");
		strcat(msg, "SAVE_MU: "); if (bSave_mu) strcat(msg, "Enabled; "); else strcat(msg, "Disabled; ");
		show_log(msg, true);
	};
	void backup_global_vars() {
		this->bPolarize = ::bPolarize; this->bDipole = ::bDipole; this->bRealEwaldSumOnly = ::bRealEwaldSumOnly;
		this->bLocalInteract = ::local_interact; this->bScaleForce = ::scale_force;
		this->bSaveVirial = ::bVirial;
	};
	void update_global_vars() {
		::bPolarize = this->bPolarize; ::bDipole = this->bDipole; ::bRealEwaldSumOnly = this->bRealEwaldSumOnly;
		::local_interact = this->bLocalInteract; ::scale_force = this->bScaleForce;
		::bVirial = this->bSaveVirial;
	};
};


class MD_PRO {
public:
	MDControlVars cvar;
	char parfname_IConstraint[256];
	float eThermal; // thermal energy of each freedom
	long nloops;

	MD_PRO() {
		strcpy(parfname_IConstraint, "\0");
		eThermal = 0.5f; // kT
	};
	void show_status() {
		char msg[256] = "\0";
		cvar.show_status();
		sprintf(msg, "thermal energy of each freedom : %5.2f kT", eThermal); show_log(msg, true);
		if (strlen(parfname_IConstraint) != 0) {show_log("Interface contraint from file : ", false); show_log(parfname_IConstraint, true);}
		else show_log("using same interface constraint", true);
	};
};

bool read_mdprocs(char *fname, ARRAY<MD_PRO> &mdpro);


class MC_SPECIE {
public:
	char name[50];
	int uid;
	MC_SPECIE() {strcpy(name, "\0"); uid = -1;};
};

class MC_SPECIES {
public:
	ARRAY<MC_SPECIE> sp;
	ARRAY<int> mtype; // 0 -- MM, 1 -- SM, 2 -- PM
	ARRAY<int> mindx; // the molecules to be randomly distributed, the molecule index of MM / SM / PM in the system
	float xf[3], xt[3]; // dimension of the cell, to randomly re-distribute the species 
	float dx, dy, dz; // dislocation distance

	ARRAY<int> order; // the order of species
	void init_order() {
		order.set_array(mtype.n);
		for (int i = 0; i < order.n; i++) order.m[i] = i;
	};
	void reorganize_order() {
		int m1, m2, nt;
		for (int i = 0; i < order.n; i++) {
			m1 = int(order.n * ranf() + 0.5); if (m1 == order.n) continue;
			m2 = int(order.n * ranf() + 0.5); if (m2 == order.n) continue;
			if (m1 != m2) { // exchange
				nt = order.m[m1];
				order.m[m1] = order.m[m2];
				order.m[m2] = nt;
			}
		}
	};
};

bool read_mc_species(char *fname, MC_SPECIES &sp);

bool read_mdprocs(char *fname, ARRAY<MD_PRO> &mdpro);

bool init_mc_species(MMOL_MD_CELL &mmcell, MC_SPECIES &sp, _bias_force_::BIAS_POINT &mdcell_bias);
void relocate_mc_species(MMOL_MD_CELL &mmcell, MC_SPECIES &sp, _bias_force_::BIAS_POINT &mdcell_bias);
