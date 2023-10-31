class GC_SVar {
public:
	double g; // fractal coefficient
	double r0; // constant length used for distance correction
	// in Lo-Palmer's technique, L-J interaction is modified as U = g * 4 * eps * [(sigma/(r + (1 - g) * r0))^12 - (sigma / (r + (1 - g) * r0))^6]
	// in another word, eps is modified to g * eps, distance is modified to (r + (1 - g) * r0)
	double rc; // correction of the distance, rc = (1 - g) * r0

	GC_SVar() {g = 0; r0 = 0; rc = 0;};
	double set(double g, double r0) {this->g = g; this->r0 = r0; rc = (1 - g) * r0; return rc;};
	double r_c() {return rc;};
};

class GC_InteractVar : public InteractVar, public GC_SVar {
public:
	GC_InteractVar() {};
};

class GC_Res {
public:
	double fg; // force on fractal coefficient
	void reset() {fg = 0;};
	GC_Res() {reset();};
	GC_Res(double v) {reset();};
	void operator = (const GC_Res &var) {
		this->fg = var.fg;
	};
	void operator = (double v) {
		fg = v;
	};
	GC_Res operator + (GC_Res &var) {
		GC_Res res;
		res.fg = fg + var.fg;
		return res;
	};
	void operator += (GC_Res &var) {
		fg += var.fg;
	};
};

class GC_InteractRes : public InteractRes, public GC_Res {
public:
	void reset() {InteractRes::reset(); GC_Res::reset();};
	GC_InteractRes() {reset();};
	GC_InteractRes(double v) {reset();};
	void operator = (const GC_InteractRes &var) {
		*((InteractRes*)this) = *((InteractRes*)(&var));
		*((GC_Res*)this) = *((GC_Res*)(&var));
	};
	void operator = (double v) {
		*((InteractRes*)this) = v;
		*((GC_Res*)this) = v;
	};
	GC_InteractRes operator + (GC_InteractRes &var) {
		GC_InteractRes res;
		*((InteractRes*)(&res)) = *((InteractRes*)this);
		*((InteractRes*)(&res)) += *((InteractRes*)(&var));
		*((GC_Res*)(&res)) = *((GC_Res*)this);
		*((GC_Res*)(&res)) += *((GC_Res*)(&var));
		return res;
	};
	void operator += (GC_InteractRes &var) {
		*((InteractRes*)this) += *((InteractRes*)(&var));
		*((GC_Res*)this) += *((GC_Res*)(&var));
	};
};

void GC_LJ_Interact_Correct(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, GC_InteractVar &var, GC_InteractRes &res);


class GC_LVar {
public:
	double g; // fractal coefficient
	double r0; // constant length used for distance correction
	// similiar to Lo-Palmer's modified short-range L-J interaction, electrostatic interaction is modified as:
	// R ==> R + (1 - g) * erfc(R/r0)
	double r_c, f_c; // r_c is the additional distance to R, (1-g) * erfc(R/r0); f_c is the correction on force, d{r_c} / dR
	GC_LVar() {g = 0; r0 = 1; r_c = 0; f_c = 0;};
	double set(double g, double r0) {this->g = g; this->r0 = r0;};
	void set_dist(double R, bool bRcOnly = false) {
		double r = R / r0;
		if (r > 5) {r_c = 0; f_c = 0; return;}
		int i;
		double v;
		ERFC(v, r, i);
		r_c = (1 - g) * v / r0;
		if (bRcOnly) return;
		EXP2(v, r, i)
		f_c = (g - 1) * v * 2 * INV_SQRT_PI;
	};
};


#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

namespace _EwaldSum_real_ {
class GC_SPME_VAR : public SPME_VAR, public GC_LVar {
public:
	GC_SPME_VAR() {};
};

void GC_Electrostatic_Interact_Correct(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, GC_SPME_VAR &var, GC_InteractRes &res);
} // end of namespace _EwaldSum_real_

#endif // _IGNORE_ELECTROSTATIC_INTERACTION_
