namespace _pmf_constraint_ {

// type of PMF_CONSTRAINT
#define PMF_CONSTRAINT_0     1
#define PMF_CONSTRAINT_1     2
#define PMF_CONSTRAINT_2     3
#define PMF_CONSTRAINT_3     4

class PMF_CONSTRAINT {
public:
	// constraint method 1: -- constraint on two PMs, the distance between the two PMs is limited with constraint potential 
	//                         -- 1/2 * kc * dx^2
	// constraint method 2: -- constraint on two PMs, the distance between the two PMs is limited with constraint potential 
	//                         -- 1/2 * kc * dx^2 + (a * dx)
	// constraint method 3: -- constraint on two PM, with one PM position fixed, another one limited with constraint potential 
	//                         -- 1/2 * kc * dx^2
	// constraint method 4: -- constraint on two PM, with one PM position fixed, another one limited with constraint potential 
	//                         -- 1/2 * kc * dx^2 + (a * dx)
	int cModel; 
	ARRAY<int> v_int; // integer variables
	ARRAY<float> v_f; // float variable

	float dr_max, r0;
	EF_DPROF ef_prof;

	void release_variables() {v_int.release(); v_f.release();};
	void set_model(int model) {
		cModel = model;
		switch (model) {
		case PMF_CONSTRAINT_0:
			v_int.SetArray(2);
			v_f.SetArray(1); // kc
			break;
		case PMF_CONSTRAINT_1:
			v_int.SetArray(2);
			v_f.SetArray(2); // kc, a
			break;
		case PMF_CONSTRAINT_2:
			v_int.SetArray(2);
			v_f.SetArray(1); // kc
			break;
		case PMF_CONSTRAINT_3:
			v_int.SetArray(2);
			v_f.SetArray(2); // kc, a
			break;
		default:
			show_log("Unknown constraint model", true);
		}
	};
	void construct_constraint_efprofile();
	void ef(float R, double &U, double &f);
};

// Ion 0 and Ion 1 index of PM
#define I0  0
#define I1  1
// float variable
#define KC       0
#define A        1

bool read_pmf_constraint(char *par_fname, PMF_CONSTRAINT &cpmf);

} // end of namespace _pmf_constraint_
