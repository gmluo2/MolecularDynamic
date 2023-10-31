// given inertial cordinate of clusters
// Jacabino matrix, inertia mass matrix, inertia spatial matrix
// torque, a, b for each cluster
// to calculate the acceleration of torque angle 
// ksi is the coefficient of NEIMO-Hoover dynamic
void Hoover_Dyn(MMOLECULE &mm, float ksi);

// calculate the matrix defined in NEIMO procedure
void NEIMO_LINEAR_POLYMER(MMOLECULE &mm);

// calculate the matrix defined in NEIMO procedure
void NEIMO_MATRIX(MMOLECULE &mm, bool new_config, bool bAcceleration);

// calculate free base acceleration 
void NEIMO_CalFreeBaseAccel(MMOLECULE &mm, bool bCalMatrix);

// calculate the last step in NEIMO procedure -- cluster acceleration
void NEIMO_CalClusterAcceleration(MMOLECULE &mm, bool new_config);

// Hoover motion will induce correction on the acceleration of child cluster
// to make sure the inertia moment of whole molecule conservation,
// the base cluster will have correction on its acceleration
void NEIMO_CalHooverForce(MMOLECULE &mm, float ksi, VECTOR<6> &force_mc);
void distribute_MassCenter_Hoover_force(MMOLECULE &mm, int nThreads = 1); // distributing (adding) Hoover force at mass center to each cluster

void NEIMO(MMOLECULE &mm, bool bCalMatrix);

// same
void NEIMO_CalHooverForce_1(MMOLECULE &mm, VECTOR<6> &force_mc); // ksi = 1, exactly same as NEIMO_CalInternalHooverForce
void distribute_MassCenter_InternalKinetic_Hoover_force(MMOLECULE &mm); //exact same as distribute_MassCenter_Hoover_force
void NEIMO_GivenHoover(MMOLECULE &mm, float ksi, bool new_config);

void NEIMO_FixedBase(MMOLECULE &mm, float ksi, bool bCalMatrix);

// this is a combination of functions NEIMO_MATRIX, and NEIMO_CalFreeBaseAccel
// about the part for calculation NON-related to the force and/or acceleration
void NEIMO_MOMENT_MATRIX(MMOLECULE &mm);


class NEIMO_THREAD_VARS : public _THREAD_VARS {
public:
	MMOLECULE *mm;
	NEIMO_THREAD_VARS() {mm = NULL;};
	void set_pars(MMOLECULE *m, int iThread) {this->mm = m; _THREAD_VARS::set_pars(iThread);};
	~NEIMO_THREAD_VARS() {mm = NULL;};
};

void Start_NEIMO_MATRIX_CALC(MMOLECULE &mm, NEIMO_THREAD_VARS &var);

void NEIMO_HingeForce(MMOLECULE &mm);
