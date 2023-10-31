void max_diff(VECTOR<6> *v1, VECTOR<6> *v2, int n, double &wmax, double &umax);

// initilize the local cordinates and calculate the inertia moment of each cluster and the whole molecule
void CalMolInertiaMoment(MMOLECULE &mm, int nThreads = 1);

// calculate Coriolis acceleration term and gyroscopic spatial force
void CalMol_ab(MMOLECULE &mm);

//adjusting free base velocity so that the total spacial moment is consistent
void MakeFreeBaseVelocityConsistent(MMOLECULE &mm);

// get self-consistent velocity and accleration based on Leap-Frog MD alghgorithm, Ethermal_trans and Ethermal_rot are in unit kT
bool LeapFrog_SelfConsistent_velocity_acceleration(bool bHooverDyn, MMOLECULE &mm, LFVERLET_MM_MD& mdpar, int md_loop, int &nconvg_loops, VECTOR<6> *vect1, VECTOR<6> *vect2, double *dv1, double *dv2, bool local_relax, char *msg_show, char *msg_log, int n_msglen, bool init_md = false);

bool MM_SpeedVerlet(bool local_relax, bool speed_verlet, MMOLECULE &mm, LFVERLET_MM_MD& mdpar, VECTOR<6> *vect1, VECTOR<6> *vect2, double &vmax_diff, double &wmax_diff, char *show_msg, char *log_msg, int n_msglen);

bool MM_SpeedVerlet_Simplified(bool speed_verlet, MMOLECULE &mm, LFVERLET_MM_MD& mdpar, VECTOR<6> *vect1, VECTOR<6> *vect2, double &vmax_diff, double &wmax_diff);

// Initilize the clusters' velocity
void InitVelocity(MMOLECULE &mm);

void MD_POLYMER(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, CMM_CELL3D<BASIC_CLUSTER> *cmm, int interact_cal_method, int nNormalMD);

void show_loop_infor(int nloop, int nsaved, int nLoops, float Ek, float Ep);

void show_molecule_struct(MMOLECULE &mm, char *fname);

void show_all_molecules();

void FIXEDBASE_MD_LPOLYMER(MMOLECULE &mm, LFVERLET_MM_MD& mdpar, CMM_CELL3D<BASIC_CLUSTER> *cmm, int interact_cal_method);





template <class MM_JOB, class CMM> class MM_CMM_JOB : public JOB_THREAD_VARS<MM_JOB> {
public:
	CMM* cmm;
	bool status;
	MM_CMM_JOB() {cmm = NULL; status = true;};
	~MM_CMM_JOB() {cmm = NULL;};
};

template <class MM_JOB, class CMM> class MM_CMM_INTERACT_JOB : public JOB_THREAD_VARS<MM_JOB> {
public:
	CMM *cmm;
	InteractVar var;
	EInteractVar evar;
	InteractRes LJ_res, eres;

	MM_CMM_INTERACT_JOB() {cmm = NULL; evar.esv1.bEwald = false;}; // NOT Ewald Sum, but the direct interaction
	~MM_CMM_INTERACT_JOB() {cmm = NULL;};
};

#define MM_INTERACT_JOB   MM_CMM_INTERACT_JOB< MACROMOL_THREAD_VARS<MMOLECULE>, CMM_CELL3D<BASIC_CLUSTER> >

void MM_ClusterAllocate2CMM(MM_CMM_JOB<MACROMOL_THREAD_VARS<MMOLECULE>, CMM_CELL3D<BASIC_CLUSTER> > *job);

void MM_FreeCellInteract(MM_INTERACT_JOB *job);
void MM_ClusterTotalForce(MACROMOL_THREAD_VARS<MMOLECULE> *job);
