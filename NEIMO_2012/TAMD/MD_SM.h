bool SM_LeapFrog_velocity_acceleration(bool bHooverDyn, SMOLECULE &m, LFVERLET_SM_MD& mdpar, int md_loop, int &nconvg_loops, char *msg_show, char *msg_log, int n_msglen);
bool SM_SpeedVerlet(bool local_relax, bool speed_verlet, SMOLECULE &m, LFVERLET_SM_MD& mdpar, double& vmax_diff, double& wmax_diff, char *msg_show, char *msg_log, int n_msglen);
bool SM_SpeedVerlet_Simplified(bool speed_verlet, SMOLECULE &m, LFVERLET_SM_MD& mdpar, double& vmax_diff, double& wmax_diff);

bool PM_LeapFrog_velocity_acceleration(bool bHooverDyn, PMOLECULE &m, LFVERLET_PM_MD& mdpar, int md_loop, int &nconvg_loops, char *msg_show, char *msg_log, int n_msglen);
bool PM_SpeedVerlet(bool local_relax, bool speed_verlet, PMOLECULE &m, LFVERLET_PM_MD& mdpar, double& vmax_diff, char *msg_show, char *msg_log, int n_msglen);
bool PM_SpeedVerlet_Simplified(bool speed_verlet, PMOLECULE &m, LFVERLET_PM_MD& mdpar, double& vmax_diff);
