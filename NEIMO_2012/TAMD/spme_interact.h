#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0

namespace _EwaldSum_real_ {
class SPME_VAR : public InteractVar {
public:
	bool bEfieldOnly;
	_EwaldSum_real_::EwaldSumRealVars1 esv1;
	bool bTholePolarize;
	_Thole_polarize_::TholeCorrTensor tv;

	SPME_VAR() {bEfieldOnly = false; bTholePolarize = false;};
	void operator = (SPME_VAR &var) {
		this->bEfieldOnly = var.bEfieldOnly;
		(*(InteractVar*)this) = (*(InteractVar*)(&var));
		esv1.cp_EwaldSum_vars(*(EwaldSumRealVars0*)(&var.esv1));
		this->bTholePolarize = var.bTholePolarize;
		this->tv.iTholeModel = var.tv.iTholeModel;
		this->tv.xcut_exp = var.tv.xcut_exp;
	};
};

// interaction on patom1 from patom2, 
// dr -- vector from patom2 to patom1
// R -- distance |dr|
// R -- |dr|^2
// E -- electric field from patom2
// F -- force from patom2 to patom1
// U -- interaction energy
void EwaldReal(EwaldSumRealVars1 &esv, MATOM *patom1, MATOM *patom2, VECTOR3& dr, double R, double R2, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a, bool bExcludeInducedDipole);
// calculate the electric field only about the real part of Ewald Sum
void EwaldReal_E(EwaldSumRealVars1 &esv, MATOM *patom1, MATOM *patom2, VECTOR3& dr, double R, double R2, VECTOR3 &E, bool bEwaldKickoff, double a);


void MTP_Interact(EwaldSumRealVars1 &esv, MATOM *patom1, MATOM *patom2, VECTOR3& dr, double R, double R2, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a, bool bExcludeInducedDipole);
void MTP_E(EwaldSumRealVars1 &esv, MATOM *patom2, VECTOR3& dr, double R, double R2, VECTOR3 &E, bool bEwaldKickoff, double a);




void CLUSTER_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, SPME_VAR &var, InteractRes &res);
void EFIELD_SPME_REAL_CMM(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, SPME_VAR &var, InteractRes &res);
void SingleMM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, SPME_VAR &var, InteractRes &res);
void MM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE *mm, int nMM, SPME_VAR &var, InteractRes &res);
void SM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE *sm, int nSM, SPME_VAR &var, InteractRes &res);
void PM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, PMOLECULE *pm, int nPM, SPME_VAR &var, InteractRes &res);

void EwaldSum_SurfaceDipole(_EwaldSum_real_::SURF_DIPOLE &mu_surf, EwaldSumRealVars0 *esv, InteractRes &res);

#ifdef _VMD_CELL_
	void VMM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE **mm, int nMM, SPME_VAR &var, InteractRes &res);
	void VSM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE **sm, int nSM, SPME_VAR &var, InteractRes &res);
	void VPM_EFIELD_SPME_REAL(CMM_CELL3D<BASIC_CLUSTER> &cell, PMOLECULE **pm, int nPM, SPME_VAR &var, InteractRes &res);
#endif

} // end of namespace _EwaldSum_real

namespace _spme_ {

	void LocalEF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res);

	void LocalEF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res);

	class SPME_REAL_VARS {
	public:
		MMOL_MD_CELL *mdcell;
		CMM_CELL3D<BASIC_CLUSTER> *cmm_cell;
		SPME_VAR* var;
		InteractRes* res;
		int nThreads;
		_THREAD_VARS thread_var;
	};

	void Asyn_SPME_EF_Real(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, int nThreads, SPME_VAR& var, InteractRes& res, SPME_REAL_VARS &svar);

// calculate E & F for each atom
void SPME_EF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);
// the function to calculate the interaction with SPME method, the final force will be accumulate to the cluster force & torque
void SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);

// calculate E & F for each atom
void SPME_EF(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);
// the function to calculate the interaction with SPME method, the final force will be accumulate to the cluster force & torque
void SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);
// polarization and the interaction with SPME method, the final force will be accumulate to the cluster force & torque
bool Polarize_SPME_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);
// excluding the interactions between induced dipoles ?
bool Polarize_SPME_Interact2(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, VSPME<_VMUatom>& vspme_induced, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);

#ifdef _VMD_CELL_
	void SPME_EF(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);
	void SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VQatom>& vspme, int nThreads, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);
	void SPME_EF(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);
	void SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);
	bool Polarize_SPME_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, _EwaldSum_real_::SPME_VAR &var, InteractRes &res);
#endif
} // end of namespace _spme_

#endif

void CLUSTER_TOTAL_FORCE(BASIC_CLUSTER *pc);
void SingleMM_TOTAL_FORCE(MMOLECULE *mm, int nc1, int nc2);
void MM_TOTAL_FORCE(MMOLECULE *mm, int nMM);
void SM_TOTAL_FORCE(SMOLECULE *sm, int nSM);
void PM_TOTAL_FORCE(PMOLECULE *pm, int nPM);

#ifdef _VMD_CELL_
	void VMM_TOTAL_FORCE(MMOLECULE **mm, int nMM);
	void VSM_TOTAL_FORCE(SMOLECULE **sm, int nSM);
	void VPM_TOTAL_FORCE(PMOLECULE **pm, int nPM);
#endif

void CLUSTER_LJ_INTERACT(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, InteractVar &var, InteractRes &res);
void SingleMM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE &mm, int n1, int n2, InteractVar &var, InteractRes &res);
void MM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE *mm, int nMM, InteractVar &var, InteractRes &res);
void SM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE *sm, int nSM, InteractVar &var, InteractRes &res);
void PM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, PMOLECULE *pm, int nPM, InteractVar &var, InteractRes &res);
void LJ_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, int nThreads, InteractVar &var, InteractRes &res);


class LJ_AsynVARS {
public:
	MMOL_MD_CELL *mdcell;
	CMM_CELL3D<BASIC_CLUSTER> *cmm_cell;
	InteractVar *var;
	InteractRes* res;
	int nThreads;
	_THREAD_VARS thread_var;
};

void Asyn_LJ_Interact(MMOL_MD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, int nThreads, InteractVar &var, InteractRes &res, LJ_AsynVARS &av);

#ifdef _VMD_CELL_
	void VMM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, MMOLECULE **mm, int nMM, InteractVar &var, InteractRes &res);
	void VSM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, SMOLECULE **sm, int nSM, InteractVar &var, InteractRes &res);
	void VPM_LJ(CMM_CELL3D<BASIC_CLUSTER> &cell, PMOLECULE **pm, int nPM, InteractVar &var, InteractRes &res);
	void LJ_Interact(MMOL_VMD_CELL *mdcell, CMM_CELL3D<BASIC_CLUSTER> *cmm_cell, int nThreads, InteractVar &var, InteractRes &res);
#endif
