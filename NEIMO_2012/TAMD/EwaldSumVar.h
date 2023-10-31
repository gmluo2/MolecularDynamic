namespace _EwaldSum_real_ {

#ifndef _KAPPA
#define _KAPPA     3.2
#endif


// the variable which is used for Ewald-Sum for real term calculation
class EwaldSumRealVars0 { // for charge - charge interaction
public:
	VECTOR3 mT1, mET1;
	double* dr;
	double R, R2, R3;
	double kR;
	double fexp2, f0, f1;

	bool bEwald;
	double kappa, kappa3;
	short iSurface; // 0 -- cubic use 4 * PI / 3V,  1 -- slab use 4PI/V
	short iSurfaceBoundary; // cubic or 2d cell 0/1
	short iExIndDipole; // excluding induced dipole -- 0/1
	double inv_3V; // 4 * PI / (3 * V)
	double inv_V; // 4 * PI / V

	EwaldSumRealVars0() {bEwald = true; inv_3V = 0; iSurface = 0; inv_V = 0; iExIndDipole = 0;};
	void init_vars(VECTOR3& dr, double R, double R2);
	void init_vars(double* dr, double R, double R2) {VECTOR3 dR(dr); init_vars(dR, R, R2);};
	void init_EwaldSum(float V, float rcut_Ewd) {
		inv_V = 4 * PI / V;
		inv_3V = inv_V / 3;
		kappa = _KAPPA / rcut_Ewd;
		kappa3 = kappa * kappa * kappa;
	};
	void cp_EwaldSum_vars(EwaldSumRealVars0 &var) {
		bEwald = var.bEwald;
		iSurfaceBoundary = var.iSurfaceBoundary;
		iSurface = var.iSurface;
		iExIndDipole = var.iExIndDipole;
		inv_3V = var.inv_3V;
		inv_V = var.inv_V;
		kappa = var.kappa;
		kappa3 = var.kappa3;
	};
};

class EwaldSumRealVars1 : public EwaldSumRealVars0 {
public:
	// more variable for dipole interaction
	DMATRIX<3> mT2, mET2, mD2;
	CMATRIX<3> mT3, mET3;
	//QMATRIX<3> mT4, mET4;

	double R4, R5, R7;//, R9;
	double kR2;
	
	double f2, f3, f4;//, f5, f6;

	EwaldSumRealVars1() {};
	void init_vars(VECTOR3& dr, double R, double R2);
	void init_vars(double* dr, double R, double R2) {VECTOR3 dR(dr); init_vars(dR, R, R2);};
	void init_vars_E(VECTOR3& dr, double R, double R2); // for E calculation only
	void init_vars_E(double *dr, double R, double R2) {VECTOR3 dR(dr); init_vars_E(dR, R, R2);}; // for E calculation only
};

class EwaldSumRealVars2 : public EwaldSumRealVars1 {
public:
	// more variable for dipole-quardupole interaction
	QMATRIX<3> mT4, mET4;

	double R9, f5, f6;

	EwaldSumRealVars2() {};
	void init_vars(VECTOR3& dr, double R, double R2);
};

// q1 - q2 of 1
void MTP_Interact(EwaldSumRealVars0 &esv, double c1, double c2, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a);
// q1, mu1 -- q2, mu2 of 1
void MTP_Interact(EwaldSumRealVars1 &esv, double c1, double *mu1, double c2, double *mu2, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a);
// q1, mu1 -- q2, mu2, Q of 1
void MTP_Interact(EwaldSumRealVars2 &esv, double c1, double *mu1, double c2, double *mu2, DMATRIX<3>& Q, VECTOR3 &E, VECTOR3 &F, double& U, bool bEwaldKickoff, double a);

// E -- c2, mu2 and/or Q
void cMTP_E(EwaldSumRealVars0 &esv, double c2, VECTOR3 &E, bool bEwaldKickoff, double a);
void MTP_E(EwaldSumRealVars1 &esv, double c2, double *mu2, VECTOR3 &E, bool bEwaldKickoff, double a);
void MTP_E(EwaldSumRealVars2 &esv, double c2, double *mu2, DMATRIX<3>& Q, VECTOR3 &E, bool bEwaldKickoff, double a);

// this is the total dipole of a cell
class SURF_DIPOLE {
public:
	float q;
	VECTOR3 qr, mu0; // total dipole moment of the whole cell, SUM(q * r) and SUM(mu0) respectively
	VECTOR3 mu; // total dipole moment of the whole cell : SUM(q * r + mu)
	VECTOR3 mu_ind; // total induced dipole moment of the whole cell: SUM(mu_ind)
	SURF_DIPOLE() {q = 0;};
	SURF_DIPOLE(double v) {reset();};
	void reset() {q = 0; V3zero(qr) V3zero(mu0) V3zero(mu) V3zero(mu_ind)};
	void reset_z() {
		//qr.v[2] = 0; mu0.v[2] = 0; mu.v[2] = 0; mu_ind.v[2] = 0;
	};
	void operator = (double v) {
		reset();
	};
	void operator += (SURF_DIPOLE &d) {
		q += d.q; V3plusV3(qr, d.qr, qr) V3plusV3(mu0, d.mu0, mu0) V3plusV3(mu, d.mu, mu) V3plusV3(mu_ind, d.mu_ind, mu_ind)
	};
};


} // end of namespace _EwaldSum_real_



// With Thole polarization model, the ion has charge distributed around the center point
// this distribution will induce the charge-charge and charge-dipole interactions different from point charge/dipole interaction
// Additional correction within a short range is required for charge-charge, charge-dipole interactions
namespace _Thole_polarize_ {
extern float xcut_exp;
class TholeCorrTensor {
public:
	int iTholeModel; // 0 -- r-exponential screening, exp(-a*u); 1 -- r^3-exponential screening, exp(-a * u^3)
	double f;
	VECTOR3 mT1;
	DMATRIX<3> mT2;
	CMATRIX<3> mT3;
	double xcut_exp;

	// Correction: [in format of paper by Andres Aguado et al, J. Chem. Phys., 119, 7471 (2003)]
	// U = f * q1 * q2;
	// E & F are from 2 to 1, with dr = r1 - r2 (from r2 to 1)
	// E = -q2 * mT1 + mT2 * mu2;  
	// F = -q1 * q2 * mT1 + q1 * mT2 * mu2 - q2 * mT2 * mu1

	TholeCorrTensor() {
		iTholeModel = 0; 
		f = 0; this->xcut_exp = ::_Thole_polarize_::xcut_exp;
	};

	bool init_vars_0(double a, double R, double R2, double *dr, bool bF = true);
	bool init_vars_1(double a, double R, double R2, double *dr, bool bF = true);
	bool init_vars(double a, double R, double R2, double *dr, bool bF = true) {
		if (iTholeModel == 0) return init_vars_0(a, R, R2, dr, bF);
		else if (iTholeModel == 1) return init_vars_1(a, R, R2, dr, bF);
		else return false;
	};
};

bool TholeCorr_E(TholeCorrTensor &tt, double a, double q2, double *mu2, double R, double R2, double *dr, double *E);
bool TholeCorr(TholeCorrTensor &tt, double a, double q1, double *mu1, double q2, double *mu2, double R, double R2, double *dr, double *E, double *F, double &U);

} // end of namespace _Thole_Polarize_

