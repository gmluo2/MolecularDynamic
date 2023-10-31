// JOB TYPE
#define UNKNOWN_JOB                -1
// MD_MODE
#define SINGLE_MM_FREE_CELL_MD     0
#define MULTI_MM_FREE_CELL_MD      1
#define MULTI_MM_PERIDIC_CELL_MD   2
#define SINGLE_CGMM_FREE_CELL_MD   3
#define MULTI_CGMM_FREE_CELL_MD    4
#define MULTI_CGMM_PERIDIC_CELL_MD 5
#define GCMD                       6
// GET MD PROCEDURE
#define GET_SINGLE_MM_MD_PROC      11
#define SHOW_SINGLE_MM_MD_PROC     12
#define GET_MULTI_MM_MD_PROC       13
#define SHOW_MULTI_MM_MD_PROC      14

#define GET_SINGLE_CGMM_MD_PROC    15
#define SHOW_SINGLE_CGMM_MD_PROC   16
#define GET_MULTI_CGMM_MD_PROC     17
#define SHOW_MULTI_CGMM_MD_PROC    18

#define GET_GCMD_PROC              19
#define SHOW_GCMD_PROC             20

// MOLECULE OPERATION IN DLL
#define MOLECULE_OPERATION         21


#if _MAIN_PROG_ == 1
short iFormat_LJ = 1; // 0 -- Standard; 1 -- Alternative, used in Amber
short MM_SpeedVerlet_method = 1;  // 0 -- mass center spatial momentum; 1 -- speed of the base cluster
// unit transform : atom mass unit, Angs. & fs ==> J, sec
// mass matrix -- rotational moment
float U_InertiaMoment = 1.66053886f; // x 10^-17 J; I * 1/dt^2 -- (u*Angs^2) * (1/fs^2)
float U_ElectFr = 0.230726f; // x 10^-17 J; e * e / (4PI*eps0 * r^2) * r -- (electron^2) / (Angs)

float U_MassAccel = 1.66053886f; // x 10^-17 J / Angs ; m * a -- (u) * (Angs / fs^2)
float U_ElectForce = 0.230726f; // x 10^-17 J / Angs; e * e /(4PI*eps0 * r^2) -- (electron^2) / (Angs^2)

float T = 25;
float kT = (float)(1.38065e-6) * (273.15f + T); // x 10^-17 J @ 25C  == 0.41164 x 10^-20 J;  1 kT * 1mol = 2.478 kJ

float unit_Ek_kT = U_InertiaMoment / kT;

float U_eps = 0.001f / 6.022f; // x 10^-17 J; eps in unit kJ / mol
float U_LJForce = U_eps; // x 10^-17 J / Angs.; LJ force = eps / Angs.

// 1 e^2/(1 Ang) = U_ElectForce / k = 1.6711404 x 10^5 K = 167114.04 K
// 1kJ/mol = U_eps / k = 120.27508 K

float U_EextFr = (float)1.6e-6; // x 10^-17 J; external E in unit kV/mm, charge has unit e
float U_EextForce = (float)1.6e-6; // x 10^-17 J / Angs; external E in unit kV/mm, charge has unit e

float eUnit_estat_kT = U_ElectFr / kT;
float eUnit_ext_kT = U_EextFr / kT;
float eUnit_LJ_kT = U_eps / kT;
float fUnit_estat_mass = U_ElectForce / U_MassAccel;
float fUnit_ext_mass = U_EextForce / U_MassAccel;
float fUnit_ext_estat = U_EextFr / U_ElectFr;
float fUnit_LJ_mass = U_eps / U_MassAccel;

float fUnit_atm_u = 6.101935e-9; // 1 atmosphere = 101325 Pa (N/m^2) = 6.101935 x 10^-9 atomic_mass_unit (u)/ (fs^2 * A)
// for idea gas, each molecule has volumn 4.08841 x 10^4 Angs^3 at 1 atmosphere and 300K, 34.45^3 Angs^3

float Eext = 0; // external E, in unit kV/mm
bool bEext = false; // apply external electric field

float r2_bound = 4; // maximum distance^2 of bound C-C is about 1.54

char errmsg[2560] = "\0";
char atompar_fname[256] = "\0";
char refresh_mol_xyz_struct_fname[256] = "\0";

FLEX_MATRIX<LJ_PARS> LJPARS;
//TORSION_DATABASE torsion_db;
DIHEDRAL_DATABASE dh_db;

CLUSTER_IMPSOLV_PARS_DATABASE ESolv_db;
//FLEX_MATRIX<CLUSTER_IMPSOLV_PARS> ESolvPars;
float rSolvent = 1; // radius of solvent molecule
float dSolv_max = 0; // the maximum radius of cluster + the maximum radius of cluster + rSolvent

bool COMMAND_MD_STOP = false;

bool bRandomForce = false;
float randf = 0, fwhm_randf = 0;

bool scale_force = false;
float force_max = 0.1f, torque_max = 0.1f;
float force2_max = 0.01f, torque2_max = 0.01f;

bool bImplicitSolvent = false;
bool bPolarize = false;
bool bDipole = false;
bool bExcludeInducedDipole = false; // exclude the interactions between induced dipole?
float fp = 0.75f; // polarized dipole increasing speed, for polarization iteration
float muConvg = 0.0003f; // convergent value of atomic polarized dipole


bool local_interact = false;
bool local_relax = false;
float dta_max_relax = 0.002f;
float Ek0Internal = 0.5f; // in kT of each freedom for internal kinetic energy of macro-molecule

float HingeLength = 1.54f;
float CC_HingeLength = 1.54f;
float CO_HingeLength = 1.43f;

float max_ksi_Hoover = 0.01f;

float min_cmm_size = 20;

LOG mlog;

//float d2_ignore = 1.0f;  //1.0f

float eps_dc = 1; // the relative dielectric constant for electro-static calculation
float d_14 = 6, d2_14 = d_14 * d_14; // the maximum distance of the 1-4 sites -- 4 * 1.5 = 6

float size_cluster = 3;
float rcut_Ewd = 9, r2cut_Ewd = rcut_Ewd * rcut_Ewd, rcut_LJ = 15, r2cut_LJ = rcut_LJ * rcut_LJ;
float rcut_LJ_norm = 2.5f; // normalized cutoff distance for L-J interaction, rc = rcut_LJ_norm * sigma

short iSurfaceBoundary = 0; // extract the surface dipole effect on the simulation

namespace _EwaldSum_recp_ {
#if _EWALD_SUM == _MultiPole_EWALD_SUM
	// Multipole Ewald-Sum parameters
	VECTOR3 h[NH][NH][NH];
	double H2[NH][NH][NH];
	double Ah[NH][NH][NH];
	float AhMin = (float)1e-8;
#endif

} // end of namespace _EwaldSum_recp_

int iPolarize = 0; // polarization model -- point-charge/dipole
namespace _Thole_polarize_ {
	bool bTholePolarize = false;
	short iTholeModel = 0;  // 0 -- exponential screening, 1 -- exponential-3 (r^3) screening
	float a = 2.1304f;
	float xcut_exp = 15; // cut of the exponential y = exp(-x) x = 15, y = 3e-7; x = 12; y = 1e-5; x = 10, y = 2.7e-5
} // end of _Thole_polarize_

int nmax_loops_polarize = 8;

bool bRealEwaldSumOnly = false;

// SPME Ewald-Sum parameters
int indx_bspline = 4;
int dim_spme[3] = {32, 32, 32};
bool bExpandSPMECellZ = false;
// end of SPME Ewald-Sum parameters

int interact_cal_method = 0;

//Brownian friction force -- Drag 
bool bBrownian = false;
float ita_Brownian = 0.f;  // linear coeff
float cd_drag_trans = 0.f;  // drag force at high velocity, non-linear part
float cd_drag_rot = 0.f;  // drag force at high velocity, non-linear part

// for MPI -- the range of clusters calculated by this computer
//            for single macromolecule -- cluster range
//            for multi-macromolecule -- macromolecule
#if _MPI_ == 1
int mpi_id = 0, nMPI = 0; // rank id in MPI, and the size
int *n1MM = NULL, *n2MM = NULL;
int *n1SM = NULL, *n2SM = NULL;
//int *n1CMM = NULL, *n2CMM = NULL;
char proc_name[20][50];
#endif

// some interaction profiles
EF_DPROF LJ12_6;

bool bCGdipole = true; // adding dipole moment to CG_BOUND
float rcut_intramol = 10, r2cut_intramol = 100; // r2cut for intramol-interaction in CGMM

bool bVirial = false; // calculate virial coefficient / Stress Tensor ?
short method_Virial = 0; // cluster virial -- 0, molecular virial -- 1

char seperator[6] = " ,;\t|";

// add a spring force along z axis on the mass center of the slab
bool bZSpring = false;
float k_zspring[2] = {0, 0};
int mode_NetForce; // 0 -- no action; 1 -- compensating netforce along x/y/z axises, 2 -- compensating netforce along z-axis
float ita_mc = 0; // we add a friction force on the mass center, to avoid high speed along z axis

// calculate short-range force, electrostatic force, ... in parallel, and the forces are to be saved in different memory
bool bMTForce = true;
#ifndef _NFS
	#define _NFS    2
#endif

double timestep_md;

namespace _bias_force_ {
#define NPAR_BIAS 5
	float par[NPAR_BIAS];

	bool biasOn = false;

#define BIAS_FORCE 0
}

#else
extern short iFormat_LJ;
extern short MM_SpeedVerlet_method;
extern float U_InertiaMoment;
extern float U_ElectFr;
extern float U_MassAccel;
extern float U_ElectForce;
extern float T;
extern float kT;
extern float U_eps;
extern float U_LJForce;
extern float U_EextFr;
extern float U_EextForce;
extern float unit_Ek_kT;
extern float eUnit_estat_kT;
extern float eUnit_ext_kT;
extern float fUnit_estat_mass;
extern float fUnit_ext_mass;
extern float fUnit_ext_estat;
extern float fUnit_atm_u; // from atmosphere to u/(fs^2*A) -- pressure

extern float eUnit_LJ_kT;
extern float fUnit_LJ_mass;


extern float Eext;
extern bool bEext;

extern float r2_bound;

extern char errmsg[2560];
extern char refresh_mol_xyz_struct_fname[256];

extern ATOM_COMM_PARS_DATABASE atompar_db;
extern FLEX_MATRIX<LJ_PARS> LJPARS;
//extern TORSION_DATABASE torsion_db;
extern DIHEDRAL_DATABASE dh_db;

extern CLUSTER_IMPSOLV_PARS_DATABASE ESolv_db;
//extern FLEX_MATRIX<CLUSTER_IMPSOLV_PARS> ESolvPars;
extern float rSolvent;
extern float dSolv_max; // the maximum radius of cluster + the maximum radius of cluster + rSolvent

extern bool COMMAND_MD_STOP;

extern bool bRandomForce;
extern float randf, fwhm_randf;

extern bool scale_force;
extern float force_max, torque_max;
extern float force2_max, torque2_max;

extern bool bImplicitSolvent;
extern bool bPolarize;
extern bool bDipole;
extern bool bExcludeInducedDipole; // exclude the interactions between induced dipole?
extern float fp; // polarized dipole increasing speed, for polarization iteration
extern float muConvg; // convergent value of atomic polarized dipole


extern bool local_interact;
extern bool local_relax;
extern float dta_max_relax;
extern float Ek0Internal; // in kT of each freedom for internal kinetic energy of macro-molecule

extern LOG mlog;

extern float HingeLength;
extern float CC_HingeLength;
extern float CO_HingeLength;

extern float max_ksi_Hoover;

extern float min_cmm_size;

//extern float d2_ignore;

extern float eps_dc;
extern float d_14, d2_14; // the maximum distance of the 1-4 sites -- 4 * 1.5 = 6
extern float size_cluster;
extern float rcut_Ewd, r2cut_Ewd, rcut_LJ, r2cut_LJ, rcut_LJ_norm;

extern short iSurfaceBoundary;

namespace _EwaldSum_recp_ {

#if _EWALD_SUM == _MultiPole_EWALD_SUM

	extern VECTOR3 h[NH][NH][NH];
	extern double H2[NH][NH][NH];
	extern double Ah[NH][NH][NH];
	extern float AhMin;
#endif

} // end of namespace _EwaldSum_recp_ 

extern int iPolarize; // polarization model -- point-charge/dipole
namespace _Thole_polarize_ {
	extern bool bTholePolarize;
	extern short iTholeModel;  // 0 -- exponential screening, 1 -- exponential-3 (r^3) screening
	extern float a; // = 2.1304;
	extern float xcut_exp; // cut of the exponential exp(-x)
} // end of _Thole_polarize_

extern int nmax_loops_polarize;


// SPME Ewald-Sum parameters
extern int indx_bspline;
extern int dim_spme[3];
extern bool bExpandSPMECellZ;
// end of SPME Ewald-Sum parameters

extern int interact_cal_method;
extern bool bRealEwaldSumOnly;

extern bool bBrownian;
extern float ita_Brownian;
extern float cd_drag_trans;
extern float cd_drag_rot;

#if _MPI_ == 1
extern int mpi_id, nMPI; // rank id in MPI, and the size
extern int *n1MM, *n2MM;
extern int *n1SM, *n2SM;
//extern int *n1CMM, *n2CMM;
extern char proc_name[20][50];
#endif

extern bool bCGdipole;
extern float rcut_intramol, r2cut_intramol; // r2cut for intramol-interaction in CGMM

extern bool bVirial; // calculate virial coefficient / Stress Tensor ?
extern short method_Virial;

extern char seperator[6];

// add a spring force along z axis on the mass center of the slab
extern bool bZSpring;
extern float k_zspring[2];
extern int mode_NetForce; // 0 -- no action; 1 -- compensating netforce along x/y/z axises, 2 -- compensating netforce along z-axis
extern float ita_mc;

// calculate short-range force, electrostatic force, ... in parallel, and the forces are to be saved in different memory
extern bool bMTForce;
#ifndef _NFS
	#define _NFS    2
#endif

extern double timestep_md;

namespace _bias_force_ {
#define NPAR_BIAS 5
	extern float par[NPAR_BIAS];
#define BIAS_FORCE 0

	extern bool biasOn;
}


#endif

