namespace _evatom_ {
class _VQatom {
public:
	bool eNeutral;
	short mMP;
	float *q; // charge of the real atom
	VECTOR3 *r; // location of the real atom
	VECTOR3 *E; // variable for electric field of the real atom
	VECTOR3 *F; // variable for the force of the real atom

	VECTOR3 u; // the normalized cordinate in SPME

	VECTOR3 E_recp; // electric field induced from the dipole moment
	VECTOR3 F_recp; // the force from the receiprical space

	ARRAY<VECTOR3> aE; // additional electric field term
	ARRAY<VECTOR3> aF; // additional force term

	_VQatom() {eNeutral = false; mMP = 0; q = NULL; r = NULL; E = NULL; F = NULL;};
	~_VQatom() {q = NULL; r = NULL; E = NULL; F = NULL;};
	void reset_local_EF(bool bF = true) {
		int i;
		V3zero(E_recp) if (bF) {V3zero(F_recp)}
		for (i = 0; i < aE.n; i++) {
			V3zero(aE.m[i]) if (bF) {V3zero(aF.m[i])}
		}
	};
};

class _VMUatom : public _VQatom {
public:
	VECTOR3 *mu; // dipole moment

	_VMUatom() {mu = NULL;};
	~_VMUatom() {mu = NULL;};
};

template <class atom> class _VE_CELL {
public:
	ARRAY<atom> va;
	float xl[3], xd[3]; // xl -- left bottom corner; xd -- width of the cell
	float V; // volumn
	double inv_V; //invV = 4 * PI / V; !!!!!!  NOT 4 * PI / (3 * V)
	double kappa;

	void reset_local_EF(bool bF = true) {
		int i, j;
		for (i = 0; i < va.n; i++) {
			V3zero(va.m[i].E_recp) if (bF) {V3zero(va.m[i].F_recp)}
			for (j = 0; j < va.m[i].aE.n; j++) {V3zero(va.m[i].aE.m[j]) if (bF) {V3zero(va.m[i].aF.m[j])}}
		}
	};
	void reset_local_E() {
		int i, j;
		for (i = 0; i < va.n; i++) {
			V3zero(va.m[i].E_recp) 
			for (j = 0; j < va.m[i].aE.n; j++) {V3zero(va.m[i].aE.m[j])}
		}
	};
	void reset_local_F() {
		int i, j;
		for (i = 0; i < va.n; i++) {
			V3zero(va.m[i].F_recp)
			for (j = 0; j < va.m[i].aF.n; j++) {V3zero(va.m[i].aF.m[j])}
		}
	};
	void dump_EF(bool bF = true, bool bCleanBuff = true, bool bE = true) {
		int i, j;
		double v;
		for (i = 0; i < va.n; i++) {
			if (bE) {
				V3plusV3((*(va.m[i].E)), va.m[i].E_recp, (*(va.m[i].E)))
				for (j = 0; j < va.m[i].aE.n; j++) {
					V3plusV3((*(va.m[i].E)), va.m[i].aE.m[j], (*(va.m[i].E)))
				}
				if (bCleanBuff) {
					V3zero(va.m[i].E_recp);
					for (j = 0; j < va.m[i].aE.n; j++) {
						V3zero(va.m[i].aE.m[j])
					}
				}
			}
			if (bF) {
				//V3plusV3((*(va.m[i].F)), va.m[i].F_recp, (*(va.m[i].F)))
				v = va.m[i].F_recp.v[0];
				for (j = 0; j < va.m[i].aF.n; j++) v += va.m[i].aF.m[j].v[0];
				va.m[i].F->v[0] += v * fUnit_estat_mass;

				v = va.m[i].F_recp.v[1];
				for (j = 0; j < va.m[i].aF.n; j++) v += va.m[i].aF.m[j].v[1];
				va.m[i].F->v[1] += v * fUnit_estat_mass;

				v = va.m[i].F_recp.v[2];
				for (j = 0; j < va.m[i].aF.n; j++) v += va.m[i].aF.m[j].v[2];
				va.m[i].F->v[2] += v * fUnit_estat_mass;

				if (bCleanBuff) {
					V3zero(va.m[i].F_recp)
					for (j = 0; j < va.m[i].aF.n; j++) {V3zero(va.m[i].aF.m[j])}
				}
			}
		}
	};

public:
	int nCluster; // nClusters = va.n. This variable is used for multi-thread funciton MacroMoleculeOperate....

	void set_atoms(int natoms) {
		va.set_array(natoms); nCluster = natoms;
	};
};

} // end of namespace _evatom_

namespace _EwaldSum_recp_ {
using namespace _evatom_;

// follow calculation of the reciprocal term of EwaldSum are only valid for cubic cell
double EwaldSum_Recip(_VE_CELL<_VQatom> &vc);

} // end of namespace _EwaldSum_recp_

//#include "EwaldSumVar.h"

