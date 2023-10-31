namespace _atomic_cmm_ {
// assuming vspme.va is same as that of ARRAY<MATOM*> atom, in another word, same as that of AtCELL *cell

	void GPU_LJ_INTERACT_ATOM(JobIndx *job, SR_AtCELL *cell);
	void GPU_LJ_INTERACT_BOUND14(JobIndx *job, SR_AtCELL *cell);
	void GPU_EwdSumReal_INTERACT_ATOM(JobIndx *job, EAtCELL *cell);
	void GPU_EwdSumReal_INTERACT_InCluster(JobIndx *job, EAtCELL *cell);
	void GPU_EwdSumReal_INTERACT_BOUND(JobIndx *job, EAtCELL *cell);

template <class _VATOM, class VSPME> void _reset_global_efield(JobIndx* job, VSPME* vspme) {
	_VATOM *atom = vspme->va.m;
	double *E;
	int ia;
	for (ia = job->n1; ia <= job->n2; ia++) {
		E = atom[ia].E->v;
		memset(E, 0, SIZE_V3);
	}
}

inline void reset_global_efield_q_2d(JobIndx* job, _spme_2d_::VSPME<_VQatom>* vspme) {
	_reset_global_efield<_VQatom, _spme_2d_::VSPME<_VQatom> >(job, vspme);
};

inline void reset_global_efield_mu_2d(JobIndx* job, _spme_2d_::VSPME<_VMUatom>* vspme) {
	_reset_global_efield<_VMUatom, _spme_2d_::VSPME<_VMUatom> >(job, vspme);
};

inline void reset_global_efield_q_3d(JobIndx* job, _spme_::VSPME<_VQatom>* vspme) {
	_reset_global_efield<_VQatom, _spme_::VSPME<_VQatom> >(job, vspme);
};

inline void reset_global_efield_mu_3d(JobIndx* job, _spme_::VSPME<_VMUatom>* vspme) {
	_reset_global_efield<_VMUatom, _spme_::VSPME<_VMUatom> >(job, vspme);
};


template <class _VATOM, class VSPME> void _reset_global_eforce(JobIndx* job, VSPME* vspme) {
	_VATOM *atom = vspme->va.m;
	double *F;
	int ia;
	for (ia = job->n1; ia <= job->n2; ia++) {
		F = atom[ia].F->v;
		memset(F, 0, SIZE_V3);
	}
}


inline void reset_global_force_q_2d(JobIndx* job, _spme_2d_::VSPME<_VQatom>* vspme) {
	_reset_global_eforce<_VQatom, _spme_2d_::VSPME<_VQatom> >(job, vspme);
};

inline void reset_global_force_mu_2d(JobIndx* job, _spme_2d_::VSPME<_VMUatom>* vspme) {
	_reset_global_eforce<_VMUatom, _spme_2d_::VSPME<_VMUatom> >(job, vspme);
};

inline void reset_global_force_q_3d(JobIndx* job, _spme_::VSPME<_VQatom>* vspme) {
	_reset_global_eforce<_VQatom, _spme_::VSPME<_VQatom> >(job, vspme);
};

inline void reset_global_force_mu_3d(JobIndx* job, _spme_::VSPME<_VMUatom>* vspme) {
	_reset_global_eforce<_VMUatom, _spme_::VSPME<_VMUatom> >(job, vspme);
};


template <class _VATOM, class VSPME> void collect_efield_sr(JobIndx* job, EAtCELL *cell, VSPME* vspme) {
	_VATOM *atom = vspme->va.m;
	double *E;
	int ia, i3a;
	for (ia = job->n1; ia <= job->n2; ia++) {
		i3a = x_idx(ia);
		E = atom[ia].E->v;
		E[0] += cell->E.m[i3a];
		E[1] += cell->E.m[i3a + 1];
		E[2] += cell->E.m[i3a + 2];
	}
}

inline void collect_efield_sr_q_2d(JobIndx* job, EAtCELL *cell, _spme_2d_::VSPME<_VQatom>* vspme) {
	collect_efield_sr<_VQatom, _spme_2d_::VSPME<_VQatom> >(job, cell, vspme);
};

inline void collect_efield_sr_mu_2d(JobIndx* job, EAtCELL *cell, _spme_2d_::VSPME<_VMUatom>* vspme) {
	collect_efield_sr<_VMUatom, _spme_2d_::VSPME<_VMUatom> >(job, cell, vspme);
};

inline void collect_efield_sr_q_3d(JobIndx* job, EAtCELL *cell, _spme_::VSPME<_VQatom>* vspme) {
	collect_efield_sr<_VQatom, _spme_::VSPME<_VQatom> >(job, cell, vspme);
};

inline void collect_efield_sr_mu_3d(JobIndx* job, EAtCELL *cell, _spme_::VSPME<_VMUatom>* vspme) {
	collect_efield_sr<_VMUatom, _spme_::VSPME<_VMUatom> >(job, cell, vspme);
};

template <class _VATOM, class VSPME> void collect_efresult_sr(JobIndx* job, EAtCELL *cell, VSPME* vspme) {
	_VATOM *atom = vspme->va.m;
	double *F;
	int ia, i3a;
	for (ia = job->n1; ia <= job->n2; ia++) {
		i3a = x_idx(ia);
		F = atom[ia].F->v;
		/*
		F[0] += cell->Fsr.m[i3a]; // has unit in mass
		F[1] += cell->Fsr.m[i3a + 1];
		F[2] += cell->Fsr.m[i3a + 2];
		*/

		F[0] += cell->Fe.m[i3a] * ::fUnit_estat_mass;
		F[1] += cell->Fe.m[i3a + 1] * ::fUnit_estat_mass;
		F[2] += cell->Fe.m[i3a + 2] * ::fUnit_estat_mass;
	}
}


inline void collect_fresult_sr_q_2d(JobIndx* job, EAtCELL *cell, _spme_2d_::VSPME<_VQatom>* vspme) {
	collect_efresult_sr<_VQatom, _spme_2d_::VSPME<_VQatom> >(job, cell, vspme);
};

inline void collect_fresult_sr_mu_2d(JobIndx* job, EAtCELL *cell, _spme_2d_::VSPME<_VMUatom>* vspme) {
	collect_efresult_sr<_VMUatom, _spme_2d_::VSPME<_VMUatom> >(job, cell, vspme);
};

inline void collect_fresult_sr_q_3d(JobIndx* job, EAtCELL *cell, _spme_::VSPME<_VQatom>* vspme) {
	collect_efresult_sr<_VQatom, _spme_::VSPME<_VQatom> >(job, cell, vspme);
};

inline void collect_fresult_sr_mu_3d(JobIndx* job, EAtCELL *cell, _spme_::VSPME<_VMUatom>* vspme) {
	collect_efresult_sr<_VMUatom, _spme_::VSPME<_VMUatom> >(job, cell, vspme);
};

void collect_Uvirial(JobIndx* job, EAtCELL *cell, SVECTOR<double, 4> *res);

void collect_sr_fresult(JobIndx* job, SR_AtCELL *cell);
void collect_Uvirial_sr(JobIndx* job, SR_AtCELL *cell, SVECTOR<double, 4> *res);

void SPME_EF(EAtCELL &cell, _spme_2d_::VSPME<_VQatom>& vspme, int nThreads, _spme_2d_::SPME_VAR& var, InteractRes& res);
void SPME_EF(EAtCELL &cell, _spme_2d_::VSPME<_VMUatom>& vspme, int nThreads, _spme_2d_::SPME_VAR& var, InteractRes& res, bool bInitCoordination);
bool Polarize_SPME_Interact(MMOL_MD_CELL &mdcell, EAtCELL &cell, _spme_2d_::VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, _spme_2d_::SPME_VAR &var, InteractRes &res);
bool Polarize_SPME_Interact2(MMOL_MD_CELL &mdcell, EAtCELL &cell, _spme_2d_::VSPME<_VMUatom>& vspme, _spme_2d_::VSPME<_VMUatom>& vspme_induced, int nThreads, bool bPolarize, _spme_2d_::SPME_VAR &var, InteractRes &res);

void SPME_EF(EAtCELL &cell, _spme_::VSPME<_VQatom>& vspme, int nThreads, SPME_VAR& var, InteractRes& res);
void SPME_EF(EAtCELL &cell, _spme_::VSPME<_VMUatom>& vspme, int nThreads, SPME_VAR &var, InteractRes &res, bool bInitCoordination);
bool Polarize_SPME_Interact(MMOL_MD_CELL &mdcell, EAtCELL &cell, _spme_::VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);
bool Polarize_SPME_Interact2(MMOL_MD_CELL &mdcell, EAtCELL &cell, _spme_::VSPME<_VMUatom>& vspme, _spme_::VSPME<_VMUatom>& vspme_induced, int nThreads, bool bPolarize, SPME_VAR &var, InteractRes &res);

void LJ_Interact(SR_AtCELL *cell, int nThreads, InteractRes &res);

void LocalF(EAtCELL &cell, _spme_::VSPME<_VQatom> &vspme, int nThreads, SPME_VAR& var, InteractRes& res);
void LocalF(EAtCELL &cell, _spme_::VSPME<_VMUatom> &vspme, int nThreads, SPME_VAR& var, InteractRes& res);

void LocalF(EAtCELL &cell, _spme_2d_::VSPME<_VQatom> &vspme, int nThreads, SPME_VAR& var, InteractRes& res);
void LocalF(EAtCELL &cell, _spme_2d_::VSPME<_VMUatom> &vspme, int nThreads, SPME_VAR& var, InteractRes& res);

// asynchrotron way for LJ interact
class AsynLJVar {
public:
	SR_AtCELL *cell;
	int nThreads;
	InteractRes *res;
};
void collect_LJ_res(SR_AtCELL *cell, int nThreads, InteractRes &res);
void cpuLJ(AsynLJVar *var);
#ifdef __GPU__
void cudaLJ(AsynLJVar *var);
#endif

} // end of namespace


