namespace _coarse_grain_ {

void gen_add_random_force_cgmm(CG_MMOLECULE &mm, _rand_::RAND_ARRAY<int> &randa);

// interactions calculated with CELL MULTIPOLE METHOD
void FreeCellCMMInteraction_SingleCGMM(CMM_CELL3D<CG_CLUSTER> &cell, CG_MMOLECULE &mm, double &V, bool &OVERLAP, bool hard_sphere, int n1, int n2); // V in kT

void FreeCellCMMInteraction_CGSM(CMM_CELL3D<CG_CLUSTER> &cell, CGSM* sm, int nSM, double &V, bool &OVERLAP); // V in kT

void SolventGaussForce(CG_MMOLECULE *cgmm, int nc1, int nc2);

void scale_force_cgcluster(CG_MMOLECULE *mm, int nc1, int nc2);
void scale_force_cgsm(CGSM *sm, int nSM);

void CG_Brownian_force(CG_MMOLECULE *mm, int nc1, int nc2);

void cgmm_IntraMolInteract(CG_MMOLECULE *mm, int ncf, int nct, double &V);

void CMM_CGImplicitSolvationForce(CMM_CELL3D<CG_CLUSTER> &cell, int n1Cell, int n2Cell);
void CGMM_ImplicitSolvationForce(CMM_CELL3D<CG_CLUSTER> *cell, CG_MMOLECULE *mm, int nc1, int nc2);

} // end of namespace _coarse_grain_ 
