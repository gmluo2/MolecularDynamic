namespace _coarse_grain_ {

bool construct_defined_cgcluster(char *unit, int nb, CG_CLUSTER* p);
bool ConstructSingleChainPolymer(CG_MMOLECULE &mm, char *polymer);
bool ConstructSingleChainPolymerFromParamFile(CG_MMOLECULE &mm, char *fname, bool save = false);

} // end of namespace _coarse_grain_ 
