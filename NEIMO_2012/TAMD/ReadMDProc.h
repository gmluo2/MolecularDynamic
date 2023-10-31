bool search(ifstream &in, char *title, char *buffer);
bool search(ifstream &in, char *title, int indx, char *buffer);
bool read_item(char *fname, char *title, char *item);
bool ReadMoleculeStruct(ifstream &in, MMOLECULE &mm);
bool ReadMoleculeStructInfo(ifstream &in, MMOLECULE &mm); // do not rotate the cluster actually, but reset the torsion angle value
bool ReadSolidMoleculeStruct(ifstream &in, SMOLECULE &m);
bool ReadPointMoleculeStruct(ifstream &in, PMOLECULE &m);
bool ConstructMoleculeStruct(MMOLECULE *m, char *fname);

bool search(ifstream &in, char *title, char *buffer, char* stop_line);
bool search(ifstream &in, char *title, int indx, char *buffer, char *stop_line);
