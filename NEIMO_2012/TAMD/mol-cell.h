
template <class MM, class SM, class PM> class CELL_FILE {
public:
	ARRAY<MM> mm;
	ARRAY<CHAIN<VECTOR3>*> mm_pos;

	ARRAY<SM> sm;
	ARRAY<CHAIN<VECTOR3>*> sm_pos;

	ARRAY<PM> pm;
	ARRAY<CHAIN<VECTOR3>*> pm_pos;

	float xyz[3]; // dimension of the cell

	CELL_FILE() {xyz[0] = 0; xyz[1] = 0; xyz[2] = 0;};
	void reset() {
		int i = 0;
		for (i = 0; i < mm.n; i++) mm.m[i].reset();
		mm.release();
		for (i = 0; i < sm.n; i++) sm.m[i].reset();
		sm.release();
		for (i = 0; i < pm.n; i++) pm.m[i].reset();
		pm.release();
		for (i = 0; i < mm_pos.n; i++) release_chain<VECTOR3>(&(mm_pos.m[i]), true);
		mm_pos.release();
		for (i = 0; i < sm_pos.n; i++) release_chain<VECTOR3>(&(sm_pos.m[i]), true);
		sm_pos.release();
		for (i = 0; i < pm_pos.n; i++) release_chain<VECTOR3>(&(pm_pos.m[i]), true);
		pm_pos.release();
	};
	void set_molecules(int nmm, int nsm, int npm) {
		reset();
		int i = 0;
		if (nmm > 0) {
			mm.SetArray(nmm); mm_pos.SetArray(nmm);
			for (i = 0; i < nmm; i++) mm_pos.m[i] = NULL;
		}
		if (nsm > 0) {
			sm.SetArray(nsm); sm_pos.SetArray(nsm);
			for (i = 0; i < nsm; i++) sm_pos.m[i] = NULL;
		}
		if (npm > 0) {
			pm.SetArray(npm); pm_pos.SetArray(npm);
			for (i = 0; i < npm; i++) pm_pos.m[i] = NULL;
		}
	};
	~CELL_FILE() {reset();};
};


template <class MM, class SM, class PM, class MD_CELL> bool ConstructMDCell(CELL_FILE<MM, SM, PM> &cf, MD_CELL &mdcell) {
	int nMM = 0, nSM = 0, nPM = 0, n = 0, ncur = 0;
	for (n = 0; n < cf.mm.n; n++) {
		nMM += number<VECTOR3>(cf.mm_pos.m[n]);
	}
	for (n = 0; n < cf.sm.n; n++) {
		nSM += number<VECTOR3>(cf.sm_pos.m[n]);
	}
	for (n = 0; n < cf.pm.n; n++) {
		nPM += number<VECTOR3>(cf.pm_pos.m[n]);
	}
	mdcell.release();
	mdcell.h[0] = cf.xyz[0] / 2;
	mdcell.h[1] = cf.xyz[1] / 2;
	mdcell.h[2] = cf.xyz[2] / 2;
	mdcell.SetMol(nMM, nSM, nPM);

	ncur = 0;
	int nsame = 0;
	MM *pmm = NULL;
	SM *psm = NULL;
	PM *ppm = NULL;
	CHAIN<VECTOR3> *pos = NULL;
	for (n = 0; n < cf.mm.n; n++) {
		pmm = cf.mm.m + n;
		pos = cf.mm_pos.m[n];
		while (pos != NULL) {
			if (!cp(mdcell.mm.m + ncur, pmm)) {
				sprintf(errmsg, "failure to copy macro-molecule to cell # %d", ncur);
				show_msg(errmsg); mdcell.release(); return false;
			}
			mdcell.mm.m[ncur].shiftMol(*(pos->p));
			ncur++;
			pos = pos->next;
		}
	}
	ncur = 0;
	for (n = 0; n < cf.sm.n; n++) {
		psm = cf.sm.m + n;
		pos = cf.sm_pos.m[n];
		while (pos != NULL) {
			if (!cp(mdcell.sm.m + ncur, psm)) {
				sprintf(errmsg, "failure to copy simple-molecule to cell # %d", ncur);
				show_msg(errmsg); mdcell.release(); return false;
			}
			mdcell.sm.m[ncur].shiftMol(*(pos->p));
			ncur++;
			pos = pos->next;
		}
	}
	ncur = 0;
	for (n = 0; n < cf.pm.n; n++) {
		ppm = cf.pm.m + n;
		pos = cf.pm_pos.m[n];
		while (pos != NULL) {
			if (!cp(mdcell.pm.m + ncur, ppm)) {
				sprintf(errmsg, "failure to copy point-molecule to cell # %d", ncur);
				show_msg(errmsg); mdcell.release(); return false;
			}
			mdcell.pm.m[ncur].shiftMol(*(pos->p));
			ncur++;
			pos = pos->next;
		}
	}
	return true;
};

template <class MM> MM* attach_mol_list(LIST<MM> **mlist, MM* m) { //, int iMol
	char errmsg[2560] = "\0";
	MM *m2 = new MM;
	if (!cp(m2, m)) { // -- default: setup tree, molecule will be copied, the tree will be setup
		sprintf(errmsg, "failure to copy molecule"); show_msg(errmsg); strcpy(errmsg, "\0");
		return NULL;
	}
	if (*mlist == NULL) {*mlist = new LIST<MM>; (*mlist)->p = m2;}
	else (*mlist)->add(m2, true);
	return m2;
};

template <class MM> bool ReadMultiMolecules2List(ifstream &in, char *item, LIST<MM> **mlist, int& nMol, void *construct_mol) { //int &ncur_mol, 
	typedef bool (*MCONSTRUCT)(MM*, char*);
	MM mm, *m = NULL;
	char buffer[256] = "\0", mol_template_fname[250] = "\0";
	int iMol = 0; nMol = 0; 
	float dx = 0, dy = 0, dz = 0;
	VECTOR3 *ds = NULL;

	//strcpy(item, "[MACRO-MOLECULE]");
	if (!search(in, item, buffer)) {
		sprintf(errmsg, "Can not find %s", item);
		return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%s %d", mol_template_fname, &nMol) != 2) {
		sprintf(errmsg, "can not get molecule template filename and the number of the molecule from %s", buffer);
		return false;
	}
	if (nMol <= 0) return true;

	iMol = 0;
	ds = new VECTOR3[nMol];
	if (search(in, "[SHIFT]", buffer)) {
		while(iMol < nMol) {
			if (in.eof()) {
				sprintf(errmsg, "get %d of molecule-relative-position only, less than %d", iMol, nMol);
				show_log(errmsg, true);
				//delete[] ds; return false;
				break;
			}
			in.getline(buffer, 250);
			if (sscanf(buffer, "%f %f %f", &dx, &dy, &dz) != 3) {
				sprintf(errmsg, "can not get molecule-relative-position from : %s", buffer);
				show_log(errmsg, true);
				//delete[] ds; return false;
				break;
			}
			ds[iMol].v[0] = dx; ds[iMol].v[1] = dy; ds[iMol].v[2] = dz;
			iMol++;
		}
	}
	else {
		//sprintf(errmsg, "Can not fine [SHIFT]"); delete[] ds; return false;
	}
	//if (!ConstructMoleculeStruct(&mm, mol_template_fname)) {
	if (!((MCONSTRUCT)construct_mol)(&mm, mol_template_fname)) {
		sprintf(buffer, "Can not construct molecule with parameter file %s", mol_template_fname);
		strcat(errmsg, buffer);
		delete[] ds; return false;
	}
	iMol = 0;
	while (iMol < nMol) {
		//m = attach_mol(&mm, ncur_mol + iMol);
		m = attach_mol_list<MM>(mlist, &mm);
		if (m == NULL) {
			sprintf(errmsg, "failure to create a new molecule");
			//ncur_mol += iMol;
			delete[] ds; return false;
		}
		m->shiftMol(ds[iMol]);
		iMol++; 
	}
	delete[] ds; 
	//ncur_mol += nMol;
	return true;
};

template <class MM> bool ReadMultiMolecules2Array(ifstream &in, char *item, ARRAY<MM> &marray, void* construct_mol()) { //int &ncur_mol, 
	typedef bool (*MCONSTRUCT)(MM*, char*);
	MM mm, *m = NULL;
	char buffer[256] = "\0", mol_template_fname[250] = "\0";
	int iMol = 0, nMol = 0; 
	float dx = 0, dy = 0, dz = 0;
	VECTOR3 *ds = NULL;

	//strcpy(item, "[MACRO-MOLECULE]");
	if (!search(in, item, buffer)) {
		sprintf(errmsg, "Can not find %s", item);
		return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%s %d", mol_template_fname, &nMol) != 2) {
		sprintf(errmsg, "can not get molecule template filename and the number of the molecule from %s", buffer);
		return false;
	}
	if (nMol <= 0) return true;

	ds = new VECTOR3[nMol];
	if (search(in, "[SHIFT]", buffer)) {
		iMol = 0;
		while(iMol < nMol) {
			if (in.eof()) {
				sprintf(errmsg, "get %d of molecule-relative-position only, less than %d", iMol, nMol);
				show_log(errmsg, true);
				//delete[] ds; return false;
				break;
			}
			in.getline(buffer, 250);
			if (sscanf(buffer, "%f %f %f", &dx, &dy, &dz) != 3) {
				sprintf(errmsg, "can not get molecule-relative-position from : %s", buffer);
				show_log(errmsg, true);
				//delete[] ds; return false;
				break;
			}
			ds[iMol].v[0] = dx; ds[iMol].v[1] = dy; ds[iMol].v[2] = dz;
			iMol++;
		}
	}
	else {
		//sprintf(errmsg, "Can not fine [SHIFT]"); delete[] ds; return false;
	}
	//if (!ConstructMoleculeStruct(&mm, mol_template_fname)) {
	if (!((MCONSTRUCT)construct_mol)(&mm, mol_template_fname)) {
		sprintf(buffer, "Can not construct molecule with parameter file %s", mol_template_fname);
		strcat(errmsg, buffer);
		delete[] ds; return false;
	}
	marray.SetArray(nMol);
	iMol = 0;
	for (iMol = 0; iMol < nMol; iMol++) {
		cp(marray.m + iMol, &mm);
		marray.m[iMol].shiftMol(ds[iMol]);
	}
	delete[] ds; 
	return true;
};

template <class MM> bool ReadDefinedMultiMolecules(ifstream &in, char *item, void *construct_mol, MM& mm, CHAIN<VECTOR3> **pos, int &nMol, bool &bShifted) {
	typedef bool (*MCONSTRUCT)(MM*, char*);
	mm.reset();
	if (*pos != NULL) release_chain<VECTOR3>(pos, true);
	nMol = 0;

	char buffer[256] = "\0", mol_template_fname[250] = "\0";
	int iMol = 0, mType = 0, n;
	float dx = 0, dy = 0, dz = 0;
	VECTOR3 *ds = NULL;

	bShifted = true;

	//strcpy(item, "[MACRO-MOLECULE]");
	if (!search(in, item, buffer)) {
		sprintf(errmsg, "Can not find %s", item);
		return false;
	}
	in.getline(buffer, 250);
	if ((n = sscanf(buffer, "%s %d %d", mol_template_fname, &nMol, &mType)) < 2) {
		sprintf(errmsg, "can not get molecule template filename and the number of the molecule from %s", buffer);
		return false;
	}
	if (n < 3) mType = 0; 
	if (nMol <= 0) {bShifted = false; return true;}

	file_pos fpos1 = current_fpos(in);
	ds = NULL;
	if (search(in, "[SHIFT]", buffer)) {
		iMol = 0;
		while(iMol < nMol) {
			ds = new VECTOR3;
			if (iMol == 0) {*pos = new CHAIN<VECTOR3>; (*pos)->p = ds;}
			else (*pos)->attach2tail(ds);

			if (in.eof()) {
				sprintf(errmsg, "get %d of molecule-relative-position only, less than %d", iMol, nMol);
				show_log(errmsg, true);
				bShifted = false;
				//release_chain<VECTOR3>(pos, true); ds = NULL;
				//return false;
				iMol++; continue;
			}
			in.getline(buffer, 250);
			if (sscanf(buffer, "%f %f %f", &dx, &dy, &dz) != 3) {
				sprintf(errmsg, "can not get molecule-relative-position from : %s", buffer);
				show_log(errmsg, true); bShifted = false;
				//release_chain<VECTOR3>(pos, true); ds = NULL;
				//return false;
				iMol++; continue;
			}
			else {ds->v[0] = dx; ds->v[1] = dy; ds->v[2] = dz;}
			iMol++;
		}
	}
	else {
		//sprintf(errmsg, "Can not fine [SHIFT]"); ds = NULL; return false;
		iMol = 0; bShifted = false;
		while(iMol < nMol) {
			ds = new VECTOR3;
			if (iMol == 0) {*pos = new CHAIN<VECTOR3>; (*pos)->p = ds;}
			else (*pos)->attach2tail(ds);
			iMol++;
		}
	}
	set_fpos(in, fpos1);
	//if (!ConstructMoleculeStruct(&mm, mol_template_fname)) {
	if (!((MCONSTRUCT)construct_mol)(&mm, mol_template_fname)) {
		sprintf(buffer, "Can not construct molecule with parameter file %s", mol_template_fname);
		strcat(errmsg, buffer);
		release_chain<VECTOR3>(pos, true); ds = NULL;
		mm.reset();
		return false;
	}
	mm.mType = mType;
	ds = NULL;
	return true;
};

template<class MM> bool ReadMultiMol2List_FromFile(LIST<MM> **mlist, void* construct_mol, char *fname, char *item) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(errmsg, "can not open file %s", fname); return false;
	}
	int nMolTypes = 0, i = 0, nmol = 0, nm = 0;
	char buffer[256] = "\0";
	if (!search(in, "[MOLECULE_TYPES]", buffer)) {
		sprintf(errmsg, "can not find [MOLECULE_TYPES] in %s", fname);
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nMolTypes) != 1) {
		sprintf(errmsg, "can not get number of molecule types from %s", buffer);
		in.close(); return false;
	}

	if (*mlist != NULL) release_list<MM>(mlist, true);
	i = 0; nmol = 0;
	while (i < nMolTypes) {
		if (!ReadMultiMolecules2List<MM>(in, item, mlist, nm, construct_mol)) {in.close(); return false;}
		nmol += nm;
		i++;
	}
	in.close(); return true;
};
/*
template<class MM, class MD_CELL> bool ReadMultiMol2Array(ARRAY<MM> &marray, (bool*)construct_mol(MM*, char*), char *fname) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) {
		sprintf(errmsg, "can not open file %s", fname); return false;
	}
	int nMolTypes = 0, i = 0, nmol = 0, nm = 0;
	char buffer[256] = "\0", item[100] = "\0";
	strcpy(item, "[MOLECULE_TYPES]");
	if (!search(in, item, buffer)) {
		sprintf(errmsg, "can not find %s in %s", item, fname);
		in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &nMolTypes) != 1) {
		sprintf(errmsg, "can not get number of molecule types from %s", buffer);
		in.close(); return false;
	}

	if (*mlist != NULL) release_list<MM>(mlist, true);
	i = 0; nmol = 0;
	while (i < nMolTypes) {
		if (!ReadMultiMolecules<MM>(in, mlist, nm, construct_mol)) {in.close(); return false;}
		nmol += nm;
		i++;
	}
	in.close(); return true;
}
*/

bool ReadDefinedCellSize(char *fname, float &x, float &y, float &z);


bool GenSolventMolPos(VECTOR3& u, int nx, int ny, int nz, CHAIN<VECTOR3> **pos_chain, int& nSolv);

template <class MM, class SM, class PM> bool ReadDefinedCell(char *fname, CELL_FILE<MM, SM, PM> &cf, void* construct_mm, void* construct_sm, void* construct_pm) {
	typedef bool (*SMCONSTRUCT)(SM*, char*);
	typedef bool (*MMCONSTRUCT)(MM*, char*);
	typedef bool (*PMCONSTRUCT)(PM*, char*);

	char title[50] = "\0", buffer[256] = "\0", buff[50] = "\0";
	ifstream check1, check2;
	cf.reset();

	check1.open(fname);
	if (!check1.is_open()) {sprintf(errmsg, "can not open file %s", fname); show_msg(errmsg); return false;}
	
	int nPM = 0, nSM = 0, nMM = 0;
	int solvent_type = -1;
	int n = 0, mol_type = 0;
	char mol_template[250] = "\0";

	int nSolvent = 0;
	float ux, uy, uz, radius = 1;
	VECTOR3 u;

	strcpy(title, "[SOLVENT]");
	if (search(check1, title, buffer)) {
		check1.getline(buffer, 250);
		if (sscanf(buffer, "%d %s %d", &mol_type, mol_template, &n) == 3) {
			if (n <= 0) {
				sprintf(errmsg, "molecule number can not < 0 : %s", buffer); show_msg(errmsg, true);
				check1.close(); return false;
			}
			if (strcmp(mol_template, "NULL") == 0 || strcmp(mol_template, "null") == 0) {
				solvent_type = -1; // unknown
			}
			else solvent_type = mol_type; // simple molecule

			switch (solvent_type) {
			case 0:
				nSM++;
				break;
			case 1:
				nMM++;
				break;
			case 2:
				nPM++;
				break;
			case -1:
				break;
			default: // unknown
				sprintf(errmsg, "UNKNOWN SOLVENT TYPE %d or Solvent is not defined", solvent_type);
				show_msg(errmsg); strcpy(errmsg, "\0");
			}
		}
		else {
			sprintf(errmsg, "UNKNOWN SOLVENT DEFINITION [type, mol, number], see : %s", buffer);
			show_msg(errmsg); strcpy(errmsg, "\0"); check1.close(); return false;
		}
	}

	char *bf = buffer;
	NextSection(&bf, " \t,;"); NextSection(&bf, " \t,;");
	int nx = 0, ny = 0, nz = 0;
	if (sscanf(bf, "%d %d %d", &nx, &ny, &nz) != 3) {ny = nx; nz = nx;}
	check1.getline(buffer, 250);
	if (sscanf(buffer, "%f %f %f", &ux, &uy, &uz) != 3) {
		sprintf(errmsg, "Can not get [a, b, c] of cubic cell from %s", buffer); show_msg(errmsg, true); 
		check1.close(); cf.reset(); return false;
	}
	u.v[0] = ux; u.v[1] = uy; u.v[2] = uz; 
	cf.xyz[0] = float(nx) * ux;
	cf.xyz[1] = float(ny) * uy;
	cf.xyz[2] = float(nz) * uz;

	int n1 = 0, n2 = 0, n3 = 0, nitem;
	strcpy(title, "[MOLECULE_TYPES]");
	if (search(check1, title, buffer)) {
		check1.getline(buffer, 250);
		if ((nitem = sscanf(buffer, "%d %d %d", &n1, &n2, &n3)) < 2) {
			sprintf(errmsg, "No types of molecule is given");
			show_msg(errmsg); check1.close(); return false;
		}
		if (nitem < 3) n3 = 0;
	}
	nMM += n1; nSM += n2; nPM += n3;
	check1.close();
	cf.set_molecules(nMM, nSM, nPM);

	ifstream in1, in2, in3, insolv;
	ARRAY<int> nm_mm, nm_sm, nm_pm;
	nm_mm.SetArray(nMM);
	nm_sm.SetArray(nSM);
	nm_pm.SetArray(nPM);
	int i = 0;
	bool status = false, bShifted = false;
	CHAIN<VECTOR3> *ch_pos = NULL;
	in1.open(fname);
	for (n = 0; n < n1; n++) {
		nm_mm.m[n] = 0;
		//status = ReadDefinedMacroMolecule(in1, cf.mm.m[n], cf.mm_pos.m + n, nm_mm.m[n]);
		status = ReadDefinedMultiMolecules<MM>(in1, "[MACRO-MOLECULE]", construct_mm, cf.mm.m[n], cf.mm_pos.m + n , nm_mm.m[n], bShifted);
		if (!status) {in1.close(); cf.reset(); return false;}
		if (!bShifted) {
			ch_pos = cf.mm_pos.m[n];
			while (ch_pos != NULL) {
				ch_pos->p->v[0] = cf.xyz[0] * (ranf() - 0.5);
				ch_pos->p->v[1] = cf.xyz[1] * (ranf() - 0.5);
				ch_pos->p->v[2] = cf.xyz[2] * (ranf() - 0.5);
				ch_pos = ch_pos->next;
			}
		}
	}
	in1.close();
	in2.open(fname);
	for (n = 0; n < n2; n++) {
		nm_sm.m[n] = 0;
		//status = ReadDefinedSimpleMolecule(in2, cf.sm.m[n], cf.sm_pos.m + n, nm_sm.m[n]);
		status = ReadDefinedMultiMolecules<SM>(in2, "[SIMPLE-MOLECULE]", construct_sm, cf.sm.m[n], cf.sm_pos.m + n, nm_sm.m[n], bShifted);
		if (!status) {in2.close(); cf.reset(); return false;}
		if (!bShifted) {
			ch_pos = cf.sm_pos.m[n];
			while (ch_pos != NULL) {
				ch_pos->p->v[0] = cf.xyz[0] * (ranf() - 0.5);
				ch_pos->p->v[1] = cf.xyz[1] * (ranf() - 0.5);
				ch_pos->p->v[2] = cf.xyz[2] * (ranf() - 0.5);
				ch_pos = ch_pos->next;
			}
		}
	}
	in2.close();
	in3.open(fname);
	for (n = 0; n < n3; n++) {
		nm_pm.m[n] = 0;
		//status = ReadDefinedPointMolecule(in3, cf.pm.m[n], cf.pm_pos.m + n, nm_pm.m[n]);
		status = ReadDefinedMultiMolecules<PM>(in3, "[POINT-MOLECULE]", construct_pm, cf.pm.m[n], cf.pm_pos.m + n, nm_pm.m[n], bShifted);
		if (!status) {in3.close(); cf.reset(); return false;}
		if (!bShifted) {
			ch_pos = cf.pm_pos.m[n];
			while (ch_pos != NULL) {
				ch_pos->p->v[0] = cf.xyz[0] * (ranf() - 0.5);
				ch_pos->p->v[1] = cf.xyz[1] * (ranf() - 0.5);
				ch_pos->p->v[2] = cf.xyz[2] * (ranf() - 0.5);
				ch_pos = ch_pos->next;
			}
		}
	}
	in3.close();

	
	insolv.open(fname);
	strcpy(title, "[SOLVENT]");
	search(insolv, title, buffer);
	insolv.getline(buffer, 250);
	sscanf(buffer, "%d %s", &mol_type, mol_template);
	insolv.getline(buffer, 250);
	/*
	if (sscanf(buffer, "%f %f %f", &ux, &uy, &uz) != 3) {
		sprintf(errmsg, "Can not get [a, b, c] of cubic cell from %s", buffer); show_msg(errmsg, true); 
		insolv.close(); cf.reset(); return false;
	}
	u.v[0] = ux; u.v[1] = uy; u.v[2] = uz; 
	cf.xyz[0] = float(ndim) * ux;
	cf.xyz[1] = float(ndim) * uy;
	cf.xyz[2] = float(ndim) * uz;
	*/

	if (solvent_type == 0) {
		//if (!ConstructSimpleMolecule(cf.sm.m + nSM-1, mol_template)) {
		if (!((SMCONSTRUCT)construct_sm)(cf.sm.m + nSM-1, mol_template)) {
				sprintf(errmsg, "failure to construct solvent molecule %s", mol_template);
				show_msg(errmsg); 
				insolv.close(); cf.reset(); return false;
		}

		status = GenSolventMolPos(u, nx, ny, nz, &(cf.sm_pos.m[nSM-1]), nSolvent);
		
		if (!status) {
			sprintf(errmsg, "failure to construct positions for solvent molecule %s", mol_template);
			show_msg(errmsg); 
			insolv.close(); cf.reset(); return false;
		}
	}
	else if (solvent_type == 1) {// macro-molecule
		if (!((MMCONSTRUCT)construct_mm)(cf.mm.m + nMM-1, mol_template)) {
				sprintf(errmsg, "failure to construct solvent molecule %s", mol_template);
				show_msg(errmsg); 
				insolv.close(); cf.reset(); return false;
		}

		status = GenSolventMolPos(u, nx, ny, nz, &(cf.mm_pos.m[nMM-1]), nSolvent);
		
		if (!status) {
			sprintf(errmsg, "failure to construct positions for solvent molecule %s", mol_template);
			show_msg(errmsg); 
			insolv.close(); cf.reset(); return false;
		}
	}
	else if (solvent_type == 2) {// point-molecule
		if (!((PMCONSTRUCT)construct_pm)(cf.pm.m + nPM-1, mol_template)) {
				sprintf(errmsg, "failure to construct solvent molecule %s", mol_template);
				show_msg(errmsg); 
				insolv.close(); cf.reset(); return false;
		}

		status = GenSolventMolPos(u, nx, ny, nz, &(cf.pm_pos.m[nPM-1]), nSolvent);
		
		if (!status) {
			sprintf(errmsg, "failure to construct positions for solvent molecule %s", mol_template);
			show_msg(errmsg); 
			insolv.close(); cf.reset(); return false;
		}
	}
	insolv.close();
	return true;
};

template <class MM, class SM, class PM, class MD_CELL> bool ConstructDefinedMDCell(char *fname, MD_CELL& mdcell, void *construct_mm, void *construct_sm, void *construct_pm, void *read_mmstruct, void *read_smstruct, void *read_pmstruct) {
	CELL_FILE<MM, SM, PM> cf;
	if (!ReadDefinedCell<MM, SM, PM>(fname, cf, construct_mm, construct_sm, construct_pm)) {cf.reset(); return false;}
	bool status = ConstructMDCell<MM, SM, PM, MD_CELL>(cf, mdcell);
	cf.reset();
	if (!status) {mdcell.release(); return false;}

	// read configuration for each macro-molecule ?
	ifstream in;
	char buffer[256] = "\0", buff[250] = "\0";
	in.open(fname); // can not be wrong, it was opened
	if (!search(in, "[MOL_STRUCT]", buffer)) {in.close(); return status;}
	
	int nm = 0;

	typedef bool (*READ_MMSTRUCT)(ifstream&, MM&);
	typedef bool (*READ_SMSTRUCT)(ifstream&, SM&);
	typedef bool (*READ_PMSTRUCT)(ifstream&, PM&);

	while (mdcell.mm.n > 0 && !in.eof()) {
		if (!search(in, "MM", buffer)) break;
		if (sscanf(buffer, "%s %d", buff, &nm) != 2) break;
		if (nm >= mdcell.mm.n) break; // strange, undefined macro-molecule
		//if (!ReadMoleculeStruct(in, mdcell.mm.m[nm])) {
		if (!((READ_MMSTRUCT)read_mmstruct)(in, mdcell.mm.m[nm])) {
			// shall we do something ?
		}
		if (nm == mdcell.mm.n - 1) break;
	}
	rewind(in);
	if (!search(in, "[MOL_STRUCT]", buffer)) {in.close(); return status;}
	while (mdcell.sm.n > 0 && !in.eof()) {
		if (!search(in, "SM", buffer)) break;
		if (sscanf(buffer, "%s %d", buff, &nm) != 2) break;
		if (nm >= mdcell.sm.n) break; // strange, undefined macro-molecule
		//if (!ReadSolidMoleculeStruct(in, mdcell.sm.m[nm])) {
		if (!((READ_SMSTRUCT)read_smstruct)(in, mdcell.sm.m[nm])) {
			// shall we do something ?
		}
		if (nm == mdcell.sm.n - 1) break;
	}
	rewind(in);
	if (!search(in, "[MOL_STRUCT]", buffer)) {in.close(); return status;}
	while (mdcell.pm.n > 0 && !in.eof()) {
		if (!search(in, "PM", buffer)) break;
		if (sscanf(buffer, "%s %d", buff, &nm) != 2) break;
		if (nm >= mdcell.pm.n) break; // strange, undefined macro-molecule
		//if (!ReadPointMoleculeStruct(in, mdcell.pm.m[nm])) {
		if (!((READ_PMSTRUCT)read_pmstruct)(in, mdcell.pm.m[nm])) {
			// shall we do something ?
		}
		if (nm == mdcell.pm.n - 1) break;
	}
	in.close();
	return status;
}
