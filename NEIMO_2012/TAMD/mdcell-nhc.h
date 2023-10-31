template <class MM, class SM, class PM, class MD_CELL> bool ReadMDCellNHC(char *fname, MD_CELL& mdcell, float dt) {
	char buffer[256] = "\0", buff[250] = "\0", msg[256] = "\0";
	if (!read_msg(fname, "[NOSE-HOOVER-CHAIN]", buffer)) return false;
	int nNHC = 0;
	if (sscanf(buffer, "%d", &nNHC) != 1) {
		sprintf(msg, "can not get the number of Nose-Hoover Chains from %s in file %s", buffer, fname); show_msg(msg);
		return false;
	}
	if (nNHC < 0) nNHC = 0;
	mdcell.set_nhc(nNHC + 1); // additional one is for default, defined in MD.INI
	if (nNHC == 0) return true;

	int iNHC = 0, n = 0, i, iflag = 0, mID = 0, imol = 0;
	float tao = 0, max_ksi = 0, MAX_KSI = 0.8/dt;
	ifstream in;
	in.open(fname);
	if (!in.is_open()) return false; // can not be wrong
	if (!search(in, "[MM-NHC]", buffer)) {
		sprintf(msg, "[MM-NHC] is not defined in %s", fname); show_msg(msg);
	}
	else {
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d %d", &n, &iflag) != 2) {
			sprintf(msg, "[MM-NHC] in %s is not defined properly, see %s", fname, buffer); show_msg(msg);
			in.close(); return false;
		}
		for (i = 0; i < n; i++) {
			sprintf(buff, "[MMNHC%d]", i);
			if (!search(in, buff, buffer)) {
				sprintf(msg, "%s is not defined in %s", buff, fname); show_msg(msg);
				in.close(); return false;
			}
			in.getline(buffer, 250);
			if (sscanf(buffer, "%d", mID) != 1) {
				sprintf(msg, "can not get molecule ID after %s in file %s", buff, fname); show_msg(msg);
				in.close(); return false;
			}
			in.getline(buffer, 250);
			if (sscanf(buffer, "%f %f", &tao, &max_ksi) != 2) {
				sprintf(msg, "can not get NHC infor. [tao, max_ksi] after %s in file %s", buff, fname); show_msg(msg);
				in.close(); return false;
			}
			if (max_ksi > MAX_KSI) max_ksi = MAX_KSI;
			mdcell.hdyn.m[iNHC].hID = mID;
			mdcell.hdyn.m[iNHC].set_vars(tao, dt, max_ksi);
			for (imol = 0; imol < mdcell.mm.n; imol++) {
				if (mdcell.mm.m[imol].mID == mID)
			}
		}
	}

	in.close();
	return status;
}
