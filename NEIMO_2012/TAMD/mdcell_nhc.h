template <class MM, class SM, class PM, class MD_CELL> bool ReadMDCellNHC(char *fname, MD_CELL& mdcell, float dt) {
	char buffer[256] = "\0", buff[250] = "\0", msg[256] = "\0";
	ifstream in;
	bool bExplicitNHC = true;
	if (!read_msg(fname, buffer, "[NOSE-HOOVER-CHAIN]", true)) {
		show_log("NHC.INI is not given, or [NOSE-HOOVER-CHAIN] is not defined in NHC.INI", true);
		show_log("only one Nose-Hoover Chain is used.", true);
		bExplicitNHC = false;
		//return false;
	}

	int iflag_MM = 0, iflag_SM = 0, iflag_PM = 0;
	int nNHC_MM = 0, nNHC_SM = 0, nNHC_PM = 0, iflag_Independent_Frame = 0;
	int iNHC = 0, n = 0, i, mType = 0, imol = 0;
	float tao = 0, max_ksi = 0, MAX_KSI = 0.8/dt;
	float ftao = 0, fmax_ksi = 0;

	int nNHC = 0;
	if (bExplicitNHC && sscanf(buffer, "%d", &nNHC) != 1) {
		sprintf(msg, "can not get the number of Nose-Hoover Chains from %s in file %s", buffer, fname); show_msg(msg);
		return false;
	}
	if (nNHC <= 0) {
		nNHC = 0;
		mdcell.set_hdyn(nNHC + 1); // additional one is for default, defined in MD.INI
		mdcell.hdyn.m[0].hID = -1; // the last one is the default, hID = -1
		goto _check_default_NHC;;
		//return true;
	}

	//if (!bExplicitNHC) goto _check_default_NHC;

	iflag_MM = 0; iflag_SM = 0; iflag_PM = 0;
	if (!read_msg(fname, buffer, "[MOL-INDEPENDENT-NHC]", true)) return false;
	if (sscanf(buffer, "%d %d %d", &iflag_MM, &iflag_SM, &iflag_PM) != 3) {
		sprintf(msg, "molecular-independent Nose-Hoover Chains is not given properly [see: %s] in file %s", buffer, fname); show_msg(msg);
		return false;
	}
	nNHC_MM = 0; nNHC_SM = 0; nNHC_PM = 0; iflag_Independent_Frame = 0;
	if (!read_msg(fname, buffer, "[MM-NHC]", true)) return false;
	if (sscanf(buffer, "%d %d", &nNHC_MM, &iflag_Independent_Frame) != 2) {
		sprintf(msg, "number of Nose-Hoover Chains for MM & Frame-Independent are not given properly [see: %s] in file %s", buffer, fname); show_msg(msg);
		return false;
	}
	if (!read_msg(fname, buffer, "[SM-NHC]", true)) return false;
	if (sscanf(buffer, "%d", &nNHC_SM) != 1) {
		sprintf(msg, "number of Nose-Hoover Chains for SM is not given properly [see: %s] in file %s", buffer, fname); show_msg(msg);
		return false;
	}
	if (!read_msg(fname, buffer, "[PM-NHC]", true)) return false;
	if (sscanf(buffer, "%d", &nNHC_PM) != 1) {
		sprintf(msg, "number of Nose-Hoover Chains for PM is not given properly [see: %s] in file %s", buffer, fname); show_msg(msg);
		return false;
	}

	nNHC = 0;
	if (iflag_MM > 0) nNHC += mdcell.mm.n * nNHC_MM;
	else nNHC += nNHC_MM;
	if (iflag_SM > 0) nNHC += mdcell.sm.n * nNHC_SM;
	else nNHC += nNHC_SM;
	if (iflag_PM > 0) nNHC += mdcell.pm.n * nNHC_PM;
	else nNHC += nNHC_PM;

	mdcell.set_hdyn(nNHC + 1); // additional one is for default, defined in MD.INI
	mdcell.hdyn.m[nNHC].hID = -1; // the last one is the default, hID = -1

	iNHC = 0; n = 0; mType = 0; imol = 0;

	if (iflag_MM > 0) {
		for (imol = 0; imol < mdcell.lfmd.n; imol++) {
			mdcell.lfmd.m[imol].set_HDyn(NULL, NULL, mdcell.hdyn.m + iNHC);
			iNHC += 1;
		}
	}
	else iNHC += nNHC_MM;
	if (iflag_SM > 0) {
		for (imol = 0; imol < mdcell.lfsmd.n; imol++) {
			mdcell.lfsmd.m[imol].set_HDyn(mdcell.hdyn.m + iNHC);
			iNHC += 1;
		}
	}
	else iNHC += nNHC_SM;
	if (iflag_PM > 0) {
		for (imol = 0; imol < mdcell.lfpmd.n; imol++) {
			mdcell.lfpmd.m[imol].set_HDyn(mdcell.hdyn.m + iNHC);
			iNHC += 1;
		}
	}
	else iNHC += nNHC_PM;

	in.open(fname);
	if (!in.is_open()) return false; // can not be wrong

	iNHC = 0; // the last one is the default
	for (i = 0; i < nNHC_MM; i++) {
		sprintf(buff, "[MMNHC%d]", i);
		if (!search(in, buff, buffer)) {
			sprintf(msg, "%s is not defined in %s", buff, fname); show_msg(msg);
			in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d", &mType) != 1) {
			sprintf(msg, "can not get molecule ID after %s in file %s", buff, fname); show_msg(msg);
			in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f %f", &tao, &max_ksi) != 2) {
			sprintf(msg, "can not get NHC infor. [tao, max_ksi] after %s in file %s", buff, fname); show_msg(msg);
			in.close(); return false;
		}
		if (max_ksi > MAX_KSI) max_ksi = MAX_KSI;

		if (iflag_MM > 0) { // each molecule has its own NHC
			for (imol = 0; imol < mdcell.mm.n; imol++) {
				if (mdcell.mm.m[imol].mType == mType) {
					mdcell.lfmd.m[imol].iHDyn->set_vars(tao, dt, max_ksi);
				}
			}
		}
		else {
			if (iNHC + i >= mdcell.hdyn.n) continue;
			mdcell.hdyn.m[iNHC + i].set_vars(tao, dt, max_ksi);
			for (imol = 0; imol < mdcell.mm.n; imol++) {
				if (mdcell.mm.m[imol].mType == mType) {
					mdcell.lfmd.m[imol].set_HDyn(NULL, NULL, mdcell.hdyn.m + iNHC + i);
				}
			}
		}
	}
	if (iflag_MM > 0) iNHC += mdcell.mm.n;
	else iNHC += nNHC_MM;

	for (i = 0; i < nNHC_SM; i++) {
		sprintf(buff, "[SMNHC%d]", i);
		if (!search(in, buff, buffer)) {
			sprintf(msg, "%s is not defined in %s", buff, fname); show_msg(msg);
			in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d", &mType) != 1) {
			sprintf(msg, "can not get molecule ID after %s in file %s", buff, fname); show_msg(msg);
			in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f %f", &tao, &max_ksi) != 2) {
			sprintf(msg, "can not get NHC infor. [tao, max_ksi] after %s in file %s", buff, fname); show_msg(msg);
			in.close(); return false;
		}
		if (max_ksi > MAX_KSI) max_ksi = MAX_KSI;

		if (iflag_SM > 0) { // each molecule has its own NHC
			for (imol = 0; imol < mdcell.sm.n; imol++) {
				if (mdcell.sm.m[imol].mType == mType) {
					mdcell.lfsmd.m[imol].hdyn->set_vars(tao, dt, max_ksi);
				}
			}
		}
		else {
			if (iNHC + i >= mdcell.hdyn.n) continue;
			mdcell.hdyn.m[iNHC + i].set_vars(tao, dt, max_ksi);
			for (imol = 0; imol < mdcell.sm.n; imol++) {
				if (mdcell.sm.m[imol].mType == mType) {
					mdcell.lfsmd.m[imol].set_HDyn(mdcell.hdyn.m + iNHC + i);
				}
			}
		}
	}

	for (i = 0; i < nNHC_PM; i++) {
		sprintf(buff, "[PMNHC%d]", i);
		if (!search(in, buff, buffer)) {
			sprintf(msg, "%s is not defined in %s", buff, fname); show_msg(msg);
			in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%d", &mType) != 1) {
			sprintf(msg, "can not get molecule ID after %s in file %s", buff, fname); show_msg(msg);
			in.close(); return false;
		}
		in.getline(buffer, 250);
		if (sscanf(buffer, "%f %f", &tao, &max_ksi) != 2) {
			sprintf(msg, "can not get NHC infor. [tao, max_ksi] after %s in file %s", buff, fname); show_msg(msg);
			in.close(); return false;
		}
		if (max_ksi > MAX_KSI) max_ksi = MAX_KSI;

		if (iflag_PM > 0) { // each molecule has its own NHC
			for (imol = 0; imol < mdcell.pm.n; imol++) {
				if (mdcell.pm.m[imol].mType == mType) {
					mdcell.lfpmd.m[imol].hdyn->set_vars(tao, dt, max_ksi);
				}
			}
		}
		else {
			if (iNHC + i >= mdcell.hdyn.n) continue;
			mdcell.hdyn.m[iNHC + i].set_vars(tao, dt, max_ksi);
			for (imol = 0; imol < mdcell.pm.n; imol++) {
				if (mdcell.pm.m[imol].mType == mType) {
					mdcell.lfpmd.m[imol].set_HDyn(mdcell.hdyn.m + iNHC + i);
				}
			}
		}
	}
	if (iflag_PM > 0) iNHC += mdcell.pm.n;
	else iNHC += nNHC_PM;

	in.close();

_check_default_NHC:
	// there might be other molecules have no NHC, connect to the default one
	int iDefaultHDyn = mdcell.hdyn.n - 1;
	for (imol = 0; imol < mdcell.lfmd.n; imol++) {
		if (mdcell.lfmd.m[imol].iHDyn == NULL) {
			mdcell.lfmd.m[imol].set_HDyn(NULL, NULL, mdcell.hdyn.m + iDefaultHDyn);
		}
	}
	for (imol = 0; imol < mdcell.lfsmd.n; imol++) {
		if (mdcell.lfsmd.m[imol].hdyn == NULL) {
			mdcell.lfsmd.m[imol].set_HDyn(mdcell.hdyn.m + iDefaultHDyn);
		}
	}
	for (imol = 0; imol < mdcell.lfpmd.n; imol++) {
		if (mdcell.lfpmd.m[imol].hdyn == NULL) {
			mdcell.lfpmd.m[imol].set_HDyn(mdcell.hdyn.m + iDefaultHDyn);
		}
	}

	return true;
}
