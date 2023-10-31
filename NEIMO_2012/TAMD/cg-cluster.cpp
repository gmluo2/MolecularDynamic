#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
#include "show.h"
#include "vector.h"
#include "Matrix.h"
#include "nhc.h"
#include "bound.h"

#include "ZMatrix.h"

#include "MM.h"
#include "cg-mm.h"
#include "Mol.h"
#include "MD.h"
#include "read.h"
#include "var.h"
#include "cluster.h"
#include "cg-cluster.h"

namespace _coarse_grain_ {
/*********************************************************************
      functions to construct coarse-grained macromolecules
*********************************************************************/

bool ParseSingleChainPolymer(char *polymer, CG_CLUSTER **cluster, int &ncs, int **nperiod) {
	char buffer[256] = "\0";
	char chead[250] = "\0", ctail[250] = "\0", cperiod[250] = "\0";
	int nblock = 0, n = 0, nbase = 0;

	char *str = polymer;
	while (get_period(&str, buffer) && strlen(buffer) > 0) {
		n++; if (str == NULL || strlen(str) == 0) break;
		get_period_number(&str, nblock);
		if (str == NULL || strlen(str) == 0) break;
	}
	ncs = n; nblock = n;
	if (nblock == 0) return false;
	if (*cluster != NULL) delete[] *cluster;
	*cluster = new CG_CLUSTER[nblock];
	if (*nperiod != NULL) {delete[] *nperiod; *nperiod = NULL;}
	if (nblock > 2) *nperiod = new int[nblock - 2];

	str = polymer;
	get_period(&str, chead);
	if (!construct_defined_cgcluster(chead, 1, *cluster)) { // this is head
		delete[] *cluster; *cluster = NULL; 
		if (*nperiod != NULL) {delete[] *nperiod; *nperiod = NULL;}
		return false;
	}

	n = 0;
	while (n < nblock - 2) {
		if (get_period(&str, cperiod) && get_period_number(&str, (*nperiod)[n])) {
			// here true will not affect the monomer because it is not head or tail
			if (!construct_defined_cgcluster(cperiod, 2, (*cluster) + n + 1)) {
				delete[] *cluster; *cluster = NULL; 
				if (*nperiod != NULL) {delete[] *nperiod; *nperiod = NULL;}
				return false;
			}
			n++;
		}
	}

	if (!get_period(&str, ctail) || !construct_defined_cgcluster(ctail, 1, *cluster + nblock - 1)) { // this is tail
		delete[] *cluster; *cluster = NULL; 
		if (*nperiod != NULL) {delete[] *nperiod; *nperiod = NULL;}
		return false;
	}

	return true;
}

// unit is composed with cluster c in order, start with head and end with tail
// using C-C as hing bound to connect the clusters
void SingleChainPolymer(CG_MMOLECULE &mm, CG_CLUSTER* head, CG_CLUSTER* cluster, int ncs, int *nperiod, CG_CLUSTER* tail) {
	int i = 0, j = 0, k = 0, n = 1;

	int nCluster = 2; // head + tail
	for (n = 0; n < ncs; n++) nCluster += nperiod[n];
	CG_CLUSTER *pc = NULL, *parent = NULL;
	int indx_cluster = 0;
	int iHinge1 = 0, iHinge2 = 0;

	mm.reset();
	mm.nCluster = nCluster;
	mm.cluster = new CG_CLUSTER[nCluster];

	indx_cluster = 0;
	pc = mm.cluster + indx_cluster;

	cp(pc, head); // copy head to macromolecules

	for (i = 0; i < ncs; i++) {
		for (j = 0; j < nperiod[i]; j++) {
			indx_cluster++; // a new monomer in macromolecule
			parent = pc;
			pc = mm.cluster + indx_cluster;
			cp(pc, cluster + i); // copy the cluster to molecule
			// setup relationship -- no parent-child relationship
			//cluster[m[i].nc[0]].set_parent_hinge(m[i].iHinge[0]);
			//setup_parent_child_relationship(cluster + m[i].nc[0]); 
			iHinge1 = (indx_cluster == 1 ? 0 : 1);
			iHinge2 = 0;
			ConnectClusters(parent, iHinge1, pc, iHinge2, SINGLE_BOUND);
		}
	}
	// tail 
	indx_cluster++; parent = pc; pc = mm.cluster + indx_cluster;
	cp(pc, tail); // copy the cluster to molecule
	// connect the cluster to its parent
	iHinge1 = (indx_cluster == 1 ? 0 : 1); iHinge2 = 0;
	ConnectClusters(parent, iHinge1, pc, iHinge2, SINGLE_BOUND);
}

bool ConstructSingleChainPolymer(CG_MMOLECULE &mm, char *polymer) {
	mm.reset();
	CG_CLUSTER *cluster = NULL;
	int ncs = 0, nblock = 0, *nperiod = NULL;
	if (!ParseSingleChainPolymer(polymer, &cluster, ncs, &nperiod)) {
		sprintf(errmsg, "polymer %s can not be recognized");
		if (cluster != NULL) {delete[] cluster; cluster = NULL;}
		if (nperiod != NULL) {delete[] nperiod; nperiod = NULL;}
		return false;
	}
	nblock = ncs - 2; // kick off head and tail
	SingleChainPolymer(mm, cluster, cluster + 1, nblock, nperiod, cluster + ncs - 1);
	if (cluster != NULL) {delete[] cluster; cluster = NULL;}
	if (nperiod != NULL) {delete[] nperiod; nperiod = NULL;}

	return true;
}

bool ConstructSingleChainPolymerFromParamFile(CG_MMOLECULE &mm, char *fname, bool save) {
	char buffer[256] = "\0", mol[250] = "\0";

	if (!read_msg(fname, mol, "[MOLECULE]")) return false;

	bool status = ConstructSingleChainPolymer(mm, mol);
	if (!status) return false;

	//if (!mm.linkPars(atompar_db, errmsg)) {show_infor(errmsg); return false;}

	// set mass center to {0, 0, 0}
	VECTOR<3> cm;
	mm.SetMassCenter(cm.v, true);

	char template_fname[256] = "\0";
	if (save) {
		if (read_msg(fname, template_fname, "[MOL_TEMPLATE_FILENAME]")) {
			{
				ofstream out;
				out.open(template_fname);
				if (out.is_open()) {
					mm.xyz_save(out); out.close();
				}
				else {
					sprintf(buffer, "can not open file %s, failure to save molecule structure", template_fname);
					show_msg(buffer);
				}
			}
		}
	}
	return true;
}

} // end of namespace _coarse_grain_
