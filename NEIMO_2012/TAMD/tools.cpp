#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "ranlib.h"

#include "def.h"
#include "show.h"
#include "vector.h"
#include "bound.h"
#include "ZMatrix.h"

#include "MM.h"

extern bool search(ifstream &in, char *title, char *buffer);
extern bool search(ifstream &in, char *title, int indx, char *buffer);
extern bool read_item(char *fname, char *title, char *item);

extern bool read_cluster(char *cname, BASIC_CLUSTER *pc);
extern bool read_monomer(char *monomer, MONOMER *m);

extern LOG mlog;
extern char errmsg[2560];

bool combine_clusters(BASIC_CLUSTER* c1, BASIC_CLUSTER* c2, BASIC_CLUSTER *d) {
	d->reset();
	int nAtoms = c1->nAtoms + c2->nAtoms;
	if (nAtoms <= 0) return true;
	int nHinges = c1->nHinges + c2->nHinges;

	d->nAtoms = nAtoms;
	d->atom = new MATOM[d->nAtoms];
	int i = 0, j = 0, n = 0, nb = 0;
	for (i = 0; i < c1->nAtoms; i++) cp(d->atom + i, c1->atom + i);
	for (i = 0; i < c2->nAtoms; i++) cp(d->atom + c1->nAtoms + i, c2->atom + i);

	int n1atom_connect = -1, n2atom_connect = -1, nb1_connect = -1, nb2_connect = -1;
	// forming the bound 
	for (i = 0; i < c1->nAtoms; i++) {
		for (j = 0; j < d->atom[i].nb; j++) {
			if (c1->atom[i].ba[j] != NULL) {
				n = c1->iAtom((MATOM*)c1->atom[i].ba[j]);
				if (n < 0) { // this bound could be in other cluster
					n1atom_connect = i; // this atom connected to another cluster
					nb1_connect = j;
					//sprintf(errmsg, "ERROR : Program is not consistent. Can not find atom defined in cluster");
					//return false;
				}
				else {
					nb = c1->atom[n].ibound(c1->atom + i);
					if (nb < 0) {
						sprintf(errmsg, "ERROR: Program is not consisten to find out the bound connection in defined cluster");
						return false;
					}
					if (!bound_connect<MATOM>(d->atom + i, j, d->atom + n, nb, c1->atom[i].bv[j], c1->atom[i].btype[j])) return false;
				}
			}
		}
	}
	for (i = 0; i < c2->nAtoms; i++) {
		for (j = 0; j < d->atom[i + c1->nAtoms].nb; j++) {
			if (c2->atom[i].ba[j] != NULL) {
				n = c2->iAtom((MATOM*)c2->atom[i].ba[j]);
				if (n < 0) { // this bound could be in other cluster
					n2atom_connect = i;  // this atom connected to another cluster
					nb2_connect = j;
					//sprintf(errmsg, "ERROR : Program is not consistent. Can not find atom defined in cluster");
					//return false;
				}
				else {
					nb = c2->atom[n].ibound(c2->atom + i);
					if (nb < 0) {
						sprintf(errmsg, "ERROR: Program is not consisten to find out the bound connection in defined cluster");
						return false;
					}
					if (!bound_connect<MATOM>(d->atom + c1->nAtoms + i, j, d->atom + c1->nAtoms + n, nb, c2->atom[i].bv[j], c2->atom[i].btype[j])) return false;
				}
			}
		}
	}
	if (n1atom_connect >= 0 && n2atom_connect >= 0) { // c1 & c2 are connected
		if (!bound_connect<MATOM>(d->atom + n1atom_connect, nb1_connect, d->atom + c1->nAtoms + n2atom_connect, nb2_connect, c2->atom[n2atom_connect].bv[nb2_connect], c2->atom[n2atom_connect].btype[nb2_connect])) return false;
		nHinges -= 2; //each cluster has a hinge used for connection
	}

	d->set_hinges(nHinges);

	n = 0;
	for (i = 0; i < c1->nHinges; i++) {
		if (c1->hinge[i].plink == c2) continue;
		d->hinge[n].vHinge = c1->hinge[i].vHinge;
		if (c1->hinge[i].indx >= 0) d->setup_hinge(n, c1->hinge[i].indx, c1->hinge[i].ibound);
		n++;
		// the connection can not be copied
	}
	for (i = 0; i < c2->nHinges; i++) {
		if (c2->hinge[i].plink == c1) continue;
		d->hinge[n].vHinge = c2->hinge[i].vHinge;
		if (c2->hinge[i].indx >= 0) d->setup_hinge(n, c2->hinge[i].indx + c1->nAtoms, c2->hinge[i].ibound);
		n++;
		// the connection can not be copied
	}
	d->d = c1->d + c2->d;
	d->eNeutral = c1->eNeutral;
	d->deuterized = c1->deuterized & c2->deuterized;

	//short mIndx; // the molecule index in the whole system, where this cluster is
	d->cID = 0;
	d->cgIndx = 0; // index of coarse-grained cluster where this cluster is in
	d->monIndx = 0; // monomer index inside the macro-molecule

	return true;
}

bool construct_tree(BASIC_CLUSTER *c, TREE<MATOM> &atree) {
	atree.release_tree(false);
	if (c->nAtoms == 0) return false;

	MATOM *patom = NULL;
	int natom = 0;
	
	LIST<MATOM> *plist = NULL, *ptmp = NULL;
	LIST< LIST<MATOM> > *plevel = NULL, *last_level = NULL;

	plevel = atree.new_level();
	plist = new LIST<MATOM>; plist->p = c->atom;
	plevel->add(plist, false);

	last_level = plevel; 
	while (last_level != NULL && last_level->p != NULL) {
		plist = last_level->p;
		plevel = NULL;
		while (plist != NULL) {
			patom = plist->p;
			plist = plist->next;
			if (patom == NULL) continue;
			for (natom = 0; natom < patom->nb; natom++) {
				if (patom->ba[natom] == NULL) continue;
				if (in_tree((MATOM*)patom->ba[natom], atree)) continue;
				if (plevel == NULL) {
					plevel = atree.new_level();
					ptmp = new LIST<MATOM>; ptmp->p = (MATOM*)patom->ba[natom];
					plevel->p = ptmp;
				}
				else plevel->p->add((MATOM*)patom->ba[natom], false);
			}
		}
		last_level = plevel;
	}
	return true;
}

bool search_connected_atoms(MATOM *patom, LIST< LIST<MATOM> >* plevel, MATOM **batom, MATOM **aatom, MATOM **tatom, bool global_search = false) {
	LIST< LIST<MATOM> > *parent = plevel->pre, *gparent = NULL, *ggparent = NULL, *child = NULL;
	if (parent != NULL) gparent = parent->pre;
	if (gparent != NULL) ggparent = gparent->pre;
	if (plevel != NULL) child = plevel->next;
	LIST<MATOM> *plist = NULL;
	MATOM *p = NULL;
	bool bound = false;
	int i = 0;

	*batom = NULL; *aatom = NULL; *tatom = NULL;

	if (parent == NULL) return false; // bound atom has to be in parent level, that is how the tree is constructed
	plist = parent->p;
	while (plist != NULL) {
		BOUND(plist->p, patom, i, bound)
		if (bound) {*batom = plist->p; break;}
		plist = plist->next;
	}
	if (*batom == NULL) return false;

	// look for the atom connected to batom for bound-angle, it must be in the levels [parent +/- 1]
	if (gparent != NULL) {
		plist = gparent->p;
		while (plist != NULL) {
			BOUND(plist->p, (*batom), i, bound)
			if (bound) {*aatom = plist->p; break;}
			plist = plist->next;
		}
	}
	if (*aatom == NULL) {
		plist = parent->p;
		while (plist != NULL) {
			if (plist->p != *batom) {
				BOUND(plist->p, (*batom), i, bound)
				if (bound) {*aatom = plist->p; break;}
			}
			plist = plist->next;
		}
	}
	if (*aatom == NULL) {
		plist = plevel->p;
		while (plist != NULL) {
			if (plist->p == patom) return false; // unfortranately, the atom in this level and after patom, is not zmatrixed yet
			if (plist->p != *batom) {
				BOUND(plist->p, (*batom), i, bound)
				if (bound) {*aatom = plist->p; break;}
			}
			plist = plist->next;
		}
	}
	if (*aatom == NULL) return false;

	// look for tatom -- torsion-angle-bounded atom
	// it must be in [ggparent, gparnet, parent, plevel, child]
	if (ggparent != NULL) {
		plist = ggparent->p;
		while (plist != NULL) {
			BOUND(plist->p, (*aatom), i, bound)
			if (bound) {*tatom = plist->p; break;}
			// the atom can not bound to bound atom, too far away
			plist = plist->next;
		}
	}
	if (*tatom == NULL) {
		if (gparent != NULL) {
			plist = gparent->p;
			while (plist != NULL) {
				if (plist->p != *batom && plist->p != *aatom) {
					BOUND(plist->p, (*aatom), i, bound)
					if (bound) {*tatom = plist->p; break;}
					if (!bound) {
						BOUND(plist->p, (*batom), i, bound)
						if (bound) {*tatom = plist->p; break;}
					}
				}
				plist = plist->next;
			}
		}
	}
	if (*tatom == NULL) {
		plist = parent->p;
		while (plist != NULL) {
			if (plist->p != *batom && plist->p != *aatom) {
				BOUND(plist->p, (*aatom), i, bound)
				if (bound) {*tatom = plist->p; break;}
				if (!bound) {
					BOUND(plist->p, (*batom), i, bound)
					if (bound) {*tatom = plist->p; break;}
				}
			}
			plist = plist->next;
		}
	}
	if (*tatom == NULL) {
		plist = plevel->p;
		while (plist != NULL) {
			if (plist->p == patom) {
				if (global_search) {plist = plist->next; continue;}
				else return false; // unfortranately, the atom in this level and after patom, is not zmatrixed yet
			}
			if (plist->p != *batom && plist->p != *aatom) {
				BOUND(plist->p, (*aatom), i, bound)
				if (bound) {*tatom = plist->p; break;}
				if (!bound) {
					BOUND(plist->p, (*batom), i, bound)
					if (bound) {*tatom = plist->p; break;}
				}
			}
			plist = plist->next;
		}
	}
	if (*tatom == NULL && global_search) {
		plist = child->p;
		while (plist != NULL) {
			if (plist->p != *batom && plist->p != *aatom) {
				BOUND(plist->p, (*aatom), i, bound)
				if (bound) {*tatom = plist->p; break;}
				if (!bound) {
					BOUND(plist->p, (*batom), i, bound)
					if (bound) {*tatom = plist->p; break;}
				}
			}
			plist = plist->next;
		}
	}

	if (*tatom == NULL) return false;

	return true;
}

extern "C" bool combine_clusters(char *par_fname) {
	::mlog.init("cb.log");

	char buffer[256] = "\0", msg[256] = "\0", ofname[256] = "\0";
	char c1_name[256] = "\0", c2_name[256] = "\0", c_name[256] = "\0";
	int iHinge1 = 0, iHinge2 = 0;
	float ta = 0;

	if (!read_item(par_fname, "[CONNECTION]", buffer)) return false;
	if (sscanf(buffer, "%s %d %s %d %f", c1_name, &iHinge1, c2_name, &iHinge2, &ta) != 5) {
		sprintf(msg, "connection parameter format : [cluster1] [hinge1] [cluster2] [hinge2] [torsion angle], see %s", buffer);
		show_msg(msg); return false;
	}
	if (!read_item(par_fname, "[OUTPUT]", buffer)) return false;
	if (sscanf(buffer, "%s %s", ofname, c_name) != 2) {
		sprintf(msg, "output parameter format : [output filename] [output cluster name], see %s", buffer);
		show_msg(msg); return false;
	}

	CLUSTER c1, c2, c;
	if (!read_cluster(c1_name, &c1)) return false;
	if (!read_cluster(c2_name, &c2)) return false;
	// connect c2 [child] to c1 [parent] with vcHinge [in c1]
	if (!ConnectClusters(&c1, iHinge1, &c2, iHinge2, SINGLE_BOUND)) return false;
	c2.parent = (BASIC_CLUSTER*)(&c1);
	c2.n2Parent = iHinge2;

	float dta = ta * fAngle2Radian - c2.km.ta;
	VECTOR3 Olink = *(c2.hinge[iHinge2].Olink);
	tweak_cluster_and_children(&c2, Olink, c2.hinge[iHinge2].vHinge, dta, true);

	if (!combine_clusters(&c1, &c2, &c)) return false;

	TREE<MATOM> atree;
	if (!construct_tree(&c, atree)) return false;

	ofstream out;
	out.open(ofname);
	if (!out.is_open()) {sprintf(msg, "can not open file %s", ofname); show_msg(msg); return false;}
	sprintf(msg, "#### %s, hinge %d ==> %s, hinge %d, with torsion angle %f", c2_name, iHinge2, c1_name, iHinge1, ta); out<<msg<<endl;
	sprintf(msg, "[%s] 0", c_name); out<<msg<<endl;
	sprintf(msg, "%d  %d  0", c.nAtoms, c.nHinges); out<<msg<<endl;
	out<<"ZMATRIX"<<endl;

	float blen = 0, bangle = 0;
	int bvalence = 0, btype = 0;
	VECTOR3 *v0 = NULL, *vb = NULL;
	VECTOR3 bvect, avect;

	LIST< LIST<MATOM> > *plevel = NULL;
	LIST<MATOM> *plist = NULL;
	MATOM *patom = NULL, *batom = NULL, *aatom = NULL, *tatom = NULL;
	patom = atree.tl->p->p;
	out<<patom->par->atom<<endl;

	bool status = true;
	int noutput = 1, nbatom = 0, naatom = 0, ntatom = 0, nb = 0;

#define iBOUND(p1, p2, nb) for (nb = 0; nb < p1->nb; nb++) {if (p1->ba[nb] == p2) break;}

	plevel = atree.tl->next;
	while (plevel != NULL) {
		plist = plevel->p;
		while (plist != NULL && plist->p != NULL) {
			if (!(status = search_connected_atoms(plist->p, plevel, &batom, &aatom, &tatom))) {
				if (noutput == 1 && batom != NULL) status = true;
				else if (noutput == 2 && batom != NULL && aatom != NULL) status = true;
			}
			if (!status) {
				sprintf(msg, "can not get enough bounded atoms for %dth atom -- %s", noutput, plist->p->par->atom);
				show_msg(msg); atree.release_tree(false); out.close(); return false;
			}
			VECT3(plist->p->r, batom->r, bvect) blen = bvect.abs();
			nbatom = c.iAtom(batom); iBOUND(plist->p, batom, nb)
			if (noutput == 1) {out<<plist->p->par->atom<<"  "<<blen<<"  "<<nbatom<<"  "<<int(plist->p->bv[nb])<<"  "<<int(plist->p->btype[nb])<<endl;}
			else if (noutput == 2) {
				naatom = c.iAtom(aatom); VECT3(aatom->r, batom->r, avect) bangle = angle(avect, bvect) * fRadian2Angle;
				out<<plist->p->par->atom<<"  "<<blen<<"  "<<nbatom<<"  "<<bangle<<"  "<<naatom<<"  "<<int(plist->p->bv[nb])<<"  "<<int(plist->p->btype[nb])<<endl;
			}
			else {
				naatom = c.iAtom(aatom); VECT3(aatom->r, batom->r, avect) bangle = angle(avect, bvect) * fRadian2Angle;
				ntatom = c.iAtom(tatom); ta = TorsionAngle(plist->p->r, batom->r, aatom->r, tatom->r) * fRadian2Angle;
				out<<plist->p->par->atom<<"  "<<blen<<"  "<<nbatom<<"  "<<bangle<<"  "<<naatom<<"  "<<ta<<"  "<<ntatom<<"  "<<int(plist->p->bv[nb])<<"  "<<int(plist->p->btype[nb])<<endl;
			}
			noutput++;
			plist = plist->next;
		}
		plevel = plevel->next;
	}

	// hinges
	if (c.nHinges == 0) {
		out.close(); atree.release_tree(false); return true;
	}

	// virtual atoms for hinges' positions
	MATOM *vatom = new MATOM[c.nHinges];
	ATOM_COMM_PARS vpar; // par for virtual atoms
	strcpy(vpar.atom, "nul");

#define vatom_indx(patom, i, indx) indx = -1; for (i = 0; i < c.nHinges; i++) {if (patom == vatom + i) {indx = i; break;}}

	int iatom = 0, i;
	for (iatom = 0; iatom < c.nHinges; iatom++) {
		vatom[iatom].par = &vpar;
		vatom[iatom].r = c.atom[c.hinge[iatom].indx].r + c.hinge[iatom].vHinge;
		vatom[iatom].set_bounds(1);
		vatom[iatom].ba[0] = c.atom + c.hinge[iatom].indx;
		vatom[iatom].bv[0] = SINGLE_BOUND;
		vatom[iatom].btype[0] = NORMAL_BOUND;
		c.atom[c.hinge[iatom].indx].ba[c.hinge[iatom].ibound] = vatom + iatom;

		// location of hinge bounded atom in atree
		if (!search<MATOM>(c.atom + c.hinge[iatom].indx, atree, &plevel, &plist)) {
			sprintf(msg, "STRANGE : can not find defined atom in the tree"); show_msg(msg);
			goto _COMBINE_CLUSTER_OVER;
		}
		// add this virtual atom to the next level of the level where bounded atom is
		plevel = plevel->next;
		if (plevel == NULL) {
			plevel = atree.new_level();
			plist = new LIST<MATOM>; plist->p = vatom + iatom; plevel->p = plist;
		}
		else {
			plevel->p->add(vatom + iatom, false);
		}

		// output the virtual atom position
		patom = vatom + iatom;
		if (!(status = search_connected_atoms(patom, plevel, &batom, &aatom, &tatom, true))) {
			if (noutput == 1 && batom != NULL) status = true;
			else if (noutput == 2 && batom != NULL && aatom != NULL) status = true;
		}
		if (!status) {
			sprintf(msg, "can not get enough bounded atoms for %dth hinge", iatom);
			show_msg(msg); atree.release_tree(false); goto _COMBINE_CLUSTER_OVER;
		}
		VECT3(patom->r, batom->r, bvect) blen = bvect.abs();
		if (noutput == 1) {out<<patom->par->atom<<"  "<<blen<<"  "<<c.hinge[iatom].indx<<"  "<<int(SINGLE_BOUND)<<"  "<<int(NORMAL_BOUND)<<endl;}
		else if (noutput == 2) {
			naatom = c.iAtom(aatom); if (naatom < 0) {vatom_indx(aatom, i, naatom) naatom += c.nAtoms;} 
			VECT3(aatom->r, batom->r, avect) bangle = angle(avect, bvect) * fRadian2Angle;
			out<<patom->par->atom<<"  "<<blen<<"  "<<c.hinge[iatom].indx<<"  "<<bangle<<"  "<<naatom<<"  "<<int(SINGLE_BOUND)<<"  "<<int(NORMAL_BOUND)<<endl;
		}
		else {
			naatom = c.iAtom(aatom); if (naatom < 0) {vatom_indx(aatom, i, naatom) naatom += c.nAtoms;}
			VECT3(aatom->r, batom->r, avect) bangle = angle(avect, bvect) * fRadian2Angle;
			ntatom = c.iAtom(tatom); if (ntatom < 0) {vatom_indx(tatom, i, ntatom) ntatom += c.nAtoms;}
			ta = TorsionAngle(patom->r, batom->r, aatom->r, tatom->r) * fRadian2Angle;
			out<<patom->par->atom<<"  "<<blen<<"  "<<c.hinge[iatom].indx<<"  "<<bangle<<"  "<<naatom<<"  "<<ta<<"  "<<ntatom<<"  "<<int(SINGLE_BOUND)<<"  "<<int(NORMAL_BOUND)<<endl;
		}
		noutput++;
	}

_COMBINE_CLUSTER_OVER:

	for (iatom = 0; iatom < c.nHinges; iatom++) {
		vatom[iatom].par = NULL;
		c.atom[c.hinge[iatom].indx].ba[c.hinge[iatom].ibound] = NULL;
	}

	out.close();
	atree.release_tree(false);
	if (vatom != NULL) delete[] vatom;
	return true;

#undef iBOUND
#undef vatom_indx
}
