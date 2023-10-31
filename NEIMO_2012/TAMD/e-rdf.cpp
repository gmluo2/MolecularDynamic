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
#include "Matrix.h"
#include "nhc.h"

#include "ZMatrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
using namespace _cmm_3d_;
#include "cluster.h"

#include "var.h"

#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"
#include "Interaction.h"
#include "Interact1.h"
#include "md_cell.h"
#include "cg-md-cell.h"

#define _READ_MM_MD_ 1
#include "read.h"

#include "ReadMDProc.h"
#include "mol-cell.h"
#include "show_md.h"

#if _SYS_ == _WINDOWS_SYS_ 
#include "..\export1.h"
#elif _SYS_ == _LINUX_SYS_
#include "export1.h"
#endif

using namespace _coarse_grain_;

extern "C"  bool construct_macromolecule(char *fname);
extern "C" bool init_read_mdproc(char *mdpar_fname);
extern bool serial_single_mol_md_proc(long isave);
extern bool serial_multi_mol_md_proc(long isave);
extern "C" bool serial_md_proc(long isave);
extern "C" void clear_read_mdproc();
extern "C" void set_md_mode(int mode);
extern "C" void set_job(int job);
extern void molecule_check_periodic_cell(MMOL_MD_CELL &mmcell);
extern char struct_fname[256];
extern MMOL_MD_CELL mm_md_cell;
extern int MD_MODE;
extern READ_MD_PROC_FILE rmdproc;

extern CGMM_MD_CELL cgmm_md_cell;

extern bool get_monomer_uid(char *monomer, int &uid);

extern bool read_given_id(char *fname, char* title, ARRAY<short>& id_cluster);

#include "distribution.h"
#include "rdf-distribt.h"

#include "e-rdf.h"

namespace _erdf_ {
using namespace _distribution_;
using namespace _MM_rdf_;

// calculate the interaction energy, Electric field, gradient electric field on each atom and LJ force
double LJ_INTERACT(VIRTUAL_MONOMER<CLUSTER> *vm1, VIRTUAL_MONOMER<CLUSTER> *vm2, CMM_CELL3D<BASIC_CLUSTER> *cmm) {
	double UTotal = 0;
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 dr, u;
	int i, j, k, na, na1, nc, nc1;
	VECTOR3 r0, r1, dr_cmm;
	double t;

	CHAIN< CMM_IMAGINE<BASIC_CLUSTER> > *ch_imag_cluster = NULL;
	CMM_IMAGINE<BASIC_CLUSTER> *cmm_cluster = NULL;
	BASIC_CLUSTER *pc1 = NULL;
	MATOM *patom1 = NULL;

	double R, R2;

	bool bound = false, bound14 = false;

	double U_LJ = 0;
	LJ_PARS *pLJ = NULL;

	for (nc = 0; nc < vm1->g.nUnit; nc++) {
		pc = (BASIC_CLUSTER*)(vm1->g.u[nc]);
		for (nc1 = 0; nc1 < vm2->g.nUnit; nc1++) {
			pc1 = (BASIC_CLUSTER*)(vm2->g.u[nc1]);
			//if ((cmm_cluster = in_imagine_rship<BASIC_CLUSTER>(pc->LJRship.ch, pc1)) == NULL) continue;
			if ((cmm_cluster = pc->LJRship.get_rship(pc1)) == NULL) continue;
			dr_cmm.v[0] = cmm_cluster->nx * cmm->xd[0];
			dr_cmm.v[1] = cmm_cluster->ny * cmm->xd[1];
			dr_cmm.v[2] = cmm_cluster->nz * cmm->xd[2];
			for (na = 0; na < pc->nAtoms; na++) {
				patom = pc->atom + na;
				if (pc->bCMMShift) {V3plusV3(patom->r, pc->dr_cmm, r0)}
				else {V32V3(patom->r, r0)}

				// LJ interaction
				for (na1 = 0; na1 < pc1->nAtoms; na1++) {
					patom1 = pc1->atom + na1;
					if (pc1->bCMMShift) {V3plusV3(patom1->r, pc1->dr_cmm, r1)}
					else {V32V3(patom1->r, r1)}
					VECT3(r0, r1, dr)
					V3plusV3(dr, dr_cmm, dr)
					//ABS2(dr, i, R2)
					R2 = dr.v[0] * dr.v[0]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[1] * dr.v[1]; if (R2 > r2cut_LJ) continue;
					R2 += dr.v[2] * dr.v[2]; if (R2 > r2cut_LJ) continue;

					if (R2 < d2_14 && pc->mIndx == pc1->mIndx) {
						BOUND(patom, patom1, i, bound) // bound atom
						if (!bound) {BOUND_SAME_ATOM(patom, patom1, i, j, bound)} // bounding to same atom ?
						if (!bound) {BOUND_CLOSE_ATOMS(patom, patom1, i, j, k, bound14)} // 1-4 interaction ?
						else bound14 = false;
					}
					else {bound = false; bound14 = false;}
					if (bound) continue; // ignore bound interaction

					if (R2 < _R2min_Interact) continue; // the two atoms are too close

					pLJ = &(LJPARS.m[patom->aindx].m[patom1->aindx]);
					if (pLJ != NULL && pLJ->epsLJ > 0.0001) {
						//LJ_V3FORCE(pLJ->epsLJ, pLJ->rLJ2, dr, R2, t1, force, t2)
						R = sqrt(R2);
						LJ_V(pLJ->epsLJ, pLJ->rLJ, R, i, t, U_LJ) // LJ12_6 profile has eps unitkJ/mol ==> kT
						if (bound14) { // 1-4 interaction, LJ /= 2
							U_LJ *= A14_LJ;
						}
						//U_LJ += U_LJ * 0.5; // all clusters are in same system
						UTotal += U_LJ; // will be corrected later
					}
				}
			}
		}
	}
	return UTotal;
}

// this function is to calculate the implicit solvation on vm1, with vm2 at a specific location
// important : implicit-solvation is not pairwise interaction
double ImplicitSolvation(VIRTUAL_MONOMER<CLUSTER> *vm1, VIRTUAL_MONOMER<CLUSTER> *vm2, CMM_CELL3D<BASIC_CLUSTER> *cmm) {
	double UTotal = 0, U = 0;
	BASIC_CLUSTER *pc = NULL;
	int nc, nc1;
	BASIC_CLUSTER *pc1 = NULL;

	for (nc = 0; nc < vm1->g.nUnit; nc++) {
		pc = (BASIC_CLUSTER*)(vm1->g.u[nc]);
		for (nc1 = 0; nc1 < vm2->g.nUnit; nc1++) {
			pc1 = (BASIC_CLUSTER*)(vm2->g.u[nc1]);
			U = ClusterImplicitSolvationEnergy(pc, pc1, cmm);
			// because ClusterImplicitSolvationEnergy calculate the implicit-solvation energy
			// of pc, with pc1 at a specific position, it averaged the implicit-solvation interacton
			// from the enviromental clusters, including those in vm2
			if (FABS(U) > 1e-6) break;
		}
		UTotal += U;
		// it is not correct to sum the implicit-solvation energy of all the clusters
		// in the vm, to say the sum is the total implicit-solvation energy of vm
		// but do we have better idear?
	}
	return UTotal;
}

void VM_INTERACT(VIRTUAL_MONOMER<CLUSTER> *vm1, VIRTUAL_MONOMER<CLUSTER> *vm2, CMM_CELL3D<BASIC_CLUSTER> *cmm, double& U_LJ, double &U_ImpSolv) {
	U_LJ = LJ_INTERACT(vm1, vm2, cmm);
	U_ImpSolv = ImplicitSolvation(vm1, vm2, cmm);
}

void single_mm_neighbor_erdf_func(NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER> *erdf_thread_var, int nrdf) {
	return single_vm_neighbor_erdf_func<MMOLECULE, CLUSTER, BASIC_CLUSTER>(erdf_thread_var, nrdf);
}

void multi_mm_neighbor_erdf_func(MMOL_MD_CELL *md_cell, ARRAY< MM_VM<MMOLECULE, CLUSTER> >* vmarray, NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER> *erdf_thread_var, int nrdf) {
	return multi_vm_neighbor_erdf_func<MMOL_MD_CELL, MMOLECULE, CLUSTER, BASIC_CLUSTER>(md_cell, vmarray, erdf_thread_var, nrdf);
}

void single_mm_non_neighbor_erdf_func(NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER> *erdf_thread_var, int nrdf) {
	return single_vm_non_neighbor_erdf_func<MMOLECULE, CLUSTER, BASIC_CLUSTER>(erdf_thread_var, nrdf);
}

void multi_mm_non_neighbor_erdf_func(MMOL_MD_CELL *md_cell, ARRAY< MM_VM<MMOLECULE, CLUSTER> >* vmarray, NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER> *erdf_thread_var, int nrdf) {
	return multi_vm_non_neighbor_erdf_func<MMOL_MD_CELL, MMOLECULE, CLUSTER, BASIC_CLUSTER>(md_cell, vmarray, erdf_thread_var, nrdf);
}

void raw_accumulate_rdf(NEIGHBOR_ERDF_VAR *rdf, NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER>* nrdf) {
	int irdf;
	/*
	for (irdf = 0; irdf < rdf->sdt.ndt; irdf++) {
		raw_distribt_combine<float, long, float, long>(rdf->sdt.dt[irdf], nrdf->sdt.dt[irdf]);
		rdf->sdt.dt[irdf].vnormal += nrdf->sdt.dt[irdf].vnormal;
	}
	*/
	for (irdf = 0; irdf < rdf->eLJ_dt.ndt; irdf++) {
		raw_distribt_combine<float, float, float, float>(rdf->eLJ_dt.dt[irdf], nrdf->eLJ_dt.dt[irdf]);
		rdf->eLJ_dt.dt[irdf].vnormal += nrdf->eLJ_dt.dt[irdf].vnormal;
	}
	for (irdf = 0; irdf < rdf->eImpSolv_dt.ndt; irdf++) {
		raw_distribt_combine<float, float, float, float>(rdf->eImpSolv_dt.dt[irdf], nrdf->eImpSolv_dt.dt[irdf]);
		rdf->eImpSolv_dt.dt[irdf].vnormal += nrdf->eImpSolv_dt.dt[irdf].vnormal;
	}
	for (irdf = 0; irdf < rdf->edt.ndt; irdf++) {
		raw_distribt_combine<float, float, float, float>(rdf->edt.dt[irdf], nrdf->edt.dt[irdf]);
		rdf->edt.dt[irdf].vnormal += nrdf->edt.dt[irdf].vnormal;
	}
}

void normalize_rdf(NEIGHBOR_ERDF_VAR *rdf) {
	int irdf, i;
	DISTRIBT<float, float> *pdt;
	for (irdf = 0; irdf < rdf->edt.ndt; irdf++) {
		for (i = 0; i < rdf->edt.dt[irdf].y.N; i++) {
			/*
			//rdf->edt.dt[irdf].y.v[i] /= rdf->sdt.dt[irdf].vnormal;
			if (rdf->sdt.dt[irdf].y.v[i] > 0) {
				rdf->eLJ_dt.dt[irdf].y.v[i] /= rdf->sdt.dt[irdf].y.v[i];
				rdf->eImpSolv_dt.dt[irdf].y.v[i] /= rdf->sdt.dt[irdf].y.v[i];
				rdf->edt.dt[irdf].y.v[i] /= rdf->sdt.dt[irdf].y.v[i];
			}
			*/
			pdt = rdf->eLJ_dt.dt + irdf;
			if (pdt->sn.v[i] > 0) {
				pdt->eb.v[i] = sqrt((double)(pdt->sn.v[i])) * pdt->y.v[i] / pdt->sn.v[i];
				pdt->y.v[i] /= pdt->sn.v[i];
				pdt->eb.v[i] /= pdt->sn.v[i];
			}
			pdt = rdf->eImpSolv_dt.dt + irdf;
			if (pdt->sn.v[i] > 0) {
				pdt->eb.v[i] = sqrt((double)(pdt->sn.v[i])) * pdt->y.v[i] / pdt->sn.v[i];
				pdt->y.v[i] /= pdt->sn.v[i];
				pdt->eb.v[i] /= pdt->sn.v[i];
			}
			pdt = rdf->edt.dt + irdf;
			if (pdt->sn.v[i] > 0) {
				pdt->eb.v[i] = sqrt((double)(pdt->sn.v[i])) * pdt->y.v[i] / pdt->sn.v[i];
				pdt->y.v[i] /= pdt->sn.v[i];
				pdt->eb.v[i] /= pdt->sn.v[i];
			}
		}
	}
}

void save_rdf(NEIGHBOR_ERDF_VAR *rdf_var, char *title, char *msg_title) {
	char ofname[256] = "\0", msg[256] = "\0";
	ofstream *out = NULL;

	int i, n;
	for (i = 0; i < rdf_var->edt.ndt; i++) {
		sprintf(ofname, "%s-%d-LJ.dat", title, i);
		out = new ofstream;
		out->open(ofname);
		/*
		for (n = 0; n < rdf_var->eLJ_dt.dt[i].N; n++) {
			(*out)<<rdf_var->eLJ_dt.dt[i].x[n]<<" "<<rdf_var->eLJ_dt.dt[i].y[n];
			if (rdf_var->sdt.dt[i].y[n] > 0) 
				(*out)<<" "<<FABS(rdf_var->eLJ_dt.dt[i].y[n]) / sqrt((double)(rdf_var->sdt.dt[i].y[n]))<<endl;
			else (*out)<<" 0"<<endl;
		}
		*/
		rdf_var->eLJ_dt.dt[i].save_x_y_eb_s(*out);
		out->close(); delete out; out = NULL;
		sprintf(msg, "LJ interaction : %s # %d ==> %s", msg_title, i, ofname);
		show_log(msg, true);

		sprintf(ofname, "%s-%d-ImpSolv.dat", title, i);
		out = new ofstream;
		out->open(ofname);
		/*
		for (n = 0; n < rdf_var->eImpSolv_dt.dt[i].N; n++) {
			(*out)<<rdf_var->eImpSolv_dt.dt[i].x[n]<<" "<<rdf_var->eImpSolv_dt.dt[i].y[n];
			if (rdf_var->sdt.dt[i].y[n] > 0) 
				(*out)<<" "<<FABS(rdf_var->eImpSolv_dt.dt[i].y[n]) / sqrt((double)(rdf_var->sdt.dt[i].y[n]))<<endl;
			else (*out)<<" 0"<<endl;
		}
		*/
		rdf_var->eImpSolv_dt.dt[i].save_x_y_eb_s(*out);
		out->close(); delete out; out = NULL;
		sprintf(msg, "Implicit-Solvation interaction : %s # %d ==> %s", msg_title, i, ofname);
		show_log(msg, true);

		sprintf(ofname, "%s-%d-e.dat", title, i);
		out = new ofstream;
		out->open(ofname);
		/*
		for (n = 0; n < rdf_var->edt.dt[i].N; n++) {
			(*out)<<rdf_var->edt.dt[i].x[n]<<" "<<rdf_var->edt.dt[i].y[n];
			if (rdf_var->sdt.dt[i].y[n] > 0) 
				(*out)<<" "<<FABS(rdf_var->edt.dt[i].y[n]) / sqrt((double)(rdf_var->sdt.dt[i].y[n]))<<endl;
			else (*out)<<" 0"<<endl;
		}
		*/
		rdf_var->edt.dt[i].save_x_y_eb_s(*out);
		out->close(); delete out; out = NULL;
		sprintf(msg, "total interaction : %s # %d ==> %s", msg_title, i, ofname);
		show_log(msg, true);
	}
	return;
}


void init_neighbor_erdf_thread_var(NEIGHBOR_ERDF_VAR* erdf_var, NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER>* erdf_thread_var) {
	init_vm_neighbor_rdf_thread_var<MMOLECULE, CLUSTER>((NEIGHBOR_RDF_VAR*)erdf_var, (NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER>*)erdf_thread_var);
	
	int i = 0;
	//erdf_thread_var->sdt.set(erdf_var->sdt.ndt);
	erdf_thread_var->eLJ_dt.set(erdf_var->eLJ_dt.ndt);
	erdf_thread_var->eImpSolv_dt.set(erdf_var->eImpSolv_dt.ndt);
	erdf_thread_var->edt.set(erdf_var->edt.ndt);
	/*
	for (i = 0; i < erdf_thread_var->sdt.ndt; i++) {
		erdf_thread_var->sdt.dt[i].set_distribt_x(erdf_var->sdt.dt[i].xf, erdf_var->sdt.dt[i].xt, erdf_var->sdt.dt[i].x.N - 1);
	}
	*/
	for (i = 0; i < erdf_thread_var->eLJ_dt.ndt; i++) {
		erdf_thread_var->eLJ_dt.dt[i].set_distribt_x(erdf_var->eLJ_dt.dt[i].xf, erdf_var->eLJ_dt.dt[i].xt, erdf_var->eLJ_dt.dt[i].x.N - 1);
	}
	for (i = 0; i < erdf_thread_var->eImpSolv_dt.ndt; i++) {
		erdf_thread_var->eImpSolv_dt.dt[i].set_distribt_x(erdf_var->eImpSolv_dt.dt[i].xf, erdf_var->eImpSolv_dt.dt[i].xt, erdf_var->eImpSolv_dt.dt[i].x.N - 1);
	}
	for (i = 0; i < erdf_thread_var->edt.ndt; i++) {
		erdf_thread_var->edt.dt[i].set_distribt_x(erdf_var->edt.dt[i].xf, erdf_var->edt.dt[i].xt, erdf_var->edt.dt[i].x.N - 1);
	}
}

void init_all_neighbors_erdf_thread_var(NEIGHBOR_ERDF_VAR* erdf_var, NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER>* erdf_thread_var) {
	init_vm_all_neighbors_rdf_thread_var<MMOLECULE, CLUSTER>((NEIGHBOR_RDF_VAR*)erdf_var, (NEIGHBOR_RDF_THREAD_VAR<MMOLECULE, CLUSTER>*)erdf_thread_var);

	int i;
	//erdf_thread_var->sdt.set_dt(erdf_var->sdt.ndt);
	erdf_thread_var->eLJ_dt.set(erdf_var->eLJ_dt.ndt);
	erdf_thread_var->eImpSolv_dt.set(erdf_var->eImpSolv_dt.ndt);
	erdf_thread_var->edt.set(erdf_var->edt.ndt);
	/*
	for (i = 0; i < erdf_thread_var->sdt.ndt; i++) {
		erdf_thread_var->sdt.dt[i].set_distribt_x(erdf_var->sdt.dt[i].xf, erdf_var->sdt.dt[i].xt, erdf_var->sdt.dt[i].x.N - 1);
	}
	*/
	for (i = 0; i < erdf_thread_var->eLJ_dt.ndt; i++) {
		erdf_thread_var->eLJ_dt.dt[i].set_distribt_x(erdf_var->eLJ_dt.dt[i].xf, erdf_var->eLJ_dt.dt[i].xt, erdf_var->eLJ_dt.dt[i].x.N - 1);
	}
	for (i = 0; i < erdf_thread_var->eImpSolv_dt.ndt; i++) {
		erdf_thread_var->eImpSolv_dt.dt[i].set_distribt_x(erdf_var->eImpSolv_dt.dt[i].xf, erdf_var->eImpSolv_dt.dt[i].xt, erdf_var->eImpSolv_dt.dt[i].x.N - 1);
	}
	for (i = 0; i < erdf_thread_var->edt.ndt; i++) {
		erdf_thread_var->edt.dt[i].set_distribt_x(erdf_var->edt.dt[i].xf, erdf_var->edt.dt[i].xt, erdf_var->edt.dt[i].x.N - 1);
	}
}

bool read_neighbor_erdf_pars(char *parfname, NEIGHBOR_ERDF_VAR** erdf_var, int &nrdf) {
	NEIGHBOR_RDF_VAR *rdf_var = NULL;
	if (!read_neighbor_rdf_pars(parfname, &rdf_var, nrdf) || rdf_var == NULL) return false;
	*erdf_var = new NEIGHBOR_ERDF_VAR[nrdf];
	int i, j;
	for (i = 0; i < nrdf; i++) {
		(*erdf_var)[i].cm_uid = rdf_var[i].cm_uid;
		(*erdf_var)[i].nm_uid = rdf_var[i].nm_uid;
		(*erdf_var)[i].nth_neighbor = rdf_var[i].nth_neighbor;
		//(*erdf_var)[i].sdt.set_dt(rdf_var[i].dt.ndt);
		(*erdf_var)[i].eLJ_dt.set(rdf_var[i].dt.ndt);
		(*erdf_var)[i].eImpSolv_dt.set(rdf_var[i].dt.ndt);
		(*erdf_var)[i].edt.set(rdf_var[i].dt.ndt);
		/*
		for (j = 0; j < (*erdf_var)[i].sdt.ndt; j++) {
			(*erdf_var)[i].sdt.dt[j].set_distribt_x(rdf_var[i].dt.dt[j].xf, rdf_var[i].dt.dt[j].xt, rdf_var[i].dt.dt[j].x.N - 1);
		}
		*/
		for (j = 0; j < (*erdf_var)[i].edt.ndt; j++) {
			(*erdf_var)[i].eLJ_dt.dt[j].set_distribt_x(rdf_var[i].dt.dt[j].xf, rdf_var[i].dt.dt[j].xt, rdf_var[i].dt.dt[j].x.N - 1);
			(*erdf_var)[i].eImpSolv_dt.dt[j].set_distribt_x(rdf_var[i].dt.dt[j].xf, rdf_var[i].dt.dt[j].xt, rdf_var[i].dt.dt[j].x.N - 1);
			(*erdf_var)[i].edt.dt[j].set_distribt_x(rdf_var[i].dt.dt[j].xf, rdf_var[i].dt.dt[j].xt, rdf_var[i].dt.dt[j].x.N - 1);
		}
	}
	delete[] rdf_var; rdf_var = NULL;
	return true;
}

bool read_non_neighbor_erdf_pars(char *parfname, NEIGHBOR_ERDF_VAR** erdf_var, int &nrdf) {
	NEIGHBOR_RDF_VAR *rdf_var = NULL;
	if (!read_non_neighbor_rdf_pars(parfname, &rdf_var, nrdf) || rdf_var == NULL) return false;
	*erdf_var = new NEIGHBOR_ERDF_VAR[nrdf];
	int i, j;
	for (i = 0; i < nrdf; i++) {
		(*erdf_var)[i].cm_uid = rdf_var[i].cm_uid;
		(*erdf_var)[i].nm_uid = rdf_var[i].nm_uid;
		(*erdf_var)[i].nth_neighbor = rdf_var[i].nth_neighbor;
		//(*erdf_var)[i].sdt.set(rdf_var[i].dt.ndt);
		(*erdf_var)[i].eLJ_dt.set(rdf_var[i].dt.ndt);
		(*erdf_var)[i].eImpSolv_dt.set(rdf_var[i].dt.ndt);
		(*erdf_var)[i].edt.set(rdf_var[i].dt.ndt);
		/*
		for (j = 0; j < (*erdf_var)[i].sdt.ndt; j++) {
			(*erdf_var)[i].sdt.dt[j].set_distribt_x(rdf_var[i].dt.dt[j].xf, rdf_var[i].dt.dt[j].xt, rdf_var[i].dt.dt[j].x.N - 1);
		}
		*/
		for (j = 0; j < (*erdf_var)[i].edt.ndt; j++) {
			(*erdf_var)[i].eLJ_dt.dt[j].set_distribt_x(rdf_var[i].dt.dt[j].xf, rdf_var[i].dt.dt[j].xt, rdf_var[i].dt.dt[j].x.N - 1);
			(*erdf_var)[i].eImpSolv_dt.dt[j].set_distribt_x(rdf_var[i].dt.dt[j].xf, rdf_var[i].dt.dt[j].xt, rdf_var[i].dt.dt[j].x.N - 1);
			(*erdf_var)[i].edt.dt[j].set_distribt_x(rdf_var[i].dt.dt[j].xf, rdf_var[i].dt.dt[j].xt, rdf_var[i].dt.dt[j].x.N - 1);
		}
	}
	delete[] rdf_var; rdf_var = NULL;
	return true;
}

} // end of namespace _erdf_

extern bool ReadMoleculeStruct(ifstream &in, MMOLECULE &mm);
extern bool serial_multi_mol_md_proc(long isave);

extern "C" bool single_mm_neighbor_erdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;
	using namespace _erdf_;
	
	::COMMAND_MD_STOP = false;
	mlog.init("erdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_ERDF_VAR *neighbor_erdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_neighbor_erdf_pars(parfname, &neighbor_erdf_var, nrdf) || neighbor_erdf_var == NULL) return false;

	erdf_distribt_single_mm_multi_md<MMOLECULE, CLUSTER, BASIC_CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_ERDF_VAR, NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER> >(
		parfname, (void*)(&ReadMoleculeStruct), (void*)(&single_mm_cal_Position), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_erdf_var, nrdf, (void*)(&init_neighbor_erdf_thread_var), (void*)(&single_mm_neighbor_erdf_func), (void*)(&VM_INTERACT), 
		(void*)(&_erdf_::raw_accumulate_rdf), (void*)(&_erdf_::normalize_rdf), (void*)(&_erdf_::save_rdf));

	mlog.close();

	if (neighbor_erdf_var != NULL) {delete[] neighbor_erdf_var; neighbor_erdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_mm_neighbor_erdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;
	using namespace _erdf_;

	::COMMAND_MD_STOP = false;
	mlog.init("erdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_ERDF_VAR *neighbor_erdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_neighbor_erdf_pars(parfname, &neighbor_erdf_var, nrdf) || neighbor_erdf_var == NULL) return false;

	erdf_distribt_multi_mm_multi_md<MMOL_MD_CELL, MMOLECULE, CLUSTER, BASIC_CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_ERDF_VAR, NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER> >(
		parfname, ::mm_md_cell, (void*)(&serial_multi_mol_md_proc), (void*)(&multi_mm_cal_Position_NonCell), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_erdf_var, nrdf, (void*)(&init_neighbor_erdf_thread_var), (void*)(&multi_mm_neighbor_erdf_func), (void*)(&VM_INTERACT), 
		(void*)(&_erdf_::raw_accumulate_rdf), (void*)(&_erdf_::normalize_rdf), (void*)(&_erdf_::save_rdf));

	mlog.close();

	if (neighbor_erdf_var != NULL) {delete[] neighbor_erdf_var; neighbor_erdf_var = NULL;} nrdf = 0;
	return true;
}


extern "C" bool single_mm_non_neighbor_erdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;
	using namespace _erdf_;

	::COMMAND_MD_STOP = false;
	mlog.init("erdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_ERDF_VAR *neighbor_erdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_erdf_pars(parfname, &neighbor_erdf_var, nrdf) || neighbor_erdf_var == NULL) return false;

	erdf_distribt_single_mm_multi_md<MMOLECULE, CLUSTER, BASIC_CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_ERDF_VAR, NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER> >(
		parfname, (void*)(&ReadMoleculeStruct), (void*)(&single_mm_cal_Position), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_erdf_var, nrdf, (void*)(&init_all_neighbors_erdf_thread_var), (void*)(&single_mm_non_neighbor_erdf_func), (void*)(&VM_INTERACT),
		(void*)(&_erdf_::raw_accumulate_rdf), (void*)(&_erdf_::normalize_rdf), (void*)(&_erdf_::save_rdf));

	mlog.close();

	if (neighbor_erdf_var != NULL) {delete[] neighbor_erdf_var; neighbor_erdf_var = NULL;} nrdf = 0;
	return true;
}

extern "C" bool multi_mm_non_neighbor_erdf_single_chain_block_copolymer(char *parfname) {
	using namespace _distribution_;
	using namespace _MM_rdf_;
	using namespace _erdf_;

	::COMMAND_MD_STOP = false;
	mlog.init("erdf_neighbor.log", false);
	VM_VAR<MMOLECULE> vm_var;
	NEIGHBOR_ERDF_VAR *neighbor_erdf_var = NULL;
	int nrdf = 0;

	if (!read_vm_var(parfname, vm_var)) return false;
	if (!read_non_neighbor_erdf_pars(parfname, &neighbor_erdf_var, nrdf) || neighbor_erdf_var == NULL) return false;

	erdf_distribt_multi_mm_multi_md<MMOL_MD_CELL, MMOLECULE, CLUSTER, BASIC_CLUSTER, VM_VAR<MMOLECULE>, NEIGHBOR_ERDF_VAR, NEIGHBOR_ERDF_THREAD_VAR<MMOLECULE, CLUSTER, BASIC_CLUSTER> >(
		parfname, ::mm_md_cell, (void*)(&serial_multi_mol_md_proc), (void*)(&multi_mm_cal_Position_InCell), vm_var, (void*)(&construct_monomer_single_chain_block_copolymer), 
		neighbor_erdf_var, nrdf, (void*)(&init_all_neighbors_erdf_thread_var), (void*)(&multi_mm_non_neighbor_erdf_func), (void*)(&VM_INTERACT),
		(void*)(&_erdf_::raw_accumulate_rdf), (void*)(&_erdf_::normalize_rdf), (void*)(&_erdf_::save_rdf));

	mlog.close();

	if (neighbor_erdf_var != NULL) {delete[] neighbor_erdf_var; neighbor_erdf_var = NULL;} nrdf = 0;
	return true;
}
