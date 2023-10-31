#include "project.h"

#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "def.h"
#include "show.h"
#include "ranlib.h"
#include "vector.h"
#include "bound.h"
#include "Matrix.h"
#include "nhc.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "CMM_2d.h"
using namespace _cmm_3d_;

#include "var.h"

#include "Interaction.h"
#include "Interact1.h"
#include "Interact2.h"
#include "Interact3.h"
#include "NEIMO.h"
#include "MD.h"
#include "cg-mm.h"
#include "cg-md.h"
#include "MD_POLYMER.h"
#include "MD_SM.h"
#include "md_cell.h"
#include "complex.h"
#include "fftw3.h"
#include "EwaldSum.h"
#include "spme.h"
#include "spme_interact.h"
#include "spme_2d.h"
#include "spme_interact_2d.h"
extern BSplineFunc<2> bsp;

using namespace _cmm_3d_;

#include "interface.h"
using namespace _interface_constraint_;

#include "mdvar.h"

#include "time.h"

using namespace _evatom_;
using namespace _EwaldSum_real_;

#include "atom_sr.h"
#include "atom_sr_interact.h"

#include "gpu_vector.h"
#include "gpu_spme.h"
#include "gpu_spme_2d.h"
#include "gpu_interact.h"
#include "gpu_interact_2d.h"

namespace _gpu_md_ {
#ifdef __GPU__
	//bool init_GPU_AtCELL(_atomic_cmm_::_AtCELL *acell, GPU_AtCELL *gpu_acell);
	bool init_GPU_EAtCELL(_atomic_cmm_::EAtCELL *acell, GPU_EAtCELL *gpu_acell);
	bool init_GPU_SR_AtCELL(_atomic_cmm_::SR_AtCELL *acell, GPU_SR_AtCELL *gpu_acell);
#endif
}

namespace _atomic_cmm_ {
#ifdef __GPU__
	void _AtCell_use_gpu(_AtCELL &acell) {
		acell.aID.release(); acell.aID.use_GPU_MapMem(true); 
		acell.cIndx.release(); acell.cIndx.use_GPU_MapMem(true);
		acell.cluster.release(); acell.cluster.use_GPU_MapMem(true); 
		acell.r.release(); acell.r.use_GPU_MapMem(true);
		acell.ir.release(); acell.ir.use_GPU_MapMem(true);
		acell.cmm.release(); acell.cmm.use_GPU_MapMem(true);

		acell.bound.release(); acell.bound.use_GPU_MapMem(true);
		acell.bound_14.release(); acell.bound_14.use_GPU_MapMem(true);

		acell.b_relship.release(); acell.b_relship.use_GPU_MapMem(true);
		acell.b14_relship.release(); acell.b14_relship.use_GPU_MapMem(true);

		acell.abound.release(); acell.abound.use_GPU_MapMem(true);

		acell.ivar.release(); acell.ivar.use_GPU_MapMem(true);
		acell.dvar.release(); acell.dvar.use_GPU_MapMem(true);

		acell.b_relship.release(); acell.b_relship.use_GPU_MapMem(true);
		acell.b14_relship.release(); acell.b14_relship.use_GPU_MapMem(true);
	};

	void EAtCell_use_gpu(EAtCELL &acell) {
		// all these results needs to be GPU-mapped memory
#define use_gpu(c) c.release(); c.use_GPU_MapMem(true);
		use_gpu(acell.c);
		use_gpu(acell.mu);

		use_gpu(acell.Fe);
		use_gpu(acell.U_estat);
		use_gpu(acell.E);
		use_gpu(acell.ve);

		use_gpu(acell.Fe_b);
		use_gpu(acell.Ub_estat);
		use_gpu(acell.E_b);
		use_gpu(acell.ve_b);

		use_gpu(acell.Fe_14);
		use_gpu(acell.U14_estat);
		use_gpu(acell.E_14);
		use_gpu(acell.ve_14);
#undef use_gpu
	}

	void SRAtCell_use_gpu(SR_AtCELL &acell) {
		// all these results needs to be GPU-mapped memory
#define use_gpu(c) c.release(); c.use_GPU_MapMem(true);
		use_gpu(acell.sr_type);
		use_gpu(acell.sr_par);

		use_gpu(acell.bShocked);

		use_gpu(acell.Fsr);
		use_gpu(acell.U_sr);
		use_gpu(acell.vsr);

		use_gpu(acell.Fsr_b);
		use_gpu(acell.Ub_sr);
		use_gpu(acell.vsr_b);

		use_gpu(acell.Fsr_14);
		use_gpu(acell.U14_sr);
		use_gpu(acell.vsr_14);
#undef use_gpu
}

#endif // __GPU__ for gpu_init_cell

#ifndef dcureset
#define dcureset(d, n) if (d != NULL && n > 0) cuMemsetD8((CUdeviceptr)(d), 0, n * SIZE_DOUBLE)
	void EAtCELL::reset_efield() {
		memset(this->E.m, 0, this->E.n * SIZE_DOUBLE); 
		memset(this->E_b.m, 0, this->E_b.n * SIZE_DOUBLE); 
		memset(this->E_14.m, 0, this->E_14.n * SIZE_DOUBLE);
#ifdef __GPU__
		//dcureset(gpu_cell.E_dev, this->E.n);
		_gpu_md_::reset_efield(&(this->gpu_cell), this->gpu_cell_dev);
		dcureset(gpu_cell.E_b_dev, this->E_b.n);
		dcureset(gpu_cell.E_14_dev, this->E_14.n);
#endif
	};
	void EAtCELL::reset_eforce() {
		memset(this->Fe.m, 0, this->Fe.n * SIZE_DOUBLE); 
		memset(this->Fe_b.m, 0, this->Fe_b.n * SIZE_DOUBLE);
		memset(this->Fe_14.m, 0, this->Fe_14.n * SIZE_DOUBLE);
		memset(this->ve.m, 0, this->ve.n * SIZE_DOUBLE);
		memset(this->ve_b.m, 0, this->ve_b.n * SIZE_DOUBLE);
		memset(this->ve_14.m, 0, this->ve_14.n * SIZE_DOUBLE);
		memset(this->U_estat.m, 0, this->U_estat.n * SIZE_DOUBLE); 
		memset(this->Ub_estat.m, 0, this->Ub_estat.n * SIZE_DOUBLE); 
		memset(this->U14_estat.m, 0, this->U14_estat.n * SIZE_DOUBLE);
#ifdef __GPU__
		dcureset(gpu_cell.Fe_dev, this->Fe.n);
		dcureset(gpu_cell.Fe_b_dev, this->Fe_b.n);
		dcureset(gpu_cell.Fe_14_dev, this->Fe_14.n);
		dcureset(gpu_cell.ve_dev, this->ve.n);
		dcureset(gpu_cell.ve_b_dev, this->ve_b.n);
		dcureset(gpu_cell.ve_14_dev, this->ve_14.n);
		dcureset(gpu_cell.U_estat_dev, this->U_estat.n);
		dcureset(gpu_cell.Ub_estat_dev, this->Ub_estat.n);
		dcureset(gpu_cell.U14_estat_dev, this->U14_estat.n);
#endif
	};
	void SR_AtCELL::reset_sr() {
		memset(this->bShocked.m, 0, this->bShocked.n * SIZE_CHAR);
		memset(this->Fsr.m, 0, this->Fsr.n * SIZE_DOUBLE);
		memset(this->Fsr_b.m, 0, this->Fsr_b.n * SIZE_DOUBLE);
		memset(this->Fsr_14.m, 0, this->Fsr_14.n * SIZE_DOUBLE);
		memset(this->vsr.m, 0, this->vsr.n * SIZE_DOUBLE);
		memset(this->vsr_b.m, 0, this->vsr_b.n * SIZE_DOUBLE);
		memset(this->vsr_14.m, 0, this->vsr_14.n * SIZE_DOUBLE);
		memset(this->U_sr.m, 0, this->U_sr.n * SIZE_DOUBLE);
		memset(this->Ub_sr.m, 0, this->Ub_sr.n * SIZE_DOUBLE);
		memset(this->U14_sr.m, 0, this->U14_sr.n * SIZE_DOUBLE);
#ifdef __GPU__
		// not necessary to reset the value here, and they will be reset in the related GPU function
		/*
		dcureset(gpu_cell.bShocked_dev, this->bShocked.n);
		//dcureset(gpu_cell.Fsr_dev, this->Fsr.n);
		dcureset(gpu_cell.Fsr_b_dev, this->Fsr_b.n);
		dcureset(gpu_cell.Fsr_14_dev, this->Fsr_14.n);
		//dcureset(gpu_cell.vsr_dev, this->vsr.n);
		dcureset(gpu_cell.vsr_b_dev, this->vsr_b.n);
		dcureset(gpu_cell.vsr_14_dev, this->vsr_14.n);
		//dcureset(gpu_cell.U_sr_dev, this->U_sr.n);
		dcureset(gpu_cell.Ub_sr_dev, this->Ub_sr.n);
		dcureset(gpu_cell.U14_sr_dev, this->U14_sr.n);
		*/
#endif
	};
#undef dcureset
#endif

#ifdef __GPU__
	bool EAtCELL::check_neighbor_memory (bool b_rcut_changed) {
#ifndef cuda_release
#define cuda_release(v) if (v != NULL) {cudaFree(v); v = NULL;}
		char buffer[256] = "\0";

		if (!bUseGPU || !bCalUseGPU) return true;
		int n = 0;
		//if (::bPolarize) {
			if (b_rcut_changed) {
				gpu_cell.rcut = this->rcut;
				cuda_release(gpu_cell.srnb_dev); cuda_release(gpu_cell.dr_dev);
			}
			if (this->Na == 0) return true;
			if (gpu_cell.n_srnb == 0 || gpu_cell.srnb_dev == NULL) {
				n = int(gpu_cell.rcut * gpu_cell.rcut * gpu_cell.rcut * 4 * PI / 3); //assuming each atom volume is 1
				if (n == 0) {
					show_log("ERROR: cut-out distance for short-range interaction is not given.", true);
					show_log("ERROR: failure to define neighbor numbers for each atom.", true);
					return false;
				}

				show_log("", true); show_log("---- neighbor memory initialization ----", true);
				size_t size_free, size_total, size2_free;
				cudaMemGetInfo(&size_free, &size_total);
				sprintf(errmsg, "device memory avaliable: %d k", int(size_free / 1024)); show_log(errmsg, true);

				if (gpu_cell.n_srnb > 1000) gpu_cell.n_srnb = 1000;
				if (n != gpu_cell.n_srnb) {
					gpu_cell.n_srnb = n;
					sprintf(buffer, "neighbor memory on device to be used: %d ", (gpu_cell.n_srnb + 1) * this->Na * sizeof(int)); show_log(buffer, true);
					if (cudaMalloc(&(gpu_cell.srnb_dev), sizeof(int) * (gpu_cell.n_srnb + 1) * this->Na) != cudaSuccess) {return false;}
					if (cudaMalloc(&(gpu_cell.dr_dev), sizeof(DOUBLE) * 4 * gpu_cell.n_srnb * Na) != cudaSuccess) {return false;}
					cudaMemcpy(gpu_cell_dev, &gpu_cell, sizeof(_gpu_md_::GPU_EAtCELL), cudaMemcpyHostToDevice);
					show_log("initialized neighbor memory on device", true);

					cudaMemGetInfo(&size2_free, &size_total);
					sprintf(errmsg, "device memory avaliable: %d k", int(size2_free / 1024)); show_log(errmsg, true);
					if (size2_free == 0) {
						show_log("No free memory on device", true); return false;
					}
					sprintf(buffer, "used memory : %d k", (size_free - size2_free) / 1024); show_log(buffer, true);
					show_log("---- neighbor memory initialization ----\n", true);

				}
			}
			return true;
		//}
		//else return true;
	}
#undef cuda_release
#endif
#endif  //__GPU__

// to initialize the cell, require estimated parameters:
//     Nb -- maximum number of bounds of each atom
//     Nca -- maximum number of atoms in each cluster
//     NC -- maximum number of atoms in each cell

// other parameters: Na -- total number of atoms
//                   Nc -- total number of cluster
//                   Nx, Ny, Nz -- dimension of CMM cell
//                   cx[3], cw[3] -- starting position of CMM, full-width of each small cell
	void init_AtCell(/*MMOL_MD_CELL &mdcell, */_AtCELL &acell, int NC, int Nx, int Ny, int Nz, float *cx, float *cw) {
#ifdef __GPU__
		_AtCell_use_gpu(acell);
#endif
		// clusters -- not required for GPU 
		int Na = acell.satom.n, Nc = acell.scluster.n;
		int Nca = 0, i;
		for (i = 0; i < acell.scluster.n; i++) {if (acell.scluster.m[i].m->nAtoms > Nca) Nca = acell.scluster.m[i].m->nAtoms;}
		acell.cmm_c.set_array(Nx * Ny * Nz);
		for (i = 0; i < acell.cmm_c.n; i++) acell.cmm_c.m[i].set_array(6);
		acell.ixc.set_array(Nc); acell.iyc.set_array(Nc); acell.izc.set_array(Nc); 


		// following parameters / arrays are to be used in GPU

		/*
		int Na; // number of atoms
		ARRAY<int> aID; // type of atoms -- relating to L-J or other short-range interaction, and the charge, polarization ...
		                  // dimension -- Na
		ARRAY<int> cIndx; // which cluster it belongs, the cluster is considered as a solid structure
		                    // dimension -- Na
	
		int Nc; // number of cluster
		int Nca; // maximum number of atom in each cluster
		ARRAY<int> cluster; // the definition of cluster, dimension -- (Nca + 1) * Nc
	     
	
		ARRAY<double> r; // positions of atoms, with dimension Na * 3
		ARRAY<double> mu; // dipole of atom, with dimension Na * 3
		ARRAY<int> ir; // the cubic index in CMM, with dimension Na * 3
		*/
		acell.Na = Na; acell.Nc = acell.scluster.n; acell.Nca = Nca;
		acell.aID.set_array(Na); acell.cIndx.set_array(Na); 
		acell.cluster.set_array(acell.scluster.n * (Nca + 1));
		acell.r.set_array(Na * 3);
		acell.ir.set_array(Na * 3);
/*
		int NC; // maximum number of atoms in each cell, NC can not be changed when the CELL is constructed !
		int Nx, Ny, Nz; // dimension of CMM cell
		double cx[3]; // starting position of CMM, most corner position, [0][0][0]
		double cw[3]; // the width of the cell
		ARRAY<int> cmm; // the atoms inside each unit cell in CMM, the atoms are indexed 0, 1 ... for each unit cell, the format is: 0 -- [real number of atoms], 1 ~ NC+1 -- atoms
			            // so the whole dimension is (NC+1) * Nx * Ny * Nz
*/
		acell.Nx = Nx; acell.Ny = Ny; acell.Nz = Nz;
		acell.NC = NC;
		acell.cmm.set_array((NC + 1) * Nx * Ny * Nz);
		acell.Nyz = Ny * Nz;

/*
		double sr_cut, ewr_cut; // the truncation distance for short-range interaction, real part Ewald-Sum
		int Newd; // maximum number of atoms within truncation distance
		ARRAY<int> EwdRel; // neighbor atoms around each atom within the truncation distance
		// the interactions between the atoms in one cluster will be excluded
		// also the interaction between the atoms bounded directly defined in [bound] will be excluded
		// dimension -- (Newd + 1) * Na

		int Newd_ex; // the other pair-wise electrostatic interaction to be excluded, e.g. 1-4
		ARRAY<int> EwdRel_exclude; // dimension 3 * Newd_ex; each bound has parameters -- a1_indx, a2_indx, c -- the fraction of the interaction to be excluded 
*/
		acell.cx[0] = cx[0]; acell.cx[1] = cx[1]; acell.cx[2] = cx[2];
		acell.cw[0] = cw[0]; acell.cw[1] = cw[1]; acell.cw[2] = cw[2];

		// the exclusive 1-4 pairwise interactions in MMOLECULE
		int nb14 = 0, j, n, nb = 0, ib;
		BASIC_CLUSTER *c = NULL, *parent = NULL;
		MATOM *patom1, *patom2, *patom3, *patom4;
		HINGE< _BASIC_CLUSTER<MATOM> > *hinge;
		int ic, ia1, ia4, im;

		for (ic = 0; ic < acell.scluster.n; ic++) {
			if (acell.scluster.m[ic].m->mType == 0) {// this cluster is within a MMOLECULE
				c = acell.scluster.m[ic].m;
				im = c->mIndx; // mIndx start from MMOLECULE
				if (c->n2Parent >= 0 && acell.scluster.get_indxAtom((BASIC_CLUSTER*)c->parent) != NULL) { 
					hinge = c->hinge + c->n2Parent;
					patom2 = c->atom + hinge->indx;
					patom3 = c->parent->atom + hinge->indx_link;
					nb14 += (patom2->nb - 1) * (patom3->nb - 1);
					nb += 1; // only hinge-bound is considered, the bound inside a cluster will be considered independently, 1-2 bound
					// we have to consider the 1-3 bounds also because of the hinge
					nb += (patom2->nb - 1) + (patom3->nb - 1); // 2's neighbors with 3, 2 with 3's neighbors, 1-3 bound
				}
			}
		}
/*
		int Nb; // number of bounds
		ARRAY<int> bound; // the bounds [atom1_index, atom2_index]
			              // so, its dimension -- 2 * Nb
*/
		acell.Nb = nb;
		acell.bound.set_array(nb * 2);

		acell.Nb14 = nb14;
		acell.bound_14.set_array(nb14 * 2);
		ib = 0; n = 0;
		for (ic = 0; ic < acell.scluster.n; ic++) {
			if (acell.scluster.m[ic].m->mType == 0) {// this cluster is within a MMOLECULE
				c = acell.scluster.m[ic].m;
				im = c->mIndx; // mIndx start from MMOLECULE
				parent = (BASIC_CLUSTER*)c->parent;
				if (c->n2Parent >= 0 && acell.scluster.get_indxAtom(parent) != NULL) { 
					hinge = c->hinge + c->n2Parent;
					patom2 = c->atom + hinge->indx;
					patom3 = c->parent->atom + hinge->indx_link;
					// 1-2 bound
					acell.bound.m[ib] = acell.satom.get_indxAtom(patom2)->indx;
					acell.bound.m[ib + 1] = acell.satom.get_indxAtom(patom3)->indx;
					ib += 2;
					// 1-3 bound
					for (ia1 = 0; ia1 < patom2->nb; ia1++) {
						if (patom2->ba[ia1] == patom3) continue;
						acell.bound.m[ib] = acell.satom.get_indxAtom(patom2->ba[ia1])->indx;
						acell.bound.m[ib + 1] = acell.satom.get_indxAtom(patom3)->indx;
						ib += 2;
					}
					for (ia1 = 0; ia1 < patom3->nb; ia1++) {
						if (patom3->ba[ia1] == patom2) continue;
						acell.bound.m[ib] = acell.satom.get_indxAtom(patom3->ba[ia1])->indx;
						acell.bound.m[ib + 1] = acell.satom.get_indxAtom(patom2)->indx;
						ib += 2;
					}
					// 1-4 bound
					for (ia1 = 0; ia1 < patom2->nb; ia1++) {
						patom1 = patom2->ba[ia1];
						if (patom1 == patom3) continue;
						for (ia4 = 0; ia4 < patom3->nb; ia4++) {
							patom4 = patom3->ba[ia4];
							if (patom4 == patom2) continue;
							acell.bound_14.m[n] = acell.satom.get_indxAtom(patom1)->indx;
							acell.bound_14.m[n + 1] = acell.satom.get_indxAtom(patom4)->indx;
							n += 2;
						}
					}
				}
			}
		}
		
		//b_relship_t and b14_relship_t, only those has bound / bound14
		acell.Na_b = 0; acell.Na_b14 = 0;
		acell.Nm_b = 0; acell.Nm_b14 = 0;
		acell.b_relship_t.set_array(Na); acell.b14_relship_t.set_array(Na);
		for (i = 0; i < Na; i++) {
			n = 0;
			for (j = 0; j < acell.Nb; j++) {
				if (acell.bound.m[j*2] == i || acell.bound.m[j*2+1] == i) n++;
			}
			acell.b_relship_t.m[i].set_array(n);
			if (n > 0) {
				n = 0;
				for (j = 0; j < acell.Nb; j++) {
					if (acell.bound.m[j*2] == i || acell.bound.m[j*2+1] == i) {acell.b_relship_t.m[i].m[n] = j; n++;}
				}
			}

			n = 0;
			for (j = 0; j < acell.Nb14; j++) {
				if (acell.bound_14.m[j*2] == i || acell.bound_14.m[j*2+1] == i) n++;
			}
			acell.b14_relship_t.m[i].set_array(n);
			if (n > 0) {
				n = 0;
				for (j = 0; j < acell.Nb14; j++) {
					if (acell.bound_14.m[j*2] == i || acell.bound_14.m[j*2+1] == i) {acell.b14_relship_t.m[i].m[n] = j;  n++;}
				}
			}
		}
		// Nm_b
		for (i = 0; i < acell.b_relship_t.n; i++) {
			if (acell.b_relship_t.m[i].n > 0) acell.Na_b += 1;
			if (acell.b_relship_t.m[i].n > acell.Nm_b) acell.Nm_b = acell.b_relship_t.m[i].n;
		}

		// all 1-2, 1-3 bound for each atom, and the 1-4 bound for each atom
		acell.abound.set_array((acell.Nm_b + 1) * acell.Na);
		memset(acell.abound.m, 0, acell.abound.n * sizeof(int));

		if (acell.Na_b > 0) {
			acell.b_relship.set_array(acell.Na_b * (acell.Nm_b + 2));
			n = 0;
			for (i = 0; i < acell.b_relship_t.n; i++) {
				if (acell.b_relship_t.m[i].n > 0) {
					acell.b_relship.m[n * (acell.Nm_b + 2)] = i; // atom index
					acell.b_relship.m[n * (acell.Nm_b + 2) + 1] = acell.b_relship_t.m[i].n; // number of bound
					memcpy(acell.b_relship.m + n * (acell.Nm_b + 2) + 2, acell.b_relship_t.m[i].m, acell.b_relship_t.m[i].n * sizeof(int)); // the list of bounds

					// list all bound atoms
					for (j = 0; j < acell.b_relship_t.m[i].n; j++) {
						acell.abound.m[i * (acell.Nm_b + 1)] += 1;
						ib = acell.b_relship_t.m[i].m[j];
						if (acell.bound.m[ib * 2] == i) acell.abound.m[i * (acell.Nm_b + 1) + j] = acell.bound.m[ib * 2 + 1];
						else if (acell.bound.m[ib * 2 + 1] == i) acell.abound.m[i * (acell.Nm_b + 1) + j] = acell.bound.m[ib * 2];
						else { // for debug purpose
							show_msg("Error: Code is not consistent for bound-checking !", true);
						}
					}

					n++;
				}
			}
		}

		//Nm_b14
		for (i = 0; i < acell.b14_relship_t.n; i++) {
			if (acell.b14_relship_t.m[i].n > 0) acell.Na_b14 += 1;
			if (acell.b14_relship_t.m[i].n > acell.Nm_b14) acell.Nm_b14 = acell.b14_relship_t.m[i].n;
		}
		if (acell.Na_b14 > 0) {
			acell.b14_relship.set_array(acell.Na_b14 * (acell.Nm_b14 + 2));
			n = 0;
			for (i = 0; i < acell.b14_relship_t.n; i++) {
				if (acell.b14_relship_t.m[i].n > 0) {
					acell.b14_relship.m[n * (acell.Nm_b14 + 2)] = i; // atom index
					acell.b14_relship.m[n * (acell.Nm_b14 + 2) + 1] = acell.b14_relship_t.m[i].n; // number of bound
					memcpy(acell.b14_relship.m + n * (acell.Nm_b14 + 2) + 2, acell.b14_relship_t.m[i].m, acell.b14_relship_t.m[i].n * sizeof(int)); // the list of bounds

					n++;
				}
			}
		}


		n = 0;
		for (i = 0; i < acell.satom.n; i++) acell.aID.m[i] = acell.satom.m[i].m->aindx;
		for (i = 0; i < acell.scluster.n; i++) {
			acell.cluster.m[i * (Nca + 1)] = acell.scluster.m[i].m->nAtoms;
			for (j = 0; j < acell.scluster.m[i].m->nAtoms; j++) {
				acell.cIndx.m[n] = i; // i == acell.scluster.m[i].indx
				acell.cluster.m[i * (Nca + 1) + j + 1] = n; // assumming that satom is indexed with the order scluster, as that shown in init_EAtCell or init_SRAtCell
				n++;
			}
		}

		// 10 parameters
		acell.ivar.set_array(10);
		acell.ivar.m[0] = acell.Na;
		acell.ivar.m[1] = acell.Nb;
		acell.ivar.m[2] = acell.Nb14;
		acell.ivar.m[3] = acell.NC;
		acell.ivar.m[4] = acell.Nc;
		acell.ivar.m[5] = acell.Nca;
		acell.ivar.m[6] = acell.Nx;
		acell.ivar.m[7] = acell.Ny;
		acell.ivar.m[8] = acell.Nz;
		acell.ivar.m[9] = acell.Nyz;

		// 9 double parameters
		acell.dvar.set_array(9);
		memcpy(acell.dvar.m, acell.cx, SIZE_V3);
		memcpy(acell.dvar.m + 3, acell.cw, SIZE_V3);
		acell.dvar.m[6] = acell.c14;
		acell.dvar.m[7] = acell.rmax_c;
		acell.dvar.m[8] = acell.rcut;
	};

	void init_EAtCell_ResVars(EAtCELL &acell) {
		//Results:
		/*
		ARRAY<double> Ex, Ey, Ez; // the electrostatic field of each atoms, dimension -- Na
		ARRAY<double> Fx_e, Fy_e, Fz_e; // force from electrostatic interaction, dimension -- Na
		*/
		int Na = acell.Na, nb = acell.Nb;
		acell.c.set_array(Na); acell.mu.set_array(Na * 3);

		acell.Fe.set_array(Na * 3);
		acell.E.set_array(Na * 3);
		acell.ve.set_array(Na * 3);
		acell.U_estat.set_array(Na);

		// temperary result:
		/*
		ARRAY<double> Ex_ex, Ey_ex, Ez_Ex; // electrostatic field from the exclusion part, e.g. 1-4, actually the field is on the first atom defined in the EwdRel_exclude, a1_indx of each pair
		                                   // the field on another atom, a2_indx, should be negative values, or inverse direction
			                               // dimension -- Newd_ex
		ARRAY<double> Fx_ex, Fy_ex, Fz_Ex; // electrostatic force from the exclusion part, e.g. 1-4, actually the field is on the first atom defined in the EwdRel_exclude, a1_indx of each pair
		                                   // the field on another atom, a2_indx, should be negative values, or inverse direction
			                               // dimension -- Newd_ex
	   */
		acell.E_14.set_array(acell.Nb14 * 3);
		acell.Fe_14.set_array(acell.Nb14 * 3);
		acell.ve_14.set_array(acell.Nb14 * 3);
		acell.U14_estat.set_array(acell.Nb14);

		acell.E_b.set_array(nb * 3);
		acell.Fe_b.set_array(nb * 3);
		acell.ve_b.set_array(nb * 3);
		acell.Ub_estat.set_array(nb);
	}

	void init_SRAtCell_ResVars(SR_AtCELL &acell) {
		int Na = acell.Na, nb = acell.Nb;
		acell.bShocked.set_array(Na);
		//ARRAY<double> Fx_sr, Fy_sr, Fz_sr; // the force from short-range interaction for each atoms, dimension -- Na
		acell.Fsr.set_array(Na * 3);
		acell.vsr.set_array(Na * 3);
		acell.U_sr.set_array(Na);

		acell.Fsr_b.set_array(nb * 3);
		acell.vsr_b.set_array(nb * 3);
		acell.Ub_sr.set_array(nb);
	}
	
	void init_EAtCell(MMOL_MD_CELL &mdcell, EAtCELL &cell, int Nx, int Ny, int Nz, float *cx, float *cw) {
#ifdef __GPU__
		EAtCell_use_gpu(cell);
#endif

		cell.scluster.set_array(mdcell.ccluster.n); cell.satom.set_array(mdcell.catom.n);
		int ia, ic;
		for (ic = 0; ic < mdcell.ccluster.n; ic++) {cell.scluster.m[ic].m = mdcell.ccluster.m[ic]; cell.scluster.m[ic].indx = ic;}
		for (ia = 0; ia < mdcell.catom.n; ia++) {cell.satom.m[ia].m = mdcell.catom.m[ia]; cell.satom.m[ia].indx = ia;}

		int NC = int(cw[0] * cw[1] * cw[2]) + 10;
		init_AtCell((*(_AtCELL*)&cell), NC, Nx, Ny, Nz, cx, cw);
		init_EAtCell_ResVars(cell);

		if (!GPU_init_EAtCell(&cell)) {
			show_log("failure to setup cell on GPU for Ewald-Sum real part electrostatic interaction", true);
		}

#ifdef __GPU__
		if (!cell.init_gpu_cell()) {
			show_msg("failure to create structure on GPU for electrostatic short-range interaction calculation.", true);
		}
#endif
	};

	void init_SRAtCell(MMOL_MD_CELL &mdcell, SR_AtCELL &cell, int Nx, int Ny, int Nz, float *cx, float *cw) {
#ifdef __GPU__
		SRAtCell_use_gpu(cell);
#endif

		cell.scluster.set_array(mdcell.bcluster.n); cell.satom.set_array(mdcell.atom.n);
		int ia, ic;
		for (ic = 0; ic < mdcell.bcluster.n; ic++) {cell.scluster.m[ic].m = mdcell.bcluster.m[ic]; cell.scluster.m[ic].indx = ic;}
		for (ia = 0; ia < mdcell.atom.n; ia++) {cell.satom.m[ia].m = mdcell.atom.m[ia]; cell.satom.m[ia].indx = ia;}

		int NC = int(cw[0] * cw[1] * cw[2]) + 10;
		init_AtCell((*(_AtCELL*)&cell), NC, Nx, Ny, Nz, cx, cw);
		init_SRAtCell_ResVars(cell);

		cell.sr_type.set_array(cell.Na);
		if (::iFormat_LJ == 0 || ::iFormat_LJ == 1) cell.n_srpars = 3;
		else {
			show_log("GPU warnning: need to update the code for unknown short-range interaction parameters", true);
			cell.n_srpars = 3; //?
		}
		cell.sr_par.set_array(cell.Na * cell.n_srpars);

		ATOM_COMM_PARS *apar = NULL;
		for (ia = 0; ia < cell.Na; ia++) {
			cell.sr_type.m[ia] = 0x00;
			if (mdcell.atom.m[ia]->par->bLJ) {
				cell.sr_type.m[ia] += 0x01; // last bit
				apar = atompar_db.search_par_indx(cell.aID.m[ia]);
				if (apar == NULL) {
					memset(cell.sr_par.m + ia * cell.n_srpars, 0, cell.n_srpars * sizeof(double));
					show_log("GPU warnning: need to update the code for unknown short-range interaction parameters", true);
				}
				else {
					cell.sr_par.m[ia * cell.n_srpars] = apar->epsLJ;
					cell.sr_par.m[ia * cell.n_srpars + 1] = apar->rLJ;
					cell.sr_par.m[ia * cell.n_srpars + 2] = apar->r_min;
				}
			}
			if (mdcell.atom.m[ia]->par->bSR) {
				cell.sr_type.m[ia] += 0x02; // last 2nd bit
				show_msg("check the code whether given SR is implemented !", true);
			}
		}

		if (!GPU_init_srAtCell(&cell)) {
			show_log("failure to setup cell on GPU for local short-range interaction", true);
		}

#ifdef __GPU__
		if (!cell.init_gpu_cell()) {
			show_msg("failure to create structure on GPU for local short-range interaction calculation.", true);
		}
#endif
	};

	void _AtCELL::distribute_cluster() {
		int ic, N, i;
		for (i = 0; i < cmm_c.n; i++) cmm_c.m[i].release();
		for (ic = 0; ic < Nc; ic++) {
			N = ixc.m[ic] * Nyz + iyc.m[ic] * Nz + izc.m[ic];
			cmm_c.m[N].attach(ic);
		}
	}

	void cmm_index_cluster(JobIndx *job, _AtCELL *cell) {
		int i;
		int ix, iy, iz;
		BASIC_CLUSTER *c;
		//MATOM *atom;
		for (i = job->n1; i <= job->n2; i++) {
			c = cell->scluster.m[i].m;
			ix = int((c->vp->v[0] - cell->cx[0]) / cell->cw[0]); TRUNCATE_RANGE(ix, 0, (cell->Nx - 1))
			iy = int((c->vp->v[1] - cell->cx[1]) / cell->cw[1]); TRUNCATE_RANGE(iy, 0, (cell->Ny - 1))
			iz = int((c->vp->v[2] - cell->cx[2]) / cell->cw[2]); TRUNCATE_RANGE(iz, 0, (cell->Nz - 1))
			cell->ixc.m[i] = ix; cell->iyc.m[i] = iy; cell->izc.m[i] = iz;
		}
	}

	void cmm_index_atom(JobIndx *job, _AtCELL *acell) {
		int i;
		int ix, iy, iz;
		MATOM *atom;
		for (i = job->n1; i <= job->n2; i++) {
			atom = acell->satom.m[i].m;
			ix = int((atom->rg.v[0] - acell->cx[0]) / acell->cw[0]); TRUNCATE_RANGE(ix, 0, (acell->Nx - 1))
			iy = int((atom->rg.v[1] - acell->cx[1]) / acell->cw[1]); TRUNCATE_RANGE(iy, 0, (acell->Ny - 1))
			iz = int((atom->rg.v[2] - acell->cx[2]) / acell->cw[2]); TRUNCATE_RANGE(iz, 0, (acell->Nz - 1))
			acell->ir.m[i * 3] = ix; acell->ir.m[i * 3 + 1] = iy; acell->ir.m[i * 3 + 2] = iz;
		}
	}

	void _distribute_atom_cpu(JobIndx *job, _AtCELL *acell) {
		char msg[256] = "\0";
	
		int ix, iy, iz, ic, ia;
		int idx, idy, idz;
		int nx, ny, nz;
		int N, Nm;
		int Ncmm, nCm;
		IndxAtom<BASIC_CLUSTER> *c;
		//IndxAtom<MATOM> *atom;
		int aIndx, am;

		int nx_max = int(acell->rmax_c / acell->cw[0]) + 1;
		int ny_max = int(acell->rmax_c / acell->cw[1]) + 1;
		int nz_max = int(acell->rmax_c / acell->cw[2]) + 1;

		for (ix = job->n1; ix <= job->n2; ix++) {
			for (iy = 0; iy < acell->Ny; iy++) {
				for (iz = 0; iz < acell->Nz; iz++) {
					N = ix * acell->Nyz + iy * acell->Nz + iz; // index of the memory of the cell
					Nm = N * (acell->NC + 1);
					acell->cmm.m[Nm] = 0; // reset the amount of atoms inside this unit cell

					for (idx = -nx_max; idx <= nx_max; idx++) {
						nx = ix + idx; if (nx < 0) continue; else if (nx >= acell->Nx) continue;
						for (idy = -ny_max; idy <= ny_max; idy++) {
							ny = iy + idy; if (ny < 0) continue; else if (ny >= acell->Ny) continue;
							for (idz = -nz_max; idz <= nz_max; idz++) {
								nz = iz + idz; if (nz < 0) continue; else if (nz >= acell->Nz) continue;
								Ncmm = nx * acell->Nyz + ny * acell->Nz + nz;
								for (ic = 0; ic < acell->cmm_c.m[Ncmm].n; ic++) {
									c = acell->scluster.m + acell->cmm_c.m[Ncmm].m[ic];
									//for (ia = 0; ia < c->m->nAtoms; ia++) {
										//aIndx = acell->satom.get_indxAtom(c->m->atom + ia)->indx;
									nCm = (acell->Nca + 1) * c->indx;
									for (ia = 0; ia < acell->cluster.m[nCm]; ia++) {
										aIndx = acell->cluster.m[nCm + ia + 1];

										if (acell->ir.m[aIndx * 3] == ix && acell->ir.m[aIndx * 3 + 1] == iy && acell->ir.m[aIndx * 3 + 2] == iz) { // the atom's index match this cell
											am = acell->cmm.m[Nm]; // number of atom inside this cell
											if (am < acell->NC - 1) {
												acell->cmm.m[Nm + am + 1] = aIndx;
												acell->cmm.m[Nm]++;
											}
											else {
												sprintf(msg, "ERROR: CMM [%d, %d, %d] is overflowed of %d", ix, iy, iz, acell->NC); show_log(msg, true);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}


	void _update_atom(JobIndx *job, _AtCELL *cell) {
		int ia;
		IndxAtom<MATOM> *atom = cell->satom.m;
		for (ia = job->n1; ia <= job->n2; ia++) {
			memcpy(cell->r.m + 3 * ia, atom[ia].m->rg.v, SIZE_V3);
			//cell->r.m[ia * 3] = atom[ia].m->rg.v[0];
			//cell->r.m[ia * 3 + 1] = atom[ia].m->rg.v[1];
			//cell->r.m[ia * 3 + 2] = atom[ia].m->rg.v[2];
		}
	}

	void _update_atomic_charge_dipole(JobIndx *job, EAtCELL *cell) {
		int ia;
		IndxAtom<MATOM> *atom = cell->satom.m;
		for (ia = job->n1; ia <= job->n2; ia++) {
			cell->c.m[ia] = atom[ia].m->c;
			memcpy(cell->mu.m + ia * 3, atom[ia].m->mu.v, SIZE_V3);
			//cell->mu.m[ia * 3] = atom[ia].m->mu.v[0];
			//cell->mu.m[ia * 3 + 1] = atom[ia].m->mu.v[1];
			//cell->mu.m[ia * 3 + 1] = atom[ia].m->mu.v[2];
		}
	}

	void _update_atoms(_AtCELL *cell) {
		JobIndx job[MAX_THREADS];
		for (int i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
		int nThreads = -1; assignJobIndx1(job, nThreads, 0, cell->satom.n - 1, 400);
		MTOperate1<JobIndx, _AtCELL>((void*)(&_update_atom), job, nThreads, cell, true);
	}

	void update_atomic_charge_dipole(EAtCELL *cell) {
		JobIndx job[MAX_THREADS];
		for (int i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);
		int nThreads = -1; assignJobIndx1(job, nThreads, 0, cell->satom.n - 1, 400);
		MTOperate1<JobIndx, EAtCELL>((void*)(&_update_atomic_charge_dipole), job, nThreads, cell, true);
#ifdef __GPU__
		if (!cell->bUseGPU) return;
		if (cell->gpu_cell_dev == NULL) cell->init_gpu_cell();
		_gpu_md_::update_charge_dipole(&(cell->gpu_cell), cell->gpu_cell_dev);
		//cudaMemcpy(cell->gpu_cell.mu_dev, cell->gpu_cell.mu_map, cell->mu.n * sizeof(double), cudaMemcpyDeviceToDevice);
		//cudaMemcpy(cell->gpu_cell.c_dev, cell->gpu_cell.c_map, cell->c.n * sizeof(double), cudaMemcpyDeviceToDevice);
#endif
	}


	bool GPU_init_srAtCell(SR_AtCELL *acell) {
#ifdef __GPU__
		return _gpu_md_::init_GPU_SR_AtCELL(acell, &acell->gpu_cell);
#else
		return true;
#endif
	}

	bool GPU_init_EAtCell(EAtCELL *acell) {
#ifdef __GPU__
		return _gpu_md_::init_GPU_EAtCELL(acell, &acell->gpu_cell);
#else
		return true;
#endif
	}


	void _distribute_atoms_cpu(_AtCELL *acell) {
		int i;
		JobIndx job[MAX_THREADS];
		for (i = 0; i < MAX_THREADS; i++) job[i].set_pars(i, job + i);

		int nThreads = -1; assignJobIndx1(job, nThreads, 0, acell->satom.n - 1, 400);
		MTOperate1<JobIndx, _AtCELL>((void*)(&_update_atom), job, nThreads, acell, true);

		nThreads = -1; assignJobIndx1(job, nThreads, 0, acell->scluster.n - 1, 400);
		MTOperate1<JobIndx, _AtCELL>((void*)(&cmm_index_cluster), job, nThreads, acell, true);

		acell->distribute_cluster();

		nThreads = -1; assignJobIndx1(job, nThreads, 0, acell->satom.n - 1, 400);
		MTOperate1<JobIndx, _AtCELL>((void*)(&cmm_index_atom), job, nThreads, acell, true);

		nThreads = -1; assignJobIndx1(job, nThreads, 0, acell->Nx - 1, 1);
		MTOperate1<JobIndx, _AtCELL>((void*)(&_distribute_atom_cpu), job, nThreads, acell, true);
	}

	void update_atoms(EAtCELL *acell) { // only the coordination
		_update_atoms((_AtCELL*)acell);
#ifdef __GPU__
		if (acell->gpu_cell_dev == NULL) acell->init_gpu_cell();
		if (!acell->bCalUseGPU) { // although GPU exist, interaction is calculated with CPU, so we have to distribute atoms in the memory of CPU
			_distribute_atoms_cpu((_AtCELL*)acell);
		}
		else if (acell->gpu_cell_dev == NULL) {
			show_log("Error report: GPU is not initialized for short-range electrostatic interaction, check the code update_atoms(EAtCELL*)", true);
		}
		else if (!_gpu_md_::relocate_atom((_gpu_md_::GPU_AtCELL*)(&(acell->gpu_cell)), (_gpu_md_::GPU_AtCELL*)(acell->gpu_cell_dev), false)) { // not necessary to update the arrya on host
			show_log("cubic of GPU-EAtCell overfolowed!", true);
		}
		//cudaMemcpy(acell->gpu_cell.r_dev, acell->gpu_cell.r_map, acell->r.n * sizeof(double), cudaMemcpyDeviceToDevice);
		//cudaMemcpy(acell->gpu_cell.ir_dev, acell->gpu_cell.ir_map, acell->ir.n * sizeof(int), cudaMemcpyDeviceToDevice);
		//cudaMemcpy(acell->gpu_cell.cmm_dev, acell->gpu_cell.cmm_map, acell->cmm.n * sizeof(int), cudaMemcpyDeviceToDevice);
#endif
	};

	void update_atoms(SR_AtCELL *acell) { // only the coordination
		_update_atoms((_AtCELL*)acell);
#ifdef __GPU__
		if (acell->gpu_cell_dev == NULL) acell->init_gpu_cell();
		if (!acell->bCalUseGPU) { // although GPU exist, interaction is calculated with CPU, so we have to distribute atoms in the memory of CPU
			_distribute_atoms_cpu((_AtCELL*)acell);
		}
		// use r_map to update device coordination and location in cell
		else if (acell->gpu_cell_dev == NULL) {
			show_log("Error report: GPU is not initialized for SR interaction, check the code update_atoms(SR_AtCELL*)", true);
		}
		else if (!_gpu_md_::relocate_atom((_gpu_md_::GPU_AtCELL*)(&(acell->gpu_cell)), (_gpu_md_::GPU_AtCELL*)(acell->gpu_cell_dev), !acell->bUseGPU)) { // not necessary to update the arrya on host
			show_log("cubic of GPU-SR_AtCell overfolowed!", true);
		}
		//cudaMemcpy(acell->gpu_cell.r_dev, acell->gpu_cell.r_map, acell->r.n * sizeof(double), cudaMemcpyDeviceToDevice);
		//cudaMemcpy(acell->gpu_cell.ir_dev, acell->gpu_cell.ir_map, acell->ir.n * sizeof(int), cudaMemcpyDeviceToDevice);
		//cudaMemcpy(acell->gpu_cell.cmm_dev, acell->gpu_cell.cmm_map, acell->cmm.n * sizeof(int), cudaMemcpyDeviceToDevice);
#endif
	};

	void distribute_atoms(EAtCELL *acell) {
#ifndef __GPU__
		_distribute_atoms_cpu((_AtCELL*)acell);
#else
		if (acell->bCalUseGPU) {
			update_atoms(acell);
			// will be transfered by GPU
			//cudaMemcpy(acell->gpu_cell.r_dev, acell->gpu_cell.r_map, acell->r.n * sizeof(double), cudaMemcpyDeviceToDevice);
			//cudaMemcpy(acell->gpu_cell.ir_dev, acell->gpu_cell.ir_map, acell->ir.n * sizeof(int), cudaMemcpyDeviceToDevice);
			//cudaMemcpy(acell->gpu_cell.cmm_dev, acell->gpu_cell.cmm_map, acell->cmm.n * sizeof(int), cudaMemcpyDeviceToDevice);
		}
		else _distribute_atoms_cpu((_AtCELL*)acell);
#endif
	};

	void distribute_atoms(SR_AtCELL *acell) {
#ifndef __GPU__
		_distribute_atoms_cpu((_AtCELL*)acell);
#else
		if (acell->bCalUseGPU) {
			update_atoms(acell);
			// will be transfered by GPU
			//cudaMemcpy(acell->gpu_cell.r_dev, acell->gpu_cell.r_map, acell->r.n * sizeof(double), cudaMemcpyDeviceToDevice);
			//cudaMemcpy(acell->gpu_cell.ir_dev, acell->gpu_cell.ir_map, acell->ir.n * sizeof(int), cudaMemcpyDeviceToDevice);
			//cudaMemcpy(acell->gpu_cell.cmm_dev, acell->gpu_cell.cmm_map, acell->cmm.n * sizeof(int), cudaMemcpyDeviceToDevice);
		}
		else _distribute_atoms_cpu((_AtCELL*)acell);
#endif
	};



	void check_distribution(_AtCELL *acell) {
		char msg[256] = "\0";
		int na = 0, i, N;
		for (i = 0; i < acell->Nx * acell->Nyz; i++) {
			N = (acell->NC + 1) * i;
			na += acell->cmm.m[N];
		}
		sprintf(msg,"CMM cell has %d atoms, expected: %d atoms", na, acell->Na); show_log(msg, true);
	}


	void GPU_collecting_efield(JobIndx *job, EAtCELL *cell) {
		int ia, i, ib;
		for (ia = job->n1; ia <= job->n2; ia++) {
			// excluding the electric field from bounded atom
			for (i = 0; i < cell->b_relship_t.m[ia].n; i++) {
				ib = cell->b_relship_t.m[ia].m[i];
				if (cell->bound.m[b_idx1(ib)] == ia) { // the first atom
					cell->E.m[x_idx(ia)] -= cell->E_b.m[x_idx(ib)]; 
					cell->E.m[y_idx(ia)] -= cell->E_b.m[y_idx(ib)]; 
					cell->E.m[z_idx(ia)] -= cell->E_b.m[z_idx(ib)]; 
				}
				else if (cell->bound.m[b_idx2(ib)] == ia) { // the second atom, electric field is on the first atom
					cell->E.m[x_idx(ia)] += cell->E_b.m[x_idx(ib)]; 
					cell->E.m[y_idx(ia)] += cell->E_b.m[y_idx(ib)]; 
					cell->E.m[z_idx(ia)] += cell->E_b.m[z_idx(ib)]; 
				}
			}

			// excluding the electric field from 1-4 bounded atom
			for (i = 0; i < cell->b14_relship_t.m[ia].n; i++) {
				ib = cell->b14_relship_t.m[ia].m[i];
				if (cell->bound_14.m[b_idx1(ib)] == ia) { // the first atom
					cell->E.m[x_idx(ia)] -= cell->E_14.m[x_idx(ib)] * cell->c14; 
					cell->E.m[y_idx(ia)] -= cell->E_14.m[y_idx(ib)] * cell->c14; 
					cell->E.m[z_idx(ia)] -= cell->E_14.m[z_idx(ib)] * cell->c14; 
				}
				else if (cell->bound_14.m[b_idx2(ib)] == ia) { // the second atom, electric field is on the first atom
					cell->E.m[x_idx(ia)] += cell->E_14.m[x_idx(ib)] * cell->c14; 
					cell->E.m[y_idx(ia)] += cell->E_14.m[y_idx(ib)] * cell->c14; 
					cell->E.m[z_idx(ia)] += cell->E_14.m[z_idx(ib)] * cell->c14; 
				}
			}
		}
	}

	void GPU_collecting_eforce(JobIndx *job, EAtCELL *cell) {
		int ia, i, ib;
		for (ia = job->n1; ia <= job->n2; ia++) {
			if (cell->bCalE == 2 || cell->bCalE == 3) {
				// excluding the electric force from bounded atom
				for (i = 0; i < cell->b_relship_t.m[ia].n; i++) {
					ib = cell->b_relship_t.m[ia].m[i];
					if (cell->bound.m[b_idx1(ib)] == ia) { // the first atom
						cell->Fe.m[x_idx(ia)] -= cell->Fe_b.m[x_idx(ib)]; 
						cell->Fe.m[y_idx(ia)] -= cell->Fe_b.m[y_idx(ib)]; 
						cell->Fe.m[z_idx(ia)] -= cell->Fe_b.m[z_idx(ib)]; 

						// collecting energy and virial term at this case only
						cell->U_estat.m[ia] -= cell->Ub_estat.m[ib];
						cell->ve.m[x_idx(ia)] -= cell->ve_b.m[x_idx(ib)];
						cell->ve.m[y_idx(ia)] -= cell->ve_b.m[y_idx(ib)];
						cell->ve.m[z_idx(ia)] -= cell->ve_b.m[z_idx(ib)];
					}
					else if (cell->bound.m[b_idx2(ib)] == ia) { // the second atom, electric field is on the first atom
						cell->Fe.m[x_idx(ia)] += cell->Fe_b.m[x_idx(ib)]; 
						cell->Fe.m[y_idx(ia)] += cell->Fe_b.m[y_idx(ib)]; 
						cell->Fe.m[z_idx(ia)] += cell->Fe_b.m[z_idx(ib)]; 
					}
				}

				// excluding the electric force from 1-4 bounded atom
				for (i = 0; i < cell->b14_relship_t.m[ia].n; i++) {
					ib = cell->b14_relship_t.m[ia].m[i];
					if (cell->bound_14.m[b_idx1(ib)] == ia) { // the first atom
						cell->Fe.m[x_idx(ia)] -= cell->Fe_14.m[x_idx(ib)] * cell->c14; 
						cell->Fe.m[y_idx(ia)] -= cell->Fe_14.m[y_idx(ib)] * cell->c14; 
						cell->Fe.m[z_idx(ia)] -= cell->Fe_14.m[z_idx(ib)] * cell->c14; 

						// collecting energy and virial term at this case only
						cell->U_estat.m[ia] -= cell->U14_estat.m[ib] * cell->c14;
						cell->ve.m[x_idx(ia)] -= cell->ve_14.m[x_idx(ib)] * cell->c14;
						cell->ve.m[y_idx(ia)] -= cell->ve_14.m[y_idx(ib)] * cell->c14;
						cell->ve.m[z_idx(ia)] -= cell->ve_14.m[z_idx(ib)] * cell->c14;
					}
					else if (cell->bound_14.m[b_idx2(ib)] == ia) { // the second atom, electric field is on the first atom
						cell->Fe.m[x_idx(ia)] += cell->Fe_14.m[x_idx(ib)] * cell->c14; 
						cell->Fe.m[y_idx(ia)] += cell->Fe_14.m[y_idx(ib)] * cell->c14; 
						cell->Fe.m[z_idx(ia)] += cell->Fe_14.m[z_idx(ib)] * cell->c14; 
					}
				}
			}
		}
	}


	void GPU_collecting_srforce(JobIndx *job, SR_AtCELL *cell) {
		int ia, i, ib;
		for (ia = job->n1; ia <= job->n2; ia++) {
			if (cell->bCalSR == 1) {
				// excluding the sr force from bounded atom
				for (i = 0; i < cell->b_relship_t.m[ia].n; i++) {
					ib = cell->b_relship_t.m[ia].m[i];
					if (cell->bound.m[b_idx1(ib)] == ia) { // the first atom
						cell->Fsr.m[x_idx(ia)] -= cell->Fsr_b.m[x_idx(ib)]; 
						cell->Fsr.m[y_idx(ia)] -= cell->Fsr_b.m[y_idx(ib)]; 
						cell->Fsr.m[z_idx(ia)] -= cell->Fsr_b.m[z_idx(ib)]; 

						// collecting energy and virial term at this case only
						cell->U_sr.m[ia] -= cell->Ub_sr.m[ib];
						cell->vsr.m[x_idx(ia)] -= cell->vsr_b.m[x_idx(ib)];
						cell->vsr.m[y_idx(ia)] -= cell->vsr_b.m[y_idx(ib)];
						cell->vsr.m[z_idx(ia)] -= cell->vsr_b.m[z_idx(ib)];
					}
					else if (cell->bound.m[b_idx2(ib)] == ia) { // the second atom, sr force is on the first atom
						cell->Fsr.m[x_idx(ia)] += cell->Fsr_b.m[x_idx(ib)]; 
						cell->Fsr.m[y_idx(ia)] += cell->Fsr_b.m[y_idx(ib)]; 
						cell->Fsr.m[z_idx(ia)] += cell->Fsr_b.m[z_idx(ib)]; 
					}
				}

				// excluding the electric field from 1-4 bounded atom
				for (i = 0; i < cell->b14_relship_t.m[ia].n; i++) {
					ib = cell->b14_relship_t.m[ia].m[i];
					if (cell->bound_14.m[b_idx1(ib)] == ia) { // the first atom
						cell->Fsr.m[x_idx(ia)] -= cell->Fsr_14.m[x_idx(ib)] * cell->c14; 
						cell->Fsr.m[y_idx(ia)] -= cell->Fsr_14.m[y_idx(ib)] * cell->c14; 
						cell->Fsr.m[z_idx(ia)] -= cell->Fsr_14.m[z_idx(ib)] * cell->c14; 

						// collecting energy and virial term at this case only
						cell->U_sr.m[ia] -= cell->U14_sr.m[ib] * cell->c14;
						cell->vsr.m[x_idx(ia)] -= cell->vsr_14.m[x_idx(ib)] * cell->c14;
						cell->vsr.m[y_idx(ia)] -= cell->vsr_14.m[y_idx(ib)] * cell->c14;
						cell->vsr.m[z_idx(ia)] -= cell->vsr_14.m[z_idx(ib)] * cell->c14;
					}
					else if (cell->bound_14.m[b_idx2(ib)] == ia) { // the second atom, sr force is on the first atom
						cell->Fsr.m[x_idx(ia)] += cell->Fsr_14.m[x_idx(ib)] * cell->c14; 
						cell->Fsr.m[y_idx(ia)] += cell->Fsr_14.m[y_idx(ib)] * cell->c14; 
						cell->Fsr.m[z_idx(ia)] += cell->Fsr_14.m[z_idx(ib)] * cell->c14; 
					}
				}
			}
		}
	}


} // end of namespace


#ifdef __GPU__
namespace _gpu_md_ {

	bool init_GPU_AtCELL(_atomic_cmm_::_AtCELL *acell, GPU_AtCELL *gpu_acell) { // gpu_acell is a host variable, although its pointer is pointing to GPU-memory
		gpu_acell->Nc = acell->Nc;
		gpu_acell->Nca = acell->Nca;
		gpu_acell->Na = acell->Na;
		gpu_acell->Nb = acell->Nb;
		gpu_acell->Nb14 = acell->Nb14;
		gpu_acell->Na_b = acell->Na_b;
		gpu_acell->Na_b14 = acell->Na_b14;
		gpu_acell->Nm_b = acell->Nm_b;
		gpu_acell->Nm_b14 = acell->Nm_b14;
		gpu_acell->c14 = acell->c14;
		gpu_acell->NC = acell->NC;
		gpu_acell->Nx = acell->Nx;
		gpu_acell->Ny = acell->Ny;
		gpu_acell->Nz = acell->Nz;
		gpu_acell->Nyz = acell->Nyz;
		memcpy(gpu_acell->cx, acell->cx, 3 * sizeof(double));
		memcpy(gpu_acell->cw, acell->cw, 3 * sizeof(double));
		gpu_acell->rmax_c = acell->rmax_c;
		gpu_acell->rcut = acell->rcut;

	/*
	int *ivar_map, *dvar_map;
	int *ivar_dev, *dvar_dev;

	int Na; // number of atoms
	int *aID_map, *aID_dev;
	//ARRAY<int> aID; // type of atoms -- relating to L-J or other short-range interaction, and the charge, polarization ...
	                  // dimension -- Na
	int *cIndx_map, *cIndx_dev;
	//ARRAY<int> cIndx; // which cluster it belongs, the cluster is considered as a solid structure
	                    // dimension -- Na
	
	// bound and bound14 are from the hinges in macromolecules
	// bound includes 1-2 bound and 1-3 bound
	int Nb; // number of bounds
	int *bound_map, *bound_dev;
	//ARRAY<int> bound; // the bounds [atom1_index, atom2_index]
	                  // so, its dimension -- 2 * Nb
	int Nb14; // 1-4 neighbors
	int *bound_14_map, *bound_14_dev;
	//ARRAY<int> bound_14; // dimension -- 2 * Nb14; each bound has parameters -- a1_indx, a2_indx, 
	double c14;

	// relationship of bound and bound_14 with the atom needs to be setup
	// this relationship is not required for GPU calculation
	int Na_b, Na_b14; // the atoms related to bound and bound_14
	int Nm_b, Nm_b14; // maximum bnumber of bound, and maximum number of 1-4 bound, for each atom
	int *b_relship_map, *b_relship_dev;
	int *b14_relship_map, *b14_relship_dev; // they have dimension Na_b * (Nm_b+2) and Na_b14 * (Nm_b14 + 2) respectively
	                                   // bounds of each atom is saved in format: [bound_atom_indx] [number bound] [bound_index in bound / bound_14]
	// end here

	//ARRAY< ARRAY_FLEX<int> > b_relship_t, b14_relship_t; // temparory relationship for all atoms, dimension Na
	 
	double *r_dev, *r_map;
	//ARRAY<double> x, y, z; // positions of atoms, with dimension Na * 3
	int *ir_map, *ir_dev;
	//ARRAY<int> ix, iy, iz; // the cubic index in CMM, with dimension Na * 3

	int NC; // maximum number of atoms in each cell, NC can not be changed when the CELL is constructed !
	int Nx, Ny, Nz; // dimension of CMM cell
	int *cmm_map, *cmm_dev; // dimension (NC+1) * Nx * Ny * Nz
	
	*/

		gpu_acell->release();
		if (cudaMalloc(&(gpu_acell->cluster), acell->cluster.n * sizeof(int)) != cudaSuccess) {gpu_acell->cluster = NULL; gpu_acell->release(); return false;}
		if (cudaMemcpy(gpu_acell->cluster, acell->cluster.m_dev, acell->cluster.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
			show_msg("failure to copy cluster information to GPU", true);
		}

		if (cudaMalloc(&(gpu_acell->aID_dev), acell->Na * sizeof(int)) != cudaSuccess) {gpu_acell->aID_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->aID_map = acell->aID.m_dev;
		if (cudaMemcpy(gpu_acell->aID_dev, gpu_acell->aID_map, acell->aID.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
			show_msg("failure to copy atomic-id information to GPU", true);
		}

		if (cudaMalloc(&(gpu_acell->cIndx_dev), acell->Na * sizeof(int)) != cudaSuccess) {gpu_acell->cIndx_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->cIndx_map = acell->cIndx.m_dev;
		if (cudaMemcpy(gpu_acell->cIndx_dev, gpu_acell->cIndx_map, acell->cIndx.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
			show_msg("failure to copy atomic cluster-id information to GPU", true);
		}

		if (cudaMalloc(&(gpu_acell->bound_dev), acell->bound.n * sizeof(int)) != cudaSuccess) {gpu_acell->bound_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->bound_map = acell->bound.m_dev;
		if (acell->bound.n > 0) {
			if (cudaMemcpy(gpu_acell->bound_dev, gpu_acell->bound_map, acell->bound.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy bound information to GPU", true);
			}
		}
		
		if (cudaMalloc(&(gpu_acell->bound_14_dev), acell->bound_14.n * sizeof(int)) != cudaSuccess) {gpu_acell->bound_14_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->bound_14_map = acell->bound_14.m_dev;
		if (acell->bound_14.n > 0) {
			if (cudaMemcpy(gpu_acell->bound_14_dev, gpu_acell->bound_14_map, acell->bound_14.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy 1-4 bound information to GPU", true);
			}
		}

		if (cudaMalloc(&(gpu_acell->b_relship_dev), acell->b_relship.n * sizeof(int)) != cudaSuccess) {gpu_acell->b_relship_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->b_relship_map = acell->b_relship.m_dev;
		if (acell->b_relship.n > 0) {
			if (cudaMemcpy(gpu_acell->b_relship_dev, gpu_acell->b_relship_map, acell->b_relship.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy bound-relationship information to GPU", true);
			}
		}
		
		if (cudaMalloc(&(gpu_acell->b14_relship_dev), acell->b14_relship.n * sizeof(int)) != cudaSuccess) {gpu_acell->b14_relship_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->b14_relship_map = acell->b14_relship.m_dev;
		if (acell->b14_relship.n > 0) {
			if (cudaMemcpy(gpu_acell->b14_relship_dev, gpu_acell->b14_relship_map, acell->b14_relship.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy 1-4 bound-relationship information to GPU", true);
			}
		}

		if (cudaMalloc(&(gpu_acell->abound_dev), acell->abound.n * sizeof(int)) != cudaSuccess) {gpu_acell->abound_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->abound_map = acell->abound.m_dev;
		if (acell->abound.n > 0) {
			if (cudaMemcpy(gpu_acell->abound_dev, gpu_acell->abound_map, acell->abound.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy abound information to GPU", true);
			}
		}

		if (cudaMalloc(&(gpu_acell->r_dev), acell->r.n * sizeof(double)) != cudaSuccess) {gpu_acell->r_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->r_map = acell->r.m_dev;
		if (acell->r.n > 0) {
			if (cudaMemcpy(gpu_acell->r_dev, gpu_acell->r_map, acell->r.n * sizeof(double), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy atomic coordinations to GPU", true);
			}
		}

		if (cudaMalloc(&(gpu_acell->ir_dev), acell->ir.n * sizeof(int)) != cudaSuccess) {gpu_acell->ir_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->ir_map = acell->ir.m_dev;
		if (acell->ir.n > 0) {
			if (cudaMemcpy(gpu_acell->ir_dev, gpu_acell->ir_map, acell->ir.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy atomic coordinations to GPU", true);
			}
		}

		if (cudaMalloc(&(gpu_acell->cmm_dev), acell->cmm.n * sizeof(int)) != cudaSuccess) {gpu_acell->cmm_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->cmm_map = acell->cmm.m_dev;
		if (acell->cmm.n > 0) {
			if (cudaMemcpy(gpu_acell->cmm_dev, gpu_acell->cmm_map, acell->cmm.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy CMM information to GPU", true);
			}
		}

		if (cudaMalloc(&(gpu_acell->ivar_dev), acell->ivar.n * sizeof(int)) != cudaSuccess) {gpu_acell->ivar_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->ivar_map = acell->ivar.m_dev;
		if (acell->ivar.n > 0) {
			if (cudaMemcpy(gpu_acell->ivar_dev, gpu_acell->ivar_map, acell->ivar.n * sizeof(int), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy variables to GPU", true);
			}
		}

		if (cudaMalloc(&(gpu_acell->dvar_dev), acell->dvar.n * sizeof(double)) != cudaSuccess) {gpu_acell->dvar_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->dvar_map = acell->dvar.m_dev;
		if (acell->dvar.n > 0) {
			if (cudaMemcpy(gpu_acell->dvar_dev, gpu_acell->dvar_map, acell->dvar.n * sizeof(double), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy variables to GPU", true);
			}
		}

		return true;
	};

	bool init_GPU_EAtCELL(_atomic_cmm_::EAtCELL *acell, GPU_EAtCELL *gpu_acell) { // gpu_acell is a host variable, although its pointer is pointing to GPU-memory
		bool status = init_GPU_AtCELL((_atomic_cmm_::_AtCELL*)acell, (GPU_AtCELL*)gpu_acell);
		if (!status) return false;

		gpu_acell->bDipole = acell->bDipole;
		gpu_acell->bCalE = acell->bCalE;

		gpu_acell->mgvar.init_gvar(::T);

		if (cudaMalloc(&(gpu_acell->c_dev), acell->c.n * sizeof(double)) != cudaSuccess) {gpu_acell->c_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->c_map = acell->c.m_dev;
		if (acell->c.n > 0) {
			if (cudaMemcpy(gpu_acell->c_dev, gpu_acell->c_map, acell->c.n * sizeof(double), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy charge to GPU", true);
			}
		}
		if (cudaMalloc(&(gpu_acell->mu_dev), acell->mu.n * sizeof(double)) != cudaSuccess) {gpu_acell->mu_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->mu_map = acell->mu.m_dev;
		if (acell->mu.n > 0) {
			if (cudaMemcpy(gpu_acell->mu_dev, gpu_acell->mu_map, acell->mu.n * sizeof(double), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy dipole to GPU", true);
			}
		}
		if (cudaMalloc(&(gpu_acell->E_dev), acell->E.n * sizeof(double)) != cudaSuccess) {gpu_acell->E_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->E_map = acell->E.m_dev;
		if (cudaMalloc(&(gpu_acell->Fe_dev), acell->Fe.n * sizeof(double)) != cudaSuccess) {gpu_acell->Fe_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->Fe_map = acell->Fe.m_dev;
		if (cudaMalloc(&(gpu_acell->ve_dev), acell->ve.n * sizeof(double)) != cudaSuccess) {gpu_acell->ve_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->ve_map = acell->ve.m_dev;
		if (cudaMalloc(&(gpu_acell->U_estat_dev), acell->U_estat.n * sizeof(double)) != cudaSuccess) {gpu_acell->U_estat_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->U_estat_map = acell->U_estat.m_dev;

		if (cudaMalloc(&(gpu_acell->E_14_dev), acell->E_14.n * sizeof(double)) != cudaSuccess) {gpu_acell->E_14_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->E_14_map = acell->E_14.m_dev;
		if (cudaMalloc(&(gpu_acell->Fe_14_dev), acell->Fe_14.n * sizeof(double)) != cudaSuccess) {gpu_acell->Fe_14_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->Fe_14_map = acell->Fe_14.m_dev;
		if (cudaMalloc(&(gpu_acell->ve_14_dev), acell->ve_14.n * sizeof(double)) != cudaSuccess) {gpu_acell->ve_14_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->ve_14_map = acell->ve_14.m_dev;
		if (cudaMalloc(&(gpu_acell->U14_estat_dev), acell->U14_estat.n * sizeof(double)) != cudaSuccess) {gpu_acell->U14_estat_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->U14_estat_map = acell->U14_estat.m_dev;

		if (cudaMalloc(&(gpu_acell->E_b_dev), acell->E_b.n * sizeof(double)) != cudaSuccess) {gpu_acell->E_b_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->E_b_map = acell->E_b.m_dev;
		if (cudaMalloc(&(gpu_acell->Fe_b_dev), acell->Fe_b.n * sizeof(double)) != cudaSuccess) {gpu_acell->Fe_b_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->Fe_b_map = acell->Fe_b.m_dev;
		if (cudaMalloc(&(gpu_acell->ve_b_dev), acell->ve_b.n * sizeof(double)) != cudaSuccess) {gpu_acell->ve_b_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->ve_b_map = acell->ve_b.m_dev;
		if (cudaMalloc(&(gpu_acell->Ub_estat_dev), acell->Ub_estat.n * sizeof(double)) != cudaSuccess) {gpu_acell->Ub_estat_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->Ub_estat_map = acell->Ub_estat.m_dev;

		return true;
	};

	bool init_GPU_SR_AtCELL(_atomic_cmm_::SR_AtCELL *acell, GPU_SR_AtCELL *gpu_acell) { // gpu_acell is a host variable, although its pointer is pointing to GPU-memory
		bool status = init_GPU_AtCELL((_atomic_cmm_::_AtCELL*)acell, (GPU_AtCELL*)gpu_acell);
		if (!status) return false;

		gpu_acell->bCalSR = acell->bCalSR;
		gpu_acell->mgvar.init_gvar(::T);

		gpu_acell->npars = acell->n_srpars;
		gpu_acell->sr_mode = ::iFormat_LJ;

		if (cudaMalloc(&(gpu_acell->sr_type_dev), acell->sr_type.n * sizeof(char)) != cudaSuccess) {gpu_acell->sr_type_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->sr_type_map = acell->sr_type.m_dev;
		if (acell->sr_type.n > 0) {
			if (cudaMemcpy(gpu_acell->sr_type_dev, gpu_acell->sr_type_map, acell->sr_type.n * sizeof(char), cudaMemcpyDeviceToDevice) != cudaSuccess) {
				show_msg("failure to copy short-range interaction type to GPU", true);
			}
		}
		
		if (cudaMalloc(&gpu_acell->sr_par_dev, acell->sr_par.n * sizeof(double)) != cudaSuccess) {gpu_acell->sr_par_dev = NULL; gpu_acell->release(); return false;}
		if (acell->sr_par.n > 0) {
			if (cudaMemcpy(gpu_acell->sr_par_dev, acell->sr_par.m, acell->sr_par.n * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
				show_msg("failure to copy atomic short-range interaction parameters to GPU", true);
			}
		}

		if (cudaMalloc(&(gpu_acell->bShocked_dev), acell->bShocked.n * sizeof(char)) != cudaSuccess) {gpu_acell->bShocked_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->bShocked_map = acell->bShocked.m_dev;

		if (cudaMalloc(&(gpu_acell->Fsr_dev), acell->Fsr.n * sizeof(double)) != cudaSuccess) {gpu_acell->Fsr_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->Fsr_map = acell->Fsr.m_dev;
		if (cudaMalloc(&(gpu_acell->vsr_dev), acell->vsr.n * sizeof(double)) != cudaSuccess) {gpu_acell->vsr_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->vsr_map = acell->vsr.m_dev;
		if (cudaMalloc(&(gpu_acell->U_sr_dev), acell->U_sr.n * sizeof(double)) != cudaSuccess) {gpu_acell->U_sr_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->U_sr_map = acell->U_sr.m_dev;

		if (cudaMalloc(&(gpu_acell->Fsr_14_dev), acell->Fsr_14.n * sizeof(double)) != cudaSuccess) {gpu_acell->Fsr_14_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->Fsr_14_map = acell->Fsr_14.m_dev;
		if (cudaMalloc(&(gpu_acell->vsr_14_dev), acell->vsr_14.n * sizeof(double)) != cudaSuccess) {gpu_acell->vsr_14_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->vsr_14_map = acell->vsr_14.m_dev;
		if (cudaMalloc(&(gpu_acell->U14_sr_dev), acell->U14_sr.n * sizeof(double)) != cudaSuccess) {gpu_acell->U14_sr_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->U14_sr_map = acell->U14_sr.m_dev;

		if (cudaMalloc(&(gpu_acell->Fsr_b_dev), acell->Fsr_b.n * sizeof(double)) != cudaSuccess) {gpu_acell->Fsr_b_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->Fsr_b_map = acell->Fsr_b.m_dev;
		if (cudaMalloc(&(gpu_acell->vsr_b_dev), acell->vsr_b.n * sizeof(double)) != cudaSuccess) {gpu_acell->vsr_b_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->vsr_b_map = acell->vsr_b.m_dev;
		if (cudaMalloc(&(gpu_acell->Ub_sr_dev), acell->Ub_sr.n * sizeof(double)) != cudaSuccess) {gpu_acell->Ub_sr_dev = NULL; gpu_acell->release(); return false;}
		gpu_acell->Ub_sr_map = acell->Ub_sr.m_dev;

		return true;
	};

} // end of namespace _gpu_md_

#endif // end of __GPU__

namespace _atomic_cmm_ {

	// Interface constraint
	void InterfaceConstraints_AtCell_thread(_atomic_cmm_::_AtCELL &cell, int ix1, int ix2, double &Ut) {
		int ix, iy, iz, i, j, ia, ia3, nz1, nz2, ic;
		IConstraint *p_ic = NULL;

		Ut = 0;
		double U, f0, z, zf, zt;
		VECTOR3 dr, f, torq;
		int nca, *cubic;
		double *r = cell.r.m;
		BASIC_CLUSTER *pc = NULL;
		CLUSTER *pcluster = NULL;

		for (ix = ix1; ix <= ix2; ix++) {
			if (ix < 0) continue;
			else if (ix >= cell.Nx) break;
			for (iy = 0; iy < cell.Ny; iy++) { 
				for (j = 0; j < 2; j++) {
					if (_interface_constraint_::barrier[j].zf < 0) {
						nz1 = 0; nz2 = int((_interface_constraint_::barrier[j].zt - cell.cx[2]) / cell.cw[2]) + 1; if (nz2 >= cell.Nz) nz2 = cell.Nz - 1;
						zf = _interface_constraint_::barrier[j].zf; zt = _interface_constraint_::barrier[j].zt;
					}
					else {
						nz2 = cell.Nz - 1; nz1 = int((_interface_constraint_::barrier[j].zf - cell.cx[2]) / cell.cw[2]) - 1; if (nz1 < 0) nz1 = 0;
						zf = _interface_constraint_::barrier[j].zf; zt = _interface_constraint_::barrier[j].zt;
					}
						//nz1 = 0; nz2 = cell.Nz - 1;
					for (iz = nz1; iz <= nz2; iz++) {
						cubic = cell.cubic(ix, iy, iz);
						nca = cubic[0]; if (nca <= 0) continue;
						for (i = 0; i < nca; i++) {
							ia = cubic[i+1];
							if (ia < 0 || ia >= cell.satom.n) continue;
							ia3 = ia * 3; 
							z = r[ia3 + 2];
							if (z > zt || z < zf) continue;
							pc = cell.scluster.m[cell.cIndx.m[ia]].m; ic = pc->cTypeIndx;
							if (ic >= _interface_constraint_::clusterIC_db.dim_x()) continue;
							p_ic = _interface_constraint_::clusterIC_db.m[ic].m + j; // j interface
							if (p_ic->ef_prof.n <= 0) continue;
							p_ic->ef(z, U, f0);
							cell.satom.m[ia].m->F.v[2] += f0;
							Ut += U;
							/*
 							f.v[2] = f0;
							pc->dyn->fc.v[5] += f0;

							if (pc->parent != NULL) {
								pcluster = (CLUSTER*)pc;
								dr.v[0] = pcluster->cm0.v[0]; // cm0 is center of cluster relative to its Op
								dr.v[1] = pcluster->cm0.v[0];
								dr.v[2] = pcluster->cm0.v[0];
								V3PV3(dr, f, torq)
								pc->dyn->fc.v[0] += torq.v[0];
								pc->dyn->fc.v[1] += torq.v[1];
								pc->dyn->fc.v[2] += torq.v[2];
							}
							*/
						}
					}
				}
			}
		}
	}

	void InterfaceConstraints_AtCell(_atomic_cmm_::_AtCELL *cell, int nThreads, double &Ut) {
		Ut = 0;
		if (nThreads > cell->Nx) nThreads = cell->Nx;
		MRunSys3<_atomic_cmm_::_AtCELL, double>((void*)(&InterfaceConstraints_AtCell_thread), *cell, 0, cell->Nx - 1, nThreads, Ut);
	}

} // end of namespace _atomic_cmm_


#if _CHECK_TIME_STAMP_ == 1
#include "time.h"
#endif

namespace _atomic_cmm_ {
	extern void SPME_EF(_atomic_cmm_::EAtCELL &cell, _spme_2d_::VSPME<_VQatom>& vspme, int nThreads, _spme_2d_::SPME_VAR& var, InteractRes& res);
	extern void SPME_EF(_atomic_cmm_::EAtCELL &cell, _spme_2d_::VSPME<_VMUatom>& vspme, int nThreads, _spme_2d_::SPME_VAR& var, InteractRes& res);
	extern bool Polarize_SPME_Interact(MMOL_MD_CELL &mdcell, _atomic_cmm_::EAtCELL &cell, _spme_2d_::VSPME<_VMUatom>& vspme, int nThreads, bool bPolarize, _spme_2d_::SPME_VAR &var, InteractRes &res);
}

extern double EwaldSum_Recip_2d_Q(MMOL_MD_CELL &mdcell);

namespace _spme_2d_{
	extern void ExplicitEInteract_2d_U(CMM_CELL3D<BASIC_CLUSTER> &cell, _spme_2d_::E2dInteractVar& evar, InteractRes &res);
	extern void ExplicitEInteract_2d_ExcludeInducedDiploe(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, E2dInteractVar& evar, InteractRes &res);
}

void test1_SPME_EwaldSum_2d() {
#ifdef __GPU__
	cudaSetDevice(_CudaDevice_); cudaDeviceReset();
#endif
	int nThreads = MAX_THREADS;
	char msg[256] = "\0";

	using namespace _atomic_cmm_;
	using namespace _spme_2d_;

	bool bExplicit = true; // calculate interaction explicitly
	::bPolarize = true; // polarizable water molecule


	::eps_dc = 1;
	MMOL_MD_CELL mdcell;
	//int nmol = 2;
	//int nmol = 8;
	int nmol = 10;
	//int nmol = 20;
	//int nmol = 52;
	//int nmol = 104;
	//int nmol = 300;
	//int nmol = 1024;
	mdcell.SetMol(0, nmol, 0);
	CMM_CELL3D<BASIC_CLUSTER> cmm;
	::mlog.init("test.log", false);

	::rcut_Ewd = 6;
	//float h[3] = {16, 16, 16};
	//float h[3] = {8, 8, 5};
	float h[3] = {5, 5, 5};
	//float h[3] = {4, 4, 4};
	//float h[3] = {32, 32, 16};
	
	//int ncell[3] = {5, 5, 5};
	int ncell[3] = {10, 10, 10};
	float corner[3] = {-h[0], -h[1], -h[2] * 2};
	float w_cell[3] = {h[0] * 2 / ncell[0], h[1] * 2 / ncell[1], h[2] * 4 / ncell[2]};
	
	int NC = 36;

	cmm.set_cell(5, 5, 5);
	cmm.set_cell_range(true, -h[0], h[0], -h[1], h[1], -h[2], h[2]);
	mdcell.h[0] = h[0]; mdcell.h[1] = h[1]; mdcell.h[2] = h[2];  mdcell.period = true;
	cmm.set_cell_pos();
	int ncs = 10;
	CMM_array_set_storage<BASIC_CLUSTER>(&cmm, ncs);

	int imol = 0, ia, i;
	SMOLECULE cw;
	VECTOR3 dr, axis;
	axis.v[2] = 1;
	//axis.v[0] = 1;

	extern bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol);
	extern bool cp(SMOLECULE *d, SMOLECULE *s, bool setup = true);
	extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
	extern void init_LJ_Pars();
	extern void construct_LJ12_6_prof(int mode);
	extern double cal_charge(MMOL_MD_CELL &mmcell);

	extern void get_bclusters(MMOL_MD_CELL &mdcell, ARRAY<BASIC_CLUSTER*>& bcluster);

	ConstructSimpleMolecule(&cw, "SPC");
	cw.c->cal_geometry_radius();
	dr.v[0] = -cw.r.v[0]; dr.v[1] = -cw.r.v[1]; dr.v[2] = -cw.r.v[2];
	cw.shiftMol(dr);
	::iFormat_LJ = 1;  // Amber
	init_LJ_Pars();
	construct_LJ12_6_prof(::iFormat_LJ);
	for (ia = 0; ia < cw.c->nAtoms; ia++) cw.c->atom[ia].par = NULL;

	double Ep = 0;
	
	//cmm_init_subcell<BASIC_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size

	float hgap = h[2] * 0.9, dhs = h[2] - hgap;
/*
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		//dr.v[2] = (h[2] * 0.8) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (ranf() > 0.5 ? 1: -1) * (hgap + dhs * ranf());
		mdcell.sm.m[imol].shiftMol(dr);
		//rand_vect3(axis.v);
		//rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI, true);

		rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI* 0.2, true);

		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
*/	
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (hgap + dhs * ranf());
		mdcell.sm.m[imol].shiftMol(dr);
		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();

		imol++;
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[2] = -dr.v[2];
		mdcell.sm.m[imol].shiftMol(dr);
		axis.v[0] = 1; axis.v[1] = 0; axis.v[2] = 0; // along x axis
		rotate_SM(mdcell.sm.m[imol], axis, PI, true); // rotate along x axis 180 degree
		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
	get_atoms(mdcell, mdcell.atom);
	get_bclusters(mdcell, mdcell.bcluster);
	mdcell.check_eneutral();

	mdcell.set_atom_multipole((bPolarize ? 1 : 0));
	mdcell.setMolIndx();
	molecule_check_periodic_cell(mdcell);
	check_atom_rg(mdcell, nThreads);

	EAtCELL acell;
	acell.rmax_c = 1; // maximum radius of the cluster
	init_EAtCell(mdcell, acell, ncell[0], ncell[1], ncell[2], corner, w_cell);

#if _CHECK_TIME_STAMP_
	mt.start();
#endif

	distribute_atoms(&acell);
	update_atomic_charge_dipole(&acell);

#if _CHECK_TIME_STAMP_
	long dt = mt.elapse();
	sprintf(msg, "used time: %d ms", dt); show_log(msg, true);
#endif

	check_distribution((_AtCELL*)&acell);

	cmm_init_distribute_cluster(mdcell.sm.m, mdcell.sm.n, cmm, true); // add simple solid molecule into cmm
	cmm_check_cluster<BASIC_CLUSTER>(cmm);
	cmm_check<BASIC_CLUSTER>(cmm, MAX_THREADS);
	_cmm_2d_::SMCheckRelshipInCell(mdcell.sm.m, mdcell.sm.n, cmm, MAX_THREADS);

	/*
	{
		for (int imol = 0; imol < mdcell.sm.n; imol++) {
			sprintf(errmsg, "m %d [%f, %f, %f]", imol, mdcell.sm.m[imol].r0.v[0], mdcell.sm.m[imol].r0.v[1], mdcell.sm.m[imol].r0.v[2]);
			show_log(errmsg, true);
			sprintf(errmsg, "m %d has %d of neighbors", imol, number<CMM_IMAGINE<BASIC_CLUSTER> >(mdcell.sm.m[imol].c->spmeRship.ch));
			show_log(errmsg, true);
		}
	}
	*/

	// calculat the total dipole of the whole cell
	cal_dipole(mdcell, 1, mdcell.mu_surf);
	acell.mu_surf = mdcell.mu_surf; cmm.mu_surf = mdcell.mu_surf;

	acell.bCalE = 2; // 1 -- Efield, 2 -- force & U
	acell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

	// SPME parameters
	extern BSplineFunc<2> bsp;

	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 32;
	//::dim_spme[0] = 20; ::dim_spme[1] = 20; ::dim_spme[2] = 20;
	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 64;
	::dim_spme[0] = 64; ::dim_spme[1] = 64; ::dim_spme[2] = 64;
	//::dim_spme[0] = 128; ::dim_spme[1] = 128; ::dim_spme[2] = 128;

	//float rspline = 4;
	//int nSpline = int(0.5 * rspline / h[0] * dim_spme[0] + 0.5);
	//init_BSplineFunc(bsp, nSpline);
	init_BSplineFunc(bsp, 6);


	VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	vspme_q.bST = true; vspme_q.bST_diagonal = true;
	VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	vspme_mu.bST = true; vspme_mu.bST_diagonal = true;
	if (bPolarize) {
		mdcell.init_polarization_buff();
		mdcell.init_dipole_hist();
		mdcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		//vspme_mu.bsp = &(bsp);
		init_spme(&mdcell, vspme_mu); // init the viral atoms
		vspme_mu.set_BSplineFunc(&bsp);
		//vspme_mu.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);
		//vspme_mu.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true);
		vspme_mu.set_cell(-h[0], -h[1], -h[2], h[0] * 2, h[1] * 2, h[2] * 2, dim_spme[0], dim_spme[1], dim_spme[2], true);
		//vspme_mu.xl[0] = cmm.xl[0]; vspme_mu.xl[1] = cmm.xl[1]; vspme_mu.xl[2] = cmm.xl[2];
		//vspme_mu.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();

		vspme_mu.init_k();
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();
	}
	else { // charge only
		mdcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		//vspme_q.bsp = &(bsp);
		vspme_q.set_BSplineFunc(&bsp);
		init_spme(&mdcell, vspme_q); // init the viral atoms
		//vspme_q.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);
		//vspme_q.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true);
		vspme_q.set_cell(-h[0], -h[1], -h[2], h[0] * 2, h[1] * 2, h[2] * 2, dim_spme[0], dim_spme[1], dim_spme[2], true);

		//vspme_q.xl[0] = cmm.xl[0]; vspme_q.xl[1] = cmm.xl[1]; vspme_q.xl[2] = cmm.xl[2];
		//vspme_q.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();

		vspme_q.init_k();
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}

	_spme_2d_::SPME_VAR spme_var;
	spme_var.set_virial(true, true);

	spme_var.spme_3d_var.esv1.bEwald = true; 
	//spme_var.spme_3d_var.esv1.init_EwaldSum(mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8, ::rcut_Ewd);
	if (bPolarize) {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_mu.V, ::rcut_Ewd);
		spme_var.spme_3d_var.esv1.iSurface = 1; // slab
	}
	else {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_q.V, ::rcut_Ewd);
		spme_var.spme_3d_var.esv1.iSurface = 1; // slab
	}
	//spme_var.spme_3d_var.bVirial = spme_var.bVirial; spme_var.spme_3d_var.bST_diagonal = spme_var.bST_diagonal;
	spme_var.spme_3d_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.spme_3d_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	spme_var.spme_3d_var.bEfieldOnly = (acell.bCalE == 1 ? true : false);

	spme_var.spme_3d_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0); 
	spme_var.spme_3d_var.esv1.iSurfaceBoundary = ::iSurfaceBoundary; 

	double kappa = spme_var.spme_3d_var.esv1.kappa;
	spme_var.K0var.bEfieldOnly = (acell.bCalE == 1 ? true : false);
	spme_var.K0var.bVirial = spme_var.bVirial;
	spme_var.K0var.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	spme_var.K0var.set(kappa, cmm.xd[0] * cmm.xd[1]);

	acell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&spme_var.spme_3d_var.esv1));
	acell.k0var = spme_var.K0var;
	acell.rcut = ::rcut_Ewd;

#ifdef __GPU__
	if (acell.gpu_cell.srnb_dev == NULL && acell.bCalUseGPU) {
		if (!acell.check_neighbor_memory(true)) return;
		gpu_reset_esr_neighbors(&acell);
		check_neighbors(&(acell.gpu_cell), acell.gpu_cell_dev);
	}
#endif

	{
		//double Ushape = _spme_2d_::ShapeEnergy_2d(cmm);
		//sprintf(errmsg, "shape energy: %f [kT]", Ushape); show_log(errmsg, true);
	}

	//V3zero(cmm.mu_surf.mu) cmm.mu_surf.q = 0; V3zero(cmm.mu_surf.qr) V3zero(cmm.mu_surf.mu0)
	InteractRes ires;
/*	
	if (bPolarize) {
		if (!_atomic_cmm_::Polarize_SPME_Interact(mdcell, acell, vspme_mu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else _atomic_cmm_::SPME_EF(acell, vspme_q, nThreads, spme_var, ires);

	Ep = ires.U;
*/

#ifdef __GPU__

	ires.reset();
	_gpu_md_::GPU_SPME_2d spme_gpu;
	spme_gpu.atom.set_array(mdcell.atom.n); _gpu_md_::init_atom(mdcell.atom, spme_gpu.atom);
	spme_gpu.init_memory(mdcell.atom.n); 
	if (!_gpu_md_::setup_gpu_spme(&spme_gpu, &::bsp, true, cmm.xl, cmm.xd, ::dim_spme, spme_var.spme_3d_var.esv1.kappa, mdcell.atom.n, (bPolarize ? 1: 0))) {
			show_msg("GPU failure to initialize the memory for SPME", true); return;
	}
	// if the cell is expanded automatically, we need to reset the cell size in real EwaldSum
	// always do it after setup_gpu_spme(...) for 2d GPU-EwaldSum
	spme_var.spme_3d_var.esv1.init_EwaldSum(spme_gpu.vspme_host.V, ::rcut_Ewd);
	acell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&spme_var.spme_3d_var.esv1));

	size_t size_free = 0, size_total = 0;
	cudaMemGetInfo(&size_free, &size_total);
	if (size_free == 0) {
		show_msg("GPU failure to initialize the memory for SPME", true); return;
	}


	memcpy(spme_gpu.r.m, acell.r.m, acell.r.n * sizeof(double));
	//_gpu_md_::update_atomic_coordinate(&spme_gpu, (MAX_THREADS < 4 ? MAX_THREADS : 4));

	memcpy(spme_gpu.q.m, acell.c.m, acell.c.n * sizeof(double));
	_gpu_md_::_spme_::update_q_device_spme(&(spme_gpu.vspme_host));

	memcpy(spme_gpu.mu.m, acell.mu.m, acell.mu.n * sizeof(double));
	_gpu_md_::_spme_::update_dipole_device_spme(&(spme_gpu.vspme_host));

	bool bstatus = _gpu_dstruct_::_gpu_util_::check_gpu();
	if (bPolarize) {
		if (!_atomic_cmm_::Polarize_SPME_Interact(mdcell, acell, spme_gpu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else _atomic_cmm_::SPME_EF(acell, spme_gpu, nThreads, spme_var, ires);

	Ep = ires.U;

	_gpu_md_::_spme_::release_device_spme(&(spme_gpu.vspme_host), spme_gpu.vspme_dev);
	spme_gpu.vspme_dev = NULL;
	spme_gpu.status.release(); spme_gpu.E.release(); spme_gpu.F.release(); spme_gpu.mu.release(); spme_gpu.q.release();
	spme_gpu.release_plans();

#endif //#ifdef __GPU__


	sprintf(errmsg, "total dipole: [%f, %f, %f]; from q: [%f, %f, %f], from induced dipole: [%f, %f, %f]", cmm.mu_surf.mu.v[0], cmm.mu_surf.mu.v[1], cmm.mu_surf.mu.v[2], cmm.mu_surf.qr.v[0], cmm.mu_surf.qr.v[1], cmm.mu_surf.qr.v[2], cmm.mu_surf.mu_ind.v[0], cmm.mu_surf.mu_ind.v[1], cmm.mu_surf.mu_ind.v[2]); show_log(errmsg, true);

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);
	
	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-SPME-EwaldSum.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);


	{
	ofstream out;
	out.open("test-atomPos.dat");
	out<<"atom pos. : "<<endl;
	int ia;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"[ SM "<<imol<<" ] : "<<mdcell.sm.m[imol].r.v[0]<<"  "<<mdcell.sm.m[imol].r.v[1]<<"  "<<mdcell.sm.m[imol].r.v[2]<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<mdcell.sm.m[imol].c->atom[ia].r.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[2]<<endl;
		}
		out<<endl;
	}
	out.close();
	}



	{
	if (!bPolarize) {
		Ep = EwaldSum_Recip_2d_Q(mdcell);
		sprintf(errmsg, "Restrict Ewald Sum [reciprocal] : %f kT", Ep); show_log(errmsg, true); 

		init_normalized_coordinates(vspme_q, nThreads);
		cal_Q_spme(vspme_q, nThreads); // SPME Q matrix
		double U = vspme_q.cal_U_recp(); //calculate the energy from the image part 	

		sprintf(errmsg, "SPME Ewald Sum [reciprocal] : %f kT", U); show_log(errmsg, true);
	}
	}


	if (bPolarize) vspme_mu.release_plans();
	else vspme_q.release_plans();

	if (!bExplicit) return;

	// explicit calculation
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
		V3zero(mdcell.sm.m[imol].c->atom[ia].E)
		V3zero(mdcell.sm.m[imol].c->atom[ia].F)
		V3zero(mdcell.sm.m[imol].c->atom[ia].mu)
		V3zero(mdcell.sm.m[imol].c->atom[ia].ind_mu)
		for (i = 0; i < mdcell.sm.m[imol].c->atom[ia].f.n; i++) {V3zero(mdcell.sm.m[imol].c->atom[ia].f.m[i]);}
		for (i = 0; i < mdcell.mu_hist.n; i++) {V3zero(mdcell.mu_hist.m[i].mu[0]); V3zero(mdcell.mu_hist.m[i].mu[1]); V3zero(mdcell.mu_hist.m[i].mu[2]);}
	}

	E2dInteractVar evar_2d;
	evar_2d.bVirial = true; evar_2d.bST_diagonal = true; 
	evar_2d.bEfieldOnly = false; 
	evar_2d.esv1.bEwald = false; // direct interaction, not EwaldSum
	evar_2d.esv1.iSurface = 1; //slab
	evar_2d.esv1.iSurfaceBoundary = 0; // normal 
	evar_2d.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	InteractRes tres;

	cmm.mu_surf = mdcell.mu_surf;

	ires.reset();
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		tres.reset();
		_spme_2d_::ExplicitEInteract_2d(cmm, mdcell.sm.m[imol].c, evar_2d, tres);
		ires += tres;
	}
	Ep = ires.U;

	tres.reset();
	_spme_2d_::ExplicitEInteract_2d_U(cmm, evar_2d, tres);
	Ep += tres.U;

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);

	sprintf(errmsg, "explicit calculation : total energy %f [kT], with dipole-dipole interaction: %f [kT]", Ep, tres.U); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-explicit.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	
	::mlog.close();
#ifdef __GPU__
	acell.gpu_cell.release(); acell.release_device(); cudaDeviceReset();
#endif
}




void test1_SPME_EwaldSum_2d_Exclude_InducedDipole() {
#ifdef __GPU__
	cudaSetDevice(_CudaDevice_); cudaDeviceReset();
#endif

	int nThreads = MAX_THREADS;
	char msg[256] = "\0";

	using namespace _atomic_cmm_;
	using namespace _spme_2d_;

	bool bExplicit = true; // calculate interaction explicitly


	::eps_dc = 1;
	MMOL_MD_CELL mdcell;
	//int nmol = 2;
	//int nmol = 8;
	int nmol = 10;
	//int nmol = 20;
	//int nmol = 52;
	//int nmol = 104;
	//int nmol = 1024;
	mdcell.SetMol(0, nmol, 0);
	CMM_CELL3D<BASIC_CLUSTER> cmm;
	::mlog.init("test.log", false);

	::rcut_Ewd = 6;
	//float h[3] = {16, 16, 16};
	//float h[3] = {8, 8, 5};
	float h[3] = {5, 5, 5};
	//float h[3] = {4, 4, 4};
	//float h[3] = {32, 32, 16};
	
	//int ncell[3] = {5, 5, 5};
	int ncell[3] = {10, 10, 10};
	float corner[3] = {-h[0], -h[1], -h[2] * 2};
	float w_cell[3] = {h[0] * 2 / ncell[0], h[1] * 2 / ncell[1], h[2] * 4 / ncell[2]};
	
	int NC = 36;

	cmm.set_cell(5, 5, 5);
	cmm.set_cell_range(true, -h[0], h[0], -h[1], h[1], -h[2], h[2]);
	mdcell.h[0] = h[0]; mdcell.h[1] = h[1]; mdcell.h[2] = h[2];  mdcell.period = true;
	cmm.set_cell_pos();
	int ncs = 10;
	CMM_array_set_storage<BASIC_CLUSTER>(&cmm, ncs);

	int imol = 0, ia, i;
	SMOLECULE cw;
	VECTOR3 dr, axis;
	axis.v[2] = 1;
	//axis.v[0] = 1;

	extern bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol);
	extern bool cp(SMOLECULE *d, SMOLECULE *s, bool setup = true);
	extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
	extern void init_LJ_Pars();
	extern void construct_LJ12_6_prof(int mode);
	extern double cal_charge(MMOL_MD_CELL &mmcell);

	extern void get_bclusters(MMOL_MD_CELL &mdcell, ARRAY<BASIC_CLUSTER*>& bcluster);

	::bPolarize = true; // polarizable water molecule
	::bExcludeInducedDipole = true; // excluding dipole-dipole interaction

	ConstructSimpleMolecule(&cw, "SPC");
	cw.c->cal_geometry_radius();
	dr.v[0] = -cw.r.v[0]; dr.v[1] = -cw.r.v[1]; dr.v[2] = -cw.r.v[2];
	cw.shiftMol(dr);
	::iFormat_LJ = 1;  // Amber
	init_LJ_Pars();
	construct_LJ12_6_prof(::iFormat_LJ);
	for (ia = 0; ia < cw.c->nAtoms; ia++) cw.c->atom[ia].par = NULL;

	double Ep = 0;
	
	//cmm_init_subcell<BASIC_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size

	float hgap = h[2] * 0.9, dhs = h[2] - hgap;
/*
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		//dr.v[2] = (h[2] * 0.8) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (ranf() > 0.5 ? 1: -1) * (hgap + dhs * ranf());
		mdcell.sm.m[imol].shiftMol(dr);
		//rand_vect3(axis.v);
		//rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI, true);

		rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI* 0.2, true);

		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
*/	
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (hgap + dhs * ranf());
		mdcell.sm.m[imol].shiftMol(dr);
		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();

		imol++;
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[2] = -dr.v[2];
		mdcell.sm.m[imol].shiftMol(dr);
		axis.v[0] = 1; axis.v[1] = 0; axis.v[2] = 0; // along x axis
		rotate_SM(mdcell.sm.m[imol], axis, PI, true); // rotate along x axis 180 degree
		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
	get_atoms(mdcell, mdcell.atom);
	get_bclusters(mdcell, mdcell.bcluster);
	mdcell.check_eneutral();

	mdcell.set_atom_multipole((bPolarize ? 1 : 0));
	mdcell.setMolIndx();
	molecule_check_periodic_cell(mdcell);
	check_atom_rg(mdcell, nThreads);

	EAtCELL acell;
	acell.rmax_c = 1; // maximum radius of the cluster
	init_EAtCell(mdcell, acell, ncell[0], ncell[1], ncell[2], corner, w_cell);

#if _CHECK_TIME_STAMP_
	mt.start();
#endif

	distribute_atoms(&acell);
	update_atomic_charge_dipole(&acell);

#if _CHECK_TIME_STAMP_
	long dt = mt.elapse();
	sprintf(msg, "used time: %d ms", dt); show_log(msg, true);
#endif

	check_distribution((_AtCELL*)&acell);


	cmm_init_distribute_cluster(mdcell.sm.m, mdcell.sm.n, cmm, true); // add simple solid molecule into cmm
	cmm_check_cluster<BASIC_CLUSTER>(cmm);
	cmm_check<BASIC_CLUSTER>(cmm, MAX_THREADS);
	_cmm_2d_::SMCheckRelshipInCell(mdcell.sm.m, mdcell.sm.n, cmm, MAX_THREADS);

	/*
	{
		for (int imol = 0; imol < mdcell.sm.n; imol++) {
			sprintf(errmsg, "m %d [%f, %f, %f]", imol, mdcell.sm.m[imol].r0.v[0], mdcell.sm.m[imol].r0.v[1], mdcell.sm.m[imol].r0.v[2]);
			show_log(errmsg, true);
			sprintf(errmsg, "m %d has %d of neighbors", imol, number<CMM_IMAGINE<BASIC_CLUSTER> >(mdcell.sm.m[imol].c->spmeRship.ch));
			show_log(errmsg, true);
		}
	}
	*/

	// calculat the total dipole of the whole cell
	cal_dipole(mdcell, 1, mdcell.mu_surf);
	acell.mu_surf = mdcell.mu_surf; cmm.mu_surf = mdcell.mu_surf;

	acell.bCalE = 2; // 1 -- Efield, 2 -- force & U
	acell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

	// SPME parameters
	extern BSplineFunc<2> bsp;

	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 32;
	//::dim_spme[0] = 20; ::dim_spme[1] = 20; ::dim_spme[2] = 20;
	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 64;
	//::dim_spme[0] = 64; ::dim_spme[1] = 64; ::dim_spme[2] = 64;
	::dim_spme[0] = 128; ::dim_spme[1] = 128; ::dim_spme[2] = 128;

	//float rspline = 4;
	//int nSpline = int(0.5 * rspline / h[0] * dim_spme[0] + 0.5);
	//init_BSplineFunc(bsp, nSpline);
	init_BSplineFunc(bsp, 6);


	VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	vspme_q.bST = true; vspme_q.bST_diagonal = true;
	VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	vspme_mu.bST = true; vspme_mu.bST_diagonal = true;

	float q0 = 0;
	_spme_2d_::VSPME<_VMUatom> vspme_mu_induced; vspme_mu_induced.mMP = 1; // with dipole
	vspme_mu_induced.bST = true; vspme_mu_induced.bST_diagonal = true;

	if (bPolarize) {
		mdcell.init_polarization_buff();
		mdcell.init_dipole_hist();
		mdcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		//vspme_mu.bsp = &(bsp);
		init_spme(&mdcell, vspme_mu); // init the viral atoms
		vspme_mu.set_BSplineFunc(&bsp);
		//vspme_mu.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);
		//vspme_mu.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true);
		vspme_mu.set_cell(-h[0], -h[1], -h[2], h[0] * 2, h[1] * 2, h[2] * 2, dim_spme[0], dim_spme[1], dim_spme[2], true);
		//vspme_mu.xl[0] = cmm.xl[0]; vspme_mu.xl[1] = cmm.xl[1]; vspme_mu.xl[2] = cmm.xl[2];
		//vspme_mu.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();

		vspme_mu.init_k();
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();

		init_spme_induced(&mdcell, vspme_mu_induced, &(q0)); // init the viral atoms
		vspme_mu_induced.set_BSplineFunc(&bsp);
		vspme_mu_induced.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true); // important, after init();
		vspme_mu_induced.init_k(); vspme_mu_induced.init_b(); vspme_mu_induced.init_C(); vspme_mu_induced.init_vars();
	}
	else { // charge only
		mdcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		//vspme_q.bsp = &(bsp);
		vspme_q.set_BSplineFunc(&bsp);
		init_spme(&mdcell, vspme_q); // init the viral atoms
		//vspme_q.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2], false);
		//vspme_q.set_cell(cmm.xl[0], cmm.xl[1], cmm.xl[2], cmm.xd[0], cmm.xd[1], cmm.xd[2], dim_spme[0], dim_spme[1], dim_spme[2], true);
		vspme_q.set_cell(-h[0], -h[1], -h[2], h[0] * 2, h[1] * 2, h[2] * 2, dim_spme[0], dim_spme[1], dim_spme[2], true);

		//vspme_q.xl[0] = cmm.xl[0]; vspme_q.xl[1] = cmm.xl[1]; vspme_q.xl[2] = cmm.xl[2];
		//vspme_q.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();

		vspme_q.init_k();
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}

	_spme_2d_::SPME_VAR spme_var;
	spme_var.set_virial(true, true);

	spme_var.spme_3d_var.esv1.bEwald = true; 
	//spme_var.spme_3d_var.esv1.init_EwaldSum(mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8, ::rcut_Ewd);
	if (bPolarize) {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_mu.V, ::rcut_Ewd);
		spme_var.spme_3d_var.esv1.iSurface = 1; // slab
	}
	else {
		spme_var.spme_3d_var.esv1.init_EwaldSum(vspme_q.V, ::rcut_Ewd);
		spme_var.spme_3d_var.esv1.iSurface = 1; // slab
	}
	//spme_var.spme_3d_var.bVirial = spme_var.bVirial; spme_var.spme_3d_var.bST_diagonal = spme_var.bST_diagonal;
	spme_var.spme_3d_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.spme_3d_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	spme_var.spme_3d_var.bEfieldOnly = (acell.bCalE == 1 ? true : false);

	spme_var.spme_3d_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0); 
	spme_var.spme_3d_var.esv1.iSurfaceBoundary = ::iSurfaceBoundary; 

	double kappa = spme_var.spme_3d_var.esv1.kappa;
	spme_var.K0var.bEfieldOnly = (acell.bCalE == 1 ? true : false);
	spme_var.K0var.bVirial = spme_var.bVirial;
	spme_var.K0var.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	spme_var.K0var.set(kappa, cmm.xd[0] * cmm.xd[1]);

	acell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&spme_var.spme_3d_var.esv1));
	acell.k0var = spme_var.K0var;
	acell.rcut = ::rcut_Ewd;

#ifdef __GPU__
	if (acell.gpu_cell.srnb_dev == NULL && acell.bCalUseGPU) {
		if (!acell.check_neighbor_memory(true)) return;
		gpu_reset_esr_neighbors(&acell);
		check_neighbors(&(acell.gpu_cell), acell.gpu_cell_dev);
	}
#endif // __GPU__
	{
		//double Ushape = _spme_2d_::ShapeEnergy_2d(cmm);
		//sprintf(errmsg, "shape energy: %f [kT]", Ushape); show_log(errmsg, true);
	}

	//V3zero(cmm.mu_surf.mu) cmm.mu_surf.q = 0; V3zero(cmm.mu_surf.qr) V3zero(cmm.mu_surf.mu0)
	InteractRes ires;
/*
	if (bPolarize) {
		if (::bExcludeInducedDipole) {
			_atomic_cmm_::Polarize_SPME_Interact2(mdcell, acell, vspme_mu, vspme_mu_induced, nThreads, true, spme_var, ires);
		}
		else if (!_atomic_cmm_::Polarize_SPME_Interact(mdcell, acell, vspme_mu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else _atomic_cmm_::SPME_EF(acell, vspme_q, nThreads, spme_var, ires);

	Ep = ires.U;
*/


#ifdef __GPU__

	ires.reset();
	_gpu_md_::GPU_SPME_2d spme_gpu;
	spme_gpu.atom.set_array(mdcell.atom.n); _gpu_md_::init_atom(mdcell.atom, spme_gpu.atom);
	spme_gpu.init_memory(mdcell.atom.n); 
	if (!_gpu_md_::setup_gpu_spme(&spme_gpu, &::bsp, true, cmm.xl, cmm.xd, ::dim_spme, spme_var.spme_3d_var.esv1.kappa, mdcell.atom.n, (bPolarize ? 1: 0))) {
			show_msg("GPU failure to initialize the memory for SPME", true); return;
	}
	// if the cell is expanded automatically, we need to reset the cell size in real EwaldSum
	// always do it after setup_gpu_spme(...) for 2d GPU-EwaldSum
	spme_var.spme_3d_var.esv1.init_EwaldSum(spme_gpu.vspme_host.V, ::rcut_Ewd);
	acell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&spme_var.spme_3d_var.esv1));

	size_t size_free = 0, size_total = 0;
	cudaMemGetInfo(&size_free, &size_total);
	if (size_free == 0) {
		show_msg("GPU failure to initialize the memory for SPME", true); return;
	}


	memcpy(spme_gpu.r.m, acell.r.m, acell.r.n * sizeof(double));
	//_gpu_md_::update_atomic_coordinate(&spme_gpu, (MAX_THREADS < 4 ? MAX_THREADS : 4));

	memcpy(spme_gpu.q.m, acell.c.m, acell.c.n * sizeof(double));
	_gpu_md_::_spme_::update_q_device_spme(&(spme_gpu.vspme_host));

	memcpy(spme_gpu.mu.m, acell.mu.m, acell.mu.n * sizeof(double));
	_gpu_md_::_spme_::update_dipole_device_spme(&(spme_gpu.vspme_host));

	bool bstatus = _gpu_dstruct_::_gpu_util_::check_gpu();

	if (bPolarize) {
		if (::bExcludeInducedDipole) {
			_atomic_cmm_::Polarize_SPME_Interact2(mdcell, acell, vspme_mu, vspme_mu_induced, nThreads, true, spme_var, ires);
		}
		else if (!_atomic_cmm_::Polarize_SPME_Interact(mdcell, acell, spme_gpu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else _atomic_cmm_::SPME_EF(acell, spme_gpu, nThreads, spme_var, ires);

	Ep = ires.U;

	_gpu_md_::_spme_::release_device_spme(&(spme_gpu.vspme_host), spme_gpu.vspme_dev);
	spme_gpu.vspme_dev = NULL;
	spme_gpu.status.release(); spme_gpu.E.release(); spme_gpu.F.release(); spme_gpu.mu.release(); spme_gpu.q.release();
	spme_gpu.release_plans();

#endif //#ifdef __GPU__


	sprintf(errmsg, "total dipole: [%f, %f, %f]; from q: [%f, %f, %f], from induced dipole: [%f, %f, %f]", cmm.mu_surf.mu.v[0], cmm.mu_surf.mu.v[1], cmm.mu_surf.mu.v[2], cmm.mu_surf.qr.v[0], cmm.mu_surf.qr.v[1], cmm.mu_surf.qr.v[2], cmm.mu_surf.mu_ind.v[0], cmm.mu_surf.mu_ind.v[1], cmm.mu_surf.mu_ind.v[2]); show_log(errmsg, true);

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);
	
	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-SPME-EwaldSum.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);


	{
	ofstream out;
	out.open("test-atomPos.dat");
	out<<"atom pos. : "<<endl;
	int ia, i;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"[ SM "<<imol<<" ] : "<<mdcell.sm.m[imol].r.v[0]<<"  "<<mdcell.sm.m[imol].r.v[1]<<"  "<<mdcell.sm.m[imol].r.v[2]<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<mdcell.sm.m[imol].c->atom[ia].r.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[2]<<endl;
		}
		out<<endl;
	}
	out.close();
	}



	{
	if (!bPolarize) {
		Ep = EwaldSum_Recip_2d_Q(mdcell);
		sprintf(errmsg, "Restrict Ewald Sum [reciprocal] : %f kT", Ep); show_log(errmsg, true); 

		init_normalized_coordinates(vspme_q, nThreads);
		cal_Q_spme(vspme_q, nThreads); // SPME Q matrix
		double U = vspme_q.cal_U_recp(); //calculate the energy from the image part 	

		sprintf(errmsg, "SPME Ewald Sum [reciprocal] : %f kT", U); show_log(errmsg, true);
	}
	}


	if (bPolarize) vspme_mu.release_plans();
	else vspme_q.release_plans();

	if (!bExplicit) return;

	// explicit calculation
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
		V3zero(mdcell.sm.m[imol].c->atom[ia].E)
		V3zero(mdcell.sm.m[imol].c->atom[ia].F)
		V3zero(mdcell.sm.m[imol].c->atom[ia].mu)
		V3zero(mdcell.sm.m[imol].c->atom[ia].ind_mu)
		for (i = 0; i < mdcell.sm.m[imol].c->atom[ia].f.n; i++) {V3zero(mdcell.sm.m[imol].c->atom[ia].f.m[i]);}
		for (i = 0; i < mdcell.mu_hist.n; i++) {V3zero(mdcell.mu_hist.m[i].mu[0]); V3zero(mdcell.mu_hist.m[i].mu[1]); V3zero(mdcell.mu_hist.m[i].mu[2]);}
	}

	E2dInteractVar evar_2d;
	evar_2d.bVirial = true; evar_2d.bST_diagonal = true; 
	evar_2d.bEfieldOnly = false; 
	evar_2d.esv1.bEwald = false; // direct interaction, not EwaldSum
	evar_2d.esv1.iSurface = 1; //slab
	evar_2d.esv1.iSurfaceBoundary = 0; //normal
	evar_2d.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	InteractRes tres;

	cmm.mu_surf = mdcell.mu_surf;

	ires.reset();
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		tres.reset();
		if (::bExcludeInducedDipole) _spme_2d_::ExplicitEInteract_2d_ExcludeInducedDiploe(cmm, mdcell.sm.m[imol].c, evar_2d, tres);
		else _spme_2d_::ExplicitEInteract_2d(cmm, mdcell.sm.m[imol].c, evar_2d, tres);
		
		ires += tres;
	}
	Ep = ires.U;

	tres.reset();
	_spme_2d_::ExplicitEInteract_2d_U(cmm, evar_2d, tres);
	Ep += tres.U;

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);

	sprintf(errmsg, "explicit calculation : total energy %f [kT], with dipole-dipole interaction: %f [kT]", Ep, tres.U); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-explicit.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	
	::mlog.close();
#ifdef __GPU__
	acell.gpu_cell.release(); acell.release_device(); cudaDeviceReset();
#endif
}



// 3d

extern void ExplicitEInteract(CMM_CELL3D<BASIC_CLUSTER> &cell, BASIC_CLUSTER *pc, EInteractVar& evar, InteractRes &res);

void test1_SPME_EwaldSum_3d() {
#ifdef __GPU__
	if (cudaSetDevice(_CudaDevice_) != cudaSuccess) {
		cout<<"GPU is not enabled!"<<endl; return;
	}
	cudaDeviceReset();
#endif

	int nThreads = MAX_THREADS;
	char msg[256] = "\0";

	using namespace _atomic_cmm_;
	using namespace _spme_;

	bool bExplicit = false; // calculate interaction explicitly
	::bPolarize = true; // polarizable water molecule


	::eps_dc = 1;
	MMOL_MD_CELL mdcell;
	//int nmol = 2;
	//int nmol = 8;
	int nmol = 20;
	//int nmol = 20;
	//int nmol = 52;
	//int nmol = 104;
	//int nmol = 300;
	//int nmol = 1024;
	mdcell.SetMol(0, nmol, 0);
	CMM_CELL3D<BASIC_CLUSTER> cmm;
	::mlog.init("test.log", false);

	::rcut_Ewd = 6;
	//float h[3] = {16, 16, 16};
	float h[3] = {8, 8, 8};
	//float h[3] = {5, 5, 5};
	//float h[3] = {4, 4, 4};
	//float h[3] = {32, 32, 16};
	
	//int ncell[3] = {5, 5, 5};
	int ncell[3] = {10, 10, 10};
	float corner[3] = {-h[0], -h[1], -h[2]};
	float w_cell[3] = {h[0] * 2 / ncell[0], h[1] * 2 / ncell[1], h[2] * 2 / ncell[2]};
	
	int NC = 36;

	cmm.set_cell(ncell[0], ncell[1], ncell[2]);
	cmm.set_cell_range(true, -h[0], h[0], -h[1], h[1], -h[2], h[2]);
	mdcell.h[0] = h[0]; mdcell.h[1] = h[1]; mdcell.h[2] = h[2];  mdcell.period = true;
	cmm.set_cell_pos();
	int ncs = 10;
	CMM_array_set_storage<BASIC_CLUSTER>(&cmm, ncs);

	int imol = 0, ia, i;
	SMOLECULE cw;
	VECTOR3 dr, axis;
	axis.v[2] = 1;
	//axis.v[0] = 1;

	extern bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol);
	extern bool cp(SMOLECULE *d, SMOLECULE *s, bool setup = true);
	extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
	extern void init_LJ_Pars();
	extern void construct_LJ12_6_prof(int mode);
	extern double cal_charge(MMOL_MD_CELL &mmcell);

	extern void get_bclusters(MMOL_MD_CELL &mdcell, ARRAY<BASIC_CLUSTER*>& bcluster);


	ConstructSimpleMolecule(&cw, "SPC");
	cw.c->cal_geometry_radius();
	dr.v[0] = -cw.r.v[0]; dr.v[1] = -cw.r.v[1]; dr.v[2] = -cw.r.v[2];
	cw.shiftMol(dr);
	::iFormat_LJ = 1;  // Amber
	init_LJ_Pars();
	construct_LJ12_6_prof(::iFormat_LJ);
	for (ia = 0; ia < cw.c->nAtoms; ia++) cw.c->atom[ia].par = NULL;

	double Ep = 0;
	
	//cmm_init_subcell<BASIC_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (h[2]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		mdcell.sm.m[imol].shiftMol(dr);
		rand_vect3(axis.v);
		rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI, true);

		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
	get_atoms(mdcell, mdcell.atom);
	get_bclusters(mdcell, mdcell.bcluster);
	mdcell.check_eneutral();

	mdcell.set_atom_multipole((bPolarize ? 1 : 0));
	mdcell.setMolIndx();
	molecule_check_periodic_cell(mdcell);
	check_atom_rg(mdcell, nThreads);

	EAtCELL acell;
#ifdef __GPU__
	acell.bCalUseGPU = bGPURealEwaldSum;
#endif
	acell.rmax_c = 1; // maximum radius of the cluster
	init_EAtCell(mdcell, acell, ncell[0], ncell[1], ncell[2], corner, w_cell);

#if _CHECK_TIME_STAMP_
	mt.start();
#endif

	distribute_atoms(&acell);
	update_atomic_charge_dipole(&acell);

#if _CHECK_TIME_STAMP_
	long dt = mt.elapse();
	sprintf(msg, "used time: %d ms", dt); show_log(msg, true);
#endif

	check_distribution((_AtCELL*)&acell);

	cmm_init_distribute_cluster(mdcell.sm.m, mdcell.sm.n, cmm, true); // add simple solid molecule into cmm
	cmm_check_cluster<BASIC_CLUSTER>(cmm);
	cmm_check<BASIC_CLUSTER>(cmm, MAX_THREADS);
	_cmm_3d_::SMCheckRelshipInCell(mdcell.sm.m, mdcell.sm.n, cmm, MAX_THREADS);

	/*
	{
		for (int imol = 0; imol < mdcell.sm.n; imol++) {
			sprintf(errmsg, "m %d [%f, %f, %f]", imol, mdcell.sm.m[imol].r0.v[0], mdcell.sm.m[imol].r0.v[1], mdcell.sm.m[imol].r0.v[2]);
			show_log(errmsg, true);
			sprintf(errmsg, "m %d has %d of neighbors", imol, number<CMM_IMAGINE<BASIC_CLUSTER> >(mdcell.sm.m[imol].c->spmeRship.ch));
			show_log(errmsg, true);
		}
	}
	*/

	// calculat the total dipole of the whole cell
	cal_dipole(mdcell, 1, mdcell.mu_surf);
	acell.mu_surf = mdcell.mu_surf; cmm.mu_surf = mdcell.mu_surf;

	acell.bCalE = 2; // 1 -- Efield, 2 -- force & U
	acell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

	// SPME parameters
	extern BSplineFunc<2> bsp;

	//::dim_spme[0] = 8; ::dim_spme[1] = 8; ::dim_spme[2] = 8;
	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 32;
	//::dim_spme[0] = 20; ::dim_spme[1] = 20; ::dim_spme[2] = 20;
	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 64;
	::dim_spme[0] = 64; ::dim_spme[1] = 64; ::dim_spme[2] = 64;
	//::dim_spme[0] = 128; ::dim_spme[1] = 128; ::dim_spme[2] = 128;

	//float rspline = 4;
	//int nSpline = int(0.5 * rspline / h[0] * dim_spme[0] + 0.5);
	//init_BSplineFunc(bsp, nSpline);
	int nBSP = 6;
	init_BSplineFunc(bsp, nBSP);


	VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	vspme_q.bST = true; vspme_q.bST_diagonal = true;
	VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	vspme_mu.bST = true; vspme_mu.bST_diagonal = true;

	if (bPolarize) {
		mdcell.init_polarization_buff();
		mdcell.init_dipole_hist();
		mdcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		vspme_mu.bsp = &(bsp);
		init_spme(&mdcell, vspme_mu); // init the viral atoms
		vspme_mu.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_mu.xl[0] = cmm.xl[0]; vspme_mu.xl[1] = cmm.xl[1]; vspme_mu.xl[2] = cmm.xl[2];
		vspme_mu.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();
		vspme_mu.init_k();
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();
	}
	else { // charge only
		mdcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		vspme_q.bsp = &(bsp);
		init_spme(&mdcell, vspme_q); // init the viral atoms
		vspme_q.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_q.xl[0] = cmm.xl[0]; vspme_q.xl[1] = cmm.xl[1]; vspme_q.xl[2] = cmm.xl[2];
		vspme_q.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();
		vspme_q.init_k();
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}

	SPME_VAR spme_var;
	spme_var.bVirial = true; spme_var.bST_diagonal = true;

	spme_var.esv1.bEwald = true; 
	//spme_var.spme_3d_var.esv1.init_EwaldSum(mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8, ::rcut_Ewd);
	if (bPolarize) {
		spme_var.esv1.init_EwaldSum(vspme_mu.V, ::rcut_Ewd);
		spme_var.esv1.iSurface = 0; // 3d
	}
	else {
		spme_var.esv1.init_EwaldSum(vspme_q.V, ::rcut_Ewd);
		spme_var.esv1.iSurface = 0; // 3d
	}
	spme_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	spme_var.bEfieldOnly = (acell.bCalE == 1 ? true : false);

	spme_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0); 
	spme_var.esv1.iSurfaceBoundary = ::iSurfaceBoundary; 

	double kappa = spme_var.esv1.kappa;
	acell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&spme_var.esv1));
	acell.rcut = ::rcut_Ewd;

#ifdef __GPU__
	if (acell.gpu_cell.srnb_dev == NULL && acell.bCalUseGPU) {
		if (!acell.check_neighbor_memory(true)) return;
		gpu_reset_esr_neighbors(&acell);
		check_neighbors(&(acell.gpu_cell), acell.gpu_cell_dev);
	}
#endif

	{
		//double Ushape = _spme_2d_::ShapeEnergy_2d(cmm);
		//sprintf(errmsg, "shape energy: %f [kT]", Ushape); show_log(errmsg, true);
	}


	//V3zero(cmm.mu_surf.mu) cmm.mu_surf.q = 0; V3zero(cmm.mu_surf.qr) V3zero(cmm.mu_surf.mu0)
	InteractRes ires;

	if (bPolarize) {
		if (!_atomic_cmm_::Polarize_SPME_Interact(mdcell, acell, vspme_mu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else _atomic_cmm_::SPME_EF(acell, vspme_q, nThreads, spme_var, ires);

	Ep = ires.U;

	{
	sprintf(errmsg, "total dipole: [%f, %f, %f]; from q: [%f, %f, %f], from induced dipole: [%f, %f, %f]", cmm.mu_surf.mu.v[0], cmm.mu_surf.mu.v[1], cmm.mu_surf.mu.v[2], cmm.mu_surf.qr.v[0], cmm.mu_surf.qr.v[1], cmm.mu_surf.qr.v[2], cmm.mu_surf.mu_ind.v[0], cmm.mu_surf.mu_ind.v[1], cmm.mu_surf.mu_ind.v[2]); show_log(errmsg, true);

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);

	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-SPME-EwaldSum_cpu.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);
	}


	ires.reset();

#ifdef __GPU__

	show_log("", true); show_log("GPU calculation:", true);

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		V6zero(mdcell.sm.m[imol].c->dyn->fc)
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			V3zero(mdcell.sm.m[imol].c->atom[ia].E)
			V3zero(mdcell.sm.m[imol].c->atom[ia].F)
			V3zero(mdcell.sm.m[imol].c->atom[ia].mu)
			V3zero(mdcell.sm.m[imol].c->atom[ia].ind_mu)
			for (i = 0; i < mdcell.sm.m[imol].c->atom[ia].f.n; i++) {V3zero(mdcell.sm.m[imol].c->atom[ia].f.m[i]);}
			for (i = 0; i < mdcell.mu_hist.n; i++) {V3zero(mdcell.mu_hist.m[i].mu[0]); V3zero(mdcell.mu_hist.m[i].mu[1]); V3zero(mdcell.mu_hist.m[i].mu[2]);}
		}
	}
	_atomic_cmm_::update_atomic_charge_dipole(&acell); // update charge and dipole on cell

	{
	size_t size_free, size_total;
	cudaMemGetInfo(&size_free, &size_total);
	sprintf(errmsg, "device memory avaliable: %d k", int(size_free / 1024)); show_log(errmsg, true);
	if (size_free == 0) {
		show_log("No free memory on device", true); return;
	}
	}

	_gpu_md_::GPU_SPME_3d spme_gpu;
	spme_gpu.atom.set_array(mdcell.atom.n); _gpu_md_::init_atom(mdcell.atom, spme_gpu.atom);
	spme_gpu.init_memory(mdcell.atom.n); 
	if (!_gpu_md_::setup_gpu_spme(&spme_gpu, &::bsp, cmm.xl, cmm.xd, ::dim_spme, spme_var.esv1.kappa, mdcell.atom.n, (bPolarize ? 1: 0))) {
			show_msg("GPU failure to initialize the memory for SPME", true); return;
	}

	size_t size_free = 0, size_total = 0;
	cudaMemGetInfo(&size_free, &size_total);
	if (size_free == 0) {
		show_msg("GPU failure to initialize the memory for SPME", true); return;
	}


	memcpy(spme_gpu.r.m, acell.r.m, acell.r.n * sizeof(double));
	//_gpu_md_::update_atomic_coordinate(&spme_gpu, (MAX_THREADS < 4 ? MAX_THREADS : 4));

	memcpy(spme_gpu.q.m, acell.c.m, acell.c.n * sizeof(double));
	_gpu_md_::_spme_::update_q_device_spme(&(spme_gpu.vspme_host));

	memcpy(spme_gpu.mu.m, acell.mu.m, acell.mu.n * sizeof(double));
	_gpu_md_::_spme_::update_dipole_device_spme(&(spme_gpu.vspme_host));

	bool bstatus = _gpu_dstruct_::_gpu_util_::check_gpu();
	if (bPolarize) {
		if (!_atomic_cmm_::Polarize_SPME_Interact(mdcell, acell, spme_gpu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else _atomic_cmm_::SPME_EF(acell, spme_gpu, nThreads, spme_var, ires);

	Ep = ires.U;

	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);

	_gpu_md_::_spme_::release_device_spme(&(spme_gpu.vspme_host), spme_gpu.vspme_dev);
	spme_gpu.vspme_dev = NULL;
	spme_gpu.status.release(); spme_gpu.E.release(); spme_gpu.F.release(); spme_gpu.mu.release(); spme_gpu.q.release();
	spme_gpu.release_plans();

	sprintf(errmsg, "total dipole: [%f, %f, %f]; from q: [%f, %f, %f], from induced dipole: [%f, %f, %f]", cmm.mu_surf.mu.v[0], cmm.mu_surf.mu.v[1], cmm.mu_surf.mu.v[2], cmm.mu_surf.qr.v[0], cmm.mu_surf.qr.v[1], cmm.mu_surf.qr.v[2], cmm.mu_surf.mu_ind.v[0], cmm.mu_surf.mu_ind.v[1], cmm.mu_surf.mu_ind.v[2]); show_log(errmsg, true);

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);
	
	{
	ofstream out;
	out.open("test-SPME-EwaldSum_gpu.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);

#endif //#ifdef __GPU__

	{
	ofstream out;
	out.open("test-atomPos.dat");
	out<<"atom pos. : "<<endl;
	int ia;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"[ SM "<<imol<<" ] : "<<mdcell.sm.m[imol].r.v[0]<<"  "<<mdcell.sm.m[imol].r.v[1]<<"  "<<mdcell.sm.m[imol].r.v[2]<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<mdcell.sm.m[imol].c->atom[ia].r.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[2]<<endl;
		}
		out<<endl;
	}
	out.close();
	}

	if (bPolarize) vspme_mu.release_plans();
	else vspme_q.release_plans();

	if (!bExplicit) return;

	// explicit calculation
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		V6zero(mdcell.sm.m[imol].c->dyn->fc)
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			V3zero(mdcell.sm.m[imol].c->atom[ia].E)
			V3zero(mdcell.sm.m[imol].c->atom[ia].F)
			for (i = 0; i < mdcell.sm.m[imol].c->atom[ia].f.n; i++) {V3zero(mdcell.sm.m[imol].c->atom[ia].f.m[i]);}
		}
	}

	EInteractVar evar;
	evar.bVirial = true; evar.bST_diagonal = true; 
	evar.bEfieldOnly = false; 
	evar.esv1.bEwald = false; // direct interaction, not EwaldSum
	evar.esv1.iSurface = 0; // cubic
	evar.esv1.iSurfaceBoundary = 0; //normal
	evar.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0);
	InteractRes tres;

	cmm.mu_surf = mdcell.mu_surf;

	ires.reset();
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		tres.reset();
		ExplicitEInteract(cmm, mdcell.sm.m[imol].c, evar, tres);
		ires += tres;
	}
	Ep = ires.U;

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);

	sprintf(errmsg, "explicit calculation : total energy %f [kT], with dipole-dipole interaction: %f [kT]", Ep, tres.U); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-explicit.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	
	::mlog.close();

#ifdef __GPU__
	acell.gpu_cell.release(); acell.release_device(); cudaDeviceReset();
#endif
}


void test1_SPME_EwaldSum_3d_Exclude_InducedDipole() {
#ifdef __GPU__
	if (cudaSetDevice(_CudaDevice_) != cudaSuccess) {
		cout<<"GPU is not enabled!"<<endl; return;
	}
	cudaDeviceReset();
#endif

	int nThreads = MAX_THREADS;
	char msg[256] = "\0";

	using namespace _atomic_cmm_;

	bool bExplicit = false; // calculate interaction explicitly

	::bPolarize = true; // polarizable water molecule
	::bExcludeInducedDipole = true; // excluding dipole-dipole interaction


	::eps_dc = 1;
	MMOL_MD_CELL mdcell;
	//int nmol = 2;
	//int nmol = 8;
	//int nmol = 10;
	int nmol = 20;
	//int nmol = 52;
	//int nmol = 512;
	//int nmol = 1024;
	mdcell.SetMol(0, nmol, 0);
	CMM_CELL3D<BASIC_CLUSTER> cmm;
	::mlog.init("test.log", false);

	//float h[3] = {16, 16, 16};
	//float h[3] = {8, 8, 8};
	//float h[3] = {5, 5, 5};
	//float h[3] = {4, 4, 4};
	//float h[3] = {32, 32, 16};
	float h[3] = {64, 64, 64};

	::rcut_Ewd = h[0] / 2; // full width is 2 * h[0] ...
	
	//int ncell[3] = {5, 5, 5};
	int ncell[3] = {10, 10, 10};
	//int ncell[3] = {30, 30, 30};
	float corner[3] = {-h[0], -h[1], -h[2]};
	float w_cell[3] = {h[0] * 2 / ncell[0], h[1] * 2 / ncell[1], h[2] * 2 / ncell[2]};
	
	int NC = 36;

	cmm.set_cell(20, 20, 20);
	cmm.set_cell_range(true, -h[0], h[0], -h[1], h[1], -h[2], h[2]);
	mdcell.h[0] = h[0]; mdcell.h[1] = h[1]; mdcell.h[2] = h[2];  mdcell.period = true;
	cmm.set_cell_pos();
	int ncs = 100;
	CMM_array_set_storage<BASIC_CLUSTER>(&cmm, ncs);

	int imol = 0, ia, i;
	SMOLECULE cw;
	VECTOR3 dr, axis;
	axis.v[2] = 1;
	//axis.v[0] = 1;

	extern bool ConstructSimpleMolecule(SMOLECULE* sm, char *mol);
	extern bool cp(SMOLECULE *d, SMOLECULE *s, bool setup = true);
	extern void cal_dipole(MMOL_MD_CELL &mmcell, int nThreads, SURF_DIPOLE &mu_surf);
	extern void init_LJ_Pars();
	extern void construct_LJ12_6_prof(int mode);
	extern double cal_charge(MMOL_MD_CELL &mmcell);

	extern void get_bclusters(MMOL_MD_CELL &mdcell, ARRAY<BASIC_CLUSTER*>& bcluster);

	ConstructSimpleMolecule(&cw, "SPC");
	cw.c->cal_geometry_radius();
	dr.v[0] = -cw.r.v[0]; dr.v[1] = -cw.r.v[1]; dr.v[2] = -cw.r.v[2];
	cw.shiftMol(dr);
	::iFormat_LJ = 1;  // Amber
	init_LJ_Pars();
	construct_LJ12_6_prof(::iFormat_LJ);
	for (ia = 0; ia < cw.c->nAtoms; ia++) cw.c->atom[ia].par = NULL;

	double Ep = 0;
	
	//cmm_init_subcell<BASIC_CLUSTER>(cmm); // setup the sub-cells for CMM in different levels when the cell size is bigger than min_cmm_size

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		mdcell.sm.m[imol].c = new BASIC_CLUSTER;
		cp(mdcell.sm.m + imol, &cw, true);
		dr.v[0] = (h[0]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[1] = (h[1]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		dr.v[2] = (h[2]) * (ranf() > 0.5 ? 1 : -1) * ranf();
		mdcell.sm.m[imol].shiftMol(dr);
		rand_vect3(axis.v);
		rotate_SM(mdcell.sm.m[imol], axis, (ranf() > 0.5 ? 1 : -1) * ranf() * PI, true);

		mdcell.sm.m[imol].c->setup_kin_dyn();
		mdcell.sm.m[imol].c->setup_cluster_fatom();
	}
	get_atoms(mdcell, mdcell.atom);
	get_bclusters(mdcell, mdcell.bcluster);
	mdcell.check_eneutral();

	mdcell.set_atom_multipole((bPolarize ? 1 : 0));
	mdcell.setMolIndx();
	molecule_check_periodic_cell(mdcell);
	check_atom_rg(mdcell, nThreads);

	EAtCELL acell;
	acell.rmax_c = 1; // maximum radius of the cluster
	init_EAtCell(mdcell, acell, ncell[0], ncell[1], ncell[2], corner, w_cell);

#if _CHECK_TIME_STAMP_
	mt.start();
#endif

	distribute_atoms(&acell);
	update_atomic_charge_dipole(&acell);

#if _CHECK_TIME_STAMP_
	long dt = mt.elapse();
	sprintf(msg, "used time: %d ms", dt); show_log(msg, true);
#endif

	check_distribution((_AtCELL*)&acell);


	cmm_init_distribute_cluster(mdcell.sm.m, mdcell.sm.n, cmm, true); // add simple solid molecule into cmm
	cmm_check_cluster<BASIC_CLUSTER>(cmm);
	cmm_check<BASIC_CLUSTER>(cmm, MAX_THREADS);
	_cmm_3d_::SMCheckRelshipInCell(mdcell.sm.m, mdcell.sm.n, cmm, MAX_THREADS);

	/*
	{
		for (int imol = 0; imol < mdcell.sm.n; imol++) {
			sprintf(errmsg, "m %d [%f, %f, %f]", imol, mdcell.sm.m[imol].r0.v[0], mdcell.sm.m[imol].r0.v[1], mdcell.sm.m[imol].r0.v[2]);
			show_log(errmsg, true);
			sprintf(errmsg, "m %d has %d of neighbors", imol, number<CMM_IMAGINE<BASIC_CLUSTER> >(mdcell.sm.m[imol].c->spmeRship.ch));
			show_log(errmsg, true);
		}
	}
	*/

	// calculat the total dipole of the whole cell
	cal_dipole(mdcell, 1, mdcell.mu_surf);
	acell.mu_surf = mdcell.mu_surf; cmm.mu_surf = mdcell.mu_surf;

	acell.bCalE = 2; // 1 -- Efield, 2 -- force & U
	acell.bDipole = ((::bPolarize || ::bDipole) ? 1 : 0);

	// SPME parameters
	extern BSplineFunc<2> bsp;

	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 32;
	//::dim_spme[0] = 20; ::dim_spme[1] = 20; ::dim_spme[2] = 20;
	//::dim_spme[0] = 32; ::dim_spme[1] = 32; ::dim_spme[2] = 64;
	::dim_spme[0] = 64; ::dim_spme[1] = 64; ::dim_spme[2] = 64;
	//::dim_spme[0] = 128; ::dim_spme[1] = 128; ::dim_spme[2] = 128;

	//float rspline = 4;
	//int nSpline = int(0.5 * rspline / h[0] * dim_spme[0] + 0.5);
	//init_BSplineFunc(bsp, nSpline);
	init_BSplineFunc(bsp, 6);


	_spme_::VSPME<_VQatom> vspme_q; vspme_q.mMP = 0; // charge only
	vspme_q.bST = true; vspme_q.bST_diagonal = true;
	_spme_::VSPME<_VMUatom> vspme_mu; vspme_mu.mMP = 1; // with dipole
	vspme_mu.bST = true; vspme_mu.bST_diagonal = true;

	float q0 = 0;
	_spme_::VSPME<_VMUatom> vspme_mu_induced; vspme_mu_induced.mMP = 1; // with dipole
	vspme_mu_induced.bST = true; vspme_mu_induced.bST_diagonal = true;

	if (bPolarize) {
		mdcell.init_polarization_buff();
		mdcell.init_dipole_hist();
		mdcell.set_atom_multipole(1); // atoms are all polarizable, if they are not defined neutral

		vspme_mu.bsp = &(bsp);
		init_spme(&mdcell, vspme_mu); // init the viral atoms
		vspme_mu.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_mu.xl[0] = cmm.xl[0]; vspme_mu.xl[1] = cmm.xl[1]; vspme_mu.xl[2] = cmm.xl[2];
		vspme_mu.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();
		vspme_mu.init_k();
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();

		vspme_mu.init_k();
		vspme_mu.init_b(); vspme_mu.init_C(); vspme_mu.init_vars();


		vspme_mu_induced.bsp = &(bsp);
		init_spme_induced(&mdcell, vspme_mu_induced, &(q0)); // init the viral atoms
		vspme_mu_induced.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_mu_induced.xl[0] = cmm.xl[0]; vspme_mu.xl[1] = cmm.xl[1]; vspme_mu.xl[2] = cmm.xl[2];
		vspme_mu_induced.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();

		//init_spme_induced(&mdcell, vspme_mu_induced, &(q0)); // init the viral atoms
		vspme_mu_induced.init_k(); vspme_mu_induced.init_b(); vspme_mu_induced.init_C(); vspme_mu_induced.init_vars();
	}
	else { // charge only
		mdcell.set_atom_multipole(0); // atoms are all polarizable, if they are not defined neutral

		vspme_q.bsp = &(bsp);
		init_spme(&mdcell, vspme_q); // init the viral atoms
		vspme_q.init(bsp.indx, ::dim_spme[0], ::dim_spme[1], ::dim_spme[2]);
		vspme_q.xl[0] = cmm.xl[0]; vspme_q.xl[1] = cmm.xl[1]; vspme_q.xl[2] = cmm.xl[2];
		vspme_q.cell_dim(cmm.xd[0], cmm.xd[1], cmm.xd[2]); // important, after init();
		vspme_q.init_k();
		vspme_q.init_b(); vspme_q.init_C(); vspme_q.init_vars();
	}

	SPME_VAR spme_var;
	spme_var.bVirial = true; spme_var.bST_diagonal = true; 

	spme_var.esv1.bEwald = true; 
	//spme_var.spme_3d_var.esv1.init_EwaldSum(mdcell.h[0] * mdcell.h[1] * mdcell.h[2] * 8, ::rcut_Ewd);
	if (bPolarize) {
		spme_var.esv1.init_EwaldSum(vspme_mu.V, ::rcut_Ewd);
		spme_var.esv1.iSurface = 0; // 3d
	}
	else {
		spme_var.esv1.init_EwaldSum(vspme_q.V, ::rcut_Ewd);
		spme_var.esv1.iSurface = 0; // 3d
	}
	spme_var.bTholePolarize = ::_Thole_polarize_::bTholePolarize;
	spme_var.tv.iTholeModel = ::_Thole_polarize_::iTholeModel;
	spme_var.bEfieldOnly = (acell.bCalE == 1 ? true : false);

	spme_var.esv1.iExIndDipole = (::bExcludeInducedDipole ? 1 : 0); 
	spme_var.esv1.iSurfaceBoundary = ::iSurfaceBoundary; 

	double kappa = spme_var.esv1.kappa;

	acell.esv.cp_EwaldSum_vars(*(_EwaldSum_real_::EwaldSumRealVars0*)(&spme_var.esv1));
	acell.rcut = ::rcut_Ewd;

#ifdef __GPU__
	if (acell.gpu_cell.srnb_dev == NULL && acell.bCalUseGPU) {
		if (!acell.check_neighbor_memory(true)) return;
		gpu_reset_esr_neighbors(&acell);
		check_neighbors(&(acell.gpu_cell), acell.gpu_cell_dev);
	}
#endif

#if _CHECK_TIME_STAMP_
	TIME time;
	{
		//double Ushape = _spme_2d_::ShapeEnergy_2d(cmm);
		//sprintf(errmsg, "shape energy: %f [kT]", Ushape); show_log(errmsg, true);
	}
#endif
	//V3zero(cmm.mu_surf.mu) cmm.mu_surf.q = 0; V3zero(cmm.mu_surf.qr) V3zero(cmm.mu_surf.mu0)
	InteractRes ires;

#if _CHECK_TIME_STAMP_
	time.start();
#endif
	if (bPolarize) {
		if (::bExcludeInducedDipole) {
			_atomic_cmm_::Polarize_SPME_Interact2(mdcell, acell, vspme_mu, vspme_mu_induced, nThreads, true, spme_var, ires);
		}
		else if (!_atomic_cmm_::Polarize_SPME_Interact(mdcell, acell, vspme_mu, nThreads, true, spme_var, ires)) {
			sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
		}
	}
	else _atomic_cmm_::SPME_EF(acell, vspme_q, nThreads, spme_var, ires);

	Ep = ires.U;
	sprintf(errmsg, "CPU-SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);

#if _CHECK_TIME_STAMP_
	dt = time.elapse();
	sprintf(msg, "CPU SPME used time : %d ms", dt); show_infor(msg, true);
#endif

	sprintf(errmsg, "total dipole: [%f, %f, %f]; from q: [%f, %f, %f], from induced dipole: [%f, %f, %f]", cmm.mu_surf.mu.v[0], cmm.mu_surf.mu.v[1], cmm.mu_surf.mu.v[2], cmm.mu_surf.qr.v[0], cmm.mu_surf.qr.v[1], cmm.mu_surf.qr.v[2], cmm.mu_surf.mu_ind.v[0], cmm.mu_surf.mu_ind.v[1], cmm.mu_surf.mu_ind.v[2]); show_log(errmsg, true);

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);
	
	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-CPU-SPME-EwaldSum.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);

	show_log("reset force & dipole!", true);

	for (imol = 0; imol < mdcell.sm.n; imol++) {
		V6zero(mdcell.sm.m[imol].c->dyn->fc)
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			V3zero(mdcell.sm.m[imol].c->atom[ia].mu)
			V3zero(mdcell.sm.m[imol].c->atom[ia].ind_mu)
			if (mdcell.sm.m[imol].c->atom[ia].mu_hist != NULL) {
				for (i = 0; i < 3; i++) {
					V3zero(mdcell.sm.m[imol].c->atom[ia].mu_hist->mu[i])
				}
			}
			V3zero(mdcell.sm.m[imol].c->atom[ia].E)
			V3zero(mdcell.sm.m[imol].c->atom[ia].F)
			for (i = 0; i < mdcell.sm.m[imol].c->atom[ia].f.n; i++) {V3zero(mdcell.sm.m[imol].c->atom[ia].f.m[i]);}
		}
	}


	show_log("CPU-SPME is done !", true);

	//******************************
	//***** using GPU for SPME *****
	//******************************
	ires.reset();
#ifdef __GPU__

	_atomic_cmm_::update_atoms(&acell);
	// calculat the total dipole of the whole cell
	cal_dipole(mdcell, 1, mdcell.mu_surf);
	acell.mu_surf = mdcell.mu_surf; cmm.mu_surf = mdcell.mu_surf;
	_atomic_cmm_::update_atomic_charge_dipole(&acell); // update charge and dipole on cell

	show_log("GPU ?", true);

	size_t size_free = 0, size_total = 0;
	cudaMemGetInfo(&size_free, &size_total);
	if (size_free == 0) {
		show_msg("GPU failure to initialize the memory for SPME", true); return;
	}

	show_log("initialize GPU configuration", true);

	_gpu_md_::GPU_SPME_3d spme_gpu;
	spme_gpu.atom.set_array(mdcell.atom.n); _gpu_md_::init_atom(mdcell.atom, spme_gpu.atom);
	spme_gpu.init_memory(mdcell.atom.n); 
	if (!_gpu_md_::setup_gpu_spme(&spme_gpu, &::bsp, cmm.xl, cmm.xd, ::dim_spme, spme_var.esv1.kappa, mdcell.atom.n, (bPolarize ? 1: 0))) {
			show_msg("GPU failure to initialize the memory for SPME", true); return;
	}

	memcpy(spme_gpu.r.m, acell.r.m, acell.r.n * sizeof(double));
	//_gpu_md_::update_atomic_coordinate(&spme_gpu, (MAX_THREADS < 4 ? MAX_THREADS : 4));

	memcpy(spme_gpu.q.m, acell.c.m, acell.c.n * sizeof(double));
	_gpu_md_::_spme_::update_q_device_spme(&(spme_gpu.vspme_host));

	if (bPolarize) memcpy(spme_gpu.mu.m, acell.mu.m, acell.mu.n * sizeof(double));
	else memset(spme_gpu.mu.m, 0, spme_gpu.mu.n * sizeof(double));
	_gpu_md_::_spme_::update_dipole_device_spme(&(spme_gpu.vspme_host));

	show_log("start GPU-interaction calculation", true);

#if _CHECK_TIME_STAMP_
	time.start();
#endif
	bool bstatus = _gpu_dstruct_::_gpu_util_::check_gpu();
	if (bPolarize) {
		if (::bExcludeInducedDipole) {
			if (!_atomic_cmm_::Polarize_SPME_Interact2(mdcell, acell, spme_gpu, nThreads, true, spme_var, ires)) {
				sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
			}
		}
		else {
			if (!_atomic_cmm_::Polarize_SPME_Interact(mdcell, acell, spme_gpu, nThreads, true, spme_var, ires)) {
				sprintf(errmsg, "polarized dipole is not convergent to 0.01D"); show_log(errmsg, true);
			}
		}
	}
	else _atomic_cmm_::SPME_EF(acell, spme_gpu, nThreads, spme_var, ires, true);


	Ep = ires.U;
	sprintf(errmsg, "GPU-SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);

#if _CHECK_TIME_STAMP_
	dt = time.elapse();
	sprintf(msg, "GPU SPME used time : %d ms", dt); show_infor(msg, true);
#endif


	_gpu_md_::_spme_::release_device_spme(&(spme_gpu.vspme_host), spme_gpu.vspme_dev);
	spme_gpu.vspme_dev = NULL;
	spme_gpu.status.release(); spme_gpu.E.release(); spme_gpu.F.release(); spme_gpu.mu.release(); spme_gpu.q.release();
	spme_gpu.release_plans();

#endif //#ifdef __GPU__

	sprintf(errmsg, "total dipole: [%f, %f, %f]; from q: [%f, %f, %f], from induced dipole: [%f, %f, %f]", cmm.mu_surf.mu.v[0], cmm.mu_surf.mu.v[1], cmm.mu_surf.mu.v[2], cmm.mu_surf.qr.v[0], cmm.mu_surf.qr.v[1], cmm.mu_surf.qr.v[2], cmm.mu_surf.mu_ind.v[0], cmm.mu_surf.mu_ind.v[1], cmm.mu_surf.mu_ind.v[2]); show_log(errmsg, true);

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);
	
	sprintf(errmsg, "SPME Ewald-Sum : total energy %f [kT]", Ep); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-GPU-SPME-EwaldSum.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	sprintf(errmsg, "SPME efield: [%f, %f, %f]", mdcell.sm.m[0].c->atom[0].E.v[0], mdcell.sm.m[0].c->atom[0].E.v[1], mdcell.sm.m[0].c->atom[0].E.v[2]);
	show_log(errmsg, true);


	{
	ofstream out;
	out.open("test-atomPos.dat");
	out<<"atom pos. : "<<endl;
	int ia;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"[ SM "<<imol<<" ] : "<<mdcell.sm.m[imol].r.v[0]<<"  "<<mdcell.sm.m[imol].r.v[1]<<"  "<<mdcell.sm.m[imol].r.v[2]<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<mdcell.sm.m[imol].c->atom[ia].r.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].r.v[2]<<endl;
		}
		out<<endl;
	}
	out.close();
	}


	if (bPolarize) vspme_mu.release_plans();
	else vspme_q.release_plans();

	if (!bExplicit) return;

	// explicit calculation
	BASIC_CLUSTER *pc = NULL;
	MATOM *patom = NULL;
	VECTOR3 F, torque;

	for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
		V3zero(mdcell.sm.m[imol].c->atom[ia].E)
		V3zero(mdcell.sm.m[imol].c->atom[ia].F)
		V3zero(mdcell.sm.m[imol].c->atom[ia].mu)
		V3zero(mdcell.sm.m[imol].c->atom[ia].ind_mu)
		for (i = 0; i < mdcell.sm.m[imol].c->atom[ia].f.n; i++) {V3zero(mdcell.sm.m[imol].c->atom[ia].f.m[i]);}
		for (i = 0; i < mdcell.mu_hist.n; i++) {V3zero(mdcell.mu_hist.m[i].mu[0]); V3zero(mdcell.mu_hist.m[i].mu[1]); V3zero(mdcell.mu_hist.m[i].mu[2]);}
	}

	EInteractVar evar;
	evar.bVirial = true; evar.bST_diagonal = true; 
	evar.bEfieldOnly = false; 
	evar.esv1.bEwald = false; // direct interaction, not EwaldSum
	evar.esv1.iSurface = 0; //3d
	evar.esv1.iSurfaceBoundary = 0; // normal
	evar.esv1.iExIndDipole = spme_var.esv1.iExIndDipole;
	InteractRes tres;

	cmm.mu_surf = mdcell.mu_surf;

	ires.reset();
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		tres.reset();
		ExplicitEInteract(cmm, mdcell.sm.m[imol].c, evar, tres);
		
		ires += tres;
	}
	Ep = ires.U;

	if (mdcell.mm.m != NULL) MOperate<MMOLECULE>((void*)(&MM_TOTAL_FORCE), mdcell.mm.m, mdcell.mm.n, nThreads);
	if (mdcell.sm.m != NULL) MOperate<SMOLECULE>((void*)(&SM_TOTAL_FORCE), mdcell.sm.m, mdcell.sm.n, nThreads);
	if (mdcell.pm.m != NULL) MOperate<PMOLECULE>((void*)(&PM_TOTAL_FORCE), mdcell.pm.m, mdcell.pm.n, nThreads);

	sprintf(errmsg, "explicit calculation : total energy %f [kT], with dipole-dipole interaction: %f [kT]", Ep, tres.U); show_log(errmsg, true);
	{
	ofstream out;
	out.open("test-explicit.dat");
	out<<"cluster force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<imol<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[0]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[1]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[2]<<"      ";
		out<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[3]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[4]<<"  "<<mdcell.sm.m[imol].c->dyn->fc.v[5]<<endl;
	}

	out<<endl<<"dipole : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].mu.v[2]<<endl;
		}
	}

	out<<endl<<"electric field : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].E.v[2]<<endl;
		}
	}

	out<<endl<<"atomic force : "<<endl;
	for (imol = 0; imol < mdcell.sm.n; imol++) {
		out<<"SM : "<<imol<<endl;
		for (ia = 0; ia < mdcell.sm.m[imol].c->nAtoms; ia++) {
			out<<ia<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[0]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[1]<<"  "<<mdcell.sm.m[imol].c->atom[ia].F.v[2]<<endl;
		}
	}

	out<<endl<<"virial : "<<endl;
	out<<"Total virial : "<<ires.STtr<<endl;
	out<<"Virial diagonal terms : "<<ires.STxx<<", "<<ires.STyy<<", "<<ires.STzz<<endl;

	out.close();
	}
	
	::mlog.close();
}

