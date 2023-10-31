#ifdef __GPU__
#include "gpu_sr.h"
//#include "gpu_esr.h"
#endif

namespace _atomic_cmm_ {
// NC is the maximum number of atoms inside each CMM cell

// to initialize the cell, require estimated parameters:
//     Nca -- maximum number of atoms in each cluster
//     NC -- maximum number of atoms in each cell

// other parameters: Na -- total number of atoms
//                   Nc -- total number of cluster
//                   Nx, Ny, Nz -- dimension of CMM cell
//                   cx[3], cw[3] -- starting position of CMM, full-width of each small cell

// distributing the atoms in cluster, into a grid-CMM
	template <class _Atom> class IndxAtom {
	public:
		_Atom *m;
		int indx;
	};

	template <class _Atom> class IndxArray : public ARRAY< IndxAtom<_Atom> > {
	public:
		IndxAtom<_Atom> *get_indxAtom(_Atom *a) {
			int i;
			for (i = 0; i < ARRAY< IndxAtom<_Atom> >::n; i++) {
				if (ARRAY< IndxAtom<_Atom> >::m[i].m == a) return ARRAY< IndxAtom<_Atom> >::m + i;
			}
			return NULL;
		};
	};


	class _AtCELL {
	public:
		IndxArray< BASIC_CLUSTER > scluster; // the clusters to be distributed to the CMM
		IndxArray< MATOM > satom; // the atoms defined in the cluster

		// clusters -- not required for GPU 
		int Nc; // number of cluster
		int Nca; // maximum number of atom in each cluster
		ARRAY<int> cluster; // the definition of cluster, dimension -- (Nca + 1) * Nc
		ARRAY<int> ixc, iyc, izc; // index of the cell for each cluster, in dimension Nc; the index is based on the coordination vp of the cluster
		ARRAY< ARRAY_FLEX<int> > cmm_c; // dimension Nx * Ny * Nz, each array ARRAY_FLEX<int> containes the clusters in the unit cell

		void distribute_cluster(); // distribute the cluster into each unit cell of cmm_c


		// following parameters / arrays are to be used in GPU

		ARRAY<int> ivar;
		ARRAY<double> dvar;

		int Na; // number of atoms
		ARRAY<int> aID; // type of atoms -- relating to L-J or other short-range interaction, and the charge, polarization ...
	                  // dimension -- Na
		ARRAY<int> cIndx; // which cluster it belongs, the cluster is considered as a solid structure
	                    // dimension -- Na
	
		// bound and bound14 are from the hinges in macromolecules
		// bound includes 1-2 bound and 1-3 bound
		// and bounds a-b and b-a are one bound
		int Nb; // number of bounds
		ARRAY<int> bound; // the bounds [atom1_index, atom2_index]
	                  // so, its dimension -- 2 * Nb
		int Nb14; // 1-4 neighbors
		ARRAY<int> bound_14; // dimension -- 2 * Nb14; each bound has parameters -- a1_indx, a2_indx, 
		double c14;
		//double c14_e; // coefficient for 1-4 electrostatic interaction
		//double c14_sr; // coefficient for 1-4 short-range interaction

		// for atoms with Non-Zero Bound
		// relationship of bound and bound_14 with the atom needs to be setup
		// this relationship is not required for GPU calculation
		int Na_b, Na_b14; // the atoms related to bound and bound_14
		int Nm_b, Nm_b14; // maximum bnumber of bound, and maximum number of 1-4 bound, for each atom
		ARRAY<int> b_relship, b14_relship; // they have dimension Na_b * (Nm_b+2) and Na_b14 * (Nm_b14 + 2) respectively
	                                   // bounds of each atom is saved in format: [bound_atom_indx] [number bound] [bound_index in bound / bound_14]
		// end here

		ARRAY< ARRAY_FLEX<int> > b_relship_t, b14_relship_t; // temparory relationship for all atoms, dimension Na

		// for short-range interaction, the repulsion force increase dramatically when two atoms are too close
		// if we calculate all interactions first and then subtract bound interaction, the accuracy on the long-range interactions will be changed since they are weak
		// so, it will be better to avoid the calculation of bound interaction, and we need an array to keep the bound atoms for each atom, between the joint clusters
		// This is required by GPU
		// for all atoms, even with Zero-Bound, format : [number of bound atoms] [bound atom indx ...]
		ARRAY<int> abound; // 1-2 and 1-3 bound for each atom, each atom has possible number of bound Nm_b, and the first element is to save the number of neighbor for the atom
		                                  // so the dimension of array is (Nm_b + 1) * Na
	                   
	
		ARRAY<double> r; // positions of atoms, with dimension Na * 3
		ARRAY<int> ir; // the cubic index in CMM, with dimension Na * 3

		int NC; // maximum number of atoms in each cell, NC can not be changed when the CELL is constructed !
		int Nx, Ny, Nz; // dimension of CMM cell
		int Nyz; // Ny * Nz
		double cx[3]; // starting position of CMM, most corner position, [0][0][0]
		double cw[3]; // the width of the cell
		double rmax_c; // maxium radius of the cluster
		ARRAY<int> cmm; // the atoms inside each unit cell in CMM, the atoms are indexed 0, 1 ... for each unit cell, the format is: 0 -- [real number of atoms], 1 ~ NC+1 -- atoms
	                // so the whole dimension is (NC+1) * Nx * Ny * Nz

		double rcut; // the truncation distance for short-range interaction, real part Ewald-Sum

		inline int* cubic(int ix, int iy, int iz) {
			int ic = ix * Nyz + iy * Nz + iz;
			return cmm.m + ic * (NC + 1);
		};
	};

#ifdef __GPU__
	bool init_gpu_cell(_AtCELL *acell, _gpu_md_::GPU_AtCELL* gpu_cell);
#endif

#ifndef x_idx 
#define x_idx(i) i * 3
#define y_idx(i) i * 3 + 1
#define z_idx(i) i * 3 + 2
#endif

#ifndef b_idx1
#define b_idx1(i) i * 2
#define b_idx2(i) i * 2 + 1
#endif


	class EAtCELL : public _AtCELL {
	public:
#ifdef __GPU__
#define cuda_release(v) if (v != NULL) {cudaFree(v); v = NULL;}
		int igpuDev;
		bool bUseGPU; // using GPU, e.g. the memory of GPU and even interaction calculation ?
		bool bCalUseGPU; // using GPU for interaction calculation ?
		_gpu_md_::GPU_EAtCELL gpu_cell;
		_gpu_md_::GPU_EAtCELL *gpu_cell_dev;
		bool init_gpu_cell() {
			if (!bUseGPU) return true;
			if (gpu_cell_dev == NULL) {
				if (cudaSetDevice(_CudaDevice_) != cudaSuccess) {
					show_infor("GPU is not enabled!", true); return false;
				}
				if (cudaMalloc(&gpu_cell_dev, sizeof(_gpu_md_::GPU_EAtCELL)) != cudaSuccess) {gpu_cell_dev = NULL; return false;}
				else {
					cudaMemcpy(gpu_cell_dev, &gpu_cell, sizeof(_gpu_md_::GPU_EAtCELL), cudaMemcpyHostToDevice);
					return true;
				}
			}
			else return true;
		};
		bool check_neighbor_memory(bool b_rcut_changed);
		void release_device() {
			gpu_cell.release();
			if (gpu_cell_dev != NULL) {cudaFree(gpu_cell_dev); gpu_cell_dev = NULL;};
		};

		EAtCELL() {
			bUseGPU = true; bCalUseGPU = true;
			gpu_cell_dev = NULL; _AtCELL::rmax_c = 1;
			igpuDev = _CudaDevice_;
		};
		~EAtCELL() {release_device();};
		
		void set_EwaldSumControlVars_3d() {
			gpu_cell.esv.bEwald = esv.bEwald;
			gpu_cell.esv.iExIndDipole = esv.iExIndDipole;
			gpu_cell.esv.iSurface = esv.iSurface;
			gpu_cell.esv.iSurfaceBoundary = esv.iSurfaceBoundary;
			gpu_cell.esv.kappa = esv.kappa;
			gpu_cell.esv.V = 4 * PI / esv.inv_V;
		};
		void set_EwaldSumControlVars_2d() {
			set_EwaldSumControlVars_3d();
			gpu_cell.esv.bEwald = esv.bEwald;
			gpu_cell.esv.iExIndDipole = esv.iExIndDipole;
			gpu_cell.esv.iSurface = esv.iSurface;
			gpu_cell.esv.iSurfaceBoundary = esv.iSurfaceBoundary;
			gpu_cell.esv.kappa = esv.kappa;
			
			gpu_cell.esv.A = k0var.A;
		};

		bool gpu_prior_cal(int dim_EwaldSum) {
			if (!bUseGPU || !bCalUseGPU) return true;
			init_gpu_cell();
			if (dim_EwaldSum == 2) set_EwaldSumControlVars_2d(); // 2d -- a control parameter could be changed
			else if (dim_EwaldSum == 3) set_EwaldSumControlVars_3d(); // 3d -- a control parameter could be changed
			else {show_msg("wrong dimension for Ewald-Sum", true); return false;}
			memcpy(gpu_cell.mu_surf, mu_surf.mu.v, SIZE_V3); // surface dipole
			memcpy(gpu_cell.ind_mu_surf, mu_surf.mu_ind.v, SIZE_V3); // induced dipole
			gpu_cell.bCalE = this->bCalE;
			gpu_cell.bDipole = this->bDipole;
			bool status = (fabs(gpu_cell.rcut - this->rcut) < 1 ? true : false);
			if (!status) {
				gpu_cell.rcut = this->rcut;
				if (!check_neighbor_memory(true)) return false;
			}
			if (cudaMemcpy(gpu_cell_dev, &gpu_cell, sizeof(_gpu_md_::GPU_EAtCELL), cudaMemcpyHostToDevice) != cudaSuccess) return false;
			if (!status) {
			{
				char buffer[256] = "\0";
				sprintf(buffer, "final number of neighbor buffer: %d", gpu_cell.n_srnb); show_log(buffer, true);
			}}
			return true;
		};
#undef cuda_release
#endif

		SURF_DIPOLE mu_surf; // surface dipole

		ARRAY<double> c; // charge with dimension Na
		ARRAY<double> mu; // dipole with dimension Na * 3

		//Results:
		ARRAY<double> E; // the electrostatic field of each atoms, dimension -- Na * 3
		ARRAY<double> Fe; // force from electrostatic interaction, dimension -- Na

		// temperary result:
		ARRAY<double> E_14; // electrostatic field from the exclusion part, e.g. 1-4, actually the field is on the first atom defined in the EwdRel_exclude, a1_indx of each pair
			                // the field on another atom, a2_indx, should be negative values, or inverse direction
				            // dimension -- Nb14 * 3
		ARRAY<double> Fe_14; // electrostatic force from the exclusion part, e.g. 1-4, actually the field is on the first atom defined in the EwdRel_exclude, a1_indx of each pair
		                     // the field on another atom, a2_indx, should be negative values, or inverse direction
			                 // dimension -- Nb14 * 3
		ARRAY<double> E_b; // electrostatic field between bounded atoms, actually the field is on the first atom defined in the bound
		                   // the field on another atom, a2_indx, should be negative values, or inverse direction
			               // dimension -- Nb * 3
		ARRAY<double> Fe_b; // electrostatic force from the exclusion part, e.g. 1-4, actually the field is on the first atom the bound
		                    // the field on another atom, a2_indx, should be negative values, or inverse direction
			                // dimension -- Nb * 3


		ARRAY<double> ve, U_estat;  // virial coefficient, dimension -- Na * 3
		ARRAY<double> ve_b, Ub_estat; // virial coefficient from bound electrostatic interaction which will be excluded, dimension -- Nb * 3
	
		ARRAY<double> ve_14, U14_estat; // virial coefficient from electrostatic interaction which will be excluded, dimension -- Nb14 * 3

		// some parameters for calculation controlling
		// ## calculating electrostatic interaction 
		char bDipole; // using dipole of atom ? [0/1]
		char bCalE; // 0 -- do not do calculation
		            // 1 -- calculate electric field
			        // 2 -- calculating force, and energy
				    // 3 -- calculating electric field, force and energy. But is this necessary?
		// Ewald-Sum parameters
		// ...

		void reset_efield();
		void reset_eforce();

		// real part of EwaldSum
		_EwaldSum_real_::EwaldSumRealVars1 esv;
		// term with k=0 in reciprocal space for 2d Ewald-Sum
		_EwaldSum_2d_::EwaldSumRecpK0Var k0var;
	};

	class SR_AtCELL : public _AtCELL {
	public:

#ifdef __GPU__
		int igpuDev;
		bool bUseGPU; // using GPU, e.g. the memory of GPU and even interaction calculation ?
		bool bCalUseGPU; // using GPU for interaction calculation ?
		_gpu_md_::GPU_SR_AtCELL gpu_cell;
		_gpu_md_::GPU_SR_AtCELL *gpu_cell_dev;
		bool init_gpu_cell() {
			if (!bUseGPU || !bCalUseGPU) return true;
			if (gpu_cell_dev == NULL) {
				if (cudaMalloc(&gpu_cell_dev, sizeof(_gpu_md_::GPU_SR_AtCELL)) != cudaSuccess) {gpu_cell_dev = NULL; return false;}
				else {
					cudaMemcpy(gpu_cell_dev, &gpu_cell, sizeof(_gpu_md_::GPU_SR_AtCELL), cudaMemcpyHostToDevice); 
					//cudaMemcpy(&gpu_cell, gpu_cell_dev, sizeof(_gpu_md_::GPU_SR_AtCELL), cudaMemcpyDeviceToHost);
					return true;
				}
			}
			else return true;
		};
		void release_device() {
			gpu_cell.release();
			if (gpu_cell_dev != NULL) {cudaFree(gpu_cell_dev); gpu_cell_dev = NULL;};
		};

		SR_AtCELL() {
			bUseGPU = true; bCalUseGPU = true;
			gpu_cell_dev = NULL; _AtCELL::rmax_c = 1; igpuDev = _CudaDevice_;
		};
		~SR_AtCELL() {release_device();};

		bool gpu_prior_cal() {
			if (!bUseGPU) return true;
			gpu_cell.bCalSR = this->bCalSR;
			init_gpu_cell();
			if (cudaMemcpy(gpu_cell_dev, &gpu_cell, sizeof(_gpu_md_::GPU_SR_AtCELL), cudaMemcpyHostToDevice) == cudaSuccess) {
				//cudaMemcpy(&gpu_cell, gpu_cell_dev, sizeof(_gpu_md_::GPU_SR_AtCELL), cudaMemcpyDeviceToHost);
				return true;
			}
			else return false;
		};
#endif

		ARRAY<char> sr_type; // type of short-range interaction, dimension -- Na
			                 // 0x01 -- L-J interaction, 0x02 -- given profile
		
		int n_srpars; // number of short-range interaction parameters for each atom, in sr_par
		ARRAY<double> sr_par; // short-range interaction parameters, e.g. for Lennard-Jones, eps, r_LJ is required for each atom
							  // dimension -- n_srpars * Na

		ARRAY<char> bShocked; // whether the atom is too close to another atom to be shocked, return result, with dimension Na

		ARRAY<double> Fsr; // the force from short-range interaction for each atoms, dimension -- Na * 3
		ARRAY<double> Fsr_b; // short-range force from the exclusion part, bounded atom, actually the field is on the first atom the bound
			                 // the field on another atom, a2_indx, should be negative values, or inverse direction
				             // dimension -- Nb * 3
		ARRAY<double> Fsr_14; // short-range force from the exclusion part, bounded atom, actually the field is on the first atom the bound
		                      // the field on another atom, a2_indx, should be negative values, or inverse direction
			                  // dimension -- Nb14 * 3
		ARRAY<double> vsr, U_sr;     // virial coefficient, dimension -- Na * 3
		ARRAY<double> vsr_b, Ub_sr; // virial coefficient from bound short-range interaction which will be excluded, dimension -- Nb * 3
		ARRAY<double> vsr_14, U14_sr; // virial coefficient from electrostatic interaction which will be excluded, dimension -- Nb14 * 3

		// some parameters for calculation controlling
		char bCalSR; // calculating short-range interaction ? [0/1]

		void reset_sr();
	};

	void init_EAtCell(MMOL_MD_CELL &mdcell, EAtCELL &cell, int Nx, int Ny, int Nz, float *cx, float *cw);
	void init_SRAtCell(MMOL_MD_CELL &mdcell, SR_AtCELL &acell, int Nx, int Ny, int Nz, float *cx, float *cw);

	//void _distribute_atoms_cpu(_AtCELL *acell);
	//void _update_atoms(_AtCELL *cell); // coordination only, keep the distribution
	void distribute_atoms(EAtCELL *acell);
	void distribute_atoms(SR_AtCELL *acell);
	void update_atoms(EAtCELL *cell); // coordination only, keep the distribution
	void update_atoms(SR_AtCELL *cell); // coordination only, keep the distribution

	void check_distribution(_AtCELL *acell);

	void update_atomic_charge_dipole(EAtCELL *cell);

	void GPU_update_atomic_charge_dipole(EAtCELL *cell);

	void GPU_collecting_efield(JobIndx *job, EAtCELL *cell);
	void GPU_collecting_eforce(JobIndx *job, EAtCELL *cell);

	void GPU_collecting_srforce(JobIndx *job, SR_AtCELL *cell);

	bool GPU_init_srAtCell(SR_AtCELL *acell);
	bool GPU_init_EAtCell(EAtCELL *acell);

	void InterfaceConstraints_AtCell(_atomic_cmm_::_AtCELL *cell, int nThreads, double &Ut);

}; // end of namespace
