#ifdef __GPU__

namespace _gpu_md_ {
#ifndef SIZE_V3
#define SIZE_V3      24
#endif

#define r0LJ_max   5  // maximum rLJ / r, or minumum r = rLJ * r0LJ_min

#define U_InertiaMoment   1.66053886                  // x 10^-17 J; I * 1/dt^2 -- (u*Angs^2) * (1/fs^2)
#define U_ElectFr         0.230726                    // x 10^-17 J; e * e / (4PI*eps0 * r^2) * r -- (electron^2) / (Angs)
#define U_MassAccel       1.66053886                  // x 10^-17 J / Angs ; m * a -- (u) * (Angs / fs^2)
#define U_ElectForce      0.230726                    // x 10^-17 J / Angs; e * e /(4PI*eps0 * r^2) -- (electron^2) / (Angs^2)

#define U_eps      1.6605778811026237130521421454666e-4
#define U_LJForce  1.6605778811026237130521421454666e-4
//float U_eps = 0.001f / 6.022f; // x 10^-17 J; eps in unit kJ / mol
//float U_LJForce = U_eps; // x 10^-17 J / Angs.; LJ force = eps / Angs.

// 1 e^2/(1 Ang) = U_ElectForce / k = 1.6711404 x 10^5 K = 167114.04 K
// 1kJ/mol = U_eps / k = 120.27508 K

#define U_EextFr     1.6e-6   // x 10^-17 J; external E in unit kV/mm, charge has unit e
#define U_EextForce  1.6e-6   // x 10^-17 J / Angs; external E in unit kV/mm, charge has unit e

	// some global variables will be used in GPU-code, so they need to be transfered from CPU-program
	class gvar {
	public:
		double kT, unit_Ek_kT;
		double eUnit_estat_kT, eUnit_ext_kT, eUnit_LJ_kT;
		double fUnit_estat_mass, fUnit_ext_mass, fUnit_ext_estat, fUnit_LJ_mass;
		void init_gvar(double T) {
			kT = 1.38065e-6 * (273.15 + T);
			unit_Ek_kT = U_InertiaMoment / kT;
			eUnit_estat_kT = U_ElectFr / kT;
			eUnit_ext_kT = U_EextFr / kT;
			eUnit_LJ_kT = U_eps / kT;
			fUnit_estat_mass = U_ElectForce / U_MassAccel;
			fUnit_ext_mass = U_EextForce / U_MassAccel;
			fUnit_ext_estat = U_EextFr / U_ElectFr;
			fUnit_LJ_mass = U_eps / U_MassAccel;
		};
	};

/*
float kT = (float)(1.38065e-6) * (273.15f + 25); // x 10^-17 J @ 25C  == 0.41164 x 10^-20 J;  1 kT * 1mol = 2.478 kJ

float unit_Ek_kT = U_InertiaMoment / kT;

float eUnit_estat_kT = U_ElectFr / kT;
float eUnit_ext_kT = U_EextFr / kT;
float eUnit_LJ_kT = U_eps / kT;
float fUnit_estat_mass = U_ElectForce / U_MassAccel;
float fUnit_ext_mass = U_EextForce / U_MassAccel;
float fUnit_ext_estat = U_EextFr / U_ElectFr;
float fUnit_LJ_mass = U_eps / U_MassAccel;
*/

//float fUnit_atm_u = 6.101935e-9; // 1 atmosphere = 101325 Pa (N/m^2) = 6.101935 x 10^-9 atomic_mass_unit (u)/ (fs^2 * A)
// for idea gas, each molecule has volumn 4.08841 x 10^4 Angs^3 at 1 atmosphere and 300K, 34.45^3 Angs^3

//float Eext = 0; // external E, in unit kV/mm
//bool bEext = false; // apply external electric field



#ifndef x_idx 
#define x_idx(i) i * 3
#define y_idx(i) i * 3 + 1
#define z_idx(i) i * 3 + 2
#endif

#ifndef b_idx1
#define b_idx1(i) i * 2
#define b_idx2(i) i * 2 + 1
#endif

	class EwaldSumReal_ControlVar {
	public:
		double kappa, V; // in 3d, V is the volumn
		bool bEwald; // calculate with EwaldSum, or direct interaction ?
		short iSurface; // 0 -- cubic use 4 * PI / 3V,  1 -- slab use 4PI/V
		short iSurfaceBoundary; // cubic or 2d cell 0/1
		short iExIndDipole; // excluding induced dipole -- 0/1

		double A; // in 2d, A is the surface area
	};


// this class, especially the memory pointer, is to be used at GPU, for the device function running on GPU
// because we can not write GPU memory like a memory at host, a structure has to be copied to GPU with cudaMemcpy
// so, we have to create a structure at host, but with memory-pointer member, pointing to the memory at GPU
// then, create a class at GPU, use cudaMemcpy to copy the host-one to GPU-one

// and because cudaFree(...) is a C-function, I believe that it does trigger the destructor of the class
// it just release the GPU-memory of its member-variable, not the pointer-related GPU-memory
// so, we have to use the destructor at host, to release the related GPU-memory pointed by the pointer
	class GPU_AtCELL { // memory pointer with _map is the page-locked mapped memory on host, with _dev is the memory on GPU device
	public:
		// clusters -- not required for GPU 
		int Nc; // number of cluster
		int Nca; // maximum number of atom in each cluster

		// **** cluster is nt on GPU ****
		int *cluster; // the definition of cluster, dimension -- (Nca + 1) * Nc
		//ARRAY<int> ixc, iyc, izc; // index of the cell for each cluster, in dimension Nc; the index is based on the coordination vp of the cluster
		//ARRAY< ARRAY_FLEX<int> > cmm_c; // dimension Nx * Ny * Nz, each array ARRAY_FLEX<int> containes the clusters in the unit cell

		// following parameters / arrays are to be used in GPU

		//ARRAY<int> ivar;
		//ARRAY<double> dvar;
		int *ivar_dev, *ivar_map;
		double *dvar_dev, *dvar_map;

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

		// for short-range interaction, the repulsion force increase dramatically when two atoms are too close
		// if we calculate all interactions first and then subtract bound interaction, the accuracy on the long-range interactions will be changed since they are weak
		// so, it will be better to avoid the calculation of bound interaction, and we need an array to keep the bound atoms for each atom, between the joint clusters
		// This is required by GPU
		// for all atoms, even with Zero-Bound, format : [number of bound atoms] [bound atom indx ...]
		int *abound_map, *abound_dev;
		//ARRAY<int> abound; // 1-2 and 1-3 bound for each atom, each atom has possible number of bound Nm_b, and the first element is to save the number of neighbor for the atom
		                                     // so the dimension of array is (Nm_b + 1) * Na

		//ARRAY< ARRAY_FLEX<int> > b_relship_t, b14_relship_t; // temparory relationship for all atoms, dimension Na
	 
		double *r_map, *r_dev;
		//ARRAY<double> r; // positions of atoms, with dimension Na * 3
		int *ir_map, *ir_dev;
		//ARRAY<int> ir; // the cubic index in CMM, with dimension Na * 3

		int NC; // maximum number of atoms in each cell, NC can not be changed when the CELL is constructed !
		int Nx, Ny, Nz; // dimension of CMM cell
		int Nyz; // Ny * Nz
		double cx[3]; // starting position of CMM, most corner position, [0][0][0]
		double cw[3]; // the width of the cell
		double rmax_c; // maxium radius of the cluster
		int *cmm_map, *cmm_dev;
		//ARRAY<int> cmm; // the atoms inside each unit cell in CMM, the atoms are indexed 0, 1 ... for each unit cell, the format is: 0 -- [real number of atoms], 1 ~ NC+1 -- atoms
			            // so the whole dimension is (NC+1) * Nx * Ny * Nz

		double rcut; // the truncation distance for short-range interaction, real part Ewald-Sum

		GPU_AtCELL() {
			ivar_map = NULL; ivar_dev = NULL;
			dvar_map = NULL; dvar_dev = NULL;
			aID_map = NULL; aID_dev = NULL;
			cIndx_map = NULL; cIndx_dev = NULL;
			bound_map = NULL; bound_dev = NULL;
			bound_14_map = NULL; bound_14_dev = NULL;
			b_relship_map = NULL; b_relship_dev = NULL;
			b14_relship_map = NULL; b14_relship_dev = NULL;
			abound_map = NULL; abound_dev = NULL;
			r_map = NULL; r_dev = NULL;
			ir_map = NULL; ir_dev = NULL;
			cmm_map = NULL; cmm_dev = NULL;
		};
		void release() {
#define cuda_release(v) if (v != NULL) {cudaFree(v); v = NULL;}
			cuda_release(ivar_dev);// cuda_release(ivar_map);
			cuda_release(dvar_dev);// cuda_release(dvar_map);
			cuda_release(aID_dev);// cuda_release(aID_map);
			cuda_release(cIndx_dev);// cuda_release(cIndx_map);
			cuda_release(bound_dev);// cuda_release(bound_map);
			cuda_release(bound_14_dev);// cuda_release(bound_14_map);
			cuda_release(b_relship_dev);// cuda_release(b_relship_map);
			cuda_release(b14_relship_dev);// cuda_release(b14_relship_map);
			cuda_release(abound_dev); // cuda_release(abound_map);
			cuda_release(r_dev);// cuda_release(r_map);
			cuda_release(ir_dev);// cuda_release(ir_map);
			cuda_release(cmm_dev);// cuda_release(cmm_map);
#undef cuda_release
		};
		~GPU_AtCELL() {release();};
		__device__ __inline__ int* device_cubic(int ix, int iy, int iz) {
			int ic = ix * Nyz + iy * Nz + iz;
			return cmm_dev + ic * (NC + 1);
		};
		__host__ __device__ __inline__ int* host_cubic(int ix, int iy, int iz) {
			int ic = ix * Nyz + iy * Nz + iz;
			return cmm_map + ic * (NC + 1);
		};

		__device__ __inline__ int* device_abound(int ia) {
			return abound_dev + ia * (Nm_b + 1);
		};
		__host__ __device__ __inline__ int* host_abound(int ia) {
			return abound_map + ia * (Nm_b + 1);
		};
	};


	class GPU_EAtCELL : public GPU_AtCELL {
	public:
		double mu_surf[3]; // surface dipole, which including the charge effect
		double ind_mu_surf[3]; // surface dipole from atomic induced dipole

		double *c_dev, *c_map; // charge with dimension Na
		double *mu_dev, *mu_map; // dipole with dimension Na * 3

		//Results:
		double *E_dev, *E_map; // the electrostatic field of each atoms, dimension -- Na * 3
		double *Fe_dev, *Fe_map; // force from electrostatic interaction, dimension -- Na * 3

		// temperary result:
		double *E_14_dev, *E_14_map; // electrostatic field from the exclusion part, e.g. 1-4, actually the field is on the first atom defined in the EwdRel_exclude, a1_indx of each pair
	                    // the field on another atom, a2_indx, should be negative values, or inverse direction
	                    // dimension -- Nb14 * 3
		double *Fe_14_dev, *Fe_14_map; // electrostatic force from the exclusion part, e.g. 1-4, actually the field is on the first atom defined in the EwdRel_exclude, a1_indx of each pair
	                     // the field on another atom, a2_indx, should be negative values, or inverse direction
	                     // dimension -- Nb14 * 3
		double *E_b_dev, *E_b_map; // electrostatic field between bounded atoms, actually the field is on the first atom defined in the bound
	                   // the field on another atom, a2_indx, should be negative values, or inverse direction
	                   // dimension -- Nb * 3
		double *Fe_b_dev, *Fe_b_map; // electrostatic force from the exclusion part, e.g. 1-4, actually the field is on the first atom the bound
	                    // the field on another atom, a2_indx, should be negative values, or inverse direction
	                    // dimension -- Nb * 3


		double *ve_dev, *ve_map;  // virial coefficient, dimension -- Na * 3
		double *U_estat_dev, *U_estat_map; // dimension -- Na
		double *ve_b_dev, *ve_b_map; // virial coefficient from bound electrostatic interaction which will be excluded, dimension -- Nb * 3
		double *Ub_estat_dev, *Ub_estat_map; // dimension -- Na
	
		double *ve_14_dev, *ve_14_map; // virial coefficient from electrostatic interaction which will be excluded, dimension -- Nb14 * 3
		double *U14_estat_dev, *U14_estat_map; // dimension -- Na

		// some parameters for calculation controlling
		// ## calculating electrostatic interaction 
		char bDipole; // using dipole of atom ? [0/1]
		char bCalE; // 0 -- do not do calculation
			        // 1 -- calculate electric field
				    // 2 -- calculating force, and energy
					// 3 -- calculating electric field, force and energy. But is this necessary?
		// Ewald-Sum parameters
		// ...
		int iExIndDipole;

		// real part of EwaldSum
		EwaldSumReal_ControlVar esv;
		gvar mgvar;

		// neighbors of each atom
		// the memory locates on device in case using GPU, or NULL in case using CPU
		bool bSRInitialized; // whether neighbors are initialized
		int n_srnb; // maximum neighbor atoms in short-range interaction for each atom, which is estimated from the cut-out distance
		int *srnb_dev; // the neighbors of all atoms, dimension (n_srnb + 1) * Na. For each atom, first element is the number of neighbors, then the index of atoms
		DOUBLE *dr_dev; // the relative coordinations from the center atom -- dr[3], and |R| -- the distance, dimension -- n_srnb * 4 * Na

#ifdef __GPU__
		__host__ __device__ int srIndx(int aIndx) {return aIndx * (n_srnb + 1);};
		__host__ __device__ int drIndx(int aIndx) {return aIndx * n_srnb * 4;};
#endif

		GPU_EAtCELL() {
			c_dev = NULL; c_map = NULL;
			mu_dev = NULL; mu_map = NULL;
			E_dev = NULL; E_map = NULL;
			Fe_dev = NULL; Fe_map = NULL;
			E_14_dev = NULL; E_14_map = NULL;
			Fe_14_dev = NULL; Fe_14_map = NULL;
			E_b_dev = NULL; E_b_map = NULL;
			Fe_b_dev = NULL; Fe_b_map = NULL;
			ve_dev = NULL; ve_map = NULL;
			U_estat_dev = NULL; U_estat_map = NULL;
			ve_b_dev = NULL; ve_b_map = NULL;
			Ub_estat_dev = NULL; Ub_estat_map = NULL;
			ve_14_dev = NULL; ve_14_map = NULL;
			U14_estat_dev = NULL; U14_estat_map = NULL;

			srnb_dev = NULL; dr_dev = NULL; n_srnb = 0; bSRInitialized = false;
		};

		void release() {
#define cuda_release(v) if (v != NULL) {cudaFree(v); v = NULL;}
			cuda_release(c_dev);
			cuda_release(mu_dev);
			cuda_release(E_dev);
			cuda_release(Fe_dev);
			cuda_release(ve_dev);
			cuda_release(U_estat_dev);

			cuda_release(E_14_dev);
			cuda_release(Fe_14_dev);
			cuda_release(ve_14_dev);
			cuda_release(U14_estat_dev);

			cuda_release(E_b_dev);
			cuda_release(Fe_b_dev);
			cuda_release(ve_b_dev);
			cuda_release(Ub_estat_dev);

			cuda_release(srnb_dev); cuda_release(dr_dev); n_srnb = 0; bSRInitialized = false;

			GPU_AtCELL::release();
#undef cuda_release
		};
	};

	class GPU_SR_AtCELL : public GPU_AtCELL {
	public:
		char *sr_type_dev, *sr_type_map; // type of short-range interaction, dimension -- Na
			                 // 0x01 -- L-J interaction, 0x02 -- given profile

		// number of parameters for short-range interaction, of each atom, e.g. eps and rLJ for Lennard-Jones interaction
		// or, pair-wise interaction would require something to describe it
		// in general, GPU calculation will get better performance, if it is not 
		int sr_mode; // 0 -- standard LJ, 1 -- alternative LJ
		int npars; 
		double *sr_par_dev;

		char *bShocked_dev, *bShocked_map; // whether the atom is too close to other atom, be shocked, dimension -- Na // 0x01 -- true, 0x00 -- false
		double *Fsr_dev, *Fsr_map; // the force from short-range interaction for each atoms, dimension -- Na * 3
		double *Fsr_b_dev, *Fsr_b_map; // short-range force from the exclusion part, bounded atom, actually the field is on the first atom the bound
			                 // the field on another atom, a2_indx, should be negative values, or inverse direction
				             // dimension -- Nb * 3
		double *Fsr_14_dev, *Fsr_14_map; // short-range force from the exclusion part, bounded atom, actually the field is on the first atom the bound
		                      // the field on another atom, a2_indx, should be negative values, or inverse direction
			                  // dimension -- Nb14 * 3
		double *vsr_dev, *vsr_map;     // virial coefficient, dimension -- Na * 3
		double *U_sr_dev, *U_sr_map; // dimension -- Na
		double *vsr_b_dev, *vsr_b_map; // virial coefficient from bound short-range interaction which will be excluded, dimension -- Nb * 3
		double *Ub_sr_dev, *Ub_sr_map; // dimension -- Nb
		double *vsr_14_dev, *vsr_14_map; // virial coefficient from electrostatic interaction which will be excluded, dimension -- Nb14 * 3
		double *U14_sr_dev, *U14_sr_map; // Nb14

		// some parameters for calculation controlling
		char bCalSR; // calculating short-range interaction ? [0/1]

		// host member function
		void reset_sr();

		gvar mgvar;

		GPU_SR_AtCELL() {
			sr_type_dev = NULL; sr_type_map = NULL;
			sr_par_dev = NULL;

			bShocked_dev = NULL; bShocked_map = NULL;
			Fsr_dev = NULL; Fsr_map = NULL;
			Fsr_b_dev = NULL; Fsr_b_map = NULL;
			Fsr_14_dev = NULL; Fsr_14_map = NULL;
			vsr_dev = NULL; vsr_map = NULL;
			U_sr_dev = NULL; U_sr_map = NULL;
			vsr_b_dev = NULL; vsr_b_map = NULL;
			Ub_sr_dev = NULL; Ub_sr_map = NULL;
			vsr_14_dev = NULL; vsr_14_map = NULL;
			U14_sr_dev = NULL; U14_sr_map = NULL;
		};

		void release() {
#define cuda_release(v) if (v != NULL) {cudaFree(v); v = NULL;}
			cuda_release(sr_type_dev); sr_type_map = NULL;
			cuda_release(sr_par_dev); sr_par_dev = NULL;

			cuda_release(bShocked_dev); bShocked_map = NULL;
			cuda_release(Fsr_dev); Fsr_map = NULL;
			cuda_release(vsr_dev); vsr_map = NULL;
			cuda_release(U_sr_dev); U_sr_map = NULL;

			cuda_release(Fsr_14_dev); Fsr_14_map = NULL;
			cuda_release(vsr_14_dev); vsr_14_map = NULL;
			cuda_release(U14_sr_dev); U14_sr_map = NULL;

			cuda_release(Fsr_b_dev); Fsr_b_map = NULL;
			cuda_release(vsr_b_dev); vsr_b_map = NULL;
			cuda_release(Ub_sr_dev); Ub_sr_map = NULL;

			GPU_AtCELL::release();
#undef cuda_release
		};
	};

	__host__ bool relocate_atom(GPU_AtCELL *acell_host, GPU_AtCELL *acell_dev, bool bUpdateHost);
	__host__ void update_charge_dipole(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev);
	__host__ void reset_efield(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev);

	__host__ void check_neighbors(GPU_EAtCELL *acell_host, GPU_EAtCELL *acell_dev);
	__host__ void update_neighbors_status(GPU_EAtCELL *acell_dev, bool status);

	__host__ void GPU_LJ_interact(GPU_SR_AtCELL *acell_host, GPU_SR_AtCELL *acell_dev);

} // end of namespace _gpu_md_


#endif // __GPU__
