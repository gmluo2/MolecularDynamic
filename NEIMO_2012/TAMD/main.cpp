#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
#include "show.h"
#define _MAIN_PROG_ 1
#include "def.h"
#include "vector.h"
#include "bound.h"
#include "ZMatrix.h"

#include "MM.h"
#include "Mol.h"

#include "var.h"

#if _SYS_ == _WINDOWS_SYS_
#include "..\export.h"
#elif _SYS_ == _LINUX_SYS_
#include "export.h"
#endif

extern "C" bool dipole_distribution(char *parfname);
extern "C" bool orientation_distribution(char *parfname);
extern "C" bool relative_orientation_distribution(char *parfname);
extern "C" bool dipole_local_orientation(char *parfname);
extern "C" bool cluster_torsion_angle(char *parfname);
extern "C" bool MultiSaveConfigXYZ(char *parfname);

extern "C" bool single_mm_neighbor_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_mm_neighbor_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool single_mm_non_neighbor_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_mm_non_neighbor_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_mm_non_neighbor_sm_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_mm_nsm_rdf_single_chain_block_copolymer(char *parfname);

extern "C" bool single_cgmm_neighbor_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_cgmm_neighbor_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool single_cgmm_non_neighbor_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_cgmm_non_neighbor_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_cgmm_non_neighbor_sm_rdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_cgmm_nsm_rdf_single_chain_block_copolymer(char *parfname);

extern "C" bool single_cgmm_neighbor_crdf_single_chain_block_copolymer(char *parfname);

extern "C" bool single_mm_neighbor_erdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_mm_neighbor_erdf_single_chain_block_copolymer(char *parfname);
extern "C" bool single_mm_non_neighbor_erdf_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_mm_non_neighbor_erdf_single_chain_block_copolymer(char *parfname);

extern "C" bool single_mm_dipole_distrbt_single_chain_block_copolymer(char *parfname);
extern "C" bool multi_mm_dipole_distrbt_single_chain_block_copolymer(char *parfname);

extern "C" bool atomic_2d_correlation_3d_single_thread(char *parfname);
extern "C" bool cluster_spatial_distribution(char *parfname);
extern "C" bool velocity_auto_correlation(char *parfname);

extern bool search(ifstream &in, char *title, char *buffer);
extern void show_msg(char *msg, bool endline);

bool get_md_cell_type(char *mdpar_fname, int& cell_type, char *errmsg) {
	char buffer[256] = "\0";
	ifstream in;
	in.open(mdpar_fname);
	if (!in.is_open()) {sprintf(errmsg, "can not open file %s", mdpar_fname); return false;}
	if (!search(in, "[MD_CELL]", buffer)) {
		sprintf(errmsg, "can not find [MD_CELL] item in %s", mdpar_fname); in.close(); return false;
	}
	in.getline(buffer, 250);
	if (sscanf(buffer, "%d", &cell_type) != 1) {
		sprintf(errmsg, "can not find cell type from %s", buffer); in.close(); return false;
	}
	in.close();
	return true;
}

void show_usage() {
	char buffer[256] = "\0";
	show_msg("NEIMO -help");
	show_msg("NEIMO -dipole [par file]");
	show_msg("NEIMO -mort [par file]");
	show_msg("NEIMO -rort [par file]");
	show_msg("NEIMO -local_dipole [par file]");
	show_msg("NEIMO -ta [par file]"); 
	show_msg("NEIMO -xyz [par file]");
	show_msg("NEIMO -sm_rdf_neighbor [par file]");
	show_msg("NEIMO -sm_rdf_non_neighbor [par file]");
	show_msg("NEIMO -mm_rdf_neighbor [par file]");
	show_msg("NEIMO -mm_rdf_non_neighbor [par file]");
	show_msg("NEIMO -mm_rdf_non_neighbor_sm [par file]");
	show_msg("NEIMO -mm_rdf_non_neighbot_nsm [par file]");
	show_msg("NEIMO -scgm_rdf_neighbor [par file]");
	show_msg("NEIMO -scgm_rdf_non_neighbor [par file]");
	show_msg("NEIMO -mcgm_rdf_neighbor [par file]");
	show_msg("NEIMO -mcgm_rdf_non_neighbor [par file]");
	show_msg("NEIMO -mcgm_rdf_non_neighbor_sm [par file]");
	show_msg("NEIMO -mcgm_rdf_non_neighbor_nsm [par file]");
	show_msg("NEIMO -scgm_crdf_neighbor [par file]");
	show_msg("NEIMO -sm_erdf_neighbor [par file]");
	show_msg("NEIMO -sm_erdf_non_neighbor [par file]");
	show_msg("NEIMO -mm_erdf_neighbor [par file]");
	show_msg("NEIMO -mm_erdf_non_neighbor [par file]");
	show_msg("NEIMO -sm_vm_dipole_distrbt [par file]");
	show_msg("NEIMO -mm_vm_dipole_distrbt [par file]");
	show_msg("NEIMO : use MD.INI for MD simulation");
}

#if _MPI_ == 1
#ifdef SEEK_SET
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#endif
#include "mpi.h"
#endif

extern "C" int run(int argc, char **argv) {
	if (argc >= 2 && (strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "-h") == 0)) {
		show_usage(); return 0;
	}

#if _MPI_ == 1
	int i = 0, id = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &::nMPI);
	MPI_Comm_rank(MPI_COMM_WORLD, &::mpi_id);
	MPI_Get_processor_name(::proc_name[::mpi_id], &i);
	for (i = 0; i < ::nMPI; i++) {
		id = ::mpi_id;
		MPI_Bcast((void*)(&id), 1, MPI_INT, i, MPI_COMM_WORLD);
		MPI_Bcast((void*)(::proc_name[id]), 45, MPI_CHAR, id, MPI_COMM_WORLD);
	}
#endif

	if (!init_env()) {
		show_msg(::errmsg, true); //return 0;
	}

	if (argc > 2) {
		if (strcmp(argv[1], "-dipole") == 0) {
			dipole_distribution(argv[2]);
		}
		else if (strcmp(argv[1], "-mort") == 0) {
			orientation_distribution(argv[2]);
		}
		else if (strcmp(argv[1], "-rort") == 0) {
			relative_orientation_distribution(argv[2]);
		}
		else if (strcmp(argv[1], "-local_dipole") == 0) {
			dipole_local_orientation(argv[2]);
		}
		else if (strcmp(argv[1], "-ta") == 0) {
			cluster_torsion_angle(argv[2]);
		}
		else if (strcmp(argv[1], "-xyz") == 0) {
			MultiSaveConfigXYZ(argv[2]);
		}
		else if (strcmp(argv[1], "-sm_rdf_neighbor") == 0) {
			single_mm_neighbor_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-sm_rdf_non_neighbor") == 0) {
			single_mm_non_neighbor_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mm_rdf_neighbor") == 0) {
			multi_mm_neighbor_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mm_rdf_non_neighbor") == 0) {
			multi_mm_non_neighbor_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mm_rdf_non_neighbor_sm") == 0) {
			multi_mm_non_neighbor_sm_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mm_rdf_non_neighbor_nsm") == 0) {
			multi_mm_nsm_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-scgm_rdf_neighbor") == 0) {
			single_cgmm_neighbor_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-scgm_rdf_non_neighbor") == 0) {
			single_cgmm_non_neighbor_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mcgm_rdf_neighbor") == 0) {
			multi_cgmm_neighbor_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mcgm_rdf_non_neighbor") == 0) {
			multi_cgmm_non_neighbor_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mcgm_rdf_non_neighbor_sm") == 0) {
			multi_cgmm_non_neighbor_sm_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mcgm_rdf_non_neighbor_nsm") == 0) {
			multi_cgmm_nsm_rdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-scgm_crdf_neighbor") == 0) {
			single_cgmm_neighbor_crdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-sm_erdf_neighbor") == 0) {
			single_mm_neighbor_erdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-sm_erdf_non_neighbor") == 0) {
			single_mm_non_neighbor_erdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mm_erdf_neighbor") == 0) {
			multi_mm_neighbor_erdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mm_erdf_non_neighbor") == 0) {
			multi_mm_non_neighbor_erdf_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-sm_vm_dipole_distrbt") == 0) {
			single_mm_dipole_distrbt_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-mm_vm_dipole_distrbt") == 0) {
			multi_mm_dipole_distrbt_single_chain_block_copolymer(argv[2]);
		}
		else if (strcmp(argv[1], "-cor_2d") == 0) {
			atomic_2d_correlation_3d_single_thread(argv[2]);
		}
		else if (strcmp(argv[1], "-cluster_distribt") == 0) {
			cluster_spatial_distribution(argv[2]);
		}
		else if (strcmp(argv[1], "-velocity_auto_cor") == 0) {
			velocity_auto_correlation(argv[2]);
		}
		clear_env(); return 1;
	}

	char mdpar_fname[256] = "MD.INI";
	int Job = UNKNOWN_JOB;
	if (!get_md_cell_type(mdpar_fname, Job, ::errmsg)) {
		show_msg(::errmsg); clear_env(); return 0;
	}

	if (Job <= 10) {
		set_job(Job);
		set_md_mode(Job);
		set_show_infor(0, 0, 40); // save xyz every 10 loops
	}
	else {
		sprintf(errmsg, "Molecule Dynamic # %d is not realized yet", Job); 
		show_msg(errmsg); clear_env(); return 0;
	}

#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 1
	sprintf(errmsg, "electrostatic interaction between atoms is ignored");
	show_log(errmsg, true);
#endif

	if (!construct_macromolecule(mdpar_fname)) {clear_env(); return 0;}
	MD_SIM(mdpar_fname);

	clear_env();
#if _MPI_ == 1
	MPI_Finalize();
#endif
	return 1;
}

bool check_data_type() {
	bool status = true;
	if (sizeof(double) != SIZE_DOUBLE) {show_log("Important: SIZE_DOUBLE is not right!", true); status = false;}
	if (sizeof(int) != SIZE_INT) {show_log("Important: SIZE_INT is not right!", true); status = false;}
	if (sizeof(float) != SIZE_FLOAT) {show_log("Important: SIZE_FLOAT is not right!", true); status = false;}
	return status;
}

#if _SYS_ == _LINUX_SYS_
int main(int argc, char **argv) {
	if (!check_data_type()) return 0;
	return run(argc, argv);
}
#endif
