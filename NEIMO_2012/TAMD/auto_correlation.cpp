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
extern void molecule_check_periodic_cell(CGMM_MD_CELL &mmcell);

extern bool get_monomer_uid(char *monomer, int &uid);

extern bool read_given_id(char *fname, char* title, ARRAY<short>& id_cluster);

#include "distribution.h"
#include "auto_correlation.h"

extern bool ReadMoleculeStruct(ifstream &in, MMOLECULE &mm, bool full = false);
extern bool ReadSolidMoleculeStruct(ifstream &in, SMOLECULE &m, bool full = false);
extern bool ReadPointMoleculeStruct(ifstream &in, PMOLECULE &m, bool full = false);

extern void calTransSpatialMomentum(MMOLECULE &mm, VECTOR3 &P);

namespace _velocity_correlation_ {

	using namespace _distribution_;

	bool read_velocity(ifstream &in, int isave, MMOL_MD_CELL &mdcell, int mtype, int mindx, int cindx, SVECTOR<float, 3> &v) {
		char mol_title[250] = "\0", buffer[256] = "\0", msg[256] = "\0";
		if (!search(in, "[MD]", isave, buffer)) return false;
		if (mtype == 0) strcpy(mol_title, "MM");
		else if (mtype == 1) strcpy(mol_title, "SM");
		else if (mtype == 2) strcpy(mol_title, "PM");
		else {
			sprintf(msg, "Unknown molecule type [0 -- MM, 1 -- SM, 2 -- PM]: %d", mtype); show_log(msg, true);
			return false;
		}

		if (!search(in, mol_title, mindx, buffer, "[MD]")) {
			sprintf(msg, "can not find %s %d", mol_title, mindx); show_log(msg, true);
			return false;
		}

		VECTOR3 P;
		if (mtype == 0) {
			if (!ReadMoleculeStruct(in, mdcell.mm.m[mindx], true)) return false;
			calClusterVelocity0(mdcell.mm.m[mindx]);
			calTransSpatialMomentum(mdcell.mm.m[mindx], P);
			v.v[0] = mdcell.mm.m[mindx].cluster[cindx].mcDyn->v.v[0];
			v.v[1] = mdcell.mm.m[mindx].cluster[cindx].mcDyn->v.v[1];
			v.v[2] = mdcell.mm.m[mindx].cluster[cindx].mcDyn->v.v[2];
		}
		else if (mtype == 1) {
			if (!ReadSolidMoleculeStruct(in, mdcell.sm.m[mindx], true)) return false;
			v.v[0] = mdcell.sm.m[mindx].c->bkm->V0.v[3];
			v.v[1] = mdcell.sm.m[mindx].c->bkm->V0.v[4];
			v.v[2] = mdcell.sm.m[mindx].c->bkm->V0.v[5];
		}
		else if (mtype == 2) {
			if (!ReadPointMoleculeStruct(in, mdcell.pm.m[mindx], true)) return false;
			v.v[0] = mdcell.pm.m[mindx].v.v[0];
			v.v[1] = mdcell.pm.m[mindx].v.v[1];
			v.v[2] = mdcell.pm.m[mindx].v.v[2];
		}
		return true;
	}

	bool read_vhist(ifstream &in, int &isave, MMOL_MD_CELL &mdcell, int mtype, int mindx, int cindx, CHAIN< SVECTOR<float, 3> > **ch, CHAIN< SVECTOR<float, 3> > **tch) {
		CHAIN< SVECTOR<float, 3> > *vch = NULL;
		rewind(in);
		SVECTOR<float, 3> v;
		int n = isave;

		while (1) {
			if (in.eof()) break;
			if (!read_velocity(in, isave, mdcell, mtype, mindx, cindx, v)) break;
			vch = new CHAIN< SVECTOR<float, 3> >; vch->p = new SVECTOR<float, 3>;
			memcpy(vch->p->v, v.v, sizeof(float) * 3);

			if (*ch == NULL) *ch = vch;

			if (*tch == NULL) *tch = vch;
			else {
				(*tch)->next = vch; *tch = vch;
			}
			n++; isave++;
		}
		if (n == 0) return false;
		else {isave = n; return true;}
	}

	bool read_vhist(char *ftitle, int istart, MMOL_MD_CELL &mdcell, int mtype, int mindx, int cindx, BUFFER< SVECTOR<float, 3> > &vhist) {
		char fname[256] = "\0";
		int isave = istart, ifile = 0;
		vhist.release();
		CHAIN< SVECTOR<float, 3> > *ch = NULL, *tch = NULL;
		ifstream *in = NULL;

		while (1) {
			sprintf(fname, "%s-%d.kin", ftitle, ifile);
			in = new ifstream;
			in->open(fname);
			if (!in->is_open()) {delete in; in = NULL; break;}
			if (!read_vhist(*in, isave, mdcell, mtype, mindx, cindx, &ch, &tch)) {
				in->close(); delete in; in = NULL;
				break;
			}
			in->close(); delete in; in = NULL;
			ifile +=1;
		}
		if (isave == istart) {
			show_log("NO md procedure was read", true);
			return false;
		}
		int nconf = isave - istart;
		vhist.init(nconf);
		tch = ch;
		SVECTOR<float, 3> *bf = NULL;
		for (int i = 0; i < nconf; i++) {
			bf = vhist.get_data(i);
			if (bf == NULL) {
				show_log("strange: vhist buffer is not avaliable", true); break;
			}
			if (tch == NULL) {
				show_log("strange: original vhist buffer is not avaliable", true); break;
			}
			memcpy(bf->v, tch->p->v, sizeof(float) * 3);
			tch = tch->next;
		}
		release_chain< SVECTOR<float, 3> >(&ch, true); tch = NULL;
		return true;
	}

	void auto_correlation_float3(BUFFER< SVECTOR<float, 3> > &vhist, DISTRIBT<int, double> &dt) {
		int i, j;
		SVECTOR<float, 3> *v1 = NULL, *v2 = NULL;
		float v_cor = 0;
		char msg[256] = "\0";

		for (i = 0; i < vhist.n_size; i++) {
			v1 = vhist.get_data(i);
			if (::COMMAND_MD_STOP) break;
			for (j = 0; j < vhist.n_size; j++) {
				if ((j - i) > dt.xt) break;
				v2 = vhist.get_data(j);
				v_cor = v2->v[0] * v1->v[0] + v2->v[1] * v1->v[1] + v2->v[2] * v1->v[2];
				dt.sample(j - i, v_cor, true);
			}

			if (i / 1000 * 1000 == i) {
				sprintf(msg, " %d", i); show_log(msg, false);
			}
		}
		sprintf(msg, " %d", i); show_log(msg, true);
	}

} // end of namespace _velocity_correlation_

bool get_md_config(char *parfname, char *mdpar_fname, int &istart) {
		char buffer[256] = "\0", msg[256] = "\0";
		if (!read_item(parfname, "[MD_CONFIG]", buffer)) {
			sprintf(msg, "failure to get MD parameter filename from %s, defined after [MD_CONFIG]", parfname); return false;
		}
		if (sscanf(buffer, "%s %d", mdpar_fname, &istart) == 2 && strlen(mdpar_fname) > 0) return true;
		else return false;
	};

extern "C" bool velocity_auto_correlation(char *parfname) {
	using namespace _velocity_correlation_;
	::mlog.init("velocity_auto_corr.log");
	::COMMAND_MD_STOP = false;

		char mdpar_fname[256] = "\0", kin_ftitle[256] = "\0";
		char msg[256] = "\0", msg1[256] = "\0", buffer[256] = "\0";

		int istart = 0, i;
		if (!get_md_config(parfname, mdpar_fname, istart)) return false;
		if (!read_item(mdpar_fname, "[MD_TITLE]", buffer)) return false;
		if (sscanf(buffer, "%s", kin_ftitle) != 1 || strlen(kin_ftitle) == 0) {
			sprintf(msg, "MD output file's title is not given in %s", mdpar_fname); show_msg(msg);
			return false;
		}


		if (!read_item(mdpar_fname, "[MD_CELL]", msg)) return false;
		int MD_MODE = 0;
		if (sscanf(msg, "%d", &MD_MODE) != 1 || (MD_MODE != MULTI_MM_FREE_CELL_MD && MD_MODE != MULTI_MM_PERIDIC_CELL_MD
			&& MD_MODE != MULTI_CGMM_FREE_CELL_MD && MD_MODE != MULTI_CGMM_PERIDIC_CELL_MD)) { // MultiMol ?
			sprintf(msg, "md is not working on multi-molecule");
			show_msg(msg); return false;
		}

		set_md_mode(MD_MODE);
		if (MD_MODE == MULTI_MM_FREE_CELL_MD || MD_MODE == MULTI_MM_PERIDIC_CELL_MD) set_job(GET_MULTI_MM_MD_PROC);
		else if (MD_MODE == MULTI_CGMM_FREE_CELL_MD || MD_MODE == MULTI_CGMM_PERIDIC_CELL_MD) set_job(GET_MULTI_CGMM_MD_PROC);

		bool status = construct_macromolecule(mdpar_fname);
		if (!status) {sprintf(msg, "failure to initilize multi-mol system"); show_msg(msg); return false;}
		for (i = 0; i < mm_md_cell.mm.n; i++) {mm_md_cell.mm.m[i].setup_kin_dyn(); mm_md_cell.mm.m[i].setup_mcDyn();}
		for (i = 0; i < mm_md_cell.sm.n; i++) mm_md_cell.sm.m[i].c->setup_kin_dyn();
		for (i = 0; i < mm_md_cell.pm.n; i++) mm_md_cell.pm.m[i].setup_dyn();


		char output_ftitle[256] = "\0", output_fname[256] = "\0";
		if (read_item(parfname, "[OUTPUT]", buffer)) {
			if (sscanf(buffer, "%s", output_ftitle) != 1 || strlen(output_ftitle) == 0) {
				sprintf(msg, "file's title is not given after [OUTPUT] in %s", output_fname);
				show_log(msg, true); return false;
			}
		}
		else {
			sprintf(msg, "file's title is not given after [OUTPUT] in %s", output_fname);
			show_log(msg, true); return false;
		}

		ofstream *out = NULL;
		ifstream in;
		in.open(parfname);
		if (!search(in, "[AUTO_CORRELATIONS]", buffer)) {
			sprintf(msg, "[AUTO_CORRELATIONS] is not define in %s", parfname); show_msg(msg); in.close(); return false;
		}
		int mtype = 0, mindx = 0, cindx = 0;
		int t0 = 0, t1 = 10, idistrbt = 0;
		BUFFER< SVECTOR<float, 3> > vhist(5000);
		DISTRIBT<int, double> dt;
		while (1) {
			sprintf(output_fname, "%s-%d.vcor.dat", output_ftitle, idistrbt);
			in.getline(buffer, 250);
			if (strlen(buffer) == 0) break;
			if (sscanf(buffer, "%d %d %d %d", &mtype, &mindx, &cindx, &t1) != 4) {
				sprintf(msg, "Unknown format for velocity auto-correlation: %s", buffer); show_msg(msg);
				show_msg("        format: [mol_type: 0 -- MM, 1 -- SM, 2 -- PM] [mol index] [cluster index]  [step max]"); 
				in.close(); return false;
			}
			dt.set_distribt_x(0, t1, t1);
			dt.reset_distribt();
			if (read_vhist(kin_ftitle, istart, mm_md_cell, mtype, mindx, cindx, vhist)) {
				sprintf(msg, "%d saved md steps were read!", vhist.n_size); show_log(msg, true);
				show_log("calculating velocity auto-correlation ...", true);

				auto_correlation_float3(vhist, dt);
				out = new ofstream;
				out->open(output_fname);
				if (!out->is_open()) {
					sprintf(msg, "failure to open file %s to ouput distribution %d", idistrbt); show_msg(msg);
					in.close(); return false;
				}
				for (i = 0; i < dt.x.N; i++) {
					if (dt.sn.v[i] > 0) {
						dt.eb.v[i] = sqrt(float(dt.sn.v[i])) / dt.sn.v[i];
						dt.y.v[i] /= dt.sn.v[i];
						dt.eb.v[i] *= dt.y.v[i];
					}
					(*out)<<i<<"  "<<dt.y.v[i]<<"  "<<dt.eb.v[i]<<endl;
				}
				out->close(); delete out; out = NULL;
				sprintf(msg, "correlation %d is saved in %s", idistrbt, output_fname);
				show_log(msg, true);
			}
			idistrbt +=1;
		}
		in.close();
		::mlog.close();
		return true;
	}


