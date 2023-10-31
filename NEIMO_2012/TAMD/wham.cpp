#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "def.h"
#include "show.h"
#include "dstruct.h"
#include "vector.h"

#include "prof.h"
#include "read.h"

#include "wham.h"

namespace _WHAM_ {
	bool WHAM::init_unbaised_prof() {
		// initialize the backup buffer for F of the windows
		Fb.set(fi.n);

		// initialize the whole profile
		f.release();
		int N = 0, n;
		int iprof, ipt;

		char msg[256] = "\0";

		float dx0 = fi.m[0].x.v[1] - fi.m[0].x.v[0];
		float dx, ddx;
		
		for (iprof = 0; iprof < fi.n; iprof++) {
			for (ipt = 0; ipt < fi.m[iprof].x.N - 1; ipt++) {
				dx = fi.m[iprof].x.v[ipt + 1] - fi.m[iprof].x.v[ipt];
				ddx = dx / dx0 - 1; ddx = (ddx >= 0 ? ddx : -ddx);
				if (ddx > 0.01f) {
					sprintf(msg, "biased window %d has different step [%f] of coordinate @ point %d from standard step [%f]", iprof, dx, ipt, dx0);
					show_log(msg, true); return false;
				}
			}
		}

		// assuming the profiles are ordered with the coordinate
		float xf = fi.m[0].x.v[0], xt = fi.m[fi.n - 1].x.t[0];
		N = int((xt - xf) / dx0 + 0.5f);
		f.set(N + 1);
		for (ipt = 0; ipt < f.x.N; ipt++) {
			f.x.v[ipt] = xf + ipt * dx0; f.s.v[ipt] = 0; f.y.v[ipt].v = 0; 
		}

		CHAIN<PtIndx> *pt = NULL, *tail = NULL, *tch = NULL;
		int i;
		bool status = true;
		for (i = 0; i < f.x.N; i++) {
			f.y.v[i].indx.release(); n = 0;
			for (iprof = 0; iprof < fi.n; iprof++) {
				for (ipt = 0; ipt < fi.m[iprof].x.N; ipt++) {
					dx = fi.m[iprof].x.v[ipt] - f.x.v[i]; dx = dx / dx0; ddx = FABS(dx);
					if (ddx < 0.1) {
						tch = new CHAIN<PtIndx>; tch->p = new PtIndx;
						tch->p->iw = iprof; tch->p->ipt = ipt;
						if (tail == NULL) pt = tch;
						else tail->next = tch;
						tail = tch;
						n++;
					}
					if (dx > 2) break; // the further points will be larger than f.x.v[n] because it is assumed that the profile has x ordered
				}
			}
			if (pt == NULL) {
				status = false; // a point has 
				sprintf(msg, "no point exists in the biased distribution when x = %f", f.x.v[i]); show_log(msg, true);
			}
			else {
				f.y.v[i].indx.SetArray(n); 
				tch = pt; ipt = 0;
				while (tch != NULL) {
					f.y.v[i].indx.m[ipt].iw = tch->p->iw;
					f.y.v[i].indx.m[ipt].ipt = tch->p->ipt;
					tch = tch->next; ipt++;
				}
				release_chain<PtIndx>(&pt, true);
			}
			pt = NULL; tail = NULL; tch = NULL;
		}
		return status;
	}

	void WHAM::construct_unbaised_prof() {
		int iw, ipt, i, j;
		PtIndx *indx = NULL;
		float v_nc = 0;
		for (i = 0; i < f.x.N; i++) {
			f.y.v[i].v = 0; f.s.v[i] = 0;
			if (f.y.v[i].indx.n == 0) continue;
			v_nc = 0;
			for (j = 0; j < f.y.v[i].indx.n; j++) {
				indx = f.y.v[i].indx.m + j;
				iw = indx->iw; ipt = indx->ipt;
				v_nc += exp(-fi.m[iw].wb.v[ipt] + fi.m[iw].F) * fi.m[iw].Nc;
			}

			for (j = 0; j < f.y.v[i].indx.n; j++) {
				indx = f.y.v[i].indx.m + j;
				iw = indx->iw; ipt = indx->ipt;
				f.y.v[i].v += fi.m[iw].y.v[ipt] * exp(-fi.m[iw].wb.v[ipt] + fi.m[iw].F) * fi.m[iw].Nc / v_nc;
			}
		}
	}

	void WHAM::backup_F() {
		int iw;
		for (iw = 0; iw < fi.n; iw++) Fb.v[iw] = fi.m[iw].F;
	}

	void WHAM::calculate_F() {
		int iw, ipt;
		float v_nc;
		int nf = 0, nt = 0;
		float xf0 = f.x.v[0], dx0 = f.x.v[1] - f.x.v[0];
		for (iw = 0; iw < fi.n; iw++) {
			nf = int((fi.m[iw].x.v[0] - xf0) / dx0);
			nt = int((fi.m[iw].x.t[0] - xf0) / dx0 + 0.9);
			v_nc = 0;
			for (ipt = 0; ipt < fi.m[iw].wb.N; ipt++) {
				v_nc += exp(-fi.m[iw].wb.v[ipt]) * f.y.v[ipt + nf].v;
			}
			fi.m[iw].F = -log(v_nc);
		}

		for (iw = 0; iw < fi.n; iw++) {
			fi.m[iw].F -= fi.m[fi.n - 1].F;
		}
	}

	float WHAM::max_diff_F() {
		float dvmax = 0, dv;
		int iw;
		dvmax = fi.m[0].F - Fb.v[0]; dvmax = FABS(dvmax);
		for (iw = 1; iw < fi.n; iw++) {
			dv = fi.m[iw].F - Fb.v[iw]; dv = FABS(dv);
			if (dvmax < dv) dvmax = dv;
		}
		return dvmax;
	}

	bool IterateWHAM(WHAM &wham, float dF_max, int max_loops) {
		bool status = true;
		float dF = 0;
		int iloop = 0, i;
		float c = 0.6f;
		for (i = 0; i < wham.fi.n; i++) wham.fi.m[i].F -= wham.fi.m[wham.fi.n - 1].F;
		while (1) {
			wham.backup_F();
			wham.construct_unbaised_prof();
			wham.calculate_F();
			dF = wham.max_diff_F();
			if (dF < dF_max) {status = true; break;}
			else {
				for (i = 0; i < wham.fi.n; i++) wham.fi.m[i].F = wham.Fb.v[i] + c * (wham.fi.m[i].F - wham.Fb.v[i]);
				// just make sure the last curve has F = 0, to avoid possible induced error on F
				for (i = 0; i < wham.fi.n; i++) wham.fi.m[i].F -= wham.fi.m[wham.fi.n - 1].F; 
			}
			iloop++;
			if (iloop > max_loops) {status = false; break;}
		}
		return status;
	}

	bool wham_read(char *parfname, WHAM &wham) {
		char buffer[256] = "\0", fname[256] = "\0", buff[256] = "\0", msg[256] = "\0";
		bool status;

		float xf = 0, xt = 0;

		ARRAY<float*> v;
		ARRAY<int> icol;
		int nCol = 4, i, NC = 0, nitem, npts = 0, iw;
		v.SetArray(nCol); icol.SetArray(nCol);
		for (i = 0; i < nCol; i++) {v.m[i] = NULL; icol.m[i] = i;}

		ifstream in;
		in.open(parfname);
		if (!in.is_open()) {
			sprintf(msg, "failure to open file %s", parfname); show_log(msg, true); return false;
		}
		if (!search(in, "[BIASED_DISTRIBUTION]", buffer)) {
			sprintf(msg, "[BIASED_DISTRIBUTION] is not defined in %s", parfname); show_log(msg, true); in.close(); return false;
		}
		nitem = 0;
		while (!in.eof()) {
			in.getline(buffer, 250);
			if (sscanf(buffer, "%s %f %f %d", fname, &xf, &xt, &NC) == 4) nitem++;
			else break;
		}
		if (nitem == 0) {
			sprintf(msg, "NO biased distribution is given in %s", parfname); show_log(msg, true); in.close(); return false;
		}
		wham.fi.SetArray(nitem);
		rewind(in);
		search(in, "[BIASED_DISTRIBUTION]", buffer);

		int npt, nf, nt;

		for(iw = 0; iw < wham.fi.n; iw++) {
			in.getline(buffer, 250);
			nitem = sscanf(buffer, "%s %f %f %d", fname, &xf, &xt, &NC);
			if (nitem != 4) continue;
			if (!read_profile(fname, 0, v, icol, npts) || npts == 0) {
				sprintf(msg, "failure to read distribution from %", fname); show_log(msg, true);
				status = false;
			}
			else status = true;
			nf = 0; nt = 0;
			for (i = 0; i < npts; i++) {
				if ((v.m[0])[i] > xf) {
					if (i == 0) nf = i;
					else if ((v.m[0])[i-1] <= xf) nf = i;
				}
				if ((v.m[0])[i] <= xt) {
					if (i == npts - 1) nt = i;
					else if ((v.m[0])[i+1] > xt) nt = i;
				}
				else break;
			}
			npt = nt - nf + 1;
			if (npt <= 1) {
				sprintf(msg, "no valid point for distribution from %s", fname); show_msg(msg, true);
				status = false;
			}

			if (status) {
				wham.fi.m[iw].set(npt);
				for (i = 0; i < npt; i++) {
					wham.fi.m[iw].x.v[i] = (v.m[0])[i + nf]; 
					wham.fi.m[iw].y.v[i] = (v.m[1])[i + nf];
					wham.fi.m[iw].s.v[i] = (v.m[2])[i + nf]; 
					wham.fi.m[iw].wb.v[i] = (v.m[3])[i + nf];
				}

				wham.fi.m[iw].Nc = NC; wham.fi.m[iw].F = 0;
			}
			for (i = 0; i < v.n; i++) {
				if (v.m[i] != NULL) delete[] v.m[i]; v.m[i] = NULL;
			}
			
			if (!status) break;
		}
		in.close();
		if (!status) {return false;}

		bool lstatus, rstatus;
		// status was true
		for (iw = 0; iw < wham.fi.n; iw++) {
			lstatus = false; rstatus = false;
			if (iw >= 1) {
				if (wham.fi.m[iw - 1].x.t[0] < wham.fi.m[iw].x.v[0]) lstatus = false;
				else lstatus = true;
			}
			else lstatus = true;

			if (iw < wham.fi.n - 1) {
				if (wham.fi.m[iw].x.t[0] < wham.fi.m[iw+1].x.v[0]) rstatus = false;
				else rstatus = true;
			}
			else rstatus = true;

			if (!lstatus) {
				sprintf(msg, "distribution %d has no jointed point to the previous one", iw); show_log(msg, true); status = false;
			}
			if (!rstatus) {
				sprintf(msg, "distribution %d has no jointed point to the next one", iw); show_log(msg, true); status = false;
			}
		}
		return status;
	}
}

extern LOG mlog;
extern "C" bool wham_pmf(char *parfname) {
	using namespace _WHAM_;
	mlog.init("wham.log");

	char msg[256] = "\0", fname[256] = "\0", title[256] = "\0";

	char ofname[256] = "\0";
	if (!read_item(parfname, "[OUTPUT_TITLE]", title)) {
		sprintf(msg, "[OUTPUT_TITLE] or the filename is not given in %s", parfname); show_log(msg, true);
		mlog.close(); return false;
	}

	WHAM wham;
	bool status = wham_read(parfname, wham);
	if (!status) {mlog.close(); return false;}
	wham.init_unbaised_prof();
	status = IterateWHAM(wham, 0.002f, 100);

	if (status) {
		sprintf(ofname, "%s.distrbt", title);
		ofstream out;
		out.open(ofname);
		if (!out.is_open()) {
			sprintf(msg, "failure to open %s", ofname); show_log(msg, true);
			mlog.close(); return false;
		}
		int i;
		for (i = 0; i < wham.f.x.N; i++) out<<wham.f.x.v[i]<<"  "<<wham.f.y.v[i].v<<"  "<<wham.f.s.v[i]<<endl;
		out.close();

		sprintf(ofname, "%s.pmf", title);
		ofstream out1;
		out1.open(ofname);
		if (!out1.is_open()) {
			sprintf(msg, "failure to open %s", ofname); show_log(msg, true);
			mlog.close(); return false;
		}
		float F0 = -log((wham.f.y.t)->v);
		for (i = 0; i < wham.f.x.N; i++) out1<<wham.f.x.v[i]<<"  "<<-log(wham.f.y.v[i].v) - F0<<"  "<<wham.f.s.v[i]<<endl;
		out1.close();

		sprintf(fname, "%s.F", title);
		ofstream out2;
		out2.open(fname);
		for (i = 0; i < wham.fi.n; i++) out2<<i<<"  "<<wham.fi.m[i].F<<endl;
		out2.close();
	}
	return status;
}

