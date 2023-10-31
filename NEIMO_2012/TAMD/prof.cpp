#include "project.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

#include "ranlib.h"

#include "prof.h"
#include "def.h"
#include "show.h"
#include "vector.h"
#include "ZMatrix.h"

char prof_path[256] = "\0";
//DPROF erfc_prof, exp2_prof;
DPROF ErfProf, ErfcProf, Exp2Prof, ExpProf, IncGammaProf_2_3;
DPROF SinProf, CosProf;

double ErfIntegrate(double xmax, double resolution)
{
	double f = 0;
	double y1 = 0, dy = 0;
	int npts = 1;
	double step = xmax;
	double xmin = 0;
	
	double y2 = y1 / 2 + (1 + exp(-xmax * xmax)) * step / 2;
	double x = xmin;
	int nloops = 0;
	while (nloops < 10 || FABS(y2 - y1) >= resolution * y2) {
		y1 = y2;
		npts *= 2;
		step /= 2;
		y2 = y1 / 2;
		dy = 0;
		x = xmin + step;
		while (x < xmax) {
			dy += exp(-x * x);
			x += step + step;
		}
		dy *= step;
		y2 += dy;
		nloops++;
	}
	return y2;
}

double erf(double x)
{
	double y = 0;
	double resolution = 0.000001;

	double x2 = x * x, x3 = x * x2, x5 = x3 * x2, x7 = x5 * x2, x9 = x7 * x2;
	if (FABS(x) < 0.15) {
		y = x - x3 / 3 + x5 / 10 - x7 / 42 + x9 / 216;
		y *= 2 * INV_SQRT_PI;
		return y;
	}
	else if (x > 5 || x < -5) {
		y = exp(-x * x) / x * INV_SQRT_PI;
		return 1 - y;
	}

	if (x == 0) y = 0;
	else if (x > 6) return 1;
	else if (x < -6) return -1;
	else if (x > 0) y = ErfIntegrate(x, resolution);
	else y = -ErfIntegrate(-x, resolution);
	//y *= 1.1283792;  // y * 2/sqrt(PI)
	y *= 2 * INV_SQRT_PI;  // y * 2/sqrt(PI)
	return y;
}

double erfc(double x)
{
	double y = 0;
	/*
	if (x > 5 || x < -5) {
		y = exp(-x * x) / x * INV_SQRT_PI;
		return y;
	}
	else return 1 - erf(x);
	*/
	return 1 - erf(x);
}

// erfc function profile -- erfc_prof
void construct_erfc(DPROF &erfc_prof) { // erfc(x)
	double x;
	int Nw = 8;
	int npts = Nw * 500;
	erfc_prof.set(npts, 1.0 / 500);
	int i = 0;
	for (i = 0; i < erfc_prof.N; i++) {
		x = erfc_prof.x[i];
		erfc_prof.y[i] = erfc(x);
	}
	erfc_prof.y[erfc_prof.N - 1] = 0; // set the infinit value 0

	double y = 0;
	ofstream out;
	out.open("erfc.dat");
	/*
	for (x = 0; x < 10; x += 0.003) {
		ERFC(y, x, i)
		out<<x<<" "<<y<<endl;
	}
	out.close();
	*/
	for (i = 0; i < erfc_prof.N; i++) {
		out<<erfc_prof.x[i]<<" "<<erfc_prof.y[i]<<endl;
	}
	out.close();
}

void construct_erf(DPROF &erf_prof) { // erf(x)
	double x;
	int Nw = 8;
	int npts = Nw * 500;
	erf_prof.set(npts, 1.0 / 500);
	int i = 0;
	for (i = 0; i < erf_prof.N; i++) {
		x = erf_prof.x[i];
		erf_prof.y[i] = erf(x);
	}
	erf_prof.y[erf_prof.N - 1] = 1; // set the infinit value 1

	double y = 0;
	ofstream out;
	out.open("erf.dat");
	/*
	for (x = 0; x < 10; x += 0.003) {
		ERF(y, x, i)
		out<<x<<" "<<y<<endl;
	}
	out.close();
	*/

	for (i = 0; i < erf_prof.N; i++) {
		out<<erf_prof.x[i]<<" "<<erf_prof.y[i]<<endl;
	}
	out.close();
}

// exp^2 function profile -- exp2_prof
void construct_exp2(DPROF &exp2_prof) {
	double x;
	int Nw = 8;
	int npts = Nw * 500;
	exp2_prof.set(npts, 1.0 / 500);
	int i = 0;
	for (i = 0; i < exp2_prof.N; i++) {
		x = exp2_prof.x[i];
		exp2_prof.y[i] = exp(-x * x);
	}
	exp2_prof.y[exp2_prof.N - 1] = 0; // set the infinit value 0
}
/*
// erfc function profile -- erfc_prof
void construct_erfc(DPROF &erfc_prof, double ksi) { // erfc(ksi * x)
	double x;
	int Nw = 8;
	int npts = Nw * 50;
	erfc_prof.set(npts, 1.0 / ksi / 50);
	int i = 0;
	for (i = 0; i < erfc_prof.npts; i++) {
		x = ksi * erfc_prof.x[i];
		erfc_prof.y[i] = erfc(x);
	}
	erfc_prof.y[erfc_prof.npts - 1] = 0; // set the infinit value 0
}

// exp(-r * r) function profile -- exp2_prof
void construct_exp2(DPROF &exp2_prof, double ksi) {
	double x;
	int Nw = 8;
	int npts = Nw * 50;
	exp2_prof.set(npts, 1.0 / ksi / 50);
	int i = 0;
	for (i = 0; i < exp2_prof.npts; i++) {
		x = ksi * exp2_prof.x[i];
		exp2_prof.y[i] = exp(-x * x);
	}
	exp2_prof.y[exp2_prof.npts - 1] = 0; // set the infinit value 0
}
*/
// exp(-r) function profile -- exp_prof
void construct_exp(DPROF &exp_prof) {
	double x;
	int Nw = 64;
	int npts = Nw * 20;
	exp_prof.set(npts, 1.0 / 20);
	int i = 0;
	for (i = 0; i < exp_prof.N; i++) {
		x = exp_prof.x[i];
		exp_prof.y[i] = exp(-x);
	}
	exp_prof.y[exp_prof.N - 1] = 0; // set the infinit value 0
}

extern bool read_profile(char *fname, int nIgnoreLines, float **x, int ix, float **y, int iy, int& npts);
bool construct_profile_from_file(DPROF &prof, char *fname) {
	char prof_fname[256] = "\0", msg[256] = "\0";
	if (strlen(prof_path) == 0) strcpy(prof_fname, fname);
	else sprintf(prof_fname, "%s%s", ::prof_path, fname);
	prof.release();

	float *x = NULL, *y = NULL;
	int npts = 0;
	bool status = read_profile(prof_fname, 0, &x, 0, &y, 1, npts);
	if (!status) {
		sprintf(msg, "failure to read profile from %s", prof_fname); show_log(msg, true); 
		if (x != NULL) {delete[] x; x = NULL;}
		if (y != NULL) {delete[] y; y = NULL;}
		prof.release(); return false;
	}
	if (npts < 2) {
		sprintf(msg, "more than 2 points is required for a profile. check file %s", prof_fname); show_log(msg, true); 
		if (x != NULL) {delete[] x; x = NULL;}
		if (y != NULL) {delete[] y; y = NULL;}
		prof.release(); return false;
	}
	float dx = (x[npts - 1] - x[0]) / (npts - 1);
	prof.set(npts + 1, dx, x[0]);
	for (int i = 0; i < npts; i++) {prof.x[i] = x[i]; prof.y[i] = y[i];}
	prof.dx = dx; prof.x[npts] = prof.x[npts - 1] + dx; prof.y[npts] = prof.y[npts - 1];
	if (x != NULL) {delete[] x; x = NULL;}
	if (y != NULL) {delete[] y; y = NULL;}
	return true;
}

void construct_sin_cos(DPROF &SinProf, DPROF &CosProf) {
	int npts = 360, i;
	double x;
	SinProf.set(npts + 1, PI2 / npts);
	CosProf.set(npts + 1, PI2 / npts);
	for (i = 0; i < npts + 1; i++) {
		x = SinProf.x[i];
		SinProf.y[i] = sin(x);
		CosProf.y[i] = cos(x);
	}
}

void construct_profiles() {
	if (ErfProf.x == NULL) {
		if (!construct_profile_from_file(ErfProf, "erf.dat")) {
			construct_erf(ErfProf);
			show_log("calculate error function", true);
		}
		/*
		{
			ofstream out;
			out.open("erf_test.dat");
			int i, ti;
			float x = 0, y, dx = 0.01;
			for (i = 0; i < 1000; i++) {
				x = i * dx;
				ERF(y, x, ti)
				out<<x<<"  "<<y<<endl;
			}
			out.close();
		}
		*/
	}

	if (ErfcProf.x == NULL) {
		if (!construct_profile_from_file(ErfcProf, "erfc.dat")) {
			construct_erfc(ErfcProf);
			show_log("calculate complementary error function", true);
		}
	}

	if (ExpProf.x == NULL) construct_exp(ExpProf);
	if (Exp2Prof.x == NULL) construct_exp2(Exp2Prof);
	if (SinProf.x == NULL || CosProf.x == NULL) construct_sin_cos(SinProf, CosProf);

	if (IncGammaProf_2_3.x == NULL) {
		if (!construct_profile_from_file(IncGammaProf_2_3, "gamma_Inc_2_3.dat")) {
			show_log("failure to construct profile for Incomplete Gamma Function Gamma(a, x) with a = 2/3", true);
			show_log("Please check if this function is required in simulation.", true);
			show_log("So far, this function is required for Thole polarization with exponential screening model exp(-a * u^3)", true);
			show_log("", true);
		}
	}
}

