#include "project.h"

#include <iostream>
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

#include "ZMatrix.h"

#include "EwaldSumVar.h"

#include "MM.h"
#include "Mol.h"
#include "CMM.h"
#include "cluster.h"
#include "Interaction.h"
#include "MD.h"
#include "MD_SM.h"

#include "var.h"

double max_diff(VECTOR3 &v1, VECTOR3 &v2) {
	int nc = 0, i = 0;
	double diff = 0;
	VECTOR<3> d;
	V3minusV3(v1, v2, d)
	MAX(d, i, diff)
	diff = FABS(diff);
	return diff;
}
//#ifndef _DISABLE_
bool SM_LeapFrog_velocity_acceleration(bool bHooverDyn, SMOLECULE &m, LFVERLET_SM_MD& mdpar, int md_loop, int &nconvg_loops, char *msg_show, char *msg_log, int n_msglen) {
	char msg[256] = "\0";
	int i = 0;
	double T = 0, Ek = 0, Ek0 = 0, dEk = 0;
	double Ethermal = 3; // translation + rotation ; in unit kT

	VECTOR3 v_new, v_old, w_new, w_old;
	double vmax_diff = 0, wmax_diff = 0;
	double vConvg = 0.001 / mdpar.dt, wConvg = 5e-5 / mdpar.dt;
	
	bool convergent = false;
	int nmax_convg = 20;

	VECTOR<6> mv0, malpha0;

	double ksi_Hoover = mdpar.hdyn->nhc->ksi.v[0];

	double ds_max = 0.5f, domega_max = 0.04, dmax = 0, dguess = 0;

	int nconvg = 0;
	while (!convergent) {
		V2wv(w_old, v_old, m.c->bkm->V0)

		Ek0 = calKineticEnergy(m); // kinetic energy 
		if (bHooverDyn) {
			T = Ek0 / Ethermal;
			mdpar.hdyn->nhc->set_dT(T - 1);
			mdpar.hdyn->nhc->verlet_NHC();
			ksi_Hoover = mdpar.hdyn->nhc->ksi.v[0];
		}
		else ksi_Hoover = 0;

		if (bHooverDyn && max_ksi_Hoover > 1e-8 && T > 4 /*ksi_Hoover / max_ksi_Hoover >= 0.98*/) {
			// decrease the speed is not efficient enough
			// so scale the kinetic energy back to thermal energy
			// and calculate the a, b since speed changed
			InitVelocity(m);
			Speed_Kinetic(mdpar, m);
			ksi_Hoover = 0; mdpar.hdyn->nhc->reset_NHC();
			md_loop = 0; // to trigger the case that do not get consistent velocity and accleration
			continue;
		}

		calAccel(m, (float)ksi_Hoover);

		if (nconvg == 0) {// backup tpp
			for (i = 0; i < 6; i++) {
				mv0.v[i] = m.c->bkm->V0.v[i]; malpha0.v[i] = m.c->bkm->alpha.v[i];
			}
		}

		if (md_loop == 0) Speed_Kinetic(mdpar, m);
		else {
			// get new speed at current position and +1/2 time step
			Speed_Verlet(mdpar, m); // reset the speed in mm to the new speed
		}

		V2wv(w_new, v_new, m.c->bkm->V0);

		vmax_diff = max_diff(v_new, v_old);
		wmax_diff = max_diff(w_new, w_old);
		Ek = calKineticEnergy(m); // kinetic energy 

		if (vmax_diff < vConvg && wmax_diff < wConvg) convergent = true;
		else {
			dEk = Ek - Ethermal; dEk /= Ethermal;
			if (dEk > 1) {
				sprintf(msg, "Kinetic energy jumps from %f to %f kT with iteration %d at loop %d", Ek0, Ek, nconvg, md_loop);
				if (n_msglen - strlen(msg_show) > 150) {strcat(msg_show, "\n"); strcat(msg_show, msg);}
				//show_infor(msg); //mlog.show(msg);
				if (nconvg > 0) {
					// using the first iteration acceleration, see below
					nconvg = nmax_convg + 1;
				}
				else {
					sprintf(msg, "ERROR : Unable to handle this problem! ");
					if (n_msglen - strlen(msg_show) > 100) {strcat(msg_show, "\n"); strcat(msg_show, msg);}
					//mlog.show(msg); //show_infor(msg); 
					nconvg = nmax_convg + 1; 
				}
				break;
			}
		}

		if (md_loop == 0) {convergent = true; break;}

		nconvg++;
		if (nconvg > nmax_convg) {convergent = true; break;}
	}

	if (!convergent) {
		sprintf(msg, "failure to get consistent speed & acceleration at loop %d \n", md_loop);
		if (n_msglen - strlen(msg_show) > 100) {strcat(msg_show, "\n"); strcat(msg_show, msg);}
		//show_infor(msg); strcpy(msg, "\0");
		// use the acceleration speed based on initial speed
		for (i = 0; i < 6; i++) {
			m.c->bkm->V0.v[i] = mv0.v[i];
			m.c->bkm->alpha.v[i] = malpha0.v[i];
			dmax = sign(malpha0.v[i]) * (i < 3 ? domega_max : ds_max);
			dguess = mv0.v[i] * mdpar.dt + 0.5 * malpha0.v[i] * mdpar.dt * mdpar.dt;
			if ((dmax > 0 && dguess > dmax) || (dmax < 0 && dguess < dmax)) {
				m.c->bkm->alpha.v[i] = 2 * (dmax - mv0.v[i] * mdpar.dt) / (mdpar.dt * mdpar.dt);
			}
		}

		Speed_Verlet(mdpar, m);
		Ek = ScaleKineticEnergy(m, Ethermal);
		Speed_Kinetic(mdpar, m);
		ksi_Hoover = 0; mdpar.hdyn->nhc->reset_NHC();
		convergent = false;
	}
	m.Ek = (float)Ek;

	if (bHooverDyn) mdpar.hdyn->nhc->accept_nextstep(); 

	return convergent;
}

bool SM_SpeedVerlet(bool local_relax, bool speed_verlet, SMOLECULE &m, LFVERLET_SM_MD& mdpar, double& vmax_diff, double& wmax_diff, char *msg_show, char *msg_log, int n_msglen) {
	char msg[256] = "\0";
	int i = 0;
	double T = 0;
	double Ethermal = 3; // translation + rotation ; in unit kT
	float Ek0 = m.Ek;

	VECTOR3 v_new, v_old, w_new, w_old;
	
	double ksi_Hoover = mdpar.hdyn->nhc->ksi.v[0];
	bool speed_scaled = false;

	if (local_relax) {
		V6zero(m.c->bkm->V0)  V6zero(m.c->bkm->alpha) 
	}

	V2wv(w_old, v_old, m.c->bkm->V0)

	//Ek0 = calKineticEnergy(m); // kinetic energy 
	T = m.Ek / Ethermal;

	if (T > 3) {
		// decrease the speed is not efficient enough
		// so scale the kinetic energy back to thermal energy
		// and calculate the a, b since speed changed
		speed_scaled = true;
		InitVelocity(m);
		//m.Ek = (float)ScaleKineticEnergy(m, Ethermal);
		ksi_Hoover = 0;
	}

	calAccel(m, (float)ksi_Hoover);

	if (local_relax || speed_scaled || !speed_verlet) Speed_Kinetic(mdpar, m);
	else {
		// get new speed at current position and +1/2 time step
		Speed_Verlet(mdpar, m); // reset the speed in mm to the new speed
	}

	V2wv(w_new, v_new, m.c->bkm->V0);

	vmax_diff = max_diff(v_new, v_old);
	wmax_diff = max_diff(w_new, w_old);
	calKineticEnergy(m); // kinetic energy 

	T = (m.Ek - Ek0) / Ethermal;
	if (T > 3) {
		sprintf(msg, "Kinetic energy jumps from %f to %f kT", Ek0, m.Ek);
		if (n_msglen - strlen(msg_show) > 150) {strcat(msg_show, "\n"); strcat(msg_show, msg);}
		//show_infor(msg); //mlog.show(msg);
		m.Ek = (float)ScaleKineticEnergy(m, Ethermal);
		Speed_Kinetic(mdpar, m);
		calKineticEnergy(m);
		return false;
	}

	bool status = true;
	if (local_relax || speed_scaled) status = false;
	return true;
}
//#endif // _DISABLE_

bool SM_SpeedVerlet_Simplified(bool speed_verlet, SMOLECULE &m, LFVERLET_SM_MD& mdpar, double& vmax_diff, double& wmax_diff) {
	VECTOR3 v_new, v_old, w_new, w_old;

	V2wv(w_old, v_old, m.c->bkm->V0)

	m.calTotalForce();
	calAccel_f(m);

	if (speed_verlet) { // get new speed at current position and +1/2 time step
		Speed_Verlet(mdpar, m); // reset the speed in mm to the new speed
	}
	else Speed_Kinetic(mdpar, m);

	V2wv(w_new, v_new, m.c->bkm->V0);

	vmax_diff = max_diff(v_new, v_old);
	wmax_diff = max_diff(w_new, w_old);
	calKineticEnergy(m); // kinetic energy 

	return true;
}

// FUNCTIONS FOR PMOLECULES
//#ifndef _DISABLE_
bool PM_LeapFrog_velocity_acceleration(bool bHooverDyn, PMOLECULE &m, LFVERLET_PM_MD& mdpar, int md_loop, int &nconvg_loops, char *msg_show, char *msg_log, int n_msglen) {
	char msg[256] = "\0";
	int i = 0;
	double T = 0, Ek = 0, Ek0 = 0, dEk = 0;
	double Ethermal = 1.5; // translation + rotation ; in unit kT

	VECTOR3 v_new, v_old;
	double vmax_diff = 0;
	double vConvg = 0.001 / mdpar.dt;
	double ds_max = 0.2, dguess, dmax;
	
	bool convergent = false;
	int nmax_convg = 20;

	VECTOR3 mv0, malpha0;

	double ksi_Hoover = mdpar.hdyn->nhc->ksi.v[0];

	int nconvg = 0;
	while (!convergent) {
		V32V3(m.v, v_old)

		Ek0 = calKineticEnergy(m); // kinetic energy 
		if (bHooverDyn) {
			T = Ek0 / Ethermal;
			mdpar.hdyn->nhc->set_dT(T - 1);
			mdpar.hdyn->nhc->verlet_NHC();
			ksi_Hoover = mdpar.hdyn->nhc->ksi.v[0];
		}
		else ksi_Hoover = 0;

		if (bHooverDyn && max_ksi_Hoover > 1e-8 && T > 4 /*ksi_Hoover / max_ksi_Hoover >= 0.98*/) {
			// decrease the speed is not efficient enough
			// so scale the kinetic energy back to thermal energy
			// and calculate the a, b since speed changed
			ScaleKineticEnergy(m, Ethermal);
			Speed_Kinetic(mdpar, m);
			ksi_Hoover = 0; mdpar.hdyn->nhc->reset_NHC();
			md_loop = 0; // to trigger the case that do not get consistent velocity and accleration
			continue;
		}

		calAccel(m, (float)ksi_Hoover);

		if (nconvg == 0) {// backup tpp
			V32V3(m.v, mv0) V32V3(m.alpha, malpha0)
		}

		if (md_loop == 0) Speed_Kinetic(mdpar, m);
		else {
			// get new speed at current position and +1/2 time step
			Speed_Verlet(mdpar, m, true); // reset the speed in mm to the new speed
		}

		V3minusV3(m.v, v_old, v_new)
		vmax_diff = v_new.abs();
		Ek = calKineticEnergy(m); // kinetic energy 

		if (vmax_diff < vConvg) convergent = true;
		else {
			dEk = Ek - Ethermal; dEk /= Ethermal;
			if (dEk > 1) {
				sprintf(msg, "Kinetic energy jumps from %f to %f kT with iteration %d at loop %d", Ek0, Ek, nconvg, md_loop);
				if (n_msglen - strlen(msg_show) > 150) {strcat(msg_show, "\n"); strcat(msg_show, msg);}
				//show_infor(msg); //mlog.show(msg);
				if (nconvg > 0) {
					// using the first iteration acceleration, see below
					nconvg = nmax_convg + 1;
				}
				else {
					sprintf(msg, "ERROR : Unable to handle this problem! ");
					if (n_msglen - strlen(msg_show) > 100) {strcat(msg_show, "\n"); strcat(msg_show, msg);}
					//mlog.show(msg); //show_infor(msg); 
					nconvg = nmax_convg + 1; 
				}
				break;
			}
		}

		if (md_loop == 0) {convergent = true; break;}

		nconvg++;
		if (nconvg > nmax_convg) {convergent = true; break;}
	}

	if (!convergent) {
		sprintf(msg, "failure to get consistent speed & acceleration at loop %d \n", md_loop);
		if (n_msglen - strlen(msg_show) > 100) {strcat(msg_show, "\n"); strcat(msg_show, msg);}
		//show_infor(msg); strcpy(msg, "\0");
		// use the acceleration speed based on initial speed
		for (i = 0; i < 3; i++) {
			m.v.v[i] = mv0.v[i];
			m.alpha.v[i] = malpha0.v[i];
			dmax = sign(malpha0.v[i]) * ds_max;
			dguess = mv0.v[i] * mdpar.dt + 0.5 * malpha0.v[i] * mdpar.dt * mdpar.dt;
			if ((dmax > 0 && dguess > dmax) || (dmax < 0 && dguess < dmax)) {
				m.alpha.v[i] = 2 * (dmax - mv0.v[i] * mdpar.dt) / (mdpar.dt * mdpar.dt);
			}
		}

		Speed_Verlet(mdpar, m, true);
		Ek = ScaleKineticEnergy(m, Ethermal);
		Speed_Kinetic(mdpar, m);
		ksi_Hoover = 0; 
		if (bHooverDyn) mdpar.hdyn->nhc->reset_NHC();
	}
	m.Ek = (float)Ek;

	if (bHooverDyn) mdpar.hdyn->nhc->accept_nextstep();

	return convergent;
}


bool PM_SpeedVerlet(bool local_relax, bool speed_verlet, PMOLECULE &m, LFVERLET_PM_MD& mdpar, double& vmax_diff, char *msg_show, char *msg_log, int n_msglen) {
	char msg[256] = "\0";
	double T = 0;
	double Ethermal = 1.5; // translation; in unit kT
	float Ek0 = m.Ek;

	VECTOR3 v_new, v_old;
	
	double ksi_Hoover = mdpar.hdyn->nhc->ksi.v[0];
	bool speed_scaled = false;

	if (local_relax) {
		V6zero(m.v)  V6zero(m.alpha) 
	}
	V32V3(m.v, v_old)

	T = m.Ek / Ethermal;

	if (T > 3) {
		// decrease the speed is not efficient enough
		// so scale the kinetic energy back to thermal energy
		// and calculate the a, b since speed changed
		speed_scaled = true;
		m.Ek = (float)ScaleKineticEnergy(m, Ethermal);
		ksi_Hoover = 0;
	}

	calAccel(m, (float)ksi_Hoover);

	if (local_relax || speed_scaled || !speed_verlet) Speed_Kinetic(mdpar, m);
	else {
		// get new speed at current position and +1/2 time step
		Speed_Verlet(mdpar, m, true); // reset the speed in mm to the new speed
	}

	V32V3(m.v, v_new)
	vmax_diff = max_diff(v_new, v_old);
	calKineticEnergy(m); // kinetic energy 

	T = (m.Ek - Ek0) / Ethermal;
	if (T > 3) {
		sprintf(msg, "Kinetic energy jumps from %f to %f kT", Ek0, m.Ek);
		if (n_msglen - strlen(msg_show) > 150) {strcat(msg_show, "\n"); strcat(msg_show, msg);}
		//show_infor(msg); //mlog.show(msg);
		m.Ek = (float)ScaleKineticEnergy(m, Ethermal);
		Speed_Kinetic(mdpar, m);
		calKineticEnergy(m);
		return false;
	}

	bool status = true;
	if (local_relax || speed_scaled) status = false;
	return true;
}
//#endif // _DISABLE_

bool PM_SpeedVerlet_Simplified(bool speed_verlet, PMOLECULE &m, LFVERLET_PM_MD& mdpar, double& vmax_diff) {
	VECTOR3 v_new, v_old;
	
	V32V3(m.v, v_old)

	m.calTotalForce();
	calAccel_f(m);

	if (speed_verlet) {// get new speed at current position and +1/2 time step
		Speed_Verlet(mdpar, m, true); // reset the speed in mm to the new speed
	}
	Speed_Kinetic(mdpar, m);

	V32V3(m.v, v_new)
	vmax_diff = max_diff(v_new, v_old);
	calKineticEnergy(m); // kinetic energy 

	return true;
}
