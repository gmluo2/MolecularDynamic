// Nose-Hoover Chain, defined by G. J. Martyna, Michael L. Klein and Mark Tuckerman in J. Chem. Phys., 97, 2635 (1992)
template <int M> class NHC : public NR_NonLinearEqs<M, NHC<M> > { // M -- length of Nose-Hoover thermostats chain 
public:
	DVECTOR<M> vksi, vksi_m, ksi, ksi_ph, ksi_mh;
	float max_ksi; // maximum value of ksi
	int N; // total freedoms
	double inv_taos, dt, dts;
	// dT is normalized temperature difference dT = (Ek - 3/2 * NkT) / (3/2 * NKT); Ek is the total kinetic energy of the system with N freedoms
	double dT;
	double inv_taos_dT;

	void set_vars(float tao, float dt, float max_ksi) { // N is the total freedomes of the system
		this->inv_taos = 1. / (tao * tao);
		this->dt = dt; this->dts = dt * dt;
		this->max_ksi = max_ksi;
	};
	void set_tao(float tao) {this->inv_taos = 1. / (tao * tao);};
	void set_freedom(int NF) {this->N = NF;};
	void set_dT(double dT) {
		if (FABS(dT) > 0.8) dT =sign(dT) * 0.8;
		this->dT = dT; this->inv_taos_dT = inv_taos * dT;
	};
	void JacobianMatrix_b();
	bool verlet_NHC(float vksi_tol = 1e-8, int MaxIter = 10);
	void accept_nextstep() {
		int i;
		for (i = 0; i < M; i++) {vksi_m.v[i] = vksi.v[i]; ksi_mh.v[i] = ksi_ph.v[i];}
	};
	void reset_NHC() {
		int i;
		for (i = 0; i < M; i++) {
			vksi_m.v[i] = 0; vksi.v[i] = 0; ksi.v[i] = 0; ksi_mh.v[i] = 0; ksi_ph.v[i] = 0;
		}
	};

	NHC();/* {
		N = 0; max_ksi = 0; 
		NR_NonLinearEqs<M, NHC<M> >::m_dclass = this;
		//NR_NonLinearEqs<M, NHC<M> >::func_JacobianMatrix_b = NULL;
		NR_NonLinearEqs<M, NHC<M> >::func_JacobianMatrix_b = (void*)(&mNHC_JacobianMatrix_b);
	};*/ // moved to the end of this file
	void set_func(void *func) {NR_NonLinearEqs<M, NHC<M> >::func_JacobianMatrix_b = func;};
	void cal_equations();
};

// function defined in NR_NonLinearEqs to calculate Jacobian Matrix and b -- -F(vt)
template <int M> void NHC<M>::JacobianMatrix_b() {
// dT is normalized temperature difference dT = (Ek - 3/2 * NkT) / (3/2 * NKT); Ek is the total kinetic energy of the system with N freedoms
	int i, j, ip, im;
	int mM = M - 1;
	for (i = 0; i < M; i++) {
		for (j = 0; j < M; j++) NR_NonLinearEqs<M, NHC<M> >::J.m[i][j] = 0;
		if (i == 0) {
			NR_NonLinearEqs<M, NHC<M> >::J.m[i][0] = 4;
			NR_NonLinearEqs<M, NHC<M> >::b.v[i] = (inv_taos_dT - NR_NonLinearEqs<M, NHC<M> >::vt.v[0]) * 4;
			if (M > 1) {
				NR_NonLinearEqs<M, NHC<M> >::J.m[i][0] += NR_NonLinearEqs<M, NHC<M> >::vt.v[1] * dts + ksi_mh.v[1] * dt;
				NR_NonLinearEqs<M, NHC<M> >::J.m[i][1] += ksi_mh.v[0] * dt;
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] -= (ksi_mh.v[0] + NR_NonLinearEqs<M, NHC<M> >::vt.v[0] * dt) * (ksi_mh.v[1] + NR_NonLinearEqs<M, NHC<M> >::vt.v[1] * dt);
			}
		}
		else if (i == mM) {
			im = i - 1;
			NR_NonLinearEqs<M, NHC<M> >::J.m[i][i] = 4;
			NR_NonLinearEqs<M, NHC<M> >::b.v[i] = -(inv_taos + NR_NonLinearEqs<M, NHC<M> >::vt.v[i]) * 4;
			if (im == 0) {
				NR_NonLinearEqs<M, NHC<M> >::J.m[i][im] -= 2 * (NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dts + ksi_mh.v[im] * dt) * N;
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] += (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * N;
			}
			else {
				NR_NonLinearEqs<M, NHC<M> >::J.m[i][mM-1] -= (NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dts + ksi_mh.v[im] * dt) * 2;
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] += (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt);
			}
		}
		else {
			im = i - 1; ip = i + 1;
			NR_NonLinearEqs<M, NHC<M> >::J.m[i][i] = 4;
			NR_NonLinearEqs<M, NHC<M> >::b.v[i] = -(ksi_mh.v[i] + NR_NonLinearEqs<M, NHC<M> >::vt.v[i] * dt) * (ksi_mh.v[ip] + NR_NonLinearEqs<M, NHC<M> >::vt.v[ip] * dt) - (inv_taos + NR_NonLinearEqs<M, NHC<M> >::vt.v[i]) * 4;
			NR_NonLinearEqs<M, NHC<M> >::J.m[i][i] += NR_NonLinearEqs<M, NHC<M> >::vt.v[ip] * dts + ksi_mh.v[ip] * dt;
			NR_NonLinearEqs<M, NHC<M> >::J.m[i][ip] += NR_NonLinearEqs<M, NHC<M> >::vt.v[i] * dts + ksi_mh.v[i] * dt;
			if (im == 0) {
				NR_NonLinearEqs<M, NHC<M> >::J.m[i][im] -= 2 * (NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dts + ksi_mh.v[im] * dt) * N;
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] += (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * N;
			}
			else {
				NR_NonLinearEqs<M, NHC<M> >::J.m[i][im] -= 2 * (NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dts + ksi_mh.v[im] * dt);
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] += (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt);
			}
		}
	}
};

template <int M> void NHC<M>::cal_equations() {
// dT is normalized temperature difference dT = (Ek - 3/2 * NkT) / (3/2 * NKT); Ek is the total kinetic energy of the system with N freedoms
	int i, j, ip, im;
	int mM = M - 1;
	for (i = 0; i < M; i++) {
		if (i == 0) {
			NR_NonLinearEqs<M, NHC<M> >::b.v[i] = (inv_taos_dT - NR_NonLinearEqs<M, NHC<M> >::vt.v[0]) * 4;
			if (M > 1) {
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] -= (ksi_mh.v[0] + NR_NonLinearEqs<M, NHC<M> >::vt.v[0] * dt) * (ksi_mh.v[1] + NR_NonLinearEqs<M, NHC<M> >::vt.v[1] * dt);
			}
		}
		else if (i == mM) {
			im = i - 1;
			NR_NonLinearEqs<M, NHC<M> >::b.v[i] = -(inv_taos + NR_NonLinearEqs<M, NHC<M> >::vt.v[i]) * 4;
			if (im == 0) {
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] += (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * N;
			}
			else {
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] += (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt);
			}
		}
		else {
			im = i - 1; ip = i + 1;
			NR_NonLinearEqs<M, NHC<M> >::b.v[i] = -(ksi_mh.v[i] + NR_NonLinearEqs<M, NHC<M> >::vt.v[i] * dt) * (ksi_mh.v[ip] + NR_NonLinearEqs<M, NHC<M> >::vt.v[ip] * dt) - (inv_taos + NR_NonLinearEqs<M, NHC<M> >::vt.v[i]) * 4;
			if (im == 0) {
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] += (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * N;
			}
			else {
				NR_NonLinearEqs<M, NHC<M> >::b.v[i] += (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt) * (ksi_mh.v[im] + NR_NonLinearEqs<M, NHC<M> >::vt.v[im] * dt);
			}
		}
	}
};

template <int M> bool NHC<M>::verlet_NHC(float vksi_tol, int MaxIter) {
	int i;
	bool status = true;
	if (vksi_tol > 0.001f / dt)  vksi_tol = 0.001f / dt;
	if (M == 1) {
		vksi.v[0] = inv_taos_dT;
		ksi_ph.v[0] = ksi_mh.v[0] + vksi.v[0] * dt;
		ksi.v[0] = (ksi_mh.v[0] + ksi_ph.v[0]) * 0.5;
		if (FABS(ksi.v[0]) > max_ksi) {
			ksi.v[0] = sign(ksi.v[0]) * max_ksi;
			ksi_ph.v[0] = ksi.v[0] * 1.5 - ksi_mh.v[0] * 0.5;
		}
	}
	else {
		vksi.v[0] = inv_taos_dT;
		for (i = 1; i < M; i++) vksi.v[i] = vksi_m.v[i];
		status = NR_NonLinearEqs<M, NHC<M> >::NRIterate(vksi, vksi, 0.8, vksi_tol, 3, MaxIter);
		for (i = 0; i < M; i++) {
			ksi_ph.v[i] = ksi_mh.v[i] + vksi.v[i] * dt;
			ksi.v[i] = (ksi_mh.v[i] + ksi_ph.v[i]) * 0.5;
			if (FABS(ksi.v[i]) > max_ksi) {
				ksi.v[i] = sign(ksi.v[i]) * max_ksi;
				ksi_ph.v[i] = ksi.v[i] * 1.5 - ksi_mh.v[i] * 0.5;
			}
		}
	}

	// check whether the equations are solved correctly, monitor the array b
	//cal_equations();

	return status;
}

inline void mNHC_JacobianMatrix_b(NHC<MNHC> *mNHC) {
	return mNHC->JacobianMatrix_b();
};

template <int M> NHC<M>::NHC() {
	N = 0; max_ksi = 0; 
	NR_NonLinearEqs<M, NHC<M> >::m_dclass = this;
	//NR_NonLinearEqs<M, NHC<M> >::func_JacobianMatrix_b = NULL;
	NR_NonLinearEqs<M, NHC<M> >::func_JacobianMatrix_b = (void*)(&mNHC_JacobianMatrix_b);
};

