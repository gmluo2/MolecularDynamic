#if _SYS_ == _LINUX_SYS_
template <class T> void release_chain_1(CHAIN<T> **ch, bool release_context = false) {
	::release_chain<T>(ch, release_context);
}
#endif

extern bool local_relax;
extern float rSolvent, dSolv_max;

template <class atom> class CMM_CELL : public FLEX_BUFF<atom> {
public:
	VECTOR<3> r0; // center position of the cell
#if _EWALD_SUM == _MultiPole_EWALD_SUM
	// electric multipole
	ELECTRO_MULTIPOLE *emp;
#endif
	bool eNeutral;

	CMM_CELL() {
		eNeutral = true;
	#if _EWALD_SUM == _MultiPole_EWALD_SUM
		emp = NULL;
	#if _IGNORE_ELECTROSTATIC_INTERACTION_ == 0
		emp = new ELECTRO_MULTIPOLE;
	#endif
	#endif  //_EWALD_SUM
	};
	void release_chain() {
		FLEX_BUFF<atom>::release();
	};
	void setup_emp() {
	#if _EWALD_SUM == _MultiPole_EWALD_SUM
		if (emp == NULL) emp = new ELECTRO_MULTIPOLE;
	#endif
	};
	~CMM_CELL() {
		release_chain();
	#if _EWALD_SUM == _MultiPole_EWALD_SUM
		if (emp != NULL) {delete emp; emp = NULL;}
	#endif
	};
	void calElectroMultipole(bool bEachCluster, bool bEwaldSum) {
	#if _EWALD_SUM == _MultiPole_EWALD_SUM
	#error this function needs to be finished to calculate the multipole with cluster-array, ca
		int i = 0, j = 0, nc[3];
		CHAIN<atom>* pch = NULL;
		atom *pc = NULL;

		VECTOR<3> dr;

		ELECTRO_MULTIPOLE *mp = this->emp, *cmp = NULL;

		mp->q = 0;
		V3zero(mp->mu)
		Mzero(mp->Q, i, j)
		//CMzero(mp->O, i, j, k)

		for (i = 0; i < 3; i++) mp->r0.v[i] = this->r0.v[i];
	
		pch = this->ch;
		while (pch != NULL) {
			pc = pch->p;
			if (pc->eNeutral) {pch = pch->next; continue;}

			if (bEachCluster) pc->calElectroMultipole(bEwaldSum); // calculate multipole of the cluster
			cmp = pc->emp;
			mp->q += pc->emp->q;
			for (i = 0; i < 3; i++) {
				dr.v[i] = cmp->r0.v[i] - mp->r0.v[i];
				mp->mu.v[i] += cmp->q * dr.v[i] + cmp->mu.v[i];
			}
			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
				// we kicked off the term in Q, - 1/2 * r^2 * delta_i,j, since it is useless in electric field and force calculation
				mp->Q.m[i][j] += cmp->Q.m[i][j] + 1.5 * (cmp->mu.v[i] * dr.v[j] 
					+ cmp->mu.v[j] * dr.v[i] + cmp->q * dr.v[i] * dr.v[j]);
			}
			pch = pch->next;
		}

		this->eNeutral = mp->neutral();
		if (this->eNeutral) return; // no charge or multipole 

//extern double kappa, kappa3;
//extern float r_Ewd, rcut_Ewd, r2cut_Ewd, rcut_LJ, r2cut_LJ;
extern VECTOR3 h[NH][NH][NH];
//extern double H2[NH][NH][NH];
//extern double Ah[NH][NH][NH];
//extern float AhMin;

		// the terms for Ewald Sum
		double hR = 0, t1 = 0, t2 = 0;
		if (bEwaldSum) {
			for (nc[0] = 0; nc[0] < NH; nc[0]++) {for (nc[1] = 0; nc[1] < NH; nc[1]++) {for (nc[2] = 0; nc[2] < NH; nc[2]++) {
				scale_uv3(mp->r0, h[nc[0]][nc[1]][nc[2]], hR)
				PERIOD(hR, 0, PI2, t1) COS(mp->fcos.m[nc[0]][nc[1]][nc[2]], t1, i, t1) SIN(mp->fsin.m[nc[0]][nc[1]][nc[2]], t1, i, t2)
				scale_uv3(mp->mu, h[nc[0]][nc[1]][nc[2]], mp->muh.m[nc[0]][nc[1]][nc[2]])
				scale_VMV(h[nc[0]][nc[1]][nc[2]], mp->Q, h[nc[0]][nc[1]][nc[2]], mp->Qh.m[nc[0]][nc[1]][nc[2]], i, j, 3)
			}}}
		}
		#endif
	};
};

#define eCell(c, i, j, k) (*(c.cell + k + j * c.nz + i * c.nyz))
#define cell_n(c, i, j, k, indx) indx = k + j * c.nz + i * c.nyz;
#define cell_indx(c, nc, i, j, k) i = nc / c.nyz; j = nc - i * c.nyz; j /= c.nz; k = nc - i * c.nyz - j * c.nz;

template <class atom> class CMM_CELL3D {
public:
	bool bChain;
	FLEX_ASYM_CMATRIX< CMM_CELL<atom> > bcell; // the deepest (basic) cell structure 
	int nx, ny, nz; // number of cell in one dimension
	int nyz; // n * n
	ARRAY< CMM_CELL<atom>* > acell;

	float hw[3]; // half width of the CMM_CELL in all dimensions
	float fw[3]; // width of the CMM_CELL in all dimensions
	float fw_max; // maximum value of the width of CMM_CELL in 3-d
	bool periodic_cell;
	float xl[3], xr[3], xd[3]; // cell range of each axis [xl, xr]; xd = xr - xl

	int nCheckCluster;

	_EwaldSum_real_::SURF_DIPOLE mu_surf;

	CMM_CELL3D() {
		bChain = true;
		nx = 0; ny = 0; nz = 0; nyz = 0; 
		nCheckCluster = 3;
		int i = 0;
		for (i = 0; i < 3; i++) {hw[i] = 0; xl[i] = 0; xr[i] = 0; xd[i] = 0;}
		periodic_cell = false; mu_surf.reset();
	};
	void set_cell(int nx, int ny, int nz) {
		int i = 0, ix = 0, iy = 0, iz = 0;
		if (acell.m != NULL) {for (i = 0; i < acell.n; i++) acell.m[i] = NULL; acell.release();}
		bcell.release();
		this->nx = 0; this->ny = 0; this->nz = 0; this->nyz = 0;
		if (nx <= 0 || ny <= 0 || nz <= 0) return;
		this->nx = nx; this->ny = ny; this->nz = nz; this->nyz = ny * nz; 
		bcell.set(nx, ny, nz);
		acell.SetArray(nx * ny * nz);
		i = 0;
		for (ix = 0; ix < bcell.nx; ix++) {for (iy = 0; iy < bcell.ny; iy++) {for (iz = 0; iz < bcell.nz; iz++) {
			acell.m[i] = &(bcell.m[ix].m[iy].m[iz]); i++;
		}}}
	};
	void set_cell_range(bool periodic_cell, float xl, float xr, float yl, float yr, float zl, float zr) {
		this->periodic_cell = periodic_cell;
		this->xl[0] = xl; this->xr[0] = xr; this->xd[0] = xr - xl;
		this->xl[1] = yl; this->xr[1] = yr; this->xd[1] = yr - yl;
		this->xl[2] = zl; this->xr[2] = zr; this->xd[2] = zr - zl;
		this->hw[0] = this->xd[0] / (this->nx + this->nx);
		this->hw[1] = this->xd[1] / (this->ny + this->ny);
		this->hw[2] = this->xd[2] / (this->nz + this->nz);
		this->fw[0] = this->hw[0] + this->hw[0];
		this->fw[1] = this->hw[1] + this->hw[1];
		this->fw[2] = this->hw[2] + this->hw[2];
		fw_max = this->fw[0];
		fw_max = (fw_max < this->fw[1] ? this->fw[1] : fw_max);
		fw_max = (fw_max < this->fw[2] ? this->fw[2] : fw_max);
	};
	void set_cell_pos() {
		int i, j, k;
		CMM_CELL<atom> *pcell = NULL;
		float xc[3] = {0, 0, 0}; // corner position
		for (i = 0; i < 3; i++) xc[i] = this->xl[i] + this->hw[i];
		for (i = 0; i < nx; i++) {for (j = 0; j < ny; j++) {for (k = 0; k < nz; k++) {
			pcell = &(bcell.m[i].m[j].m[k]);
			pcell->r0.v[0] = i * this->fw[0] + xc[0];
			pcell->r0.v[1] = j * this->fw[1] + xc[1];
			pcell->r0.v[2] = k * this->fw[2] + xc[2];
		}}}
	};
	void reset_cell_chain() {
		int nc = 0; for (nc = 0; nc < acell.n; nc++) acell.m[nc]->release_chain();
	};
	void reset_basic_cell_matrix() {
		int i, j, k;
		for (i = 0; i < acell.n; i++) acell.m[i] = NULL; acell.release();
		for (i = 0; i < nx; i++) {for (j = 0; j < ny; j++) {for (k = 0; k < nz; k++) {
			bcell.m[i].m[j].m[k].release_chain();
		}}}
		bcell.release();
		this->nx = 0; this->ny = 0; this->nz = 0; this->nyz = 0;
	};
	~CMM_CELL3D() {
		reset_basic_cell_matrix();
	};
	int number_atoms() {
		int ic, nc = 0;
		for (ic = 0; ic < acell.n; ic++) nc += acell.m[ic]->length();
		return nc;
	};
};

#define CMMCELL_r2INDX(r, c, i, j, k, dr) dr.v[0] = r.v[0] - c.xl[0]; dr.v[1] = r.v[1] - c.xl[1]; dr.v[2] = r.v[2] - c.xl[2]; \
	i = int(dr.v[0] / c.fw[0]); j = int(dr.v[1] / c.fw[1]); k = int(dr.v[2] / c.fw[2]); \
	if (i < 0) i = 0; if (j < 0) j = 0; if (k < 0) k = 0; \
	if (i >= c.nx) i = c.nx - 1; if (j >= c.ny) j = c.ny - 1; if (k >= c.nz) k = c.nz - 1;

#define InCell(r, r0, dr, hw, i, status) status = true; \
	for (i = 0; i < dr.n; i++) {\
		dr.v[i] = r.v[i] - r0.v[i]; dr.v[i] = FABS(dr.v[i]); \
		if (dr.v[i] >= hw[i]) {status = false; break;}\
	}

#define PeriodCellCloestDist2(dr, cell, R2, u, nx, ny, nz) u.v[0] = FABS(dr.v[0]); u.v[1] = FABS(dr.v[1]); u.v[2] = FABS(dr.v[2]); \
	if (u.v[0] > cell.xd[0] * 0.5) {u.v[0] = cell.xd[0] - u.v[0]; nx = (dr.v[0] > 0 ? -1 : 1);} else nx = 0;\
	if (u.v[1] > cell.xd[1] * 0.5) {u.v[1] = cell.xd[1] - u.v[1]; ny = (dr.v[1] > 0 ? -1 : 1);} else ny = 0;\
	if (u.v[2] > cell.xd[2] * 0.5) {u.v[2] = cell.xd[2] - u.v[2]; nz = (dr.v[2] > 0 ? -1 : 1);} else nz = 0;\
	R2 = u.v[0] * u.v[0] + u.v[1] * u.v[1] + u.v[2] * u.v[2];

template <class atom> void cmm_check_cluster(CMM_CELL3D<atom> &cell) {
	int nc = 0, n_new = 0;
	int ni = 0, nj = 0, nk = 0;
	VECTOR<3> dr, rc;
	int vi = 0;
	bool in_cell;

	CMM_CELL<atom> *pcell = NULL;
	atom* pc = NULL;
	ITERATOR<atom, CMM_CELL<atom> > iterator;
	for (nc = 0; nc < cell.acell.n; nc++) { 
		pcell = cell.acell.m[nc]; iterator.set(pcell); iterator.start_iterate();
		while (!iterator.eoi()) {
			pc = iterator.current();
			V32V3((*(pc->vp)), rc) // use vp position
			InCell(rc, pcell->r0, dr, cell.hw, ni, in_cell)
			if (!in_cell) {
				CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
				cell_n(cell, ni, nj, nk, n_new)
				if (n_new >= cell.acell.n) n_new = cell.acell.n - 1;
				if (n_new != nc) {
					// current cluster will be moved away to new cell
					// so plist go to next cluster first
					cell.acell.m[nc]->detach((atom*)pc);
					cell.acell.m[n_new]->attach((atom*)pc);
					iterator.next_iterate();
					continue;
				}
			}
			iterator.next_iterate();
		}
	}
	//for (nc = 0; nc < cell.acell.n; nc++) subcell_check_cluster(cell.acell.m[nc]);
};

template <class atom> void cmm_init_subcell(CMM_CELL3D<atom> &cell) {
	// subcell is constructed for Multipole method for Ewald Sum
#if _EWALD_SUM == _MultiPole_EWALD_SUM
	int nc = 0;
	for (nc = 0; nc < cell.acell.n; nc++) {
		if (cell.acell.m[nc].emp != NULL) cell.acell.m[nc].emp->d = cell.hw[0];
	}
#endif
};
/*
template <class atom> void subcell_check_cluster(CMM_CELL<atom> &cell) {
	if (cell.subCell == NULL) return;
	int ni = 0, nj = 0, nk = 0;
	int nc[3];
	VECTOR3 rc, dr;
	CHAIN<atom> *ch = NULL;
	CMM_CELL<atom> *c = NULL;
	atom *pc = NULL;

	int N = cell.subCell->n;

	for (ni = 0; ni < N; ni++) {for (nj = 0; nj < N; nj++) {for (nk = 0; nk < N; nk++) {
		c = &(cell.subCell->cell[ni][nj][nk]);
		c->release_chain();
	}}}

	ch = cell.ch;
	while (ch != NULL) {
		pc = ch->p;
		//for (i = 0; i < 3; i++) rc.v[i] = pc->Op->v[i];
		//V3plusV3((*(pc->vp)), pc->dr_cmm, rc)
		V32V3((*(pc->vp)), rc)
		SUBCELL_r2INDX(rc, (*(cell.subCell)), nc[0], nc[1], nc[2], dr)
		if (nc[0] > N) nc[0] = N; if (nc[0] < 0) nc[0] = 0;
		if (nc[1] > N) nc[1] = N; if (nc[1] < 0) nc[1] = 0;
		if (nc[2] > N) nc[2] = N; if (nc[2] < 0) nc[2] = 0;
		cell.subCell->cell[nc[0]][nc[1]][nc[2]].attach(pc);

		ch = ch->next;
	}
	for (ni = 0; ni < N; ni++) {for (nj = 0; nj < N; nj++) {for (nk = 0; nk < N; nk++) {
		subcell_check_cluster(cell.subCell->cell[ni][nj][nk]);
	}}}
};
*/

#if _EWALD_SUM == _MultiPole_EWALD_SUM
template <class atom> class CMM_MULTIPOLES_THREAD_VARS : public _THREAD_VARS {
public:
	CMM_CELL3D<atom> *pcell;
	bool bEachCluster;
	bool bEwaldSum; // calculating the parameters for EwaldSum ?
	int n1Cell, n2Cell;

	CMM_MULTIPOLES_THREAD_VARS() {
		pcell = NULL; bEachCluster = false; bEwaldSum = false; n1Cell = 0; n2Cell = 0;
	};
	void set_pars(int iThread, CMM_CELL3D<atom>* pcell, bool bEachCluster, bool bEwaldSum, int n1Cell, int n2Cell) {
		this->pcell = pcell; this->n1Cell = n1Cell; this->n2Cell = n2Cell;
		this->bEachCluster = bEachCluster; this->bEwaldSum = bEwaldSum;
		_THREAD_VARS::set_pars(iThread);
	};
	~CMM_MULTIPOLES_THREAD_VARS() {};
};

//void calColumbMultipoles(CMM_CELL3D<BASIC_CLUSTER> &cell, bool bEachCluster, bool bEwaldSum, int n1Cell, int n2Cell) {
template <class atom> void calColumbMultipoles(CMM_CELL3D<atom> &cell, bool bEachCluster, bool bEwaldSum, int n1Cell, int n2Cell) {
	int ncell;
	for (ncell = n1Cell; ncell <= n2Cell; ncell++) {
		if (ncell >= cell.acell.n) break;
		cell.acell.m[ncell]->calElectroMultipole(bEachCluster, bEwaldSum);
	}
}

// suppose the multipoles of all the cells are calculated
template <class atom> void calCellElectroMultipole(CMM_CELL3D<atom> &cell) {
	int nCell, i;
	VECTOR3 dr;
	V3zero(cell.mu)
	for (nCell = 0; nCell < cell.acell.n; nCell++) {
		for (i = 0; i < 3; i++) {
			cell.mu.v[i] += cell.acell.m[nCell]->emp->mu.v[i];
			if (cell.acell.m[nCell]->emp->q != 0) {
				dr.v[i] = cell.acell.m[nCell]->emp->r0.v[i];
				cell.mu.v[i] += cell.acell.m[nCell]->emp->q * dr.v[i];
			}
		}
	}
}

#if _SYS_ == _WINDOWS_SYS_ 
template <class atom> int calColumbMultipolesThread(LPVOID *var) {
#elif  _SYS_ == _LINUX_SYS_
template <class atom> void* calColumbMultipolesThread(void *var) {
#endif
	CMM_MULTIPOLES_THREAD_VARS<atom> *cmm_var = (CMM_MULTIPOLES_THREAD_VARS<atom>*)var;
	calColumbMultipoles<atom>(*(cmm_var->pcell), cmm_var->bEachCluster, cmm_var->bEwaldSum, cmm_var->n1Cell, cmm_var->n2Cell);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(cmm_var->hMMThread);
	AfxEndThread(1, TRUE);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

//void calColumbMultipoles(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, bool bEachCluster, bool bEwaldSum, int nThreads);
template <class atom> void calColumbMultipoles(CMM_CELL3D<atom> &cell, int n1Cell, int n2Cell, bool bEachCluster, bool bEwaldSum, int nThreads) {
	if (n2Cell >= cell.nCells) n2Cell = cell.nCells - 1;
	int nc = 0;
	CMM_MULTIPOLES_THREAD_VARS<atom> *thread_vars = new CMM_MULTIPOLES_THREAD_VARS<atom>[nThreads];
	int iThread = 0;
	int n1 = 0, n2 = 0;
	int nstep = (n2Cell - n1Cell + 1) / nThreads; n1 = n1Cell;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = n2Cell;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &cell, bEachCluster, bEwaldSum, n1, n2);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)calColumbMultipolesThread<atom>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &calColumbMultipolesThread<atom>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].pcell = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
	
	calCellElectroMultipole<atom>(cell); // the dipole of the whole cell
};
#endif // Multipole Ewald_Sum

//**********************************************************************************
//****  MULTI-THREAD FUNCTIONS FOR CELL OPERATION WITHOUT ADDITIONAL PARAMETERS ****
//**********************************************************************************

//typedef void (*CELL_OPERATOR)(CMM_CELL3D<BASIC_CLUSTER>&, int, int);

template <class atom> class CELL_THREAD_VARS : public _THREAD_VARS {
public:
	typedef void (*CELL_OPERATOR)(CMM_CELL3D<atom>&, int, int);

	CMM_CELL3D<atom> *cell;
	int n1Cell, n2Cell; // the interaction of the cell index in the range [n1Cell, n2Cell]

	//CELL_OPERATOR CellOperator;
	void* CellOperator;

	CELL_THREAD_VARS() {
		cell = NULL; CellOperator = NULL; n1Cell = 0; n2Cell = 0;
	};
	void set_pars(int i, CMM_CELL3D<atom> *c, int n1, int n2) {
		cell = c; n1Cell = n1; n2Cell = n2;
		_THREAD_VARS::set_pars(i);
	};
	//void set_func(CELL_OPERATOR func) {this->CellOperator = (CELL_OPERATOR)func;};
	void set_func(void* func) {this->CellOperator = (void*)func;};
	~CELL_THREAD_VARS() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class atom> int CELL_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class atom> void* CELL_Thread(void *void_par) { // thread-function
#endif
	CELL_THREAD_VARS<atom> *par = (CELL_THREAD_VARS<atom>*)void_par;
	//((CELL_THREAD_VARS<atom>::CELL_OPERATOR)(par->CellOperator))(*(par->cell), par->n1Cell, par->n2Cell);
	typedef void (*CELL_OPERATOR)(CMM_CELL3D<atom>&, int, int);
	((CELL_OPERATOR)(par->CellOperator))(*(par->cell), par->n1Cell, par->n2Cell);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

//void CellOperate(CELL_OPERATOR func, CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, int nThreads) {
//template <class atom> void CellOperate(CELL_OPERATOR func, CMM_CELL3D<atom> &cell, int n1Cell, int n2Cell, int nThreads) {
template <class atom> void CellOperate(void* func, CMM_CELL3D<atom> &cell, int n1Cell, int n2Cell, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*CELL_OPERATOR)(CMM_CELL3D<atom>&, int, int);
		((CELL_OPERATOR)(func))(cell, n1Cell, n2Cell);
		return;
	}

	if (n2Cell >= cell.acell.n) n2Cell = cell.acell.n - 1;
	int nc = 0;
	CELL_THREAD_VARS<atom> *thread_vars = new CELL_THREAD_VARS<atom>[nThreads];
	int iThread = 0;
	int n1 = 0, n2 = 0;
	int nstep = (n2Cell - n1Cell + 1) / nThreads; n1 = n1Cell;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = n2Cell;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &cell, n1, n2);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)CELL_Thread<atom>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &CELL_Thread<atom>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].cell = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
};



//**********************************************************************************
//****  MULTI-THREAD FUNCTIONS FOR CELL OPERATION WITH -- a bool variable  ****
//**********************************************************************************

//typedef void (*CELL_OPERATOR_1)(CMM_CELL3D<BASIC_CLUSTER>&, bool, int, int);

template <class atom, class TVAR> class CELL_THREAD_VARS_1 : public _THREAD_VARS {
public:
	typedef void (*CELL_OPERATOR_1)(CMM_CELL3D<atom>&, TVAR&, int, int);

	CMM_CELL3D<atom> *cell;
	int n1Cell, n2Cell; // the interaction of the cell index in the range [n1Cell, n2Cell]
	TVAR var;

	//CELL_OPERATOR_1 CellOperator;
	void* CellOperator;

	CELL_THREAD_VARS_1() {
		cell = NULL; CellOperator = NULL; n1Cell = 0; n2Cell = 0;
	};
	void set_pars(int i, CMM_CELL3D<atom> *c, TVAR& var, int n1, int n2) {
		cell = c; n1Cell = n1; n2Cell = n2; this->var = var;
		_THREAD_VARS::set_pars(i);
	};
	//void set_func(CELL_OPERATOR_1 func) {this->CellOperator = (CELL_OPERATOR_1)func;};
	void set_func(void* func) {this->CellOperator = (void*)func;};
	~CELL_THREAD_VARS_1() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class atom, class TVAR> int CELL_Thread_1(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class atom, class TVAR> void* CELL_Thread_1(void *void_par) { // thread-function
#endif
	CELL_THREAD_VARS_1<atom, TVAR> *par = (CELL_THREAD_VARS_1<atom, TVAR>*)void_par;
	//((CELL_THREAD_VARS_1<atom, TVAR>::CELL_OPERATOR_1)(par->CellOperator))(*(par->cell), par->var, par->n1Cell, par->n2Cell);
	typedef void (*CELL_OPERATOR)(CMM_CELL3D<atom>&, TVAR, int, int);
	((CELL_OPERATOR)(par->CellOperator))(*(par->cell), par->var, par->n1Cell, par->n2Cell);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

//void CellOperate1(CELL_OPERATOR_1 func, CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, bool bvar, int nThreads);
template <class atom, class TVAR> void CellOperate1(void* func, CMM_CELL3D<atom> &cell, int n1Cell, int n2Cell, TVAR& var, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*CELL_OPERATOR)(CMM_CELL3D<atom>&, TVAR, int, int);
		((CELL_OPERATOR)func)(cell, var, n1Cell, n2Cell); return;
	}
	if (n2Cell >= cell.acell.n) n2Cell = cell.acell.n - 1;
	int nc = 0;
	CELL_THREAD_VARS_1<atom, TVAR> *thread_vars = new CELL_THREAD_VARS_1<atom, TVAR>[nThreads];
	int iThread = 0;
	int n1 = 0, n2 = 0;
	int nstep = (n2Cell - n1Cell + 1) / nThreads; n1 = n1Cell;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = n2Cell;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &cell, var, n1, n2);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)CELL_Thread_1<atom, TVAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &CELL_Thread_1<atom, TVAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].cell = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
};



//**********************************************************************************
//****  MULTI-THREAD FUNCTIONS FOR CELL OPERATION WITH -- return a result  ****
//**********************************************************************************

//typedef void (*CELL_OPERATOR_2)(CMM_CELL3D<BASIC_CLUSTER>&, int, int, double &res);

template <class atom, class TRES> class CELL_THREAD_VARS_2 : public _THREAD_VARS {
public:
	typedef void (*CELL_OPERATOR_2)(CMM_CELL3D<atom>&, int, int, TRES*);

	CMM_CELL3D<atom> *cell;
	int n1Cell, n2Cell; // the interaction of the cell index in the range [n1Cell, n2Cell]
	TRES* res;

	//CELL_OPERATOR_2 CellOperator;
	void* CellOperator;

	CELL_THREAD_VARS_2() {
		cell = NULL; CellOperator = NULL; n1Cell = 0; n2Cell = 0;
	};
	void set_pars(int i, CMM_CELL3D<atom> *c, int n1, int n2, TRES *res) {
		cell = c; n1Cell = n1; n2Cell = n2; this->res = res; 
		_THREAD_VARS::set_pars(i);
	};
	//void set_func(CELL_OPERATOR_2 func) {this->CellOperator = (CELL_OPERATOR_2)func;};
	void set_func(void* func) {this->CellOperator = (void*)func;};
	~CELL_THREAD_VARS_2() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class atom, class TRES> int CELL_Thread_2(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class atom, class TRES> void* CELL_Thread_2(void *void_par) { // thread-function
#endif
	CELL_THREAD_VARS_2<atom, TRES> *par = (CELL_THREAD_VARS_2<atom, TRES>*)void_par;
	//((CELL_THREAD_VARS_1<atom, TVAR>::CELL_OPERATOR_1)(par->CellOperator))(*(par->cell), par->var, par->n1Cell, par->n2Cell);
	typedef void (*CELL_OPERATOR)(CMM_CELL3D<atom>&, int, int, TRES &res);
	((CELL_OPERATOR)(par->CellOperator))(*(par->cell), par->n1Cell, par->n2Cell, *(par->res));
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class atom, class TRES> void CellOperate2(void* func, CMM_CELL3D<atom> &cell, int n1Cell, int n2Cell, TRES& res, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*CELL_OPERATOR)(CMM_CELL3D<atom>&, int, int, TRES &res);
		res = (TRES)0;
		((CELL_OPERATOR)func)(cell, n1Cell, n2Cell, res); return;
	}
	if (n2Cell >= cell.acell.n) n2Cell = cell.acell.n - 1;
	int nc = 0;
	CELL_THREAD_VARS_2<atom, TRES> *thread_vars = new CELL_THREAD_VARS_2<atom, TRES>[nThreads];
	TRES tres[MAX_THREADS];
	int iThread = 0;
	int n1 = 0, n2 = 0;
	int nstep = (n2Cell - n1Cell + 1) / nThreads; n1 = n1Cell;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = n2Cell;
		else {
			nstep = (n2Cell - n1 + 1) / (nThreads - iThread);
			if (nstep * (nThreads - iThread) < (n2Cell - n1 + 1)) nstep += 1;
			n2 = n1 + nstep - 1;
		}

		tres[iThread] = 0;
		thread_vars[iThread].set_pars(iThread, &cell, n1, n2, tres + iThread);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)CELL_Thread_2<atom, TRES>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &CELL_Thread_2<atom, TRES>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	res = (TRES)0;
	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].cell = NULL;
		res += tres[iThread];
	}
	delete[] thread_vars; thread_vars = NULL;
};

//**********************************************************************************
//****  MULTI-THREAD FUNCTIONS FOR CELL OPERATION WITH -- a chain  ****
//**********************************************************************************
/*
template <class atom> void kickoff_cluster_in_subcell(CMM_CELL<atom> *cell, atom *pc) {
	if (cell->ch == NULL) return; // empty cell, subcell is also empty of course

	int ni, nj, nk;
	if (cell->subCell != NULL) {
		for (ni = 0; ni < cell->subCell->n; ni++) {
			for (nj = 0; nj < cell->subCell->n; nj++) {
				for (nk = 0; nk < cell->subCell->n; nk++) {
					kickoff_cluster_in_subcell<atom>(&(cell->subCell->cell[ni][nj][nk]), pc);
				}
			}
		}
	}
	cell->detach(pc);
};
*/
template <class atom> void cmm_check_cluster_in_cell(CMM_CELL3D<atom> *cell, int n1Cell, int n2Cell, CHAIN<atom> **ch) {
	int nc = 0;
	int ni = 0;
	VECTOR<3> dr, rc;
	bool in_cell;
	CHAIN<atom> *ch_tail = NULL, *tch = NULL;
	if (ch != NULL && *ch != NULL) ch_tail = (*ch)->get_tail();

	ITERATOR<atom, CMM_CELL<atom> > it;

	CMM_CELL<atom> *pcell = NULL;
	atom* pc = NULL;
	for (nc = n1Cell; nc <= n2Cell; nc++) { 
		if (nc >= cell->acell.n) break;
		pcell = cell->acell.m[nc]; it.set(pcell); it.start_iterate();
		while (!it.eoi()) {
			pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
			V32V3((*(pc->vp)), rc)  // use vp position
			InCell(rc, pcell->r0, dr, cell->hw, ni, in_cell)
			// current cluster will be moved away to new cell
			// so go to next cluster first
			it.next_iterate();
			if (!in_cell) {
				pcell->detach(pc);

				tch = new CHAIN<atom>; tch->p = pc;
				if (ch_tail == NULL) *ch = tch;
				else ch_tail->next = tch;
				ch_tail = tch;

				//kickoff_cluster_in_subcell<atom>(cell->cell + nc, pc);
			}
		}
	}
};
/*
template <class atom> void distribute_cluster_into_subcell(CMM_CELL<atom> *cell, atom *pc, VECTOR3& rc) {
	cell->attach(pc);

	int ni, nj, nk;
	VECTOR3 dr;
	if (cell->subCell != NULL) {
		SUBCELL_r2INDX(rc, (*(cell->subCell)), ni, nj, nk, dr)
		if (ni > 2) ni = 2; if (ni < 0) ni = 0;
		if (nj > 2) nj = 2; if (nj < 0) nj = 0;
		if (nk > 2) nk = 2; if (nk < 0) nk = 0;
		distribute_cluster_into_subcell<atom>(&(cell->subCell->cell[ni][nj][nk]), pc, rc);
	}
};
*/
template <class atom> void cmm_distribute_cluster(CMM_CELL3D<atom> &cell, CHAIN<atom> *ch) {
	int nm = 0, nc = 0;
	int ni = 0, nj = 0, nk = 0;
	SVECTOR<double, 3> dr, rc;
	int vi = 0;

	CHAIN<atom> *tch = ch;
	atom* pc = NULL;

	while (tch != NULL) {
		pc = tch->p;
		V32V3((*(pc->vp)), rc)  // use vp position

		CMMCELL_r2INDX(rc, cell, ni, nj, nk, dr)
		cell_n(cell, ni, nj, nk, nc)
		if (nc >= cell.acell.n) nc = cell.acell.n - 1;

		cell.acell.m[nc]->attach(pc);

		tch = tch->next;
	}
};

template <class atom> class CHECK_CELL_THREAD_VARS : public _THREAD_VARS {
public:
	typedef void (*CHECK_CELL_OPERATOR)(CMM_CELL3D<atom>*, int, int, CHAIN<atom>**);

	CMM_CELL3D<atom> *cell;
	int n1Cell, n2Cell; // the interaction of the cell index in the range [n1Cell, n2Cell]

	CHAIN<atom> *ch;

	CHECK_CELL_OPERATOR CellOperator;
	//void* CellOperator;

	CHECK_CELL_THREAD_VARS() {
		cell = NULL; CellOperator = NULL; n1Cell = 0; n2Cell = 0; ch = NULL;
	};
	void set_pars(int i, CMM_CELL3D<atom> *c, int n1, int n2) {
		cell = c; n1Cell = n1; n2Cell = n2;
		_THREAD_VARS::set_pars(i);
	};
	void set_func(CHECK_CELL_OPERATOR func) {this->CellOperator = (CHECK_CELL_OPERATOR)func;};
	//void set_func(void* func) {this->CellOperator = (void*)func;};
	~CHECK_CELL_THREAD_VARS() {ch = NULL;};
};

#if _SYS_ == _WINDOWS_SYS_
template <class atom> int CHECK_CELL_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class atom> void* CHECK_CELL_Thread(void *void_par) { // thread-function
#endif
	CHECK_CELL_THREAD_VARS<atom> *par = (CHECK_CELL_THREAD_VARS<atom>*)void_par;
	//((CHECK_CELL_THREAD_VARS<atom>::CHECK_CELL_OPERATOR)(par->CellOperator))(par->cell, par->n1Cell, par->n2Cell, &(par->ch));
	typedef void (*CHECK_CELL_OPERATOR)(CMM_CELL3D<atom>*, int, int, CHAIN<atom>**);
	((CHECK_CELL_OPERATOR)(par->CellOperator))(par->cell, par->n1Cell, par->n2Cell, &(par->ch));
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class atom> void cmm_check_cluster_parallel(CMM_CELL3D<atom> *cell, int n1Cell, int n2Cell, CHAIN<atom> **ch, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*CHECK_CELL_OPERATOR)(CMM_CELL3D<atom>*, int, int, CHAIN<atom>**);
		cmm_check_cluster_in_cell<atom>(cell, n1Cell, n2Cell, ch); return;
	}
	if (n2Cell >= cell->acell.n) n2Cell = cell->acell.n - 1;
	int iThread = 0, nc = 0;
	int n1 = 0, n2 = 0;
	int nstep = (n2Cell - n1Cell + 1) / nThreads; n1 = n1Cell;

	CHECK_CELL_THREAD_VARS<atom> *thread_vars = new CHECK_CELL_THREAD_VARS<atom>[nThreads];

	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = n2Cell;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, cell, n1, n2);
		thread_vars[iThread].set_func(cmm_check_cluster_in_cell<atom>);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)CHECK_CELL_Thread<atom>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &CHECK_CELL_Thread<atom>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	CHAIN<atom> *tail = NULL;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].ch != NULL) {
			if (*ch == NULL) *ch = thread_vars[iThread].ch;
			else {tail = (*ch)->get_tail(); tail->next = thread_vars[iThread].ch;}
			thread_vars[iThread].ch = NULL;
			thread_vars[iThread].cell = NULL;
		}
	}
	delete[] thread_vars; thread_vars = NULL;
};

//void cmm_check(CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads);
#if _CHECK_TIME_STAMP_ == 1
#ifndef _DEF_TIME_STAMP_
#include "time.h"
#endif
#endif

template <class atom> void cmm_check(CMM_CELL3D<atom> &cell, int nThreads) {
	CHAIN<atom> *ch_re_distribt = NULL;

#if _CHECK_TIME_STAMP_ == 1
	TIME time;
	time.start();
#endif

	cmm_check_cluster_parallel<atom>(&cell, 0, cell.acell.n - 1, &ch_re_distribt, nThreads);

#if _CHECK_TIME_STAMP_ == 1
	sprintf(::errmsg, "check-cluster in CMM takes %d ms", time.glance());
	show_infor(errmsg);
#endif
	/*
	char buffer[256] = "\0";
	sprintf(buffer, "re-distributed %d clusters", number(ch_re_distribt));
	show_log(buffer, true);
	*/
	if (ch_re_distribt == NULL) return;

#if _CHECK_TIME_STAMP_ == 1
	time.start();
#endif

	cmm_distribute_cluster<atom>(cell, ch_re_distribt);

#if _CHECK_TIME_STAMP_ == 1
	sprintf(::errmsg, "re-distributing %d clusters, takes %d ms", number(ch_re_distribt), time.glance());
	show_infor(errmsg);
#endif
};


template <class MACROMOL, class atom> class MACROMOL_CELL_THREAD_VARS : public _THREAD_VARS {
public:
	typedef void (*MACROMOL_CELL_OPERATOR)(CMM_CELL3D<atom>*, MACROMOL*, int, int);

	CMM_CELL3D<atom> *cell;
	MACROMOL *mm;
	int nc1, nc2;

	void* Operator;

	MACROMOL_CELL_THREAD_VARS() {
		cell = NULL; mm = NULL; nc1 = 0; nc2 = 0; Operator = NULL;
	};
	void set_pars(int i, CMM_CELL3D<atom> *cell, MACROMOL *mm, int nc1, int nc2) {
		this->cell = cell; this->mm = mm; this->nc1 = nc1; this->nc2 = nc2;
		_THREAD_VARS::set_pars(i);
	};
	void set_func(void* func) {this->Operator = (void*)func;};
	~MACROMOL_CELL_THREAD_VARS() {cell = NULL; mm = NULL;};
};


#if _SYS_ == _WINDOWS_SYS_
template <class MACROMOL, class atom> int MACROMOL_CELL_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MACROMOL, class atom> void* MACROMOL_CELL_Thread(void *void_par) { // thread-function
#endif
	MACROMOL_CELL_THREAD_VARS<MACROMOL, atom> *par = (MACROMOL_CELL_THREAD_VARS<MACROMOL, atom>*)void_par;
	//((MACROMOL_CELL_THREAD_VARS<MACROMOL, atom>::MACROMOL_CELL_OPERATOR)(par->Operator))(par->cell, par->mm, par->nc1, par->nc2);
	typedef void (*MACROMOL_CELL_OPERATOR)(CMM_CELL3D<atom>*, MACROMOL*, int, int);
	((MACROMOL_CELL_OPERATOR)(par->Operator))(par->cell, par->mm, par->nc1, par->nc2);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
}

template <class MACROMOL, class atom> void MacroMolCellOperate(void* func, CMM_CELL3D<atom> &cell, MACROMOL &mm, int nc1, int nc2, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*MACROMOL_CELL_OPERATOR)(CMM_CELL3D<atom>*, MACROMOL*, int, int);
		((MACROMOL_CELL_OPERATOR)(func))(&cell, &mm, nc1, nc2);
		return;
	}

	if (nc2 >= mm.nCluster) nc2 = mm.nCluster - 1;
	int iThread = 0, nc = 0;
	int n1 = 0, n2 = 0;
	int nstep = (nc2 - nc1 + 1) / nThreads; n1 = nc1;

	MACROMOL_CELL_THREAD_VARS<MACROMOL, atom> *thread_vars = new MACROMOL_CELL_THREAD_VARS<MACROMOL, atom>[nThreads];

	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = nc2;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &cell, &mm, n1, n2);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MACROMOL_CELL_Thread<MACROMOL, atom>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MACROMOL_CELL_Thread<MACROMOL, atom>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	delete[] thread_vars; thread_vars = NULL;
}



template <class MACROMOL, class atom, class TVAR> class MACROMOL_CELL_THREAD_VARS_1 : public _THREAD_VARS {
public:
	typedef void (*MACROMOL_CELL_OPERATOR_1)(CMM_CELL3D<atom>*, MACROMOL*, int, int, TVAR&);

	CMM_CELL3D<atom> *cell;
	MACROMOL *mm;
	int nc1, nc2;
	TVAR var;

	//CELL_OPERATOR_1 CellOperator;
	void* CellOperator;

	MACROMOL_CELL_THREAD_VARS_1() {
		mm = NULL; cell = NULL; CellOperator = NULL; nc1 = 0; nc2 = 0;
	};
	void set_pars(int iThread, CMM_CELL3D<atom> *c, MACROMOL* mm, int n1, int n2, TVAR &var) {
		this->mm = mm; this->cell = c; nc1 = n1; nc2 = n2; this->var = var;
		_THREAD_VARS::set_pars(iThread);
	};
	//void set_func(CELL_OPERATOR_1 func) {this->CellOperator = (CELL_OPERATOR_1)func;};
	void set_func(void* func) {this->CellOperator = (void*)func;};
	~MACROMOL_CELL_THREAD_VARS_1() {};
};

#if _SYS_ == _WINDOWS_SYS_
template <class MACROMOL, class atom, class TVAR> int MACROMOL_CELL_Thread_1(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class MACROMOL, class atom, class TVAR> void* MACROMOL_CELL_Thread_1(void *void_par) { // thread-function
#endif
	MACROMOL_CELL_THREAD_VARS_1<MACROMOL, atom, TVAR> *par = (MACROMOL_CELL_THREAD_VARS_1<MACROMOL, atom, TVAR>*)void_par;
	//((MACROMOL_CELL_THREAD_VARS_1<MACROMOL, atom, TVAR>::MACROMOL_CELL_OPERATOR_1)(par->CellOperator))(par->cell, par->mm, par->nc1, par->nc2, par->var);
	typedef void (*MACROMOL_CELL_OPERATOR)(CMM_CELL3D<atom>*, MACROMOL*, int, int, TVAR&);
	((MACROMOL_CELL_OPERATOR)(par->CellOperator))(par->cell, par->mm, par->nc1, par->nc2, par->var);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

//void CellOperate1(CELL_OPERATOR_1 func, CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, bool bvar, int nThreads);
template <class MACROMOL, class atom, class TVAR> void MacroMolCellOperate1(void* func, CMM_CELL3D<atom> &cell, MACROMOL &mm, int nc1, int nc2,  TVAR var, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*MACROMOL_CELL_OPERATOR)(CMM_CELL3D<atom>*, MACROMOL*, int, int, TVAR&);
		((MACROMOL_CELL_OPERATOR)(func))(&cell, &mm, nc1, nc2, var);
		return;
	}

	if (nc2 >= mm.nCluster) nc2 = mm.nCluster - 1;
	int nc = 0;
	MACROMOL_CELL_THREAD_VARS_1<MACROMOL, atom, TVAR> *thread_vars = new MACROMOL_CELL_THREAD_VARS_1<MACROMOL, atom, TVAR>[nThreads];
	int iThread = 0;
	int n1 = 0, n2 = 0;
	int nstep = (nc2 - nc1 + 1) / nThreads; n1 = nc1;
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) n2 = nc2;
		else n2 = n1 + nstep - 1;

		thread_vars[iThread].set_pars(iThread, &cell, &mm, n1, n2, var);
		thread_vars[iThread].set_func(func);
		n1 += nstep;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MACROMOL_CELL_Thread_1<MACROMOL, atom, TVAR>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		CloseHandle(thread_vars[iThread].hMMThread); thread_vars[iThread].hMMThread = NULL;
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MACROMOL_CELL_Thread_1<MACROMOL, atom, TVAR>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	for (iThread = 0; iThread < nThreads; iThread++) {
		thread_vars[iThread].cell = NULL;
	}
	delete[] thread_vars; thread_vars = NULL;
};




template <class SM, class CMM_CLUSTER> class SM_CELL_THREAD_VARS : public _THREAD_VARS {
public:
	CMM_CELL3D<CMM_CLUSTER> *cell;
	SM *sm;
	int nSM;

	void* Operator;

	SM_CELL_THREAD_VARS() {
		cell = NULL; sm = NULL; nSM = 0; Operator = NULL;
	};
	void set_pars(int iThread, CMM_CELL3D<CMM_CLUSTER> *cell, SM *sm, int nSM) {
		this->cell = cell; this->sm = sm; this->nSM = nSM;
		_THREAD_VARS::set_pars(iThread);
	};
	void set_func(void* func) {this->Operator = (void*)func;};
	~SM_CELL_THREAD_VARS() {cell = NULL; sm = NULL; nSM = 0;};
};

#if _SYS_ == _WINDOWS_SYS_
template <class SM, class CMM_CLUSTER> int SM_CELL_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class SM, class CMM_CLUSTER> void* SM_CELL_Thread(void *void_par) { // thread-function
#endif
	typedef void (*SM_CELL_OPERATOR)(SM *sm, int, CMM_CELL3D<CMM_CLUSTER>*);
	SM_CELL_THREAD_VARS<SM, CMM_CLUSTER> *par = (SM_CELL_THREAD_VARS<SM, CMM_CLUSTER>*)void_par;
	((SM_CELL_OPERATOR)(par->Operator))(par->sm, par->nSM, par->cell);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class SM, class CMM_CLUSTER> void SMCellOperate(void* func, CMM_CELL3D<CMM_CLUSTER> &cell, SM *sm, int nSM, int nThreads) {
	if (nThreads <= 1) {
		typedef void (*SM_CELL_OPERATOR)(SM *sm, int, CMM_CELL3D<CMM_CLUSTER>*);
		((SM_CELL_OPERATOR)(func))(sm, nSM, &cell);
		return;
	}

	int nstep = nSM / nThreads;
	if (nstep < 1) {nstep = 1; nThreads = nSM;}
	int iThread = 0;
	int nstart = 0, nm = 0;

	SM_CELL_THREAD_VARS<SM, CMM_CLUSTER> *thread_vars = new SM_CELL_THREAD_VARS<SM, CMM_CLUSTER>[nThreads];

	for (iThread = 0; iThread < nThreads; iThread++) {
		if (iThread == nThreads - 1) nm = nSM - nstart;
		else nm = nstep;

		thread_vars[iThread].set_pars(iThread, &cell, sm + nstart, nm);
		thread_vars[iThread].set_func(func);
		nstart += nm;
	}

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)SM_CELL_Thread<SM, CMM_CLUSTER>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &SM_CELL_Thread<SM, CMM_CLUSTER>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
	}
#endif

	delete[] thread_vars; thread_vars = NULL;
};

namespace _cmm_3d_ {

/************************************************************************************
                        implicit solvation for clusters
************************************************************************************/
// pc0 is the center cluster, nc0[3] is the cell's (where pc0 is) index in basic_cell_matrix
// dnc[3] is the offset of cell index, where the cluster inside could be neighbors
// xd[3] is the dimension of basic_cell
template <class atom> void CheckImpSolvNeighbors_Cluster(CMM_CELL3D<atom> &cell, atom *pc0, int *nc0, int *dnc, float r0_solvent, IMAG_RELSHIP<atom> &impsolv_rship) {
	bool imagine_cell = cell.periodic_cell;
	int i, j, k, ni, nj, nk, nx, ny, nz;
	SVECTOR<double, 3> r0, dr, r1;

	if (pc0->bCMMShift) {V3plusV3(pc0->rc_solv, pc0->dr_cmm, r0)} // r0 is pc0's rc
	else {V32V3(pc0->rc_solv, r0)}
	
	float d2, rcut = 0, r2cut;
	CMM_CELL<atom> *pcell = NULL;
	atom *pc = NULL;
	int n_relationships = 0;
	ITERATOR<atom, CMM_CELL<atom> > iterator;

	for (i = nc0[0] - dnc[0]; i <= nc0[0] + dnc[0]; i++) {
		ni = i; nx = 0;
		if (ni < 0) {
			if (imagine_cell) {while (ni < 0) {ni += cell.nx; nx -= 1;}}
			else continue;
		}
		else if (ni >= cell.bcell.nx) {
			if (imagine_cell) {while (ni >= cell.bcell.nx) {ni -= cell.bcell.nx; nx += 1;}}
			else continue;
		}
		for (j = nc0[1] - dnc[1]; j <= nc0[1] + dnc[1]; j++) {
			nj = j; ny = 0;
			if (nj < 0) {
				if (imagine_cell) {while (nj < 0) {nj += cell.bcell.ny; ny -= 1;}}
				else continue;
			}
			else if (nj >= cell.bcell.ny) {
				if (imagine_cell) {while (nj >= cell.bcell.ny) {nj -= cell.bcell.ny; ny += 1;}}
				else continue;
			}
			for (k = nc0[2] - dnc[2]; k <= nc0[2] + dnc[2]; k++) {
				nk = k; nz = 0;
				if (nk < 0) {
					if (imagine_cell) {while (nk < 0) {nk += cell.bcell.nz; nz -= 1;}}
					else continue;
				}
				else if (nk >= cell.bcell.nz) {
					if (imagine_cell) {while (nk >= cell.bcell.nz) {nk -= cell.bcell.nz; nz += 1;}}
					else continue;
				}
				pcell = &(cell.bcell.m[ni].m[nj].m[nk]);
				// this the deepest level
				iterator.set(pcell); iterator.start_iterate();
				while (!iterator.eoi()) {
					pc = iterator.current(); if (pc == NULL) {iterator.next_iterate(); continue;}
					if (pc == pc0) {
						if (!imagine_cell) {
							iterator.next_iterate();
							continue;
						}
						else if (nx == 0 && ny == 0 && nz == 0) {
							iterator.next_iterate();
							continue;
						}
					}
					if (imagine_cell) {
						if (pc->bCMMShift) {V3plusV3(pc->rc_solv, pc->dr_cmm, r1)}
						else {V32V3(pc->rc_solv, r1)}
						if (nx != 0) r1.v[0] += nx * cell.xd[0];
						if (ny != 0) r1.v[1] += ny * cell.xd[1];
						if (nz != 0) r1.v[2] += nz * cell.xd[2];
					}
					else {V32V3(pc->rc_solv, r1)}
					VECT3(r0, r1, dr) // dr = r1 - r0
					V3ABS2(dr, d2)
					rcut = pc0->psolv->r0 + pc->psolv->r0 + r0_solvent + r0_solvent;
					r2cut = rcut * rcut;
					if (d2 <= r2cut) {
						if (::local_relax) {
							if (d2 + d2 < r2cut) {
								impsolv_rship.attach(pc, nx, ny, nz);
								n_relationships++;
								if (n_relationships > 20) return; // local_relax case: do not consider too much
							}
						}
						else impsolv_rship.attach(pc, nx, ny, nz);
					}
					iterator.next_iterate();
				}
			}
		}
	}
};

template <class atom> void ClusterImpSolvNeighborsInCell(CMM_CELL3D<atom> *cell, atom *pc, int *dnc) {
	// nc is the subcell index of the cluster
	// dnc is the neighbor subcell to search around
	int nc[3];
	if (pc->psolv == NULL) return;
	SVECTOR<double, 3> r0, dr;
	// position of cluster
	V32V3((*(pc->vp)), r0)
	// cell index where cluster in the basic_cell_matrix
	CMMCELL_r2INDX(r0, (*cell), nc[0], nc[1], nc[2], dr)

	pc->impsolv_rship.release_relationship();
	CheckImpSolvNeighbors_Cluster<atom>(*cell, pc, nc, dnc, ::rSolvent, (*((IMAG_RELSHIP<atom>*)&(pc->impsolv_rship))));
};

template <class atom> void CheckImpSolvNeighbors_Cell(CMM_CELL3D<atom> &cell, int nc1, int nc2) {
	float r0_solvent = rSolvent;
	float xd[3] = {cell.xd[0] / cell.bcell.nx, cell.xd[1] / cell.bcell.ny, cell.xd[2] / cell.bcell.nz};
	int nc[3], dnc[3];
	float rcut = dSolv_max;

	if (::local_relax) {dnc[0] = 1; dnc[1] = 1; dnc[2] = 1;}
	else {
		dnc[0] = int(rcut / xd[0] + 0.5f); if (dnc[0] == 0) dnc[0] = 1; if (dnc[0] * xd[0] < rcut) dnc[0] += 1;
		dnc[1] = int(rcut / xd[1] + 0.5f); if (dnc[1] == 0) dnc[1] = 1; if (dnc[1] * xd[1] < rcut) dnc[1] += 1;
		dnc[2] = int(rcut / xd[2] + 0.5f); if (dnc[2] == 0) dnc[2] = 1; if (dnc[2] * xd[2] < rcut) dnc[2] += 1;
		dnc[0] += 1; dnc[1] += 1; dnc[2] += 1;
	}

	CMM_CELL<atom> *pcell = NULL;
	atom *pc = NULL;
	int i, j, k;
	ITERATOR<atom, CMM_CELL<atom> > iterator;
	for (i = nc1; i <= nc2; i++) {
		if (i >= cell.bcell.nx) break;
		for (j = 0; j < cell.bcell.ny; j++) {for (k = 0; k < cell.bcell.nz; k++) {
			pcell = &(cell.bcell.m[i].m[j].m[k]); iterator.set(pcell);
			nc[0] = i; nc[1] = j; nc[2] = k;
			iterator.start_iterate();
			while (!iterator.eoi()) {
				pc = iterator.current();
				if (pc == NULL) {iterator.next_iterate(); continue;}
				if (pc->psolv == NULL) continue;
				pc->impsolv_rship.release_relationship();
				CheckImpSolvNeighbors_Cluster<atom>(cell, pc, nc, dnc, r0_solvent, pc->impsolv_rship);
				iterator.next_iterate();
			}
		}}
	}
	return;
};

// check the implicit-solvation neighbors of each cluster in cell, in parallel
/*
template <class atom> void CheckImpSolvNeighborsInCell(CMM_CELL3D<atom> &cell) {
	CellOperate<atom>((void*)(&CheckImpSolvNeighbors_Cell<atom>), cell, 0, cell.bcell.nx, MAX_THREADS);
};
*/

void CMMCheckImpSolvNeighbors_Cell(CMM_CELL3D<BASIC_CLUSTER> &cell, int nc1, int nc2);


template <class MACROMOL, class atom> void MacroMolClusterImpSolvNeighborsInCell(CMM_CELL3D<atom> *cell, MACROMOL *mm, int nc1, int nc2) {
	atom *pc = NULL;
	//float r0_solvent = rSolvent;
	float xd[3] = {cell->xd[0] / cell->bcell.nx, cell->xd[1] / cell->bcell.ny, cell->xd[2] / cell->bcell.nz};
	int dnc[3];
	float rcut = dSolv_max;

	SVECTOR<double, 3> r0;
	
	if (::local_relax) {dnc[0] = 1; dnc[1] = 1; dnc[2] = 1;}
	else {
		dnc[0] = int(rcut / xd[0] + 0.5f); if (dnc[0] == 0) dnc[0] = 1; if (dnc[0] * xd[0] < rcut) dnc[0] += 1;
		dnc[1] = int(rcut / xd[1] + 0.5f); if (dnc[1] == 0) dnc[1] = 1; if (dnc[1] * xd[1] < rcut) dnc[1] += 1;
		dnc[2] = int(rcut / xd[2] + 0.5f); if (dnc[2] == 0) dnc[2] = 1; if (dnc[2] * xd[2] < rcut) dnc[2] += 1;
		dnc[0] += 1; dnc[1] += 1; dnc[2] += 1;
	}

	for (int mnc = nc1; mnc <= nc2; mnc++) {
		if (mnc >= mm->nCluster) break;
		pc = (atom*)(mm->cluster + mnc);
		if (pc->psolv == NULL) continue;
		// position of cluster
		V32V3((*(pc->vp)), r0)
		
		ClusterImpSolvNeighborsInCell<atom>(cell, pc, dnc);
	}
};

void MMClusterImpSolvNeighborsInCell(CMM_CELL3D<BASIC_CLUSTER> *cell, MMOLECULE *mm, int nc1, int nc2);

void CheckClusterRelship(BASIC_CLUSTER *pc, CMM_CELL3D<BASIC_CLUSTER> &cell);
void ParallelCheckRelshipInCells(CMM_CELL3D<BASIC_CLUSTER> &cell, int n1Cell, int n2Cell, int nThreads);

// check the relationship of a cluster in imagine cells
void MMCheckRelshipInCell(MMOLECULE &mm, int nc1, int nc2, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads);

void SMRelshipInCell(SMOLECULE *sm, int nSM, CMM_CELL3D<BASIC_CLUSTER> *cell);
void SMCheckRelshipInCell(SMOLECULE *sm, int nSM, CMM_CELL3D<BASIC_CLUSTER> &cell, int nThreads);

} // end of namespace _cmm_3d_


namespace _EwaldSum_recp_ {
	void InitEwaldSumPars(float V, double &inv_3V, double &kappa, double &kappa3);
} // end of namespace _EwaldSum_recp_ 

float cluster_max_size(CMM_CELL3D<BASIC_CLUSTER> &cell);

void cmm_init_distribute_cluster(MMOLECULE &mm, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain);
void cmm_init_distribute_cluster(SMOLECULE* sm, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain);
void cmm_init_distribute_cluster(PMOLECULE* pm, int n, CMM_CELL3D<BASIC_CLUSTER> &cell, bool reset_cell_chain);

void calClusterElectroMultipole(CMM_CELL3D<BASIC_CLUSTER> &cell, bool bEwaldSum, int n1Cell, int n2Cell);


/*********************************************************************************
*********************************************************************************/
template <class TCLUSTER, class TMM, class TSM, class TPM, class CELL, class MDCELL> class MOperateVAR {
public:
	TMM *mm;
	TSM *sm;
	TPM *pm;
	int nMM, nSM, nPM;
	CELL *cell;
	MDCELL *mmcell;

	int nMainThreads, nMinorThreads;
	MOperateVAR() {mm = NULL; sm = NULL; pm = NULL; nMM = 0; nSM = 0; nPM = 0; cell = NULL; mmcell = NULL; nMainThreads = 1; nMinorThreads = 1;};
};

// IMPORTANT: TVAR needs an operator =; TRES needs two functions, initialize() and accumulate(TRES&)
template <class TCLUSTER, class TMM, class TSM, class TPM, class CELL, class MDCELL, class TVAR, class TRES> class MOperateThreadVAR : public _THREAD_VARS {
public:
	MOperateVAR<TCLUSTER, TMM, TSM, TPM, CELL, MDCELL> mv;
	TVAR v;
	TRES res;

	void *func;

	MOperateThreadVAR() {func = NULL;};
};

#if _SYS_ == _WINDOWS_SYS_
template <class TCLUSTER, class TMM, class TSM, class TPM, class CELL, class MDCELL, class TVAR, class TRES> int MOperateInCell_Thread(LPVOID *void_par) { // thread-function
#elif _SYS_ == _LINUX_SYS_
template <class TCLUSTER, class TMM, class TSM, class TPM, class CELL, class MDCELL, class TVAR, class TRES> void* MOperateInCell_Thread(void *void_par) { // thread-function
#endif
	MOperateThreadVAR<TCLUSTER, TMM, TSM, TPM, CELL, MDCELL, TVAR, TRES> *par = (MOperateThreadVAR<TCLUSTER, TMM, TSM, TPM, CELL, MDCELL, TVAR, TRES>*)void_par;
	typedef void (*MInCell_OPERATOR)(TMM*, int, TSM*, int, TPM*, int, CELL*, MDCELL*, int, TVAR&, TRES&);
	((MInCell_OPERATOR)(par->func))(par->mv.mm, par->mv.nMM, par->mv.sm, par->mv.nSM, par->mv.pm, par->mv.nPM, par->mv.cell, par->mv.mmcell, par->mv.nMinorThreads, par->v, par->res);
#if _SYS_ == _WINDOWS_SYS_
	SetEvent(par->hMMThread);
	return 1;
#elif _SYS_ == _LINUX_SYS_
	return NULL;
#endif
};

template <class TCLUSTER, class TMM, class TSM, class TPM, class CELL, class MDCELL, class TVAR, class TRES> void MOperateInCell(MOperateVAR<TCLUSTER, TMM, TSM, TPM, CELL, MDCELL> &mvar, void* func, TVAR& v, TRES& res) {
	int nMMstep = mvar.nMM / mvar.nMainThreads, nSMstep = mvar.nSM / mvar.nMainThreads, nPMstep = mvar.nPM / mvar.nMainThreads;
	int nMMThreads = 0, nSMThreads = 0, nPMThreads;
	if (nSMstep < 10) nSMstep = 10;
	if (nPMstep < 10) nPMstep = 10;

	if (mvar.nMM <= 0 || mvar.mm == NULL) {nMMstep = 1; nMMThreads = 0;}
	else {
		if (nMMstep < 1) {nMMstep = 1; nMMThreads = mvar.nMM;}
		else {
			nMMThreads = mvar.nMainThreads;
			nMMstep = mvar.nMM / mvar.nMainThreads;
			if ((nMMstep * nMMThreads + 1) < mvar.nMM) nMMstep += 1;
		}
	}
	
	if (mvar.nSM <= 0 || mvar.sm == NULL) {nSMstep = 1; nSMThreads = 0;}
	else {
		if (mvar.nSM < 10) {
			nSMstep = 10; nSMThreads = 1;
		}
		else {
			nSMThreads = mvar.nMainThreads;
			nSMstep = mvar.nSM / mvar.nMainThreads;
			if ((nSMstep * nSMThreads + 1) < mvar.nSM) nSMstep += 1;
		}
	}

	if (mvar.nPM <= 0 || mvar.pm == NULL) {nPMstep = 1; nPMThreads = 0;}
	else {
		if (nPMstep < 10) {nPMstep = 10; nPMThreads = 1;}
		else {
			nPMThreads = mvar.nMainThreads;
			nPMstep = mvar.nPM / mvar.nMainThreads;
			if ((nPMstep * nPMThreads + 1) < mvar.nPM) nPMstep += 1;
		}
	}

	int iThread = 0, nThreads = (nMMThreads > nSMThreads ? nMMThreads : nSMThreads);
	if (nThreads <= 0) return;
	int nMMstart = 0, nMMleft = mvar.nMM;
	int nSMstart = 0, nSMleft = mvar.nSM;
	int nPMstart = 0, nPMleft = mvar.nPM;

	MOperateThreadVAR<TCLUSTER, TMM, TSM, TPM, CELL, MDCELL, TVAR, TRES> *thread_vars = NULL;
	thread_vars = new MOperateThreadVAR<TCLUSTER, TMM, TSM, TPM, CELL, MDCELL, TVAR, TRES>[nThreads];

	for (iThread = 0; iThread < nThreads; iThread++) {
		if (mvar.mm != NULL) {
			if (nMMstart < mvar.nMM) {
				thread_vars[iThread].mv.mm = mvar.mm + nMMstart;
				if (iThread == nThreads - 1) {
					thread_vars[iThread].mv.nMM = nMMleft;
				}
				else {
					thread_vars[iThread].mv.nMM = (nMMstep <= nMMleft ? nMMstep : nMMleft);
				}
			}
			else {
				thread_vars[iThread].mv.mm = NULL;
				thread_vars[iThread].mv.nMM = 0;
			}
			nMMstart += thread_vars[iThread].mv.nMM;
			nMMleft -= thread_vars[iThread].mv.nMM;
		}
		else {
			thread_vars[iThread].mv.mm = NULL;
			thread_vars[iThread].mv.nMM = 0;
		}

		if (mvar.sm != NULL) {
			if (nSMstart < mvar.nSM) {
				thread_vars[iThread].mv.sm = mvar.sm + nSMstart;
				if (iThread == nThreads - 1) {
					thread_vars[iThread].mv.nSM = nSMleft;
				}
				else {
					thread_vars[iThread].mv.nSM = (nSMstep <= nSMleft ? nSMstep : nSMleft);
				}
			}
			else {
				thread_vars[iThread].mv.sm = NULL;
				thread_vars[iThread].mv.nSM = 0;
			}
			nSMstart += thread_vars[iThread].mv.nSM;
			nSMleft -= thread_vars[iThread].mv.nSM;
		}
		else {
			thread_vars[iThread].mv.sm = NULL;
			thread_vars[iThread].mv.nSM = 0;
		}

		if (mvar.pm != NULL) {
			if (nPMstart < mvar.nPM) {
				thread_vars[iThread].mv.pm = mvar.pm + nPMstart;
				if (iThread == nThreads - 1) {
					thread_vars[iThread].mv.nPM = nPMleft;
				}
				else {
					thread_vars[iThread].mv.nPM = (nPMstep <= nPMleft ? nPMstep : nPMleft);
				}
			}
			else {
				thread_vars[iThread].mv.pm = NULL;
				thread_vars[iThread].mv.nPM = 0;
			}
			nPMstart += thread_vars[iThread].mv.nPM;
			nPMleft -= thread_vars[iThread].mv.nPM;
		}
		else {
			thread_vars[iThread].mv.pm = NULL;
			thread_vars[iThread].mv.nPM = 0;
		}
		
		thread_vars[iThread].mv.cell = mvar.cell;
		thread_vars[iThread].mv.mmcell = mvar.mmcell;
		thread_vars[iThread].v = v;
		thread_vars[iThread].mv.nMainThreads = mvar.nMainThreads;
		thread_vars[iThread].mv.nMinorThreads = mvar.nMinorThreads;

		thread_vars[iThread].set_pars(iThread);
		thread_vars[iThread].func = func;
		thread_vars[iThread].res.initialize();
	}

	res.initialize();

#if _SYS_ == _WINDOWS_SYS_ 
	for (iThread = 0; iThread < nThreads; iThread++) {
		ResetEvent(thread_vars[iThread].hMMThread);
		thread_vars[iThread].thread = AfxBeginThread((AFX_THREADPROC)MOperateInCell_Thread<TCLUSTER, TMM, TSM, TPM, CELL, MDCELL, TVAR, TRES>, (LPVOID)(thread_vars + iThread), THREAD_PRIORITY_IDLE); // THREAD_PRIORITY_IDLE
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		WaitForSingleObject(thread_vars[iThread].hMMThread, INFINITE);
		res.accumulate(thread_vars[iThread].res);
	}
#elif _SYS_ == _LINUX_SYS_
	for (iThread = 0; iThread < nThreads; iThread++) {
		pthread_create(&(thread_vars[iThread].thread), NULL, &MOperateInCell_Thread<TCLUSTER, TMM, TSM, TPM, CELL, MDCELL, TVAR, TRES>, (void *)(thread_vars + iThread));
	}
	for (iThread = 0; iThread < nThreads; iThread++) {
		if (thread_vars[iThread].thread != 0) pthread_join(thread_vars[iThread].thread, NULL); // wait the thread over
		res.accumulate(thread_vars[iThread].res);
	}
#endif

	delete[] thread_vars; thread_vars = NULL;
};






/**************************************************************************************************************************
	following functions are specially designed for the CMM using ARRAY to store the cluster inside each CMM_CELL
**************************************************************************************************************************/
template <class atom> bool combine_AtomAllocatedCMM_array(CMM_CELL3D<atom> *cell, int ncell, CMM_CELL3D<atom> *res_cell, bool bReplaceResCell = true, bool bResetCell = true) {
	char msg[256] = "\0";
	bool status = true;
	int ix, iy, iz, icell, ic;
	CMM_CELL<atom> *res = NULL, *pcell = NULL;
	ITERATOR<atom, CMM_CELL<atom> > it;
	atom *pc = NULL;
	int size_atom = sizeof(atom*);
	for (ix = 0; ix < res_cell->bcell.nx; ix++) { for (iy = 0; iy < res_cell->bcell.ny; iy++) {for (iz = 0; iz < res_cell->bcell.nz; iz++) {
		res = &(res_cell->bcell.m[ix].m[iy].m[iz]); if (bReplaceResCell) res->release();
		for (icell = 0; icell < ncell; icell++) {
			pcell = &(cell[icell].bcell.m[ix].m[iy].m[iz]); 
			if (pcell->bChain) {
				it.set(pcell); it.start_iterate();
				for (ic = 0; ic < pcell->length(); ic++) {
					pc = it.current(); if (pc == NULL) {it.next_iterate(); continue;}
					if (!res->attach(pc)) {
						status = false; sprintf(msg, "subcell [%d, %d, %d] leaks", ix, iy, iz); show_log(msg, true); break;
					}
					it.next_iterate();
				}
			}
			else { // in case using array 
				if (pcell->nbs == 0) continue;
				res->require_memory(pcell->nbs);
				memcpy(res->ba.m + res->nbs, pcell->ba.m, pcell->nbs * size_atom);
				res->nbs += pcell->nbs;
			}
			
			if (bResetCell) pcell->release();
		}
	}}}
	return status;
};

template <class atom, class MM> bool MM_AllocateCMM_array(MM *mm, CMM_CELL3D<atom> *cell) {
	char msg[256] = "\0";
	bool status = true;
	int ix = 0, iy = 0, iz = 0, ic;
	atom *pc = NULL;
	CMM_CELL<atom> *pcell = NULL;
	SVECTOR<double, 3> rc, dr;
	for (ic = 0; ic < mm->nCluster; ic++) {
		pc = (atom*)(mm->cluster + ic);
		V32V3((*(pc->vp)), rc) // use vp position
		CMMCELL_r2INDX(rc, (*cell), ix, iy, iz, dr)
		pcell = &(cell->bcell.m[ix].m[iy].m[iz]);
		if (!pcell->attach(pc)) {
			sprintf(msg, "subcell [%d, %d, %d] leaks", ix, iy, iz); show_log(msg, true); status = false;
		}
	}
	return status;
};

template <class atom, class MM> bool MM_ClusterAllocateCMM_array(CMM_CELL3D<atom> *cell, MM *mm, int nc1, int nc2) {
	char msg[256] = "\0";
	bool status = true;
	int ix = 0, iy = 0, iz = 0, ic;
	atom *pc = NULL;
	CMM_CELL<atom> *pcell = NULL;
	SVECTOR<double, 3> rc, dr;
	for (ic = nc1; ic <= nc2; ic++) {
		if (ic >= mm->nCluster) break;
		pc = (atom*)(mm->cluster + ic);
		V32V3((*(pc->vp)), rc) // use vp position
		CMMCELL_r2INDX(rc, (*cell), ix, iy, iz, dr)
		pcell = &(cell->bcell.m[ix].m[iy].m[iz]);
		if (!pcell->attach(pc)) {
			sprintf(msg, "subcell [%d, %d, %d] leaks", ix, iy, iz); show_log(msg, true); status = false;
		}
	}
	return status;
};

template <class atom> void CMM_array_reset_storage(CMM_CELL3D<atom> *cell) {
	int ix, iy, iz;
	for (ix = 0; ix < cell->bcell.nx; ix++) { for (iy = 0; iy < cell->bcell.ny; iy++) {for (iz = 0; iz < cell->bcell.nz; iz++) {
		cell->bcell.m[ix].m[iy].m[iz].release();
	}}}
};

template <class atom> void CMM_array_set_storage(CMM_CELL3D<atom> *cell, int nmax) {
	int ix, iy, iz;
	cell->bChain= false;
	for (ix = 0; ix < cell->bcell.nx; ix++) { for (iy = 0; iy < cell->bcell.ny; iy++) {for (iz = 0; iz < cell->bcell.nz; iz++) {
		cell->bcell.m[ix].m[iy].m[iz].use_array(nmax);
	}}}
};

template <class atom> void CMM_construct(CMM_CELL3D<atom> *dest, CMM_CELL3D<atom> *source) {
	dest->set_cell(source->nx, source->ny, source->nz);
	dest->set_cell_range(source->periodic_cell, source->xl[0], source->xr[0], source->xl[1], source->xr[1], source->xl[2], source->xr[2]);
	dest->set_cell_pos();
};

template <class atom> void CMM_array_construct_for_storage(CMM_CELL3D<atom> *dest, CMM_CELL3D<atom> *source, int nmax) {
	CMM_construct<atom>(dest, source);
	dest->bChain = false;
	for (int i = 0; i < dest->acell.n; i++) dest->acell.m[i]->use_array(nmax);
};

bool MM_AllocateCluster_CMM_array(MMOLECULE *mm, int nmol, CMM_CELL3D<BASIC_CLUSTER> *cell);
bool SM_AllocateCluster_CMM_array(SMOLECULE *sm, int nmol, CMM_CELL3D<BASIC_CLUSTER> *cell);
bool PM_AllocateCluster_CMM_array(PMOLECULE *pm, int nmol, CMM_CELL3D<BASIC_CLUSTER> *cell);
bool AllocateCluster_CMM_array(ARRAY<BASIC_CLUSTER*> &bc, CMM_CELL3D<BASIC_CLUSTER> *cell);



/***********************************************************************************************
	MULTI-THREAD Operator on Multi-Clusters, with CMM
***********************************************************************************************/

template <class atom> class CMM_AtomJob {
public:
	ARRAY<atom*> atoms;
	CMM_CELL3D<atom> *cmm;

	CMM_AtomJob() {
		this->cmm = NULL;
	};
	void release() {
		this->cmm = NULL;
		if (atoms.m != NULL) {
			memset(atoms.m, 0, atoms.n * sizeof(atom*));
			atoms.release();
		}
	};
	~CMM_AtomJob() {release();};
};

template <class atom> void AssignAtomJob_CMM(ARRAY<atom*> &atoms, CMM_AtomJob<atom> *job, int nJobs) {
	int i = 0, ia;
	int n1, n2, nstep, n;
	for (i = 0; i < nJobs; i++) {job[i].release();}
	if (atoms.n > 0) {
		n1 = 0; nstep = atoms.n / nJobs;
		if (nstep == 0 || (nstep * nJobs + 1) < atoms.n) nstep += 1;
		for (i = 0; i < nJobs; i++) {
			n2 = n1 + nstep - 1; if (i == nJobs - 1 || n2 >=atoms.n) n2 = atoms.n - 1;
			n = n2 - n1 + 1; if (n <= 0) break;
			job[i].atoms.SetArray(n);
			memcpy(job[i].atoms.m, atoms.m + n1, n * sizeof(atom*));
			n1 = n2 + 1;
		}
	}
};

void AllocateCluster_CMM(JOB_THREAD_VARS<CMM_AtomJob<BASIC_CLUSTER> > *job);
