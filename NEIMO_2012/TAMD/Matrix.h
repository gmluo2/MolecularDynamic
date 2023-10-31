
// Gauss-Jordan elimination with full pivoting

/********************************************************************
Linear equation solution by Gauss-Jordan elimination. 
a[][] is the input matrix; b[] is the right-hand side vector.
On output, a[][] is replaced by its matrix inverse, 
and b[] is replaced b the corresponding set of solution vectors.
********************************************************************/

extern char errmsg[2560];

// Linear equation solution by Gauss-Jordan elimination
// Output -- a is replaced by its matrix inverse supposely, but actually, it does not
//           and b is replaced by the corresponding set of solution vector
template <int N> bool GaussJ(DMATRIX<N> &a, DVECTOR<N> &b) // matrix a -- n x n , b -- n 
{
	if (a.n != b.n) {
		sprintf(errmsg, "A & B matrixes do not match");
		return false;
	}
	int n = a.n, m = 1.0;

	int *indxc = NULL, *indxr = NULL, *ipiv = NULL;
	int i = 0, irow = 0, icol = 0, j = 0, k = 0, l = 0, ll = 0;
	double big = 0, dum = 0, pivinv = 0, temp = 0;

	indxc = new int[n]; for (i = 0; i < n; i++) indxc[i] = -1;
	indxr = new int[n]; for (i = 0; i < n; i++) indxr[i] = -1;
	ipiv = new int[n];  for (j = 0; j < n; j++) ipiv[j] = -1;
	for (i = 0; i < n; i++) {
		big = 0;
		// find the pivot element
		for (j = 0; j < n; j++) {
			if (ipiv[j] != 0) {
				for (k = 0; k < n; k++) {
					if (ipiv[k] == -1) { // k column is not checked yet
						if (fabs(a.m[j][k]) > big) {
							big = fabs(a.m[j][k]);
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 0) {
						// this column was labeled several times
						sprintf(errmsg, "GaussJ : Singular Matrix - 1");
						delete[] indxc; indxc = NULL;
						delete[] indxr; indxr = NULL;
						delete[] ipiv; ipiv = NULL;
						return false;
					}
				}
			}
		}
		(ipiv[icol])++; // label that this row is checked
		/****************************************************************************************
		we now have the pivot element, so we interchange rows, if needed, to put the pivot
		element on the diagonal. The columns are not physically interchanged, only labeled:
		indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
		indxr[i] is the row in which that pivot element was originally located. If indxr[i] !=
		indxc[i] there is an implied column interchanged. With this form of bookkeeping, the
		solution b will end up in the correct order, and the inverse matrix will be scrambled
		by columns.
		****************************************************************************************/
		if (irow != icol) {
			// move the pivot element to diagonal by interchanging the rows
			//for (l = 0; l < n; l++) SWAP(a[irow][l], a[icol][l])
			//SWAP(b[irow], b[icol]);
			for (l = 0; l < n; l++) SWAP(a.m[irow][l], a.m[icol][l], temp)
			SWAP(b.v[irow], b.v[icol], temp)
		}
		indxr[i] = irow;
		indxc[i] = icol;
		//if (a[icol][icol] == 0) {sprintf(errmsg, "GaussJ : Singular Matrix - 2"); return false;}
		if (fabs(a.m[icol][icol]) <= 1e-10) {
			sprintf(errmsg, "GaussJ : Singular Matrix - 2"); 
			delete[] indxc; indxc = NULL;
			delete[] indxr; indxr = NULL;
			delete[] ipiv; ipiv = NULL;
			return false;
		}
		pivinv = 1.0 / a.m[icol][icol]; //pivinv = 1 / a[icol][icol];
		a.m[icol][icol] = 1.0; //a[icol][icol] = 1;
		for (l = 0; l < n; l++) a.m[icol][l] *= pivinv; //a[icol][l] *= pivinv;
		b.v[icol] *= pivinv; //for (l = 0; l < m; l++) b[icol][l] *= pivinv;
		for (ll = 0; ll < n; ll++) {
			if (ll != icol) {
				dum = a.m[ll][icol]; //dum = a[ll][icol];
				//if (dum == 0) continue; // obviously, this line is fine
				a.m[ll][icol] = 0; //a[ll][icol] = 0;
				for (l = 0; l < n; l++) a.m[ll][l] -= a.m[icol][l] * dum; //a[ll][l] -= a[icol][l] * dum;
				b.v[ll] -= b.v[icol] * dum; //b[ll] -= b[icol] * dum;
			}
		}
	}
	
	for (l = n - 1; l >= 0; l--) {
		if (indxr[l] != indxc[l]) 
			for (k = 0; k < n; k++) SWAP(a.m[k][indxr[l]], a.m[k][indxc[l]], temp); //SWAP(a[k][indxr[l]], a[k][indxc[l]])
	}
	
	delete[] indxc; indxc = NULL;
	delete[] indxr; indxr = NULL;
	delete[] ipiv; ipiv = NULL;
	return true;
}

// LOW UP DECOMPOSE
// Given matrix a (n x n), this routine replaces it by the LU decomposition of a rowwise
// permutation of itself. a is input. a is output, arranged for LU decomposed way (diagnoal elements belongs to the UP matrix);
// indx[0 .. n-1] is an output vector that records the row permutation effected by the partital
// pivoting; d is output as +/-1 depending on whether the number of row interchanges was even
// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
// or invert a matrix
template <int N> bool ludcmp(DMATRIX<N> &a, int *indx, int& d, double *vv)
{
	int n = a.n;
	int i = 0, imax = 0, j = 0, k = 0;
	double big, dum, sum, temp;
	//double *vv = new double[n];

	d = 1;
	for (i = 0; i < n; i++) {
		big = 0;
		for (j = 0; j < n; j++) {temp = fabs(a.m[i][j]); if (temp > big) big = temp;}
		if (big < 1.0e-40) {
			sprintf(errmsg, "Singular matrix in routine ludcmp"); 
			return false;
		}

		vv[i] = 1.f / big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a.m[i][j];
			for (k = 0; k < i; k++) sum -= a.m[i][k] * a.m[k][j];
			a.m[i][j] = sum;
		}
		big = 0;
		for (i = j; i < n; i++) {
			sum = a.m[i][j];
			for (k = 0; k < j; k++) sum -= a.m[i][k] * a.m[k][j];
			a.m[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = a.m[imax][k];
				a.m[imax][k] = a.m[j][k];
				a.m[j][k] = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a.m[j][j] == 0) a.m[j][j] = 1.0e-30; // to avoid the singular
		if (j != n) {
			dum = 1.0 / a.m[j][j];
			for (i = j + 1; i < n; i++) a.m[i][j] *= dum;
		}
	}
	//delete[] vv;
	return true;
}

// BACK SUBSITUTING
// Solves the set of n linear equation A * X = B. Here A (n x n) is input, not original matrix,
// but rather as its LU decomposition, determined by the routine ludcmp. indx[0 .. n-1] is input
// as the permutation vector returned by ludcmp. b[0 .. n-1] is input as the right-hand side vector
// B, and returns with the solution vector X. a, indx are not modified by this routine and can be
// left in place for successive calls with different right-hand sides b. This routine takes
// into account the possibility that b will begin with many zero elements, so it is efficient for
// use in matrix inversion.
template <int N> bool lubksb(DMATRIX<N> &a, int *indx, double *b)
{
	int n = a.n;
	int i, ii = -1, ip, j;
	double sum = 0;

	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0) for (j = ii; j <= i - 1; j++) sum -= a.m[i][j] * b[j];
		else if (sum != 0) ii = i;
		b[i] = sum;
	}
	for (i = n - 1; i >= 0; i--) {
		sum = b[i];
		for (j = i + 1; j < n; j++) sum -= a.m[i][j] * b[j];
		b[i] = sum / a.m[i][i];
	}
	return true;
}

// Solve linear set of equation A * X = B
// input a -- MATRIX A, b -- array B; X will be output to b
template <int N> bool SolveEquation(DMATRIX<N> &a, double *b)
{
	int *indx = new int[a.n];
	double *vv = new double[a.n];
	int d = 0;
	if (!ludcmp<N>(a, indx, d, vv)) {delete[] indx; delete[] vv; return false;}
	if (!lubksb<N>(a, indx, b)) {delete[] indx; delete[] vv; return false;}
	delete[] indx; delete[] vv;
	return true;
}

// Inverse matrix a
// input a; output y -- inverse matrix of a
// Inverse matrix a
// input a; output y -- inverse matrix of a
template <int N> bool DInverseMatrix(DMATRIX<N> &a, DMATRIX<N> &y, int *indx, double *vv, double *c)
{
	int n = a.n, d = 0;
	int i, j;

	// normalize the coefficient to maximum value, modified on 06/02/2014
	double vmax = 0, v;
	for (i = 0; i < N; i++) {
		if (i == 0) vmax = FABS(a.m[i][i]);
		else {
			v = FABS(a.m[i][i]);
			vmax = (vmax < v ? v : vmax);
		}
	}
	if (vmax == 0) return false;
	v = 1.0 / vmax;
	for (i = 0; i < N; i++) for (j = 0; j < N; j++) a.m[i][j] *= v;
	
	//int *indx = new int[n];
	//T *c = new T[n];
	//T vv = new T[n];

	if (!ludcmp<N>(a, indx, d, vv)) {
		//delete[] indx; delete[] c; delete[] vv;
		return false;
	}
	for (j = 0; j < n; j++) {
		//for (i = 0; i < n; i++) c[i] = 0;
		memset(c, 0, n * SIZE_DOUBLE);
		c[j] = 1;
		lubksb<N>(a, indx, c);
		for (i = 0; i < n; i++) y.m[i][j] = c[i];
	}
	//delete[] indx; delete[] c; delete[] vv;

	// renormalized the result, modified on 06/02/2014
	for (i = 0; i < N; i++) for (j = 0; j < N; j++) y.m[i][j] *= v; // re-normalize the inverse matrix

	return true;
}

template <int N> bool InverseDMatrix(DMATRIX<N> &a, DMATRIX<N> &y) {
	int indx[N];
	double vv[N], c[N];

	bool status = DInverseMatrix<N>(a, y, indx, vv, c);

	return status;
}

template <int N> bool InverseMatrix(MATRIX<N> &a, MATRIX<N> &y) {
	int indx[N];
	double vv[N], c[N];
	DMATRIX<N> m1, m2;
	int i, j;

	//for (i = 0; i < N; i++) for (j = 0; j < N; j++) m1.m[i][j] = a.m[i][j];
	memcpy(m1.m[0], a.m[0], N * N * SIZE_DOUBLE);
	bool status = DInverseMatrix<N>(m1, m2, indx, vv, c);
	memcpy(y.m[0], m2.m[0], N * N * SIZE_DOUBLE);

	return status;
}

// Newton-Raphson method for Nonlinear Systems of Equations
template <int N, class DERIVED_CLASS> class NR_NonLinearEqs {
public:
	DMATRIX<N> J; // Jacobian matrix J
	DMATRIX<N> invJ; // inverse J
	DVECTOR<N> b; // -F(v)
	DVECTOR<N> vt; // temporary variable for the result
	DVECTOR<N> dv;

	DERIVED_CLASS *m_dclass;

	// "A Member Function Template Shall not be Virtual" -- defined in ANSI/ISO Standard
	//virtual void JacobianMatrix_b(); // construct the Jacobian Matrix J and b with given vt;
	void *func_JacobianMatrix_b;
	NR_NonLinearEqs() {m_dclass = NULL; func_JacobianMatrix_b = NULL;};
	bool NRIterate(DVECTOR<N>& v_init, DVECTOR<N> &vres, double vstep, double vtol, int Nmin, int Nmax); // calculate the equations with maximum loops Nmax, and the tolerate resolution
};

template <int N, class DERIVED_CLASS> bool NR_NonLinearEqs<N, DERIVED_CLASS>::NRIterate(DVECTOR<N> &v_init, DVECTOR<N> &vres, double vstep, double vtol, int Nmin, int Nmax) {
	bool status = true;
	int i, j, iloop = 0;
	double dvmax = vtol + vtol;
	//for (i = 0; i < N; i++) vt.v[i] = v_init.v[i];
	memcpy(vt.v, v_init.v, N * SIZE_DOUBLE);

	typedef void (*FUNC)(DERIVED_CLASS *m);

	while (iloop < Nmax && (iloop < Nmin || dvmax > vtol)) {
		((FUNC)func_JacobianMatrix_b)(m_dclass);
		status = InverseDMatrix<N>(J, invJ);
		if (!status) return false;
		//dv = invJ * b;
		MpV(invJ, b, dv, i, j, N)
		MAX(dv, i, dvmax); dvmax = FABS(dvmax);
		for (i = 0; i < N; i++) vt.v[i] = vt.v[i] + vstep * dv.v[i];
		if (dvmax <= vtol && iloop >= Nmin) {
			//for (i = 0; i < N; i++) vres.v[i] = vt.v[i];
			memcpy(vres.v, vt.v, N * SIZE_DOUBLE);
			return true;
		}
		iloop+=1;
	}
	//for (i = 0; i < N; i++) vres.v[i] = vt.v[i];
	memcpy(vres.v, vt.v, N * SIZE_DOUBLE);
	return false;
};

