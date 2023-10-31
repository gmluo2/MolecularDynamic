#define PI 3.141592653589793100
#define PI2 6.283185307179586200
#define Angle2Radian  0.017453292519943295  // (PI / 180)
#define Radian2Angle  57.295779513082323    // (180 / PI)

#define fPI 3.1415926535897931f
#define fPI2 6.2831853071795862f
#define fAngle2Radian  0.017453292519943295f  // (PI / 180)
#define fRadian2Angle  57.295779513082323f    // (180 / PI)

#define LOOP1(i, N) for (i = 0; i < N; i++) {
#define END_LOOP1 }
#define LOOP2(i, j, N) for (i = 0; i < N; i++) {for (j = 0; j < N; j++) {
#define END_LOOP2 }}

template <class T, unsigned short N> class _VECTOR_BASE {
public:
	short int n;
	T v[N];
	_VECTOR_BASE() {
		n = N;
		memset(v, 0, n * sizeof(T));
	};
	_VECTOR_BASE(T vd) {
		n = N;
		if (vd == 0) {memset(v, 0, n * sizeof(T)); return;}
		for (int i = 0; i < N; i++) v[i] = vd;
	};
	void Show(char* str = NULL) {
		char buffer[256] = "[\0", buff[60] = "\0";
		int i = 0;
		for (i = 0; i < n; i++) {
			if (i != 0) strcat(buffer, ",  ");
			sprintf(buff, "%f", v[i]);
			strcat(buffer, buff);
		}
		strcat(buffer, "]");
		if (str != NULL) strcpy(str, buffer);
		else show_msg(buffer);
	};
	void operator = (_VECTOR_BASE<T, N>& v) {
		memcpy(this->v, v.v, N * sizeof(T));
	};
	void operator += (_VECTOR_BASE<T, N>& v) {
		int i = 0;
		for (i = 0; i < n; i++) this->v[i] += v.v[i];
	};
	void operator -= (_VECTOR_BASE<T, N>& v) {
		int i = 0;
		for (i = 0; i < n; i++) this->v[i] -= v.v[i];
	};
	void operator *= (T v) {
		int i = 0;
		for (i = 0; i < n; i++) this->v[i] *= v;
	};
	void operator /= (T v) {
		int i = 0;
		for (i = 0; i < n; i++) this->v[i] /= v;
	};
	_VECTOR_BASE<T, N> operator + (_VECTOR_BASE<T, N>& v) {
		_VECTOR_BASE<T, N> r;
		int i = 0;
		for (i = 0; i < n; i++) r.v[i] = this->v[i] + v.v[i];
		return r;
	};
	_VECTOR_BASE<T, N> operator - (_VECTOR_BASE<T, N>& v) {
		_VECTOR_BASE<T, N> r;
		int i = 0;
		for (i = 0; i < n; i++) r.v[i] = this->v[i] - v.v[i];
		return r;
	};
	_VECTOR_BASE<T, N> operator * (T v) {
		_VECTOR_BASE<T, N> r;
		int i = 0;
		for (i = 0; i < n; i++) r.v[i] = this->v[i] * v;
		return r;
	};
	_VECTOR_BASE<T, N> operator / (T v) {
		_VECTOR_BASE<T, N> r;
		int i = 0;
		for (i = 0; i < n; i++) r.v[i] = this->v[i] / v;
		return r;
	};
	T operator * (_VECTOR_BASE<T, N>& u) {
		int i;
		T res = 0;
		for (i = 0; i < N; i++) res += this->v[i] * u.v[i];
		return res;
	};
	T abs() {
		int i;
		T res = 0;
		for (i = 0; i < N; i++) res += v[i] * v[i];
		return (T)sqrt(res);
	};
};

class VECTOR3 : public _VECTOR_BASE<double, 3> { //three dimention
public:
	VECTOR3() {memset(v, 0, SIZE_V3);};
	VECTOR3(double vd) {v[0] = vd; v[1] = vd; v[2] = vd;};
	VECTOR3(double *vd) {memcpy(_VECTOR_BASE<double, 3>::v, vd, SIZE_V3);};
	VECTOR3(float x, float y, float z) {
		v[0] = x; v[1] = y; v[2] = z;
	};
	void Reset() {
		memset(v, 0, SIZE_V3);
	};
	bool operator == (VECTOR3& vector) {
		//return (v[0] == vector.v[0] && v[1] == vector.v[1] && v[2] == vector.v[2]);
		return (memcmp(this->v, vector.v, SIZE_V3) == 0 ? true : false);
	};
	bool operator != (VECTOR3& vector) {
		/*
		bool status = true;
		if (v[0] == vector.v[0] && v[1] == vector.v[1] && v[2] == vector.v[2]) status = true;
		else status = false;
		return (!status);
		*/
		return (memcmp(this->v, vector.v, SIZE_V3) == 0 ? false : true);
	};
	void operator = (const VECTOR3& vector) {
		//v[0] = vector.v[0]; v[1] = vector.v[1]; v[2] = vector.v[2];
		memcpy(this->v, vector.v, SIZE_V3);
		return;
	};
	VECTOR3 operator + (VECTOR3& v) {
		VECTOR3 r;
		int i = 0;
		for (i = 0; i < 3; i++) r.v[i] = this->v[i] + v.v[i];
		return r;
	};
	VECTOR3 operator - (VECTOR3& v) {
		VECTOR3 r;
		int i = 0;
		for (i = 0; i < 3; i++) r.v[i] = this->v[i] - v.v[i];
		return r;
	};
	void operator += (VECTOR3& vector) {
		v[0] += vector.v[0];
		v[1] += vector.v[1];
		v[2] += vector.v[2];
	};
	void operator -= (VECTOR3& vector) {
		v[0] -= vector.v[0];
		v[1] -= vector.v[1];
		v[2] -= vector.v[2];
	};
	void operator *= (double x) {
		v[0] = v[0] * x;
		v[1] = v[1] * x;
		v[2] = v[2] * x;
	};
	void operator /= (double x) {
		v[0] /= x;
		v[1] /= x;
		v[2] /= x;
	};
	VECTOR3 operator * (double x) {
		VECTOR3 result;
		result.v[0] = v[0] * x;
		result.v[1] = v[1] * x;
		result.v[2] = v[2] * x;
		return result;
	};
	VECTOR3 operator / (double x) {
		VECTOR3 result;
		result.v[0] = v[0] / x;
		result.v[1] = v[1] / x;
		result.v[2] = v[2] / x;
		return result;
	};
	VECTOR3 operator ^ (VECTOR3& v1) {
		VECTOR3 result;
		result.v[0] = v[1] * v1.v[2] - v[2] * v1.v[1];
		result.v[1] = v[2] * v1.v[0] - v[0] * v1.v[2];
		result.v[2] = v[0] * v1.v[1] - v[1] * v1.v[0];
		return result;
	};
	void operator ^= (VECTOR3& v1) {
		double x = v[1] * v1.v[2] - v[2] * v1.v[1];
		double y = v[2] * v1.v[0] - v[0] * v1.v[2];
		double z = v[0] * v1.v[1] - v[1] * v1.v[0];
		v[0] = x;
		v[1] = y;
		v[2] = z;
		return;
	};
	double operator * (VECTOR3& v1) {
		double f = 0;
		f += this->v[0] * v1.v[0];
		f += this->v[1] * v1.v[1];
		f += this->v[2] * v1.v[2];
		return f;
	};
	void Abs() {
		if (v[0] < 0) v[0] = -v[0];
		if (v[1] < 0) v[1] = -v[1];
		if (v[2] < 0) v[2] = -v[2];
	};
	void stretch(double length) {
		double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		len = length / len;
		v[0] *= len; v[1] *= len; v[2] *= len;
	};
	double theta() { // in radian, [0, PI]
		double r = this->abs();
		double x = v[2] / r;
		if (FABS(x) >= 1) return (sign(x) == 1 ? 0 : PI);
		return (acos(x));
	};
	double phi() { // in radian, [0, 2 PI]
		double x = 0, phi = 0;
		if (FABS(v[0]) < 1e-8) {
			show_infor("strange: v[0] is very small for phi", true);
			if (FABS(v[1]) < 1e-8) phi = 0; 
			else if (v[1] > 0) phi = PI/2;
			else if (v[1] < 0) phi = PI*1.5; // 3/2 PI
		}
		else {
			x = v[1] / v[0];
			if (x >= 0 && v[0] > 0) phi = atan(x);
			else if (x >= 0 && v[0] < 0) phi = PI + atan(x);
			else if (x < 0 && v[0] < 0) phi = PI + atan(x);
			else phi = PI2 + atan(x); 
		}
		return phi;
	};
};

template <class T, unsigned short N> class SVECTOR { // simplified vector just for data component, no constructor and member functions
public:
	T v[N];
};

inline void cp(VECTOR3 &u, VECTOR3 &v, VECTOR3 &res) {
	res.v[0] = u.v[1] * v.v[2] - u.v[2] * v.v[1]; // uy * vz - uz * vy
	res.v[1] = u.v[2] * v.v[0] - u.v[0] * v.v[2]; // uz * vx - ux * vz
	res.v[2] = u.v[0] * v.v[1] - u.v[1] * v.v[0]; // ux * vy - uy * vx
}

template <unsigned short N> class VECTOR : public _VECTOR_BASE<double, N> {
public:
	VECTOR() {};
	VECTOR(double vd) {
		if (vd == 0) {memset(_VECTOR_BASE<double, N>::v, 0, N * SIZE_DOUBLE); return;}
		for (int i = 0; i < N; i++) this->v[i] = vd;
	};

	void operator = (const VECTOR<N> &v1) {
		//for (int i = 0; i < N; i++) this->v[i] = v1.v[i];
		memcpy(_VECTOR_BASE<double, N>::v, v1.v, N * SIZE_DOUBLE);
	};
};

template <unsigned short N> class DVECTOR : public _VECTOR_BASE<double, N> {
public:
	DVECTOR() {};
};

template <class T, unsigned short N> class _MATRIX_BASE { // N x N matrix
public:
	char n;
	T m[N][N];
	_MATRIX_BASE(int n = N) {
		this->n = N;
		memset(m[0], 0, N * N  * sizeof(T));
	};
	void SetZero() {
		/*
		int i = 0, j = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) m[i][j] = 0;
		*/
		memset(m[0], 0, N * N  * sizeof(T));
	};
	_MATRIX_BASE(T v) {
		n = N;
		memset(m[0], 0, N * N  * sizeof(T));
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			m[i][i] = v; // m[i][i] = v;
		}
	};
	_MATRIX_BASE(T *v) { // v is N dimensional
		n = N;
		memset(m[0], 0, N * N  * sizeof(T));
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			m[i][i] = v[i]; // m[k][k] = v[k];
		}
	};
	_MATRIX_BASE(T* u, T* v) { // u & v are N dimensional
		n = N;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				m[i][j] = u[i] * v[j];
			}
		}
	};
	void operator = (const _MATRIX_BASE<T, N>& matrix) {
		if (this->n != matrix.n) return;
		memcpy(m[0], matrix.m, N * N * sizeof(T));
	};
	void operator += (_MATRIX_BASE<T, N>& matrix) {
		if (this->n != matrix.n) return;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) this->m[i][j] += matrix.m[i][j];
	};
	void operator -= (_MATRIX_BASE<T, N>& matrix) {
		if (this->n != matrix.n) return;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) this->m[i][j] -= matrix.m[i][j];
	};
	void operator *= (T v) {
		int i = 0, j = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) this->m[i][j] *= v;
	};
	void operator /= (T v) {
		int i = 0, j = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) this->m[i][j] /= v;
	};
	_MATRIX_BASE<T, N> operator + (_MATRIX_BASE<T, N> &m1) {
		_MATRIX_BASE<T, N> matrix;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) matrix.m[i][j] = this->m[i][j] + m1.m[i][j];
		}
		return matrix;
	};
	_MATRIX_BASE<T, N> operator - (_MATRIX_BASE<T, N> &m1) {
		_MATRIX_BASE<T, N> matrix;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) matrix.m[i][j] = this->m[i][j] - m1.m[i][j];
		}
		return matrix;
	};
	_MATRIX_BASE<T, N> operator* (T v) {
		_MATRIX_BASE<T, N> matrix;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) matrix.m[i][j] = this->m[i][j] * v;
		}
		return matrix;
	};
	_MATRIX_BASE<T, N> operator / (T v) {
		_MATRIX_BASE<T, N> matrix;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) matrix.m[i][j] = this->m[i][j] / v;
		}
		return matrix;
	};
	void Show(char *title = NULL, char *str = NULL) {
		char buffer[2560] = "\0", buff[250] = "\0";
		int i, j;
		if (title != NULL) sprintf(buffer, "%s : \n", title);
		else sprintf(buffer, "MATRIX : \n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {sprintf(buff, "  %f", m[i][j]); strcat(buffer, buff);}
			strcat(buffer, "\n");
		}
		if (str != NULL) strcpy(str, buffer);
		else show_msg(buffer);
	};
};
/*
template <int N> class MATRIX : public _MATRIX_BASE<double, N> {
public:
	MATRIX() {};
	void operator = (double v) {
		int i = 0, j = 0;
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				this->m[i][j] = v;
			}
		}
	};
};
*/

template <unsigned short N> class DMATRIX : public _MATRIX_BASE<double, N> {
public:
	DMATRIX() {};
	DMATRIX(double v) {
		if (v == 0) {
			memset(&(this->m[0][0]), 0, N * N * SIZE_DOUBLE);
			return;
		}
		int i = 0, j = 0;
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				this->m[i][j] = v;
			}
		}
	};
};

template <unsigned short N> class MATRIX : public DMATRIX<N> {
public:
	MATRIX() {};
	MATRIX(double v) {
		int i, j;
		if (v == 0) {
			memset(&(this->m[0][0]), 0, N * N * SIZE_DOUBLE); return;
		}
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				this->m[i][j] = v;
			}
		}
	};
};

class R : public _MATRIX_BASE<double, 3> { // rotation matrix
public:
	int axis;
	double arc_angle;
	R() {
		axis = 0;
		arc_angle = 0.0;
		memset(m[0], 0, SIZE_M3);
	};
	R(int axis, double angle, bool Arc = false) {
		int i = 0, j = 0;
		this->axis = axis;
		if (Arc) {this->arc_angle = angle;}
		else {this->arc_angle = angle * Angle2Radian;}
		switch (axis) {
		case 1: // x axis
			m[0][0] = 1.0;
			m[0][1] = 0.0;
			m[0][2] = 0.0;
			m[1][0] = 0.0;
			m[1][1] = cos(arc_angle);
			m[1][2] = -sin(arc_angle);
			m[2][0] = 0.0;
			m[2][1] = -m[1][2];
			m[2][2] = m[1][1];
			break;
		case 2: // y axis
			m[0][0] = cos(arc_angle);
			m[0][1] = 0.0;
			m[0][2] = sin(arc_angle);
			m[1][0] = 0.0;
			m[1][1] = 1.0;
			m[1][2] = 0.0;
			m[2][0] = -m[0][2];
			m[2][1] = 0.0;
			m[2][2] = m[0][0];
			break;
		case 3: // z axis
			m[0][0] = cos(arc_angle);
			m[0][1] = -sin(arc_angle);
			m[0][2] = 0.0;
			m[1][0] = -m[0][1];
			m[1][1] = m[0][0];
			m[1][2] = 0.0;
			m[2][0] = 0.0;
			m[2][1] = 0.0;
			m[2][2] = 1.0;
			break;
		default: //undefined
			/*
			for (i = 0; i < 3; i++) {
				for(j = 0; j < 3; j++) {
					m[i][j] = 0.0;
				}
			}
			*/
			memset(m[0], 0, SIZE_M3);
			break;
		}
	};

	friend R Inverse(R& r2) {
		return R(r2.axis, -r2.arc_angle, true);
	};

	void operator = (R& r2) {
		memcpy(m[0], r2.m[0], SIZE_M3);
	};

	R operator * (R& r2) {
		R result;
		int i, j, k;
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 3; k++) {
					result.m[i][j] += m[i][k] * r2.m[k][j];
				}
			}
		}
		return result;
	};

	VECTOR3 operator * (VECTOR3& v) {
		VECTOR3 result;
		int i, j;
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				result.v[i] += m[i][j] * v.v[j];
			}
		}
		return result;
	};

	void rotate(VECTOR3 &v, VECTOR3 &result) {
		int i, j;
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				result.v[i] += m[i][j] * v.v[j];
			}
		}
	};
};

void rotate(VECTOR3 &axis, double omega, VECTOR3 &v, VECTOR3 &result, bool arc);
void rotate(double theta, double phi, double omega, VECTOR3 &v, VECTOR3 &result, bool arc);
double angle(VECTOR3 &u, VECTOR3 &v); //angle between u & v, in arc degree

// create the MATRIX for rotation
void RMatrix(VECTOR3 &axis, double omega, R &r, bool arc);

// Product of N x N MATRIXs 
template<unsigned short N> void product(MATRIX<N> &m1, MATRIX<N> &m2, MATRIX<N> &res) {
	int i = 0, j = 0, k = 0;
	LOOP2(i, j, N)
		res.m[i][j] = 0;
		LOOP1(k, N)
			res.m[i][j] += m1.m[i][k] * m2.m[k][j];
		END_LOOP1
	END_LOOP2
};

template <unsigned short N> class CMATRIX { // N x N x N matrix
public:
	short int n;
	double m[N][N][N];
	CMATRIX(int n = N) {
		this->n = N;
		/*
		int i = 0, j = 0, k = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++) 
			m[i][j][k] = 0;
			*/
		memset(m[0][0], 0, N * N * N * SIZE_DOUBLE);
	};
	void SetZero() {
		/*
		int i = 0, j = 0, k = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++) m[i][j][k] = 0;
		*/
		memset(m[0][0], 0, N * N * N * SIZE_DOUBLE);
	};
	void operator = (const CMATRIX<N> &matrix) {
		if (matrix.m == NULL) {SetZero(); return;}
		if (this->n != matrix.n) return;
		memcpy(m[0][0], matrix.m[0][0], N * N * N * SIZE_DOUBLE);
	};
	void Show(char *title = NULL) {
		char buffer[2560] = "\0", buff[250] = "\0";
		int i, j, k;
		if (title != NULL) sprintf(buffer, "%s : \n", title);
		else sprintf(buffer, "CUBIC MATRIX : \n");
		for (i = 0; i < n; i++) {
			sprintf(buff, "LAYER %d\n", i); strcat(buffer, buff);
			for (j = 0; j < n; j++) {
				for (k = 0; k < n; k++) {
					sprintf(buff, "  %f", m[i][j][k]); strcat(buffer, buff);
				}
			}
			strcat(buffer, "\n");
		}
		show_msg(buffer);
	};
};

template <unsigned short N> class QMATRIX { // N x N x N x N matrix
public:
	short int n;
	double m[N][N][N][N];
	QMATRIX(int n = N) {
		this->n = N;
		/*
		int i = 0, j = 0, k = 0, l = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++) for (l = 0; l < n; l++)
			m[i][j][k][l] = 0;
			*/
		memset(m[0][0][0], 0, N * N * N * N * SIZE_DOUBLE);
	};
	void SetZero() {
		/*
		int i = 0, j = 0, k = 0, l = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++)  for (l = 0; l < n; l++) 
			m[i][j][k][l] = 0;
			*/
		memset(m[0][0][0], 0, N * N * N * N * SIZE_DOUBLE);
	};
	void operator = (const QMATRIX<N>& matrix) {
		if (this->m == NULL) {SetZero(); return;}
		if (this->n != matrix.n) return;
		memcpy(m[0][0][0], matrix.m[0][0][0], N * N * N * N * SIZE_DOUBLE);
	};
	void Show(char *title = NULL) {
		char buffer[20560] = "\0", buff[250] = "\0";
		int i, j, k, l;
		if (title != NULL) sprintf(buffer, "%s : \n", title);
		else sprintf(buffer, "Quadru-MATRIX : \n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				sprintf(buff, "LAYER [%d, %d]\n", i, j); strcat(buffer, buff);
				for (k = 0; k < n; k++) {
					for (l = 0; l < n; l++) {
						sprintf(buff, "  %f", m[i][j][k][l]); strcat(buffer, buff);
					}
				}
				strcat(buffer, "\n");
			}
		}
		show_msg(buffer);
	};
};

//#define CMzero(cm, i, j, k) for (i = 0; i < cm.n; i++) for (j = 0; j < cm.n; j++) for (k = 0; k < cm.n; k++) cm.m[i][j][k] = 0;
#define CMzero(cm) memset(cm.m[0][0], 0, cm.n * cm.n * cm.n * SIZE_DOUBLE);

//#define QMzero(qm, i, j, k, l) for (i = 0; i < qm.n; i++) for (j = 0; j < qm.n; j++) for (k = 0; k < qm.n; k++) for (l = 0; l < qm.n; l++) qm.m[i][j][k][l] = 0;
#define QMzero(qm) memset(qm.m[0][0][0], 0, qm.n * qm.n * qm.n * qm.n * SIZE_DOUBLE);

// random direction
void rand_vect3(double *v);
double rand_circle_radius();
void rand_circle_vect(double *v);
double rand_sphere_shell_radius();
void project(double *u, double *axis, double *p);
