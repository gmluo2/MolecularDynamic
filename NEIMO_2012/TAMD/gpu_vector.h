#ifdef __GPU__
namespace _gpu_dstruct_ {

#ifndef PI
#define PI 3.141592653589793100
#define PI2 6.283185307179586200
#define Angle2Radian  0.017453292519943295  // (PI / 180)
#define Radian2Angle  57.295779513082323    // (180 / PI)

#define fPI 3.1415926535897931f
#define fPI2 6.2831853071795862f
#define fAngle2Radian  0.017453292519943295f  // (PI / 180)
#define fRadian2Angle  57.295779513082323f    // (180 / PI)
#endif

#ifndef LOOP1
#define LOOP1(i, N) for (i = 0; i < N; i++) {
#define END_LOOP1 }
#define LOOP2(i, j, N) for (i = 0; i < N; i++) {for (j = 0; j < N; j++) {
#define END_LOOP2 }}
#endif

#ifndef SIZE_V3
#define SIZE_V3      24
#endif

#ifndef sign
#define sign(x) (x > 1 ? 1 : -1)
#endif

template <class T, unsigned short N> class _VECTOR_BASE {
public:
	short int n;
	T v[N];
	__device__ _VECTOR_BASE() {
		n = N;
		//memset(v, 0, n * sizeof(T));
		for (int i = 0; i < N; i++) v[i] = 0;
	};
	__device__ _VECTOR_BASE(T vd) {
		n = N;
		for (int i = 0; i < N; i++) v[i] = vd;
	};
	__device__ void operator = (_VECTOR_BASE<T, N>& v) const {
		for (int i = 0; i < N; i++) this->v[i] = v.v[i];
		//memcpy(this->v, v.v, N * sizeof(T));
		//cuMemcpy((CUdeviceptr)(this->v), (CUdeviceptr)(v.v), N * sizeof(T));
	};
	__device__ void operator += (_VECTOR_BASE<T, N>& v) const {
		int i = 0;
		for (i = 0; i < n; i++) this->v[i] += v.v[i];
	};
	__device__ void operator -= (_VECTOR_BASE<T, N>& v) const {
		int i = 0;
		for (i = 0; i < n; i++) this->v[i] -= v.v[i];
	};
	__device__ void operator *= (T v) const {
		int i = 0;
		for (i = 0; i < n; i++) this->v[i] *= v;
	};
	__device__ void operator /= (T v) const {
		int i = 0;
		for (i = 0; i < n; i++) this->v[i] /= v;
	};
	__device__ _VECTOR_BASE<T, N> operator + (_VECTOR_BASE<T, N>& v) const {
		_VECTOR_BASE<T, N> r;
		int i = 0;
		for (i = 0; i < n; i++) r.v[i] = this->v[i] + v.v[i];
		return r;
	};
	__device__ _VECTOR_BASE<T, N> operator - (_VECTOR_BASE<T, N>& v) const {
		_VECTOR_BASE<T, N> r;
		int i = 0;
		for (i = 0; i < n; i++) r.v[i] = this->v[i] - v.v[i];
		return r;
	};
	__device__ _VECTOR_BASE<T, N> operator * (T v) const {
		_VECTOR_BASE<T, N> r;
		int i = 0;
		for (i = 0; i < n; i++) r.v[i] = this->v[i] * v;
		return r;
	};
	__device__ _VECTOR_BASE<T, N> operator / (T v) const {
		_VECTOR_BASE<T, N> r;
		int i = 0;
		for (i = 0; i < n; i++) r.v[i] = this->v[i] / v;
		return r;
	};
	__device__ T operator * (_VECTOR_BASE<T, N>& u) const {
		int i;
		T res = 0;
		for (i = 0; i < N; i++) res += this->v[i] * u.v[i];
		return res;
	};
	__device__ T abs() {
		int i;
		T res = 0;
		for (i = 0; i < N; i++) res += v[i] * v[i];
		return (T)sqrt(res);
	};
};

template <class T> class VECTOR3: public _VECTOR_BASE<T, 3> { //three dimention
private:
	//int sizeV3;
public:
	__device__ VECTOR3() {
		//sizeV3 = sizeof(T) * 3;
		//memset(_VECTOR_BASE<T, 3>::v, 0, sizeV3);
		_VECTOR_BASE<T, 3>::v[0] = 0; _VECTOR_BASE<T, 3>::v[1] = 0; _VECTOR_BASE<T, 3>::v[2] = 0;
	};
	__device__ VECTOR3(T vd) {
		//sizeV3 = sizeof(T) * 3;
		_VECTOR_BASE<T, 3>::v[0] = vd; _VECTOR_BASE<T, 3>::v[1] = vd; _VECTOR_BASE<T, 3>::v[2] = vd;
	};
	__device__ VECTOR3(T *vd) {
		//sizeV3 = sizeof(T) * 3;
		//memcpy(_VECTOR_BASE<T, 3>::v, vd, sizeV3);
		_VECTOR_BASE<T, 3>::v[0] = vd[0];
		_VECTOR_BASE<T, 3>::v[1] = vd[1];
		_VECTOR_BASE<T, 3>::v[2] = vd[2];
	};
	__device__ VECTOR3(float x, float y, float z) {
		//sizeV3 = sizeof(T) * 3;
		_VECTOR_BASE<T, 3>::v[0] = x; _VECTOR_BASE<T, 3>::v[1] = y; _VECTOR_BASE<T, 3>::v[2] = z;
	};
	__device__ void Reset() {
		//memset(_VECTOR_BASE<T, 3>::v, 0, sizeV3);
		_VECTOR_BASE<T, 3>::v[0] = 0; _VECTOR_BASE<T, 3>::v[1] = 0; _VECTOR_BASE<T, 3>::v[2] = 0;
	};
	__device__ bool operator == (VECTOR3<T>& vector) {
		//return (v[0] == vector.v[0] && v[1] == vector.v[1] && v[2] == vector.v[2]);
		return ((this->v[0] == vector.v[0] && this->v[1] == vector.v[1] && this->v[2] == vector.v[2]) ? true : false);
	};
	__device__ bool operator != (VECTOR3<T>& vector) {
		return ((this->v[0] != vector.v[0] || this->v[1] != vector.v[1] || this->v[2] != vector.v[2]) ? true : false);
	};
	__device__ void operator = (const VECTOR3<T>& vector) {
		this->v[0] = vector.v[0]; this->v[1] = vector.v[1]; this->v[2] = vector.v[2];
		return;
	};
	__device__ VECTOR3<T> operator + (VECTOR3<T>& v) {
		VECTOR3<T> r(this->v[0] + v.v[0], this->v[1] + v.v[1], this->v[2] + v.v[2]);
		return r;
	};
	__device__ VECTOR3<T> operator - (VECTOR3<T>& v) {
		VECTOR3<T> r(this->v[0] - v.v[0], this->v[1] - v.v[1], this->v[2] - v.v[2]);
		return r;
	};
	__device__ void operator += (VECTOR3<T>& vector) {
		_VECTOR_BASE<T, 3>::v[0] += vector.v[0];
		_VECTOR_BASE<T, 3>::v[1] += vector.v[1];
		_VECTOR_BASE<T, 3>::v[2] += vector.v[2];
	};
	__device__ void operator -= (VECTOR3<T>& vector) {
		_VECTOR_BASE<T, 3>::v[0] -= vector.v[0];
		_VECTOR_BASE<T, 3>::v[1] -= vector.v[1];
		_VECTOR_BASE<T, 3>::v[2] -= vector.v[2];
	};
	__device__ void operator *= (T x) {
		_VECTOR_BASE<T, 3>::v[0] = _VECTOR_BASE<T, 3>::v[0] * x;
		_VECTOR_BASE<T, 3>::v[1] = _VECTOR_BASE<T, 3>::v[1] * x;
		_VECTOR_BASE<T, 3>::v[2] = _VECTOR_BASE<T, 3>::v[2] * x;
	};
	__device__ void operator /= (T x) {
		_VECTOR_BASE<T, 3>::v[0] /= x;
		_VECTOR_BASE<T, 3>::v[1] /= x;
		_VECTOR_BASE<T, 3>::v[2] /= x;
	};
	__device__ VECTOR3<T> operator * (T x) {
		VECTOR3<T> result(_VECTOR_BASE<T, 3>::v[0] * x, _VECTOR_BASE<T, 3>::v[1] * x, _VECTOR_BASE<T, 3>::v[2] * x);
		return result;
	};
	__device__ VECTOR3<T> operator / (T x) {
		VECTOR3<T> result(_VECTOR_BASE<T, 3>::v[0] / x, _VECTOR_BASE<T, 3>::v[1] / x, _VECTOR_BASE<T, 3>::v[2] / x);
		return result;
	};
	__device__ VECTOR3<T> operator ^ (VECTOR3<T>& v1) {
		VECTOR3<T> result;
		result.v[0] = _VECTOR_BASE<T, 3>::v[1] * v1.v[2] - _VECTOR_BASE<T, 3>::v[2] * v1.v[1];
		result.v[1] = _VECTOR_BASE<T, 3>::v[2] * v1.v[0] - _VECTOR_BASE<T, 3>::v[0] * v1.v[2];
		result.v[2] = _VECTOR_BASE<T, 3>::v[0] * v1.v[1] - _VECTOR_BASE<T, 3>::v[1] * v1.v[0];
		return result;
	};
	__device__ void operator ^= (VECTOR3<T>& v1) {
		double x = _VECTOR_BASE<T, 3>::v[1] * v1.v[2] - _VECTOR_BASE<T, 3>::v[2] * v1.v[1];
		double y = _VECTOR_BASE<T, 3>::v[2] * v1.v[0] - _VECTOR_BASE<T, 3>::v[0] * v1.v[2];
		double z = _VECTOR_BASE<T, 3>::v[0] * v1.v[1] - _VECTOR_BASE<T, 3>::v[1] * v1.v[0];
		_VECTOR_BASE<T, 3>::v[0] = x;
		_VECTOR_BASE<T, 3>::v[1] = y;
		_VECTOR_BASE<T, 3>::v[2] = z;
		return;
	};
	__device__ double operator * (VECTOR3<T>& v1) {
		T f = 0;
		f += this->v[0] * v1.v[0];
		f += this->v[1] * v1.v[1];
		f += this->v[2] * v1.v[2];
		return f;
	};
	__device__ void Abs() {
		if (this->v[0] < 0) this->v[0] = -this->v[0];
		if (this->v[1] < 0) this->v[1] = -this->v[1];
		if (this->v[2] < 0) this->v[2] = -this->v[2];
	};
	__device__ void stretch(T length) {
		T len = sqrt(this->v[0] * this->v[0] + this->v[1] * this->v[1] + this->v[2] * this->v[2]);
		len = length / len;
		this->v[0] *= len; this->v[1] *= len; this->v[2] *= len;
	};
	__device__ T theta() { // in radian
		T r = this->abs();
		T x = _VECTOR_BASE<T, 3>::v[2] / r;
		if (fabs(x) >= 1) return (sign(x) == 1 ? 0 : PI);
		return (acos(x));
	};
	__device__ T phi() { // in radian
		T x = 0, phi = 0;
		if (fabs(this->v[0]) < 0.00002f) {
			if (fabs(this->v[1]) < 0.00002f) phi = 0;
			else if (this->v[1] > 0) phi = PI/2;
			else if (this->v[1] < 0) phi = PI*1.5; // 3/2 PI
		}
		else {
			x = this->v[1] / this->v[0];
			if (x >= 0 && this->v[0] > 0) phi = atan(x);
			else if (x >= 0 && this->v[0] < 0) phi = PI + atan(x);
			else if (x < 0 && this->v[0] < 0) phi = PI + atan(x);
			else phi = PI2 + atan(x); 
		}
		return phi;
	};
};

template <class T, unsigned short N> class SVECTOR { // simplified vector just for data component, no constructor and member functions
public:
	T v[N];
};

template <class T> __device__ void cp(VECTOR3<T> &u, VECTOR3<T> &v, VECTOR3<T> &res) {
	res.v[0] = u.v[1] * v.v[2] - u.v[2] * v.v[1]; // uy * vz - uz * vy
	res.v[1] = u.v[2] * v.v[0] - u.v[0] * v.v[2]; // uz * vx - ux * vz
	res.v[2] = u.v[0] * v.v[1] - u.v[1] * v.v[0]; // ux * vy - uy * vx
}


template <class T, unsigned short N> class VECTOR : public _VECTOR_BASE<T, N> {
private:
	int sizeV;
public:
	__device__ VECTOR() {sizeV = sizeof(T);};
	__device__ VECTOR(T vd) {
		for (int i = 0; i < N; i++) this->v[i] = vd;
	};

	__device__ void operator = (const VECTOR<T, N> &v1) {
		for (int i = 0; i < N; i++) this->v[i] = v1.v[i];
	};
};

template <unsigned short N> class DVECTOR : public _VECTOR_BASE<double, N> {
public:
	__device__ DVECTOR() {};
};

template <class T, unsigned short N> class _MATRIX_BASE { // N x N matrix
public:
	char n;
	T m[N][N];
	__device__ _MATRIX_BASE() {
		this->n = N;
		//memset(m[0], 0, N * N  * sizeof(T));
		int i, j;
		for (i = 0; i < N; i++) {for (j = 0; j < N; j++) m[i][j] = 0;}
	};
	__device__ void SetZero() {
		int i, j;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) m[i][j] = 0;

		//memset(m[0], 0, N * N  * sizeof(T));
	};
	__device__ _MATRIX_BASE(T v) {
		n = N;
		//memset(m[0], 0, N * N  * sizeof(T));
		int i, j;
		for (i = 0; i < n; i++) { for (j = 0; j < N; j++) m[i][i] = 0;}
		for (i = 0; i < n; i++) m[i][i] = v; // m[i][i] = v;
	};
	__device__ _MATRIX_BASE(T *v) { // v is N dimensional
		n = N;
		//memset(m[0], 0, N * N  * sizeof(T));
		int i, j;
		for (i = 0; i < n; i++) { for (j = 0; j < N; j++) m[i][i] = 0;}
		for (i = 0; i < n; i++) m[i][i] = v[i]; // m[k][k] = v[k];
	};
	__device__ _MATRIX_BASE(T* u, T* v) { // u & v are N dimensional
		n = N;
		int i, j;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				m[i][j] = u[i] * v[j];
			}
		}
	};
	__device__ void operator = (const _MATRIX_BASE<T, N>& matrix) {
		if (this->n != matrix.n) return;
		//memcpy(m[0], matrix.m, N * N * sizeof(T));
		int i, j;
		for (i = 0; i < n; i++) { for (j = 0; j < N; j++) m[i][i] = matrix.m[i][j];}
	};
	__device__ void operator += (_MATRIX_BASE<T, N>& matrix) {
		if (this->n != matrix.n) return;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) this->m[i][j] += matrix.m[i][j];
	};
	__device__ void operator -= (_MATRIX_BASE<T, N>& matrix) {
		if (this->n != matrix.n) return;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) this->m[i][j] -= matrix.m[i][j];
	};
	__device__ void operator *= (T v) {
		int i = 0, j = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) this->m[i][j] *= v;
	};
	__device__ void operator /= (T v) {
		int i = 0, j = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) this->m[i][j] /= v;
	};
	__device__ _MATRIX_BASE<T, N> operator + (_MATRIX_BASE<T, N> &m1) {
		_MATRIX_BASE<T, N> matrix;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) matrix.m[i][j] = this->m[i][j] + m1.m[i][j];
		}
		return matrix;
	};
	__device__ _MATRIX_BASE<T, N> operator - (_MATRIX_BASE<T, N> &m1) {
		_MATRIX_BASE<T, N> matrix;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) matrix.m[i][j] = this->m[i][j] - m1.m[i][j];
		}
		return matrix;
	};
	__device__ _MATRIX_BASE<T, N> operator* (T v) {
		_MATRIX_BASE<T, N> matrix;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) matrix.m[i][j] = this->m[i][j] * v;
		}
		return matrix;
	};
	__device__ _MATRIX_BASE<T, N> operator / (T v) {
		_MATRIX_BASE<T, N> matrix;
		int i = 0, j = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) matrix.m[i][j] = this->m[i][j] / v;
		}
		return matrix;
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
	__device__ DMATRIX() {};
	__device__ DMATRIX(double v) {
		int i, j;
		for (i = 0; i < N; i++) {for (j = 0; j < N; j++) this->m[i][j] = v;}
	};
};

template <class T, unsigned short N> class MATRIX : public _MATRIX_BASE<T, N> {
public:
	__device__ MATRIX() {};
	__device__ MATRIX(T v) {
		int i, j;
		for (i = 0; i < N; i++) {for (j = 0; j < N; j++) this->m[i][j] = v;}
	};
};

class R : public _MATRIX_BASE<double, 3> { // rotation matrix
public:
	int axis;
	double arc_angle;
	__device__ R() {
		axis = 0;
		arc_angle = 0.0;
		//memset(m[0], 0, SIZE_M3);
		m[0][0] = 0; m[0][1] = 0; m[0][2] = 0;
		m[1][0] = 0; m[1][1] = 0; m[1][2] = 0;
		m[2][0] = 0; m[2][1] = 0; m[2][2] = 0;
	};
	__device__ R(int axis, double angle, bool Arc = false) {
		//int i = 0, j = 0;
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
			m[0][0] = 0; m[0][1] = 0; m[0][2] = 0;
			m[1][0] = 0; m[1][1] = 0; m[1][2] = 0;
			m[2][0] = 0; m[2][1] = 0; m[2][2] = 0;
			//memset(m[0], 0, SIZE_M3);
			break;
		}
	};

	__device__ friend R Inverse(R& r2) {
		return R(r2.axis, -r2.arc_angle, true);
	};

	__device__ void operator = (R& r2) {
		//memcpy(m[0], r2.m[0], SIZE_M3);
		m[0][0] = r2.m[0][0]; m[0][1] = r2.m[0][1]; m[0][2] = r2.m[0][2];
		m[1][0] = r2.m[1][0]; m[1][1] = r2.m[1][1]; m[1][2] = r2.m[1][2];
		m[2][0] = r2.m[2][0]; m[2][1] = r2.m[2][1]; m[2][2] = r2.m[2][2];
	};

	__device__ R operator * (R& r2) {
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
};

template <class T> __device__ VECTOR3<T> rotate(R &r, VECTOR3<T>& v) {
	VECTOR3<T> result;
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			result.v[i] += r.m[i][j] * v.v[j];
		}
	}
	return result;
};

template <class T> __device__ void rotate(R &r, VECTOR3<T> &v, VECTOR3<T> &result) {
	int i, j;
	for (i = 0; i < 3; i++) {
		result.v[i] = 0;
		for (j = 0; j < 3; j++) {
			result.v[i] += r.m[i][j] * v.v[j];
		}
	}
};


template <class T> __device__ void rotate(VECTOR3<T> &axis, T omega, VECTOR3<T> &v, VECTOR3<T> &result, bool arc)
{
	T theta, phi;
	theta = axis.theta();
	phi = axis.phi();
	result = rotate<T>(R(3, -phi, true), v);
	result = rotate<T>(R(2, -theta, true), result);
	result = rotate<T>(R(3, omega, arc), result);
	result = rotate<T>(R(2, theta, true), result);
	result = rotate<T>(R(3, phi, true), result);
}

template <class T> __device__ void rotate(T theta, T phi, T omega, VECTOR3<T> &v, VECTOR3<T> &result, bool arc)
{
	result = rotate<T>(R(3, -phi, arc), v);
	result = rotate<T>(R(2, -theta, arc), result);
	result = rotate<T>(R(3, omega, arc), result);
	result = rotate<T>(R(2, theta, arc), result);
	result = rotate<T>(R(3, phi, arc), result);
}

//angle between u & v, in arc degree
template <class T> __device__ T angle(VECTOR3<T> &u, VECTOR3<T> &v) {
	double uv = (u.v[0] * u.v[0] + u.v[1] * u.v[1] + u.v[2] * u.v[2]) *
		(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]);
	if (uv == 0) return 0;
	uv = sqrt(uv);
	double c = (u.v[0] * v.v[0] + u.v[1] * v.v[1] + u.v[2] * v.v[2]);
	c /= uv;
	if (c >= 1.0) return 0; //sometime, this could happen accidently, slightly > 1
	else if (c <= -1.0) return PI; //sometime, this could happen accidently, slightly < -1
	else return (T)acos(c);
}

template <class T> __device__ void RotMatrix(int axis, T angle, bool arc, MATRIX<T, 3> &r, MATRIX<T, 3> *inv_r) {
	int i, j;
	double arc_angle = angle;
	if (!arc) arc_angle = angle * Angle2Radian;
	switch (axis) {
	case 1: // x axis
			r.m[0][0] = 1.0;
			r.m[0][1] = 0.0;
			r.m[0][2] = 0.0;
			r.m[1][0] = 0.0;
			r.m[1][1] = cos(arc_angle);
			r.m[1][2] = -sin(arc_angle);
			r.m[2][0] = 0.0;
			r.m[2][1] = -r.m[1][2];
			r.m[2][2] = r.m[1][1];
			if (inv_r != NULL) {
				inv_r->m[0][0] = 1.0;
				inv_r->m[0][1] = 0.0;
				inv_r->m[0][2] = 0.0;
				inv_r->m[1][0] = 0.0;
				inv_r->m[1][1] = r.m[1][0]; // cos(-angle) = cos(angle)
				inv_r->m[1][2] = -r.m[1][2]; // sin(-angle) = -sin(angle)
				inv_r->m[2][0] = 0.0;
				inv_r->m[2][1] = -inv_r->m[1][2];
				inv_r->m[2][2] = inv_r->m[1][1];
			}
			break;
	case 2: // y axis
			r.m[0][0] = cos(arc_angle);
			r.m[0][1] = 0.0;
			r.m[0][2] = sin(arc_angle);
			r.m[1][0] = 0.0;
			r.m[1][1] = 1.0;
			r.m[1][2] = 0.0;
			r.m[2][0] = -r.m[0][2];
			r.m[2][1] = 0.0;
			r.m[2][2] = r.m[0][0];
			if (inv_r != NULL) {
				inv_r->m[0][0] = r.m[0][0]; // cos(-angle)
				inv_r->m[0][1] = 0.0;
				inv_r->m[0][2] = -r.m[0][2]; //sin(-angle);
				inv_r->m[1][0] = 0.0;
				inv_r->m[1][1] = 1.0;
				inv_r->m[1][2] = 0.0;
				inv_r->m[2][0] = -inv_r->m[0][2];
				inv_r->m[2][1] = 0.0;
				inv_r->m[2][2] = inv_r->m[0][0];
			}
			break;
	case 3: // z axis
			r.m[0][0] = cos(arc_angle);
			r.m[0][1] = -sin(arc_angle);
			r.m[0][2] = 0.0;
			r.m[1][0] = -r.m[0][1];
			r.m[1][1] = r.m[0][0];
			r.m[1][2] = 0.0;
			r.m[2][0] = 0.0;
			r.m[2][1] = 0.0;
			r.m[2][2] = 1.0;
			if (inv_r != NULL) {
				inv_r->m[0][0] = r.m[0][0]; //cos(-arc_angle);
				inv_r->m[0][1] = -r.m[0][1]; //sin(-arc_angle);
				inv_r->m[0][2] = 0.0;
				inv_r->m[1][0] = -inv_r->m[0][1];
				inv_r->m[1][1] = inv_r->m[0][0];
				inv_r->m[1][2] = 0.0;
				inv_r->m[2][0] = 0.0;
				inv_r->m[2][1] = 0.0;
				inv_r->m[2][2] = 1.0;
			}
			break;
	default: //undefined
			for (i = 0; i < 3; i++) {
				for(j = 0; j < 3; j++) {
					r.m[i][j] = 0.0;
				}
			}
			break;
	}
}

// create the MATRIX for rotation
template <class T> __device__ void RMatrix(VECTOR3<T> &axis, T omega, R &r, bool arc)
{
	int i, j, k;
	if (FABS(omega) < 1e-8) {I3Matrix(r)  return;}
	if (axis.v[0] == 0 && axis.v[1] == 0 && axis.v[2] == 0) {I3Matrix(r) return;}
	double theta, phi;
	theta = axis.theta();
	phi = axis.phi();
	/*
	result = R(3, -phi) * v;
	result = R(2, -theta) * result;
	result = R(3, omega) * result;
	result = R(2, theta) * result;
	result = R(3, phi) * result;
	*/
	/*
	r = R(3, phi, true) * R(2, theta, true) * R(3, omega, arc) * R(2, -theta, true) * R(3, -phi, true);
	*/


	MATRIX<T, 3> R3_phi, R2_theta, R3_omega, inv_R2_theta, inv_R3_phi, m1, m2, m3;
	RotMatrix<T>(3, phi, true, R3_phi, &inv_R3_phi);
	RotMatrix<T>(2, theta, true, R2_theta, &inv_R2_theta);
	RotMatrix<T>(3, omega, arc, R3_omega, NULL);

	MpM(R3_phi, R2_theta, m1, i, j, k, 3)
	MpM(m1, R3_omega, m2, i, j, k, 3)
	MpM(inv_R2_theta, inv_R3_phi, m3, i, j, k, 3)
	MpM(m2, m3, r, i, j, k, 3)

}

// Product of N x N MATRIXs 
template<class T, unsigned short N> __device__ void product(MATRIX<T, N> &m1, MATRIX<T, N> &m2, MATRIX<T, N> &res) {
	int i = 0, j = 0, k = 0;
	LOOP2(i, j, N)
		res.m[i][j] = 0;
		LOOP1(k, N)
			res.m[i][j] += m1.m[i][k] * m2.m[k][j];
		END_LOOP1
	END_LOOP2
};

template <class T, unsigned short N> class CMATRIX { // N x N x N matrix
public:
	short int n;
	T m[N][N][N];
	__device__ CMATRIX(int n = N) {
		this->n = N;

		int i = 0, j = 0, k = 0;
		for (i = 0; i < n; i++) {for (j = 0; j < n; j++) {for (k = 0; k < n; k++) m[i][j][k] = 0;}}

		//memset(m[0][0], 0, N * N * N * sizeof(T));
	};
	__device__ void SetZero() {

		int i = 0, j = 0, k = 0;
		for (i = 0; i < n; i++) {for (j = 0; j < n; j++) {for (k = 0; k < n; k++) m[i][j][k] = 0;}}

		//memset(m[0][0], 0, N * N * N * SIZE_DOUBLE);
	};
	__device__ void operator = (const CMATRIX<T, N> &matrix) {
		if (matrix.m == NULL) {SetZero(); return;}
		if (this->n != matrix.n) return;
		int i = 0, j = 0, k = 0;
		for (i = 0; i < n; i++) {for (j = 0; j < n; j++) {for (k = 0; k < n; k++) m[i][j][k] = matrix.m[i][j][k];}}
		//memcpy(m[0][0], matrix.m[0][0], N * N * N * sizeof(T));
	};
};

template <class T, unsigned short N> class QMATRIX { // N x N x N x N matrix
public:
	short int n;
	T m[N][N][N][N];
	__device__ QMATRIX(int n = N) {
		this->n = N;
		/*
		int i = 0, j = 0, k = 0, l = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++) for (l = 0; l < n; l++)
			m[i][j][k][l] = 0;
			*/
		memset(m[0][0][0], 0, N * N * N * N * SIZE_DOUBLE);
	};
	__device__ void SetZero() {
		/*
		int i = 0, j = 0, k = 0, l = 0;
		for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++)  for (l = 0; l < n; l++) 
			m[i][j][k][l] = 0;
			*/
		memset(m[0][0][0], 0, N * N * N * N * SIZE_DOUBLE);
	};
	__device__ void operator = (const QMATRIX<T, N>& matrix) {
		if (this->m == NULL) {SetZero(); return;}
		if (this->n != matrix.n) return;
		memcpy(m[0][0][0], matrix.m[0][0][0], N * N * N * N * SIZE_DOUBLE);
	};
};

//#define CMzero(cm, i, j, k) for (i = 0; i < cm.n; i++) for (j = 0; j < cm.n; j++) for (k = 0; k < cm.n; k++) cm.m[i][j][k] = 0;
//#define CMzero(cm) memset(cm.m[0][0], 0, cm.n * cm.n * cm.n * SIZE_DOUBLE);

//#define QMzero(qm, i, j, k, l) for (i = 0; i < qm.n; i++) for (j = 0; j < qm.n; j++) for (k = 0; k < qm.n; k++) for (l = 0; l < qm.n; l++) qm.m[i][j][k][l] = 0;
//#define QMzero(qm) memset(qm.m[0][0][0], 0, qm.n * qm.n * qm.n * qm.n * SIZE_DOUBLE);

// random vector
//void rand_vect3(double *v);

template <class T> class ARRAY {
public:
	bool bExt; // external controlled array, e.g. the memory is allocated from host
	int n;
	T *m;

	__host__ void external(bool bExt) {this->bExt = bExt;};
	__host__ void release_host() {
		if (!bExt) return;
		if (m != NULL) {cudaFree(m); m = NULL;}
		n = 0;
	};
	__host__ bool set_array_host(int n) {
		bExt = true;
		if (cudaMalloc(&(m), n * sizeof(T)) != cudaSuccess) return false; this->n = n;
		cuMemsetD8((CUdeviceptr)m, 0, n * sizeof(T));
		return true;
	};

	__device__ ARRAY(bool bExt = false) {n = 0; m = NULL; this->bExt = bExt;};
	__device__ void release() {
		if (bExt) return;
		if (m != NULL) {
			delete[] m; 
			m = NULL;
		}
		n = 0;
	};
	__host__ __device__ void _init() {
		if (bExt) return;
		m = NULL; n = 0;
	};

	__device__ bool SetArray(int n) {
		if (bExt) return true;
		if (this->n == n) return true;
		release(); if (n <= 0) return true;
		m = new T[n]; 
		//m = (T*)malloc(n * sizeof(T));
		this->n = n; return true;
	};

	__device__ bool set_array(int n) {return SetArray(n);};
	__device__ bool exist(T v) {
		for (int i = 0; i < n; i++) {
			if (m[i] == v) return true;
		}
		return false;
	};
	__device__ int indx(T v) {
		for (int i = 0; i < n; i++) {
			if (m[i] == v) return i;
		}
		return -1;
	};
	__device__ void release_pointer() {
		m = NULL; n = 0;
	};
	/*
	__device__ void expand(int nIncrease) {
		int n1 = n;
		int n2 = n + nIncrease;
		T *bf = m;
		release_pointer();
		set_array(n2);
		if (bf == NULL) return;
		//for (i = 0; i < n1; i++) m[i] = bf[i];
		memcpy(m, bf, n1 * sizeof(T));
		delete[] bf;
		//free(bf);
	};
	__device__ void remove_one(int indx) {
		if (n == 0 || indx < 0) indx += n;
		if (indx < 0 || indx >= n) return;
		T* bf = m;
		int n2 = n - 1;
		if (n2 <= 0) {release(); return;}
		else {
			release_pointer();
			set_array(n2);
			if (indx > 0) memcpy(m, bf, indx * sizeof(T));
			if (indx < n) memcpy(m + indx, bf + indx + 1, (n - indx - 1) * sizeof(T));
			delete bf; bf = NULL;
			return;
		}
	};
	*/
	__device__ ~ARRAY() {release();};
};
/*
// T is the type without pointer inside
template <class T> class ARRAY_FLEX {
public:
	int n, nm; // n is the dimension used, nm is the memory size
	T *m;

	__device__ void _init() {
		m = NULL; n = 0; m = 0;
	};
	__device__ ARRAY_FLEX() {
		n = 0; m = NULL; nm = 0;
	};
	__device__ void release(bool release_memory = false) {
		if (release_memory) {
			if (m != NULL) {delete[] m; m = NULL;} n = 0; nm = 0;
		}
		else n = 0;
	};

	__device__ void SetArray(int nm) {
		if (this->nm == nm) return;
		release(true); if (nm <= 0) return;
		m = new T[nm]; this->nm = nm; n = 0;
	};

	__device__ void set_array(int nm) {SetArray(nm);};
	__device__ void release_pointer() {
		m = NULL; nm = 0; n = 0;
	};

	__device__ void expand(int nIncrease) {
		int n1 = nm;
		int n2 = nm + nIncrease;
		T *bf = m;
		release_pointer();
		set_array(n2);
		if (bf == NULL) return;
		memcpy(m, bf, n1 * sizeof(T));
		delete[] bf;
	};
	__device__ void attach(T v) {
		if (n == nm) expand(5);
		memcpy(m + n, &v, sizeof(T)); 
		n++;
	};
	__device__ ~ARRAY_FLEX() {release(true);};
};


template <class T> class FLEX_ARRAY {
public:
	ARRAY<T*> a;

	__device__ T* get(int indx) {
		if (indx < 0) indx += a.n;
		if (indx < 0) return NULL;
		else if (indx < a.n) return a.m[indx];
		else return NULL;
	};
	__device__ void attach(T* v) {
		a.expand(1);
		a.m[a.n - 1] = v;
	};
	__device__ T** add(int n) {
		int n0 = a.n;
		a.expand(n);
		T* bf = new T[n];
		for (int i = 0; i < n; i++) a.m[n0 + i] = bf + i;
		return a.m + n0;
	};
	__device__ T** attach(T* v, int n) {
		int n0 = a.n;
		a.expand(n);
		for (int i = 0; i < n; i++) a.m[n0 + i] = v + i;
		return a.m + n0;
	};
	__device__ int indx(T* v) {
		for (int i = 0; i < a.n; i++) {
			if (a.m[i] == v) return i;
		}
		return -1;
	};
	__device__ void remove(int indx) {
		return a.remove_one(indx);
	};
	__device__ void del(T* v) {
		int vindx = indx(v);
		if (vindx < 0) return;
		this->remove(vindx); 
	};
	__device__ void release(bool release_context) {
		int i;
		if (release_context) {
			for (i = 0; i < a.n; i++) {delete a.m[i]; a.m[i] = NULL;}
		}
		else {
			for (i = 0; i < a.n; i++) a.m[i] = NULL;
		}
		a.release();
	};
	__device__ ~FLEX_ARRAY() {release(false);};
};
*/

template <class T> class FLEX_MATRIX {
public:
	ARRAY<T> mb; // whole memory block of the matrix
	int nx, ny;

	__host__ __device__ void _init() {mb._init(); nx = 0; ny = 0;};
	FLEX_MATRIX() {nx = 0; ny = 0;};
	__host__ void set_host(int nx, int ny) {
		mb.set_array_host(nx * ny); this->nx = nx; this->ny = ny;
	};
	__host__ void release_host() {
		mb.release_host(); nx = 0; ny = 0;
	};
	__device__ void release() {
		mb.release(); nx = 0; ny = 0;
	};
	__device__ T* mblock() {return mb.m;};
	__device__ int mem_size() {return mb.n;};
	__device__ int dim_x() {return nx;};
	__device__ int dim_y() {return ny;};
	__device__ void set(int nx, int ny) {
		this->nx = nx; this->ny = ny;
	};
	__host__ __device__ __inline__ T* e(int ix, int iy) {
		return mb.m + ix * ny + iy;
	};
	__host__ __device__ __inline__ T* el(int ix) {
		return mb.m + ix * ny;
	};

	__device__ ~FLEX_MATRIX() {mb.release();};
	__device__ void operator = (FLEX_MATRIX<T> &fm) {
		int i;
		if (fm.mb.m != NULL && this->mb.m != NULL) {
			//memcpy(this->mb.m, fm.mb.m, this->mb.n * sizeof(T));
			for (i = 0; i < this->mb.n; i++) this->mb.m[i] = fm.mb.m[i];
		}
	};
	__device__ void operator = (const int v) {
		if (mb.m == NULL) return;
		//memset(this->mb.m, v, this->mb.n * sizeof(T)); return;
		int i;
		for (i = 0; i < this->mb.n; i++) this->mb.m[i] = v;
	};
	__device__ void operator += (FLEX_MATRIX<T> &fm) {
		int i;
		if (this->mb.n != fm.mb.n || fm.mb.m != NULL || this->mb.m != NULL) return;
		T* fmb = fm.mblock(), mb = this->mb.m;
		for (i = 0; i < mb.n; i++) mb[i] += fmb[i];
	};
	__device__ void operator -= (FLEX_MATRIX<T> &fm) {
		int i;
		if (this->mb.n != fm.mb.n || fm.mb.m != NULL || this->mb.m != NULL) return;
		T* fmb = fm.mblock(), mb = this->mb.m;
		for (i = 0; i < mb.n; i++) mb[i] -= fmb[i];
	};
};

template <class T> class FLEX_CMATRIX { // flexible cubic matrix
public:
	ARRAY<T> mb; // memory block
	int nx, ny, nz;
	int nyz;

	FLEX_CMATRIX() {nx = 0; ny = 0; nz = 0; nyz = 0;};

	__host__ __device__ void _init() {nx = 0; ny = 0; nz = 0; nyz = 0; mb._init();};
	__host__ void release_host() {
		mb.release_host(); nx = 0; ny = 0; nz = 0; nyz = 0;
	};
	__device__ void release() {mb.release(); nx = 0; ny = 0; nz = 0; nyz = 0;};
	__device__ T* mblock() {return mb.m;};
	__device__ int mem_size() {return mb.n;};
	__device__ int dim_x() {return nx;};
	__device__ int dim_y() {return ny;};
	__device__ int dim_z() {return nz;};
	__host__ __device__ __inline__ T* e(int ix, int iy, int iz) {
		return mb.m + ix * nyz + iy * nz + iz;
	};
	__host__ __device__ __inline__ T* el(int ix, int iy) {
		return mb.m + ix * nyz + iy * nz;
	};
	__host__ __device__ __inline__ T* elayer(int ix) {
		return mb.m + ix * nyz;
	};

	__device__ void set(int nx, int ny, int nz) {this->nx = nx; this->ny = ny; this->nz = nz; this->nyz = ny * nz;};
	__host__ void set_host(int nx, int ny, int nz) {
		mb.set_array_host(nx * ny * nz); this->nx = nx; this->ny = ny; this->nz = nz; this->nyz = ny * nz;
	};
	__device__ ~FLEX_CMATRIX() {mb.release();};
	__device__ void operator = (FLEX_CMATRIX<T> &fm) {
		if (fm.mb.m == NULL || this->mb.m == NULL || this->mb.n != fm.mb.n) return;
		//memcpy(mb.m, fm.mblock(), mb.n * sizeof(T));
		int i;
		for (i = 0; i < this->mb.n; i++) this->mb.m[i] = fm.mb.m[i];
	};
	__device__ void operator += (FLEX_CMATRIX<T> &fm) {
		int i;
		if (fm.mb.m == NULL || this->mb.m == NULL || this->mb.n != fm.mb.n) return;
		T* fmb = fm.mb.m;
		for (i = 0; i < mb.n; i++) mb.m[i] += fmb[i];
	};
	__device__ void operator -= (FLEX_CMATRIX<T> &fm) {
		int i;
		if (fm.mb.m == NULL || this->mb.m == NULL || this->mb.n != fm.mb.n) return;
		T* fmb = fm.mblock();
		for (i = 0; i < mb.n; i++) mb.m[i] -= fmb[i];
	};
};



	class complex {
	public:
		double x, y;
		__host__ __device__ complex() {x = 0.; y = 0.;};
		__host__ __device__ complex(double x, double y) {this->x = x; this->y = y;};
		__host__ __device__ complex (double x) {this->x = x; y = 0.;};

		__host__ __device__ __inline__ void operator = (const complex &c) {
			this->x = c.x;
			this->y = c.y;
		};
		__host__ __device__ __inline__ void operator = (const cuDoubleComplex &c) {
			this->x = c.x; this->y = c.y;
		};

		__host__ __device__ __inline__ void operator = (double d) {
			this->x = d;
			this->y = 0.;
		};

		__host__ __device__ __inline__ complex operator + (const complex &c) {
			complex t;
			t.x = this->x + c.x;
			t.y = this->y + c.y;
			return t;
		};

		__host__ __device__ __inline__ complex operator + (double d) {
			complex t;
			t.x = this->x + d;
			t.y = this->y;
			return t;
		};

		__host__ __device__ __inline__ complex operator - (const complex &c) {
			complex t;
			t.x = this->x - c.x;
			t.y = this->y - c.y;
			return t;
		};

		__host__ __device__ __inline__ complex operator - (double d) {
			complex t;
			t.x = this->x - d;
			t.y = this->y;
			return t;
		};

		__host__ __device__ __inline__ complex operator * (const complex &c) {
			complex t;
			t.x = this->x * c.x - this->y * c.y;
			t.y = this->x * c.y + this->y * c.x;
			return t;
		};

		__host__ __device__ __inline__ complex operator * (double d) {
			complex t;
			t.x = this->x * d;
			t.y = this->y * d;
			return t;
		};

		__host__ __device__ __inline__ complex operator / (const complex &c) {
			complex t;
			t = (*this) * complex(c.x, -c.y) / (c.x * c.x + c.y * c.y);
			return t;
		};

		__host__ __device__ __inline__ complex operator / (double d) {
			complex t;
			t = complex(this->x / d, this->y / d);
			return t;
		};

		__host__ __device__ __inline__ void operator += (const complex &c)	{
			this->x += c.x;
			this->y += c.y;
		};

		__host__ __device__ __inline__ void operator += (double d) {
			this->x += d;
		};

		__host__ __device__ __inline__ void operator -= (const complex &c) {
			this->x -= c.x;
			this->y -= c.y;
		};
		__host__ __device__ __inline__ void operator -= (double d) {
			this->x -= d;
		};

		__host__ __device__ __inline__ void operator *= (const complex &c) {
			DOUBLE x = this->x * c.x - this->y * c.y;
			DOUBLE y = this->x * c.y + this->y * c.x;
			this->x = x; this->y = y;
		};

		__host__ __device__ __inline__ void operator *= (double d) {
			this->x *= d;
			this->y *= d;
		};

		__host__ __device__ __inline__ void operator /= (const complex &c) {
			complex t = (*this) / c;
			(*this) = t;
		};

		__host__ __device__ __inline__ void operator /= (double d) {
			this->x /= d;
			this->y /= d;
		};
	};

	__host__ __device__ static __inline__ double cabs(const complex &c) {
		return sqrt(c.x * c.x + c.y * c.y);
	};

	__host__ __device__ static __inline__ double cabs2(const complex &c) {
		return (c.x * c.x + c.y * c.y);
	};

		//#define PI 3.1415926535847
		//#define PI 3.1415926
	__host__ __device__ static __inline__ double cphase(const complex &c) {
		// phase -PI ~ PI  -- VERY IMPORTANT in reflectivity calculation
		if (c.y == 0 && c.x >= 0) return 0.;
		else if (c.y == 0 && c.x < 0) return PI;
		else if (c.x == 0 && c.y > 0) return PI/2.;
		else if (c.x == 0 && c.y < 0) return -PI / 2.;
		else if (c.x > 0 && c.y > 0) return atan(c.y / c.x);
		else if (c.x < 0 && c.y > 0) return PI + atan(c.y / c.x);
		else if (c.x < 0 && c.y < 0) return -PI + atan(c.y / c.x);
		else return atan(c.y / c.x);
	};

	__host__ __device__ static __inline__ complex conjugate(const complex &c) {
		return complex(c.x, -c.y);
	};

	__host__ __device__ __inline__ complex clog(const complex &c) {
		double v = log(sqrt(c.x * c.x + c.y * c.y));
		return complex(v, cphase(c));
	};

	__host__ __device__ __inline__ complex cexp(const complex &c) {
		complex t = complex(cos(c.y), sin(c.y)) * exp(c.x);
		return t;
	};

	__host__ __device__ __inline__ complex cpow(const complex &c, double exp) {
		complex t;
		double angle = cphase(c) * exp;
		t = complex(cos(angle), sin(angle)) * pow(cabs(c), exp);
		return t;
	};

	__host__ __device__ __inline__ complex csqrt(const complex &c) {
		complex t = cpow(c, 0.5);
		return t;
	};

// sqrt can give two results, +/-
// this function returns the result with positive imaginary part 
// this is for the x-ray reflectivity
	__host__ __device__ __inline__ complex csqrt_py(const complex &c) {
		complex t = cpow(c, 0.5);
		if (t.y < 0) {t.x = -t.x; t.y = -t.y;}
		return t;
	};


namespace _gpu_util_ {
	__host__ __device__ __inline__ void gpu_grid(int nJob, int dimBlock, unsigned int &ngrid, int gridmax = -1) {
		ngrid = nJob / dimBlock;
		if (ngrid == 0) ngrid = 1;
		if (ngrid * dimBlock < nJob) ngrid += 1;
		if (gridmax > 0 && ngrid > gridmax) ngrid = gridmax;
	};
	__host__ __device__ __inline__ void gpu_grid(int nJob, int dimBlock, int &ngrid, int gridmax = -1) {
		ngrid = nJob / dimBlock;
		if (ngrid == 0) ngrid = 1;
		if (ngrid * dimBlock < nJob) ngrid += 1;
		if (gridmax > 0 && ngrid > gridmax) ngrid = gridmax;
	};

	__host__ __device__ __inline__ int gpuJob_nEachThread(int nJob, int dimGrid, int dimBlock) {
		int nEachThread = nJob / (dimGrid * dimBlock);
		if (nEachThread == 0) nEachThread = 1;
		else if (nEachThread * dimGrid * dimBlock < nJob) nEachThread += 1;
		return nEachThread;
	};




	template <class T> class HostArray {
	public:
		int n;
		T *m;
		bool bGPU_MapMem;
		T *m_dev;

		__host__ void _init() {m = NULL; n = 0; bGPU_MapMem = false; m_dev = NULL;};

		__host__ HostArray() {
			n = 0; m = NULL; bGPU_MapMem = false;
			m_dev = NULL;
		};

		__host__ void use_GPU_MapMem(bool bGPU_MapMem) {release(); this->bGPU_MapMem = bGPU_MapMem; };

		__host__ void release() {
			if (this->bGPU_MapMem) {
				if (m_dev != NULL) {
					//cudaFree(m_dev); m_dev = NULL; // can not free this from device, although it is mapped to device
					m_dev = NULL;
				}
				if (m != NULL) {
					cudaFreeHost(m); m = NULL;
				}
			}
			else {
				if (m != NULL) {
					delete[] m; m = NULL;
				}
			}
			n = 0;
		};

		__host__ bool SetArray(int n) {
			if (this->n == n) return true;
			release(); if (n <= 0) return true;

			if (this->bGPU_MapMem) {
				if (cudaHostAlloc((void**)(&m), n * sizeof(T), cudaHostAllocMapped) != cudaSuccess) {
					//show_msg("failure to create GPU-mapped host-memory"); 
					m = NULL; n = 0; return false;
				}
				if (cudaHostGetDevicePointer((void**)&m_dev, m, 0) != cudaSuccess) {
					//show_msg("failure to get device pointer to mapped memory"); 
					release(); return false;
				}
				this->n = n;
				memset(m, 0, n * sizeof(T));
				return true;
			}
			else {
				m = new T[n]; this->n = n; return true;
			}
			
			return true;
		};

		__host__ bool set_array(int n) {return SetArray(n);};
		__host__ void release_pointer() {
			m = NULL; n = 0;
			m_dev = NULL;
		};
		~HostArray() {release();};
	};

	/*
	template <class T> __global__ void sum_array_kernel(T *a_dev, int n, T *s_block) {
		__shared__ T *sum;
		if (threadIdx.x == 0) {sum = new T[blockDim.x]; memset(sum, 0, blockDim.x * sizeof(T));}
		__syncthreads();

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = gpuJob_nEachThread(n, gridDim.x, blockDim.x);
		int i, ia;
		
		for (i = 0; i < nEachThread; i++) { // each thread
			ia = tid * nEachThread + i;
			if (ia > n) break;
			sum[threadIdx.x] += a_dev[ia];
		}
		__syncthreads();
		if (threadIdx.x == 0) {// each block
			s_block[blockIdx.x] = 0;
			for (i = 0; i < blockDim.x; i++) s_block[blockIdx.x] += sum[i];
			delete[] sum; sum = NULL;
		}
		__syncthreads();
	};

	template <class T> __global__ void sum_array_cuml_kernel(T *a1_dev, T *a2_dev, int n, T *s_block) {
		__shared__ T *sum;
		if (threadIdx.x == 0) {sum = new T[blockDim.x]; memset(sum, 0, blockDim.x * sizeof(T));}
		__syncthreads();

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		int nEachThread = gpuJob_nEachThread(n, gridDim.x, blockDim.x);
		int i, ia;
		
		for (i = 0; i < nEachThread; i++) { // each thread
			ia = tid * nEachThread + i;
			if (ia > n) break;
			sum[threadIdx.x] += a1_dev[ia] * a2_dev[ia];
		}
		__syncthreads();
		if (threadIdx.x == 0) {// each block
			s_block[blockIdx.x] = 0;
			for (i = 0; i < blockDim.x; i++) s_block[blockIdx.x] += sum[i];
			delete[] sum; sum = NULL;
		}
		__syncthreads();
	};

	template <class T> __host__ T sum_array(T *a_dev, int n) {
		int dimBlock = 128;
		int dimGrid = n / dimBlock;
		if (dimGrid == 0) dimGrid = 1;
		else if (dimGrid * dimBlock < n) dimGrid += 1;
		if (dimGrid > 32) dimGrid = 32; // do not have too many blocks
		HostArray<T> s_block; s_block.use_GPU_MapMem(true); s_block.set_array(dimGrid);

		cudaStream_t cuStream = NULL;
		cudaStreamCreate(&cuStream);

		sum_array_kernel<T><<< dimGrid, dimBlock, 256, cuStream >>>(a_dev, n, s_block.m_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		T sum = 0;
		for (int i = 0; i < s_block.n; i++) sum += s_block.m[i];
		s_block.release();
		return sum;
	};

	template <class T> __host__ T sum_array_cuml(T *a1_dev, T *a2_dev, int n) {
		int dimBlock = 128;
		int dimGrid = n / dimBlock;
		if (dimGrid == 0) dimGrid = 1;
		else if (dimGrid * dimBlock < n) dimGrid += 1;
		if (dimGrid > 32) dimGrid = 32; // do not have too many blocks
		HostArray<T> s_block; s_block.use_GPU_MapMem(true); 
		if (!s_block.set_array(dimGrid)) {
			//show_msg("failure to initialize GPU-mapped array!");
			return 0;
		}

		cudaStream_t cuStream = NULL;
		cudaStreamCreate(&cuStream);

		sum_array_cuml_kernel<T><<< dimGrid, dimBlock, 256, cuStream >>>(a1_dev, a2_dev, n, s_block.m_dev);
		cudaStreamSynchronize(cuStream);
		cudaStreamDestroy(cuStream);
		T sum = 0;
		for (int i = 0; i < s_block.n; i++) sum += s_block.m[i];
		s_block.release();
		return sum;
	};
	*/

	__host__ double sum_array(double *a_dev, int n, double *buf_dev, double *buf_host, int dim_buf);
	__host__ float sum_float_array(float *a_dev, int n, float *buf_dev, float *buf_host, int dim_buf);
	__host__ double sum_array_cuml(double *a1_dev, double *a2_dev, int n, double *buf_dev, double *buf_host, int dim_buf);
	__host__ float sum_float_array_cuml(float *a1_dev, float *a2_dev, int n, float *buf_dev, float *buf_host, int dim_buf);
	__host__ float sum_float_double_array_cuml(float *a1_dev, double *a2_dev, int n, float *buf_dev, float *buf_host, int dim_buf);

	__host__ __inline__ bool check_gpu() {
		size_t size_free, size_total;
		cudaMemGetInfo(&size_free, &size_total);
		return (size_free == 0 ? false : true);
	};

};

__device__ double atomicAdd(double *address, double val);


} // end of _gpu_dstruct_

#endif // __GPU__
