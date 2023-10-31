class complex {
public:
	double x, y;
	complex() {x = 0.; y = 0.;};
	complex(double x, double y) {this->x = x; this->y = y;};
	complex (double x) {this->x = x; y = 0.;};

	void operator = (const complex& c);
	void operator = (double d);
	complex operator + (const complex& c);
	complex operator + (double d);
	complex operator - (const complex& c);
	complex operator - (double d);
	complex operator * (const complex& c);
	complex operator * (double d);
	complex operator / (const complex& c);
	complex operator / (double d);
	void operator += (const complex& c);
	void operator += (double d);
	void operator -= (const complex& c);
	void operator -= (double d);
	void operator *= (const complex& c);
	void operator *= (double d);
	void operator /= (const complex& c);
	void operator /= (double d);
};

double abs(const complex& c);
double abs2(const complex& c);
double phase(const complex& c);  // in rad. 0 ~ 2PI
complex pow(const complex& c, double exp);
complex sqrt(const complex& c);
complex sqrt_py(const complex& c);
complex exp(const complex& c);
complex log(const complex& c);
complex sin(const complex& c);
complex cos(const complex& c);
complex tan(const complex& c);
complex conjugate(const complex& c);

