//#include <math.h>
#include <iostream>
#include <cmath>
using namespace std;
#include "complex.h"

double abs(const complex& c)
{
	return sqrt(c.x * c.x + c.y * c.y);
}

double abs2(const complex& c)
{
	return (c.x * c.x + c.y * c.y);
}

#define PI 3.1415926535847
//#define PI 3.1415926
double phase(const complex& c)
{
	// phase -PI ~ PI  -- VERY IMPORTANT in reflectivity calculation
	if (c.y == 0 && c.x >= 0) return 0.;
	else if (c.y == 0 && c.x < 0) return PI;
	else if (c.x == 0 && c.y > 0) return PI/2.;
	else if (c.x == 0 && c.y < 0) return -PI / 2.;
	else if (c.x > 0 && c.y > 0) return atan(c.y / c.x);
	else if (c.x < 0 && c.y > 0) return PI + atan(c.y / c.x);
	else if (c.x < 0 && c.y < 0) return -PI + atan(c.y / c.x);
	else return atan(c.y / c.x);
}

complex conjugate(const complex &c)
{
	return complex(c.x, -c.y);
}

complex log(const complex &c) 
{
	complex t;
	t = complex(log(abs(c)), phase(c));
	return t;
}

complex exp(const complex &c)
{
	complex t = complex(cos(c.y), sin(c.y)) * exp(c.x);
	return t;
}

complex pow(const complex &c, double exp)
{
	complex t;
	double angle = phase(c) * exp;
	t = complex(cos(angle), sin(angle)) * pow(abs(c), exp);
	return t;
}

complex sqrt(const complex &c)
{
	complex t = pow(c, 0.5);
	return t;
}

// sqrt can give two results, +/-
// this function returns the result with positive imaginary part 
// this is for the x-ray reflectivity
complex sqrt_py(const complex &c)
{
	complex t = pow(c, 0.5);
	if (t.y < 0) {t.x = -t.x; t.y = -t.y;}
	return t;
}

void complex::operator = (const complex& c)
{
	this->x = c.x;
	this->y = c.y;
}

void complex::operator = (double d)
{
	this->x = d;
	this->y = 0.;
}

complex complex::operator + (const complex &c)
{
	complex t;
	t.x = this->x + c.x;
	t.y = this->y + c.y;
	return t;
}

complex complex::operator + (double d)
{
	complex t;
	t.x = this->x + d;
	t.y = this->y;
	return t;
}

complex complex::operator - (const complex &c)
{
	complex t;
	t.x = this->x - c.x;
	t.y = this->y - c.y;
	return t;
}

complex complex::operator - (double d)
{
	complex t;
	t.x = this->x - d;
	t.y = this->y;
	return t;
}

complex complex::operator * (const complex &c)
{
	complex t;
	t.x = this->x * c.x - this->y * c.y;
	t.y = this->x * c.y + this->y * c.x;
	return t;
}

complex complex::operator * (double d)
{
	complex t;
	t.x = this->x * d;
	t.y = this->y * d;
	return t;
}

complex complex::operator / (const complex &c)
{
	complex t;
	t = (*this) * conjugate(c) / abs2(c);
	return t;
}

complex complex::operator / (double d)
{
	complex t;
	t = complex(this->x / d, this->y / d);
	return t;
}

void complex::operator += (const complex &c)
{
	this->x += c.x;
	this->y += c.y;
}

void complex::operator += (double d)
{
	this->x += d;
}

void complex::operator -= (const complex &c)
{
	this->x -= c.x;
	this->y -= c.y;
}
void complex::operator -= (double d)
{
	this->x -= d;
}

void complex::operator *= (const complex &c)
{
	double x = this->x * c.x - this->y * c.y;
	double y = this->x * c.y + this->y * c.x;
	this->x = x; this->y = y;
}

void complex::operator *= (double d)
{
	this->x *= d;
	this->y *= d;
}

void complex::operator /= (const complex &c)
{
	complex t = (*this) / c;
	(*this) = t;
}

void complex::operator /= (double d)
{
	this->x /= d;
	this->y /= d;
}



