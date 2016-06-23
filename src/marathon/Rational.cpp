#ifndef RATIONAL_CPP_
#define RATIONAL_CPP_

#include "../../include/marathon/Rational.h"

namespace marathon {

// Rational: just a wrapper around boost Rational data type
rational::rational() :
		_r(0) {

}

rational::rational(const rational& o) :
		_r(o._r) {
}

rational::rational(boost::multiprecision::cpp_rational r) :
		_r(r) {
}

rational::rational(int n) {
	_r = boost::multiprecision::cpp_rational(n);
}

rational::rational(int num, int denom) {
	_r = boost::multiprecision::cpp_rational(num, denom);
}

bool rational::operator==(const rational& o) const {
	return _r == o._r;
}

bool rational::operator!=(const rational& o) const {
	return _r != o._r;
}

void rational::operator+=(const rational& o) {
	_r += o._r;
}

void rational::operator-=(const rational& o) {
	_r -= o._r;
}

void rational::operator*=(const rational& o) {
	_r *= o._r;
}

rational rational::operator*(const rational& o) const {
	boost::multiprecision::cpp_rational res = _r * o._r;
	return rational(res);
}

rational rational::operator-(const rational& o) const {
	boost::multiprecision::cpp_rational res = _r - o._r;
	return rational(res);
}

rational rational::operator+(const rational& o) const {
	boost::multiprecision::cpp_rational res = _r + o._r;
	return rational(res);
}

rational rational::operator/(const rational& o) const {
	boost::multiprecision::cpp_rational res = _r / o._r;
	return rational(res);
}

void rational::operator/=(const rational& o) {
	_r /= o._r;
}

std::string rational::to_string_fraction() const {
	return _r.str();
}

std::string rational::to_string_dec_float(int precision) const {
	boost::multiprecision::cpp_dec_float_100 x(_r);
	return x.str(precision);
}


bool rational::operator<(const rational& o) const {
	return _r < o._r;
}

bool rational::operator>(const rational& o) const {
	return _r > o._r;
}

void rational::operator=(const rational& o) {
	_r = o._r;
}

std::ostream& operator<<(std::ostream& out, const rational& r) {
	out << r.to_string_fraction();
	return out;
}

}

#endif
