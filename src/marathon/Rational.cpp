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

rational::rational(long n) {
	_r = boost::multiprecision::cpp_rational(n);
}

rational::rational(long num, long denom) {
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

std::string rational::to_string() const {
	return _r.str();
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
	out << r.to_string();
	return out;
}

}

#endif
