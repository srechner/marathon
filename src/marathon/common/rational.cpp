#ifndef RATIONAL_CPP_
#define RATIONAL_CPP_

#include "../../../include/marathon/common/rational.h"


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

bool rational::operator==(const rational& o) {
	return _r == o._r;
}

bool rational::operator!=(const rational& o) {
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

rational rational::operator/(const rational& o) const {
	boost::multiprecision::cpp_rational res = _r / o._r;
	return rational(res);
}

void rational::operator/=(const rational& o) {
	_r /= o._r;
}

void rational::stream_to(std::ostream& os) const {
	os << _r;
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
	r.stream_to(out);
	return out;
}

}

#endif
