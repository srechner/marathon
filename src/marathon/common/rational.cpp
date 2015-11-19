#ifndef RATIONAL_CPP_
#define RATIONAL_CPP_

#include "../../../include/marathon/common/rational.h"


namespace marathon {

// Rational: just a wrapper around boost Rational data type
Rational::Rational() :
		_r(0) {

}

Rational::Rational(const Rational& o) :
		_r(o._r) {
}

Rational::Rational(boost::multiprecision::cpp_rational r) :
		_r(r) {
}

Rational::Rational(int n) {
	_r = boost::multiprecision::cpp_rational(n);
}

Rational::Rational(int num, int denom) {
	_r = boost::multiprecision::cpp_rational(num, denom);
}

bool Rational::operator==(const Rational& o) {
	return _r == o._r;
}

bool Rational::operator!=(const Rational& o) {
	return _r != o._r;
}

void Rational::operator+=(const Rational& o) {
	_r += o._r;
}

void Rational::operator-=(const Rational& o) {
	_r -= o._r;
}

void Rational::operator*=(const Rational& o) {
	_r *= o._r;
}

Rational Rational::operator*(const Rational& o) const {
	boost::multiprecision::cpp_rational res = _r * o._r;
	return Rational(res);
}

Rational Rational::operator-(const Rational& o) const {
	boost::multiprecision::cpp_rational res = _r - o._r;
	return Rational(res);
}

Rational Rational::operator/(const Rational& o) const {
	boost::multiprecision::cpp_rational res = _r / o._r;
	return Rational(res);
}

void Rational::operator/=(const Rational& o) {
	_r /= o._r;
}

void Rational::stream_to(std::ostream& os) const {
	os << _r;
}

bool Rational::operator<(const Rational& o) const {
	return _r < o._r;
}

bool Rational::operator>(const Rational& o) const {
	return _r > o._r;
}

void Rational::operator=(const Rational& o) {
	_r = o._r;
}

std::ostream& operator<<(std::ostream& out, const Rational& r) {
	r.stream_to(out);
	return out;
}

}

#endif
