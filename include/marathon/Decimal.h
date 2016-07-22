/*
 * Rational.h
 *
 *  Created on: Mar 27, 2015
 *      Author: rechner
 */

#ifndef FLOAT_H
#define FLOAT_H

#include <ostream>
#include "Rational.h"
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace marathon {

/**
 * just a wrapper around boost cpp_dec_float data type
 */
template<int precision>
class decimal {

private:
	boost::multiprecision::number<boost::multiprecision::cpp_dec_float<precision>> _d;

public:
	decimal() :
			_d(0) {

	}

	decimal(const rational& o) :
			_d(o._r) {
	}

	void operator=(const decimal<precision>& o) {
		_d = o._d;
	}

	void operator+=(const decimal<precision>& o) {
		_d += o._d;
	}

	void operator-=(const decimal<precision>& o) {
		_d -= o._d;
	}

	void operator*=(const decimal<precision>& o) {
		_d *= o._d;
	}

	void operator/=(const decimal<precision>& o) {
		_d /= o._d;
	}

	decimal<precision> operator*(const decimal<precision>& o) const {
		return _d * o._d;
	}

	decimal<precision> operator-(const decimal<precision>& o) const {
		return _d - o._d;
	}

	decimal<precision> operator+(const decimal<precision>& o) const {
		return _d + o._d;
	}

	decimal<precision> operator/(const decimal<precision>& o) const {
		return _d / o._d;
	}

	bool operator==(const decimal<precision>& o) const {
		return _d == o._d;
	}

	bool operator!=(const decimal<precision>& o) const {
		return _d != o._d;
	}

	bool operator<(const decimal<precision>& o) const {
		return _d < o._d;
	}

	bool operator>(const decimal<precision>& o) const {
		return _d > o._d;
	}

	std::string to_string(int digits = precision) const {
		return _d.str(digits);
	}
};

template<int precision>
std::ostream& operator<<(std::ostream& out, const decimal<precision>& r) {
	out << r.to_string();
	return out;
}

}

#endif /* FLOAT_H */
