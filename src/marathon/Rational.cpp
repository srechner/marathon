/*
 * Rational.cpp
 *
 * Created on: Mar 27, 2015
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2016, Steffen Rechner
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef RATIONAL_CPP_
#define RATIONAL_CPP_

#include "../../include/marathon/Rational.h"

namespace marathon {

// Rational: just a wrapper around boost Rational data type
	rational::rational() :
			_r(0) {
	}

	rational::rational(const rational &o) :
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

	rational::rational(integer num, integer denom) {
		_r = boost::multiprecision::cpp_rational(num, denom);
	}

	bool rational::operator==(const rational &o) const {
		return _r == o._r;
	}

	bool rational::operator!=(const rational &o) const {
		return _r != o._r;
	}

	void rational::operator+=(const rational &o) {
		_r += o._r;
	}

	void rational::operator-=(const rational &o) {
		_r -= o._r;
	}

	void rational::operator*=(const rational &o) {
		_r *= o._r;
	}

	rational rational::operator*(const rational &o) const {
		boost::multiprecision::cpp_rational res = _r * o._r;
		return rational(res);
	}

	rational rational::operator-(const rational &o) const {
		boost::multiprecision::cpp_rational res = _r - o._r;
		return rational(res);
	}

	rational rational::operator+(const rational &o) const {
		boost::multiprecision::cpp_rational res = _r + o._r;
		return rational(res);
	}

	rational rational::operator/(const rational &o) const {
		boost::multiprecision::cpp_rational res = _r / o._r;
		return rational(res);
	}

	void rational::operator/=(const rational &o) {
		_r /= o._r;
	}

	std::string rational::to_string() const {
		return _r.str();
	}

	bool rational::operator<(const rational &o) const {
		return _r < o._r;
	}

	bool rational::operator>(const rational &o) const {
		return _r > o._r;
	}

	void rational::operator=(const rational &o) {
		_r = o._r;
	}

	std::ostream &operator<<(std::ostream &out, const rational &r) {
		out << r.to_string();
		return out;
	}

}

#endif
