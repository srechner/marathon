/*
 * rational.h
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

#ifndef RATIONAL_H_
#define RATIONAL_H_

#include "integer.h"
#include <ostream>
#include "boost/multiprecision/cpp_int.hpp"

namespace marathon {

	/**
	 * just a wrapper around boost Rational data type
	 */
	class Rational {

	private:
		boost::multiprecision::cpp_rational _r;

	public:
		Rational() :
				_r(0) {
		}

		Rational(const Rational &o) :
				_r(o._r) {
		}

		Rational(boost::multiprecision::cpp_rational r) :
				_r(r) {
		}

		Rational(long n) {
			_r = boost::multiprecision::cpp_rational(n);
		}

		Rational(long num, long denom) {
			_r = boost::multiprecision::cpp_rational(num, denom);
		}

		Rational(const Integer& n) {
			_r = boost::multiprecision::cpp_rational(n);
		}

		Rational(const Integer& num, const Integer& denom) {
			_r = boost::multiprecision::cpp_rational(num, denom);
		}

		bool operator==(const Rational &o) const {
			return _r == o._r;
		}

		bool operator!=(const Rational &o) const {
			return _r != o._r;
		}

		void operator+=(const Rational &o) {
			_r += o._r;
		}

		void operator-=(const Rational &o) {
			_r -= o._r;
		}

		void operator*=(const Rational &o) {
			_r *= o._r;
		}

		Rational operator*(const Rational &o) const {
			boost::multiprecision::cpp_rational res = _r * o._r;
			return Rational(res);
		}

		Rational operator-(const Rational &o) const {
			boost::multiprecision::cpp_rational res = _r - o._r;
			return Rational(res);
		}

		Rational operator+(const Rational &o) const {
			boost::multiprecision::cpp_rational res = _r + o._r;
			return Rational(res);
		}

		Rational operator/(const Rational &o) const {
			boost::multiprecision::cpp_rational res = _r / o._r;
			return Rational(res);
		}

		Rational operator/(const Integer &o) const {
			return operator/(Rational(o,1));
		}

		Rational operator/(const size_t &o) const {
			return operator/(Rational(o,1));
		}

		void operator/=(const Rational &o) {
			_r /= o._r;
		}

		std::string to_string() const {
			return _r.str();
		}

		bool operator<(const Rational &o) const {
			return _r < o._r;
		}

		bool operator<=(const Rational &o) const {
			return _r <= o._r;
		}

		bool operator>=(const Rational &o) const {
			return _r >= o._r;
		}

		bool operator<=(const double &c) const {
			return _r.convert_to<double>() <= c;
		}

		bool operator>=(const double &c) const {
			return _r.convert_to<double>() >= c;
		}

		bool operator<(const double &c) const {
			return _r.convert_to<double>() < c;
		}

		bool operator>(const double &c) const {
			return _r.convert_to<double>() > c;
		}

		bool operator>(const Rational &o) const {
			return _r > o._r;
		}

		void operator=(const Rational &o) {
			_r = o._r;
		}

		template<typename T>
		T convert_to() const {
			return _r.convert_to<T>();
		}

		Integer getNumerator() const {
			boost::multiprecision::cpp_int num = numerator(_r);
			return Integer(num);
		}

		Integer getDenominator() const {
			boost::multiprecision::cpp_int denom = denominator(_r);
			return Integer(denom);
		}
	};

	inline
	std::ostream &operator<<(std::ostream &out, const Rational &r) {
		out << r.to_string();
		return out;
	}

	template<>
	inline
	Rational Rational::convert_to<Rational>() const {
		return _r;
	}
}


#endif /* RATIONAL_H_ */
