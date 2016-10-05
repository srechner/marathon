/*
 * Rational.h
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

#include "Types.h"
#include <ostream>
#include "boost/multiprecision/cpp_int.hpp"

namespace marathon {

	//	typedef boost::multiprecision::cpp_rational rational;

	/**
	 * just a wrapper around boost rational data type
	 */
	class rational {

		template<int precision>
		friend
		class decimal;

	private:
		boost::multiprecision::cpp_rational _r;

	public:
		rational() :
				_r(0) {
		}

		rational(const rational &o) :
				_r(o._r) {
		}

		rational(boost::multiprecision::cpp_rational r) :
				_r(r) {
		}

		rational(long n) {
			_r = boost::multiprecision::cpp_rational(n);
		}

		rational(long num, long denom) {
			_r = boost::multiprecision::cpp_rational(num, denom);
		}

		rational(integer num, integer denom) {
			_r = boost::multiprecision::cpp_rational(num, denom);
		}

		bool operator==(const rational &o) const {
			return _r == o._r;
		}

		bool operator!=(const rational &o) const {
			return _r != o._r;
		}

		void operator+=(const rational &o) {
			_r += o._r;
		}

		void operator-=(const rational &o) {
			_r -= o._r;
		}

		void operator*=(const rational &o) {
			_r *= o._r;
		}

		rational operator*(const rational &o) const {
			boost::multiprecision::cpp_rational res = _r * o._r;
			return rational(res);
		}

		rational operator-(const rational &o) const {
			boost::multiprecision::cpp_rational res = _r - o._r;
			return rational(res);
		}

		rational operator+(const rational &o) const {
			boost::multiprecision::cpp_rational res = _r + o._r;
			return rational(res);
		}

		rational operator/(const rational &o) const {
			boost::multiprecision::cpp_rational res = _r / o._r;
			return rational(res);
		}

		void operator/=(const rational &o) {
			_r /= o._r;
		}

		std::string to_string() const {
			return _r.str();
		}

		bool operator<(const rational &o) const {
			return _r < o._r;
		}

		bool operator>(const rational &o) const {
			return _r > o._r;
		}

		void operator=(const rational &o) {
			_r = o._r;
		}

		template<typename T>
		T convert_to() const {
			return _r.convert_to<T>();
		}
	};

	inline
	std::ostream &operator<<(std::ostream &out, const rational &r) {
		out << r.to_string();
		return out;
	}
}

#endif /* RATIONAL_H_ */
