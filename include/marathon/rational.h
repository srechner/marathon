/*
 * Rational.h
 *
 *  Created on: Mar 27, 2015
 *      Author: rechner
 */

#ifndef RATIONAL_H_
#define RATIONAL_H_

#include <ostream>
#include "boost/multiprecision/cpp_int.hpp"

//  just a wrapper around boost rational data type

namespace marathon {

class Rational {

private:
	boost::multiprecision::cpp_rational _r;

public:
	Rational();
	Rational(const Rational& o);
	Rational(boost::multiprecision::cpp_rational r);
	Rational(int n);
	Rational(int num, int denom);

	void operator=(const Rational& o);

	bool operator==(const Rational& o);
	bool operator!=(const Rational& o);
	void operator+=(const Rational& o);
	void operator-=(const Rational& o);
	void operator*=(const Rational& o);
	void operator/=(const Rational& o);

	Rational operator*(const Rational& o) const;
	Rational operator-(const Rational& o) const;
	Rational operator/(const Rational& o) const;

	bool operator<(const Rational& o) const;
	bool operator>(const Rational& o) const;

	void stream_to(std::ostream& os) const;

	template<typename T>
	T convert_to() const {
		return _r.convert_to<T>();
	}
};

std::ostream& operator<<(std::ostream& out, const Rational& r);

}

#endif /* RATIONAL_H_ */
