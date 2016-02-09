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

namespace marathon {

/**
 * just a wrapper around boost rational data type
 */
class rational {

private:
	boost::multiprecision::cpp_rational _r;

public:
	rational();
	rational(const rational& o);
	rational(boost::multiprecision::cpp_rational r);
	rational(int n);
	rational(int num, int denom);

	void operator=(const rational& o);

	bool operator==(const rational& o);
	bool operator!=(const rational& o);
	void operator+=(const rational& o);
	void operator-=(const rational& o);
	void operator*=(const rational& o);
	void operator/=(const rational& o);

	rational operator*(const rational& o) const;
	rational operator-(const rational& o) const;
	rational operator/(const rational& o) const;

	bool operator<(const rational& o) const;
	bool operator>(const rational& o) const;

	void stream_to(std::ostream& os) const;

	template<typename T>
	T convert_to() const {
		return _r.convert_to<T>();
	}
};

std::ostream& operator<<(std::ostream& out, const rational& r);

}

#endif /* RATIONAL_H_ */
