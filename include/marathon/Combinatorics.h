/*
 * Combinatorics.h
 *
 *  Created on: May 11, 2016
 *      Author: rechner
 */

#ifndef INCLUDE_MARATHON_COMBINATORICS_H_
#define INCLUDE_MARATHON_COMBINATORICS_H_

#include "Rational.h"

namespace marathon {
namespace combinatorics {

// variables for computation of binomial coefficients
extern int _nrow, _ncol;	// number of rows and columns in _bimom table
extern rational* _binom;	// table with computed binomonial coefficients

/**
 * Computes n choose k.
 */
rational binom(const int n, const int k);

/**
 * Computes n!.
 */
rational faculty(const int n);

}
}

#endif /* INCLUDE_MARATHON_COMBINATORICS_H_ */
