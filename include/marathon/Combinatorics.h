/*
 * Combinatorics.h
 *
 * Created on: May 11, 2016
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

#ifndef INCLUDE_MARATHON_COMBINATORICS_H_
#define INCLUDE_MARATHON_COMBINATORICS_H_

#include "Types.h"
#include "Rational.h"

namespace marathon {

	class Combinatorics {

	private:

		// declare friend functions
		friend void marathon::init();
		friend void marathon::cleanup();

		// variables for computation of binomial coefficients
		static int _nrow, _ncol;    // number of rows and columns in _bimom table
		static integer *_binom;     // table with computed binomonial coefficients

	public:

		/**
		 * Initialize the Combinatorics Component.
		 */
		static void init();

		/**
		 * Cleanup the Combinatorics Component.
		 */
		static void cleanup();


		/**
		 * Computes n choose k.
		 */
		static integer binom(const int n, const int k);

		/**
		 * Computes n!.
		 */
		static integer factorial(const int n);

	};
}

#endif /* INCLUDE_MARATHON_COMBINATORICS_H_ */
