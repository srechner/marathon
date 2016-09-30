/*
 * Eigenvalues.h
 *
 * Created on: Mar 24, 2016
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

#ifndef INCLUDE_MARATHON_EIGENVALUES_H_
#define INCLUDE_MARATHON_EIGENVALUES_H_

#include "StateGraph.h"

namespace marathon {

    namespace Eigenvalue {

        /**
         * Options for computation of eigenvalues.
         */
        enum eigenvalue_t {
            // Eigenvalue options
            _2ndLargestMagnitude,
            _2ndLargestAlgebraic,
        };

        /**
         * Computes the eigenvalue with second largest magnitute of the
         * transition matrix of mc.
         * @param mc State Graph Representation of Markov Chain.
         * @param which Which Eigenvalue to compute. Options are: 2nd largest in magnitute and 2nd largest algebraic.
         * @return The corresponding Eigenvalue.
         */
        template<typename T>
        T eigenvalue(const StateGraph *mc, eigenvalue_t which);
    }
}

#endif /* INCLUDE_MARATHON_EIGENVALUES_H_ */
