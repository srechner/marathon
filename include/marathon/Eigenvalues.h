/*
 * Eigenvalues.h
 *
 *  Created on: Mar 24, 2016
 *      Author: rechner
 */

#ifndef INCLUDE_MARATHON_EIGENVALUES_H_
#define INCLUDE_MARATHON_EIGENVALUES_H_

#include "StateGraph.h"

namespace marathon {
    namespace eigenvalue {

        /**
         * Options for computation of eigenvalues.
         */
        enum eigenvalue_t {
            // eigenvalue options
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
