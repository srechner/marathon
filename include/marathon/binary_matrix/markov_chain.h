/*
 * Created on: Apr 4, 2018
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

#ifndef MARATHON_BINARY_MATRIX_MARKOV_CHAIN_H
#define MARATHON_BINARY_MATRIX_MARKOV_CHAIN_H

#include "marathon/markov_chain.h"

namespace marathon {
    namespace binary_matrix {
        /**
         * Base class for all Markov chains that generate bipartite graphs / binary matrices.
         */
        class MarkovChain : public ::marathon::MarkovChain {

        protected:

            BinaryMatrix _currentState;  // current state
            const size_t _nrow;         // number of rows
            const size_t _ncol;         // number of columns

        public:

            /**
             * Create a Markov chain object with a initial state.
             * @param bin Initial state.
             */
            MarkovChain(BinaryMatrix bin) :
                    _currentState(std::move(bin)),
                    _nrow(_currentState.getNumRows()),
                    _ncol(_currentState.getNumCols()) {

            }

            /**
             * Set the current state.
             * @param s Binary matrix.
             */
            virtual void setCurrentState(const State &s) override {

                // try to cast s to the correct type
                auto m = dynamic_cast<const BinaryMatrix *>(&s);
                if (m == nullptr)
                    throw std::runtime_error("Error! State is not a binary matrix!");

                _currentState = *m;
            }

            /**
             * Return the current state of the Markov chain.
             * @return
             */
            virtual const BinaryMatrix &getCurrentState() const override {
                return _currentState;
            }

        };
    }
}


#endif //MARATHON_BINARY_MATRIX_MARKOV_CHAIN_MARGIN_FIXED_H
