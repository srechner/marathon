/*
 * Created on: Jan 18, 2017
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2017, Steffen Rechner
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

#ifndef MARATHON_BINARY_MATRIX_MARKOV_CHAIN_MARGIN_INTERVAL_H
#define MARATHON_BINARY_MATRIX_MARKOV_CHAIN_MARGIN_INTERVAL_H

#include "instance.h"
#include "marathon/binary_matrix/markov_chain.h"
#include "realize.h"
#include "count.h"

#include <sstream>

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * Base class for all Markov chains that generate bipartite graphs / binary matrices.
             */
            class MarkovChain : public ::marathon::binary_matrix::MarkovChain {

            protected:

                const interval_margin::Instance _inst;          // lower and upper matrix margins
                fixed_margin::Instance _current_margin;         // row and column sums of current state

            public:

                /**
                 * Create a Markov chain.
                 * @param inst Four-Tuple of integer vectors.
                 */
                MarkovChain(Instance inst) :
                        ::marathon::binary_matrix::MarkovChain(realize_fast(inst)),
                        _inst(std::move(inst)),
                        _current_margin(fixed_margin::Instance(_currentState)) {
                }

                /**
                 * Create a Markov chain.
                 * Use the binary matrix as initial state.
                 * @param inst Four-Tuple of integer vectors.
                 * @param m Binary matrix.
                 */
                MarkovChain(Instance inst, BinaryMatrix m) :
                        ::marathon::binary_matrix::MarkovChain(std::move(m)),
                        _inst(std::move(inst)),
                        _current_margin(fixed_margin::Instance(_currentState)) {

                }

                /**
                 * Create a Markov chain.
                 * @param inst String-encoded Instance.
                 * Instances for this chain have the form "l1-u1,l2-u2,l3-u3;l4-u4,l5-u5", where li is the ith lower
                 * bound and ui is the ith upper bound. For convenience, if li = ui, the string 'li-ui' can be replaced
                 * by 'li'. The semicolon separates the row sums from the column sums.
                 */
                explicit MarkovChain(const std::string &inst) : MarkovChain(Instance(inst)) {

                }

                /**
                 * Create a Markov chain.
                 * @param rowsum_lower Sequence of lower bounds for row sums.
                 * @param rowsum_upper Sequence of upper bounds for row sums.
                 * @param colsum_lower Sequence of lower bounds for column sums.
                 * @param colsum_upper Sequence of upper bounds for column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns.
                 */
                MarkovChain(
                        const int *rowsum_lower,
                        const int *rowsum_upper,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        size_t nrow,
                        size_t ncol
                ) : MarkovChain(Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper, nrow, ncol)) {

                }

                /**
                 * Create a Markov chain.
                 * @param rowsum_lower Sequence of lower bounds for row sums.
                 * @param rowsum_upper Sequence of upper bounds for row sums.
                 * @param colsum_lower Sequence of lower bounds for column sums.
                 * @param colsum_upper Sequence of upper bounds for column sums.
                 */
                MarkovChain(
                        const std::vector<int> &rowsum_lower,
                        const std::vector<int> &rowsum_upper,
                        const std::vector<int> &colsum_lower,
                        const std::vector<int> &colsum_upper
                ) : MarkovChain(
                        &rowsum_lower[0], &rowsum_upper[0], &colsum_lower[0], &colsum_upper[0],
                        rowsum_lower.size(), colsum_lower.size()
                ) {

                }

                /**
                 * Set the current state.
                 * @param s Binary matrix.
                 */
                virtual void setCurrentState(const State &s) override {
                    auto m = dynamic_cast<const BinaryMatrix *>(&s);
                    if (m == nullptr)
                        throw std::runtime_error("Error! State is not a binary matrix!");
                    _currentState = *m;
                    _current_margin = fixed_margin::Instance(_currentState);
                }
            };
        }
    }
}


#endif //MARATHON_BINARY_MATRIX_MARKOV_CHAIN_MARGIN_INTERVAL_H
