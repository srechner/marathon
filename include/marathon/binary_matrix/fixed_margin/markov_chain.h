/*
 * Created on: Nov 1, 2016
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

#ifndef MARATHON_BINARY_MATRIX_MARKOV_CHAIN_MARGIN_FIXED_H
#define MARATHON_BINARY_MATRIX_MARKOV_CHAIN_MARGIN_FIXED_H

#include "realize.h"
#include "count.h"
#include "enumerate.h"
#include "marathon/binary_matrix/markov_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Base class for all Markov chains that generate bipartite graphs / binary matrices.
             */
            class MarkovChain : public ::marathon::binary_matrix::MarkovChain {

                friend class KannanPath;

            protected:

                const Instance _inst;        // instance description

            public:

                /**
                 * Create a Markov chain for the given instance.
                 * @param inst Row and Column sums.
                 */
                explicit MarkovChain(Instance inst) :
                        ::marathon::binary_matrix::MarkovChain(realize(inst)),
                        _inst(std::move(inst)) {
                }

                /**
                 * Create a Markov chain from a certain binary matrix.
                 * @param m Binary matrix.
                 */
                explicit MarkovChain(BinaryMatrix m) :
                        ::marathon::binary_matrix::MarkovChain(std::move(m)),
                        _inst(Instance(_currentState)) {

                }

                /**
                 * Create a Markov chain for the given instance.
                 * Use a specified state as the initial state.
                 * @param inst Row and column sums.
                 * @param bin BinaryMatrix used as initial state.
                 */
                MarkovChain(Instance inst, BinaryMatrix bin) :
                        ::marathon::binary_matrix::MarkovChain(std::move(bin)),
                        _inst(std::move(inst)) {
                }


                /**
                 * Create a Markov chain instance based on a string-encoded sequence pair.
                 * Instances have the form "2,2,2;1,2,1,2".
                 * The semicolon separates the row sums from the column sums.
                 */
                explicit MarkovChain(const std::string &inst) :
                        MarkovChain(Instance(inst)) {

                }

                /**
                 * Create a Markov chain instance based on the given row and column sums.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns.
                 */
                MarkovChain(
                        const int *rowsum,
                        const int *colsum,
                        size_t nrow,
                        size_t ncol
                ) : MarkovChain(Instance(rowsum, colsum, nrow, ncol)) {

                }

                /**
                 * Create a Markov chain instance based on the given row and column sums.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 */
                MarkovChain(
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) : MarkovChain(&rowsum[0], &colsum[0], rowsum.size(), colsum.size()) {

                }

                /**
                 * Return the instance
                 * @return Return the number of rows.
                 */
                const Instance &getInstance() const {
                    return _inst;
                }

            };
        }
    }
}


#endif //MARATHON_BINARY_MATRIX_MARKOV_CHAIN_MARGIN_FIXED_H
