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
#include "marathon/markov_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Base class for all Markov chains that generate bipartite graphs / binary matrices.
             */
            class MarkovChain : public ::marathon::MarkovChain {

                friend class KannanPath;

            protected:

                const Instance inst;

            public:

                /**
                 * Create a Markov chain for the given instance.
                 * @param inst Row and Column sums.
                 */
                explicit MarkovChain(const Instance &seq) :
                        ::marathon::MarkovChain(realize(seq)),
                        inst(seq) {
                }

                /**
                 * Create a Markov chain for the given instance.
                 * Use a specified state as initial state.
                 * @param inst Row and Column sums.
                 * @param bin BinaryMatrix used as initial state.
                 */
                MarkovChain(const Instance &inst, const BinaryMatrix &bin) :
                        ::marathon::MarkovChain(&bin),
                        inst(inst) {
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
                        const int nrow,
                        const int ncol
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
                 * Create a Markov chain instance based where m is used as initial state.
                 * @param m Binary matrix used as initial state
                 */
                explicit MarkovChain(const BinaryMatrix &m) :
                        ::marathon::MarkovChain(&m),
                        inst(Instance(m)) {
                }

                /**
                 * Create a Markov chain instance as a copy of another one.
                 * @param mc Markov chain object.
                 */
                explicit MarkovChain(const MarkovChain &mc) :
                        ::marathon::MarkovChain(mc.getCurrentState()),
                        inst(mc.inst) {

                }


                /**
                 * Return the instance
                 * @return Return the number of rows.
                 */
                const Instance &getInstance() const {
                    return inst;
                }


                /**
                 * Return the current state of the Markov chain.
                 * @return
                 */
                virtual const BinaryMatrix *getCurrentState() const override {
                    return (BinaryMatrix *) currentState;
                }

                /**
                 * Randomize the current state of the Markov chain.
                 * @param steps Number of steps.
                 * @return The current state after randomization.
                 */
                virtual const BinaryMatrix *randomize(int steps) override {
                    return (BinaryMatrix *) marathon::MarkovChain::randomize(steps);
                }

                /**
                 * Create a copy of the Markov chain object.
                 * @return Copy of the Markov chain.
                 */
                virtual MarkovChain* copy() const override = 0;
            };
        }
    }
}


#endif //MARATHON_BINARY_MATRIX_MARKOV_CHAIN_MARGIN_FIXED_H
