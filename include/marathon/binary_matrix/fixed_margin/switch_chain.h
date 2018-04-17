/*
 * Created on: Nov 24, 2014
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

#ifndef MARATHON_BINARY_MATRIX_FIXED_MARGIN_SWITCH_CHAIN_H
#define MARATHON_BINARY_MATRIX_FIXED_MARGIN_SWITCH_CHAIN_H

#include <vector>
#include <stack>
#include "marathon/binary_matrix/fixed_margin/markov_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Implements the Markov chain defined by
             *
             * R. Kannan, P. Tetali, and S. Vempala.
             * Simple Markov-chain algorithms for generating bipartite graphs and tournaments.
             * Random Structures & Algorithms 14 (1997), 293–308.
             */
            class SwitchChain : public MarkovChain {

            public:

                /**
                 * Create a Markov chain for the given instance.
                 */
                explicit SwitchChain(Instance seq) : MarkovChain(std::move(seq)) {

                }

                /**
                 * Create a Markov chain for the given instance.
                 * Use a specified state as initial state.
                 * @param inst Row and Column sums.
                 * @param bin BinaryMatrix used as initial state.
                 */
                SwitchChain(Instance inst, BinaryMatrix bin)
                        : MarkovChain(std::move(inst), std::move(bin)) {

                }

                /**
                 * Create a Markov chain instance based on a string-encoded sequence pair.
                 * Instances have the form "2,2,2;1,2,1,2".
                 * The semicolon separates the row sums from the column sums.
                 */
                explicit SwitchChain(const std::string &inst) : MarkovChain(inst) {

                }

                /**
				 * Create a Markov chain instance based on the given row and column sums.
				 * @param rowsum Sequence of row sums.
				 * @param colsum Sequence of column sums.
				 */
                SwitchChain(
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) : MarkovChain(rowsum, colsum) {

                }

                /**
				 * Create a Markov chain instance based on the given row and column sums.
				 * @param rowsum Sequence of row sums.
				 * @param colsum Sequence of column sums.
				 * @param nrow Number of rows.
				 * @param ncol Number of columns.
				 */
                SwitchChain(
                        const int *rowsum,
                        const int *colsum,
                        size_t nrow,
                        size_t ncol
                ) : MarkovChain(rowsum, colsum, nrow, ncol) {

                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s State object.
                 * @param process Function object evaluated for each adjacent state.
                 */
                void adjacentStates(
                        const State &x,
                        const std::function<void(const State &, const marathon::Rational &)> &process
                ) const override {

                    const size_t nrow = _inst.getNumRows();
                    const size_t ncol = _inst.getNumCols();

                    // create a copy of x
                    BinaryMatrix s(static_cast<const BinaryMatrix &> (x));

                    // each state has proposal prob. of p
                    const Rational p(4, nrow * (nrow - 1) * ncol * (ncol - 1));

                    Rational loop(0);

                    // Definition of Kannan, Tetali, Vempala
                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {
                            for (int k = i + 1; k < nrow; k++) {
                                for (int l = j + 1; l < ncol; l++) {

                                    /**
                                     * A switch is possible in the following situations
                                     *
                                     * a) alternating cycle ( (i,j)=1, (k,j)=0, (k,l)=1, (i,l)=0 )
                                     *         i   k
                                     *       1 2 3 4 5
                                     *     -----------
                                     * j 1 | x 1 x 0 x
                                     * l 2 | x 0 x 1 x
                                     *   3 | x x x x x
                                     *
                                     * b) symmetric case (attention: not regarded in paper!)
                                     *         i   k
                                     *       1 2 3 4 5
                                     *     -----------
                                     * j 1 | x 0 x 1 x
                                     * l 2 | x 1 x 0 x
                                     *   3 | x x x x x
                                     *
                                     */

                                    if (s.isCheckerBoardUnit(i, j, k, l)) {

                                        // switch the cycle
                                        s.flipSubmatrix(i, j, k, l);

                                        // process adjacent state
                                        process(s, p);

                                        // undo modifications
                                        s.flipSubmatrix(i, j, k, l);

                                    } else {
                                        // loop
                                        loop += p;
                                    }
                                }
                            }
                        }
                    }
                    // create loop
                    process(s, loop);
                }

                /**
                 * Randomize the current state of the Markov chain.
                 */
                virtual void step() override {

                    const int nrow = (int) _inst.getNumRows();
                    const int ncol = (int) _inst.getNumCols();

                    // select four random integers i,j,k,l

                    // select two different row indices
                    int i = rg.nextInt(nrow);
                    int k = rg.nextInt(nrow);
                    while (i == k)
                        k = rg.nextInt(nrow);

                    // select two different column indices
                    int j = rg.nextInt(ncol);
                    int l = rg.nextInt(ncol);
                    while (j == l)
                        l = rg.nextInt(ncol);

                    // if i,j,k,l makes a switchable cycle
                    if (_currentState.isCheckerBoardUnit(i, j, k, l)) {
                        _currentState.flipSubmatrix(i, j, k, l);      // switch the cycle
                    }
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual std::unique_ptr<marathon::MarkovChain> copy() const override {
                    return std::make_unique<SwitchChain>(_inst, _currentState);
                }
            };
        }
    }
}


#endif /* MARATHON_BINARY_MATRIX_FIXED_MARGIN_SWITCH_CHAIN_H */
