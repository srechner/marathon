/*
 * Created on: Nov 05, 2017
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

#ifndef MARATHON_BINARY_MATRIX_FIXED_MARGIN_EDGE_SWITCH_CHAIN_H
#define MARATHON_BINARY_MATRIX_FIXED_MARGIN_EDGE_SWITCH_CHAIN_H

#include <vector>
#include <stack>
#include "marathon/binary_matrix/fixed_margin/markov_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Implements a variant of the Markov chain defined by Kannan et al.
             */
            class EdgeSwitchChain : public MarkovChain {

            protected:

                const size_t total;                         // number of edges
                std::vector<std::pair<uint, uint>> edges;   // edge array


                /**
                 * Create an edge array form a binary matrix.
                 * @param M
                 * @return
                 */
                std::vector<std::pair<uint, uint>> initEdges(const BinaryMatrix &M) const {

                    std::vector<std::pair<uint, uint>> edges(M.getTotal());

                    const size_t nrow = M.getNumRows();
                    const size_t ncol = M.getNumCols();
                    int k = 0;
                    for (uint i = 0; i < nrow; i++) {
                        for (uint j = 0; j < ncol; j++) {
                            if (M.get(i, j)) {
                                edges[k].first = i;
                                edges[k].second = j;
                                k++;
                            }
                        }
                    }
                    return edges;
                }

            public:

                /**
                 * Create a Markov chain for the given instance.
                 * @param seq Bi-graphical pair of integer vectors.
                 */
                explicit EdgeSwitchChain(Instance seq)
                        : MarkovChain(std::move(seq)),
                          total(_inst.getTotal()),
                          edges(initEdges(_currentState)) {

                }

                /**
                 * Create a Markov chain for the given instance.
                 * Use a specified state as initial state.
                 * @param inst Row and Column sums.
                 * @param bin BinaryMatrix used as initial state.
                 */
                EdgeSwitchChain(Instance inst, BinaryMatrix bin)
                        : MarkovChain(std::move(inst), std::move(bin)),
                          total(_currentState.getTotal()),
                          edges(initEdges(_currentState)) {
                }

                /**
                 * Create a Markov chain instance based on a string-encoded sequence pair.
                 * Instances have the form "2,2,2;1,2,1,2".
                 * The semicolon separates the row sums from the column sums.
                 */
                explicit EdgeSwitchChain(const std::string &inst)
                        : EdgeSwitchChain(Instance(inst)) {

                }

                /**
				 * Create a Markov chain instance based on the given row and column sums.
				 * @param rowsum Sequence of row sums.
				 * @param colsum Sequence of column sums.
				 */
                EdgeSwitchChain(
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) : EdgeSwitchChain(&rowsum[0], &colsum[0], rowsum.size(), colsum.size()) {

                }

                /**
				 * Create a Markov chain instance based on the given row and column sums.
				 * @param rowsum Sequence of row sums.
				 * @param colsum Sequence of column sums.
				 * @param nrow Number of rows.
				 * @param ncol Number of columns.
				 */
                EdgeSwitchChain(
                        const int *rowsum,
                        const int *colsum,
                        size_t nrow,
                        size_t ncol
                ) : EdgeSwitchChain(Instance(rowsum, colsum, nrow, ncol)) {

                }

                /**
                 * Randomize the current state of the Markov chain.
                 */
                virtual void step() override {

                    const int nrow = (int) _inst.getNumRows();
                    const int ncol = (int) _inst.getNumCols();

                    // select two edges
                    int a = rg.nextInt(total);     // a = {i,j}
                    int b = rg.nextInt(total);     // b = {k,l}
                    while (a == b)
                        b = rg.nextInt(total);

                    int i = edges[a].first;
                    int j = edges[a].second;
                    int k = edges[b].first;
                    int l = edges[b].second;

                    assert(_currentState.get(i, j));
                    assert(_currentState.get(k, l));

                    // if i,j,k,l is switchable
                    if (!_currentState.get(i, l) && !_currentState.get(k, j)) {

                        // switch
                        _currentState.set(i, j, 0);
                        _currentState.set(k, l, 0);
                        _currentState.set(i, l, 1);
                        _currentState.set(k, j, 1);

                        // adjust edge array
                        edges[a] = std::make_pair(i, l);
                        edges[b] = std::make_pair(k, j);
                    }
                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s
                 * @param process
                 */
                void adjacentStates(
                        const State &x,
                        const std::function<void(const State &, const marathon::Rational &)> &process
                ) const override {

                    const int nrow = (int) _inst.getNumRows();
                    const int ncol = (int) _inst.getNumCols();

                    // create a copy of x
                    BinaryMatrix s(static_cast<const BinaryMatrix &>(x));

                    auto edges = initEdges(s);

                    // each choice has a prob. of p
                    const Rational p(2, total * (total - 1));

                    Rational loop(0);

                    // select two edges
                    for (int a = 0; a < total; a++) {
                        for (int b = a + 1; b < total; b++) {

                            int i = edges[a].first;
                            int j = edges[a].second;
                            int k = edges[b].first;
                            int l = edges[b].second;

                            assert(s.get(i, j));
                            assert(s.get(k, l));

                            // if i,j,k,l is switchable
                            if (!s.get(i, l) && !s.get(k, j)) {

                                // switch
                                s.set(i, j, 0);
                                s.set(k, l, 0);
                                s.set(i, l, 1);
                                s.set(k, j, 1);

                                process(s, p);

                                // undo modification
                                s.set(i, j, 1);
                                s.set(k, l, 1);
                                s.set(i, l, 0);
                                s.set(k, j, 0);

                            } else {
                                loop += p;
                            }
                        }
                    }

                    // create loop
                    process(x, loop);
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual std::unique_ptr<marathon::MarkovChain> copy() const override {
                    return std::make_unique<EdgeSwitchChain>(_inst, _currentState);
                }

            };
        }
    }
}


#endif /* MARATHON_BINARY_MATRIX_FIXED_MARGIN_EDGE_SWITCH_CHAIN_H */
