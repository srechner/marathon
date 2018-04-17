/*
 * Created on: Oct 19, 2016
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

#ifndef _SWITCH_CHAIN_INTELLIGENT_H_
#define _SWITCH_CHAIN_INTELLIGENT_H_

#include "marathon/binary_matrix/fixed_margin/switch_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Implements a variant of the Markov chain defined by Kannan et al.
             *
             * R. Kannan, P. Tetali, and S. Vempala.
             * Simple Markov-chain algorithms for generating bipartite graphs and tournaments.
             * Random Structures & Algorithms 14 (1997), 293–308.
             */
            class SwitchChainIntelligent : public MarkovChain {

            protected:

                // temporary memory
                int *X;
                int *Y;

            public:

                /**
                 * Create a Markov chain.
                 * @param inst Lower and upper bounds on row and column sums.
                 */
                SwitchChainIntelligent(const Instance &inst) : MarkovChain(inst) {
                    const int ncol = (int) inst.getNumCols();
                    X = new int[ncol];
                    Y = new int[ncol];
                }


                /**
                 * Create a Markov chain.
                 * @param inst String-encoded lower and upper bounds on row and column sums.
                 * Instances have the form "2,2,2;1,2,1,2".
                 * The semicolon separates the row sums from the column sums.
                 */
                SwitchChainIntelligent(const std::string &inst)
                        : SwitchChainIntelligent(Instance(inst)) {

                }

                /**
                 * Create a Markov chain as a copy of another one.
                 * @param mc Markov chain.
                 */
                SwitchChainIntelligent(const SwitchChainIntelligent &mc) : MarkovChain(mc) {

                }

                virtual ~SwitchChainIntelligent() {
                    delete[] X;
                    delete[] Y;
                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s
                 * @param process
                 */
                void adjacentStates(
                        const State &state,
                        const std::function<void(const State &, const marathon::Rational &)> &process
                ) const override {

                    // create a copy of x
                    BinaryMatrix A(static_cast<const BinaryMatrix &>(state));

                    // propability of 'lazy loop'
                    const Rational lazy(1, 100);

                    const int nrow = _inst.getNumRows();
                    const int ncol = _inst.getNumCols();

                    // collects all the loop propability
                    Rational loop = lazy;

                    // for all combinations of rows
                    for (int i = 0; i < nrow; i++) {
                        for (int j = i + 1; j < nrow; j++) {

                            // determine the set of column indices k for which A_ik != A_jk
                            int x = 0;
                            int y = 0;
                            for (int k = 0; k < ncol; k++) {
                                const bool Aik = A.get(i, k);
                                const bool Ajk = A.get(j, k);
                                if (Aik != Ajk) {
                                    if (Aik) {
                                        X[x] = k;
                                        x++;
                                    } else {
                                        Y[y] = k;
                                        y++;
                                    }
                                }
                            }

                            // no switch possible?
                            if (x == 0 || y == 0) {
                                loop += (Rational(1) - lazy) * Rational(2, nrow * (nrow - 1));
                            } else {

                                // determine proposal propability
                                const Rational p = (Rational(1) - lazy) * Rational(2, nrow * (nrow - 1) * x * y);

                                // select two random column indices that make a switchable cycle
                                for (int m = 0; m < x; m++) {
                                    for (int n = 0; n < y; n++) {

                                        // switch entries
                                        A.flip(i, X[m]);
                                        A.flip(i, Y[n]);
                                        A.flip(j, X[m]);
                                        A.flip(j, Y[n]);

                                        // process adjacent state
                                        process(A, p);

                                        // undo changes
                                        A.flip(i, X[m]);
                                        A.flip(i, Y[n]);
                                        A.flip(j, X[m]);
                                        A.flip(j, Y[n]);
                                    }
                                }
                            }
                        }
                    }

                    // process loop propability
                    process(A, loop);
                }


                /**
                 * Randomize the current state of the Markov chain.
                 * @param steps Number of steps.
                 * @return The current state after randomization.
                 */
                virtual void step() override {

                    const int nrow = (int) _inst.getNumRows();
                    const int ncol = (int) _inst.getNumCols();

                    // with propability of 1/2, remain in current state
                    if (rg.nextDouble() < 0.5)
                        return;

                    // select two different random rows
                    int i = rg.nextInt(nrow);
                    int k = rg.nextInt(nrow);
                    while (i == k)
                        k = rg.nextInt(nrow);

                    // determine the set of column indices j for which A_ij != A_kk
                    int x = 0;
                    int y = 0;
                    for (int j = 0; j < ncol; j++) {
                        bool Aij = _currentState.get(i, j);
                        bool Akj = _currentState.get(k, j);
                        if (Aij != Akj) {
                            if (Aij) {
                                X[x] = j;
                                x++;
                            } else {
                                Y[y] = j;
                                y++;
                            }
                        }
                    }

                    // select two random column indices that make a switchable cycle
                    const int m = rg.nextInt(x);
                    const int n = rg.nextInt(y);

                    // switch entries
                    _currentState.flipSubmatrix(i, X[m], k, Y[n]);
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual std::unique_ptr<marathon::MarkovChain> copy() const {
                    return std::make_unique<SwitchChainIntelligent>(_inst);
                }

                /**
                 * Return an identifier for the Markov chain.
                 * @return
                 */
                virtual std::string name() const {
                    return "switchI";
                }
            };
        }
    }
}


#endif /* _SWITCH_CHAIN_INTELLIGENT_H_ */
