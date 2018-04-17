/*
 * curveball.h
 *
 * Created on: Aug 24, 2015
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

#ifndef CURVEBALL_H_
#define CURVEBALL_H_

#include "marathon/binary_matrix/fixed_margin/markov_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Markov chain defined in
             *
             *   G. Strona, D. Nappo, F. Boccacci, S. Fattorini, and J. San-Miguel-Ayanz.
             *   A fast and unbiased procedure to randomize ecological binary matrices with
             *   fixed row and column totals. Nature communications 5 (2014).
             *   doi: 10.1038/ncomms5114.
             *
             * Alternatively described in
             *
             *   C. J. Carstens.
             *   Proof of uniform sampling of binary matrices with fixed row sums and column sums
             *   for the fast curveball algorithm. Physical Review E 91 (2015), 042812.
             *   doi: 10.1103/physreve.91.042812.
             */
            class Curveball : public MarkovChain {

            protected:

                // auxiliary array
                std::vector<int> tmp1;

            public:

                /**
                 * Create a Markov chain.
                 * @param inst Row and column sums.
                 */
                explicit Curveball(Instance inst) : MarkovChain(std::move(inst)) {
                    tmp1.resize(_currentState.getNumCols());
                }

                /**
                * Create a Markov chain.
                * @param m Binary matrix used as initial state
                */
                explicit Curveball(BinaryMatrix m) : MarkovChain(std::move(m)) {
                    tmp1.resize(_currentState.getNumCols());
                }

                /**
                 * Create a Markov chain for the given instance.
                 * Use a specified state as initial state.
                 * @param inst Row and Column sums.
                 * @param bin BinaryMatrix used as initial state.
                 */
                Curveball(Instance inst, BinaryMatrix bin)
                        : MarkovChain(std::move(inst), std::move(bin)) {
                    tmp1.resize(_currentState.getNumCols());
                }

                /**
                 * Create a Markov chain.
                 * @param inst String-encoded instance.
                 * Instances have the form "2,2,2;1,2,1,2".
                 * The semicolon separates the row sums from the column sums.
                 */
                explicit Curveball(const std::string &inst) : Curveball(Instance(inst)) {

                }

                /**
                * Create a Markov chain.
                * @param rowsum Sequence of row sums.
                * @param colsum Sequence of column sums.
                */
                Curveball(
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) : Curveball(Instance(rowsum, colsum)) {

                }

                /**
                 * Create a Markov chain.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns.
                 */
                Curveball(
                        const int *rowsum,
                        const int *colsum,
                        size_t nrow,
                        size_t ncol
                ) : Curveball(Instance(rowsum, colsum, nrow, ncol)) {

                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s
                 * @param process
                 */
                virtual void adjacentStates(
                        const State &x,
                        const std::function<void(const State &, const marathon::Rational &)> &f
                ) const override {

                    const int nrow = (int) _inst.getNumRows();
                    const int ncol = (int) _inst.getNumCols();

                    // convert state reference
                    const BinaryMatrix &X = static_cast<const BinaryMatrix &>(x);

                    // Make a copy of the current state
                    BinaryMatrix A(X);

                    // auxiliary array
                    std::vector<int> tmp1(ncol);
                    std::vector<int> tmp2(ncol);

                    /**
                     * Definition of Strona et. al: A fast and unbiased procedure
                     * to randomize ecological binary matrices with fixed row and
                     * column totals.
                     */

                    // randomly select two row indices
                    for (int i = 0; i < nrow; i++) {
                        for (int k = i + 1; k < nrow; k++) {

                            // select the indices that occur in on the Ai and Ak, but not in both
                            int a = 0;
                            int b = 0;

                            // for each column position
                            for (int j = 0; j < ncol; j++) {

                                int A_ij = X.get(i, j);
                                int A_kj = X.get(k, j);

                                if (A_ij != A_kj) {

                                    // store index j in temporary array
                                    tmp1[a + b] = j;

                                    if (A_ij)
                                        a++;
                                    else
                                        b++;
                                }
                            }

                            // calculate the probability of this choice
                            const Integer num_subsets = binom(a + b, a);
                            const Integer num_row_sel = binom(nrow, 2);
                            const Rational p(1, num_subsets * num_row_sel);

                            // simulate all combinations of choosing a out of (a+b) columns
                            CombinationGenerator<int> cg(&tmp1[0], &tmp2[0], a + b, a);
                            do {

                                // set A[i,j]=0 and A[k,j]=1 for j in tmp1
                                for (int l = 0; l < a + b; l++) {
                                    int j = tmp1[l];
                                    A.set(i, j, 0);
                                    A.set(k, j, 1);
                                }

                                // for each selected column j: set A[i,j]=1 and A[k,j]=0
                                for (int l = 0; l < a; l++) {
                                    int j = tmp2[l];
                                    A.set(i, j, 1);
                                    A.set(k, j, 0);
                                }

                                // process adjacent state
                                f(A, p);
                                assert(_inst.isValid(A));

                            } while (cg.next());

                            // restore original state of rows i and k
                            for (int l = 0; l < a + b; l++) {
                                const int j = tmp1[l];
                                A.set(i, j, X.get(i, j));
                                A.set(k, j, X.get(k, j));
                            }
                        }
                    }
                }

                /**
                 * Randomize the current state of the Markov chain.
                 */
                virtual void step() override {

                    const size_t nrow = _inst.getNumRows();
                    const size_t ncol = _inst.getNumCols();

                    if (nrow == 1 || ncol == 1)
                        return;

                    // randomly select two row indices
                    int i = rg.nextInt(nrow);
                    int k = rg.nextInt(nrow);
                    while (i == k)
                        k = rg.nextInt(nrow);

                    // select the indices that occur in on the Ai and Ak, but not in both
                    int a = 0;
                    int b = 0;

                    // for each column position
                    for (int j = 0; j < ncol; j++) {

                        bool A_ij = _currentState.get(i, j);
                        bool A_kj = _currentState.get(k, j);

                        if (A_ij != A_kj) {

                            // store index j in temporary array
                            tmp1[a + b] = j;

                            if (A_ij)
                                a++;
                            else
                                b++;
                        }
                    }

                    // randomly select a out of (a+b) elements to row i
                    rg.shuffle<int>(&tmp1[0], a + b);

                    for (int l = 0; l < a; l++) {
                        int j = tmp1[l];
                        _currentState.set(i, j, 1);
                        _currentState.set(k, j, 0);
                    }

                    // the remaining elements to go row j
                    for (int l = a; l < a + b; l++) {
                        int j = tmp1[l];
                        _currentState.set(i, j, 0);
                        _currentState.set(k, j, 1);
                    }
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual std::unique_ptr<marathon::MarkovChain> copy() const override {
                    return std::make_unique<Curveball>(_inst, _currentState);
                }
            };
        }
    }
}


#endif /* CURVEBALL_H_ */
