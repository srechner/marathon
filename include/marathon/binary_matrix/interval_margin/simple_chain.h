/*
 * Created on: Jan 17, 2017
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *         Linda Strowick <linda.strowick@student.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2017, Steffen Rechner, Linda Strowick
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

#ifndef _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_SIMPLE_CHAIN_H
#define _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_SIMPLE_CHAIN_H

#include <vector>
#include <stack>
#include "marathon/binary_matrix/interval_margin/markov_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * Simple Markov chain defined in
             *
             * Steffen Rechner, Linda Strowick, Matthias Müller-Hannemann.
             * Uniform sampling of bipartite graphs with degrees in prescribed intervals.
             * Journal of Complex Networks (2017). DOI: 10.1093/comnet/cnx059
             */
            class SimpleChain : public MarkovChain {

                friend class SwitchChainIntelligent;

            protected:

                // the probabilities of selecting each kind of operation type
                const Rational p_switch_rat = Rational(1, 3);
                const Rational p_shift_rat = Rational(1, 3);
                const Rational p_flip_rat = Rational(1, 3);


                void simulateSwitch(
                        const BinaryMatrix *s,
                        const std::function<void(const State *, const Rational &)> &process,
                        const int *rowsum,
                        const int *colsum
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    Rational loop(0);

                    // simulate switch operations
                    const Rational p_switch =
                            p_switch_rat * Rational(4, ncol * (ncol - 1) * nrow * (nrow - 1));
                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {
                            for (int k = i + 1; k < nrow; k++) {
                                for (int l = j + 1; l < ncol; l++) {

                                    // if Switch is possible
                                    if (s->isCheckerBoardUnit(i, j, k, l)) {

                                        BinaryMatrix s2(*s);

                                        // switch the cycle
                                        s2.flipSubmatrix(i, j, k, l);

                                        // process state
                                        process(&s2, p_switch);
                                    } else {
                                        loop += p_switch;
                                    }
                                }
                            }
                        }
                    }

                    process(s, loop);
                }

                void simulateFlip(
                        const BinaryMatrix *s,
                        const std::function<void(const State *, const Rational &)> &process,
                        const int *rowsum,
                        const int *colsum
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    Rational loop(0);

                    // simulate flip operations
                    const Rational p_flip = p_flip_rat * Rational(1, ncol * nrow);
                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {

                            //Flip 1 -> 0
                            if (s->get(i, j)) {

                                // flip is allowed by lower and upper margins?
                                if (rowsum[i] > inst.rowsum_lower[i] && colsum[j] > inst.colsum_lower[j]) {
                                    BinaryMatrix s2(*s);
                                    s2.flip(i, j);
                                    process(&s2, p_flip);
                                } else {
                                    loop += p_flip;
                                }
                            } else {   //Flip 0 -> 1

                                // flip is allowed by lower and upper margins?
                                if (rowsum[i] < inst.rowsum_upper[i] && colsum[j] < inst.colsum_upper[j]) {
                                    BinaryMatrix s2(*s);
                                    s2.flip(i, j);
                                    process(&s2, p_flip);
                                } else {
                                    loop += p_flip;
                                }
                            }
                        }
                    }

                    process(s, loop);
                }

                void simulateShift(
                        const BinaryMatrix *s,
                        const std::function<void(const State *, const Rational &)> &process,
                        const int *rowsum,
                        const int *colsum
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    Rational loop(0);

                    // simulate shifts
                    const Rational p_shift = p_shift_rat * Rational(1, nrow * ncol * (nrow + ncol - 2));
                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {

                            // simulate choosing k as row index
                            for (int k = 0; k < nrow; k++) {
                                if (k != i) {

                                    // if we can shift a one from (i,j) to (k,j)
                                    if (s->get(i, j) && !s->get(k, j) &&
                                        rowsum[i] > inst.rowsum_lower[i] && rowsum[k] < inst.rowsum_upper[k]) {

                                        BinaryMatrix s2(*s);

                                        // shift a one from (i,j) to (k,j)
                                        s2.flip(i, j);
                                        s2.flip(k, j);
                                        process(&s2, p_shift);
                                    } // if we can shift a one from (k,j) to (i,j)
                                    else if (!s->get(i, j) && s->get(k, j) &&
                                             rowsum[k] > inst.rowsum_lower[k] && rowsum[i] < inst.rowsum_upper[i]) {

                                        BinaryMatrix s2(*s);

                                        // shift a one from (k,j) to (i,j)
                                        s2.flip(i, j);
                                        s2.flip(k, j);
                                        process(&s2, p_shift);
                                    } else {
                                        loop += p_shift;
                                    }
                                }
                            }

                            // simulate choosing k as row index
                            for (int k = 0; k < ncol; k++) {
                                if (k != j) {

                                    // if we can shift a one from (i,j) to (i,k)
                                    if (s->get(i, j) && !s->get(i, k) &&
                                        colsum[j] > inst.colsum_lower[j] && colsum[k] < inst.colsum_upper[k]) {

                                        BinaryMatrix s2(*s);

                                        // shift a one from (i,j) to (i,k)
                                        s2.flip(i, j);
                                        s2.flip(i, k);
                                        process(&s2, p_shift);
                                    } // if we can shift a one from (i,k) to (i,j)
                                    else if (!s->get(i, j) && s->get(i, k) &&
                                             colsum[k] > inst.colsum_lower[k] && colsum[j] < inst.colsum_upper[j]) {

                                        BinaryMatrix s2(*s);

                                        // shift a one from (i,k) to (i,j)
                                        s2.flip(i, j);
                                        s2.flip(i, k);
                                        process(&s2, p_shift);
                                    } else {
                                        loop += p_shift;
                                    }
                                }
                            }
                        }
                    }

                    process(s, loop);
                }


                inline void applySwitch(BinaryMatrix &s) {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    /** apply switch operation **/
                    count_switch++;

                    // choose two random row indices i < j
                    int i = rg.nextInt(nrow);
                    int k = rg.nextInt(nrow);
                    while (i == k)
                        k = rg.nextInt(nrow);
                    if (i > k)
                        std::swap(i, k);

                    // choose two random column indices i < j
                    int j = rg.nextInt(ncol);
                    int l = rg.nextInt(ncol);
                    while (j == l) {
                        l = rg.nextInt(ncol);
                    }
                    if (j > l)
                        std::swap(j, l);

                    // check if edges are flippable
                    if (s.isCheckerBoardUnit(i, j, k, l)) {
                        s.flipSubmatrix(i, j, k, l); // switch!
                    } else {
                        // loop
                        loop_switch++;
                    }
                }

                void applyShift(BinaryMatrix& s) {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    /** apply a shift operation **/
                    count_shift++;

                    // select a random row index i
                    int i = rg.nextInt(nrow);

                    // select a random column index j
                    int j = rg.nextInt(ncol);

                    // select a random index k (row or column)
                    int k = rg.nextInt(nrow + ncol);
                    while (k == i || k - nrow == j)
                        k = rg.nextInt(nrow + ncol);

                    if (k < nrow) {      // if k is a row index

                        // if we can shift a one from (i,j) to (k,j)
                        if (s.get(i, j) && !s.get(k, j) &&
                            current_margin.rowsum[i] > inst.rowsum_lower[i] &&
                            current_margin.rowsum[k] < inst.rowsum_upper[k]) {

                            // shift a one from (i,j) to (k,j)
                            s.flip(i, j);
                            s.flip(k, j);
                            current_margin.rowsum[i]--;
                            current_margin.rowsum[k]++;
                        }
                            // if we can shift a one from (k,j) to (i,j)
                        else if (!s.get(i, j) && s.get(k, j) &&
                                 current_margin.rowsum[k] > inst.rowsum_lower[k] &&
                                 current_margin.rowsum[i] < inst.rowsum_upper[i]) {

                            // shift a one from (k,j) to (i,j)
                            s.flip(i, j);
                            s.flip(k, j);
                            current_margin.rowsum[k]--;
                            current_margin.rowsum[i]++;
                        } else {
                            // loop
                            loop_shift++;
                        }

                    } else {            // if k is a column index

                        k -= nrow;

                        // if we can shift a one from (i,j) to (i,k)
                        if (s.get(i, j) && !s.get(i, k) &&
                            current_margin.colsum[j] > inst.colsum_lower[j] &&
                            current_margin.colsum[k] < inst.colsum_upper[k]) {

                            // shift a one from (i,j) to (i,k)
                            s.flip(i, j);
                            s.flip(i, k);
                            current_margin.colsum[j]--;
                            current_margin.colsum[k]++;
                        } // if we can shift a one from (i,k) to (i,j)
                        else if (!s.get(i, j) && s.get(i, k) &&
                                 current_margin.colsum[k] > inst.colsum_lower[k] &&
                                 current_margin.colsum[j] < inst.colsum_upper[j]) {

                            // shift a one from (i,k) to (i,j)
                            s.flip(i, j);
                            s.flip(i, k);
                            current_margin.colsum[k]--;
                            current_margin.colsum[j]++;
                        } else {
                            // loop
                            loop_shift++;
                        }
                    }
                }


                inline void applyFlip(BinaryMatrix& s) {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    /** apply flip operation **/

                    count_flip++;

                    // select a random row index i and a random column index j
                    int i = rg.nextInt(nrow);
                    int j = rg.nextInt(ncol);

                    if (s.get(i, j)) {
                        // Flip 1 -> 0

                        // flip is allowed by lower and upper margins?
                        if (current_margin.rowsum[i] > inst.rowsum_lower[i] &&
                            current_margin.colsum[j] > inst.colsum_lower[j]) {
                            s.flip(i, j);
                            current_margin.rowsum[i]--;
                            current_margin.colsum[j]--;
                        } else {
                            // loop
                            loop_flip++;
                        }
                    } else {
                        //Flip 0 -> 1

                        // flip is allowed by lower and upper margins?
                        if (current_margin.rowsum[i] < inst.rowsum_upper[i] &&
                            current_margin.colsum[j] < inst.colsum_upper[j]) {
                            s.flip(i, j);
                            current_margin.rowsum[i]++;
                            current_margin.colsum[j]++;
                        } else {
                            // loop
                            loop_flip++;
                        }
                    }
                }

            public:

                // comparison with double is much faster
                double p_switch = p_switch_rat.convert_to<double>();
                double p_flip = p_flip_rat.convert_to<double>();
                double p_shift = p_shift_rat.convert_to<double>();

                // number of every operation
                size_t count_switch = 0;
                size_t count_flip = 0;
                size_t count_shift = 0;

                // loop count of every operation
                size_t loop_switch = 0;
                size_t loop_flip = 0;
                size_t loop_shift = 0;

                /**
                 * Create a Markov chain.
                 * @param inst Four-Tuple of integer vectors.
                 */
                explicit SimpleChain(const Instance &inst) : MarkovChain(inst) {

                }

                /**
                 * Create a Markov chain. Use a certain binary matrix as initial state.
                 * @param inst Four-Tuple of integer vectors.
                 * @param bin Binary Matrix.
                 */
                SimpleChain(const Instance &inst, const BinaryMatrix &bin)
                        : MarkovChain(inst, bin) {

                }

                /**
                 * Create a Markov chain.
                 * @param inst String-encoded instance.
                 * Instances for this chain have the form "l1-u1,l2-u2,l3-u3;l4-u4,l5-u5", where li is the ith lower
                 * bound and ui is the ith upper bound. For convenience, if li = ui, the string 'li-ui' can be replaced
                 * by 'li'. The semicolon separates the row sums from the column sums.
                 */
                explicit SimpleChain(const std::string &inst)
                        : SimpleChain(Instance(inst)) {

                }

                /**
				 * Create a Markov chain.
				 * @param rowsum_lower Sequence of lower bounds for row sums.
				 * @param rowsum_upper Sequence of upper bounds for row sums.
				 * @param colsum_lower Sequence of lower bounds for column sums.
				 * @param colsum_upper Sequence of upper bounds for column sums.
				 */
                SimpleChain(
                        const std::vector<int> &rowsum_lower,
                        const std::vector<int> &rowsum_upper,
                        const std::vector<int> &colsum_lower,
                        const std::vector<int> &colsum_upper
                ) : MarkovChain(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper) {

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
                SimpleChain(
                        const int *rowsum_lower,
                        const int *rowsum_upper,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        const int nrow,
                        const int ncol
                ) : MarkovChain(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper, nrow, ncol) {

                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s
                 * @param process
                 */
                virtual
                void adjacentStates(
                        const State *x,
                        const std::function<void(const State *, const marathon::Rational &)> &process
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    // create a copy of x
                    const BinaryMatrix *s = (const BinaryMatrix *) x;

                    // create temporary array of row and column sums
                    int *rowsum = new int[nrow];
                    int *colsum = new int[ncol];
                    memset(rowsum, 0, nrow * sizeof(int));
                    memset(colsum, 0, ncol * sizeof(int));
                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {
                            int sij = s->get(i, j);
                            rowsum[i] += sij;
                            colsum[j] += sij;
                        }
                    }

                    simulateSwitch(s, process, rowsum, colsum);
                    simulateFlip(s, process, rowsum, colsum);
                    simulateShift(s, process, rowsum, colsum);
                }

                /**
                 * Randomize the state s by applying a single transition.
                 */
                virtual void step() override {

                    // decide whether to use flip, switch, or shift
                    double p = rg.nextDouble();

                    if (p < p_switch)
                        applySwitch(currentState);
                    else if (p < p_switch + p_shift)
                        applyShift(currentState);
                    else if (p < p_switch + p_shift + p_flip)
                        applyFlip(currentState);
                    else
                        throw std::runtime_error("Error at SwitchShiftFlipSimple::step()");
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual std::unique_ptr<marathon::MarkovChain> copy() const override {
                    return std::make_unique<SimpleChain>(inst, currentState);
                }
            };
        }
    }
}


#endif /* _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_SIMPLE_CHAIN_H */
