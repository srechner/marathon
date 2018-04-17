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

            protected:

                // the probabilities of selecting each kind of operation type
                const Rational p_switch_rat = Rational(1, 3);
                const Rational p_shift_rat = Rational(1, 3);
                const Rational p_flip_rat = Rational(1, 3);


                void simulateSwitch(
                        const BinaryMatrix &s,
                        const std::function<void(const State &, const Rational &)> &process,
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) const {

                    Rational loop(0);

                    // simulate switch operations
                    const Rational p_switch =
                            p_switch_rat * Rational(4, _ncol * (_ncol - 1) * _nrow * (_nrow - 1));
                    for (size_t i = 0; i < _nrow; i++) {
                        for (size_t j = 0; j < _ncol; j++) {
                            for (size_t k = i + 1; k < _nrow; k++) {
                                for (size_t l = j + 1; l < _ncol; l++) {

                                    // if Switch is possible
                                    if (s.isCheckerBoardUnit(i, j, k, l)) {

                                        BinaryMatrix s2(s);

                                        // switch the cycle
                                        s2.flipSubmatrix(i, j, k, l);

                                        // process state
                                        process(s2, p_switch);
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
                        const BinaryMatrix &s,
                        const std::function<void(const State &, const Rational &)> &process,
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) const {

                    Rational loop(0);

                    // simulate flip operations
                    const Rational p_flip = p_flip_rat * Rational(1, _ncol * _nrow);
                    for (size_t i = 0; i < _nrow; i++) {
                        for (size_t j = 0; j < _ncol; j++) {

                            //Flip 1 -> 0
                            if (s.get(i, j)) {

                                // flip is allowed by lower and upper margins?
                                if (rowsum[i] > _inst._rowsum_lower[i] && colsum[j] > _inst._colsum_lower[j]) {
                                    BinaryMatrix s2(s);
                                    s2.flip(i, j);
                                    process(s2, p_flip);
                                } else {
                                    loop += p_flip;
                                }
                            } else {   //Flip 0 -> 1

                                // flip is allowed by lower and upper margins?
                                if (rowsum[i] < _inst._rowsum_upper[i] && colsum[j] < _inst._colsum_upper[j]) {
                                    BinaryMatrix s2(s);
                                    s2.flip(i, j);
                                    process(s2, p_flip);
                                } else {
                                    loop += p_flip;
                                }
                            }
                        }
                    }

                    process(s, loop);
                }

                void simulateShift(
                        const BinaryMatrix &s,
                        const std::function<void(const State &, const Rational &)> &process,
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) const {

                    Rational loop(0);

                    // simulate shifts
                    const Rational p_shift = p_shift_rat * Rational(1, _nrow * _ncol * (_nrow + _ncol - 2));
                    for (size_t i = 0; i < _nrow; i++) {
                        for (size_t j = 0; j < _ncol; j++) {

                            // simulate choosing k as row index
                            for (size_t k = 0; k < _nrow; k++) {
                                if (k != i) {

                                    // if we can shift a one from (i,j) to (k,j)
                                    if (s.get(i, j) && !s.get(k, j) &&
                                        rowsum[i] > _inst._rowsum_lower[i] && rowsum[k] < _inst._rowsum_upper[k]) {

                                        BinaryMatrix s2(s);

                                        // shift a one from (i,j) to (k,j)
                                        s2.flip(i, j);
                                        s2.flip(k, j);
                                        process(s2, p_shift);
                                    } // if we can shift a one from (k,j) to (i,j)
                                    else if (!s.get(i, j) && s.get(k, j) &&
                                             rowsum[k] > _inst._rowsum_lower[k] && rowsum[i] < _inst._rowsum_upper[i]) {

                                        BinaryMatrix s2(s);

                                        // shift a one from (k,j) to (i,j)
                                        s2.flip(i, j);
                                        s2.flip(k, j);
                                        process(s2, p_shift);
                                    } else {
                                        loop += p_shift;
                                    }
                                }
                            }

                            // simulate choosing k as row index
                            for (size_t k = 0; k < _ncol; k++) {
                                if (k != j) {

                                    // if we can shift a one from (i,j) to (i,k)
                                    if (s.get(i, j) && !s.get(i, k) &&
                                        colsum[j] > _inst._colsum_lower[j] && colsum[k] < _inst._colsum_upper[k]) {

                                        BinaryMatrix s2(s);

                                        // shift a one from (i,j) to (i,k)
                                        s2.flip(i, j);
                                        s2.flip(i, k);
                                        process(s2, p_shift);
                                    } // if we can shift a one from (i,k) to (i,j)
                                    else if (!s.get(i, j) && s.get(i, k) &&
                                             colsum[k] > _inst._colsum_lower[k] && colsum[j] < _inst._colsum_upper[j]) {

                                        BinaryMatrix s2(s);

                                        // shift a one from (i,k) to (i,j)
                                        s2.flip(i, j);
                                        s2.flip(i, k);
                                        process(s2, p_shift);
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

                    /** apply switch operation **/
                    count_switch++;

                    // choose two random row indices i < j
                    size_t i = rg.nextInt(_nrow);
                    size_t k = rg.nextInt(_nrow);
                    while (i == k)
                        k = rg.nextInt(_nrow);
                    if (i > k)
                        std::swap(i, k);

                    // choose two random column indices i < j
                    size_t j = rg.nextInt(_ncol);
                    size_t l = rg.nextInt(_ncol);
                    while (j == l) {
                        l = rg.nextInt(_ncol);
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

                void applyShift(BinaryMatrix &s) {

                    /** apply a shift operation **/
                    count_shift++;

                    // select a random row index i
                    size_t i = rg.nextInt(_nrow);

                    // select a random column index j
                    size_t j = rg.nextInt(_ncol);

                    // select a random index k (row or column)
                    size_t k = rg.nextInt(_nrow + _ncol);
                    while (k == i || k - _nrow == j)
                        k = rg.nextInt(_nrow + _ncol);

                    if (k < _nrow) {      // if k is a row index

                        // if we can shift a one from (i,j) to (k,j)
                        if (s.get(i, j) && !s.get(k, j) &&
                            _current_margin._rowsum[i] > _inst._rowsum_lower[i] &&
                            _current_margin._rowsum[k] < _inst._rowsum_upper[k]) {

                            // shift a one from (i,j) to (k,j)
                            s.flip(i, j);
                            s.flip(k, j);
                            _current_margin._rowsum[i]--;
                            _current_margin._rowsum[k]++;
                        }
                            // if we can shift a one from (k,j) to (i,j)
                        else if (!s.get(i, j) && s.get(k, j) &&
                                 _current_margin._rowsum[k] > _inst._rowsum_lower[k] &&
                                 _current_margin._rowsum[i] < _inst._rowsum_upper[i]) {

                            // shift a one from (k,j) to (i,j)
                            s.flip(i, j);
                            s.flip(k, j);
                            _current_margin._rowsum[k]--;
                            _current_margin._rowsum[i]++;
                        } else {
                            // loop
                            loop_shift++;
                        }

                    } else {            // if k is a column index

                        k -= _nrow;

                        // if we can shift a one from (i,j) to (i,k)
                        if (s.get(i, j) && !s.get(i, k) &&
                            _current_margin._colsum[j] > _inst._colsum_lower[j] &&
                            _current_margin._colsum[k] < _inst._colsum_upper[k]) {

                            // shift a one from (i,j) to (i,k)
                            s.flip(i, j);
                            s.flip(i, k);
                            _current_margin._colsum[j]--;
                            _current_margin._colsum[k]++;
                        } // if we can shift a one from (i,k) to (i,j)
                        else if (!s.get(i, j) && s.get(i, k) &&
                                 _current_margin._colsum[k] > _inst._colsum_lower[k] &&
                                 _current_margin._colsum[j] < _inst._colsum_upper[j]) {

                            // shift a one from (i,k) to (i,j)
                            s.flip(i, j);
                            s.flip(i, k);
                            _current_margin._colsum[k]--;
                            _current_margin._colsum[j]++;
                        } else {
                            // loop
                            loop_shift++;
                        }
                    }
                }


                inline void applyFlip(BinaryMatrix &s) {

                    /** apply flip operation **/

                    count_flip++;

                    // select a random row index i and a random column index j
                    size_t i = rg.nextInt(_nrow);
                    size_t j = rg.nextInt(_ncol);

                    if (s.get(i, j)) {
                        // Flip 1 -> 0

                        // flip is allowed by lower and upper margins?
                        if (_current_margin._rowsum[i] > _inst._rowsum_lower[i] &&
                            _current_margin._colsum[j] > _inst._colsum_lower[j]) {
                            s.flip(i, j);
                            _current_margin._rowsum[i]--;
                            _current_margin._colsum[j]--;
                        } else {
                            // loop
                            loop_flip++;
                        }
                    } else {
                        //Flip 0 -> 1

                        // flip is allowed by lower and upper margins?
                        if (_current_margin._rowsum[i] < _inst._rowsum_upper[i] &&
                            _current_margin._colsum[j] < _inst._colsum_upper[j]) {
                            s.flip(i, j);
                            _current_margin._rowsum[i]++;
                            _current_margin._colsum[j]++;
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
                explicit SimpleChain(Instance inst) : MarkovChain(std::move(inst)) {

                }

                /**
                 * Create a Markov chain. Use a certain binary matrix as initial state.
                 * @param inst Four-Tuple of integer vectors.
                 * @param bin Binary Matrix.
                 */
                SimpleChain(Instance inst, BinaryMatrix bin)
                        : MarkovChain(std::move(inst), std::move(bin)) {

                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s
                 * @param process
                 */
                virtual
                void adjacentStates(
                        const State &x,
                        const std::function<void(const State &, const marathon::Rational &)> &process
                ) const override {

                    // create a copy of x
                    const BinaryMatrix &s = static_cast<const BinaryMatrix &>(x);

                    // create temporary array of row and column sums
                    std::vector<int> rowsum(_nrow);
                    std::vector<int> colsum(_ncol);
                    for (size_t i = 0; i < _nrow; i++) {
                        for (size_t j = 0; j < _ncol; j++) {
                            bool sij = s.get(i, j);
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
                        applySwitch(_currentState);
                    else if (p < p_switch + p_shift)
                        applyShift(_currentState);
                    else if (p < p_switch + p_shift + p_flip)
                        applyFlip(_currentState);
                    else
                        throw std::runtime_error("Error at SwitchShiftFlipSimple::step()");
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual std::unique_ptr<marathon::MarkovChain> copy() const override {
                    return std::make_unique<SimpleChain>(_inst, _currentState);
                }
            };
        }
    }
}


#endif /* _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_SIMPLE_CHAIN_H */
