/*
 * Created on: Jan 17, 2017
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

#ifndef _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_INFORMED_CHAIN_H
#define _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_INFORMED_CHAIN_H

#include <vector>
#include <stack>
#include "marathon/binary_matrix/interval_margin/markov_chain.h"
#include "marathon/combination_generator.h"

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * Informed Markov chain defined in
             *
             * Steffen Rechner, Linda Strowick, Matthias Müller-Hannemann.
             * Uniform sampling of bipartite graphs with degrees in prescribed intervals.
             * Journal of Complex Networks (2017). DOI: 10.1093/comnet/cnx059
             */
            class InformedChain : public MarkovChain {

            protected:

                std::vector<size_t> tmp1, tmp2, tmp3;    // auxiliary arrays

                // probability of selecting a horizontal operation
                const double q_horizontal;
                const double q_vertical;
                const Rational q_horizontal_rat;
                const Rational q_vertical_rat;

                // probability of selecting each kind of operation type
                const Rational p_trade_rat = Rational(1, 3);
                const Rational p_multishift_rat = Rational(1, 3);
                const Rational p_multiflip_rat = Rational(1, 3);

                /**
                 * Apply horizontal Trade.
                 */
                inline
                void applyTradeHorizontal(BinaryMatrix &s) {

                    /* apply horizontal trade */
                    count_trade++;

                    // randomly select two row indices
                    size_t i = rg.nextInt(_nrow);
                    size_t j = rg.nextInt(_nrow);
                    while (i == j)
                        j = rg.nextInt(_nrow);

                    // select the indices that occur in on the s_i and s_j, but not in both
                    size_t a = 0;
                    size_t b = 0;

                    // for each column position
                    for (size_t k = 0; k < _ncol; k++) {

                        bool s_ik = s.get(i, k);
                        bool s_jk = s.get(j, k);

                        if (s_ik != s_jk) {

                            // store index k in temporary array
                            tmp1[a + b] = k;

                            if (s_ik)
                                a++;
                            else
                                b++;
                        }
                    }

                    // if there is nothing to choose from
                    if (a + b == 0) {
                        loop_trade++;
                        return;
                    }

                    // randomly select a out of (a+b) elements to row i
                    rg.shuffle<size_t>(&tmp1[0], a + b);

                    // does the selection of a elements produce a loop?
                    bool allsame = true;
                    for (size_t k = 0; k < a; k++) {
                        if (!s.get(i, tmp1[k])) {
                            allsame = false;
                            break;
                        }
                    }

                    if (allsame) {
                        loop_trade++;
                        return;
                    }

                    for (size_t k = 0; k < a; k++) {
                        s.set(i, tmp1[k], true);
                        s.set(j, tmp1[k], false);
                    }

                    // the remaining elements to go row j
                    for (size_t k = a; k < a + b; k++) {
                        s.set(i, tmp1[k], false);
                        s.set(j, tmp1[k], true);
                    }
                }

                /**
                 * Apply horizontal Multi-Shift
                 * @param s
                 */
                inline void applyMultiShiftHorizontal(BinaryMatrix &s) {

                    /* apply horizontal multi-shift */
                    count_multishift++;

                    // randomly select a row index
                    size_t i = rg.nextInt(_nrow);

                    // select the positions k of a) the entries s_ik=1 that can be switched to zero
                    // and b) the entries s_ik=0 that can be switched to one
                    size_t a = 0;
                    size_t b = 0;

                    // for each column position
                    for (size_t k = 0; k < _ncol; k++) {
                        bool s_ik = s.get(i, k);

                        if (s_ik && _current_margin._colsum[k] > _inst._colsum_lower[k]) {
                            tmp1[a + b] = k;
                            a++;
                        } else if (!s_ik && _current_margin._colsum[k] < _inst._colsum_upper[k]) {
                            tmp1[a + b] = k;
                            b++;
                        }
                    }

                    // if there is nothing to choose from
                    if (a + b == 0) {
                        loop_multishift++;
                        return;
                    }

                    // select a out of (a+b) positions
                    rg.shuffle<size_t>(&tmp1[0], a + b);

                    // does the selection of a elements produce a loop?
                    bool allsame = true;
                    for (size_t k = 0; k < a; k++) {
                        if (!s.get(i, tmp1[k])) {
                            allsame = false;
                            break;
                        }
                    }

                    if (allsame) {
                        loop_multishift++;
                        return;
                    }

                    // set positive entries
                    for (size_t k = 0; k < a; k++) {
                        size_t j = tmp1[k];
                        _current_margin._colsum[j] += s.get(i, j) ? 0 : 1;
                        s.set(i, j, true);
                    }

                    // set negative entries
                    for (size_t k = a; k < a + b; k++) {
                        size_t j = tmp1[k];
                        _current_margin._colsum[j] += s.get(i, j) ? -1 : 0;
                        s.set(i, tmp1[k], false);
                    }
                }

                inline void applyTradeVertical(BinaryMatrix &s) {

                    /* apply vertical trade */

                    count_trade++;

                    // randomly select two column indices
                    size_t i = rg.nextInt(_ncol);
                    size_t j = rg.nextInt(_ncol);
                    while (i == j)
                        j = rg.nextInt(_ncol);

                    // select the indices that occur in on the s_i and s_j, but not in both
                    size_t a = 0;
                    size_t b = 0;

                    // for each row position
                    for (size_t k = 0; k < _nrow; k++) {

                        bool s_ki = s.get(k, i);
                        bool s_kj = s.get(k, j);

                        if (s_ki != s_kj) {

                            // store index k
                            tmp1[a + b] = k;

                            if (s_ki)
                                a++;
                            else
                                b++;
                        }
                    }

                    // if there is nothing to choose from
                    if (a + b == 0) {
                        loop_trade++;
                        return;
                    }

                    // randomly select a out of (a+b) positions
                    rg.shuffle<size_t>(&tmp1[0], a + b);

                    // does the selection of a elements produce a loop?
                    bool allsame = true;
                    for (size_t k = 0; k < a; k++) {
                        if (!s.get(tmp1[k], i)) {
                            allsame = false;
                            break;
                        }
                    }

                    if (allsame) {
                        loop_trade++;
                        return;
                    }

                    for (size_t k = 0; k < a; k++) {
                        s.set(tmp1[k], i, true);
                        s.set(tmp1[k], j, false);
                    }

                    // the remaining elements to go column j
                    for (size_t k = a; k < a + b; k++) {
                        s.set(tmp1[k], i, false);
                        s.set(tmp1[k], j, true);
                    }
                }

                /**
                 * Vertical Multi-Shift to s.
                 * @param s Binary Matrix.
                 */
                inline void applyMultiShiftVertical(BinaryMatrix &s) {

                    /* apply vertical multi-shift */

                    count_multishift++;

                    // randomly select a column index
                    size_t j = rg.nextInt(_ncol);

                    // select the positions i of a) the entries sij=1 that can be switched to zero
                    // and b) the entries sij=0 that can be switched to one
                    size_t a = 0;
                    size_t b = 0;

                    // for each row position
                    for (size_t i = 0; i < _nrow; i++) {
                        bool sij = s.get(i, j);

                        if (sij && _current_margin._rowsum[i] > _inst._rowsum_lower[i]) {
                            tmp1[a + b] = i;
                            a++;
                        } else if (!sij && _current_margin._rowsum[i] < _inst._rowsum_upper[i]) {
                            tmp1[a + b] = i;
                            b++;
                        }
                    }

                    // if there is nothing to choose from
                    if (a + b == 0) {
                        loop_multishift++;
                        return;
                    }

                    // select a out of (a+b) positions
                    rg.shuffle<size_t>(&tmp1[0], a + b);

                    // does the selection of a elements produce a loop?
                    bool allsame = true;
                    for (size_t k = 0; k < a; k++) {
                        size_t i = tmp1[k];
                        if (!s.get(i, j)) {
                            allsame = false;
                            break;
                        }
                    }

                    if (allsame) {
                        loop_multishift++;
                        return;
                    }

                    // set positive entries
                    for (size_t k = 0; k < a; k++) {
                        size_t i = tmp1[k];
                        _current_margin._rowsum[i] += s.get(i, j) ? 0 : 1;
                        s.set(i, j, true);
                    }

                    // set negative entries
                    for (size_t k = a; k < a + b; k++) {
                        size_t i = tmp1[k];
                        _current_margin._rowsum[i] += s.get(i, j) ? -1 : 0;
                        s.set(i, j, false);
                    }
                }


                /**
                 * Apply a horizontal Multi-Flip operation to s.
                 * @param s Binary Matrix.
                 */
                inline void applyMultiFlipHorizontal(BinaryMatrix &s) {

                    /* apply horizontal multi-flip */

                    count_multiflip++;

                    // randomly select a row index
                    size_t i = rg.nextInt(_nrow);

                    // select the positions j of a) the entries sij=1 that can be switched to zero
                    // and b) the entries sij=0 that can be switched to one
                    size_t a = 0;
                    size_t b = 0;

                    // for each column position
                    for (size_t j = 0; j < _ncol; j++) {
                        bool sij = s.get(i, j);
                        if (sij && _current_margin._colsum[j] > _inst._colsum_lower[j]) {
                            tmp1[a + b] = j;
                            a++;
                        } else if (!sij && _current_margin._colsum[j] < _inst._colsum_upper[j]) {
                            tmp1[a + b] = j;
                            b++;
                        }
                    }

                    // determine maximal number of positions to flip
                    int range = _inst._rowsum_upper[i] - _inst._rowsum_lower[i];
                    int t = std::min(range, (int) (a + b));

                    // if there is nothing to choose from
                    if (t == 0) {
                        loop_multiflip++;
                        //printf("typeA\n");
                        return;
                    }

                    // choose a random subset from tmp1 of size t
                    rg.select<size_t>(&tmp2[0], &tmp1[0], a + b, t);

                    // choose a subset of tmp2 uniformly at random
                    size_t k = rg.subset<size_t>(&tmp3[0], &tmp2[0], t);

                    // if nothing to flip
                    if (k == 0) {
                        loop_multiflip++;
                        //printf("typeB\n");
                        return;
                    }

                    // determine by how much rowsum[i] would be changed by a flip of the k selected entries
                    int diff = 0;
                    for (size_t l = 0; l < k; l++)
                        diff += s.get(i, tmp3[l]) ? -1 : 1;

                    // if entries can be flipped without violating the lower and upper bounds
                    if (_inst._rowsum_lower[i] <= _current_margin._rowsum[i] + diff &&
                        _current_margin._rowsum[i] + diff <= _inst._rowsum_upper[i]) {

                        // flip selected entries
                        for (size_t l = 0; l < k; l++) {
                            size_t j = tmp3[l];
                            _current_margin._colsum[j] += s.get(i, j) ? -1 : 1;
                            s.flip(i, j);
                        }
                        _current_margin._rowsum[i] += diff;
                    } else {
                        loop_multiflip++;
                        //printf("typeC (k=%2i, diff=%2i, range=%2i, a=%2i, b=%2i)\n", k, diff, range, a, b);
                        return;
                    }
                }


                /**
                 * Apply a vertical Multi-Flip operation to s.
                 * @param s Binary Matrix.
                 */
                inline void applyMultiFlipVertical(BinaryMatrix &s) {

                    /* apply vertical multi-flip */

                    count_multiflip++;

                    // randomly select a column index
                    size_t j = rg.nextInt(_ncol);

                    // select the positions i of a) the entries sij=1 that can be switched to zero
                    // and b) the entries sij=0 that can be switched to one
                    size_t a = 0;
                    size_t b = 0;

                    // for each column position
                    for (size_t i = 0; i < _nrow; i++) {
                        bool sij = s.get(i, j);
                        if (sij && _current_margin._rowsum[i] > _inst._rowsum_lower[i]) {
                            tmp1[a + b] = i;
                            a++;
                        } else if (!sij && _current_margin._rowsum[i] < _inst._rowsum_upper[i]) {
                            tmp1[a + b] = i;
                            b++;
                        }
                    }

                    // determine maximal number of positions to flip
                    int range = _inst._colsum_upper[j] - _inst._colsum_lower[j];
                    int t = std::min(range, (int) (a + b));

                    // if there is nothing to choose from
                    if (t == 0) {
                        loop_multiflip++;
                        //printf("typeA\n");
                        return;
                    }

                    // choose a random subset from tmp1 of size t
                    rg.select<size_t>(&tmp2[0], &tmp1[0], a + b, t);

                    // choose a subset of tmp2 uniformly at random
                    size_t k = rg.subset<size_t>(&tmp3[0], &tmp2[0], t);

                    if (k == 0) {
                        loop_multiflip++;
                        //printf("typeB\n");
                        return;
                    }

                    // determine by how much colsum[j] would be changed by a flip of the k selected entries
                    int diff = 0;
                    for (int l = 0; l < k; l++)
                        diff += s.get(tmp3[l], j) ? -1 : 1;

                    // if entries can be flipped without violating the lower and upper bounds
                    if (_inst._colsum_lower[j] <= _current_margin._colsum[j] + diff &&
                        _current_margin._colsum[j] + diff <= _inst._colsum_upper[j]) {

                        // flip selected entries
                        for (size_t l = 0; l < k; l++) {
                            size_t i = tmp3[l];
                            _current_margin._rowsum[i] += s.get(i, j) ? -1 : 1;
                            s.flip(i, j);
                        }
                        _current_margin._colsum[j] += diff;
                    } else {
                        loop_multiflip++;
                        //printf("typeC (k=%2i, diff=%2i, range=%2i, a=%2i, b=%2i)\n", k, diff, range, a, b);
                        return;
                    }
                }


                /**
                 * Simulate horizontal trade operations.
                 * @param s
                 * @param Rational
                 * @return
                 */
                void simulateTradeHorizontal(
                        const BinaryMatrix &s,
                        const std::function<void(const State &, const Rational &)> &process,
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) const {

                    /* simulate horizontal trade */

                    Rational loop(0);   // loop probability

                    // temporary memory
                    std::vector<size_t> tmp1(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp2(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp3(std::max(_nrow, _ncol));

                    // randomly select two row indices
                    for (size_t i = 0; i < _nrow; i++) {
                        for (size_t j = i + 1; j < _nrow; j++) {

                            // select the indices that occur in on the s_i and s_j, but not in both
                            size_t a = 0;
                            size_t b = 0;

                            // for each column position
                            for (size_t k = 0; k < _ncol; k++) {

                                bool s_ik = s.get(i, k);
                                bool s_jk = s.get(j, k);

                                if (s_ik != s_jk) {

                                    // store index k at X array
                                    tmp1[a + b] = k;

                                    if (s_ik)
                                        a++;
                                    else
                                        b++;
                                }
                            }

                            // if there is nothing to choose from
                            if (a + b == 0) {
                                loop += q_horizontal_rat * p_trade_rat *
                                        Rational(2, _nrow * (_nrow - 1));
                                continue;
                            }

                            // simulate randomly select a out of (a+b) elements
                            CombinationGenerator<size_t> cg(&tmp1[0], &tmp2[0], a + b, a);

                            // the probability for each choice of the loop
                            const Rational p =
                                    (q_horizontal_rat * p_trade_rat *
                                     Rational(2, _nrow * (_nrow - 1))) / binom(a + b, a);
                            do {

                                // create copy of s
                                BinaryMatrix s2(s);

                                for (int k = 0; k < a + b; k++) {
                                    s2.set(i, tmp1[k], false);
                                    s2.set(j, tmp1[k], true);
                                }

                                for (int k = 0; k < a; k++) {
                                    s2.set(i, tmp2[k], true);
                                    s2.set(j, tmp2[k], false);
                                }

                                process(s2, p);

                            } while (cg.next());
                        }
                    }

                    // process loop
                    process(s, loop);
                }


                /**
                 * Simulate horizontal multi-shift operations.
                 * @param s
                 * @param Rational
                 * @return
                 */
                void simulateMultiShiftHorizontal(
                        const BinaryMatrix &s,
                        const std::function<void(const State &, const Rational &)> &process,
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) const {

                    /* simulate horizontal multi-shift */

                    Rational loop(0);   // loop probability

                    // temporary memory
                    std::vector<size_t> tmp1(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp2(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp3(std::max(_nrow, _ncol));

                    // randomly select a row index
                    for (size_t i = 0; i < _nrow; i++) {

                        // select the positions k of a) the entries s_ik=1 that can be switched to zero
                        // and b) the entries s_ik=0 that can be switched to one
                        size_t a = 0;
                        size_t b = 0;

                        // for each column position
                        for (size_t k = 0; k < _ncol; k++) {
                            bool s_ik = s.get(i, k);

                            if (s_ik && colsum[k] > _inst._colsum_lower[k]) {
                                tmp1[a + b] = k;
                                a++;
                            } else if (!s_ik && colsum[k] < _inst._colsum_upper[k]) {
                                tmp1[a + b] = k;
                                b++;
                            }
                        }

                        // if there is nothing to choose from
                        if (a + b == 0) {
                            loop += q_horizontal_rat * p_multishift_rat * Rational(1, _nrow);
                            continue;
                        }

                        // simulate randomly select a out of (a+b) elements
                        CombinationGenerator<size_t> cg(&tmp1[0], &tmp2[0], a + b, a);

                        // the probability for each choice of the loop
                        const Rational p = q_horizontal_rat * p_multishift_rat * Rational(1, _nrow) /
                                           binom(a + b, a);

                        do {
                            // create copy of s
                            BinaryMatrix s2(s);

                            for (size_t k = 0; k < a + b; k++)
                                s2.set(i, tmp1[k], false);

                            for (size_t k = 0; k < a; k++)
                                s2.set(i, tmp2[k], true);

                            process(s2, p);

                        } while (cg.next());

                    }

                    // process loop
                    process(s, loop);
                }

                /**
                 * Simulate vertical trade operations.
                 * @param s
                 * @param Rational
                 * @return
                 */
                void simulateTradeVertical(
                        const BinaryMatrix &s,
                        const std::function<void(const State &, const Rational &)> &process,
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) const {

                    /* simulate vertical trade */

                    Rational loop(0);   // loop probability

                    // temporary memory
                    std::vector<size_t> tmp1(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp2(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp3(std::max(_nrow, _ncol));

                    // randomly select two column indices
                    for (size_t i = 0; i < _ncol; i++) {
                        for (size_t j = i + 1; j < _ncol; j++) {

                            // select the indices that occur in on the s_i and s_j, but not in both
                            size_t a = 0;
                            size_t b = 0;

                            // for each row position
                            for (size_t k = 0; k < _nrow; k++) {

                                bool s_ki = s.get(k, i);
                                bool s_kj = s.get(k, j);

                                if (s_ki != s_kj) {

                                    // store index k at X array
                                    tmp1[a + b] = k;

                                    if (s_ki)
                                        a++;
                                    else
                                        b++;
                                }
                            }

                            // if there is nothing to choose from
                            if (a + b == 0) {
                                loop += q_vertical_rat * p_trade_rat * Rational(2, _ncol * (_ncol - 1));
                                continue;
                            }

                            // simulate randomly select a out of (a+b) elements
                            CombinationGenerator<size_t> cg(&tmp1[0], &tmp2[0], a + b, a);

                            // the probability for each choice of the loop
                            const Rational p =
                                    q_vertical_rat * p_trade_rat * Rational(2, _ncol * (_ncol - 1)) /
                                    binom(a + b, a);
                            do {
                                // create copy of s
                                BinaryMatrix s2(s);

                                for (size_t k = 0; k < a + b; k++) {
                                    s2.set(tmp1[k], i, false);
                                    s2.set(tmp1[k], j, true);
                                }

                                for (size_t k = 0; k < a; k++) {
                                    s2.set(tmp2[k], i, true);
                                    s2.set(tmp2[k], j, false);
                                }

                                process(s2, p);

                            } while (cg.next());
                        }
                    }

                    // process loop
                    process(s, loop);
                }

                /**
                 * Simulate vertical multi-shift operation.
                 * @param s
                 * @param Rational
                 * @return
                 */
                void simulateMultiShiftVertical(
                        const BinaryMatrix &s,
                        const std::function<void(const State &, const Rational &)> &process,
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) const {

                    /* simulate vertical multi-shift */

                    Rational loop(0);   // loop probability

                    // temporary memory
                    std::vector<size_t> tmp1(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp2(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp3(std::max(_nrow, _ncol));

                    // randomly select a column index
                    for (size_t j = 0; j < _ncol; j++) {

                        // select the positions i of a) the entries s_ij=1 that can be switched to zero
                        // and b) the entries s_ij=0 that can be switched to one
                        size_t a = 0;
                        size_t b = 0;

                        // for each row position
                        for (size_t i = 0; i < _nrow; i++) {
                            bool sij = s.get(i, j);

                            if (sij && rowsum[i] > _inst._rowsum_lower[i]) {
                                tmp1[a + b] = i;
                                a++;
                            } else if (!sij && rowsum[i] < _inst._rowsum_upper[i]) {
                                tmp1[a + b] = i;
                                b++;
                            }
                        }

                        // if there is nothing to choose from
                        if (a + b == 0) {
                            loop += q_vertical_rat * p_multishift_rat * Rational(1, _ncol);
                            continue;
                        }

                        // simulate randomly select a out of (a+b) elements
                        CombinationGenerator<size_t> cg(&tmp1[0], &tmp2[0], a + b, a);

                        // the probability for each choice of the loop
                        const Rational p = q_vertical_rat * p_multishift_rat * Rational(1, _ncol) / binom(a + b, a);
                        do {
                            // create copy of s
                            BinaryMatrix s2(s);

                            for (size_t k = 0; k < a + b; k++)
                                s2.set(tmp1[k], j, false);

                            for (size_t k = 0; k < a; k++)
                                s2.set(tmp2[k], j, true);

                            process(s2, p);

                        } while (cg.next());
                    }

                    // process loop
                    process(s, loop);
                }


                /**
                 * Simulate horizontal multi-flip operation.
                 * @param s
                 * @param Rational
                 * @return
                 */
                void simulateMultiFlipHorizontal(
                        const BinaryMatrix &s,
                        const std::function<void(const State &, const Rational &)> &process,
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) const {

                    /* simulate horizontal multi-flip */

                    Rational loop(0);   // loop probability

                    // temporary memory
                    std::vector<size_t> tmp1(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp2(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp3(std::max(_nrow, _ncol));

                    // randomly select a row index
                    for (size_t i = 0; i < _nrow; i++) {

                        // select the positions j of a) the entries sij=1 that can be switched to zero
                        // and b) the entries sij=0 that can be switched to one
                        size_t a = 0;
                        size_t b = 0;

                        // for each column position
                        for (size_t j = 0; j < _ncol; j++) {
                            bool sij = s.get(i, j);
                            if (sij && colsum[j] > _inst._colsum_lower[j]) {
                                tmp1[a + b] = j;
                                a++;
                            } else if (!sij && colsum[j] < _inst._colsum_upper[j]) {
                                tmp1[a + b] = j;
                                b++;
                            }
                        }

                        // determine maximal number of positions to flip
                        int range = _inst._rowsum_upper[i] - _inst._rowsum_lower[i];
                        int t = std::min(range, (int) (a + b));

                        // if there is nothing to choose from
                        if (t == 0) {
                            loop += q_horizontal_rat * p_multiflip_rat * Rational(1, _nrow);
                            continue;
                        }

                        // simulate randomly selecting t out of (a+b) elements
                        CombinationGenerator<size_t> cg(&tmp1[0], &tmp2[0], a + b, t);
                        do {

                            // simulate selecting random subsets from tmp2

                            // each of the 2^t subsets has the same probability p of being selected
                            const Rational p =
                                    q_horizontal_rat * p_multiflip_rat * Rational(1, _nrow) /
                                    (binom(a + b, t) * pow(2, t));

                            // simulate selecting subset of size k
                            for (size_t k = 0; k <= t; k++) {

                                // simulate selecting k objects out of t
                                CombinationGenerator<size_t> cg_sub(&tmp2[0], &tmp3[0], t, k);
                                do {

                                    // determine by how much rowsum[i] would be changed by a flip of the k selected entries
                                    int diff = 0;
                                    for (size_t l = 0; l < k; l++)
                                        diff += s.get(i, tmp3[l]) ? -1 : 1;

                                    // if entries can be flipped without violating the lower and upper bounds
                                    if (_inst._rowsum_lower[i] <= rowsum[i] + diff &&
                                        rowsum[i] + diff <= _inst._rowsum_upper[i]) {

                                        // create copy of s
                                        BinaryMatrix s2(s);

                                        // flip selected entries
                                        for (size_t l = 0; l < k; l++) {
                                            size_t j = tmp3[l];
                                            s2.flip(i, j);
                                        }

                                        process(s2, p);

                                    } else {
                                        loop += p;
                                    }
                                } while (cg_sub.next());
                            }

                        } while (cg.next());
                    }

                    // process loop
                    process(s, loop);
                }


                /**
                 * Simulate vertical multi-flip operation.
                 * @param s
                 * @param Rational
                 * @return
                 */
                void simulateMultiFlipVertical(
                        const BinaryMatrix &s,
                        const std::function<void(const State &, const Rational &)> &process,
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) const {

                    /* simulate vertical multi-fliop */

                    Rational loop(0);   // loop probability

                    // temporary memory
                    std::vector<size_t> tmp1(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp2(std::max(_nrow, _ncol));
                    std::vector<size_t> tmp3(std::max(_nrow, _ncol));

                    // randomly select a column index
                    for (size_t j = 0; j < _ncol; j++) {

                        // select the positions i of a) the entries sij=1 that can be switched to zero
                        // and b) the entries sij=0 that can be switched to one
                        size_t a = 0;
                        size_t b = 0;

                        // for each row position
                        for (size_t i = 0; i < _nrow; i++) {
                            bool sij = s.get(i, j);
                            if (sij && rowsum[i] > _inst._rowsum_lower[i]) {
                                tmp1[a + b] = i;
                                a++;
                            } else if (!sij && rowsum[i] < _inst._rowsum_upper[i]) {
                                tmp1[a + b] = i;
                                b++;
                            }
                        }

                        // determine maximal number of positions to flip
                        int range = _inst._colsum_upper[j] - _inst._colsum_lower[j];
                        int t = std::min(range, (int) (a + b));

                        // if there is nothing to choose from
                        if (t == 0) {
                            loop += q_vertical_rat * p_multiflip_rat * Rational(1, _ncol);
                            continue;
                        }

                        // simulate randomly selecting t out of (a+b) elements
                        CombinationGenerator<size_t> cg(&tmp1[0], &tmp2[0], a + b, t);
                        do {

                            // simulate selecting random subsets from tmp2

                            // each of the 2^t subsets has the same probability p of being selected
                            const Rational p =
                                    q_vertical_rat * p_multiflip_rat * Rational(1, _ncol) /
                                    (binom(a + b, t) * pow(2, t));

                            // simulate selecting subset of size k
                            for (size_t k = 0; k <= t; k++) {

                                // simulate selecting k objects out of t
                                CombinationGenerator<size_t> cg_sub(&tmp2[0], &tmp3[0], t, k);
                                do {

                                    // determine by how much colsum[j] would be changed by a flip of the k selected entries
                                    int diff = 0;
                                    for (size_t l = 0; l < k; l++) {
                                        size_t i = tmp3[l];
                                        diff += s.get(i, j) ? -1 : 1;
                                    }

                                    // if entries can be flipped without violating the lower and upper bounds
                                    if (_inst._colsum_lower[j] <= colsum[j] + diff &&
                                        colsum[j] + diff <= _inst._colsum_upper[j]) {

                                        // create copy of s
                                        BinaryMatrix s2(s);

                                        // flip selected entries
                                        for (size_t l = 0; l < k; l++) {
                                            size_t i = tmp3[l];
                                            s2.flip(i, j);
                                        }

                                        process(s2, p);

                                    } else {
                                        loop += p;
                                    }
                                } while (cg_sub.next());
                            }

                        } while (cg.next());
                    }

                    // process loop
                    process(s, loop);
                }

                void init() {
                    tmp1.resize(std::max(_nrow, _ncol));
                    tmp2.resize(std::max(_nrow, _ncol));
                    tmp3.resize(std::max(_nrow, _ncol));
                    p_trade = p_trade_rat.convert_to<double>();
                    p_multishift = p_multishift_rat.convert_to<double>();
                    p_multiflip = p_multiflip_rat.convert_to<double>();
                }

            public:

                // convert to double for faster comparison
                double p_trade;
                double p_multishift;
                double p_multiflip;

                // number of every operation
                size_t count_trade = 0;
                size_t count_multiflip = 0;
                size_t count_multishift = 0;

                // loop count of every operation
                size_t loop_trade = 0;
                size_t loop_multiflip = 0;
                size_t loop_multishift = 0;

                /**
                 * Create a Markov chain.
                 * @param inst Lower and upper bounds on the row and column sums.
                 */
                explicit InformedChain(Instance inst) :
                        MarkovChain(std::move(inst)),
                        q_horizontal(_ncol / (double) (_nrow + _ncol)),
                        q_vertical(_nrow / (double) (_nrow + _ncol)),
                        q_horizontal_rat(Rational(_ncol, _nrow + _ncol)),
                        q_vertical_rat(Rational(_nrow, _nrow + _ncol)) {
                    init();
                }

                /**
                 * Create a Markov chain. Use a certain binary matrix as initial state.
                 * @param inst Lower and upper bounds on the row and column sums.
                 * @param bin Binary matrix.
                 */
                InformedChain(Instance inst, BinaryMatrix bin)
                        : MarkovChain(std::move(inst), std::move(bin)),
                          q_horizontal(_ncol / (double) (_nrow + _ncol)),
                          q_vertical(_nrow / (double) (_nrow + _ncol)),
                          q_horizontal_rat(Rational(_ncol, _nrow + _ncol)),
                          q_vertical_rat(Rational(_nrow, _nrow + _ncol)) {
                    init();
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

                    // create a copy of x
                    BinaryMatrix s(static_cast<const BinaryMatrix &>(x));

                    // create temporary array of row and column sums
                    std::vector<int> rowsum(_nrow);
                    std::vector<int> colsum(_ncol);
                    for (size_t i = 0; i < _nrow; i++) {
                        for (size_t j = 0; j < _ncol; j++) {
                            if (s.get(i, j)) {
                                rowsum[i]++;
                                colsum[j]++;
                            }
                        }
                    }

                    simulateTradeHorizontal(s, process, rowsum, colsum);
                    simulateTradeVertical(s, process, rowsum, colsum);
                    simulateMultiShiftHorizontal(s, process, rowsum, colsum);
                    simulateMultiShiftVertical(s, process, rowsum, colsum);
                    simulateMultiFlipHorizontal(s, process, rowsum, colsum);
                    simulateMultiFlipVertical(s, process, rowsum, colsum);
                }

                /**
                 * Randomize the current state of the Markov chain.
                 */
                virtual void step() override {

                    // select two random numbers from [0,1)
                    const double p = rg.nextDouble();
                    const double q = rg.nextDouble();

                    // which operation will be applied?
                    if (p < p_trade) {                                   // apply trade

                        if (q < q_horizontal) {
                            applyTradeHorizontal(_currentState);
                        } else if (q < q_horizontal + q_vertical) {
                            applyTradeVertical(_currentState);
                        } else {
                            throw std::runtime_error("Error while applying trade.");
                        }

                    } else if (p < p_trade + p_multishift) {      // apply multi-shift

                        if (q < q_horizontal) {
                            applyMultiShiftHorizontal(_currentState);
                        } else if (q < q_horizontal + q_vertical) {
                            applyMultiShiftVertical(_currentState);
                        } else {
                            throw std::runtime_error("Error while applying multi-shift.");
                        }

                    } else if (p < p_trade + p_multishift +
                                   p_multiflip) {                                                       // apply multi-flip

                        if (q < q_horizontal) {
                            applyMultiFlipHorizontal(_currentState);
                        } else if (q < q_horizontal + q_vertical) {
                            applyMultiFlipVertical(_currentState);
                        } else {
                            throw std::runtime_error("Error while applying multi-flip.");
                        }
                    } else {
                        throw std::runtime_error("Error while applying operation.");
                    }
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual std::unique_ptr<marathon::MarkovChain> copy() const override {
                    return std::make_unique<InformedChain>(_inst, _currentState);
                }

            };
        }
    }
}

#endif /* _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_INFORMED_CHAIN_H */
