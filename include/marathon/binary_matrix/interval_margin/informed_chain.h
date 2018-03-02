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

                // temporary memory
                int *tmp1;
                int *tmp2;
                int *tmp3;

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
                void applyTradeHorizontal(BinaryMatrix *s) {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    /* apply horizontal trade */
                    count_trade++;

                    // randomly select two row indices
                    int i = rg.nextInt(nrow);
                    int j = rg.nextInt(nrow);
                    while (i == j)
                        j = rg.nextInt(nrow);

                    // select the indices that occur in on the s_i and s_j, but not in both
                    int a = 0;
                    int b = 0;

                    // for each column position
                    for (int k = 0; k < ncol; k++) {

                        bool s_ik = s->get(i, k);
                        bool s_jk = s->get(j, k);

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
                    rg.shuffle<int>(tmp1, a + b);

                    // does the selection of a elements produce a loop?
                    bool allsame = true;
                    for (int k = 0; k < a; k++) {
                        if (!s->get(i, tmp1[k])) {
                            allsame = false;
                            break;
                        }
                    }

                    if (allsame) {
                        loop_trade++;
                        return;
                    }

                    for (int k = 0; k < a; k++) {
                        s->set(i, tmp1[k], true);
                        s->set(j, tmp1[k], false);
                    }

                    // the remaining elements to go row j
                    for (int k = a; k < a + b; k++) {
                        s->set(i, tmp1[k], false);
                        s->set(j, tmp1[k], true);
                    }
                }

                /**
                 * Apply horizontal Multi-Shift
                 * @param s
                 */
                inline void applyMultiShiftHorizontal(BinaryMatrix *s) {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    /* apply horizontal multi-shift */
                    count_multishift++;

                    // randomly select a row index
                    int i = rg.nextInt(nrow);

                    // select the positions k of a) the entries s_ik=1 that can be switched to zero
                    // and b) the entries s_ik=0 that can be switched to one
                    int a = 0;
                    int b = 0;

                    // for each column position
                    for (int k = 0; k < ncol; k++) {
                        bool s_ik = s->get(i, k);

                        if (s_ik && current_margin.colsum[k] > inst.colsum_lower[k]) {
                            tmp1[a + b] = k;
                            a++;
                        } else if (!s_ik && current_margin.colsum[k] < inst.colsum_upper[k]) {
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
                    rg.shuffle<int>(tmp1, a + b);

                    // does the selection of a elements produce a loop?
                    bool allsame = true;
                    for (int k = 0; k < a; k++) {
                        if (!s->get(i, tmp1[k])) {
                            allsame = false;
                            break;
                        }
                    }

                    if (allsame) {
                        loop_multishift++;
                        return;
                    }

                    // set positive entries
                    for (int k = 0; k < a; k++) {
                        int j = tmp1[k];
                        current_margin.colsum[j] += s->get(i, j) ? 0 : 1;
                        s->set(i, j, true);
                    }

                    // set negative entries
                    for (int k = a; k < a + b; k++) {
                        int j = tmp1[k];
                        current_margin.colsum[j] += s->get(i, j) ? -1 : 0;
                        s->set(i, tmp1[k], false);
                    }
                }

                inline void applyTradeVertical(BinaryMatrix *s) {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    /* apply vertical trade */

                    count_trade++;

                    // randomly select two column indices
                    int i = rg.nextInt(ncol);
                    int j = rg.nextInt(ncol);
                    while (i == j)
                        j = rg.nextInt(ncol);

                    // select the indices that occur in on the s_i and s_j, but not in both
                    int a = 0;
                    int b = 0;

                    // for each row position
                    for (int k = 0; k < nrow; k++) {

                        bool s_ki = s->get(k, i);
                        bool s_kj = s->get(k, j);

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
                    rg.shuffle<int>(tmp1, a + b);

                    // does the selection of a elements produce a loop?
                    bool allsame = true;
                    for (int k = 0; k < a; k++) {
                        if (!s->get(tmp1[k], i)) {
                            allsame = false;
                            break;
                        }
                    }

                    if (allsame) {
                        loop_trade++;
                        return;
                    }


                    for (int k = 0; k < a; k++) {
                        s->set(tmp1[k], i, true);
                        s->set(tmp1[k], j, false);
                    }

                    // the remaining elements to go column j
                    for (int k = a; k < a + b; k++) {
                        s->set(tmp1[k], i, false);
                        s->set(tmp1[k], j, true);
                    }
                }

                /**
                 * Vertical Multi-Shift to s.
                 * @param s Binary Matrix.
                 */
                inline void applyMultiShiftVertical(BinaryMatrix *s) {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    /* apply vertical multi-shift */

                    count_multishift++;

                    // randomly select a column index
                    int j = rg.nextInt(ncol);

                    // select the positions i of a) the entries sij=1 that can be switched to zero
                    // and b) the entries sij=0 that can be switched to one
                    int a = 0;
                    int b = 0;

                    // for each row position
                    for (int i = 0; i < nrow; i++) {
                        bool sij = s->get(i, j);

                        if (sij && current_margin.rowsum[i] > inst.rowsum_lower[i]) {
                            tmp1[a + b] = i;
                            a++;
                        } else if (!sij && current_margin.rowsum[i] < inst.rowsum_upper[i]) {
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
                    rg.shuffle<int>(tmp1, a + b);

                    // does the selection of a elements produce a loop?
                    bool allsame = true;
                    for (int k = 0; k < a; k++) {
                        int i = tmp1[k];
                        if (!s->get(i, j)) {
                            allsame = false;
                            break;
                        }
                    }

                    if (allsame) {
                        loop_multishift++;
                        return;
                    }

                    // set positive entries
                    for (int k = 0; k < a; k++) {
                        int i = tmp1[k];
                        current_margin.rowsum[i] += s->get(i, j) ? 0 : 1;
                        s->set(i, j, true);
                    }

                    // set negative entries
                    for (int k = a; k < a + b; k++) {
                        int i = tmp1[k];
                        current_margin.rowsum[i] += s->get(i, j) ? -1 : 0;
                        s->set(i, j, false);
                    }
                }


                /**
                 * Apply a horizontal Multi-Flip operation to s.
                 * @param s Binary Matrix.
                 */
                inline void applyMultiFlipHorizontal(BinaryMatrix *s) {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    /* apply horizontal multi-flip */

                    count_multiflip++;

                    // randomly select a row index
                    int i = rg.nextInt(nrow);

                    // select the positions j of a) the entries sij=1 that can be switched to zero
                    // and b) the entries sij=0 that can be switched to one
                    int a = 0;
                    int b = 0;

                    // for each column position
                    for (int j = 0; j < ncol; j++) {
                        bool sij = s->get(i, j);
                        if (sij && current_margin.colsum[j] > inst.colsum_lower[j]) {
                            tmp1[a + b] = j;
                            a++;
                        } else if (!sij && current_margin.colsum[j] < inst.colsum_upper[j]) {
                            tmp1[a + b] = j;
                            b++;
                        }
                    }

                    // determine maximal number of positions to flip
                    int range = inst.rowsum_upper[i] - inst.rowsum_lower[i];
                    int t = std::min(range, a + b);

                    // if there is nothing to choose from
                    if (t == 0) {
                        loop_multiflip++;
                        //printf("typeA\n");
                        return;
                    }

                    // choose a random subset from tmp1 of size t
                    rg.select<int>(tmp2, tmp1, a + b, t);

                    // choose a subset of tmp2 uniformly at random
                    int k = rg.subset<int>(tmp3, tmp2, t);

                    // if nothing to flip
                    if(k == 0) {
                        loop_multiflip++;
                        //printf("typeB\n");
                        return;
                    }

                    // determine by how much rowsum[i] would be changed by a flip of the k selected entries
                    int diff = 0;
                    for (int l = 0; l < k; l++)
                        diff += s->get(i, tmp3[l]) ? -1 : 1;

                    // if entries can be flipped without violating the lower and upper bounds
                    if (inst.rowsum_lower[i] <= current_margin.rowsum[i] + diff &&
                        current_margin.rowsum[i] + diff <= inst.rowsum_upper[i]) {

                        // flip selected entries
                        for (int l = 0; l < k; l++) {
                            int j = tmp3[l];
                            current_margin.colsum[j] += s->get(i, j) ? -1 : 1;
                            s->flip(i, j);
                        }
                        current_margin.rowsum[i] += diff;
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
                inline void applyMultiFlipVertical(BinaryMatrix *s) {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    /* apply vertical multi-flip */

                    count_multiflip++;

                    // randomly select a column index
                    int j = rg.nextInt(ncol);

                    // select the positions i of a) the entries sij=1 that can be switched to zero
                    // and b) the entries sij=0 that can be switched to one
                    int a = 0;
                    int b = 0;

                    // for each column position
                    for (int i = 0; i < nrow; i++) {
                        bool sij = s->get(i, j);
                        if (sij && current_margin.rowsum[i] > inst.rowsum_lower[i]) {
                            tmp1[a + b] = i;
                            a++;
                        } else if (!sij && current_margin.rowsum[i] < inst.rowsum_upper[i]) {
                            tmp1[a + b] = i;
                            b++;
                        }
                    }

                    // determine maximal number of positions to flip
                    int range = inst.colsum_upper[j] - inst.colsum_lower[j];
                    int t = std::min(range, a + b);

                    // if there is nothing to choose from
                    if (t == 0) {
                        loop_multiflip++;
                        //printf("typeA\n");
                        return;
                    }

                    // choose a random subset from tmp1 of size t
                    rg.select<int>(tmp2, tmp1, a + b, t);

                    // choose a subset of tmp2 uniformly at random
                    int k = rg.subset<int>(tmp3, tmp2, t);

                    if(k == 0) {
                        loop_multiflip++;
                        //printf("typeB\n");
                        return;
                    }

                    // determine by how much colsum[j] would be changed by a flip of the k selected entries
                    int diff = 0;
                    for (int l = 0; l < k; l++)
                        diff += s->get(tmp3[l], j) ? -1 : 1;

                    // if entries can be flipped without violating the lower and upper bounds
                    if (inst.colsum_lower[j] <= current_margin.colsum[j] + diff &&
                        current_margin.colsum[j] + diff <= inst.colsum_upper[j]) {

                        // flip selected entries
                        for (int l = 0; l < k; l++) {
                            int i = tmp3[l];
                            current_margin.rowsum[i] += s->get(i, j) ? -1 : 1;
                            s->flip(i, j);
                        }
                        current_margin.colsum[j] += diff;
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
                        const BinaryMatrix *s,
                        const std::function<void(const State *, const Rational &)> &process,
                        const int *rowsum,
                        const int *colsum
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    Rational loop(0);

                    /* simulate horizontal trade */

                    // randomly select two row indices
                    for (int i = 0; i < nrow; i++) {
                        for (int j = i + 1; j < nrow; j++) {

                            // select the indices that occur in on the s_i and s_j, but not in both
                            int a = 0;
                            int b = 0;

                            // for each column position
                            for (int k = 0; k < ncol; k++) {

                                bool s_ik = s->get(i, k);
                                bool s_jk = s->get(j, k);

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
                                        Rational(2, nrow * (nrow - 1));
                                continue;
                            }

                            // simulate randomly select a out of (a+b) elements
                            CombinationGenerator<int> cg(tmp1, tmp2, a + b, a);

                            // the probability for each choice of the loop
                            const Rational p =
                                    (q_horizontal_rat * p_trade_rat *
                                    Rational(2, nrow * (nrow - 1))) / binom(a + b, a);
                            do {

                                // create copy of s
                                BinaryMatrix s2(*s);

                                for (int k = 0; k < a + b; k++) {
                                    s2.set(i, tmp1[k], false);
                                    s2.set(j, tmp1[k], true);
                                }

                                for (int k = 0; k < a; k++) {
                                    s2.set(i, tmp2[k], true);
                                    s2.set(j, tmp2[k], false);
                                }

                                process(&s2, p);

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
                        const BinaryMatrix *s,
                        const std::function<void(const State *, const Rational &)> &process,
                        const int *rowsum,
                        const int *colsum
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    Rational loop(0);

                    /* simulate horizontal multi-shift */

                    // randomly select a row index
                    for (int i = 0; i < nrow; i++) {

                        // select the positions k of a) the entries s_ik=1 that can be switched to zero
                        // and b) the entries s_ik=0 that can be switched to one
                        int a = 0;
                        int b = 0;

                        // for each column position
                        for (int k = 0; k < ncol; k++) {
                            bool s_ik = s->get(i, k);

                            if (s_ik && colsum[k] > inst.colsum_lower[k]) {
                                tmp1[a + b] = k;
                                a++;
                            } else if (!s_ik && colsum[k] < inst.colsum_upper[k]) {
                                tmp1[a + b] = k;
                                b++;
                            }
                        }

                        // if there is nothing to choose from
                        if (a + b == 0) {
                            loop += q_horizontal_rat * p_multishift_rat * Rational(1, nrow);
                            continue;
                        }

                        // simulate randomly select a out of (a+b) elements
                        CombinationGenerator<int> cg(tmp1, tmp2, a + b, a);

                        // the probability for each choice of the loop
                        const Rational p = q_horizontal_rat * p_multishift_rat * Rational(1, nrow) /
                                           binom(a + b, a);

                        do {
                            // create copy of s
                            BinaryMatrix s2(*s);

                            for (int k = 0; k < a + b; k++)
                                s2.set(i, tmp1[k], false);

                            for (int k = 0; k < a; k++)
                                s2.set(i, tmp2[k], true);

                            process(&s2, p);

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
                        const BinaryMatrix *s,
                        const std::function<void(const State *, const Rational &)> &process,
                        const int *rowsum,
                        const int *colsum
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    Rational loop(0);

                    /* simulate vertical trade */

                    // randomly select two column indices
                    for (int i = 0; i < ncol; i++) {
                        for (int j = i + 1; j < ncol; j++) {

                            // select the indices that occur in on the s_i and s_j, but not in both
                            int a = 0;
                            int b = 0;

                            // for each row position
                            for (int k = 0; k < nrow; k++) {

                                bool s_ki = s->get(k, i);
                                bool s_kj = s->get(k, j);

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
                                loop += q_vertical_rat * p_trade_rat *
                                        Rational(2, ncol * (ncol - 1));
                                continue;
                            }

                            // simulate randomly select a out of (a+b) elements
                            CombinationGenerator<int> cg(tmp1, tmp2, a + b, a);

                            // the probability for each choice of the loop
                            const Rational p =
                                    q_vertical_rat * p_trade_rat * Rational(2, ncol * (ncol - 1)) /
                                    binom(a + b, a);
                            do {
                                // create copy of s
                                BinaryMatrix s2(*s);

                                for (int k = 0; k < a + b; k++) {
                                    s2.set(tmp1[k], i, false);
                                    s2.set(tmp1[k], j, true);
                                }

                                for (int k = 0; k < a; k++) {
                                    s2.set(tmp2[k], i, true);
                                    s2.set(tmp2[k], j, false);
                                }

                                process(&s2, p);

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
                        const BinaryMatrix *s,
                        const std::function<void(const State *, const Rational &)> &process,
                        const int *rowsum,
                        const int *colsum
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    Rational loop(0);

                    /* simulate vertical multi-shift */

                    // randomly select a column index
                    for (int j = 0; j < ncol; j++) {

                        // select the positions i of a) the entries s_ij=1 that can be switched to zero
                        // and b) the entries s_ij=0 that can be switched to one
                        int a = 0;
                        int b = 0;

                        // for each row position
                        for (int i = 0; i < nrow; i++) {
                            bool sij = s->get(i, j);

                            if (sij && rowsum[i] > inst.rowsum_lower[i]) {
                                tmp1[a + b] = i;
                                a++;
                            } else if (!sij && rowsum[i] < inst.rowsum_upper[i]) {
                                tmp1[a + b] = i;
                                b++;
                            }
                        }

                        // if there is nothing to choose from
                        if (a + b == 0) {
                            loop += q_vertical_rat * p_multishift_rat * Rational(1, ncol);
                            continue;
                        }

                        // simulate randomly select a out of (a+b) elements
                        CombinationGenerator<int> cg(tmp1, tmp2, a + b, a);

                        // the probability for each choice of the loop
                        const Rational p = q_vertical_rat * p_multishift_rat * Rational(1, ncol) /
                                           binom(a + b, a);
                        do {
                            // create copy of s
                            BinaryMatrix s2(*s);

                            for (int k = 0; k < a + b; k++)
                                s2.set(tmp1[k], j, false);

                            for (int k = 0; k < a; k++)
                                s2.set(tmp2[k], j, true);

                            process(&s2, p);

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
                        const BinaryMatrix *s,
                        const std::function<void(const State *, const Rational &)> &process,
                        const int *rowsum,
                        const int *colsum
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    Rational loop(0);

                    /* simulate horizontal multi-flip */

                    // randomly select a row index
                    for (int i = 0; i < nrow; i++) {

                        // select the positions j of a) the entries sij=1 that can be switched to zero
                        // and b) the entries sij=0 that can be switched to one
                        int a = 0;
                        int b = 0;

                        // for each column position
                        for (int j = 0; j < ncol; j++) {
                            bool sij = s->get(i, j);
                            if (sij && colsum[j] > inst.colsum_lower[j]) {
                                tmp1[a + b] = j;
                                a++;
                            } else if (!sij && colsum[j] < inst.colsum_upper[j]) {
                                tmp1[a + b] = j;
                                b++;
                            }
                        }

                        // determine maximal number of positions to flip
                        int range = inst.rowsum_upper[i] - inst.rowsum_lower[i];
                        int t = std::min(range, a + b);

                        // if there is nothing to choose from
                        if (t == 0) {
                            loop += q_horizontal_rat * p_multiflip_rat * Rational(1, nrow);
                            continue;
                        }

                        // simulate randomly selecting t out of (a+b) elements
                        CombinationGenerator<int> cg(tmp1, tmp2, a + b, t);
                        do {

                            // simulate selecting random subsets from tmp2

                            // each of the 2^t subsets has the same probability p of being selected
                            const Rational p =
                                    q_horizontal_rat * p_multiflip_rat * Rational(1, nrow) /
                                    (binom(a + b, t) * pow(2, t));

                            // simulate selecting subset of size k
                            for (int k = 0; k <= t; k++) {

                                // simulate selecting k objects out of t
                                CombinationGenerator<int> cg_sub(tmp2, tmp3, t, k);
                                do {

                                    // determine by how much rowsum[i] would be changed by a flip of the k selected entries
                                    int diff = 0;
                                    for (int l = 0; l < k; l++)
                                        diff += s->get(i, tmp3[l]) ? -1 : 1;

                                    // if entries can be flipped without violating the lower and upper bounds
                                    if (inst.rowsum_lower[i] <= rowsum[i] + diff &&
                                        rowsum[i] + diff <= inst.rowsum_upper[i]) {

                                        // create copy of s
                                        BinaryMatrix s2(*s);

                                        // flip selected entries
                                        for (int l = 0; l < k; l++) {
                                            int j = tmp3[l];
                                            s2.flip(i, j);
                                        }

                                        process(&s2, p);

                                    } else {
                                        loop += p;
                                    }
                                } while(cg_sub.next());
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
                        const BinaryMatrix *s,
                        const std::function<void(const State *, const Rational &)> &process,
                        const int *rowsum,
                        const int *colsum
                ) const {

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    Rational loop(0);

                    /* simulate vertical multi-fliop */

                    // randomly select a column index
                    for (int j = 0; j < ncol; j++) {

                        // select the positions i of a) the entries sij=1 that can be switched to zero
                        // and b) the entries sij=0 that can be switched to one
                        int a = 0;
                        int b = 0;

                        // for each row position
                        for (int i = 0; i < nrow; i++) {
                            bool sij = s->get(i, j);
                            if (sij && rowsum[i] > inst.rowsum_lower[i]) {
                                tmp1[a + b] = i;
                                a++;
                            } else if (!sij && rowsum[i] < inst.rowsum_upper[i]) {
                                tmp1[a + b] = i;
                                b++;
                            }
                        }

                        // determine maximal number of positions to flip
                        int range = inst.colsum_upper[j] - inst.colsum_lower[j];
                        int t = std::min(range, a + b);

                        // if there is nothing to choose from
                        if (t == 0) {
                            loop += q_vertical_rat * p_multiflip_rat * Rational(1, ncol);
                            continue;
                        }

                        // simulate randomly selecting t out of (a+b) elements
                        CombinationGenerator<int> cg(tmp1, tmp2, a + b, t);
                        do {

                            // simulate selecting random subsets from tmp2

                            // each of the 2^t subsets has the same probability p of being selected
                            const Rational p =
                                    q_vertical_rat * p_multiflip_rat * Rational(1, ncol) /
                                    (binom(a + b, t) * pow(2, t));

                            // simulate selecting subset of size k
                            for (int k = 0; k <= t; k++) {

                                // simulate selecting k objects out of t
                                CombinationGenerator<int> cg_sub(tmp2, tmp3, t, k);
                                do {

                                    // determine by how much colsum[j] would be changed by a flip of the k selected entries
                                    int diff = 0;
                                    for (int l = 0; l < k; l++) {
                                        int i = tmp3[l];
                                        diff += s->get(i, j) ? -1 : 1;
                                    }

                                    // if entries can be flipped without violating the lower and upper bounds
                                    if (inst.colsum_lower[j] <= colsum[j] + diff &&
                                        colsum[j] + diff <= inst.colsum_upper[j]) {

                                        // create copy of s
                                        BinaryMatrix s2(*s);

                                        // flip selected entries
                                        for (int l = 0; l < k; l++) {
                                            int i = tmp3[l];
                                            s2.flip(i, j);
                                        }

                                        process(&s2, p);

                                    } else {
                                        loop += p;
                                    }
                                } while(cg_sub.next());
                            }

                        } while (cg.next());
                    }

                    // process loop
                    process(s, loop);
                }

                void init() {
                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();
                    tmp1 = new int[std::max(nrow, ncol)];
                    tmp2 = new int[std::max(nrow, ncol)];
                    tmp3 = new int[std::max(nrow, ncol)];

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
                explicit InformedChain(const Instance &inst) :
                        MarkovChain(inst),
                        q_horizontal_rat(Rational(inst.getNumCols(), inst.getNumRows() + inst.getNumCols())),
                        q_vertical_rat(Rational(inst.getNumRows(), inst.getNumRows() + inst.getNumCols())),
                        q_horizontal(inst.getNumCols() / (double) (inst.getNumRows() + inst.getNumCols())),
                        q_vertical(inst.getNumRows() / (double) (inst.getNumRows() + inst.getNumCols())) {
                    init();
                }

                /**
                 * Create a Markov chain. Use a certain binary matrix as initial state.
                 * @param inst Lower and upper bounds on the row and column sums.
                 * @param bin Binary matrix.
                 */
                InformedChain(const Instance &inst, const BinaryMatrix &bin)
                        : MarkovChain(inst, bin),
                        q_horizontal_rat(Rational(inst.getNumCols(), inst.getNumRows() + inst.getNumCols())),
                        q_vertical_rat(Rational(inst.getNumRows(), inst.getNumRows() + inst.getNumCols())),
                        q_horizontal(inst.getNumCols() / (double) (inst.getNumRows() + inst.getNumCols())),
                        q_vertical(inst.getNumRows() / (double) (inst.getNumRows() + inst.getNumCols())) {
                    init();
                }


                /**
                 * Create a Markov chain.
                 * @param inst String-encoded instance of lower and upper bounds on the row and column sums.
                 * Instances for this chain have the form "l1-u1,l2-u2,l3-u3;l4-u4,l5-u5", where li is the ith lower
                 * bound and ui is the ith upper bound. For convenience, if li = ui, the string 'li-ui' can be replaced
                 * by 'li'. The semicolon separates the row sums from the column sums.
                 */
                explicit InformedChain(const std::string &inst)
                        : InformedChain(Instance(inst)) {

                }


                /**
                 * Create a Markov chain.
                 * @param rowsum_lower Sequence of lower bounds for row sums.
                 * @param rowsum_upper Sequence of upper bounds for row sums.
                 * @param colsum_lower Sequence of lower bounds for column sums.
                 * @param colsum_upper Sequence of upper bounds for column sums.
                 */
                InformedChain(
                        const std::vector<int> &rowsum_lower,
                        const std::vector<int> &rowsum_upper,
                        const std::vector<int> &colsum_lower,
                        const std::vector<int> &colsum_upper
                ) : InformedChain(Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper)) {

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
                InformedChain(
                        const int *rowsum_lower,
                        const int *rowsum_upper,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        const int nrow,
                        const int ncol
                ) : InformedChain(
                        Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper, nrow, ncol)) {

                }

                virtual ~InformedChain() {
                    delete[] tmp1;
                    delete[] tmp2;
                    delete[] tmp3;
                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s
                 * @param process
                 */
                void adjacentStates(
                        const State *x,
                        const std::function<void(const State *, const marathon::Rational &)> &process
                ) const override {

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
                            if(s->get(i,j)) {
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

                    /* clean up */
                    delete[] rowsum;
                    delete[] colsum;
                }

                /**
                 * Randomize the current state of the Markov chain.
                 */
                virtual void step() override {

                    BinaryMatrix *s = (BinaryMatrix *) currentState;

                    // select two random numbers from [0,1)
                    double p = rg.nextDouble();
                    double q = rg.nextDouble();

                    // which operation will be applied?
                    if (p < p_trade) {                                   // apply trade

                        if (q < q_horizontal) {
                            applyTradeHorizontal(s);
                        } else if (q < q_horizontal + q_vertical) {
                            applyTradeVertical(s);
                        } else {
                            throw std::runtime_error("Error while applying trade.");
                        }

                    } else if (p < p_trade + p_multishift) {      // apply multi-shift

                        if (q < q_horizontal) {
                            applyMultiShiftHorizontal(s);
                        } else if (q < q_horizontal + q_vertical) {
                            applyMultiShiftVertical(s);
                        } else {
                            throw std::runtime_error("Error while applying multi-shift.");
                        }

                    } else if (p < p_trade + p_multishift + p_multiflip) {                                                       // apply multi-flip

                        if (q < q_horizontal) {
                            applyMultiFlipHorizontal(s);
                        } else if(q < q_horizontal + q_vertical) {
                            applyMultiFlipVertical(s);
                        } else {
                            throw std::runtime_error("Error while applying multi-flip.");
                        }
                    }
                    else {
                        throw std::runtime_error("Error while applying operation.");
                    }
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual InformedChain *copy() const override {
                    auto mc = new InformedChain(inst);
                    mc->setCurrentState(this->getCurrentState());
                    return mc;
                }

            };
        }
    }
}

#endif /* _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_INFORMED_CHAIN_H */
