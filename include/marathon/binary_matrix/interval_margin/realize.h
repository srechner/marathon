/*
 * Created on: May 4, 2017
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


#ifndef MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_REALIZATION_H
#define MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_REALIZATION_H

#include "marathon/binary_matrix/binary_matrix.h"
#include "marathon/binary_matrix/interval_margin/instance.h"
#include "marathon/binary_matrix/fixed_margin/realize.h"

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * Test the realisability of a bipartite degree sequence that is described
             * by the lower and upper bounds.
             * @param inst Four-Tuple of integer vectors.
             * @return True, if a graphical sequence exists that does not violate the upper
             * and lower bounds, or False, if no such sequence exists.
             */
            inline
            bool isRealizable(const Instance &seq) {

                const int nrow = seq._rowsum_lower.size();
                const int ncol = seq._colsum_lower.size();

                // are all row sums in the right range?
                for (int i = 0; i < nrow; i++) {
                    if (seq._rowsum_upper[i] < seq._rowsum_lower[i] || seq._rowsum_lower[i] > ncol) {
                        return false;
                    }
                }

                // are all column sums in the right range?
                for (int j = 0; j < ncol; j++) {
                    if (seq._colsum_upper[j] < seq._colsum_lower[j] || seq._colsum_lower[j] > nrow) {
                        return false;
                    }
                }

                // create sorted sequence of lower row sum
                int *rowsum_lower_sorted = new int[nrow];
                for (int i = 0; i < nrow; i++)
                    rowsum_lower_sorted[i] = std::max(seq._rowsum_lower[i], 0);
                std::sort(rowsum_lower_sorted, rowsum_lower_sorted + nrow, [](const int a, const int b) {
                    return a > b;
                });

                // create sorted sequence of lower column sum
                int *colsum_lower_sorted = new int[ncol];
                for (int j = 0; j < ncol; j++)
                    colsum_lower_sorted[j] = std::max(seq._colsum_lower[j], 0);
                std::sort(colsum_lower_sorted, colsum_lower_sorted + ncol, [](const int a, const int b) {
                    return a > b;
                });

                // create conjugates sequence of upper rowsum
                int *rowsum_upper_conjugated = new int[ncol];
                conjugate(rowsum_upper_conjugated, &seq._rowsum_upper[0], ncol, nrow);

                // create conjugates sequence of upper colsum
                int *colsum_upper_conjugated = new int[nrow];
                conjugate(colsum_upper_conjugated, &seq._colsum_upper[0], nrow, ncol);

                // the sequence tuple is realizable when
                // a) the conjugated upper row sums dominate the lower column sums
                bool realizable = true;
                int S1 = 0;
                int S2 = 0;
                for (int k = 0; k < ncol; k++) {
                    S1 += rowsum_upper_conjugated[k];
                    S2 += colsum_lower_sorted[k];
                    if (S1 < S2) {
                        realizable = false;
                        break;
                    }
                }

                // b) the conjugated upper column sums dominate the lower row sums
                if (realizable) {
                    S1 = 0;
                    S2 = 0;
                    for (int k = 0; k < nrow; k++) {
                        S1 += colsum_upper_conjugated[k];
                        S2 += rowsum_lower_sorted[k];
                        if (S1 < S2) {
                            realizable = false;
                            break;
                        }
                    }
                }

                delete[] rowsum_lower_sorted;
                delete[] colsum_lower_sorted;
                delete[] rowsum_upper_conjugated;
                delete[] colsum_upper_conjugated;

                return realizable;
            }


            /**
             * Test the realisability of a bipartite degree sequence that is described
             * by the lower and upper bounds.
             * @param rowsum_lower Lower bound for row sums.
             * @param rowsum_upper Upper bound for row sums.
             * @param colsum_lower Lower bound for column sums.
             * @param colsum_upper Upper bound for column sums.
             * @param nrow Number of rows.
             * @param ncol Number of columns.
             * @return True, if a graphical sequence exists that does not violate the upper
             * and lower bounds, or False, if no such sequence exists.
             */
            inline
            bool isRealizable(
                    const int *rowsum_lower,
                    const int *rowsum_upper,
                    const int *colsum_lower,
                    const int *colsum_upper,
                    const int nrow,
                    const int ncol
            ) {
                return isRealizable(Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper, nrow, ncol));
            }

            namespace detail {

                struct A {
                    size_t index;
                    int value;
                };

                /**
                 * Rise the lower column sums to its maximum, such that my_colsum_lower[i] <= my_colsum_upper[i]
                 * for each i and my_rowsum_upper majorize my_colsum_lower.
                 * @param my_rowsum_lower
                 * @param my_rowsum_upper
                 * @param my_colsum_lower
                 * @param my_colsum_upper
                 * @param nrow
                 * @param ncol
                 * @param col_order
                 */
                int rise(
                        int *my_rowsum_lower,
                        int *my_rowsum_upper,
                        int *my_colsum_lower,
                        int *my_colsum_upper,
                        const int nrow,
                        const int ncol,
                        A *col_order
                ) {

                    // create conjugates sequence of lower and upper rowsum
                    int *rowsum_upper_conjugated = new int[ncol];
                    memset(rowsum_upper_conjugated, 0, ncol * sizeof(int));
                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < my_rowsum_upper[i]; j++)
                            rowsum_upper_conjugated[j]++;
                    }

                    // create arrays of accumulated sums
                    int *sum_rowsum_upper_conj = new int[ncol];
                    int *sum_colsum_lower = new int[ncol];
                    sum_rowsum_upper_conj[0] = rowsum_upper_conjugated[0];
                    sum_colsum_lower[0] = my_colsum_lower[0];
                    for (int j = 1; j < ncol; j++) {
                        sum_rowsum_upper_conj[j] = sum_rowsum_upper_conj[j - 1] + rowsum_upper_conjugated[j];
                        sum_colsum_lower[j] = sum_colsum_lower[j - 1] + my_colsum_lower[j];
                    }

                    // choose k as the largest position such that sum_rowsum_upper_conj[k] == sum_colsum_lower[k]
                    int k = -1;
                    for (int j = ncol - 1; j >= 0; j--) {
                        if (sum_rowsum_upper_conj[j] == sum_colsum_lower[j]) {
                            k = j;
                            break;
                        }
                    }

                    for (int j = k + 1; j < ncol && sum_rowsum_upper_conj[ncol - 1] > sum_colsum_lower[ncol - 1]; j++) {

                        // can my_colsum_lower[j] be increased?
                        if (my_colsum_lower[j] < my_colsum_upper[j]) {

                            // j is the smallest index > k such that the my_colsum_lower[j] < my_colsum_upper[j]

                            // determine maximal amount of which my_colsum_lower[j] can be increased
                            int min_slack = sum_rowsum_upper_conj[j] - sum_colsum_lower[j];
                            for (int l = j + 1; l < ncol; l++)
                                min_slack = std::min(min_slack, sum_rowsum_upper_conj[l] - sum_colsum_lower[l]);

                            // iteratively increase my_colsum_lower[j] and shift it to the left
                            int l = j;
                            while (min_slack > 0 && my_colsum_lower[l] < my_colsum_upper[l]) {

                                // find maximal position i such that my_colsum_lower[i] == my_colsum_lower[l]
                                int i = l;
                                while (i > 0 && my_colsum_lower[i - 1] == my_colsum_lower[l]) {
                                    min_slack = std::min(min_slack,
                                                         sum_rowsum_upper_conj[i - 1] - sum_colsum_lower[i - 1]);
                                    i--;
                                }

                                // determine the largest amount of which my_colsum_lower[l] can be increased
                                int diff = std::min(min_slack, my_colsum_upper[l] - my_colsum_lower[l]);

                                // the new value of my_colsum_lower[i] must not exceed my_colsum_lower[i-1]
                                if (i > 0)
                                    diff = std::min(diff, my_colsum_lower[i - 1] - my_colsum_lower[i]);

                                // increase my_colsum_lower[i] by diff
                                my_colsum_lower[i] += diff;

                                // swap the columns i and l
                                std::swap(my_colsum_upper[i], my_colsum_upper[l]);
                                std::swap(col_order[i], col_order[l]);

                                // update value of slack
                                min_slack -= diff;

                                l = i;
                            }

                            // update accumulated column sums todo: improve!
                            sum_colsum_lower[0] = my_colsum_lower[0];
                            for (int l = 1; l < ncol; l++)
                                sum_colsum_lower[l] = sum_colsum_lower[l - 1] + my_colsum_lower[l];

                            // advance k to the largest position such that sum_rowsum_upper_conj[k] == sum_colsum_lower[k] todo: improve!
                            for (int l = 0; l < ncol; l++) {
                                if (sum_rowsum_upper_conj[l] == sum_colsum_lower[l]) {
                                    k = l;
                                    j = k;
                                }
                            }
                        }
                    }

                    int left = sum_rowsum_upper_conj[ncol - 1] - sum_colsum_lower[ncol - 1];

                    delete[] rowsum_upper_conjugated;
                    delete[] sum_colsum_lower;
                    delete[] sum_rowsum_upper_conj;

                    return left;
                }


                /**
                 * Transform a realizable four-tuple into a bi-graphical pair of
                 * integer vectors.
                 * @param Lower and upper bounds on row and column sums.
                 * @return Bi-graphical pair of integer vectors.
                 */
                fixed_margin::Instance
                construct_bigraphical(const interval_margin::Instance &margin) {

                    /**************************************************************************
                     * Algorithm presented in
                     *
                     *    S. Rechner. An Optimal Realization Algorithm for Bipartite Graphs
                     *    with Degrees in Prescribed Intervals. arXiv preprint (2017).
                     *    arXiv:1708.05520v1 [cs.DS].
                     *************************************************************************/

                    /**************************************************************************
                     * Step Zero: Preparations
                     *************************************************************************/

                    const int nrow = (int) margin._rowsum_lower.size();
                    const int ncol = (int) margin._colsum_lower.size();

                    struct A {
                        int index;
                        int lower;
                        int upper;
                    };


                    A *rows = new A[nrow];
                    A *cols = new A[ncol];
                    for (int i = 0; i < nrow; i++) {
                        rows[i].index = i;
                        rows[i].lower = std::max(margin._rowsum_lower[i], 0);
                        rows[i].upper = std::min(margin._rowsum_upper[i], ncol);
                    }
                    for (int j = 0; j < ncol; j++) {
                        cols[j].index = j;
                        cols[j].lower = std::max(margin._colsum_lower[j], 0);
                        cols[j].upper = std::min(margin._colsum_upper[j], nrow);
                    }

                    // sort rows by lower bound
                    std::sort(rows, rows + nrow,
                              [](const A &a, const A &b) {
                                  return a.lower == b.lower ? a.upper > b.upper : a.lower > b.lower;
                              });

                    // sort columns by lower bound
                    std::sort(cols, cols + ncol,
                              [](const A &a, const A &b) {
                                  return a.lower == b.lower ? a.upper > b.upper : a.lower > b.lower;
                              });

                    // create sequences of row sums and column sums that ultimately will contain the bi-graphical sequence
                    fixed_margin::Instance res(nrow, ncol);

                    // create conjugate sequence to lower column sums
                    int *colsum_lower_conj = new int[nrow];
                    marathon::binary_matrix::conjugate(colsum_lower_conj, &margin._colsum_lower[0], nrow, ncol);

                    // first and last occurrence of each value in cols.lower
                    int *colsum_lower_first = new int[nrow + 1];
                    for (int i = 0; i <= nrow; i++)
                        colsum_lower_first[i] = INT_MAX;
                    for (int j = ncol - 1; j >= 0; j--)
                        colsum_lower_first[cols[j].lower] = j;

                    /***************************************************************************
                     * Step One: Rise lower column sums.
                     **************************************************************************/

                    // determine amount by which the lower column sums must be increased
                    int rowsum_lower_agg = std::accumulate(rows, rows + nrow, 0,
                                                           [](const int res, const A &a) { return res + a.lower; });
                    int colsum_lower_conj_agg = std::accumulate(colsum_lower_conj, colsum_lower_conj + nrow, 0);
                    int delta1 = 0;
                    for (int i = nrow - 1; i >= 0; i--) {
                        const int slack = rowsum_lower_agg - colsum_lower_conj_agg;
                        delta1 = std::max(delta1, slack);
                        rowsum_lower_agg -= rows[i].lower;
                        colsum_lower_conj_agg -= colsum_lower_conj[i];
                    }

                    int x = delta1;

                    // k will be the right-most position such that cols[k].upper can be increased further
                    int k = ncol - 1;
                    while (delta1 > 0) {

                        // proceed to right-most position such that cols[k].upper can be increased further
                        while (cols[k].lower == cols[k].upper)
                            k--;

                        const int x = cols[k].lower;

                        // find smallest position f such that rows[f].lower = rows[k].lower
                        const int f = colsum_lower_first[x];
                        assert(f != INT_MAX);

                        // switch columns
                        std::swap(cols[f], cols[k]);

                        // increment cols[f].lower
                        cols[f].lower++;

                        // update first array
                        if (f < ncol && cols[f + 1].lower == x)
                            colsum_lower_first[x]++;
                        else
                            colsum_lower_first[x] = INT_MAX;
                        if (colsum_lower_first[x + 1] == INT_MAX)
                            colsum_lower_first[x + 1] = f;

                        delta1--;
                    }

                    // copy final column sums
                    for (int j = 0; j < ncol; j++)
                        res._colsum[j] = cols[j].lower;

                    // create conjugate sequence of column sums
                    std::vector<int> colsum_conj = conjugate(res._colsum, nrow);

                    // create conjugate sequence to upper row sums
                    std::vector<int> rowsum_upper_conj = conjugate(margin._rowsum_upper, ncol);

                    // at this point there is a binary matrix whose column sums are described by colsum
                    assert(marathon::binary_matrix::isDominating(colsum_conj, margin._rowsum_lower));
                    assert(marathon::binary_matrix::isDominating(rowsum_upper_conj, res._colsum));
                    assert(marathon::binary_matrix::interval_margin::isRealizable(
                            &margin._rowsum_lower[0], &margin._rowsum_upper[0], &res._colsum[0], &res._colsum[0], nrow,
                            ncol));

                    /**************************************************************************
                     * Step Two: Rise the lower row sums
                     *************************************************************************/

                    // determine amount by which rows[].lower must be increased
                    int delta2 = std::accumulate(res._colsum.begin(), res._colsum.end(), 0) -
                                 std::accumulate(rows, rows + nrow, 0,
                                                 [](const int res, const A &a) { return res + a.lower; });
                    assert(delta2 >= 0);

                    int y = delta2;

                    // first occurrence of each rows[
                    int *rowsum_lower_first = new int[ncol + 1];
                    for (int i = 0; i <= ncol; i++)
                        rowsum_lower_first[i] = INT_MAX;
                    for (int i = nrow - 1; i >= 0; i--)
                        rowsum_lower_first[rows[i].lower] = i;

                    // right-most position such that rows[k].lower can be increased
                    k = nrow - 1;
                    while (delta2 > 0) {

                        // proceed to right-most position such that rows[k].lower can be increased
                        while (rows[k].lower == rows[k].upper)
                            k--;

                        const int x = rows[k].lower;

                        // determine right-most position f such that rows[f].lower == rows[k].lower
                        const int f = rowsum_lower_first[x];
                        assert(f != INT_MAX);

                        // switch rows
                        std::swap(rows[f], rows[k]);

                        // increment rows[f].lower
                        rows[f].lower++;

                        // update last array
                        if (f < nrow && rows[f + 1].lower == x)
                            rowsum_lower_first[x]++;
                        else
                            rowsum_lower_first[x] = INT_MAX;
                        if (rowsum_lower_first[x + 1] == INT_MAX)
                            rowsum_lower_first[x + 1] = f;

                        delta2--;
                    }

                    // restore original row and column order
                    for (int i = 0; i < nrow; i++) {
                        res._rowsum[rows[i].index] = rows[i].lower;
                        res._rowindex[i] = i;
                    }
                    for (int j = 0; j < ncol; j++) {
                        res._colsum[cols[j].index] = cols[j].lower;
                        res._colindex[j] = j;
                    }

                    // at this point, rowsum and colsum are realizable
                    assert(marathon::binary_matrix::fixed_margin::isRealizable(res));

                    // clean up
                    delete[] rowsum_lower_first;
                    delete[] colsum_lower_conj;
                    delete[] colsum_lower_first;
                    delete[] rows;
                    delete[] cols;

                    return res;
                }
            }

            /**
             * Construct a binary matrix of size nrow times ncol
             * whose row and column sums lie in the intervals
             * prescribed the lower and upper sums.
             * @param margin Lower and upper bounds on row and Column sums.
             */
            marathon::binary_matrix::BinaryMatrix
            realize_fast(const Instance &margin) {

                // check realizability
                if (!isRealizable(margin))
                    throw std::runtime_error("Error! Four-tuple of integer vectors is not realisable!");

                // transform four-tuple into a bi-graphical pair of integer vectors
                auto bigraphical = detail::construct_bigraphical(margin);

                // determine a realization
                return marathon::binary_matrix::fixed_margin::realize(bigraphical);
            }

            /**
             * Construct a binary matrix of size nrow times ncol
             * whose row and column sums lie in the intervals
             * prescribed the lower and upper sums.
             * @param rowsum_lower Lower bounds on each row sum.
             * @param rowsum_upper Upper bounds on each row sum.
             * @param colsum_lower Lower bounds on each column sum.
             * @param colsum_upper Upper bounds on each column sum.
             */
            marathon::binary_matrix::BinaryMatrix
            realize_fast(
                    const int *rowsum_lower,
                    const int *rowsum_upper,
                    const int *colsum_lower,
                    const int *colsum_upper,
                    const int nrow,
                    const int ncol
            ) {
                return realize_fast(Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper, nrow, ncol));
            }


            /**
             * Construct a random binary matrix whose margins lie in prescribed intervals.
             * @param margin Lower and upper bounds on row and column sums.
             */
            marathon::binary_matrix::BinaryMatrix
            realize_fast_random(const Instance &margin) {

                // check realizability
                if (!isRealizable(margin))
                    throw std::runtime_error("Error! Four-tuple of integer vectors is not realisable!");

                // transform four-tuple into a bi-graphical pair of integer vectors
                auto bigraphical = detail::construct_bigraphical(margin);

                // determine a random realization
                auto bin = marathon::binary_matrix::fixed_margin::realize_random(bigraphical);

                return bin;

            }

            /**
             * Construct a random binary matrix of size nrow times ncol
             * whose row and column sums lie in the intervals
             * prescribed the lower and upper sums.
             * @param rowsum_lower Lower bounds on each row sum.
             * @param rowsum_upper Upper bounds on each row sum.
             * @param colsum_lower Lower bounds on each column sum.
             * @param colsum_upper Upper bounds on each column sum.
             */
            marathon::binary_matrix::BinaryMatrix
            realize_fast_random(
                    const int *rowsum_lower,
                    const int *rowsum_upper,
                    const int *colsum_lower,
                    const int *colsum_upper,
                    const int nrow,
                    const int ncol
            ) {
                return realize_fast_random(
                        Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper, nrow, ncol));
            }


            /**
             * Construct a binary matrix of size nrow times ncol
             * whose row and column sums lie in the intervals
             * prescribed the lower and upper sums.
             * @param rowsum_lower Lower bounds on each row sum.
             * @param rowsum_upper Upper bounds on each row sum.
             * @param colsum_lower Lower bounds on each column sum.
             * @param colsum_upper Upper bounds on each column sum.
             * @param nrow The number of rows.
             * @param ncol The number of columns.
             */
            inline
            marathon::binary_matrix::BinaryMatrix
            realize_slow(const Instance &seq) {

                const int nrow = seq._rowsum_upper.size();
                const int ncol = seq._colsum_lower.size();

                if (!isRealizable(seq))
                    throw std::runtime_error("Error! Four-tuple of integer vectors is not realizable!");

                /**************************************************************************
                 * Transform the lower column sum and the upper row sum to a realizable
                 * sequence. The algorithm has been derived from the proof of theorem 1.6
                 * of Manfred Schocker: 'On Graphs with Degrees in Prescribed Intervals'
                 *************************************************************************/

                int *my_rowsum_lower = new int[nrow];
                int *my_rowsum_upper = new int[nrow];
                int *my_colsum_lower = new int[ncol];
                int *my_colsum_upper = new int[ncol];

                // create non-increasing order of row sums
                detail::A *row_order = new detail::A[nrow];
                for (int i = 0; i < nrow; i++)
                    row_order[i] = {seq._rowindex[i], seq._rowsum_lower[i]};
                std::sort(row_order, row_order + nrow, [](const detail::A &a, const detail::A &b) {
                    return a.value > b.value;
                });
                for (int i = 0; i < nrow; i++) {
                    my_rowsum_lower[i] = seq._rowsum_lower[row_order[i].index];
                    my_rowsum_upper[i] = seq._rowsum_upper[row_order[i].index];
                }

                // create non-increasing order of column sums
                detail::A *col_order = new detail::A[ncol];
                for (int j = 0; j < ncol; j++)
                    col_order[j] = {seq._colindex[j], seq._colsum_lower[j]};
                std::sort(col_order, col_order + ncol, [](const detail::A &a, const detail::A &b) {
                    return a.value > b.value;
                });
                for (int j = 0; j < ncol; j++) {
                    my_colsum_lower[j] = seq._colsum_lower[col_order[j].index];
                    my_colsum_upper[j] = seq._colsum_upper[col_order[j].index];
                }

                // trim row sums to fit bounds
                for (int i = 0; i < nrow; i++) {
                    my_rowsum_lower[i] = std::max(my_rowsum_lower[i], 0);
                    my_rowsum_upper[i] = std::min(my_rowsum_upper[i], ncol);
                }

                // trim column sums to fit bounds
                for (int j = 0; j < ncol; j++) {
                    my_colsum_lower[j] = std::max(my_colsum_lower[j], 0);
                    my_colsum_upper[j] = std::min(my_colsum_upper[j], nrow);
                }

                /**************************************************************************
                 * Phase One: Rise the lower column sums.
                 *************************************************************************/
                detail::rise(my_rowsum_lower, my_rowsum_upper, my_colsum_lower, my_colsum_upper, nrow, ncol,
                             col_order);

                /**************************************************************************
                 * Phase Two: Rise the lower row sums.
                 *************************************************************************/
                detail::rise(my_colsum_lower, my_colsum_upper, my_rowsum_lower, my_rowsum_upper, ncol, nrow,
                             row_order);

                /**************************************************************
                 * Realize mofified sequence via H-H-algorithm.
                 *************************************************************/
                marathon::binary_matrix::BinaryMatrix bin_tmp =
                        marathon::binary_matrix::fixed_margin::realize(
                                my_rowsum_lower,
                                my_colsum_lower,
                                nrow,
                                ncol
                        );

                /**************************************************************
                 * Re-shuffle rows and columns to original order.
                 *************************************************************/
                marathon::binary_matrix::BinaryMatrix bin(nrow, ncol);
                for (int i = 0; i < nrow; i++) {
                    for (int j = 0; j < ncol; j++) {
                        int a = row_order[i].index;
                        int b = col_order[j].index;
                        bin.set(a, b, bin_tmp.get(i, j));
                    }
                }

                delete[] my_rowsum_lower;
                delete[] my_rowsum_upper;
                delete[] my_colsum_lower;
                delete[] my_colsum_upper;

                return bin;

            }

            /**
             * Construct a binary matrix of size nrow times ncol
             * whose row and column sums lie in the intervals
             * prescribed the lower and upper sums.
             * @param rowsum_lower Lower bounds on each row sum.
             * @param rowsum_upper Upper bounds on each row sum.
             * @param colsum_lower Lower bounds on each column sum.
             * @param colsum_upper Upper bounds on each column sum.
             * @param nrow The number of rows.
             * @param ncol The number of columns.
             */
            inline
            marathon::binary_matrix::BinaryMatrix
            realize_slow(
                    const int *rowsum_lower,
                    const int *rowsum_upper,
                    const int *colsum_lower,
                    const int *colsum_upper,
                    const int nrow,
                    const int ncol
            ) {
                return realize_slow(Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper, nrow, ncol));
            }


            /**
             * Construct a binary matrix of size nrow times ncol
             * whose row and column sums lie in the intervals
             * prescribed the lower and upper sums.
             * @param rowsum_lower Lower bounds on each row sum.
             * @param rowsum_upper Upper bounds on each row sum.
             * @param colsum_lower Lower bounds on each column sum.
             * @param colsum_upper Upper bounds on each column sum.
             */
            inline
            marathon::binary_matrix::BinaryMatrix
            realize_slow(
                    const std::vector<int> &rowsum_lower,
                    const std::vector<int> &rowsum_upper,
                    const std::vector<int> &colsum_lower,
                    const std::vector<int> &colsum_upper
            ) {

                return realize_slow(&rowsum_lower[0], &rowsum_upper[0],
                                    &colsum_lower[0], &colsum_upper[0],
                                    rowsum_lower.size(), colsum_lower.size());
            }
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_REALIZATION_H
