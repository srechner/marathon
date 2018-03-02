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


#ifndef MARATHON_FIXED_MARGIN_REALIZATION_H
#define MARATHON_FIXED_MARGIN_REALIZATION_H

#include "marathon/basic_random.h"
#include "marathon/binary_matrix/binary_matrix.h"
#include "marathon/binary_matrix/common.h"
#include "instance.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Return true, if a binary matrix exists that possesses the given row and column sum.
             * @param rowsum Sequence of row sums.
             * @param colsum Sequence of column sums.
             * @param nrow Number of rows.
             * @param ncol Number of columns.
             * @return True, if (rowsum, colsum) is a realizable sequence pair. False, otherwise.
             */
            inline
            bool isRealizable(const Instance& m) {

                // the total sum of each sequence
                int sum_rowsum = 0;
                int sum_colsum = 0;

                // are all row sums in the right range?
                for (int i = 0; i < m.getNumRows(); i++) {
                    sum_rowsum += m.rowsum[i];
                    if (m.rowsum[i] < 0 || m.rowsum[i] > m.getNumCols()) {
                        return false;
                    }
                }

                // are all column sums in the right range?
                for (int j = 0; j < m.getNumCols(); j++) {
                    sum_colsum += m.colsum[j];
                    if (m.colsum[j] < 0 || m.colsum[j] > m.getNumRows()) {
                        return false;
                    }
                }

                // do they have the same total?
                if (sum_rowsum != sum_colsum)
                    return false;

                // compute conjugated row sums
                int *rowsum_conjugated = new int[m.getNumCols()];
                conjugate(rowsum_conjugated, &m.rowsum[0], m.getNumCols(), m.getNumRows());

                // sort column sums descendingly
                int *colsum_sorted = new int[m.getNumCols()];
                memcpy(colsum_sorted, &m.colsum[0], m.getNumCols() * sizeof(int));
                std::sort(colsum_sorted, colsum_sorted + m.getNumCols(), [](const int a, const int b) { return a > b; });

                // the degree sequence is realizable when the conjugated rowsum dominate the colsum
                bool realizable = true;
                int S1 = 0;
                int S2 = 0;
                for (int k = 0; k < m.getNumCols(); k++) {
                    S1 += rowsum_conjugated[k];
                    S2 += colsum_sorted[k];
                    if (S1 < S2) {
                        realizable = false;
                        break;
                    }
                }

                delete[] rowsum_conjugated;
                delete[] colsum_sorted;

                return realizable;
            }

            /**
             * Return true, if a binary matrix exists that possesses the given row and column sum.
             * @param rowsum Sequence of row sums.
             * @param colsum Sequence of column sums.
             * @param nrow Number of rows.
             * @param ncol Number of columns.
             * @return True, if (rowsum, colsum) is a realizable sequence pair. False, otherwise.
             */
            inline
            bool isRealizable(
                    const int *rowsum,
                    const int *colsum,
                    const int nrow,
                    const int ncol
            ) {
                return isRealizable(Instance(rowsum, colsum, nrow, ncol));
            }

            /**
             * Return true, if a binary matrix exists that possesses the given row and column sum.
             * @param rowsum Sequence of row sums.
             * @param colsum Sequence of column sums.
             * @return True, if (rowsum, colsum) is a realizable sequence pair. False, otherwise.
             */
            inline
            bool isRealizable(
                    const std::vector<int> &rowsum,
                    const std::vector<int> &colsum
            ) {
                return isRealizable(&rowsum[0], &colsum[0], rowsum.size(), colsum.size());
            }

            /**
             * Construct a binary matrix of size nrow times ncol
             * whose row and column sums match the given integer
             * sequence.
             * @param seq Sequence of row and column sums.
             */
            inline
            BinaryMatrix *
            realize(const Instance &seq) {

                const int nrow = (int) seq.getNumRows();
                const int ncol = (int) seq.getNumCols();

                /**************************************************************
                 * Running Time: O(nrow + ncol*log(ncol) + total),
                 * where total = sum(rowsum) = sum(colsum).
                 *************************************************************/

                /**************************************************************
                 * Check realizability.
                 *************************************************************/

                struct A {
                    int index;
                    int degree;
                };

                // sorted vector of index-degree-pairs
                A *column = new A[ncol];
                for (int j = 0; j < ncol; j++) {
                    column[j] = {seq.colindex[j], seq.colsum[j]};
                }

                // sort column pairs descendingly by degree
                std::sort(column, column + ncol, [](const A &a, const A &b) {
                    return a.degree > b.degree;
                });

                // create conjugate sequence of row sums
                int *rowsum_conjugated = new int[ncol];
                conjugate(rowsum_conjugated, &seq.rowsum[0], ncol, nrow);

                int rowsum_conjugated_sum = 0;
                int colsum_sum = 0;
                for (int j = 0; j < ncol; j++) {
                    rowsum_conjugated_sum += rowsum_conjugated[j];
                    colsum_sum += column[j].degree;
                    if (rowsum_conjugated_sum < colsum_sum) {
                        // sequence is not realizable
                        delete[] column;
                        delete[] rowsum_conjugated;
                        return nullptr;
                    }
                }

                /******************************************************************
                 * Apply Ryser's algorithm.
                 *****************************************************************/

                BinaryMatrix *bip = new BinaryMatrix(nrow, ncol);

                // last[k] = max { j : column[j].degree == k }
                int *last = new int[nrow + 1];
                for (int k = 0; k < nrow + 1; k++)
                    last[k] = -1;
                for (int j = 0; j < ncol; j++) {
                    last[column[j].degree] = j;
                }

                // for each row
                for (int x = 0; x < nrow; x++) {

                    // index of current row
                    const int i = seq.rowindex[x];
                    const int di = seq.rowsum[x];

                    // distribute di ones in row i
                    for (int k = di - 1; k >= 0; k--) {

                        // index j is node with kth-largest degree
                        int j = column[k].index;
                        int dj = column[k].degree;
                        assert(dj != 0);

                        bip->set(i, j, 1);
                        column[k].degree--;

                        // find largest index l such that column[l].degree == d
                        int l = last[dj];

                        // repair sorting of column vector by swapping elements k and l
                        std::swap(column[l], column[k]);

                        // update last array
                        if (l > 0 && column[l - 1].degree == dj)
                            last[dj] = l - 1;
                        else
                            last[dj] = -1;
                        if (last[dj - 1] == -1)
                            last[dj - 1] = l;
                    }
                }

                delete[] column;
                delete[] rowsum_conjugated;
                delete[] last;

                return bip;
            }


            /**
             * Construct a binary matrix of size nrow times ncol
             * whose row and column sums match the given integer
             * sequence.
             * @param rowsum Sequence of row sums.
             * @param colsum Sequence of column sums.
             * @param nrow The number of rows.
             * @param ncol The number of columns.
             */
            inline
            BinaryMatrix *
            realize(
                    const int *rowsum,
                    const int *colsum,
                    const int nrow,
                    const int ncol
            ) {
                return realize(Instance(rowsum, colsum, nrow, ncol));
            }


            /**
             * Construct a binary matrix of size rowsum.size() times colsum.size()
             * whose row and column sums match the given integer sequences.
             * @param rowsum Sequence of row sums.
             * @param colsum Sequence of column sums.
             */
            inline
            BinaryMatrix *
            realize(
                    const std::vector<int> &rowsum,
                    const std::vector<int> &colsum
            ) {
                return realize(&rowsum[0], &colsum[0], rowsum.size(), colsum.size());
            }


            /**
             * Construct a random binary matrix of size nrow times ncol
             * whose row and column sums match the given integer
             * sequence.
             * @param m Row and column sums.
             */
            inline
            BinaryMatrix *
            realize_random(const Instance& m) {

                const int nrow = m.getNumRows();
                const int ncol = m.getNumCols();

                /**************************************************************
                 * Running Time: O(nrow + ncol*log(ncol) + total),
                 * where total = sum(rowsum) = sum(colsum).
                 *************************************************************/

                /**************************************************************
                 * Check realizability.
                 *************************************************************/

                struct A {
                    int index;
                    int degree;
                };

                // sorted vector of index-degree-pairs
                A *cols = new A[ncol];
                for (int j = 0; j < ncol; j++)
                    cols[j] = {m.colindex[j], m.colsum[j]};

                A *rows = new A[nrow];
                for (int i = 0; i < nrow; i++)
                    rows[i] = {m.rowindex[i], m.rowsum[i]};

                // sort columns descendingly by degree
                std::sort(cols, cols + ncol, [](const A &a, const A &b) {
                    return a.degree > b.degree;
                });

                // sort rows descendingly by degree
                std::sort(rows, rows + nrow, [](const A &a, const A &b) {
                    return a.degree > b.degree;
                });

                // check range of degrees
                if (rows[0].degree > ncol || rows[nrow - 1].degree < 0
                    || cols[0].degree > nrow || cols[ncol - 1].degree < 0)
                    return nullptr;

                // create conjugate sequence of row sums
                int *rowsum_conjugated = new int[ncol];
                conjugate(rowsum_conjugated, &m.rowsum[0], ncol, nrow);

                int rowsum_conjugated_sum = 0;
                int colsum_sum = 0;
                for (int j = 0; j < ncol; j++) {
                    rowsum_conjugated_sum += rowsum_conjugated[j];
                    colsum_sum += cols[j].degree;
                    if (rowsum_conjugated_sum < colsum_sum) {
                        // sequence is not realizable
                        delete[] cols;
                        delete[] rows;
                        delete[] rowsum_conjugated;
                        return nullptr;
                    }
                }

                /******************************************************************
                 * Apply Ryser's algorithm with random vertex selection.
                 *****************************************************************/

                // matrix realization
                BinaryMatrix *bip = new BinaryMatrix(nrow, ncol);

                // RNG
                marathon::BasicRandom r;

                // determine first and last occurrence of each degree
                int *row_first = new int[ncol + 1];
                int *row_last = new int[ncol + 1];
                int *col_first = new int[nrow + 1];
                int *col_last = new int[nrow + 1];

                // initialize row_first and row_last
                for (int k = 0; k < ncol + 1; k++)
                    row_first[k] = row_last[k] = nrow;

                // initialize col_first and col_last
                for (int k = 0; k < nrow + 1; k++)
                    col_first[k] = col_last[k] = ncol;

                // row_first[k] = min { i : rows[i].degree == k }
                for (int i = nrow - 1; i >= 0; i--)
                    row_first[rows[i].degree] = i;

                // row_last[k] = max { i : rows[i].degree == k }
                for (int i = 0; i < nrow; i++)
                    row_last[rows[i].degree] = i;

                // col_first[k] = min { j : column[j].degree == k }
                for (int j = ncol - 1; j >= 0; j--)
                    col_first[cols[j].degree] = j;

                // col_last[k] = max { j : column[j].degree == k }
                for (int j = 0; j < ncol; j++)
                    col_last[cols[j].degree] = j;

                // while there is a row or column with positive degree
                while (rows[0].degree > 0 || cols[0].degree > 0) {

                    // determine number of rows and columns with positive degree
                    int num_rows_left = row_first[0];
                    int num_cols_left = col_first[0];

                    // randomly select a row or column with positive degree
                    int x = r.nextInt(num_rows_left + num_cols_left);

                    // if we selected a row
                    if (x < num_rows_left) {

                        assert(rows[x].degree > 0);

                        // calculate auxiliary variables
                        int i = rows[x].index;
                        int d = rows[x].degree;
                        int e = cols[d - 1].degree;
                        int f = col_first[e];
                        int l = col_last[e];
                        int n = l - f + 1;
                        int k = d - f;

                        /**
                         * Example:
                         *                 0  1  2  3
                         * rows.degree = [ 5, 3, 2, 1 ]
                         *                 x
                         *
                         * x = 0
                         * d = 5
                         *
                         *                 0  1  2  3  4  5, 6, 7
                         * cols.degree = [ 3, 2, 1, 1, 1, 1, 1, 0 ]
                         *                       f        d  l
                         * e = 1
                         * n = 6-2+1 = 5
                         * k = 5-2   = 3
                         */

                        // the first f columns will be used without having a choice
                        for (int q = f - 1; q >= 0; q--) {

                            // determine matrix position and set entry to one
                            int j = cols[q].index;
                            bip->set(i, j, 1);

                            // reduce degree
                            int m = cols[q].degree;
                            cols[q].degree--;

                            // update first and last
                            if (col_first[m] == col_last[m]) {
                                col_first[m] = col_last[m] = ncol;
                            } else {
                                col_last[m] = q - 1;
                            }
                            if (col_first[m - 1] == ncol) {
                                col_first[m - 1] = col_last[m - 1] = q;
                            } else {
                                col_first[m - 1] = q;
                            }
                        }

                        if (k < n) {
                            // randomly shuffle the columns cols[f..f+n]
                            r.shuffle<A>(cols + f, n);
                        }

                        // select k random columns out of n
                        for (int q = l; q >= l - k + 1; q--) {

                            // insert one
                            int j = cols[q].index;
                            bip->set(i, j, 1);
                            cols[q].degree--;

                            // update first and last
                            if (col_first[e] == col_last[e]) {
                                col_first[e] = col_last[e] = ncol;
                            } else {
                                col_last[e] = q - 1;
                            }
                            if (col_first[e - 1] == ncol) {
                                col_first[e - 1] = col_last[e - 1] = q;
                            } else {
                                col_first[e - 1] = q;
                            }
                        }


                        // restore descending order of rows
                        while (rows[x].degree > 0) {

                            // determine last row with same degree
                            int l = row_last[rows[x].degree];

                            // switch to position l
                            std::swap(rows[x], rows[l]);

                            // update first and last
                            if (row_first[rows[x].degree] == row_last[rows[x].degree]) {
                                row_first[rows[x].degree] = row_last[rows[x].degree] = nrow;
                            } else {
                                row_last[rows[x].degree]--;
                            }

                            // reduce d to minimal allowed value
                            if (l == nrow - 1) {
                                rows[l].degree = 0;
                                row_first[0] = row_last[0] = l;
                            } else {
                                rows[l].degree = rows[l + 1].degree;
                                row_first[rows[l].degree] = l;
                            }

                            x = l;
                        }


                    } else { // if we selected a column

                        // switch the roles of rows and columns
                        x -= num_rows_left;

                        assert(cols[x].degree > 0);

                        // calculate auxiliary variables
                        int i = cols[x].index;
                        int d = cols[x].degree;
                        int e = rows[d - 1].degree;
                        int f = row_first[e];
                        int l = row_last[e];
                        int n = l - f + 1;
                        int k = d - f;

                        /**
                         * Example:
                         *                 0  1  2  3
                         * cols.degree = [ 5, 3, 2, 1 ]
                         *                 x
                         *
                         * x = 0
                         * d = 5
                         *
                         *                 0  1  2  3  4  5, 6, 7
                         * rows.degree = [ 3, 2, 1, 1, 1, 1, 1, 0 ]
                         *                       f        d  l
                         * e = 1
                         * n = 6-2+1 = 5
                         * k = 5-2   = 3
                         */

                        // the first f rows will be used without having a choice
                        for (int q = f - 1; q >= 0; q--) {

                            // determine matrix position and set entry to one
                            int j = rows[q].index;
                            bip->set(j, i, 1);

                            // reduce degree
                            int m = rows[q].degree;
                            rows[q].degree--;

                            // update first and last
                            if (row_first[m] == row_last[m]) {
                                row_first[m] = row_last[m] = nrow;
                            } else {
                                row_last[m] = q - 1;
                            }
                            if (row_first[m - 1] == nrow) {
                                row_first[m - 1] = row_last[m - 1] = q;
                            } else {
                                row_first[m - 1] = q;
                            }
                        }

                        if (k < n) {
                            // randomly shuffle the rows rows[f..f+n]
                            r.shuffle<A>(rows + f, n);
                        }

                        // select k random rows out of n
                        for (int q = l; q >= l - k + 1; q--) {

                            // insert one
                            int j = rows[q].index;
                            bip->set(j, i, 1);
                            rows[q].degree--;

                            // update first and last
                            if (row_first[e] == row_last[e]) {
                                row_first[e] = row_last[e] = nrow;
                            } else {
                                row_last[e] = q - 1;
                            }
                            if (row_first[e - 1] == nrow) {
                                row_first[e - 1] = row_last[e - 1] = q;
                            } else {
                                row_first[e - 1] = q;
                            }
                        }


                        // restore descending order of columns
                        while (cols[x].degree > 0) {

                            // determine last column with same degree
                            int l = col_last[cols[x].degree];

                            // switch to position l
                            std::swap(cols[x], cols[l]);

                            // update first and last
                            if (col_first[cols[x].degree] == col_last[cols[x].degree]) {
                                col_first[cols[x].degree] = col_last[cols[x].degree] = ncol;
                            } else {
                                col_last[cols[x].degree]--;
                            }

                            // reduce d to minimal allowed value
                            if (l == ncol - 1) {
                                cols[l].degree = 0;
                                col_first[0] = col_last[0] = l;
                            } else {
                                cols[l].degree = cols[l + 1].degree;
                                col_first[cols[l].degree] = l;
                            }

                            x = l;
                        }

                    }
                }

                delete[] rows;
                delete[] cols;
                delete[] rowsum_conjugated;
                delete[] row_first;
                delete[] row_last;
                delete[] col_first;
                delete[] col_last;

                return bip;
            }

            /**
             * Construct a random binary matrix of size nrow times ncol
             * whose row and column sums match the given integer
             * sequence.
             * @param rowsum Sequence of row sums.
             * @param colsum Sequence of column sums.
             * @param nrow The number of rows.
             * @param ncol The number of columns.
             */
            inline
            BinaryMatrix *
            realize_random(
                    const int *rowsum,
                    const int *colsum,
                    const int nrow,
                    const int ncol
            ) {
                return realize_random(Instance(rowsum, colsum, nrow, ncol));
            }
        }
    }
}

#endif //MARATHON_FIXED_MARGIN_REALIZATION_H
