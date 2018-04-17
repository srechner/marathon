/*
 * Created on: Oct 06, 2017
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

#ifndef MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_ENUMERATE_H
#define MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_ENUMERATE_H

// marathon includes
#include "marathon/enumerate.h"
#include "marathon/combination_generator.h"
#include "marathon/binary_matrix/interval_margin/count.h"
#include "marathon/binary_matrix/fixed_margin/random_generator_exact.h"

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * An enumerator engine that produces all binary matrices whose row and column
             * sums are bounded by prescribed lower and upper bounds.
             */
            class Enumerator : private Counter, public marathon::Enumerator {

                BinaryMatrix _bin;       // binary matrix
                int *_ascending;         // array of ascending integers starting at zero

                /**
				 * Recursively sample a random binary matrices whose row and column sums
				 * lie in the prescribed intervals.
				 * @param rowsum_lower Sequence of lower row sums.
				 * @param rowsum_upper Sequence of upper row sums.
				 * @param rowsum_index Original order of rows.
				 * @param rowsum_upper_conj Conjugate sequence of upper row sums.
				 * @param nrow Number of rows.
				 * @param colsum_lower Sequence of lower column sums.
				 * @param colsum_upper Sequence of upper column sums.
				 * @param colsum_index Original order of columns.
				 * @param colsum_upper_conj Conjugate sequence of upper column sums.
				 * @param ncol Number of columns
				 * @param column_groups List of groups. Each group contains a number of columns
				 * with identical values of (lower, upper).
				 * @param group_id Index of the group considered in this step of recursion.
				 * @param aggr_in_row Amount this has been so far distributed in the current row.
                 * @param f Function that is evaluated for each valid binary matrix.
				 */
                void enumerate_recursive(
                        const int *rowsum_lower,
                        const int *rowsum_upper,
                        const int *rowsum_index,
                        int *rowsum_upper_conj,
                        const int nrow,
                        int *colsum_lower,
                        int *colsum_upper,
                        int *colsum_index,
                        int *colsum_upper_conj,
                        const int ncol,
                        const std::vector<Group> &column_groups,
                        const int group_id,
                        const int aggr_in_row,
                        const std::function<void(const State &)> &f
                ) {

                    // if the whole matrix has completely been processed
                    if (nrow == 0) {
                        if (colsum_lower[0] == 0)
                            f(_bin);   // process matrix
                        return;     // we are done
                    }

                    // if the current row has been completely processed.
                    if (group_id == column_groups.size()) {

                        assert(aggr_in_row >= rowsum_lower[0] && aggr_in_row <= rowsum_upper[0]);

                        // prepare columns for next step
                        std::vector<A> new_columns(ncol);
                        for (int j = 0; j < ncol; j++) {
                            new_columns[j].index = colsum_index[j];
                            new_columns[j].lower = colsum_lower[j];
                            new_columns[j].upper = std::min(nrow - 1, colsum_upper[j]);
                        }
                        std::sort(new_columns.begin(), new_columns.end(), [](const A &a, const A &b) {
                            return a.lower == b.lower ? a.upper > b.upper : a.lower > b.lower;
                        });
                        int *colsum_lower_new = new int[ncol];
                        int *colsum_upper_new = new int[ncol];
                        int *colsum_index_new = new int[ncol];
                        for (int j = 0; j < ncol; j++) {
                            colsum_lower_new[j] = new_columns[j].lower;
                            colsum_upper_new[j] = new_columns[j].upper;
                            colsum_index_new[j] = new_columns[j].index;
                        }
                        int *colsum_upper_conj_new = new int[nrow - 1];
                        conjugate(colsum_upper_conj_new, colsum_upper_new, nrow - 1, ncol);

                        // modify rowsum_upper_conj
                        for (int i = 0; i < rowsum_upper[0]; i++)
                            rowsum_upper_conj[i]--;

                        // if new sequence is bi-graphical
                        if (isDominating(colsum_upper_conj_new, rowsum_lower + 1, nrow - 1) &&
                            isDominating(rowsum_upper_conj, colsum_lower_new, ncol)) {

                            // re-group columns
                            std::vector<Group> new_groups;
                            group_columns(new_columns, new_groups);

                            // recursively process the nrow-1 times ncol submatrix.
                            enumerate_recursive(
                                    rowsum_lower + 1,
                                    rowsum_upper + 1,
                                    rowsum_index + 1,
                                    rowsum_upper_conj,
                                    nrow - 1,
                                    colsum_lower_new,
                                    colsum_upper_new,
                                    colsum_index_new,
                                    colsum_upper_conj_new,
                                    ncol,
                                    new_groups,
                                    0,
                                    0,
                                    f
                            );

                        }

                        // undo modification
                        for (int i = 0; i < rowsum_upper[0]; i++)
                            rowsum_upper_conj[i]++;

                        delete[] colsum_lower_new;
                        delete[] colsum_upper_new;
                        delete[] colsum_index_new;
                        delete[] colsum_upper_conj_new;

                        return;
                    }

                    // we backtrack all choices of distributing ones in the current group

                    const int gl = column_groups[group_id].lower;
                    const int gu = column_groups[group_id].upper;

                    // n is the number of columns in the current group
                    const int n = column_groups[group_id].size;

                    // a is the maximal number of ones that can be distributed without violating the
                    // upper bound on the current row
                    const int a = rowsum_upper[0] - aggr_in_row;

                    // b is the minimal number of ones that can be distributed to this group such
                    // that we may reach a row total of at least rowsum_lower[0]
                    const int b = rowsum_lower[0] - aggr_in_row - column_groups[group_id].getFirst();

                    const int lower = std::max(0, b);
                    const int upper = std::min(n, a);

                    // first and last indices of current group
                    const int first = column_groups[group_id].getFirst();
                    const int last = column_groups[group_id].getLast();

                    // auxiliary array used for backtracking combinations
                    int *select = new int[std::max(upper, 0)];

                    // Each value of k in the range [lower,upper] results in a set of matrices.
                    // for each valid choice
                    for (int k = lower; k <= upper; k++) {

                        // apply choice: distribute k ones in the current group
                        colsum_upper_conj[gu - 1] -= k;
                        for (int l = 0, j = last; l < k; l++, j--) {
                            colsum_lower[j] = std::max(0, gl - 1);
                            colsum_upper[j] = gu - 1;
                        }

                        // simulate all ways to select k columns out of n
                        CombinationGenerator<int> cg(_ascending, select, n, k);
                        do {

                            // for each selected column
                            for (int l = 0; l < k; l++) {

                                // index of selected column in non-increasing order
                                const int j = first + select[k - l - 1];

                                // set matrix entry
                                _bin.set(rowsum_index[0], colsum_index[j], 1);

                                // restore non-increasing order of _colsum_lower
                                std::swap(colsum_index[last - l], colsum_index[j]);
                            }

                            // recursively fill the remaining matrix
                            enumerate_recursive(
                                    rowsum_lower,
                                    rowsum_upper,
                                    rowsum_index,
                                    rowsum_upper_conj,
                                    nrow,
                                    colsum_lower,
                                    colsum_upper,
                                    colsum_index,
                                    colsum_upper_conj,
                                    ncol,
                                    column_groups,
                                    group_id + 1,
                                    aggr_in_row + k,
                                    f
                            );

                            // undo column selection
                            for (int l = 0; l < k; l++) {

                                // index of selected column in non-increasing order
                                const int j = first + select[l];

                                // restore original order
                                std::swap(colsum_index[j], colsum_index[last - k + l + 1]);

                                // reset matrix entry
                                _bin.set(rowsum_index[0], colsum_index[j], 0);
                            }

                        } while (cg.next());

                        // undo modification
                        colsum_upper_conj[gu - 1] += k;
                        for (int l = 0, j = last; l < k; l++, j--) {
                            colsum_lower[j] = gl;
                            colsum_upper[j] = gu;
                        }
                    }

                    delete[] select;
                }

            public:

                /**
				 * Create a enumerator for binary matrices whose row and column
				 * sums lie in the prescribed intervals.
				 * @param inst Four-tuple of integer vectors.
				 */
                Enumerator(Instance inst)
                        : Counter(std::move(inst)) {

                    _bin = BinaryMatrix(_nrow, _ncol);

                    // create ascending array of integers
                    _ascending = new int[_ncol];
                    for (int j = 0; j < _ncol; j++)
                        _ascending[j] = j;
                }

                /**
				 * Create an enumerator for binary matrices whose row and column
				 * sums lie in the prescribed intervals.
				 * @param rowsum_lower Lower bounds on the row sums.
				 * @param rowsum_upper Upper bounds on the row sums.
				 * @param colsum_lower Lower bounds on the column sums.
				 * @param colsum_upper Upper bounds on the column sums.
				 */
                Enumerator(
                        const std::vector<int> &rowsum_lower,
                        const std::vector<int> &rowsum_upper,
                        const std::vector<int> &colsum_lower,
                        const std::vector<int> &colsum_upper
                ) : Enumerator(
                        Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper)) {

                }

                /**
				 * Create a enumerator for binary matrices whose row and column
				 * sums lie in the prescribed intervals.
				 * @param rowsum_lower Lower bounds on the row sums.
				 * @param rowsum_upper Upper bounds on the row sums.
				 * @param colsum_lower Lower bounds on the column sums.
				 * @param colsum_upper Upper bounds on the column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns.
				 */
                Enumerator(
                        const int *rowsum_lower,
                        const int *rowsum_upper,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        size_t nrow,
                        size_t ncol
                ) : Enumerator(
                        std::vector<int>(rowsum_lower, rowsum_lower + nrow),
                        std::vector<int>(rowsum_upper, rowsum_upper + nrow),
                        std::vector<int>(colsum_lower, colsum_lower + ncol),
                        std::vector<int>(colsum_upper, colsum_upper + ncol)
                ) {

                }

                virtual ~Enumerator() override {
                    delete _ascending;
                }

                /**
                 * Enumerate all binary matrices of size nrow times ncol whose
                 * row and column sums are bounded by given integer vectors.
                 * @param f Function that is evaluated for each binary matrix.
                 */
                void enumerate(const std::function<void(const State &)> f) override {

                    int *rowsum_lower = new int[_nrow];
                    int *rowsum_upper = new int[_nrow];
                    int *rowsum_index = new int[_nrow];
                    int *colsum_lower = new int[_ncol];
                    int *colsum_upper = new int[_ncol];
                    int *colsum_index = new int[_ncol];

                    // create row sequence sorted descendingly by lower bound
                    std::vector<A> rows(_nrow);
                    for (int i = 0; i < _nrow; i++) {
                        rows[i].index = i;
                        rows[i].lower = std::max(_seq._rowsum_lower[i], 0);
                        rows[i].upper = std::min(_seq._rowsum_upper[i], _ncol);
                    }
                    std::sort(rows.begin(), rows.end(), [](const A &a, const A &b) {
                        return a.lower == b.lower ? a.upper > b.upper : a.lower > b.lower;
                    });

                    for (int i = 0; i < _nrow; i++) {
                        rowsum_lower[i] = rows[i].lower;
                        rowsum_upper[i] = rows[i].upper;
                        rowsum_index[i] = rows[i].index;
                    }

                    // create column sequence sorted descendingly by lower bound
                    std::vector<A> columns(_ncol);
                    for (int i = 0; i < _ncol; i++) {
                        columns[i].index = i;
                        columns[i].lower = std::max(_seq._colsum_lower[i], 0);
                        columns[i].upper = std::min(_seq._colsum_upper[i], _nrow);
                    }
                    std::sort(columns.begin(), columns.end(), [](const A &a, const A &b) {
                        return a.lower == b.lower ? a.upper > b.upper : a.lower > b.lower;
                    });

                    for (int j = 0; j < _ncol; j++) {
                        colsum_lower[j] = columns[j].lower;
                        colsum_upper[j] = columns[j].upper;
                        colsum_index[j] = columns[j].index;
                    }

                    // create conjuate sequences of upper column sum
                    int *colsum_upper_conj = new int[_nrow];
                    conjugate(colsum_upper_conj, colsum_upper, _nrow, _ncol);

                    // create conjugate sequence of upper row sums
                    int *rowsum_upper_conj = new int[_ncol];
                    conjugate(rowsum_upper_conj, rowsum_upper, _ncol, _nrow);

                    // group columns by identical pairs of lower and upper bounds
                    std::vector<Group> groups;
                    group_columns(columns, groups);

                    // clear variables
                    _bin.clear();

                    // select the binary matrix with number target by traversing the enumeration tree
                    enumerate_recursive(
                            rowsum_lower,
                            rowsum_upper,
                            rowsum_index,
                            rowsum_upper_conj,
                            _nrow,
                            colsum_lower,
                            colsum_upper,
                            colsum_index,
                            colsum_upper_conj,
                            _ncol,
                            groups,
                            0,
                            0,
                            f);

                    delete[] rowsum_lower;
                    delete[] rowsum_upper;
                    delete[] rowsum_index;
                    delete[] rowsum_upper_conj;
                    delete[] colsum_lower;
                    delete[] colsum_upper;
                    delete[] colsum_index;
                    delete[] colsum_upper_conj;
                }
            };
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_RANDOM_GENERATOR_EXACT_H
