/*
 * Created on: Nov 02, 2016
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


#ifndef MARATHON_BINARY_MATRIX_FIXED_MARGIN_ENUMERATE_H
#define MARATHON_BINARY_MATRIX_FIXED_MARGIN_ENUMERATE_H

#include <cstdlib>
#include <vector>
#include <numeric>
#include <algorithm>

// marathon includes
#include "marathon/enumerate.h"
#include "marathon/binary_matrix/binary_matrix.h"
#include "marathon/binary_matrix/fixed_margin/instance.h"
#include "marathon/combination_generator.h"
#include "count.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Enumerator engine for the set of binary matrices with prescribed margins.
             */
            class Enumerator : private Counter, public marathon::Enumerator {

            protected:

                // auxiliary data structure
                struct A {
                    size_t index;
                    int value;
                };

                // invariant class members
                size_t *_roworder;          // permutation of row indices
                A *_columns;                // sorted column margins and indices
                int *_colsum_first;         // first position of each value of colsum
                int *_colsum_last;          // last position of each value of colsum
                int *_ascending;            // array of ascending integers starting at zero

                // class members modified during execution
                BinaryMatrix _bin;          // binary matrix modified during execution


                /**
                 * (Help Function) Enumerate all binary matrices of size nrow times ncol which agree with the given
                 * rowsums and conjugated column sums. Works by recursive backtracking with clever cutting techniques.
                 * @param rowsum Integer array with row sums that are to be realized.
                 * @param colsum_conj Conjugate Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param x Consider columns i with colsum[i]=x.
                 * @param k_aggr Number of ones already distributed in current row.
                 * @param rowsum_aggr Aggregated row sums up to the current row.
                 * @param colsum_conj_aggr Aggregated conjugated column sums.
                 */
                void enumerate_recursive(
                        const int *rowsum,
                        const size_t *roworder,
                        A *columns,
                        int *colsum_conj,
                        int *colsum_first,
                        int *colsum_last,
                        const int nrow,
                        const int x,
                        const int k_aggr,
                        const int rowsum_aggr,
                        const int colsum_conj_aggr,
                        const std::function<void(const State &)> &f
                ) {

                    // if the current row has been completely processed
                    if (k_aggr == rowsum[0]) {

                        // if the whole matrix has been processed
                        if (nrow == 1) {
                            // we are done
                            f(_bin);
                            return;
                        } else {

                            // the matrix has not yet been completely processed
                            // recursively process the nrow-1 times ncol submatrix
                            enumerate_recursive(
                                    &rowsum[1], &roworder[1],
                                    columns, colsum_conj, colsum_first, colsum_last,
                                    nrow - 1, 1, 0, 0, 0, f);
                            return;
                        }
                    }

                    // Consider the set of entries X = {j : colsum[j] = x}.

                    // n is the number of entries j such that colsum[j] = x.
                    const int n = colsum_conj[x - 1] - colsum_conj[x];

                    // b is the minimal number of entries in X that can be set to one
                    // so that we may reach a total of rowsum[0]
                    const int b = rowsum[0] - k_aggr - colsum_conj[x];

                    // s is the maximal amount of which colsum_conj[x-1] can be decreased
                    // without violating the Gale-Ryser Condition
                    const int s = (colsum_conj_aggr + colsum_conj[x - 1]) - (rowsum_aggr + rowsum[x]) - k_aggr;

                    // upper and lower bounds
                    const int lower = std::max(0, b);
                    const int upper = std::min({rowsum[0], n, s});

                    // create auxiliary array
                    int *selection = new int[std::max(upper, 0)];

                    // Each value of k in the range [lower,upper] results in a set of matrices.
                    for (int k = lower; k <= upper; k++) {

                        /*
                         * Distribute k ones to the n columns of the current row that currently have value x.
                         * First element with value x is colsum_first[x].
                         * Last element with value x is colsum_last[x].
                         */

                        colsum_conj[x - 1] -= k;

                        // smallest column index with column sum x
                        const int first = colsum_first[x];

                        // simulate all combinations of distributing k elements to n slots
                        marathon::CombinationGenerator<int> cg(_ascending, selection, n, k);
                        do {

                            // for each selected column
                            for (int i = k - 1; i >= 0; i--) {

                                // j is the index of the selected column in decreasing order
                                const int j = first + selection[i];
                                assert(columns[j].value == x);

                                // set matrix element t one (restore original order)
                                _bin.set(roworder[0], columns[j].index, 1);
                                columns[j].value--;

                                // to maintain the non-increasing order of columns, switch column j with
                                // the last column of value x
                                const int last = colsum_last[x];    // last column index with column sum x
                                assert(last != -1);

                                std::swap(columns[j], columns[last]);

                                // update last array
                                if (first < last)
                                    colsum_last[x]--;
                                else
                                    colsum_last[x] = -1;
                                if (colsum_last[x - 1] == -1)
                                    colsum_last[x - 1] = last;

                                // update first array
                                if (first == last)
                                    colsum_first[x] = -1;
                                colsum_first[x - 1] = last;
                            }

                            // recursively fill the remaining matrix
                            enumerate_recursive(
                                    rowsum,
                                    roworder,
                                    columns,
                                    colsum_conj,
                                    colsum_first,
                                    colsum_last,
                                    nrow,
                                    x + 1,
                                    k_aggr + k,
                                    rowsum_aggr + rowsum[x],
                                    colsum_conj_aggr + colsum_conj[x - 1] + k,
                                    f
                            );

                            // undo changes

                            // for each selected column in reverse order
                            for (int i = 0; i < k; i++) {

                                // j is the index of the previously selected column
                                const int j = first + selection[i];

                                // column j was previously stored at position k
                                const int k = colsum_first[x - 1];
                                assert(colsum_first[x - 1] != -1);
                                assert(columns[k].value == x - 1);

                                // reset matrix entry
                                _bin.set(roworder[0], columns[k].index, 0);
                                columns[k].value++;

                                // restore original order
                                std::swap(columns[j], columns[k]);

                                // update first and last array at position x-1
                                if (colsum_first[x - 1] == colsum_last[x - 1])
                                    colsum_first[x - 1] = colsum_last[x - 1] = -1;
                                else
                                    colsum_first[x - 1]++;

                                // update first and last array at position x
                                if (colsum_first[x] == -1)
                                    colsum_first[x] = colsum_last[x] = j;
                                else
                                    colsum_last[x]++;
                            }

                        } while (cg.next());


                        // undo changes
                        colsum_conj[x - 1] += k;
                    }

                    delete[] selection;
                }

            public:

                /**
                 * Create a Enumerator for the set of binary matrices
                 * with exactly prescribed margins.
                 * @param seq Margin of row and column sums.
                 */
                Enumerator(Instance inst) : Counter(std::move(inst)) {

                    /**
                     * Derived from Miller and Harrison:
                     * Exact Sampling and Counting for Fixed-Margin Matrices.
                     * The Annals of Statistics, 2013.
                     */

                    if (!_realizable)
                        return;

                    // create ascending order of rows
                    A *row = new A[_nrow];
                    for (int i = 0; i < _nrow; i++)
                        row[i] = {_inst._rowindex[i], _inst._rowsum[i]};
                    std::sort(row, row + _nrow, [](const A &a, const A &b) { return a.value > b.value; });
                    _roworder = new size_t[_nrow];
                    for (size_t i = 0; i < _nrow; i++)
                        _roworder[i] = row[i].index;
                    delete[] row;

                    // create descending order of columns
                    _columns = new A[_ncol];
                    for (int j = 0; j < _ncol; j++)
                        _columns[j] = {_inst._colindex[j], _inst._colsum[j]};
                    std::sort(_columns, _columns + _ncol,
                              [](const A &a, const A &b) { return a.value > b.value; });

                    // set up array of first and last occurrence of each value of colsum
                    _colsum_first = new int[_nrow + 1];
                    _colsum_last = new int[_nrow + 1];
                    for (int i = 0; i <= _nrow; i++) {
                        _colsum_first[i] = -1;    // initialize with dummy value
                        _colsum_last[i] = -1;     // initialize with dummy value
                    }
                    for (int j = _ncol - 1; j >= 0; j--)
                        _colsum_first[_columns[j].value] = j;
                    for (int j = 0; j < _ncol; j++)
                        _colsum_last[_columns[j].value] = j;

                    // create binary matrix object
                    _bin = BinaryMatrix(_nrow, _ncol);

                    // create ascending array of integers
                    _ascending = new int[_ncol];
                    for (int j = 0; j < _ncol; j++)
                        _ascending[j] = j;

                }


                virtual ~Enumerator() override {
                    if (!_realizable)
                        return;

                    delete[] _roworder;
                    delete[] _columns;
                    delete[] _colsum_first;
                    delete[] _colsum_last;
                    delete _ascending;
                }

                /**
                 * Enumerate all binary matrices of size nrow times ncol whose
                 * row and column sums match the given integer vectors.
                 * @param f Function that is evaluated for each binary matrix.
                 */
                void enumerate(const std::function<void(const State &)> f) override {

                    if (!_realizable)
                        return;

                    // reset variables
                    A *columns = new A[_ncol];
                    int *rowsum = new int[_nrow + 2];
                    int *colsum_first = new int[_nrow + 1];
                    int *colsum_last = new int[_nrow + 1];
                    int *colsum_conj = new int[_nrow + 1];
                    memcpy(columns, _columns, _ncol * sizeof(A));
                    memcpy(rowsum, &_rowsum[0], (_nrow + 2) * sizeof(int));
                    memcpy(colsum_first, _colsum_first, (_nrow + 1) * sizeof(int));
                    memcpy(colsum_last, _colsum_last, (_nrow + 1) * sizeof(int));
                    memcpy(colsum_conj, &_colsum_conj[0], (_nrow + 1) * sizeof(int));

                    // select the binary matrix with number target by traversing the enumeration tree
                    enumerate_recursive(
                            rowsum,
                            _roworder,
                            columns,
                            colsum_conj,
                            colsum_first,
                            colsum_last,
                            _nrow,
                            1,
                            0,
                            0,
                            0,
                            f
                    );

                    delete[] columns;
                    delete[] rowsum;
                    delete[] colsum_first;
                    delete[] colsum_last;
                    delete[] colsum_conj;
                }
            };
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_FIXED_MARGIN_ENUMERATE_H
