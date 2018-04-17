/*
 * Created on: Sep 20, 2016
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

#ifndef MARATHON_BINARY_MATRIX_FIXED_MARGIN_COUNT_H
#define MARATHON_BINARY_MATRIX_FIXED_MARGIN_COUNT_H

#include <numeric>
#include <algorithm>
#include <mutex>
#include <boost/unordered_map.hpp>

// marathon includes
#include "marathon/count.h"
#include "decompose.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Class for counting the number of binary matrices with prescribed margins.
             */
            class Counter : public marathon::Counter {

            protected:

                // invariant class members
                const Instance _inst;            // instance description
                const size_t _nrow;            // number of rows
                const size_t _ncol;            // number of columns
                std::vector<int> _rowsum;    // sequence of row sums
                std::vector<int> _colsum_conj;  // conjugated sequence of column sums
                bool _realizable;            // is sequence realizable at all?

                // class members that will be altered during the execution
                boost::unordered_map<marathon::Integer, marathon::Integer> tmp;

                /**
                 * Store the value in the hash table for later use.
                 * The sequence pair (rowsum, colsum_conj) is used as key.
                 * @param rowsum Sequence of row sums.
                 * @param nrow Number of rows.
                 * @param x Number of matrices with prescribed row sums and column sums.
                 */
                void
                storeToTable(
                        const int *rowsum,
                        const int *colsum_conj,
                        const int nrow,
                        const marathon::Integer &value
                ) {
                    // compress key
                    const int base = (int) _ncol + 1;
                    const marathon::Integer key1 = evaluate_polynomial(rowsum, nrow, base);
                    const marathon::Integer key = evaluate_polynomial(colsum_conj, nrow, base, key1);

                    // store solution in hash map
                    tmp[key] = value;
                }

                /**
                 * Load a precomputed value from hash table.
                 * The sequence pair (rowsum, colsum_conj) is used as key.
                 * @param rowsum Sequence of row sums.
                 * @param nrow Number of rows.
                 * @param res (Output paramter) The value of the precomputed element will be stored in this variable.
                 * @return True, if a precomputed value for the key is already stored in table or false, otherwise.
                 */
                bool loadFromTable(
                        const int *rowsum,
                        const int *colsum_conj,
                        const int nrow,
                        marathon::Integer &res
                ) {

                    // compress key
                    const int base = (int) _ncol + 1;
                    const marathon::Integer key1 = evaluate_polynomial(rowsum, nrow, base);
                    const marathon::Integer key = evaluate_polynomial(colsum_conj, nrow, base, key1);

                    // crawl table for result
                    const auto it = tmp.find(key);
                    const auto end = tmp.end();

                    // solution found?
                    if (it != end) {
                        // use precomputed value
                        res = it->second;
                        return true;
                    }

                    // the result is not in table
                    return false;
                }

                /**
                 * (Help Function) Count the number of binary matrices of size nrow times ncol which agree with the given
                 * rowsums and conjugated column sums. Works by recursive backtracking with clever cutting techniques.
                 * @param rowsum Integer array with row sums that are to be realized.
                 * @param colsum_conj Conjugate Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param x Consider columns i with colsum[i]=x.
                 * @param k_aggr Number of ones already distributed in current row.
                 * @param rowsum_aggr Aggregated row sums up to the current row.
                 * @param colsum_conj_aggr Aggregated conjugated column sums.

                 */
                marathon::Integer
                count_recursive(
                        const int *rowsum,
                        int *colsum_conj,
                        const int nrow,
                        const int x,
                        const int k_aggr,
                        const int rowsum_aggr,
                        const int colsum_conj_aggr
                ) {

                    // if the current row has been completely processed
                    if (k_aggr == rowsum[0]) {

                        // if the whole matrix has been processed
                        if (nrow == 1) {
                            return 1;
                        } else {

                            // the matrix has not yet been completely processed
                            // recursively process the nrow-1 times ncol submatrix

                            // try to find whether the solution of the subproblem
                            // is already stored in the hash table
                            marathon::Integer res;
                            bool found = loadFromTable(&rowsum[1], colsum_conj, nrow - 1, res);

                            // if solution is not already known
                            if (!found) {

                                // recursively count the nrow-1 times ncol submatrix
                                res = count_recursive(&rowsum[1], colsum_conj, nrow - 1, 1, 0, 0, 0);

                                // store solution in hash map
                                storeToTable(&rowsum[1], colsum_conj, nrow - 1, res);
                            }

                            return res;
                        }
                    } else {

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

                        // for each valid choice of distributing ones to the set X
                        marathon::Integer res = 0;
                        for (int k = lower; k <= upper; k++) {

                            // virtually distribute k ones in the current row
                            colsum_conj[x - 1] -= k;

                            // recursively count the number of matrices that result from this choice
                            marathon::Integer val = count_recursive(
                                    rowsum,
                                    colsum_conj,
                                    nrow,
                                    x + 1,
                                    k_aggr + k,
                                    rowsum_aggr + rowsum[x],
                                    colsum_conj_aggr + colsum_conj[x - 1] + k
                            );

                            // multiply by the number of ways to distribute k ones to n entries
                            val *= marathon::binom(n, k);

                            // increase the total number of matrices by the matrices resulting from this choice
                            res += val;

                            // undo modification
                            colsum_conj[x - 1] += k;
                        }

                        return res;
                    }
                }

            public:

                /**
                 * Create a Counter for the number of binary matrices
                 * with exactly prescribed margins.
                 * @param seq Margin of row and column sums.
                 */
                Counter(Instance seq) :
                        _inst(std::move(seq)),
                        _nrow(_inst.getNumRows()),
                        _ncol(_inst.getNumCols()) {

                    /**
                     * Idea by Miller and Harrison:
                     * Exact Sampling and Counting for Fixed-Margin Matrices.
                     * The Annals of Statistics, 2013.
                     */

                    // check realizability
                    _realizable = isRealizable(_inst);

                    if (!_realizable)
                        return;

                    // descendingly sort row sums
                    _rowsum = _inst._rowsum;
                    std::sort(_rowsum.begin(), _rowsum.end(), [](int a, int b) { return a > b; });

                    // add two dummy elements at the end of the vector
                    _rowsum.push_back(0);
                    _rowsum.push_back(0);

                    // conjugate column sums
                    _colsum_conj = conjugate(_inst._colsum, _nrow + 1);
                }

                /**
                 * Create a Counter for the number of binary matrices
                 * with exactly prescribed margins.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param ncol Sequence of columns.
                 */
                Counter(
                        const int *rowsum,
                        const int *colsum,
                        size_t nrow,
                        size_t ncol
                ) : Counter(Instance(rowsum, colsum, nrow, ncol)) {

                }

                /**
                 * Create a Counter for the number of binary matrices
                 * with exactly prescribed margins.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 */
                Counter(
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) : Counter(&rowsum[0], &colsum[0], rowsum.size(), colsum.size()) {

                }


                /**
                 * Count the number of binary matrices of size nrow times ncol whose
                 * row and column sums match the given integer vectors.
                 * @return The number of binary matrices that agree with the given row and column sums.
                 */
                marathon::Integer count() override {

                    marathon::Integer res(0);

                    if (!_realizable)
                        return res;

                    // try to find solution the current subproblem in the hash table
                    bool found = loadFromTable(&_rowsum[0], &_colsum_conj[0], _nrow, res);

                    // solution is not already known
                    if (!found) {

                        // create new working arrays
                        int *rowsum = new int[_nrow + 2];
                        int *colsum_conj = new int[_nrow + 1];
                        memcpy(rowsum, &_rowsum[0], (_nrow + 2) * sizeof(int));
                        memcpy(colsum_conj, &_colsum_conj[0], (_nrow + 1) * sizeof(int));

                        // start calculations
                        res = count_recursive(rowsum, colsum_conj, _nrow, 1, 0, 0, 0);
                        assert(res > 0);

                        // store solution in hash map
                        storeToTable(&_rowsum[0], &_colsum_conj[0], _nrow, res);

                        // clean up
                        delete[] rowsum;
                        delete[] colsum_conj;
                    }

                    return res;
                }
            };
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_FIXED_MARGIN_COUNT_H
