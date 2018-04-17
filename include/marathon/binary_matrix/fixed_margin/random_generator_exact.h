/*
 * Created on: May 17, 2017
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

#ifndef MARATHON_BINARY_MATRIX_FIXED_MARGIN_RANDOM_GENERATOR_EXACT_H
#define MARATHON_BINARY_MATRIX_FIXED_MARGIN_RANDOM_GENERATOR_EXACT_H

// marathon includes
#include "marathon/random_device.h"
#include "marathon/binary_matrix/random_generator.h"
#include "marathon/binary_matrix/fixed_margin/count.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * A generator for binary matrices with prescribed row sums
             * and column sums.
             * This generator DOES NOT use sequence decomposition.
             */
            class RandomGeneratorExact :
                    public marathon::binary_matrix::RandomGenerator,
                    private marathon::binary_matrix::fixed_margin::Counter {

                // auxiliary data structure
                struct A {
                    size_t index;
                    int value;
                };

                // invariant class members
                int *_roworder;             // permutation of row indices
                A *_columns;                // sorted column margins and indices
                int *_colsum_first;         // first position of each value of colsum
                int *_colsum_last;          // last position of each value of colsum

                // class members modified during execution
                RandomDevice _rg;           // RandomDevice Generator
                int *_aux;                  // auxiliary array
                Integer _mul;               // multiplier used to avoid division
                Integer _target;            // we want to select the matrix with this number
                Integer _num_matrices;      // number of matrices with prescribed margins
                BinaryMatrix *_bin;         // binary matrix created during execution

                /**
                 * (Help function) Recursively sample a binary matrix.
                 * @param nrow Number of rows of sub matrix.
                 * @param rowsum Sequence of row sums.
                 * @param roworder Original row indices.
                 * @param pos
                 * @param running_sum
                 * @param running_rowsum
                 * @param running_colsum_conj
                 */
                void exact_sample_recursive(
                        const int *rowsum,
                        const int *roworder,
                        A *columns,
                        int *colsum_conj,
                        int *colsum_first,
                        int *colsum_last,
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
                            // we are done
                            return;
                        } else {

                            // the matrix has not yet been completely processed
                            // recursively process the nrow-1 times ncol submatrix
                            exact_sample_recursive(
                                    &rowsum[1], &roworder[1],
                                    columns, colsum_conj, colsum_first, colsum_last,
                                    nrow - 1, 1, 0, 0, 0);
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

                    /*
                     * Each value of k in the range [lower,upper] results in a number of matrices val(k).
                     * We seek the smallest k such that (val(k_lower) + ... + val(k)) > target.
                     */
                    int k = lower;

                    // virtually distribute k ones in the current row
                    colsum_conj[x - 1] -= k;

                    // if there are multiple values of k to choose from
                    if (lower != upper) {

                        marathon::Integer skipped = 0;

                        // determine how many matrices can be generated by each choice of k
                        for (; k <= upper; k++) {

                            // count the number of matrices that result from this choice
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

                            //std::cout << " k=" << k << ": val = " << val << std::endl;

                            // smallest k found!
                            if ((skipped + val) * _mul > _target)
                                break;

                            // increase the total number of matrices by the matrices resulting from this choice
                            skipped += val;

                            // distribute an additional one to the current row
                            colsum_conj[x - 1]--;
                        }

                        // reduce t by the number of matrices that result from choices smaller than k
                        _target -= skipped * _mul;
                    }

                    assert(k <= upper);

                    // std::cout << "reduce " << k << " colums with colsum " << pos << " by one" << std::endl;

                    /*
                     * k is the smallest value such that (val(l) + ... + val(k)) > target.
                     * Distribute k ones to the n columns of the current row that currently
                     * have value pos.
                     * First element with value pos is colsum_first[pos].
                     * Last element with value pos is colsum_last[pos].
                     */

                    // randomly select k columns out of n
                    _rg.select(_aux, n, k);

                    // for each selected column
                    for (int i = k - 1; i >= 0; i--) {

                        const int first = colsum_first[x];
                        const int last = colsum_last[x];
                        assert(first != -1);
                        assert(last != -1);

                        // j is the index of the selected column in decreasing order
                        const int j = first + _aux[i];
                        assert(columns[j].value == x);

                        // set matrix element t one (restore original order)
                        _bin->set(roworder[0], columns[j].index, 1);
                        columns[j].value--;

                        // to keep the sorted order of columns, switch column j with
                        // the last column of value pos
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

                    /*
                     * The number of matrices in the next level of recursion is
                     * scaled up by binom(n,k), thus the matrix with number t in
                     * the current level of recursion corresponds to the matrix
                     * with number t / binom(n,k).
                     */
                    _mul *= marathon::binom(n, k);

                    // recursively fill the remaining matrix
                    exact_sample_recursive(
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
                            colsum_conj_aggr + colsum_conj[x - 1] + k
                    );
                }

            public:

                /**
                 * Create a generator that produces (exactly!) uniformly distributed
                 * random binary matrices whose row and column sums match prescribed sequences.
                 * (This generator does not use sequence decomposition.)
                 * @param Row and column sums.
                 */
                RandomGeneratorExact(Instance margin) : Counter(std::move(margin)) {

                    /**
                     * Idea by Miller and Harrison:
                     * Exact Sampling and Counting for Fixed-Margin Matrices.
                     * The Annals of Statistics, 2013.
                     */

                    if (!_realizable)
                        return;

                    // count the number binary matrices
                    _num_matrices = count();
                    assert(_num_matrices > 0);

                    // create ascending order of rows
                    A *row = new A[_nrow];
                    for (int i = 0; i < _nrow; i++)
                        row[i] = {_inst._rowindex[i], _inst._rowsum[i]};
                    std::sort(row, row + _nrow, [](const A &a, const A &b) { return a.value > b.value; });
                    _roworder = new int[_nrow];
                    for (int i = 0; i < _nrow; i++)
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

                    // auxiliary array
                    _aux = new int[_ncol + 1];

                    // create binary matrix object
                    _bin = new BinaryMatrix(_nrow, _ncol);

                }

                /**
                 * Create a generator that produces (exactly!) uniformly distributed
                 * random binary matrices whose row and column sums match prescribed sequences.
                 * (This generator does not use sequence decomposition.)
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns
                 */
                RandomGeneratorExact(
                        const int *rowsum,
                        const int *colsum,
                        size_t nrow,
                        size_t ncol
                ) : RandomGeneratorExact(Instance(rowsum, colsum, nrow, ncol)) {


                }

                /**
                 * Create a generator that produces (exactly!) uniformly distributed
                 * random binary matrices whose row and column sums match prescribed sequences.
                 * (This generator does not use sequence decomposition.)
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 */
                RandomGeneratorExact(
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) : RandomGeneratorExact(&rowsum[0], &colsum[0], rowsum.size(), colsum.size()) {

                }

                virtual ~RandomGeneratorExact() {

                    if (!_realizable)
                        return;

                    delete[] _roworder;
                    delete[] _columns;
                    delete[] _colsum_first;
                    delete[] _colsum_last;
                    delete[] _aux;
                    delete _bin;
                }

                /**
                 * Return a binary matrix of size nrow times ncol whose
                 * row and column sums match the given integer vectors.
                 * The sample is distributed (exactly!) uniform from the set
                 * of all binary matrices with the prescribed row sums and
                 * column sums.
                 * @return Random binary matrix.
                 */
                const BinaryMatrix &next() override {

                    if (!_realizable)
                        throw std::runtime_error("Error! Vector pair not bi-graphical!");

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
                    _bin->clear();
                    _mul = 1;

                    // select a random number in range [0,num_matrices)
                    _target = _rg.nextInt(_num_matrices);

                    // select the binary matrix with number target by traversing the enumeration tree
                    exact_sample_recursive(
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
                            0
                    );

                    delete[] columns;
                    delete[] rowsum;
                    delete[] colsum_first;
                    delete[] colsum_last;
                    delete[] colsum_conj;

                    return *_bin;
                }


                /**
                 * Create an independent copy of an existing random generator.
                 * @return Copy of this random generator.
                 */
                std::unique_ptr<marathon::RandomGenerator> copy() const {
                    // todo: improve performance by avoiding the redundant construction of auxiliary tables
                    return std::make_unique<RandomGeneratorExact>(_inst);
                }
            };
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_FIXED_MARGIN_RANDOM_GENERATOR_EXACT_H
