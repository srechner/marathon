/*
 * Created on: May 19, 2017
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

#ifndef MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_RANDOM_GENERATOR_EXACT_H
#define MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_RANDOM_GENERATOR_EXACT_H

// marathon includes
#include "marathon/random_device.h"
#include "marathon/binary_matrix/random_generator.h"
#include "marathon/binary_matrix/interval_margin/count.h"
#include "marathon/binary_matrix/fixed_margin/random_generator_exact.h"

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * A random generator for binary matrices whose row sums and column sums
             * lie in prescribed intervals. This random generator produces (exactly!)
             * uniform distributed samples.
             */
            class RandomGeneratorExact :
                    public marathon::binary_matrix::RandomGenerator,
                    private Counter {

                Integer _num_matrices;   // total number of matrices
                Integer _target;         // we want to select the matrix with number target
                Integer _mul;            // to avoid division
                RandomDevice _rg;              // random generator
                BinaryMatrix *_bin;      // binary matrix

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
				 */
                void sample_recursive(
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
                        const int aggr_in_row
                ) {

                    // if the whole matrix has completely been processed
                    if (nrow == 0) {
                        return; // we are done
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
                        for (int j = 0; j < ncol; j++) {
                            colsum_lower[j] = new_columns[j].lower;
                            colsum_upper[j] = new_columns[j].upper;
                            colsum_index[j] = new_columns[j].index;
                        }
                        conjugate(colsum_upper_conj, colsum_upper, nrow, ncol);

                        // modify rowsum_upper_conj
                        for (int i = 0; i < rowsum_upper[0]; i++)
                            rowsum_upper_conj[i]--;

                        // is the sequence tuple realizable at all?
                        assert(isDominating(colsum_upper_conj, rowsum_lower + 1, nrow - 1) &&
                               isDominating(rowsum_upper_conj, colsum_lower, ncol));

                        // re-group columns
                        std::vector<Group> new_groups;
                        group_columns(new_columns, new_groups);

                        // recursively process the nrow-1 times ncol submatrix.
                        sample_recursive(
                                rowsum_lower + 1,
                                rowsum_upper + 1,
                                rowsum_index + 1,
                                rowsum_upper_conj,
                                nrow - 1,
                                colsum_lower,
                                colsum_upper,
                                colsum_index,
                                colsum_upper_conj,
                                ncol,
                                new_groups,
                                0,
                                0);

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

                    // Each value of k in the range [lower,upper] results in a number of matrices val(k).
                    // We seek the smallest k such that (val(k_lower) + ... + val(k)) > target.
                    int k = lower;

                    Integer skipped = 0;

                    // for each valid choice
                    for (; k <= upper; k++) {

                        // apply choice: distribute k ones in the current group
                        colsum_upper_conj[gu - 1] -= k;
                        const int last = column_groups[group_id].getLast();
                        for (int j = last; j > last - k; j--) {
                            colsum_lower[j] = std::max(0, gl - 1);
                            colsum_upper[j] = gu - 1;
                        }

                        // recursively calculate the number of matrices that result from this choice
                        Integer val = count_recursive(
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
                                aggr_in_row + k
                        );

                        // multiply by the number of ways to distribute k ones to n entries
                        val *= marathon::binom(n, k);

                        //printf("%ssample_recursive(%i,%i,%i): k=%i, val=%s\n", S, nrow, group_id, aggr_in_row, k,
                        //       val.convert_to<std::string>().c_str());

                        // smallest k found!
                        if ((skipped + val) * _mul > _target) {
                            break;
                        }

                        // increase the total number of matrices by the matrices resulting from this choice
                        skipped += val;

                        // undo modification
                        colsum_upper_conj[gu - 1] += k;
                        for (int j = last; j > last - k; j--) {
                            colsum_lower[j] = gl;
                            colsum_upper[j] = gu;
                        }
                    }

                    assert(k <= upper);

                    // reduce t by the number of matrices that result from choices smaller than k
                    _target -= skipped * _mul;

                    // std::cout << "reduce " << k << " colums with colsum " << pos << " by one" << std::endl;

                    /*
                     * k is the smallest value such that (val(l) + ... + val(k)) > target.
                     * Distribute k ones to the n columns of the current row that currently
                     * have value pos.
                     * First element with value pos is colsum_first[pos].
                     * Last element with value pos is colsum_last[pos].
                     */

                    // randomly select k columns out of n
                    const int first = column_groups[group_id].getFirst();
                    const int last = column_groups[group_id].getLast();
                    _rg.shuffle<int>(colsum_index + first, n);

                    // for each selected column
                    for (int j = last; j > last - k; j--)
                        _bin->set(rowsum_index[0], colsum_index[j], 1);

                    /*
					 * The number of matrices in the next level of recursion is
					 * scaled up by binom(n,k), thus the matrix with number t in
					 * the current level of recursion corresponds to the matrix
					 * with number t / binom(n,k).
					 */
                    _mul *= marathon::binom(n, k);

                    // recursively fill the remaining matrix
                    sample_recursive(
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
                            aggr_in_row + k
                    );
                }

            public:

                /**
				 * Create a random generator for binary matrices whose row and column
				 * sums lie in the prescribed intervals.
				 * @param rowsum_lower Lower bounds on the row sums.
				 * @param rowsum_upper Upper bounds on the row sums.
				 * @param colsum_lower Lower bounds on the column sums.
				 * @param colsum_upper Upper bounds on the column sums.
				 */
                RandomGeneratorExact(const Instance &inst)
                        : Counter(inst) {
                    _bin = new BinaryMatrix(_nrow, _ncol);
                    _num_matrices = count();
                    assert(_num_matrices > 0);
                }

                /**
				 * Create a random generator for binary matrices whose row and column
				 * sums lie in the prescribed intervals.
				 * @param rowsum_lower Lower bounds on the row sums.
				 * @param rowsum_upper Upper bounds on the row sums.
				 * @param colsum_lower Lower bounds on the column sums.
				 * @param colsum_upper Upper bounds on the column sums.
				 */
                RandomGeneratorExact(
                        const std::vector<int> &rowsum_lower,
                        const std::vector<int> &rowsum_upper,
                        const std::vector<int> &colsum_lower,
                        const std::vector<int> &colsum_upper
                ) : RandomGeneratorExact(
                        Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper)) {

                }

                /**
				 * Create a random generator for binary matrices whose row and column
				 * sums lie in the prescribed intervals.
				 * @param rowsum_lower Lower bounds on the row sums.
				 * @param rowsum_upper Upper bounds on the row sums.
				 * @param colsum_lower Lower bounds on the column sums.
				 * @param colsum_upper Upper bounds on the column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns.
				 */
                RandomGeneratorExact(
                        const int *rowsum_lower,
                        const int *rowsum_upper,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        const int nrow,
                        const int ncol
                ) : RandomGeneratorExact(
                        std::vector<int>(rowsum_lower, rowsum_lower + nrow),
                        std::vector<int>(rowsum_upper, rowsum_upper + nrow),
                        std::vector<int>(colsum_lower, colsum_lower + ncol),
                        std::vector<int>(colsum_upper, colsum_upper + ncol)
                ) {

                }

#ifdef USE_BOOST_SERIALIZATION
                /**
                 * Create a random generator based on the information previously dumped to file.
                 * @param ifs Input file stream.
                 */
                explicit RandomGeneratorExact(std::ifstream &ifs)
                        : Counter(ifs) {
                    _bin = new BinaryMatrix(_nrow, _ncol);
                    _num_matrices = count();
                    assert(_num_matrices > Integer(0));
                }
#endif


                virtual ~RandomGeneratorExact() {
                    delete _bin;
                }

                /**
				 * Return a binary matrix whose row and column sums match the
                 * given integer vectors. The sample is assured to be distributed
                 * (exactly!) uniform.
				 * @return Random binary matrix.
				 */
                const BinaryMatrix &next() override {

                    //std::cout << "=========================" << std::endl;

                    if (_num_matrices == 0)
                        throw std::runtime_error("Error! Four-tuple of integer vectors is not realizable!");

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

                    // select a random number in range [0,num_matrices)
                    _target = _rg.nextInt(_num_matrices);

                    //std::cout << "num_matrices=" << _num_matrices << " target=" << _target << std::endl;

                    // clear variables
                    _bin->clear();
                    _mul = 1;

                    // select the binary matrix with number target by traversing the enumeration tree
                    sample_recursive(
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
                            0);

                    delete[] rowsum_lower;
                    delete[] rowsum_upper;
                    delete[] rowsum_index;
                    delete[] rowsum_upper_conj;
                    delete[] colsum_lower;
                    delete[] colsum_upper;
                    delete[] colsum_index;
                    delete[] colsum_upper_conj;

                    return *_bin;

                }

                /**
                 * Create an independent copy of the random generator.
                 * @return Copy of this random generator.
                 */
                std::unique_ptr<marathon::RandomGenerator> copy() const {
                    // todo: improve performance by avoiding the redundant construction of auxiliary tables
                    return std::make_unique<RandomGeneratorExact>(_seq);
                }
            };
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_RANDOM_GENERATOR_EXACT_H
