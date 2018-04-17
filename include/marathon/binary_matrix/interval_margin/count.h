/*
 * Created on: May 04, 2017
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

#ifndef MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_COUNT_H
#define MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_COUNT_H

#include "marathon/binary_matrix/common.h"
#include "marathon/binary_matrix/interval_margin/instance.h"
#include "marathon/binary_matrix/interval_margin/realize.h"
#include "marathon/binary_matrix/fixed_margin/decompose.h"
#include "marathon/binary_matrix/fixed_margin/count.h"

#ifdef USE_BOOST_SERIALIZATION
#include <fstream>
#include <boost/serialization/boost_unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * A counter for the number of binary matrices, whose row and column
             * sums lie in intervals prescribed by integer sequences.
             */
            class Counter : public marathon::Counter {

            protected:

                /* Auxiliary Classes */
                struct A {
                    int index;
                    int lower;
                    int upper;
                };

                struct Group {
                    int lower;
                    int upper;
                    int last;
                    int size;

                    Group() {

                    }

                    Group(const int lower, const int upper, const int last)
                            : lower(lower), upper(upper), last(last), size(1) {

                    }

                    int getLast() const {
                        return last;
                    }

                    int getFirst() const {
                        return last - size + 1;
                    }
                };

                // class members
                Instance _seq;
                int _nrow, _ncol;
                boost::unordered_map<Integer, Integer> _tmp;

                /**
                 * Store the key-value-pair in the hash table.
                 * @param nrow Number of rows. Part of the key.
                 * @param colsum_lower Sequence of lower column sums. Part of the key.
                 * @param colsum_upper Sequence of upper column sums. Part of the key.
                 * @param ncol Number of columns
                 * @param value Number of matrices for the sub problem.
                 */
                void store(
                        const int nrow,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        const int ncol,
                        Integer &value
                ) {
                    const int base = _nrow + 1;
                    auto key0 = evaluate_polynomial(colsum_lower, ncol, base, nrow);
                    auto key = evaluate_polynomial(colsum_upper, ncol, base, key0);
                    _tmp[key] = value;
                }

                /**
                 * Load the value for a given key from the hash table.
                 * @param nrow Number of rows. Part of the key.
                 * @param colsum_lower Sequence of lower column sums. Part of the key.
                 * @param colsum_upper Sequence of upper column sums. Part of the key.
                 * @param ncol Number of columns
                 * @param value (Output parameter) If found, the value is stored here.
                 * @return True, if the key-value-pair is stored in table. False, otherwise.
                 */
                bool load(
                        const int nrow,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        const int ncol,
                        Integer &value
                ) const {

                    const int base = _nrow + 1;
                    auto key0 = evaluate_polynomial(colsum_lower, ncol, base, nrow);
                    auto key = evaluate_polynomial(colsum_upper, ncol, base, key0);

                    auto it = _tmp.find(key);
                    auto end = _tmp.end();

                    if (it == end)
                        return false;
                    value = it->second;
                    return true;
                }

                /**
                 * Divide the column sequences into groups of identical (lower,upper)-pairs.
                 * @param columns Sequence of lower and upper column sums.
                 * @param groups (Output parameter) List of column groups.
                 */
                void group_columns(
                        const std::vector<A> &columns,
                        std::vector<Group> &groups
                ) {
                    // skip columns with non-positive upper bound
                    int j = _ncol - 1;
                    while (j >= 0 && columns[j].upper <= 0)
                        j--;

                    if (j < 0)
                        return;

                    Group currentGroup(columns[j].lower, columns[j].upper, j);
                    j--;
                    while (j >= 0) {
                        if (columns[j].lower == currentGroup.lower &&
                            columns[j].upper == currentGroup.upper) {
                            // extend current group
                            currentGroup.size++;
                        } else {
                            // create new group
                            groups.push_back(currentGroup);
                            currentGroup = Group(columns[j].lower, columns[j].upper, j);
                        }
                        j--;
                    }
                    // push last group
                    groups.push_back(currentGroup);
                }


                /**
                 * Recursively count the number of matrices whose row and column sums
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
                 * @return The number of matrices prescribed by the intervals.
                 */
                Integer count_recursive(
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

                    /*char S[1000];
                    memset(S, 0, 1000);
                    int ii = 0;
                    for (int i = 0; i < (_rowsum_lower.size() - nrow) * ncol + group_id; i++, ii++)
                        S[ii] = ' ';
                    printf("%scount_recursive(%i,%i,%i)\n", S, nrow, group_id, aggr_in_row);
                    printf("%saggr_in_row                  = %i\n", S, aggr_in_row);
                    printf("%snrow                         = %i\n", S, nrow);
                    printf("%sgroup_id                     = %i\n", S, group_id);
                    printf("%scolumns_lower                =", S);
                    for (int i = 0; i < ncol; i++)
                        printf(" %i", _colsum_lower[i]);
                    printf("\n");
                    printf("%scolumns_upper                =", S);
                    for (int i = 0; i < ncol; i++)
                        printf(" %i", _colsum_upper[i]);
                    printf("\n");
                    printf("%scolumn_groups                = ", S);
                    for (int i = 0; i < column_groups.size(); i++)
                        printf("[%i, %i, %i, %i] ", column_groups[i].lower, column_groups[i].upper,
                               column_groups[i].size, column_groups[i].getFirst());
                    printf("\n");
                    printf("%scolsum_upper_conj            =", S);
                    for (int i = 0; i < nrow; i++)
                        printf(" %i", colsum_upper_conj[i]);
                    printf("\n");
                    printf("%srowsum_lower                 =", S);
                    for (int i = 0; i < nrow; i++)
                        printf(" %i", rowsum_lower[i]);
                    printf("\n");
                    printf("%srowsum_upper_conj            =", S);
                    for (int i = 0; i < ncol; i++)
                        printf(" %i", rowsum_upper_conj[i]);
                    printf("\n");*/

                    // if the whole matrix has completely been processed
                    if (nrow == 0) {

                        // are all lower bounds on the column sums fulfilled?
                        return colsum_lower[0] > 0 ? 0 : 1;
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

                        Integer res;

                        // is the sequence tuple realizable at all?
                        if (!isDominating(colsum_upper_conj, rowsum_lower + 1, nrow - 1) ||
                            !isDominating(rowsum_upper_conj, colsum_lower, ncol)) {
                            res = 0;
                        } else {

                            // try to find solution in hash table
                            bool found = load(nrow - 1, colsum_lower_new, colsum_upper_new, ncol, res);
                            if (!found) {

                                // re-group columns
                                std::vector<Group> new_groups;
                                group_columns(new_columns, new_groups);

                                // recursively process the nrow-1 times ncol submatrix.
                                res = count_recursive(
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
                                        0);

                                // store result
                                store(nrow - 1, colsum_lower_new, colsum_upper_new, ncol, res);
                            }
                        }

                        // undo modification
                        for (int i = 0; i < rowsum_upper[0]; i++)
                            rowsum_upper_conj[i]++;

                        delete[] colsum_lower_new;
                        delete[] colsum_upper_new;
                        delete[] colsum_index_new;
                        delete[] colsum_upper_conj_new;

                        return res;
                    }

                    // we backtrack all choices of distributing ones in the current group
                    Integer res = 0;

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

                    /*printf("%sgl                           = %i\n", S, gl);
                    printf("%sgu                           = %i\n", S, gu);
                    printf("%sn                            = %i\n", S, n);
                    printf("%sa                            = %i\n", S, a);
                    printf("%sb                            = %i\n", S, b);
                    printf("%slower                        = %i\n", S, lower);
                    printf("%supper                        = %i\n", S, upper);*/

                    // for each valid choice
                    for (int k = lower; k <= upper; k++) {

                        // apply choice: distribute k ones in the current group
                        colsum_upper_conj[gu - 1] -= k;
                        const int last = column_groups[group_id].getLast();
                        for (int l = 0, j = last; l < k; l++, j--) {
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

                        //printf("%scount_recursive(%i,%i,%i): k=%i, val=%s\n", S, nrow, group_id, aggr_in_row, k,
                        //      val.convert_to<std::string>().c_str());

                        // increase the total number of matrices by the matrices resulting from this choice
                        res += val;

                        // undo modification
                        colsum_upper_conj[gu - 1] += k;
                        for (int l = 0, j = last; l < k; l++, j--) {
                            colsum_lower[j] = gl;
                            colsum_upper[j] = gu;
                        }
                    }

                    //printf("%scount_recursive(%i,%i,%i): res = %s\n", S, nrow, group_id, aggr_in_row,
                    //       res.convert_to<std::string>().c_str());
                    return res;
                }

            public:

                /**
				 * Create a Counter for the number of binary matrices whose row and column
				 * sums lie in the prescribed intervals.
				 * @param seq Four-Tuple of integer vectors.
				 */
                Counter(Instance seq) :
                        _seq(std::move(seq)),
                        _nrow((int) _seq._rowsum_upper.size()),
                        _ncol((int) _seq._colsum_lower.size()) {

                }

                /**
                 * Create a Counter for the number of binary matrices whose row and column
                 * sums lie in the prescribed intervals.
                 * @param rowsum_lower Lower bounds on the row sums.
                 * @param rowsum_upper Upper bounds on the row sums.
                 * @param colsum_lower Lower bounds on the column sums.
                 * @param colsum_upper Upper bounds on the column sums.
                 */
                Counter(
                        const std::vector<int> &rowsum_lower,
                        const std::vector<int> &rowsum_upper,
                        const std::vector<int> &colsum_lower,
                        const std::vector<int> &colsum_upper
                ) : Counter(Instance(&rowsum_lower[0], &rowsum_upper[0],
                                     &colsum_lower[0], &colsum_upper[0],
                                     (int) rowsum_lower.size(), (int) colsum_lower.size())) {

                }

#ifdef USE_BOOST_SERIALIZATION
                /**
                 * Create a counter for the number of binary matrices from a
                 * previously dumped instance.
                 * @param ifs Input file stream.
                 */
                Counter(std::ifstream &ifs) {

                    // create and open an archive for input
                    boost::archive::text_iarchive ia(ifs);

                    ia >> _seq.rowsum_lower;
                    ia >> _seq.rowsum_upper;
                    ia >> _seq._colsum_lower;
                    ia >> _seq._colsum_upper;

                    _nrow = (int) _seq.rowsum_lower.size();
                    _ncol = (int) _seq._colsum_lower.size();

                    ia >> _tmp;
                }
#endif

                /**
                 * Count the number of binary matrices whose row and column sums lie in the
                 * prescribed intervals.
                 * @return Number of matrices.
                 */
                Integer count() override {

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

                    // try to find pre-computed result in hash table
                    Integer res;
                    bool found = load(_nrow, colsum_lower, colsum_upper, _ncol, res);
                    if (!found) {

                        // group columns by identical pairs of lower and upper bounds
                        std::vector<Group> groups;
                        group_columns(columns, groups);

                        // start recursive counting
                        res = count_recursive(
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

                        // store in table
                        store(_nrow, colsum_lower, colsum_upper, _ncol, res);
                    }

                    delete[] rowsum_lower;
                    delete[] rowsum_upper;
                    delete[] rowsum_index;
                    delete[] rowsum_upper_conj;
                    delete[] colsum_lower;
                    delete[] colsum_upper;
                    delete[] colsum_index;
                    delete[] colsum_upper_conj;

                    return res;
                }

#ifdef USE_BOOST_SERIALIZATION
                /**
                 * Dump the content of the hash table to disk.
                 * @param ofs Output stream.
                 */
                void dump(std::ofstream ofs) {

                    // create archive
                    boost::archive::text_oarchive oa(ofs);

                    // dump attributes
                    oa << _seq.rowsum_lower;
                    oa << _seq.rowsum_upper;
                    oa << _seq._colsum_lower;
                    oa << _seq._colsum_upper;

                    // dump hash table
                    oa << _tmp;
                }
#endif
            };
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_COUNT_H
