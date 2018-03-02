/*
 * Created on: Jun 06, 2017
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

#ifndef MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_MARGIN_H
#define MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_MARGIN_H

#include <numeric>
#include <algorithm>
#include <mutex>
#include <boost/unordered_map.hpp>

#include "../binary_matrix.h"

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * Lower and upper bounds on matrix margins.
             */
            struct Instance {

                std::vector<int> rowsum_lower;
                std::vector<int> rowsum_upper;
                std::vector<int> colsum_lower;
                std::vector<int> colsum_upper;
                std::vector<int> rowindex;
                std::vector<int> colindex;

                /**
                 * Dummy constructor.
                 */
                Instance() {

                }

                /**
                 * Define margins for a matrix of size nrow times ncol.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns
                 * @return
                 */
                Instance(const int nrow, const int ncol) {
                    rowsum_lower.resize(nrow);
                    rowsum_upper.resize(nrow);
                    colsum_lower.resize(ncol);
                    colsum_upper.resize(ncol);
                    rowindex.resize(nrow);
                    colindex.resize(ncol);
                    for (int i = 0; i < nrow; i++)
                        rowindex[i] = i;
                    for (int j = 0; j < ncol; j++)
                        colindex[j] = j;
                }

                /**
                 * Define margins for a matrix of size nrow times ncol.
                 * @param rowsum_lower Lower bounds on row sums.
                 * @param rowsum_upper Upper bounds on row sums.
                 * @param colsum_lower Lower bounds on column sums.
                 * @param colsum_upper Upper bounds on column sums.
                 * @return
                 */
                Instance(
                        const std::vector<int>& rowsum_lower,
                        const std::vector<int>& rowsum_upper,
                        const std::vector<int>& colsum_lower,
                        const std::vector<int>& colsum_upper
                ) : rowsum_lower(rowsum_lower),
                    rowsum_upper(rowsum_upper),
                    colsum_lower(colsum_lower),
                    colsum_upper(colsum_upper) {

                    const int nrow = rowsum_lower.size();
                    const int ncol = colsum_lower.size();

                    rowindex.resize(nrow);
                    colindex.resize(ncol);

                    for(int i=0; i<nrow; i++)
                        rowindex[i] = i;
                    for(int j=0; j<ncol; j++)
                        colindex[j] = j;
                }

                /**
                 * Define margins for a matrix of size nrow times ncol.
                 * @param rowsum_lower Lower bounds on row sums.
                 * @param rowsum_upper Upper bounds on row sums.
                 * @param colsum_lower Lower bounds on column sums.
                 * @param colsum_upper Upper bounds on column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns
                 * @return
                 */
                Instance(
                        const int *rowsum_lower,
                        const int *rowsum_upper,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        const int nrow,
                        const int ncol
                ) : rowsum_lower(rowsum_lower, rowsum_lower + nrow),
                    rowsum_upper(rowsum_upper, rowsum_upper + nrow),
                    colsum_lower(colsum_lower, colsum_lower + ncol),
                    colsum_upper(colsum_upper, colsum_upper + ncol) {

                    rowindex.resize(nrow);
                    colindex.resize(ncol);

                    for (int i = 0; i < nrow; i++)
                        rowindex[i] = i;
                    for (int j = 0; j < ncol; j++)
                        colindex[j] = j;
                }

                /**
                 * Define margins for a matrix of size nrow times ncol.
                 * @param m Margins that will be copied.
                 * @return
                 */
                Instance(const Instance &m)
                        : rowsum_lower(m.rowsum_lower),
                          rowsum_upper(m.rowsum_upper),
                          colsum_lower(m.colsum_lower),
                          colsum_upper(m.colsum_upper),
                          rowindex(m.rowindex),
                          colindex(m.colindex) {
                }


                /**
                 * Define margins for a matrix of size nrow times ncol.
                 * @param inst String encoded margins. Instances have the form "(l-u)*;(l-u)*",
                 * where each l is a lower, and each u is an upper bound on a certain row or column sum.
                 * For convenience, if li = ui, the string 'li-ui' can be replaced by 'li'.
                 * The semicolon separates the row sums from the column sums.
                 *
                 * For example, the instance encoding "1-2,2,2-3;0-2,0-1,1-1,1-3" corresponds to
                 *
                 *    lower bounds on row sums:    (1,2,2)
                 *    upper bounds on row sums:    (2,2,3)
                 *    lower bounds on column sums: (0,0,1,1)
                 *    upper bounds on column sums: (2,1,1,3)
                 *
                 * @return Instance representation.
                 */
                Instance(const std::string &inst) {

                    // convert to C string
                    const char *str = inst.c_str();

                    int k = 0;
                    int currentLower = 0;
                    int currentUpper = -1;
                    int *currentNumber = &currentLower;

                    std::vector<int> *current_vec_lower = &rowsum_lower;
                    std::vector<int> *current_vec_upper = &rowsum_upper;

                    // parse string
                    while (str[k] != '\0') {

                        switch (str[k]) {
                            case '0':
                            case '1':
                            case '2':
                            case '3':
                            case '4':
                            case '5':
                            case '6':
                            case '7':
                            case '8':
                            case '9':
                                *currentNumber *= 10;
                                *currentNumber += (int) (str[k] - 48);
                                break;
                            case '-':
                                currentUpper = 0;
                                currentNumber = &currentUpper;
                                break;
                            case ',':
                                current_vec_lower->push_back(currentLower);
                                current_vec_upper->push_back(currentUpper == -1 ? currentLower : currentUpper);
                                currentLower = 0;
                                currentUpper = -1;
                                currentNumber = &currentLower;
                                break;
                            case ';':
                                current_vec_lower->push_back(currentLower);
                                current_vec_upper->push_back(currentUpper == -1 ? currentLower : currentUpper);
                                currentLower = 0;
                                currentUpper = -1;
                                currentNumber = &currentLower;
                                current_vec_lower = &colsum_lower;
                                current_vec_upper = &colsum_upper;
                                break;
                            default:
                                throw std::runtime_error("Expection: malformed instance encoding!");
                        }
                        k++;
                    }

                    current_vec_lower->push_back(currentLower);
                    current_vec_upper->push_back(currentUpper == -1 ? currentLower : currentUpper);

                    const int nrow = rowsum_lower.size();
                    const int ncol = colsum_lower.size();

                    rowindex.resize(nrow);
                    colindex.resize(ncol);

                    for (int i = 0; i < nrow; i++)
                        rowindex[i] = i;
                    for (int j = 0; j < ncol; j++)
                        colindex[j] = j;

                }

                /**
                 * Return the number of rows.
                 * @return Number of rows
                 */
                const size_t getNumRows() const {
                    return rowsum_lower.size();
                }


                /**
                 * Return the number of columns.
                 * @return Number of columns
                 */
                const size_t getNumCols() const {
                    return colsum_lower.size();
                }


                /**
                 * Return a string.
                 * @return
                 */
                const std::string toString() const {

                    const int nrow = rowsum_lower.size();
                    const int ncol = colsum_lower.size();

                    std::stringstream ss;
                    ss << "rowindex              =";
                    for (int i = 0; i < nrow; i++)
                        ss << " " << rowindex[i];
                    ss << "\n";
                    ss << "rowsum_lower          =";
                    for (int i = 0; i < nrow; i++)
                        ss << " " << rowsum_lower[i];
                    ss << "\n";
                    ss << "rowsum_upper          =";
                    for (int i = 0; i < nrow; i++)
                        ss << " " << rowsum_upper[i];
                    ss << "\n";
                    ss << "colindex              =";
                    for (int j = 0; j < ncol; j++)
                        ss << " " << colindex[j];
                    ss << "\n";
                    ss << "colsum_lower          =";
                    for (int j = 0; j < ncol; j++)
                        ss << " " << colsum_lower[j];
                    ss << "\n";
                    ss << "colsum_upper          =";
                    for (int j = 0; j < ncol; j++)
                        ss << " " << colsum_lower[j];
                    return ss.str();
                }


                /**
                 * Verify whether is valid with respect to the lower and upper row and column sums.
                 * @param s State object.
                 * @return True, if s is valid or False, otherwise.
                 */
                bool isValid(const BinaryMatrix& bin) const {

                    const int nrow = (int) rowsum_lower.size();
                    const int ncol = (int) colsum_lower.size();

                    if (bin.getNumRows() != nrow || bin.getNumCols() != ncol)
                        return false;

                    int *rowsum = new int[nrow];
                    int *colsum = new int[ncol];
                    memset(rowsum, 0, nrow * sizeof(int));
                    memset(colsum, 0, ncol * sizeof(int));

                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {
                            if (bin.get(i, j)) {
                                rowsum[i]++;
                                colsum[j]++;
                            }
                        }
                    }

                    bool valid = true;
                    for (int i = 0; i < nrow; i++) {
                        if (rowsum[i] < rowsum_lower[i] || rowsum[i] > rowsum_upper[i]) {
                            valid = false;
                            break;
                        }
                    }

                    if (valid) {
                        for (int j = 0; j < ncol; j++) {
                            if (colsum[j] < colsum_lower[j] || colsum[j] > colsum_upper[j]) {
                                valid = false;
                                break;
                            }
                        }
                    }

                    delete[] rowsum;
                    delete[] colsum;

                    return valid;
                }

            };
        }
    }
}


/* Special Instances Hardcoded */

const marathon::binary_matrix::interval_margin::Instance darwin_margin_p1(
        {14,13,14,10,12,2,10,1,10,11,6,2,17},
        {15,14,15,11,13,3,11,2,11,12,7,3,18},
        {4,4,11,10,10,8, 9,10,8, 9,3,10,4,7, 9,3,3},
        {5,5,12,11,11,9,10,11,9,10,4,11,5,8,10,4,4}
);

const marathon::binary_matrix::interval_margin::Instance darwin_margin_pm1(
        {13,12,13, 9,11,1, 9,0, 9,10,5,1,16},
        {15,14,15,11,13,3,11,2,11,12,7,3,18},
        {3,3,10, 9, 9,7, 8, 9,7, 8,2, 9,3,6, 8,2,2},
        {5,5,12,11,11,9,10,11,9,10,4,11,5,8,10,4,4}
);

const marathon::binary_matrix::interval_margin::Instance darwin_margin_pm2(
        {12,11,12, 8,10,0, 8,0, 8, 9,4,0,15},
        {16,15,16,12,14,4,12,3,12,13,8,4,19},
        {2,2, 9, 8, 8, 6, 7, 8, 6, 7,1, 8,2,5, 7,1,1},
        {6,6,13,12,12,10,11,12,10,11,5,12,6,9,11,5,5}
);

#endif //MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_MARGIN_H
