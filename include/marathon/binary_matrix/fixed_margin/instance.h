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

#ifndef MARATHON_BINARY_MATRIX_FIXED_MARGIN_MARGIN_H
#define MARATHON_BINARY_MATRIX_FIXED_MARGIN_MARGIN_H

#include <numeric>
#include <algorithm>
#include <mutex>
#include <boost/unordered_map.hpp>

#include "../binary_matrix.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Pair of row and column sums.
             */
            struct Instance {

                std::vector<int> _rowsum;
                std::vector<int> _colsum;
                std::vector<size_t> _rowindex;
                std::vector<size_t> _colindex;

                // Dummy constructor
                Instance() = default;

                /**
                 * Define margins for a matrix of size nrow times ncol.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns.
                 */
                Instance(size_t nrow, size_t ncol) {
                    _rowsum.resize(nrow);
                    _colsum.resize(ncol);
                    _rowindex.resize(nrow);
                    _colindex.resize(ncol);
                    for (size_t i = 0; i < nrow; i++)
                        _rowindex[i] = i;
                    for (size_t j = 0; j < ncol; j++)
                        _colindex[j] = j;
                }

                /**
                 * Define margins for a matrix of size nrow times ncol.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns.
                 */
                Instance(
                        const int *rowsum,
                        const int *colsum,
                        size_t nrow,
                        size_t ncol
                ) {

                    this->_rowsum.assign(rowsum, rowsum + nrow);
                    this->_colsum.assign(colsum, colsum + ncol);
                    this->_rowindex.resize(nrow);
                    this->_colindex.resize(ncol);

                    for (size_t i = 0; i < nrow; i++)
                        this->_rowindex[i] = i;
                    for (size_t j = 0; j < ncol; j++)
                        this->_colindex[j] = j;
                }

                /**
                 * Define margins for a matrix of size nrow times ncol.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 */
                Instance(
                        std::vector<int> rowsum,
                        std::vector<int> colsum
                ) : _rowsum(std::move(rowsum)),
                    _colsum(std::move(colsum)) {

                    const size_t nrow = _rowsum.size();
                    const size_t ncol = _colsum.size();

                    _rowindex.resize(nrow);
                    _colindex.resize(ncol);

                    for (size_t i = 0; i < nrow; i++)
                        _rowindex[i] = i;
                    for (size_t j = 0; j < ncol; j++)
                        _colindex[j] = j;

                }

                /**
                * Define margins for a matrix of size nrow times ncol.
                * @param String-encoded margins. Instances have the form "2,2,2;1,2,1,2".
                * The semicolon separates the row sums from the column sums.
                */
                Instance(const std::string &inst) {

                    // convert to C string
                    const char *str = inst.c_str();

                    int k = 0;
                    int currentNumber = 0;

                    std::vector<int> *current_vec = &_rowsum;

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
                                currentNumber *= 10;
                                currentNumber += (int) (str[k] - 48);
                                break;
                            case ',':
                                current_vec->push_back(currentNumber);
                                currentNumber = 0;
                                break;
                            case ';':
                                current_vec->push_back(currentNumber);
                                currentNumber = 0;
                                current_vec = &_colsum;
                                break;
                            default:
                                throw std::runtime_error("Error while parsing instance.");
                        }
                        k++;
                    }

                    current_vec->push_back(currentNumber);

                    _rowindex.resize(_rowsum.size());
                    _colindex.resize(_colsum.size());

                    for (size_t i = 0; i < _rowsum.size(); i++)
                        _rowindex[i] = i;
                    for (size_t j = 0; j < _colsum.size(); j++)
                        _colindex[j] = j;
                }

                /**
                * Define margins for a matrix of size nrow times ncol.
                * @param Binary matrix.
                */
                Instance(const BinaryMatrix &bin) {

                    const size_t nrow = bin.getNumRows();
                    const size_t ncol = bin.getNumCols();

                    _rowsum.resize(nrow);
                    _colsum.resize(ncol);
                    _rowindex.resize(nrow);
                    _colindex.resize(ncol);

                    // determine row and column sums of currentState
                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {
                            int sij = bin.get(i, j);
                            if (sij) {
                                _rowsum[i] += sij;
                                _colsum[j] += sij;
                            }
                        }
                    }

                    for (size_t i = 0; i < _rowsum.size(); i++)
                        _rowindex[i] = i;
                    for (size_t j = 0; j < _colsum.size(); j++)
                        _colindex[j] = j;
                }


                /**
                 * Return a string representation.
                 * @return
                 */
                const std::string toString() const {
                    std::stringstream ss;
                    ss << "_rowindex              =";
                    for (int i = 0; i < _rowindex.size(); i++)
                        ss << " " << _rowindex[i];
                    ss << "\n";
                    ss << "rowsum                =";
                    for (int i = 0; i < _rowsum.size(); i++)
                        ss << " " << _rowsum[i];
                    ss << "\n";
                    ss << "_colindex              =";
                    for (int j = 0; j < _colindex.size(); j++)
                        ss << " " << _colindex[j];
                    ss << "\n";
                    ss << "colsum                =";
                    for (int j = 0; j < _colsum.size(); j++)
                        ss << " " << _colsum[j];
                    return ss.str();
                }

                /**
                 * Return the number of rows.
                 * @return Number of rows.
                 */
                const size_t getNumRows() const {
                    return _rowsum.size();
                }

                /**
                 * Return the number of columns.
                 * @return Number of columns.
                 */
                const size_t getNumCols() const {
                    return _colsum.size();
                }

                /**
                * Determine the sum of all row sum entries or column sum entries.
                * @return Row total.
                */
                const uint getTotal() const {
                    return std::accumulate(_rowsum.begin(), _rowsum.end(), 0u);
                }

                /**
                 * Verify whether the row and column sums of a binary matrix match the prescribed ones.
                 * @param s Binary Matrix.
                 * @return True, if s is valid or False, otherwise.
                 */
                virtual bool isValid(const BinaryMatrix &bin) const {

                    const int nrow = (int) getNumRows();
                    const int ncol = (int) getNumCols();

                    if (bin.getNumRows() != nrow || bin.getNumCols() != ncol)
                        return false;

                    int *rowsum_tmp = new int[nrow];
                    int *colsum_tmp = new int[ncol];
                    memset(rowsum_tmp, 0, nrow * sizeof(int));
                    memset(colsum_tmp, 0, ncol * sizeof(int));

                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {
                            if (bin.get(i, j)) {
                                rowsum_tmp[i]++;
                                colsum_tmp[j]++;
                            }
                        }
                    }

                    bool valid = true;
                    for (int i = 0; i < nrow; i++) {
                        if (rowsum_tmp[i] != this->_rowsum[i]) {
                            valid = false;
                            break;
                        }
                    }

                    if (valid) {
                        for (int j = 0; j < ncol; j++) {
                            if (colsum_tmp[j] != this->_colsum[j]) {
                                valid = false;
                                break;
                            }
                        }
                    }

                    delete[] rowsum_tmp;
                    delete[] colsum_tmp;

                    return valid;
                }
            };
        }
    }
}

/* Special Instances Hardcoded */

const marathon::binary_matrix::fixed_margin::Instance darwin_margin(darwin_matrix);

#endif //MARATHON_BINARY_MATRIX_FIXED_MARGIN_MARGIN_H
