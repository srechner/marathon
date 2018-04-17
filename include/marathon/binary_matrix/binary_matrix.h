/*
 * Created on: Nov 24, 2014
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

#ifndef MARATHON_BINARY_MATRIX_H
#define MARATHON_BINARY_MATRIX_H

#include <cstdlib>
#include <iostream>

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS  // so dynamic_bitset becomes hashable

#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

#include "marathon/state.h"

namespace marathon {
    namespace binary_matrix {

        /**
         * This class represents a 0-1-matrix of size nrow times ncol.
         * It can be interpreted as biadjacency matrix of a bipartite graph.
         */
        class BinaryMatrix : public State {

        private:

            size_t _nrow, _ncol;            // number of rows and columns
            boost::dynamic_bitset<> _bits;  // the matrix is stored as a single row in a bitset.

            /**
             * Return the internal index p of the element (i,j).
             * @param i Row index.
             * @param j Column index.
             * @return Internal position of element (i,j).
             */
            size_t coord_transform(size_t i, size_t j) const {
                const size_t p = i * _ncol + j;
                assert(p < _bits.size());
                /*if (p >= _bits.size()) {
                    std::stringstream ss;
                    ss << "marathon::Exception: Invalid Matrix Element: ( " << i << "," << j << ")";
                    throw std::runtime_error(ss.str());
                }*/
                return p;
            }

        public:

            /**
             * Create a binary matrix of size 0 times 0.
             */
            BinaryMatrix() : _nrow(0), _ncol(0) {

            }

            /**
             * Create a new binary matrix with nrow rows and ncol columns.
             * If given, initialize the binary matrix with the values stored in bits.
             * @param nrow Number of rows.
             * @param ncol Number of columns.
             */
            BinaryMatrix(
                    size_t nrow,
                    size_t ncol
            ) : _nrow(nrow), _ncol(ncol),
                _bits(boost::dynamic_bitset<>(nrow * ncol)) {

            }

            /**
             * Create a new binary matrix with nrow rows and ncol columns.
             * If given, initialize the binary matrix with the values stored in bits.
             * @param nrow Number of rows.
             * @param ncol Number of columns.
             * @param bits Vector of booleans that is interpreted as a nrow times ncol that has been
             * flattened to a single row.
             */
            BinaryMatrix(size_t nrow, size_t ncol, const std::vector<bool> &bits)
                    : BinaryMatrix(nrow, ncol) {

                for (size_t i = 0; i < nrow * ncol; i++) {
                    _bits.set(i, bits[i]);
                }
            }

            /**
             * Create a binary matrix from a binary string.
             * @param nrow Number of rows.
             * @param ncol Number of columns.
             * @param str Binary string iterpreted as a nrow times ncol matrix flattened to a single row.
             */
            BinaryMatrix(size_t nrow, size_t ncol, const std::string &str) :
                    BinaryMatrix(nrow, ncol) {

                for (size_t i = 0; i < nrow * ncol; i++) {
                    if (str[i] == '0')
                        _bits.set(i, 0);
                    else if (str[i] == '1')
                        _bits.set(i, 1);
                    else {
                        std::stringstream ss;
                        ss << "Error! Unexpected symbol: " << str[i];
                        throw std::runtime_error(ss.str());
                    }

                }
            }

            /**
             * Return a reference to the internal matrix representation.
             * @return Internal matrix representation.
             */
            const boost::dynamic_bitset<> &getBitset() const {
                return _bits;
            }

            /**
             * Return the number of nonzero entries.
             * @return Number of nonzero entries.
             */
            size_t getTotal() const {
                return _bits.count();
            }

            /**
             * Return the number of rows.
             * @return Return the number of rows.
             */
            size_t getNumRows() const {
                return _nrow;
            }

            /**
             * @return Return the number of columns.
             */
            size_t getNumCols() const {
                return _ncol;
            }

            /**
             * Return the value of matrix entry (i,j).
             * @param i Row index.
             * @param j Column index.
             * @return Value of entry (i,j).
             */
            bool get(size_t i, size_t j) const {
                const size_t p = coord_transform(i, j);
                return _bits[p];
            }

            /**
             * Set matrix entry (i,j) to b.
             * @param i Row index.
             * @param j Column index.
             * @param b Zero or One.
             */
            void set(size_t i, size_t j, bool b) {
                const size_t p = coord_transform(i, j);
                _bits[p] = b;
            }

            /**
             * Set all elements to zero.
             */
            void clear() {
                _bits.reset();
            }

            /**
             * Flip the matrix entry (i,j).
             * @param i Row index.
             * @param j Column index.
             */
            void flip(size_t i, size_t j) {
                const size_t p = coord_transform(i, j);
                _bits[p].flip();
            }

            /**
             * Return a fingerprint of the binary matrix.
             * @return Hash value of this matrix.
             */
            size_t hashValue() const override {
                return boost::hash_value(_bits.m_bits);
            }

            /**
             * Compare operator.
             * @param x State pointer that will be interpreted as a binary matrix pointer.
             * @return 0, if x is identical to this matrix, 1 or -1, otherwise.
             */
            int compare(const State &x) const override {

                // try to convert state to binary matrix
                auto d = dynamic_cast<const BinaryMatrix *>(&x);

                if(d == nullptr)
                    return -1;

                if (operator==(*d))
                    return 0;
                else {
                    if (operator<(*d))
                        return -1;
                }
                return 1;
            }

            /**
             * Return a compact string representation of the binary matrix.
             * @return A string of zeros and ones. The rows of the matrix are concatenated to a single line.
             */
            std::string toString() const override {
                std::stringstream ss;
                for (unsigned int i = 0; i < _bits.size(); i++)
                    ss << _bits[i];
                return ss.str();
            }

            /**
             * Create nicely formatted matrix representation.
             * @param name Name of the matrix.
             * @return
             */
            std::string fancyString() const {
                std::stringstream ss;
                for (int i = 0; i < _nrow; i++) {
                    for (int j = 0; j < _ncol; j++) {
                        ss << "  " << get(i, j);
                    }
                    ss << std::endl;
                }
                return ss.str();
            }

            /**
             * Return a deep copy of the matrix.
             * @return Copy of this matrix.
             */
            std::unique_ptr<State> copy() const override {
                return std::make_unique<BinaryMatrix>(*this);
            }

            /**
             * Tests if (i,j), (i,l), (l,j), and (k,l) is a submatrix of form
             *    (0 1)   or   (1 0)
             *    (1 0)        (0 1).
             * @param i Row index.
             * @param j Column index.
             * @param k Row index.
             * @param l Column index.
             * @return
             */
            bool isCheckerBoardUnit(uint i, uint j, uint k, uint l) const {

                const bool a = get(i, j);
                const bool b = get(i, l);
                const bool c = get(k, j);
                const bool d = get(k, l);

                return (a == d) && (b == c) && (a == !b);
            }

            /**
             * Flip the entries (i,j), (i,l), (k,j), and (k,l).
             * @param i Row index.
             * @param j Column index.
             * @param k Row index.
             * @param l Column index.
             */
            void flipSubmatrix(uint i, uint j, uint k, uint l) {
                flip(i, j);
                flip(i, l);
                flip(k, j);
                flip(k, l);
            }

            /**
             * Comparison operator with respect to lexicographical order.
             * @param s Binary matrix.
             * @return
             */
            bool operator<(const BinaryMatrix &s) const {
                return getNumRows() == s.getNumRows()
                       && getNumCols() == s.getNumCols()
                       && _bits < s._bits;
            }

            /**
             * Equal operator. Two matrices are equal if they have the same size
             * and if their elements are identical.
             * @param s Binary matrix.
             * @return
             */
            bool operator==(const BinaryMatrix &s) const {
                return getNumRows() == s.getNumRows()
                       && getNumCols() == s.getNumCols()
                       && _bits == s._bits;
            }

        };

    }
}

// overload standard hash function of BinaryMatrix objects
namespace std {
    template<>
    struct hash<marathon::binary_matrix::BinaryMatrix> {
        size_t operator()(const marathon::binary_matrix::BinaryMatrix &x) const {
            return x.hashValue();
        }
    };
}


// Special hard coded instances
const marathon::binary_matrix::BinaryMatrix darwin_matrix(
        13, 17, "00111111110111111"
                "11111111110101100"
                "11111111111101100"
                "00111001010110111"
                "11101111110101100"
                "00000000001010000"
                "00111111100101100"
                "00000000000100000"
                "00111111110100100"
                "00111111110101100"
                "00111011010000000"
                "00110000000000000"
                "11111111111111111");

#endif /* MARATHON_BINARY_MATRIX_H */
