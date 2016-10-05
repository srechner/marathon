/*
 * BinaryMatrix.h
 *
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

#ifndef STATE_H_
#define STATE_H_

#include <cstdlib>
#include <iostream>

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS  // so dynamic_bitset becomes hashable

#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

#include "marathon/State.h"

namespace marathon {
	namespace bipgraph {

		/**
		 * This class represents a 0-1-matrix of size nrow times ncol.
		 * It can be interpreted as biadjacency matrix of a bipartite graph.
		 */
		class BinaryMatrix : public State {

		public:

			int nrow, ncol;        // number of rows and columns
			boost::dynamic_bitset<> M;// the matrix is stored as a single row in a bitset.

			BinaryMatrix() :
					nrow(0), ncol(0) {

			}

			BinaryMatrix(int nrows, int ncols, const bool *bits = nullptr) :
					nrow(nrows), ncol(ncols) {
				if (bits == nullptr) {
					M = boost::dynamic_bitset<>(nrows * ncols);
				} else {
					for (int i = 0; i < nrows * ncols; i++) {
						M.push_back(bits[i]);
					}
				}
			}

			BinaryMatrix(int nrows, int ncols, const std::string &str) :
					nrow(nrows), ncol(ncols) {
				M = boost::dynamic_bitset<>(str);
			}

			BinaryMatrix(const BinaryMatrix &s) :
					nrow(s.nrow), ncol(s.ncol) {
				M = boost::dynamic_bitset<>(s.M);
			}

			virtual ~BinaryMatrix() {

			}

			int getNumRows() const {
				return nrow;
			}

			int getNumColumns() const {
				return ncol;
			}

			bool get(int u, int v) const {
				int p = u * ncol + v;
				if (p >= M.m_num_bits) {
					std::cerr << "has_edge::Error! " << u << " " << v << " " << M.m_num_bits
					          << std::endl;
				}
				return M[p];
			}

			void set(int u, int v, bool b) {
				int p = u * ncol + v;
				if (p >= M.m_num_bits) {
					std::cerr << "set_edge::Error! " << u << " " << v << " " << M.m_num_bits
					          << std::endl;
				}
				M[p] = b;
			}

			void flip(int u, int v) {
				int p = u * ncol + v;
				if (p >= M.m_num_bits) {
					std::cerr << "flip_edge::Error! " << u << " " << v << " "
					          << M.m_num_bits << std::endl;
				}
				M[p].flip();
			}

			size_t hashValue() const {
				return boost::hash_value(M.m_bits);
			}

			int compare(const State *x) const {
				const BinaryMatrix *d = (const BinaryMatrix *) x;
				if (operator==(*d))
					return 0;
				else {
					if (operator<(*d))
						return -1;
				}
				return 1;
			}

			std::string toString() const {
				std::stringstream ss;
				for (unsigned int i = 0; i < M.size(); i++)
					ss << M[i];
				return ss.str();
			}

			State *copy() const {
				return new BinaryMatrix(*this);
			}

			bool is_switchable(int u1, int v1, int u2, int v2) const {

				/* translate node labels */
				if (u1 > v1) {
					// rotate
					uint tmp = u1;
					u1 = v1;
					v1 = u2;
					u2 = v2;
					v2 = tmp;
				}

				v1 -= nrow;
				v2 -= nrow;

#ifdef DEBUG
				std::cout << "is switchable " << u1 << " " << v1 << " " << u2 << " " << v2
<< std::endl;

std::cout << "(" << u1 << "," << v1 << "): " << get(u1, v1)
<< std::endl;
std::cout << "(" << u1 << "," << v2 << "): " << get(u1, v2)
<< std::endl;
std::cout << "(" << u2 << "," << v1 << "): " << get(u2, v1)
<< std::endl;
std::cout << "(" << u2 << "," << v2 << "): " << get(u2, v2)
<< std::endl;
#endif

				bool a = get(u1, v1);
				bool b = get(u1, v2);
				bool c = get(u2, v1);
				bool d = get(u2, v2);

				return (a == d) && (b == c) && (a == !b);
			}

			void switch_4_cycle(int u1, int v1, int u2, int v2) {

				if (u1 > v1) {
					// rotate
					uint tmp = u1;
					u1 = v1;
					v1 = u2;
					u2 = v2;
					v2 = tmp;
				}

				// translate node labels
				v1 -= nrow;
				v2 -= nrow;

#ifdef DEBUG
				std::cout << "switch cycle [ " << u1 << " " << v1 << " " << u2 << " " << v2
<< " ]" << std::endl;
#endif

				flip(u1, v1);
				flip(u1, v2);
				flip(u2, v1);
				flip(u2, v2);
			}

			void get_row(int u, boost::dynamic_bitset<> &row) const {
				row.resize(ncol);
				for (int v = 0; v < ncol; v++) {
					row[v] = get(u, v);
				}
			}

			void operator=(const BinaryMatrix &s) {
				nrow = s.nrow;
				ncol = s.ncol;
				M = boost::dynamic_bitset<>(s.M);
			}

			bool operator<(const BinaryMatrix &s) const {
				return nrow == s.nrow && ncol == s.ncol && M < s.M;
			}

			bool operator==(const BinaryMatrix &s) const {
				return nrow == s.nrow && ncol == s.ncol && M == s.M;
			}

		};

	}
}

#endif /* STATE_H_ */
