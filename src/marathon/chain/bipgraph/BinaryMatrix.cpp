/*
 * BinaryMatrix.cpp
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

// so dynamic_bitset becomes hashable
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

// project includes
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include "../../../../include/marathon/chain/bipgraph/BinaryMatrix.h"

//#define DEBUG

namespace marathon {
namespace chain {
namespace bipgraph {

BinaryMatrix::BinaryMatrix() :
		nrow(0), ncol(0) {

}

BinaryMatrix::BinaryMatrix(int nrows, int ncols, const bool* bits) :
		nrow(nrows), ncol(ncols) {
	if (bits == nullptr) {
		M = boost::dynamic_bitset<>(nrows * ncols);
	} else {
		for (int i = 0; i < nrows * ncols; i++) {
			M.push_back(bits[i]);
		}
	}
}

BinaryMatrix::BinaryMatrix(int nrows, int ncols, const std::string& str) :
		nrow(nrows), ncol(ncols) {
	M = boost::dynamic_bitset<>(str);
}

BinaryMatrix::BinaryMatrix(const BinaryMatrix& s) :
		nrow(s.nrow), ncol(s.ncol) {
	M = boost::dynamic_bitset<>(s.M);
}

BinaryMatrix::~BinaryMatrix() {

}

int BinaryMatrix::getNumRows() const {
	return nrow;
}

int BinaryMatrix::getNumColumns() const {
	return ncol;
}

void BinaryMatrix::set(int u, int v, bool b) {
	int p = u * ncol + v;
	if (p >= M.m_num_bits) {
		std::cerr << "set_edge::Error! " << u << " " << v << " " << M.m_num_bits
				<< std::endl;
	}
	M[p] = b;
}

bool BinaryMatrix::get(int u, int v) const {
	int p = u * ncol + v;
	if (p >= M.m_num_bits) {
		std::cerr << "has_edge::Error! " << u << " " << v << " " << M.m_num_bits
				<< std::endl;
	}
	return M[p];
}

void BinaryMatrix::flip(int u, int v) {
	int p = u * ncol + v;
	if (p >= M.m_num_bits) {
		std::cerr << "flip_edge::Error! " << u << " " << v << " "
				<< M.m_num_bits << std::endl;
	}
	M[p].flip();
}

bool BinaryMatrix::is_switchable(int u1, int v1, int u2, int v2) const {

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

void BinaryMatrix::switch_4_cycle(int u1, int v1, int u2, int v2) {

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

void BinaryMatrix::operator=(const BinaryMatrix& s) {
	nrow = s.nrow;
	ncol = s.ncol;
	M = boost::dynamic_bitset<>(s.M);
}

bool BinaryMatrix::operator==(const BinaryMatrix &s) const {
	return nrow == s.nrow && ncol == s.ncol && M == s.M;
}

bool BinaryMatrix::operator<(const BinaryMatrix &s) const {
	return nrow == s.nrow && ncol == s.ncol && M < s.M;
}

void BinaryMatrix::get_row(int u, boost::dynamic_bitset<>& row) const {
	row.resize(ncol);
	for (int v = 0; v < ncol; v++) {
		row[v] = get(u, v);
	}
}

size_t BinaryMatrix::hashValue() const {
	return boost::hash_value(M.m_bits);
}

int BinaryMatrix::compare(const State* x) const {
	const BinaryMatrix* d = (const BinaryMatrix*) x;
	if (operator==(*d))
		return 0;
	else {
		if (operator <(*d))
			return -1;
	}
	return 1;
}

State* BinaryMatrix::copy() const {
	return new BinaryMatrix(*this);
}

std::string BinaryMatrix::toString() const {
	std::stringstream ss;
	for (unsigned int i = 0; i < M.size(); i++)
		ss << M[i];
	return ss.str();
}

}
}
}

