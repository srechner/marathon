/*
 * HavelHakimi.cpp
 *
 * Created on: Sep 7, 2015
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

#include <iostream>
#include <algorithm>
#include <cstring>

#include "../../../../include/marathon/chain/bipgraph/HavelHakimi.h"

namespace marathon {
namespace chain {
namespace bipgraph {

struct A {
	int id;
	int deg;
	int numForbidden;

	A(int id, int deg, int numForbidden = 0) :
			id(id), deg(deg), numForbidden(numForbidden) {
	}
};

BinaryMatrix* HavelHakimiBipartite(const std::vector<int>& u,
		const std::vector<int>& v) {

	int m = u.size();
	int n = v.size();

	bool M[u.size() * v.size()];
	memset(M, 0, m * n * sizeof(bool));

	// sum of node degrees
	int sum[2] = { 0, 0 };

	// check if node degrees sum up to same number
	for (std::vector<int>::const_iterator it = u.begin(); it != u.end(); ++it)
		sum[0] += *it;
	for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it)
		sum[1] += *it;
	if (sum[0] != sum[1]) {
		//std::cerr << "Error 1! Degree Sequence not graphical!" << std::endl;
		return nullptr;
	}

	// check if degrees have valid size
	for (std::vector<int>::const_iterator it = u.begin(); it != u.end(); ++it) {
		if (*it > v.size()) {
			//std::cerr << "Error 2! Degree Sequence not graphical!" << std::endl;
			return nullptr;
		}
	}

	for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it) {
		if (*it > u.size()) {
			//std::cerr << "Error 3! Degree Sequence not graphical!" << std::endl;
			return nullptr;
		}
	}

	std::vector<A> vv;
	for (int i = 0; i < v.size(); i++) {
		vv.push_back(A(i, v[i]));
	}

	for (int i = 0; i < m; i++) {

		// sort by degree
		std::sort(vv.begin(), vv.end(), [](const A& a, const A& b) {
			return a.deg < b.deg;
		});

		for (int j = n - u[i]; j < n; j++) {
			if (vv[j].deg == 0) {
				//std::cerr << "Error 4! Degree Sequence not graphical!"
				//		<< std::endl;
				return nullptr;
			} else {
				M[i * n + vv[j].id] = 1;
				vv[j].deg--;
			}
		}
	}

	BinaryMatrix* bip = new BinaryMatrix(m, n, M);

	return bip;
}

BinaryMatrix* HavelHakimiBipartiteForbidden(const std::vector<int>& u,
		const std::vector<int>& v, bool* const forbidden) {

	int nrow = u.size();
	int ncol = v.size();

	// sum of node degrees
	int sum[2] = { 0, 0 };

	// binary matrix
	bool* M = new bool[nrow * ncol];
	memset(M, 0, nrow * ncol * sizeof(bool));

	// determine the number of forbidden entries in each row and col
	int* forbidden_rowsum = new int[nrow];
	int* forbidden_colsum = new int[ncol];
	memset(forbidden_rowsum, 0, nrow * sizeof(int));
	memset(forbidden_colsum, 0, ncol * sizeof(int));
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			if (forbidden[i * ncol + j]) {
				forbidden_rowsum[i]++;
				forbidden_colsum[j]++;
			}
		}
	}

	// check if node degrees sum up to same number
	if (sum[0] != sum[1]) {
		std::cerr << "Error: sums are not equal!" << std::endl;
		return nullptr;
	}

	// check if degrees have valid size
	for (int i = 0; i < nrow; i++) {
		if (u[i] > ncol - forbidden_rowsum[i]) {
			std::cerr << "Error: not enough ones in row " << i << "!"
					<< std::endl;
			return nullptr;
		}
	}
	for (int j = 0; j < ncol; j++) {
		if (v[j] > nrow - forbidden_colsum[j]) {
			std::cerr << "Error: not enough ones in col " << j << "!"
					<< std::endl;
			return nullptr;
		}
	}

	// start realizing

	std::vector<A> uu;
	for (int i = 0; i < u.size(); i++) {
		// determine number of forbidden elements in row i
		int numForbidden = 0;
		for (int j = 0; j < ncol; j++) {
			if (forbidden[i * ncol + j]) {
				numForbidden++;
			}
		}
		uu.push_back(A(i, u[i], numForbidden));
	}

	// descendlingly sort uu by number of forbidden elements
	std::sort(uu.begin(), uu.end(), [](const A& a, const A& b) {
		return a.numForbidden > b.numForbidden;
	});

	std::vector<A> vv;
	for (int i = 0; i < v.size(); i++) {
		vv.push_back(A(i, v[i]));
	}

	// for each row
	for (int i = 0; i < nrow; i++) {

		// sort by degree
		std::sort(vv.begin(), vv.end(), [](const A& a, const A& b) {
			return a.deg < b.deg;
		});

		// for each one that has to be distributed in line i
		for (int j = ncol - 1; uu[i].deg > 0; j--) {

			// if (i,vv[j].id) is not a forbidden entry
			if (!forbidden[uu[i].id * ncol + vv[j].id]) {
				if (vv[j].deg == 0) {
					return nullptr;
				} else {
					M[uu[i].id * ncol + vv[j].id] = 1;
					vv[j].deg--;
					uu[i].deg--;
				}
			}
		}
	}

	// assert that nodes have zero degree now
	for (auto a : uu) {
		assert(a.deg == 0);
	}
	for (auto a : vv) {
		assert(a.deg == 0);
	}

	BinaryMatrix* bip = new BinaryMatrix(nrow, ncol, M);

	delete[] M;
	delete[] forbidden_rowsum;
	delete[] forbidden_colsum;

	return bip;
}

}
}
}
