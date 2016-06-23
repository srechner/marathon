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

	A(int id, int deg) :
			id(id), deg(deg) {
	}

	bool operator<(const A& a) const {
		return deg < a.deg;
	}
};

DenseBipartiteGraph* HavelHakimiBipartite(const std::vector<int>& u,
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
		std::sort(vv.begin(), vv.end());

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

	DenseBipartiteGraph* bip = new DenseBipartiteGraph(m, n, M);

	return bip;
}

DenseBipartiteGraph* HavelHakimiBipartiteForbidden(const std::vector<int>& u,
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
	memset(forbidden_colsum, 0, ncol* sizeof(int));
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

	std::vector<A> vv;
	for (int i = 0; i < v.size(); i++) {
		vv.push_back(A(i, v[i]));
	}

	// for each row
	for (int i = 0; i < nrow; i++) {

		std::sort(vv.begin(), vv.end());

		int ui = u[i];

		// for each one that has to be distributed in line i
		for (int j = ncol - 1; ui > 0; j--) {

			// if (i,vv[j].id) is not a forbidden entry
			if (!forbidden[i * ncol + vv[j].id]) {
				if (vv[j].deg == 0) {
					return nullptr;
				} else {
					M[i * ncol + vv[j].id] = 1;
					vv[j].deg--;
					ui--;
				}
			}
		}
	}

	// assert that nodes have zero degree now
	int sumv = 0;
	for (auto a : vv) {
		sumv += a.deg;
	}
	assert(sumv == 0);

	DenseBipartiteGraph* bip = new DenseBipartiteGraph(nrow, ncol, M);

	delete[] M;
	delete[] forbidden_rowsum;
	delete[] forbidden_colsum;

	return bip;
}

}
}
}
