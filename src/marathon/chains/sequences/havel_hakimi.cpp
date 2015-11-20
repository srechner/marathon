#include "../../../../include/marathon/chains/sequences/havel_hakimi.h"

#include <iostream>
#include <algorithm>
#include <cstring>

namespace marathon {

namespace chain {

namespace sequence {

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

bool HavelHakimiBipartite(const std::vector<int>& u,
		const std::vector<int>& v, bool* M) {

	int m = u.size();
	int n = v.size();

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
		return false;
	}

	// check if degrees have valid size
	for (std::vector<int>::const_iterator it = u.begin(); it != u.end();
			++it) {
		if (*it > v.size()) {
			//std::cerr << "Error 2! Degree Sequence not graphical!" << std::endl;
			return false;
		}
	}
	for (std::vector<int>::const_iterator it = v.begin(); it != v.end();
			++it) {
		if (*it > u.size()) {
			//std::cerr << "Error 3! Degree Sequence not graphical!" << std::endl;
			return false;
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
				return false;
			} else {
				M[i * n + vv[j].id] = 1;
				vv[j].deg--;
			}
		}
	}

	return true;
}

}

}

}
