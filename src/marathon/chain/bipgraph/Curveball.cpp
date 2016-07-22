/*
 * swap.cpp
 *
 *  Created on: Jun 4, 2015
 *      Author: rechner
 */

#include <iostream>
#include <map>

#include "../../../../include/marathon/chain/bipgraph/Curveball.h"
#include "../../../../include/marathon/Random.h"
#include "../../../../include/marathon/Combinatorics.h"

namespace marathon {
namespace chain {
namespace bipgraph {

Curveball::Curveball(const std::string& line) :
		SwitchChain(line) {

}

Curveball::Curveball(const std::vector<int>& u, const std::vector<int>& v) :
		SwitchChain(u, v) {

}

Curveball::Curveball(int* const rowsum, int* const colsum, const int nrow,
		const int ncol) :
		SwitchChain(rowsum, colsum, nrow, ncol) {

}

void Curveball::computeNeighbours(const State* x,
		std::vector<std::pair<State*, rational>>& neighbors) const {

	// Declare Variables
	const BinaryMatrix* s = (const BinaryMatrix*) x;
	const int nrows = s->getNumRows();
	boost::dynamic_bitset<>* A;

	/**
	 * Definition of Strona et. al: A fast and unbiased procedure
	 * to randomize ecological binary matrices with fixed row and
	 * column totals.
	 *
	 */

	// a) transform matrix into adjacency list (nothing to do)
	A = new boost::dynamic_bitset<>[nrows];
	for (int i = 0; i < nrows; i++)
		s->get_row(i, A[i]);

	boost::dynamic_bitset<> Aij, Aji, Bi, Bj, X;

	// b) select two different row indices
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < nrows; j++) {

			if (i != j) {

				/*std::cout << "i=" << i << ", j=" << j << std::endl;
				 std::cout << "A[i]=" << A[i] << std::endl;
				 std::cout << "A[j]=" << A[j] << std::endl;*/

				// c) compute intersections of edges
				Aij = A[i] - A[j];
				Aji = A[j] - A[i];
				int d = Aij.count();

				/*std::cout << "A[i]-A[j]=" << Aij << std::endl;
				 std::cout << "A[j]-A[i]=" << Aji << std::endl;*/

				// d) remove A[i]-A[j] from A[i], which are d elements being removed...
				// and fill with the same number of elements from A[i]^A[j].
				// Therefore, backtrack all subsets of A[i]^A[j] of size d
				X = A[i] ^ A[j];

				//std::cout << "A[i]^A[j]=" << X << std::endl;

				// determine number of subsets
				const rational numSubsets = ::marathon::combinatorics::binom(
						X.count(), Aij.count());

				// proposal probability of each adjacent state
				const rational p = rational(1)
						/ (rational(nrows) * rational(nrows - 1) * numSubsets);

				// prepare line i and j of adjacent matrix
				Bi = A[i] - Aij;
				Bj = A[j] - Aji;

				// recursively enumerate all subsets and add a adjacent state for each subset
				backtrackSubsets(X, d, p, Bi, Bj, X, s, i, j, neighbors);
			}
		}
	}

	delete[] A;
}

void Curveball::backtrackSubsets(boost::dynamic_bitset<>& population, int n,
		const rational& p, const boost::dynamic_bitset<>& Bi,
		const boost::dynamic_bitset<>& Bj, const boost::dynamic_bitset<>& X,
		const BinaryMatrix* s, const int i, const int j,
		std::vector<std::pair<State*, rational>>& myneighbours) const {

	boost::dynamic_bitset<> tmp(population.size(), 0);
	backtrackSubsetsRecursive(population, tmp, 0, 0, n, Bi, Bj, X, s, i, j, p,
			myneighbours);

}

void Curveball::backtrackSubsetsRecursive(
		const boost::dynamic_bitset<>& population, boost::dynamic_bitset<>& sub,
		int pos, int m, int n, const boost::dynamic_bitset<>& Bi,
		const boost::dynamic_bitset<>& Bj, const boost::dynamic_bitset<>& X,
		const BinaryMatrix* s, const int i, const int j,
		const rational& p,
		std::vector<std::pair<State*, rational>>& myneighbours) const {

	// stop of recursion
	if (m == n) {

		//std::cout << "  " << sub << std::endl;

		const int ncols = s->getNumColumns();

		// fill Bi with elements from sub
		const boost::dynamic_bitset<> Ci = Bi | sub;

		// Fill Bj with remaining elements
		const boost::dynamic_bitset<> Cj = Bj | (X - sub);

		BinaryMatrix *s2 = new BinaryMatrix(*s);
		for (int k = 0; k < ncols; k++) {
			s2->set(i, k, Ci[k]);
			s2->set(j, k, Cj[k]);
		}
		myneighbours.push_back(std::make_pair(s2, p));
		return;
	} else if (pos >= population.size()) {
		return;
	} else {
		while (pos < population.size()) {
			if (population[pos]) {
				sub[pos] = 1;
				backtrackSubsetsRecursive(population, sub, pos + 1, m + 1, n,
						Bi, Bj, X, s, i, j, p, myneighbours);
				sub[pos] = 0;
			}
			pos++;
		}
	}
}

void Curveball::randomize(State* x, const uint32_t t) const {

	BinaryMatrix* s = (BinaryMatrix*) x;

	const int nrows = u.size();
	const int ncols = v.size();

	if (nrows == 1 || ncols == 1) {
		return;
	}

	for (int l = 0; l < t; l++) {

		// randomly select two row indices
		int tmp[ncols];
		::marathon::random::select(nrows, 2, tmp);
		int i = tmp[0];
		int j = tmp[1];

		if (i >= nrows || j >= nrows) {
			std::cerr << "error! " << i << " " << j << std::endl;
		}

		// select the indices that occur in on the Mi and Mj, but not in both
		int X[ncols];
		int a = 0;
		int b = 0;

		// for each column position
		for (int k = 0; k < ncols; k++) {

			bool A_ik = s->get(i, k);
			bool A_jk = s->get(j, k);

			if (A_ik != A_jk) {

				// store index k at X array
				X[a + b] = k;

				if (A_ik)
					a++;
				else
					b++;
			}
		}

		// select a out of (a+b) positions
		bool Y[a + b];
		memset(Y, false, (a + b) * sizeof(bool));
		::marathon::random::select(a + b, a, tmp);
		for (int k = 0; k < a; k++)
			Y[tmp[k]] = true;

		// for each integer in X
		for (int k = 0; k < a + b; k++) {
			s->set(i, X[k], Y[k]);
			s->set(j, X[k], !Y[k]);
		}
	}
}

marathon::rational Curveball::loopProbability(const State* x) const {

	const BinaryMatrix* s = (const BinaryMatrix*) x;

	rational loop(0);

	const int nrows = u.size();
	const int ncols = v.size();

	// select two different row indices
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < nrows; j++) {

			if (i != j) {

				// count number of edges that are in A_i but not in A_j (and vice versa)
				int cij = 0;
				int cji = 0;

				// for each column position
				for (int k = 0; k < ncols; k++) {

					bool A_ik = s->get(i, k);
					bool A_jk = s->get(j, k);

					if (A_ik != A_jk) {

						if (A_ik)
							cij++;
						else
							cji++;
					}
				}

				// determine number of subsets
				const rational numSubsets = ::marathon::combinatorics::binom(
						cij + cji, cij);

				// proposal probability of each adjacent state
				const rational p = rational(1)
						/ (rational(nrows) * rational(nrows - 1) * numSubsets);

				// exactly one subset results in a loop
				loop += p;
			}
		}
	}

	return loop;
}

}
}
}
