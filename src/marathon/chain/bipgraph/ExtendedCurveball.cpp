/*
 * curveball2.cpp
 *
 *  Created on: Sep 7, 2015
 *      Author: rechner
 */

#include "../../../../include/marathon/chain/bipgraph/ExtendedCurveball.h"

namespace marathon {
namespace chain {
namespace bipgraph {

ExtendedCurveball::ExtendedCurveball(const std::string& line) :
		Curveball(line) {
}

void ExtendedCurveball::computeNeighbours(const State* x,
		std::vector<std::pair<State*, rational>>& neighbors) const {

	// Declare Variables
	int nrows, ncols, nTwoPartitions;
	boost::dynamic_bitset<>* A;

	const DenseBipartiteGraph* s = (const DenseBipartiteGraph*) x;

	nrows = s->get_nrows();
	ncols = s->get_ncols();

	/**
	 * Modification of the Curveball algorithm by Annabell Berger.
	 */

	// a) transform matrix into adjacency list (nothing to do)
	A = new boost::dynamic_bitset<>[nrows];
	for (int i = 0; i < nrows; i++)
		s->get_row(i, A[i]);

	// b) for each 2-partition of the set {0..n-1}
	std::vector<std::vector<std::pair<int, int>> > twoPartitons;
	backtrackTwoPartitions(nrows, twoPartitons);
	nTwoPartitions = twoPartitons.size();

	//std::cout << "found " << nTwoPartitions << " two-partitions" << std::endl;

	// for each two-partition P do in parallel:
#pragma omp parallel
	{
		// thread-local variables
		int d, i, j, k, l, nSubSets;
		boost::dynamic_bitset<> Aij, Aji, Bi, Bj, X;
		std::list<boost::dynamic_bitset<>> subsets;
		std::vector<std::pair<State*, rational>> myneighbors;

#pragma omp for
		for (l = 0; l < nTwoPartitions; l++) {

			std::vector<std::pair<int, int> > P = twoPartitons[l];

			// b) for each (i,j) in P
			for (std::pair<int, int> ij : P) {

				i = ij.first;
				j = ij.second;

				/*			std::cout << " i=" << i << ", j=" << j << std::endl;
				 std::cout << " A[i]=" << A[i] << std::endl;
				 std::cout << " A[j]=" << A[j] << std::endl;*/

				// c) compute intersections of edges
				Aij = A[i] - A[j];
				Aji = A[j] - A[i];
				d = Aij.count();

				//std::cout << " A[i]-A[j]=" << Aij << std::endl;
				//std::cout << " A[j]-A[i]=" << Aji << std::endl;

				// d) remove A[i]-A[j] from A[i], which are d elements being removed...
				// and fill with the same number of elements from A[i]^A[j].
				// Therefore, backtrack all subsets of A[i]^A[j] of size d
				X = A[i] ^ A[j];

				//std::cout << " A[i]^A[j]=" << X << std::endl;

				subsets.clear();
				backtrackSubsets(X, d, subsets);
				nSubSets = subsets.size();

				const rational p(1, nTwoPartitions * nSubSets * P.size());

				/*std::cout << " found " << nSubSets << " subsets of size " << d
				 << std::endl;*/

				// for each subset
				for (boost::dynamic_bitset<> sub : subsets) {

					//std::cout << "  " << sub << std::endl;

					// fill Bi with elements from sub
					Bi = A[i] - Aij;
					Bi |= sub;

					// Fill Bj with remaining elements
					Bj = A[j] - Aji;
					Bj |= X - sub;

					DenseBipartiteGraph *s2 = new DenseBipartiteGraph(*s);
					for (k = 0; k < ncols; k++) {
						s2->set_edge(i, k, Bi[k]);
						s2->set_edge(j, k, Bj[k]);
					}

					/*std::cout << " new neighbour: " << s2
					 << " with probability of " << p << std::endl;*/

					myneighbors.push_back(std::make_pair(s2, p));
				}
			}
		}
		// copy thread-local variables to global neighbour list
#pragma omp critical
		{
			neighbors.insert(neighbors.end(), myneighbors.begin(),
					myneighbors.end());
		}

	}

	delete[] A;
}

void ExtendedCurveball::backtrackTwoPartitions(int n,
		std::vector<std::vector<std::pair<int, int> > >& partitions) const {

// allocate working memory
	int* a = new int[n];
	int* b = new int[n];
	memset(a, 0, n * sizeof(int));
	memset(b, 0, n * sizeof(int));
	backtrackTwoPartitionsRecursive(n, a, b, 0, 0, partitions);
	delete[] a;
	delete[] b;
}

// just another recursive backtracking procedure
void ExtendedCurveball::backtrackTwoPartitionsRecursive(int n, int* a, int* b, int k,
		int pos,
		std::vector<std::vector<std::pair<int, int> > >& partitions) const {
	if (n - 2 * k <= 1) {
		// found new two-partition
		std::vector<std::pair<int, int>> tmp;

		// convert to vector of pairs
		for (int i = 0; i < 2 * k; i += 2) {
			tmp.push_back(std::make_pair(b[i], b[i + 1]));
		}

		// store solution to result list
		partitions.push_back(tmp);
	} else {
		int i = pos;
		while (i < n) {
			if (a[i] == 0) {
				a[i] = k + 1;
				b[2 * k] = i;
				int j = i + 1;
				while (j < n) {
					if (a[j] == 0) {
						a[j] = k + 1;
						b[2 * k + 1] = j;
						backtrackTwoPartitionsRecursive(n, a, b, k + 1, i + 1,
								partitions);
						a[j] = 0;
					}
					j++;
				}
				a[i] = 0;
			}
			i++;
		}
	}
}

}

}

}
