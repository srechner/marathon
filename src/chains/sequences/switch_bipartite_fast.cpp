/*
 * chain_swap_bipartite_fast.cpp
 *
 *  Created on: Jun 4, 2015
 *      Author: rechner
 */

#include <iostream>
#include "../../../include/marathon/chains/sequences/switch_chain_bipartite_berger.h"

namespace marathon {

namespace chain {

namespace sequence {

SwitchBipartiteFast::SwitchBipartiteFast(std::string line) :
		SwitchBipartite(line) {
}

void SwitchBipartiteFast::computeNeighbours(const DenseBipartiteGraph& s,
		boost::unordered_map<DenseBipartiteGraph, Rational>& neighbors) const {
	uint i, j, k, l, m, n;
	DenseBipartiteGraph s2;
	int numNonAdjacentEdgePairs = 0;

	/**
	 * Definition of Annabell Berger
	 *
	 * Let M be a set of all non-adjacent pairs of edges of s,
	 * extended by the artificial edge pair (a,a) which is necessary
	 * to get a ergodic chain.
	 *
	 * Choose a pair (e,e') of edges from M uniformly and random.
	 * Let these pairs be e=(i,j) and e'=(k,l) with i!=k and j!=l.
	 * In Practice this would be done by an efficient data strutcure.
	 *
	 * If the edges (i,l) and (k,j) are missing, switch the cycle.
	 */

	m = s.get_nrows();
	n = s.get_ncols();

	/***************************************************************
	 * First Step: Count number of "natural" non-adjacent edge pairs
	 **************************************************************/

	// choose edge e=(i,j)
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {

			// Choose edge e'=(k,l)
			for (k = i + 1; k < m; k++) {
				for (l = 0; l < n; l++) {

					// if e and e' are actual edges
					if (s.has_edge(i, j) && s.has_edge(k, l)) {

						// if e and e' are non adjacent
						if (l != j)
							numNonAdjacentEdgePairs++;

					}
				}
			}
		}
	}

	// add an additional edge a, which is non-adjacent to all
	numNonAdjacentEdgePairs += sum;

	/*****************************************************
	 * Second Step: Select two non-adjacent edge pairs
	 ****************************************************/

	// choose edge e=(i,j)
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {

			// choose second edge e'=(k,l), different from e
			for (k = i + 1; k < m; k++) {
				for (l = 0; l < n; l++) {

					// if e and e' are non adjacent
					if (l != j) {

						// if e and e' are actual edges
						if (s.has_edge(i, j) && s.has_edge(k, l)) {

							s2 = DenseBipartiteGraph(s);

							// alternating cycle ( (i,j)=1, (k,j)=0, (k,l)=1, (i,l)=0 )
							if (!s.has_edge(i, l) && !s.has_edge(k, j)) {

								// switch the cycle
								s2.flip_edge(i, j); // (i,j) = 0
								s2.flip_edge(k, l); // (k,l) = 0
								s2.flip_edge(i, l); // (i,l) = 1
								s2.flip_edge(k, j); // (k,j) = 1
								//std::cout << "a -> " << s2 << " (" << s << ")" << std::endl;
							}

							neighbors[s2] += Rational(1, numNonAdjacentEdgePairs);
						}
					}
				}
			}
		}
	}

	// choose second edge as artificial edge
	neighbors[s] += Rational(sum, numNonAdjacentEdgePairs);
}

}

}

}
