/*
 * chain_swap_bipartite_fast.cpp
 *
 *  Created on: Jun 4, 2015
 *      Author: rechner
 */

#include "../../../include/marathon/chains/sequences/switch_chain_bipartite_fast.h"
#include <iostream>

namespace marathon {

namespace chain {

namespace sequence {

SwapChainBipartiteFast::SwapChainBipartiteFast(std::string line) :
		SwapChainBipartite(line) {

}

void SwapChainBipartiteFast::computeNeighbors(const DenseBipartiteGraph& s,
		std::vector<DenseBipartiteGraph>& neighbors) const {
	uint i, j, k, l, m, n;
	DenseBipartiteGraph s2;

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

	m = s.get_m();
	n = s.get_n();

	// choose edge e=(i,j)
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {

			// Choose edge e'=(k,l)
			for (k = i + 1; k < m; k++) {
				for (l = 0; l < n; l++) {

					if (l != j) {

						// only non-adjacent pairs of edges are considered
						if (s.has_edge(i, j) && s.has_edge(k, l)) {

							//std::cout << "Fall 1";

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

							neighbors.push_back(s2);
						}
					}
				}
			}
		}
	}

	/**
	 * Consider artificial edge pair (a,a).
	 * This automatically results to a single self-loop.
	 */
	s2 = DenseBipartiteGraph(s);
	neighbors.push_back(s2);
}

}

}

}
