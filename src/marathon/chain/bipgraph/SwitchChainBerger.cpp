/*
 * chain_swap_bipartite_fast.cpp
 *
 *  Created on: Jun 4, 2015
 *      Author: rechner
 */

#include <iostream>

#include "../../../../include/marathon/chain/bipgraph/SwitchChainBerger.h"

namespace marathon {
namespace chain {
namespace bipgraph {

SwitchChainBerger::SwitchChainBerger(const std::string& input) :
		SwitchChain(input) {
}

void SwitchChainBerger::computeNeighbours(const State* x,
		std::vector<std::pair<State*, rational>>& neighbors) const {

	uint i, j, k, l, m, n;
	DenseBipartiteGraph s2;
	int numNonAdjacentEdgePairs = 0;
	rational loop(0);

	const DenseBipartiteGraph* s = (const DenseBipartiteGraph*) x;

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

	m = s->get_nrows();
	n = s->get_ncols();

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
					if (s->has_edge(i, j) && s->has_edge(k, l)) {

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

	const rational p(1, numNonAdjacentEdgePairs);

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
						if (s->has_edge(i, j) && s->has_edge(k, l)) {

							DenseBipartiteGraph* s2 = new DenseBipartiteGraph(
									*s);

							// alternating cycle ( (i,j)=1, (k,j)=0, (k,l)=1, (i,l)=0 )
							if (!s->has_edge(i, l) && !s->has_edge(k, j)) {

								// switch the cycle
								s2->flip_edge(i, j); // (i,j) = 0
								s2->flip_edge(k, l); // (k,l) = 0
								s2->flip_edge(i, l); // (i,l) = 1
								s2->flip_edge(k, j); // (k,j) = 1
								//std::cout << "a -> " << s2 << " (" << s << ")" << std::endl;
							}

							neighbors.push_back(std::make_pair(s2, p));
						}
					}
				}
			}
		}
	}

	// choose second edge as artificial edge
	loop += rational(sum, numNonAdjacentEdgePairs);

	DenseBipartiteGraph *sl = new DenseBipartiteGraph(*s);
	neighbors.push_back(std::make_pair(sl, loop));
}

}
}
}
