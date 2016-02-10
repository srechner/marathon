/*
 * jerrum_sinclair_vigoda04.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: rechner
 */

#include "../../../../include/marathon/chains/matching/matching_chain_JSV04.h"

namespace marathon {

namespace chain {

namespace matching {

JerrumSinclairVigoda04::JerrumSinclairVigoda04(const std::string& input) :
		Broder86(input), num_perfect_matching(0) {
	num_near_perfect_matching = (uint*) malloc(0);
}

JerrumSinclairVigoda04::~JerrumSinclairVigoda04() {
	free(num_near_perfect_matching);
}

void JerrumSinclairVigoda04::computeNeighbours(const BipartiteMatching& s,
		std::unordered_map<BipartiteMatching, rational>& neighbors) const {

	// Variables
	const uint n = Broder86::g->getNumberOfNodes();
	const uint k = s.k;
	uint u, v, z, x, y;
	BipartiteMatching s2;
	edgelist edges;

	neighbors.clear();

	// Implementierung der Transitionen nach Jerrum, Sinclair, Vigoda 2004

	if (s.is_perfect()) {	// if s is perfect matching

		// Choose matching edge e=(u,v) and remove it
		for (u = 0; u < s.n / 2; u++) {
			v = s.mates[u];
			s2 = BipartiteMatching(s);
			s2.removeEdge(u, v);
			neighbors[s2] += rational(2, n);
		}

	} else if (s.is_near_perfect()) {	// is a near-perfect matching

		u = s.unmatched[0];
		v = s.unmatched[1];

		// choose node z
		for (z = 0; z < n; z++) {

			s2 = BipartiteMatching(s);

			// Three Cases to rule them all

			if ((z == u || z == v) && Broder86::g->hasEdge(u, v)) {
				// Case 2i
				s2.addEdge(u, v);
				neighbors[s2] += rational(1, n);
			} else if (z >= n / 2 && Broder86::g->hasEdge(u, z)
					&& s.mates[z] != n) {	// Case 2ii
				x = s.mates[z];
				s2.removeEdge(x, z);
				s2.addEdge(u, z);
				s2.unmatched[0] = x;
				s2.unmatched[1] = v;
				neighbors[s2] += rational(1, n);
			} else if (z < n / 2 && Broder86::g->hasEdge(z, v)
					&& s.mates[z] != n) {	// Case 2iii
				y = s.mates[z];
				s2.removeEdge(z, y);
				s2.addEdge(z, v);
				s2.unmatched[0] = u;
				s2.unmatched[1] = y;
				neighbors[s2] += rational(1, n);
			} else {
				// loop
				neighbors[s2] += rational(1, n);
			}
		}
	}
}

void JerrumSinclairVigoda04::computeWeights(
		std::vector<BipartiteMatching>& states,
		std::vector<rational>& weights) {

	const uint n = Broder86::g->getNumberOfNodes();
	uint u, v;
	typename std::vector<BipartiteMatching>::iterator it;

	num_perfect_matching = 0;
	num_near_perfect_matching = (uint*) realloc(num_near_perfect_matching,
			n * n / 4 * sizeof(uint));

	memset(num_near_perfect_matching, 0, n * n / 4 * sizeof(uint));

	// count the number of u-v-partitions
	for (it = states.begin(); it != states.end(); ++it) {

		u = it->unmatched[0];
		v = it->unmatched[1] - n / 2;

		if (it->is_perfect())
			num_perfect_matching++;
		else if (it->is_near_perfect())
			num_near_perfect_matching[u * (n / 2) + v]++;
	}

	weights.resize(states.size());
	for (int i = 0; i < states.size(); i++)
		weights[i] = getWeight(states[i]);
}

rational JerrumSinclairVigoda04::getWeight(const BipartiteMatching& s) const {

	if (s.is_perfect()) {
		return 1;
	} else if (s.is_near_perfect()) {
		uint u = s.unmatched[0];
		uint v = s.unmatched[1] - s.n / 2;
		rational r(num_perfect_matching,
				num_near_perfect_matching[u * s.n / 2 + v]);
		return r;
	} else {
		// cannot happen
		return -1;
	}
}

}

}

}
