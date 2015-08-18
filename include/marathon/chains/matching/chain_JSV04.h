/*
 * JS89.h
 *
 *  Created on: Nov 18, 2014
 *      Author: rechner
 */

#ifndef JSV04_CHAIN_H_
#define JSV04_CHAIN_H_

#include "chain_JS89.h"

#include <queue>

namespace marathon {

namespace chain {

namespace matching {

class MatchingChain04: public MatchingChain89 {

private:
	uint num_perfect_matching;
	uint *num_near_perfect_matching;

public:

	MatchingChain04(std::string line) :
			MatchingChain89(line), num_perfect_matching(0) {
		num_near_perfect_matching = (uint*) malloc(0);
	}

	~MatchingChain04() {
		free(num_near_perfect_matching);
	}

	void compute_neighbors(const BipartiteMatching& s,
			std::vector<BipartiteMatching>& neighbors) const {

		// Variables
		const uint n = MatchingChain89::g.getNumberOfNodes();
		const uint k = s.k;
		uint u, v, z, x, y;
		BipartiteMatching s2;
		edgelist edges;

		neighbors.clear();

		// Implementierung der Transitionen nach Jerrum, Sinclair, Vigoda 2004

		if (2 * k == n) {	// if s is perfect matching

			// Choose matching edge e=(u,v) and remove it
			for (u = 0; u < s.n / 2; u++) {
				v = s.mates[u];
				s2 = BipartiteMatching(s);
				s2.removeEdge(u, v);
				neighbors.push_back(s2);
			}

			for (v = s.n / 2; v < s.n; v++) {
				u = s.mates[v];
				s2 = BipartiteMatching(s);
				s2.removeEdge(u, v);
				neighbors.push_back(s2);
			}

		} else if (2 * k == n - 2) {	// is a near-perfect matching

			u = s.unmatched[0];
			v = s.unmatched[1];

			// choose node z
			for (z = 0; z < n; z++) {

				s2 = BipartiteMatching(s);

				// Three Cases to rule them all

				if ((z == u || z == v) && MatchingChain89::g.hasEdge(u, v)) {// Case 2i
					s2.addEdge(u, v);
				} else if (z >= n / 2 && MatchingChain89::g.hasEdge(u, z)
						&& s.mates[z] != n) {	// Case 2ii
					x = s.mates[z];
					s2.removeEdge(x, z);
					s2.addEdge(u, z);
					s2.unmatched[0] = x;
					s2.unmatched[1] = v;
				} else if (z < n / 2 && MatchingChain89::g.hasEdge(z, v)
						&& s.mates[z] != n) {	// Case 2iii
					y = s.mates[z];
					s2.removeEdge(z, y);
					s2.addEdge(z, v);
					s2.unmatched[0] = u;
					s2.unmatched[1] = y;
				}

				neighbors.push_back(s2);
			}
		}
	}

	void compute_weights(std::vector<BipartiteMatching>& states) {

		const uint n = MatchingChain89::g.getNumberOfNodes();
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
				num_near_perfect_matching[u * n / 2 + v]++;
		}

	}

	Rational getWeight(uint i) const {

		const BipartiteMatching s = states[i];

		if (s.is_perfect())
			return 1;
		else if (s.is_near_perfect()) {
			uint u = s.unmatched[0];
			uint v = s.unmatched[1] - s.n / 2;
			return Rational(num_perfect_matching,
					num_near_perfect_matching[u * s.n / 2 + v]);
		} else {
			// cannot happen
			return -1;
		}
	}
};

}

}

}

#endif /* JS89_H_ */
