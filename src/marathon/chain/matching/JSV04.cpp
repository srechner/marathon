/*
 * JSV04.cpp
 *
 * Created on: Aug 24, 2015
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

#include "../../../../include/marathon/chain/matching/JSV04.h"

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

void JerrumSinclairVigoda04::computeNeighbours(const State* xs,
		std::vector<std::pair<State*, rational>>& neighbors) const {

	// Variables
	const BipartiteMatching* s = (const BipartiteMatching*) xs;
	const uint n = Broder86::g->getNumberOfNodes();
	const uint k = s->k;
	uint u, v, z, x, y;
	edgelist edges;

	const rational p(2, n);
	const rational q(1, n);

	neighbors.clear();

	// Implementierung der Transitionen nach Jerrum, Sinclair, Vigoda 2004

	if (s->is_perfect()) {	// if s is perfect matching

		// Choose matching edge e=(u,v) and remove it
		for (u = 0; u < s->n / 2; u++) {
			v = s->mates[u];
			BipartiteMatching *s2 = new BipartiteMatching(*s);
			s2->removeEdge(u, v);
			neighbors.push_back(std::make_pair(s2, p));
		}

	} else if (s->is_near_perfect()) {	// is a near-perfect matching

		u = s->unmatched[0];
		v = s->unmatched[1];

		// choose node z
		for (z = 0; z < n; z++) {

			BipartiteMatching* s2 = new BipartiteMatching(*s);

			// Three Cases to rule them all

			if ((z == u || z == v) && Broder86::g->hasEdge(u, v)) {
				// Case 2i
				s2->addEdge(u, v);
			} else if (z >= n / 2 && Broder86::g->hasEdge(u, z)
					&& s->mates[z] != n) {	// Case 2ii
				x = s->mates[z];
				s2->removeEdge(x, z);
				s2->addEdge(u, z);
				s2->unmatched[0] = x;
				s2->unmatched[1] = v;
			} else if (z < n / 2 && Broder86::g->hasEdge(z, v)
					&& s->mates[z] != n) {	// Case 2iii
				y = s->mates[z];
				s2->removeEdge(z, y);
				s2->addEdge(z, v);
				s2->unmatched[0] = u;
				s2->unmatched[1] = y;
			} else {
				// loop
			}

			neighbors.push_back(std::make_pair(s2, q));
		}
	}
}

void JerrumSinclairVigoda04::computeWeights(
		const std::vector<const State*>& states,
		std::vector<rational>& weights) {

	const uint n = Broder86::g->getNumberOfNodes();

	num_perfect_matching = 0;
	num_near_perfect_matching = (uint*) realloc(num_near_perfect_matching,
			n * n / 4 * sizeof(uint));

	memset(num_near_perfect_matching, 0, n * n / 4 * sizeof(uint));

	// count_recursive the number of u-v-partitions
	for (const State* x : states) {

		const BipartiteMatching* s = (const BipartiteMatching*) x;

		int u = s->unmatched[0];
		int v = s->unmatched[1] - n / 2;

		if (s->is_perfect())
			num_perfect_matching++;
		else if (s->is_near_perfect())
			num_near_perfect_matching[u * (n / 2) + v]++;
	}

	weights.resize(states.size());
	for (int i = 0; i < states.size(); i++)
		weights[i] = getWeight(states[i]);
}

rational JerrumSinclairVigoda04::getWeight(const State* x) const {

	const BipartiteMatching* s = (const BipartiteMatching*) x;

	if (s->is_perfect()) {
		return 1;
	} else if (s->is_near_perfect()) {
		uint u = s->unmatched[0];
		uint v = s->unmatched[1] - s->n / 2;
		rational r(num_perfect_matching,
				num_near_perfect_matching[u * s->n / 2 + v]);
		return r;
	} else {
		// cannot happen
		return -1;
	}
}

}
}
}
