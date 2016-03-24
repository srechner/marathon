/*
 * JS89Chain.cpp
 *
 *  Created on: Feb 11, 2015
 *      Author: steffen
 */

#include "../../../../include/marathon/chain/matching/Broder86.h"

namespace marathon {
namespace chain {
namespace matching {

Broder86::Broder86(const std::string& inst) :
			MarkovChain(inst) {
	// construct input object from input instance
	parseInstance(inst);
}

Broder86::~Broder86() {
	if (g != nullptr)
		delete g;
}

void Broder86::parseInstance(const std::string& inst) {
	if (g != nullptr)
		delete g;
	g = new SparseBipartiteGraph(inst);
}

State* Broder86::computeArbitraryState() {

	int n = g->getNumberOfNodes();

	std::vector<int> mate;
	g->cardmax_matching(mate);

	// Count Number of Matching Edges
	int k = 0;
	for (int i = 0; i < mate.size(); i++)
		if (mate[i] != n)
			k++;

	if (k != g->getNumberOfNodes()) {
		std::cerr << "Graph doesn't have a perfect matching!" << std::endl;
		return nullptr;
	}

	int unmatched[2] = { n, n };
	State* s = new BipartiteMatching(g->getNumberOfNodes(), k / 2, unmatched,
			&mate[0]);

	return s;
}

void Broder86::computeNeighbours(const State* x,
		std::vector<std::pair<State*, rational>>& neighbors) const {

	// Variables
	const BipartiteMatching* s = (const BipartiteMatching*) x;
	const int n = g->getNumberOfNodes();
	const int k = s->k;
	int u, v;
	edgelist edges;

	neighbors.clear();

	// Implementierung der Transitionen nach Vorlage JS89

	// Gehe über jede Kante
	g->getEdges(edges);

	const rational p(1, edges.size());

	//std::cout << "compute neighbors for s=" << *s << std::endl;

	// Für jede Wahl einer Kante e=(u,v) wird Transition konstruiert
	for (typename edgelist::iterator it = edges.begin(); it != edges.end();
			++it) {

		u = it->first;
		v = it->second;

		// create new state as copy of old one
		BipartiteMatching* s2 = new BipartiteMatching(*s);

		// Transition 1
		if (2 * k == n && s->mates[u] == v) {
			// remove edge (u,v)
			s2->mates[u] = n;
			s2->mates[v] = n;
			s2->k = k - 1;
			s2->unmatched[0] = u;
			s2->unmatched[1] = v;
		} else if (2 * k + 2 == n) {

			// Transition 2
			if (s->mates[u] == n && s->mates[v] == n) {
				// add edge (u,v)
				s2->mates[u] = v;
				s2->mates[v] = u;
				s2->k = k + 1;
				s2->unmatched[0] = n;
				s2->unmatched[1] = n;
			}
			// Transition 3a
			else if (s->mates[u] != n && s->mates[v] == n) {
				// remove edge (u, mate[u])
				// add edge (u,v)
				int w = s->mates[u];
				s2->mates[w] = n;
				s2->mates[u] = v;
				s2->mates[v] = u;
				// w=mate[u] (>= n/2) becomes unmatched node in bipartition group 1
				s2->unmatched[1] = w;
			}
			// Transition 3b
			else if (s->mates[u] == n && s->mates[v] != n) {
				// remove edge (v, mate[v])
				// add edge (u,v)
				int w = s->mates[v];
				s2->mates[w] = n;
				s2->mates[u] = v;
				s2->mates[v] = u;
				// w=mate[v] (< n/2) becomes unmatched node in bipartition group 0
				s2->unmatched[0] = w;
			} else {
				// stay in s
			}
		}
		// sonst
		else {
			// Verbleibe im aktuellen Zustand
		}

		neighbors.push_back(std::make_pair(s2, p));
	}
}

}
}
}
