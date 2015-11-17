/*
 * JS89Chain.cpp
 *
 *  Created on: Feb 11, 2015
 *      Author: steffen
 */

#include "../../../include/marathon/chains/matching/matching_chain_JS89.h"

namespace marathon {

namespace chain {

namespace matching {

Broder86::Broder86(std::string line) :
		g(SparseBipartiteGraph(line)) {

}

Broder86::~Broder86() {
}

bool Broder86::computeArbitraryState(BipartiteMatching& s) const {

	int n = g.getNumberOfNodes();

	std::vector<int> mate;
	g.cardmax_matching(mate);

	// Count Number of Matching Edges
	int k = 0;
	for (int i = 0; i < mate.size(); i++)
		if (mate[i] != n)
			k++;

	if (k != g.getNumberOfNodes()) {
		std::cerr << "Graph doesn't have a perfect matching!" << std::endl;
		return false;
	}

	int unmatched[2] = { n, n };
	s = BipartiteMatching(g.getNumberOfNodes(), k / 2, unmatched, &mate[0]);

	return true;
}

void Broder86::computeNeighbours(const BipartiteMatching& s,
		boost::unordered_map<BipartiteMatching, Rational>& neighbors) const {

	// Variables
	const int n = g.getNumberOfNodes();
	const int k = s.k;
	int u, v;
	BipartiteMatching s2;
	edgelist edges;

	neighbors.clear();

	// Implementierung der Transitionen nach Vorlage JS89

	// Gehe über jede Kante
	g.getEdges(edges);

	//std::cout << "compute neighbors for s=" << s << std::endl;

	// Für jede Wahl einer Kante e=(u,v) wird Transition konstruiert
	for (typename edgelist::iterator it = edges.begin(); it != edges.end();
			++it) {

		u = it->first;
		v = it->second;

		// create new state as copy of old one
		s2 = BipartiteMatching(s);

		// Transition 1
		if (2 * k == n && s.mates[u] == v) {
			// remove edge (u,v)
			s2.mates[u] = n;
			s2.mates[v] = n;
			s2.k = k - 1;
			s2.unmatched[0] = u;
			s2.unmatched[1] = v;
		} else if (2 * k + 2 == n) {

			// Transition 2
			if (s.mates[u] == n && s.mates[v] == n) {
				// add edge (u,v)
				s2.mates[u] = v;
				s2.mates[v] = u;
				s2.k = k + 1;
				s2.unmatched[0] = n;
				s2.unmatched[1] = n;
			}
			// Transition 3a
			else if (s.mates[u] != n && s.mates[v] == n) {
				// remove edge (u, mate[u])
				// add edge (u,v)
				int w = s.mates[u];
				s2.mates[w] = n;
				s2.mates[u] = v;
				s2.mates[v] = u;
				// w=mate[u] (>= n/2) becomes unmatched node in bipartition group 1
				s2.unmatched[1] = w;
			}
			// Transition 3b
			else if (s.mates[u] == n && s.mates[v] != n) {
				// remove edge (v, mate[v])
				// add edge (u,v)
				int w = s.mates[v];
				s2.mates[w] = n;
				s2.mates[u] = v;
				s2.mates[v] = u;
				// w=mate[v] (< n/2) becomes unmatched node in bipartition group 0
				s2.unmatched[0] = w;
			} else {
				// stay in s
			}
		}
		// sonst
		else {
			// Verbleibe im aktuellen Zustand
		}

		neighbors[s2] += Rational(1, edges.size());

		//std::cout << std::endl;
	}
}

void Broder86::canonicalPath(int s, int t, std::list<int>& path) const {

	path.clear();

	// TODO: use indices instead of states

	// Implementation after JS89
	std::list<int> init_seqment;
	std::list<int> main_seqment;
	std::list<int> final_seqment;

	const int n = g.getNumberOfNodes();
	int i, a, sw;
	BipartiteMatching u, v;
	std::queue<BipartiteMatching> q;
	boost::unordered_map<BipartiteMatching, Rational> neighbors;
	boost::unordered_map<BipartiteMatching, BipartiteMatching> prev;
	BipartiteMatching null;

	if (s == t)
		return;

	BipartiteMatching start_states[2] = { states[s], states[t] };

	std::list<int>* segments[2] = { &init_seqment, &final_seqment };

	// Initial Segment and Final Segment by BFS
	for (i = 0; i < 2; i++) {

		// clear queue and prevs
		while (!q.empty())
			q.pop();
		prev.clear();

		// start BFS with s respective t
		q.push(start_states[i]);
		prev[start_states[i]] = null;
		while (!q.empty()) {
			u = q.front();
			q.pop();

			if (2 * u.k == n) {	// Path to Perfect Matching
				// Reconstruct Path and return
				do {
					segments[i]->push_front(indices.find(u)->second);
					u = prev[u];
				} while (!(u == null));
				break;
			} else {
				// for all neighbors of u
				computeNeighbours(u, neighbors);
				for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
					v = it->first;

					if (prev.find(v) == prev.end()) {	// not visited yet
						prev[v] = u;
						q.push(v);
					}
				}
			}
		}
	}

	// clean up

	final_seqment.reverse();

	// Main Segment
	BipartiteMatching I = states[init_seqment.back()];
	BipartiteMatching F = states[final_seqment.front()];
	BipartiteMatching s2 = I;

	BipartiteMatching x[2] = { I, F };

	std::vector<int> cycle;
	std::vector<int>::iterator it2;

	bool unrolled[n];
	memset(unrolled, 0, n * sizeof(bool));

	// unroll symmetric difference of I and F Cycle for Cycle
	for (i = 0; i < n; i++) {

		// if some node lies on a cycle which isn't unrolled yet
		if (I.mates[i] != F.mates[i] && !unrolled[i]) {

			// Collect Cycle Nodes and store in list
			cycle.clear();
			a = i;

			sw = 0;	// switch variable
			do {
				cycle.push_back(a);
				unrolled[a] = true;
				a = x[sw].mates[a];
				sw = 1 - sw;
			} while (a != i);

			// unroll cycle

			// first step: remove (u0, v0)

			int u0 = cycle[0];
			int v0 = cycle[1];
			s2.removeEdge(u0, v0);
			main_seqment.push_back(indices.find(s2)->second);

			// replace each edge (u_j, v_j) by (u_j, v_j-1)
			for (int j = 1; j < cycle.size() / 2; j++) {
				int u_j = cycle[2 * j];
				int v_j = cycle[2 * j + 1];
				int v_jj = cycle[2 * j - 1];
				s2.mates[u_j] = v_jj;
				s2.mates[v_jj] = u_j;
				s2.mates[v_j] = n;
				s2.unmatched[1] = v_j;
				main_seqment.push_back(indices.find(s2)->second);
			}

			// last step: add (u0, v_last)
			s2.addEdge(u0, cycle.back());
			main_seqment.push_back(indices.find(s2)->second);
		}
	}

	/*if (s == 0 && t == 12) {
	 std::cout << "init" << std::endl;
	 for (std::list<int>::iterator it = init_seqment.begin();
	 it != init_seqment.end(); ++it) {
	 std::cout << *it << " " << states[*it] << std::endl;
	 }

	 std::cout << "main" << std::endl;
	 for (std::list<int>::iterator it = main_seqment.begin();
	 it != main_seqment.end(); ++it) {
	 std::cout << *it << " " << states[*it] << std::endl;
	 }

	 std::cout << "final" << std::endl;
	 for (std::list<int>::iterator it = final_seqment.begin();
	 it != final_seqment.end(); ++it) {
	 std::cout << *it << " " << states[*it] << std::endl;
	 }
	 }*/

	path.insert(path.end(), init_seqment.begin(), init_seqment.end());
	path.insert(path.end(), main_seqment.begin(), main_seqment.end());
	path.insert(path.end(), ++final_seqment.begin(), final_seqment.end());
}

}

}

}
