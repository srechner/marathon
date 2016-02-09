/*
 * JS89Chain.cpp
 *
 *  Created on: Feb 11, 2015
 *      Author: steffen
 */

#include "../../../../include/marathon/chains/matching/matching_chain_JS89.h"

namespace marathon {

namespace chain {

namespace matching {

Broder86::Broder86(const std::string& inst) {
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

bool Broder86::computeArbitraryState(BipartiteMatching& s) {

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
		return false;
	}

	int unmatched[2] = { n, n };
	s = BipartiteMatching(g->getNumberOfNodes(), k / 2, unmatched, &mate[0]);

	return true;
}

void Broder86::computeNeighbours(const BipartiteMatching& s,
		std::unordered_map<BipartiteMatching, rational>& neighbors) const {

	// Variables
	const int n = g->getNumberOfNodes();
	const int k = s.k;
	int u, v;
	BipartiteMatching s2;
	edgelist edges;

	neighbors.clear();

	// Implementierung der Transitionen nach Vorlage JS89

	// Gehe über jede Kante
	g->getEdges(edges);

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

		neighbors[s2] += rational(1, edges.size());

		//std::cout << std::endl;
	}
}

void Broder86::constructPath(const StateGraph* sg,
		int s, int t, std::list<int>& path) const {

	path.clear();

	// reinterpret state graph object
	const _StateGraph<BipartiteMatching>* _sg = (const _StateGraph<
			BipartiteMatching>*) sg;

	if (s == t)
		return;

	// Implementation after JS89
	std::list<int> init_seqment;
	std::list<int> main_seqment;
	std::list<int> final_seqment;

	const int n = _sg->getState(0).n;
	const int omega = sg->getNumStates();

	// working variables
	int i, c, a, sw;
	BipartiteMatching u, v;
	std::queue<int> q;

	// arrays to store the labels for BFS
	int* prev = new int[omega];
	bool* visited = new bool[omega];

	//std::unordered_map<BipartiteMatching, BipartiteMatching> prev;
	//BipartiteMatching null;

	int start_states[2] = { s, t };

	std::list<int>* segments[2] = { &init_seqment, &final_seqment };

	// Initial Segment and Final Segment by BFS
	for (i = 0; i < 2; i++) {

		// clear queue and prevs
		while (!q.empty())
			q.pop();

		// start BFS with s respective t
		memset(visited, 0, omega*sizeof(bool));
		q.push(start_states[i]);
		visited[start_states[i]] = 1;
		prev[start_states[i]] = -1;
		while (!q.empty()) {
			c = q.front();
			q.pop();
			u = _sg->getState(c);

			if (2 * u.k == n) {	// Path to Perfect Matching
				// Reconstruct Path and return
				do {
					segments[i]->push_front(c);
					c = prev[c];
				} while (c != -1);
				break;
			} else {
				// for all neighbors of u
				for(Transition* t : sg->getOutArcs(c)) {
					// if not visited yet
					if(!visited[t->v]) {
						q.push(t->v);
						prev[t->v] = c;
						visited[t->v] = 1;
					}
				}
				/*computeNeighbours(u, neighbors);
				for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
					v = it->first;

					if (prev.find(v) == prev.end()) {	// not visited yet
						prev[v] = u;
						q.push(v);
					}
				}*/
			}
		}
	}

	final_seqment.reverse();

	// Main Segment
	BipartiteMatching I = _sg->getState(init_seqment.back());
	BipartiteMatching F = _sg->getState(final_seqment.front());
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
			main_seqment.push_back(_sg->findState(s2));

			// replace each edge (u_j, v_j) by (u_j, v_j-1)
			for (int j = 1; j < cycle.size() / 2; j++) {
				int u_j = cycle[2 * j];
				int v_j = cycle[2 * j + 1];
				int v_jj = cycle[2 * j - 1];
				s2.mates[u_j] = v_jj;
				s2.mates[v_jj] = u_j;
				s2.mates[v_j] = n;
				s2.unmatched[1] = v_j;
				main_seqment.push_back(_sg->findState(s2));
			}

			// last step: add (u0, v_last)
			s2.addEdge(u0, cycle.back());
			main_seqment.push_back(_sg->findState(s2));
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

	delete[] prev;
	delete[] visited;
}

}

}

}
