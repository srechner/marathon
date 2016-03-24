/*
 * JS89CanPath.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#include "../../../../include/marathon/chain/matching/JS89CanPath.h"
#include "../../../../include/marathon/chain/matching/Broder86.h"

namespace marathon {
namespace chain {
namespace matching {

void JS89Path::construct(const StateGraph* sg, const int s, const int t,
		std::list<int>& path) const {

	const Broder86* mc = (const Broder86*) sg->getMarkovChain();
	const SparseBipartiteGraph* g = mc->g;

	path.clear();

	if (s == t)
		return;

	// Implementation after JS89
	std::list<int> init_seqment;
	std::list<int> main_seqment;
	std::list<int> final_seqment;

	const int omega = sg->getNumStates();

	// working variables
	int i, c, a, sw;
	const int n = g->getNumberOfNodes();
	const BipartiteMatching* u, *v;
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
		memset(visited, 0, omega * sizeof(bool));
		q.push(start_states[i]);
		visited[start_states[i]] = 1;
		prev[start_states[i]] = -1;
		while (!q.empty()) {
			c = q.front();
			q.pop();
			u =
					(const ::marathon::chain::matching::BipartiteMatching*) sg->getState(
							c);

			if (2 * u->k == n) {	// Path to Perfect Matching
				// Reconstruct Path and return
				do {
					segments[i]->push_front(c);
					c = prev[c];
				} while (c != -1);
				break;
			} else {
				// for all neighbors of u
				for (Transition* t : sg->getOutArcs(c)) {
					// if not visited yet
					if (!visited[t->v]) {
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
	const BipartiteMatching *I = (const BipartiteMatching*) sg->getState(
			init_seqment.back());
	const BipartiteMatching *F = (const BipartiteMatching*) sg->getState(
			final_seqment.front());
	BipartiteMatching s2(*I);

	const BipartiteMatching* x[2] = { I, F };

	std::vector<int> cycle;
	std::vector<int>::iterator it2;

	bool unrolled[n];
	memset(unrolled, 0, n * sizeof(bool));

	// unroll symmetric difference of I and F Cycle for Cycle
	for (i = 0; i < n; i++) {

		// if some node lies on a cycle which isn't unrolled yet
		if (I->mates[i] != F->mates[i] && !unrolled[i]) {

			// Collect Cycle Nodes and store in list
			cycle.clear();
			a = i;

			sw = 0;	// switch variable
			do {
				cycle.push_back(a);
				unrolled[a] = true;
				a = x[sw]->mates[a];
				sw = 1 - sw;
			} while (a != i);

			// unroll cycle

			// first step: remove (u0, v0)

			int u0 = cycle[0];
			int v0 = cycle[1];
			s2.removeEdge(u0, v0);
			main_seqment.push_back(sg->indexOf(&s2));

			// replace each edge (u_j, v_j) by (u_j, v_j-1)
			for (int j = 1; j < cycle.size() / 2; j++) {
				int u_j = cycle[2 * j];
				int v_j = cycle[2 * j + 1];
				int v_jj = cycle[2 * j - 1];
				s2.mates[u_j] = v_jj;
				s2.mates[v_jj] = u_j;
				s2.mates[v_j] = n;
				s2.unmatched[1] = v_j;
				main_seqment.push_back(sg->indexOf(&s2));
			}

			// last step: add (u0, v_last)
			s2.addEdge(u0, cycle.back());
			main_seqment.push_back(sg->indexOf(&s2));
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
