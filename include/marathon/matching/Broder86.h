/*
 * Broder86.h
 *
 * Created on: Nov 18, 2014
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

#ifndef JS89_CHAIN_H_
#define JS89_CHAIN_H_

#include <queue>

#include "marathon/MarkovChain.h"
#include "BipartiteMatching.h"
#include "SparseBipartiteGraph.h"

namespace marathon {
	namespace matching {

		class Broder86 : public MarkovChain {

			friend class JS89Path;

		protected:

			SparseBipartiteGraph *g = nullptr;

			/**
			 * Instances have the form "110101011".
			 * Such a 0-1-String is interpreted as a biadjacency matrix of a bipartite graph, flattened to a single line.
			 * Thus, the input string above corresponds to the biadjacency  matrix
			 *
			 *  1 1 0
			 *  1 0 1
			 *  0 1 1
			 *
			 *  which is the graph
			 *
			 *  u1  u2  u3
			 *  |\ / \ /|
			 *  | X   X |
			 *  |/ \ / \|
			 *  v1  v2  v3
			 * .
			 */
			void parseInstance(const std::string &inst) {
				if (g != nullptr)
					delete g;
				g = new SparseBipartiteGraph(inst);
			}

		public:

			Broder86(const std::string &inst) {
				// construct input object from input instance
				parseInstance(inst);
			}

			~Broder86() {
				if (g != nullptr)
					delete g;
			}

			virtual State *computeArbitraryState() const {

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

				int unmatched[2] = {n, n};
				State *s = new BipartiteMatching(g->getNumberOfNodes(), k / 2, unmatched,
				                                 &mate[0]);

				return s;
			}

			virtual void computeNeighbours(const State *x,
			                               std::vector<std::pair<State *, rational>> &neighbors) const {

				// Variables
				const BipartiteMatching *s = (const BipartiteMatching *) x;
				const int n = g->getNumberOfNodes();
				const int k = s->k;
				int u, v;
				edgelist edges;

				neighbors.clear();

				// Implementierung der Transitionen nach Vorlage JS89

				// Gehe über jede Kante
				g->getEdges(edges);

				const rational p((long) 1, (long) edges.size());

				//std::cout << "compute neighbors for s=" << *s << std::endl;

				// Für jede Wahl einer Kante e=(u,v) wird Transition konstruiert
				for (typename edgelist::iterator it = edges.begin(); it != edges.end();
				     ++it) {

					u = it->first;
					v = it->second;

					// create new state as copy of old one
					BipartiteMatching *s2 = new BipartiteMatching(*s);

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

		};
	}
}
#endif /* JS89_H_ */

