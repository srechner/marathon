/*
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

#include "markov_chain.h"

namespace marathon {
	namespace matching {

		/**
		 * Markov chain presented by
		 *
		 * A. Z. Broder.
		 * How hard is it to marry at random? (On the approximation of the permanent).
		 * Proceedings of the Eighteenth Annual ACM Symposium on Theory of Computing (1986), 50–58.
		 */
		class Broder86 : public MarkovChain {

		public:

			/**
			 * Create a Markov chain object.
			 * @param inst instance representation (see base class)
			 */
			Broder86(const std::string &inst)
					: MarkovChain(inst) {

			}

			/**
			 * Create a Markov chain object.
			 * @param g Bipartite graph.
			 */
			Broder86(const SparseBipartiteGraph& g) : MarkovChain(g) {

			}

			/**
			 * Create a Markov chain object with specified initial state.
			 * @param g Bipartite Graph.
			 * @param m Perfect or near-perfect matching in g.
			 */
			Broder86(const SparseBipartiteGraph& g, const BipartiteMatching& m) :
					MarkovChain(g, m) {

			}

			/**
			 * Create a Markov chain object as a copy of another one.
			 * @param mc Another Markov chain object.
			 */
			Broder86(const Broder86& mc) :
					MarkovChain(mc) {

			}

			/**
             * Create a copy of this MarkovChain.
             * @return
             */
			virtual ::marathon::MarkovChain *copy() const {
				return new ::marathon::matching::Broder86(*this);
			}


			/**
			 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
			 * For each pair (x,p) call the function f.
			 * @param x A state.
			 * @param process Callback function to call for each adjacent state.
			 */
			virtual void adjacentStates(
					const State *x,
					const std::function<void(const State *, const marathon::Rational &)> &f
			) const override {

				// Variables
				const BipartiteMatching *s = (const BipartiteMatching *) x;

				// Jerrum and Sinclair, 1989. Approximating the Permanent.

				const Rational p((long) 1, (long) _edges.size());

				// for each edge e=(u,v)
				for (auto it = _edges.begin(); it != _edges.end(); ++it) {

					int u = it->first;
					int v = it->second;

					// create new state as copy of old one
					BipartiteMatching s2(*s);

					// if matching is perfect and (u,v) can be removed from it
					if (s->is_perfect() && s->mates[u] == v) {
						// remove edge (u,v)
						s2.mates[u] = -1;
						s2.mates[v] = -1;
						s2.k--;
						s2.unmatched[0] = u;
						s2.unmatched[1] = v;
					} else if (s->is_near_perfect()) {	// matching is near-perfect

						// (u,v) can be added
						if (!s->isMatched(u) && !s->isMatched(v)) {
							// add edge (u,v)
							s2.mates[u] = v;
							s2.mates[v] = u;
							s2.k++;
							s2.unmatched[0] = -1;
							s2.unmatched[1] = -1;
						}
						else if (s->isMatched(u) && !s->isMatched(v)) {
							// remove edge (u, mate[u])
							// add edge (u,v)
							int w = s->mates[u];
							s2.mates[w] = -1;
							s2.mates[u] = v;
							s2.mates[v] = u;
							s2.unmatched[1] = w;	// w=mate[u] (>= n/2) becomes unmatched node
						}
						else if (!s->isMatched(u) && s->isMatched(v)) {
							// remove edge (v, mate[v])
							// add edge (u,v)
							int w = s->mates[v];
							s2.mates[w] = -1;
							s2.mates[u] = v;
							s2.mates[v] = u;
							s2.unmatched[0] = w;	// w=mate[v] (< n/2) becomes unmatched node
						} else {
							// stay in s
						}
					}
					else {
						// stay at current state
					}

					f(&s2, p);
				}
			}

			virtual void step() override {

				auto m = (BipartiteMatching*) getCurrentState();

				// select an edge (u,v) uniformly at random
				const int e = rg.nextInt(_edges.size());
				const int u = _edges[e].first;
				const int v = _edges[e].second;

				// if m is perfect and (u,v) can be removed
				if (m->is_perfect() && m->mates[u] == v) {
					// remove edge (u,v)
					m->mates[u] = -1;
					m->mates[v] = -1;
					m->k--;
					m->unmatched[0] = u;
					m->unmatched[1] = v;
				} else if (m->is_near_perfect()) {	// matching is near-perfect

					// if (u,v) can be added to m
					if (!m->isMatched(u) && !m->isMatched(v)) {
						// add edge (u,v) to m
						m->mates[u] = v;
						m->mates[v] = u;
						m->k++;
						m->unmatched[0] = -1;
						m->unmatched[1] = -1;
					}
					else if (m->isMatched(u) && !m->isMatched(v)) {
						// remove edge (u, mate[u]) from m
						// add edge (u,v) to m
						int w = m->mates[u];
						m->mates[w] = -1;
						m->mates[u] = v;
						m->mates[v] = u;
						m->unmatched[1] = w;	// w=mate[u] (>= n/2) becomes unmatched
					}
					else if (!m->isMatched(u) && m->isMatched(v)) {
						// remove edge (v, mate[v]) from m
						// add edge (u,v) to m
						int w = m->mates[v];
						m->mates[w] = -1;
						m->mates[u] = v;
						m->mates[v] = u;
						m->unmatched[0] = w;	// w=mate[v] (< n/2) becomes unmatched
					} else {
						// stay at current state
					}
				}
			}
		};
	}
}
#endif /* JS89_H_ */

