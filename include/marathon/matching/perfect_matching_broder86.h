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
			Broder86(SparseBipartiteGraph g) : MarkovChain(std::move(g)) {

			}

			/**
			 * Create a Markov chain object with specified initial state.
			 * @param g Bipartite Graph.
			 * @param m Perfect or near-perfect matching in g.
			 */
			Broder86(SparseBipartiteGraph g, BipartiteMatching m) :
					MarkovChain(std::move(g), std::move(m)) {

			}

			/**
             * Create a copy of this MarkovChain.
             * @return
             */
			virtual std::unique_ptr<marathon::MarkovChain> copy() const {
				return std::make_unique<Broder86>(*this);
			}


			/**
			 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
			 * For each pair (x,p) call the function f.
			 * @param x A state.
			 * @param process Callback function to call for each adjacent state.
			 */
			virtual void adjacentStates(
					const State &x,
					const std::function<void(const State &, const marathon::Rational &)> &f
			) const override {

				// Variables
				const BipartiteMatching &s = static_cast<const BipartiteMatching &>(x);

				// Jerrum and Sinclair, 1989. Approximating the Permanent.

				const Rational p((long) 1, (long) _edges.size());

				// for each edge e=(u,v)
				for (auto it = _edges.begin(); it != _edges.end(); ++it) {

					size_t u = it->first;
					size_t v = it->second;

					// create new state as copy of old one
					BipartiteMatching s2(s);

					// if matching is perfect and (u,v) can be removed from it
					if (s.is_perfect() && s.getMate(u) == v) {

						// remove edge (u,v)
						s2.removeEdge(u,v);

					} else if (s.is_near_perfect()) {	// matching is near-perfect

						// (u,v) can be added
						if (!s.isMatched(u) && !s.isMatched(v)) {

							// add edge (u,v)
							s2.addEdge(u,v);
						}
						else if (s.isMatched(u) && !s.isMatched(v)) {
							// remove edge (u, mate[u])
							// add edge (u,v)
							size_t w = s.getMate(u);
							s2.mates[w] = SIZE_MAX;
							s2.mates[u] = v;
							s2.mates[v] = u;
							s2.unmatched2 = w;			// w=mate[u] (>= n/2) becomes unmatched node
						}
						else if (!s.isMatched(u) && s.isMatched(v)) {
							// remove edge (v, mate[v])
							// add edge (u,v)
							size_t w = s.getMate(v);
							s2.mates[w] = SIZE_MAX;
							s2.mates[u] = v;
							s2.mates[v] = u;
							s2.unmatched1 = w;	// w=mate[v] (< n/2) becomes unmatched node
						} else {
							// stay in s
						}
					}
					else {
						// stay at current state
					}

					f(s2, p);
				}
			}

			virtual void step() override {

				// select an edge (u,v) uniformly at random
				const size_t e = rg.nextInt(_edges.size());
				const size_t u = _edges[e].first;
				const size_t v = _edges[e].second;

				// if m is perfect and (u,v) can be removed
				if (currentState.is_perfect() && currentState.mates[u] == v) {
					// remove edge (u,v)
					currentState.mates[u] = SIZE_MAX;
					currentState.mates[v] = SIZE_MAX;
					currentState._edges--;
					currentState.unmatched1 = u;
					currentState.unmatched2 = v;
				} else if (currentState.is_near_perfect()) {	// matching is near-perfect

					// if (u,v) can be added to m
					if (!currentState.isMatched(u) && !currentState.isMatched(v)) {
						// add edge (u,v) to m
						currentState.mates[u] = v;
						currentState.mates[v] = u;
						currentState._edges++;
						currentState.unmatched1 = SIZE_MAX;
						currentState.unmatched2 = SIZE_MAX;
					}
					else if (currentState.isMatched(u) && !currentState.isMatched(v)) {
						// remove edge (u, mate[u]) from m
						// add edge (u,v) to m
						size_t w = currentState.mates[u];
						currentState.mates[w] = SIZE_MAX;
						currentState.mates[u] = v;
						currentState.mates[v] = u;
						currentState.unmatched2 = w;	// w=mate[u] (>= n/2) becomes unmatched
					}
					else if (!currentState.isMatched(u) && currentState.isMatched(v)) {
						// remove edge (v, mate[v]) from m
						// add edge (u,v) to m
						size_t w = currentState.mates[v];
						currentState.mates[w] = SIZE_MAX;
						currentState.mates[u] = v;
						currentState.mates[v] = u;
						currentState.unmatched1 = w;	// w=mate[v] (< n/2) becomes unmatched
					} else {
						// stay at current state
					}
				}
			}
		};
	}
}
#endif /* JS89_H_ */

