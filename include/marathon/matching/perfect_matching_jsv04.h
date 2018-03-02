/*
 * JSV04.h
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

#ifndef JSV04_CHAIN_H_
#define JSV04_CHAIN_H_

#include <queue>

#include "count.h"
#include "markov_chain.h"

namespace marathon {
    namespace matching {

        /**
         * Markov chain defined by
         *
         * M. Jerrum, A. Sinclair, and E. Vigoda.
         * A polynomial-time approximation algorithm for the permanent of a matrix with nonnegative entries.
         * Journal of the ACM 51 (2004), 671–697. doi: 10.1145/1008731.1008738.
         */
        class JSVChain : public MarkovChain {

        protected:

            Counter *_cnt;      // counter for perfect and near-perfect matchings

            /**
             * Define a hash function for integer pairs.
             */
            struct hash {
                size_t operator()(const std::pair<int, int>& p) const {
                    return p.first ^ p.second;
                }
            };

            // store weight of each near perfect matchings for fast access
            std::unordered_map<std::pair<int, int>, double, hash> _weight;

            /**
             * Initialize the non-constant class members.
             */
            void init_weight() {
                const int n = (int) _g.getNumberOfNodes()/2;
                const Rational M = _cnt->countPerfect();
                for(int u=0; u<n; u++) {
                    for(int v = n; v<2*n; v++) {
                        const Rational Nuv = _cnt->countNearPerfect(u,v);
                        if(Nuv > 0) {
                            auto uv = std::make_pair(u, v);
                            _weight[uv] = (M / Nuv).convert_to<double>();
                        }
                    }
                }
            }

        public:

            /**
			 * Create a Markov chain object.
			 * @param inst instance representation (see base class)
			 */
            explicit JSVChain(const std::string &inst)
                    : MarkovChain(inst) {
                _cnt = new Counter(_g);
                init_weight();
            }

            /**
             * Create a Markov chain object.
             * @param g Bipartite graph.
             */
            explicit JSVChain(const SparseBipartiteGraph& g)
                    : MarkovChain(g) {
                _cnt = new Counter(_g);
                init_weight();
            }

            /**
             * Create a Markov chain object with specified initial state.
             * @param g Bipartite Graph.
             * @param m Perfect or near-perfect matching in g.
             */
            JSVChain(const SparseBipartiteGraph& g, const BipartiteMatching& m) :
                    MarkovChain(g, m) {
                _cnt = new Counter(_g);
                init_weight();
            }

            /**
             * Create a Markov chain object as a copy of another.
             * @param mc Another Markov chain object.
             */
            JSVChain(const JSVChain& mc) :
                   MarkovChain(mc), _weight(mc._weight) {

                _cnt = mc._cnt; // do not copy but re-use same counter object
            }

            /**
             * Create a copy of this MarkovChain.
             * @return
             */
            virtual ::marathon::MarkovChain *copy() const {
                return new ::marathon::matching::JSVChain(*this);
            }


            /**
             * Return the weight of state x.
             * @param x Perfect or near-perfect matching object.
             * @return w(x)
             */
            Rational getWeight(const State *x) const override {

                const BipartiteMatching *s = (const BipartiteMatching *) x;

                if (s->is_perfect()) {
                    return 1;
                } else if (s->is_near_perfect()) {
                    const int u = s->unmatched[0];
                    const int v = s->unmatched[1];
                    return Rational(_cnt->countPerfect(), _cnt->countNearPerfect(u, v));
                } else {
                    // cannot happen
                    throw std::runtime_error("Error! Illegal Matching!");
                }
            }

            /**
             * Enumerate the set of adjacent states t of state s and the transition probability p(s,t).
             * For each pair (t, p) evaluate the function f once.
             * @param s State.
             * @param f Function to be evaluated for each adjacent state.
             */
            void adjacentStates(
                    const State *s,
                    const std::function<void(const State *, const marathon::Rational &)> &f
            ) const override {

                // Variables
                auto x = (const BipartiteMatching *) s;
                const int n = (int) _g.getNumberOfNodes();

                const Rational p(2, n);
                const Rational q(1, n);

                Rational loop(0);

                // Transition rules presented by Jerrum, Sinclair, Vigoda 2004.

                // if x is perfect matching
                if (x->is_perfect()) {

                    // create copy of x
                    BipartiteMatching y(*x);

                    // Choose a matching edge e=(u,v)
                    for (int u = 0; 2 * u < n; u++) {

                        const int v = x->mates[u];

                        // remove (u,v)
                        y.removeEdge(u, v);

                        const Rational wx = getWeight(x);
                        const Rational wy = getWeight(&y);

                        // metr = p * min(1, w(y)/w(x))
                        Rational metr = p;
                        if (wy < wx) {
                            metr *= wy / wx;
                            loop += p - metr;
                        }

                        f(&y, metr);

                        // undo changes
                        y.addEdge(u, v);
                    }

                } else if (x->is_near_perfect()) {    // is a near-perfect matching

                    const int u = x->unmatched[0];
                    const int v = x->unmatched[1];

                    // choose node z
                    for (int z = 0; z < n; z++) {

                        // create copy
                        BipartiteMatching y(*x);

                        // Three Cases to rule them all

                        if ((z == u || z == v) && _g.hasEdge(u, v)) {
                            y.addEdge(u, v);
                        } else if (2 * z >= n && _g.hasEdge(u, z) && x->isMatched(z)) {
                            int a = x->mates[z];
                            y.removeEdge(a, z);
                            y.addEdge(u, z);
                            y.unmatched[0] = a;
                            y.unmatched[1] = v;
                        } else if (2 * z < n && _g.hasEdge(z, v) && x->isMatched(z)) {
                            int b = x->mates[z];
                            y.removeEdge(z, b);
                            y.addEdge(z, v);
                            y.unmatched[0] = u;
                            y.unmatched[1] = b;
                        } else {
                            // loop
                        }

                        const Rational wx = getWeight(x);
                        const Rational wy = getWeight(&y);

                        // metr = p * min(1, w(y)/w(x))
                        Rational metr = q;
                        if (wy < wx) {
                            metr *= wy / wx;
                            loop += q - metr;
                        }

                        f(&y, metr);
                    }
                }

                // process loop transition
                f(s, loop);
            }


            /**
             * Apply a single step of the Markov chain to the current state.
             */
            void step() override {

                auto x = (BipartiteMatching *) getCurrentState();
                const int n = (int) _g.getNumberOfNodes();

                // draw real number between zero and one uniformly at random
                double r = rg.nextDouble();

                // if x is perfect matching
                if (x->is_perfect()) {

                    // weight of state x
                    const double wx = 1.0;
                    //const double wx = getWeight(x).convert_to<double>();

                    // Choose matching edge e=(u,v) uniformly at random and remove it
                    const int u = rg.nextInt(n);
                    const int v = x->mates[u];
                    x->removeEdge(u, v);

                    // apply metropolis rule
                    auto uv = std::make_pair(x->unmatched[0], x->unmatched[1]);
                    const double wy = _weight[uv];
                    //const double wy = getWeight(x).convert_to<double>();
                    const double metr = wy;         // metr = min(1.0, wy / wx);
                    if (r >= metr) {
                        // reject transition by reversing the operation
                        x->addEdge(u, v);
                    }

                } else if (x->is_near_perfect()) {    // if x is a near-perfect matching

                    // determine unmatched nodes
                    const int u = x->unmatched[0];
                    const int v = x->unmatched[1];

                    // weight of state x
                    auto uv = std::make_pair(x->unmatched[0],x->unmatched[1]);
                    const double wx = _weight[uv];
                    //const double wx = getWeight(x).convert_to<double>();

                    // choose a node z uniformly at random
                    const int z = rg.nextInt(n);

                    // if z is an unmatched node
                    if ((z == u || z == v) && _g.hasEdge(u, v)) {
                        x->addEdge(u, v);

                        // apply metropolis rule
                        const double wy = 1.0;  // y is a perfect matching
                        const double metr = wy / wx; // metr = std::min(1.0, wy / wx);
                        if (r >= metr) {
                            // reject transition by reversing the operation
                            x->removeEdge(u, v);
                        }


                    } else if (2 * z >= n && _g.hasEdge(u, z) && x->isMatched(z)) {

                        // shift edge
                        const int a = x->mates[z];
                        x->mates[u] = z;
                        x->mates[z] = u;
                        x->mates[a] = -1;
                        x->unmatched[0] = a;

                        // apply metropolis rule
                        auto uv = std::make_pair(x->unmatched[0], x->unmatched[1]);
                        const double wy = _weight[uv];
                        const double metr = wy / wx;  // metr = std::min(1.0, wy / wx);
                        if (r >= metr) {
                            // reject transition by reversing the operation
                            x->mates[z] = a;
                            x->mates[a] = z;
                            x->mates[u] = -1;
                            x->unmatched[0] = u;
                        }

                    } else if (2 * z < n && _g.hasEdge(z, v) && x->isMatched(z)) {

                        // shift edge
                        const int b = x->mates[z];
                        x->mates[v] = z;
                        x->mates[z] = v;
                        x->mates[b] = -1;
                        x->unmatched[1] = b;

                        // apply metropolis rule
                        auto uv = std::make_pair(x->unmatched[0], x->unmatched[1]);
                        const double wy = _weight[uv];
                        const double metr = wy/wx; // metr = std::min(1.0, wy / wx);
                        if (r >= metr) {
                            // reject transition by reversing the operation
                            x->mates[z] = b;
                            x->mates[b] = z;
                            x->mates[v] = -1;
                            x->unmatched[1] = v;
                        }
                    }
                }
            }
        };
    }
}

#endif /* JS89_H_ */
