/*
 * Created on: Nov 1, 2016
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

#ifndef MARATHON_MATCHINGCHAIN_H
#define MARATHON_MATCHINGCHAIN_H

#include "marathon/markov_chain.h"
#include "marathon/matching/realize.h"

namespace marathon {
    namespace matching {

        /**
         * Markov chains that sample bipartite matchings.
         */
        class MarkovChain : public ::marathon::MarkovChain {

            friend class JS89Path;

        protected:

            SparseBipartiteGraph _g;
            edgelist _edges;

        public:

            /**
			 * Create a Markov chain object with specified initial state.
			 * @param g Bipartite graph.
			 * @param m Perfect or near-perfect matching in g.
			 */
            MarkovChain(const SparseBipartiteGraph& g, const BipartiteMatching& m)
            : _g(g) {
                _g.getEdges(_edges);
                currentState = m.copy();
            }

            /**
             * Create a Markov chain object.
             * @param g Bipartite graph.
             */
            MarkovChain(const SparseBipartiteGraph& g)
                    : _g(g) {
                _g.getEdges(_edges);
                currentState = cardmax_matching(g);
            }


            /**
             * Instances of this chain have the form "110101011".
             * Such a 0-1-String is interpreted as the bi-adjacency matrix M = (m_ij) of a bipartite graph G=(V,E).
             * The bi-adjacency matrix M is defined as m_ij = 1, if (i,j) is in E, or 0, otherwise. The rows of M
             * are concatenated to a single string.
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
            MarkovChain(const std::string &inst) :
                    _g(inst) {
                _g.getEdges(_edges);
                currentState = cardmax_matching(_g);
            }


            /**
             * Create a Markov chain as a copy of another one.
             * @param mc Another Markov chain object.
             */
            MarkovChain(const MarkovChain& mc) :
                    _g(mc._g), _edges(mc._edges) {
                // copy the current state of mc
                currentState = mc.currentState->copy();
            }

            virtual ~MarkovChain() {

            }

            /**
             * Return the current state of the Markov chain.
             * @return
             */
            virtual const BipartiteMatching *getCurrentState() const override {
                return (const BipartiteMatching*) currentState;
            }

            /**
             * Create an independent copy of the Markov chain.
             * @return Copy of this Markov chain.
             */
            virtual marathon::matching::MarkovChain* copy() const override = 0;

        };
    }
}

#endif //MARATHON_MATCHINGCHAIN_H
