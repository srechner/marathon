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

            SparseBipartiteGraph _g;            // bipartite graph
            edgelist _edges;                    // list of edges in _g
            BipartiteMatching currentState;     // current state of the Markov chain

        public:

            /**
			 * Create a Markov chain object with specified initial state.
			 * @param g Bipartite graph.
			 * @param m Perfect or near-perfect matching in g.
			 */
            MarkovChain(SparseBipartiteGraph g, BipartiteMatching m)
                    : _g(std::move(g)), _edges(_g.getEdges()), currentState(std::move(m)) {

            }

            /**
             * Create a Markov chain object.
             * @param g Bipartite graph.
             */
            MarkovChain(SparseBipartiteGraph g)
                    : _g(std::move(g)), _edges(_g.getEdges()), currentState(cardmax_matching(_g)) {

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
                    _g(inst), _edges(_g.getEdges()), currentState(cardmax_matching(_g)) {
                
            }


            /**
             * Return the current state of the Markov chain.
             * @return
             */
            virtual const BipartiteMatching &getCurrentState() const override {
                return currentState;
            }

            /**
             * Set the current state.
             * @param s Bipartite matching.
             */
            virtual void setCurrentState(const State &s) override {

                // try to cast s to the correct type
                auto m = dynamic_cast<const BipartiteMatching *>(&s);
                if (m == nullptr)
                    throw std::runtime_error("Error! State is not a binary matrix!");

                currentState = *m;
            }

        };
    }
}

#endif //MARATHON_MATCHINGCHAIN_H
