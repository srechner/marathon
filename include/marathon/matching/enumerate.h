/*
 * Created on: Jan 09, 2018
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


#ifndef MARATHON_MATCHING_ENUMERATE_H
#define MARATHON_MATCHING_ENUMERATE_H

#include <cstdlib>
#include <vector>
#include <numeric>
#include <algorithm>

// marathon includes
#include "marathon/enumerate.h"
#include "marathon/matching/count.h"

namespace marathon {
    namespace matching {

        /**
         * Class for enumerating all perfect and near-perfect matchings in a bipartite graph.
         */
        class Enumerator : public marathon::Enumerator {

        protected:

            const SparseBipartiteGraph _g;              // bipartite graph
            const int _n;                               // number of nodes in each vertex set

            void enumerate_recursive(
                    BipartiteMatching &match,
                    const int u,
                    const int unmatched1,
                    const int unmatched2,
                    const std::function<void(const State *)>& f
            ) const {

                // if each vertex of vertex set U has been considered
                if (u == _n) {

                    //
                    if(unmatched1 != -1) {
                        match.mates[unmatched1] = -1;
                        match.mates[unmatched2] = -1;
                        f(&match);
                        match.mates[unmatched1] = -2;
                        match.mates[unmatched2] = -2;
                    }
                    else {
                        f(&match);
                    }
                    return;
                }

                // if node u is already matched (are marked as permanently unmatched)
                if (match.isMatched(u)) {
                    return enumerate_recursive(match, u + 1, unmatched1, unmatched2, f);
                }

                // simulate all choices to add an edge (u,v) to the matching

                // iterate over all adjacent nodes
                std::vector<int> neighbors;
                _g.getNeighbors(u, neighbors);
                for (int v : neighbors) {

                    if (!match.isMatched(v)) {
                        // add edge (u,v) to matching
                        match.mates[u] = v;
                        match.mates[v] = u;
                        match.k++;

                        enumerate_recursive(match, u+1, unmatched1, unmatched2, f);

                        // undo changes
                        match.mates[u] = -1;
                        match.mates[v] = -1;
                        match.k--;
                    }
                }
            }

        public:

            /**
             * Create an Enumerator object.
             * @param
             */
            Enumerator(const SparseBipartiteGraph &g) :
                    _g(g),
                    _n((int) (g.getNumberOfNodes() / 2)) {

            }


            virtual ~Enumerator() override {

            }

            /**
             * Enumerate all perfect and near-perfect matchings an an bipartite graph.
             * @param f Function that is evaluated for each binary matrix.
             */
            void enumerate(const std::function<void(const State *)> f) override {

                BipartiteMatching match(2 * _n); // empty matching on 2*_n nodes

                // enumerate perfect matchings
                enumerate_recursive(match, 0, -1, -1, f);

                // enumerate near-perfect matchings
                for(int u=0; u<_n; u++){
                    for(int v=_n; v<2*_n; v++) {

                        // mark (u,v) as permanently unmatched
                        match.mates[u] = -2;
                        match.mates[v] = -2;

                        enumerate_recursive(match, 0, u, v, f);

                        // undo marking
                        match.mates[u] = -1;
                        match.mates[v] = -1;
                    }
                }

            }
        };
    }
}

#endif //MARATHON_MATCHING_ENUMERATE_H
