/*
 * Created on: Jan 08, 2018
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2018, Steffen Rechner
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

#ifndef MARATHON_MATCHING_COUNT_H
#define MARATHON_MATCHING_COUNT_H

// marathon includes
#include <marathon/state_graph.h>
#include "sparse_bipartite_graph.h"
#include "marathon/count.h"
#include "bipartite_matching.h"
#include "marathon/random_generator.h"

namespace marathon {
    namespace matching {

        /**
         * Class for counting the perfect and near-perfect matchings in an bipartite graph.
         */
        class Counter : public marathon::Counter {

        protected:

            const SparseBipartiteGraph _g;              // bipartite graph
            const size_t _n;                            // number of nodes in each vertex set

            marathon::Integer _num_perfect;             // number of perfect matchings
            marathon::Integer _num_near_perfect;        // number of near-perfect matchings
            std::vector<marathon::Integer> _num_near_perfect_uv;    // number of near-perfect matchings with (u,v) being unmatched
            std::unordered_map<SparseBipartiteGraph, marathon::Integer> _tmp;   // store temporal results

            /**
             * auxiliary function
             * @param match
             * @param value
             */
            void store_in_table(const BipartiteMatching &match, const marathon::Integer &value) {

                // create bipartite graph in which all matched nodes are removed
                boost::dynamic_bitset<> bits;
                for (int u = 0; u < _n; u++) {
                    if (!match.isMatched(u)) {
                        for (int v = _n; v < 2 * _n; v++) {
                            if (!match.isMatched(v)) {
                                //std::cout << "u=" << u << " v=" << v << " => " << _g.hasEdge(u,v) << std::endl;
                                bool b = _g.hasEdge(u, v);
                                bits.push_back(b);
                            }
                        }
                    }
                }

                SparseBipartiteGraph b(bits);
                // todo: b.canonize();

                // store entry with key b
                _tmp[b] = value;
            }

            /**
             * auxiliary function
             * @param match
             * @param value
             */
            bool load_from_table(const BipartiteMatching &match, marathon::Integer &key) const {

                // create bipartite graph in which all matched nodes are removed
                boost::dynamic_bitset<> bits;
                for (int u = 0; u < _n; u++) {
                    if (!match.isMatched(u)) {
                        for (int v = _n; v < 2 * _n; v++) {
                            if (!match.isMatched(v)) {
                                //std::cout << "u=" << u << " v=" << v << " => " << _g.hasEdge(u,v) << std::endl;
                                bool b = _g.hasEdge(u, v);
                                bits.push_back(b);
                            }
                        }
                    }
                }

                SparseBipartiteGraph b(bits);
                // todo: b.canonize();

                // try to load precomputed value
                auto it = _tmp.find(b);
                if (it != _tmp.end()) {
                    key = it->second;
                    return true;
                }

                return false;
            }


            /**
             * auxiliary function
             * @param u
             * @param v
             * @return
             */
            int translate(const int u, const int v) const {
                return u * _n + (v - _n);
            }

            /**
             * Auxiliary function used for backtracking.
             * @param match
             * @param u
             * @return
             */
            marathon::Integer count_recursive(
                    BipartiteMatching &match,
                    size_t u
            ) {

                // if each vertex of vertex set U has been considered
                if (u == _n) {
                    return Integer(1);
                }

                // if node u is already matched (are marked as permanently unmatched)
                if (match.isMatched(u)) {
                    return count_recursive(match, u + 1);
                }

                // simulate all choices to add an edge (u,v) to the matching
                marathon::Integer res(0);

                // iterate over all adjacent nodes
                std::vector<size_t> neighbors = _g.getNeighbors(u);
                for (size_t v : neighbors) {

                    if (!match.isMatched(v)) {

                        // add edge (u,v) to matching
                        match.addEdge(u, v);

                        // try to load cached result
                        marathon::Integer x;
                        bool found = load_from_table(match, x);
                        if (!found) {
                            // recursively count the number of solutions resulting from this choice
                            //std::cout << "start computing" << match << ": " << x << std::endl;
                            x = count_recursive(match, u + 1);
                        } else {
                            //std::cout << "found in table: " << match << ": " << x << std::endl;
                        }

                        res += x;

                        // undo changes
                        match.removeEdge(u, v);
                    }

                }

                // store result
                //std::cout << "store in table: " << match << ": " << res << std::endl;
                store_in_table(match, res);

                return res;
            }

        public:

            /**
             * Create a Counter for the perfect and near-perfect matchings in a bipartite graph.
             * @param g Bipartite graph.
             */
            Counter(SparseBipartiteGraph g) : _g(std::move(g)), _n((int) (_g.getNumberOfNodes() / 2)) {

                // auxiliary variable used for backtracking
                BipartiteMatching match(2 * _n); // empty matching on 2*_n nodes

                // count perfect matchings
                _num_perfect = count_recursive(match, 0);

                _num_near_perfect_uv.resize(_n * _n);

                // for each pair of unmatched nodes
                for (size_t u = 0; u < _n; u++) {
                    for (size_t v = _n; v < 2 * _n; v++) {

                        // mark (u,v) as permanently unmatched
                        match.mates[u] = SIZE_MAX - 1;
                        match.mates[v] = SIZE_MAX - 1;

                        // try to load cached result
                        marathon::Integer x;
                        bool found = load_from_table(match, x);
                        if (!found) {
                            // start counting
                            //std::cout << "start computing" << match << ": " << x << std::endl;
                            x = count_recursive(match, 0);
                        } else {
                            //std::cout << "found in table: " << match << ": " << x << std::endl;
                        }

                        _num_near_perfect_uv[translate(u, v)] += x;
                        _num_near_perfect += x;

                        // undo marking
                        match.mates[u] = SIZE_MAX;
                        match.mates[v] = SIZE_MAX;
                    }
                }

                /*std::cout << "counts:" << std::endl;
                std::cout << "perfect         : " << _num_perfect << std::endl;
                std::cout << "near perfect    : " << _num_near_perfect << std::endl;
                for (int u = 0; u < _n; u++) {
                    for (int v = _n; v < 2 * _n; v++) {
                        printf(" %2i %2i: ", u, v);
                        std::cout << _num_near_perfect_uv[translate(u, v)] << std::endl;
                    }
                }*/
            }

            /**
             * Count the number of perfect and near-perfect matchings in the bipartite graph.
             */
            marathon::Integer count() override {
                return countPerfect() + countNearPerfect();
            }

            /**
             * Count the number of perfect matchings in the bipartite graph.
             */
            marathon::Integer countPerfect() const {
                return _num_perfect;
            }

            /**
             * Count the number of near-perfect matchings in the bipartite graph.
             */
            marathon::Integer countNearPerfect() const {
                return _num_near_perfect;
            }

            /**
             * Count the number of near-perfect matchings where u and v are unmatched vertices.
             */
            marathon::Integer countNearPerfect(const int u, const int v) const {
                return _num_near_perfect_uv[translate(u, v)];
            }
        };
    }
}

#endif //MARATHON_MATCHING_COUNT_H
