/*
 * Created on: Nov 20, 2014
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

#ifndef STATE_JS89_H_
#define STATE_JS89_H_

#include <cstring>
#include <cstdlib>
#include <ostream>

#include "marathon/state.h"
#include "sparse_bipartite_graph.h"

namespace marathon {
    namespace matching {

        /**
         * Class representing a matching in a bipartite graph.
         */
        class BipartiteMatching : public State {

            friend class Broder86;
            friend class JSVChain;
            friend class Counter;
            friend class Enumerator;
            friend class JS89Path;
            friend class RandomGeneratorExact;

        public:

            size_t n;                   // number of nodes in each vertex sets
            size_t _edges;              // number of edges
            size_t unmatched1;          // index of unmatched node in first vertex set (or SIZE_MAX if none)
            size_t unmatched2;          // index of unmatched node in second vertex set (or SIZE_MAX if none)
            std::vector<size_t> mates;  // matching partner of each node (or SIZE_MAX if none)

        public:

            /**
             * Create empty matching in a bipartite graph with n nodes.
             * @param n
             */
            BipartiteMatching(size_t n = 0) :
                    n(n), _edges(0) {
                mates.resize(n, SIZE_MAX);
            }

            /**
             * Create a matching of size n.
             * @param n Number of nodes.
             * @param k Number of matching edges.
             * @param unmatched1 Index of first unmatched node (or SIZE_MAX, if none).
             * @param unmatched2 Index of second unmatched node (or SIZE_MAX, if none).
             * @param matex Vector of matching partners for each node.
             */
            BipartiteMatching(
                    size_t n, size_t k, std::vector<size_t> mates, size_t unmatched1, size_t unmatched2) :
                    n(n), _edges(k), unmatched1(unmatched1), unmatched2(unmatched2), mates(std::move(mates)) {
            }


            /**
             * Add the edge (u,v) to the matching.
             * @param u Index of first node.
             * @param v Index of second node.
             */
            void addEdge(size_t u, size_t v) {
                if (u > v)
                    addEdge(v, u);
                else {
                    mates[u] = v;
                    mates[v] = u;
                    _edges++;
                    unmatched1 = SIZE_MAX;
                    unmatched2 = SIZE_MAX;
                }
            }

            /**
             * Remove the edge (u,v) from the matching.
             * @param u Index of first node.
             * @param v Index of second node.
             */
            void removeEdge(size_t u, size_t v) {
                if (u > v)
                    removeEdge(v, u);
                else {
                    mates[u] = SIZE_MAX;
                    mates[v] = SIZE_MAX;
                    _edges--;
                    unmatched1 = u;
                    unmatched2 = v;
                }
            }


            bool operator==(const BipartiteMatching &s) const {

                if (n != s.n)
                    return false;
                if (_edges != s._edges)
                    return false;

                const int ret = memcmp(&mates[0], &s.mates[0], n * sizeof(size_t));
                return (ret == 0);
            }

            bool operator<(const BipartiteMatching &s) const {

                if (n != s.n)
                    return false;
                if (_edges != s._edges)
                    return false;

                return memcmp(&mates[0], &s.mates[0], n * sizeof(size_t)) < 0;
            }

            /**
             * Return a hash value of this matching.
             * @return Hash value.
             */
            size_t hashValue() const {
                return boost::hash_range(mates.begin(), mates.end());
            }

            /**
             * Compare the matching with another one.
             * @param x Bipartite matching.
             * @return Zero, if x equals this. Negative value, if this < x. Positive value, if x < this.
             */
            int compare(const State &x) const override {

                // try to convert state to bipartite matching
                auto b = dynamic_cast<const BipartiteMatching *>(&x);

                if (b == nullptr)
                    return -1;

                if (n != b->n || _edges != b->_edges)
                    return -1;

                // element-wise comparison of mates
                const int res = memcmp(&mates[0], &b->mates[0], this->n * sizeof(int));
                return res;
            }

            /**
             * Create a compact string representation of the matching.
             * @return
             */
            std::string toString() const {

                std::stringstream ss;

                ss << "[";
                ss.width((int) log10(n) + 1);
                for (size_t i = 0; i < n - 1; i++) {
                    ss.width((int) log10(n) + 1);
                    if (mates[i] == SIZE_MAX)
                        ss << "-";
                    else
                        ss << mates[i];
                    ss << ", ";
                }
                ss.width((int) log10(n) + 1);
                if (mates[n - 1] == SIZE_MAX)
                    ss << "-";
                else
                    ss << mates[n - 1];
                ss << "]";

                return ss.str();
            }

            /**
             * Return a copy of the current object.
             * @return
             */
            std::unique_ptr<State> copy() const {
                return std::make_unique<BipartiteMatching>(*this);
            }

            /**
             * Is this a perfect matching?
             * @return
             */
            bool is_perfect() const {
                return 2 * _edges == n;
            }

            /**
             * Is this a near-perfect matching?
             * @return
             */
            bool is_near_perfect() const {
                return 2 * _edges == n - 2;
            }

            /**
             * Is vertex v matched?
             * @param v
             * @return
             */
            bool isMatched(size_t v) const {
                return mates[v] != SIZE_MAX;
            }

            /**
             * Return the matching partner of node v.
             * @param v Node index.
             * @return Matching partner of node v, or SIZE_MAX if v is unmatched.
             */
            size_t getMate(size_t v) const {
                return mates[v];
            }

            /**
             * Return the number of matching edges.
             * @return Number of matching edges.
             */
            size_t getNumberOfEdges() const {
                return _edges;
            }

            /**
             * Return the number of vertices in each vertex sets.
             * @return Number of vertices in each vertex set.
             */
            size_t getVertexSize() const {
                return n;
            }
        };

    }
}

// overload standard hash function of BipartiteMatching objects
namespace std {
    template<>
    struct hash<marathon::matching::BipartiteMatching> {
        size_t operator()(const marathon::matching::BipartiteMatching &x) const {
            return x.hashValue();
        }
    };
}


#endif /* STATE_H_ */

