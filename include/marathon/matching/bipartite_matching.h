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
         * Class reprentation of matchings in a bipartite graph.
         */
        class BipartiteMatching : public State {

        public:

            int n, k;            // number of nodes			// number of edges
            int unmatched[2];    // indices of unmatched nodes (if any)
            int *mates;

            /**
             * Create empty matching in a bipartite graph with n nodes.
             * @param n
             */
            BipartiteMatching(int n = 0) :
                    n(n), k(0) {
                mates = (int *) malloc(n * sizeof(int));
                for (int i = 0; i < n; i++)
                    mates[i] = -1;
            }

            /**
             * Create a matching as a copy of the matching s.
             * @param s Bipartite matching.
             */
            BipartiteMatching(const BipartiteMatching &s) :
                    n(s.n), k(s.k) {
                unmatched[0] = s.unmatched[0];
                unmatched[1] = s.unmatched[1];
                mates = (int *) malloc(n * sizeof(int));
                memcpy(mates, s.mates, n * sizeof(int));
            }


            /**
             * Create a matching of size n.
             * @param n Number of nodes.
             * @param k Number of matching edges.
             * @param unmatched1 Index of first unmatched node (or -1, if none).
             * @param unmatched2 Index of second unmatched node (or -1, if none).
             * @param matching Matching object.
             */
            BipartiteMatching(
                    int n, int k, int unmatched1, int unmatched2, int *matching) :
                    n(n), k(k) {
                this->unmatched[0] = unmatched[0];
                this->unmatched[1] = unmatched[1];
                mates = (int *) malloc(n * sizeof(int));
                memcpy(this->mates, matching, n * sizeof(int));
            }

            ~BipartiteMatching() {
                free(mates);
            }

            /**
             * Add the edge (u,v) to the matching.
             * @param u Index of first node.
             * @param v Index of second node.
             */
            void addEdge(int u, int v) {
                if (u > v)
                    addEdge(v, u);
                else {
                    mates[u] = v;
                    mates[v] = u;
                    k++;
                    unmatched[0] = -1;
                    unmatched[1] = -1;
                }
            }

            /**
             * Remove the edge (u,v) from the matching.
             * @param u Index of first node.
             * @param v Index of second node.
             */
            void removeEdge(int u, int v) {
                if (u > v)
                    removeEdge(v, u);
                else {
                    mates[u] = -1;
                    mates[v] = -1;
                    k--;
                    unmatched[0] = u;
                    unmatched[1] = v;
                }
            }

/*
            void operator=(BipartiteMatching const &s) {
                n = s.n;
                k = s.k;
                unmatched[0] = s.unmatched[0];
                unmatched[1] = s.unmatched[1];
                mates = (int *) realloc(mates, n * sizeof(int));
                memcpy(mates, s.mates, n * sizeof(int));
            }*/

            bool operator==(const BipartiteMatching &s) const {
                if (n != s.n)
                    return false;
                if (k != s.k)
                    return false;
                const int ret = memcmp(mates, s.mates, n * sizeof(int));
                return (ret == 0);
            }

            bool operator<(const BipartiteMatching &s) const {
                return memcmp(mates, s.mates, n * sizeof(int)) < 0;
            }

            size_t hashValue() const {
                //return unmatched[0] * n + unmatched[1];
                return boost::hash_range(this->mates, this->mates + n);
            }

            int compare(const State *x) const override {
                const BipartiteMatching *b = (const BipartiteMatching *) x;
                const int res = memcmp(this->mates, b->mates, this->n * sizeof(int));
                return res;
            }

            /**
             * Create a compact string representation of the matching.
             * @return
             */
            std::string toString() const {

                std::stringstream ss;

                ss << "[";
                ss.width(log10(n) + 1);
                for (int i = 0; i < n - 1; i++) {
                    ss.width(log10(n) + 1);
                    if (mates[i] == -1)
                        ss << "-";
                    else
                        ss << mates[i];
                    ss << ", ";
                }
                ss.width(log10(n) + 1);
                if (mates[n - 1] == -1)
                    ss << "-";
                else
                    ss << mates[n - 1];
                ss << "]";

                return ss.str();
            }

            State *copy() const {
                return new BipartiteMatching(*this);
            }

            bool is_perfect() const {
                return 2 * k == n;
            }

            bool is_near_perfect() const {
                return 2 * k == n - 2;
            }

            bool isMatched(int v) const {
                return mates[v] != -1;
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

