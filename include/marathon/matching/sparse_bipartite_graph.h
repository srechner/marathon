/*
 * Created on: 22.06.2013
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

#ifndef BIPARTITEGRAPH_H_
#define BIPARTITEGRAPH_H_

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS  // so dynamic_bitset becomes hashable

#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_set.hpp>

namespace marathon {
    namespace matching {

        typedef std::pair<size_t, size_t> edge;
        typedef std::vector<edge> edgelist;

        /**
         * Sparse representation of a bipartite graph.
         */
        class SparseBipartiteGraph {

        private:

            std::vector<size_t> adj;      // adj[i] is the index of the first node adjacent to i
            std::vector<size_t> edges;    // end points of edges

            /**
             * Create a bitset representation of the bipartite graph.
             * @return Bitset that represents the bi-adjacency matrix of this bipartite graph.
             */
            boost::dynamic_bitset<> convert_to_bitset() const {

                const size_t n = getNumberOfNodes() / 2;

                boost::dynamic_bitset<> bits(n*n);
                for (size_t u = 0; u < n; u++) {
                    for (size_t k = adj[u]; k < adj[u + 1]; k++) {
                        size_t v = edges[k];
                        bits[u * n + (v - n)] = 1;
                    }
                }

                return bits;
            }

        public:

            /**
             * Create a bipartite graph from an instance encoding.
             *
             * @param inst Binary string of length n*n.
             *
             * Such a 0-1-String is interpreted as the bi-adjacency matrix M = (m_ij) of a bipartite graph G=(V,E).
             * The bi-adjacency matrix M is defined as m_ij = 1, if (i,j) is in E, or 0, otherwise. The rows of M
             * are concatenated to a single string.
             * Thus, the input string "110101011" corresponds to the biadjacency  matrix
             *
             *  1 1 0
             *  1 0 1
             *  0 1 1
             *
             *  which is equivalent to the graph
             *
             *  u1  u2  u3
             *  |\ / \ /|
             *  | X   X |
             *  |/ \ / \|
             *  v1  v2  v3
             *
             */
            SparseBipartiteGraph(const std::string &inst) {

                size_t n = (size_t) sqrt(inst.length());

                if (n * n != inst.length()) {
                    throw std::runtime_error("Error! Malformed input instance: " + inst);
                }

                adj.resize(2 * n + 1);
                size_t k = 0;
                for (size_t u = 0; u < n; u++) {
                    for (size_t v = n; v < 2 * n; v++) {
                        char c = inst.at(k);
                        if (c == '1') {
                            edges.push_back(v);
                            adj[u + 1] = edges.size();
                        }
                        k++;
                    }
                }

                for (size_t v = n; v < 2 * n; v++) {
                    for (size_t u = 0; u < n; u++) {
                        size_t k = u * n + (v - n);
                        char c = inst.at(k);
                        if (c == '1') {
                            edges.push_back(u);
                            adj[v + 1] =  edges.size();
                        }
                    }
                }
            }

            /**
             * Create a bipartite graph from a bitset.
             * @param bits Bitset of length n*n.
             *
             * Such a bitset is interpreted as the bi-adjacency matrix M = (m_ij) of a bipartite graph G=(V,E).
             * The bi-adjacency matrix M is defined as m_ij = 1, if (i,j) is in E, or 0, otherwise. The rows of M
             * are concatenated to a single string.
             * Thus, the bitset containing "110101011" corresponds to the biadjacency  matrix
             *
             *  1 1 0
             *  1 0 1
             *  0 1 1
             *
             *  which is equivalent to graph
             *
             *  u1  u2  u3
             *  |\ / \ /|
             *  | X   X |
             *  |/ \ / \|
             *  v1  v2  v3
             */
            SparseBipartiteGraph(const boost::dynamic_bitset<> &bits) {

                size_t n = (size_t) sqrt(bits.size());

                if (n * n != bits.size()) {
                    throw std::runtime_error("Error! Malformed input instance.");
                }

                adj.resize(2 * n + 1);
                size_t k = 0;
                for (size_t u = 0; u < n; u++) {
                    for (size_t v = n; v < 2 * n; v++) {
                        if (bits[k]) {
                            edges.push_back(v);
                            adj[u + 1] = edges.size();
                        }
                        k++;
                    }
                }

                for (size_t v = n; v < 2 * n; v++) {
                    for (size_t u = 0; u < n; u++) {
                        size_t k = u * n + (v - n);
                        if (bits[k]) {
                            edges.push_back(u);
                            adj[v + 1] = edges.size();
                        }
                    }
                }
            }

            /**
             * Is the edge (u,v) included?
             * @param u Node index.
             * @param v Node index.
             * @return True, if (u,v) is included. False, otherwise.
             */
            bool hasEdge(size_t u, size_t v) const {
                //return boost::edge(u, v, g).second;

                // apply binary search to scan adjacency list
                size_t l = adj[u];
                size_t r = adj[u + 1];
                while (l < r) {
                    const size_t m = (l + r) / 2;
                    if (edges[m] == v)
                        return true;
                    else if (edges[m] < v)
                        l = m + 1;
                    else
                        r = m;
                }
                return false;
            }

            /**
             * Return the total number of nodes.
             * @return
             */
            size_t getNumberOfNodes() const {
                return adj.size() - 1;
            }

            /**
             * Return the total number of edges.
             * @return
             */
            size_t getNumberOfEdges() const {
                return edges.size();
            }

            /**
             * Compile a vector of edges.
             * @return Vector of edges.
             */
            edgelist
            getEdges() const {

                edgelist E;

                for (size_t u = 0; u < adj.size() - 1; u++) {
                    for (size_t i = adj[u]; i < adj[u + 1]; i++) {
                        size_t v = edges[i];
                        E.push_back(std::make_pair(u, v));
                    }
                }

                return E;
            }


            /**
             * Return the set of adjacent states.
             * @param v Node index.
             * @return List of adjacent states.
             */
            std::vector<size_t>
            getNeighbors(size_t v) const {

                std::vector<size_t> neighbors;

                for (size_t i = adj[v]; i < adj[v + 1]; i++) {
                    neighbors.push_back(edges[i]);
                }

                return neighbors;
            }

            /**
             * Return a unique string representation.
             * @return
             */
            std::string toString() const {

                boost::dynamic_bitset<> bits = convert_to_bitset();

                std::string buf;
                boost::to_string(bits, buf);

                return buf;
            }

            /**
             * Return the degree of node u.
             * @param u Node index.
             * @return Number of nodes adjacent to u.
             */
            size_t getDegree(size_t u) const {
                return adj[u + 1] - adj[u];
            }

            /**
             * Return the index of the first edge including node v.
             * @param v node index.
             * @return
             */
            size_t getIndexOfFirstEdge(size_t v) const {
                return adj[v];
            }

            /**
             * Return the index of the last edge including node v.
             * @param v node index.
             * @return
             */
            size_t getIndexOfLastEdge(size_t v) const {
                return adj[v + 1] - 1;
            }

            /**
             * Determine whether two bipartite graphs are identical.
             * @param s Bipartite graph.
             * @return
             */
            bool operator==(const SparseBipartiteGraph &g) const {
                return adj == g.adj && edges == g.edges;
            }

            /**
             * Create a hash value of the bipartite graph.
             * @return hash value
             */
            size_t hashValue() const {
                return boost::hash_range(adj.begin(), adj.end()) ^
                       boost::hash_range(edges.begin(), edges.end());
            }
        };

        inline
        std::ostream &operator<<(std::ostream &out, const SparseBipartiteGraph &bip) {
            out << bip.toString();
            return out;
        }
    }
}

// overload standard hash function of BinaryMatrix objects
namespace std {
    template<>
    struct hash<marathon::matching::SparseBipartiteGraph> {
        size_t operator()(const marathon::matching::SparseBipartiteGraph &x) const {
            return x.hashValue();
        }
    };
}


#endif /* BIPARTITEGRAPH_H_ */
