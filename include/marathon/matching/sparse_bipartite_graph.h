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

        typedef std::pair<int, int> edge;
        typedef std::vector<edge> edgelist;

        /**
         * Sparse representation of a bipartite graph.
         */
        class SparseBipartiteGraph {

        private:

            std::vector<int> adj;      // adj[i] is the index of the first node adjacent to i
            std::vector<int> edges;    // end points of edges

            void convert_to_bitset(
                    boost::dynamic_bitset<> &bits) const {

                //std::vector<int> adj;

                size_t n = getNumberOfNodes() / 2;

                bits.resize(n * n);

                for (int u = 0; u < n; u++) {

                    for (int k = adj[u]; k < adj[u + 1]; k++) {
                        int v = edges[k];
                        bits[u * n + (v - n)] = 1;
                    }
                }
            }

        public:

            SparseBipartiteGraph(std::string inst) {

                int n = (int) sqrt(inst.length());

                adj.resize(2 * n + 1);
                int k = 0;
                for (int u = 0; u < n; u++) {
                    for (int v = n; v < 2 * n; v++) {
                        char c = inst.at(k);
                        if (c == '1') {
                            edges.push_back(v);
                            adj[u + 1] = (int) edges.size();
                        }
                        k++;
                    }
                }

                for (int v = n; v < 2 * n; v++) {
                    for (int u = 0; u < n; u++) {
                        int k = u * n + (v - n);
                        char c = inst.at(k);
                        if (c == '1') {
                            edges.push_back(u);
                            adj[v + 1] = (int) edges.size();
                        }
                    }
                }
            }

            SparseBipartiteGraph(const SparseBipartiteGraph &b) {
                //g = undirected_graph(b.g);
                adj = b.adj;
                edges = b.edges;
            }

            SparseBipartiteGraph(const boost::dynamic_bitset<> &bits) {

                int n = (int) sqrt(bits.size());

                adj.resize(2 * n + 1);
                int k = 0;
                for (int u = 0; u < n; u++) {
                    for (int v = n; v < 2 * n; v++) {
                        if (bits[k]) {
                            edges.push_back(v);
                            adj[u + 1] = (int) edges.size();
                        }
                        k++;
                    }
                }

                for (int v = n; v < 2 * n; v++) {
                    for (int u = 0; u < n; u++) {
                        int k = u * n + (v - n);
                        if (bits[k]) {
                            edges.push_back(u);
                            adj[v + 1] = (int) edges.size();
                        }
                    }
                }
            }

            bool hasEdge(const int u, const int v) const {
                //return boost::edge(u, v, g).second;

                // apply binary search to scan adjacency list
                int l = adj[u];
                int r = adj[u + 1];
                while (l < r) {
                    const int m = (l + r) / 2;
                    if (edges[m] == v)
                        return true;
                    else if (edges[m] < v)
                        l = m + 1;
                    else
                        r = m;
                }
                return false;
            }

            size_t getNumberOfNodes() const {
                //return boost::num_vertices(g);
                return adj.size() - 1;
            }

            size_t getNumberOfEdges() const {
                return edges.size();
            }

            void getEdges(edgelist &E) const {

                E.clear();

                for (int u = 0; u < adj.size() - 1; u++) {
                    for (int i = adj[u]; i < adj[u + 1]; i++) {
                        int v = edges[i];
                        E.push_back(std::make_pair(u, v));
                    }
                }
            }


            /**
             * Return the set of adjacent states.
             * @param v Node index.
             * @param neighbors (out parameter) List of adjacent states.
             */
            void getNeighbors(int v, std::vector<int> &neighbors) const {

                neighbors.clear();

                for (int i = adj[v]; i < adj[v + 1]; i++) {
                    neighbors.push_back(edges[i]);
                }

            }

            /**
             * Return a unique string representation.
             * @return
             */
            std::string toString() const {

                boost::dynamic_bitset<> bits;
                convert_to_bitset(bits);

                std::string buf;
                boost::to_string(bits, buf);

                return buf;
            }

            /**
             * Return the degree of node u.
             * @param u Node index.
             * @return
             */
            int getDegree(const int u) const {
                return adj[u + 1] - adj[u];
            }

            /**
             * Return the index of the first edge including node v.
             * @param v node index.
             * @return
             */
            int getIndexOfFirstEdge(const int v) const {
                return adj[v];
            }

            /**
             * Return the index of the last edge including node v.
             * @param v node index.
             * @return
             */
            int getIndexOfLastEdge(const int v) const {
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
