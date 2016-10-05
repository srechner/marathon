/*
 * SparseBipartiteGraph.h
 *
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

#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_set.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

namespace marathon {
	namespace matching {

		// undirected graph
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> undirected_graph;
		typedef std::pair<uint, uint> edge;
		typedef std::vector<edge> edgelist;

		class SparseBipartiteGraph {

		private:
			undirected_graph g;

		public:
			SparseBipartiteGraph(std::string hash) {

				int n = 2 * (int) sqrt(hash.length());

				g = undirected_graph(n);

				for (int k = 0; k < hash.length(); k++) {

					//if (hash.at(hash.length() - k - 1) == '1') {
					if (hash.at(k) == '1') {
						int u = k / (n / 2);
						int v = (k % (n / 2)) + (n / 2);
						add_edge(u, v, g);
					}
				}
			}

			SparseBipartiteGraph(const SparseBipartiteGraph &b) {
				g = undirected_graph(b.g);
			}

			SparseBipartiteGraph(size_t n) {
				g = undirected_graph(n);
			}

			~SparseBipartiteGraph() {
			}

			void addEdge(int u, int v) {
				add_edge(u, v, g);
			}

			void removeEdge(int u, int v) {
				boost::remove_edge(u, v, g);
			}

			bool hasEdge(int u, int v) const {
				return boost::edge(u, v, g).second;
			}

			unsigned int getNumberOfNodes() const {
				return boost::num_vertices(g);
			}

			unsigned int getNumberOfEdges() const {
				return boost::num_edges(g);
			}

			void getEdges(edgelist &E) const {

				E.clear();

				typename boost::graph_traits<undirected_graph>::edge_iterator ei, ei_end;

				for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
					edge e(source(*ei, g), target(*ei, g));
					E.push_back(e);
				}
			}

			void getNeighbors(int v,
			                  std::vector<int> &neighbors) const {

				neighbors.clear();

				typename boost::graph_traits<undirected_graph>::adjacency_iterator vi,
						vi_end;
				typename boost::graph_traits<undirected_graph>::vertex_descriptor vv =
						vertex(v, g);

				for (boost::tie(vi, vi_end) = boost::adjacent_vertices(vv, g); vi != vi_end;
				     ++vi) {
					neighbors.push_back(*vi);
				}

			}

			void convert_to_bitset(
					boost::dynamic_bitset<> &bits) const {

				std::vector<int> adj;

				unsigned int n = getNumberOfNodes() / 2;

				bits.resize(n * n);

				for (int u = 0; u < n; u++) {

					getNeighbors(u, adj);

					for (std::vector<int>::iterator it = adj.begin(); it != adj.end();
					     ++it) {

						int v = *it;

						unsigned int i = n * u + v - n;

						bits[n * n - i - 1] = 1;
					}
				}
			}

			std::string toString() const {

				boost::dynamic_bitset<> bits;
				convert_to_bitset(bits);

				std::string buf;
				boost::to_string(bits, buf);

				return buf;
			}


			void cardmax_matching(std::vector<int> &mates) const {

				int n = getNumberOfNodes();

				mates.clear();
				mates.resize(n);

				std::vector<boost::graph_traits<undirected_graph>::vertex_descriptor> mateys(
						n);

				bool success = checked_edmonds_maximum_cardinality_matching(g, &mateys[0]);
				assert(success);

				boost::graph_traits<undirected_graph>::vertex_iterator vi, vi_end;
				for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
					if (mateys[*vi]
					    != boost::graph_traits<undirected_graph>::null_vertex()) {
						//std::cout << "{" << *vi << ", " << mateys[*vi] << "}" << std::endl;
						mates[*vi] = mateys[*vi];
						mates[mateys[*vi]] = *vi;
					} else
						mates[*vi] = n;
				}
			}
		};

		inline
		std::ostream &operator<<(std::ostream &out, const SparseBipartiteGraph &bip) {
			out << bip.toString() << std::endl;
			return out;
		}

	}
}

#endif /* BIPARTITEGRAPH_H_ */
