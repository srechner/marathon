/*
 * Created on: Jan 12, 2018
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


#ifndef MARATHON_MATCHING_REALIZATION_H
#define MARATHON_MATCHING_REALIZATION_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include "marathon/matching/bipartite_matching.h"
#include "marathon/matching/sparse_bipartite_graph.h"

namespace marathon {
    namespace matching {

        /**
         * Create an arbitrary perfect matching in the bipartite graph g.
         * @param g Bipartite graph.
         * @return Perfect matching in g.
         */
        BipartiteMatching
        cardmax_matching(const SparseBipartiteGraph &g) {

            const size_t n = g.getNumberOfNodes();

            // transform graph into boost graph representation
            typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> undirected_graph;
            undirected_graph gg(n);
            edgelist edges = g.getEdges();
            for (edge e : edges) {
                const int u = e.first;
                const int v = e.second;
                add_edge(u, v, gg);
            }

            // calculate perfect matching
            std::vector<boost::graph_traits<undirected_graph>::vertex_descriptor> mateys(n);
            bool success = checked_edmonds_maximum_cardinality_matching(gg, &mateys[0]);
            assert(success);

            std::vector<size_t> mates(n);

            boost::graph_traits<undirected_graph>::vertex_iterator vi, vi_end;
            for (boost::tie(vi, vi_end) = vertices(gg); vi != vi_end; ++vi) {
                if (mateys[*vi] != boost::graph_traits<undirected_graph>::null_vertex()) {
                    //std::cout << "{" << *vi << ", " << mateys[*vi] << "}" << std::endl;
                    mates[*vi] = mateys[*vi];
                    mates[mateys[*vi]] = *vi;
                } else {
                    mates[*vi] = SIZE_MAX;
                }
            }


            // Count Number of Matching Edges
            size_t k = 0;
            for (size_t i = 0; i < n; i++) {
                if (mates[i] != SIZE_MAX) {
                    k++;
                }
            }

            if (k != n) {
                throw std::runtime_error("Bipartite graph doesn't have a perfect matching!");
            } else {
                return BipartiteMatching(n, k / 2, mates, SIZE_MAX, SIZE_MAX);
            }

        }
    }
}

#endif //MARATHON_MATCHING_REALIZATION_H
