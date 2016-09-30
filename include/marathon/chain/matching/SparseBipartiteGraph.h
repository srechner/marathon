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

namespace marathon {
namespace chain {
namespace matching {

// undirected graph
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> undirected_graph;
typedef std::pair<uint, uint> edge;
typedef std::vector<edge> edgelist;

class SparseBipartiteGraph {

private:
	undirected_graph g;

public:
	SparseBipartiteGraph(const SparseBipartiteGraph& b);
	SparseBipartiteGraph(size_t n);
	SparseBipartiteGraph(std::string hash);
	virtual ~SparseBipartiteGraph();

	unsigned int getNumberOfNodes() const;
	unsigned int getNumberOfEdges() const;
	void getEdges(edgelist& edges) const;

	void addEdge(int u, int v);
	bool hasEdge(int u, int v) const;
	void removeEdge(int u, int v);

	void getNeighbors(int v, std::vector<int>& neighbors) const;

	void cardmax_matching(std::vector<int>& mates) const;

	std::string toString() const;

	void convert_to_bitset(boost::dynamic_bitset<>&) const;

	friend std::ostream& operator<<(std::ostream& os, const SparseBipartiteGraph& bip);
};

}
}
}

#endif /* BIPARTITEGRAPH_H_ */
