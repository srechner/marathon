/*
 * BipartiteGraph.h
 *
 *  Created on: 22.06.2013
 *      Author: steffen
 */

#ifndef BIPARTITEGRAPH_H_
#define BIPARTITEGRAPH_H_

#include <vector>
#include <iostream>

#include "undirected_graph.h"

namespace marathon {

namespace chain {

namespace matching {

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
