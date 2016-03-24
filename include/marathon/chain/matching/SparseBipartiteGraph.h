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
