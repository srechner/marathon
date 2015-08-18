/*
 * BipartiteGraph.cpp
 *
 *  Created on: 22.06.2013
 *      Author: steffen
 */

// project includes
#include "../../../include/marathon/chains/matching/sparse_bipartite_graph.h"

// system includes
#include <cmath>
#include <boost/graph/max_cardinality_matching.hpp>

namespace marathon {

namespace chain {

namespace matching {

SparseBipartiteGraph::SparseBipartiteGraph(std::string hash) {

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

SparseBipartiteGraph::SparseBipartiteGraph(const SparseBipartiteGraph& b) {
	g = undirected_graph(b.g);
}

SparseBipartiteGraph::SparseBipartiteGraph(size_t n) {
	g = undirected_graph(n);
}

SparseBipartiteGraph::~SparseBipartiteGraph() {
}

void SparseBipartiteGraph::addEdge(int u, int v) {
	add_edge(u, v, g);
}

void SparseBipartiteGraph::removeEdge(int u, int v) {
	boost::remove_edge(u, v, g);
}

bool SparseBipartiteGraph::hasEdge(int u, int v) const {
	return boost::edge(u, v, g).second;
}

unsigned int SparseBipartiteGraph::getNumberOfNodes() const {
	return boost::num_vertices(g);
}

unsigned int SparseBipartiteGraph::getNumberOfEdges() const {
	return boost::num_edges(g);
}

void SparseBipartiteGraph::getEdges(edgelist& E) const {

	E.clear();

	typename boost::graph_traits<undirected_graph>::edge_iterator ei, ei_end;

	for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		edge e(source(*ei, g), target(*ei, g));
		E.push_back(e);
	}
}

void SparseBipartiteGraph::getNeighbors(int v,
		std::vector<int>& neighbors) const {

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

void SparseBipartiteGraph::convert_to_bitset(
		boost::dynamic_bitset<>& bits) const {

	std::vector < int > adj;

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

std::string SparseBipartiteGraph::toString() const {

	boost::dynamic_bitset<> bits;
	convert_to_bitset(bits);

	std::string buf;
	boost::to_string(bits, buf);

	return buf;
}

std::ostream& operator<<(std::ostream& out, const SparseBipartiteGraph& bip) {
	out << bip.toString() << std::endl;
	return out;
}

void SparseBipartiteGraph::cardmax_matching(std::vector<int>& mates) const {

	int n = getNumberOfNodes();

	mates.clear();
	mates.resize(n);

	std::vector < boost::graph_traits < undirected_graph > ::vertex_descriptor
			> mateys(n);

	bool success = checked_edmonds_maximum_cardinality_matching(g, &mateys[0]);
	assert(success);

	boost::graph_traits<undirected_graph>::vertex_iterator vi, vi_end;
	for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
		if (mateys[*vi]
				!= boost::graph_traits < undirected_graph > ::null_vertex()) {
			//std::cout << "{" << *vi << ", " << mateys[*vi] << "}" << std::endl;
			mates[*vi] = mateys[*vi];
			mates[mateys[*vi]] = *vi;
		} else
			mates[*vi] = n;
	}
}

}

}

}
