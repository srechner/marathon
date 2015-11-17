/*
 * State.cpp
 *
 *  Created on: Nov 24, 2014
 *      Author: rechner
 */

// project includes
#include "../../../include/marathon/chains/sequences/dense_bipartite_graph.h"

// system includes
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <functional>

//#define DEBUG

namespace marathon {

namespace chain {

namespace sequence {

DenseBipartiteGraph::DenseBipartiteGraph() :
		nrows(0), ncols(0) {

}

DenseBipartiteGraph::DenseBipartiteGraph(int nrows, int ncols, const bool* bits) :
		nrows(nrows), ncols(ncols) {
	if (bits == nullptr) {
		M = boost::dynamic_bitset<>(nrows * ncols);
	} else {
		for (int i = 0; i < nrows * ncols; i++) {
			M.push_back(bits[i]);
		}
	}
}

DenseBipartiteGraph::DenseBipartiteGraph(int nrows, int ncols,
		const std::string& str) :
		nrows(nrows), ncols(ncols) {
	M = boost::dynamic_bitset<>(str);
}

DenseBipartiteGraph::DenseBipartiteGraph(const DenseBipartiteGraph& s) :
		nrows(s.nrows), ncols(s.ncols) {
	M = boost::dynamic_bitset<>(s.M);
}

DenseBipartiteGraph::~DenseBipartiteGraph() {

}

int DenseBipartiteGraph::get_nrows() const {
	return nrows;
}

int DenseBipartiteGraph::get_ncols() const {
	return ncols;
}

void DenseBipartiteGraph::set_edge(int u, int v, bool b) {
	M[COORD_TRANSFORM(u, v, ncols)] = b;
}

bool DenseBipartiteGraph::has_edge(int u, int v) const {
	return M[COORD_TRANSFORM(u, v, ncols)];
}

void DenseBipartiteGraph::flip_edge(int u, int v) {
	M[COORD_TRANSFORM(u, v, ncols)].flip();
}

bool DenseBipartiteGraph::is_switchable(int u1, int v1, int u2, int v2) const {

	/* translate node labels */
	if (u1 > v1) {
		// rotate
		uint tmp = u1;
		u1 = v1;
		v1 = u2;
		u2 = v2;
		v2 = tmp;
	}

	v1 -= nrows;
	v2 -= nrows;

#ifdef DEBUG
	std::cout << "is switchable " << u1 << " " << v1 << " " << u2 << " " << v2
	<< std::endl;

	std::cout << "(" << u1 << "," << v1 << "): " << has_edge(u1, v1)
	<< std::endl;
	std::cout << "(" << u1 << "," << v2 << "): " << has_edge(u1, v2)
	<< std::endl;
	std::cout << "(" << u2 << "," << v1 << "): " << has_edge(u2, v1)
	<< std::endl;
	std::cout << "(" << u2 << "," << v2 << "): " << has_edge(u2, v2)
	<< std::endl;
#endif

	bool a = has_edge(u1, v1);
	bool b = has_edge(u1, v2);
	bool c = has_edge(u2, v1);
	bool d = has_edge(u2, v2);

	return (a == d) && (b == c) && (a == !b);
}

void DenseBipartiteGraph::switch_4_cycle(int u1, int v1, int u2, int v2) {

	if (u1 > v1) {
		// rotate
		uint tmp = u1;
		u1 = v1;
		v1 = u2;
		u2 = v2;
		v2 = tmp;
	}

	// translate node labels
	v1 -= nrows;
	v2 -= nrows;

#ifdef DEBUG
	std::cout << "switch cycle [ " << u1 << " " << v1 << " " << u2 << " " << v2
	<< " ]" << std::endl;
#endif

	flip_edge(u1, v1);
	flip_edge(u1, v2);
	flip_edge(u2, v1);
	flip_edge(u2, v2);
}

void DenseBipartiteGraph::operator=(const DenseBipartiteGraph& s) {
	nrows = s.nrows;
	ncols = s.ncols;
	M = boost::dynamic_bitset<>(s.M);
}

bool DenseBipartiteGraph::operator==(const DenseBipartiteGraph &s) const {
	return nrows == s.nrows && ncols == s.ncols && M == s.M;
}

std::ostream& operator<<(std::ostream& os, const DenseBipartiteGraph& bip) {

	for (unsigned int i = 0; i < bip.M.size(); i++)
		os << bip.M[i];

	return os;
}

size_t hash_value(const DenseBipartiteGraph& s) {
	return boost::hash_value(s.M.m_bits);
	//return boost::hash_range(s.M, s.M + s.m * s.ncols);
}

void DenseBipartiteGraph::get_row(int u, boost::dynamic_bitset<>& row) const {
	row.resize(ncols);
	for (int v = 0; v < ncols; v++) {
		row[v] = has_edge(u, v);
	}
}

}

}

}
