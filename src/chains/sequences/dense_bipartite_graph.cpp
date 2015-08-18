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
		m(0), n(0) {

}

DenseBipartiteGraph::DenseBipartiteGraph(int m, int n, const bool* bits) :
		m(m), n(n) {
	if (bits == nullptr) {
		M = boost::dynamic_bitset<>(m * n);
	} else {
		for (int i = 0; i < m * n; i++) {
			M.push_back(bits[i]);
		}
	}
}

DenseBipartiteGraph::DenseBipartiteGraph(const DenseBipartiteGraph& s) :
		m(s.m), n(s.n) {
	M = boost::dynamic_bitset<>(s.M);
}

DenseBipartiteGraph::~DenseBipartiteGraph() {

}

int DenseBipartiteGraph::get_m() const {
	return m;
}

int DenseBipartiteGraph::get_n() const {
	return n;
}

bool DenseBipartiteGraph::has_edge(int u, int v) const {
	return M[COORD_TRANSFORM(u, v, n)];
}

void DenseBipartiteGraph::flip_edge(int u, int v) {
	M[COORD_TRANSFORM(u, v, n)].flip();
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

	v1 -= m;
	v2 -= m;

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
	v1 -= m;
	v2 -= m;

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
	m = s.m;
	n = s.n;
	M = boost::dynamic_bitset<>(s.M);
}

bool DenseBipartiteGraph::operator==(const DenseBipartiteGraph &s) const {
	return m == s.m && n == s.n && M == s.M;
}

std::ostream& operator<<(std::ostream& os, const DenseBipartiteGraph& bip) {

	for (unsigned int i = 0; i < bip.M.size(); i++)
		os << bip.M[i];

	return os;
}

size_t hash_value(const DenseBipartiteGraph& s) {
	return boost::hash_value(s.M.m_bits);
	//return boost::hash_range(s.M, s.M + s.m * s.n);
}

}

}

}
