/*
 * State.h
 *
 *  Created on: Nov 24, 2014
 *      Author: rechner
 */

#ifndef STATE_H_
#define STATE_H_

// so dynamic_bitset becomes hashable
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include <cstdlib>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

namespace marathon {

namespace chain {

namespace sequence {

inline int COORD_TRANSFORM(int x, int y, int ld) {
	return (x*ld+y);
}

class DenseBipartiteGraph {

private:
	int m, n;
	boost::dynamic_bitset<> M;

public:

	DenseBipartiteGraph();
	DenseBipartiteGraph(int m, int n, const bool* bits = nullptr);
	DenseBipartiteGraph(const DenseBipartiteGraph& s);
	virtual ~DenseBipartiteGraph();

	int get_m() const;
	int get_n() const;
	bool has_edge(int u, int v) const;
	void flip_edge(int u, int v);
	bool is_switchable(int u1, int u2, int v1, int v2) const;
	void switch_4_cycle(int u1, int u2, int v1, int v2);

	void operator=(const DenseBipartiteGraph& s);
	bool operator==(const DenseBipartiteGraph &rhs) const;

	friend std::ostream& operator<<(std::ostream& os, const DenseBipartiteGraph& bip);
	friend size_t hash_value(const DenseBipartiteGraph& s);
};

}

}

}

#endif /* STATE_H_ */
