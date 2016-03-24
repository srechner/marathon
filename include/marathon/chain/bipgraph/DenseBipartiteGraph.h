/*
 * State.h
 *
 *  Created on: Nov 24, 2014
 *      Author: rechner
 */

#ifndef STATE_H_
#define STATE_H_

#include <cstdlib>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

#include "../../State.h"

namespace marathon {
namespace chain {
namespace bipgraph {

class DenseBipartiteGraph: public State {

public:

	static int COORD_TRANSFORM(const int x, const int y, const int ld) {
		return (x * ld + y);
	}

	int nrows, ncols;
	boost::dynamic_bitset<> M;

	DenseBipartiteGraph();
	DenseBipartiteGraph(int nrows, int ncols, const bool* bits = nullptr);
	DenseBipartiteGraph(int nrows, int ncols, const std::string& str);
	DenseBipartiteGraph(const DenseBipartiteGraph& s);
	virtual ~DenseBipartiteGraph();

	int get_nrows() const;
	int get_ncols() const;
	bool has_edge(int u, int v) const;
	void flip_edge(int u, int v);
	void set_edge(int u, int v, bool);
	bool is_switchable(int u1, int u2, int v1, int v2) const;
	void switch_4_cycle(int u1, int u2, int v1, int v2);

	void get_row(int u, boost::dynamic_bitset<>& row) const;

	size_t hash_value() const;
	int compare_to(const State* x) const;
	std::string to_string() const;

	void operator=(const DenseBipartiteGraph& s);
	bool operator<(const DenseBipartiteGraph &rhs) const;
	bool operator==(const DenseBipartiteGraph &rhs) const;

};

}
}
}

#endif /* STATE_H_ */
