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

class DenseBipartiteGraph {

private:
	int nrows, ncols;
	boost::dynamic_bitset<> M;

public:

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

	void operator=(const DenseBipartiteGraph& s);
	bool operator==(const DenseBipartiteGraph &rhs) const;

	friend std::ostream& operator<<(std::ostream& os,
			const DenseBipartiteGraph& bip);
};

inline int COORD_TRANSFORM(int x, int y, int ld) {
	return (x * ld + y);
}

}
}
}

// overload hash function for this class
namespace std {
template<>
struct hash<::marathon::chain::sequence::DenseBipartiteGraph> {
	std::size_t operator()(
			const marathon::chain::sequence::DenseBipartiteGraph& s) const {
		return s.hash_value();
	}
};
}

#endif /* STATE_H_ */
