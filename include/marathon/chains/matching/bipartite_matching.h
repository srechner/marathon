/*
 * State.h
 *
 *  Created on: Nov 20, 2014
 *      Author: rechner
 */

#ifndef STATE_JS89_H_
#define STATE_JS89_H_

#include <cstring>
#include <cstdlib>
#include <ostream>

namespace marathon {

namespace chain {

namespace matching {

class BipartiteMatching {

public:

	int n;				// number of nodes
	int k;				// number of edges
	int unmatched[2];	// indices of unmatched nodes (if any)
	int* mates;

	BipartiteMatching();

	BipartiteMatching(const BipartiteMatching& s);
	BipartiteMatching(int n, int k, int unmatched[2], int* matching);
	virtual ~BipartiteMatching();

	void addEdge(int u, int v);
	void removeEdge(int u, int v);

	void operator=(BipartiteMatching const& s);
	bool operator==(BipartiteMatching const&rhs) const;

	bool is_perfect() const;
	bool is_near_perfect() const;

	friend std::ostream& operator<<(std::ostream& os, BipartiteMatching& bip);
};

size_t hash_value(const BipartiteMatching& s);

}

}

}
#endif /* STATE_H_ */

