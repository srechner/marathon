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
	bool operator==(const BipartiteMatching &rhs) const;
	size_t hash_value() const;

	bool is_perfect() const;
	bool is_near_perfect() const;

	friend std::ostream& operator<<(std::ostream& os, BipartiteMatching& bip);
};

}
}
}

// overload hash function for this class
namespace std {
template<>
struct hash<::marathon::chain::matching::BipartiteMatching> {

	typedef ::marathon::chain::matching::BipartiteMatching argument_type;
	typedef size_t result_type;

	result_type operator ()(const argument_type& x) const {
		return x.hash_value();
	}
};
}




#endif /* STATE_H_ */

