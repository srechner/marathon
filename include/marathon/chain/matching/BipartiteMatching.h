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

#include "../../State.h"

namespace marathon {
namespace chain {
namespace matching {

class BipartiteMatching: public State {

public:

	int n, k;				// number of nodes			// number of edges
	int unmatched[2];	// indices of unmatched nodes (if any)
	int* mates;

	BipartiteMatching();

	BipartiteMatching(const BipartiteMatching& s);
	BipartiteMatching(int n, int k, int unmatched[2], int* matching);
	virtual ~BipartiteMatching();

	void addEdge(int u, int v);
	void removeEdge(int u, int v);

	void operator=(BipartiteMatching const& s);
	bool operator==(const BipartiteMatching &s) const;
	bool operator<(const BipartiteMatching& s) const;
	size_t hashValue() const;
	int compare(const State*) const;
	std::string toString() const;
	State* copy()  const;

	bool is_perfect() const;
	bool is_near_perfect() const;
};

}
}
}


#endif /* STATE_H_ */

