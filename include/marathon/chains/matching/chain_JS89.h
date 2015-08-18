/*
 * JS89.h
 *
 *  Created on: Nov 18, 2014
 *      Author: rechner
 */

#ifndef JS89_CHAIN_H_
#define JS89_CHAIN_H_

// STL includes
#include <queue>

#include "../../markov_chain.hpp"
// project includes
#include "bipartite_matching.h"
#include "sparse_bipartite_graph.h"

// sampling-core includes

namespace marathon {

namespace chain {

namespace matching {

class MatchingChain89: public MarkovChain<BipartiteMatching> {

protected:
	SparseBipartiteGraph g;

public:
	MatchingChain89(std::string line);
	virtual ~MatchingChain89();

	virtual bool computeArbitraryState(BipartiteMatching& s) const;
	virtual void computeNeighbors(const BipartiteMatching& s,
			std::vector<BipartiteMatching>& neighbors) const;

	void canonicalPath(int s, int t, std::list<int>& path) const;

};

}

}

}
#endif /* JS89_H_ */

