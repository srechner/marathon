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

// project includes
#include "../../common/markov_chain.hpp"
#include "sparse_bipartite_graph.h"
#include "bipartite_matching.h"

// sampling-core includes

namespace marathon {

namespace chain {

namespace matching {

class Broder86: public MarkovChain<BipartiteMatching> {

protected:
	SparseBipartiteGraph g;

public:
	Broder86(std::string line);
	virtual ~Broder86();

	virtual bool computeArbitraryState(BipartiteMatching& s) const;
	virtual void computeNeighbours(const BipartiteMatching& s,
			boost::unordered_map<BipartiteMatching, Rational>& neighbors) const;

	void canonicalPath(int s, int t, std::list<int>& path) const;

};

}

}

}
#endif /* JS89_H_ */

