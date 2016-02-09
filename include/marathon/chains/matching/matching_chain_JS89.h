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

	SparseBipartiteGraph *g = nullptr;

	void parseInstance(const std::string& inst);

public:

	Broder86(const std::string& instance);
	~Broder86();

	virtual bool computeArbitraryState(BipartiteMatching& s);
	virtual void computeNeighbours(const BipartiteMatching& s,
			std::unordered_map<BipartiteMatching, rational>& neighbors) const;

	void constructPath(const StateGraph* sg, int s, int t,
			std::list<int>& path) const;

};

}

}

}
#endif /* JS89_H_ */

