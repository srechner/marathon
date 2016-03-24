/*
 * Chain.h
 *
 *  Created on: Nov 24, 2014
 *      Author: rechner
 */

#ifndef _BIP_DEG_CHAIN_H_
#define _BIP_DEG_CHAIN_H_

#include "../../MarkovChain.h"
#include <vector>
#include <stack>

#include "../bipgraph/DenseBipartiteGraph.h"

namespace marathon {
namespace chain {
namespace bipgraph {

/**
 * Implements the Markov chain defined by Kannan et al.
 */
class SwitchChain: public MarkovChain {

	friend class KannanPath;

protected:

	// variables and methods for general purpose
	std::vector<int> u;
	std::vector<int> v;
	int sum;

	/**
	 * Instances have the form "2,2,2;1,2,1,2".
	 * The semicolon separates two degree sequences of both bipartition sets.
	 */
	virtual void parseInstance(const std::string& line);

public:

	SwitchChain(const std::string& inst);
	virtual ~SwitchChain();

	virtual State* computeArbitraryState();
	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;
};

}
}
}

#endif /* CHAIN_H_ */
