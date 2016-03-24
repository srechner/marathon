/*
 * chain_swap_bipartite_fast.h
 *
 *  Created on: Jun 4, 2015
 *      Author: rechner
 */

#ifndef CHAIN_SWAP_BIPARTITE_FAST_H_
#define CHAIN_SWAP_BIPARTITE_FAST_H_

#include "../bipgraph/SwitchChain.h"

namespace marathon {
namespace chain {
namespace bipgraph {

class SwitchChainBerger: public SwitchChain {

public:
	SwitchChainBerger(const std::string& input);

protected:

	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

};

}
}
}

#endif /* CHAIN_SWAP_BIPARTITE_FAST_H_ */
