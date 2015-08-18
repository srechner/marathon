/*
 * chain_swap_bipartite_fast.h
 *
 *  Created on: Jun 4, 2015
 *      Author: rechner
 */

#ifndef CHAIN_SWAP_BIPARTITE_FAST_H_
#define CHAIN_SWAP_BIPARTITE_FAST_H_

#include "switch_chain_bipartite.h"

namespace marathon {

namespace chain {

namespace sequence {

class SwapChainBipartiteFast: public SwapChainBipartite {
public:
	SwapChainBipartiteFast(std::string line);

protected:
	virtual void computeNeighbors(const DenseBipartiteGraph& s,
			std::vector<DenseBipartiteGraph>& neighbors) const;
};

}

}

}

#endif /* CHAIN_SWAP_BIPARTITE_FAST_H_ */
