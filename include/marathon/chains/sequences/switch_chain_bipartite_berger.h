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

class SwitchBipartiteFast: public SwitchBipartite {

public:
	SwitchBipartiteFast(const std::string& input);

protected:

	virtual void computeNeighbours(const DenseBipartiteGraph& s,
			std::unordered_map<DenseBipartiteGraph, rational>& neighbors) const;
};

}

}

}

#endif /* CHAIN_SWAP_BIPARTITE_FAST_H_ */
