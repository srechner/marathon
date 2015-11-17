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
	SwitchBipartiteFast(std::string line);

protected:

	virtual void computeNeighbours(const DenseBipartiteGraph& s,
			boost::unordered_map<DenseBipartiteGraph, Rational>& neighbors) const;
};

}

}

}

#endif /* CHAIN_SWAP_BIPARTITE_FAST_H_ */
