/*
 * JS89.h
 *
 *  Created on: Nov 18, 2014
 *      Author: rechner
 */

#ifndef JS89_CHAIN_H_
#define JS89_CHAIN_H_

#include <queue>

#include "../../MarkovChain.h"
#include "BipartiteMatching.h"
#include "SparseBipartiteGraph.h"

namespace marathon {
namespace chain {
namespace matching {

class Broder86: public MarkovChain {

	friend class JS89Path;

protected:

	SparseBipartiteGraph *g = nullptr;

	void parseInstance(const std::string& inst);

public:

	Broder86(const std::string& instance);
	~Broder86();

	virtual State* computeArbitraryState();
	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

};

}

}

}
#endif /* JS89_H_ */

