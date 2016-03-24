/*
 * JS89.h
 *
 *  Created on: Nov 18, 2014
 *      Author: rechner
 */

#ifndef JSV04_CHAIN_H_
#define JSV04_CHAIN_H_

#include <queue>

#include "Broder86.h"

namespace marathon {
namespace chain {
namespace matching {

class JerrumSinclairVigoda04: public Broder86 {

protected:

	uint num_perfect_matching;
	uint *num_near_perfect_matching;

	rational getWeight(const State* s) const;

public:

	JerrumSinclairVigoda04(const std::string& input);

	~JerrumSinclairVigoda04();

	void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

	void computeWeights(const std::vector<const State*>& states,
			std::vector<rational>& weights);
};

}
}
}

#endif /* JS89_H_ */
