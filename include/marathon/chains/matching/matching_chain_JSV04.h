/*
 * JS89.h
 *
 *  Created on: Nov 18, 2014
 *      Author: rechner
 */

#ifndef JSV04_CHAIN_H_
#define JSV04_CHAIN_H_

#include <queue>
#include "matching_chain_JS89.h"

namespace marathon {

namespace chain {

namespace matching {

class JerrumSinclairVigoda04: public Broder86 {

protected:

	uint num_perfect_matching;
	uint *num_near_perfect_matching;

	rational getWeight(const BipartiteMatching& s) const;

public:

	JerrumSinclairVigoda04(const std::string& input);

	~JerrumSinclairVigoda04();

	void computeNeighbours(const BipartiteMatching& s,
			std::unordered_map<BipartiteMatching, rational>& neighbors) const;

	void computeWeights(std::vector<BipartiteMatching>& states, std::vector<rational>& weights);
};

}

}

}

#endif /* JS89_H_ */
