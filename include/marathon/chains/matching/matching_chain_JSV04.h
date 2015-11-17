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

private:
	uint num_perfect_matching;
	uint *num_near_perfect_matching;

public:

	JerrumSinclairVigoda04(std::string line);

	~JerrumSinclairVigoda04();

	void computeNeighbours(const BipartiteMatching& s,
			boost::unordered_map<BipartiteMatching, Rational>& neighbors) const;

	void computeWeights(std::vector<BipartiteMatching>& states);

	Rational getWeight(int i) const;
};

}

}

}

#endif /* JS89_H_ */
