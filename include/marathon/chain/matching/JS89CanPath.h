/*
 * JS89CanPath.h
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#ifndef INCLUDE_MARATHON_CHAINS_MATCHING_JS89CANPATH_H_
#define INCLUDE_MARATHON_CHAINS_MATCHING_JS89CANPATH_H_

#include "../../PathConstructionScheme.h"
#include "../../MarkovChain.h"

namespace marathon {
namespace chain {
namespace matching {

class JS89Path: public PathConstructionScheme {

	/**
	 * Path construction Scheme for canonical paths from Jerrum and Sinclair 1989.
	 * @param sg Pointer to a state graph object.
	 * @param s Index of a state (Start).
	 * @param t Index of a state (End).
	 * @param path List of State indices that corresponds to the state.
	 */
	virtual void construct(const StateGraph* sg, const int s, const int t,
			std::list<int>& path) const;

};

}
}
}

#endif /* INCLUDE_MARATHON_CHAINS_MATCHING_JS89CANPATH_H_ */
