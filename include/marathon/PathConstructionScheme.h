/*
 * PathConstructionScheme.h
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#ifndef INCLUDE_MARATHON_PATHCONSTRUCTIONSCHEME_H_
#define INCLUDE_MARATHON_PATHCONSTRUCTIONSCHEME_H_

#include "StateGraph.h"

namespace marathon {

/**
 * A virtual base class for construction schemes of Canonical Paths.
 */
class PathConstructionScheme {

public:

	virtual ~PathConstructionScheme() {

	}

	/**
	 * Construct a path between states s and t in Graph sg.
	 * @param sg A pointer to a state graph object at which the path is embedded.
	 * @param s The index of the paths start state.
	 * @param t The index of the paths final state.
	 * @param path A list of state indices that represent the path.
	 */
	virtual void construct(const StateGraph* sg, const int s, const int t,
			std::list<int>& path) const = 0;

};

}

#endif /* INCLUDE_MARATHON_PATHCONSTRUCTIONSCHEME_H_ */
