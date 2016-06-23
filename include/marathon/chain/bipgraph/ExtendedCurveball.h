/*
 * curveball.h
 *
 *  Created on: Aug 24, 2015
 *      Author: rechner
 */

#ifndef CURVEBALL2_H_
#define CURVEBALL2_H_

#include "../bipgraph/Curveball.h"

namespace marathon {
namespace chain {
namespace bipgraph {

class ExtendedCurveball: public Curveball {

public:
	ExtendedCurveball(const std::string& line);

protected:
	void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

	/**
	 * Backtrack all two-partitions of the set {1..n} and store them in list.
	 */
	void backtrackTwoPartitions(int n,
			std::vector<std::vector<std::pair<int, int> > >& partitions) const;

	// just another recursive backtracking procedure
	void backtrackTwoPartitionsRecursive(int n, int* a, int* b, int k, int pos,
			std::vector<std::vector<std::pair<int, int> > >& partitions) const;

};

}

}

}

#endif /* CURVEBALL2_H_ */
