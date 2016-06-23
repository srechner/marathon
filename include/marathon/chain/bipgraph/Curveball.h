/*
 * curveball.h
 *
 *  Created on: Aug 24, 2015
 *      Author: rechner
 */

#ifndef CURVEBALL_H_
#define CURVEBALL_H_

#include "../bipgraph/SwitchChain.h"

namespace marathon {
namespace chain {
namespace bipgraph {

class Curveball: public SwitchChain {
public:
	Curveball(const std::string& line);
	Curveball(const std::vector<int>& u, const std::vector<int>& v);
	Curveball(int* const rowsum, int* const colsum, const int nrow, const int ncol);

	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

	void randomize(State* s, const uint32_t t=1) const;

	/**
	 * Return the loop probability of State s.
	 */
	virtual rational loopProbability(const State* s) const;

protected:

	/**
	 * Backtrack all subsets of size n from population. Each subset is stored in list.
	 */
	void backtrackSubsets(boost::dynamic_bitset<>& population, int n,
			const rational& p, const boost::dynamic_bitset<>& Bi,
			const boost::dynamic_bitset<>& Bj, const boost::dynamic_bitset<>& X,
			const DenseBipartiteGraph* s, const int i, const int j,
			std::vector<std::pair<State*, rational>>& myneighbours) const;

	// just another recursive backtracking procedure
	void backtrackSubsetsRecursive(const boost::dynamic_bitset<>& population,
			boost::dynamic_bitset<>& tmp, int pos, int m, int n,
			const boost::dynamic_bitset<>& Bi,
			const boost::dynamic_bitset<>& Bj, const boost::dynamic_bitset<>& X,
			const DenseBipartiteGraph* s, const int i, const int j,
			const rational& p,
			std::vector<std::pair<State*, rational>>& myneighbours) const;

};

}

}

}

#endif /* CURVEBALL_H_ */
