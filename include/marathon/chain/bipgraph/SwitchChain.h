/*
 * Chain.h
 *
 *  Created on: Nov 24, 2014
 *      Author: rechner
 */

#ifndef _BIP_DEG_CHAIN_H_
#define _BIP_DEG_CHAIN_H_

#include "../../MarkovChain.h"
#include <vector>
#include <stack>

#include "BinaryMatrix.h"

namespace marathon {
namespace chain {
namespace bipgraph {

/**
 * Implements the Markov chain defined by Kannan et al.
 */
class SwitchChain: public MarkovChain {

	friend class KannanPath;

public:

	// variables and methods for general purpose
	std::vector<int> u;
	std::vector<int> v;
	int sum;

	/**
	 * Instances have the form "2,2,2;1,2,1,2".
	 * The semicolon separates two degree sequences of both bipartition sets.
	 */
	virtual void parseInstance(const std::string& line);

public:

	SwitchChain(const std::string& inst);
	SwitchChain(const std::vector<int>& u, const std::vector<int>& v);
	SwitchChain(int* const rowsum, int* const colsum, const int nrow, const int ncol);
	virtual ~SwitchChain();

	virtual State* computeArbitraryState() const;


	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

	/**
	 * Return the loop probability of State s.
	 */
	virtual rational loopProbability(const State* s) const;

	/**
	 * Randomize the state s by applying a single transition.
	 */
	virtual void randomize(State* s, const uint32_t t=1) const;
};

}
}
}

#endif /* CHAIN_H_ */
