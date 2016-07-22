/*
 * LindaChain.h
 *
 *  Created on: Jul 21, 2016
 *      Author: rechner
 */

#ifndef INCLUDE_MARATHON_CHAIN_BIPGRAPH_LINDACHAIN_H_
#define INCLUDE_MARATHON_CHAIN_BIPGRAPH_LINDACHAIN_H_

#include "../../MarkovChain.h"
#include <vector>
#include <stack>

#include "BinaryMatrix.h"

namespace marathon {
namespace chain {
namespace bipgraph {

/**
 * Implements the Markov chain developed by Linda Strowik.
 */
class LindaChain: public MarkovChain {

	friend class KannanPath;

private:

	// upper and lower row and column sums
	std::vector<int> u_lower, u_upper, v_lower, v_upper;
	int nrow, ncol;		// the number of rows and columns

public:

	/**
	 * Parse string s that encodes four vectors of integers.
	 */
	LindaChain(std::string& s);

	/**
	 * Create the MarkovChain object w.r.t. the given vectors.
	 */
	LindaChain(std::vector<int>& u_lower, std::vector<int>& u_upper,
			std::vector<int>& v_lower, std::vector<int>& v_upper);

	/**
	 * Default Destructor.
	 */
	virtual ~LindaChain();

	/**
	 * Create a 0-1-matrix out whose row and column sums agree with
	 * the constraints given by u_lower, u_upper, v_lower, v_upper.
	 */
	virtual State* computeArbitraryState() const;

	/**
	 * Randomize the state s by applying t transition steps.
	 */
	virtual void randomize(State* s, const uint32_t t = 1) const;
};

}
}
}

#endif /* INCLUDE_MARATHON_CHAIN_BIPGRAPH_LINDACHAIN_H_ */
