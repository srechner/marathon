/*
 * markovchain.h
 *
 *  Created on: Mar 3, 2016
 *      Author: rechner
 */

#ifndef MARATHON_MARKOVCHAIN_H_
#define MARATHON_MARKOVCHAIN_H_

#include <climits>
#include <string>
#include <list>
#include <set>
#include <random>

#include "Rational.h"
#include "State.h"
#include "StateGraph.h"

namespace marathon {

class StateGraph;

/**
 * Virtual Markov chain base class.
 */
class MarkovChain {

protected:

public:

	/*
	 * Standard Destructor
	 */
	virtual ~MarkovChain();


	/**
	 * Return a human readable name (identifier) of the Markov chain.
	 */
	std::string getName() const;

	/* Virtual Function that has to be implemented by subclasses */

	/**
	 * Computes an arbitrary state and store it the state object s.
	 * @return A pointer to a state object or nullptr if state space is empty.
	 */
	virtual State * computeArbitraryState() const;

	/**
	 *  Compute the set of adjacent states of s with corresponding proposal probability.
	 *  @param s A pointer to the state for which its neighbours are to be computed.
	 *  @param neighbors A vector with pointers to adjacent state objects that and their proposal
	 *  probabilities.
	 */
	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const ;

	/**
	 * Return the loop probability of state s.
	 * @param s: The state whose loop probability is to determined.
	 */
	virtual rational loopProbability(const State* s) const;

	/**
	 *  Computes weights for each state.
	 *  @param states The Vector of states.
	 *  @param weights The Vector of weights. After calling the method, this vector must
	 *  have the same size as states and is filled with rationals.
	 */
	virtual void computeWeights(const std::vector<const State*>& states,
			std::vector<rational>& weights);

	/**
	 * Apply a random transition to the state. Used to simulate a random walk.
	 * @param s A pointer to a state, which is randomly modified by the method.
	 * @param t The number of steps in the walk.
	 */
	virtual void randomize(State* s, const uint32_t t=1) const;

};

}

#endif /* INCLUDE_MARATHON_MARKOVCHAIN_H_ */
