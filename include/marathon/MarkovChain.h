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

	std::string instance;		// the instance representation


public:

	/**
	 * Create A Markov Chain Object for the input s.
	 * @param s input string of the markov chain.
	 * false if is to be constructed by Omega^2 procedure.
	 */
	MarkovChain(const std::string& s, int seed = 0);

	/*
	 * Standard Destructor
	 */
	virtual ~MarkovChain();

	/**
	 * @return A reference to the string instance.
	 */
	const std::string& getInstance() const;

	/**
	 * Return a human readable name (identifier) of the Markov chain.
	 */
	std::string getName() const;

	/* Virtual Function that has to be implemented by subclasses */

	/**
	 * Computes an arbitrary state and store it the state object s.
	 * @return A pointer to a state object or nullptr if state space is empty.
	 */
	virtual State * computeArbitraryState() = 0;

	/**
	 *  Compute the set of adjacent states of s with corresponding proposal probability.
	 *  @param s A pointer to the state for which its neighbours are to be computed.
	 *  @param neighbors A vector with pointers to adjacent state objects that and their proposal
	 *  probabilities.
	 */
	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const = 0;

	/**
	 *  Computes weights for each state.
	 *  @param states The Vector of states.
	 *  @param weights The Vector of weights. After calling the method, this vector must
	 *  have the same size as states and is filled with rationals.
	 */
	virtual void computeWeights(const std::vector<const State*>& states,
			std::vector<rational>& weights);

	/**
	 * Apply a random walk and return the current state at the end of the walk.
	 *
	 * @param t The number of steps in the walk.
	 *
	 * @return A state, randomly selected from the probability distribution p^(t)_s,
	 *         where s is the state that is constructed via the computeArbitraryState
	 *         method.
	 */
	State* randomWalk(const int t) {
		marathon::State* s = computeArbitraryState();

		for (int i = 0; i < t; i++)
			randomize(s);

		return s;
	}

	/**
		 * Apply a random transition to the state. Used to simulate a random walk.
		 * @param s A pointer to a state, which is randomly modified by the method.
		 */
		virtual void randomize(State* s) const;

};

}

#endif /* INCLUDE_MARATHON_MARKOVCHAIN_H_ */
