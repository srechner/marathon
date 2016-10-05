/*
 * MarkovChain.h
 *
 * Created on: Mar 3, 2016
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2016, Steffen Rechner
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef MARATHON_MARKOVCHAIN_H_
#define MARATHON_MARKOVCHAIN_H_

#include <climits>
#include <string>
#include <list>
#include <set>
#include <random>

#include "Combinatorics.h"
#include "Random.h"
#include "Rational.h"
#include "State.h"

namespace marathon {

class StateGraph;

/**
 * Virtual Markov chain base class.
 */
class MarkovChain {

protected:

public:

	MarkovChain() {
		Combinatorics::init();
		Random::init();
	}

	/*
	 * Standard Destructor
	 */
	virtual ~MarkovChain() {
		Combinatorics::cleanup();
		Random::cleanup();
	}


	/**
	 * Return a human readable name (identifier) of the Markov chain.
	 */
	std::string getName() const {
		return "[unknown]";
	}

	/* Virtual Function that has to be implemented by subclasses */

	/**
	 * Computes an arbitrary state and store it the state object s.
	 * @return A pointer to a state object or nullptr if state space is empty.
	 */
	virtual State * computeArbitraryState() const {
		// dummy implementation
		std::cout << "marathon::Exception: computeArbitraryState is not implemented!"
		          << std::endl;
		return nullptr;
	}

	/**
	 *  Compute the set of adjacent states of s with corresponding proposal probability.
	 *  @param s A pointer to the state for which its neighbours are to be computed.
	 *  @param neighbors A vector with pointers to adjacent state objects that and their proposal
	 *  probabilities.
	 */
	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const {
		// dummy implementation
		std::cout << "marathon::Exception: computeNeighbouts is not implemented!"
		          << std::endl;
	}

	/**
	 * Return the loop probability of state s.
	 * @param s: The state whose loop probability is to determined.
	 */
	virtual rational loopProbability(const State* s) const {
		rational res(0);
		std::vector<std::pair<State*, rational>> N;
		this->computeNeighbours(s, N);
		for (auto x : N) {
			if (x.first->compare(s) == 0) {
				res += x.second;
			}
			delete x.first;
		}
		return res;
	}

	/**
	 *  Computes weights for each state.
	 *  @param states The Vector of states.
	 *  @param weights The Vector of weights. After calling the method, this vector must
	 *  have the same size as states and is filled with rationals.
	 */
	virtual void computeWeights(const std::vector<const State*>& states,
			std::vector<rational>& weights) {
		weights.clear();
		weights.reserve(states.size());
		for (int i = 0; i < states.size(); i++)
			weights.push_back(1);
	}

	/**
	 * Apply a random transition to the state. Used to simulate a random walk.
	 * @param s A pointer to a state, which is randomly modified by the method.
	 * @param t The number of steps in the walk.
	 */
	virtual void randomize(State* s, const uint32_t t=1) const {
		// dummy implementation
		std::cout << "marathon::Exception: randomize is not implemented!"
		          << std::endl;
	}

};

}

#endif /* INCLUDE_MARATHON_MARKOVCHAIN_H_ */
