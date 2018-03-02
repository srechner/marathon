/*
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

#include "marathon/integer.h"
#include "basic_random.h"
#include "rational.h"
#include "state.h"

namespace marathon {

	/**
	 * Virtual Markov chain base class.
	 */
	class MarkovChain {

	protected:

		// a pointer to the current state object during randomization
		State *currentState = nullptr;

		// random generator used for randomization
		marathon::BasicRandom rg;

	public:

		/**
		 * Standard Constructor
		 */
		MarkovChain() = default;

		/**
         * Declare a Markov chain with initial state s.
         * @param s Initial state.
         */
		MarkovChain(const State *s) {
			currentState = s->copy();
		}

		/*
		 * Standard Destructor
		 */
		virtual ~MarkovChain() {
			delete currentState;
		}

		/**
		 * Randomize the current state of the Markov chain.
		 * @param steps Number of steps.
		 * @return The current state after randomization.
		 */
		virtual const State *randomize(int steps) {
			for (int i = 0; i < steps; i++)
                step();
			return getCurrentState();
		}

		/**
		 * Return the current state of the Markov chain.
		 * @return
		 */
		virtual const State *getCurrentState() const {
			return currentState;
		}

		/**
		 * Return weight of state s.
		 * @param s State of the Markov chain.
		 * @return w(s).
		 */
		virtual Rational getWeight(const State* s) const {
			return Rational(1);
		}


		/**
		 * Set the current state of the Markov chain.
		 * A copy of state s is used from now on as current State.
		 * @param s State pointer that is used to create the current state.
		 */
		virtual
		void setCurrentState(const State *s) {
			currentState = s->copy();
		}

		/**
		 *  Generate the set of states that are adjacent to s and call the function f for each state.
		 *  @param s A pointer to the state for which its neighbours are to be computed.
		 *  @param f A std::function object that represents a function Each adjacent state is to be processed by the function f.
		 */
		virtual void adjacentStates(
				const State *s,
				const std::function<void(const State *, const Rational &)> &f
		) const {
			throw std::runtime_error("marathon::Exception: virtual method adjacentStates is not implemented!");
		}

		/**
		 * Apply one step of the Markov chain to the current state.
		 * @return The current state after randomization.
		 */
		virtual void step() {
			throw std::runtime_error("marathon::Exception: virtual method step() is not implemented!");
		}

		/**
		 * Create a copy of this MarkovChain.
		 * @return
		 */
		virtual MarkovChain *copy() const {
			throw std::runtime_error("marathon::Exception: virtual method copy() is not implemented!");
		}
	};
}

#endif /* INCLUDE_MARATHON_MARKOVCHAIN_H_ */
