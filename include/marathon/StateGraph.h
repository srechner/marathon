/*
 * StateGraph.h
 *
 * Created on: Aug 4, 2015
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

#ifndef STATEGRAPH_H
#define STATEGRAPH_H

#include "State.h"
#include "MarkovChain.h"

#include <climits>
#include <list>
#include <string>
#include <unordered_map>
#include <map>
#include "Transition.h"

namespace marathon {

	class StateGraph;

	/**
	 * State Graph representation. A State Graph is
	 * a directed, weighted graph that represents a instance
	 * of a Markov Chain for a certain input instance.
	 */
	class StateGraph {

	protected:

		MarkovChain *mc;            // pointer to markov chain object

		/* Vertices and its attributes */
		std::vector<const State *> states;            // The set of states.
		std::vector<rational> weights;                // A weight for each state
		std::unordered_map<State *, int, State::Hash, State::Equal> indices;    // State -> Index

		/* The transition set and views on it. */
		std::vector<Transition *> arcs;
		std::vector<std::vector<Transition *>> outArcs;
		std::vector<std::vector<Transition *>> inArcs;

		/* Variables used for State Graph Construction */
		int nextIndex = 0;
		std::set<int> reexpand;

		/* Helper Functions */

		/**
		 * This is a private method that is called during state graph expansion.
		 * It computes all neighbouring states of state s and insert them into the state graph repectively
		 * into the leftover structures that store the states and arcs for next expandStateGraph().
		 * @param i The index of the state that is to be expanded.
		 * @param limit The maximal number of states.
		 * @param lastStop The size of the state graph when this expansion has been triggered.
		 * @param verbose If true, additional debug information is printed.
		 * @param True, if all adjacent states could be inserted in the state graph.
		 */
		void expandState(const int i, const int limit, const int lastStop,
		                 const bool verbose) {
			if (verbose)
				std::cout << "expand state " << i << std::endl;

			const State *s = (const State *) getState(i);

			// compute states that are neighboured to s
			std::vector<std::pair<State *, rational>> neighbours;
			mc->computeNeighbours(s, neighbours);

			/**********************************************************************
			 *  eliminate duplicate states by sorting and compressing
			 **********************************************************************/
			std::sort(neighbours.begin(), neighbours.end(),
			          [](const std::pair<State *, rational> &s1, const std::pair<State *, rational> &s2) -> bool {
				          return s1.first->compare(s2.first) < 0;
			          });

			// compress
			std::vector<std::pair<State *, rational>> tmp;
			tmp.reserve(neighbours.size());
			auto it = neighbours.begin();
			tmp.push_back(*it);
			for (++it; it != neighbours.end(); ++it) {
				// if duplicate state
				if (it->first->compare(tmp.back().first) == 0) {
					tmp.back().second += it->second;
					delete it->first;
				} else {
					tmp.push_back(*it);
				}
			}
			//neighbours.swap(tmp);
			neighbours = std::move(tmp);
			tmp.clear();

			bool complete = true;

			// sum of proposal probabilites of adjacent states
			rational sum(0);

			// for all adjacent states
			for (auto it = neighbours.begin(); it != neighbours.end(); ++it) {

				// neighbour state s2 with proposal prob. kappa(s,s2)
				State *s2 = it->first;
				const rational &kappa = it->second;
				sum += kappa;

				// Look if s2 already known
				int j = indexOf(s2);    // j is the index of state s2

				// if s2 is already known
				if (j != -1) {

					if (verbose) {
						std::cout << " " << j << ": " << s2->toString() << " " << kappa
						          << " already known state" << std::endl;
					}

					// if arc (i,j) has not already included last round
					if (j >= lastStop) {
						// add a transition arc to the state graph
						addArc(i, j, kappa);
					}

					delete s2;
				}
					// if s2 is not seen so far and the maximal number of states is not reached
				else if (getNumStates() < limit) {

					// add s2 to the vector of states
					j = addState(s2);

					if (verbose) {
						std::cout << " " << j << ": " << s2->toString() << " " << kappa
						          << " new state" << std::endl;
					}

					// add a transition arc to the state graph
					addArc(i, j, kappa);
				}
					// if s2 is not seen so far and the maximal number of states is reached
				else {
					reexpand.insert(i);
					complete = false;
					delete s2;
				}
			}

			// do the proposal probabilites sum up to one?
			if (sum != rational(1)) {
				std::cerr
						<< "marathon::StateGraph::expandState: Error: Sum of proposal probabilities of state "
						<< s << " is not one!" << std::endl;
			}
		}

	public:

		/**
		 * Standard Constructor. Creates an empty State Graph.
		 * @param mc A pointer to the Markov Chain Object that defines transition rules, etc.
		 * @param limit A limit on the number of states of the graph. The graph can later on be expanded by the expand() method.
		 */
		StateGraph(MarkovChain *mc, const int limit = INT_MAX) : mc(mc) {
			expand(limit);
		}

		/**
		 * Standard Destructor. Remove everything.
		 */
		virtual ~StateGraph() {
			// delete all transitions
			for (int i = 0; i < arcs.size(); i++)
				delete arcs[i];
		}

		/**
		 * Expands an existing state graph to a given maximum of states.
		 * @param limit The maximal number of states after the expansion
		 * @param verbose Enables or disables additional debug output
		 * @return the number of states that has been added during the expansion
		 */
		void expand(const int limit = INT_MAX, const bool verbose = false) {
			const int sizeLast = getNumStates();

			// nothing to do?
			if (limit <= sizeLast)
				return;

			// if state graph is empty
			if (sizeLast == 0) {

				// Start with arbitrary State
				State *s1 = mc->computeArbitraryState();
				if (s1 == nullptr) {

					if (verbose) {
						// print warning if statespace is empty
						std::cerr << "Warning! Empty state!" << std::endl;
					}
					return;
				}

				if (verbose) {
					std::cout << "Start state is " << s1->toString() << std::endl;
				}

				// add initial state
				addState(s1);
			}

			/***********************************************************
			 * Expand states that have not been expanded completely
			 **********************************************************/

			if (verbose)
				std::cout << "re-expand states from last round" << std::endl;

			// gather the indices that have to re-expanded
			std::vector<int> tmp(reexpand.begin(), reexpand.end());
			reexpand.clear();

			// for all states that are not yet completely expanded
			for (int y : tmp) {
				// re-expand state with index y
				expandState(y, limit, sizeLast, verbose);
			}

			/***********************************************************
			 * Explore further regions of the state graph.
			 **********************************************************/

			if (verbose)
				std::cout << "explore further regions" << std::endl;

			for (int y = sizeLast; y < getNumStates(); y++) {
				// expand state with index y
				expandState(y, limit, 0, verbose);
			}

			/***********************************************************
			 * Compute weights for metropolis rule
			 **********************************************************/
			mc->computeWeights(states, weights);

			/***********************************************************
			 * Apply Metropolis Rule:
			 * P(i,j) = kappa(i,j) * min( w(j) / w(i) , 1)
			 ************************************************************/

			// gather the newly introduced loop probability
			std::vector<rational> generatedLoopProbability;
			generatedLoopProbability.reserve(getNumStates());

			// Transform Transition probabilities
			for (Transition *t : getArcs()) {

				const int i = t->u;
				const int j = t->v;
				const rational kappa = t->p;
				const rational w_i = weights[i];
				const rational w_j = weights[j];

				// apply metropolis rule
				if (w_j < w_i) {

					const rational metr = w_j / w_i;

					// P(i,j) = kappa(i,j) * min(1, w(j)/w(i))
					t->p = kappa * metr;

					// add remaining probability as loop probability
					const rational prob = (rational(1) - metr) * kappa;
					generatedLoopProbability[i] += prob;
				}
			}

			// add additional loop probability
			for (int i = 0; i < generatedLoopProbability.size(); i++) {
				addTransitionProbability(i, i, generatedLoopProbability[i]);
			}

			/***********************************************************
			 * Self-Check: Verifiy that chain is reversible
			 ***********************************************************/
			const rational Z = getZ();

			for (const Transition *t : getArcs()) {

				if (t->u >= t->v)
					continue;

				const rational stat_u = getWeight(t->u) / Z;
				const rational stat_v = getWeight(t->v) / Z;
				const rational puv = t->p;

				// find anti parallel transition arc
				for (const Transition *t2 : getInArcs(t->u)) {

					// t2 is anti parallel to t
					if (t2->v == t->v) {

						const rational pvu = t2->p;

						if (stat_u * puv != stat_v * pvu) {
							std::cerr << "Error! Chain is not reversible!" << std::endl;
							std::cerr << "P(" << t->u << "," << t->v << ")=" << puv
							          << std::endl;
							std::cerr << "P(" << t->v << "," << t->u << ")=" << pvu
							          << std::endl;
							std::cerr << "stat(" << t->u << ")=" << stat_u << std::endl;
							std::cerr << "stat(" << t->v << ")=" << stat_v << std::endl;
							std::cerr << "stat(" << t->u << ")*P(" << t->u << ","
							          << t->v << ")=" << stat_u * puv << std::endl;
							std::cerr << "stat(" << t->v << ")*P(" << t->v << ","
							          << t->u << ")=" << stat_v * pvu << std::endl;
						}

						assert(stat_u * puv == stat_v * pvu);
						break;
					}
				}
			}

			/*********************************************************************
			 * Print Information about states
			 ********************************************************************/

			if (verbose) {
				std::cout << "state size: " << getNumStates() << std::endl;
				for (int ii = 0; ii < getNumStates(); ii++) {
					const State *s = getState(ii);
					std::cout << ii << ": " << s->toString() << " " << weights[ii]
					          << std::endl;
				}

				std::cout << "transition size: " << getNumTransitions() << std::endl;
				for (const Transition *t : getArcs()) {
					std::cout << t->u << " " << t->v << " " << t->p << std::endl;
				}
			}
		}

		/**
		 * Return a pointer to the corresponding Markov Chain Object.
		 */
		MarkovChain *getMarkovChain() const {
			return mc;
		}

		/**
		 * Add a new transition to the state graph that represents a loop.
		 */
		//int addLoopArc(const int u, const rational& p);

		/**
		 * Adds a transition arc to the graph.
		 * Precondition: The state graph does not already contain an arc between state u and state v.
		 * @return Returns the index of the new transition.
		 */
		int addArc(const int u, const int v, const rational &p) {
			return addArc(new Transition(u, v, p));
		}

		/**
		 * Adds a transition arc to the graph.
		 * Precondition: The state graph does not already contain an arc between state t.u and state t.v.
		 * @return Returns the index of the new transition.
		 */
		int addArc(Transition *t) {
			// add the arc to the arc vector
			arcs.push_back(t);

			// add a pointer to the transition in u's and v's outarc/inarc array
			outArcs[t->u].push_back(t);
			inArcs[t->v].push_back(t);

			return arcs.size() - 1;
		}

		/**
		 * Return a pointer to the arc that connects u with v or nullptr, if no such arc exists.
		 */
		Transition *getArc(int u, int v) const {
			// search for transition (u,v)
			for (Transition *t : outArcs[u]) {
				assert(t->u == u);
				if (t->v == v) {
					return t;
				}
			}
			return nullptr;
		}

		/**
		 * Returns the number of states of the state graph
		 */
		size_t getNumStates() const {
			return states.size();
		}

		/**
		 * Returns the number of Transitions/Arcs of the state graph
		 */
		size_t getNumTransitions() const {
			return arcs.size();
		}

		/**
		 * Returns the transition probability P_uv for going from states[u] to states[v]
		 */
		rational getTransitionProbability(int u, int v) const {
			const Transition *t = getArc(u, v);
			if (t == nullptr)
				return 0;
			else
				return t->p;
		}

		/**
		 * Set P(u,v) to p
		 */
		void setTransitionProbability(int u, int v, rational p) {
			Transition *t = getArc(u, v);

			if (t != nullptr) {
				t->p = p;
			} else {
				// no transition found? add a new one
				addArc(new Transition(u, v, p));
			}
		}

		/**
		 * Increases P(u,v) by an amount of p.
		 */
		void addTransitionProbability(int u, int v, rational p) {
			Transition *t = getArc(u, v);

			if (t != nullptr) {
				t->p += p;
			} else {
				// no transition found? add a new one
				addArc(new Transition(u, v, p));
			}
		}

		/**
		 * Sets the weight of state[i] to p.
		 */
		void setWeight(const int i, const rational p) {
			weights[i] = p;
		}

		/**
		 * Return the weight of state i.
		 */
		rational getWeight(const int i) const {
			return weights[i];
		}

		/**
		 * Return the minimal weight of a state.
		 */
		rational getMinWeight() const {
			return *std::min(weights.begin(), weights.end());
		}

		/**
		 * Return the sum of all weights.
		 */
		rational getZ() const {
			rational Z = 0;
			for (auto w : weights)
				Z += w;
			return Z;
		}

		/**
		 * Return a vector of weights for each state.
		 */
		const std::vector<rational> &getWeights() const {
			return weights;
		}

		/**
		 * Returns a reference to the outgoing arcs of state v.
		 */
		const std::vector<Transition *> &getOutArcs(int v) const {
			return outArcs[v];
		}

		/**
		 * Returns a reference to the ingoing arcs of state v.
		 */
		const std::vector<Transition *> &getInArcs(int v) const {
			return inArcs[v];
		}

		/**
		 * Returns a reference to the vector of all arcs in the state graph.
		 */
		const std::vector<Transition *> &getArcs() const {
			return arcs;
		}

		/**
		 * Return a pointer to arc with index i.
		 * @param i The index of the arc.
		 * @return A pointer to the i'th transition.
		 */
		Transition *getArc(const int i) const {
			return (Transition *) &arcs[i];
		}

		/**
		 * Returns the number of adjacent states of state[v]
		 */
		int getNumOutArcs(int v) const {
			return this->getOutArcs(v).size();
		}

		/**
		 * Removes all States and Transitions and re-initializes the state graph.
		 */
		virtual void clear() {
			arcs.clear();
			inArcs.clear();
			outArcs.clear();
			weights.clear();
		}

		/**
		 * Add a new State to the state graph.
		 * @param s The State to insert.
		 * @return The index of the state after insertion.
		 */
		int addState(State *s) {
			// add state to the vector of states
			states.push_back(s);
			indices[s] = states.size() - 1;
			outArcs.push_back(std::vector<Transition *>());
			inArcs.push_back(std::vector<Transition *>());
			weights.push_back(1);
			return states.size() - 1;
		}

		/**
		 * Returns a reference to the State with index i.
		 */
		const State *getState(int i) const {
			return states[i];
		}

		/**
		 * Returns a reference to a vector of States.
		 */
		const std::vector<const State *> &getStates() const {
			return states;
		}

		/**
		 * Returns the index of a state or -1 if the state graph does not contain this state.
		 */
		int indexOf(const State *s) const {
			auto it = indices.find((State *) s);
			if (it != indices.end())
				return it->second;
			else
				return -1;
		}
	};

	inline
	std::ostream &operator<<(std::ostream &out, const StateGraph &sg) {

		out << "n " << sg.getNumStates() << " m " << sg.getNumTransitions() << "\n";
		for (int i = 0; i < sg.getNumStates(); i++) {
			for (Transition *t : sg.getOutArcs(i))
				out << t->u << " " << t->v << " " << t->p << "\n";
		}

		return out;
	}

}

#endif /* STATEGRAPH_H */
