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
class MarkovChain;


/**
 * State Graph representation. A State Graph is
 * a directed, weighted graph that represents a instance
 * of a Markov Chain for a certain input instance.
 */
class StateGraph {

protected:

	MarkovChain* mc;			// pointer to markov chain object

	/* Vertices and its attributes */
	std::vector<const State*> states;			// The set of states.
	std::vector<rational> weights;				// A weight for each state
	std::unordered_map<State*, int, State::Hash, State::Equal> indices;	// State -> Index

	/* The transition set and views on it. */
	std::vector<Transition*> arcs;
	std::vector<std::vector<Transition*>> outArcs;
	std::vector<std::vector<Transition*>> inArcs;

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
			const bool verbose);

public:

	/**
	 * Standard Constructor. Creates an empty State Graph.
	 * @param mc A pointer to the Markov Chain Object that defines transition rules, etc.
	 * @param limit A limit on the number of states of the graph. The graph can later on be expanded by the expand() method.
	 */
	StateGraph(MarkovChain* mc, const int limit = INT_MAX);

	/**
	 * Standard Destructor. Remove everything.
	 */
	virtual ~StateGraph();

	/**
	 * Expands an existing state graph to a given maximum of states.
	 * @param limit The maximal number of states after the expansion
	 * @param verbose Enables or disables additional debug output
	 * @return the number of states that has been added during the expansion
	 */
	void expand(const int limit = INT_MAX, const bool verbose = false);

	/**
	 * Return a pointer to the corresponding Markov Chain Object.
	 */
	MarkovChain* getMarkovChain() const;

	/**
	 * Add a new transition to the state graph that represents a loop.
	 */
	//int addLoopArc(const int u, const rational& p);

	/**
	 * Adds a transition arc to the graph.
	 * Precondition: The state graph does not already contain an arc between state u and state v.
	 * @return Returns the index of the new transition.
	 */
	int addArc(const int u, const int v, const rational& p);

	/**
	 * Adds a transition arc to the graph.
	 * Precondition: The state graph does not already contain an arc between state t.u and state t.v.
	 * @return Returns the index of the new transition.
	 */
	int addArc(Transition* t);

	/**
	 * Return a pointer to the arc that connects u with v or nullptr, if no such arc exists.
	 */
	Transition* getArc(int u, int v) const;

	/**
	 * Returns the number of states of the state graph
	 */
	size_t getNumStates() const;

	/**
	 * Returns the number of Transitions/Arcs of the state graph
	 */
	size_t getNumTransitions() const;

	/**
	 * Returns the transition probability P_uv for going from states[u] to states[v]
	 */
	rational getTransitionProbability(int u, int v) const;

	/**
	 * Set P(u,v) to p
	 */
	void setTransitionProbability(int u, int v, rational p);

	/**
	 * Increases P(u,v) by an amount of p.
	 */
	void addTransitionProbability(int u, int v, rational p);

	/**
	 * Sets the weight of state[i] to p.
	 */
	void setWeight(const int i, const rational p);

	/**
	 * Return the weight of state i.
	 */
	rational getWeight(const int i) const;

	/**
	 * Return the minimal weight of a state.
	 */
	rational getMinWeight() const;

	/**
	 * Return the sum of all weights.
	 */
	rational getZ() const;

	/**
	 * Return a vector of weights for each state.
	 */
	const std::vector<rational>& getWeights() const;

	/**
	 * Returns a reference to the outgoing arcs of state v.
	 */
	const std::vector<Transition*>& getOutArcs(int v) const;

	/**
	 * Returns a reference to the ingoing arcs of state v.
	 */
	const std::vector<Transition*>& getInArcs(int v) const;

	/**
	 * Returns a reference to the vector of all arcs in the state graph.
	 */
	const std::vector<Transition*>& getArcs() const;

	/**
	 * Return a pointer to arc with index i.
	 * @param i The index of the arc.
	 * @return A pointer to the i'th transition.
	 */
	Transition* getArc(const int i) const;

	/**
	 * Returns the number of adjacent states of state[v]
	 */
	int getNumOutArcs(int v) const;

	/**
	 * Removes all States and Transitions and re-initializes the state graph.
	 */
	virtual void clear();

	/**
	 * Add a new State to the state graph.
	 * @param s The State to insert.
	 * @return The index of the state after insertion.
	 */
	int addState(State* s);

	/**
	 * Returns a reference to the State with index i.
	 */
	const State* getState(int i) const;

	/**
	 * Returns a reference to a vector of States.
	 */
	const std::vector<const State*>& getStates() const;

	/**
	 * Returns the index of a state or -1 if the state graph does not contain this state.
	 */
	int indexOf(const State * s) const;
};

std::ostream& operator<<(std::ostream& out, const StateGraph& sg);

}

#endif /* STATEGRAPH_H */
