/*
 * markov_chain_base.h
 *
 *  Created on: Aug 4, 2015
 *      Author: rechner
 */

#ifndef STATEGRAPH_H
#define STATEGRAPH_H

#include "transition.h"
#include <list>
#include <string>
#include <unordered_map>

namespace marathon {

class StateGraph;

/*****************************************************************************
 * forward declaration of StateGraph's friend functions
 ****************************************************************************/

namespace cpu {

template<typename T>
T secondLargestEigenvalue(const StateGraph*);

template<typename T>
class DenseTransitionMatrix;

}

/**
 * A State Graph is a directed, weighted graph that represents a instance
 * of a Markov Chain for a certain input instance. A state graph
 * is implemented as a adjacency list.
 */
class StateGraph {

	/***************************************************
	 * declare friends
	 ***************************************************/

	template<typename T>
	friend T cpu::secondLargestEigenvalue(const StateGraph*);

	template<typename T>
	friend class cpu::DenseTransitionMatrix;

protected:

	int numStates;
	int numArcs;
	std::string instance;

	// The transition set and views on it.
	std::vector<Transition> arcs;
	std::vector<std::vector<Transition*>> outArcs;
	std::vector<std::vector<Transition*>> inArcs;

	// stationary distribution
	std::vector<rational> stationary_distribution;

public:

	/**
	 * Standard Constructor. Creates an empty State Graph.
	 */
	StateGraph();

	/**
	 * Standard Destructor. Remove everything.
	 */
	virtual ~StateGraph();

	/**
	 * Rerserve memory for n states and m arcs.
	 */
	virtual void resize(int n, int m);

	/*
	 * Returns a void pointer to the input instance that lies behind the state graph.
	 *
	 */
	const std::string& getInstance() const;

	/**
	 * Set the instance.
	 */
	void setInstance(const std::string& instance);

	/**
	 * Adds a transition arc to the graph.
	 * Precondition: The state graph does not already contain an arc between state u and state v.
	 */
	void addArc(int u, int v, rational p);

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
	 * Returns the stationary probability of state[i].
	 */
	rational getStationary(const int i) const;

	/**
	 * Sets the stationary probability of state[i] to p.
	 */
	void setStationary(const int i, const rational p);

	/**
	 * Returns the smallest stationary probability of all states
	 */
	rational getMinimalStationary() const;

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
	const std::vector<Transition>& getArcs() const;

	/**
	 * Returns the number of adjacent states of state[v]
	 */
	int getNumOutArcs(int v) const;

	/**
	 * Removes all States and Transitions and re-initializes the state graph.
	 */
	virtual void clear();
};

/**
 * This is the implementation of the abstract StateGraph class.
 */
template<typename State>
class _StateGraph: public StateGraph {

private:

	/* Variables */

	std::vector<State> states;					// The set of states.
	std::unordered_map<State, int> indices;		// State -> Index

public:

	/**
	 * Create Empty state graph.
	 */
	_StateGraph() :
			StateGraph() {
	}

	/**
	 * Standard Destructor. Remove everything.
	 */
	~_StateGraph() {

	}

	/**
	 * Rerserve memory for n states and m arcs.
	 */
	virtual void resize(int n, int m) {
		StateGraph::resize(n, m);
		states.resize(n);
	}

	/**
	 * Removes all States and Transitions and re-initializes the state graph.
	 */
	void clear() {
		StateGraph::clear();
		states.clear();
	}

	/* Non-inherited methods */

	/**
	 * Add a new State to the state graph.
	 * @param s The State to insert.
	 * @return The index of the state after insertion.
	 */
	int addState(const State& s) {
		// add state to the vector of states
		states[numStates] = s;
		indices[s] = numStates;
		numStates++;
		return numStates - 1;
	}

	/**
	 * Returns the State with index i.
	 */
	const State& getState(int i) const {
		return states[i];
	}

	/**
	 * Returns a reference to a vector of States.
	 */
	const std::vector<State>& getStates() const {
		return states;
	}

	/**
	 * Returns the index of a state or -1 if the state graph does not contain this state.
	 */
	int findState(const State& s) const {
		auto it = indices.find(s);
		if (it != indices.end())
			return it->second;
		else
			return -1;
	}

};

std::ostream& operator<<(std::ostream& out, const StateGraph& sg);

}

#endif /* STATEGRAPH_H */
