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
 * Adjacency Array Representation of a State Graph.
 *
 */
class StateGraph {

	/***************************************************
	 * declare friends
	 ***************************************************/

	template<typename T>
	friend T cpu::secondLargestEigenvalue(const StateGraph*);

	template<typename T>
	friend class cpu::DenseTransitionMatrix;

public:

	/* Variables */

	// number of states
	size_t numStates;

	// edge array representation of state graph
	std::vector<Transition> arcs;
	std::vector<uint> outgoing_arcs;

	// stationary distribution
	std::vector<Rational> stationary_distribution;

public:

	/* Constructor and Destructor */
	StateGraph();
	virtual ~StateGraph();

	/* Getter Methods */

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
	Rational getTransitionProbability(int u, int v) const;

	/**
	 * Set P(u,v) to p
	 */
	void setTransitionProbability(int u, int v, Rational p);

	/**
	 * Returns the stationary probability of state[i].
	 */
	Rational getStationary(int i) const;

	/**
	 * Returns the smallest stationary probability of all states
	 */
	Rational getMinimalStationary() const;

	/**
	 * Returns the index of the first outgoing transition of state[v] or getNumTransition(),
	 * if v does not have an outgoing arc.
	 **/
	int getIndexOfFirstTransition(int v) const;

	/**
	 * Returns the index of the last outgoing transition of state[v] or getNumTransition(),
	 * if v does not have an outgoing arc.
	 */
	int getIndexOfLastTransition(int v) const;

	/**
	 * Returns the number of adjacent states of state[v]
	 */
	int getNumTransitions(int v) const;

	/**
	 * Returns the Transition with given index
	 */
	Transition getTransition(int index) const;

	/**
	 * Returns the weight of state[i]
	 */
	virtual Rational getWeight(int i) const;

	/**
	 * Constructs a path between state[u] and state[v]
	 */
	virtual void canonicalPath(int u, int v, std::list<int>& path) const;

	virtual void constructStateGraph(bool verbose = false) = 0;

	virtual void printStates() const = 0;
};

}

#endif /* STATEGRAPH_H */
