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

protected:

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
	size_t getNumStates() const;
	size_t getNumArcs() const;
	Rational getTransitionProbability(int u, int v) const;
	Rational getStationary(int i) const;
	Rational getMinimalStationary() const;
	virtual Rational getWeight(int i) const;

	virtual void canonicalPath(int u, int v, std::list<int>& path) const;
	virtual void constructStatespace() = 0;
};

}

#endif /* STATEGRAPH_H */
