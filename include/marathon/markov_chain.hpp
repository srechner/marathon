/*
 * MarkovChain.h
 *
 *  Created on: Feb 11, 2015
 *      Author: steffen
 */

#ifndef MARKOVCHAIN_H_
#define MARKOVCHAIN_H_

// STL includes
#include <vector>
#include <list>
#include <queue>

// boost includes
#include "boost/unordered_map.hpp"

#include "state_graph.h"

//#define MY_DEBUG

namespace marathon {

template<class State>
class MarkovChain: public StateGraph {

protected:

	std::vector<State> states;					// vector of states
	boost::unordered_map<State, int> indices;	// State -> Index

public:

	MarkovChain() {

	}

	virtual ~MarkovChain() {

	}

	/* Virtual Methods */

	// computes arbitrary (start) state
	virtual bool computeArbitraryState(State& s) const = 0;

	// compute vector of neighbors of s
	virtual void computeNeighbors(const State& s,
			std::vector<State>& neighbors) const = 0;

	// methods for non-uniform chains
	virtual void computeWeights(std::vector<State>& states) {

	}

public:

	/**
	 * Construct State Graph
	 */
	void constructStatespace() {

		// Variables
		State s;
		int i, j;
		std::queue<int> q;
		std::vector<State> neighbors;
		typename std::vector<State>::iterator it;

		// Transition Matrix (sparse)
		boost::unordered_map<std::pair<int, int>, Rational> transition_matrix;
		typename boost::unordered_map<State, int>::iterator it2;
		typename boost::unordered_map<std::pair<int, int>, Rational>::iterator it3;

		// Start with arbitrary State
		if (!computeArbitraryState(s))
			return;

#ifdef MY_DEBUG
		std::cout << "Start state is " << s << std::endl;
#endif

		// BFS through state space
		q.push(0);
		states.push_back(s);
		indices[s] = 0;

		while (!q.empty()) {
			i = q.front();
			s = states[i];
			q.pop();

#ifdef MY_DEBUG
			std::cout << s << std::endl;
#endif

			neighbors.clear();
			computeNeighbors(s, neighbors);

			for (it = neighbors.begin(); it != neighbors.end(); ++it) {

				s = *it;
#ifdef MY_DEBUG
				std::cout << " " << s << std::flush;
#endif

				// Look if s already known
				it2 = indices.find(s);

				// s not seen so far
				if (it2 == indices.end()) {
					j = indices.size();
					indices[s] = j;
					states.push_back(s);
					q.push(j);
#ifdef MY_DEBUG
					std::cout << " new" << std::endl;
#endif
				}
				// s already known
				else {
					j = it2->second;
#ifdef MY_DEBUG
					std::cout << " skipped" << std::endl;
#endif
				}

				// add transition propability
				transition_matrix[std::make_pair(i, j)] += Rational(1,
						(int) neighbors.size());
			}
		}

		numStates = states.size();

		// compute weights for metropolis rule
		computeWeights(states);

		// apply Metropolis Rule

		// margin distribution (row-sums)
		Rational *margin = new Rational[numStates];
		for (i = 0; i < numStates; i++)
			margin[i] = 0;

		for (it3 = transition_matrix.begin(); it3 != transition_matrix.end();
				++it3) {

			i = it3->first.first;
			j = it3->first.second;
			Rational w_i = getWeight(i);
			Rational w_j = getWeight(j);
			Rational metr = 1;
			if (w_j < w_i)
				metr = w_j / w_i;

#ifdef MY_DEBUG
			std::cout << " (" << std::setw(2) << it3->first.first << ", "
			<< std::setw(2) << it3->first.second << "): " << std::setw(2)
			<< it3->second << " * " << metr << " = " << (it3->second * metr)
			<< std::endl;
#endif

			it3->second *= metr;
			margin[i] += it3->second;
		}

		// remaining probability goes to "stay-at-state" probabiltiy
		for (i = 0; i < numStates; i++) {
			std::pair<int, int> ii(i, i);
			transition_matrix[ii] += Rational(1) - margin[i];
		}

		// Compile Stationary Distrubtion
		stationary_distribution.resize(numStates);

		Rational Z = 0;
		for (i = 0; i < numStates; i++)
			Z += getWeight(i);

		for (i = 0; i < numStates; i++) {
			stationary_distribution[i] = getWeight(i) / Z;
		}

		// Transform unordered_map into edge array
		for (it3 = transition_matrix.begin(); it3 != transition_matrix.end();
				++it3) {
			int u = it3->first.first;
			int v = it3->first.second;
			Rational p_uv = it3->second;
			arcs.push_back(Transition(u, v, p_uv));
		}

		// sort edge array by start node
		sort(arcs.begin(), arcs.end(),
				TransitionComparator());

		// collect for each state the first of its arcs
		outgoing_arcs.resize(numStates + 1);
		outgoing_arcs[numStates] = arcs.size();
		for (int i = arcs.size() - 1; i >= 0; i--) {
			int u = arcs[i].u;
			outgoing_arcs[u] = i;
		}

#ifdef MY_DEBUG
		int ii = 0;
		std::cout << "state size: " << states.size() << std::endl;
		for (typename std::vector<State>::iterator it = states.begin();
				it != states.end(); ++it)
		std::cout << ii++ << ": " << *it << " " << std::endl;

		std::cout << "transition size: " << transition_matrix.size() << std::endl;
#endif

		delete[] margin;
	}

};

} /* end namespace marathon */

#endif /* MARKOVCHAIN_H_ */
