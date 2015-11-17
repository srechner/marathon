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

namespace marathon {

// provide hash function for pair<int,int>
struct PairHasher {
public:
	size_t operator()(const std::pair<int, int>& x) const {
		size_t h = x.first * 31 + x.second;
		return h;
	}
};

// provide == function for pair<int,int>
struct PairEqual {
public:
	bool operator()(const std::pair<int, int>& t1,
			const std::pair<int, int>& t2) const {
		return t1.first == t2.first && t1.second == t2.second;
	}
};

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

	/**
	 * computes arbitrary (start) state
	 */
	virtual bool computeArbitraryState(State& s) const = 0;

	/**
	 *  compute vector of adjacent states of s with corresponding proposal prob.
	 */
	virtual void computeNeighbours(const State& s,
			boost::unordered_map<State, Rational>& neighbors) const = 0;

	/**
	 *  How to compute weights for each state.
	 */
	virtual void computeWeights(std::vector<State>& states) {

	}

public:

	/**
	 * Construct State Graph
	 */
	void constructStateGraph(bool verbose) {

		// Declare Variables
		State s1, s2;
		int i, j;
		std::queue<int> q;
		boost::unordered_map<State, Rational> neighbours;
		typename boost::unordered_map<State, int>::iterator it2;

		// Transition Matrix as temporary storage(sparse)
		//boost::unordered_map<std::pair<int, int>, Rational, PairHasher, PairEqual> transition_matrix;
		std::map<std::pair<int, int>, Rational> transition_matrix;

		// Start with arbitrary State
		if (!computeArbitraryState(s1)) {
			if (verbose)
				std::cout << "Empty Statespace!" << std::endl;
			return;
		}

		if (verbose)
			std::cout << "Start state is " << s1 << std::endl;

		q.push(0);
		states.push_back(s1);
		indices[s1] = 0;

		// Run BFS through state space
		while (!q.empty()) {
			i = q.front();
			s1 = states[i];
			q.pop();

			if (verbose)
				std::cout << i << ": " << s1 << std::endl;

			neighbours.clear();
			computeNeighbours(s1, neighbours);

			// check for valid transition rules: proposal probability must sum to 1
			Rational sum = 0;
			for (auto v : neighbours) {

				// neighbour state s2 with proposal prob. kappa(s1,s2)
				s2 = v.first;
				Rational kappa = v.second;
				sum += kappa;

				// Look if s2 already known
				it2 = indices.find(s2);

				// s2 not seen so far
				if (it2 == indices.end()) {
					j = indices.size();
					indices[s2] = j;
					states.push_back(s2);
					q.push(j);

					if (verbose)
						std::cout << " " << j << ": " << s2 << " " << kappa
								<< " new" << std::endl;
				}
				// s2 already known
				else {
					j = it2->second;
					if (verbose)
						std::cout << " " << j << ": " << s2 << " " << kappa
								<< " skipped" << std::endl;
				}

				// add proposal propability kappa(i,j)
				transition_matrix[std::make_pair(i, j)] = kappa;

			}
			assert(sum == Rational(1));
		}

		numStates = states.size();

		/***********************************************************
		 * Compute weights for metropolis rule
		 **********************************************************/
		computeWeights(states);

		/*********************************************************************
		 * Print Information about states
		 ********************************************************************/

		if (verbose) {
			int ii = 0;
			std::cout << "state size: " << states.size() << std::endl;
			for (typename std::vector<State>::iterator it = states.begin();
					it != states.end(); ++it) {
				std::cout << ii << ": " << *it << " " << getWeight(ii)
						<< std::endl;
				ii++;
			}

			std::cout << "transition size: " << transition_matrix.size()
					<< std::endl;
		}

		/*******************************************************
		 * Compile Stationary Distrubtion
		 ******************************************************/
		stationary_distribution.resize(numStates);

		Rational Z = 0;
		for (i = 0; i < numStates; i++)
			Z += getWeight(i);

		for (i = 0; i < numStates; i++) {
			stationary_distribution[i] = getWeight(i) / Z;
		}

		/***********************************************************
		 * Apply Metropolis Rule:
		 * P(i,j) = kappa(i,j) * min( w(j) / w(i) , 1)
		 ************************************************************/

		// Transform Transition probabilities
		for (auto it = transition_matrix.begin(); it != transition_matrix.end();
				++it) {

			i = it->first.first;
			j = it->first.second;
			Rational kappa = it->second;
			Rational w_i = getWeight(i);
			Rational w_j = getWeight(j);

			// metr = min(w_j / w_i, 1)
			Rational metr(1);
			if (w_j < w_i)
				metr = w_j / w_i;

			if (verbose)
				std::cout << " (" << std::setw(2) << i << ", " << std::setw(2)
						<< j << "): " << std::setw(2) << kappa << " * " << metr
						<< " = " << (kappa * metr) << std::endl;

			// add remaining probability as loop probability
			if (kappa != Rational(1)) {
				std::pair<int, int> ii = std::make_pair(i, i);
				transition_matrix[ii] += (Rational(1) - metr) * kappa;
			}

			// P(i,j) = kappa(i,j) * metr
			it->second *= metr;
		}

		/***********************************************************
		 * Check that chain is reversible
		 ***********************************************************/
		for (auto it = transition_matrix.begin(); it != transition_matrix.end();
				++it) {
			int i = it->first.first;
			int j = it->first.second;
			Rational pij = it->second;
			Rational pji = transition_matrix[std::make_pair(j, i)];
			Rational stat_i = getStationary(i);
			Rational stat_j = getStationary(j);

			if (stat_i * pij != stat_j * pji) {
				std::cerr << "Error! Chain is not reversible!" << std::endl;
				std::cerr << "P(" << i << "," << j << ")=" << pij << std::endl;
				std::cerr << "P(" << j << "," << i << ")=" << pji << std::endl;
				std::cerr << "stat(" << i << ")=" << stat_i << std::endl;
				std::cerr << "stat(" << j << ")=" << stat_j << std::endl;
				std::cerr << "stat(" << i << ")*P(" << i << "," << j << ")="
						<< stat_i * pij << std::endl;
				std::cerr << "stat(" << j << ")*P(" << j << "," << i << ")="
						<< stat_j * pji << std::endl;
			}

			assert(stat_i * pij == stat_j * pji);
		}

		/********************************************************
		 * Transform unordered_map into edge array
		 *******************************************************/
		for (auto it = transition_matrix.begin(); it != transition_matrix.end();
				++it) {
			int u = it->first.first;
			int v = it->first.second;
			Rational p_uv = it->second;
			arcs.push_back(Transition(u, v, p_uv));
		}

		// sort edge array by start node
		sort(arcs.begin(), arcs.end(), TransitionComparator());

		// collect for each state the first of its arcs
		outgoing_arcs.resize(numStates + 1);
		for (i = 0; i <= numStates; i++)
			outgoing_arcs[i] = arcs.size();
		for (int i = arcs.size() - 1; i >= 0; i--) {
			int u = arcs[i].u;
			outgoing_arcs[u] = i;
		}

	}

};

} /* end namespace marathon */

#endif /* MARKOVCHAIN_H_ */
