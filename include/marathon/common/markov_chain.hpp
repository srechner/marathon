/*
 * @file: MarkovChain.h
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
#include <unordered_map>
#include <map>

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

/**
 * Parent-Class of Markov-Chain. This class just represents a abstraction
 * of a Markov Chain that does not need a template parameter.
 */
class SamplingChain {

public:

	virtual ~SamplingChain() {

	}

	/**
	 * The only method of this class is to construct a state graph out
	 * of the input instance.
	 */
	virtual StateGraph* constructStateGraph(bool verbose = false) = 0;

	virtual void constructPath(const StateGraph* sg, int i, int j, std::list<int>& path) const {

	}

};

template<typename State>
class MarkovChain: public SamplingChain {

protected:

	MarkovChain() {
	}

	virtual ~MarkovChain() {

	}

	/* Virtual Methods */

	/**
	 * computes arbitrary (start) state
	 */
	virtual bool computeArbitraryState(State& s) = 0;

	/**
	 *  Compute the set of adjacent states of s with corresponding proposal propability.
	 */
	virtual void computeNeighbours(const State& s,
			std::unordered_map<State, rational>& neighbors) const = 0;

	/**
	 *  Computes weights for each state.
	 */
	virtual void computeWeights(std::vector<State>& states, std::vector<rational>& weights) {
		weights.clear();
		for(int i=0; i<states.size(); i++)
			weights.push_back(1);
	}

public:

	/**
	 * Construct a State Graph for the given input instance.
	 */
	StateGraph* constructStateGraph(bool verbose = false) {

		// Declare Variables
		State s1, s2;
		int i, j;
		std::unordered_map<State, rational> neighbours;
		typename std::unordered_map<State, int>::iterator it2;

		// Temporary variables for construction of a state graph
		// TODO: Replace map by unordered map
		//boost::unordered_map<std::pair<int, int>, rational, PairHasher, PairEqual> transition_matrix;
		std::map<std::pair<int, int>, rational> transition_matrix;
		std::queue<int> q;

		std::vector<State> states;					// vector of states
		std::unordered_map<State, int> indices;	// State -> Index

		// Start with arbitrary State
		if (!computeArbitraryState(s1)) {
			// exception if statespace is empty
			throw std::runtime_error("empty statespace!");
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
			rational sum = 0;
			for (auto v : neighbours) {

				// neighbour state s2 with proposal prob. kappa(s1,s2)
				s2 = v.first;
				rational kappa = v.second;
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
			assert(sum == rational(1));
		}

		// create new state graph
		int n = states.size();
		_StateGraph<State>* sg = new _StateGraph<State>(n);

		// add the states to the state graph
		for (int i = 0; i < n; i++)
			sg->addState(states[i]);

		// tell the state graph how the input instance looks like
		//sg->setInstance();

		/***********************************************************
		 * Compute weights for metropolis rule
		 **********************************************************/
		std::vector<rational> weights;
		computeWeights(states, weights);

		/*********************************************************************
		 * Print Information about states
		 ********************************************************************/

		if (verbose) {
			std::cout << "state size: " << states.size() << std::endl;
			for (int ii=0; ii<states.size(); ii++) {
				std::cout << ii << ": " << states[ii] << " " << weights[ii]
						<< std::endl;
			}

			std::cout << "transition size: " << transition_matrix.size()
					<< std::endl;
		}

		/*******************************************************
		 * Compile Stationary Distrubtion
		 ******************************************************/
		rational Z = 0;
		for (i = 0; i < n; i++)
			Z += weights[i];

		for (i = 0; i < n; i++) {
			sg->setStationary(i, weights[i] / Z);
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
			rational kappa = it->second;
			rational w_i = weights[i];
			rational w_j = weights[j];

			// metr = min(w_j / w_i, 1)
			rational metr(1);
			if (w_j < w_i)
				metr = w_j / w_i;

			if (verbose)
				std::cout << " (" << std::setw(2) << i << ", " << std::setw(2)
						<< j << "): " << std::setw(2) << kappa << " * " << metr
						<< " = " << (kappa * metr) << std::endl;

			// add remaining probability as loop probability
			if (kappa != rational(1)) {
				std::pair<int, int> ii = std::make_pair(i, i);
				transition_matrix[ii] += (rational(1) - metr) * kappa;
			}

			// P(i,j) = kappa(i,j) * metr
			it->second *= metr;
		}

		/***********************************************************
		 * Self-Check: Verifiy that chain is reversible
		 ***********************************************************/
		for (auto it = transition_matrix.begin(); it != transition_matrix.end();
				++it) {
			int i = it->first.first;
			int j = it->first.second;
			rational pij = it->second;
			rational pji = transition_matrix[std::make_pair(j, i)];
			rational stat_i = sg->getStationary(i);
			rational stat_j = sg->getStationary(j);

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
		 * Update state graph.
		 *******************************************************/
		for (auto it = transition_matrix.begin(); it != transition_matrix.end();
				++it) {
			int u = it->first.first;
			int v = it->first.second;
			rational p_uv = it->second;
			sg->addArc(u, v, p_uv);
		}

		return (StateGraph*) sg;
	}
};

} /* end namespace marathon */

#endif /* MARKOVCHAIN_H_ */
