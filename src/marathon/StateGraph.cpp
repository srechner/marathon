/*
 * StateGraph.cpp
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


#include <string>
#include "marathon/StateGraph.h"

namespace marathon {

StateGraph::StateGraph(MarkovChain* mc, const int limit) :
		mc(mc) {
	expand(limit);
}

StateGraph::~StateGraph() {
	// delete all transitions
	for (int i = 0; i < arcs.size(); i++)
		delete arcs[i];
}

size_t StateGraph::getNumStates() const {
	return states.size();
}

MarkovChain* StateGraph::getMarkovChain() const {
	return mc;
}

void StateGraph::expandState(const int i, const int limit, const int lastStop,
		const bool verbose) {

	if (verbose)
		std::cout << "expand state " << i << std::endl;

	const State* s = (const State*) getState(i);

	// compute states that are neighboured to s
	std::vector<std::pair<State*, rational>> neighbours;
	mc->computeNeighbours(s, neighbours);

	/**********************************************************************
	 *  eliminate duplicate states by sorting and compressing
	 **********************************************************************/
	std::sort(neighbours.begin(), neighbours.end(),
			[](const std::pair<State*, rational>& s1, const std::pair<State*, rational>& s2) -> bool {
				return s1.first->compare(s2.first) < 0;
			});

	// compress
	std::vector<std::pair<State*, rational>> tmp;
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
	neighbours = std::move(tmp);
	tmp.clear();

	bool complete = true;

	// sum of proposal probabilites of adjacent states
	rational sum(0);

	// for all adjacent states
	for (auto it = neighbours.begin(); it != neighbours.end(); ++it) {

		// neighbour state s2 with proposal prob. kappa(s,s2)
		State* s2 = it->first;
		const rational& kappa = it->second;
		sum += kappa;

		// Look if s2 already known
		int j = indexOf(s2);	// j is the index of state s2

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

void StateGraph::expand(const int limit, const bool verbose) {

	const int sizeLast = getNumStates();

	// nothing to do?
	if (limit <= sizeLast)
		return;

	// if state graph is empty
	if (sizeLast == 0) {

		// Start with arbitrary State
		State* s1 = mc->computeArbitraryState();
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
	std::vector<rational> generatedLoopProbability(getNumStates());

	// Transform Transition probabilities
	for (Transition* t : getArcs()) {

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
	for(int i=0; i<generatedLoopProbability.size(); i++) {
		addTransitionProbability(i, i, generatedLoopProbability[i]);
	}

	/***********************************************************
	 * Self-Check: Verifiy that chain is reversible
	 ***********************************************************/
	const rational Z = getZ();

	for (const Transition* t : getArcs()) {

		if (t->u >= t->v)
			continue;

		const rational stat_u = getWeight(t->u) / Z;
		const rational stat_v = getWeight(t->v) / Z;
		const rational puv = t->p;

		// find anti parallel transition arc
		for (const Transition* t2 : getInArcs(t->u)) {

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
			const State* s = getState(ii);
			std::cout << ii << ": " << s->toString() << " " << weights[ii]
					<< std::endl;
		}

		std::cout << "transition size: " << getNumTransitions() << std::endl;
		for (const Transition* t : getArcs()) {
			std::cout << t->u << " " << t->v << " " << t->p << std::endl;
		}
	}
}

int StateGraph::addArc(const int u, const int v, const rational& p) {
	return addArc(new Transition(u, v, p));
}

int StateGraph::addArc(Transition* t) {

	// add the arc to the arc vector
	arcs.push_back(t);

	// add a pointer to the transition in u's and v's outarc/inarc array
	outArcs[t->u].push_back(t);
	inArcs[t->v].push_back(t);

	return arcs.size() - 1;
}

Transition* StateGraph::getArc(int u, int v) const {

	// search for transition (u,v)
	for (Transition* t : outArcs[u]) {
		assert(t->u == u);
		if (t->v == v) {
			return t;
		}
	}
	return nullptr;
}

size_t StateGraph::getNumTransitions() const {
	return arcs.size();
}

rational StateGraph::getTransitionProbability(int u, int v) const {
	const Transition* t = getArc(u, v);
	if (t == nullptr)
		return 0;
	else
		return t->p;
}

void StateGraph::setTransitionProbability(int u, int v, rational p) {

// search for transition (u,v)
	Transition* t = getArc(u, v);

	if (t != nullptr) {
		t->p = p;
	} else {
		// no transition found? add a new one
		addArc(new Transition(u, v, p));
	}

}

void StateGraph::addTransitionProbability(int u, int v, rational p) {

// search for transition (u,v)
	Transition* t = getArc(u, v);

	if (t != nullptr) {
		t->p += p;
	} else {
		// no transition found? add a new one
		addArc(new Transition(u, v, p));
	}
}

void StateGraph::setWeight(const int i, const rational p) {
	weights[i] = p;
}

rational StateGraph::getWeight(const int i) const {
	return weights[i];
}

const std::vector<rational>& StateGraph::getWeights() const {
	return weights;
}

rational StateGraph::getZ() const {
	rational Z = 0;
	for (auto w : weights)
		Z += w;
	return Z;
}

rational StateGraph::getMinWeight() const {
	return *std::min(weights.begin(), weights.end());
}

const std::vector<Transition*>& StateGraph::getOutArcs(int v) const {
	return outArcs[v];
}

const std::vector<Transition*>& StateGraph::getInArcs(int v) const {
	return inArcs[v];
}

const std::vector<Transition*>& StateGraph::getArcs() const {
	return arcs;
}

Transition* StateGraph::getArc(const int i) const {
	return (Transition*) &arcs[i];
}

int StateGraph::getNumOutArcs(int v) const {
	return this->getOutArcs(v).size();
}

void StateGraph::clear() {
	arcs.clear();
	inArcs.clear();
	outArcs.clear();
	weights.clear();
}

int StateGraph::addState(State* s) {
// add state to the vector of states
	states.push_back(s);
	indices[s] = states.size() - 1;
	outArcs.push_back(std::vector<Transition*>());
	inArcs.push_back(std::vector<Transition*>());
	weights.push_back(1);
	return states.size() - 1;
}

const State* StateGraph::getState(int i) const {
	return states[i];
}

const std::vector<const State*>& StateGraph::getStates() const {
	return states;
}

int StateGraph::indexOf(const State * s) const {
	auto it = indices.find((State*) s);
	if (it != indices.end())
		return it->second;
	else
		return -1;
}

std::ostream& operator<<(std::ostream& out, const StateGraph& sg) {

	out << "n " << sg.getNumStates() << " m " << sg.getNumTransitions() << "\n";
	for (int i = 0; i < sg.getNumStates(); i++) {
		for (Transition* t : sg.getOutArcs(i))
			out << t->u << " " << t->v << " " << t->p << "\n";
	}

	return out;
}

}
