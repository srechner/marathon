/*
 * markov_chain.cpp
 *
 *  Created on: Feb 15, 2016
 *      Author: steffen
 */

#include "../../include/marathon/MarkovChain.h"
#include <cstdlib>

marathon::MarkovChain::~MarkovChain() {

}

std::string marathon::MarkovChain::getName() const {
	return "[unknown]";
}

void marathon::MarkovChain::computeWeights(
		const std::vector<const State*>& states,
		std::vector<rational>& weights) {
	weights.clear();
	weights.reserve(states.size());
	for (int i = 0; i < states.size(); i++)
		weights.push_back(1);
}

void marathon::MarkovChain::randomize(State* s, const uint32_t t) const {
	// dummy implementation
	std::cout << "marathon::Exception: randomize is not implemented!"
			<< std::endl;
}

marathon::State * marathon::MarkovChain::computeArbitraryState() const {
	// dummy implementation
	std::cout << "marathon::Exception: computeArbitraryState is not implemented!"
			<< std::endl;
	return nullptr;
}

void marathon::MarkovChain::computeNeighbours(const State* s,
		std::vector<std::pair<State*, rational>>& neighbors) const {
	// dummy implementation
	std::cout << "marathon::Exception: computeNeighbouts is not implemented!"
			<< std::endl;
}

marathon::rational marathon::MarkovChain::loopProbability(
		const State* s) const {

	rational res(0);
	std::vector<std::pair<State*, rational>> N;
	this->computeNeighbours(s, N);
	for (auto x : N) {
		if (x.first->compare(s) == 0) {
			res += x.second;
		}
		delete x.first;
	}
	return res;
}
