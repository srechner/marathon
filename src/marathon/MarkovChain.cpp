/*
 * markov_chain.cpp
 *
 *  Created on: Feb 15, 2016
 *      Author: steffen
 */

#include "../../include/marathon/MarkovChain.h"
#include <cstdlib>

marathon::MarkovChain::MarkovChain(const std::string& s, int seed) :
		instance(s) {
	srand(seed);
}

marathon::MarkovChain::~MarkovChain() {

}

const std::string& marathon::MarkovChain::getInstance() const {
	return this->instance;
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


void marathon::MarkovChain::randomize(State* s) const {
	// dummy implementation
}
