/*
 * MarkovChain.cpp
 *
 * Created on: Feb 15, 2016
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

#include "../../include/marathon/MarkovChain.h"

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
