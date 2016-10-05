/*
 * RandomSamples.cpp
 *
 * Created on: Nov 24, 2015
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

// system includes
#include <iostream>

// marathon includes
#include "marathon/marathon.h"

int main(int argc, char** argv) {

	if (argc != 5) {
		std::cout << "Create N random samples using the specified Markov chain.\n"
				<< std::endl;
		std::cout << "usage: RandomSamples <chain> <instance> <N> <t>"
				<< std::endl;
		std::cout
				<< "   chain: Markov chain specifier. Allowed values are: [js89,jsv04,switch]."
				<< std::endl;
		std::cout
				<< "   instance: encoded instance string (parsable by the specified chain)."
				<< std::endl;
		std::cout << "   N: number of samples" << std::endl;
		std::cout << "   t: length of random walk" << std::endl;
		return 1;
	}

	// command line argument
	int N, t;
	std::string inst(argv[2]);
	N = atoi(argv[3]);
	t = atoi(argv[4]);

	// Abstract Markov Chain Object
	marathon::MarkovChain *mc;

	// check which chain is selected and instanciate with proper implementation
	if (strcmp(argv[1], "js89") == 0)
		mc = new marathon::matching::Broder86(inst);
	else if (strcmp(argv[1], "jsv04") == 0)
		mc = new marathon::matching::JerrumSinclairVigoda04(inst);
	else if (strcmp(argv[1], "switch") == 0)
		mc = new marathon::bipgraph::SwitchChain(inst);
	else {
		std::cerr << "unknown chain specifier: " << argv[1] << std::endl;
		return 1;
	}

	// create random samples
	marathon::State* s = mc->computeArbitraryState();

	for (int i = 0; i < N; i++) {
		mc->randomize(s, t);
		std::cout << *s << std::endl;
	}

	delete mc;

	return 0;
}
