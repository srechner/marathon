/*
 * ListStates.cpp
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
#include <climits>

// marathon includes
#include "marathon/marathon.h"

int main(int argc, char** argv) {

	if (argc != 3) {
		std::cout << "usage: states <matching|bipgraph> <instance>"
				<< std::endl;
		return 1;
	}

	// command line argument
	int limit = INT_MAX;
	std::string inst(argv[2]);
	if (argc == 4)
		limit = atoi(argv[3]);

	// declare Markov Chain Object
	marathon::MarkovChain* mc;

	// check which chain is selected
	if (strcmp(argv[1], "matching") == 0) {
		mc = new marathon::chain::matching::Broder86(inst);
	} else if (strcmp(argv[1], "bipgraph") == 0) {
		mc = new marathon::chain::bipgraph::SwitchChain(inst);
	} else {
		std::cerr << "unknown mode: " << argv[1] << std::endl;
		return 1;
	}

	// declare state graph object
	marathon::StateGraph* sg = new marathon::StateGraph(mc);

	// output all states
	for (const marathon::State* s : sg->getStates())
		std::cout << *s << std::endl;

	delete mc;
	delete sg;

	return 0;
}
