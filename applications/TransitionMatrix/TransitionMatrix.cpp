/*
 * TransitionMatrix.cpp
 *
 * Created on: Nov 24, 2014
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

#include "marathon/marathon.h"
#include "marathon/TransitionMatrixCBLAS.h"
#include "marathon/TransitionMatrixCuBLAS.h"

int main(int argc, char** argv) {

	if (argc != 3) {
		std::cout
				<< "usage: transition_matrix <js89|jsv04|swapBip|swapBipFast> <instance>"
				<< std::endl;
		return 1;
	}

	// parse command line arguments
	std::string inst(argv[2]);

	// Init library
	marathon::init();

	// Abstract Markov Chain Object
	marathon::MarkovChain *mc;

	// check which chain is selected and instanciate with proper implementation
	if (strcmp(argv[1], "js89") == 0)
		mc = new marathon::chain::matching::Broder86(inst);
	else if (strcmp(argv[1], "jsv04") == 0)
		mc = new marathon::chain::matching::JerrumSinclairVigoda04(inst);
	else if (strcmp(argv[1], "swapBip") == 0)
		mc = new marathon::chain::bipgraph::SwitchChain(inst);
	else if (strcmp(argv[1], "swapBipFast") == 0)
		mc = new marathon::chain::bipgraph::SwitchChainBerger(inst);
	else {
		std::cerr << "unknown chain specifier: " << argv[1] << std::endl;
		return 1;
	}

	// declare State Graph object
	marathon::StateGraph* sg = new marathon::StateGraph(mc);

	// print transition matrix
	marathon::TransitionMatrix<double>* P = new marathon::TransitionMatrixCBLAS<double>(sg);
	std::cout << P << std::endl;

	delete sg;
	delete mc;
	delete P;

	// finalize library
	marathon::cleanup();

	return 0;
}
