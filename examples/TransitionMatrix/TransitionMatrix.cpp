// system includes
#include <iostream>
#include <string>

#include "marathon/marathon.h"

int main(int argc, char** argv) {

	if (argc < 3 || argc > 4) {
		std::cout
				<< "usage: transition_matrix <js89|jsv04|swapBip|swapBipFast|curveball> <instance> [max-states]"
				<< std::endl;
		return 1;
	}

	// parse command line arguments
	std::string inst(argv[2]);

	int limit = INT_MAX;
	if (argc == 4)
		limit = atoi(argv[3]);

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
	else if (strcmp(argv[1], "curveball") == 0)
		mc = new marathon::chain::bipgraph::Curveball(inst);
	else {
		std::cerr << "unknown chain specifier: " << argv[1] << std::endl;
		return 1;
	}

	// declare State Graph object
	marathon::StateGraph* sg = new marathon::StateGraph(mc, limit);

	// print transition matrix
	marathon::TransitionMatrix<double>* P = new marathon::TransitionMatrixCBLAS<double>(sg);
	std::cout << P << std::endl;

	delete sg;
	delete mc;
	delete P;

	// finalize library
	marathon::finalize();

	return 0;
}
