// system includes
#include <iostream>
#include <string>

// project includes
#include "marathon/cpu/transition_matrix.h"
#include "marathon/chains/matching/matching_chain_JSV04.h"
#include "marathon/chains/matching/matching_chain_JS89.h"
#include "marathon/chains/sequences/switch_chain_bipartite.h"
#include "marathon/chains/sequences/switch_chain_bipartite_berger.h"

// marathon includes
#include "marathon/marathon.h"

int main(int argc, char** argv) {

	if (argc != 3) {
		std::cout
				<< "usage: transition_matrix <js89|jsv04|swapBip|swapBipFast> <instance>"
				<< std::endl;
		return 1;
	}

	// command line arguments
	std::string inst(argv[2]);

	// Init library
	marathon::init();

	// Abstract Markov Chain Object
	marathon::SamplingChain *mc;

	// check which chain is selected and instanciate with proper implementation
	if (strcmp(argv[1], "js89") == 0)
		mc = new marathon::chain::matching::Broder86();
	else if (strcmp(argv[1], "jsv04") == 0)
		mc = new marathon::chain::matching::JerrumSinclairVigoda04();
	else if (strcmp(argv[1], "swapBip") == 0)
		mc = new marathon::chain::sequence::SwitchBipartite();
	else if (strcmp(argv[1], "swapBipFast") == 0)
		mc = new marathon::chain::sequence::SwitchBipartiteFast();
	else {
		std::cerr << "unknown chain specifier: " << argv[1] << std::endl;
		return 1;
	}

	// State Graph Instance
	marathon::StateGraph* sg = mc->constructStateGraph(inst);

	// print transition matrix
	marathon::cpu::DenseTransitionMatrix<double> P;
	P.initFromStateGraph(sg);
	std::cout << P << std::endl;

	delete sg;
	delete mc;

	// finalize library
	marathon::finalize();

	return 0;
}
