// system includes
#include <iostream>
#include <string>

// project includes
#include "marathon/chains/matching/matching_chain_JSV04.h"
#include "marathon/chains/matching/matching_chain_JS89.h"
#include "marathon/chains/sequences/switch_chain_bipartite.h"
#include "marathon/chains/sequences/switch_chain_bipartite_berger.h"

// marathon includes
#include "marathon/marathon.h"

int main(int argc, char** argv) {

	if (argc != 4) {
		std::cout
				<< "usage: totalMixingTime <js89|jsv04|swapBip|swapBipFast> <instance> <epsilon>"
				<< std::endl;
		return 1;
	}

	// command line arguments
	std::string inst(argv[2]);
	float eps = atof(argv[3]);

	// Declare StateGraph object
	marathon::StateGraph *sg = nullptr;

	// check which chain is selected
	if (strcmp(argv[1], "js89") == 0)
		sg = new marathon::chain::matching::Broder86(inst);
	else if (strcmp(argv[1], "jsv04") == 0)
		sg = new marathon::chain::matching::JerrumSinclairVigoda04(inst);
	else if (strcmp(argv[1], "swapBip") == 0)
		sg = new marathon::chain::sequence::SwitchBipartite(inst);
	else if (strcmp(argv[1], "swapBipFast") == 0)
		sg = new marathon::chain::sequence::SwitchBipartiteFast(inst);
	else {
		std::cerr << "unknown chain specifier: " << argv[1] << std::endl;
		return 1;
	}

	// construct state graph
	sg->constructStateGraph();

	// Init library
	marathon::init();

	/**
	 * compute total mixing time:
	 * 1. try gpu method
	 * 2. if gpu is not sucessfull (e.g. not enough memory) try cpu implementation
	 */
	int t = marathon::totalMixingTime<double>(sg, eps, marathon::device_t::GPU_ONLY);
	if(t == -1)
		t = marathon::totalMixingTime<double>(sg, eps, marathon::device_t::CPU_ONLY);

	// compute spectral gap
	double lambda = marathon::secondLargestEigenvalue<double>(sg);

	// print information
	std::cout << "number of states:          " << sg->getNumStates()
			<< std::endl;
	std::cout << "number of transition arcs: " << sg->getNumTransitions()
			<< std::endl;
	std::cout << "total mixing time:         " << t << std::endl;
	std::cout << "spectral gap:              " << (1-lambda) << std::endl;

	delete sg;

	// finalize library
	marathon::finalize();

	return 0;
}
