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

	if (argc != 3) {
		std::cout << "usage: states <matching|bipgraph> <instance>"
				<< std::endl;
		return 1;
	}

	// command line arguments
	std::string inst(argv[2]);

	marathon::StateGraph *sg;

	// check which chain is selected
	if (strcmp(argv[1], "matching") == 0) {
		sg = new marathon::chain::matching::Broder86(inst);
	} else if (strcmp(argv[1], "bipgraph") == 0) {
		sg = new marathon::chain::sequence::SwitchBipartite(inst);
	} else {
		std::cerr << "unknown mode: " << argv[1] << std::endl;
		return 1;
	}

	sg->constructStateGraph();
	sg->printStates();

	delete sg;

	return 0;
}
