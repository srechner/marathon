// system includes
#include <iostream>
#include <string>

// marathon includes
#include "marathon/marathon.h"

int main(int argc, char** argv) {

	if (argc != 3 && argc != 4) {
		std::cout << "usage: states <matching|bipgraph> <instance> [max-states]"
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
	marathon::StateGraph* sg = new marathon::StateGraph(mc, limit);

	// output all states
	for (const marathon::State* s : sg->getStates())
		std::cout << *s << std::endl;

	delete mc;
	delete sg;

	return 0;
}
