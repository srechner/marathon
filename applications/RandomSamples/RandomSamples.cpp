// system includes
#include <iostream>
#include <string>

// marathon includes
#include "marathon/marathon.h"

int main(int argc, char** argv) {

	if (argc != 5) {
		std::cout << "Create N random samples using the specified Markov chain.\n"
				<< std::endl;
		std::cout << "usage: RandomSamples <chain> <instance> <N> <t>"
				<< std::endl;
		std::cout
				<< "   chain: Markov chain specifier. Allowed values are: [js89,jsv04,switch,curveball,linda]."
				<< std::endl;
		std::cout
				<< "   instance: encoded instance string (parsable by the specified chain)."
				<< std::endl;
		std::cout << "   N: number of samples" << std::endl;
		std::cout << "   t: length of random walk" << std::endl;
		return 1;
	}

	marathon::init();

	// command line argument
	int N, t;
	std::string inst(argv[2]);
	N = atoi(argv[3]);
	t = atoi(argv[4]);

	// Abstract Markov Chain Object
	marathon::MarkovChain *mc;

	// check which chain is selected and instanciate with proper implementation
	if (strcmp(argv[1], "js89") == 0)
		mc = new marathon::chain::matching::Broder86(inst);
	else if (strcmp(argv[1], "jsv04") == 0)
		mc = new marathon::chain::matching::JerrumSinclairVigoda04(inst);
	else if (strcmp(argv[1], "switch") == 0)
		mc = new marathon::chain::bipgraph::SwitchChain(inst);
	else if (strcmp(argv[1], "curveball") == 0)
		mc = new marathon::chain::bipgraph::Curveball(inst);
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

	marathon::finalize();

	return 0;
}
