// system includes
#include <iostream>
#include <string>

// project includes
#include "marathon/chains/matching/chain_JSV04.h"
#include "marathon/chains/matching/chain_JS89.h"
#include "marathon/chains/sequences/switch_chain_bipartite.h"
#include "marathon/chains/sequences/switch_chain_bipartite_fast.h"

// marathon includes
#include "marathon.h"

//using namespace std;
//using namespace marathon;

int main(int argc, char** argv) {

	if (argc != 4) {
		std::cout << "usage: totalMixingTime <js89|jsv04|swapBip|swapBipFast> <instance> <epsilon>"
				<< std::endl;
		return 1;
	}

	// command line arguments
	std::string inst(argv[2]);
	float eps = atof(argv[3]);

	marathon::StateGraph *sg = nullptr;
	
	// check which chain is selected
	if (strcmp(argv[1], "js89") == 0)
		sg = new marathon::chain::matching::MatchingChain89(inst);
	else if (strcmp(argv[1], "jsv04") == 0)
		sg = new marathon::chain::matching::MatchingChain04(inst);
	else if (strcmp(argv[1], "swapBip") == 0)
		sg = new marathon::chain::sequence::SwapChainBipartite(inst);
	else if (strcmp(argv[1], "swapBipFast") == 0)
		sg = new marathon::chain::sequence::SwapChainBipartiteFast(inst);
	else {
		std::cerr << "unknown chain specifier: " << argv[1] << std::endl;
		return 1;
	}

	// construct state graph
	sg->constructStatespace();

	// compute total mixing time
	int t = marathon::cpu::totalMixingTime<double>(sg, eps);

	std::cout << t << std::endl;
	
	delete sg;

	return 0;
}
