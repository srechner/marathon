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

	// command line argument
	std::string inst(argv[2]);

	// check which chain is selected
	if (strcmp(argv[1], "matching") == 0) {

		// declare markov chain
		typedef marathon::chain::matching::Broder86 Chain;
		typedef marathon::chain::matching::BipartiteMatching State;
		Chain mc;

		// construct state graph
		marathon::StateGraph* sg = mc.constructStateGraph(inst);

		// reinterpret state graph
		marathon::_StateGraph<State>* _sg = (marathon::_StateGraph<State>*) sg;

		// output all states
		for (auto s : _sg->getStates())
			std::cout << s << std::endl;

		delete sg;

	} else if (strcmp(argv[1], "bipgraph") == 0) {

		// declare markov chain
		typedef marathon::chain::sequence::SwitchBipartite Chain;
		typedef marathon::chain::sequence::DenseBipartiteGraph State;
		Chain mc;

		// construct state graph
		marathon::StateGraph* sg = mc.constructStateGraph(inst);

		// reinterpret state graph
		marathon::_StateGraph<State>* _sg = (marathon::_StateGraph<State>*) sg;

		// output all states
		for (auto s : _sg->getStates())
			std::cout << s << std::endl;

		delete sg;

	} else {
		std::cerr << "unknown mode: " << argv[1] << std::endl;
		return 1;
	}

	return 0;
}
