// system includes
#include <iostream>
#include <string>

#include "marathon/marathon.h"

int main(int argc, char** argv) {

	// check command line arguments
	if (argc != 4) {
		std::cout
				<< "usage: MixingBounds <js89|jsv04|swapBip|swapBipFast> <instance> <epsilon>"
				<< std::endl;
		return 1;
	}

	// parse command line arguments
	std::string inst(argv[2]);
	float eps = atof(argv[3]);

	// Init library
	marathon::init();

	// Declare Markov chain and path construction scheme objects
	marathon::MarkovChain* mc;
	marathon::PathConstructionScheme* pcs;

	// check which chain is selected
	if (strcmp(argv[1], "js89") == 0) {
		mc = new marathon::chain::matching::Broder86(inst);
		pcs = new marathon::chain::matching::JS89Path;
	} else if (strcmp(argv[1], "jsv04") == 0) {
		mc = new marathon::chain::matching::JerrumSinclairVigoda04(inst);
		pcs = new marathon::chain::matching::JS89Path;
	} else if (strcmp(argv[1], "swapBip") == 0) {
		mc = new marathon::chain::bipgraph::SwitchChain(inst);
		pcs = new marathon::chain::bipgraph::KannanPath;
	} else if (strcmp(argv[1], "swapBipFast") == 0) {
		mc = new marathon::chain::bipgraph::SwitchChainBerger(inst);
		pcs = new marathon::chain::bipgraph::KannanPath;
	} else {
		std::cerr << "unknown chain specifier: " << argv[1] << std::endl;
		return -1;
	}

	// construct state graph
	marathon::StateGraph* sg = new marathon::StateGraph(mc);

	// compute total mixing time
	int t = marathon::totalMixingTime<double>(sg, eps,
			marathon::device_t::CPU_ONLY);

	// compute spectral bound
	double lambda = marathon::eigenvalue<double>(sg,
			marathon::eigenvalue_t::_2ndLargestMagnitude);
	if (fabs(lambda) < 1e-15)
		lambda = 0.0;
	else
		lambda = fabs(lambda);

	double pimin = (sg->getMinWeight() / sg->getZ()).convert_to<double>();
	double lower_spectral = 0.5 * lambda / (1.0 - lambda) * -log(2.0 * eps);
	double upper_spectral = -log(pimin * eps) / (1.0 - lambda);

	// compute congestion bound
	const marathon::rational load = marathon::pathCongestion(sg, *pcs);
	double upper_congestion = load.convert_to<double>() * -log(pimin * eps);

	// print information
	std::cout << "number of states:          " << sg->getNumStates()
			<< std::endl;
	std::cout << "number of transition arcs: " << sg->getNumTransitions()
			<< std::endl;
	std::cout << "total mixing time:         " << t << std::endl;
	std::cout << "lower spectral bound:      " << lower_spectral << std::endl;
	std::cout << "upper spectral bound:      " << upper_spectral << std::endl;
	std::cout << "upper congestion bound:    " << upper_congestion << std::endl;

	delete mc;
	delete sg;
	delete pcs;

	// finalize library
	marathon::finalize();

	return 0;
}