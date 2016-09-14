// system includes
#include <iostream>
#include <string>

#include "marathon/marathon.h"

int main(int argc, char **argv) {

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
    marathon::MarkovChain *mc;
    marathon::PathConstructionScheme *pcs;

    // check which chain is selected
    if (strcmp(argv[1], "js89") == 0) {
        mc = new marathon::chain::matching::Broder86(inst);
        pcs = new marathon::chain::matching::JS89Path();
    } else if (strcmp(argv[1], "jsv04") == 0) {
        mc = new marathon::chain::matching::JerrumSinclairVigoda04(inst);
        pcs = new marathon::chain::matching::JS89Path();
    } else if (strcmp(argv[1], "swapBip") == 0) {
        mc = new marathon::chain::bipgraph::SwitchChain(inst);
        pcs = new marathon::chain::bipgraph::KannanPath();
    } else if (strcmp(argv[1], "swapBipFast") == 0) {
        mc = new marathon::chain::bipgraph::SwitchChainBerger(inst);
        pcs = new marathon::chain::bipgraph::KannanPath();
    } else {
        std::cerr << "unknown chain specifier: " << argv[1] << std::endl;
        return -1;
    }

    // construct state graph
    marathon::StateGraph *sg = new marathon::StateGraph(mc);

    // compute total mixing time
    int t = marathon::totalMixingTime<double>(sg, eps,
                                              marathon::device_t::CPU_ONLY);

    // compute bounds
    double lower_spectral = marathon::lowerSpectralBound(sg, eps);
    double upper_spectral = marathon::upperSpectralBound(sg, eps);
    double upper_congestion = marathon::upperPathCongestionBound(sg, *pcs, eps);

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
