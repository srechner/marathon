// system includes
#include <iostream>

// marathon includes
#include "marathon/spectral_bound.h"
#include "marathon/binary_matrix/fixed_margin/switch_chain.h"
#include "marathon/binary_matrix/fixed_margin/edge_switch_chain.h"
#include "marathon/binary_matrix/fixed_margin/curveball.h"

// auxiliary functions
#include "helper.h"

int main(int argc, char **argv) {

    // parse command line arguments
    if (!parse_arguments(argc, argv)) {
        print_help_message();
        return -1;
    }

    // create Markov chain
    std::unique_ptr<marathon::MarkovChain> mc;
    switch (chain) {
        case classical_switch:
            mc = std::make_unique<marathon::binary_matrix::fixed_margin::SwitchChain>(inst);
            break;
        case edge_switch:
            mc = std::make_unique<marathon::binary_matrix::fixed_margin::EdgeSwitchChain>(inst);
            break;
        case curveball:
            mc = std::make_unique<marathon::binary_matrix::fixed_margin::Curveball>(inst);
            break;
    }

    // construct state graph
    marathon::StateGraph sg(*mc);

    // calculate spectral bound of the total mixing time
    marathon::SpectralBoundCalculator<double> sbc(sg);

    // lower spectral bound
    double lspec = sbc.lowerSpectralBound(eps);

    // upper spectral bound
    double uspec = sbc.upperSpectralBound(eps);

    std::cout << "lower spectral bound: " << lspec << std::endl;
    std::cout << "upper spectral bound: " << uspec << std::endl;

    return 0;
}
