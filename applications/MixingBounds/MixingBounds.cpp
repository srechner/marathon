/*
 * MixingBounds.cpp
 *
 * Created on: Nov 24, 2015
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2016, Steffen Rechner
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// system includes
#include <iostream>

#include "marathon/marathon.h"
#include "marathon/chain/bipgraph/KannanCanPath.h"
#include "marathon/chain/matching/JS89CanPath.h"
#include "marathon/MixingTime.h"

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
    int t = marathon::MixingTime::totalMixingTime<double>(sg, eps,
                                              marathon::MixingTime::device_t::CPU_ONLY);

    // compute bounds
    double lower_spectral = marathon::MixingTime::lowerSpectralBound(sg, eps);
    double upper_spectral = marathon::MixingTime::upperSpectralBound(sg, eps);
    double upper_congestion = marathon::MixingTime::upperPathCongestionBound(sg, *pcs, eps);

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
    marathon::cleanup();

    return 0;
}
