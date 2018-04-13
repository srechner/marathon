//
// Created by rechner on 3/1/18.
//

#ifndef MARATHON_HELPER_H
#define MARATHON_HELPER_H

#include <iostream>

// Markov chains
enum chain_t {
    classical_switch,
    edge_switch,
    curveball,
} chain;

std::string inst;       // string encoded problem instance

/**
 * Parse command line arguments and extract parameters.
 * @param argc Number of arguments.
 * @param argv Argument list.
 * @return True, if arguments are specified correctly. False, otherwise.
 */
bool parse_arguments(int argc, char **argv) {
    if (argc != 3)
        return false;

    // parse problem type
    if (strcmp(argv[1], "classical-switch") == 0)
        chain = classical_switch;
    else if (strcmp(argv[1], "edge-switch") == 0)
        chain = edge_switch;
    else if (strcmp(argv[1], "curveball") == 0)
        chain = curveball;
    else {
        std::cerr << "Unknown CHAIN specifier: " << argv[1] << std::endl;
        return false;
    }

    inst = std::string(argv[2]);

    return true;
}

void print_help_message() {
    std::cout << "Usage: transitionMatrix CHAIN INSTANCE" << std::endl;
    std::cout << "Construct the transition matrix of a specified Markov chain." << std::endl;
    std::cout << std::endl;
    std::cout << "The parameter CHAIN must be one of the following: " << std::endl;
    std::cout << "  'classical-switch':" << std::endl;
    std::cout << "       Markov chain defined by 'Kannan et al. Simple Markov-chain  algorithms" << std::endl;
    std::cout << "       for generating bipartite graphs and tournaments. Random Structures and" << std::endl;
    std::cout << "       Algorithms 14 (1997), 293–308'." << std::endl;
    std::cout << "  'edge-switch':" << std::endl;
    std::cout << "       Variant of the Markov-chain suggested by 'Kannan et al.' based on an" << std::endl;
    std::cout << "       informed edge selection at the cost of a larger memory consumption." << std::endl;
    std::cout << "  'curveball':" << std::endl;
    std::cout << "       Markov chain defined by 'Strona et al. A fast and unbiased procedure to" << std::endl;
    std::cout << "       randomize ecological binary matrices with fixed row and column totals." << std::endl;
    std::cout << "       Nature communications 5 (2014).'" << std::endl;
    std::cout << std::endl;
    std::cout << "The parameter INSTANCE is a string-encoded input instance of the form \"r*;c*\"." << std::endl;
    std::cout << "While the i-th r defines the sum of row i, the j-th c is the sum of column j." << std::endl;
    std::cout << "For example, the instance \"2,2,2;1,2,1,2\" corresponds to the" << std::endl;
    std::cout << std::endl;
    std::cout << "        row sums:    (2,2,2)" << std::endl;
    std::cout << "        column sums: (1,2,1,2)" << std::endl;
    std::cout << std::endl;
}

#endif //MARATHON_HELPER_H
