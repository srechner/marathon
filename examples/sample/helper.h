//
// Created by rechner on 3/1/18.
//

#ifndef MARATHON_HELPER_H
#define MARATHON_HELPER_H

#include <iostream>

// sampling algorithms
enum algorithms {
    exact, simple, informed, unknown
} alg = unknown;

int N = -1;                 // number of samples
int steps = 0;              // steps for MCMC methods
std::string inst = "";      // string encoded problem instance

bool parse_arguments(int argc, char **argv) {

    // scan arguments
    for (int i = 1; i < argc; i++) {
        const std::string arg(argv[i]);

        if (arg.compare("-m") == 0) {

            // next argument is the sampling method
            if (i + 1 >= argc) {
                std::cerr << "Error: Missing argument!" << std::endl;
                return false;
            }

            std::string tmp(argv[i + 1]);
            if (tmp.compare("exact") == 0)
                alg = exact;
            else if (tmp.compare("simple") == 0)
                alg = simple;
            else if (tmp.compare("informed") == 0)
                alg = informed;
            else {
                std::cerr << "Error: Unknown method specifier: " << tmp << std::endl;
                return false;
            }

            i++;
        } else if (arg.compare("-N") == 0) {

            // next argument is the number of samples
            if (i + 1 >= argc) {
                std::cerr << "Error: Missing argument!" << std::endl;
                return false;
            }

            N = atoi(argv[i + 1]);

            i++;
        } else if (arg.compare("-s") == 0) {

            // next argument is the number of steps
            if (i + 1 >= argc) {
                std::cerr << "Error: Missing argument!" << std::endl;
                return false;
            }

            steps = atoi(argv[i + 1]);

            i++;
        } else {
            // arg is an instance
            inst = std::string(arg);
        }
    }

    // valid number of steps?
    if (steps <= 0 && alg != exact)
        return false;

    // valid number of samples?
    if (N <= 0)
        return false;

    return true;
}

void print_help_message() {
    std::cout << "Usage: sample METHOD NUMBER [OPTION] INSTANCE" << std::endl;
    std::cout
            << "Construct a NUMBER of uniformly distributed binary matrices whose row and column sums lie in the intervals specified by INSTANCE."
            << std::endl;
    std::cout << std::endl;
    std::cout << "The parameter METHOD must be one of the following:" << std::endl;
    std::cout << std::endl;
    std::cout << "   'exact'   : Use an exact sampling algorithm (not based on MCMC)." << std::endl;
    std::cout << "   'simple'  : Use the first Markov chain defined by " << std::endl;
    std::cout << "               'Rechner et al. Uniform sampling of bipartite graphs with degrees" << std::endl;
    std::cout << "               in prescribed intervals. Journal of Complex Networks (2017)'." << std::endl;
    std::cout << "   'informed': Use the second Markov chain defined by " << std::endl;
    std::cout << "               'Rechner et al. Uniform sampling of bipartite graphs with degrees" << std::endl;
    std::cout << "               in prescribed intervals. Journal of Complex Networks (2017)'." << std::endl;
    std::cout << std::endl;
    std::cout << "If METHOD is 'simple' or 'informed', the optional parameter" << std::endl << std::endl;
    std::cout << "   -s STEPS" << std::endl << std::endl;
    std::cout << "must be used to specify a random walk of length STEPS. (Larger is better.)" << std::endl;
    std::cout << std::endl;
    std::cout << "The parameter INSTANCE is a string of the form \"(l-u)*;(l-u)*\"." << std::endl;
    std::cout << "The semicolon separates the row sums from the column sums. The i-th pair on the" << std::endl;
    std::cout << "left side of the semicolon defines lower and upper bounds on the sum of row i." << std::endl;
    std::cout << "In contrast, the j-th pair on the right side of the semicolon bounds the sum of" << std::endl;
    std::cout << "column j. For convenience, a token 'l-u' can be replaced by 'l' if l=u." << std::endl;
    std::cout << "For example, the string \"1-2,2,2-3;0-2,0-1,1-1,1-3\" corresponds to" << std::endl;
    std::cout << std::endl;
    std::cout << "        lower bounds on row sums:    (1,2,2)" << std::endl;
    std::cout << "        upper bounds on row sums:    (2,2,3)" << std::endl;
    std::cout << "        lower bounds on column sums: (0,0,1,1)" << std::endl;
    std::cout << "        upper bounds on column sums: (2,1,1,3)" << std::endl;
    std::cout << std::endl;
}


#endif //MARATHON_HELPER_H
