//
// Created by rechner on 3/1/18.
//

#ifndef MARATHON_HELPER_H
#define MARATHON_HELPER_H

#include <iostream>

// problem type
enum problem_t {
    fixed_margin, interval_margin, matching
} problem;

// string encoded problem instance
std::string inst;

/**
 * Parse command line arguments and extract parameters.
 * @param argc Number of arguments.
 * @param argv Argument list.
 * @return True, if arguments are specified correctly. False, otherwise.
 */
bool parse_arguments(int argc, char **argv) {
    if (argc != 2)
        return false;

    inst = std::string(argv[1]);

    return true;
}

void print_help_message() {
    std::cout << "Usage: enumerate INSTANCE" << std::endl;
    std::cout << "Enumerate the set of binary matrices whose row and column sums lie in the" << std::endl;
    std::cout << "intervals described by the INSTANCE encoding.\n" << std::endl;
    std::cout << "The parameter INSTANCE is a string of the form \"(l-u)*;(l-u)*\"." << std::endl;
    std::cout << "The semicolon separates the row sums from the column sums. The i-th pair on the" << std::endl;
    std::cout << "left side of the semicolon defines lower and upper bounds on the sum of row i." << std::endl;
    std::cout << "In contrast, the j-th pair on the right side of the semicolon defines the sum" << std::endl;
    std::cout << "of column j. For convenience, a token 'l-u' can be replaced by 'l' if l=u." << std::endl;
    std::cout << "For example, the string \"1-2,2,2-3;0-2,0-1,1-1,1-3\" corresponds to" << std::endl;
    std::cout << std::endl;
    std::cout << "        lower bounds on row sums:    (1,2,2)" << std::endl;
    std::cout << "        upper bounds on row sums:    (2,2,3)" << std::endl;
    std::cout << "        lower bounds on column sums: (0,0,1,1)" << std::endl;
    std::cout << "        upper bounds on column sums: (2,1,1,3)" << std::endl;
    std::cout << std::endl;
}


#endif //MARATHON_HELPER_H
