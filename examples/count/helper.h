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
    if (strcmp(argv[1], "matching") == 0)
        problem = matching;
    else if (strcmp(argv[1], "fixed") == 0)
        problem = fixed_margin;
    else if (strcmp(argv[1], "interval") == 0)
        problem = interval_margin;
    else {
        std::cerr << "Unknown problem type: " << argv[1] << std::endl;
        return false;
    }

    inst = std::string(argv[2]);

    return true;
}

void print_help_message() {
    std::cout << "Usage: count TYPE INSTANCE" << std::endl;
    std::cout << "Determine the size of the sample space specified by TYPE and INSTANCE." << std::endl;
    std::cout << std::endl;
    std::cout << "The TYPE of sampling problem must be one of the following:\n" << std::endl;
    std::cout << "  'matching':" << std::endl;
    std::cout << "     Count the number of perfect and near perfect matchings in the bipartite" << std::endl;
    std::cout << "     graph that is encoded in the INSTANCE description.\n" << std::endl;
    std::cout << "     An INSTANCE has the form \"([0|1]^n)^n\", where n is a positive integer." << std::endl;
    std::cout << "     Such a 0-1 string is interpreted as the n times n bi-adjacency matrix" << std::endl;
    std::cout << "     M = (m_ij) of a bipartite graph G=(U,V,E) with |U|=|V|, such that\n" << std::endl;
    std::cout << "         m_ij = 1, if (i,j) is in E, or 0, otherwise.\n" << std::endl;
    std::cout << "     For example, the input string  '110101011' corresponds to\n" << std::endl;
    std::cout << "                             u1  u2  u3" << std::endl;
    std::cout << "             1 1 0           |\\ / \\ /|" << std::endl;
    std::cout << "        M =  1 0 1      G =  | X   X |" << std::endl;
    std::cout << "             0 1 1           |/ \\ / \\|" << std::endl;
    std::cout << "                             v1  v2  v3" << std::endl;
    std::cout << "  'fixed':" << std::endl;
    std::cout << "     Count the number of binary matrices whose row and column sums match" << std::endl;
    std::cout << "     prescribed integers.\n\n";
    std::cout << "     An INSTANCE of this TYPE has the form \"r(,r)*;c(,c)*\", where r and c are" << std::endl;
    std::cout << "     positive integers. The i-th integer on the left side of the semicolon" << std::endl;
    std::cout << "     describes the sum of row i. Likewise, the j-th integer on the right side is" << std::endl;
    std::cout << "     the sum of column j.\n" << std::endl;
    std::cout << "     For example, the instance \"2,2,2;1,2,1,2\" corresponds to" << std::endl;
    std::cout << std::endl;
    std::cout << "        row sums:    (2,2,2)" << std::endl;
    std::cout << "        column sums: (1,2,1,2)" << std::endl;
    std::cout << std::endl;
    std::cout << "  'interval':" << std::endl;
    std::cout << "     Count the number of binary matrices whose row and column sums lie in" << std::endl;
    std::cout << "     the intervals that are prescribed by the encoded INSTANCE.\n" << std::endl;
    std::cout << "     The parameter INSTANCE is a string of the form \"l-u(,l-u)*;l-u(,l-u)*\"." << std::endl;
    std::cout << "     The semicolon separates the row sums from the column sums. The i-th" << std::endl;
    std::cout << "     pair on the left side of the semicolon prescribes a lower and upper" << std::endl;
    std::cout << "     bound on the sum of row i. In contrast, the j-th pair on the right side" << std::endl;
    std::cout << "     of the semicolon defines the sum of column j." << std::endl;
    std::cout << "     For convenience, a pair 'l-u' can be replaced by 'l' if l=u.\n" << std::endl;
    std::cout << "     For example, the string \"1-2,2,2-3;0-2,0-1,1-1,1-3\" corresponds to" << std::endl;
    std::cout << std::endl;
    std::cout << "        lower bounds on row sums:    (1,2,2)" << std::endl;
    std::cout << "        upper bounds on row sums:    (2,2,3)" << std::endl;
    std::cout << "        lower bounds on column sums: (0,0,1,1)" << std::endl;
    std::cout << "        upper bounds on column sums: (2,1,1,3)" << std::endl;
    std::cout << std::endl;
}

#endif //MARATHON_HELPER_H
