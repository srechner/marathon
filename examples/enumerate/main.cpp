// system includes
#include <iostream>

// marathon includes
#include "marathon/matching/enumerate.h"
#include "marathon/binary_matrix/fixed_margin/enumerate.h"
#include "marathon/binary_matrix/interval_margin/enumerate.h"

// auxiliary functions
#include "helper.h"

int main(int argc, char **argv) {

    // parse command line arguments
    if (!parse_arguments(argc, argv)) {
        print_help_message();
        return -1;
    }

    // create enumerator object
    marathon::binary_matrix::interval_margin::Enumerator enm(inst);

    // enumerate states
    size_t num_states = 0;
    enm.enumerate([&](const marathon::State &s) {

        // convert state to binary matrix
        const auto &x = static_cast<const marathon::binary_matrix::BinaryMatrix &> (s);

        // output binary matrix
        std::stringstream ss("A");
        ss << "A" << ++num_states << " = \n\n";
        std::cout << ss.str() << x.fancyString() << std::endl;

    });

    std::cout << "Enumerated " << num_states << " matrices." << std::endl;

    return 0;
}
