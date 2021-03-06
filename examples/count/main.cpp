// system includes
#include <iostream>

// marathon includes
#include "marathon/matching/count.h"
#include "marathon/binary_matrix/fixed_margin/count.h"
#include "marathon/binary_matrix/interval_margin/count.h"

// auxiliary functions
#include "helper.h"

int main(int argc, char **argv) {

    // parse command line arguments
    if (!parse_arguments(argc, argv)) {
        print_help_message();
        return -1;
    }

    // create counter object
    std::unique_ptr<marathon::Counter> cnt;
    switch (problem) {
        case matching:
            cnt = std::make_unique<marathon::matching::Counter>(inst);
            break;
        case fixed_margin:
            cnt = std::make_unique<marathon::binary_matrix::fixed_margin::Counter>(inst);
            break;
        case interval_margin:
            cnt = std::make_unique<marathon::binary_matrix::interval_margin::Counter>(inst);
            break;
    }

    // count the number of objects
    marathon::Integer size = cnt->count();

    // output
    std::cout << size << std::endl;
    
    return 0;
}
