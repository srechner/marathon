// system includes
#include <iostream>

// marathon includes
#include <marathon/binary_matrix/interval_margin/random_generator_exact.h>
#include <marathon/binary_matrix/interval_margin/random_generator_mcmc.h>

// auxiliary functions
#include "helper.h"

int main(int argc, char **argv) {

    // parse command line arguments
    if (!parse_arguments(argc, argv)) {
        print_help_message();
        return -1;
    }

    // decode instance description
    marathon::binary_matrix::interval_margin::Instance margin(inst);

    // create random generator
    std::unique_ptr<marathon::binary_matrix::RandomGenerator> rg;

    switch (alg) {
        case exact: {
            rg = std::make_unique<marathon::binary_matrix::interval_margin::RandomGeneratorExact>(margin);
            break;
        }
        case simple: {
            auto type = marathon::binary_matrix::interval_margin::RandomGeneratorMCMC::simple;
            rg = std::make_unique<marathon::binary_matrix::interval_margin::RandomGeneratorMCMC>(margin, type, steps);
            break;
        }
        case informed: {
            auto type = marathon::binary_matrix::interval_margin::RandomGeneratorMCMC::informed;
            rg = std::make_unique<marathon::binary_matrix::interval_margin::RandomGeneratorMCMC>(margin, type, steps);
            break;
        }
    }

    // create N random samples
    for (int i = 0; i < N; i++) {

        // create a random binary matrix
        const marathon::binary_matrix::BinaryMatrix &s = rg->next();

        // process sample...

        // output string representation
        std::stringstream ss;
        ss << "A" << i << " = \n\n";
        std::cout << ss.str() << s.fancyString() << std::endl;
    }

    return 0;
}
