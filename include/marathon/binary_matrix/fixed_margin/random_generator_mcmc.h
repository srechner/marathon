/*
 * Created on: May 19, 2017
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2017, Steffen Rechner
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


#ifndef MARATHON_BINARY_MATRIX_FIXED_MARGIN_RANDOM_GENERATOR_MCMC_H
#define MARATHON_BINARY_MATRIX_FIXED_MARGIN_RANDOM_GENERATOR_MCMC_H

#include "marathon/binary_matrix/binary_matrix.h"
#include "marathon/binary_matrix/random_generator.h"
#include "marathon/binary_matrix/fixed_margin/realize.h"
#include "marathon/binary_matrix/fixed_margin/decompose.h"
#include "marathon/binary_matrix/fixed_margin/switch_chain.h"
#include "marathon/binary_matrix/fixed_margin/edge_switch_chain.h"
#include "marathon/binary_matrix/fixed_margin/curveball.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * A generator for binary matrices with prescribed row sums
             * and column sums.
             * This generator produces (approximately) uniformly
             * distributed samples.
             * The generator does not use sequence decomposition.
             */
            class RandomGeneratorMCMC :
                    public marathon::binary_matrix::RandomGenerator {

            private:

                // Markov chain
                MarkovChain *mc;
                int steps;

            public:

                /**
                 * Markov chain specifier.
                 */
                enum chain_t {
                    classical_switch, edge_switch, curveball
                };


                /**
                 * Create a generator for binary matrices with prescribed row sums and column sums.
                 * This generator uses a Markov chain Monte Carlo approach to create random samples.
                 * @param margin Row and column sums.
                 * @param chain Markov chain specifier.
                 * @param steps Length of the random walk
                 */
                RandomGeneratorMCMC(
                        const Instance &margin,
                        const chain_t chain,
                        const int steps
                ) : steps(steps) {

                    // create Markov chain object
                    switch(chain) {
                        case classical_switch:
                            mc = new SwitchChain(margin);
                            break;
                        case edge_switch:
                            mc = new EdgeSwitchChain(margin);
                            break;
                        case curveball:
                            mc = new Curveball(margin);
                            break;
                    }
                }


                ~RandomGeneratorMCMC() {
                    delete mc;
                }

                /**
                 * Return a binary matrix of size nrow times ncol whose
                 * row and column sums match the given integer vectors.
                 * The sample is distributed (approximately) uniform from the set
                 * of all binary matrices with the prescribed row sums and
                 * column sums.
                 * @return Random binary matrix.
                 */
                const BinaryMatrix *next() override {

                    // apply random walk of length <steps>
                    return mc->randomize(steps);
                }
            };
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_FIXED_MARGIN_RANDOM_GENERATOR_MCMC_H
