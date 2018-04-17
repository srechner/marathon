/*
 * Created on: May 30, 2017
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


#ifndef MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_RANDOM_GENERATOR_MCMC_H
#define MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_RANDOM_GENERATOR_MCMC_H

#include "marathon/binary_matrix/binary_matrix.h"
#include "marathon/binary_matrix/random_generator.h"
#include "marathon/binary_matrix/interval_margin/simple_chain.h"
#include "marathon/binary_matrix/interval_margin/simple_chain_dynamic.h"
#include "marathon/binary_matrix/interval_margin/informed_chain_dynamic.h"
#include "marathon/binary_matrix/interval_margin/informed_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {


            /**
             * A random generator for binary matrices whose row sums and column sums
             * lie in prescribed intervals. This random generator produces (approximately)
             * uniform distributed samples.
             */
            class RandomGeneratorMCMC :
                    public marathon::binary_matrix::RandomGenerator {

            public:

                /**
                 * Markov chain specifier.
                 */
                enum chain_t {
                    simple, informed, simple_dynamic, informed_dynamic, unknown
                };

            private:

                size_t _steps;
                std::unique_ptr<interval_margin::MarkovChain> _mc;

            public:


                /**
                 * Create a generator for binary matrices with prescribed row sums and column sums.
                 * This generator uses a Markov chain Monte Carlo approach to create random samples.
                 * @param seq Four-Tuple of integer vectors.
                 * @param chain Markov chain specifier.
                 * @param steps Length of the random walk.
                 */
                RandomGeneratorMCMC(
                        Instance inst,
                        chain_t chain,
                        size_t steps
                ) : _steps(steps) {

                    switch (chain) {
                        case simple:
                            _mc = std::make_unique<SimpleChain>(std::move(inst));
                            break;
                        case informed:
                            _mc = std::make_unique<InformedChain>(std::move(inst));
                            break;
                        case simple_dynamic:
                            _mc = std::make_unique<SimpleChainDynamic>(std::move(inst));
                            break;
                        case informed_dynamic:
                            _mc = std::make_unique<InformedChainDynamic>(std::move(inst));
                            break;
                        default:
                            throw std::runtime_error("Error! unknown chain specifier!");
                    }
                }

                /**
                 * Create a random generator with a predefined Markov chain.
                 * @param mc Markov chain.
                 * @param steps Number of steps.
                 */
                RandomGeneratorMCMC(
                        const std::unique_ptr<interval_margin::MarkovChain> &mc,
                        size_t steps)
                        : _steps(steps) {

                    auto x = static_cast<interval_margin::MarkovChain *>(mc->copy().release());
                    _mc = std::unique_ptr<interval_margin::MarkovChain>(x);
                }


                /**
                 * Return a binary matrix whose row and column sums match the
                 * given integer vectors. The sample is assured to be distributed
                 * (approximately) uniform.
                 * @return Random binary matrix.
                 */
                const BinaryMatrix &next() override {

                    // apply random walk
                    return static_cast<const BinaryMatrix &>(_mc->randomize(_steps));
                }

                /**
                 * Create an independent copy of the random generator.
                 * @return Copy of this random generator.
                 */
                std::unique_ptr<marathon::RandomGenerator> copy() const override {
                    return std::make_unique<RandomGeneratorMCMC>(_mc, _steps);
                }
            };
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_RANDOM_GENERATOR_MCMC_H
