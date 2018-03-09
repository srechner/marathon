/*
 * Created on: Jan 09, 2018
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


#ifndef MARATHON_MATCHING_RANDOM_GENERATOR_MCMC_H
#define MARATHON_MATCHING_RANDOM_GENERATOR_MCMC_H

#include "marathon/random_generator.h"
#include "marathon/matching/bipartite_matching.h"
#include "marathon/matching/perfect_matching_jsv04.h"
#include "marathon/matching/perfect_matching_broder86.h"

namespace marathon {
    namespace matching {


        /**
         * A random generator for perfect and near-perfect matchings in a bipartite graph.
         * This random generator produces (approximately) uniform distributed samples.
         */
        class RandomGeneratorMCMC :
                public marathon::RandomGenerator {

        public:

            enum chain_t {
                broder, jsv04
            };

        private:

            MarkovChain *mc;
            int steps;

        public:

            RandomGeneratorMCMC(
                    const SparseBipartiteGraph& g,
                    const chain_t method,
                    const int steps
            ) : steps(steps) {
                switch(method) {
                    case broder:
                        mc = new Broder86(g);
                        break;
                    case jsv04:
                        mc = new JSVChain(g);
                        break;
                }
            }

            /**
             * Create a random generator as the copy of another one.
             * @param rg Random generator.
             */
            RandomGeneratorMCMC(const RandomGeneratorMCMC &rg) :
                    mc(rg.mc->copy()), steps(rg.steps) {

            }

            virtual ~RandomGeneratorMCMC() {
                delete mc;
            }

            /**
             * Return a binary matrix whose row and column sums match the
             * given integer vectors. The sample is assured to be distributed
             * (approximately) uniform.
             * @return Random binary matrix.
             */
            const BipartiteMatching *next() override {

                // apply random walk
                return (const BipartiteMatching *) mc->randomize(steps);
            }

            /**
             * Create an independent copy of the random generator.
             * @return Copy of this random generator.
             */
            RandomGeneratorMCMC *copy() const override {
                return new RandomGeneratorMCMC(*this);
            }
        };
    }
}

#endif //MARATHON_MATCHING_RANDOM_GENERATOR_MCMC_H
