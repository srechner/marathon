/*
 * Created on: Mar 06, 2018
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2016, Steffen Rechner
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

#ifndef MARATHON_BINARY_MATRIX_SAMPLING_ENGINE_H
#define MARATHON_BINARY_MATRIX_SAMPLING_ENGINE_H

#include "marathon/binary_matrix/binary_matrix.h"
#include "marathon/binary_matrix/random_generator.h"
#include <vector>
#include <thread>
#include <functional>

namespace marathon {
    namespace binary_matrix {

        /*
         * Abstract base class designed for the parallel construction
         * of large numbers of binary matrices.
         */
        class SamplingEngine {

        private:

            // vector of random generators
            std::vector<std::unique_ptr<marathon::binary_matrix::RandomGenerator>> vec;

        public:

            /**
             * Create an abstract sampling engine.
             */
            SamplingEngine(const marathon::binary_matrix::RandomGenerator &rg) {

                int T = std::max(1u, std::thread::hardware_concurrency());

                // create T random generators
                vec.reserve(T);
                for (int i = 0; i < T; i++) {
                    auto x = static_cast<marathon::binary_matrix::RandomGenerator *>(rg.copy().release());
                    vec.push_back(std::unique_ptr<marathon::binary_matrix::RandomGenerator>(x));
                }
            }

            /**
             * Evaluate the function  a large number of samples
             * @tparam T Return type of function fun.
             * @param number Number of samples.
             * @param fun Function that is evaluated on each sample.
             * @return vector of samples
             */
            template<class T>
            std::vector<T>
            sample(
                    const int number,
                    const std::function<T(const BinaryMatrix &)> fun
            ) {

                std::vector<T> samples(number);

                // vector container stores threads
                std::vector<std::thread> workers;

                const int t = (int) vec.size();

                // launch threads
                for (int i = 0; i < t; i++) {
                    workers.push_back(std::thread([this, t, number, &fun, &samples, i]() {

                        // divide work into chunks of approximately equal size
                        const int begin = i * ((number + t - 1) / t);
                        const int end = std::min(number, (i + 1) * ((number + t - 1) / t));

                        // create random samples
                        for (int j = begin; j < end; j++) {
                            samples[j] = fun(vec[i]->next());
                        }

                    }));
                }

                // wait for threads to finish
                for (int i = 0; i < t; i++)
                    workers[i].join();

                return samples;
            }

            /**
             * Create a large number of random binary matrices.
             * @param number Number of samples to be produced.
             * @@return Vector of binary matrices.
             */
            virtual
            std::vector<BinaryMatrix>
            sample(const int number) {
                return sample<BinaryMatrix>(number, [](const BinaryMatrix &m) {
                    // return the identity of m
                    return m;
                });
            }
        };
    }
}


#endif //MARATHON_BINARY_MATRIX_SAMPLING_ENGINE_H
