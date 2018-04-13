/*
 * Created on: Oct 17, 2017
 * Author: Hjalmar Boulouednine <hjalmar@b-nine.de>
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

#ifndef MARATHON_STATISTICS_H
#define MARATHON_STATISTICS_H

#include <unordered_map>

#include "marathon/state.h"
#include "marathon/enumerate.h"
#include "marathon/mixing_time.h"
#include "marathon/random_generator.h"

namespace marathon {

    // create an alias for histograms
    template<typename T=double>
    using SamplingHistogram = std::unordered_map<T, T>;

    /**
     * Calculate the expected value of a given the function calculated
     * for the given MarkovChain by enumeration of all states.
     * @param mc The given MarkovChain object
     * @param f Function that maps a state to a number of type T.
     * @return the mean value computed over all states
     */
    template<typename T>
    T expectedValue(
            Enumerator *enumerator,
            const std::function<T(const State &)> f) {

        T res = 0;
        size_t num_states = 0;
        enumerator->enumerate([&](const State &s) {
            res += f(s);
            ++num_states;
        });
        return res / T(num_states);
    }


    /**
     * Calculate the sample mean of the function f using the random generator
     * to produce N samples.
     * @param rg Random generator.
     * @param f Function mapping each state to a number of type T.
     * @param N Number of samples
     * @return Sample mean computed over N samples.
     */
    template<typename T>
    T sampleMean(
            marathon::RandomGenerator *rg,
            const std::function<T(const State &)> &f,
            const int N) {

        T res = 0;
        for (int l = 0; l < N; l++) {
            const auto &s = rg->next();
            res += f(s);
        }
        return res / T(N);
    }


    /**
     * Calculate the propbability distribution of a function f by enumerating each state.
     * @param enumerator State Enumerator Object.
     * @param f Function that maps states to number of type T.
     * @return Histogram of probability.
     */
    template<typename T=double>
    const SamplingHistogram<T>
    samplingHistogram(
            marathon::Enumerator *enumerator,
            std::function<T(const State *)> f
    ) {

        // create empty sampling histogram
        SamplingHistogram<T> res;

        // enumerate states
        size_t num_states = 0;
        enumerator->enumerate([&](const State *s) {
            res[f(s)]++;
            ++num_states;
        });

        // normalize histogram
        for (auto &p : res) {
            p.second = p.second / T(num_states);
        }

        return res;
    }

    /**
     * Calculate a sampling histogram of a function f by taking N random samples.
     * @param rg Random Generator object.
     * @param f Function that maps states to number of type T.
     * @param N Number of samples.
     * @return Histogram of probability.
     */
    template<typename T=double>
    const SamplingHistogram<T>
    samplingHistogram(
            marathon::RandomGenerator *rg,
            std::function<T(const State *)> f,
            const int N
    ) {

        // create empty sampling histogram
        SamplingHistogram<T> res;

        // draw N random samples
        for (int i = 0; i < N; i++) {
            const auto x = rg->next();
            res[f(x)]++;
        }

        // normalize histogram
        for (auto &p : res) {
            p.second = p.second / T(N);
        }

        return res;
    }

}

#endif /* MARATHON_STATISTICS_H */