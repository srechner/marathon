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

#ifndef MARATHON_EMPIRICAL_MIXING_TIME_H
#define MARATHON_EMPIRICAL_MIXING_TIME_H

#include <unordered_map>

#include "marathon/state.h"
#include "marathon/mixing_time.h"
#include "statistics.h"

namespace marathon {

    template<class T>
    class EmpiricalMixingTimeCalculator {

    protected:

        const std::shared_ptr<MarkovChain> mc;
        const std::function<T(const State &)> f;

    public:

        /**
         * Create a calcuator object.
         * @param mc Markov chain object.
         * @param f Function object evaluated to determine the empirical mixing time.
         */
        EmpiricalMixingTimeCalculator(
                std::shared_ptr<MarkovChain> mc,
                std::function<T(const State &)> f
        ) : mc(std::move(mc)), f(std::move(f)) {

        }

        /**
       * Computes the number of steps until the mean value of the given metric,
       * calculated for N copies of the given chain, is close enough to the given mean
       * value.
       * @param eps The value used to decide wether the two means are close enough.
       * 0.01 was suitebale in most cases.
       * @param mean The given "true" mean value.
       * @param N The number of MarkovChain objects used du compute the mean value of
       * each step.
       * @return The number of steps until the computed mean is close to the given mean.
       */
        int meanConvergenceTime(
                double eps,
                T mean,
                int N) {

            MarkovChain **mc_population = new MarkovChain *[N];
            for (int i = 0; i < N; i++) {
                mc_population[i] = mc->copy();
            }

            const int MAX_STEPS = 10000;
            T val;
            T first_val = 0;
            int res = -1;

            for (int i = 0; i < N; ++i) {
                const State *s = mc_population[i]->getCurrentState();
                first_val += f(s);
            }
            first_val = first_val / N;

            for (int i = 1; i < MAX_STEPS; ++i) {
                for (int i = 0; i < N; ++i) {
                    mc_population[i]->randomize(1);
                    const State *s = (State *) mc_population[i]->getCurrentState();
                    val += f(s);
                }
                // compute the current mean
                val = val / N;
                // normalize the mean using the given "true" mean
                val = (val - mean) / (first_val - mean);
                if (val < eps) {
                    res = i;
                    break;
                }
            }

            // cleanup
            for (int i = 0; i < N; i++)
                delete mc_population[i];
            delete[] mc_population;
            return res;
        }


        /**
         * Computes the number of steps until the variation distance of the distribution
         * of the given metric, calculated for N copies of the given chain, and the
         * given distribution is smaller than eps.
         * @param mc MarkovChain object.
         * @param N Number of samples used to approximate each sampling distribution.
         * @param eps The value used to decide wether the two means are close enough.
         * 0.01 was suitebale in most cases.
         * @param q_limit Limiting distribution of the target function f.
         * @param stepsize Number of steps after which the sampling histograms will be evaluated.
         * @param Enable verbose output.
         * @return Empirical Mixing Time.
         */
        int empiricalMixingTime(
                const std::unordered_map<T, double> q_limit,
                const int N,
                const double eps,
                const int stepsize = 1,
                const int bins = INT_MAX,
                bool verbose = false
        ) {

            // create a population of Markov chain objects
            std::vector<std::unique_ptr<MarkovChain>> mc_population(N);
            for (int i = 0; i < N; i++)
                mc_population[i] = mc->copy();

            // sampling distribution after t steps
            std::unordered_map<T, double> q_t;

            // auxiliary array to store evaluted function values
            std::vector<T> tmp(N);

            // number of steps (will be returned)
            int t = 0;

            // calculate initial variation distance
            T s = f(mc->getCurrentState());
            q_t[s] = 1.0;
            double distance = variationDistance<T, double>(q_t, q_limit, bins);

            /*std::cout << "distribution at step " << t << std::endl;
            for(auto kv : q_t)
                std::cout << kv.first << " " << kv.second << std::endl;
            std::cout << t << "\t" << distance << std::endl;
            std::cout << " -------------- " << std::endl;*/

            // simulate Markov chains
            while (distance > eps) {

#pragma omp parallel for
                for (int i = 0; i < N; i++) {

                    // evolve Markov chain
                    const State &s = mc_population[i]->randomize(stepsize);

                    // evaluate target function
                    tmp[i] = f(s);
                }

                // construct sampling histogram
                q_t.clear();
                for (int i = 0; i < N; ++i)
                    q_t[tmp[i]]++;

                for (auto &kv : q_t)
                    kv.second /= N;

                /*std::cout << "distribution at step " << t << std::endl;
                for(auto kv : q_t)
                    std::cout << kv.first << " " << kv.second << std::endl;
                std::cout << t << "\t" << distance << std::endl;
                std::cout << " -------------- " << std::endl;*/

                if (verbose)
                    std::cout << t << "\t" << distance << std::endl;

                t += stepsize;
                distance = variationDistance<T, double>(q_limit, q_t, bins);
            }

            if (verbose)
                std::cout << t << "\t" << distance << std::endl;

            return t;
        }

    };
}

#endif /* MARATHON_CONVERGENCE_H */