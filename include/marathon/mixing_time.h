/*
 * Created on: Sep 23, 2016
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

#ifndef PROJECT_MIXINGTIME_H
#define PROJECT_MIXINGTIME_H

#include "variation_distance.h"
#include "state_graph.h"
#include "transition_matrix.h"

namespace marathon {

    /**
     * A class for the computation of total mixing time and its bounds.
     * @tparam T
     */
    template<class T=double>
    class MixingTimeCalculator {

    protected:

        const MarkovChain &mc;              // markov chain
        const StateGraph &sg;               // state graph
        const size_t omega;                 // number of states
        const std::vector<T> stationary;    // stationary distribution

        std::vector<T> initStationary() const {

            // construct stationary distribution as array
            std::vector<T> pi(omega);
            const Rational Z = sg.getNormalizingConstant();
            for (size_t i = 0; i < omega; i++)
                pi[i] = (sg.getWeight(i) / Z).convert_to<T>();
            return pi;
        }

    public:

        /**
         * Calculate a Mixing Time Calculator.
         * @param sg State Graph.
         */
        explicit MixingTimeCalculator(const StateGraph &sg) :
                mc(sg.getMarkovChain()),
                sg(sg),
                omega(sg.getNumStates()),
                stationary(initStationary()) {

        }

        /**
         * Determine the mixing time of state graph sg while starting with a probability vector p_start.
         * @tparam T One of the following: float, double, Rational.
         * @param p_start Initial probability distribution.
         * @param eps Threshold parameter.
         * @return The minimal number of steps t such that ||P^t * p_start - pi|| < eps.
         */
        uint mixingTime(const std::vector<T> &p_start, double eps) {

            // check trivial cases
            if (omega == 0)
                throw std::runtime_error("Error! Empty state space!");
            else if (omega == 1)
                return 0;

            // store a single row of the transition matrices P^t and P^t+1
            std::vector<T> curr(p_start);
            std::vector<T> next(omega);

            uint32_t t = 0;     // the exponent of P^t

            // d is variation distance between P^t[i,] and pi
            T d = variationDistance<T>(curr, stationary);

            while (d >= eps) {

                // for each entry of the probability vector
                for (size_t j = 0; j < omega; j++) {

                    // simulate a single step of matrix multiplication
                    T x(0);

                    // for each transition (i,j)
                    for (Transition *kj : sg.getInArcs(j)) {

                        const size_t k = kj->from;
                        const Rational pkj = kj->weight;

                        x += curr[k] * pkj.convert_to<T>();
                    }

                    next[j] = x;
                }

                std::swap(next, curr);

                d = variationDistance<T>(curr, stationary);
                //printf("%i: %f\n", t, d);
                t++;
            }

            return t;
        }


        /**
         * Determine the mixing time of state graph sg while starting at state i.
         * @tparam T One of the following: float, double, Rational.
         * @param i Index of initial state.
         * @param eps Threshold parameter.
         * @return The minimal number of steps t such that ||P^t[i,] - pi|| < eps.
         */
        uint mixingTime(const size_t i, const double eps) {

            if (omega == 0 || i < 0 || i >= omega)
                throw std::runtime_error("Error! Invalid state index!");

            // p_start is initialized as the i-th row of the unit matrix
            std::vector<T> p_start(omega, 0);
            p_start[i] = T(1);

            return mixingTime(p_start, eps);
        }


        /**
         * Determine the maximal mixing time of the states in range [from,to].
         * @param begin Index of first state to be considered in calculations.
         * @param end Index of first state not to be considered in calculations.
         * @param eps Distance to stationary distribution.
         * @return The smallest number of steps t until max_{i=begin..end-1}(||P^t_i-pi|| < eps)
         */
        uint mixingTime(size_t begin, size_t end, double eps) {

            // check trivial cases
            if (omega == 0)
                throw std::runtime_error("Error! State space is empty!");
            else if (omega == 1)
                return 0;

            uint t_max = 0;

            // for each row of the transition matrix P
#pragma omp parallel for
            for (size_t i = std::max(begin, 0lu); i < std::min(end, omega); i++) {

                // determine mixing time while starting with state i
                uint t = mixingTime(i, eps);

#pragma omp critical
                t_max = std::max(t, t_max);
            }

            return t_max;
        }


        /**
         * Determine the total mixing time t of a state graph.
         * The total mixing time is defined as the smallest integer t,
         * such that max_{i=1..V} ||P^t[i,] - pi|| < eps, where P is
         * the transition matrix of a Markov chain, V is its number of
         * states and pi is its stationary distribution.
         * Running Time: O(t*V^3), where V is the number of states.
         * Memory Complexity: O(V).
         * @param eps Distance to stationary distribution.
         * @return The smallest number of steps t until max_{i=1..omega}(||P^t_i-pi|| < eps)
         */
        uint totalMixingTime(const double eps) {
            return mixingTime(0, sg.getNumStates(), eps);
        }


        /**
         * Determines the total mixing time t of a state graph.
         * The total mixing time is defined as the smallest integer t,
         * such that max_{i=1..V} ||P^t[i,] - pi|| < eps, where P is
         * the transition matrix of a Markov chain, V is its number of
         * states and pi is its stationary distribution.
         * Running Time: O(t*V^3), where V is the number of states.
         * Memory Complexity: O(V^2), where V is the number of states.
         * @tparam T One of the following: float, double, Rational.
         * @param sg Pointer to a state graph object.
         * @param eps The distance to stationary distribution.
         * @return
         */
        uint totalMixingTimeDense(const double eps) {

            /**********************************************************************
             * This implementation searches the total mixing time of mc by a two
             * step procedure. Starting with matrix M=P, it squares M until
             * dist(M,pi) < eps. By doing so, the method computes the power of
             * two r, such that dist(P^r,pi) <= eps and l = r/2 such that
             * dist(P^l, pi) > eps.
             *
             * After finding the limits l and r, it uses binary search for finding
             * the smallest t such that dist(P^t,pi) < eps.
             *
             * In each step of binary search, the boundaries l and r are
             * transformed. To compute P^m with m=(l+r)/2 we first compute P^(m-l)
             * and multiply the result to P^l. By this, we need three temporary
             * matrices.
             *********************************************************************/

            // check trivial cases
            if (omega == 0)
                throw std::runtime_error("Error! Empty state space!");
            else if (omega == 1)
                return 0;

            const TransitionMatrix<T> P(sg);    // transition matrix

            // Phase 1: Search for smallest k such that ||P^(2^k) - pi|| < eps
            uint k = 0;
            TransitionMatrix<T> A(P);           // invariant A = P^(2^k)

            while (totalVariationDistance(A, stationary) >= eps) {
                A = A * A;                      // A = P^(2^(k+1))
                k++;                            // invariant restored
            }

            if (k == 0)
                return 1;

            // it holds: ||P^(2^k) - pi|| < eps <= ||P^(2^(k-1)) - pi||

            uint l = 1u << (k - 1);       // l = 2^(k-1)
            uint u = 1u << k;             // u = 2^k

            // Phase 2: Binary search

            // invariant: ||P^u - pi|| < eps <= ||P^l - pi||
            while (u - l > 1) {

                uint m = (l + u) / 2;

                A = P.pow(m);                   // A = P^m

                T d = totalVariationDistance(A, stationary);

                if (d < eps) {
                    u = m;
                } else {
                    l = m;
                }
            }

            return u;
        }
    };
}


#endif //PROJECT_MIXINGTIME_H
