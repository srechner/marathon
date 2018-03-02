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

        const MarkovChain *mc;      // markov chain
        const StateGraph *sg;       // state graph
        const size_t omega;         // number of states
        const T *stationary;        // stationary distribution

        T *initStationary() const {

            // construct stationary distribution as array
            T *pi = new T[omega];
            const Rational Z = sg->getNormalizingConstant();
            for (size_t i = 0; i < omega; i++)
                pi[i] = (sg->getWeight(i) / Z).convert_to<T>();
            return pi;
        }

    public:

        /**
         * Calculate a Mixing Time Calculator.
         * @param sg State Graph.
         */
        explicit MixingTimeCalculator(const StateGraph *sg) :
                mc(sg->getMarkovChain()),
                sg(sg),
                omega(sg->getNumStates()),
                stationary(initStationary()) {

        }

        virtual ~MixingTimeCalculator() {
            delete[] stationary;
        }


        /**
         * Determine the mixing time of state graph sg while starting with a probability vector p_start.
         * @tparam T One of the following: float, double, Rational.
         * @param p_start Initial probability distribution.
         * @param eps Threshold parameter.
         * @return The minimal number of steps t such that ||P^t * p_start - pi|| < eps.
         */
        int mixingTime(const T *p_start, const double eps) {

            // check trivial cases
            if (omega == 0)
                return -1;
            else if (omega == 1)
                return 0;

            // store a single row of the transition matrices P^t and P^t+1
            T *curr = new T[omega];
            T *next = new T[omega];

            // copy the start distribution
            for (size_t i = 0; i < omega; i++)
                curr[i] = p_start[i];

            uint32_t t = 0;     // the exponent of P^t

            // d is variation distance between P^t[i,] and pi
            T d = variationDistance<T>(curr, stationary, omega);

            while (d >= eps) {

                // for each entry of the probability vector
                for (size_t j = 0; j < omega; j++) {

                    // simulate a single step of matrix multiplication
                    T x(0);

                    // for each transition (i,j)
                    for (Transition *kj : sg->getInArcs(j)) {

                        const size_t k = kj->from;
                        const Rational pkj = kj->weight;

                        x += curr[k] * pkj.convert_to<T>();
                    }

                    next[j] = x;
                }

                std::swap(next, curr);

                d = variationDistance<T>(curr, stationary, omega);
                //printf("%i: %f\n", t, d);
                t++;
            }

            delete[] curr;
            delete[] next;

            return t;
        }


        /**
         * Determine the mixing time of state graph sg while starting at state i.
         * @tparam T One of the following: float, double, Rational.
         * @param i Index of initial state.
         * @param eps Threshold parameter.
         * @return The minimal number of steps t such that ||P^t[i,] - pi|| < eps.
         */
        int mixingTime(const int i, const double eps) {

            if (omega == 0 || i < 0 || i >= omega)
                return -1;

            // p_start is initialized as the i-th row of the unit matrix
            T *p_start = new T[omega];
            for (size_t j = 0; j < omega; j++)
                p_start[j] = T(0);
            p_start[i] = T(1);

            int t = mixingTime(p_start, eps);

            delete[] p_start;

            return t;
        }


        /**
         * Determine the maximal mixing time of the states in range [from,to].
         * @param begin Index of first state to be considered in calculations.
         * @param end Index of first state not to be considered in calculations.
         * @param eps Distance to stationary distribution.
         * @return The smallest number of steps t until max_{i=begin..end-1}(||P^t_i-pi|| < eps)
         */
        int mixingTime(const int begin, const int end, const double eps) {

            // check trivial cases
            if (omega == 0)
                return -1;
            else if (omega == 1)
                return 0;

            int t_max = 0;

            // for each row of the transition matrix P
#pragma omp parallel for
            for (int i = std::max(begin, 0); i < std::min(end, (int) omega); i++) {

                // determine mixing time while starting with state i
                int t = mixingTime(i, eps);

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
        int totalMixingTime(const double eps) {
            return mixingTime(0, sg->getNumStates(), eps);
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
        int totalMixingTimeDense(const double eps) {

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

            /* Variables */
            TransitionMatrix<T> *P = nullptr;                     // transition matrix of sg
            TransitionMatrix<T> *tmp[3] = {nullptr, nullptr, nullptr}; // working matrices

            // check trivial cases
            if (omega == 0)
                return -1;
            else if (omega == 1)
                return 0;

            // try to allocate memory
            try {
                // allocate memory
                P = new TransitionMatrix<T>(sg);
                tmp[0] = new TransitionMatrix<T>(omega);
                tmp[1] = new TransitionMatrix<T>(omega);
                tmp[2] = new TransitionMatrix<T>(omega);
            }
            catch (std::bad_alloc) {
                if (P != nullptr)
                    delete P;
                if (tmp[0] != nullptr)
                    delete P;
                if (tmp[1] != nullptr)
                    delete P;
                if (tmp[2] != nullptr)
                    delete P;
                // react to bad memory allocation
                return -1;
            }

            // Search for mixing time
            uint l = 0;
            uint r = 1;

            // tmp[0] = P
            tmp[0]->copy(P);

            // First Phase: Square tmp[0] until dist(tmp[0], pi) < eps
            T d = totalVariationDistance(tmp[0], stationary);

            while (d >= eps) {
                tmp[1]->mult(tmp[0], tmp[0]);           // tmp[1] gets tmp[0]*tmp[0]
                tmp[0]->swap(tmp[1]);                   // tmp[0] gets tmp[1]
                d = totalVariationDistance(tmp[0], stationary);
                //std::cout << "l=" << l << " r=" << r << " d=" << d << std::endl;
                l = r;
                r *= 2;
                if (r >= 100000000) {
                    l = r = UINT_MAX;
                    break;
                }
            }

            /*
             * State of the variables:
             *
             * tmp[0] = P^r
             * tmp[1] = P^l
             *
             * dist(tmp[1], pi) <= eps < dist(tmp[0], pi)
             */

            // Second Phase: Binary Search
            // Invariant: tmp[1] = P^l
            while (l < r - 1) {
                uint m = (l + r) / 2;

                tmp[2]->pow(P, m - l, tmp[0]);  // tmp[2] gets P^(m-l)
                tmp[0]->mult(tmp[1], tmp[2]);   // tmp[0] gets P^l * P^(m-l) = P^m
                d = totalVariationDistance(tmp[0], stationary);

                if (d >= eps) {
                    l = m;
                    tmp[0]->swap(tmp[1]);
                } else {
                    r = m;
                }
            }

            // free memory
            delete P;
            delete tmp[0];
            delete tmp[1];
            delete tmp[2];

            return r;
        }
    };
}


#endif //PROJECT_MIXINGTIME_H
