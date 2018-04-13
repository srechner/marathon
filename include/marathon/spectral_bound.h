/*
 * Created on: Mar 24, 2016
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

#ifndef INCLUDE_MARATHON_EIGENVALUES_H_
#define INCLUDE_MARATHON_EIGENVALUES_H_

#ifdef USE_ARMADILLO


#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include <iostream>
#include <algorithm>

#include "state_graph.h"

namespace marathon {

    template<typename T>
    class SpectralBoundCalculator {

    protected:

        const StateGraph &sg;
        T lambda_max = -1;
        Rational pimin;

    public:

        /**
         * Create a Spectral Bound Calculator.
         * @param sg
         */
        SpectralBoundCalculator(const StateGraph &sg) : sg(sg) {
            pimin = sg.getMinWeight() / sg.getNormalizingConstant();
        }


        /**
         * Computes the eigenvalue with second largest magnitute of the
         * transition matrix of mc.
         * @param mc State Graph Representation of Markov Chain.
         * @param which_eig Which Eigenvalue to compute. Options are: 2nd largest in magnitute and 2nd largest algebraic.
         * @return The corresponding Eigenvalue.
         */
        T eigenvalue() {

            size_t omega = sg.getNumStates();          // Number of states
            size_t numArcs = sg.getNumTransitions();   // Number of transitions

            // Check trivial cases
            if (omega <= 1)
                return 1;

            arma::umat locations(2, numArcs);
            arma::Col<T> values(numArcs);

            int i = 0;
            for (const Transition *t : sg.getArcs()) {

                Rational wu = sg.getWeight(t->from);
                Rational wv = sg.getWeight(t->to);
                T p = t->weight.convert_to<T>();

                // if transition matrix is not symmetric
                if(wu != wv) {
                    // symmetrize (see theory paper)
                    p *= sqrt((wu / wv).convert_to<T>());
                }

                locations(0, i) = t->from;
                locations(1, i) = t->to;
                values(i) = p;
                i++;
            }

            arma::SpMat<T> sp_A(locations, values);
            arma::Col<T> eigval; // eigenvalues

            T lambda;

            // compute eigenvalues
            if(omega == 2) {

                if(!arma::eigs_sym(eigval, sp_A, 1, "sm") || eigval.size() != 1)
                    throw std::runtime_error("Error while computing eigenvalues: Could not determine smallest eigenvalue!");

                lambda = eigval[0];
            }
            else {

                int k = 2;  // how many eigenvalues
                while (!arma::eigs_sym(eigval, sp_A, k) || eigval.size() != k) {
                    k++;
                }

                std::sort(eigval.begin(), eigval.end(), [](double x, double y)
                          { return std::abs(x) < std::abs(y); }
                );

                // eigval.size() - 1 is the position of the largest
                // eigval.size() - 2 is the position of the second largest
                lambda = eigval[eigval.size() - 2];
            }


            //if (fabs(lambda) < std::numeric_limits<T>::epsilon())
            //	lambda = 0.0;

            return lambda;
        }


        /**
         * Let 1 = x1 > x2 >= x_3 >= ... >= x_N > -1 be the real eigenvalues
         * of a transition matrix P.
         * Return max(x_2, x_N).
         * @return
         */
        T getLambdaMax() {

            if (lambda_max == -1)
                lambda_max = eigenvalue();

            return lambda_max;
        }



        /**
         * Compute the lower spectral bound of the state graph.
         * @param eps The distance to stationary distribution.
         * @return A lower bound of the total mixing time.
         */
        T lowerSpectralBound(double eps) {


            if (lambda_max == -1)
                lambda_max = eigenvalue();

            const double lambda = fabs(lambda_max);

            return 0.5 * (lambda / (1.0 - lambda)) * -log(2.0 * eps);
        }


        /**
         * Compute the upper spectral bound of the state graph.
         * @param eps The distance to stationary distribution.
         * @return A lower bound of the total mixing time.
         */
        T upperSpectralBound(double eps) {

            if (lambda_max == -1)
                lambda_max = eigenvalue();

            const double lambda = fabs(lambda_max);

            return -log(eps * pimin.convert_to<double>()) / (1.0 - lambda);
        }

    };
}

#endif

#endif /* INCLUDE_MARATHON_EIGENVALUES_H_ */
