/*
 * Eigenvalues.h
 *
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

#include "StateGraph.h"
#include "marathon/Memory.h"
#include "marathon/marathon.h"
#include "arpack++/arlssym.h"

#include <iostream>
#include <algorithm>


namespace marathon {

	namespace Eigenvalue {

		/**
		 * Options for computation of eigenvalues.
		 */
		enum eigenvalue_t {
			// Eigenvalue options
					_2ndLargestMagnitude,
			_2ndLargestAlgebraic,
		};

		/**
		 * Computes the eigenvalue with second largest magnitute of the
		 * transition matrix of mc.
		 * @param mc State Graph Representation of Markov Chain.
		 * @param which_eig Which Eigenvalue to compute. Options are: 2nd largest in magnitute and 2nd largest algebraic.
		 * @return The corresponding Eigenvalue.
		 */
		template<typename T>
		T eigenvalue(const StateGraph *sg, eigenvalue_t which_eig) {

			// Variables
			size_t omega;            // Number of states
			size_t numArcs;            // Number of transitions
			unsigned int i;

			omega = sg->getNumStates();

			// Check trivial cases
			if (omega <= 1)
				return 1;

			// count_recursive number of non-zero elements in left upper part of transition matrix
			numArcs = 0;
			for (Transition *t : sg->getArcs()) {
				if (t->u <= t->v)
					numArcs++;
			}

			// 0. Try to alloate memory
			T *p = alloc<T>(numArcs);
			if (p == nullptr) {
				throw std::runtime_error(
						"marathon::cpu::eigenvalues::exception: bad memory alloc!");
			}
			int *col = alloc<int>(numArcs);
			if (col == nullptr) {
				marathon::free(p);
				throw std::runtime_error(
						"marathon::cpu::eigenvalues::exception: bad memory alloc!");
			}
			int *row = alloc<int>(omega + 1);
			if (row == nullptr) {
				marathon::free(p);
				marathon::free(col);
				throw std::runtime_error(
						"marathon::cpu::eigenvalues::exception: bad memory alloc!");
			}

			// 1. Define Matrix
			ARluSymMatrix <T> B;
			for (i = 0; i < omega + 1; i++)
				row[i] = numArcs;

			// make a copy of the arc pointers
			std::vector < Transition * > copy(sg->getArcs());

			// sort transitions ascending
			std::sort(copy.begin(), copy.end(),
			          [](const Transition *a, const Transition *b) -> bool {
				          if (a->u == b->u)
					          return a->v < b->v;
				          else
					          return a->u < b->u;
			          });
			i = 0;
			for (const Transition *t : copy) {
				if (t->u <= t->v) {
					// symmetrize (see theory paper)
					T x = t->p.convert_to<T>();
					rational y = sg->getWeight(t->u) / sg->getWeight(t->v);
					x *= sqrt(y.convert_to<T>());
					p[i] = x;
					col[i] = t->v;
					if (i < row[t->u])
						row[t->u] = i;
					i++;
				}
			}

			/*
			std::cout << "P: ";
			for (i = 0; i < numArcs; i++)
				std::cout << p[i] << " ";
			std::cout << std::endl;

			std::cout << "row: ";
			for (i = 0; i <= omega; i++)
				std::cout << row[i] << " ";
			std::cout << std::endl;

			std::cout << "col: ";
			for (i = 0; i < numArcs; i++)
				std::cout << col[i] << " ";
			std::cout << std::endl;
			*/

			B.DefineMatrix(omega, numArcs, p, col, row, 'L');

			// 2. Define Problem

			// Symmetric Standard Eigenvalue Problem
			// Compute the two eigenvalues with largest magnitude

			ARluSymStdEig <T> prob;
			int nev;
			T tol = std::numeric_limits<T>::epsilon() * 10.0;// one power more than machine precision
			int maxit = 1000 * omega;
			char which[3];

			// decide which eigenvalues to compute
			if (omega == 2) {
				nev = 1;
				strcpy(which, "SM");
			} else if (which_eig == _2ndLargestMagnitude) {
				nev = 2;
				strcpy(which, "LM");
			} else if (which_eig == _2ndLargestAlgebraic) {
				nev = 2;
				strcpy(which, "LA");
			} else {
				throw std::runtime_error(
						"marathon::cpu::Eigenvalue::exception: unknown option: "
						+ std::to_string(which_eig));
			}

			int ncv = std::min(2 * nev, (int) omega);

			/*std::cout << "omega = " << omega << std::endl;
			std::cout << "nev   = " << nev << std::endl;
			std::cout << "ncv   = " << ncv << std::endl;
			std::cout << "maxit = " << maxit << std::endl;
			std::cout << "tol   = " << tol << std::endl;*/

			prob.DefineParameters(omega, nev, &B, &ARluSymMatrix<T>::MultMv, which, ncv, tol, maxit);

			// 3. Solve Problem
			prob.FindEigenvalues();

			/*std::cout << "found " << prob.ConvergedEigenvalues() << std::endl;
			for (int i = 0; i < prob.ConvergedEigenvalues(); i++) {
				std::cout << i << ": " << std::setprecision(std::numeric_limits<T>::digits10) << prob.Eigenvalue(i) << std::endl;
			}*/

			// free memory
			marathon::free(p);
			marathon::free(col);
			marathon::free(row);

			return prob.Eigenvalue(0);
		}

		/**
		 * Compute the lower spectral bound of the state graph.
		 * @param sg A pointer to a state graph object.
		 * @param eps The distance to stationary distribution.
		 * @return A lower bound of the total mixing time.
		 */
		template<typename T>
		T lowerSpectralBound(const StateGraph *sg, T eps) {

			const double lambda = fabs(
					::marathon::Eigenvalue::eigenvalue<T>(sg,
					                                      ::marathon::Eigenvalue::eigenvalue_t::_2ndLargestMagnitude));

			//if (lambda < std::numeric_limits<T>::epsilon())
			//	lambda = 0.0;

			return 0.5 * (lambda / (1.0 - lambda)) * -log(2.0 * eps);
		}


		/**
		 * Compute the upper spectral bound of the state graph.
		 * @param sg A pointer to a state graph object.
		 * @param eps The distance to stationary distribution.
		 * @return A lower bound of the total mixing time.
		 */
		template<typename T>
		T upperSpectralBound(const StateGraph *sg, T eps) {

			const double lambda = fabs(
					marathon::Eigenvalue::eigenvalue<T>(sg,
					                                    ::marathon::Eigenvalue::eigenvalue_t::_2ndLargestMagnitude));

			//std::cout << lambda << std::endl;

			//if (fabs(lambda) < std::numeric_limits<T>::epsilon())
			//	return std::numeric_limits<T>::infinity();

			const rational pimin = sg->getMinWeight() / sg->getZ();

			return -log(eps * pimin.convert_to<double>()) / (1.0 - lambda);
		}
	}
}

#endif /* INCLUDE_MARATHON_EIGENVALUES_H_ */
