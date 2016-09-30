/*
 * Eigenvalue.cpp
 *
 * Created on: Nov 23, 2014
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

#ifndef EIGENVALUES_H_
#define EIGENVALUES_H_

// system includes
#include <iostream>
#include <algorithm>

#include "marathon/Eigenvalue.h"
#include "marathon/Memory.h"
#include "marathon/marathon.h"
#include "arpack++/arlssym.h"

namespace marathon {
	namespace Eigenvalue {

/**
 * Computes the second largest value of the spectrum of transition matrix P.
 *
 * Let the sequence of eigenvalues be (1 >= l1 >= l2 >= ... >= lN > -1).
 * The method uses ARPACK++ library to compute the two eigenvalues with
 * largest magnitute and returns the lesser one.
 */
		template<typename T>
		T eigenvalue(const StateGraph *sg, eigenvalue_t opt) {

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
			T * p = alloc<T>(numArcs);
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
			int *row = alloc<int>(omega+1);
			if (row == nullptr) {
				marathon::free(p);
				marathon::free(col);
				throw std::runtime_error(
						"marathon::cpu::eigenvalues::exception: bad memory alloc!");
			}

			// 1. Define Matrix
			ARluSymMatrix<T> B;
			for (i = 0; i < omega + 1; i++)
				row[i] = numArcs;

			// make a copy of the arc pointers
			std::vector<Transition *> copy(sg->getArcs());

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

			ARluSymStdEig<T> prob;
			int nev;
			T tol = std::numeric_limits<T>::epsilon() * 10.0;// one power more than machine precision
			int maxit = 1000 * omega;
			char which[3];

			// decide which eigenvalues to compute
			if (omega == 2) {
				nev = 1;
				strcpy(which, "SM");
			} else if (opt == _2ndLargestMagnitude) {
				nev = 2;
				strcpy(which, "LM");
			} else if (opt == _2ndLargestAlgebraic) {
				nev = 2;
				strcpy(which, "LA");
			} else {
				throw std::runtime_error(
						"marathon::cpu::Eigenvalue::exception: unknown option: "
						+ std::to_string(opt));
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
		 * Export Template Specialization
		 */
		template float eigenvalue<float>(const StateGraph *, eigenvalue_t);

		template double eigenvalue<double>(const StateGraph *, eigenvalue_t);

	}
}

#endif /* EIGENVALUES_H_ */
