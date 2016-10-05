/*
 * MixingTime.h
 *
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

#include "StateGraph.h"
#include "PathCongestion.h"
#include "PathConstructionScheme.h"
#include "TransitionMatrixCBLAS.h"
#include "Eigenvalue.h"

#ifdef CUDA
#include "../../include/marathon/TransitionMatrixCuBLAS.h"
#include "../../include/marathon/TransitionMatrixCuBLASXt.h"
#endif

namespace marathon {

	namespace MixingTime {

		/**
		 * Options for computation of total mixing time.
		 */
		enum device_t {
			CPU_ONLY,
#ifdef CUDA
			GPU_ONLY, HYBRID
#endif
		};

		/**
		 * Total mixing time
		 *
		 * The total mixing time t is given by the smallest integer t, such that
		 * dist(P^t, pi) < epsilon, where P is the transition matrix of mc and pi is
		 * its stationary distribution.
		 *
		 * This implementation searches the total mixing time of mc by a two step procedure.
		 * Starting with matrix M=P, it squares M until dist(M,pi) < eps. By doing so,
		 * the methods computes the power of two r, such that dist(P^r,pi) <= eps and
		 * l = r/2 such that dist(P^l, pi) > eps. After finding the limits l and r, it
		 * uses binary search for finding the smallest t such that dist(P^t,pi) < eps.
		 *
		 * In each step of binary search, the boundaries l and r are transformed.
		 * To compute P^m with m=(l+r)/2 we first compute P^(m-l) and multiply the result
		 * to P^l. By this, we need three temporary matrices.
		 *
		 * @param sg A pointer to a state graph object.
		 * @param eps The distance to stationary distribution.
		 * @param dev Determines which matrix multiplication library to use.
		 *    Can be one of the following: CPU_ONLY, GPU_ONLY, HYBRID.
		 */
		template<typename T>
		int totalMixingTime(const StateGraph *sg, const T eps,
		                    device_t device=CPU_ONLY) {

			/* Variables */
			const size_t omega = sg->getNumStates();    // number of states
			TransitionMatrix<T> *P;                     // transition matrix of sg
			TransitionMatrix<T> *tmp[3];                // working matrices
			T *pi;                                      // stationary distribution

			// check trivial cases
			if (omega == 0)
				return -1;
			else if (omega == 1)
				return 0;

			// convert stationary distribution into floating point number array
			pi = new T[omega];
			const rational Z = sg->getZ();
			for (int i = 0; i < omega; i++)
				pi[i] = (sg->getWeight(i) / Z).convert_to<T>();

			// try to allocate memory
			try {
				// decide with mode to use
				switch (device) {

					case device_t::CPU_ONLY:
						// allocate memory
						P = new TransitionMatrixCBLAS<T>(sg);
						tmp[0] = new TransitionMatrixCBLAS<T>(omega);
						tmp[1] = new TransitionMatrixCBLAS<T>(omega);
						tmp[2] = new TransitionMatrixCBLAS<T>(omega);
						break;

#ifdef CUDA
					case device_t::GPU_ONLY:

						// allocate memory
						P = new TransitionMatrixCuBLAS<T>(sg);
						tmp[0] = new TransitionMatrixCuBLAS<T>(omega);
						tmp[1] = new TransitionMatrixCuBLAS<T>(omega);
						tmp[2] = new TransitionMatrixCuBLAS<T>(omega);

						// prepare stationary distribution in device memory
						T *pi_d;
						cudaMalloc((void **) &pi_d, omega * sizeof(T));
						cudaMemcpy(pi_d, pi, omega * sizeof(T), cudaMemcpyHostToDevice);
						std::swap(pi_d, pi);
						delete[] pi_d;
						break;

					case device_t::HYBRID:

						// allocate memory
						P = new TransitionMatrixCuBLASXt<T>(sg);
						tmp[0] = new TransitionMatrixCuBLASXt<T>(omega);
						tmp[1] = new TransitionMatrixCuBLASXt<T>(omega);
						tmp[2] = new TransitionMatrixCuBLASXt<T>(omega);
						break;
#endif

					default:
						std::cerr << "marathon::TotalMixingTime: Error! unknown option: "
						          << device << std::endl;
						return -1;
				}
			}
			catch(std::bad_alloc) {
				// react to bad memory allocation
				return -1;
			}

			// Search for mixing time
			uint l = 0;
			uint r = 1;

			// tmp[0] = P
			tmp[0]->copy(P);

			// First Phase: Square tmp[0] until dist(tmp[0], pi) < eps
			T d = tmp[0]->totalVariationDistance(pi);

			while (d >= eps) {
				// tmp[1] = tmp[0]*tmp[0]
				tmp[1]->mult(tmp[0], tmp[0]);
				tmp[0]->swap(tmp[1]);
				d = tmp[0]->totalVariationDistance(pi);
				//std::cout << "l=" << l << " r=" << r << " d=" << d << std::endl;
				l = r;
				r *= 2;
			}

			/*
			 * State of the variables:
			 *
			 * tmp[0] = P^r
			 * tmp[1] = P^l
			 *
			 * dist_l = dist(tmp[1], pi) <= eps < dist(tmp[0], pi) = dist_r
			 */

			// Second Phase: Binary Search
			// Invariant: tmp[1] = P^l
			while (l < r - 1) {
				uint m = (l + r) / 2;

				// tmp[2] =  P^(m-l)
				tmp[2]->pow(P, m - l);

				// tmp[0] = P^l * P^(m-l) = P^m
				tmp[0]->mult(tmp[1], tmp[2]);
				T dist_m = tmp[0]->totalVariationDistance(pi);

				if (dist_m >= eps) {
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

#ifdef CUDA
			if (device == GPU_ONLY)
				cudaFree(pi);
			else
#endif
				delete[] pi;

			return r;
		}
	}
}


#endif //PROJECT_MIXINGTIME_H
