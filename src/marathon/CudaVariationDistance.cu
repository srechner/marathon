/*
 * CudaVariationDistance.cu
 *
 * Created on: Mar 23, 2016
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

#include <thrust/reduce.h>
#include <thrust/device_ptr.h>

#include "marathon/TransitionMatrix.h"

namespace marathon {

	namespace cuda {

		template<typename T, uint blockSize>
		__global__ void variationDistanceKernel2(const T *mat, const T *pi, T *d,
		                                         const size_t n, const size_t ld) {

			__shared__ T sdata[blockSize * sizeof(T)];

			uint tid = threadIdx.x;
			uint i = blockIdx.y * blockDim.y + threadIdx.y;
			uint j = threadIdx.x;    // just one block per line

			T sum = 0.0;
			while (j < n) {
				sum += fabs(mat[i * ld + j] - pi[j]);
				j += blockSize;
			}

			sdata[tid] = sum;

			__syncthreads();

			// do reduction in shared mem

			if (blockSize >= 1024) {
				if (tid < 512) {
					sdata[tid] = sum = sum + sdata[tid + 512];
				}

				__syncthreads();
			}

			if (blockSize >= 512) {
				if (tid < 256) {
					sdata[tid] = sum = sum + sdata[tid + 256];
				}

				__syncthreads();
			}

			if (blockSize >= 256) {
				if (tid < 128) {
					sdata[tid] = sum = sum + sdata[tid + 128];
				}

				__syncthreads();
			}

			if (blockSize >= 128) {
				if (tid < 64) {
					sdata[tid] = sum = sum + sdata[tid + 64];
				}

				__syncthreads();
			}

			if (tid < 32) {
				// now that we are using warp-synchronous programming (below)
				// we need to declare our shared memory volatile so that the compiler
				// doesn't reorder stores to it and induce incorrect behavior.
				volatile T *smem = sdata;

				if (blockSize >= 64) {
					smem[tid] = sum = sum + smem[tid + 32];
				}

				if (blockSize >= 32) {
					smem[tid] = sum = sum + smem[tid + 16];
				}

				if (blockSize >= 16) {
					smem[tid] = sum = sum + smem[tid + 8];
				}

				if (blockSize >= 8) {
					smem[tid] = sum = sum + smem[tid + 4];
				}

				if (blockSize >= 4) {
					smem[tid] = sum = sum + smem[tid + 2];
				}

				if (blockSize >= 2) {
					smem[tid] = sum = sum + smem[tid + 1];
				}
			}

			// store result
			if (tid == 0)
				d[i] = sum / 2.0;
		}

		template<typename T>
		extern void cudaVariationDistance(const T *data, const size_t n,
		                                  const size_t ld, const T *pi, T *dist) {

			const uint blockSize = 512;
			const dim3 grid(1, n);
			variationDistanceKernel2<T, blockSize> << < grid, blockSize >> > (data, pi, dist,
					n, ld);
		}

		template<typename T>
		extern T cudaTotalVariationDistance(const T *data, const size_t n,
		                                    const size_t ld, const T *pi) {

			T *dist;
			cudaMalloc(&dist, n * sizeof(T));

			// compute variation distance
			cudaVariationDistance(data, n, ld, pi, dist);

			// use thrust to find maximum
			thrust::device_ptr<T> thrust_ptr(dist);
			const T res = ::thrust::reduce(thrust_ptr, thrust_ptr + n, 0.0,
			                               ::thrust::maximum<T>());

			cudaFree(dist);
			return res;
		}

// export template specializations
		template void cudaVariationDistance<float>(const float *, const size_t,
		                                           const size_t, const float *, float *);

		template void cudaVariationDistance<double>(const double *, const size_t,
		                                            const size_t, const double *, double *);

		template float cudaTotalVariationDistance<float>(const float *, const size_t,
		                                                 const size_t, const float *);

		template double cudaTotalVariationDistance<double>(const double *, const size_t,
		                                                   const size_t, const double *);

	}

}