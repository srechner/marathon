/*
 * CudaVariationDistance.cu
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

#include <thrust/reduce.h>
#include <thrust/device_ptr.h>

#include "../../include/marathon/TransitionMatrix.h"

namespace marathon {

template<typename T, uint blockSize>
__global__ void variationDistanceKernel2(const T* mat, const T * pi, T* d,
		const size_t n, const size_t ld) {

	__shared__ T sdata[blockSize * sizeof(T)];

	uint tid = threadIdx.x;
	uint i = blockIdx.y * blockDim.y + threadIdx.y;
	uint j = threadIdx.x;	// just one block per line

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
extern void cudaVariationDistance(const T* data, const size_t n,
		const size_t ld, const T * pi, T* dist) {

	const uint blockSize = 512;
	const dim3 grid(1, n);
	variationDistanceKernel2<T, blockSize> <<<grid, blockSize>>>(data, pi, dist,
			n, ld);
}

template<typename T>
extern T cudaTotalVariationDistance(const T* data, const size_t n,
		const size_t ld, const T * pi) {

	T* dist;
	cudaMalloc(&dist, n * sizeof(T));

	// compute variation distance
	cudaVariationDistance(data, n, ld, pi, dist);

	// use thrust to find maximum
	thrust::device_ptr < T > thrust_ptr(dist);
	const T res = ::thrust::reduce(thrust_ptr, thrust_ptr + n, 0.0,
			::thrust::maximum<T>());

	cudaFree(dist);
	return res;
}

// export template specializations
template void cudaVariationDistance<float>(const float*, const size_t,
		const size_t, const float *, float*);
template void cudaVariationDistance<double>(const double*, const size_t,
		const size_t, const double *, double*);
template float cudaTotalVariationDistance<float>(const float*, const size_t,
		const size_t, const float *);
template double cudaTotalVariationDistance<double>(const double*, const size_t,
		const size_t, const double *);
}
