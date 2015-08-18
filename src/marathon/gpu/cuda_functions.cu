#include "../../../include/marathon/exceptions.h"

#include <cublas_v2.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>

namespace marathon {

namespace gpu {

namespace cuda {

cublasHandle_t cublasHandle;

extern "C" void initCublas() {
	cublasCreate_v2(&cublasHandle);
}

extern "C" void finalizeCublas() {
	cublasDestroy_v2(cublasHandle);
}

extern "C" void multFloat(const float* A, const size_t ldA, const float* B,
		const size_t ldB, float* C, const size_t ldC, const size_t n) {
	cublasStatus_t err;

	const float alpha_d = 1.0;
	const float beta_d = 0.0;

	err = cublasSgemm_v2(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n,
			&alpha_d, B, ldB, A, ldA, &beta_d, C, ldC);

	if (err != CUBLAS_STATUS_SUCCESS) {
		throw CUBLAS_EXCEPTION;
	}
}

extern "C" void multDouble(const double* A, const size_t ldA, const double* B,
		const size_t ldB, double* C, const size_t ldC, const size_t n) {
	cublasStatus_t err;

	const double alpha_d = 1.0;
	const double beta_d = 0.0;

	err = cublasDgemm_v2(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n,
			&alpha_d, B, ldB, A, ldA, &beta_d, C, ldC);

	if (err != CUBLAS_STATUS_SUCCESS) {
		throw CUBLAS_EXCEPTION;
	}
}

extern "C" void allocMemory(void** ptr, size_t size) {
	cudaError_t err = cudaMalloc(ptr, size);
	if (err != cudaError_t::cudaSuccess)
		throw BAD_CUDA_MEMORY_EXCEPTION;
}

extern "C" void allocMemory2D(void** ptr, size_t* pitch, size_t width,
		size_t height) {
	cudaError_t err = cudaMallocPitch(ptr, pitch, width, height);
	if (err != cudaError_t::cudaSuccess)
		throw BAD_CUDA_MEMORY_EXCEPTION;
}

extern "C" void freeMemory(void* ptr) {
	cudaError_t err = cudaFree(ptr);
	if (err != cudaError_t::cudaSuccess)
		throw BAD_CUDA_MEMORY_EXCEPTION;
}

extern "C" void setMemory2D(void* ptr, size_t pitch, int value, size_t width,
		size_t height) {
	cudaError_t err = cudaMemset2D(ptr, pitch, value, width, height);
	if (err != cudaError_t::cudaSuccess)
		throw BAD_CUDA_MEMORY_EXCEPTION;
}

extern "C" void copy2DDeviceToDevice(void* dst, size_t dpitch, void* src,
		size_t spitch, size_t width, size_t height) {
	cudaError_t err = cudaMemcpy2D(dst, dpitch, src, spitch, width, height,
			cudaMemcpyDeviceToDevice);
	if (err != cudaError_t::cudaSuccess)
		throw BAD_CUDA_MEMORY_EXCEPTION;
}

extern "C" void copy2DDeviceToHost(void* dst, size_t dpitch, void* src,
		size_t spitch, size_t width, size_t height) {
	cudaError_t err = cudaMemcpy2D(dst, dpitch, src, spitch, width, height,
			cudaMemcpyDeviceToHost);
	if (err != cudaError_t::cudaSuccess)
		throw BAD_CUDA_MEMORY_EXCEPTION;
}

extern "C" void copy2DHostToDevice(void* dst, size_t dpitch, void* src,
		size_t spitch, size_t width, size_t height) {
	cudaError_t err = cudaMemcpy2D(dst, dpitch, src, spitch, width, height,
			cudaMemcpyHostToDevice);
	if (err != cudaError_t::cudaSuccess)
		throw BAD_CUDA_MEMORY_EXCEPTION;
}

extern "C" void copyHostToDevice(void* dst, void* src, size_t size) {
	cudaError_t err = cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice);
	if (err != cudaError_t::cudaSuccess)
		throw BAD_CUDA_MEMORY_EXCEPTION;
}

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

extern "C" void variationDistanceFloat(const float* ptr, const float * pi,
		float* dist, const size_t n, const size_t ld) {

	const uint blockSize = 512;
	const dim3 grid(1, n);
	variationDistanceKernel2<float, blockSize> <<<grid, blockSize>>>(ptr, pi,
			dist, n, ld);
}

extern "C" void variationDistanceDouble(const double* ptr, const double * pi,
		double* dist, const size_t n, const size_t ld) {
	const uint blockSize = 512;
	const dim3 grid(1, n);
	variationDistanceKernel2<double, blockSize> <<<grid, blockSize>>>(ptr, pi,
			dist, n, ld);
}

extern "C" float totalVariationDistanceFloat(const float* ptr, const float * pi,
		float* dist, const size_t n, const size_t ld) {

	variationDistanceFloat(ptr, pi, dist, n, ld);
	thrust::device_ptr<float> thrust_ptr(dist);
	return thrust::reduce(thrust_ptr, thrust_ptr + n, (float) 0.0,
			thrust::maximum<float>());
}

extern "C" double totalVariationDistanceDouble(const double* ptr,
		const double * pi, double* dist, const size_t n, const size_t ld) {
	variationDistanceDouble(ptr, pi, dist, n, ld);
	thrust::device_ptr<double> thrust_ptr(dist);
	return thrust::reduce(thrust_ptr, thrust_ptr + n, (double) 0.0,
			thrust::maximum<double>());
}

extern "C" float minVariationDistanceFloat(const float* ptr, const float * pi,
		float* dist, const size_t n, const size_t ld) {

	variationDistanceFloat(ptr, pi, dist, n, ld);
	thrust::device_ptr<float> thrust_ptr(dist);
	return thrust::reduce(thrust_ptr, thrust_ptr + n, (float) 0.0,
			thrust::minimum<float>());
}

extern "C" double minVariationDistanceDouble(const double* ptr,
		const double * pi, double* dist, const size_t n, const size_t ld) {
	variationDistanceDouble(ptr, pi, dist, n, ld);
	thrust::device_ptr<double> thrust_ptr(dist);
	return thrust::reduce(thrust_ptr, thrust_ptr + n, (double) 0.0,
			thrust::minimum<double>());
}

}

}

}
