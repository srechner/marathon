/*
 * CudaWrapper.cu
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

namespace marathon {

extern void myCudaMallocPitch(void** devPtr, size_t* pitch, size_t width,
		size_t height) {
	cudaMallocPitch(devPtr, pitch, width, height);
}

extern void myCudaMalloc(void** devPtr, size_t size) {
	cudaMalloc(devPtr, size);
}

extern void myCudaFree(void* devPtr) {
	cudaFree(devPtr);
}

extern void myCudaMemcpyHostToDevice(void* dst, const void * src,
		size_t count) {
	cudaMemcpy(dst, src, count, cudaMemcpyHostToDevice);
}

extern void myCudaMemset2D(void* devPtr, size_t pitch, int value, size_t width,
		size_t height) {
	cudaMemset2D(devPtr, pitch, value, width, height);
}

extern void myCudaMemcpy2DHostToDevice(void* dst, size_t dpitch,
		const void* src, size_t spitch, size_t width, size_t height) {
	cudaMemcpy2D(dst, dpitch, src, spitch, width, height,
			cudaMemcpyHostToDevice);
}

extern void myCudaMemcpy2DDeviceToHost(void* dst, size_t dpitch,
		const void* src, size_t spitch, size_t width, size_t height) {
	cudaMemcpy2D(dst, dpitch, src, spitch, width, height,
			cudaMemcpyDeviceToHost);
}

extern void myCudaMemcpy2DDeviceToDevice(void* dst, size_t dpitch,
		const void* src, size_t spitch, size_t width, size_t height) {
	cudaMemcpy2D(dst, dpitch, src, spitch, width, height,
			cudaMemcpyDeviceToDevice);
}

}

