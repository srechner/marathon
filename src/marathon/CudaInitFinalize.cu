/*
 * CudaInitFinalize.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

#include <cublas_v2.h>
#include <cublasXt.h>

namespace marathon {

// a handle for using cublas library
cublasHandle_t cublasHandle;
cublasXtHandle_t cublasXtHandle;

extern void cudaInit() {

	// check of cuda capable gpu is available
	int deviceCount = 0;
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

	if (error_id != cudaSuccess || deviceCount == 0)
		return;

	cublasCreate_v2(&cublasHandle);

	cublasXtCreate(&cublasXtHandle);

	// use first gpu
	// TODO: support multi gpu usage
	int devices[] = { 0 };
	int numDevices = 1;

	cublasXtDeviceSelect(cublasXtHandle, numDevices, devices);

}

extern void cudaFinalize() {
	cublasDestroy_v2(cublasHandle);
	cublasXtDestroy(cublasXtHandle);
}

}

