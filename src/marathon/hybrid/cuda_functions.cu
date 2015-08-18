#include "../../../include/marathon/exceptions.h"

#include <cublasXt.h>
#include "stdio.h"

namespace marathon {

namespace hybrid {

namespace cuda {

cublasXtHandle_t cublasXtHandle;

extern "C" void initCublasXt() {
	cublasXtCreate(&cublasXtHandle);

	// use first gpu
	// TODO: support multi gpu usage
	int devices[] = { 0 };
	int numDevices = 1;

	cublasXtDeviceSelect(cublasXtHandle, numDevices, devices);

	// set block dimension
	// TODO: add logic to choose optimal value
	//cublasXtSetBlockDim(cublasXtHandle, 6000);
}

extern "C" void finalizeCublasXt() {
	cublasXtDestroy(cublasXtHandle);
}

extern "C" void multFloatXt(const float* A, const size_t ldA, const float* B,
		const size_t ldB, float* C, const size_t ldC, const size_t n) {
	cublasStatus_t err;

	const float alpha_d = 1.0;
	const float beta_d = 0.0;

	err = cublasXtSgemm(cublasXtHandle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n,
			&alpha_d, B, ldB, A, ldA, &beta_d, C, ldC);

	if (err != CUBLAS_STATUS_SUCCESS) {
		throw CUBLAS_EXCEPTION;
	}
}

extern "C" void multDoubleXt(const double* A, const size_t ldA, const double* B,
		const size_t ldB, double* C, const size_t ldC, const size_t n) {

	cublasStatus_t err;

	const double alpha_d = 1.0;
	const double beta_d = 0.0;

	err = cublasXtDgemm(cublasXtHandle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n,
			&alpha_d, B, ldB, A, ldA, &beta_d, C, ldC);

	if (err != CUBLAS_STATUS_SUCCESS) {
		throw CUBLAS_EXCEPTION;
	}
}

}

}

}
