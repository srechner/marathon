/*
 * TransitionMatrixCuBLAS.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

#ifdef CUDA

#include "../../include/marathon/TransitionMatrixCuBLAS.h"

#include <cublas_v2.h>

namespace marathon {

// a handle for using cublas library
extern cublasHandle_t cublasHandle;

template<>
void TransitionMatrixCuBLAS<float>::mult(const TransitionMatrix<float>* A,
		const TransitionMatrix<float>* B) {

	const float alpha = 1.0;
	const float beta = 0.0;

	// use cublas
	cublasSgemm_v2(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, this->n, this->n,
			this->n, &alpha, B->getData(), A->getLeadDimension(), A->getData(),
			A->getLeadDimension(), &beta, this->getData(),
			this->getLeadDimension());
}

template<>
void TransitionMatrixCuBLAS<double>::mult(const TransitionMatrix<double>* A,
		const TransitionMatrix<double>* B) {

	const double alpha = 1.0;
	const double beta = 0.0;

	// use cublas
	cublasDgemm_v2(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, this->n, this->n,
			this->n, &alpha, B->getData(), A->getLeadDimension(), A->getData(),
			A->getLeadDimension(), &beta, this->getData(),
			this->getLeadDimension());

}

}

#endif
