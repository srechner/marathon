/*
 * TransitionMatrixCPU.h
 *
 * Created on: Mar 22, 2016
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

#ifndef INCLUDE_MARATHON_TRANSITIONMATRIXCUBLASXT_H_
#define INCLUDE_MARATHON_TRANSITIONMATRIXCUBLASXT_H_

#include "TransitionMatrixCBLAS.h"

#ifdef CUDA

#include <cuda_runtime.h>
#include <cublasXt.h>

namespace marathon {

	template<typename T>
	class TransitionMatrixCuBLASXt : public TransitionMatrixCBLAS<T> {

	protected:

		void init_handle() {

			cublasXtCreate(&handle);

			// determine the number of gpu's
			int numDevices;
			cudaGetDeviceCount(&numDevices);

			// use all GPU's
			int devices[numDevices];
			for (int i = 0; i < numDevices; i++)
				devices[i] = i;

			cublasXtDeviceSelect(handle, numDevices, devices);
		}

		// a handle for the cublasXt library
		cublasXtHandle_t handle;

	public:

		TransitionMatrixCuBLASXt(const int n) :
				TransitionMatrixCBLAS<T>(n) {
			init_handle();
		}

		TransitionMatrixCuBLASXt(const StateGraph *sg) :
				TransitionMatrixCBLAS<T>(sg) {
			init_handle();
		}

		virtual ~TransitionMatrixCuBLASXt() {
			cublasXtDestroy(handle);
		}

		/**
		 * Multiply A with B and write the result to this.
		 * @param A A pointer to matrix A. Will not be changed.
		 * @param B A pointer to matrix B. Will not be changed.
		 */
		virtual void mult(const TransitionMatrix <T> *A,
		                  const TransitionMatrix <T> *B);

	};

	template<>
	void TransitionMatrixCuBLASXt<float>::mult(const TransitionMatrix<float> *A,
	                                           const TransitionMatrix<float> *B) {

		const float alpha = 1.0;
		const float beta = 0.0;

		// use cublasXt
		cublasXtSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, this->N, this->N,
		              this->N, &alpha, B->getData(), A->getLeadDimension(), A->getData(),
		              A->getLeadDimension(), &beta, this->getData(),
		              this->getLeadDimension());
	}

	template<>
	void TransitionMatrixCuBLASXt<double>::mult(const TransitionMatrix<double> *A,
	                                            const TransitionMatrix<double> *B) {

		const double alpha = 1.0;
		const double beta = 0.0;

		// use cblas
		cublasXtDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, this->N, this->N,
		              this->N, &alpha, B->getData(), A->getLeadDimension(), A->getData(),
		              A->getLeadDimension(), &beta, this->getData(),
		              this->getLeadDimension());

	}
}

#endif /* CUDA */

#endif /* INCLUDE_MARATHON_TRANSITIONMATRIXCUBLASXT_H_ */
