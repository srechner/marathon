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

#ifndef INCLUDE_MARATHON_TRANSITIONMATRIXCUBLAS_H_
#define INCLUDE_MARATHON_TRANSITIONMATRIXCUBLAS_H_

#include "TransitionMatrixCBLAS.h"

#ifdef CUDA

#include <cuda_runtime.h>
#include <cublas_v2.h>

namespace marathon {

	/**
	 * Declare external functions.
	 */
	namespace cuda {

		template<typename T>
		extern void cudaVariationDistance(const T *data, const size_t n,
		                                  const size_t ld, const T *pi, T *dist);

		template<typename T>
		extern T cudaTotalVariationDistance(const T *data, const size_t n,
		                                    const size_t ld, const T *pi);
	}

	template<typename T>
	class TransitionMatrixCuBLAS : public TransitionMatrix<T> {

	protected:

		/**
		 * Member Variables.
		 */
		cublasHandle_t handle;      // for matrix multiplication

		/**
		 * Copy the content of matrix P to this.
		 */
		virtual void copy(const TransitionMatrix<T> *P) {

			const TransitionMatrixCuBLAS<T> *X =
					(const TransitionMatrixCuBLAS<T> *) P;
			this->N = X->N;
			this->ld = X->ld;
			this->pitch = X->pitch;
			cudaMemcpy2D(this->data, this->ld * sizeof(T),
			             X->data, X->ld * sizeof(T),
			             this->N * sizeof(T), this->N,
			             cudaMemcpyDeviceToDevice);
		}

		/**
		 * Return a pointer to a Transition matrix of subtype instance.
		 */
		virtual TransitionMatrix<T> *generateSubTypeInstance(const int n) {
			TransitionMatrix<T> *res = new TransitionMatrixCuBLAS<T>(n);
			return res;
		}

	public:

		TransitionMatrixCuBLAS(const uint32_t N) {

			cublasCreate_v2(&handle);

			this->N = N;
			// allocate aligned 2d memory
			cudaMallocPitch((void**) &this->data, &this->pitch, N * sizeof(T), N);
			this->ld = (uint32_t) (this->pitch / sizeof(T));
		}

		TransitionMatrixCuBLAS(const StateGraph *sg) :
				TransitionMatrixCuBLAS<T>(sg->getNumStates()) {

			TransitionMatrixCBLAS<T> tmp(sg);

			cudaMemcpy2D(this->data, this->pitch,
			             tmp.getData(), tmp.getLeadDimension()*sizeof(T),
			             this->N * sizeof(T), this->N,
			             cudaMemcpyHostToDevice);
		}

		virtual ~TransitionMatrixCuBLAS() {
			cudaFree(this->data);
			cublasDestroy_v2(handle);
		}

		/**
		 * Overwrite the current matrix with unity matrix.
		 */
		virtual void setEye() {

			TransitionMatrixCBLAS<T> tmp(this->N);
			tmp.setEye();

			cudaMemcpy2D(this->data, this->ld * sizeof(T),
			             tmp.getData(), tmp.getLeadDimension() * sizeof(T),
			             this->N * sizeof(T), this->N,
			             cudaMemcpyHostToDevice);
		}

		/**
		 * Overwrite the current matrix with zeroes.
		 */
		virtual void setZero() {
			cudaMemset2D(this->data, this->ld * sizeof(T), 0, this->N * sizeof(T),
			             this->N);
		}

		/**
		 * Return a string that represents the matrix.
		 */
		virtual std::string to_string() const {

			TransitionMatrixCBLAS<T> tmp(this->N);

			cudaMemcpy2D(tmp.getData(),
			             tmp.getLeadDimension() * sizeof(T), this->data,
			             this->ld * sizeof(T), this->N * sizeof(T), this->N,
			             cudaMemcpyDeviceToHost);

			return tmp.to_string();
		}

		/**
		 * Multiply A with B and write the result to this.
		 * @param A A pointer to matrix A. Will not be changed.
		 * @param B A pointer to matrix B. Will not be changed.
		 */
		virtual void mult(const TransitionMatrix<T> *A,
		                  const TransitionMatrix<T> *B);

		/**
		 * Compute the variation distance of each state to the distribution pi
		 * @param pi A pointer to a probability distribution.
		 * @param dist Out parameter.
		 */
		virtual void variationDistance(const T *pi, T *dist) const {

			// call of external function
			cuda::cudaVariationDistance(this->data, this->N, this->ld, pi, dist);
		}

		/**
		 * Compute the total variation distance to the distribution.
		 * @param pi A probability distribution.
		 */
		virtual T totalVariationDistance(const T *pi) const {

			// call of external function
			return cuda::cudaTotalVariationDistance(this->data, this->N, this->ld, pi);
		}
	};


	// Template Specialization
	template<>
	void TransitionMatrixCuBLAS<float>::mult(
			const TransitionMatrix<float> *A,
			const TransitionMatrix<float> *B) {

		const float alpha = 1.0;
		const float beta = 0.0;

		// use cublas
		cublasSgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, this->N, this->N,
		               this->N, &alpha, B->getData(), A->getLeadDimension(), A->getData(),
		               A->getLeadDimension(), &beta, this->getData(),
		               this->getLeadDimension());
	}

	template<>
	void TransitionMatrixCuBLAS<double>::mult(
			const TransitionMatrix<double> *A,
			const TransitionMatrix<double> *B) {

		const double alpha = 1.0;
		const double beta = 0.0;

		// use cublas
		cublasDgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_N, this->N, this->N,
		               this->N, &alpha, B->getData(), A->getLeadDimension(), A->getData(),
		               A->getLeadDimension(), &beta, this->getData(),
		               this->getLeadDimension());

	}

}

#endif /* CUDA */

#endif /* INCLUDE_MARATHON_TRANSITIONMATRIXCUBLAS_H_ */
