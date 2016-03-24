/*
 * TransitionMatrixCPU.h
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#ifndef INCLUDE_MARATHON_TRANSITIONMATRIXCUBLAS_H_
#define INCLUDE_MARATHON_TRANSITIONMATRIXCUBLAS_H_

#include "marathon.h"
#include "TransitionMatrix.h"
#include "TransitionMatrixCBLAS.h"
#include "CudaWrapper.h"

namespace marathon {

/**
 * Declare external functions.
 */
template<typename T>
extern void cudaVariationDistance(const T* data, const size_t n,
		const size_t ld, const T * pi, T* dist);

template<typename T>
extern T cudaTotalVariationDistance(const T* data, const size_t n,
		const size_t ld, const T * pi);

template<typename T>
class TransitionMatrixCuBLAS: public TransitionMatrix<T> {

protected:

	/**
	 * Copy the content of matrix P to this.
	 */
	virtual void copy(const TransitionMatrix<T>* P) {

		const TransitionMatrixCuBLAS<T>* X =
				(const TransitionMatrixCuBLAS<T>*) P;
		this->n = X->n;
		this->ld = X->ld;
		size_t pitch = this->ld * sizeof(T);
		myCudaMemcpy2DDeviceToDevice(this->data, this->ld * sizeof(T), X->data,
				X->ld * sizeof(T), this->n * sizeof(T), this->n);
	}

	/**
	 * Return a pointer to a Transition matrix of subtype instance.
	 */
	virtual TransitionMatrix<T>* generateSubTypeInstance(const int n) {
		TransitionMatrix<T>* res = new TransitionMatrixCuBLAS<T>(n);
		return res;
	}

public:

	TransitionMatrixCuBLAS(const int n) {
		this->n = n;

		// allocate aligned 2d memory
		size_t pitch;
		myCudaMallocPitch((void**) &this->data, &pitch, n * sizeof(T), n);
		this->ld = pitch / sizeof(T);
	}

	TransitionMatrixCuBLAS(const StateGraph* sg) :
			TransitionMatrixCuBLAS<T>(sg->getNumStates()) {

		TransitionMatrixCBLAS<T> tmp(sg);

		myCudaMemcpy2DHostToDevice(this->data, this->ld * sizeof(T),
				tmp.getData(), tmp.getLeadDimension() * sizeof(T),
				this->n * sizeof(T), this->n);
	}

	virtual ~TransitionMatrixCuBLAS() {
		myCudaFree(this->data);
	}

	/**
	 * Overwrite the current matrix with unity matrix.
	 */
	virtual void setEye() {

		TransitionMatrixCBLAS<T> tmp(this->n);
		tmp.setEye();

		myCudaMemcpy2DHostToDevice(this->data, this->ld * sizeof(T),
				tmp.getData(), tmp.getLeadDimension() * sizeof(T),
				this->n * sizeof(T), this->n);
	}

	/**
	 * Overwrite the current matrix with zeroes.
	 */
	virtual void setZero() {

		// call of external function
		myCudaMemset2D(this->data, this->ld * sizeof(T), 0, this->n * sizeof(T),
				this->n);
	}

	/**
	 * Return a string that represents the matrix.
	 */
	virtual std::string to_string() const {

		TransitionMatrixCBLAS<T> tmp(this->n);

		myCudaMemcpy2DDeviceToHost(tmp.getData(),
				tmp.getLeadDimension() * sizeof(T), this->data,
				this->ld * sizeof(T), this->n * sizeof(T), this->n);

		return tmp.to_string();
	}

	/**
	 * Multiply A with B and write the result to this.
	 * @param A A pointer to matrix A. Will not be changed.
	 * @param B A pointer to matrix B. Will not be changed.
	 */
	virtual void mult(const TransitionMatrix<T>* A,
			const TransitionMatrix<T>* B);

	/**
	 * Compute the variation distance of each state to the distribution pi
	 * @param pi A pointer to a probability distribution.
	 * @param dist Out parameter.
	 */
	virtual void variationDistance(const T* pi, T* dist) const {

		// call of external function
		cudaVariationDistance(this->data, this->n, this->ld, pi, dist);
	}

	/**
	 * Compute the total variation distance to the distribution.
	 * @param pi A probability distribution.
	 */
	virtual T totalVariationDistance(const T* pi) const {

		// call of external function
		return cudaTotalVariationDistance(this->data, this->n, this->ld, pi);
	}

};

}

#endif /* INCLUDE_MARATHON_TRANSITIONMATRIXCUBLAS_H_ */
