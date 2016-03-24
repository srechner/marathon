/*
 * TransitionMatrixCPU.h
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#ifndef INCLUDE_MARATHON_TRANSITIONMATRIXCUBLASXT_H_
#define INCLUDE_MARATHON_TRANSITIONMATRIXCUBLASXT_H_

#include "TransitionMatrixCBLAS.h"

namespace marathon {

template<typename T>
class TransitionMatrixCuBLASXt: public TransitionMatrixCBLAS<T> {

protected:

public:

	TransitionMatrixCuBLASXt(const int n) :
			TransitionMatrixCBLAS<T>(n) {

	}

	TransitionMatrixCuBLASXt(const StateGraph* sg) :
			TransitionMatrixCBLAS<T>(sg) {
	}

	virtual ~TransitionMatrixCuBLASXt() {

	}

	/**
	 * Multiply A with B and write the result to this.
	 * @param A A pointer to matrix A. Will not be changed.
	 * @param B A pointer to matrix B. Will not be changed.
	 */
	virtual void mult(const TransitionMatrix<T>* A,
			const TransitionMatrix<T>* B);

};

}

#endif /* INCLUDE_MARATHON_TRANSITIONMATRIXCUBLASXT_H_ */
