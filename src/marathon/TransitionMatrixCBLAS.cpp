/*
 * TransitionMatrixCBLAS.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

#include "../../include/marathon/TransitionMatrixCBLAS.h"
#include <cblas.h>


namespace marathon {
namespace tm {

template<>
void TransitionMatrixCBLAS<float>::mult(const TransitionMatrix<float>* A,
		const TransitionMatrix<float>* B) {

	const float alpha = 1.0;
	const float beta = 0.0;

	// use cblas
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->n, this->n,
			this->n, alpha, B->getData(), B->getLeadDimension(), A->getData(),
			A->getLeadDimension(), beta, this->data, this->ld);
}

template<>
void TransitionMatrixCBLAS<double>::mult(const TransitionMatrix<double>* A,
		const TransitionMatrix<double>* B) {

	const double alpha = 1.0;
	const double beta = 0.0;

	// use cblas
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->n, this->n,
			this->n, alpha, B->getData(), B->getLeadDimension(), A->getData(),
			A->getLeadDimension(), beta, this->data, this->ld);

}
}

}
