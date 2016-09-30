/*
 * TransitionMatrixCBLAS.cpp
 *
 * Created on: Mar 23, 2016
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

#include "../../include/marathon/TransitionMatrixCBLAS.h"

#include <cblas.h>

namespace marathon {

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
