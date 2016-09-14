/*
 * TransitionMatrixCBLAS.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

#include "../../include/marathon/config.h"

#ifdef CUDA

#include "../../include/marathon/TransitionMatrixCuBLASXt.h"
#include <cublasXt.h>

namespace marathon {

    extern cublasXtHandle_t cublasXtHandle;

    namespace tm {

        template<>
        void TransitionMatrixCuBLASXt<float>::mult(const TransitionMatrix<float> *A,
                                                   const TransitionMatrix<float> *B) {

            const float alpha = 1.0;
            const float beta = 0.0;

            // use cublasXt
            cublasXtSgemm(cublasXtHandle, CUBLAS_OP_N, CUBLAS_OP_N, this->n, this->n,
                          this->n, &alpha, B->getData(), A->getLeadDimension(), A->getData(),
                          A->getLeadDimension(), &beta, this->getData(),
                          this->getLeadDimension());
        }

        template<>
        void TransitionMatrixCuBLASXt<double>::mult(const TransitionMatrix<double> *A,
                                                    const TransitionMatrix<double> *B) {

            const double alpha = 1.0;
            const double beta = 0.0;

            // use cblas
            cublasXtDgemm(cublasXtHandle, CUBLAS_OP_N, CUBLAS_OP_N, this->n, this->n,
                          this->n, &alpha, B->getData(), A->getLeadDimension(), A->getData(),
                          A->getLeadDimension(), &beta, this->getData(),
                          this->getLeadDimension());

        }
    }

}

#endif
