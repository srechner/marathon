/*
 * Matrix.h
 *
 *  Created on: Feb 20, 2015
 *      Author: rechner
 */

#ifndef HYBRID_TRANSITION_MATRIX_H_
#define HYBRID_TRANSITION_MATRIX_H_

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cblas.h>

// project includes
#include "../exceptions.h"
#include "../state_graph.h"
#include "analyzer.h"
#include "../cpu/transition_matrix.h"

namespace marathon {

namespace hybrid {

/*********************************************************************
 * Wrapper Functions for CUDA functionality
 ********************************************************************/
namespace cuda {

extern "C" void multFloatXt(const float* A, const size_t ldA, const float* B,
		const size_t ldB, float* C, const size_t ldC, const size_t n);
extern "C" void multDoubleXt(const double* A, const size_t ldA, const double* B,
		const size_t ldB, double* C, const size_t ldC, const size_t n);

}

template<typename T>
class DenseTransitionMatrix : public ::marathon::cpu::DenseTransitionMatrix<T> {

public:

	/**
	 * This matrix becomes product of A and B
	 */
	virtual void mult(const DenseTransitionMatrix<T>& A,
			const DenseTransitionMatrix<T>& B);

};

}

}

#endif /* HYBRID_TRANSITION_MATRIX_H_ */
