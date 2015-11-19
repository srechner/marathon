/*
 * Matrix.h
 *
 *  Created on: Feb 20, 2015
 *      Author: rechner
 */

#ifndef DEVICE_TRANSITION_MATRIX_H_
#define DEVICE_TRANSITION_MATRIX_H_

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cblas.h>

// project includes
#include "../common/exceptions.h"
#include "../common/state_graph.h"
#include "analyzer.h"
#include "../cpu/transition_matrix.h"

namespace marathon {

namespace gpu {

/******************
 *
 * This is a device-side clone of Host::TransitionMatrix.
 * However, implementation differs from Host Implementation.
 *
 *****************/
template<typename T>
class DenseTransitionMatrix {

	/**************************************************************************
	 * Forward Declaration of friend functions
	 *************************************************************************/
	template<typename S>
	friend S totalVariationDistance(const DenseTransitionMatrix<S>&, const S*,
			S*);
	template<typename S>
	friend S minVariationDistance(const DenseTransitionMatrix<S>&, const S*,
			S*);
	template<typename S>
	friend void variationDistance(const DenseTransitionMatrix<S>&, const S*,
			S*);
	template<typename S>
	friend std::ostream& operator<<(std::ostream& out,
			const DenseTransitionMatrix<S>& m);

protected:

	size_t n, ld;
	T* data;

public:

	DenseTransitionMatrix();

	DenseTransitionMatrix(size_t n);

	virtual ~DenseTransitionMatrix();

	// swaps content of this and B
	void swapContent(DenseTransitionMatrix<T>& B);

	void init(size_t n);
	void initFromStateGraph(const StateGraph* mc);

	void setZero();
	void setEye();

	void copy(const gpu::DenseTransitionMatrix<T>& m);
	void copy(const cpu::DenseTransitionMatrix<T>& m);

	// matrix operations

	/**
	 * This matrix becomes product of A and B
	 */
	void mult(const DenseTransitionMatrix<T>& A,
			const DenseTransitionMatrix<T>& B);

	/**
	 * This matrix becomes A^k.
	 *
	 * The matrix tmp is used as working matrix.
	 */
	void pow(const DenseTransitionMatrix<T>& A, int k,
			DenseTransitionMatrix<T>& tmp);

	DenseTransitionMatrix& operator=(DenseTransitionMatrix<T> m);

	template<typename S>
	friend std::ostream& operator<<(std::ostream& out,
			const DenseTransitionMatrix<S>& m);

};

}/* namespace Host */

}/* namespace Sampling */

#endif /* MATRIX_CUH_ */
