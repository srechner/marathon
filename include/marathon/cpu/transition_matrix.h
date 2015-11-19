/*
 * Matrix.h
 *
 *  Created on: Feb 20, 2015
 *      Author: rechner
 */

#ifndef HOST_TRANSITION_MATRIX_H_
#define HOST_TRANSITION_MATRIX_H_

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cblas.h>

// project includes
#include "../common/exceptions.h"
#include "../common/state_graph.h"
#include "analyzer.h"

namespace marathon {

// Forward Declaration of device side Transition Matrix
namespace gpu {

template<typename T>
class DenseTransitionMatrix;

template<typename S>
std::ostream& operator<<(std::ostream& out, const DenseTransitionMatrix<S>& m);

}

namespace cpu {

template<typename T>
class DenseTransitionMatrix {

	/**************************************************************************
	 * Forward Declaration of friends
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

	template<typename S>
	friend class gpu::DenseTransitionMatrix;

	template<typename S>
	friend std::ostream& cpu::operator<<(std::ostream& out,
			const cpu::DenseTransitionMatrix<S>& m);

	template<typename S>
	friend std::ostream& ::marathon::gpu::operator<<(std::ostream& out,
			const ::marathon::gpu::DenseTransitionMatrix<S>& m);

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

	void copy(const DenseTransitionMatrix<T>& m);

	// matrix operations

	/**
	 * This matrix becomes product of A and B
	 */
	virtual void mult(const DenseTransitionMatrix<T>& A,
			const DenseTransitionMatrix<T>& B);

	/**
	 * This matrix becomes A^k.
	 *
	 * The matrix tmp is used as working matrix.
	 */
	void pow(const DenseTransitionMatrix<T>& A, int k,
			DenseTransitionMatrix<T>& tmp);

	DenseTransitionMatrix& operator=(DenseTransitionMatrix<T> m);

};

}/* namespace cpu */

}/* namespace marathon */

#endif /* MATRIX_CUH_ */
