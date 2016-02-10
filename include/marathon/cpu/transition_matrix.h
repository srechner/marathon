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

template<typename T>
std::ostream& operator<<(std::ostream& out, const DenseTransitionMatrix<T>& m);

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

	DenseTransitionMatrix() :
			n(0), ld(0), data(nullptr) {

	}

	DenseTransitionMatrix(size_t n) {
		data = nullptr;
		init(n);
	}

	virtual ~DenseTransitionMatrix() {
		if (data != nullptr)
			delete[] data;
	}

	// swaps content of this and B
	void swapContent(DenseTransitionMatrix<T>& B) {
		T* tmp_data = data;
		size_t tmp_n = n;
		size_t tmp_ld = ld;

		data = B.data;
		n = B.n;
		ld = B.ld;

		B.data = tmp_data;
		B.n = tmp_n;
		B.ld = tmp_ld;
	}

	void init(size_t n) {
		if (data != nullptr)
			delete[] data;
		this->n = n;
		this->ld = n;
		data = new T[n * ld];
		if (data == nullptr)
			throw BAD_HOST_MALLOC_EXCEPTION;
	}

	void initFromStateGraph(const StateGraph* mc) {
		size_t omega = mc->getNumStates();
		init(omega);

		// convert to dense transition matrix
		memset(data, 0, ld * n * sizeof(T));

		//#pragma omp parallel for if(omega>1000)
		for (const Transition& t : mc->getArcs()) {
			data[t.u * ld + t.v] = t.p.convert_to<T>();
		}
	}

	void setZero() {
		assert(data != nullptr);
		memset(data, 0, n * ld * sizeof(T));
	}

	void setEye() {
		setZero();
#pragma omp parallel for if(n > 1000)
		for (size_t i = 0; i < n; i++)
			data[i * ld + i] = 1.0;

	}

	void copy(const DenseTransitionMatrix<T>& m) {
		assert(n == m.n && ld == m.ld);
		assert(data != nullptr && m.data != nullptr);
		memcpy(data, m.data, n * ld * sizeof(T));
	}

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
			DenseTransitionMatrix<T>& tmp) {
		assert(n == A.n && n == tmp.n);
		assert(n == A.ld && n == tmp.ld);

		// init
		if (k == 0) {
			this->setEye();
		} else {
			this->copy(A);
		}

		// create binary representation of k
		int bin[32];
		memset(bin, 0, 32 * sizeof(int));
		int l = 31;
		while (k > 0) {
			bin[l] = k % 2;
			k >>= 1;
			l--;
		}
		l += 2;

#ifdef DEBUG
		std::cout << "bin: ";
		for (int i = 0; i < 32; i++) {
			std::cout << bin[i];
		}
		std::cout << " l=" << l << std::endl;
#endif

		// binary exponentation - Left to Right (see Knuth Seminumerical Alg. Vol. 2 page 461)
		while (l < 32) {

			// square
			tmp.mult(*this, *this);
			tmp.swapContent(*this);

			// multiply
			if (bin[l] == 1) {
				// this = this*A
				tmp.mult(*this, A);
				tmp.swapContent(*this);
			}

			l++;
		}
	}

	DenseTransitionMatrix& operator=(DenseTransitionMatrix<T> m) {
		swapContent(m);
		return *this;
	}

};

template<>
void DenseTransitionMatrix<float>::mult(const DenseTransitionMatrix<float>& A,
		const DenseTransitionMatrix<float>& B);
template<>
void DenseTransitionMatrix<double>::mult(const DenseTransitionMatrix<double>& A,
		const DenseTransitionMatrix<double>& B);

template<typename T>
std::ostream& operator<<(std::ostream& out, const DenseTransitionMatrix<T>& m) {

	out << "n=" << m.n << ", ld=" << m.ld << std::endl;

	for (size_t i = 0; i < m.n; i++) {
		for (size_t j = 0; j < m.n; j++) {
			out << " " << std::setprecision(8) << std::fixed
					<< m.data[i * m.n + j];
		}
		out << std::endl;
	}

	return out;
}

}/* namespace cpu */

}/* namespace marathon */

#endif /* HOST_TRANSITION_MATRIX_H_ */
