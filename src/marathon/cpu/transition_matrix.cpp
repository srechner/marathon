#ifndef _HOST_TRANSITION_MATRIX_CPP
#define _HOST_TRANSITION_MATRIX_CPP

#include "../../../include/marathon/cpu/transition_matrix.h"

namespace marathon {

namespace cpu {

template<typename T>
DenseTransitionMatrix<T>::DenseTransitionMatrix() :
		n(0), ld(0), data(nullptr) {

}

template<typename T>
DenseTransitionMatrix<T>::DenseTransitionMatrix(size_t n) {
	data = nullptr;
	init(n);
}

template<typename T>
DenseTransitionMatrix<T>::~DenseTransitionMatrix() {
	if (data != nullptr)
		delete[] data;
}

template<typename T>
void DenseTransitionMatrix<T>::swapContent(DenseTransitionMatrix<T>& B) {

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

template<typename T>
void DenseTransitionMatrix<T>::init(size_t n) {
	if (data != nullptr)
		delete[] data;
	this->n = n;
	this->ld = n;
	data = new T[n * ld];
	if (data == nullptr)
		throw BAD_HOST_MALLOC_EXCEPTION;
}

template<typename T>
void DenseTransitionMatrix<T>::setZero() {
	assert(data != nullptr);
	memset(data, 0, n * ld * sizeof(T));
}

template<typename T>
void DenseTransitionMatrix<T>::setEye() {
	setZero();
#pragma omp parallel for if(n > 1000)
	for (size_t i = 0; i < n; i++)
		data[i * ld + i] = 1.0;
}

template<typename T>
void DenseTransitionMatrix<T>::copy(const DenseTransitionMatrix<T>& m) {
	assert(n == m.n && ld == m.ld);
	assert(data != nullptr && m.data != nullptr);
	memcpy(data, m.data, n * ld * sizeof(T));
}

template<typename T>
void DenseTransitionMatrix<T>::initFromStateGraph(const StateGraph* mc) {

	size_t omega = mc->getNumStates();
	init(omega);

	// convert to dense transition matrix
	memset(data, 0, ld * n * sizeof(T));

#pragma omp parallel for if(omega>1000)
	for (size_t i = 0; i<mc->getNumArcs(); i++) {
		Transition uv = mc->arcs[i];
		data[uv.u * ld + uv.v] = uv.p.convert_to<T>();
	}
}

template<typename T>
void DenseTransitionMatrix<T>::pow(const DenseTransitionMatrix<T>& A, int k,
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

template<typename T>
DenseTransitionMatrix<T>& DenseTransitionMatrix<T>::operator=(
		DenseTransitionMatrix<T> m) {
	swapContent(m);
	return *this;
}

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

template<>
void DenseTransitionMatrix<float>::mult(const DenseTransitionMatrix<float>& A,
		const DenseTransitionMatrix<float>& B) {

	assert(A.n == B.n && A.n == n);
	assert(A.ld == B.ld && A.ld == ld);

	const float alpha_d = 1.0;
	const float beta_d = 0.0;

	// use cblas
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A.n, A.n, A.n,
			alpha_d, B.data, B.ld, A.data, A.ld, beta_d, data, ld);
}

template<>
void DenseTransitionMatrix<double>::mult(const DenseTransitionMatrix<double>& A,
		const DenseTransitionMatrix<double>& B) {

	assert(A.n == B.n && A.n == n);
	assert(A.ld == B.ld && A.ld == ld);

	const double alpha_d = 1.0;
	const double beta_d = 0.0;

	// use cblas
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha_d,
			B.data, B.ld, A.data, A.ld, beta_d, data, ld);
}

/**
 * Explicite Template Specification
 */

template DenseTransitionMatrix<float>::DenseTransitionMatrix();
template DenseTransitionMatrix<float>::DenseTransitionMatrix(size_t);
template DenseTransitionMatrix<float>::~DenseTransitionMatrix();
template void DenseTransitionMatrix<float>::swapContent(
		DenseTransitionMatrix<float>&);
template void DenseTransitionMatrix<float>::init(size_t);
template void DenseTransitionMatrix<float>::setZero();
template void DenseTransitionMatrix<float>::setEye();
template std::ostream& operator<<(std::ostream& out,
		const DenseTransitionMatrix<float>& m);
template DenseTransitionMatrix<float>& DenseTransitionMatrix<float>::operator=(
		DenseTransitionMatrix<float> m);
template void DenseTransitionMatrix<float>::pow(
		const DenseTransitionMatrix<float>& A, int k,
		DenseTransitionMatrix<float>& tmp);
template void DenseTransitionMatrix<float>::initFromStateGraph(
		const StateGraph*);
template void DenseTransitionMatrix<float>::copy(
		const DenseTransitionMatrix<float>&);

template DenseTransitionMatrix<double>::DenseTransitionMatrix();
template DenseTransitionMatrix<double>::DenseTransitionMatrix(size_t);
template DenseTransitionMatrix<double>::~DenseTransitionMatrix();
template void DenseTransitionMatrix<double>::swapContent(
		DenseTransitionMatrix<double>&);
template void DenseTransitionMatrix<double>::init(size_t);
template void DenseTransitionMatrix<double>::setZero();
template void DenseTransitionMatrix<double>::setEye();
template void DenseTransitionMatrix<double>::copy(
		const DenseTransitionMatrix<double>&);
template void DenseTransitionMatrix<double>::initFromStateGraph(
		const StateGraph*);
template void DenseTransitionMatrix<double>::pow(
		const DenseTransitionMatrix<double>& A, int k,
		DenseTransitionMatrix<double>& tmp);
template DenseTransitionMatrix<double>& DenseTransitionMatrix<double>::operator=(
		DenseTransitionMatrix<double> m);
template std::ostream& operator<<(std::ostream& out,
		const DenseTransitionMatrix<double>& m);

}

}

#endif
