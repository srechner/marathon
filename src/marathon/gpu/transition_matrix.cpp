#ifndef _DEVICE_TRANSITION_MATRIX_CU
#define _DEVICE_TRANSITION_MATRIX_CU

#include "../../../include/marathon/gpu/transition_matrix.h"

namespace marathon {

namespace gpu {

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
		cuda::freeMemory(data);
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
	this->n = n;
	if (data != nullptr)
		cuda::freeMemory(data);
	cuda::allocMemory2D((void**) &data, &ld, n * sizeof(T), n);
	ld /= sizeof(T);
}

template<typename T>
void DenseTransitionMatrix<T>::setZero() {
	cuda::setMemory2D(data, 0, n, ld);
}

template<typename T>
void DenseTransitionMatrix<T>::setEye() {
	::marathon::cpu::DenseTransitionMatrix<T> hostMatrix(n);
	hostMatrix.setEye();
	cuda::copy2DHostToDevice(data, ld * sizeof(T), hostMatrix.data,
			hostMatrix.ld * sizeof(T), n * sizeof(T), n);
}

template<typename T>
void DenseTransitionMatrix<T>::copy(const DenseTransitionMatrix<T>& m) {
	assert(n == m.n);
	assert(data != nullptr && m.data != nullptr);
	cuda::copy2DHostToDevice(data, ld * sizeof(T), m.data, m.ld * sizeof(T),
			n * sizeof(T), n);
}

template<typename T>
void DenseTransitionMatrix<T>::copy(const ::marathon::cpu::DenseTransitionMatrix<T>& m) {
	assert(n == m.n);
	assert(data != nullptr && m.data != nullptr);
	cuda::copy2DHostToDevice(data, ld * sizeof(T), m.data, m.ld * sizeof(T),
			n * sizeof(T), n);
}

template<typename T>
void DenseTransitionMatrix<T>::initFromStateGraph(const StateGraph* mc) {

	size_t omega = mc->getNumStates();
	init(omega);
	::marathon::cpu::DenseTransitionMatrix<T> hostMatrix(omega);
	hostMatrix.initFromStateGraph(mc);
	copy(hostMatrix);
}

template<typename T>
void DenseTransitionMatrix<T>::pow(
		const DenseTransitionMatrix<T>& A, int k,
		DenseTransitionMatrix<T>& tmp) {
	assert(n == A.n && n == tmp.n);
	assert(ld == A.ld && ld == tmp.ld);

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

	cpu::DenseTransitionMatrix<T> hostMatrix(m.n);
	cuda::copy2DDeviceToHost(hostMatrix.data, hostMatrix.ld * sizeof(T), m.data,
			m.ld * sizeof(T), m.n * sizeof(T), m.n);
	out << hostMatrix;
	return out;
}

template<>
void DenseTransitionMatrix<float>::mult(
		const DenseTransitionMatrix<float>& A,
		const DenseTransitionMatrix<float>& B) {

	cuda::multFloat(A.data, A.ld, B.data, B.ld, data, ld, n);
}

template<>
void DenseTransitionMatrix<double>::mult(
		const DenseTransitionMatrix<double>& A,
		const DenseTransitionMatrix<double>& B) {

	cuda::multDouble(A.data, A.ld, B.data, B.ld, data, ld, n);
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
template void DenseTransitionMatrix<float>::initFromStateGraph(const StateGraph*);
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
