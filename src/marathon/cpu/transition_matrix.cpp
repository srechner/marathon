#include "../../../include/marathon/cpu/transition_matrix.h"

/**
 * Explicit Template Specification
 */

namespace marathon {
namespace cpu {

/*
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
		const marathon::StateGraph*);
template void DenseTransitionMatrix<double>::pow(
		const DenseTransitionMatrix<double>& A, int k,
		DenseTransitionMatrix<double>& tmp);
template DenseTransitionMatrix<double>& DenseTransitionMatrix<double>::operator=(
		DenseTransitionMatrix<double> m);
template std::ostream& operator<<(std::ostream& out,
		const DenseTransitionMatrix<double>& m);
*/

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

template class DenseTransitionMatrix<float>;
template class DenseTransitionMatrix<double>;

}
}
