#ifndef _HOST_TRANSITION_MATRIX_CPP
#define _HOST_TRANSITION_MATRIX_CPP

#include "../../../include/marathon/hybrid/transition_matrix.h"

namespace marathon {

namespace hybrid {

template<>
void DenseTransitionMatrix<float>::mult(const DenseTransitionMatrix<float>& A,
		const DenseTransitionMatrix<float>& B) {
	cuda::multFloatXt(A.data, A.ld, B.data, B.ld, data, ld, n);
}

template<>
void DenseTransitionMatrix<double>::mult(const DenseTransitionMatrix<double>& A,
		const DenseTransitionMatrix<double>& B) {
	cuda::multDoubleXt(A.data, A.ld, B.data, B.ld, data, ld, n);
}

}

}

#endif
