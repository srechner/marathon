#ifndef VARIATION_DISTANCE_CU_
#define VARIATION_DISTANCE_CU_

#include "../../../include/marathon/gpu/analyzer.h"

#include <iostream>

namespace marathon {

namespace gpu {

void init() {
	cuda::initCublas();
}

void finalize() {
	cuda::finalizeCublas();
}

template<>
void variationDistance<float>(const DenseTransitionMatrix<float>& P,
		const float* pi, float* tmp) {
	cuda::variationDistanceFloat(P.data, pi, tmp, P.n, P.ld);
}

template<>
void variationDistance<double>(const DenseTransitionMatrix<double>& P,
		const double* pi, double* tmp) {
	cuda::variationDistanceDouble(P.data, pi, tmp, P.n, P.ld);
}

template<>
float totalVariationDistance<float>(const DenseTransitionMatrix<float>& P,
		const float* pi, float* tmp) {
	return cuda::totalVariationDistanceFloat(P.data, pi, tmp, P.n, P.ld);
}

template<>
double totalVariationDistance<double>(const DenseTransitionMatrix<double>& P,
		const double* pi, double* tmp) {
	return cuda::totalVariationDistanceDouble(P.data, pi, tmp, P.n, P.ld);
}

template<>
float minVariationDistance<float>(const DenseTransitionMatrix<float>& P,
		const float* pi, float* tmp) {
	return cuda::minVariationDistanceFloat(P.data, pi, tmp, P.n, P.ld);
}

template<>
double minVariationDistance<double>(const DenseTransitionMatrix<double>& P,
		const double* pi, double* tmp) {
	return cuda::minVariationDistanceDouble(P.data, pi, tmp, P.n, P.ld);
}

/**
 * Explicit template specialization
 */

template float totalVariationDistance<float>(
		const DenseTransitionMatrix<float>& P, const float* pi, float* tmp);
template double totalVariationDistance<double>(
		const DenseTransitionMatrix<double>& P, const double* pi, double* tmp);

}

}

#endif
