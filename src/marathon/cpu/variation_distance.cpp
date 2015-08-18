#ifndef VARIATION_DISTANCE_CU_
#define VARIATION_DISTANCE_CU_

#include <iostream>

#include "../../../include/marathon/cpu/analyzer.h"

/**
 * Host Implementation of Variation Distance helper function
 */
template<typename T>
void marathon::cpu::variationDistance(const DenseTransitionMatrix<T>& P,
		const T* pi, T* tmp) {

#pragma omp parallel for if(P.n > 1000)
	for (size_t i = 0; i < P.n; i++) {
		T sum = 0;
		size_t j;
		for (j = 0; j < P.n; j++)
			sum += fabs(P.data[i * P.ld + j] - pi[j]);
		tmp[i] = sum / 2.0;
	}
}

template<typename T>
T marathon::cpu::totalVariationDistance(const DenseTransitionMatrix<T>& P,
		const T* pi, T* tmp) {

	marathon::cpu::variationDistance < T > (P, pi, tmp);
	return *std::max_element(tmp, tmp + P.n);
}

template<typename T>
T marathon::cpu::minVariationDistance(const DenseTransitionMatrix<T>&P,
		const T* pi, T* tmp) {

	marathon::cpu::variationDistance<T>(P, pi, tmp);
	return *std::min_element(tmp, tmp + P.n);
}

/**
 * Explicit template specialization
 */

template float marathon::cpu::totalVariationDistance<float>(
		const DenseTransitionMatrix<float>& P, const float* pi, float* tmp);
template double marathon::cpu::totalVariationDistance<double>(
		const DenseTransitionMatrix<double>& P, const double* pi, double* tmp);

#endif
