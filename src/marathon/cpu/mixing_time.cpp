#include "../../../include/marathon/cpu/analyzer.h"

template<typename T>
int marathon::cpu::totalMixingTime(const StateGraph* mc, const T epsilon) {

	/* Variables */
	size_t omega;							// number of states
	DenseTransitionMatrix<T> P;				// transition matrix
	DenseTransitionMatrix<T> tmp[3];		// working memory
	T* pi;									// stationary distribution
	T* dist;								// working array for

	omega = mc->getNumStates();

	// check trivial cases
	if (omega == 0)
		return -1;
	else if (omega == 1)
		return 0;

	pi = new T[omega];
	dist = new T[omega];
	if (pi == nullptr)
		return -1;
	if (dist == nullptr) {
		delete[] pi;
		return -1;
	}

	// try to allocate memory
	try {
		tmp[0].init(omega);
		tmp[1].init(omega);
		tmp[2].init(omega);
		P.initFromStateGraph(mc);

	} catch (int n) {
		delete[] pi;
		delete[] dist;
		std::cerr << "bad malloc" << std::endl;
		return -1;
	}

	// convert stationary distribution into floating point number array
#pragma omp parallel for if(omega > 1000)
	for (size_t i = 0; i < omega; i++) {
		pi[i] = mc->getStationary(i).convert_to<T>();
	}

	// Search for mixing time
	uint l = 0;
	uint r = 1;

	// First Phase: Square tmp[0] until dist(tmp[0], pi) < eps
	tmp[0].copy(P);
	T d = totalVariationDistance<T>(tmp[0], pi, dist);

	while (d >= epsilon) {
		tmp[1].mult(tmp[0], tmp[0]);
		tmp[0].swapContent(tmp[1]);
		d = totalVariationDistance<T>(tmp[0], pi, dist);
		l = r;
		r *= 2;
	}

	/*
	 * State of the variables:
	 *
	 * tmp[0] = P^r
	 * tmp[1] = P^l
	 *
	 * dist_l = dist(tmp[1], pi) <= eps < dist(tmp[0], pi) = dist_r
	 */

	// Second Phase: Binary Search
	// Invariant: tmp[1] = P^l
	while (l < r - 1) {
		uint m = (l + r) / 2;

		// tmp[2] =  P^(m-l)
		tmp[2].pow(P, m - l, tmp[0]);

		// tmp[0] = P^l * P^(m-l) = P^m
		tmp[0].mult(tmp[1], tmp[2]);
		T dist_m = totalVariationDistance<T>(tmp[0], pi, dist);

		if (dist_m >= epsilon) {
			l = m;
			tmp[0].swapContent(tmp[1]);
		} else {
			r = m;
		}
	}

	// free memory
	delete[] pi;
	delete[] dist;
	return r;
}

/**
 * Export Template Specialization
 */
template int marathon::cpu::totalMixingTime<float>(const StateGraph*,
		const float);
template int marathon::cpu::totalMixingTime<double>(const StateGraph*,
		const double);
