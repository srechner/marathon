/*
 * TotalMixingTime.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

#include "../../include/marathon/marathon.h"
#include "../../include/marathon/TransitionMatrixCBLAS.h"
#ifdef CUDA
#include "../../include/marathon/TransitionMatrixCuBLAS.h"
#include "../../include/marathon/TransitionMatrixCuBLASXt.h"
#endif

template<typename T>
int marathon::totalMixingTime(const StateGraph* sg, const T eps,
		device_t device) {

	/* Variables */
	const size_t omega = sg->getNumStates();	// number of states
	TransitionMatrix<T> *P;						// transition matrix of sg
	TransitionMatrix<T> *tmp[3];				// working matrices
	T* pi;										// stationary distribution

	// check trivial cases
	if (omega == 0)
		return -1;
	else if (omega == 1)
		return 0;

	// convert stationary distribution into floating point number array
	pi = new T[omega];
	const rational Z = sg->getZ();
	for (int i = 0; i < omega; i++)
		pi[i] = (sg->getWeight(i) / Z).convert_to<T>();

	// decide with mode to use
	switch (device) {

	case device_t::CPU_ONLY:
		// allocate memory
		P = new TransitionMatrixCBLAS<T>(sg);
		tmp[0] = new TransitionMatrixCBLAS<T>(omega);
		tmp[1] = new TransitionMatrixCBLAS<T>(omega);
		tmp[2] = new TransitionMatrixCBLAS<T>(omega);
		break;

	case device_t::GPU_ONLY:
#ifdef CUDA
		// allocate memory
		P = new TransitionMatrixCuBLAS<T>(sg);
		tmp[0] = new TransitionMatrixCuBLAS<T>(omega);
		tmp[1] = new TransitionMatrixCuBLAS<T>(omega);
		tmp[2] = new TransitionMatrixCuBLAS<T>(omega);
		// prepare stationary distribution in device memory
		T* pi_d;
		myCudaMalloc((void**) &pi_d, omega * sizeof(T));
		myCudaMemcpyHostToDevice(pi_d, pi, omega * sizeof(T));
		std::swap(pi_d, pi);
#else
		std::cerr
				<< "marathon::TotalMixingTime: Error: Library was compiled without CUDA support. Will use CBLAS implementation."
				<< std::endl;

		// allocate memory
		P = new TransitionMatrixCBLAS<T>(sg);
		tmp[0] = new TransitionMatrixCBLAS<T>(omega);
		tmp[1] = new TransitionMatrixCBLAS<T>(omega);
		tmp[2] = new TransitionMatrixCBLAS<T>(omega);
#endif
		break;

	case device_t::HYBRID:

#ifdef CUDA
		// allocate memory
		P = new TransitionMatrixCuBLASXt<T>(sg);
		tmp[0] = new TransitionMatrixCuBLASXt<T>(omega);
		tmp[1] = new TransitionMatrixCuBLASXt<T>(omega);
		tmp[2] = new TransitionMatrixCuBLASXt<T>(omega);
#else
		std::cerr
				<< "marathon::TotalMixingTime: Error: Library was compiled without CUDA support. Will use CBLAS implementation."
				<< std::endl;

		// allocate memory
		P = new TransitionMatrixCBLAS<T>(sg);
		tmp[0] = new TransitionMatrixCBLAS<T>(omega);
		tmp[1] = new TransitionMatrixCBLAS<T>(omega);
		tmp[2] = new TransitionMatrixCBLAS<T>(omega);
#endif
		break;

	default:
		std::cerr << "marathon::TotalMixingTime: Error! unknown option: "
				<< device << std::endl;
		return 1;
		break;
	}

	// react to bad memory allocation
	if (pi == nullptr || tmp[0] == nullptr || tmp[1] == nullptr
			|| tmp[2] == nullptr)
		return -1;

	// Search for mixing time
	uint l = 0;
	uint r = 1;

	// tmp[0] = P
	tmp[0]->copy(P);

	// First Phase: Square tmp[0] until dist(tmp[0], pi) < eps
	T d = tmp[0]->totalVariationDistance(pi);

	while (d >= eps) {
		// tmp[1] = tmp[0]*tmp[0]
		tmp[1]->mult(tmp[0], tmp[0]);
		tmp[0]->swap(tmp[1]);
		d = tmp[0]->totalVariationDistance(pi);
		//std::cout << "l=" << l << " r=" << r << " d=" << d << std::endl;
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
		tmp[2]->pow(P, m - l);

		// tmp[0] = P^l * P^(m-l) = P^m
		tmp[0]->mult(tmp[1], tmp[2]);
		T dist_m = tmp[0]->totalVariationDistance(pi);

		if (dist_m >= eps) {
			l = m;
			tmp[0]->swap(tmp[1]);
		} else {
			r = m;
		}
	}

	// free memory
	delete P;
	delete tmp[0];
	delete tmp[1];
	delete tmp[2];

#ifdef CUDA
	if (device == GPU_ONLY)
	myCudaFree(pi);
	else
#endif
	delete[] pi;

	return r;
}

/**
 * Export Template Specialization
 */
template int marathon::totalMixingTime<float>(const StateGraph*, const float,
		device_t);
template int marathon::totalMixingTime<double>(const StateGraph*, const double,
		device_t);
