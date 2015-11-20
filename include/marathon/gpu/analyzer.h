/*
 * Toolkit.h
 *
 *  Created on: Feb 12, 2015
 *      Author: rechner
 */

#ifndef DEVICE_ANALYZER_H_
#define DEVICE_ANALYZER_H_

// project includes
#include "../common/rational.h"
#include "../common/state_graph.h"
#include "transition_matrix.h"

namespace marathon {

namespace gpu {

/*********************************************************************
 * Forward declaration of Transition Matrix
 ********************************************************************/
template<typename T>
class DenseTransitionMatrix;

/*********************************************************************
 * Wrapper Functions for CUDA functionality
 ********************************************************************/
namespace cuda {

extern "C" void allocMemory(void** ptr, size_t size);
extern "C" void allocMemory2D(void** ptr, size_t* pitch, size_t width,
		size_t height);
extern "C" void freeMemory(void* ptr);
extern "C" void setMemory2D(void* ptr, int c, size_t pitch, size_t height);
extern "C" void copyHostToDevice(void* dst, void* src, size_t size);
extern "C" void copy2DHostToDevice(void* dst, size_t dpitch, void* src,
		size_t spitch, size_t width, size_t height);
extern "C" void copy2DDeviceToDevice(void* dst, size_t dpitch, void* src,
		size_t spitch, size_t width, size_t height);
extern "C" void copy2DDeviceToHost(void* dst, size_t dpitch, void* src,
		size_t spitch, size_t width, size_t height);
extern "C" void variationDistanceFloat(const float* ptr, const float * pi,
		float* dist, const size_t n, const size_t ld);
extern "C" void variationDistanceDouble(const double* ptr, const double * pi,
		double* dist, const size_t n, const size_t ld);
extern "C" float totalVariationDistanceFloat(const float* ptr, const float * pi,
		float* dist, const size_t n, const size_t ld);
extern "C" double totalVariationDistanceDouble(const double* ptr,
		const double * pi, double* dist, const size_t n, const size_t ld);
extern "C" float minVariationDistanceFloat(const float* ptr, const float * pi,
		float* dist, const size_t n, const size_t ld);
extern "C" double minVariationDistanceDouble(const double* ptr,
		const double * pi, double* dist, const size_t n, const size_t ld);
extern "C" bool initCublas();
extern "C" void finalizeCublas();
extern "C" void multFloat(const float* A, const size_t ldA, const float* B,
		const size_t ldB, float* C, const size_t ldC, const size_t n);
extern "C" void multDouble(const double* A, const size_t ldA, const double* B,
		const size_t ldB, double* C, const size_t ldC, const size_t n);

}

/***********************************************************************
 * Device API for Analysis of Mixing Time Properties of Markov chains
 *
 * The main task of this program is to compute the mixing time
 * and various upper and lower bounds of specific Markov chains.
 *
 ***********************************************************************/

/**
 * Initializes the cublas library for matrix multiplication.
 * Returns true, if cublas could be initialized successfully and a CUDA capable
 * device could be detected, else returns false.
 */
bool init();

/**
 * Finalizes the cublas library.
 */
void finalize();

template<typename T>
void variationDistance(const DenseTransitionMatrix<T>& P, const T* pi, T* tmp);

template<typename T>
T minVariationDistance(const DenseTransitionMatrix<T>& P, const T* pi, T* tmp);

/*************************************************************************
 * Computes the total variation distance d(P, pi).
 ************************************************************************/
template<typename T>
T totalVariationDistance(const DenseTransitionMatrix<T>& P, const T* pi,
		T* tmp);

/**************************************************************************
 * Total mixing time
 *
 * The total mixing time t is given by the smallest integer t, such that
 * dist(P^t, pi) < epsilon, where P is the transition matrix of mc and pi is
 * its stationary distribution.
 *
 * This method searches the total mixing time of mc by performing finger search.
 * Starting with matrix M=P, it squares M until dist(M,pi) < eps. By doing so,
 * the methods computes the power of two r, such that dist(P^r,pi) <= eps and
 * l = r/2 such that dist(P^l, pi) > eps. After finding the limits l and r, it
 * uses binary search for finding the smallest t such that dist(P^t,pi) < eps.
 *
 * In each step of binary search, the boundaries l and r are transformed.
 * To compute P^m with m=(l+r)/2 we first compute P^(m-l) and multiply the result
 * to P^l. By this, we need three temporary matrices.
 *
 *****************************************************************************/
template<typename T>
int totalMixingTime(const StateGraph* mc, const T epsilon);

// mixing time distribution
template<typename T>
void mixingTimeDistribution(const StateGraph* mc, std::vector<uint>& distr,
		const T epsilon);

} /* namespace Device */

} /* namespace Sampling */

#endif /* TOOLKIT_H_ */
