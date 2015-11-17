/*
 * Toolkit.h
 *
 *  Created on: Feb 12, 2015
 *      Author: rechner
 */

#ifndef HOST_ANALYZER_H_
#define HOST_ANALYZER_H_

// project includes
#include "../rational.h"
#include "../state_graph.h"
#include "transition_matrix.h"

// STL includes
#include <list>

/**
 * Host API for Analysis of Mixing Time Properties of Markov chains
 *
 * The main task of this program is to compute the mixing time
 * and various upper and lower bounds of specific Markov chains.
 *
 */

namespace marathon {

namespace cpu {

/*********************************************************************
 * Forward declaration of Transition Matrix
 ********************************************************************/
template<typename T>
class DenseTransitionMatrix;

template<typename T>
void variationDistance(const DenseTransitionMatrix<T>& P, const T* pi, T* tmp);

template<typename T>
T minVariationDistance(const DenseTransitionMatrix<T>& P, const T* pi, T* tmp);

/**
 * Computes the total variation distance d(P, pi).
 */
template<typename T>
T totalVariationDistance(const DenseTransitionMatrix<T>& P, const T* pi,
		T* tmp);

/**
 *  Computes upper congestion bound by canonical path method.
 *  Path construction scheme can be given by function pointer.
 */
Rational pathCongestion(const StateGraph* mc);

/**
 * Computes the eigenvalue with second largest magnitute of the
 * transition matrix of mc.
 */
template<typename T>
T secondLargestEigenvalue(const StateGraph* mc);

/**
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
 */
template<typename T>
int totalMixingTime(const StateGraph* mc, const T epsilon);


/**
 * Computes the diameter of the graph, i.e. the maximal length of a shortest path
 * between some nodes of the graph. Each arc has a length of 1.
 */
int diameter(const StateGraph* G);

/**
 * Computes a histogram of the length of shortest path in G. Each arc of the graph
 * contributes one to the length of a path. For each length l for l in 0..diameter(G),
 * the number of shortest paths of G is stored in the vector count[l].
 */
void pathLengthHistogram(std::vector<long>& count, const StateGraph* G);


} /* namespace cpu */

} /* namespace marathon */

#endif /* TOOLKIT_H_ */
