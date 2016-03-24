/*
 * marathon.h
 * 
 * main header for the marathon library.
 *
 *
 *  Created on: Jun 15, 2015
 *      Author: rechner
 */

#ifndef MARATHON_H_
#define MARATHON_H_

#include "PathConstructionScheme.h"
#include "MarkovChain.h"
#include "StateGraph.h"

namespace marathon {

/*!
 * @brief The marathon library has to be initialized to use it.
 */
void init();

/*!
 * @brief It is better to close the library at the end.
 */
void finalize();

/**
 * Options for computation of total mixing time.
 */

enum device_t {
	CPU_ONLY, GPU_ONLY, HYBRID
};

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
int totalMixingTime(const StateGraph* mc, const T epsilon, device_t device =
		CPU_ONLY);

/**
 *  Computes upper congestion bound by canonical path method.
 *  Path construction scheme can be given by function pointer.
 *  @param sg A pointer to a state graph object.
 *  @param constructPath An object of a path construction scheme class.
 */
rational pathCongestion(const StateGraph* sg,
		const PathConstructionScheme& pcs);

/**
 * Options for computation of eigenvalues.
 */
enum eigenvalue_t {
	// eigenvalue options
	_2ndLargestMagnitude,
	_2ndLargestAlgebraic,
};

/**
 * Computes the eigenvalue with second largest magnitute of the
 * transition matrix of mc.
 * @param mc State Graph Representation of Markov Chain.
 * @param which Which Eigenvalue to compute. Options are: 2nd largest in magnitute and 2nd largest algebraic.
 * @return The corresponding Eigenvalue.
 */
template<typename T>
T eigenvalue(const StateGraph* mc, eigenvalue_t which);

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

}

// include Transition Matrix classes
#include "TransitionMatrixCBLAS.h"
#include "TransitionMatrixCuBLAS.h"
#include "TransitionMatrixCuBLASXt.h"

// include Markov chains
#include "chain/matching/Broder86.h"
#include "chain/matching/JSV04.h"
#include "chain/bipgraph/SwitchChain.h"
#include "chain/bipgraph/SwitchChainBerger.h"

// include Path Construction Schemes
#include "chain/matching/JS89CanPath.h"
#include "chain/bipgraph/KannanCanPath.h"

#endif /* MARATHON_H_ */
