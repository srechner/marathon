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
#include "Random.h"
#include "Combinatorics.h"

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
 * This implementation searches the total mixing time of mc by a two step procedure.
 * Starting with matrix M=P, it squares M until dist(M,pi) < eps. By doing so,
 * the methods computes the power of two r, such that dist(P^r,pi) <= eps and
 * l = r/2 such that dist(P^l, pi) > eps. After finding the limits l and r, it
 * uses binary search for finding the smallest t such that dist(P^t,pi) < eps.
 *
 * In each step of binary search, the boundaries l and r are transformed.
 * To compute P^m with m=(l+r)/2 we first compute P^(m-l) and multiply the result
 * to P^l. By this, we need three temporary matrices.
 *
 * @param sg A pointer to a state graph object.
 * @param eps The distance to stationary distribution.
 * @param dev Determines which matrix multiplication library to use.
 *    Can be one of the following: CPU_ONLY, GPU_ONLY, HYBRID.
 */
template<typename T>
int totalMixingTime(const StateGraph* sg, const T eps, device_t dev =
		CPU_ONLY);

/**
 * Computes upper congestion bound by canonical path method.
 * Path construction scheme can be given by function pointer.
 * @param sg A pointer to a state graph object.
 * @param constructPath An object of a path construction scheme class.
 * @param eps The distance to stationary distribution.
 * @return Maximum congestion of a path.
 */
template<typename T>
T upperPathCongestionBound(const StateGraph* sg,
		const PathConstructionScheme& pcs, T eps);

/**
 * Compute the lower spectral bound of the state graph.
 * @param sg A pointer to a state graph object.
 * @param eps The distance to stationary distribution.
 * @return A lower bound of the total mixing time.
 */
template<typename T>
T lowerSpectralBound(const StateGraph* sg, T eps);


/**
 * Compute the upper spectral bound of the state graph.
 * @param sg A pointer to a state graph object.
 * @param eps The distance to stationary distribution.
 * @return A lower bound of the total mixing time.
 */
template<typename T>
T upperSpectralBound(const StateGraph* sg, T eps);


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
#include "chain/bipgraph/Curveball.h"
#include "chain/bipgraph/ExtendedCurveball.h"
#include "chain/bipgraph/CurveballForbiddenEntries.h"

// include Path Construction Schemes
#include "chain/matching/JS89CanPath.h"
#include "chain/bipgraph/KannanCanPath.h"

#endif /* MARATHON_H_ */
