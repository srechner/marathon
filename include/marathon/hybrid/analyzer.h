/*
 * Toolkit.h
 *
 *  Created on: Feb 12, 2015
 *      Author: rechner
 */

#ifndef HYBRID_ANALYZER_H_
#define HYBRID_ANALYZER_H_

// project includes
#include "../rational.h"
#include "../state_graph.h"
#include "transition_matrix.h"

/**
 * Host/Device API for Analysis of Mixing Time Properties of Markov chains
 *
 * The main task of this program is to compute the mixing time
 * and various upper and lower bounds of specific Markov chains.
 *
 */

namespace marathon {

namespace hybrid {

namespace cuda {

extern "C" void initCublasXt();
extern "C" void finalizeCublasXt();

}

void init();

void finalize();

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

} /* namespace Hybrid */

} /* namespace Sampling */

#endif /* TOOLKIT_H_ */
