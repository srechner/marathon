/*
 * marathon.cpp
 *
 *  Created on: Mar 24, 2016
 *      Author: rechner
 */

#include "../../include/marathon/marathon.h"
#include "../../include/marathon/Eigenvalues.h"
#include "../../include/marathon/PathCongestion.h"

#include <cmath>

/**
 * Implement Eigenvalue functions.
 */

template<typename T>
T marathon::lowerSpectralBound(const StateGraph* sg, T eps) {

	const double lambda = fabs(
			eigenvalue::eigenvalue<T>(sg,
					eigenvalue::eigenvalue_t::_2ndLargestMagnitude));

	if (lambda < 1e-10)
		return std::numeric_limits<T>::infinity();

	return 0.5 * lambda / (1.0 - lambda) * -log(2.0 * eps);
}

template<typename T>
T marathon::upperSpectralBound(const StateGraph* sg, T eps) {

	const double lambda = fabs(
			eigenvalue::eigenvalue<T>(sg,
					eigenvalue::eigenvalue_t::_2ndLargestMagnitude));

	if (fabs(lambda) < 1e-10)
		return std::numeric_limits<T>::infinity();

	const rational pimin = sg->getMinWeight() / sg->getZ();

	return -log(eps * pimin.convert_to<double>()) / (1.0 - lambda);
}

/**
 * Implement Path Congestion functions.
 */

template<typename T>
T marathon::upperPathCongestionBound(const StateGraph* sg,
		const PathConstructionScheme& pcs, T eps) {

	const rational load = marathon::pathCongestion::pathCongestion(sg, pcs);
	const rational pimin = sg->getMinWeight() / sg->getZ();

	return load.convert_to<T>() * -log(pimin.convert_to<T>() * eps);
}


/**
 * Export Template Specializations
 */

template float marathon::lowerSpectralBound(const StateGraph*, float);
template double marathon::lowerSpectralBound(const StateGraph*, double);
template float marathon::upperSpectralBound(const StateGraph*, float);
template double marathon::upperSpectralBound(const StateGraph*, double);
template float marathon::upperPathCongestionBound(const StateGraph* sg, const PathConstructionScheme&, float);
template double marathon::upperPathCongestionBound(const StateGraph* sg, const PathConstructionScheme&, double);
