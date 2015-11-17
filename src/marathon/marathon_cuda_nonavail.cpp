/*
 * marathon.cpp
 *
 *  Created on: Nov 17, 2015
 *      Author: steffen
 */

#include "../../include/marathon.h"
#include "../../include/marathon/cpu/analyzer.h"
#include "../../include/marathon/gpu/analyzer.h"
#include "../../include/marathon/hybrid/analyzer.h"

namespace marathon {

template<typename T>
int totalMixingTime(const StateGraph* mc, const T epsilon, device_t device) {
	switch (device) {
	case CPU_ONLY:
		return cpu::totalMixingTime(mc, epsilon);
	case GPU_ONLY:
	case HYBRID:
		std::cerr
				<< "Error! marathon build without CUDA support. Will use cpu implementation."
				<< std::endl;
		return cpu::totalMixingTime(mc, epsilon);
	}
}

Rational pathCongestion(const StateGraph* mc) {
	return cpu::pathCongestion(mc);
}

template<typename T>
T secondLargestEigenvalue(const StateGraph* mc) {
	return cpu::secondLargestEigenvalue<T>(mc);
}

int diameter(const StateGraph* G) {
	return cpu::diameter(G);
}

void pathLengthHistogram(std::vector<long>& count, const StateGraph* G) {
	return cpu::pathLengthHistogram(count, G);
}

void init() {
	// do nothing
}

void finalize() {
	// do nothing
}

/** template specialization for export */
template int totalMixingTime<float>(const StateGraph*, const float, device_t);
template int totalMixingTime<double>(const StateGraph*, const double, device_t);
template float secondLargestEigenvalue(const StateGraph*);
template double secondLargestEigenvalue(const StateGraph*);

}
