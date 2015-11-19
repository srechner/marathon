/*
 * marathon.cpp
 *
 *  Created on: Nov 17, 2015
 *      Author: steffen
 */

#include "../../include/marathon/marathon.h"
#include "../../include/marathon/cpu/analyzer.h"
#include "../../include/marathon/gpu/analyzer.h"
#include "../../include/marathon/hybrid/analyzer.h"

namespace marathon {

// indicates, whether library has already be initialized
bool isInit = false;

template<typename T>
int totalMixingTime(const StateGraph* mc, const T epsilon, device_t device) {
	switch (device) {
	case CPU_ONLY:
		return cpu::totalMixingTime(mc, epsilon);
		break;
	case GPU_ONLY:
		// library needs to be initialized
		if (!isInit)
			init();
		return gpu::totalMixingTime(mc, epsilon);
		break;
	case HYBRID:
		// library needs to be initialized
		if (!isInit)
			init();

		return hybrid::totalMixingTime(mc, epsilon);
		break;
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
	marathon::gpu::init();
	marathon::hybrid::init();
	isInit = true;
}

void finalize() {
	marathon::gpu::finalize();
	marathon::hybrid::init();
	isInit = true;
}

/** template specialization for export */
template int totalMixingTime<float>(const StateGraph*, const float, device_t);
template int totalMixingTime<double>(const StateGraph*, const double, device_t);
template float secondLargestEigenvalue(const StateGraph*);
template double secondLargestEigenvalue(const StateGraph*);

}
