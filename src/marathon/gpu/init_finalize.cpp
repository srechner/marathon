/*
 * init_finalize.cpp
 *
 *  Created on: Nov 20, 2015
 *      Author: rechner
 */

#include "../../../include/marathon/gpu/analyzer.h"

namespace marathon {

namespace gpu {

bool init() {
	return cuda::initCublas();
}

void finalize() {
	cuda::finalizeCublas();
}

}

}
