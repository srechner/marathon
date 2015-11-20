/*
 * init_finalize.cpp
 *
 *  Created on: Nov 20, 2015
 *      Author: rechner
 */

#include "../../../include/marathon/hybrid/transition_matrix.h"

namespace marathon {

namespace hybrid {

bool init() {
	return cuda::initCublasXt();
}

void finalize() {
	cuda::finalizeCublasXt();
}

}

}

