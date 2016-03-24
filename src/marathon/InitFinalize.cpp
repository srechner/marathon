/*
 * InitFinalize.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

#include "../../include/marathon/marathon.h"

namespace marathon {

#ifdef CUDA
// declare external function that initialize and finalze cuda stuff
extern void cudaInit();
extern void cudaFinalize();
#endif

void init() {

#ifdef CUDA
	cudaInit();
#endif
}

void finalize() {

#ifdef CUDA
	cudaFinalize();
#endif

}

}
