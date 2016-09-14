/*
 * InitFinalize.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

#include "../../include/marathon/marathon.h"
#include "../../include/marathon/Random.h"
#include "../../include/marathon/Combinatorics.h"

namespace marathon {

#ifdef CUDA

// declare external function that initialize and finalze cuda stuff
	extern bool cudaInit();

	extern void cudaFinalize();

	bool cuda_init_success = false;
#endif

	void init() {

#ifdef CUDA
		cuda_init_success = cudaInit();
#endif

		// init RNG
		std::random_device rd;
		random::rng.seed(rd());


		// init binomial coefficients
		combinatorics::_nrow = 1;
		combinatorics::_ncol = 1;
		combinatorics::_binom = new rational[combinatorics::_nrow * combinatorics::_ncol];
		combinatorics::_binom[0] = 1;

	}

	void finalize() {

#ifdef CUDA
		if (cuda_init_success)
			cudaFinalize();
#endif

		delete[] combinatorics::_binom;
	}

}
