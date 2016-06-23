/*
 * Random.cpp
 *
 *  Created on: May 11, 2016
 *      Author: rechner
 */

#include "../../include/marathon/Random.h"

#include <random>
#include <iostream>
#include <climits>

namespace marathon {
namespace random {

std::default_random_engine rng;
std::uniform_real_distribution<double> real_dist(0.0, 1.0);
std::uniform_int_distribution<int> int_dist(0, INT_MAX);

double nextDouble() {
	const double r = real_dist(rng);
	return r;
}

int nextInt(int a, int b) {
	const int N = b - a;
	//const double r =  nextDouble();
	const int res = a + int_dist(rng) % N;
	//std::cout << "a=" << a << " b=" << b << " r=" << r << " N=" << N << " res=" << res << std::endl;
	return res;
}

int nextInt(int b) {
	return nextInt(0,b);
}

void select(int N, int n, int* selection) {

	/**************************************************************************
	 * Generate Random Combination of k out of n numbers.
	 * Use Selection Sampling (see Knuth - TAoCP Section 3.4.2 Algorithm S)
	 *************************************************************************/

	int t, m;
	double U;

	m = 0;
	t = 0;

	while (m < n) {

		U = nextDouble();	// U is uniformly distributed between 0 and 1

		// select t+1 with probability (n-m)/(N-t)
		if ((N - t) * U < (n - m)) {
			selection[m] = t;
			m++;
		}

		t++;
	}
}

}
}

