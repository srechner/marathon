/*
 * Random.cpp
 *
 * Created on: May 11, 2016
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2016, Steffen Rechner
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "../../include/marathon/Random.h"

#include <iostream>
#include <climits>

std::default_random_engine marathon::Random::rng;
std::uniform_real_distribution<double> marathon::Random::real_dist(0.0, 1.0);
std::uniform_int_distribution<int> marathon::Random::int_dist(0, INT_MAX);

void marathon::Random::init() {
	std::random_device rd;
	rng.seed(rd());
}

void marathon::Random::cleanup() {

}

double marathon::Random::nextDouble() {
	const double r = real_dist(rng);
	return r;
}

int marathon::Random::nextInt(int a, int b) {
	const int N = b - a;
	//const double r =  nextDouble();
	const int res = a + int_dist(rng) % N;
	//std::cout << "a=" << a << " b=" << b << " r=" << r << " N=" << N << " res=" << res << std::endl;
	return res;
}

int marathon::Random::nextInt(int b) {
	return nextInt(0, b);
}

void marathon::Random::select(int n, int k, int *selection) {

	/**************************************************************************
	 * Generate Random Combination of k out of n numbers.
	 * Use Selection Sampling (see Knuth - TAoCP Section 3.4.2 Algorithm S)
	 *************************************************************************/

	int t, m;
	double U;

	m = 0;
	t = 0;

	while (m < k) {

		U = nextDouble();    // U is uniformly distributed between 0 and 1

		// select t+1 with probability (n-m)/(N-t)
		if ((n - t) * U < (k - m)) {
			selection[m] = t;
			m++;
		}

		t++;
	}
}
