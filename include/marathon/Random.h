/*
 * Random.h
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

#ifndef INCLUDE_MARATHON_RANDOM_H_
#define INCLUDE_MARATHON_RANDOM_H_

#include <random>

namespace marathon {

	namespace Random {

		namespace detail {
			static std::default_random_engine rng;        // Random Number Generator
			static std::uniform_real_distribution<double> real_dist;
			static std::uniform_int_distribution<int> int_dist;
			static uint32_t num_clients = 0;
		}


		/**
		 * Initialize the Random Component.
		 */
		static
		void init() {
			if (detail::num_clients == 0) {
				std::random_device rd;
				detail::rng.seed(rd());
			}
			detail::num_clients++;
		}

		/**
		 * Cleanup the Random Component.
		 */
		static
		void cleanup() {
			detail::num_clients--;
			if(detail::num_clients== 0) {
				// to notghin
			}
		}

		/**
		 * Return a random double of the intervall [0,1).
		 */
		static
		double nextDouble() {
			const double r = detail::real_dist(detail::rng);
			return r;
		}

		/**
		 * Return a random integer of the intervall [a,b).
		 */
		static
		int nextInt(int a, int b) {
			const int N = b - a;
			//const double r =  nextDouble();
			const int res = a + detail::int_dist(detail::rng) % N;
			//std::cout << "a=" << a << " b=" << b << " r=" << r << " N=" << N << " res=" << res << std::endl;
			return res;
		}

		/**
		 * Return a random integer of the intervall [0,b).
		 */
		static
		int nextInt(int b) {
			return nextInt(0, b);
		}

		/**
		 * Randomly select k integers from the range [0..n).
		 * @param n: The integers to choose from
		 * @param k: The number of integers to choose.
		 * @param selection: An integer array of size k where the selected number are stored.
		 */
		static
		void select(int n, int k, int *selection) {

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

		/**
		 * Shuffle the array, i.e. create a random permutation.
		 */
		template<class T>
		static
		void shuffle(T *data, int size) {
			for (int i = size; i > 2; i--) {
				int r = nextInt(i);
				std::swap(data[i - 1], data[r]);
			}
		}
	}
};

#endif /* INCLUDE_MARATHON_RANDOM_H_ */
