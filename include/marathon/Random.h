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

	class Random {

	private:

		// declare friend functions
		friend void marathon::init();
		friend void marathon::cleanup();

		static std::default_random_engine rng;        // Random Number Generator
		static std::uniform_real_distribution<double> real_dist;
		static std::uniform_int_distribution<int> int_dist;

	public:

		/**
		 * Initialize the Random Component.
		 */
		static
		void init();

		/**
		 * Cleanup the Random Component.
		 */
		static
		void cleanup();

		/**
		 * Return a random double of the intervall [0,1).
		 */
		static
		double nextDouble();

		/**
		 * Return a random integer of the intervall [a,b).
		 */
		static
		int nextInt(int a, int b);

		/**
		 * Return a random integer of the intervall [0,b).
		 */
		static
		int nextInt(int b);

		/**
		 * Randomly select k integers from the range [0..n).
		 * @param n: The integers to choose from
		 * @param k: The number of integers to choose.
		 * @param selection: An integer array of size k where the selected number are stored.
		 */
		static
		void select(int n, int k, int *selection);


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

	};

};

#endif /* INCLUDE_MARATHON_RANDOM_H_ */
