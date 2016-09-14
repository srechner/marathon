/*
 * Random.h
 *
 *  Created on: May 11, 2016
 *      Author: rechner
 */

#ifndef INCLUDE_MARATHON_RANDOM_H_
#define INCLUDE_MARATHON_RANDOM_H_

#include <random>

namespace marathon {
    namespace random {

        extern std::default_random_engine rng;        // Random Number Generator
        extern std::uniform_real_distribution<double> real_dist;
        extern std::uniform_int_distribution<int> int_dist;

        /**
         * Return a random double of the intervall [0,1).
         */
        double nextDouble();

        /**
         * Return a random integer of the intervall [a,b).
         */
        int nextInt(int a, int b);

        /**
         * Return a random integer of the intervall [0,b).
         */
        int nextInt(int b);

        /**
         * Randomly select k integers from the range [0..n).
         * @param n: The integers to choose from
         * @param k: The number of integers to choose.
         * @param selection: An integer array of size k where the selected number are stored.
         */
        void select(int n, int k, int *selection);


        /**
         * Shuffle the array, i.e. create a random permutation.
         */
        template<class T>
        void shuffle(T *data, int size) {
            for (int i = size; i > 2; i--) {
                int r = nextInt(i);
                std::swap(data[i - 1], data[r]);
            }
        }
    }
};

#endif /* INCLUDE_MARATHON_RANDOM_H_ */
