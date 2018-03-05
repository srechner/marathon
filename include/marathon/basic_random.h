/*
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

#ifndef INCLUDE_MARATHON_BASIC_RANDOM_H_
#define INCLUDE_MARATHON_BASIC_RANDOM_H_

#include <random>
#include <marathon/integer.h>

namespace marathon {

    /**
     * Basic Random Functions.
     */
    class BasicRandom {

    private:

        static std::mutex m;
        static std::mt19937 rng_init;            // pseudo random number generator
        static std::random_device rd;           // hardware random number generator (expensive)

        std::mt19937 rng;                        // random number generator for this object
        std::uniform_real_distribution<double> real_dist;
        std::uniform_int_distribution<int> int_dist;

    public:

        /**
         * Create a Random Generator Object.
         */
        BasicRandom() {
            m.lock();
            rng.seed(rng_init());
            m.unlock();
        }

        /**
         * Return a uniformly distributed real number from the interval [0,1).
          */
        double nextDouble() {
            const double r = real_dist(rng);
            return r;
        }

        /**
          * Return a uniformly distributed integer from the interval [a,b).
          */
        int nextInt(int a, int b) {
            const int N = b - a;
            //const double r =  nextDouble();
            const int res = a + (int_dist(rng) % N);
            //std::cout << "a=" << a << " b=" << b << " r=" << r << " N=" << N << " res=" << res << std::endl;
            return res;
        }

        /**
         * Return a uniformly distributed integer from the interval [0,b).
         */
        int nextInt(int b) {
            const int r = nextInt(0, b);
            return r;
        }

        /**
         * Return a uniformly distributed integer from the interval [0,b).
         */
        Integer nextInt(const Integer &b) {

            // check trivial cases
            if (b <= 1)
                return 0;

            // determine number of bits to represent an integer in [0,b)
            const int k = msb(b) + 1;

            Integer res(0); // will be returned

            // create random integer with k bits until a number smaller than b is found.
            // the expected number of iterations is less than two.
            do {

                // create a random integer with k bits
                res = 0;

                const int c = 30;
                const int cc = 1 << c;

                // create k random bits in batches of size c
                int i = 0;
                while (i + c < k) {
                    const int u = nextInt(cc);
                    res = (res * cc) + u;
                    i += c;
                }

                // the remaining batch has a size of d <= c
                const int d = k - i;
                const int dd = 1 << d;
                const int u = nextInt(dd);
                res = (res * dd) + u;

            } while (res >= b);

            return res;
        }


        /**
         * Randomly select k integers from the range [0..n).
         * @param dst: An integer array of size k where the selected number are stored.
         * @param n: The number of integers to choose from.
         * @param k: The number of integers to choose.
         */
        void select(int *dst, const int n, const int k) {

            assert(k <= n);

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
                    dst[m] = t;
                    m++;
                }

                t++;
            }
        }


        /**
         * Select a subset of the set { src[0], src[1], ..., src[n-1] } of size k uniformly at random.
         * Each subset has a probability of 1/(binom(n,k)).
         * @tparam T Type of the objects.
         * @param dst Destination array of size k.
         * @param src Source array of size n.
         * @param n The number of objects to choose from.
         * @param k The number of objects to choose.
         */
        template<class T>
        void select(T *dst, const T *src, const int n, const int k) {

            /**************************************************************************
             * Generate Random Combination of k out of n numbers.
             * Use Selection Sampling (see Knuth - TAoCP Section 3.4.2 Algorithm S)
             *************************************************************************/

            size_t t, m;
            double U;

            m = 0;
            t = 0;

            while (m < k) {

                U = nextDouble();    // U is uniformly distributed between 0 and 1

                // select t+1 with probability (n-m)/(N-t)
                if ((n - t) * U < (k - m)) {
                    dst[m] = src[t];
                    m++;
                }

                t++;
            }
        }

        /**
         * Shuffle the array, i.e. create a random permutation.
         */
        template<class T>
        void shuffle(T *data, int size) {
            for (int i = size; i > 1; i--) {
                int r = nextInt(i);
                std::swap(data[i - 1], data[r]);
            }
        }


        /**
         * Select a subset of the set { src[0], src[1], ..., src[n-1] } uniformly at random.
         * Each subset is selected with probability 1/(2^n).
         * @param dst Destination array of size n.
         * @param src Source array of size n.
         * @param n Number of elements.
         * @return Size of the random subset.
         */
        template<class T>
        int subset(T *dst, const T *src, const int n) {

            int sz = 0;        // size of the subset
            int i = 0;         // index of the current element

            const int k = 30;     // k random bits are generated in one step

            // for each element of src
            while (i < n) {

                // the number of bits generated in each step
                int e = std::min(k, n - i);

                // create a random integer with e bits
                int x = nextInt(1 << e);      // 0 <= x < 2^e

                // for each bit
                for (int l = 0; l < e; l++) {

                    // if x mod 2 == 1
                    if (x & 1) {
                        // subset will contain src[i]
                        dst[sz] = src[i];
                        sz++;
                    }
                    x >>= 1;                // x = x/2
                    i++;
                }
            }

            return sz;
        }
    };
};

// initialize random devices
std::mutex marathon::BasicRandom::m;
std::random_device marathon::BasicRandom::rd;
std::mt19937 marathon::BasicRandom::rng_init(rd());

#endif /* INCLUDE_MARATHON_BASIC_RANDOM_H_ */
