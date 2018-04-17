/*
 * Created on: Dec 17, 2016
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

#ifndef INCLUDE_MARATHON_COMBINATION_GENERATOR_H_
#define INCLUDE_MARATHON_COMBINATION_GENERATOR_H_

#include "marathon/integer.h"

namespace marathon {

    /**
     * Engine for the systematic enumeration of combinations.
     * @tparam T Data type of objects.
     */
    template<class T>
    class CombinationGenerator {

    private:

        const size_t n;         // the number of items
        const size_t k;         // the number of items to choose
        const T *src;           // array of n objects
        T *dst;                 // array of k objects (is reordered after each next)
        std::vector<size_t> c;  // temporary memory
        size_t j;               // smallest index s.t. c[j+1] > j
        size_t x;

    public:

        /**
         * Create a combination generator that generates all combinations of k out of n items.
         * @param src An array of n objects.
         * @param dst An array of k objects. This array is reodered after each call of next().
         * @param n The total number of objects.
         * @param k The number of objects to choose.
         * @return
         */
        CombinationGenerator(const T *src, T *dst, size_t n, size_t k) :
                n(n), k(k), src(src), dst(dst) {
            c.resize(k+3);
            reset();
        }

        /**
        * Reset the generator to its initial configuration.
        */
        void reset() {
            // initialize
            for (j = 1; j <= k; j++)
                c[j] = j - 1;
            c[k + 1] = n;
            c[k + 2] = 0;

            // generate output
            for (size_t i = 1; i <= k; i++)
                dst[i - 1] = src[c[i]];
        }

        /**
         * Generate the next combination.
         * @return True, if a new combination has been generated, false if all combination have been generated.
         */
        bool next() {

            /***********************************************************************************
             * Generate all combinations of k out of n items with Algorihtm L of TAoCP -
             * Combinatorial Algorithms Part I
             **********************************************************************************/

            j = 1;
            while (c[j] + 1 == c[j + 1]) {
                c[j] = j - 1;
                j++;
            }

            if (j > k)
                return false;

            c[j] = c[j] + 1;

            // generate output
            for (size_t i = 1; i <= k; i++)
                dst[i - 1] = src[c[i]];

            return true;

#ifdef false
            /***********************************************************************************
             * Generate all combinations of k out of n items with Algorihtm T of TAoCP -
             * Combinatorial Algorithms Part I
             **********************************************************************************/

            // j is the smallest index such that c[j+1] > j

            if (j > 0) {
                x = j;
                c[j] = x;                // increase c[j]
                j--;
            }
            else if (c[1] + 1 < c[2]) {   // easy case
                c[1]++;
            }
            else {                       // find j
                j = 1;
                do {
                    j++;
                    c[j - 1] = j - 2;
                    x = c[j] + 1;
                } while (x == c[j + 1]);

                if (j > k) {
                    return false;
                }
                else {                  // increase c[j]
                    c[j] = x;
                    j--;
                }
            }

            // generate output
            for (int i = 1; i <= k; i++)
                dst[i - 1] = src[c[i]];

            return true;
#endif
        }
    };
}


#endif /* INCLUDE_MARATHON_PERMUTATION_GENERATOR_H_ */
