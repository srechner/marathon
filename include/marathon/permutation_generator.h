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

#ifndef INCLUDE_MARATHON_PERMUTATION_GENERATOR_H_
#define INCLUDE_MARATHON_PERMUTATION_GENERATOR_H_

#include "marathon/integer.h"

namespace marathon {

    /**
     * Engine for the systematic generation of permuations.
     * @tparam T Data type of objects.
     */
    template<class T>
    class PermutationGenerator {

    private:

        const int n;            // the number of items
        const T *src;           // array of n objects
        T *dst;                 // array of n objects

    public:

        /**
         * Create a perumutation generator for n objects.
         * Let src_1, ..., src_n be objects. The generator sequentially generates
         * all permutations of the set { a1, ..., an }, visiting them in lexicographic order.
         * Each permutation is stored in the array dst.
         * @param src An array of n objects
         * @param dst An array of n objects
         * @param n The number of objects
         * @return
         */
        PermutationGenerator(const T *src, T *dst, const int n) : n(n), src(src), dst(dst) {
            reset();
        }

        ~PermutationGenerator() {

        }

        /**
         * Reset the generator to its initial configuration.
         */
        void reset() {
            // copy objects from src to dst
            memcpy(dst, src, n * sizeof(T));

            // initially sort data ascendingly by value
            std::sort(dst, dst + n);
        }

        /**
         * Generate next permutation
         * @return True, if a new permutation has been generated, false if all permutations have been generated.
         */
        bool next() {

            /***********************************************************************************
             * Generate lexicographic permutations with Algorihtm L of TAoCP -
             * Combinatorial Algorithms Part I, page 319
             **********************************************************************************/

            // 1. find the largest j such that aj can be increased
            int j = n - 2;
            while (dst[j] >= dst[j + 1]) {
                if (j == 0)
                    return false;
                else
                    j--;
            }

            // 2. increase aj by the smallest feasible amount
            int l = n - 1;
            while (dst[j] >= dst[l])
                l--;
            std::swap(dst[j], dst[l]);

            // 3. find the lexicographically least way to extend the new a1 ... aj to a complete pattern.
            int k = j + 1;
            l = n - 1;
            while (k < l) {
                std::swap(dst[k], dst[l]);
                k++;
                l--;
            }

            return true;
        }
    };
}

#endif /* INCLUDE_MARATHON_PERMUTATION_GENERATOR_H_ */
