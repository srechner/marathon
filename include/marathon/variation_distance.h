/*
 * Created on: Oct 19, 2017
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

#ifndef MARATHON_VARIATION_DISTANCE_H_
#define MARATHON_VARIATION_DISTANCE_H_

#include "transition_matrix.h"

namespace marathon {

    /**
	 * Calculate the variation distance of the two vectors p1 and p2.
	 * @tparam T One of the following: float, double, Rational.
	 * @param p1 Vector of length n.
	 * @param p2 Vector of length n.
	 * @return Variation Distance between p1 and p2.
	 */
    template<class T>
    T variationDistance(const std::vector<T> &p1, const std::vector<T> &p2) {

        if (p1.size() != p2.size())
            throw std::runtime_error("Error! Vectors of unequal length!");

        T sum(0);
        for (size_t j = 0; j < p1.size(); j++) {
            T x = p1[j] - p2[j];
            if (x >= T(0))
                sum += x;
            else
                sum -= x;
        }

        return sum / T(2);
    }

    /**
     * Calculate the total variation distance of matrix P
     */
    template<typename T>
    T totalVariationDistance(const TransitionMatrix<T> &P, const std::vector<T> &pi) {

        const size_t omega = P.getDimension();

        T max(0);

#pragma omp parallel for if(omega > 100)
        for (size_t i = 0; i < omega; i++) {

            T sum(0);
            for (size_t j = 0; j < omega; j++) {
                T pij = P.get(i, j);
                T x = pij - pi[j];
                if (x >= 0)
                    sum += x;
                else
                    sum -= x;
            }

#pragma omp critical
            max = sum > max ? sum : max;
        }

        return max / T(2);
    }


    /**
     * Computes the variation distance of two given normalized histograms.
     * @param hist0 The first histogram.
     * @param hist1 The second histogram.
     * @param bins Number of bins.
     * @return The variation distance value.
     */
    template<class T1, typename T2 = double>
    T2 variationDistance(
            const std::unordered_map<T1, T2> &hist0,
            const std::unordered_map<T1, T2> &hist1,
            const int bins = INT_MAX
    ) {

        std::unordered_map<T1, T2> hist_temp(hist0);
        hist_temp.insert(hist1.begin(), hist1.end());

        // determine number of bins
        const int k = std::min((int) hist_temp.size(), bins);

        // error handling
        if (k == 0)
            return -1;
        else if (k == 1)
            return 0;

        // determine minimum and maximum
        T2 min, max;
        min = max = hist_temp.begin()->first;
        for (const auto &kv : hist_temp) {
            if (kv.first < min)
                min = kv.first;
            if (kv.first > max)
                max = kv.first;
        }

        // Distribution each observation to one of k bins of same size.

        // constant required for calculation of bin id
        const T2 z = (max - min) * T2(10001) / T2(10000);

        // build vectors of length k
        std::vector<T2> vec0(k);
        std::vector<T2> vec1(k);

        for (const auto &p : hist_temp) {

            // calculate bin in which p must be inserted
            const int bin_id = (p.first - min) * k / z;

            //std::cout << "bin(" << p.first << ")=" << bin_id << std::endl;

            auto v0 = hist0.find(p.first);
            auto v1 = hist1.find(p.first);
            if (v0 != hist0.end()) {
                vec0[bin_id] += v0->second;
            }
            if (v1 != hist1.end()) {
                vec1[bin_id] += v1->second;
            }
        }

        /*std::cout << "bin0:" << std::endl;
        for (int i = 0; i < k; i++)
            std::cout << vec0[i] << std::endl;
        std::cout << "bin1:" << std::endl;
        for (int i = 0; i < k; i++)
            std::cout << vec1[i] << std::endl;*/

        return marathon::variationDistance<T2>(vec0, vec1);
    }
}


#endif //MARATHON_VARIATION_DISTANCE_H_
