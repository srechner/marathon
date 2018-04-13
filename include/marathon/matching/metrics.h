/*
 * Created on: Jan 08, 2018
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

#ifndef MATCHING_METRICS_H_
#define MATCHING_METRICS_H_

#include "bipartite_matching.h"

namespace marathon {
    namespace matching {

        /**
         * Calculate the size of the symmetric difference of both matchings.
         * @param bin1 A bipartite matching.
         * @param bin2 A bipartite matching.
         * @return Size of symmetric difference.
         */
        size_t hammingDistance(const BipartiteMatching &m1, const BipartiteMatching &m2) {
            const size_t n = m1.getVertexSize();
            size_t d = 0;
            for (size_t i = 0; i < n; i++) {
                if (m1.getMate(i) != m2.getMate(i)) {
                    d++;
                }
            }
            return d;
        }
    }
}


#endif /* MATCHING_METRICS_H_ */
