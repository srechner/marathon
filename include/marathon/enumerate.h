/*
 * Created on: Oct 19, 2017
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2017, Steffen Rechner
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

#ifndef MARATHON_ENUMERATE_H
#define MARATHON_ENUMERATE_H

#include <numeric>
#include <algorithm>
#include <mutex>
#include <boost/unordered_map.hpp>

#include "marathon/state.h"

namespace marathon {

    /**
     * Abstract base class designed for counting the number of states of a Markov chain.
     */
    class Enumerator {

    public:

        virtual ~Enumerator() = default;

        /**
         * Enumerate all states. For each state, evoke the function f.
         * @param f Function that is evaluated for each state.
         */
        virtual void enumerate(const std::function<void(const State &)> f) = 0;
    };
}

#endif //MARATHON_BINARY_MATRIX_FIXED_MARGIN_COUNT_H
