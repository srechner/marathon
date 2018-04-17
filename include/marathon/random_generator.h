/*
 * Created on: May 19, 2017
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

#ifndef MARATHON_RANDOM_GENERATOR_H
#define MARATHON_RANDOM_GENERATOR_H

#include "omp.h"
#include <functional>
#include "marathon/state.h"
#include "marathon/rational.h"

namespace marathon {

    /**
     * Base class for Random Generators.
     */
    class RandomGenerator {

    public:

        /**
         * Return a random state.
         * @return Random state.
         */
        virtual const State &next() = 0;


        /**
         * Create a copy of the random generator.
         * @return A random generator that is identical to the current one.
         */
        virtual std::unique_ptr<RandomGenerator> copy() const = 0;

    };
}


#endif /* MARATHON_RANDOM_GENERATOR_H */
