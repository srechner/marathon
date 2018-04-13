/*
 * Created on: Oct 20, 2017
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

#ifndef MARATHON_BINARY_MATRIX_RANDOM_GENERATOR_H
#define MARATHON_BINARY_MATRIX_RANDOM_GENERATOR_H

// marathon includes
#include "marathon/random_generator.h"
#include "marathon/binary_matrix/binary_matrix.h"

namespace marathon {
    namespace binary_matrix {

        /**
         * An abstract base class for random generators producing
         * binary matrices.
         */
        class RandomGenerator : public marathon::RandomGenerator {

        public:

            /**
             * Return a random binary matrix.
             * @return Random binary matrix.
             */
            virtual const BinaryMatrix& next() override = 0;

        };
    }
}

#endif //MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_RANDOM_GENERATOR_EXACT_H
