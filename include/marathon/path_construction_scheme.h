/*
 * PathConstructionScheme.h
 *
 * Created on: Mar 22, 2016
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

#ifndef INCLUDE_MARATHON_PATHCONSTRUCTIONSCHEME_H_
#define INCLUDE_MARATHON_PATHCONSTRUCTIONSCHEME_H_

#include "state_graph.h"

namespace marathon {

    /**
     * A virtual base class for construction schemes of Canonical Paths.
     */
    class PathConstructionScheme {

    public:

        /**
         * Construct a path between states s and t in Graph sg.
         * @param sg A pointer to a state graph object at which the path is embedded.
         * @param s The index of the paths start state.
         * @param t The index of the paths final state.
         * @param path A list of state indices that represent the path.
         */
        virtual std::list<int> construct(const StateGraph &sg, const int s, const int t) const = 0;

    };

}

#endif /* INCLUDE_MARATHON_PATHCONSTRUCTIONSCHEME_H_ */
