/*
 * Diameter.h
 *
 * Created on: Sep 23, 2016
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

#ifndef PROJECT_DIAMETER_H
#define PROJECT_DIAMETER_H

#include "StateGraph.h"

namespace marathon {

	class Diameter {

		/**
		 * Computes the diameter of the graph, i.e. the maximal length of a shortest path
		 * between some nodes of the graph. Each arc has a length of 1.
		 */
		static
		int
		diameter(const StateGraph *G);

		/**
		 * Computes a histogram of the length of shortest path in G. Each arc of the graph
		 * contributes one to the length of a path. For each length l for l in 0..diameter(G),
		 * the number of shortest paths of G is stored in the vector count[l].
		 */
		static
		void
		pathLengthHistogram(std::vector<long> &count, const StateGraph *G);

	};

}


#endif //PROJECT_DIAMETER_H
