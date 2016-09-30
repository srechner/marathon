/*
 * HavalHakimi.h
 *
 * Created on: Dec 5, 2014
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

#ifndef HAVAL_HAKIMI_H_
#define HAVAL_HAKIMI_H_

#include <vector>
#include <cstdlib>

#include "BinaryMatrix.h"

namespace marathon {
namespace chain {
namespace bipgraph {

/**
 * Construct a bipartite graph with vertex degrees sequence u and v.
 */
BinaryMatrix* HavelHakimiBipartite(const std::vector<int>& u,
		const std::vector<int>& v);

/**
 * Construct a bipartite graph with vertex degrees sequence u and v
 * where all entries that are marked as 'forbidden' are set to zero.
 * @param u Degree Sequence of first bipartition set.
 * @param v Degree Sequence of second bipartition set.
 * @param forbidden binary matrix of size |u|*|v|.
 */
BinaryMatrix* HavelHakimiBipartiteForbidden(const std::vector<int>& u,
		const std::vector<int>& v, bool* const forbidden);

}
}
}

#endif /* HAVAL_HAKIMI_H_ */
