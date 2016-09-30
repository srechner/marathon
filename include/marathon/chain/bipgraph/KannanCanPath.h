/*
 * KannanCanPath.h
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

#ifndef INCLUDE_MARATHON_CHAINS_SEQUENCES_KANNANCANPATH_H_
#define INCLUDE_MARATHON_CHAINS_SEQUENCES_KANNANCANPATH_H_

#include "../../PathConstructionScheme.h"
#include "BinaryMatrix.h"

namespace marathon {
namespace chain {
namespace bipgraph {

class KannanPath: public PathConstructionScheme {

protected:

	int next_red_edge(int col, bool* red_edges, int m, int n) const;
	int next_blue_edge(int row, bool* blue_edges, int m, int n) const;
	void trace_cycle(bool* blue_edges, bool* red_edges, int m, int n, int i,
			int j, std::vector<int>& cycle) const;
	void splice_cycle(std::vector<int> cycle,
			std::list<std::vector<int> >& cycles, const int m,
			const int n) const;
	void cycle_decomposition(
			const BinaryMatrix& x,
			const BinaryMatrix& y,
			std::list<std::vector<int> >& cycles) const;
	struct cycle_comparator {
		bool operator()(const std::vector<int>& c1,
				const std::vector<int>& c2) {
			return c1.size() < c2.size();
		}
	};

public:

	virtual void construct(const StateGraph* sg, const int s, const int t,
			std::list<int>& path) const;

};

}
}
}

#endif /* INCLUDE_MARATHON_CHAINS_SEQUENCES_KANNANCANPATH_H_ */
