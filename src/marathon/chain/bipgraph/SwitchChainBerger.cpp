/*
 * SwitchChainBerger.cpp
 *
 * Created on: Jun 4, 2015
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

#include <iostream>

#include "../../../../include/marathon/chain/bipgraph/SwitchChainBerger.h"

namespace marathon {
namespace chain {
namespace bipgraph {

SwitchChainBerger::SwitchChainBerger(const std::string& input) :
		SwitchChain(input) {
}

void SwitchChainBerger::computeNeighbours(const State* x,
		std::vector<std::pair<State*, rational>>& neighbors) const {

	uint i, j, k, l, m, n;
	BinaryMatrix s2;
	int numNonAdjacentEdgePairs = 0;
	rational loop(0);

	const BinaryMatrix* s = (const BinaryMatrix*) x;

	/**
	 * Definition of Annabell Berger
	 *
	 * Let M be a set of all non-adjacent pairs of edges of s,
	 * extended by the artificial edge pair (a,a) which is necessary
	 * to get a ergodic chain.
	 *
	 * Choose a pair (e,e') of edges from M uniformly and random.
	 * Let these pairs be e=(i,j) and e'=(k,l) with i!=k and j!=l.
	 * In Practice this would be done by an efficient data strutcure.
	 *
	 * If the edges (i,l) and (k,j) are missing, switch the cycle.
	 */

	m = s->getNumRows();
	n = s->getNumColumns();

	/***************************************************************
	 * First Step: Count number of "natural" non-adjacent edge pairs
	 **************************************************************/

	// choose edge e=(i,j)
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {

			// Choose edge e'=(k,l)
			for (k = i + 1; k < m; k++) {
				for (l = 0; l < n; l++) {

					// if e and e' are actual edges
					if (s->get(i, j) && s->get(k, l)) {

						// if e and e' are non adjacent
						if (l != j)
							numNonAdjacentEdgePairs++;

					}
				}
			}
		}
	}

	// add an additional edge a, which is non-adjacent to all
	numNonAdjacentEdgePairs += sum;

	const rational p(1, numNonAdjacentEdgePairs);

	/*****************************************************
	 * Second Step: Select two non-adjacent edge pairs
	 ****************************************************/

	// choose edge e=(i,j)
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {

			// choose second edge e'=(k,l), different from e
			for (k = i + 1; k < m; k++) {
				for (l = 0; l < n; l++) {

					// if e and e' are non adjacent
					if (l != j) {

						// if e and e' are actual edges
						if (s->get(i, j) && s->get(k, l)) {

							BinaryMatrix* s2 = new BinaryMatrix(
									*s);

							// alternating cycle ( (i,j)=1, (k,j)=0, (k,l)=1, (i,l)=0 )
							if (!s->get(i, l) && !s->get(k, j)) {

								// switch the cycle
								s2->flip(i, j); // (i,j) = 0
								s2->flip(k, l); // (k,l) = 0
								s2->flip(i, l); // (i,l) = 1
								s2->flip(k, j); // (k,j) = 1
								//std::cout << "a -> " << s2 << " (" << s << ")" << std::endl;
							}

							neighbors.push_back(std::make_pair(s2, p));
						}
					}
				}
			}
		}
	}

	// choose second edge as artificial edge
	loop += rational(sum, numNonAdjacentEdgePairs);

	BinaryMatrix *sl = new BinaryMatrix(*s);
	neighbors.push_back(std::make_pair(sl, loop));
}

}
}
}
