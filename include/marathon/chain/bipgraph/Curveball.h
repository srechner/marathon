/*
 * Curveball.h
 *
 * Created on: Aug 24, 2015
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

#ifndef CURVEBALL_H_
#define CURVEBALL_H_

#include "../bipgraph/SwitchChain.h"

namespace marathon {
namespace chain {
namespace bipgraph {

class Curveball: public SwitchChain {
public:
	Curveball(const std::string& line);
	Curveball(const std::vector<int>& u, const std::vector<int>& v);
	Curveball(int* const rowsum, int* const colsum, const int nrow, const int ncol);

	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

	void randomize(State* s, const uint32_t t=1) const;

	/**
	 * Return the loop probability of State s.
	 */
	virtual rational loopProbability(const State* s) const;

protected:

	/**
	 * Backtrack all subsets of size n from population. Each subset is stored in list.
	 */
	void backtrackSubsets(boost::dynamic_bitset<>& population, int n,
			const rational& p, const boost::dynamic_bitset<>& Bi,
			const boost::dynamic_bitset<>& Bj, const boost::dynamic_bitset<>& X,
			const BinaryMatrix* s, const int i, const int j,
			std::vector<std::pair<State*, rational>>& myneighbours) const;

	// just another recursive backtracking procedure
	void backtrackSubsetsRecursive(const boost::dynamic_bitset<>& population,
			boost::dynamic_bitset<>& tmp, int pos, int m, int n,
			const boost::dynamic_bitset<>& Bi,
			const boost::dynamic_bitset<>& Bj, const boost::dynamic_bitset<>& X,
			const BinaryMatrix* s, const int i, const int j,
			const rational& p,
			std::vector<std::pair<State*, rational>>& myneighbours) const;

};

}

}

}

#endif /* CURVEBALL_H_ */
