/*
 * SwitchChain.h
 *
 * Created on: Nov 24, 2014
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

#ifndef _BIP_DEG_CHAIN_H_
#define _BIP_DEG_CHAIN_H_

#include "../../MarkovChain.h"
#include <vector>
#include <stack>

#include "BinaryMatrix.h"

namespace marathon {
namespace chain {
namespace bipgraph {

/**
 * Implements the Markov chain defined by Kannan et al.
 */
class SwitchChain: public MarkovChain {

	friend class KannanPath;

public:

	// variables and methods for general purpose
	std::vector<int> u;
	std::vector<int> v;
	int sum;

	/**
	 * Instances have the form "2,2,2;1,2,1,2".
	 * The semicolon separates two degree sequences of both bipartition sets.
	 */
	virtual void parseInstance(const std::string& line);

public:

	SwitchChain(const std::string& inst);
	SwitchChain(const std::vector<int>& u, const std::vector<int>& v);
	SwitchChain(int* const rowsum, int* const colsum, const int nrow, const int ncol);
	virtual ~SwitchChain();

	virtual State* computeArbitraryState() const;


	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

	/**
	 * Return the loop probability of State s.
	 */
	virtual rational loopProbability(const State* s) const;

	/**
	 * Randomize the state s by applying a single transition.
	 */
	virtual void randomize(State* s, const uint32_t t=1) const;
};

}
}
}

#endif /* CHAIN_H_ */
