/*
 * Broder86.h
 *
 * Created on: Nov 18, 2014
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

#ifndef JS89_CHAIN_H_
#define JS89_CHAIN_H_

#include <queue>

#include "../../MarkovChain.h"
#include "BipartiteMatching.h"
#include "SparseBipartiteGraph.h"

namespace marathon {
namespace chain {
namespace matching {

class Broder86: public MarkovChain {

	friend class JS89Path;

protected:

	SparseBipartiteGraph *g = nullptr;

	/**
	 * Instances have the form "110101011".
	 * Such a 0-1-String is interpreted as a biadjacency matrix of a bipartite graph, flattened to a single line.
	 * Thus, the input string above corresponds to the biadjacency  matrix
	 *
	 *  1 1 0
	 *  1 0 1
	 *  0 1 1
	 *
	 *  which is the graph
	 *
	 *  u1  u2  u3
	 *  |\ / \ /|
	 *  | X   X |
	 *  |/ \ / \|
	 *  v1  v2  v3
	 * .
	 */
	void parseInstance(const std::string& inst);

public:

	Broder86(const std::string& instance);
	~Broder86();

	virtual State* computeArbitraryState() const;
	virtual void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

};

}

}

}
#endif /* JS89_H_ */

