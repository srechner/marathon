/*
 * BipartiteMatching.h
 *
 * Created on: Nov 20, 2014
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

#ifndef STATE_JS89_H_
#define STATE_JS89_H_

#include <cstring>
#include <cstdlib>
#include <ostream>

#include "../../State.h"

namespace marathon {
namespace chain {
namespace matching {

class BipartiteMatching: public State {

public:

	int n, k;				// number of nodes			// number of edges
	int unmatched[2];	// indices of unmatched nodes (if any)
	int* mates;

	BipartiteMatching();

	BipartiteMatching(const BipartiteMatching& s);
	BipartiteMatching(int n, int k, int unmatched[2], int* matching);
	virtual ~BipartiteMatching();

	void addEdge(int u, int v);
	void removeEdge(int u, int v);

	void operator=(BipartiteMatching const& s);
	bool operator==(const BipartiteMatching &s) const;
	bool operator<(const BipartiteMatching& s) const;
	size_t hashValue() const;
	int compare(const State*) const;
	std::string toString() const;
	State* copy()  const;

	bool is_perfect() const;
	bool is_near_perfect() const;
};

}
}
}


#endif /* STATE_H_ */

