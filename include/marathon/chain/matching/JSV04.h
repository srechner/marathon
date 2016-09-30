/*
 * JSV04.h
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

#ifndef JSV04_CHAIN_H_
#define JSV04_CHAIN_H_

#include <queue>

#include "Broder86.h"

namespace marathon {
namespace chain {
namespace matching {

class JerrumSinclairVigoda04: public Broder86 {

protected:

	uint num_perfect_matching;
	uint *num_near_perfect_matching;

	rational getWeight(const State* s) const;

public:

	JerrumSinclairVigoda04(const std::string& input);

	~JerrumSinclairVigoda04();

	void computeNeighbours(const State* s,
			std::vector<std::pair<State*, rational>>& neighbors) const;

	void computeWeights(const std::vector<const State*>& states,
			std::vector<rational>& weights);
};

}
}
}

#endif /* JS89_H_ */
