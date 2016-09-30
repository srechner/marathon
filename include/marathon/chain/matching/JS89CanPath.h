/*
 * JS89CanPath.h
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

#ifndef INCLUDE_MARATHON_CHAINS_MATCHING_JS89CANPATH_H_
#define INCLUDE_MARATHON_CHAINS_MATCHING_JS89CANPATH_H_

#include "../../PathConstructionScheme.h"
#include "../../MarkovChain.h"

namespace marathon {
namespace chain {
namespace matching {

class JS89Path: public PathConstructionScheme {

	/**
	 * Path construction Scheme for canonical paths from Jerrum and Sinclair 1989.
	 * @param sg Pointer to a state graph object.
	 * @param s Index of a state (Start).
	 * @param t Index of a state (End).
	 * @param path List of State indices that corresponds to the state.
	 */
	virtual void construct(const StateGraph* sg, const int s, const int t,
			std::list<int>& path) const;

};

}
}
}

#endif /* INCLUDE_MARATHON_CHAINS_MATCHING_JS89CANPATH_H_ */
