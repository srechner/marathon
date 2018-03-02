/*
 * Transition.h
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
#ifndef TRANSITION_H_
#define TRANSITION_H_

#include <cstdlib>
#include <vector>

#include "rational.h"

namespace marathon {

	/**
	 * Transition Arc Representation of State Graph
	 */
	class Transition {

	public:

		uint from, to;            	// from, to
		Rational weight;            // P(u,v)

		Transition() :
				from((uint) -1), to((uint) -1), weight(0) {
		}

		Transition(uint u, uint v, Rational p) :
				from(u), to(v), weight(p) {
		}

		virtual ~Transition() {

		}
	};

	struct TransitionComparator {
		bool operator()(const Transition &a, const Transition &b) {
			return a.from == b.from ? a.to < b.to : a.from < b.from;
		}
	};

}

#endif /* TRANSITION_H_ */

