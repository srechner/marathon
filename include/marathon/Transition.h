/*
 * transition.h
 *
 *  Created on: Jun 4, 2015
 *      Author: rechner
 */

#ifndef TRANSITION_H_
#define TRANSITION_H_

#include <cstdlib>
#include <vector>

#include "Rational.h"

namespace marathon {

/**
 * Transition Arc Representation of State Graph
 */
class Transition {

public:

	uint u, v;			// from, to
	rational p;			// P(u,v)

	Transition();
	Transition(uint u, uint v, rational p);
	virtual ~Transition();
};

struct TransitionComparator {
	bool operator()(const Transition& a, const Transition& b);
};

}

#endif /* TRANSITION_H_ */

