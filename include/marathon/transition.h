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
#include "rational.h"

namespace marathon {

/**
 * Transition Arc Representation of State Graph
 */
class Transition {

public:

	uint u, v;		// from, to
	Rational p;		// P(u,v)

	Transition();
	Transition(uint u, uint v, Rational p);

	virtual ~Transition();

	typedef std::vector<Transition>::iterator transition_iterator;
};

struct TransitionComparator {
	bool operator()(const Transition& a, const Transition& b);
};

}

#endif /* TRANSITION_H_ */

