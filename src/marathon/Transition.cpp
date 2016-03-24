/*
 * transition.cpp
 *
 *  Created on: Jun 13, 2015
 *      Author: rechner
 */

#include "../../include/marathon/Transition.h"

namespace marathon {

Transition::Transition() :
		u((uint) -1), v((uint) -1), p(0) {
}

Transition::Transition(uint u, uint v, rational p) :
		u(u), v(v), p(p) {
}

Transition::~Transition() {

}



bool TransitionComparator::operator()(const Transition& a,
		const Transition& b) {
	return a.u == b.u ? a.v < b.v : a.u < b.u;
}

}
