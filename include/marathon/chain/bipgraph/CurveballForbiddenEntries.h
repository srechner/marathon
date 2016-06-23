/*
 * curveball.h
 *
 *  Created on: Aug 24, 2015
 *      Author: rechner
 */

#ifndef CURVEBALLFRFC_H_
#define CURVEBALLFRFC_H_

#include "../../MarkovChain.h"

namespace marathon {
namespace chain {
namespace bipgraph {

/**
 * Sampling 0-1-matrices where the entries of row and col have to be zero.
 */
class CurveballForbiddenEntries: public MarkovChain {

protected:

	const std::vector<int> u;
	const std::vector<int> v;
	bool* forbidden;

public:

	/**
	 * Creates random realizations that realize (u,v) and where the first k columns in
	 * row 0 are set to zero.
	 */
	CurveballForbiddenEntries(const std::vector<int>& u, const std::vector<int>& v,
			bool* const forbidden);

	virtual State * computeArbitraryState() const;

	void randomize(State* s, const uint32_t t = 1) const;

protected:

};

}
}
}

#endif /* CURVEBALLFRFC_H_ */
