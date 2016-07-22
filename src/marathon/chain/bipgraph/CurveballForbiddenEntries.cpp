/*
 * swap.cpp
 *
 *  Created on: Jun 4, 2015
 *      Author: rechner
 */

#include "../../../../include/marathon/chain/bipgraph/CurveballForbiddenEntries.h"

#include <iostream>
#include <map>

#include "../../../../include/marathon/Random.h"
#include "../../../../include/marathon/Combinatorics.h"
#include "../../../../include/marathon/chain/bipgraph/HavelHakimi.h"

namespace marathon {
namespace chain {
namespace bipgraph {

CurveballForbiddenEntries::CurveballForbiddenEntries(const std::vector<int>& u,
		const std::vector<int>& v, bool* const forbidden) :
		u(u), v(v), forbidden(forbidden) {

}

State* CurveballForbiddenEntries::computeArbitraryState() const {

	// not a valid instance
	if (u.size() < 1 || v.size() < 1) {
		return nullptr;
	}

	// try to find a realization
	return HavelHakimiBipartiteForbidden(u, v, forbidden);
}

void CurveballForbiddenEntries::randomize(State* x, const uint32_t t) const {

	BinaryMatrix* s = (BinaryMatrix*) x;

	const int nrow = u.size();
	const int ncol = v.size();

	if (nrow == 1 || ncol == 1) {
		return;
	}

	for (int l = 0; l < t; l++) {

		// randomly select two row indices
		int tmp[ncol];
		::marathon::random::select(nrow, 2, tmp);
		int i = tmp[0];
		int j = tmp[1];

		if (i >= nrow || j >= nrow) {
			std::cerr << "error! " << i << " " << j << std::endl;
		}

		// select the indices that occur in on the Mi and Mj, but not in both
		int X[ncol];
		int a = 0;
		int b = 0;

		// for each column position
		for (int k = 0; k < ncol; k++) {

			if (!forbidden[i * ncol + k] && !forbidden[j * ncol + k]) {

				bool A_ik = s->get(i, k);
				bool A_jk = s->get(j, k);

				if (A_ik != A_jk) {

					// store index k at X array
					X[a + b] = k;

					if (A_ik)
						a++;
					else
						b++;
				}
			}
		}

		// select a out of (a+b) positions
		bool Y[a + b];
		memset(Y, false, (a + b) * sizeof(bool));
		::marathon::random::select(a + b, a, tmp);
		for (int k = 0; k < a; k++)
			Y[tmp[k]] = true;

		// for each integer in X
		for (int k = 0; k < a + b; k++) {
			s->set(i, X[k], Y[k]);
			s->set(j, X[k], !Y[k]);
		}
	}
}

}
}
}
