/*
 * haval_hakimi.h
 *
 *  Created on: Dec 5, 2014
 *      Author: rechner
 */

#ifndef HAVAL_HAKIMI_H_
#define HAVAL_HAKIMI_H_

#include <vector>
#include <cstdlib>

#include "BinaryMatrix.h"

namespace marathon {
namespace chain {
namespace bipgraph {

/**
 * Construct a bipartite graph with vertex degrees sequence u and v.
 */
BinaryMatrix* HavelHakimiBipartite(const std::vector<int>& u,
		const std::vector<int>& v);

/**
 * Construct a bipartite graph with vertex degrees sequence u and v
 * where all entries that are marked as 'forbidden' are set to zero.
 * @param u Degree Sequence of first bipartition set.
 * @param v Degree Sequence of second bipartition set.
 * @param forbidden binary matrix of size |u|*|v|.
 */
BinaryMatrix* HavelHakimiBipartiteForbidden(const std::vector<int>& u,
		const std::vector<int>& v, bool* const forbidden);

}
}
}

#endif /* HAVAL_HAKIMI_H_ */
