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

namespace marathon {
namespace chain {
namespace bipgraph {

bool HavelHakimiBipartite(const std::vector<int>& u,
		const std::vector<int>& v, bool* M);

bool isRealizable(const int* rowsum, const int* colsum, const int nrow, const int ncol);

}
}
}

#endif /* HAVAL_HAKIMI_H_ */
