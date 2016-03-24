/*
 * PathCongestion.h
 *
 *  Created on: Mar 24, 2016
 *      Author: rechner
 */

#ifndef INCLUDE_MARATHON_PATHCONGESTION_H_
#define INCLUDE_MARATHON_PATHCONGESTION_H_

#include "marathon.h"

namespace marathon {
namespace pathCongestion {

rational pathCongestion(const StateGraph* sg,
		const PathConstructionScheme& pcs);

}
}



#endif /* INCLUDE_MARATHON_PATHCONGESTION_H_ */
