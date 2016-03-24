/*
 * KannanPath.h
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#ifndef INCLUDE_MARATHON_CHAINS_SEQUENCES_KANNANCANPATH_H_
#define INCLUDE_MARATHON_CHAINS_SEQUENCES_KANNANCANPATH_H_

#include "DenseBipartiteGraph.h"
#include "../../PathConstructionScheme.h"

namespace marathon {
namespace chain {
namespace bipgraph {

class KannanPath: public PathConstructionScheme {

protected:

	int next_red_edge(int col, bool* red_edges, int m, int n) const;
	int next_blue_edge(int row, bool* blue_edges, int m, int n) const;
	void trace_cycle(bool* blue_edges, bool* red_edges, int m, int n, int i,
			int j, std::vector<int>& cycle) const;
	void splice_cycle(std::vector<int> cycle,
			std::list<std::vector<int> >& cycles, const int m,
			const int n) const;
	void cycle_decomposition(
			const DenseBipartiteGraph& x,
			const DenseBipartiteGraph& y,
			std::list<std::vector<int> >& cycles) const;
	struct cycle_comparator {
		bool operator()(const std::vector<int>& c1,
				const std::vector<int>& c2) {
			return c1.size() < c2.size();
		}
	};

public:

	virtual void construct(const StateGraph* sg, const int s, const int t,
			std::list<int>& path) const;

};

}
}
}

#endif /* INCLUDE_MARATHON_CHAINS_SEQUENCES_KANNANCANPATH_H_ */
