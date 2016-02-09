/*
 * Chain.h
 *
 *  Created on: Nov 24, 2014
 *      Author: rechner
 */

#ifndef _BIP_DEG_CHAIN_H_
#define _BIP_DEG_CHAIN_H_

#include "../../common/markov_chain.hpp"
#include "dense_bipartite_graph.h"

#include <vector>

namespace marathon {

namespace chain {

namespace sequence {

/**
 * Implements the Markov chain defined by Kannan et al.
 */
class SwitchBipartite: public MarkovChain<DenseBipartiteGraph> {

protected:
	std::vector<int> u;
	std::vector<int> v;
	int sum;

public:

	SwitchBipartite(const std::string& inst);
	virtual bool computeArbitraryState(DenseBipartiteGraph& s);
	virtual void computeNeighbours(const DenseBipartiteGraph& s,
			std::unordered_map<DenseBipartiteGraph, rational>& neighbors) const;

	virtual void constructPath(const StateGraph* sg, int from, int to,
			std::list<int>& path) const;

	virtual void parseInstance(const std::string& line);

protected:

	struct cycle_comparator;

	// helper functions for canonical path construction
	int next_red_edge(int col, bool* red_edges, int m, int n) const;
	int next_blue_edge(int row, bool* blue_edges, int m, int n) const;
	void splice_cycle(std::vector<int> cycle,
			std::list<std::vector<int> >& cycles) const;
	void trace_cycle(bool* blue_edges, bool* red_edges, int m, int n,
			int i, int j, std::vector<int>& cycle) const;
	void cycle_decomposition(const DenseBipartiteGraph& x,
			const DenseBipartiteGraph& y, std::list<std::vector<int> >& cycles) const;

};

}

}

}

#endif /* CHAIN_H_ */
