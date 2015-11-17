/*
 * Chain.h
 *
 *  Created on: Nov 24, 2014
 *      Author: rechner
 */

#ifndef _BIP_DEG_CHAIN_H_
#define _BIP_DEG_CHAIN_H_

#include "../../markov_chain.hpp"
#include "dense_bipartite_graph.h"

#include <vector>

namespace marathon {

namespace chain {

namespace sequence {

/**
 * Implements the Markov chain defined by Kannan et al.
 */
class SwitchBipartite: public MarkovChain<DenseBipartiteGraph> {

public:

	SwitchBipartite(const std::string& line);
	virtual bool computeArbitraryState(DenseBipartiteGraph& s) const;
	virtual void computeNeighbours(const DenseBipartiteGraph& s,
			boost::unordered_map<DenseBipartiteGraph, Rational>& neighbors) const;
	virtual void canonicalPath(int from, int to, std::list<int>& path) const;

protected:

	std::vector<int> u;
	std::vector<int> v;
	int sum;

	// helper functions for canonical path construction
	int next_red_edge(int col, bool* red_edges, int m, int n) const;
	int next_blue_edge(int row, bool* blue_edges, int m, int n) const;
	void splice_cycle(std::vector<int> cycle,
			std::list<std::vector<int> >& cycles) const;
	void trace_cycle(bool* blue_edges, bool* red_edges, int m, int n, int i,
			int j, std::vector<int>& cycle) const;
	void cycle_decomposition(const DenseBipartiteGraph& x,
			const DenseBipartiteGraph& y,
			std::list<std::vector<int> >& cycles) const;
	struct cycle_comparator;

};

}

}

}

#endif /* CHAIN_H_ */
