/*
 * KannanPath.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#include "../../../../include/marathon/chain/bipgraph/KannanCanPath.h"

namespace marathon {
namespace chain {
namespace bipgraph {

int KannanPath::next_red_edge(int col, bool* red_edges, int m, int n) const {
	for (int i = 0; i < m; i++) {
		if (red_edges[i * n + col])
			return i;
	}
	// no edge found in column
	return -1;
}

int KannanPath::next_blue_edge(int row, bool* blue_edges, int m, int n) const {
	for (int j = 0; j < n; j++) {
		if (blue_edges[row * n + j])
			return j;
	}
	// no blue edge found in row
	return -1;
}

void KannanPath::trace_cycle(bool* blue_edges, bool* red_edges, int m, int n,
		int i, int j, std::vector<int>& cycle) const {

	while (j != -1) {

		// add (i,j) to cycle
		cycle.push_back(i);
		cycle.push_back(m + j);

		// remove blue edge (i,j)
		blue_edges[i * n + j] = false;

		i = next_red_edge(j, red_edges, m, n);

		// remove red edge (i,j)
		red_edges[i * n + j] = false;

		j = next_blue_edge(i, blue_edges, m, n);
	};
}

void KannanPath::splice_cycle(std::vector<int> cycle,
		std::list<std::vector<int> >& cycles, const int m, const int n) const {

#ifdef DEBUG
	std::cout << "splice cycle [ ";
	for (std::vector<int>::iterator it = cycle.begin(); it != cycle.end();
			++it) {
		std::cout << *it << " ";
	}
	std::cout << "]" << std::endl;
#endif

	std::vector<int> c;
	bool removed[cycle.size()];
	int first_seen[n + m];
	int i, j;

	memset(removed, 0, cycle.size() * sizeof(bool));

	for (i = 0; i < n + m; i++)
		first_seen[i] = n + m;

	for (i = 0; i < cycle.size(); i++) {

#ifdef DEBUG
		std::cout << i << ": " << cycle[i] << std::endl;
#endif

		// smaller cycle detected
		if (first_seen[cycle[i]] != n + m) {

#ifdef DEBUG
			std::cout << "smaller cycle detected" << std::endl;
#endif

			// extract smaller cycle and store in cycle list
			c.clear();
			for (j = first_seen[cycle[i]]; j != i; j++) {
				if (!removed[j]) {
					c.push_back(cycle[j]);
					first_seen[cycle[j]] = n + m;
					removed[j] = true;
				}
			}
			cycles.push_back(c);
		}
		first_seen[cycle[i]] = i;
	}

	// not removed vertices
	c.clear();
	for (i = 0; i < cycle.size(); i++) {
		if (!removed[i])
			c.push_back(cycle[i]);
	}
	cycles.push_back(c);
}

void KannanPath::cycle_decomposition(
		const BinaryMatrix& x,
		const BinaryMatrix& y,
		std::list<std::vector<int> >& cycles) const {

	const int nrow = x.getNumRows();
	const int ncol = x.getNumColumns();

	bool red[nrow * ncol];
	bool blue[nrow * ncol];

	std::vector<int> cycle;
	std::vector<int>::iterator cit;
	std::list<std::vector<int> >::iterator it;

	memset(red, 0, nrow * ncol * sizeof(bool));
	memset(blue, 0, nrow * ncol * sizeof(bool));

	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			if (x.get(i, j) && !y.get(i, j))
				blue[i*ncol+j] = true;
			else if (!x.get(i, j) && y.get(i, j))
				red[i*ncol+j] = true;
		}
	}

	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			if (blue[i*ncol+j]) {
				// start of alternating Cycle in x found
				cycle.clear();
				// trace cycle
				trace_cycle(blue, red, nrow, ncol, i, j, cycle);
				// try to splice cycles into smaller ones
				splice_cycle(cycle, cycles, nrow, ncol);
			}
		}
	}
}

void KannanPath::construct(const StateGraph* sg, const int s, const int t,
		std::list<int>& path) const {

#ifdef DEBUG
	std::cout << "from=" << states[from] << std::endl;
	std::cout << "to  =" << states[to] << std::endl;
#endif

	// start with from
	path.push_back(s);

	// continously modify u
	BinaryMatrix u = *((BinaryMatrix*) sg->getState(s));
	BinaryMatrix v = *((const BinaryMatrix*) sg->getState(t));

	std::list<std::vector<int> > cycles;
	std::list<std::vector<int> >::iterator it;

	int i, l;
	std::vector<int> w;
	std::vector<int>::iterator it2;

	// decompose symmetric difference of u and v into cycles
	cycle_decomposition(u, v, cycles);

#ifdef DEBUG
	std::cout << "found " << cycles.size() << " cycles" << std::endl;
#endif

	// sort cycles by length
	cycles.sort(cycle_comparator());

	while (!cycles.empty()) {

		// deal with smallest cycle
		w = cycles.front();
		cycles.pop_front();

		assert(w.size() % 2 == 0);
		assert(w.size() >= 4);
#ifdef DEBUG
		std::cout << "[";
		for (it2 = w.begin(); it2 != w.end(); ++it2) {
			std::cout << " " << *it2;
		}
		std::cout << " ]" << std::endl;
#endif

		l = w.size();

		/*
		 * find first vertex w[i] s.t. the cycle
		 *     (w[0], w[1], ..., w[i], w[i+1], w[l-i-2], w[l-i-1], ... , w[l-1])
		 * is alternating in u
		 */
		i = 0;
		while (!u.is_switchable(w[i], w[i + 1], w[l - i - 2], w[l - i - 1])) {
			i++;
		}

#ifdef DEBUG
		std::cout << "switch cycle: [ " << w[i] << " " << w[i + 1] << " "
		<< w[l - i - 2] << " " << w[l - i - 1] << " ]" << std::endl;
#endif

		/**
		 * (w[i], w[i+1], w[l-i-2], w[l-i-1]) is switchable
		 *    switch it!
		 */
		u.switch_4_cycle(w[i], w[i + 1], w[l - i - 2], w[l - i - 1]);
		int x = sg->indexOf(&u);
		path.push_back(x);

		/**
		 * two smaller alternating cycles are created:
		 *  1: (w[0], ... , w[i], w[l-i-1], ... , w[l-1])
		 *  2: (w[i+1], ... , w[l-i-2])
		 *
		 *  Erase first cycle by series of switches
		 */

		for (int j = i; j > 0; j--) {
			u.switch_4_cycle(w[j - 1], w[j], w[l - j - 1], w[l - j]);
			int x = sg->indexOf(&u);
			path.push_back(x);
		}

		/**
		 * We are left with the second cycle which is smaller
		 * Continue with this cycle
		 */

		// if second cycle has at least 4 edges
		if (l - 2 * i - 2 >= 4) {
			w.clear();
			for (int j = i + 1; j <= l - i - 2; j++) {
				w.push_back(w[j]);
			}
			cycles.push_front(w);
		}
	}
}

}
}
}
