/*
 * State.h
 *
 *  Created on: Nov 24, 2014
 *      Author: rechner
 */

#ifndef STATE_H_
#define STATE_H_

#include <cstdlib>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

#include "../../State.h"

namespace marathon {
namespace chain {
namespace bipgraph {

/**
 * This class represents a 0-1-matrix of size nrow times ncol.
 * It can be interpreted as biadjacency matrix of a bipartite graph.
 */
class BinaryMatrix: public State {

public:

	int nrow, ncol;		// number of rows and columns
	boost::dynamic_bitset<> M;// the matrix is stored as a single row in a bitset.

	BinaryMatrix();
	BinaryMatrix(int nrows, int ncols, const bool* bits = nullptr);
	BinaryMatrix(int nrows, int ncols, const std::string& str);
	BinaryMatrix(const BinaryMatrix& s);
	virtual ~BinaryMatrix();

	int getNumRows() const;
	int getNumColumns() const;
	bool get(int u, int v) const;
	void set(int u, int v, bool);
	void flip(int u, int v);

	size_t hashValue() const;
	int compare(const State* x) const;
	std::string toString() const;
	State* copy() const;

	bool is_switchable(int u1, int u2, int v1, int v2) const;
	void switch_4_cycle(int u1, int u2, int v1, int v2);

	void get_row(int u, boost::dynamic_bitset<>& row) const;

	void operator=(const BinaryMatrix& s);
	bool operator<(const BinaryMatrix &rhs) const;
	bool operator==(const BinaryMatrix &rhs) const;

};

}
}
}

#endif /* STATE_H_ */
