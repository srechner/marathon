/*
 * types.h
 *
 *  Created on: Feb 11, 2015
 *      Author: steffen
 */

#ifndef TYPES_H_
#define TYPES_H_

// system includes
#include <cstdlib>
#include <ostream>
// boost includes

#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_set.hpp>
#include <boost/graph/adjacency_list.hpp>

// project includes

namespace marathon {

namespace chain {

namespace matching {

// undirected graph
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> undirected_graph;
typedef std::pair<uint, uint> edge;
typedef std::vector<edge> edgelist;


}

}

}

#endif /* TYPES_H_ */
