// project includes
#include "../../../../include/marathon/chain/bipgraph/LindaChain.h"
#include "../../../../include/marathon/Random.h"

//#define DEBUG

namespace marathon {
namespace chain {
namespace bipgraph {

LindaChain::LindaChain(std::string& s) {
	// Todo: Implement constructor.
	nrow = 0;
	ncol = 0;
}

LindaChain::LindaChain(std::vector<int>& u_lower, std::vector<int>& u_upper,
		std::vector<int>& v_lower, std::vector<int>& v_upper) :
		u_lower(u_lower), u_upper(u_upper), v_lower(v_lower), v_upper(v_upper), nrow(
				u_lower.size()), ncol(v_lower.size()) {

	// check length of row and colum sums
	if (u_lower.size() != u_upper.size() || v_lower.size() != v_upper.size()) {
		std::cerr
				<< "marathon::LindaChain::Error: length of degree sequences does not match!"
				<< std::endl;
	}

}

LindaChain::~LindaChain() {

}

State* LindaChain::computeArbitraryState() const {

	// Create Empty 0-1-matrix of size m times n.
	BinaryMatrix* g = new BinaryMatrix(nrow, ncol);

	// Todo: Set entries of g to zero or one by solving the max flow problem.
	// use g->set_edge(i, j, true) to set entries of g.

	return g;
}

/**
 * Modify the state x by application of t random transitions.
 */
void LindaChain::randomize(State* x, const uint32_t t) const {

	// interpret x as BinaryMatrix object
	BinaryMatrix* s = (BinaryMatrix*) x;

	// do t steps
	for (int x = 0; x < t; x++) {

		// Todo: Transitionsregeln implementieren

		/**
		 *	Beispielfunktionen die du benutzen kannst:
		 *
		 *	// ganze Zahl gleichverteilte im Bereich [0..nrow-1] erzeugen
		 *  int x = ::marathon::random::nextInt(nrow);
		 *
		 *  // Gleitkommazahl gleichverteilt im Bereich [0..1) erzeugen
		 *  double y = ::marathon::random::nextDouble();
		 *
		 *  // Zufällige Auswahl von k ganzen Zahlen aus dem Bereich [0..nrow-1] (ohne Wiederholung):
		 *  int auswahl[k];
		 *  ::marathon::random::select(nrow, k, auswahl);
		 *
		 *  // Zufälliges Permutieren eines Arrays:
		 *  int array[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		 *  ::marathon::random::shuffle<int>(array, 9);
		 *
		 *
		 *  // die einträge der Matrix s können mit
		 *  // s->flip(i,j)  bzw. s->set(i,v, b) geändert werden
		 */

	}
}

}
}
}
