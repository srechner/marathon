#include "../../include/marathon/state_graph.h"

namespace marathon {

StateGraph::StateGraph() :
		numStates(0) {

}

StateGraph::~StateGraph() {

}

size_t StateGraph::getNumStates() const {
	return numStates;
}

size_t StateGraph::getNumArcs() const {
	return arcs.size();
}

Rational StateGraph::getStationary(int i) const {
	return stationary_distribution[i];
}

Rational StateGraph::getMinimalStationary() const {
	return *std::min_element(stationary_distribution.begin(),
			stationary_distribution.end());
}

Rational StateGraph::getTransitionProbability(int u, int v) const {
	Rational r = 0;
	for (size_t i = outgoing_arcs[u]; i < outgoing_arcs[u + 1]; i++) {
		Transition t = arcs[i];
		assert(t.u == u);
		if (t.v == v) {
			r = t.p;
			break;
		}
	}
	return r;
}

Rational StateGraph::getWeight(int i) const {
	return Rational(1);
}

void StateGraph::canonicalPath(int u, int v, std::list<int>& path) const {

}

}
