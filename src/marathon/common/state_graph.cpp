#include "../../../include/marathon/common/state_graph.h"

namespace marathon {

StateGraph::StateGraph() :
		numStates(0) {

}

StateGraph::~StateGraph() {

}

size_t StateGraph::getNumStates() const {
	return numStates;
}

size_t StateGraph::getNumTransitions() const {
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

void StateGraph::setTransitionProbability(int u, int v, Rational p) {
	// search for transition (u,v)
	for(int i=getIndexOfFirstTransition(u); i<= getIndexOfLastTransition(u); i++) {
		Transition& t = arcs[i];
		if(t.v == v) {
			t.p = p;
			break;
		}
	}
}

Rational StateGraph::getWeight(int i) const {
	return Rational(1);
}

void StateGraph::canonicalPath(int u, int v, std::list<int>& path) const {

}

int StateGraph::getIndexOfFirstTransition(int v) const {
	return outgoing_arcs[v];
}

int StateGraph::getIndexOfLastTransition(int v) const {
	return outgoing_arcs[v + 1] - 1;
}

Transition StateGraph::getTransition(int index) const {
	return arcs[index];
}

int StateGraph::getNumTransitions(int v) const {
	if (outgoing_arcs[v] == arcs.size())
		return 0;
	return outgoing_arcs[v + 1] - outgoing_arcs[v];
}

}
