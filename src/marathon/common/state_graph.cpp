#include "../../../include/marathon/common/state_graph.h"

#include <string>

namespace marathon {

StateGraph::StateGraph(int n) : numStates(0),
		instance("") {

	inArcs = new std::vector<Transition*>[n];
	outArcs = new std::vector<Transition*>[n];
	stationary_distribution.reserve(n);

}

StateGraph::~StateGraph() {

	if (inArcs != nullptr)
		delete[] inArcs;

	if (outArcs != nullptr)
		delete[] outArcs;
}

size_t StateGraph::getNumStates() const {
		return numStates;
}

void StateGraph::addArc(int u, int v, rational p) {
	Transition* t = new Transition(u, v, p);
	arcs.push_back(t);
	outArcs[u].push_back(t);
	inArcs[v].push_back(t);
}

const std::string& StateGraph::getInstance() const {
	return instance;
}

void StateGraph::setInstance(const std::string& instance) {
	this->instance = instance;
}

size_t StateGraph::getNumTransitions() const {
	return arcs.size();
}

rational StateGraph::getTransitionProbability(int u, int v) const {
	// search for transition (u,v)
	rational r = 0;
	for (Transition* t : outArcs[u]) {
		assert(t->u == u);
		if (t->v == v) {
			return t->p;
		}
	}
	return 0;
}

void StateGraph::setTransitionProbability(int u, int v, rational p) {
	// search for transition (u,v)
	for (Transition* t : outArcs[u]) {
		if (t->v == v) {
			t->p = p;
			break;
		}
	}
}

rational StateGraph::getStationary(const int i) const {
	return stationary_distribution[i];
}

void StateGraph::setStationary(const int i, const rational p) {
	stationary_distribution[i] = p;
}

rational StateGraph::getMinimalStationary() const {
	return *std::min_element(stationary_distribution.begin(),
			stationary_distribution.end());
}

const std::vector<Transition*>& StateGraph::getOutArcs(int v) const {
	return outArcs[v];
}

const std::vector<Transition*>& StateGraph::getInArcs(int v) const {
	return inArcs[v];
}

const std::vector<Transition*>& StateGraph::getArcs() const {
	return arcs;
}

int StateGraph::getNumOutArcs(int v) const {
	return this->getOutArcs(v).size();
}
void StateGraph::clear() {

	arcs.clear();

	for (size_t i = 0; i < arcs.size(); i++)
		delete arcs[i];

	if (inArcs != nullptr)
		delete[] inArcs;
	if (outArcs != nullptr)
		delete[] outArcs;
	stationary_distribution.clear();
}

std::ostream& operator<<(std::ostream& out, const StateGraph& sg) {

	out << "n " << sg.getNumStates() << " m " << sg.getNumTransitions() << "\n";
	for (int i = 0; i < sg.getNumStates(); i++) {
		for (Transition* t : sg.getOutArcs(i))
			out << t->u << " " << t->v << " " << t->p << "\n";
	}

	return out;
}

}
