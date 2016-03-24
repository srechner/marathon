/*
 * MarkovChainDense.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

/*
 * markov_chain.cpp
 *
 *  Created on: Feb 15, 2016
 *      Author: steffen
 */

#include "../../include/marathon/MarkovChainDense.h"

marathon::MarkovChainDense::MarkovChainDense(const std::string& s, bool sparse) :
		instance(s), sparse(sparse) {

}

marathon::MarkovChainDense::~MarkovChainDense() {

}

const std::string& marathon::MarkovChainDense::getInstance() const {
	return this->instance;
}

void marathon::MarkovChainDense::constructPath(const StateGraph* sg, int i,
		int j, std::list<int>& path) const {

}

void marathon::MarkovChainDense::computeWeights(
		const std::vector<const State*>& states,
		std::vector<rational>& weights) {
	weights.clear();
	weights.reserve(states.size());
	for (int i = 0; i < states.size(); i++)
		weights.push_back(1);
}

int marathon::MarkovChainDense::expandStateGraph(StateGraph* sg,
		const int maxStates, const bool verbose) {
	return expandDense(sg, maxStates, verbose);
}

int marathon::MarkovChainDense::expandDense(StateGraph* sg, const int maxStates,
		const bool verbose) {

	const int sizeLast = sg->getNumStates();

	// further explore the state space

	enumerateStateSpace(maxStates, verbose);

	for (int i = sizeLast; i < states.size(); i++)
		sg->addState(states[i]);

	const int sizeNew = sg->getNumStates();

	// compute new weights
	computeWeights(sg->states, sg->weights);

	/*************************************************************
	 * Compute Transition probability between old and new states.
	 ************************************************************/

	// for each old state
#pragma omp parallel for if(!verbose)
	for (int i = 0; i < sizeLast; i++) {

		// sum of all transition probabilities of state i
		rational sum_i = 0;

		// for each new state
		for (int j = sizeLast; j < sizeNew; j++) {

			// apply metropolis rule
			const rational metr = metropolis(sg, i, j);

			// add transition arc
			if (metr != 0) {
#pragma omp critical
				sg->addArc(i, j, metr);
				sum_i += metr;
			}
		}

		// reduce loop probability of state i by the new assigned probability
		Transition* loop_i = sg->getLoop(i);
		loop_i->p -= sum_i;
	}

	/*************************************************************
	 * Compute Transition probability between new states and all states.
	 ************************************************************/

	// for each new state
#pragma omp parallel for if(!verbose)
	for (int i = sizeLast; i < sizeNew; i++) {

		// sum of all transition probabilities of state i
		rational sum_i = 0;

		// for each state (old and new)
		for (int j = 0; j < sizeNew; j++) {

			// if non-loop probability
			if (i != j) {

				// apply metropolis rule
				const rational metr = metropolis(sg, i, j);

				// add transition arc
				if (metr != 0) {
#pragma omp critical
					sg->addArc(i, j, metr);
					sum_i += metr;
				}
			}
		}

		// set loop probability
#pragma omp critical
		sg->addArc(i, i, rational(1) - sum_i);
	}

	return sizeNew - sizeLast;
}

marathon::rational marathon::MarkovChainDense::metropolis(const StateGraph* sg,
		const int i, const int j) {
	const rational wi = sg->getWeight(i);
	const rational wj = sg->getWeight(j);
	rational metr = 1;
	if (wj < wi)
		metr = wj / wi;
	metr *= proposalProbability(sg->getState(i), sg->getState(j));
	return metr;
}


SwitchBipartite::StackItem::StackItem(std::vector<int>& u, std::vector<int>& v,
		int sum) :
		i(0), j(0), total(sum), nrow(u.size()), ncol(v.size()) {
	bits = new bool[nrow * ncol];
	rowsum = new int[nrow];
	colsum = new int[ncol];
	memset(bits, 0, nrow * ncol * sizeof(bool));
	for (int i = 0; i < nrow; i++)
		rowsum[i] = u[i];
	for (int j = 0; j < ncol; j++)
		colsum[j] = v[j];
}

SwitchBipartite::StackItem::StackItem(const StackItem& s) :
		i(s.i), j(s.j), total(s.total), nrow(s.nrow), ncol(s.ncol) {
	bits = new bool[nrow * ncol];
	rowsum = new int[nrow];
	colsum = new int[ncol];
	memcpy(bits, s.bits, nrow * ncol * sizeof(bool));
	memcpy(rowsum, s.rowsum, nrow * sizeof(int));
	memcpy(colsum, s.colsum, ncol * sizeof(int));
}

SwitchBipartite::StackItem::~StackItem() {
	delete[] bits;
	delete[] rowsum;
	delete[] colsum;
}

std::string SwitchBipartite::StackItem::to_string() const {

	std::stringstream ss;
	const int p = i * ncol + j;

	for (int k = 0; k < p; k++)
		ss << " ";
	ss << "i=" << i << ", j=" << j << ", total=" << total << ", rowsum=";
	for (int k = 0; k < nrow; k++)
		ss << rowsum[k] << ",";
	ss << " colsum=";
	for (int k = 0; k < ncol; k++)
		ss << colsum[k] << ",";
	ss << " bits=";
	for (int k = 0; k < ncol * nrow; k++)
		ss << (int) bits[k];

	return ss.str();
}

void SwitchBipartite::expandStack(std::stack<StackItem>& stack, const int limit,
		const bool verbose) {

	// pop item from stack
	const StackItem s = stack.top();
	stack.pop();

	const int p = s.i * s.ncol + s.j;

	if (verbose)
		std::cout << s << std::flush;

	// found a solution?
	if (s.total == 0) {
#pragma omp critical
		{
			if (states.size() < limit) {
				states.push_back(
						new DenseBipartiteGraph(s.nrow, s.ncol, s.bits));
				if (verbose)
					printf("*\n");
			} else {
				if (verbose)
					printf("#\n");
				stack.push(s);
			}
		}
	}
	// not a solution? keep backtracking!
	else {

		if (verbose) {
			printf("\n");
		}

		// make a copy of the stack item
		StackItem t(s);

		// forward checking: if there are enough columns for the remaining ones
		if (t.ncol - t.j > t.rowsum[t.i]) {

			// prepare variables for recursive call
			t.bits[p] = 0;	// bits[i,j]=0

			if (s.i < s.nrow && s.j == s.ncol - 1) {
				t.i = s.i + 1;
				t.j = 0;
				stack.push(t);
			} else if (s.i < s.nrow && s.j < s.ncol - 1) {
				t.j = s.j + 1;
				stack.push(t);
			}
		}

		// forward checking: column j has enough ones left
		if (s.rowsum[s.i] > 0 && s.colsum[s.j] > 0) {

			// prepare variable for recursive call
			t.bits[p] = 1;
			t.rowsum[s.i] = s.rowsum[s.i] - 1;
			t.colsum[s.j] = s.colsum[s.j] - 1;
			t.total = s.total - 1;

			// forward checking: skip non-realizable sequence
			if (isRealizable(t.rowsum, t.colsum, t.nrow, t.ncol)) {
				if (s.i < s.nrow && s.j == s.ncol - 1) {
					t.i = s.i + 1;
					t.j = 0;
					stack.push(t);
				} else if (s.i < s.nrow && s.j < s.ncol - 1) {
					t.j = s.j + 1;
					stack.push(t);
				}
			}
		}
	}
}

void SwitchBipartite::enumerateStateSpace(const int limit, const bool verbose) {

	// first call of enumerateStateSpace?
	if (states.empty()) {

		// init call stack
		callstack.push(StackItem(u, v, sum));

		// construct start items for parallel search
		while (callstack.size() < numLocalStacks && states.size() < limit
				&& !callstack.empty()) {

			// keep stacking items
			expandStack(callstack, limit, verbose);
		}

		// transfer stack items into local stacks
		for (int i = 0; !callstack.empty(); i++, callstack.pop()) {
			const StackItem& s = callstack.top();
			local_stacks[i].push(StackItem(s));
		}
	}

	// for each local call stack do in parallel
#pragma omp parallel for if(!verbose)
	for (int k = 0; k < numLocalStacks; k++) {
		// work local stack
		while (states.size() < limit && !local_stacks[k].empty()) {
			expandStack(local_stacks[k], limit, verbose);
		}
	}
}

rational SwitchBipartite::proposalProbability(const State* u,
		const State* v) const {

	// convert pointers
	const DenseBipartiteGraph* _u = (const DenseBipartiteGraph*) u;
	const DenseBipartiteGraph* _v = (const DenseBipartiteGraph*) v;

	const long m = _u->get_nrows();
	const long n = _u->get_ncols();

	// determine symmetric difference of bipartite graphs
	boost::dynamic_bitset<> symdif = _u->M ^ _v->M;

	// check whether symdif is a cycle of length four
	if (symdif.count() == 4) {
		return rational(4, m * (m + 1) * n * (n + 1));
	} else
		return 0;
}

