/*
 * MarkovChainDense.h
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#ifndef INCLUDE_MARATHON_MARKOVCHAINDENSE_H_
#define INCLUDE_MARATHON_MARKOVCHAINDENSE_H_

#include <climits>
#include <string>
#include <list>
#include <set>

#include "Rational.h"
#include "State.h"
#include "StateGraph.h"

namespace marathon {

class StateGraph;

/**
 * Virtual MarkovChain Base Class.
 */
class MarkovChainDense {

	friend class StateGraph;

protected:

	std::string instance;		// the instance representation
	bool sparse;			// whether the construction is to be sparse or dense

	/* Variables used for State Graph Construction */
	std::vector<const State*> states;

	/* Private Helper Functions */

	/**
	 * Apply metropolis rule to transition(i,j).
	 */
	rational metropolis(const StateGraph* sg, const int i, const int j);

	int expandDense(StateGraph* sg, const int maxStates, const bool verbose);

	/* Virtual Function that has to be implemented by subclasses */

	/**
	 *  Computes weights for each state.
	 *  @param states The Vector of states.
	 *  @param weights The Vector of weights. After calling the method, this vector must
	 *  have the same size as states and is filled with rationals.
	 */
	virtual void computeWeights(const std::vector<const State*>& states,
			std::vector<rational>& weights);

	/**
	 * Enumerate the state space up a limit.
	 * @param limit The maximal number of states.
	 */
	virtual void enumerateStateSpace(const int limit,
			const bool verbose = false);

	/**
	 * Return the proposal probability P(u,v) for going from state u to state v.
	 */
	virtual rational proposalProbability(const State* u, const State* v) const;

public:

	/**
	 * Create A Markov Chain Object for the input s.
	 * @param s input string of the markov chain.
	 * @param sparse true, if state graph is to be constructed by graph enumeration,
	 * false if is to be constructed by Omega^2 procedure.
	 */
	MarkovChainDense(const std::string& s, bool sparse = true);

	/*
	 * Standard Destructor
	 */
	virtual ~MarkovChainDense();

	/**
	 * @return A reference to the string instance.
	 */
	const std::string& getInstance() const;

	/**
	 * Expands an existing state graph to a maximum of maxStates states.
	 * @param maxStates The maximal number of states after the expansion
	 * @param Enables or disables additional debug output
	 * @return the number of states that has been added during the expansion
	 */
	int expandStateGraph(StateGraph* sg, int maxStates = INT_MAX, bool verbose =
			false);

	/**
	 * Method Stub for computing a canonical path in State Graph sg connecting
	 * state i to state j. The path is represented by a list of indices of the
	 * states and are returned in the path object.
	 */
	virtual void constructPath(const StateGraph* sg, int i, int j,
			std::list<int>& path) const;

	// variables needed for state space enumeration
	struct StackItem {
		int i, j, total;
		const int nrow, ncol;
		bool* bits;
		int* rowsum;
		int* colsum;

		StackItem(std::vector<int>& u, std::vector<int>& v, int sum);
		StackItem(const StackItem& s);
		~StackItem();
		std::string to_string() const;
		friend inline std::ostream& operator<<(std::ostream& out,
				const StackItem& s) {
			out << s.to_string();
			return out;
		}
	};
	std::stack<StackItem> callstack;
	std::stack<StackItem>* local_stacks;
	const int numLocalStacks = 64;
	void expandStack(std::stack<StackItem>& stack, const int limit,
			const bool verbose);

};

}

#endif /* INCLUDE_MARATHON_MARKOVCHAINDENSE_H_ */
