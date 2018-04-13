/*
 * Created on: Aug 4, 2015
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2016, Steffen Rechner
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef STATEGRAPH_H
#define STATEGRAPH_H

#include "state.h"
#include "markov_chain.h"

#include <climits>
#include <list>
#include <string>
#include <unordered_map>
#include <map>
#include "transition.h"

namespace marathon {

    /**
     * State Graph representation. A State Graph is
     * a directed, weighted graph that represents a instance
     * of a Markov Chain for a certain input instance.
     * This class implements a sparse state graph.
     */
    class StateGraph {

    protected:

        const MarkovChain &mc;                                                  // pointer to markov chain object

        /* Vertices and its attributes */
        std::vector<std::unique_ptr<State>> states;                             // The set of states.
        std::unordered_map<State *, size_t, State::Hash, State::Equal> indices;    // State -> Index

        /* The transition set and views on it. */
        std::vector<Transition *> arcs;
        std::vector<std::vector<Transition *>> outArcs;
        std::vector<std::vector<Transition *>> inArcs;

        /* Variables used for State Graph Construction */
        std::set<size_t> reexpand;

        /**
         * Self-Check: Verifiy that chain is reversible.
         */
        bool reversibilityCheck() const {

            const Rational Z = getNormalizingConstant();

            for (const Transition *t : getArcs()) {

                if (t->from >= t->to)
                    continue;

                const Rational stat_u = getWeight(t->from) / Z;
                const Rational stat_v = getWeight(t->to) / Z;
                const Rational puv = t->weight;

                // find anti parallel transition arc
                for (const Transition *t2 : getInArcs(t->from)) {

                    // t2 is anti parallel to t
                    if (t2->to == t->to) {

                        const Rational pvu = t2->weight;

                        if (stat_u * puv != stat_v * pvu) {
                            std::cerr << "Error! Chain is not reversible!" << std::endl;
                            std::cerr << "P(" << t->from << "," << t->to << ")=" << puv
                                      << std::endl;
                            std::cerr << "P(" << t->to << "," << t->from << ")=" << pvu
                                      << std::endl;
                            std::cerr << "stat(" << t->from << ")=" << stat_u << std::endl;
                            std::cerr << "stat(" << t->to << ")=" << stat_v << std::endl;
                            std::cerr << "stat(" << t->from << ")*P(" << t->from << ","
                                      << t->to << ")=" << stat_u * puv << std::endl;
                            std::cerr << "stat(" << t->to << ")*P(" << t->to << ","
                                      << t->from << ")=" << stat_v * pvu << std::endl;

                            return false;

                        }
                    }
                }
            }

            return true;
        }

        /**
         * Print information about the state graph.
         */
        void printInformation() const {

            std::cout << "state size: " << getNumStates() << std::endl;
            for (int ii = 0; ii < getNumStates(); ii++) {
                const State &s = getState(ii);
                std::cout << ii << ": " << s.toString() << " " << getWeight(ii)
                          << std::endl;
            }

            std::cout << "transition size: " << getNumTransitions() << std::endl;
            for (const Transition *t : getArcs()) {
                std::cout << t->from << " " << t->to << " " << t->weight << std::endl;
            }
        }


        /**
         * Helper Functions
         */

        /**
         * This is a private method that is called during state graph expansion.
         * It computes all neighbouring states of state s and insert them into the state graph repectively
         * into the leftover structures that store the states and arcs for next expandStateGraph().
         * @param i The index of the state that is to be expanded.
         * @param limit The maximal number of states.
         * @param lastStop The size of the state graph when this expansion has been triggered.
         * @param verbose If true, additional debug information is printed.
         * @param True, if all adjacent states could be inserted in the state graph.
         */
        void expandState(const int i, const int limit, const int lastStop,
                         const bool verbose) {
            if (verbose)
                std::cout << "expand state " << i << std::endl;

            const State &s = getState(i);

            // sum of proposal probabilites of adjacent states
            Rational sum(0);

            // for each state that is adjacent to s
            mc.adjacentStates(s, [this, i, limit, lastStop, verbose, &s, &sum]
                    (const State &s2, const marathon::Rational &p) {

                if (p == 0)
                    return;

                // sum proposal propability
                sum += p;

                // Determine the index of s2
                int j = indexOf(s2);    // j is the index of state s2

                // if s2 is already known
                if (j != -1) {

                    if (verbose) {
                        std::cout << " " << j << ": " << s2.toString() << " " << p
                                  << " already known state" << std::endl;
                    }

                    // if arc (i,j) has not already included last round
                    if (j >= lastStop) {

                        // increase transition probability
                        Transition *t = getArc(i, j);
                        if (t == nullptr) {
                            addArc(i, j, p);
                        } else {
                            t->weight += p;
                        }
                    }
                }
                    // if s2 is not seen so far and the maximal number of states is not reached
                else if (getNumStates() < limit) {

                    // add s2 to the vector of states
                    j = addState(s2);

                    if (verbose) {
                        std::cout << " " << j << ": " << s2.toString() << " " << p
                                  << " new state" << std::endl;
                    }

                    // add a transition arc to the state graph
                    addArc(i, j, p);
                }
                    // if s2 is not seen so far and the maximal number of states is reached
                else {
                    reexpand.insert(i);
                }
            });

            // do the proposal probabilites sum up to one?
            if (sum != Rational(1)) {
                std::stringstream ss;
                ss << "Expection in marathon::StateGraph::expandState: Sum of transition probabilities of state "
                   << i << "(" << s << ") is " << sum << "!= 1!";
                throw std::runtime_error(ss.str());
            }
        }

    public:

        /**
         * Standard Constructor. Creates an empty State Graph.
         * @param mc A pointer to the Markov Chain Object that defines transition rules, etc.
         * @param limit A limit on the number of states of the graph. The graph can later on be expanded by the expand() method.
         * @param verbose Enable debug output during state graph construction.
         */
        StateGraph(
                const MarkovChain &mc,
                const int limit = INT_MAX,
                const bool verbose = false
        ) : mc(mc) {
            expand(limit, verbose);
        }

        /**
         * Standard Destructor. Remove arcs and states.
         */
        virtual ~StateGraph() {
            // delete all transitions
            for (int i = 0; i < arcs.size(); i++)
                delete arcs[i];
        }

        /**
         * @return Returns a reference to the corresponding Markov Chain Object.
         */
        const MarkovChain &getMarkovChain() const {
            return mc;
        }

        /**
         * Adds a transition arc to the graph.
         * Precondition: The state graph does not already contain an arc between state u and state v.
         * @return Returns the index of the new transition.
         */
        int addArc(const int u, const int v, const Rational &p) {
            return addArc(new Transition(u, v, p));
        }

        /**
         * Adds a transition arc to the graph.
         * Precondition: The state graph does not already contain an arc between state t.u and state t.v.
         * @return Returns the index of the new transition.
         */
        int addArc(Transition *t) {
            // add the arc to the arc vector
            arcs.push_back(t);

            // add a pointer to the transition in u's and v's outarc/inarc array
            outArcs[t->from].push_back(t);
            inArcs[t->to].push_back(t);

            return arcs.size() - 1;
        }

        /**
         * @return Return a pointer to the arc that connects u with v or nullptr, if no such arc exists.
         */
        Transition *getArc(int u, int v) const {
            // search for transition (u,v)
            for (Transition *t : outArcs[u]) {
                assert(t->from == u);
                if (t->to == v) {
                    return t;
                }
            }
            return nullptr;
        }

        /**
         * Returns the number of states of the state graph
         */
        size_t getNumStates() const {
            return states.size();
        }

        /**
         * Returns the number of Transitions/Arcs of the state graph
         */
        size_t getNumTransitions() const {
            return arcs.size();
        }

        /**
         * Returns the transition probability P_uv for going from states[u] to states[v]
         */
        Rational getTransitionProbability(int u, int v) const {
            const Transition *t = getArc(u, v);
            if (t == nullptr)
                return 0;
            else
                return t->weight;
        }

        /**
         * Set P(u,v) to p
         */
        void setTransitionProbability(int u, int v, Rational p) {
            Transition *t = getArc(u, v);

            if (t != nullptr) {
                t->weight = p;
            } else {
                // no transition found? add a new one
                addArc(new Transition(u, v, p));
            }
        }

        /**
         * Increases P(u,v) by an amount of p.
         */
        void addTransitionProbability(int u, int v, Rational p) {
            Transition *t = getArc(u, v);

            if (t != nullptr) {
                t->weight += p;
            } else {
                // no transition found? add a new one
                addArc(new Transition(u, v, p));
            }
        }


        /**
         * Return the weight of state i.
         */
        Rational getWeight(const int i) const {
            return mc.getWeight(*(states[i]));
        }

        /**
         * Return the minimal weight of a state.
         */
        Rational getMinWeight() const {
            Rational min_weight = getWeight(0);
            for (int i = 1; i < states.size(); i++) {
                Rational wi = getWeight(i);
                if (wi < min_weight)
                    min_weight = wi;
            }
            return min_weight;
        }

        /**
         * Return the sum of all weights.
         */
        Rational getNormalizingConstant() const {
            Rational Z = 0;
            for (const auto &s : states)
                Z += mc.getWeight(*s);
            return Z;
        }


        /**
         * Returns a reference to the outgoing arcs of state v.
         */
        const std::vector<Transition *> &getOutArcs(int v) const {
            return outArcs[v];
        }

        /**
         * Returns a reference to the ingoing arcs of state v.
         */
        const std::vector<Transition *> &getInArcs(int v) const {
            return inArcs[v];
        }

        /**
         * Returns a reference to the vector of all arcs in the state graph.
         */
        const std::vector<Transition *> &getArcs() const {
            return arcs;
        }

        /**
         * Return a pointer to arc with index i.
         * @param i The index of the arc.
         * @return A pointer to the i'th transition.
         */
        Transition *getArc(const int i) const {
            return (Transition *) &arcs[i];
        }

        /**
         * Returns the number of adjacent states of state[v]
         */
        int getNumOutArcs(int v) const {
            return this->getOutArcs(v).size();
        }

        /**
         * Removes all States and Transitions and re-initializes the state graph.
         */
        virtual void clear() {
            arcs.clear();
            inArcs.clear();
            outArcs.clear();
        }

        /**
         * Add a new State to the state graph.
         * @param s The State to insert.
         * @return The index of the state after insertion.
         */
        size_t addState(const State &s) {

            // create a copy
            std::unique_ptr<State> c = s.copy();

            const size_t index = states.size();

            // add state to the vector of states
            indices[c.get()] = index;
            //indices.insert(std::make_pair(c.get(), index));
            states.push_back(std::move(c));

            // prepare vector of transitions
            outArcs.push_back(std::vector<Transition *>());
            inArcs.push_back(std::vector<Transition *>());
            return states.size() - 1;
        }

        /**
         * Returns a reference to the State with index i.
         */
        const State &getState(int i) const {
            return *(states[i]);
        }

        /**
         * Returns a reference to a vector of States.
         */
        const std::vector<std::unique_ptr<State>> &getStates() const {
            return states;
        }

        /**
         * Returns the index of a state, or SIZE_MAX if the state graph does not contain this state.
         */
        size_t indexOf(const State &s) const {
            auto it = indices.find(const_cast<State*>(&s));
            if (it != indices.end())
                return it->second;
            else
                return SIZE_MAX;
        }

        /**
     * Expands an existing state graph to a given maximum of states.
     * @param limit The maximal number of states after the expansion
     * @param verbose Enables or disables additional debug output
     * @return the number of states that has been added during the expansion
     */
        virtual
        void expand(const size_t limit = SIZE_MAX, const bool verbose = false) {

            const size_t sizeLast = getNumStates();

            // nothing to do?
            if (limit <= sizeLast)
                return;

            // if state graph is empty
            if (sizeLast == 0) {

                // Start with arbitrary State
                const State &s1 = mc.getCurrentState();

                if (verbose) {
                    std::cout << "Start state is " << s1.toString() << std::endl;
                }

                // add initial state
                addState(s1);
            }

            /***********************************************************
             * Expand states that have not been expanded completely
             **********************************************************/

            if (verbose)
                std::cout << "re-expand states from last round" << std::endl;

            // gather the indices that have to re-expanded
            std::vector<int> tmp(reexpand.begin(), reexpand.end());
            reexpand.clear();

            // for all states that are not yet completely expanded
            for (int y : tmp) {
                // re-expand state with index y
                expandState(y, limit, sizeLast, verbose);
            }

            /***********************************************************
             * Explore further regions of the state graph.
             **********************************************************/

            if (verbose)
                std::cout << "explore further regions" << std::endl;

            for (int y = sizeLast; y < getNumStates(); y++) {
                // expand state with index y
                expandState(y, limit, 0, verbose);
            }

            /***********************************************************
             * Self-Check: Verifiy that chain is reversible
             ***********************************************************/
            assert(reversibilityCheck());

            /*********************************************************************
             * Print Information about states
             ********************************************************************/
            if (verbose)
                printInformation();
        }

    };

    inline
    std::ostream &operator<<(std::ostream &out, const StateGraph &sg) {

        out << "n " << sg.getNumStates() << " m " << sg.getNumTransitions() << "\n";
        for (int i = 0; i < sg.getNumStates(); i++) {
            for (Transition *t : sg.getOutArcs(i))
                out << t->from << " " << t->to << " " << t->weight << "\n";
        }

        return out;
    }

}

#endif /* STATEGRAPH_H */
