/*
 * Created on: Mar 22, 2016
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

#ifndef INCLUDE_MARATHON_CHAINS_MATCHING_JS89CANPATH_H_
#define INCLUDE_MARATHON_CHAINS_MATCHING_JS89CANPATH_H_

#include "marathon/path_construction_scheme.h"
#include "marathon/markov_chain.h"
#include "perfect_matching_broder86.h"

namespace marathon {
    namespace matching {

        /**
         * Path construction scheme presented by
         *
         * M. Jerrum and A. Sinclair. Approximating the permanent.
         * SIAM Journal on Computing 18 (1989), 1149–1178.
         * doi: 10.1137/0218077.
         */
        class JS89Path : public PathConstructionScheme {

            /**
             * Path construction Scheme for canonical paths from Jerrum and Sinclair 1989.
             * @param sg Pointer to a state graph object.
             * @param s Index of a state (Start).
             * @param t Index of a state (End).
             * @return List of State indices that corresponds to the state.
             */
            virtual std::list<int> construct(const StateGraph &sg, const int s, const int t) const override {

                const Broder86 &mc = static_cast<const Broder86 &>(sg.getMarkovChain());

                std::list<int> path;    // will be returned

                if (s == t)
                    return path;

                // Implementation after JS89
                std::list<int> init_seqment;
                std::list<int> main_seqment;
                std::list<int> final_seqment;

                const size_t omega = sg.getNumStates();

                // working variables
                int i, c, a, sw;
                const size_t n = mc._g.getNumberOfNodes();
                const BipartiteMatching *u, *v;
                std::queue<int> q;

                // arrays to store the labels for BFS
                int *prev = new int[omega];
                bool *visited = new bool[omega];

                //std::unordered_map<BipartiteMatching, BipartiteMatching> prev;
                //BipartiteMatching null;

                int start_states[2] = {s, t};

                std::list<int> *segments[2] = {&init_seqment, &final_seqment};

                // Initial Segment and Final Segment by BFS
                for (i = 0; i < 2; i++) {

                    // clear queue and prevs
                    while (!q.empty())
                        q.pop();

                    // start BFS with s respective t
                    memset(visited, 0, omega * sizeof(bool));
                    q.push(start_states[i]);
                    visited[start_states[i]] = 1;
                    prev[start_states[i]] = -1;
                    while (!q.empty()) {
                        c = q.front();
                        q.pop();
                        const BipartiteMatching &u = static_cast<const BipartiteMatching &>(sg.getState(c));

                        if (2 * u._edges == n) {    // Path to Perfect Matching
                            // Reconstruct Path and return
                            do {
                                segments[i]->push_front(c);
                                c = prev[c];
                            } while (c != -1);
                            break;
                        } else {
                            // for all neighbors of u
                            for (Transition *t : sg.getOutArcs(c)) {
                                // if not visited yet
                                if (!visited[t->to]) {
                                    q.push(t->to);
                                    prev[t->to] = c;
                                    visited[t->to] = 1;
                                }
                            }
                        }
                    }
                }

                final_seqment.reverse();

                // Main Segment
                const int ii = init_seqment.back();
                const int ff = final_seqment.front();
                const BipartiteMatching &I = static_cast<const BipartiteMatching &>(sg.getState(ii));
                const BipartiteMatching &F = static_cast<const BipartiteMatching &>(sg.getState(ff));
                BipartiteMatching s2(I);

                //std::cout << "I = " << I->toString() << std::endl;
                //std::cout << "F = " << F->toString() << std::endl;

                const BipartiteMatching *x[2] = {&I, &F};

                std::vector<size_t> cycle;
                std::vector<size_t>::iterator it2;

                std::vector<bool> unrolled(n);

                // unroll symmetric difference of I and F Cycle for Cycle
                for (i = 0; i < n; i++) {

                    // if some node lies on a cycle which isn't unrolled yet
                    if (I.mates[i] != F.mates[i] && !unrolled[i]) {

                        // Collect Cycle Nodes and store in list
                        cycle.clear();
                        a = i;

                        sw = 0;    // switch variable
                        do {
                            cycle.push_back(a);
                            unrolled[a] = true;
                            a = x[sw]->getMate(a);
                            sw = 1 - sw;
                        } while (a != i);

                        // unroll cycle

                        // first step: remove (u0, v0)

                        int u0 = cycle[0];
                        int v0 = cycle[1];
                        s2.removeEdge(u0, v0);
                        int x = sg.indexOf(s2);
                        assert(x != -1);
                        main_seqment.push_back(x);

                        // replace each edge (u_j, v_j) by (u_j, v_j-1)
                        for (int j = 1; j < cycle.size() / 2; j++) {
                            size_t u_j = cycle[2 * j];
                            size_t v_j = cycle[2 * j + 1];
                            size_t v_jj = cycle[2 * j - 1];
                            s2.mates[u_j] = v_jj;
                            s2.mates[v_jj] = u_j;
                            s2.mates[v_j] = -1;
                            s2.unmatched2 = v_j;
                            x = sg.indexOf(s2);
                            assert(x != -1);
                            main_seqment.push_back(x);
                        }

                        // last step: add (u0, v_last)
                        s2.addEdge(u0, cycle.back());
                        x = sg.indexOf(s2);
                        assert(x != -1);
                        main_seqment.push_back(x);
                    }
                }

                /*std::cout << "init" << std::endl;
                for (std::list<int>::iterator it = init_seqment.begin();
                     it != init_seqment.end(); ++it) {
                    std::cout << *it << " " << sg->getState(*it) << std::endl;
                }

                std::cout << "main" << std::endl;
                for (std::list<int>::iterator it = main_seqment.begin();
                     it != main_seqment.end(); ++it) {
                    std::cout << *it << " " << sg->getState(*it) << std::endl;
                }

                std::cout << "final" << std::endl;
                for (std::list<int>::iterator it = final_seqment.begin();
                     it != final_seqment.end(); ++it) {
                    std::cout << *it << " " << sg->getState(*it) << std::endl;
                }*/

                path.insert(path.end(), init_seqment.begin(), init_seqment.end());
                path.insert(path.end(), main_seqment.begin(), main_seqment.end());
                path.insert(path.end(), ++final_seqment.begin(), final_seqment.end());

                delete[] prev;
                delete[] visited;

                return path;
            }

        };

    }
}

#endif /* INCLUDE_MARATHON_CHAINS_MATCHING_JS89CANPATH_H_ */
