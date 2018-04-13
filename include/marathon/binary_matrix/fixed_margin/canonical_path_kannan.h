/*
 * KannanCanPath.h
 *
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

#ifndef INCLUDE_MARATHON_CHAINS_SEQUENCES_KANNANCANPATH_H_
#define INCLUDE_MARATHON_CHAINS_SEQUENCES_KANNANCANPATH_H_

#include "marathon/path_construction_scheme.h"
#include "marathon/binary_matrix/binary_matrix.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Path construction scheme presented by
             *
             *   R. Kannan, P. Tetali, and S. Vempala.
             *   Simple Markov-chain algorithms for generating bipartite graphs and tournaments.
             *   Random Structures & Algorithms 14 (1997), 293–308.
             *   doi: 10. 1002 /(sici )1098 - 2418(199907) 14:4<293::aid-rsa1>3.0.co;2-g.
             */
            class KannanPath : public PathConstructionScheme {

            protected:

                int next_red_edge(int col, bool *red_edges, int m, int n) const {
                    for (int i = 0; i < m; i++) {
                        if (red_edges[i * n + col])
                            return i;
                    }
                    // no edge found in column
                    return -1;
                }

                int next_blue_edge(int row, bool *blue_edges, int m, int n) const {
                    for (int j = 0; j < n; j++) {
                        if (blue_edges[row * n + j])
                            return j;
                    }
                    // no blue edge found in row
                    return -1;
                }

                void trace_cycle(bool *blue_edges, bool *red_edges, int m, int n,
                                 int i, int j, std::vector<int> &cycle) const {

                    while (j != -1) {

                        // add (i,j) to cycle
                        cycle.push_back(i);
                        cycle.push_back(m + j);

                        // remove blue edge (i,j)
                        blue_edges[i * n + j] = false;

                        i = next_red_edge(j, red_edges, m, n);

                        // remove red edge (i,j)
                        red_edges[i * n + j] = false;

                        j = next_blue_edge(i, blue_edges, m, n);
                    };
                }

                void splice_cycle(std::vector<int> cycle,
                                  std::list<std::vector<int>> &cycles, const int m, const int n) const {

#ifdef DEBUG
                    std::cout << "splice cycle [ ";
                for (std::vector<int>::iterator it = cycle.begin(); it != cycle.end();
                        ++it) {
                    std::cout << *it << " ";
                }
                std::cout << "]" << std::endl;
#endif

                    std::vector<int> c;
                    bool removed[cycle.size()];
                    int first_seen[n + m];
                    int i, j;

                    memset(removed, 0, cycle.size() * sizeof(bool));

                    for (i = 0; i < n + m; i++)
                        first_seen[i] = n + m;

                    for (i = 0; i < cycle.size(); i++) {

#ifdef DEBUG
                        std::cout << i << ": " << cycle[i] << std::endl;
#endif

                        // smaller cycle detected
                        if (first_seen[cycle[i]] != n + m) {

#ifdef DEBUG
                            std::cout << "smaller cycle detected" << std::endl;
#endif

                            // extract smaller cycle and store in cycle list
                            c.clear();
                            for (j = first_seen[cycle[i]]; j != i; j++) {
                                if (!removed[j]) {
                                    c.push_back(cycle[j]);
                                    first_seen[cycle[j]] = n + m;
                                    removed[j] = true;
                                }
                            }
                            cycles.push_back(c);
                        }
                        first_seen[cycle[i]] = i;
                    }

                    // not removed vertices
                    c.clear();
                    for (i = 0; i < cycle.size(); i++) {
                        if (!removed[i])
                            c.push_back(cycle[i]);
                    }
                    cycles.push_back(c);
                }

                void cycle_decomposition(
                        const BinaryMatrix &x,
                        const BinaryMatrix &y,
                        std::list<std::vector<int>> &cycles) const {

                    const int nrow = x.getNumRows();
                    const int ncol = x.getNumCols();

                    bool red[nrow * ncol];
                    bool blue[nrow * ncol];

                    std::vector<int> cycle;
                    std::vector<int>::iterator cit;
                    std::list<std::vector<int> >::iterator
                            it;

                    memset(red, 0, nrow * ncol * sizeof(bool));
                    memset(blue, 0, nrow * ncol * sizeof(bool));

                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {
                            if (x.get(i, j) && !y.get(i, j))
                                blue[i * ncol + j] = true;
                            else if (!x.get(i, j) && y.get(i, j))
                                red[i * ncol + j] = true;
                        }
                    }

                    for (int i = 0; i < nrow; i++) {
                        for (int j = 0; j < ncol; j++) {
                            if (blue[i * ncol + j]) {
                                // start of alternating Cycle in x found
                                cycle.clear();
                                // trace cycle
                                trace_cycle(blue, red, nrow, ncol, i, j, cycle);
                                // try to splice cycles into smaller ones
                                splice_cycle(cycle, cycles, nrow, ncol);
                            }
                        }
                    }
                }

            public:

                virtual std::list<int> construct(const StateGraph &sg, const int s, const int t) const override {

#ifdef DEBUG
                    std::cout << "from=" << sg->getState(s) << std::endl;
    std::cout << "to  =" << sg->getState(t) << std::endl;
#endif

                    std::list<int> path;
                    path.push_back(s);  // start with from

                    // continously modify u
                    BinaryMatrix u(static_cast<const BinaryMatrix &>(sg.getState(s)));
                    BinaryMatrix v(static_cast<const BinaryMatrix &>(sg.getState(t)));

                    std::list<std::vector<int>> cycles;
                    std::list<std::vector<int> >::iterator it;

                    int i, l;
                    // decompose symmetric difference of u and v into cycles
                    cycle_decomposition(u, v, cycles);

#ifdef DEBUG
                    std::cout << "found " << cycles.size() << " cycles" << std::endl;
#endif

                    // sort cycles by length
                    cycles.sort([](const std::vector<int> &c1, const std::vector<int> &c2) -> bool {
                        return c1.size() < c2.size();
                    });

                    while (!cycles.empty()) {

                        // deal with smallest cycle
                        auto cycle = cycles.front();
                        cycles.pop_front();

                        assert(cycle.size() % 2 == 0);
                        assert(cycle.size() >= 4);
#ifdef DEBUG
                        std::cout << "[";
        for (auto it2 = cycle.begin(); it2 != cycle.end(); ++it2) {
            std::cout << " " << *it2;
        }
        std::cout << " ]" << std::endl;
#endif

                        l = cycle.size();

                        /*
                         * find first vertex w[i] s.t. the cycle
                         *     (w[0], w[1], ..., w[i], w[i+1], w[l-i-2], w[l-i-1], ... , w[l-1])
                         * is alternating in u
                         */
                        i = 0;
                        while (true) {

                            int u1 = cycle[i];
                            int v1 = cycle[i + 1];
                            int u2 = cycle[l - i - 2];
                            int v2 = cycle[l - i - 1];

                            /* translate node labels */
                            if (u1 > v1) {
                                // rotate
                                uint tmp = u1;
                                u1 = v1;
                                v1 = u2;
                                u2 = v2;
                                v2 = tmp;
                            }

                            v1 -= u.getNumRows();
                            v2 -= u.getNumRows();

                            if (u.isCheckerBoardUnit(u1, v1, u2, v2))
                                break;

                            i++;
                        }

#ifdef DEBUG
                        std::cout << "switch cycle: [ " << cycle[i] << " " << cycle[i + 1] << " "
        << cycle[l - i - 2] << " " << cycle[l - i - 1] << " ]" << std::endl;
#endif

                        /**
                         * (w[i], w[i+1], w[l-i-2], w[l-i-1]) is switchable
                         *    switch it!
                         */

                        int u1 = cycle[i];
                        int v1 = cycle[i + 1];
                        int u2 = cycle[l - i - 2];
                        int v2 = cycle[l - i - 1];

                        if (u1 > v1) {
                            // rotate
                            uint tmp = u1;
                            u1 = v1;
                            v1 = u2;
                            u2 = v2;
                            v2 = tmp;
                        }

                        // translate node labels
                        v1 -= u.getNumRows();
                        v2 -= u.getNumRows();

                        assert(u.isCheckerBoardUnit(u1, v1, u2, v2));
                        u.flipSubmatrix(u1, v1, u2, v2);
                        int x = sg.indexOf(u);
                        assert(x != -1);
                        path.push_back(x);

                        /**
                         * two smaller alternating cycles are created:
                         *  1: (w[0], ... , w[i], w[l-i-1], ... , w[l-1])
                         *  2: (w[i+1], ... , w[l-i-2])
                         *
                         *  Erase first cycle by series of switches
                         */

                        for (int j = i; j > 0; j--) {

                            int u1 = cycle[j - 1];
                            int v1 = cycle[j];
                            int u2 = cycle[l - j - 1];
                            int v2 = cycle[l - j];

                            if (u1 > v1) {
                                // rotate
                                uint tmp = u1;
                                u1 = v1;
                                v1 = u2;
                                u2 = v2;
                                v2 = tmp;
                            }

                            // translate node labels
                            v1 -= u.getNumRows();
                            v2 -= u.getNumRows();

                            u.flipSubmatrix(u1, v1, u2, v2);
                            int x = sg.indexOf(u);
                            assert(x != -1);
                            path.push_back(x);
                        }

                        /**
                         * We are left with the second cycle which is smaller
                         * Continue with this cycle
                         */

                        // if second cycle has at least 4 edges
                        if (l - 2 * i - 2 >= 4) {
                            cycle.clear();
                            for (int j = i + 1; j <= l - i - 2; j++) {
                                cycle.push_back(cycle[j]);
                            }
                            cycles.push_front(cycle);
                        }
                    }

                    return path;
                }

            };
        }
    }
}

#endif /* INCLUDE_MARATHON_CHAINS_SEQUENCES_KANNANCANPATH_H_ */
