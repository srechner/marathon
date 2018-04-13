/*
 * PathCongestion.h
 *
 * Created on: Mar 24, 2016
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

#ifndef INCLUDE_MARATHON_PATHCONGESTION_H_
#define INCLUDE_MARATHON_PATHCONGESTION_H_

//#define DEBUG

#include <unordered_map>
#include "marathon/rational.h"
#include "marathon/state_graph.h"
#include "marathon/path_construction_scheme.h"

namespace marathon {

    template<typename T>
    class CongestionBoundCalculator {

    protected:

        const StateGraph &sg;               // State Graph
        const PathConstructionScheme &pcs;  // Path construction scheme
        const size_t omega;                 // number of states
        const Rational pimin;               // minimal stationary distribution
        Rational load;                      // maximal congestion

        struct pair_hash {
            inline std::size_t operator()(const std::pair<int, int> &v) const {
                return v.first * 31 + v.second;
            }
        };

        /**
         * Calculate the maximal congestion.
         * @return
         */
        Rational pathCongestion() {

            const MarkovChain &mc = sg.getMarkovChain();

            Rational max_congestion(0);
            std::unordered_map<std::pair<int, int>, Rational, pair_hash> congestion;

            if (omega <= 1)
                return Rational(0);

            // calculate normalizing constant
            Rational Z = sg.getNormalizingConstant();

#ifndef DEBUG
#pragma omp parallel
#endif
            {
                typename std::unordered_map<std::pair<int, int>, Rational, pair_hash> local_congestion;
                typename std::unordered_map<std::pair<int, int>, Rational, pair_hash>::iterator c_it;
                typename std::list<int>::const_iterator p_it;
                int u, v;

#ifndef DEBUG
#pragma omp for schedule(dynamic)
#endif
                for (int i = 0; i < omega; i++) {
                    for (int j = 0; j < omega; j++) {

                        if (i != j) {

                            // compute path from i to j
                            std::list<int> path = pcs.construct(sg, i, j);

#ifdef DEBUG
                            std::cout << "path from state " << i << " to state " << j
                                  << ":" << std::endl;
                        for (p_it = path.begin(); p_it != path.end(); ++p_it) {
                            std::cout << *p_it << " ";
                        }
                        std::cout << std::endl;
#endif

                            Rational x = sg.getWeight(i) * sg.getWeight(j)
                                         * (path.size() - 1) / (Z * Z);

                            assert(path.front() == i);
                            assert(path.back() == j);

                            if (!path.empty()) {

                                // walk along the path
                                u = i;
                                p_it = path.begin();
                                for (++p_it; p_it != path.end(); ++p_it) {
                                    v = *p_it;
                                    assert(v != -1);
                                    // augment congestion of (u,v) by pi[i]*pi[j]*|path| (number of transitions)
                                    std::pair<uint, uint> uv(u, v);
#ifdef DEBUG
                                    std::cout << "augment congestion of transition ("
                    << u << "," << v << ") by " << x
                    << std::endl;
#endif
                                    local_congestion[uv] += x;
                                    u = v;
                                }
                            }

#ifdef DEBUG
                            for (auto it = local_congestion.begin();
                    it != local_congestion.end(); ++it) {
                std::cout << "(" << it->first.first << ","
                << it->first.second << "): " << it->second
                << std::endl;
            }
#endif
                        }
                    }
                }

#ifndef DEBUG
#pragma omp critical
#endif
                {
                    for (c_it = local_congestion.begin();
                         c_it != local_congestion.end(); ++c_it) {
                        std::pair<int, int> uv = c_it->first;
                        assert(uv.first != -1);
                        assert(uv.second != -1);
                        Rational c_uv = c_it->second;
                        congestion[uv] += c_uv;
                    }
                }
            }

            int max_u, max_v;
            // compute maximal congestion
            for (auto c_it = congestion.begin(); c_it != congestion.end(); ++c_it) {
                std::pair<int, int> uv = c_it->first;
                int u = uv.first;
                Rational c = c_it->second;
                Rational P_uv = sg.getTransitionProbability(uv.first, uv.second);
                Rational Q_uv = sg.getWeight(u) * P_uv / Z;
                assert(P_uv != Rational(-1));
                c /= Q_uv;
                if (c > max_congestion) {
                    max_congestion = c;
                    max_u = uv.first;
                    max_v = uv.second;
                }
            }

#ifdef DEBUG
            std::cout << "max congestion is on edge (" << max_u << "," << max_v
<< ") with congestion " << congestion[std::make_pair(max_u, max_v)]
<< " and normalized congestion " << max_congestion << std::endl;
#endif

            return max_congestion;
        }

    public:

        CongestionBoundCalculator(
                const StateGraph &sg,
                const PathConstructionScheme &pcs)
                : sg(sg), pcs(pcs), omega(sg.getNumStates()),
                  pimin(sg.getMinWeight() / sg.getNormalizingConstant()),
                  load(-1) {

        }


        /**
         * Computes upper congestion bound by canonical path method.
         * Path construction scheme can be given by function pointer.
         * @param eps The distance to stationary distribution.
         * @return Maximum congestion of a path.
         */
        T upperPathCongestionBound(T eps) {

            if (load == -1)
                load = pathCongestion();

            return load.convert_to<T>() * -log(pimin.convert_to<T>() * eps);
        }

    };

}


#endif /* INCLUDE_MARATHON_PATHCONGESTION_H_ */
