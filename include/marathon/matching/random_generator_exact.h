/*
 * Created on: Jan 11, 2018
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2017, Steffen Rechner
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

#ifndef MARATHON_MATCHING_RANDOM_GENERATOR_EXACT_H
#define MARATHON_MATCHING_RANDOM_GENERATOR_EXACT_H

// marathon includes
#include "marathon/random_device.h"
#include "marathon/random_generator.h"
#include "marathon/matching/enumerate.h"

namespace marathon {
    namespace matching {

        /**
         * A random generator that produces uniformly distributed perfect
         * and near-perfect matchings of a bipartite graph.
         */
        class RandomGeneratorExact :
                public marathon::RandomGenerator,
                private marathon::matching::Counter {

        public:

            /**
             * Options for the random generator.
             * uniform: Create uniformly distributed random samples.
             * jsv04: Create random samples according to the distribution
             *        that is specified in Jerrum et al. (2004).
             */
            enum distribution_t {
                uniform, jsv04
            };

        protected:

            const distribution_t dist;

            marathon::RandomDevice rg;               // random generator
            std::vector<std::pair<int, int>> Nuv;    // list of node pairs with non-empty sample sets
            marathon::Integer Z;                     // total number of samples
            BipartiteMatching match;                 // will be returned

            void exact_sample_recursive(
                    const size_t u,
                    const size_t unmatched1,
                    const size_t unmatched2,
                    const Integer &target
            ) {

                //std::cout << "u=" << u << " unmatched=" << unmatched1 << " unmatched2=" << unmatched2 << " target=" << target << std::endl;

                // if each vertex of vertex set U has been considered
                if (u == _n) {

                    // prepare matching
                    if (unmatched1 != SIZE_MAX) {
                        match.mates[unmatched1] = SIZE_MAX;
                        match.mates[unmatched2] = SIZE_MAX;
                    }

                    // stop recursion
                    return;
                }

                // if node u is already matched (are marked as permanently unmatched)
                if (match.isMatched(u)) {
                    exact_sample_recursive(u + 1, unmatched1, unmatched2, target);
                    return;
                }

                // simulate all choices to add an edge (u,v) to the matching

                marathon::Integer skipped(0);

                // iterate over all adjacent nodes
                std::vector<size_t> neighbors = _g.getNeighbors(u);

                // Each choice of k results in a number of val(k) matchings.
                // We seek the smallest k such that val(0)+val(1)+...+val(k) > target
                size_t k = 0;
                for (; k < neighbors.size(); k++) {

                    const size_t v = neighbors[k];

                    if (!match.isMatched(v)) {

                        // add edge (u,v) to matching
                        match.mates[u] = v;
                        match.mates[v] = u;
                        match._edges++;

                        // count the number of matrices that result from this choice
                        Integer val = count_recursive(match, u + 1);

                        //std::cout << "k=" << k << ": val=" << val << std::endl;

                        // found the right sub-tree?
                        if (skipped + val > target) {
                            break;
                        }

                        skipped += val;

                        // undo changes
                        match.mates[u] = SIZE_MAX;
                        match.mates[v] = SIZE_MAX;
                        match._edges--;
                    }
                }

                //std::cout << "k=" << k << " neighbors.size()=" << neighbors.size() << std::endl;
                assert(k < neighbors.size());

                const size_t v = neighbors[k];

                // reduce the target by the number of elements skipped and continue sampling
                exact_sample_recursive(u + 1, unmatched1, unmatched2, target - skipped);
            }


            /**
             * Return a uniformly distributed random sample.
             * @return
             */
            const BipartiteMatching &next_uniform() {

                // prepare Matching
                match._edges = 0;
                for (size_t i = 0; i < 2 * _n; i++) {
                    match.mates[i] = SIZE_MAX;
                }

                // determine random number
                const Integer target = rg.nextInt(Z);

                //std::cout << "Z=" << Z << std::endl;
                //std::cout << "target=" << target << std::endl;

                Integer skipped(0);

                // count perfect matchings
                Integer val;
                bool found = load_from_table(match, val);
                assert(val);

                //std::cout << "num perfect = " << val << std::endl;

                // if target is an index of a perfect matching
                if (skipped + val > target) {
                    exact_sample_recursive(0, SIZE_MAX, SIZE_MAX, target);
                    return match;
                }

                skipped += val;

                // target is an index of a near-perfect matching
                for (size_t u = 0; u < _n; u++) {
                    for (size_t v = _n; v < 2 * _n; v++) {

                        // mark (u,v) as permanently unmatched
                        match.mates[u] = SIZE_MAX - 1;
                        match.mates[v] = SIZE_MAX - 1;

                        found = load_from_table(match, val);
                        assert(found);

                        //std::cout << "unmatched u=" << u << " v=" << v << ": " << match->toString() << ": skipped=" << skipped << " val=" << val << std::endl;

                        // if target is an index of a perfect matching
                        if (skipped + val > target) {
                            exact_sample_recursive(0, u, v, target - skipped);
                            return match;
                        }

                        skipped += val;

                        // undo marking
                        match.mates[u] = SIZE_MAX;
                        match.mates[v] = SIZE_MAX;

                    }
                }
            }

            /**
             * Return a random sample distributed according to the distribution
             * defined by Jerrum et al. 2004.
             * @return
             */
            const BipartiteMatching &next_jsv04() {

                // prepare Matching
                match._edges = 0;
                for (size_t i = 0; i < 2 * _n; i++) {
                    match.mates[i] = SIZE_MAX;
                }

                // First select one of (Nuv.size()+1) subsets of samples uniformly at random
                int k = rg.nextInt(Nuv.size() + 1);

                if (k == Nuv.size()) {
                    // draw a perfect matching uniformy at random
                    const Integer Z = countPerfect();
                    const Integer target = rg.nextInt(Z);
                    exact_sample_recursive(0, SIZE_MAX, SIZE_MAX, target);
                    return match;
                } else {

                    // determine unmatched nodes
                    int u = Nuv[k].first;
                    int v = Nuv[k].second;

                    // prepare matching
                    match.mates[u] = SIZE_MAX - 1;
                    match.mates[v] = SIZE_MAX - 1;

                    // draw a near-perfect matching uniformy at random from the set
                    // of near-perfect matchings with u and v being unmatch
                    const Integer Z = countNearPerfect(u, v);
                    const Integer target = rg.nextInt(Z);
                    exact_sample_recursive(0, u, v, target);
                    return match;
                }
            }


        public:

            /**
             * Create a random generator that produces perfect and near-perfect matchings of a bipartite graph
             * according to either the uniform or the jsv04 distribution.
             * @param g Bipartite graph.
             * @param dist Specifies from which distribution samples are produced. Allowed values are: uniform, jsv04
             */
            RandomGeneratorExact(
                    SparseBipartiteGraph g,
                    distribution_t dist = uniform
            ) : Counter(std::move(g)), dist(dist) {

                // todo: use binary search

                match = BipartiteMatching(_g.getNumberOfNodes());
                Z = count();

                // create list of non-empty near-perfect matching sets
                for (size_t u = 0; u < _n; u++) {
                    for (size_t v = _n; v < 2 * _n; v++) {
                        Integer x = countNearPerfect(u, v);
                        if (x > 0) {
                            Nuv.push_back(std::make_pair(u, v));
                        }
                    }
                }
            }

            /**
             * Return a random perfect or near-perfect matching.
             * @return Random binary matrix.
             */
            const BipartiteMatching &next() override {

                switch (dist) {
                    case uniform:
                        return next_uniform();
                    case jsv04:
                        return next_jsv04();
                    default:
                        throw std::runtime_error("Error! Unknown distribution type!");
                }

            }

            /**
             * Return an independent copy of the random generator.
             * @return Copy of the random generator.
             */
            std::unique_ptr<RandomGenerator> copy() const override {
                return std::make_unique<RandomGeneratorExact>(*this);
            }
        };
    }
}

#endif //MARATHON_MATCHING_RANDOM_GENERATOR_EXACT_H
