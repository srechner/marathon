/*
 * Created on: Feb 17, 2017
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

#ifndef _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_INFORMED_CHAIN_DYNAMIC_H
#define _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_INFORMED_CHAIN_DYNAMIC_H

#include <vector>
#include <stack>
#include "marathon/binary_matrix/interval_margin/informed_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * Informed Markov chain with dynamic probability adjustment defined in
             *
             * Steffen Rechner, Linda Strowick, Matthias Müller-Hannemann.
             * Uniform sampling of bipartite graphs with degrees in prescribed intervals.
             * Journal of Complex Networks (2017). DOI: 10.1093/comnet/cnx059
             */
            class InformedChainDynamic : public InformedChain {

            protected:

                uint32_t steps_until_next_adjustment;
                const uint32_t adjustment_rate = 100;

            public:

                /**
                 * Create a Markov chain.
                 * @param inst Four-Tuple of lower and upper bounds.
                 */
                explicit InformedChainDynamic(const Instance &inst)
                        : InformedChain(inst) {
                    reset();
                }

                /**
                 * Create a Markov chain. Use a certain binary matrix as initial state.
                 * @param inst Four-Tuple of integer vectors.
                 * @param bin Binary Matrix.
                 */
                explicit InformedChainDynamic(const Instance &inst, const BinaryMatrix &bin)
                        : InformedChain(inst, bin) {
                    reset();
                }


                /**
                 * Create a Markov chain.
                 * @param inst String-encoded instance of lower and upper bounds on the row and column sums.
                 * Instances for this chain have the form "l1-u1,l2-u2,l3-u3;l4-u4,l5-u5", where li is the ith lower
                 * bound and ui is the ith upper bound. For convenience, if li = ui, the string 'li-ui' can be replaced
                 * by 'li'. The semicolon separates the row sums from the column sums.
                 */
                explicit InformedChainDynamic(const std::string &inst) :
                        InformedChainDynamic(Instance(inst)) {

                }


                /**
				 * Create a Markov chain.
				 * @param rowsum_lower Sequence of lower bounds for row sums.
				 * @param rowsum_upper Sequence of upper bounds for row sums.
				 * @param colsum_lower Sequence of lower bounds for column sums.
				 * @param colsum_upper Sequence of upper bounds for column sums.
				 */
                InformedChainDynamic(
                        const std::vector<int> &rowsum_lower,
                        const std::vector<int> &rowsum_upper,
                        const std::vector<int> &colsum_lower,
                        const std::vector<int> &colsum_upper
                ) : InformedChainDynamic(
                        Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper)) {

                }

                /**
				 * Create a Markov chain.
				 * @param rowsum_lower Sequence of lower bounds for row sums.
				 * @param rowsum_upper Sequence of upper bounds for row sums.
				 * @param colsum_lower Sequence of lower bounds for column sums.
				 * @param colsum_upper Sequence of upper bounds for column sums.
				 * @param nrow Number of rows.
				 * @param ncol Number of columns.
				 */
                InformedChainDynamic(
                        const int *rowsum_lower,
                        const int *rowsum_upper,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        const int nrow,
                        const int ncol
                ) : InformedChainDynamic(
                        Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper, nrow, ncol)) {

                }

                /**
                 * Reset the probability for selecting on operation to default and the loop count to zero.
                 */
                void reset() {

                    p_trade = p_trade_rat.convert_to<double>();
                    p_multiflip = p_multiflip_rat.convert_to<double>();
                    p_multishift = p_multishift_rat.convert_to<double>();

                    count_trade = 0;
                    count_multiflip = 0;
                    count_multishift = 0;

                    loop_trade = 0;
                    loop_multiflip = 0;
                    loop_multishift = 0;

                    steps_until_next_adjustment = adjustment_rate;
                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s
                 * @param process
                 */
                void adjacentStates(
                        const State &x,
                        const std::function<void(const State &, const marathon::Rational &)> &process
                ) const override {
                    throw std::runtime_error("Error! Method adjacentStates is not reasonable for time-inhomogeneous Markov chains.");
                }

                /**
                 * Randomize the state s by applying a single transition.
                 */
                void step() override {

                    steps_until_next_adjustment--;

                    if (steps_until_next_adjustment == 0) {

                        // adjust probabilites

                        if (count_trade > 0 && count_multiflip > 0 && count_multishift > 0) {

                            // determine success rates for each operation
                            double success_trade = (double) (count_trade - loop_trade) / (double) count_trade;
                            double success_multiflip =
                                    (double) (count_multiflip - loop_multiflip) / (double) count_multiflip;
                            double success_multishift =
                                    (double) (count_multishift - loop_multishift) / (double) count_multishift;

                            /*printf("success rates:\n");
                            printf("success_trade=%f\n", success_trade);
                            printf("success_multiflip=%f\n", success_multiflip);
                            printf("success_multishift=%f\n", success_multishift);*/

                            // each operation gets probability that is proportional to its success rate
                            double total = success_trade + success_multiflip + success_multishift;

                            if (total > 0.01) {

                                double p_select_trade_new = success_trade / total;
                                double p_select_multiflip_new = success_multiflip / total;
                                double p_select_multishift_new = success_multishift / total;

                                // determine the new selection probability as mean between old and new value
                                p_trade = (p_trade + p_select_trade_new) / 2.0;
                                p_multiflip = (p_multiflip + p_select_multiflip_new) / 2.0;
                                p_multishift = (p_multishift + p_select_multishift_new) / 2.0;

                                /*printf("probabilites after adjustment:\n");
                                printf("p_trade=%f\n", p_trade);
                                printf("p_multiflip=%f\n", p_multiflip);
                                printf("p_multishift=%f\n", p_multishift);*/
                            }
                        }

                        steps_until_next_adjustment = adjustment_rate;
                    }

                    // run a steps of the simple switch-flip-shift chain
                    InformedChain::step();
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual std::unique_ptr<marathon::MarkovChain> copy() const override {
                    return std::make_unique<InformedChainDynamic>(_inst, _currentState);
                }
            };
        }
    }
}


#endif /* _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_INFORMED_CHAIN_DYNAMIC_H */
