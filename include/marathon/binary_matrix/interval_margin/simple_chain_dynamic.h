/*
 * Created on: Feb 17, 2017
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2017, Steffen Rechner, Linda Strowick
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

#ifndef _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_SIMPLE_CHAIN_DYNAMIC_H
#define _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_SIMPLE_CHAIN_DYNAMIC_H

#include <vector>
#include <stack>
#include "marathon/binary_matrix/interval_margin/simple_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace interval_margin {

            /**
             * Simple Markov chain with dynamic probability adjustment defined in
             *
             * Steffen Rechner, Linda Strowick, Matthias Müller-Hannemann.
             * Uniform sampling of bipartite graphs with degrees in prescribed intervals.
             * Journal of Complex Networks (2017). DOI: 10.1093/comnet/cnx059             *
             */
            class SimpleChainDynamic : public SimpleChain {

            protected:

                uint32_t steps_until_next_adjustment;
                const uint32_t adjustment_rate = 100;

            public:

                /**
                 * Create a Markov chain.
                 * @param m Four-Tuple of lower and upper bounds.
                 */
                explicit SimpleChainDynamic(const Instance &m) : SimpleChain(m) {
                    reset();
                }

                /**
                 * Create a Markov chain. Use a certain binary matrix as initial state.
                 * @param inst Four-Tuple of integer vectors.
                 * @param bin Binary Matrix.
                 */
                SimpleChainDynamic(const Instance &inst, const BinaryMatrix &bin)
                        : SimpleChain(inst, bin) {
                    reset();
                }


                /**
                 * Create a Markov chain.
                 * @param inst string-encoded instance. Instances for this chain have the
                 * form "l1-u1,l2-u2,l3-u3;l4-u4,l5-u5", where li is the ith lower
                 * bound and ui is the ith upper bound. For convenience, if li = ui,
                 * the string 'li-ui' can be replaced by 'li'.
                 * The semicolon separates the row sums from the column sums.
                 */
                explicit SimpleChainDynamic(const std::string &inst) :
                        SimpleChainDynamic(Instance(inst)) {

                }


                /**
                 * Create a Markov chain.
                 * @param rowsum_lower Sequence of lower bounds for row sums.
                 * @param rowsum_upper Sequence of upper bounds for row sums.
                 * @param colsum_lower Sequence of lower bounds for column sums.
                 * @param colsum_upper Sequence of upper bounds for column sums.
                 */
                SimpleChainDynamic(
                        const std::vector<int> &rowsum_lower,
                        const std::vector<int> &rowsum_upper,
                        const std::vector<int> &colsum_lower,
                        const std::vector<int> &colsum_upper
                ) : SimpleChainDynamic(
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
                SimpleChainDynamic(
                        const int *rowsum_lower,
                        const int *rowsum_upper,
                        const int *colsum_lower,
                        const int *colsum_upper,
                        const int nrow,
                        const int ncol
                ) : SimpleChainDynamic(Instance(rowsum_lower, rowsum_upper, colsum_lower, colsum_upper, nrow, ncol)) {

                }

                /**
                 * Reset the probability for selecting on operation to default and the loop count to zero.
                 */
                void reset() {

                    p_switch = p_switch_rat.convert_to<double>();
                    p_flip = p_flip_rat.convert_to<double>();
                    p_shift = p_shift_rat.convert_to<double>();

                    count_switch = 0;
                    count_flip = 0;
                    count_shift = 0;

                    loop_switch = 0;
                    loop_flip = 0;
                    loop_shift = 0;

                    steps_until_next_adjustment = adjustment_rate;
                }

                /**
                  * Set the current state of the Markov chain.
                  * A copy of state s is used from now on as current State.
                  * @param s State pointer that is used to create the current state.
                  */
                void setCurrentState(const State *s) override {
                    SimpleChain::setCurrentState(s);
                    reset();
                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s
                 * @param process
                 */
                void adjacentStates(
                        const State *x,
                        const std::function<void(const State *, const marathon::Rational &)> &process
                ) const override {
                    std::cerr << "adjacentStates method not reasonable for this Markov chain" << std::endl;
                }

                /**
                 * Randomize the state s by applying a single transition.
                 */
                void step() override {

                    steps_until_next_adjustment--;

                    if (steps_until_next_adjustment == 0) {

                        // adjust probabilites
                        if (count_switch > 0 && count_flip > 0 && count_shift > 0) {

                            // determine success rates for each operation
                            double success_switch = (double) (count_switch - loop_switch) / (double) count_switch;
                            double success_flip = (double) (count_flip - loop_flip) / (double) count_flip;
                            double success_shift = (double) (count_shift - loop_shift) / (double) count_shift;

                            /*printf("success rates:\n");
                            printf("success_switch=%f\n", success_switch);
                            printf("success_flip=%f\n", success_flip);
                            printf("success_shift=%f\n", success_shift);*/

                            // each operation gets probability that is proportional to its success rate
                            double total = success_switch + success_flip + success_shift;

                            if (total > 0.01) {

                                double p_select_switch_new = success_switch / total;
                                double p_select_flip_new = success_flip / total;
                                double p_select_shift_new = success_shift / total;

                                // determine the new selection probability as mean between old and new value
                                p_switch = (p_switch + p_select_switch_new) / 2.0;
                                p_flip = (p_flip + p_select_flip_new) / 2.0;
                                p_shift = (p_shift + p_select_shift_new) / 2.0;

                                /*printf("probabilites after adjustment:\n");
                                printf("p_switch=%f\n", p_switch);
                                printf("p_flip=%f\n", p_flip);
                                printf("p_shift=%f\n", p_shift);*/
                            }
                        }

                        steps_until_next_adjustment = adjustment_rate;
                    }

                    SimpleChain::step();
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual SimpleChainDynamic *copy() const override {
                    auto mc = new SimpleChainDynamic(inst);
                    mc->setCurrentState(this->getCurrentState());
                    return mc;
                }
            };
        }
    }
}


#endif /* _MARATHON_BINARY_MATRIX_INTERVAL_MARGIN_SIMPLE_CHAIN_DYNAMIC_H */
