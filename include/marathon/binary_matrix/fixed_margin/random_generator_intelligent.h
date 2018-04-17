/*
 * Created on: May 14, 2017
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


#ifndef MARATHON_BINARY_MATRIX_FIXED_MARGIN_RANDOM_GENERATOR_INTELLIGENT_H
#define MARATHON_BINARY_MATRIX_FIXED_MARGIN_RANDOM_GENERATOR_INTELLIGENT_H

#include "marathon/binary_matrix/binary_matrix.h"
#include "marathon/random_generator.h"
#include "marathon/binary_matrix/fixed_margin/random_generator_mcmc.h"
#include "marathon/binary_matrix/fixed_margin/random_generator_exact.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            class RandomGeneratorIntelligent :
                    public marathon::binary_matrix::RandomGenerator {

            private:

                // lists of sample generators and primitive sequences
                std::vector<Instance> primitives;
                std::vector<RandomGenerator *> random_generators;
                std::vector<bool> inverse;

                // binary matrix
                BinaryMatrix *bin;

            public:

                /**
                 * Create a generator that produces (approximately) uniformly distributed
                 * random binary matrices whose row and column sums match prescribed sequences.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns
                 */
                RandomGeneratorIntelligent(
                        const int *rowsum,
                        const int *colsum,
                        const int nrow,
                        const int ncol
                ) {

                    // check realizability
                    if (!isRealizable(rowsum, colsum, nrow, ncol))
                        return;

                    // decompose sequence into primitive sequences
                    decompose(rowsum, colsum, nrow, ncol, [this](
                                      const marathon::binary_matrix::fixed_margin::Instance &seq) {

                                  Instance copy(seq);

                                  // if there are more ones than zeroes in the table
                                  const int total = std::accumulate(copy._rowsum, copy._rowsum + copy.getNumRows, 0);
                                  if (2 * total > copy.getNumRows * copy.getNumCols) {
                                      // inverse the sequence
                                      for (int i = 0; i < copy.getNumRows; i++)
                                          copy._rowsum[i] = copy.getNumCols - copy._rowsum[i];
                                      for (int j = 0; j < copy.getNumCols; j++)
                                          copy._colsum[j] = copy.getNumRows - copy._colsum[j];
                                      inverse.push_back(true);
                                  } else {
                                      inverse.push_back(false);
                                  }

                                  // add sequence
                                  primitives.push_back(copy);

                                  /**
                                   * Create a primitive exact sampler for each primitive sequence.
                                   * For small sequences use exact sampling, for large sequences use
                                   * approximate sampling.
                                   */
                                  RandomGenerator *rg;

                                  // if sequence is small enough
                                  if (copy.getNumRows == 1 || copy.getNumRows * copy.getNumCols < 100) {
                                      // use exact sampling
                                      //rg = new RandomGeneratorExact(
                                      //       copy.rowsum, copy.colsum, copy.nrow, copy.ncol);
                                  } else {
                                      // use approximate sampling
                                      rg = new RandomGeneratorMCMC(
                                              copy._rowsum, copy._colsum, copy.getNumRows, copy.getNumCols);
                                  }

                                  // add random generator
                                  random_generators.push_back(rg);
                              }
                    );

                    // create binary matrix object
                    bin = new BinaryMatrix(nrow, ncol);
                }

                /**
                 * Create a generator that produces (exactly!) uniformly distributed
                 * random binary matrices whose row and column sums match prescribed sequences.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 */
                RandomGeneratorIntelligent(
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) : RandomGeneratorIntelligent(&rowsum[0], &colsum[0], rowsum.size(), colsum.size()) {

                }

                ~RandomGeneratorIntelligent() {
                    for (int i = 0; i < random_generators.size(); i++)
                        delete random_generators[i];
                    delete bin;
                }

                /**
                 * Return a binary matrix of size nrow times ncol whose
                 * row and column sums match the given integer vectors.
                 * The sample is distributed (approximately) uniform from the set
                 * of all binary matrices with the prescribed row sums and
                 * column sums.
                 * @return Random binary matrix.
                 */
                const BinaryMatrix *next() override {

                    // if sequence is not realizable
                    if (primitives.size() == 0)
                        return nullptr;

                    // special case if sequence is already primitive:
                    if (primitives.size() == 1)
                        return random_generators[0]->next();

                    // independently sample primitive sequences
                    for (int l = 0; l < primitives.size(); l++) {
                        const BinaryMatrix *sub = random_generators[l]->next();
                        const Instance &seq = primitives[l];

                        // translate submatrix into final matrix
                        for (int i = 0; i < seq.getNumRows; i++) {
                            const int ii = seq._rowindex[i];
                            for (int j = 0; j < seq.getNumCols; j++) {
                                const int jj = seq._colindex[j];
                                bool b = sub->get(i, j);
                                bin->set(ii, jj, b ^ inverse[l]);
                            }
                        }
                    }

                    return bin;
                }
            };
        }
    }
}

#endif //MARATHON_BINARY_MATRIX_FIXED_MARGIN_RANDOM_GENERATOR_INTELLIGENT_H
