/*
 * Created on: Feb 06, 2017
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

#ifndef BINARY_MATRIX_METRICS_H_
#define BINARY_MATRIX_METRICS_H_

#ifdef USE_ARMADILLO
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#endif

#include "binary_matrix.h"
#include "marathon/rational.h"
#include "marathon/integer.h"

namespace marathon {
    namespace binary_matrix {

        /**
         * Return the number of 2x2 submatrices of the form
         * (0 1)    or   (1 0)
         * (1 0)         (0 1)
         * normalized by the total number of 2x2 submatrices.
         * @param bin Binary Matrix.
         * @return The rate of checkerboard units.
         */
        template<class T=double>
        T checkerBoardRate(const BinaryMatrix &bin) {

            const size_t nrow = bin.getNumRows();
            const size_t ncol = bin.getNumCols();

            size_t num = 0;
#pragma omp parallel
            {
                size_t my_num = 0;
#pragma omp for
                for (int i1 = 0; i1 < nrow; i1++) {
                    for (int i2 = i1 + 1; i2 < nrow; i2++) {
                        for (int j1 = 0; j1 < ncol; j1++) {
                            for (int j2 = j1 + 1; j2 < ncol; j2++) {
                                if (bin.isCheckerBoardUnit(i1, j1, i2, j2)) {
                                    my_num++;
                                }
                            }
                        }
                    }
                }
#pragma omp critical
                num += my_num;
            }

            return T(num * 4) / T(nrow * (nrow - 1) * ncol * (ncol - 1));
        }

        /**
         * Calculate the nestedness defined by Robert and Stone, 1990.
         *
         *   A. Roberts and L. Stone.
         *   Island-sharing by archipelago species. Oecologia 83 (1990), 560–567. doi: 10.1007/bf00317210.
         *
         * @param bin Binary Matrix.
         * @return S2 Nestedness score.
         */
        template<class T=double>
        T nestednessS2(const BinaryMatrix &bin) {

            const size_t nrow = bin.getNumRows();
            const size_t ncol = bin.getNumCols();

            size_t S = 0;
            for (int i = 0; i < nrow; i++) {
                for (int j = i + 1; j < nrow; j++) {

                    int sij = 0;
                    for (int k = 0; k < ncol; k++) {
                        sij += bin.get(i, k) * bin.get(j, k);
                    }
                    S += sij * sij;
                }
            }

            return T(S * 2) / T(nrow * (nrow - 1));
        }

        /**
         * Calculate the nestedness defined by Patterson and Attmar.
         *
         *   B. D. Patterson and W. Atmar.
         *   Nested subsets and the structure of insular mammalian faunas and archipelagos.
         *   Biological Journal of the Linnean Society 28 (1986), 65–82.
         *   doi: 10.1111/j.1095-8312.1986.tb01749.x.
         *
         * @param bin Binary Matrix.
         * @return Nestedness score.
         */
        template<class T=double>
        T nestedSubset(const BinaryMatrix &bin) {

            const size_t nrow = bin.getNumRows();
            const size_t ncol = bin.getNumCols();

            // determine column sums
            int *colsum = new int[ncol];
            memset(colsum, 0, ncol * sizeof(int));
            for (int i = 0; i < nrow; i++) {
                for (int j = 0; j < ncol; j++) {
                    colsum[j] += bin.get(i, j);
                }
            }

            // for each row: determine minimal column sum
            int *m = new int[nrow];
            for (int i = 0; i < nrow; i++) {
                m[i] = nrow;
                for (int j = 0; j < ncol; j++) {
                    if (bin.get(i, j)) {
                        m[i] = std::min(m[i], colsum[j]);
                    }
                }
            }

            // calculate nestedness
            size_t S = 0;
            for (int i = 0; i < nrow; i++) {
                for (int j = 0; j < ncol; j++) {
                    if (!bin.get(i, j) && colsum[j] > m[i]) {
                        S++;
                    }
                }
            }

            delete[] colsum;
            delete[] m;

            return T(S);
        }

        /**
         * Calculate the number of entries that differ from each other.
         * @param bin1 A binary matrix.
         * @param bin2 A binary matrix.
         * @return Hamming distance between matrices.
         */
        size_t hammingDistance(const BinaryMatrix &bin1, const BinaryMatrix &bin2) {
            boost::dynamic_bitset<> x = bin1.getBitset() ^ bin2.getBitset();
            return x.count();
        }


        /**
         * Calculate the nestedness measure based on overlap and decreasing fill.
         *
         *   M. Almeida-Neto, P. Guimarães, P. R. Guimarães, R. D. Loyola, and W. Ulrich.
         *   A consistent metric for nestedness analysis in ecological systems: reconciling
         *   concept and measurement.
         *   Oikos 117 (2008), 1227–1239. doi: 10.1111/j.0030-1299.2008.16644.x.
         *
         * @param bin A binary matrix.
         * @return NODF
         */
        template<class T=double>
        T NODF(const BinaryMatrix &bin) {

            const int nrow = bin.getNumRows();
            const int ncol = bin.getNumCols();

            // determine row and column sums
            int *rowsum = new int[nrow];
            int *colsum = new int[ncol];
            memset(rowsum, 0, nrow * sizeof(int));
            memset(colsum, 0, ncol * sizeof(int));
            for (int i = 0; i < nrow; i++) {
                for (int j = 0; j < ncol; j++) {
                    if (bin.get(i, j)) {
                        rowsum[i]++;
                        colsum[j]++;
                    }
                }
            }

            /*
             * score_col is the sum of scores from pairwise comparison of
             * each column j1 against all columns j2 to its right.
             * If the degree of column j1 is smaller or equal to the degree
             * of column j2, then the score is zero.
             * If they have different degrees, the score is the percentage
             * of elements in the second column which also appear in the first
             * column.
             * score_rows is found similarly for pairwise comparisons of each
             * row against all rows below it.
             * The sum of score_col and score_row is then normalized by the total
             * number of pairwise comparisons.
             */

            T score_col = 0;
            T score_row = 0;

            // for each pair of columns
            for (int j1 = 0; j1 < ncol; j1++) {
                for (int j2 = j1 + 1; j2 < ncol; j2++) {

                    // if the degree of column j1 is smaller or equal to the degree of column j2
                    if (colsum[j1] <= colsum[j2]) {
                        // score_col += 0
                    } else {
                        // determine the percentage of elements in the
                        // second column which also appear in the first
                        int score = 0;
                        for (int i = 0; i < nrow; i++) {
                            if (bin.get(i, j2) && bin.get(i, j1)) {
                                score++;
                            }
                        }
                        score_col += T(score) / T(colsum[j2]);
                    }
                }
            }

            // for each pair of rows
            for (int i1 = 0; i1 < nrow; i1++) {
                for (int i2 = i1 + 1; i2 < nrow; i2++) {

                    // if both rows have the same degree
                    if (rowsum[i1] <= rowsum[i2]) {
                        // score_row += 0
                    } else {
                        // determine the percentage of elements in the
                        // second row which also appear in the first
                        int score = 0;
                        for (int j = 0; j < ncol; j++) {
                            if (bin.get(i1, j) && bin.get(i2, j)) {
                                score++;
                            }
                        }
                        score_row += T(score) / T(rowsum[i2]);
                    }
                }
            }

            delete[] rowsum;
            delete[] colsum;

            T nodf = T(100 * 2) * (score_row + score_col) / T(nrow * (nrow - 1) + ncol * (ncol - 1));

            return nodf;
        }


#ifdef USE_ARMADILLO


        /**
         * Calculate the spectral radius of the symmetric adjacency matrix
         * defined by expanding the bi-adjacency matrix bin.
         *
         *   P. P. Staniczenko, J. C. Kopp, and S. Allesina. The ghost of nestedness in
         *   ecological networks. Nature Communications 4 (Jan. 2013), 1391.
         *   doi: 10.1038/ncomms2422.
         *
         * @tparam T Can be float or double.
         * @param bin A binary matrix.
         * @return The absolute value of the maximum real eigenvalue
         * from the adjacency matrix defined by bin.
         */
        template<class T=double>
        T spectralRadius(const BinaryMatrix &bin) {

            const int nrow = bin.getNumRows();
            const int ncol = bin.getNumCols();
            const int total = bin.getTotal();

            // transform binary matrix into armadillo sparse matrix format
            arma::umat locations(2, 2*total);
            arma::Col<T> values(2*total);
            int k = 0;
            for(int i=0; i<nrow; i++) {
                for(int j=0; j<ncol; j++) {
                    if(bin.get(i,j)) {
                        locations(0, k) = i;
                        locations(1, k) = nrow + j;
                        values(k) = T(1);
                        k++;

                        locations(0, k) = nrow + j;
                        locations(1, k) = i;
                        values(k) = T(1);
                        k++;
                    }
                }
            }

            // define sparse matrix
            arma::SpMat<T> A(locations, values, nrow+ncol, nrow+ncol);

            // find largest eigenvalues
            arma::Col<T> eigval;

            k=1;                    // number of eigenvalues
            const int K_MAX = 10;   // maximal number of eigenvalues

            // while not successful calculate the largest k eigenvalues
            while(!arma::eigs_sym(eigval, A, k) || eigval.size() != k) {
                k++;
                if(k > K_MAX) {
                    throw std::runtime_error("Error while calculating eigenvalues!");
                }
            }

            return fabs(eigval[k-1]);
        }

#endif
    }
}


#endif /* BINARY_MATRIX_METRICS_H_ */
