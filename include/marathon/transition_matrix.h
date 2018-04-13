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

#ifndef INCLUDE_MARATHON_TRANSITIONMATRIX_H_
#define INCLUDE_MARATHON_TRANSITIONMATRIX_H_

#include "state_graph.h"

#ifdef USE_ARMADILLO

#include <armadillo>

#endif

#ifdef USE_BLAS

#include <cblas.h>

#endif

namespace marathon {

    /**
     * Virtual Base Class for Transition Matrix.
     */
    template<class T=double>
    class TransitionMatrix {

    protected:

        size_t N;               // number of rows and columns
        size_t ld;              // lead dimension (upper bound on n)
        std::vector<T> data;    // actual data array

    public:

        /**
         * Standard Constructor. Create uninitialized transition matrix of size N times N.
         * @param N number of rows or columns
         */
        TransitionMatrix(const size_t N) :
                N(N),
                ld(((N + 255) / 256) * 256) // lead dimension is next mulitple of 256
        {
            data.resize(N * ld, 0);
        }

        /**
         * Constructor. Create Transition Matrix from State Graph.
         * @param sg Pointer to state graph object.
         */
        TransitionMatrix(const StateGraph &sg) :
                TransitionMatrix(sg.getNumStates()) {

            for (const Transition *t : sg.getArcs()) {
                this->data[t->from * ld + t->to] = t->weight.convert_to<T>();
            }
        }

        /**
         *  Return size of the matrix.
         */
        size_t getDimension() const {
            return N;
        }

        /**
         * Return lead dimension of the matrix.
         */
        size_t getLeadDimension() const {
            return ld;
        }

        /**
         * Return a pointer to the data.
         */
        const std::vector<T> &getData() const {
            return data;
        }

        /**
         * Return P[i,j].
         * @param i row index
         * @param j column index
         * @return P[i,j]
         */
        T get(size_t i, size_t j) const {
            return data[i * ld + j];
        }

        /**
         * Set P[i,j] to x.
         * @param i row index.
         * @param j column index.
         * @param x value of type T
         */
        void set(size_t i, size_t j, T x) {
            data[i * ld + j] = x;
        }

        /**
         * Overwrite the current matrix with zeroes.
         */
        virtual void clear() {
            data.resize(N * ld, T(0));
        }

        /**
         * Compute P^k.
         * @param P A pointer to a Transition Matrix.
         * @param k Exponent.
         * @return P^k
         */
        TransitionMatrix<T> pow(uint k) const {

            // init matrix
            if (k == 0) {
                return eye(N);
            }

            // create binary representation of k
            int bin[32];
            memset(bin, 0, 32 * sizeof(int));
            int l = 31;
            int kk = k;
            while (kk > 0) {
                bin[l] = kk % 2;
                kk >>= 1;
                l--;
            }
            l += 2;

#ifdef DEBUG
            std::cout << "bin: ";
            for (int i = 0; i < 32; i++) {
                std::cout << bin[i];
            }
            std::cout << " l=" << l << std::endl;
#endif


            TransitionMatrix<T> A(*this); // will be returned

            // binary exponentation - Left to Right (see Don. Knuth: Seminumerical Alg. Vol. 2 page 461)
            while (l < 32) {

                // square
                A = A * A;

                // multiply
                if (bin[l] == 1)
                    A = A * *this;

                l++;
            }

            return A;
        }

        /**
         * Matrix multiplication.
         * @param P Transition matrix.
         * @return P * this
         */
        TransitionMatrix<T> operator*(const TransitionMatrix<T> &P) const {

            TransitionMatrix<T> X(N);  // will be returned

#pragma omp parallel for
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    T p_ij = 0;
                    for (size_t k = 0; k < N; k++) {
                        p_ij += this->get(i, k) * P.get(k, j);
                    }
                    X.set(i, j, p_ij);
                }
            }

            return X;
        }

        /**
         * Return a string that represents the matrix.
         */
        virtual std::string to_string() const {

            std::stringstream ss;
            ss << "\n";
            for (size_t i = 0; i < this->N; i++) {
                ss << "  ";
                for (size_t j = 0; j < this->N - 1; j++) {
                    ss << std::setprecision(std::numeric_limits<T>::digits10) << std::fixed
                       << this->data[i * this->ld + j] << "  ";
                }
                ss << std::setprecision(std::numeric_limits<T>::digits10) << std::fixed
                   << this->data[i * this->ld + this->N - 1];
                ss << "\n";
            }

            return ss.str();

        }

        /**
         * To output into streams.
         */
        friend inline std::ostream &operator<<(std::ostream &out,
                                               const TransitionMatrix<T> &s) {
            out << s.to_string();
            return out;
        }

        /**
         * Return the identity matrix with N rows and columns.
         * @param N Number of rows and columns.
         * @return Identity matrix.
         */
        static TransitionMatrix<T> eye(size_t N) {
            TransitionMatrix<T> P(N);
            for (size_t i = 0; i < N; i++)
                P.set(i,i,1);
            return P;
        }

    };


    /***********************************************************************
     * template specializations
     **********************************************************************/

#ifdef USE_BLAS

    template<>
    TransitionMatrix<float> TransitionMatrix<float>::operator*(const TransitionMatrix<float> &P) const {

        const float alpha = 1.0;
        const float beta = 0.0;

        TransitionMatrix<float> X(N);

        // use cblas
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, alpha,
                    &P.data[0], P.ld, &data[0], ld, beta, &X.data[0], X.ld);

        return X;
    }

    template<>
    TransitionMatrix<double> TransitionMatrix<double>::operator*(const TransitionMatrix<double> &P) const {

        const double alpha = 1.0;
        const double beta = 0.0;

        TransitionMatrix<double> X(N);

        // use cblas
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, alpha,
                    &P.data[0], P.ld, &data[0], ld, beta, &X.data[0], X.ld);

        return X;
    }

#endif

}

#endif /* INCLUDE_MARATHON_TRANSITIONMATRIX_H_ */
