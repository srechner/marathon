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

        size_t N;       // number of rows and columns
        size_t ld;      // lead dimension (upper bound on n)
        T *data;          // actual data array

    public:

        /**
         * Standard Constructor. Create uninitialized transition matrix of size N times N.
         * @param N number of rows or columns
         */
        TransitionMatrix(const size_t N) :
                N(N),
                ld(((N + 255) / 256) * 256) // lead dimension is next mulitple of 256
        {
            data = new T[N * ld];
        }

        /**
         * Copy constructor
         * @param t Transition matrix
         */
        TransitionMatrix(const TransitionMatrix &t) : TransitionMatrix(t.getDimension()) {
            copy(&t);
        }

        /**
         * Constructor. Create Transition Matrix from State Graph.
         * @param sg Pointer to state graph object.
         */
        TransitionMatrix(const StateGraph *sg) :
                TransitionMatrix(sg->getNumStates()) {

            setZero();

            for (const Transition *t : sg->getArcs()) {
                this->data[t->from * ld + t->to] = t->weight.convert_to<T>();
            }
        }

        /**
         * Standard Destructor
         */
        virtual ~TransitionMatrix() {
            delete[] data;
        }

        /**
         *  Return size of the matrix.
         */
        uint32_t getDimension() const {
            return N;
        }

        /**
         * Return lead dimension of the matrix.
         */
        uint32_t getLeadDimension() const {
            return ld;
        }

        /**
         * Return a pointer to the data.
         */
        T *getData() const {
            return data;
        }

        /**
         * Return P[i,j].
         * @param i row index
         * @param j column index
         * @return P[i,j]
         */
        T get(const size_t i, const size_t j) const {
            return data[i * ld + j];
        }

        /**
         * Set P[i,j] to x.
         * @param i row index.
         * @param j column index.
         * @param x value of type T
         */
        void set(const size_t i, const size_t j, const T x) {
            data[i * ld + j] = x;
        }

        /**
         * Overwrite the current matrix with unity matrix.
         */
        virtual void setEye() {
            setZero();
            for (int i = 0; i < this->N; i++) {
                this->data[i * this->ld + i] = 1;
            }
        }

        /**
         * Overwrite the current matrix with zeroes.
         */
        virtual void setZero() {
            // todo: add template specialization
            for (size_t i = 0; i < getDimension(); i++) {
                for (size_t j = 0; j < getDimension(); j++) {
                    set(i, j, T(0));
                }
            }
        }


        /**
         * Copy the content of P to this matrix.
         * @param T Pointer to transition matrix.
         */
        virtual void copy(const TransitionMatrix *P) {
            printf("generic copy\n");
            for (size_t i = 0; i < getDimension(); i++) {
                for (size_t j = 0; j < getDimension(); j++) {
                    set(i, j, P->get(i, j));
                }
            }
        }

        /**
         * Multiply A with B and write the result to this.
         * @param A A pointer to matrix A. Will not be changed.
         * @param B A pointer to matrix B. Will not be changed.
         */
        virtual void mult(const TransitionMatrix<T> *A,
                          const TransitionMatrix<T> *B) {

            printf("generic mult\n");
#pragma omp parallel for
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    T p_ij = 0;
                    for (size_t k = 0; k < N; k++) {
                        p_ij += A->get(i, k) * B->get(k, j);
                    }
                    set(i, j, p_ij);
                }
            }
        }

        /**
         * Compute P^k and write the result to this.
         * @param P A pointer to a Transition Matrix.
         * @param k Exponent.
         * @param tmp Transition matrix used for temporary memory.
         */
        void pow(const TransitionMatrix<T> *P, const int k, TransitionMatrix<T> *tmp = nullptr) {

            const int omega = P->N;

            // create temporary matrix of the same subtype as calling instance.
            bool tmp_given = tmp != nullptr;
            if (!tmp_given) {
                tmp = new TransitionMatrix<T>(omega);
            }

            // init matrix
            if (k == 0) {
                this->setEye();
            } else {
                this->copy(P);
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

            // binary exponentation - Left to Right (see Knuth Seminumerical Alg. Vol. 2 page 461)
            while (l < 32) {

                // square
                tmp->mult(this, this);
                std::swap(this->data, tmp->data);

                // multiply
                if (bin[l] == 1) {
                    // this = this*P
                    tmp->mult(this, P);
                    std::swap(this->data, tmp->data);
                }

                l++;
            }

            if (!tmp_given)
                delete tmp;
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
         * Swap the content of the Matrix with another matrix.
         *
         */
        void swap(TransitionMatrix<T> *P) {
            std::swap(N, P->N);
            std::swap(ld, P->ld);
            std::swap(data, P->data);
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
         * To output into streams.
         */
        friend inline std::ostream &operator<<(std::ostream &out,
                                               const TransitionMatrix<T> *s) {
            out << s->to_string();
            return out;
        }
    };

    /***********************************************************************
     * template specializations
     **********************************************************************/

    template<>
    void TransitionMatrix<float>::copy(const TransitionMatrix<float> *P) {
        memcpy(data, P->data, N * ld * sizeof(float));
    }

    template<>
    void TransitionMatrix<double>::copy(const TransitionMatrix<double> *P) {
        memcpy(data, P->data, N * ld * sizeof(double));
    }

    template<>
    void TransitionMatrix<float>::setZero() {
        size_t bytes = N * ld * sizeof(float);
        memset(data, 0, bytes);
    }

    template<>
    void TransitionMatrix<double>::setZero() {
        size_t bytes = N * ld * sizeof(double);
        memset(data, 0, bytes);
    }

#ifdef USE_BLAS

    template<>
    void TransitionMatrix<float>::mult(const TransitionMatrix<float> *A,
                                       const TransitionMatrix<float> *B) {

        const float alpha = 1.0;
        const float beta = 0.0;

        // use cblas
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->N, this->N,
                    this->N, alpha, B->getData(), B->getLeadDimension(), A->getData(),
                    A->getLeadDimension(), beta, this->data, this->ld);

    }

    template<>
    void TransitionMatrix<double>::mult(const TransitionMatrix<double> *A,
                                        const TransitionMatrix<double> *B) {

        const double alpha = 1.0;
        const double beta = 0.0;

        // use cblas
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->N, this->N,
                    this->N, alpha, B->getData(), B->getLeadDimension(), A->getData(),
                    A->getLeadDimension(), beta, this->data, this->ld);

    }

#endif

}

#endif /* INCLUDE_MARATHON_TRANSITIONMATRIX_H_ */
