/*
 * TransitionMatrix.h
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#ifndef INCLUDE_MARATHON_TRANSITIONMATRIX_H_
#define INCLUDE_MARATHON_TRANSITIONMATRIX_H_

#include "StateGraph.h"

namespace marathon {

    /**
     * Virtual Base Class for Transition Matrix.
     */
    template<typename T=double>
    class TransitionMatrix {

    protected:

        size_t n;        // number of rows and columns
        size_t ld;        // lead dimension (upper bound on n)
        T *data;        // actual data

        /**
         * Return a pointer to a Transition matrix of an appropriate subtype.
         */
        virtual TransitionMatrix<T> *generateSubTypeInstance(const int n) = 0;

    public:

        virtual ~TransitionMatrix() {

        }

        /**
         *  Return size of the matrix.
         */
        size_t getN() const {
            return n;
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
        T *getData() const {
            return data;
        }

        /**
         * Copy the content of matrix P to this.
         */
        virtual void copy(const TransitionMatrix<T> *P) = 0;

        /**
         * Overwrite the current matrix with unity matrix.
         */
        virtual void setEye() = 0;

        /**
         * Overwrite the current matrix with zeroes.
         */
        virtual void setZero() = 0;

        /**
         * Multiply A with B and write the result to this.
         * @param A A pointer to matrix A. Will not be changed.
         * @param B A pointer to matrix B. Will not be changed.
         */
        virtual void mult(const TransitionMatrix<T> *A,
                          const TransitionMatrix<T> *B) = 0;

        /**
         * Compute P^k and write the result to this.
         * @param P A pointer to a Transition Matrix.
         * @param k Exponent.
         */
        void pow(const TransitionMatrix<T> *P, const int k) {

            const int omega = P->n;

            // create temporary matrix of the same subtype as calling instance.
            TransitionMatrix<T> *tmp = generateSubTypeInstance(omega);

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
        }

        /**
         * Return a string that represents the matrix.
         */
        virtual std::string to_string() const = 0;

        /**
         * Compute the variation distance of each state to the distribution pi
         * @param pi A pointer to a probability distribution.
         * @param dist Out parameter.
         *
         */
        virtual void variationDistance(const T *pi, T *dist) const = 0;

        /**
         * Compute the total variation distance to the distribution.
         * @param pi A probability distribution.
         */
        virtual T totalVariationDistance(const T *pi) const = 0;

        /**
         * Swap the content of the Matrix with another matrix.
         *
         */
        void swap(TransitionMatrix<T> *P) {
            std::swap(n, P->n);
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

}

#endif /* INCLUDE_MARATHON_TRANSITIONMATRIX_H_ */
