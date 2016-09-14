/*
 * TransitionMatrixCPU.h
 *
 *  Created on: Mar 22, 2016
 *      Author: steffen
 */

#ifndef INCLUDE_MARATHON_TRANSITIONMATRIXCBLAS_H_
#define INCLUDE_MARATHON_TRANSITIONMATRIXCBLAS_H_

#include "TransitionMatrix.h"

namespace marathon {
	namespace tm {

		template<typename T>
		class TransitionMatrixCBLAS : public TransitionMatrix<T> {

		protected:

			/**
			 * Copy the content of matrix P to this.
			 */
			virtual void copy(const TransitionMatrix<T> *P) {
				const TransitionMatrixCBLAS<T> *X = (const TransitionMatrixCBLAS<T> *) P;
				this->n = X->n;
				this->ld = X->ld;
				memcpy(this->data, X->data, this->ld * this->n * sizeof(T));
			}

			/**
			 * Return a pointer to a Transition matrix of subtype instance.
			 */
			virtual TransitionMatrix<T> *generateSubTypeInstance(const int n) {
				TransitionMatrix<T> *res = new TransitionMatrixCBLAS<T>(n);
				return res;
			}

		public:

			TransitionMatrixCBLAS(const int n) {
				this->n = n;
				this->ld = n;
				this->data = new T[this->n * this->ld];
			}

			TransitionMatrixCBLAS(const StateGraph *sg) :
					TransitionMatrixCBLAS<T>(sg->getNumStates()) {

				setZero();

				for (const Transition *t : sg->getArcs()) {
					this->data[t->u * this->ld + t->v] = t->p.convert_to<T>();
				}
			}

			virtual ~TransitionMatrixCBLAS() {
				delete[] this->data;
			}

			/**
			 * Overwrite the current matrix with unity matrix.
			 */
			virtual void setEye() {
				setZero();
				for (int i = 0; i < this->n; i++) {
					this->data[i * this->ld + i] = 1;
				}
			}

			/**
			 * Overwrite the current matrix with zeroes.
			 */
			virtual void setZero() {
				size_t bytes = this->n * this->ld * sizeof(T);
				memset(this->data, 0, bytes);
			}

			/**
			 * Return a string that represents the matrix.
			 */
			virtual std::string to_string() const {
				std::stringstream ss;
				//ss << "[ ";
				for (size_t i = 0; i < this->n; i++) {
					ss << " ";
					for (size_t j = 0; j < this->n - 1; j++) {
						ss << std::setprecision(std::numeric_limits<T>::digits10) << std::fixed
						   << this->data[i * this->ld + j] << " ";
					}
					ss << std::setprecision(std::numeric_limits<T>::digits10) << std::fixed
					   << this->data[i * this->ld + this->n - 1];

					/*if (i < this->n - 1)
						ss << ";\n";
					else
						ss << " ]";*/
					ss << "\n";
				}

				return ss.str();
			}

			/**
			 * Multiply A with B and write the result to this.
			 * @param A A pointer to matrix A. Will not be changed.
			 * @param B A pointer to matrix B. Will not be changed.
			 */

			virtual void mult(const TransitionMatrix<T> *A,
			                  const TransitionMatrix<T> *B);

			/**
			 * Compute the variation distance of each state to the distribution pi
			 * @param pi A pointer to a probability distribution.
			 * @param dist Out parameter.
			 */
			virtual void variationDistance(const T *pi, T *dist) const {

#pragma omp parallel for if(this->n > 1000)
				for (int i = 0; i < this->n; i++) {
					T sum = 0;
					for (int j = 0; j < this->n; j++)
						sum += fabs(this->data[i * this->ld + j] - pi[j]);
					dist[i] = sum / 2.0;
				}
			}

			/**
			 * Compute the total variation distance to the distribution.
			 * @param pi A probability distribution.
			 */
			virtual T totalVariationDistance(const T *pi) const {
				T *tmp = new T[this->n];
				variationDistance(pi, tmp);
				delete[] tmp;
				return *std::max_element(tmp, tmp + this->n);
			}

		};

//typedef TransitionMatrixCBLAS TransitionMatrixCBLAS<double>;

	}
}

#endif /* INCLUDE_MARATHON_TRANSITIONMATRIXCBLAS_H_ */
