/*
 * Created on: Sep 16, 2016
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
#ifndef MARATHON_INTEGER_H
#define MARATHON_INTEGER_H

#include <vector>
#include <mutex>
#include <boost/multiprecision/cpp_int.hpp>

namespace marathon {

	/**
	 * A wrapper for boost cpp_int type
	 */
	typedef boost::multiprecision::cpp_int Integer;

	namespace detail {

		// variables for computation of binomial coefficients
		static std::vector<std::vector<Integer>> binom_table = {{1}};
		static volatile bool under_construction = false;
		static std::mutex binom_table_mutex;

		/**
		 * Resize the table of binomial coefficients
		 * @param new_nrow New number of rows.
		 * @param new_ncol New number of columns.
		 */
		inline
		void resize(size_t new_nrow, size_t new_ncol) {

			// lock table
			std::lock_guard<std::mutex> lock(detail::binom_table_mutex);

			size_t old_nrow = binom_table.size();
			size_t old_ncol = binom_table[0].size();

			if (new_nrow < old_nrow && new_ncol < old_ncol)
				return;

			under_construction = true;

			// printf("enlarge table from (%i times %i) to (%i times %i)\n", old_nrow, old_ncol, new_nrow, new_ncol);

			/********************************************************
			 *  Fill the missing parts RU and BOTTOM row by row.
			 *
			 *          0              old_ncol       new_ncol
			 *        0 +--------------|--------------+
			 *          |      LU      |      RU      |
			 * old_nrow |--------------|--------------|
			 *          |            BOTTOM           |
			 * new_nrow +-----------------------------+
			 *
			 ******************************************************/

			// fill RU part
			binom_table[0].reserve(new_ncol);
			for (size_t j = old_ncol; j < new_ncol; j++)
				binom_table[0].push_back(0);
			for (size_t i = 1; i < old_nrow; i++) {
				binom_table[i].reserve(new_ncol);
				for (size_t j = old_ncol; j < new_ncol; j++) {
					Integer bij = binom_table[i - 1][j] + binom_table[i - 1][j - 1];
					binom_table[i].push_back(bij);
				}
			}

			// fill BOTTOM part
			binom_table.reserve(new_nrow);
			for (size_t i = old_nrow; i < new_nrow; i++) {
				binom_table.push_back({1}); // binom(i,0)=1
				binom_table[i].reserve(new_ncol);
				for (size_t j = 1; j < new_ncol; j++) {
					Integer bij = binom_table[i - 1][j] + binom_table[i - 1][j - 1];
					binom_table[i].push_back(bij);
				}
			}

			/*printf("after resize:\n");
			for(int i=0; i<binom_table.size(); i++) {
				for(int j=0; j<binom_table[i].size(); j++) {
					printf("%i ", binom_table[i][j].convert_to<int>());
				}
				printf("\n");
			}*/

			under_construction = false;
		}
	}

	/**
	 * Computes n choose k.
	 */
	inline
	Integer binom(size_t n, size_t k) {

		// binom(n,k) = binom(n, n-k)
		if (k > n - k)
			return binom(n, n - k);

		// spin while table is under construction
		while(detail::under_construction);

		const size_t nrow = detail::binom_table.size();
		const size_t ncol = detail::binom_table[0].size();

		// if the result is not yet computed? extend the table!
		if (n >= nrow || k >= ncol) {

			// the new table sizes (at least +50% rows)
			const size_t new_nrow = std::max(n+1, (nrow * 3) / 2);
			const size_t new_ncol = (new_nrow / 2) + 1;
			assert(new_nrow > nrow);

			detail::resize(new_nrow, new_ncol);
		}

		return detail::binom_table[n][k];
	}

	/**
	 * Calculate the mulitinomial coefficient (k_0 + k_1 + ... + k_{n-1})! / (k_0! * k_1! * ... * k_{n-1}!)
	 * @param k Multinomial coefficients.
	 * @param n Number of coefficients.
	 * @return Multinomial coefficient (k_0 + k_1 + ... + k_{n-1})! / (k_0! * k_1! * ... * k_{n-1}!).
	 */
	inline
	marathon::Integer multinomial(
	        const int* k,
	        const int n
	) {
		marathon::Integer res = 1;
		int sum_k = 0;

		for(int i=0; i<n; i++) {
			sum_k += k[i];
			res *= binom(sum_k, k[i]);
		}

		return res;
	}


	/**
	 * Computes n!.
	 */
	inline
	Integer factorial(size_t n) {
		Integer res(1);
		for (size_t i = 2; i <= n; i++) {
			res = res * Integer(i);
		}
		return res;
	}

	/**
	 * Calculate b^e.
	 * @param b Base.
	 * @param e Exponent.
	 * @return b^e
	 */
	inline
	Integer pow(size_t b, size_t e) {

		Integer a = b;
		size_t f = e;
		Integer res(1);

		while (f > 0) {
			if (f & 1)
				res *= a;
			a = a * a;
			f = f >> 1;
		}

		return res;
	}

	/**
	 * Evaluate the polynomial P(X) = start + coef[0] * base^0 + coef[1] * base^1 + ... + coef[n-1] * base^{n-1}
	 * @param coef Array of coefficients.
	 * @param n Number of coefficients.
	 * @param base Base to which the polynomialial as evaluated.
	 * @param start Initial value.
	 * @return P(X) = start + coef[0] * base^0 + coef[1] * base^1 + ... + coef[n-1] * base^{n-1}
	 */
	inline
	marathon::Integer
	evaluate_polynomial(
			const int* coef,
			const int n,
			const int base,
	        const marathon::Integer& start = 0
	) {
		marathon::Integer res = start;
		for (int i=0; i<n; i++) {
			res = res * base + coef[i];
		}
		return res;
	}

}

#endif //MARATHON_INTEGER_H
