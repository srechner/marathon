/*
 * BipartiteMatching.h
 *
 * Created on: Nov 20, 2014
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

#ifndef STATE_JS89_H_
#define STATE_JS89_H_

#include <cstring>
#include <cstdlib>
#include <ostream>

#include "marathon/State.h"

namespace marathon {
	namespace matching {

		class BipartiteMatching : public State {

		public:

			int n, k;                // number of nodes			// number of edges
			int unmatched[2];    // indices of unmatched nodes (if any)
			int *mates;

			BipartiteMatching() :
					n(0), k(0) {
				mates = (int *) malloc(n * sizeof(int));
			}

			BipartiteMatching(const BipartiteMatching &s) :
					n(s.n), k(s.k) {
				unmatched[0] = s.unmatched[0];
				unmatched[1] = s.unmatched[1];
				mates = (int *) malloc(n * sizeof(int));
				memcpy(mates, s.mates, n * sizeof(int));
			}

			BipartiteMatching(int n, int k, int unmatched[2],
			                  int *matching) :
					n(n), k(k) {
				this->unmatched[0] = unmatched[0];
				this->unmatched[1] = unmatched[1];
				mates = (int *) malloc(n * sizeof(int));
				memcpy(this->mates, matching, n * sizeof(int));
			}

			~BipartiteMatching() {
				free(mates);
			}

			void addEdge(int u, int v) {
				mates[u] = v;
				mates[v] = u;
				k++;
				unmatched[0] = n;
				unmatched[1] = n;
			}

			void removeEdge(int u, int v) {
				mates[u] = n;
				mates[v] = n;
				k--;
				unmatched[0] = u;
				unmatched[1] = v;
			}

			void operator=(BipartiteMatching const &s) {
				n = s.n;
				k = s.k;
				unmatched[0] = s.unmatched[0];
				unmatched[1] = s.unmatched[1];
				mates = (int *) realloc(mates, n * sizeof(int));
				memcpy(mates, s.mates, n * sizeof(int));
			}

			bool operator==(const BipartiteMatching &s) const {
				if (n != s.n)
					return false;
				const int ret = memcmp(mates, s.mates, n * sizeof(int));
				return (ret == 0);
			}

			bool operator<(const BipartiteMatching &s) const {
				return memcmp(mates, s.mates, n * sizeof(int));
			}

			size_t hashValue() const {
				return unmatched[0] * n + unmatched[1];
			}

			int compare(const State *x) const {
				const BipartiteMatching *b = (const BipartiteMatching *) x;
				const int res = memcmp(this->mates, b->mates, this->n * sizeof(int));
				return res;
			}

			std::string toString() const {

				std::stringstream ss;

				ss << "[";
				ss.width(log10(n) + 1);
				for (int i = 0; i < n - 1; i++) {
					ss.width(log10(n) + 1);
					if (mates[i] == n)
						ss << "-";
					else
						ss << mates[i];
					ss << ", ";
				}
				ss.width(log10(n) + 1);
				if (mates[n - 1] == n)
					ss << "-";
				else
					ss << mates[n - 1];
				ss << "]";

				return ss.str();
			}

			State *copy() const {
				return new BipartiteMatching(*this);
			}

			bool is_perfect() const {
				return 2 * k == n;
			}

			bool is_near_perfect() const {
				return 2 * k == n - 2;
			}
		};

	}
}


#endif /* STATE_H_ */
