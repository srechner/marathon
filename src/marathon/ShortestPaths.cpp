/*
 * ShortestPaths.cpp
 *
 * Created on: Sep 28, 2015
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

#include "marathon/Diameter.h"

#include <queue>

void marathon::Diameter::pathLengthHistogram(std::vector<long>& count,
		const StateGraph* G) {

	int n = G->getNumStates();

	// cnt[i] is number of times a length of i is observed
	long* cnt = new long[n + 1];
	memset(cnt, 0, (n + 1) * sizeof(long));
	int maxlen = 0;

	// for each state
#pragma omp parallel
	{
		// thread local variables
		int thread_max = 0;
		long* cnt_local = new long[n + 1];
		memset(cnt_local, 0, (n + 1) * sizeof(long));
		int* len = new int[n];

#pragma omp for
		for (int s = 0; s < n; s++) {

			// init temporary length array
			for (int i = 0; i < n; i++)
				len[i] = n + 1;	// not reachable

			// start breadth first search
			std::queue<int> q;
			q.push(s);
			len[s] = 0;

			while (!q.empty()) {

				int v = q.front();

				// iterate over neighbours of v
				for (Transition* t : G->getOutArcs(v)) {

					int l = len[v] + 1;

					if (l < len[t->v]) {
						len[t->v] = l;
						q.push(t->v);
					}
				}

				q.pop();
			}

			for (int i = 0; i < n; i++) {
				if (len[i] != n + 1) {
					cnt_local[len[i]]++;
					if (len[i] > thread_max)
						thread_max = len[i];
				} else {
					// should not happen if MC is irreducible
				}
			}
		}

#pragma omp critical
		{
			if (thread_max > maxlen)
				maxlen = thread_max;

			for (int i = 0; i <= n; i++)
				cnt[i] += cnt_local[i];
		}

		delete[] len;
		delete[] cnt_local;
	}

	count.assign(cnt, cnt + maxlen + 1);

	delete[] cnt;
}

int marathon::Diameter::diameter(const StateGraph* G) {
	std::vector<long> count;
	pathLengthHistogram(count, G);
	return count.size() - 1;
}

