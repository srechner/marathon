/*
 * shortest_paths.cpp
 *
 *  Created on: Sep 28, 2015
 *      Author: rechner
 */

#include "../../include/marathon/marathon.h"

#include <queue>
#include <iostream>

void marathon::pathLengthHistogram(std::vector<long>& count,
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

int marathon::diameter(const StateGraph* G) {
	std::vector<long> count;
	pathLengthHistogram(count, G);
	return count.size() - 1;
}

