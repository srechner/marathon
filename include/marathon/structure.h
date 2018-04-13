/*
 * Created on: Sep 23, 2016
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

#ifndef PROJECT_DIAMETER_H
#define PROJECT_DIAMETER_H

#include "state_graph.h"

namespace marathon {

    /**
     * Determine whether or not the given state graph is biparite
     * @param G Pointer to a state graph.
     * @return True, if G is bipartite, False, if not.
     */
    bool isBipartite(const StateGraph *G) {

        const size_t N = G->getNumStates();

        // try to color the states red and blue.
        int *color = new int[N];
        for (int i = 0; i < N; i++)
            color[i] = -1;

        // start graph traversal
        std::stack<int> s;
        s.push(0);
        color[0] = 0;

        bool cont = true;                   // used for early termination
        while (!s.empty() && cont) {

            int i = s.top();                // index of current state
            int c = color[i];               // color of state i
            s.pop();

            // for all adjacent states
            for (const auto& t : G->getOutArcs(i)) {

                int j = t->to;               // j is adjacent to i

                if (color[j] == 0) {         // j is not yet visited
                    // give j the opposite color of i and add j to stack
                    color[j] = 1 - c;
                    s.push(j);
                } else if (color[j] == c) {    // adjacent state with same color detected
                    cont = false;
                    break;
                }
            }
        }

        delete[] color;

        return cont;
    }


    /**
     * Computes a histogram of the length of shortest path in G. Each arc of the graph
     * contributes one to the length of a path. For each length l for l in 0..diameter(G),
     * the number of shortest paths of G is stored in the vector count[l].
     */
    void pathLengthHistogram(std::vector<long> &count, const StateGraph &G) {
        int n = G.getNumStates();

        // cnt[i] is number of times a length of i is observed
        long *cnt = new long[n + 1];
        memset(cnt, 0, (n + 1) * sizeof(long));
        int maxlen = 0;

        // for each state
#pragma omp parallel
        {
            // thread local variables
            int thread_max = 0;
            long *cnt_local = new long[n + 1];
            memset(cnt_local, 0, (n + 1) * sizeof(long));
            int *len = new int[n];

#pragma omp for
            for (int s = 0; s < n; s++) {

                // init temporary length array
                for (int i = 0; i < n; i++)
                    len[i] = n + 1;    // not reachable

                // start breadth first search
                std::queue<int> q;
                q.push(s);
                len[s] = 0;

                while (!q.empty()) {

                    int v = q.front();

                    // iterate over neighbours of v
                    for (Transition *t : G.getOutArcs(v)) {

                        int l = len[v] + 1;

                        if (l < len[t->to]) {
                            len[t->to] = l;
                            q.push(t->to);
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


    /**
     * Calculate the shortest path connecting the node 'from' with 'to'.
     * The length of the path is measured by its number of arcs.
     * @param sg State graph.
     * @param from Node index.
     * @param to Node index.
     * @return Length of shortest path.
     */
    int distance(const StateGraph& sg, const int from, const int to) {

        // run bfs
        std::queue< std::pair<int,int> > q;
        std::vector<bool> visited(sg.getNumStates());

        // insert start vertex with distance 0
        q.push(std::make_pair(from, 0));
        visited[from] = true;

        while(!q.empty()) {

            // extract node
            const auto& vd = q.front();
            const int v = vd.first;
            const int d = vd.second;
            q.pop();

            // early termination
            if(v == to) {
                return d;
            }

            // for all outgoing arcs (v,w)
            for(const auto& t : sg.getOutArcs(v)) {
                const int w = t->to;
                if(!visited[w]) {
                    visited[w] = true;
                    q.push(std::make_pair(w,d+1));
                }
            }
        }

        // cannot happen
        return -1;
    }

    /**
     * Computes the diameter of the graph, i.e. the maximal length of a shortest path
     * between some nodes of the graph. Each arc has a length of 1.
     */
    int diameter(const StateGraph &G) {
        std::vector<long> count;
        pathLengthHistogram(count, G);
        return count.size() - 1;
    }

}


#endif //PROJECT_DIAMETER_H
