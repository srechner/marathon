/*
 * State.cpp
 *
 *  Created on: Nov 20, 2014
 *      Author: rechner
 */

// project includes
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "../../../../include/marathon/chains/matching/bipartite_matching.h"

namespace marathon {

namespace chain {

namespace matching {

BipartiteMatching::BipartiteMatching() :
		n(0), k(0) {
	mates = (int*) malloc(n * sizeof(int));
}

BipartiteMatching::BipartiteMatching(const BipartiteMatching& s) :
		n(s.n), k(s.k) {
	unmatched[0] = s.unmatched[0];
	unmatched[1] = s.unmatched[1];
	mates = (int*) malloc(n * sizeof(int));
	memcpy(mates, s.mates, n * sizeof(int));
}

BipartiteMatching::BipartiteMatching(int n, int k, int unmatched[2],
		int* matching) :
		n(n), k(k) {
	this->unmatched[0] = unmatched[0];
	this->unmatched[1] = unmatched[1];
	mates = (int*) malloc(n * sizeof(int));
	memcpy(this->mates, matching, n * sizeof(int));
}

BipartiteMatching::~BipartiteMatching() {
	free(mates);
}

void BipartiteMatching::operator=(BipartiteMatching const& s) {
	n = s.n;
	k = s.k;
	unmatched[0] = s.unmatched[0];
	unmatched[1] = s.unmatched[1];
	mates = (int*) realloc(mates, n * sizeof(int));
	memcpy(mates, s.mates, n * sizeof(int));
}

bool BipartiteMatching::operator==(const BipartiteMatching &s) const {
	if (n != s.n)
		return false;
	int ret = memcmp(mates, s.mates, n * sizeof(int));
	return (ret == 0);
}

std::ostream& operator<<(std::ostream& out, BipartiteMatching& s) {

	out << "[";
	out.width(log10(s.n) + 1);
	for (int i = 0; i < s.n - 1; i++) {
		out.width(log10(s.n) + 1);
		if (s.mates[i] == s.n)
			out << "-";
		else
			out << s.mates[i];
		out << ", ";
	}
	out.width(log10(s.n) + 1);
	if (s.mates[s.n - 1] == s.n)
		out << "-";
	else
		out << s.mates[s.n - 1];
	out << "]";

	return out;
}

bool BipartiteMatching::is_perfect() const {
	return 2 * k == n;
}

bool BipartiteMatching::is_near_perfect() const {
	return 2 * k == n - 2;
}

void BipartiteMatching::addEdge(int u, int v) {
	mates[u] = v;
	mates[v] = u;
	k++;
	unmatched[0] = n;
	unmatched[1] = n;
}

void BipartiteMatching::removeEdge(int u, int v) {
	mates[u] = n;
	mates[v] = n;
	k--;
	unmatched[0] = u;
	unmatched[1] = v;
}

size_t BipartiteMatching::hash_value() const {
	return unmatched[0] * n + unmatched[1];
}

}
}
}

