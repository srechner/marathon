// project includes
#include "../../../../include/marathon/chain/bipgraph/HavelHakimi.h"
#include "../../../../include/marathon/chain/bipgraph/SwitchChain.h"
#include "../../../../include/marathon/Random.h"

//#define DEBUG

namespace marathon {
namespace chain {
namespace bipgraph {

SwitchChain::SwitchChain(const std::string& line) :
		sum(0) {
	parseInstance(line);
}

SwitchChain::SwitchChain(const std::vector<int>& u, const std::vector<int>& v) :
		u(u), v(v), sum(0) {
	for (int i = 0; i < u.size(); i++) {
		sum += u[i];
	}
}

SwitchChain::SwitchChain(int* const rowsum, int* const colsum, const int nrow,
		const int ncol) {

	u = std::vector<int>(rowsum, rowsum + nrow);
	v = std::vector<int>(colsum, colsum + ncol);

	for (int i = 0; i < u.size(); i++) {
		sum += u[i];
	}
}

SwitchChain::~SwitchChain() {
}

void SwitchChain::parseInstance(const std::string& line) {
	std::string copy(line);

	// split string at ";"
	std::vector<std::string> vec;

	std::string delimiter = ";";
	size_t pos = 0;
	std::string token;
	while ((pos = copy.find(delimiter)) != std::string::npos) {
		token = copy.substr(0, pos);
		vec.push_back(token);
		copy.erase(0, pos + delimiter.length());
	}
	vec.push_back(copy);

	if (vec.size() != 2) {
		std::cerr << "invalid syntax: " << line << std::endl;
		return;
	}

	delimiter = ",";
	pos = 0;

	copy = vec[0];
	while ((pos = copy.find(delimiter)) != std::string::npos) {
		token = copy.substr(0, pos);
		u.push_back(atoi(token.c_str()));
		copy.erase(0, pos + delimiter.length());
	}
	u.push_back(atoi(copy.c_str()));

	copy = vec[1];
	while ((pos = copy.find(delimiter)) != std::string::npos) {
		token = copy.substr(0, pos);
		v.push_back(atoi(token.c_str()));
		copy.erase(0, pos + delimiter.length());
	}
	v.push_back(atoi(copy.c_str()));

	for (auto it = v.begin(); it != v.end(); ++it)
		sum += *it;

#ifdef DEBUG
	std::cout << "u=[";
	for (std::vector<int>::iterator it = u.begin(); it != u.end(); ++it) {
		std::cout << *it << " ";
	}
	std::cout << "]" << std::endl;

	std::cout << "v=[";
	for (std::vector<int>::iterator it = v.begin(); it != v.end(); ++it) {
		std::cout << *it << " ";
	}
	std::cout << "]" << std::endl;
#endif

}

State* SwitchChain::computeArbitraryState() const {

	// not a valid instance
	if (u.size() < 1 || v.size() < 1) {
		return nullptr;
	}

	return HavelHakimiBipartite(u, v);
}

void SwitchChain::computeNeighbours(const State* x,
		std::vector<std::pair<State*, rational>>& neighbors) const {

	const BinaryMatrix* s = (const BinaryMatrix*) x;

	const long m = u.size();
	const long n = v.size();

	// each state has proposal prob. of 4 / (m*(m+1)*n*(n+1))
	const rational p(4, m * (m + 1) * n * (n + 1));

	rational loop(0);

	// Definition of Kannan, Tetali, Vempala
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = i; k < m; k++) {
				for (int l = j; l < n; l++) {

					/**
					 * A switch is possible in the following situations
					 *
					 * a) alternating cycle ( (i,j)=1, (k,j)=0, (k,l)=1, (i,l)=0 )
					 *         i   k
					 *       1 2 3 4 5
					 *     -----------
					 * j 1 | x 1 x 0 x
					 * l 2 | x 0 x 1 x
					 *   3 | x x x x x
					 *
					 * b) symmetric case (attention: not regarded in paper!)
					 *         i   k
					 *       1 2 3 4 5
					 *     -----------
					 * j 1 | x 0 x 1 x
					 * l 2 | x 1 x 0 x
					 *   3 | x x x x x
					 *
					 */

					// situation a
					if (s->get(i, j) && s->get(k, l)
							&& !s->get(i, l) && !s->get(k, j)) {

						// make a copy
						BinaryMatrix *s2 = new BinaryMatrix(*s);

						// switch the cycle
						s2->flip(i, j); // (i,j) = 0
						s2->flip(k, l); // (k,l) = 0
						s2->flip(i, l); // (i,l) = 1
						s2->flip(k, j); // (k,j) = 1

						neighbors.push_back(std::make_pair(s2, p));
					}
					// situation b
					else if (!s->get(i, j) && !s->get(k, l)
							&& s->get(i, l) && s->get(k, j)) {

						// make a copy
						BinaryMatrix *s2 = new BinaryMatrix(*s);

						// switch the cycle
						s2->flip(i, j); // (i,j) = 1
						s2->flip(k, l); // (k,l) = 1
						s2->flip(i, l); // (i,l) = 0
						s2->flip(k, j); // (k,j) = 0

						neighbors.push_back(std::make_pair(s2, p));
					} else {
						// loop
						loop += p;
					}
				}
			}
		}
	}
	// create loop
	BinaryMatrix *s2 = new BinaryMatrix(*s);
	neighbors.push_back(std::make_pair(s2, loop));
}

void SwitchChain::randomize(State* x, const uint32_t t) const {

	BinaryMatrix* s = (BinaryMatrix*) x;

	const int nrows = u.size();
	const int ncols = v.size();

	for (int tt = 0; tt < t; tt++) {

		// select four random integers i,j,k,l

		// 0 <= i <= k < nrows
		int i = ::marathon::random::nextInt(nrows);
		int k = ::marathon::random::nextInt(nrows);
		if (i > k)
			std::swap(i, k);

		// 0 <= j <= l < ncols
		int j = ::marathon::random::nextInt(ncols);
		int l = ::marathon::random::nextInt(ncols);
		if (j > l)
			std::swap(j, l);

		//std::cout << i << " " << k << " " << j << " " << l << std::endl;

		// check if edges are flippable
		const bool ij = s->get(i, j);
		const bool kl = s->get(k, l);
		const bool il = s->get(i, l);
		const bool kj = s->get(k, j);

		// if i,j,k,l makes a switchable cycle
		if (ij == kl && il == kj && ij != kj) {

			// switch the cycle
			s->flip(i, j);
			s->flip(k, l);
			s->flip(i, l);
			s->flip(k, j);
		}
	}
}

rational SwitchChain::loopProbability(const State* x) const {

	const BinaryMatrix* s = (const BinaryMatrix*) x;

	const long m = u.size();
	const long n = v.size();

	// each state has proposal prob. of 4 / (m*(m+1)*n*(n+1))
	const rational p(4, m * (m + 1) * n * (n + 1));

	rational loop(0);

	// Definition of Kannan, Tetali, Vempala
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = i; k < m; k++) {
				for (int l = j; l < n; l++) {

					/**
					 * A switch is possible in the following situations
					 *
					 * a)
					 *         i   k
					 *       1 2 3 4 5
					 *     -----------
					 * j 1 | x 1 x 0 x
					 * l 2 | x 0 x 1 x
					 *   3 | x x x x x
					 *
					 * b)
					 *         i   k
					 *       1 2 3 4 5
					 *     -----------
					 * j 1 | x 0 x 1 x
					 * l 2 | x 1 x 0 x
					 *   3 | x x x x x
					 *
					 */
					if (s->get(i, j) != s->get(k, l)
							|| s->get(i, l) != s->get(k, j)) {
						loop += p;
					}
				}
			}
		}
	}

	return loop;
}

}
}
}
