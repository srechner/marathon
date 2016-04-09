// project includes
#include "../../../../include/marathon/chain/bipgraph/HavelHakimi.h"
#include "../../../../include/marathon/chain/bipgraph/SwitchChain.h"

//#define DEBUG

namespace marathon {
namespace chain {
namespace bipgraph {

SwitchChain::SwitchChain(const std::string& line, int seed) :
		MarkovChain(line, seed), sum(0) {
	parseInstance(line);
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

State* SwitchChain::computeArbitraryState() {

	// not a valid instance
	if (u.size() < 2 || v.size() < 2) {
		return nullptr;
	}

	bool M[u.size() * v.size()];
	if (HavelHakimiBipartite(u, v, M) != 0) {
		return new DenseBipartiteGraph(u.size(), v.size(), M);
	} else
		return nullptr;
}

void SwitchChain::computeNeighbours(const State* x,
		std::vector<std::pair<State*, rational>>& neighbors) const {

	const DenseBipartiteGraph* s = (const DenseBipartiteGraph*) x;

	const long m = s->get_nrows();
	const long n = s->get_ncols();

	// each state has proposal prob. of 4 / (m*(m+1)*n*(n+1))
	const rational p(4, m * (m + 1) * n * (n + 1));

	// construct neighbours in parallel
#pragma omp parallel
	{
		// each thread has its own set of neighbours
		std::vector<std::pair<State*, rational>> myneighbours;
		rational myloop(0);

#pragma omp for
		// Definition of Kannan, Tetali, Vempala
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = i; k < m; k++) {
					for (int l = j; l < n; l++) {

						// alternating cycle ( (i,j)=1, (k,j)=0, (k,l)=1, (i,l)=0 )
						if (s->has_edge(i, j) && s->has_edge(k, l)
								&& !s->has_edge(i, l) && !s->has_edge(k, j)) {

							// make a copy
							DenseBipartiteGraph *s2 = new DenseBipartiteGraph(
									*s);

							// switch the cycle
							s2->flip_edge(i, j); // (i,j) = 0
							s2->flip_edge(k, l); // (k,l) = 0
							s2->flip_edge(i, l); // (i,l) = 1
							s2->flip_edge(k, j); // (k,j) = 1

							myneighbours.push_back(std::make_pair(s2, p));
						}
						// symmetric case (attention: not regarded in paper!)
						else if (!s->has_edge(i, j) && !s->has_edge(k, l)
								&& s->has_edge(i, l) && s->has_edge(k, j)) {

							// make a copy
							DenseBipartiteGraph *s2 = new DenseBipartiteGraph(
									*s);

							// switch the cycle
							s2->flip_edge(i, j); // (i,j) = 1
							s2->flip_edge(k, l); // (k,l) = 1
							s2->flip_edge(i, l); // (i,l) = 0
							s2->flip_edge(k, j); // (k,j) = 0

							myneighbours.push_back(std::make_pair(s2, p));
						} else {
							// loop
							myloop += p;
						}
					}
				}
			}
		}

		// copy thread own neighbors into global neighbour array
#pragma omp critical
		{
			DenseBipartiteGraph *s2 = new DenseBipartiteGraph(*s);
			neighbors.push_back(std::make_pair(s2, myloop));
			neighbors.insert(neighbors.end(), myneighbours.begin(),
					myneighbours.end());
		}
	}
}

void SwitchChain::randomize(State* x) const {

	DenseBipartiteGraph* s = (DenseBipartiteGraph*) x;

	const int nrows = u.size();
	const int ncols = v.size();

	// select four random integers i,j,k,l

	// 0 <= i <= k < nrows
	int i = rand() % nrows;
	int k = rand() % nrows;
	if (i > k)
		std::swap(i, k);

	// 0 <= j <= l < ncols
	int j = rand() % ncols;
	int l = rand() % ncols;
	if (j > l)
		std::swap(j, l);

	// check if edges are flippable
	const bool ij = s->has_edge(i, j);
	const bool kl = s->has_edge(k, l);
	const bool il = s->has_edge(i, l);
	const bool kj = s->has_edge(k, j);

	// if i,j,k,l makes a switchable cycle
	if (ij == kl && il == kj && ij != kj) {

		// switch the cycle
		s->flip_edge(i, j);
		s->flip_edge(k, l);
		s->flip_edge(i, l);
		s->flip_edge(k, j);
	}

}

}
}
}
