// project includes
#include <boost/unordered_map.hpp>

#include "../../../include/marathon/cpu/analyzer.h"

namespace marathon {

namespace cpu {

//#define DEBUG

/**
 * Computes the congestion bound of state graph with given path construction scheme.
 */
Rational pathCongestion(const StateGraph* mc) {

	int omega = mc->getNumStates();

	Rational max_congestion;
	boost::unordered_map<std::pair<int, int>, Rational> congestion;

	if (omega <= 1)
		return Rational();

#pragma omp parallel
	{
		boost::unordered_map<std::pair<int, int>, Rational> local_congestion;
		boost::unordered_map<std::pair<int, int>, Rational>::iterator c_it;
		std::list<int> path;
		typename std::list<int>::const_iterator p_it;
		int u, v;

#pragma omp for schedule(dynamic)
		for (int i = 0; i < omega; i++) {
			for (int j = 0; j < omega; j++) {

				if (i != j) {

					// compute path from i to j
					path.clear();
					mc->canonicalPath(i, j, path);

#ifdef DEBUG
#pragma omp critical
					{
						for (p_it = path.begin(); p_it != path.end(); ++p_it) {
							std::cout << *p_it << " ";
						}
						std::cout << std::endl;
					}
#endif

					Rational x = mc->getStationary(i) * mc->getStationary(j)
							* (path.size() - 1);

					assert(path.front() == i);
					assert(path.back() == j);

					if (!path.empty()) {
						// walk along the path
						u = i;
						p_it = path.begin();
						for (++p_it; p_it != path.end(); ++p_it) {
							v = *p_it;
							// augment congestion of (u,v) by pi[i]*pi[j]*|path| (number of transitions)
							std::pair<uint, uint> uv(u, v);
							local_congestion[uv] += x;
							u = v;
						}
					}
				}
			}
		}

#pragma omp critical
		{
			for (c_it = local_congestion.begin();
					c_it != local_congestion.end(); ++c_it) {
				std::pair<int, int> uv = c_it->first;
				Rational c_uv = c_it->second;
				congestion[uv] += c_uv;
			}
		}
	}

	// compute maximum
	boost::unordered_map<std::pair<int, int>, Rational>::iterator c_it;
	for (c_it = congestion.begin(); c_it != congestion.end(); ++c_it) {
		std::pair<int, int> uv = c_it->first;
		int u = uv.first;
		Rational c = c_it->second;
		Rational P_uv = mc->getTransitionProbability(uv.first, uv.second);
		assert(P_uv != Rational(-1));
		c /= mc->getStationary(u) * P_uv;
		if (c > max_congestion) {
			max_congestion = c;
		}
	}

	//std::cout << "end congestion" << std::endl;

	return max_congestion;
}

}

}
