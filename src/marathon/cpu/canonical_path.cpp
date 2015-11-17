// project includes
#include <boost/unordered_map.hpp>
//#include <map>

#include "../../../include/marathon/cpu/analyzer.h"

namespace marathon {

namespace cpu {

//#define DEBUG

/**
 * Computes the congestion bound of state graph with given path construction scheme.
 */
Rational pathCongestion(const StateGraph* mc) {

	int omega = mc->getNumStates();

	Rational max_congestion(0);
	boost::unordered_map<std::pair<int, int>, Rational> congestion;

	if (omega <= 1)
		return Rational();

#ifndef DEBUG
#pragma omp parallel
#endif
	{
		boost::unordered_map<std::pair<int, int>, Rational> local_congestion;
		boost::unordered_map<std::pair<int, int>, Rational>::iterator c_it;
		std::list<int> path;
		typename std::list<int>::const_iterator p_it;
		int u, v;

#ifndef DEBUG
#pragma omp for schedule(dynamic)
#endif
		for (int i = 0; i < omega; i++) {
			for (int j = 0; j < omega; j++) {

				if (i != j) {

					// compute path from i to j
					path.clear();
					mc->canonicalPath(i, j, path);

#ifdef DEBUG
					std::cout << "path from state " << i << " to state " << j
							<< ":" << std::endl;
					for (p_it = path.begin(); p_it != path.end(); ++p_it) {
						std::cout << *p_it << " ";
					}
					std::cout << std::endl;
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
#ifdef DEBUG
							std::cout << "augment congestion of transition ("
									<< u << "," << v << ") by " << x
									<< std::endl;
#endif
							local_congestion[uv] += x;
							u = v;
						}
					}

#ifdef DEBUG
					for (auto it = local_congestion.begin();
							it != local_congestion.end(); ++it) {
						std::cout << "(" << it->first.first << ","
								<< it->first.second << "): " << it->second
								<< std::endl;
					}
#endif
				}
			}
		}

#ifndef DEBUG
#pragma omp critical
#endif
		{
			for (c_it = local_congestion.begin();
					c_it != local_congestion.end(); ++c_it) {
				std::pair<int, int> uv = c_it->first;
				Rational c_uv = c_it->second;
				congestion[uv] += c_uv;
			}
		}
	}

	int max_u, max_v;
	// compute maximal congestion
	for (auto c_it = congestion.begin(); c_it != congestion.end(); ++c_it) {
		std::pair<int, int> uv = c_it->first;
		int u = uv.first;
		Rational c = c_it->second;
		Rational P_uv = mc->getTransitionProbability(uv.first, uv.second);
		Rational Q_uv = mc->getStationary(u) * P_uv;
		assert(P_uv != Rational(-1));
		c /= Q_uv;
		if (c > max_congestion) {
			max_congestion = c;
			max_u = uv.first;
			max_v = uv.second;
		}
	}

#ifdef DEBUG
	std::cout << "max congestion is on edge (" << max_u << "," << max_v
			<< ") with congestion " << congestion[std::make_pair(max_u, max_v)]
			<< " and normalized congestion " << max_congestion << std::endl;
#endif

	return max_congestion;
}

}

}
