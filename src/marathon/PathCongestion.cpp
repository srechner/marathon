// project includes
#include <boost/unordered_map.hpp>
//#include <map>

#include "../../include/marathon/marathon.h"

namespace marathon {
namespace pathCongestion {

//#define DEBUG

rational pathCongestion(const StateGraph* sg, const PathConstructionScheme& pcs) {

	const int omega = sg->getNumStates();
	const MarkovChain* mc = sg->getMarkovChain();

	rational max_congestion(0);
	//boost::unordered_map<std::pair<int, int>, rational> congestion;
	// todo: use unordered map
	std::map<std::pair<int, int>, rational> congestion;

	if (omega <= 1)
		return rational();

	// calculate normalizing constant
	rational Z = sg->getZ();

#ifndef DEBUG
#pragma omp parallel
#endif
	{
		boost::unordered_map<std::pair<int, int>, rational> local_congestion;
		boost::unordered_map<std::pair<int, int>, rational>::iterator c_it;
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
					pcs.construct(sg, i, j, path);
					//mc->constructPath(sg, i, j, path);

#ifdef DEBUG
					std::cout << "path from state " << i << " to state " << j
					<< ":" << std::endl;
					for (p_it = path.begin(); p_it != path.end(); ++p_it) {
						std::cout << *p_it << " ";
					}
					std::cout << std::endl;
#endif

					rational x = sg->getWeight(i) * sg->getWeight(j)
							* (path.size() - 1) / (Z * Z);

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
				rational c_uv = c_it->second;
				congestion[uv] += c_uv;
			}
		}
	}

	int max_u, max_v;
	// compute maximal congestion
	for (auto c_it = congestion.begin(); c_it != congestion.end(); ++c_it) {
		std::pair<int, int> uv = c_it->first;
		int u = uv.first;
		rational c = c_it->second;
		rational P_uv = sg->getTransitionProbability(uv.first, uv.second);
		rational Q_uv = sg->getWeight(u) * P_uv / Z;
		assert(P_uv != rational(-1));
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

