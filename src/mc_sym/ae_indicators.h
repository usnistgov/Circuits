#pragma once

#include "bf.h"
#include "ae_invariants.h"
#include "ae_group.h"
#include "ae_local_connections.h"
#include "ae_helper.h"
#include "logger.h"

namespace ae {

	using namespace bfl;

	template<uint64_t N>
	bf_vector<N> get_indicators(const bf_tt<N> &f)
	{
		bf_vector<N> indicators;

		auto walsh = f.walsh_spectrum();

		auto abs_dist = walsh.absolute_distribution();

		indicators.reserve(abs_dist.size());

		for (auto &m : abs_dist)
		{
			bf_tt<N> fi(init::zero);

			for (uint64_t i = 0; i < bf_tt<N>::TableSize; i++)
				if (m.first == std::abs(walsh[i]))
					fi.set(i);

			indicators.push_back(fi);
		}

		return indicators;
	}

	template<uint64_t N>
	bf_set<N> get_indicators(const bf_tt<N> &f, uint64_t distance) {

		bf_set<N> indicators;

		for (auto const &fc : local_connections(f, distance))
		{
			auto inds = get_indicators(fc);

			std::copy(std::begin(inds), std::end(inds), std::inserter(indicators, indicators.end()));
		}

		return indicators;
	}

	template<uint64_t N>
	uint64_t indicator_max_distance(const bf_vector<N> &fvec)
	{
		assert(fvec.size() >= 2);

		const auto vecsize = fvec.size();

		for (uint64_t distance = 0; ; distance++)
		{
			LOG_TRACE << "ind checking for distance = " << distance << LOG_ENDL;

			std::vector<std::set<affine_invariants>> aivec;

			aivec.reserve(vecsize);

			bool all_distinct{ true };

			for (auto g : fvec)
			{
				std::set<affine_invariants> aiset;

				const auto h = get_local_connections(g, distance);

				LOG_DEBUG << "#local connections : " << h.size() << LOG_ENDL;

#pragma omp parallel for
				for (int lcIndex = 0; lcIndex < h.size(); lcIndex++)
				{
					auto inds = get_indicators(h[lcIndex]);

					std::set<affine_invariants> aisettemp;

					for (const auto &h2 : inds)
						aisettemp.insert(calculate_affine_invariants(h2));

#pragma omp critical
					{
						std::copy(std::begin(aisettemp), std::end(aisettemp), std::inserter(aiset, std::begin(aiset)));
						//for (const auto &h2 : inds)
						//	aiset.insert(calculate_affine_invariants(h2));
					}
				}

				LOG_DEBUG << "function has " << aiset.size() << " distinct indicator signatures" << LOG_ENDL;

				for (uint64_t i = 0; i < aivec.size(); i++)
					if (aiset == aivec[i])
					{
						all_distinct = false;
						break;
					}

				if (!all_distinct)
					break;

				aivec.push_back(aiset);
			}

			if(all_distinct)
				return distance;
		}
	}

}