#pragma once

#include <numeric>

#include "bf.h"
#include "combinatorics.h"
#include "ae_group.h"
#include "benchmark.h"
#include "ae_helper.h"
#include "logger.h"
#include "utils.h"
#include "combinatorics.h"

namespace ae {

	using namespace bfl;
	using namespace combinatorics;

	template<uint64_t N>
	bf_vector<N> get_local_connections(const bf_tt<N> &f, uint64_t distance)
	{
		assert(distance <= bf_tt<N>::TableSize);

		// the only function at distance zero is the function itself
		if (distance == 0)
			return { f };

		bf_vector<N> connections;

		// TODO : reserve space for C(2^n, dist) functions

		for (auto &positions : subset(bf_tt<N>::TableSize, distance))
			connections.push_back(flip_entries(f, positions));

		return connections;
	}

	template<uint64_t N>
	class LocalConnections {

	public:

		template<uint64_t N>
		struct Iterator : public std::iterator<std::forward_iterator_tag, bf_tt<N>>
		{
			Iterator(const bf_tt<N> &f, subset::Iterator it)
				: f_(f), it_(it) {
			}

			bool operator != (const Iterator<N> &right) const {
				return it_ != right.it_;
			}

			Iterator<N> operator++() {
				++it_;
				return *this;
			}

			bf_tt<N> operator*() const {
				return flip_entries(f_, *it_);
			}

		private:

			bf_tt<N> f_;
			subset::Iterator it_;
		};

		LocalConnections(const bf_tt<N> &f, uint64_t distance)
			: f_(f), s_(bf_tt<N>::TableSize, distance) {
		}

		Iterator<N> begin()  {
			return Iterator<N>(f_, s_.begin());
		}

		Iterator<N> end()  {
			return Iterator<N>(f_, s_.end());
		}

	private:
		bf_tt<N> f_;
		subset s_;
	};

	template<uint64_t N>
	LocalConnections<N> local_connections(const bf_tt<N> &f, uint64_t distance) {

		return LocalConnections<N>(f, distance);
	}

	template<uint64_t N>
	uint64_t local_connection_max_distance(const bf_vector<N> &fvec)
	{
		assert(fvec.size() >= 2);

		const auto vecsize = fvec.size();
			
		for (uint64_t distance = 1; ; distance++)
		{
			LOG_DEBUG << "lc checking for distance = " << distance << LOG_ENDL;

			std::vector<std::set<affine_invariants>> aivec(vecsize);

//			aivec.reserve(vecsize);

			//bool all_distinct{ true };

			//for (auto g : fvec)

#pragma omp parallel for
			for(int x = 0; x < fvec.size(); x++)
			{
				//std::set<affine_invariants> aiset;

				auto &g = fvec[x];

				for (auto h : local_connections(g, distance))
					//aiset.insert(calculate_affine_invariants(h));
					aivec[x].insert(calculate_affine_invariants(h));

				//for (uint64_t i = 0; i < aivec.size(); i++)
				//	if (aiset == aivec[i])
				//	{
				//		all_distinct = false;
				//		break;
				//	}

				//if (!all_distinct)
				//	break;

				//aivec.push_back(aiset);
			}

			//if (all_distinct)
			if(all_distinct(aivec))
				return distance;
		}
	}

}

