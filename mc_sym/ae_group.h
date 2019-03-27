#pragma once

#include <algorithm>
#include <vector>
#include <set>
#include <unordered_set>
#include <omp.h>

#include "bf.h"
#include "ae_invariants.h"
#include "ae_transformation.h"
#include "benchmark.h"
#include "logger.h"

namespace ae {

	using namespace bfl;

	template<uint64_t N>
	struct group
	{
		group() {
		}

		group(const affine_invariants &a) : ai(a) {
		}

		group(affine_invariants &&a) 
			: ai(std::move(a)) {
		}

		group(const group<N> &rhs) 
			:	ai(rhs.ai), 
				functions(rhs.functions) {
		}

		group(group<N> &&rhs)
			:	ai(std::move(rhs.ai)), 
				functions(std::move(rhs.functions)) {
		}

		group<N>& operator=(const group<N> &rhs) {
			ai			= rhs.ai;
			functions	= rhs.functions;
			return *this;
		}

		group<N>& operator=(group<N> &&rhs) {
			ai			= std::move(rhs.ai);
			functions	= std::move(rhs.functions);
			return *this;
		}


		affine_invariants ai;
		mutable bf_vector<N> functions;
	};

	template<uint64_t N>
	using group_vector = std::vector<group<N>>;

	template<uint64_t N>
	using group_set = std::set<group<N>>;

	template<uint64_t N>
	bool operator==(const group<N> &gr1, const group<N> &gr2) {
		return gr1.ai == gr2.ai;
	}

	template<uint64_t N>
	bool operator<(const group<N> &gr1, const group<N> &gr2) {
		return gr1.ai < gr2.ai;
	}

	template<uint64_t N>
	auto find_group(const bf_tt<N> &f, const group_vector<N> &groups) -> decltype(std::begin(groups))
	{
		auto ai = calculate_affine_invariants(f);

		auto p = std::lower_bound(std::begin(groups), std::end(groups), group<N>(ai));

		assert(p != std::end(groups) && p->ai == ai);

		return p;
	}

	template<uint64_t N>
	uint64_t get_group_index(const bf_tt<N> &f, const group_vector<N> &groups)
	{
		auto gr = find_group(f, groups);

		auto index = std::distance(std::begin(groups), gr);

		return index;
	}

	template<uint64_t N, typename Cont>
	std::set<uint64_t> get_group_indices(const Cont &cont, const group_vector<N> &groups)
	{
		std::set<uint64_t> indices;

		for (auto f : cont)
			indices.insert(get_group_index(f, groups));

		return indices;
	}

	template<uint64_t N>
	[[deprecated("inefficient")]] group_vector<N> generate_groups(const bf_vector<N> &funcs)
	{
		benchmark bm("generate_groups", benchmark::report::end);

		group_vector<N> groups;

		for (const auto &f : funcs)
		{
			auto ai = calculate_affine_invariants(f);

			// check whether the signature is a new one or an existing one
			auto p = std::find(std::begin(groups), std::end(groups), group<N>(ai));

			if (p == std::end(groups))
			{
				// signature not found, create a new group
				group<N> gr;
				gr.ai = ai;
				gr.functions.push_back(f);
				groups.push_back(gr);
			}
			else
			{
				// add to the group with the same signature
				p->functions.push_back(f);
			}
		}

		// the groups must be sorted with respect to affine invariants in order to utilize binary search
		std::sort(std::begin(groups), std::end(groups));

		return groups;
	}

	template<uint64_t N, typename SetCont>
	inline bool add_function_to_groupset(SetCont &cont, const bf_tt<N> &f) {

		auto gr = group<N>(calculate_affine_invariants(f));

		auto p = cont.insert(gr);

		// independent of whether there was an insertion or not,
		// p points to the group with the correct signature
		p.first->functions.push_back(f);

		return p.second;
	}

	template<uint64_t N>
	group_vector<N> generate_groups_fast(const bf_vector<N> &funcs)
	{
		//benchmark bm("generate_groups_fast", benchmark::report::end);

		group_set<N> gs;

		for (const auto &f : funcs)
			add_function_to_groupset(gs, f);

		return group_vector<N>(gs.begin(), gs.end());
	}

	template<uint64_t N, typename SetCont>
	group_set<N> merge_groupsets(const std::vector<SetCont> &cont) {

		group_set<N> gs;

		for (auto &gsl : cont)
		{
			for (auto &gr : gsl)
			{
				auto p = gs.insert(gr);

				if (!p.second)
					std::copy(std::begin(gr.functions), std::end(gr.functions), std::back_inserter(p.first->functions));
			}
		}

		return gs;
	}

	template<uint64_t N>
	group_vector<N> generate_groups_parallel(const bf_vector<N> &funcs)
	{
		benchmark bm("generate_groups_parallel", benchmark::report::end);

		std::vector<group_set<N>> gsv(omp_get_max_threads());

#pragma omp parallel for
		for(int fIndex = 0; fIndex < funcs.size(); fIndex++)
			add_function_to_groupset(gsv[omp_get_thread_num()], funcs[fIndex]);

		auto gs = merge_groupsets<N>(gsv);

		return group_vector<N>(gs.begin(), gs.end());
	}

}

namespace std {

	template<uint64_t N>
	struct hash<ae::group<N>>
	{
		size_t operator()(const ae::group<N> &gr) const {

			//! Note: functions do not contribute to hash because we want the 
			//! groups to have the same hash only if they have the same signature
			return hash<ae::affine_invariants>()(gr.ai);
		}
	};
}