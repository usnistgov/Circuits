#pragma once

#include <array>
#include <vector>
#include "bf.h"
#include "ae.h"
#include "topology.h"
#include "topology_slp.h"
#include "slp.h"

using namespace bfl;
using namespace ae;

template<uint64_t N, uint64_t K>
struct topology_eval_data
{
	bf_tt<N> fout;
	bf_tt<N> frep;
	std::array<bf_tt<N>, 2 * K> inputs;
};


template<uint64_t N>
uint64_t find_mc(const bf_tt<N> &f, const class_distinguisher_base<N> &cs, const std::vector<circuit_data<N>> &sd)
{
	uint64_t mc{ 0 };

	auto rep = cs.get_representative(f);

	circuit_data<N> s;

	s.frep = rep;

	auto p = std::lower_bound(std::begin(sd), std::end(sd), s);

	if (p->frep != rep)
		std::cout << "could not find the slp" << std::endl;
	else
		mc = p->topology.gateCount();

	return mc;
}


template<uint64_t N>
std::map<affine_invariants, uint64_t> compute_mc_map(const class_distinguisher_base<N> &cs, const std::vector<circuit_data<N>> &sd)
{
	std::map<affine_invariants, uint64_t> mcmap;

	auto reps = cs.representatives();

	auto groups = cs.groups();

	uint64_t countSame{ 0 }, countDifferent{ 0 };

	for (auto &gr : groups)
	{
		std::set<uint64_t> mcset;

		for (auto f : gr.functions)
		{
			auto mc = find_mc(f, cs, sd);
			mcset.insert(mc);
		}

		if (mcset.size() > 1)
		{
			for (auto x : mcset)
				LOG_TRACE << x << ' ';
			LOG_TRACE << LOG_ENDL;

			countDifferent++;
		}
		else
		{
			mcmap[gr.ai] = *mcset.begin();
			countSame++;
		}
	}

	LOG_DEBUG << "same : " << countSame << LOG_ENDL;
	LOG_DEBUG << "different : " << countDifferent << LOG_ENDL;

	return mcmap;
}

template<uint64_t N>
uint64_t find_mc(const bf_tt<N> f, const class_distinguisher_base<N> &cs, const std::vector<circuit_data<N>> &sd, const std::map<affine_invariants, uint64_t> &mcmap)
{
	auto ai = calculate_affine_invariants(f);

	auto p = mcmap.find(ai);

	if (p != mcmap.end())
		return p->second;
	else
		return find_mc(f, cs, sd);
}


template<uint64_t N>
std::map<bf_tt<N>, uint64_t> load_mc_data(const std::string &filename)
{
	std::map<bf_tt<N>, uint64_t> ret;

	std::ifstream ifile(filename);

	if (!ifile.is_open())
		return ret;

	while (!ifile.eof())
	{
		std::string line;

		std::getline(ifile, line);

		std::istringstream stream(line);

		if (!line.empty())
		{
			uint64_t mc;
			std::string anf;

			stream >> mc >> anf;

			ret[bf_tt<N>(anf)] = mc;
		}
	}

	return ret;
}
