#pragma once

#include "bf.h"
#include "ae.h"
#include "topology.h"
#include "topology_slp.h"
#include "logger.h"
#include "ae_distinguisher.h"

using namespace bfl;
using namespace ae;

template<uint64_t N>
struct circuit_data
{
	bf_tt<N> frep;
	Topology topology;
	bf_vector<N> inputs;
	bf_tt<N> ffinal;
	uint64_t outputMask;

	bf_tt<N> compute() const {

		TopologySLP<N> slp(topology);

		auto f = slp.evaluateMask(inputs, outputMask);

		f ^= ffinal;

		return f;
	}

	void print(bfl::mode m = bfl::mode::anf) const {

		LOG_ASSERT_EQ(inputs.size(), 2 * topology.gateCount());

		std::cout << str(frep, m) << '\n';
		std::cout << topology.str() << '\n';
		std::cout << outputMask << '\n';
		for (int k = 0; k < inputs.size(); k++)
			std::cout << str(inputs[k], m) << '\n';
		std::cout << str(ffinal, m) << '\n';
		std::cout << '\n';
	}
};

template<uint64_t N>
bool operator<(const circuit_data<N> &sd1, const circuit_data<N> &sd2)
{
	return sd1.frep < sd2.frep;
}

template<uint64_t N>
std::vector<circuit_data<N>> load_circuit_data(const std::string &fileName)
{
	std::vector<circuit_data<N>> circdata;

	LOG_INFO << "loading circuit data '" << fileName << "'" << LOG_ENDL;

	std::ifstream ifile(fileName);

	if (ifile.is_open())
	{
		std::string line;

		while (std::getline(ifile, line))
		{
			if (line.empty())
				continue;

			circuit_data<N> sd;

			sd.frep.assign(line);

			std::getline(ifile, line);
			sd.topology.set(line);

			std::getline(ifile, line);
			sd.outputMask = std::strtoul(line.data(), nullptr, 10);

			for (int i = 0; i < 2 * sd.topology.gateCount(); i++)
			{
				std::getline(ifile, line);
				sd.inputs.push_back(bf_tt<N>(line));
			}

			std::getline(ifile, line);
			sd.ffinal.assign(line);

			circdata.push_back(sd);
		}
	}

	std::sort(std::begin(circdata), std::end(circdata));

	LOG_INFO << circdata.size() << " circuits loaded" << LOG_ENDL;

	LOG_ERROR_IF(circdata.empty()) << "no circuits loaded" << LOG_ENDL;

	return circdata;
}

template<uint64_t N>
std::vector<circuit_data<N>> load_circuit_data() {

	static_assert(N >= 2 && N <= 6, "");

	return load_circuit_data<N>("../data/circuits/n" + std::to_string(N) + "_circ.txt");
}

template<uint64_t N>
circuit_data<N> build_slp(const bf_tt<N> &f, class_distinguisher_base<N> &cs, const std::vector<circuit_data<N>> &sd)
{
	auto rep = cs.get_representative(f);

	circuit_data<N> s;

	s.frep = rep;

	auto p = std::lower_bound(std::begin(sd), std::end(sd), s);

	if (p == std::end(sd) || p->frep != rep)
	{
		LOG_ERROR << "could not find the slp" << LOG_ENDL;
		LOG_ERROR << "f   : " << f << LOG_ENDL;
		LOG_ERROR << "rep : " << rep << LOG_ENDL;
	}
	else
	{
		s = *p;

		affine_mapping<N> am;

		if (is_affine_equivalent(rep, f, am))
		{
			s.frep = f;

			for (auto &g : s.inputs)
				g = transform(g, am.A, am.a);

			s.ffinal = transform(s.ffinal, am.A, am.a);

			s.ffinal ^= am.b;
		}
	}

	return s;
}
