#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <algorithm>
#include <fstream>

#include "bf_base.h"
#include "bf_tt.h"
#include "bf_anf.h"
#include "ae_transformation.h"
#include "matrix.h"
#include "benchmark.h"

namespace bfl
{
	enum class sort { no, yes };

	enum class mode {anf, hex, bin};

	template<uint64_t N>
	bf_vector<N> load_functions(const std::string &fileName, mode m = mode::anf, sort s = sort::no)
	{
		std::ifstream file(fileName);

		bf_vector<N> functions;

		if (file)
		{
			progress_bar pb("loading functions from " + fileName, file_size(file));

			std::string line;

			while (!file.eof())
			{
				if (file >> line)
				{
					bf_tt<N> f;

					if (m == mode::anf)
						f.assign(line);
					else if (m == mode::hex)
						f.table().set_hex(line);

					functions.push_back(f);
				}

				pb.update(line.size() + 2);
			}

			if (s == sort::yes)
				std::sort(std::begin(functions), std::end(functions));
		}

		LOG_INFO << "loaded " << functions.size() << " functions from '" << fileName << LOG_ENDL;

		return functions;
	}

	template<uint64_t N, uint64_t K>
	bf_tuple_vector<N, K> load_function_tuples(const std::string &fileName)
	{
		bf_tuple_vector<N, K> result;

		auto bfvec = load_functions<N>(fileName);

		for (int i = 0; i < bfvec.size() / K; i++)
		{
			bf_tuple<N, K> arr;

			for (int j = 0; j < K; j++)
				arr[j] = bfvec[K * i + j];

			result.push_back(arr);
		}

		return result;
	}

	template<uint64_t N, typename Iter>
	bool write_functions(Iter first, Iter last, const std::string &filename, mode m = mode::anf)
	{
		std::ofstream file(filename);

		if (!file.is_open())
			return false;

		while(first != last)
		{
			if (m == mode::anf)
				file << (*first).anf().str() << std::endl;
			else if (m == mode::hex)
				file << (*first).table().hex_str() << std::endl;

			++first;
		}

		return true;
	}

	template<uint64_t N>
	bool write_functions(const bf_vector<N> &vec, const std::string &filename, mode m = mode::anf)
	{
		return write_functions(vec.begin(), vec.end(), filename, m);
	}

	template<uint64_t N, typename T>
	bf_tt<N> flip_entries(const bf_tt<N> &f, const std::vector<T> &entries)
	{
		bf_tt<N> g(f);

		for (auto x : entries)
			g.flip(x);

		return g;
	}

	template<uint64_t N, class FwIter>
	void remove_affine_terms(FwIter first, FwIter last)
	{
		while (first != last)
		{
			(*first).remove_affine_terms();
			++first;
		}
	}

	template<uint64_t N>
	std::array<uint64_t, N + 1> degree_distribution(const bf_tt<N> &f)
	{
		std::array<uint64_t, N + 1> dist;

		std::fill(std::begin(dist), std::end(dist), 0);

		auto g = f.anf();

		for (uint64_t x = 0; x < bf_tt<N>::TableSize; x++)
			if (g.get(x))
				dist[popcnt(x)]++;

		return dist;
	}

	template<uint64_t N>
	uint64_t dimension(const bf_tt<N> &f)
	{
		auto max_values = f.autocorrelation_spectrum().max_value_count();

		LOG_ASSERT(popcnt(max_values) == 1);

		return N - msb_index(max_values);
	}

	template<uint64_t N>
	bool is_reducible(const bf_tt<N> &f) {

		return (dimension(f) < N);
	}

	template<uint64_t M, uint64_t N>
	bf_tt<M> extend(const bf_tt<N> &f) {

		return bf_tt<M>(f);
	}

	template<uint64_t M, uint64_t N>
	bf_tt<M> shrink(const bf_tt<N> &f)
	{
		static_assert(M < N);

		bf_tt<M> g(init::zero);

		for (uint64_t x = 0; x < power_of_2(M); x++)
			g.set(x, f.get(x + power_of_2(N) - power_of_2(M)));

		return g;
	}

	template<uint64_t N>
	bf_tt<N> shrink_dimension(const bf_anf<N> &f)
	{
		auto mask = f.dependent_variable_mask();

		auto dim = popcnt(mask);

		if (dim == N)
			return f.tt();

		uint64_t unused_var = N - 1;

		matrix<N> m;

		for (uint64_t i = 0; i < dim; i++)
		{
			if ((mask >> i) & 1)
			{
				// if xi is used do not remap it
				m[i] = (1ULL << i);
			}
			else
			{
				// find the largest used variable
				while (!((mask >> unused_var) & 1))
					unused_var--;

				// swap the variables i and unused var (i is unused)
				m[i] = (1ULL << unused_var);
				m[unused_var] = (1ULL << i);

				unused_var--;
			}
		}

		return transform(f.tt(), m);
	}

	template<uint64_t N>
	bf_tt<N> shrink_dimension(const bf_tt<N> &f)
	{
		return shrink_dimension(f.anf());
	}

	template<uint64_t N>
	bf_tt<N - 1> reduce_dimension(const bf_tt<N> &f)
	{
		bf_anf<N> mask(init::zero);

		mask.set(power_of_2(N - 1));

		auto p = decompose(f, mask.tt());

		return p.second;
	}


	template<uint64_t N>
	uint64_t distance(const bf_tt<N> &f, const bf_tt<N> &g)
	{
		bf_tt<N> h(f);
		h ^= g;
		return h.weight();
	}

	template<uint64_t N>
	bf_tt<N> linear_function(uint64_t mask)
	{
		assert(mask < power_of_2(N));

		bf_anf<N> f(init::zero);

		for (uint64_t i = 0; i < N; i++)
			if ((mask >> i) & 1)
				f.add_monomial(1ULL << i);

		return f.tt();
	}

	template<uint64_t N>
	bool is_affine(const bf_tt<N> &f) {

		const auto &affines = bf_tt<N>::get_affine_functions();

		return std::find(std::begin(affines), std::end(affines), f) != std::end(affines);
	}

	template<uint64_t N>
	std::pair<bf_tt<N - 1>, bf_tt<N - 1>> decompose(const bf_tt<N> &f, const bf_tt<N> &fmask)
	{
		assert(fmask.is_balanced());

		bf_tt<N - 1> f0(init::zero);
		bf_tt<N - 1> f1(init::zero);

		uint64_t f0Index = 0;
		uint64_t f1Index = 0;

		for (uint64_t x = 0; x < power_of_2(N); x++)
		{
			if (fmask.get(x))
				f1.set(f1Index++, f.get(x));
			else
				f0.set(f0Index++, f.get(x));
		}

		return std::make_pair(f0, f1);
	}

	template<int N, int M, class ElemType = uint8_t>
	bf_vector<N> compute_vectorial_functions(const ElemType *sbox)
	{
		bf_vector<N> vec(M);

		for (auto &f : vec)
			f.clear();

		for (int x = 0; x < power_of_2(N); x++)
		{
			for (int i = 0; i < M; i++)
			{
				vec[i].set(x, (sbox[x] >> i) & 1);
			}
		}

		return vec;
	}

	template<int N, int M>
	std::vector<uint64_t> compute_lookup_table(const bf_vector<N> &vec)
	{
		std::vector<uint64_t> table;

		table.reserve(power_of_2(N));

		for (uint64_t x = 0; x < power_of_2(N); x++)
		{
			uint64_t val{ 0 };

			for (uint64_t m = 0; m < M; m++)
				val ^= vec[m].get(x) << m;

			table.push_back(val);
		}

		return table;
	}

	template<uint64_t N>
	std::array<uint8_t, power_of_2(N)> to_bytes(const bf_tt<N> &f)
	{
		std::array<uint8_t, power_of_2(N)> ret;

		for (uint64_t i = 0; i < power_of_2(N); i++)
			ret[i] = static_cast<uint8_t>(f.get(i));

		return ret;
	}

	template<uint64_t N>
	std::ostream& operator<<(std::ostream& os, const bf_tt<N> &f)
	{
		return os << f.anf().str();
	}

	template<uint64_t N>
	std::ostream& operator<<(std::ostream& os, const bf_anf<N> &f)
	{
		return os << f.str();
	}

	template<uint64_t N>
	std::string str(const bf_tt<N> &f, bfl::mode m)
	{
		if (m == bfl::mode::anf)
			return f.anf().str();
		else if (m == bfl::mode::hex)
			return f.table().hex_str();
		else if (m == bfl::mode::bin)
			return f.table().bin_str();

		return {};
	}

	template<typename Cont>
	void print_bf(	const Cont& cont,
					bfl::mode m = bfl::mode::anf,
					bool bfinalNewLine = true,
					std::ostream& os = std::cout)
	{
		for (const auto &f : cont)
			os << str(f, m) << std::endl;

		if (bfinalNewLine)
			os << std::endl;
	}


}

