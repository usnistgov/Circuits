#pragma once

#include <array>
#include <vector>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "constant_time.h"
#include "utils.h"


namespace bfl {


	using spectrum_value_t		= int32_t;
	using spectrum_size_t		= std::make_unsigned<spectrum_value_t>::type;
	using spectrum_distribution = std::vector<std::pair<spectrum_value_t, spectrum_size_t>>;

	// TODO: What is the maximum number of distinct walsh spectrum entries for n > 6 ?
	// TODO: Make sure the maximum array size is large enough to store the distribution
	//using spectrum_distribution = dyn_array<std::pair<spectrum_value_t, spectrum_size_t>, 10>;

	bool operator==(const spectrum_distribution &d1, const spectrum_distribution &d2);

	bool operator<(const spectrum_distribution &d1, const spectrum_distribution &d2);

	int compare(const spectrum_distribution &d1, const spectrum_distribution &d2);

	template<uint64_t N>
	class spectrum
	{
	public:

		// auto-correlation calculation squares the walsh entries, doubling the size required to store the result
		static_assert((sizeof(spectrum_value_t) * 8) >= 2 * (N + 2), "insufficient size to store spectrum values");
		static_assert((sizeof(spectrum_size_t) * 8) >= (N + 1), "insufficient size to store spectrum counts");

		static const uint64_t TableSize = power_of_2(N);

		enum class zeros {exclude , include};

		spectrum() {}

		~spectrum() {}

		spectrum_distribution absolute_distribution(zeros z = zeros::include) const {

			//return calculate_distribution([](spectrum_value_t val) { return std::abs(val); }, z);

			spectrum_distribution sd;

			// Optimization hint to avoid multiple reallocations
			sd.reserve(N + (N / 2));

			for (uint64_t i = 0; i < TableSize; i++)
			{
				// Ignoring zero vales in the distribution improves performance
				// There cannot be two spectrums that only differ in the number of zeros
				if (z == zeros::exclude && _spectrum[i] == 0)
					continue;

				auto value = ct::abs(_spectrum[i]);

				auto p = std::find_if(std::begin(sd), std::end(sd), [=](std::pair<spectrum_value_t, spectrum_size_t> &x) {
					return x.first == value;
				});

				if (p != std::end(sd))
					p->second++;
				else
					sd.push_back(std::make_pair(value, 1));
			}

			// spectrum_distribution comparison is based on the entries being sorted
			std::sort(std::begin(sd), std::end(sd));

			return sd;
		}

		spectrum_value_t operator[](uint64_t index) const {
			return _spectrum[index];
		}

		spectrum_value_t& operator[](uint64_t index) {
			return _spectrum[index];
		}

		void print() const {
			for (auto w : _spectrum)
				std::cout << std::setw(4) << w << " ";
			std::cout << std::endl;
		}

		uint64_t zero_count() const {

			return std::count_if(_spectrum.begin(), _spectrum.end(), [](spectrum_value_t x) {
				return x == 0;
			});
		}

		uint64_t max_value_count() const {

			return std::count_if(_spectrum.begin(), _spectrum.end(), [](spectrum_value_t x) {
				return std::abs(x) == power_of_2(N);
			});
		}

		bool operator == (const spectrum<N> &rhs) const {
			return (_spectrum == rhs._spectrum);
		}

	protected:

		template<class Fn>
		[[deprecated]] spectrum_distribution calculate_distribution(Fn f, zeros z) const {

			spectrum_distribution sd;

			// Optimization hint to avoid multiple reallocations
			sd.reserve(N + 2);

			for (uint64_t i = 0; i < TableSize; i++)
			{
				// Ignoring zero vales in the distribution improves performance
				// There cannot be two spectrums that only differ in the number of zeros
				if (z == zeros::exclude && _spectrum[i] == 0)
					continue;

				auto value = f(_spectrum[i]);

				auto p = std::find_if(std::begin(sd), std::end(sd), [=](std::pair<spectrum_value_t, spectrum_size_t> &x) {
					return x.first == value;
				});

				if (p != std::end(sd))
					p->second++;
				else
					sd.push_back(std::make_pair(value, 1));
			}

			// comparison operator assumes that the entries are sorted
			std::sort(std::begin(sd), std::end(sd));

			return sd;
		}

		std::array<spectrum_value_t, TableSize> _spectrum;
	};

	template<uint64_t N, uint64_t BlockCount>
	struct fast_walsh_transform {

		static const uint64_t BlockSize = (spectrum<N>::TableSize / BlockCount);

		static void apply(spectrum<N> &arr) {

			for (uint64_t i = 0; i < BlockCount; i++)
			{
				for (uint64_t j = 0; j < BlockSize / 2; j++)
				{
					auto sum  = arr[i * BlockSize + j] + arr[i * BlockSize + BlockSize / 2 + j];
					auto diff = arr[i * BlockSize + j] - arr[i * BlockSize + BlockSize / 2 + j];

					arr[i * BlockSize + j] = sum;
					arr[i * BlockSize + BlockSize / 2 + j] = diff;
				}
			}

			fast_walsh_transform<N, BlockCount / 2>::apply(arr);
		}
	};

	template<uint64_t N>
	struct fast_walsh_transform<N, 0> {
		static void apply(spectrum<N>&) {
		}
	};

	template<uint64_t N>
	void compute_walsh_spectrums()
	{
		const uint64_t TableSize = power_of_2(N);

		using spec_t = std::array<int32_t, TableSize>;

		std::vector<spec_t> spec;


		for (uint64_t x = 0; x < power_of_2(TableSize); x++)
		{
			bf_tt<N> f(init::zero);

			for (uint64_t i = 0; i < TableSize; i++)
				f.set(i, get_bit(x, i));

			auto walsh = f.walsh_spectrum();

			spec.push_back(spec_t());

			for (uint64_t i = 0; i < walsh.TableSize; i++)
				spec[x][i] = walsh[i];
		}


		for (auto &s : spec)
		{
			for (auto x : s)
				std::cout << std::setw(2) << x << " ";
			std::cout << std::endl;
		}
	}

	template<uint64_t N>
	std::vector<spectrum<N>> load_walsh_spectrums(const std::string &filename) {

		const uint64_t TableSize = power_of_2(N);

		std::vector<spectrum<N>> spec;

		std::ifstream ifile(filename);

		if (ifile.is_open())
		{
			spec.reserve(power_of_2(TableSize));

			for (uint64_t i = 0; i < power_of_2(TableSize); i++)
			{
				spectrum<N> s;

				for (uint64_t j = 0; j < TableSize; j++)
					ifile >> s[j];

				spec.push_back(s);
			}
		}

		return spec;
	}

}

namespace std {

	template<>
	struct hash<bfl::spectrum_distribution> {

		size_t operator()(const bfl::spectrum_distribution &sd) const {

			size_t digest{ 0 };

			for (const auto &e : sd)
				digest ^= std::hash<bfl::spectrum_value_t>()(e.first) ^
				std::hash<bfl::spectrum_size_t>()(e.second);

			return digest;
		}
	};
}
