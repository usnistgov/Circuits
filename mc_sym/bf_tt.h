#pragma once

#include <cstdint>
#include <mutex>
#include "bf_monomial.h"
#include "bf_base.h"
#include "spectrum.h"

namespace bfl {

	template<uint64_t N>
	class bf_tt : public bf_base<N> {

	public:

		friend class bf_anf<N>;

		using table_t = typename bf_base<N>::table_t;

		static const uint64_t MaxNonlinearity = (power_of_2(N - 1) - power_of_2((N / 2) - 1));

		bf_tt() {}

		explicit bf_tt(init in) : bf_base<N>(in) {}

		bf_tt(const bf_tt<N> &rhs) {
			this->table_ = rhs.table_;
		}

		template<uint64_t M>
		bf_tt(const bf_tt<M> &rhs) {
			static_assert(N > M, "");
			// TODO: make this more efficient by copying the truth table
			assign(rhs.anf().str());
		}

		bf_tt(const bf_anf<N> &rhs) {
			this->table_ = rhs.table_;
			this->table_.mobius_transform();
		}

		bf_tt(const table_t &tab) {
			this->table_ = tab;
		}

		bf_tt(const std::string &anf, ordering ord = ordering::rtl1) {
			assign(anf, ord);
		}

		bf_tt<N>& operator=(const bf_anf<N> &rhs) {
			this->table_ = rhs.table_;
			this->table_.mobius_transform();
			return *this;
		}

		bf_tt<N>& operator=(const std::string &anf) {
			assign(anf);
			return *this;
		}

		bf_anf<N> anf() const {
			table_t newtab(this->table_);
			newtab.mobius_transform();
			return bf_anf<N>(newtab);
		}

		uint64_t weight() const {
			return this->table_.weight();
		}

		void assign(const std::string &anf, ordering ord = ordering::rtl1) {
			bf_anf<N> fanf(anf, ord);
			*this = fanf;
		}

		spectrum<N> walsh_spectrum() const {

			spectrum<N> s = polarity_truth_table();

			fast_walsh_transform<N, bf_tt<N>::TableSize / 2>::apply(s);

			return s;
		}

		uint64_t nonlinearity() const {

			auto walsh = walsh_spectrum();

			return nonlinearity(walsh);
		}

		uint64_t nonlinearity(const spectrum<N> &walsh) const {

			uint64_t maxWalsh{ 0 };

			for (uint64_t i = 0; i < spectrum<N>::TableSize; i++)
			{
				uint64_t value = std::abs(walsh[i]);

				if (value > maxWalsh)
					maxWalsh = value;
			}

			return power_of_2(N - 1) - (maxWalsh / 2);
		}

		bool is_bent() const {

			if (is_odd(N)) 
				return false;

			return nonlinearity() == MaxNonlinearity;
		}

		bool is_balanced() const {

			return weight() == power_of_2(N - 1);
		}

		spectrum<N> autocorrelation_spectrum() const {

			auto walsh = walsh_spectrum();

			return autocorrelation_spectrum(walsh);
		}

		spectrum<N> autocorrelation_spectrum(const spectrum<N> &walsh) const {
			
			spectrum<N> acs;

			for (uint64_t i = 0; i < spectrum<N>::TableSize; i++)
				acs[i] = walsh[i] * walsh[i];

			fast_walsh_transform<N, spectrum<N>::TableSize / 2>::apply(acs);

			for (uint64_t i = 0; i < spectrum<N>::TableSize; i++)
				acs[i] >>= N;

			return acs;
		}

		void remove_affine_terms() {

			bf_anf<N> tmp(*this);

			tmp.remove_affine_terms();

			*this = tmp;
		}

		void flip(uint64_t index) {
			this->table_.flip(index);
		}

		static const bf_vector<N>& get_linear_functions() {

			static std::mutex m;
			std::lock_guard<std::mutex> lock(m);

			if (linear_functions_.empty())
				linear_functions_ = generate_linear_functions(false);

			return linear_functions_;
		}

		static const bf_vector<N>& get_affine_functions() {

			static std::mutex m;
			std::lock_guard<std::mutex> lock(m);

			if (affine_functions_.empty())
				affine_functions_ = generate_linear_functions(true);

			return affine_functions_;
		}

		static const bf_vector<N>& get_coordinate_functions() {

			if (coordinate_functions_.empty())
				coordinate_functions_ = generate_coordinate_functions();

			return coordinate_functions_;
		}

		void remove_variables(uint64_t vars) {

			auto f = anf();

			for (uint64_t x = 0; x < table_t::TableSize; x++)
			{
				if (f.get(x) && (x & vars))
				{
					f.set(x, 0);
					f.flip((x & (~vars)) & (power_of_2(N) - 1));
				}
			}

			*this = f;
		}

	private:

		static bf_vector<N> generate_linear_functions(bool affine = false) {

			std::vector<bf_tt<N>> lin;

			for (uint64_t mask = 0; mask < power_of_2(N); mask++)
			{
				bf_anf<N> f(init::zero);

				for (uint64_t m = 0; m < N; m++)
					if ((mask >> m) & 1)
						f.add_monomial(1ULL << m);

				lin.push_back(f.tt());

				if (affine)
				{
					f.set(0);
					lin.push_back(f.tt());
				}
			}

			return lin;
		}

		static bf_vector<N> generate_coordinate_functions() {

			std::vector<bf_tt<N>> coords;

			for (uint64_t m = 0; m < N; m++)
			{
				bf_anf<N> f(init::zero);

				f.add_monomial(1ULL << m);

				coords.push_back(f.tt());
			}


			return coords;
		}

		spectrum<N> polarity_truth_table() const {

			spectrum<N> s;

			for (uint64_t a = 0; a < bf_base<N>::TableSize; a++)
				s[a] = static_cast<spectrum_value_t>(1 - 2 * this->table_.get(a));

			return s;
		}

		// TODO: Obsolete : replaced by template based fast_walsh_transform()
		static void walsh_transform(spectrum<N> &s) {

			uint64_t blockLength = 1;

			for (uint64_t a = 0; a < N; a++)
			{
				blockLength <<= 1;

				for (uint64_t b = 0; b < spectrum<N>::TableSize / blockLength; b++)
				{
					for (uint64_t c = 0; c < blockLength / 2; c++)
					{
						auto sum	= s[b * blockLength + c] + s[b * blockLength + blockLength / 2 + c];
						auto diff	= s[b * blockLength + c] - s[b * blockLength + blockLength / 2 + c];

						s[b * blockLength + c] = sum;
						s[b * blockLength + blockLength / 2 + c] = diff;
					}
				}
			}
		}

		// TODO: generate those functions on the fly as opposed to caching
		//		 Or, use an iterator mechanism to traverse all functions
		static bf_vector<N> linear_functions_;
		static bf_vector<N> affine_functions_;
		static bf_vector<N> coordinate_functions_;
	};

	template<uint64_t N>
	bf_vector<N> bf_tt<N>::linear_functions_;

	template<uint64_t N>
	bf_vector<N> bf_tt<N>::affine_functions_;

	template<uint64_t N>
	bf_vector<N> bf_tt<N>::coordinate_functions_;

	template<uint64_t N>
	bf_tt<N> operator&(const bf_tt<N>& lhs, const bf_tt<N> &rhs) {

		// TODO : check which implementation is faster for tables of size 1
		//bf_tt<N> res(lhs);
		//res &= rhs;
		//return res;

		bf_tt<N> res;
		auto size = bf_tt<N>::table_t::WordCount;
		auto &tsrc1 = lhs.table();
		auto &tsrc2 = rhs.table();
		auto &tdest = res.table();
		for (uint64_t i = 0; i < size; i++)
			tdest[i] = tsrc1[i] & tsrc2[i];
		return res;
	}

	template<uint64_t N>
	bf_tt<N> operator^(const bf_tt<N>& lhs, const bf_tt<N> &rhs) {
		// TODO : check which implementation is faster for tables of size 1
		//bf_tt<N> res(lhs);
		//res ^= rhs;
		//return res;

		bf_tt<N> res;
		auto size = bf_tt<N>::table_t::WordCount;
		auto &tsrc1 = lhs.table();
		auto &tsrc2 = rhs.table();
		auto &tdest = res.table();
		for (uint64_t i = 0; i < size; i++)
			tdest[i] = tsrc1[i] ^ tsrc2[i];
		return res;
	}
}

namespace std {

	template<uint64_t N>
	struct hash<bfl::bf_tt<N>> {

		size_t operator()(const bfl::bf_tt<N> &f) const {

			size_t digest{ 0 };

			// TODO : improve the hash function
			for (uint64_t i = 0; i < bfl::bf_tt<N>::table_t::WordCount; i++)
				digest ^= (hash<bfl::bf_tt<N>::table_t::word_t>()(f.table().word(i)) << 1);

			return digest;
		}
	};

	template<uint64_t N>
	struct hash<bfl::bf_vector<N>> {

		size_t operator()(const bfl::bf_vector<N> &fvec) const {

			size_t digest{ 0 };

			// TODO : improve the hash function
			for (const auto &f : fvec)
				digest ^= (hash<bfl::bf_tt<N>>()(f) << 1);

			return digest;
		}
	};
}