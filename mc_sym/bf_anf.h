#pragma once

#include "bf_base.h"
#include "bf_monomial.h"

namespace bfl {

	template<uint64_t N>
	class bf_anf : public bf_base<N> {

	public:

		friend class bf_tt<N>;

		using table_t = typename bf_base<N>::table_t;

		bf_anf() {}

		bf_anf(init in) : bf_base<N>(in) {}

		bf_anf(const bf_anf<N> &rhs) {
			this->table_ = rhs.table_;
		}

		bf_anf(const table_t &tab) {
			this->table_ = tab;
		}

		bf_anf(const bf_tt<N> &rhs) {
			this->table_ = rhs.table_;
			this->table_.mobius_transform();
		}

		bf_anf(const std::string &anf, ordering ord = ordering::rtl1) {
			assign(anf, ord);
		}

		bf_anf<N>& operator=(const bf_tt<N> &rhs) {
			this->table_ = rhs.table_;
			this->table_.mobius_transform();
			return *this;
		}

		bf_anf<N>& operator=(const std::string &anf) {
			assign(anf);
			return *this;
		}

		bf_tt<N> tt() const {
			table_t newtab(this->table_);
			newtab.mobius_transform();
			return bf_tt<N>(newtab);
		}

		std::string str(ordering ord = ordering::rtl1) const {

			if (this->table_.is_zero())
				return { "0" };

			std::string s;

			for (uint64_t i = 0; i < table_t::TableSize; i++)
			{
				if (this->table_.get(i))
				{
					if (!s.empty())
						s.append("+");

					s += monomial<N>::str(i, ord);
				}
			}

			return s;
		}

		uint64_t monomial_count() const {
			return this->table_.weight();
		}

		uint64_t degree() const {

			// note: constant functions are assigned algebraic degree 1
			uint64_t deg{ 1 };

			for (uint64_t i = 1; i < table_t::TableSize; i++)
			{
				if (this->table_.get(i))
				{
					auto mondeg = popcnt(i);

					if (mondeg > deg)
						deg = mondeg;
				}
			}

			return deg;
		}

		bool has_terms_with_degree(uint64_t d) const {

			for (uint64_t i = 0; i < table_t::TableSize; i++)
			{
				if (this->table_.get(i) && popcnt(i) == d)
					return true;
			}

			return false;
		}

		uint64_t dependent_variable_count() const {

			auto vars = dependent_variable_mask();

			return popcnt(vars);
		}

		uint64_t dependent_variable_mask() const {

			uint64_t mask{ 0 };

			for (uint64_t i = 0; i < table_t::TableSize; i++)
				if (this->table_.get(i))
					mask |= i;

			return mask;
		}

		bool is_complete() const {

			return (dependent_variable_count() == N);
		}

		void add_monomial(uint64_t index) {
			this->table_.flip(index);
		}

		void add_monomial(const std::string &mon) {
			add_monomial(monomial<N>::index(mon));
		}

		void assign(const std::string &anf, ordering ord = ordering::rtl1) {

			this->table_.clear();

			std::size_t start = 0, end;

			do
			{
				end = anf.find_first_of('+', start);

				std::string monstr = anf.substr(start, end == std::string::npos ? (anf.size() - start) : (end - start));

				if (monstr != "0")
					this->table_.set(monomial<N>::index(monstr, ord));

				if (end == std::string::npos)
					break;

				start = end + 1;

			} while (start < anf.size());
		}

		void remove_affine_terms() {
			*this &= this->affine_term_mask_;
		}

		bf_pair<N> factor(uint64_t mon) const {

			bf_anf<N> f1(init::zero), f2(init::zero);

			for (uint64_t i = 0; i < table_t::TableSize; i++)
				if (this->table_.get(i))
				{
					if ((i & mon) == mon)
						f1.set(i & (~mon));
					else
						f2.set(i);
				}

			return bf_pair<N>{f1.tt(), f2.tt()};
		}

		uint64_t twin_pair_count() const {

			std::array<uint64_t, N> masks;

			for (uint64_t i = 0; i < N; i++)
				masks[i] = (~(1 << i)) & ((1 << N) - 1);

			for (uint64_t x = 0; x < table_t::TableSize; x++)
			{
				if (this->table_.get(x))
				{
					for (uint64_t j = 0; j < N; j++)
						if ((x >> j) & 1)
							masks[j] &= x;
				}
			}

			uint64_t twinCount{ 0 };

			for (uint64_t i = 0; i < N - 1; i++)
			{
				if (masks[i])
				{
					for (uint64_t j = i + 1; j < N; j++)
					{
						if (((masks[i] >> j) & 1) && ((masks[j] >> i) & 1))
							twinCount++;
					}
				}
			}

			return twinCount;
		}

		bool has_twins() const {

			return (twin_pair_count() > 0);
		}

	private:

		static bf_anf<N> get_affine_term_mask() {

			bf_anf<N> mask(init::zero);

			mask.set(0);

			for (uint64_t i = 0; i < N; i++)
				mask.set(1ULL << i);

			mask.table().complement();

			return mask;
		}

		static bf_anf<N> affine_term_mask_;
	};

	template<uint64_t N>
	bf_anf<N> bf_anf<N>::affine_term_mask_ = bf_anf<N>::get_affine_term_mask();

}
