#pragma once

#include "bf_tables.h"
#include <random>
#include <set>

namespace bfl {

	//! Boolean function representation types
	enum class rep { tt, anf };

	/// Initialization flags
	enum class init { none, zero, random };

	template<uint64_t N>
	class bf_base {

	public:

		using table_t = table_u64<N>;

		static const uint64_t Variables = N;

		static const uint64_t TableSize = table_t::TableSize;

		void clear() {
			table_.clear();
		}

		void random() {

			table_.clear();

			for (uint64_t i = 0; i < TableSize; i++)
				table_.set(i, parity(rd_()));
		}

		bool is_zero() const {
			return table_.is_zero();
		}

		table_t& table() {
			return table_;
		}

		const table_t& table() const {
			return table_;
		}

		void complement() {
			this->table_.complement();
		}

		uint64_t table_size() const {
			return table_t::TableSize;
		}

		void set(uint64_t index, uint64_t val = 1) {
			table_.set(index, val);
		}

		uint64_t get(uint64_t index) const {
			return table_.get(index);
		}

		bool operator==(const bf_base<N> &rhs) const {
			return (this->table_ == rhs.table_);
		}

		bool operator!=(const bf_base<N> &rhs) const {
			return (this->table_ != rhs.table_);
		}

		bool operator<(const bf_base<N> &rhs) const {
			return this->table_ < rhs.table_;
		}

		void operator^=(const bf_base<N> &rhs) {
			this->table_ ^= rhs.table_;
		}

		void operator&=(const bf_base<N> &rhs) {
			this->table_ &= rhs.table_;
		}

		void operator|=(const bf_base<N> &rhs) {
			this->table_ |= rhs.table_;
		}

		void flip(uint64_t index) {
			this->table_.flip(index);
		}

	protected:

		bf_base() {
			this->table_.clear();
		}

		bf_base(init in) {
			switch (in) 
			{
				case init::zero:
					clear();
					break;
				case init::random:
					random();
					break;
				case init::none:
					break;
			};
		}

		table_t table_;

		static std::random_device rd_;
	};

	template<uint64_t N>
	std::random_device bf_base<N>::rd_;

	template<uint64_t N>
	class bf_tt;

	template<uint64_t N>
	class bf_anf;

	template<uint64_t N, uint64_t D>
	using bf_array = std::array<bf_tt<N>, D>;

	template<uint64_t N>
	using bf_vector = std::vector<bf_tt<N>>;

	//template<uint64_t N>
	//class bf_vector : public std::vector<bf_tt<N>> {
	//public:
	//	bf_vector() = default;
	//	bf_vector(size_t initial) {
	//		reserve(initial);
	//	}
	//};

	template<uint64_t N, uint64_t K>
	using bf_tuple = std::array<bf_tt<N>, K>;

	template<uint64_t N, uint64_t K>
	using bf_tuple_vector = std::vector<bf_tuple<N, K>>;

	template<uint64_t N>
	using bf_set = std::set<bf_tt<N>>;

	template<uint64_t N>
	using bf_pair = std::pair<bf_tt<N>, bf_tt<N>>;
}
