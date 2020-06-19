#pragma once

#include <cassert>
#include <array>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <type_traits>

namespace bfl {

	template<uint64_t N, typename T>
	struct table_bit {

		static_assert((N >= 3) && (N <= 12), "Supported range for N is: 3 <= N <= 12");

		static_assert(std::is_integral<T>::value && !std::is_signed<T>::value, "Unsigned integral type required");

		using word_t = T;

		static const uint64_t	Inputs = N;

		static const uint64_t	TableSize = 1ULL << Inputs;

		static const uint64_t	TableBytes = TableSize / 8;

		static const uint64_t	WordBytes = sizeof(word_t);

		static const uint64_t	WordBits = 1;

		static const uint64_t	WordCount = TableSize;

		table_bit() {
		}

		table_bit(const table_bit<N, T> &rhs) {
			data_ = rhs.data_;
		}

		table_bit<N, T>& operator=(const table_bit<N, T> &rhs) {
			data_ = rhs.data_;
			return *this;
		}

		void clear() {
			for (uint64_t i = 0; i < WordCount; i++)
				data_[i] = 0;
		}

		uint64_t get(uint64_t index) const {
			assert(index < TableSize);
			return data_[index];
		}

		void set(uint64_t index, uint64_t val = 1) {
			assert(index < TableSize);
			data_[index] = static_cast<word_t>(val & 1);
		}

		void clear(uint64_t index) {
			assert(index < TableSize);
			data_[index] = 0;
		}

		word_t word(uint64_t index) const {
			assert(index < WordCount);
			return data_[index];
		}

		bool is_zero() const {
			for (uint64_t i = 0; i < WordCount; i++)
				if (data_[i] != 0)
					return false;
			return true;
		}

		void complement() {
			for (uint64_t i = 0; i < WordCount; i++)
				data_[i] ^= static_cast<word_t>(1);
		}

		void flip(uint64_t index) {
			assert(index < TableSize);
			data_[index] ^= static_cast<word_t>(1);
		}

		void operator^=(const table_bit<N, T> &rhs) {
			for (uint64_t i = 0; i < WordCount; i++)
				data_[i] ^= rhs.data_[i];
		}

		void operator&=(const table_bit<N, T> &rhs) {
			for (uint64_t i = 0; i < WordCount; i++)
				data_[i] &= rhs.data_[i];
		}

		void operator|=(const table_bit<N, T> &rhs) {
			for (uint64_t i = 0; i < WordCount; i++)
				data_[i] |= rhs.data_[i];
		}

		bool operator==(const table_bit<N, T> &rhs) const {
			for (uint64_t i = 0; i < WordCount; i++)
				if (data_[i] != rhs.data_[i])
					return false;
			return true;
		}

		bool operator<(const table_bit<N, T> &rhs) const {
			for (uint64_t i = 0; i < WordCount; i++)
				if (data_[i] < rhs.data_[i])
					return true;
				else if (data_[i] > rhs.data_[i])
					return false;
			return false;
		}

		uint64_t weight() const {
			uint64_t w{ 0 };
			for (uint64_t i = 0; i < WordCount; i++)
				w += data_[i];
			return w;
		}

		void mobius_transform() {

			mobius_block<T, WordCount, WordCount / 2>::transform(data_);
		}

		std::string bin_str() const {

			std::string s;
			s.reserve(TableSize);

			for (uint64_t i = 0; i < TableSize; i++)
				s.push_back(static_cast<char>('0' + get(TableSize - 1 - i)));

			return s;
		}

		std::string hex_str() const {

			std::ostringstream os;

			auto b = bytes();

			std::for_each(std::begin(b), std::end(b), [&](const word_t &w) {
				os << std::hex << std::setw(2) << std::setfill('0') << w;
			});

			return os.str();
		}

		std::array<uint8_t, TableBytes> bytes() const {

			std::array<uint8_t, TableBytes> r;

			//for (uint64_t i = 0; i < TableBytes; i++)
			//	r[TableBytes - 1 - i] = word_byte(data_[(i / WordBytes)], i % WordBytes);

			return r;
		}

	private:
		std::array<word_t, WordCount> data_;
	};

	template<uint64_t N>
	using table_bit8 = table_bit<N, uint8_t>;

	template<uint64_t N>
	using table_bit16 = table_bit<N, uint16_t>;

	template<uint64_t N>
	using table_bit32 = table_bit<N, uint32_t>;

	template<uint64_t N>
	using table_bit64 = table_bit<N, uint64_t>;

}