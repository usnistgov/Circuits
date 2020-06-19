#pragma once

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <type_traits>
#include "utils.h"


namespace bfl {

	namespace  {

		template<typename T>
		struct word_dimension {
			static const uint64_t value = log_of_2<8 * sizeof(T)>::value;
		};

		template<typename T, uint64_t N, uint64_t BlockCount>
		struct mobius_block {

			static const uint64_t BlockSize = (N / BlockCount);

			template<typename Arr>
			static void transform(Arr &arr) {

				for (uint64_t i = 0; i < BlockCount; i++)
				{
					for (uint64_t j = 0; j < BlockSize / 2; j++)
						arr[i * BlockSize + BlockSize / 2 + j] ^= arr[i * BlockSize + j];
				}

				mobius_block<T, N, BlockCount / 2>::transform(arr);
			}
		};

		template<typename T, uint64_t N>
		struct mobius_block<T, N, 0> {
			template<typename Arr>
			static void transform(Arr&) {
			}
		};
	}

	template<uint64_t N, typename T>
	struct table_block {

		static_assert((N >= 1) && (N <= 31), "Unsupported dimension");

		static_assert(std::is_integral<T>::value && std::is_unsigned<T>::value, "Unsigned integral type required");

		using word_t = T;

		static const uint64_t	Inputs		=	N;

		static const uint64_t	TableSize	=	1ULL << Inputs;

		static const uint64_t	TableBytes	=	(TableSize + 7) / 8;

		static const uint64_t	WordBytes	=	word_bytes<word_t>();

		static const uint64_t	WordBits	=	word_bits<word_t>();

		static const uint64_t	WordCount	=	(TableSize + WordBits - 1) / WordBits;

		static const word_t		WordMask	=	(TableSize >= WordBits) ?	full_mask<word_t>() :
																			partial_mask<word_t, TableSize % WordBits>();

		static const uint64_t	WordDimension = word_dimension<word_t>::value;

		static const uint64_t	MobiusMask[6];

		//using data_t = std::array<word_t, WordCount>;
		using data_t = std::vector<word_t>;

		table_block()
			 : data_(WordCount) 
		{
		}

		table_block(const table_block<N, T> &rhs) {
			data_ = rhs.data_;
		}

		table_block<N, T>& operator=(const table_block<N, T> &rhs) {
			data_ = rhs.data_;
			return *this;
		}

		void clear() {
			for (uint64_t i = 0; i < WordCount; i++)
				data_[i] = 0;
		}

		uint64_t get(uint64_t index) const {
			return (data_[index / WordBits] >> (index % WordBits) & 1);
		}

		void set(uint64_t index, uint64_t val = 1) {
			data_[index / WordBits] &= ~(static_cast<word_t>(1) << (index % WordBits));
			data_[index / WordBits] |= static_cast<word_t>(val) << (index % WordBits);
		}

		void clear(uint64_t index) {
			data_[index / WordBits] &= ~(static_cast<word_t>(1) << (index % WordBits));
		}

		// TODO : unify access functions
		word_t word(uint64_t index) const {
			return data_[index];
		}

		word_t operator[](uint64_t index) const {
			return data_[index];
		}

		word_t& operator[](uint64_t index) {
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
				data_[i] = ~data_[i];
			data_[WordCount - 1] &= WordMask;
		}

		void flip(uint64_t index) {
			data_[index / WordBits] ^= static_cast<word_t>(1) << (index % WordBits);
		}

		void operator^=(const table_block<N, T> &rhs) {
			for (uint64_t i = 0; i < WordCount; i++)
				data_[i] ^= rhs.data_[i];
		}

		void operator&=(const table_block<N, T> &rhs) {
			for (uint64_t i = 0; i < WordCount; i++)
				data_[i] &= rhs.data_[i];
		}

		void operator|=(const table_block<N, T> &rhs) {
			for (uint64_t i = 0; i < WordCount; i++)
				data_[i] |= rhs.data_[i];
		}

		bool operator==(const table_block<N, T> &rhs) const {
			for (uint64_t i = 0; i < WordCount; i++)
				if (data_[i] != rhs.data_[i])
					return false;
			return true;
		}

		bool operator!=(const table_block<N, T> &rhs) const {
			for (uint64_t i = 0; i < WordCount; i++)
				if (data_[i] != rhs.data_[i])
					return true;
			return false;
		}

		bool operator<(const table_block<N, T> &rhs) const {
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
				w += popcnt(data_[i]);
			return w;
		}

		void mobius_transform() {

			static_assert(WordDimension <= 6, "Dimensions greater than 6 are not supported");

			for (uint64_t w = 0; w < WordCount; w++)
				for (uint64_t i = 0; i < WordDimension; i++)
				{
					// TODO: erroneus compilation in clang release mode
					//data_[w] ^= (data_[w] & static_cast<word_t>(MobiusMask[i])) << (static_cast<word_t>(1) << i);

					word_t tmp = (data_[w] & static_cast<word_t>(MobiusMask[i]));
					tmp <<= (static_cast<word_t>(1) << i);
					data_[w] ^= tmp;					
				}

			data_[WordCount - 1] &= WordMask;

			mobius_block<T, WordCount, WordCount / 2>::transform(data_);
		}

		std::string bin_str() const {

			std::string s(TableSize, '0');

			for (uint64_t i = 0; i < TableSize; i++)
				s[i] = (static_cast<char>('0' + get(TableSize - 1 - i)));

			return s;
		}

		void set_hex(const std::string &hexstr) {
			
			clear();

			auto b = hex_to_bytes(hexstr);

			for (int64_t i = b.size() - 1; i >= 0; i--)
			{
				auto index = b.size() - 1 - i;
				data_[index / WordBytes] ^= static_cast<word_t>(b[i]) << (8 * (index % WordBytes));
			}
		}

		std::string hex_str() const {

			std::ostringstream os;

			auto b = bytes();

			for(auto w : b) {
				os << std::hex << std::setw(2) << std::setfill('0') << static_cast<uint64_t>(w);
			}

			return os.str();
		}

		std::array<uint8_t, TableBytes> bytes() const {

			std::array<uint8_t, TableBytes> r;

			// TODO : change the ordering of the bytes in the result
			for (uint64_t i = 0; i < TableBytes; i++)
				r[TableBytes - 1 - i] = word_byte(data_[(i / WordBytes)], i % WordBytes);
			
			return r;
		}

	private:
		data_t data_;
	};

	template<uint64_t N, typename T>
	const uint64_t table_block<N, T>::MobiusMask[6] = {
		0x5555555555555555,
		0x3333333333333333,
		0x0f0f0f0f0f0f0f0f,
		0x00ff00ff00ff00ff,
		0x0000ffff0000ffff,
		0x00000000ffffffff
	};

	template<uint64_t N>
	using table_u8 = table_block<N, uint8_t>;

	template<uint64_t N>
	using table_u16 = table_block<N, uint16_t>;

	template<uint64_t N>
	using table_u32 = table_block<N, uint32_t>;

	template<uint64_t N>
	using table_u64 = table_block<N, uint64_t>;
}

