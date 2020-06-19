#pragma once

#include <cstdint>
#include <iterator>
#include <vector>
#include <fstream>
#include <algorithm>		// std::reverse
#include "combinatorics.h"
#include "platform.h"

#define UNUSED(v) (void)(v)

/// Lookup table of Hamming weights for 8-bit values
extern uint64_t _ByteWeight[];

template<typename T>
uint64_t popcnt(T x)
{
	uint64_t weight{ 0 };

	for (int i = 0; i < sizeof(T); i++)
	{
		weight += _ByteWeight[x & 0xff];
		x >>= 8;
	}

	return weight;
}

#if defined(ENABLE_HARDWARE_POPCNT)
template<>
uint64_t popcnt<uint32_t>(uint32_t x);

template<>
uint64_t popcnt<uint64_t>(uint64_t x);
#endif // ENABLE_HARDWARE_POPCNT


/// Calculates the least significant bit index (lsb) of the supplied value
uint64_t lsb_index(uint64_t a);

/// Calculates the most significant bit index (msb) of the supplied value
uint64_t msb_index(uint64_t a);

/// returns the least significant bit of the input
template<typename T>
T lsb(const T t)  {
	return t & (static_cast<T>(1));
}

/// checks if the input is an odd number
template<typename T>
bool is_odd(const T t) {
	return lsb(t) == 1;
}

/// checks if the input is an even number
template<typename T>
bool is_even(const T t) {
	return lsb(t) == 0;
}

/// calculates the parity of the input
template<typename T>
uint64_t parity(const T x) {
	return (popcnt(x) & 1);
}

/// checks if the input has odd parity
template<typename T>
bool is_odd_parity(const T t) {
	return is_odd(parity(t));
}

/// checks if the input has even parity
template<typename T>
bool is_even_parity(const T t) {
	return is_even(parity(t));
}

/// computes the integer 2^x
constexpr uint64_t power_of_2(uint64_t x) {
	return (1ULL << x);
}

constexpr bool is_power_of_2(uint64_t t) {
	return (t & (t - 1)) == 0;
}

template<typename T>
T log2_ceil(T x) {

	return msb_index(x) + (popcnt(x & (x - 1)) > 0);
}

template<typename T>
constexpr bool select(uint64_t s, const T t1, const T t2) {
	return (s == 0 ? t2 : t1);
}

template<typename T>
uint8_t word_byte(T t, uint64_t byteIndex) {
	return (t >> (8 * byteIndex)) & 0xff;
}

template<typename T>
constexpr uint64_t word_bytes() {
	return sizeof(T);
}

template<typename T>
constexpr uint64_t word_bits() {
	return sizeof(T) * 8;
}

template<typename T, size_t N>
constexpr size_t block_count() {
	return (N + word_bits<T>() - 1) / word_bits<T>();
}

///
template<typename T>
std::string format_integer(T t, char sep = ',') {

	std::string digits;

	// The maximum size of the formatted numbers for each type size
	//  8-bits ->  3 chars
	// 16-bits ->  6 chars
	// 32-bits -> 13 chars
	// 64-bits -> 26 chars
	// (#bits / 3) is a good approximation for the upper bound on size
	digits.reserve(word_bits<T>() / 3);

	if (std::is_signed<T>::value && (t < 0)) {
		t *= -1;
		digits.push_back('-');
	}

	do
	{
		digits.push_back('0' + (t % 10));

		if ((digits.size() + 1) % 4 == 0)
			digits.push_back(sep);

		t /= 10;

	} while (t != 0);

	if (digits.back() == sep)
		digits.pop_back();

	std::reverse(digits.begin(), digits.end());

	return digits;
}


/// create a mask of all 1s for type T
template<typename T>
constexpr T full_mask() {
	return static_cast<T>(-1);
}

// create a mask of L 1s in least significant bits for type T
template<typename T, uint64_t L>
constexpr T partial_mask() {
	static_assert(L < (8 * sizeof(T)), "Partial mask invalid length.");
	return ((L == 0) ? 0 : ((static_cast<T>(1) << (L)) - 1));
	//return bf::select(L, (static_cast<T>(1) << L) - 1, T{});
	//return (full_mask<T>() >> ((8 * sizeof(T) - L)) * (L != 0));
	//return full_mask<T>() >> kronecker_delta(8 * sizeof(T) - L);
}

template<typename T, uint64_t L>
struct lsb_mask;

template<typename T>
struct lsb_mask<T, 0> {
	static const T value = T{};
};

template<typename T, uint64_t L>
struct lsb_mask {
	static const T value = (lsb_mask<T, L - 1>::value << 1) ^ T { 1 };
};

template<uint64_t>
struct log_of_2;

template<>
struct log_of_2<1> {
	static const uint64_t value = 0;
};

template<uint64_t N>
struct  log_of_2 {
	static_assert(N % 2 == 0, "N must be a power of 2.");
	static const uint64_t value = 1 + log_of_2<N / 2>::value;
};

template<uint64_t kb>
struct KB {
	constexpr static uint64_t bits	= power_of_2(13) * kb;
	constexpr static uint64_t bytes = power_of_2(10) * kb;
};

template<uint64_t mb>
struct MB {
	constexpr static uint64_t bits   = power_of_2(23) * mb;
	constexpr static uint64_t bytes  = power_of_2(20) * mb;
	constexpr static uint64_t kbytes = power_of_2(10) * mb;
};

//! checks whether the inputs have value 1 in the same bit positions
template<typename T>
bool is_disjoint(const T t1, const T t2)
{
	return (t1 & t2) == 0;
}

//! returns the i'th bit of the variable t
template<typename T>
T get_bit(const T &t, uint64_t i)
{
	return (t >> i) & static_cast<T>(1);
}

template<typename T>
T unit_vector(uint64_t i) 
{
	return static_cast<T>(1) << i;
}

template<typename T>
void set_bit(T &t, uint64_t pos)
{
	t |= unit_vector<T>(pos);
}

template<typename T>
void clear_bit(T &t, uint64_t pos)
{
	t &= (full_mask<T>() ^ unit_vector<T>(pos));
}

template<typename T>
void flip_bit(T &t, uint64_t pos)
{
	t ^= unit_vector<T>(pos);
}

//
// STL Helpers
//
template<typename FwIter, typename Func>
void for_all_pairs(FwIter first, FwIter last, Func f)
{
	if (first != last)
	{
		FwIter trailer = first;
		++first;

		for (; first != last; ++first, ++trailer)
			for (FwIter it = first; it != last; ++it)
				f(*trailer, *it);
	}
}

template<typename Cont>
uint64_t distinct_pairs(const Cont& cont)
{
	uint64_t result{ 0 };

	for_all_pairs(cont.begin(), cont.end(), [&](const typename Cont::const_reference e1,
		const typename Cont::const_reference e2) {
		result += (e1 != e2);
	});

	return result;
}

template<typename Cont>
bool all_distinct(const Cont& cont) {

	return distinct_pairs(cont) == combinatorics::choose(cont.size(), 2);
}

template<typename Cont>
void print_container(const Cont& cont, char sep, bool newline = true) {

	for (const auto &x : cont)
		std::cout << x << sep;

	if (newline)
		std::cout << std::endl;
}


std::vector<uint8_t> hex_to_bytes(const std::string &hexstr);

uint64_t file_size(std::ifstream &file);

class mask
{
public:

	enum class type { all, nonzero };

	mask(uint64_t size, type t = type::all) : size_(size), t_(t) {}

	class Iterator : public std::iterator<std::forward_iterator_tag, std::vector<uint64_t>>
	{
	public:

		Iterator(uint64_t size, uint64_t value) : size_(size), value_(value) {

			vec_.reserve(size);
		}

		bool operator!=(const Iterator &rhs) const {

			return value_ != rhs.value_;
		}

		const std::vector<uint64_t>& operator*() const {

			vec_.clear();

			for (uint64_t i = 0; i < 8 * sizeof(uint64_t); i++)
				if (get_bit(value_, i))
					vec_.push_back(i);

			return vec_;
		}

		Iterator operator++() {
			value_++;
			return *this;
		}

	private:
		uint64_t size_;
		uint64_t value_;
		mutable std::vector<uint64_t> vec_;
	};

	Iterator begin() {
		return Iterator(size_, (t_ == type::nonzero) ? 1 : 0);
	}

	Iterator end() {
		return Iterator(size_, 1ULL << size_);
	}

private:

	uint64_t size_;
	type t_;
};
