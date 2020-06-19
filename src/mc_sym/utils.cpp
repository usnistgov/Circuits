#include "utils.h"
#include <cassert>
#include <immintrin.h>

uint64_t _ByteWeight[] = {
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

#if defined(ENABLE_HARDWARE_POPCNT)
template<>
uint64_t popcnt<uint32_t>(uint32_t x) {
	return _mm_popcnt_u32(x);
}

template<>
uint64_t popcnt<uint64_t>(uint64_t x) {
	return _mm_popcnt_u64(x);
}
#endif // ENABLE_HARDWARE_POPCNT

uint64_t lsb_index(uint64_t a)
{
//	if (a == 0)
//		return static_cast<uint64_t>(-1);
//
//#ifdef _MSC_VER
//	unsigned long res;
//#else
//	uint32_t res;
//#endif
//
//	_BitScanForward64(&res, a);
//	return res;

	for (uint64_t i = 0; i < word_bits<uint64_t>(); i++)
		if (get_bit(a, i))
			return i;

	return static_cast<uint64_t>(-1);
}


uint64_t msb_index(uint64_t a)
{
//	if(a == 0)
//		return static_cast<uint64_t>(-1);
//
//#ifdef _MSC_VER
//	unsigned long res;
//#else
//	uint32_t res;
//#endif
//
//	_BitScanReverse64(&res, a);
//	return res;

	for (int64_t i = word_bits<uint64_t>() - 1; i >= 0; i--)
		if (get_bit(a, i))
			return i;

	return static_cast<uint64_t>(-1);
}


struct hex_digit
{
	static bool to_uint8(const char ch, uint8_t &val) {

		if (is_decimal(ch))
			val = ch - '0';
		else if (is_lowercase_hexadecimal(ch))
			val = ch - 'a' + 10;
		else if (is_uppercase_hexadecimal(ch))
			val = ch - 'A' + 10;
		else
			return false;

		return true;
	}

	static bool is_decimal(const char ch) {
		return (ch >= '0' && ch <= '9');
	}

	static bool is_lowercase_hexadecimal(const char ch) {
		return	(ch >= 'a' && ch <= 'f');
	}

	static bool is_uppercase_hexadecimal(const char ch) {
		return	(ch >= 'A' && ch <= 'F');
	}
};


std::vector<uint8_t> hex_to_bytes(const std::string &hexstr) {

	std::vector<uint8_t> result;

	result.reserve((hexstr.size() + 1) / 2);

	uint64_t	pos{ 0 };
	uint8_t		byte{ 0 };

	// if there are odd number of hex digits, pretend that there is a leading zero digit
	if (is_odd(hexstr.size()))
		pos++;

	for (auto ch : hexstr)
	{
		byte <<= 4;

		uint8_t value;

		auto valid = hex_digit::to_uint8(ch, value);

		if (valid)
			byte ^= value;
		else 
			throw std::invalid_argument(std::string("invalid hexadecimal character : ") + ch);

		pos++;

		if (is_even(pos))
		{
			result.push_back(byte);
			byte = 0;
		}
	}

	assert(is_even(pos));

	return result;
}

uint64_t file_size(std::ifstream &file)
{
	const auto curpos = file.tellg();

	file.seekg(0, std::ios_base::end);

	const auto size = file.tellg();

	file.seekg(curpos);

	return size;
}
