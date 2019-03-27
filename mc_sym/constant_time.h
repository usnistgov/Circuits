#pragma once

#include <cassert>
#include "utils.h"

namespace ct
{
	template<typename T>
	T extend_lsb(T t) {

		return ~((t & 1) - 1);
	}

	template<typename T>
	T msb(T t) {

		return (t >> (word_bits<T>() - 1));
	}

	template<typename T>
	T neg(T t) {

		return (~t) + static_cast<T>(1);
	}

	template<typename T>
	auto abs(T t) -> typename std::make_unsigned<T>::type {

		static_assert(std::is_signed<T>::value, "");

		//
		// http://graphics.stanford.edu/~seander/bithacks.html#IntegerAbs
		//
		T mask{ t >> (word_bits<T>() - 1) };

		return  (t + mask) ^ mask;
	}

	template<typename T>
	void swap_if(T &a, T &b, T cond) {

		T sum	{ a ^ b };
		T mask	{ extend_lsb<T>(cond) };

		a ^= sum & mask;
		b ^= sum & mask;
	}

	//! xor b to a if cond is true
	template<typename T>
	void xor_if(T &a, T b, T cond) {

		T mask = extend_lsb<T>(cond);

		a ^= (b & mask);
	}

	template<typename T>
	void set_if(T &a, T b, T cond) {

		a ^= select<T>(a ^ b, static_cast<T>(0), cond);
	}

	/// if cond is true return a, otherwise return b
	/// res = c.(a+b)+b
	template<typename T>
	T select(T a, T b, T cond) {

		return (extend_lsb<T>(cond) & (a ^ b)) ^ b;
	}

	template<typename T>
	T is_nonzero(T t) {

		auto n{ neg<T>(t) };
		auto p{ msb<T>(t & n) };

		return msb<T>(t ^ n) ^ p;
	}

	template<typename T>
	T is_zero(T t) {

		return is_nonzero<T>(t) ^ static_cast<T>(1);
	}

	template<typename T>
	T is_all_ones(T t) {

		return is_zero(t + static_cast<T>(1));
	}

	template<typename T>
	T is_greater(T a, T b) {

		return msb(b - a);
	}

	template<typename T>
	T is_smaller(T a, T b) {

		return msb(a - b);
	}

	template<typename T>
	T is_not_equal(T a, T b) {

		return (is_smaller<T>(a, b) | is_greater<T>(a, b));
	}

	template<typename T>
	T is_equal(T a, T b) {

		return is_not_equal<T>(a, b) ^ static_cast<T>(1);
	}

	template<typename T>
	T is_power_of_two(T a) {

		return is_zero<T>(a & (a - 1));
	}
}

