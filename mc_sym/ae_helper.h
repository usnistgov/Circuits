#pragma once

#include <cstdint>
#include <vector>
#include <numeric>

#include "bf.h"

namespace ae {

	using namespace bfl;

	template<uint64_t N>
	bf_vector<N> load_representatives(bfl::sort srt = bfl::sort::no) {

		static_assert(N >= 2 && N <= 6, "");

		return load_functions<N>("../data/classes/n" + std::to_string(N) + "_classes.txt", mode::anf, srt);
	}

	class elimination_list
	{
	public:
		elimination_list(uint64_t itemCount) : _vec(itemCount) {
			std::iota(_vec.begin(), _vec.end(), 0);
		}

		uint64_t count() const {
			return _vec.size();
		}

		bool remove(uint64_t index) {

			auto backIndex = count() - 1;

			if (index != backIndex)
				//std::swap(_vec[index], _vec[backIndex]);
				_vec[index] = _vec[backIndex];

			_vec.pop_back();

			return (count() == 1);
		}

		uint64_t get(uint64_t index) const {
			return _vec[index];
		}

		uint64_t head() const {
			return _vec[0];
		}

	private:
		std::vector<uint64_t> _vec;
	};
}


