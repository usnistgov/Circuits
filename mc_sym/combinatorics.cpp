#include "combinatorics.h"
#include <algorithm>

namespace combinatorics {


	uint64_t choose(uint64_t n, uint64_t k) {

		assert(k <= n);

		if (k == 0 || n <= 1)
			return 1;

		k = std::min(k, n - k);

		//if ((n > 50) && std::abs((static_cast<int64_t>(n / 2) - static_cast<int64_t>(k)) < 10))
		//	return combin(n - 1, k - 1) + combin(n - 1, k);

		uint64_t result{ 1 };

		for (uint64_t i = 1; i <= k; i++) {
			result *= n - (i - 1);
			result /= i;
		}

		return result;
	}

	uint64_t factorial(uint64_t n) {

		uint64_t res{ 1 };

		for (uint64_t i = 2; i <= n; i++)
			res *= i;

		return res;
	}


	subset::subset(uint64_t objectCount, uint64_t chooseCount)
		: _objectCount(objectCount),
		_chooseCount(chooseCount),
		_items(chooseCount)
	{
		assert(chooseCount <= objectCount);

		init();
	
		_terminateIndex = choose(objectCount, chooseCount);
	}

	subset::~subset()
	{
	}

	void subset::init()
	{
		for (uint64_t i = 0; i < _chooseCount; i++)
			_items[i] = i;

		_stateIndex = 0;
	}

	uint64_t subset::count() const {

		return (_terminateIndex - _stateIndex);

		//if (_chooseCount == 0)
		//	return 1;

		//uint64_t value = _objectCount;

		//for (uint64_t i = 2; i <= _chooseCount; i++) {
		//	value *= _objectCount - (i - 1);
		//	value /= i;
		//}

		//return value;
	}

	bool subset::next()
	{
		if (_chooseCount == 0)
			return false;

		_items[_chooseCount - 1]++;

		for (int64_t a = _chooseCount - 1; a >= 0; a--)
		{
			if (_items[a] >= (_objectCount - ((_chooseCount - 1) - a)))
			{
				if (a == 0)
				{
					// End of search states
					init();
					return false;
				}

				_items[a - 1]++;

				for (uint64_t b = a; b < _chooseCount; b++)
					_items[b] = _items[b - 1] + 1;
			}
			else
				break;
		}

		_stateIndex++;

		return true;
	}

	//////
	permutation::permutation(uint64_t uItems) : _state(uItems) {
		init();
	}

	void permutation::init()
	{
		std::iota(_state.begin(), _state.end(), 0);
	}

	bool permutation::next()
	{
		uint64_t a, t;

		const auto items = _state.size();

		// From "Algorithms and Programming", Alexander Shen, page 34.
		for (int64_t k = static_cast<int64_t>(items - 2); k >= 0; k--)
		{
			// find an element that is smaller than the next
			if (_state[k] < _state[k + 1])
			{
				for (t = k + 1; t < items; t++)
				{
					if (_state[k] > _state[t])
						break;
				}

				// exchange
				std::swap(_state[k], _state[t - 1]);

				// reverse
				for (a = 0; a < (items - (k + 1)) / 2; a++)
				{
					std::swap(_state[k + a + 1], _state[items - 1 - a]);
				}

				return true;
			}
		}

		return false;
	}

}