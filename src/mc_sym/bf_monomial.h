#pragma once

#include <cstdint>
#include <string>

namespace bfl {


	// Variable naming
	enum class ordering { rtl1, ltr1, rtl0, ltr0 };

	template<uint64_t N>
	struct monomial {

		static uint64_t index(const std::string &monstr, ordering ord = ordering::rtl1) {

			uint64_t mon{ 0 };

			if (monstr != "0" && monstr != "1")
			{
				int64_t index{ -1 };

				for (auto c : monstr)
				{
					//if (c == ' ' || c == '\t')
					//{
					//	// ignore whitespace
					//}
					//else 
					if (isdigit(c))
					{
						if (index == -1)
							index = (c - '0');
						else
							index = 10 * index + (c - '0');
					}
					else if(index != -1)
					{
						mon ^= (static_cast<uint64_t>(1) << char_to_index(index, ord));
						index = -1;
					}
				}

				if(index != -1)
					mon ^= (static_cast<uint64_t>(1) << char_to_index(index, ord));
			}

			return mon;
		}

		static std::string str(uint64_t mon, ordering ord = ordering::rtl1)
		{
			if (mon == 0)
				return { "1" };

			std::string s;

			for (uint64_t i = 0; i < N; i++)
			{
				if ((mon >> i) & 1)
				{
					s.append("x");
					s.append(std::to_string(index_to_char(i, ord)));
				}
			}

			return s;
		}

	private:

		static uint64_t char_to_index(uint64_t index, ordering ord) {

			if (ord == ordering::rtl1)
				return index - 1;
			else if (ord == ordering::rtl0)
				return index;
			else if (ord == ordering::ltr1)
				return N - 1 - (index - 1);
			else
				return N - 1 - (index);
		}

		static uint64_t index_to_char(uint64_t index, ordering ord) {

			if (ord == ordering::rtl1)
				return index + 1;
			else if (ord == ordering::rtl0)
				return index;
			else if (ord == ordering::ltr1)
				return N - index;
			else
				return N - (index + 1);
		}

	};

	template<uint64_t N>
	std::vector<uint64_t> get_monomials(uint64_t degree) {

		std::vector<uint64_t> ret;

		for (uint64_t x = 0; x < power_of_2(N); x++) {
			
			if(popcnt(x) == degree)
				ret.push_back(x);
		}

		return ret;
	}


}
