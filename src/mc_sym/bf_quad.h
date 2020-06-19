#pragma once

#include "bf.h"
#include "combinatorics.h"

namespace bfl {

	using namespace combinatorics;

	template<uint64_t N, typename T = uint64_t>
	class bf_quad {

	public:

		explicit bf_quad(T terms = 0) : _terms(terms) {
			assert(monomial_count() <= word_bits<T>());
		}

		bf_quad(const bf_tt<N> &f) : bf_quad() {

			auto  &ind = monomial_indices();

			auto g = f.anf();

			for (uint64_t i = 0; i < bf_tt<N>::TableSize; i++)
			{
				if (g.get(i))
				{
					auto p = std::find(ind.begin(), ind.end(), i);

					if (p == std::end(ind))
						LOG_WARN << "bfquad : discarding term " << monomial<N>::str(i) << LOG_ENDL;
					else
					{
						uint64_t index = std::distance(ind.begin(), p);

						if (index >= word_bits<T>())
						{
							LOG_WARN << "bfquad : index out of range" << index << LOG_ENDL;
							_terms = 0;
							break;
						}

						_terms |= unit_vector<T>(index);
					}
				}
			}					
		}

		bf_tt<N> to_bf() const {

			bf_anf<N> f(init::zero);

			auto &ind = monomial_indices();
			for (uint64_t i = 0; i < monomial_count(); i++)
				if (get_bit(_terms, i))
					f.set(ind[i]);

			return f.tt();
		}

		T as_integer() const {
			return _terms;
		}

		bf_quad<N> operator^=(bf_quad<N> rhs) {
			return bf_quad<N>(_terms ^ rhs._terms);
		}

		static constexpr uint64_t monomial_count() {
			return N * (N - 1) / 2;
		}

		static const std::vector<uint64_t>& monomial_indices() {

			static std::vector<uint64_t> indices;

			if (indices.empty())
			{
				indices.reserve(choose(N, 2));

				for (auto sub : subset(N, 2)) {
					T term{};
					for (auto x : sub)
						term |= unit_vector<T>(x);
					indices.push_back(term);
				}
			}

			return indices;
		}

	private:
		T _terms;
	};

	template<uint64_t N, typename T>
	bool operator==(bf_quad<N, T> a, bf_quad<N, T> b) {
		return a.as_integer() == b.as_integer();
	}

	template<uint64_t N, typename T>
	bool operator<(bf_quad<N, T> a, bf_quad<N, T> b) {
		return a.as_integer() < b.as_integer();
	}
}
