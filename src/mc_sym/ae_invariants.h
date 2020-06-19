#pragma once

#include "bf.h"
#include "spectrum.h"

namespace ae {

	using namespace bfl;

	struct affine_invariants
	{
		uint64_t degree;
		spectrum_distribution walsh_dist;
		spectrum_distribution ac_dist;

		affine_invariants() = default;

		~affine_invariants() = default;

		template<uint64_t N>
		affine_invariants(const bf_tt<N> &f) {
			*this = calculate_affine_invariants(f);
		}

		affine_invariants(const affine_invariants &rhs) 
			:	degree(rhs.degree), 
				walsh_dist(rhs.walsh_dist), 
				ac_dist(rhs.ac_dist) {
		}

		affine_invariants& operator=(const affine_invariants &rhs) {		
			degree		= rhs.degree;
			walsh_dist	= rhs.walsh_dist;
			ac_dist		= rhs.ac_dist;
			return *this;
		}

		affine_invariants(affine_invariants &&rhs)
			:	degree(rhs.degree), 
				walsh_dist(std::move(rhs.walsh_dist)), 
				ac_dist(std::move(rhs.ac_dist)) {
		}

		affine_invariants& operator=(affine_invariants &&rhs) {
			degree		= rhs.degree;
			walsh_dist	= std::move(rhs.walsh_dist);
			ac_dist		= std::move(rhs.ac_dist);
			return *this;
		}

		void print() const {
			std::cout << degree << " : ";
			for (auto m : ac_dist)
				std::cout << m.first << '-' << m.second << ' ';
			std::cout << " : ";
			for (auto m : walsh_dist)
				std::cout << m.first << '-' << m.second << ' ';
			std::cout << std::endl;
		}

		void print_compact() const {

			std::cout << degree << std::endl;

			std::cout << ac_dist.size() << ' ';
			for (const auto &m : ac_dist)
				std::cout << static_cast<int64_t>(m.first) << ' ' << static_cast<int64_t>(m.second) << ' ';
			std::cout << std::endl;

			std::cout << walsh_dist.size() << ' ';
			for (const auto &m : walsh_dist)
				std::cout << static_cast<int64_t>(m.first) << ' ' << static_cast<int64_t>(m.second) << ' ';
			std::cout << std::endl;
		}
	};

	bool operator==(const affine_invariants &ai1, const affine_invariants &ai2);

	bool operator!=(const affine_invariants &ai1, const affine_invariants &ai2);

	bool operator<(const affine_invariants &ai1, const affine_invariants &ai2);

	template<uint64_t N>
	affine_invariants calculate_affine_invariants(const bf_tt<N> &f)
	{
		affine_invariants ai;

		ai.degree = f.anf().degree();

		auto walsh = f.walsh_spectrum();

		ai.walsh_dist = walsh.absolute_distribution(spectrum<N>::zeros::exclude);

		ai.ac_dist = f.autocorrelation_spectrum(walsh).absolute_distribution(spectrum<N>::zeros::exclude);

		return ai;
	}

	template<uint64_t N>
	bool compare_affine_invariants(const bf_tt<N> &f, const affine_invariants &ai)
	{
		if (f.anf().degree() != ai.degree)
			return false;

		auto walshF = f.walsh_spectrum();

		if (walshF.absolute_distribution(spectrum<N>::zeros::exclude) != ai.walsh_dist)
			return false;

		if (f.autocorrelation_spectrum(walshF).absolute_distribution(spectrum<N>::zeros::exclude) != ai.ac_dist)
			return false;

		return true;
	}

	template<uint64_t N>
	bool compare_affine_invariants(const bf_tt<N> &f, const bf_tt<N> &g)
	{
		if (f.anf().degree() != g.anf().degree())
			return false;

		auto walshF = f.walsh_spectrum();
		auto walshG = g.walsh_spectrum();

		if (walshF.absolute_distribution(spectrum<N>::zeros::exclude) !=
			walshG.absolute_distribution(spectrum<N>::zeros::exclude))
			return false;

		if (f.autocorrelation_spectrum(walshF).absolute_distribution(spectrum<N>::zeros::exclude) !=
			g.autocorrelation_spectrum(walshG).absolute_distribution(spectrum<N>::zeros::exclude))
			return false;

		return true;
	}
}

namespace std {

	template<>
	struct hash<ae::affine_invariants> {

		size_t operator()(const ae::affine_invariants &ai) const {
			
			return	hash<uint64_t>()(ai.degree) ^
					hash<bfl::spectrum_distribution>()(ai.walsh_dist) ^
					hash<bfl::spectrum_distribution>()(ai.ac_dist);
		}
	};
}

