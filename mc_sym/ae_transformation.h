#pragma once

#include <array>
#include <numeric>

#include "bf.h"
#include "bf_utils.h"
#include "matrix.h"
#include "combinatorics.h"
#include "ae_invariants.h"
#include "logger.h"
#include "utils.h"
#include "constant_time.h"

namespace ae {

	using namespace bfl;
	using namespace combinatorics;

	template<uint64_t N>
	struct affine_mapping
	{
		enum class init{uninitialized, zero, identity};

		affine_mapping(init in = init::uninitialized){

			if (in == init::identity) {

				A = matrix<N>::ae_identity();
				a = 0;
				b.clear();
			}
			else if (in == init::zero) {

				A.clear();
				a = 0;
				b.clear();
			}
		}

		matrix<N>	A;
		uint64_t	a;
		bf_tt<N>	b;
	};

	template<uint64_t N>
	using perm = std::array<uint64_t, power_of_2(N)>;

	template<uint64_t N>
	using perm_vector = std::vector<perm<N>>;

	template<uint64_t N>
	perm<N> construct_permutation(const matrix<N> &A, uint64_t a = 0)
	{
		const uint64_t dimension = N;
		const uint64_t tableSize = power_of_2(dimension);

		perm<N> permutation;

		for (uint64_t i = 0; i < tableSize; i++)
		{
			uint64_t index = 0;

			for (int j = 0; j < dimension; j++)
			{
				//index <<= 1;
				index ^= (popcnt(A[j] & i) & 1) << j;
			}

			index ^= a;

			permutation[i] = index;
		}

		return permutation;
	}

	template<uint64_t N>
	perm<N> construct_permutation(const affine_mapping<N> &am)
	{
		return construct_permutation(am.A, am.a);
	}

	template<uint64_t N>
	bf_tt<N> transform(const bf_tt<N> &f, const perm<N> &perm)
	{
		bf_tt<N> g;

		for (uint64_t i = 0; i < bf_tt<N>::TableSize; i++)
			g.set(i, f.get(perm[i]));

		return g;
	}

	template<uint64_t N>
	bf_tt<N> transform(const bf_tt<N> &f, const matrix<N> &A, uint64_t a = 0)
	{
		auto perm = construct_permutation(A, a);

		return transform(f, perm);
	}

	template<uint64_t N>
	bf_tt<N> transform_inner(const bf_tt<N> &f, const affine_mapping<N> &am)
	{
		return transform(f, am.A, am.a);
	}

	template<uint64_t N>
	bf_tt<N> transform(const bf_tt<N> &f, const affine_mapping<N> &am)
	{
		auto g = transform(f, am.A, am.a);

		g ^= am.b;

		return g;
	}

	template<uint64_t N, class FwIter, class OutIter>
	void transform(FwIter first, FwIter last, OutIter dest, const perm<N> &perm)
	{
		while (first != last)
			*dest++ = transform(*first++, perm);
	}

	template<uint64_t N, class FwIter, class OutIter>
	void transform(FwIter first, FwIter last, OutIter dest, const affine_mapping<N> &am)
	{
		auto perm = construct_permutation(am);

		while (first != last)
		{
			*dest = transform(*first, perm);
			*dest ^= am.b;

			++dest;
			++first;
		}
	}

	template<uint64_t N>
	bool verify_affine_transformation(const bf_tt<N> &f, const bf_tt<N> &g, const affine_mapping<N> &am)
	{
		return (transform(f, am) == g);
	}

	template<uint64_t N>
	perm_vector<N> get_all_transformations()
	{
		perm_vector<N> ret;

		invertible_matrix_generator<N> img;

		do
		{
			for (uint64_t a = 0; a < power_of_2(N); a++)
				ret.push_back(construct_permutation(img.get_matrix(), a));

		} while (img.next());

		return ret;
	}

	template<uint64_t N, uint64_t K>
	perm_vector<N> get_self_mappings(const bf_tuple<N, K> &tuple, const perm_vector<N> &transformations)
	{
		perm_vector<N> ret;

		bf_tuple<N, K> original(tuple);

		bfl::remove_affine_terms<N>(std::begin(original), std::end(original));

		for (auto const &t : transformations)
		{
			bf_tuple<N, K> transformed;

			transform<N>(std::begin(original), std::end(original), std::begin(transformed), t);

			bfl::remove_affine_terms<N>(std::begin(transformed), std::end(transformed));

			if (transformed == original)
				ret.push_back(t);
		}

		return ret;
	}

	template<uint64_t N>
	class dimension_reducer
	{
	public:
		dimension_reducer() {
			init();
		}

		std::pair<bf_tt<N>, bool> reduce(const bf_tt<N> &f) const {

			auto dimOriginal = f.anf().dependent_variable_count();

			auto g = f;

			for (uint64_t i = 0; i < N; i++)
			{
				for (auto &perm : permvec_)
				{
					auto h = transform(g, perm).anf();

					h.remove_affine_terms();

					auto dim = h.dependent_variable_count();

					if (dim < dimOriginal)
						return { h.tt(), true };
				}

				g = transform(g, permrotate_);
			}

			return { bf_tt<N>{}, false };
		}

		bool is_reducible(const bf_tt<N> &f) const {

			auto p = reduce(f);

			return p.second;
		}

	private:

		void init() {

			permvec_.reserve(power_of_2(N - 1) - 1);

			matrix<N> m = matrix<N>::ae_identity();

			for (uint64_t mask = 1; mask < power_of_2(N - 1); mask++)
			{
				m.setColumn(N - 1, mask | power_of_2(N - 1));

				permvec_.push_back(construct_permutation(m));
			}

			matrix<N> mrotate;
			for (uint64_t i = 0; i < N; i++)
				mrotate.setRow(i, (1ULL << ((i + 1) % N)));

			permrotate_ = construct_permutation(mrotate, 0);

		}

		std::array<uint64_t, power_of_2(N)> permrotate_;
		std::vector<std::array<uint64_t, power_of_2(N)>> permvec_;
	};

	template<uint64_t N>
	bf_tt<N> swap_variables(const bf_tt<N> &f, uint64_t i, uint64_t j)
	{
		LOG_VERIFY(i >= 0 && i <= (N - 1));
		LOG_VERIFY(j >= 0 && j <= (N - 1));

		matrix<N> A = matrix<N>::ae_identity();

		A[i] = unit_vector<uint64_t>(j);
		A[j] = unit_vector<uint64_t>(i);

		return transform(f, A);
	}

	template<uint64_t N>
	bf_tt<N> rotate_variables_left(const bf_tt<N> &f, uint64_t r = 1)
	{
		matrix<N> mrotate;

		for (uint64_t i = 0; i < N; i++)
			mrotate.setRow(i, unit_vector<uint64_t>(((i + r) % N)));

		return transform(f, mrotate);
	}

	template<uint64_t N>
	bf_tt<N> rotate_variables_right(const bf_tt<N> &f, uint64_t r = 1)
	{
		return rotate_variables_left(f, N - (r % N));
	}

	namespace
	{
		template<uint64_t N>
		class weight_dist
		{
		public:
			weight_dist(const bf_tt<N> &f)
			{
				index_base_[0] = 0;
				for (uint64_t i = 1; i <= N; i++)
					index_base_[i] = index_base_[i - 1] + power_of_2(i);

				std::fill(std::begin(dist_), std::end(dist_), 0);

				for (uint64_t i = 0; i < power_of_2(N); i++)
				{
					for (uint64_t l = 0; l < N; l++)
						dist_[index(l, i >> (N - 1 - l))] += f.get(i);
				}			
			}

			uint64_t get(uint64_t level, uint64_t part) const {

				return dist_[index(level, part)];
			}

			bool is_trivial_a0(uint64_t level) const {

				const uint64_t Blocks = power_of_2(level);

				const uint64_t BlockSize = power_of_2(N) / Blocks;

				for (uint64_t b = 0; b < 2 * Blocks; b++)
				{
					auto value1 = get(level, b);

					if (!(value1 == 0 || value1 == BlockSize / 2))
						return false;
				}

				return true;
			}

			bool is_trivial_a1(uint64_t level) const {

				const uint64_t Blocks = power_of_2(level);

				const uint64_t BlockSize = power_of_2(N) / Blocks;

				for (uint64_t b = 0; b < Blocks; b++)
				{
					auto value1 = get(level, 2 * b + 0);
					auto value2 = get(level, 2 * b + 1);

					if (!(value1 == 0 || value1 == BlockSize / 2) ||
						!(value2 == 0 || value2 == BlockSize / 2)
						)
						return false;
				}

				return true;
			}

		private:

			uint64_t index(uint64_t level, uint64_t part) const {

				return index_base_[level] + part;
			}

			std::array<uint64_t, N + 1> index_base_;
			std::array<uint64_t, power_of_2(N + 1) - 2> dist_;
		};

		template<uint64_t N>
		class independent_subspace
		{
		public:

			static const uint64_t Elements = power_of_2(N);

			independent_subspace() {

				span_.reserve(Elements);
				span_.push_back(0);

				std::fill(std::begin(used_), std::end(used_), false);
			}

			void add_row(uint64_t row) {

				uint64_t prevSize = span_.size();

				for (uint64_t i = 0; i < prevSize; i++)
				{
					uint64_t newElement = span_[i] ^ row;

					span_.push_back(newElement);

					used_[newElement] = true;
				}
			}

			void remove_last_row() {

				for (uint64_t i = span_.size() / 2; i < span_.size(); i++)
					used_[span_[i]] = false;

				span_.resize(span_.size() / 2);
			}

			std::vector<uint64_t> get_vectors() {

				std::vector<uint64_t> vecs;

				vecs.reserve(Elements - span_.size());

				for (uint64_t i = 1; i < Elements; i++)
					if (!used_[i])
						vecs.push_back(i);

				return vecs;
			}

			uint64_t get_next_vector(uint64_t &prev) {

				while (prev < Elements)
				{
					if (!used_[prev])
						return prev++;
					else
						prev++;
				}

				return 0;
			}

		private:

			std::vector<uint64_t> span_;
			std::array<bool, Elements> used_;
		};

		template<typename Iter>
		bool split_indices(Iter first, Iter last, uint64_t value) {

			auto mid = first + (last + 1 - first) / 2;

			while (first < mid && last >= mid)
			{
				auto cond = parity(*first & value);

				conditional_swap(cond, *first, *last);

				first += (cond ^ 1);
				last -= cond;
			}

			return true;
		}
	}


	template<uint64_t N, uint64_t M>
	class find_transformation;

	template<uint64_t N>
	class find_transformation<N, 0>
	{
	public:
		static bool find(const std::array<uint64_t, power_of_2(N)> &transformation,
			independent_subspace<N> &sub,
			const weight_dist<N> &fwdist,
			const std::array<uint8_t, power_of_2(N)> &g,
			matrix<N> &A, uint64_t &a)
		{
			UNUSED(transformation);
			UNUSED(sub);
			UNUSED(fwdist);
			UNUSED(g);
			UNUSED(A);
			UNUSED(a);

			//for (uint64_t i = 0; i < power_of_2(N); i++)
			//{
			//	if (g[transformation[i]] != f[i])
			//		return false;
			//}

			return true;
		}
	};

	template<uint64_t N, uint64_t M>
	class find_transformation
	{
	public:

		static const uint64_t Level = (N - M);

		static const uint64_t Blocks = power_of_2(Level);

		static const uint64_t BlockSize = power_of_2(N) / Blocks;

		static bool find(
			const std::array<uint64_t, power_of_2(N)> &transformation,
			independent_subspace<N> &sub,
			const weight_dist<N> &fwdist,
			const std::array<uint8_t, power_of_2(N)> &g,
			matrix<N> &A, uint64_t &a)
		{
			uint64_t value;

			uint64_t vecIndex = 1;

			while((value = sub.get_next_vector(vecIndex)) != 0)
			{
				LOG_TRACE_IF((N - M) < 2) << "setting row " << (M - 1) << " to " << value << LOG_ENDL;

				A.setRow(M - 1, value);

				auto tcopy = transformation;

				split_indices(tcopy, value);

				// check for a_i = 0
				bool valid{ true };

				for (uint64_t b = 0; b < 2 * Blocks; b++)
				{
					uint64_t wg{ 0 };

					for (uint64_t j = 0; j < BlockSize / 2; j++)
						wg += g[tcopy[b * (BlockSize / 2) + j]];

					if (wg != fwdist.get(Level, b))
					{
						valid = false;
						break;
					}
				}

				if (valid)
				{
					a &= ~(1ULL << (M - 1));

					sub.add_row(value);

					if (find_transformation<N, M - 1>::find(tcopy, sub, fwdist, g, A, a))
						return true;
					
					sub.remove_last_row();
				}

				// check for a_i = 1
				valid = true;

				for (uint64_t b = 0; b < Blocks; b++)
				{
					uint64_t wg0{ 0 }, wg1{ 0 };

					for (uint64_t j = 0; j < BlockSize / 2; j++)
					{
						wg0 += g[tcopy[(2 * b + 0) * (BlockSize / 2) + j]];
						wg1 += g[tcopy[(2 * b + 1) * (BlockSize / 2) + j]];

						std::swap(tcopy[(2 * b + 0) * (BlockSize / 2) + j], tcopy[(2 * b + 1) * (BlockSize / 2) + j]);
					}

					if (wg0 != fwdist.get(Level, 2 * b + 1) || wg1 != fwdist.get(Level, 2 * b + 0))
					{
						valid = false;
						break;
					}
				}

				if (valid)
				{
					a |= (1ULL << (M - 1));

					sub.add_row(value);

					if (find_transformation<N, M - 1>::find(tcopy, sub, fwdist, g, A, a))
						return true;
					
					sub.remove_last_row();
				}
			}

			return false;
		}

	private:

		static void split_indices(std::array<uint64_t, power_of_2(N)> &dest, uint64_t value) {

			for (uint64_t b = 0; b < Blocks; b++)
			{
				auto first = b * BlockSize;
				auto last = (b + 1) * BlockSize - 1;

				auto mid = first + (last + 1 - first) / 2;

				while (first < mid && last >= mid)
				{
					//auto cond = parity(dest[first] & value);

					//ct::swap_if(dest[first], dest[last], cond);

					//first += (cond ^ 1);
					//last -= cond;

					if (parity(dest[first] & value))
					{
						std::swap(dest[first], dest[last]);
						last--;
					}
					else
						first++;
				}
			}
		}
	};


	//! Search for a vector 'a' that transforms f to g
	template<uint64_t N>
	bool is_affine_equivalent_a(const bf_tt<N> &f, const bf_tt<N> &g, affine_mapping<N> &am) {

		auto A = matrix<N>::ae_identity();

		auto gn(g);
		gn.remove_affine_terms();

		for (uint64_t a = 0; a < power_of_2(N); a++) {

			auto h = transform(f, A, a);

			auto h2(h);
			h2.remove_affine_terms();

			if (h2 == g)
			{
				am.A = A;
				am.a = a;
				am.b = h2 ^ h;
				return true;
			}
		}

		return false;
	}

	template<uint64_t N>
	bool is_affine_equivalent(const bf_tt<N> &f, const bf_tt<N> &g, affine_mapping<N> &am)
	{
		bool retVal{ false };

		//if (is_affine_equivalent_a(f, g, am))
		//{
		//	retVal = true;
		//}
		//else
		{
			auto fw = f.weight();

			weight_dist<N> fwdist(f);

			bf_vector<N> lcandidate;

			for (auto &l : bf_tt<N>::get_affine_functions())
			{
				if (bf_tt<N>(g ^ l).weight() == fw)
					lcandidate.push_back(l);
			}

			LOG_TRACE << "f : " << f.anf().str() << LOG_ENDL;
			LOG_TRACE << "g : " << g.anf().str() << LOG_ENDL;
			LOG_TRACE << "affine candidates : " << lcandidate.size() << LOG_ENDL;
			for (auto l : lcandidate)
			{
				LOG_TRACE << l << LOG_ENDL;
			}

			for (auto &l : lcandidate)
			{
				LOG_TRACE << "processing affine function : " << l << LOG_ENDL;

				bf_tt<N> g1(g ^ l);

				std::array<uint64_t, power_of_2(N)> transformation;

				std::iota(std::begin(transformation), std::end(transformation), 0);

				auto garr = to_bytes(g1);

				affine_mapping<N> amlocal(affine_mapping<N>::init::zero);

				independent_subspace<N> sub;

				auto eq = find_transformation<N, N>::find(transformation, sub, fwdist, garr, amlocal.A, amlocal.a);

				if (eq)
				{
					am = amlocal;
					am.b = l;
					retVal = true;
					break;
				}

				LOG_TRACE << eq << " : " << l << LOG_ENDL;
			}
		}

		LOG_ASSERT_IF(retVal, verify_affine_transformation(f, g, am));

		return retVal;
	}

	template<uint64_t N>
	bool is_affine_equivalent_parallel(const bf_tt<N> &f, const bf_tt<N> &g, affine_mapping<N> &am)
	{
		auto fw = f.weight();

		bool retVal{ false };

		weight_dist<N> fwdist(f);

		bf_vector<N> lcandidate;

		for (auto &l : bf_tt<N>::get_affine_functions())
		{
			bf_tt<N> g1(g);

			g1 ^= l;

			if (g1.weight() == fw)
				lcandidate.push_back(l);
		}

		LOG_TRACE << "affine candidates : " << lcandidate.size() << LOG_ENDL;

#pragma omp parallel for schedule(dynamic)
		for (int lIndex = 0; lIndex < lcandidate.size(); lIndex++)
		{
			if (retVal)
				continue;

			auto l = lcandidate[lIndex];

#pragma omp critical
			LOG_TRACE << "processing affine function : " << l << LOG_ENDL;

			bf_tt<N> g1(g);

			g1 ^= l;

			std::array<uint64_t, power_of_2(N)> transformation;

			std::iota(std::begin(transformation), std::end(transformation), 0);

			auto garr = to_bytes(g1);

			affine_mapping<N> amlocal(affine_mapping<N>::init::zero);

			independent_subspace<N> sub;

			auto eq = find_transformation<N, N>::find(transformation, sub, fwdist, garr, amlocal.A, amlocal.a);

#pragma omp critical
			LOG_TRACE << eq << " : " << l << LOG_ENDL;

			if (eq)
			{
#pragma omp critical
				{
					if (!retVal)
					{
						am = amlocal;
						am.b = l;
						retVal = true;
						//break;
					}

				}
			}
		}

		return retVal;
	}

	template<uint64_t N>
	bool is_affine_equivalent(const bf_tt<N> &f, const bf_tt<N> &g)
	{
		affine_mapping<N> am;

		return is_affine_equivalent(f, g, am);
	}

	template<uint64_t N>
	bool is_affine_equivalent_wsig(const bf_tt<N> &f, const bf_tt<N> &g)
	{
		if (!compare_signatures(f, g))
			return false;

		return is_affine_equivalent(f, g);
	}

	template<uint64_t N>
	[[deprecated]] bool find_affine_mapping(const bf_tt<N> &f, const bf_tt<N> &g, affine_mapping<N> &am)
	{
		bf_tt<N> gr(g);
		gr.remove_affine_terms();

		bf_tt<N> g0(g);
		g0.flip(0);

		auto aig0 = calculate_affine_invariants(g0);

		std::vector<ae::affine_invariants> aivec;

		aivec.reserve(power_of_2(N));

		for (uint64_t i = 0; i < power_of_2(N); i++)
		{
			bf_tt<N> fi(f);
			fi.flip(i);
			aivec.push_back(calculate_affine_invariants(fi));
		}

		std::vector<uint64_t> a_candidates;

		for (uint64_t i = 0; i < power_of_2(N); i++)
		{
			if (aivec[i] == aig0)
				a_candidates.push_back(i);
		}

		std::vector<std::vector<uint64_t>> d_candidates(N);

		for (uint64_t j = 0; j < N; j++)
		{
			bf_tt<N> gj(g);

			gj.flip(1ULL << j);

			auto ai = calculate_affine_invariants(gj);

			for (uint64_t i = 0; i < power_of_2(N); i++)
				if (ai == aivec[i])
					d_candidates[j].push_back(i);
		}

		if (a_candidates.empty())
			return false;

		for (auto &v : d_candidates)
			if (v.empty())
				return false;

		enumeration e_dk(N);
		for (uint64_t i = 0; i < N; i++)
			e_dk.set_limit(i, d_candidates[i].size());

		uint64_t search_space = a_candidates.size();
		for (uint64_t i = 0; i < N; i++)
			search_space *= d_candidates[i].size();

		auto sslog = log2(search_space);

		LOG_DEBUG << "find_mapping search space: " << sslog << LOG_ENDL;

		matrix<N> A;

		for (auto a : a_candidates)
		{
			for (auto &dk : e_dk)
			{
				bool column_used{ false };

				for (uint64_t i = 0; i < N && !column_used; i++)
				{
					for (uint64_t j = 0; j < i; j++)
					{
						if (d_candidates[i][dk[i]] == d_candidates[j][dk[j]])
						{
							column_used = true;
							break;
						}
					}

					A[i] = a ^ d_candidates[i][dk[i]];
				}

				if (column_used)
					continue;

				auto t = transform(f, A, a);

				auto tr = t;

				tr.remove_affine_terms();

				if (tr == gr)
				{
					am.A = A;
					am.a = a;
					am.b = t;
					am.b ^= g;
					return true;
				}
			}
		}

		return false;
	}
}