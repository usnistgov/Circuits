#pragma once

#include "bf.h"
#include "matrix.h"
#include "ae.h"
#include "ae_transformation.h"
#include <array>
#include <vector>
#include <cassert>


namespace bfl {

	namespace sym {

		using namespace ae;

		template<uint64_t N>
		bf_tt<N> elementary(uint64_t k) {

			assert(k <= N);

			bf_anf<N> e;

			for (uint64_t x = 0; x < e.table_size(); x++)
			{
				if (popcnt(x) == k)
					e.set(x);
			}

			return e.tt();
		}

		template<uint64_t N, typename... Args>
		bf_tt<N> elementary(uint64_t k, Args... args)
		{
			assert(k <= N);

			auto f = elementary<N>(k);
			f ^= elementary<N>(args...);
			return f;
		}

		template<uint64_t N>
		bf_vector<N> get_elementary_symmetric_functions() {

			bf_vector<N> ret(N + 1);

			for (uint64_t k = 0; k <= N; k++)
				ret[k] = elementary<N>(k);

			return ret;
		}

		template<uint64_t N>
		bf_vector<N> get_all_symmetric_functions() {

			bf_vector<N> ret;
			ret.reserve(power_of_2(N + 1));

			auto sigma = get_elementary_symmetric_functions<N>();

			for (uint64_t mask = 0; mask < power_of_2(N + 1); mask++)
			{
				bf_tt<N> f(init::zero);

				for (uint64_t i = 0; i < (N + 1); i++)
					if (get_bit(mask, i))
						f ^= sigma[i];

				ret.push_back(f);
			}

			return ret;
		}

		template<uint64_t N>
		bf_vector<N> get_nonlinear_symmetric_functions() {

			auto elementaries = sym::get_elementary_symmetric_functions<N>();

			bf_vector<N> allsym;

			for (uint64_t mask = 1; mask < power_of_2(N - 1); mask++)
			{
				bf_tt<N> f(init::zero);

				for (uint64_t m = 0; m < (N - 1); m++)
					if ((mask >> m) & 1)
						f ^= elementaries[2 + m];

				allsym.push_back(f);
			}

			return allsym;
		}

		template<uint64_t N>
		bf_tt<N> get_symmetric_mask_degree(uint64_t mask)
		{
			bf_tt<N> f(init::zero);

			for (uint64_t i = 0; i <= N; i++)
				if (get_bit(mask, i))
					f ^= elementary<N>(i);

			return f;
		}

		template<uint64_t N>
		bf_tt<N> get_symmetric_mask_weight(uint64_t mask)
		{
			bf_tt<N> f;

			for (uint64_t x = 0; x < power_of_2(N); x++)
				f.set(x, get_bit(mask, popcnt(x)));

			return f;
		}

		template<uint64_t N>
		bf_tt<N> counting(uint64_t k)
		{
			assert(k <= N);

			bf_tt<N> e(init::zero);

			uint64_t kfac{ 1 };
			for (uint64_t i = 2; i <= k; i++)
				kfac *= i;

			for (uint64_t i = k; i <= N; i++)
			{
				uint64_t comb{ 1 };

				for (uint64_t j = 0; j < k; j++)
					comb *= (i - j);

				comb /= kfac;

				if (comb & 1)
					e ^= elementary<N>(i);
			}

			return e;
		}

		template<uint64_t N>
		bf_tt<N> threshold(uint64_t k)
		{
			assert(k <= N);

			bf_tt<N> e(init::zero);

			uint64_t kfac{ 1 };
			for (uint64_t i = 2; i < k; i++)
				kfac *= i;

			for (uint64_t i = k; i <= N; i++)
			{
				uint64_t comb{ 1 };

				for (uint64_t j = 1; j < i; j++)
					comb *= (i - j);

				comb /= kfac;

				if (comb & 1)
					e ^= elementary<N>(i);
			}

			return e;
		}

		template<uint64_t N>
		bool is_symmetric(const bf_anf<N> &f) {

			std::array<uint64_t, N + 1> degreeDist;

			std::fill(std::begin(degreeDist), std::end(degreeDist), 0);

			for (uint64_t x = 1; x < f.table_size(); x++)
			{
				if (f.get(x))
					degreeDist[popcnt(x)]++;
			}

			uint64_t termCount{ N };
			for (uint64_t w = 1; w < N; w++)
			{
				if ((degreeDist[w] > 0) && (degreeDist[w] != termCount))
					return false;

				termCount *= (N - w);
				termCount /= (w + 1);
			}

			return true;
		}

		template<uint64_t N>
		bool is_symmetric(const bf_tt<N> &f) {

			std::array<uint64_t, N + 1> output;

			output[0] = f.get(0);

			for (uint64_t w = 1; w <= N; w++)
				output[w] = f.get(power_of_2(w) - 1);

			for (uint64_t x = 0; x < power_of_2(N); x++)
				if (f.get(x) != output[popcnt(x)])
					return false;

			return true;

			//return is_symmetric(f.anf());
		}

		template<uint64_t N>
		std::vector<uint64_t> elementary_decomposition(const bf_anf<N> &f) {

			assert(is_symmetric(f));

			if (!is_symmetric(f))
				return {};

			std::vector<uint64_t> decomp;

			for (uint64_t d = 0; d <= N; d++)
				if (f.has_terms_with_degree(d))
					decomp.push_back(d);

			return decomp;
		}

		template<uint64_t N>
		std::vector<uint64_t> elementary_decomposition(const bf_tt<N> &f) {

			return elementary_decomposition(f.anf());
		}

		template<uint64_t N>
		std::string elementary_decomposition_str(const bf_tt<N> &f) {

			if (!is_symmetric(f))
				return { "<not symmetric>" };
			
			if (f.is_zero())
				return { "<0>" };

			std::ostringstream os;

			auto decomp = elementary_decomposition(f);
			
			bool firstTerm{ true };

			for (auto x : decomp)
			{
				if (!firstTerm)
					os << '+';
				else
					firstTerm = false;

				os << char(228) << "(" << N << "," << x << ")";
			}

			return os.str();
		}

		template<uint64_t N>
		std::vector<uint64_t> counting_decomposition(const bf_tt<N> &f) {

			std::set<uint64_t> wset;

			// handle zero weight separately
			if (f.get(0)) 
				wset.insert(0);

			// check only one truth table entry for each weight
			for (uint64_t x = 1; x <= N; x++)
				if (f.get(power_of_2(x) - 1))
					wset.insert(x);

			return std::vector<uint64_t>(wset.begin(), wset.end());		
		}

		template<uint64_t N>
		std::string counting_decomposition_str(const bf_tt<N> &f) {

			if (!is_symmetric(f))
				return { "<not symmetric>" };

			if (f.is_zero())
				return { "<0>" };

			std::ostringstream os;

			auto decomp = counting_decomposition(f);

			bool firstTerm{ true };

			for (auto x : decomp)
			{
				if (!firstTerm)
					os << '+';
				else
					firstTerm = false;

				os << "E(" << N << "," << x << ")";
			}

			return os.str();
		}

		template<uint64_t N>
		std::pair<bool, bf_tt<N - 1>> sum3_reduce(const bf_tt<N> &f)
		{
			static_assert(N > 3, "");

			std::pair<bool, bf_tt<N - 1>> res;

			for (uint64_t x = 0; x < power_of_2(N); x++)
			{
				uint64_t ms{ 0 };
				for (uint64_t i = 0; i < 3; i++)
				{
					ms <<= 1;
					ms |= get_bit(x, N - 1 - i);
				}

				uint64_t sum = popcnt(ms);

				uint64_t val = f.get(x);

				uint64_t pos = ((sum << (N - 3)) | (x & partial_mask<uint64_t, (N - 3)>()));

				// check for a contradiction of reducability
				if (!(ms == 0 || ms == 1 || ms == 3 || ms == 7))
				{
					if (val != res.second.get(pos))
					{
						res.first = false;
						return res;
					}
				}
				else
					res.second.set(pos, val);
			}

			res.first = true;

			return res;
		}

		template<uint64_t N>
		affine_mapping<N> get_twin_transformation() {

			affine_mapping<N> am;

			uint64_t blocks = N / 3;

			am.A = matrix<N>::ae_identity();
			am.a = 0;
			am.b.clear();

			for (uint64_t i = 0; i < blocks; i++)
			{
				am.A[2 * i + 0] = (0x1ULL << (2 * i)) ^ (1ULL << (N - 1 - i));
				am.A[2 * i + 1] = (0x3ULL << (2 * i));
				am.A[N - 1 - i] = (0x3ULL << (2 * i)) ^ (1ULL << (N - 1 - i));

				am.a ^= 1ULL << (N - 1 - i);
			}

			return am;
		}

		template<uint64_t N>
		affine_mapping<N> get_twin_transformation_inverse() {

			affine_mapping<N> am;

			uint64_t blocks = N / 3;

			am.A = matrix<N>::ae_identity();
			am.a = 0;
			am.b.clear();

			for (uint64_t i = 0; i < blocks; i++)
			{
				am.A[2 * i + 0] = (0x3ULL << (2 * i)) ^ (1ULL << (N - 1 - i));
				am.A[2 * i + 1] = (0x2ULL << (2 * i)) ^ (1ULL << (N - 1 - i));
				am.A[N - 1 - i] = (0x3ULL << (2 * i));

				am.a ^= (0x3ULL << (2 * i)) ^ (1ULL << (N - 1 - i));
			}

			return am;
		}

		template<uint64_t N>
		bf_tt<N> symmetric_to_twin(const bf_tt<N> &f)
		{
			auto am = get_twin_transformation<N>();

			auto tw = transform(f, am);

			//tw.remove_affine_terms();

			return tw;
		}
	
		template<uint64_t N>
		bf_tt<N> twin_to_symmetric(const bf_tt<N> &f)
		{
			auto am = get_twin_transformation_inverse<N>();

			auto tw = transform(f, am);

			//tw.remove_affine_terms();

			return tw;
		}


		// TODO:
		enum class fn_tag {elementary, counting, threshold};

		template<uint64_t N>
		class bf_sym
		{
		public:

			static_assert(N < 64, "");

			bf_sym() : _sigma_mask(0) {}

			bf_sym(uint64_t mask) : _sigma_mask(mask) {

				LOG_ASSERT(msb_index(mask) <= N);
			}

			uint64_t get_mask() const {
				return _sigma_mask;
			}

			uint64_t degree() const {

				// constant functions are assigned degree 1 for consistency with bflib
				if (_sigma_mask == 0)
					return 1;

				return msb_index(_sigma_mask);
			}

			void print() const {
				for (uint64_t i = 0; i <= N; i++)
					if (get_bit(_sigma_mask, i))
						std::cout << i << ' ';
				std::cout << std::endl;
			}

			void transform() {

				_sigma_mask = do_transform(_sigma_mask);
				_sigma_mask &= partial_mask<uint64_t, N + 1>();
			}

			bf_tt<N> to_bf() const {

				auto e_decomp = do_transform(_sigma_mask);

				bf_tt<N> f;

				for (uint64_t x = 0; x < power_of_2(N); x++)
					f.set(x, get_bit(e_decomp, popcnt(x)));

				return f;
			}

			std::vector<uint64_t> get_elementary_decomposition() const {

				std::vector<uint64_t> v;

				for (uint64_t i = 0; i <= N; i++)
					if (get_bit(_sigma_mask, i))
						v.push_back(i);

				return v;
			}

			std::vector<uint64_t> get_counting_decomposition() const {

				std::vector<uint64_t> v;

				auto e_decomp = do_transform(_sigma_mask);

				for (uint64_t i = 0; i <= N; i++)
					if (get_bit(e_decomp, i))
						v.push_back(i);

				return v;
			}

			void operator^=(const bf_sym<N> &rhs) {
				_sigma_mask ^= rhs._sigma_mask;
			}


		private:

			static uint64_t do_transform(uint64_t m) {

				for (uint64_t i = 0; i < 6; i++)
					m = m ^ ((m & sigma_e_transform[i]) << power_of_2(i));

				return m;
			}

			uint64_t _sigma_mask;

			static const uint64_t sigma_e_transform[];
		};

		template<uint64_t N>
		const uint64_t bf_sym<N>::sigma_e_transform[] = {

			0x5555555555555555,
			0x3333333333333333,
			0x0f0f0f0f0f0f0f0f,
			0x00ff00ff00ff00ff,
			0x0000ffff0000ffff,
			0x00000000ffffffff
		};

		template<uint64_t N>
		void test_bf_sym()
		{
			benchmark bm("Testing N = " + std::to_string(N));

			bf_vector<N> elems;

			for (uint64_t k = 0; k <= N; k++)
				elems.push_back(sym::elementary<N>(k));

			//	progress_bar pb("processing", power_of_2(N + 1));

			uint64_t tested{ 0 };

#pragma omp parallel for
			for (int mask = 0; mask < power_of_2(N - 1); mask++)
			{
				bf_tt<N> fsym(init::zero);

				for (uint64_t i = 0; i < N - 1; i++)
					if (get_bit(mask, i))
						fsym ^= elems[i + 2];

				bf_sym<N> gsym(mask << 2);

				auto g = gsym.to_bf();

				auto sigma_decomp = sym::elementary_decomposition(fsym);
				auto counting_decomp = sym::counting_decomposition(fsym);

				LOG_ASSERT(fsym == g);
				LOG_ASSERT(sigma_decomp == gsym.get_elementary_decomposition());
				LOG_ASSERT(counting_decomp == gsym.get_counting_decomposition());

#pragma omp critical
				{
					std::cout << (++tested) << " : ";
					for (auto x : sigma_decomp)
						std::cout << x << ' ';
					std::cout << std::endl;
				}

				//#pragma omp critical
				//		pb.update();
			}
		}

	}
}