
#include "bf.h"
#include "bf_sym.h"
#include "mc.h"
#include "slp.h"
#include "utils.h"
#include "benchmark.h"

using namespace bfl;
using namespace bfl::sym;
using namespace ae;


uint64_t weight_bits(uint64_t n)
{
	for (uint64_t s = 1; s < 64; s++)
	{
		if ((n + 1) <= power_of_2(s))
			return s;
	}

	return 0;
}

struct bit_configuration
{
	bit_configuration(uint64_t N) : _N(N), digits(weight_bits(N)) {
		digits[0] = N;
	}

	uint64_t reduce_all() {

		uint64_t reductions{ 0 };

		for(uint64_t i = 0; i < digits.size() - 1; i++)
			reductions += reduce_column(i);

		return reductions;
	}

	uint64_t reduce_upto_degree(uint64_t degree) {

		uint64_t reductions{ 0 };

		auto pos = degree_max_col(degree);

		for (uint64_t i = 0; i < pos; i++)
			reductions += reduce_column(i);

		return reductions;
	}

	uint64_t reduce_until_varcount(uint64_t n, uint64_t degree = 0) {

		uint64_t reductions{ 0 };

		auto vars = bits();

		uint64_t maxColIndex = (degree == 0) ? digits.size() : degree_max_col(degree);

		for (uint64_t col = 0; maxColIndex; col++)
		{
			while ((vars > n) && reduce_column_by_one(col))
			{
				reductions++;
				vars--;
			}
		}

		return reductions;
	}

	uint64_t reduce_column(uint64_t col) {

		uint64_t triplets{ 0 };

		while (reduce_column_by_one(col))
			triplets++;

		return triplets;
	}

	bool reduce_lossy(uint64_t col) {

		if (digits[col] == 0 || (col == digits.size() - 1))
			return false;

		--digits[col];
		++digits[col + 1];

		return true;
	}

	bool reduce_column_by_one(uint64_t col) {

		if (digits[col] < 3)
			return false;

		digits[col] -= 2;
		digits[col + 1] += 1;

		return true;
	}

	uint64_t bits() const {
		return std::accumulate(digits.begin(), digits.end(), 0ULL);
	}

	uint64_t width() const {
		return std::count_if(digits.begin(), digits.end(), [](uint64_t c) {
			return c != 0;
		});
	}

	void print() const {

		uint64_t maxval = *(std::max_element(digits.begin(), digits.end()));

		auto hwbits = digits.size();

		for (uint64_t y = maxval; y > 0; y--)
		{
			for (uint64_t x = 0; x < hwbits; x++)
			{
				if (digits[hwbits - 1 - x] >= y)
					std::cout << '*';
				else if (y == 1)
					std::cout << '.';
				else
					std::cout << ' ';
			}

			std::cout << std::endl;
		}
	}

	uint64_t degree_max_col(uint64_t degree) {

		uint64_t digit_pos{ 0 };

		while ((digit_pos < digits.size()) && (power_of_2(digit_pos + 1) <= degree))
			digit_pos++;

		return digit_pos;
	}

	std::vector<std::vector<uint64_t>> weight_indices() const {

		std::vector<uint64_t> coef;

		uint64_t currentCoef{ 1 };
		for (auto x : digits)
		{
			for (uint64_t i = 0; i < x; i++)
				coef.push_back(currentCoef);

			currentCoef *= 2;
		}

		std::vector<std::vector<uint64_t>> wi(_N + 1);

		for (uint64_t i = 0; i < power_of_2(coef.size()); i++)
		{
			uint64_t weight{ 0 };

			for (uint64_t j = 0; j < coef.size(); j++)
				if (get_bit(i, j))
					weight += coef[j];
	
			wi[weight].push_back(i);
		}

		return wi;
	}

	std::vector<uint64_t> weight_coefficients() const {

		std::vector<uint64_t> wc;

		uint64_t current_coefficient{ 1 };
		for (auto x : digits)
		{
			for (uint64_t i = 0; i < x; i++)
				wc.push_back(current_coefficient);

			current_coefficient *= 2;
		}

		return wc;
	}

	uint64_t max_weight() const {

		auto wc = weight_coefficients();

		return std::accumulate(std::begin(wc), std::end(wc), 0ULL);
	}

	std::vector<uint64_t> block2_indices() const {

		std::vector<uint64_t> indices;

		bool inblock{ false };
		for (uint64_t i = 0; i < digits.size(); i++)
		{
			if (digits[i] == 2)
			{
				if (!inblock)
				{
					indices.push_back(i);
					inblock = true;
				}
			}
			else
				inblock = false;


		}

		return indices;
	}

	std::vector<uint64_t> digits;
	uint64_t _N;
};


void print_bit_configurations(uint64_t N, uint64_t cols)
{
	std::vector<bit_configuration> bcv;

	bit_configuration bc(N);

	bcv.push_back(bc);

	for (uint64_t i = 0; i < cols; i++)
	{
		bc.reduce_column(i);
		bcv.push_back(bc);
	}

	auto hwbits = weight_bits(N);

	std::cout << "N = " << N << std::endl;

	for (uint64_t j = N; j > 0; j--)
	{
		for (uint64_t i = 0; i < bcv.size(); i++)
		{
			for (uint64_t k = 0; k < cols + 1; k++)
				if (bcv[bcv.size() - 1 - i].digits[cols - k] >= j)
					std::cout << '*';
				else
					std::cout << ' ';

			std::cout << "   ";
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;
}


template<uint64_t N>
circuit_data<N> generate_slp(const bf_tt<N> &f)
{
	static auto cs = class_distinguisher<N>();
	static auto sd = load_circuit_data<N>();

	auto slp = build_slp(f, *cs, sd);

	return slp;
}


template<uint64_t N, uint64_t M>
void produce_slp(const bf_vector<N> &hwbits_serialized, 
				 const std::vector<std::vector<uint64_t>>& vecmap, 
				 const uint64_t reductions,
				 const uint64_t minDegree)
{
	// initialize class distinguisher and circuit data
	generate_slp(bf_tt<M>{});

	std::map<uint64_t, uint64_t> mcmap;

	uint64_t mc_tight{ 0 };

	progress_bar pb("producing slp's", power_of_2(N + 1));

#pragma omp parallel for //schedule(dynamic)
	for (int wmask = 1; wmask < power_of_2(N + 1); wmask += 1)
	{
#pragma omp critical
		pb.update();

		//bf_tt<N> fsym(init::zero);
		bf_tt<M> g(init::zero);


		bf_sym<N> fs(wmask);

		//fs.transform();

		auto Sigmadecomp = fs.get_elementary_decomposition();

		if (Sigmadecomp.empty() || Sigmadecomp[0] == 0 || Sigmadecomp[0] == 1)
			continue;

		auto degree = Sigmadecomp.back();

		if (degree < minDegree)
			continue;

		auto Edecomp = fs.get_counting_decomposition();

		for (auto e : Edecomp)
			for (auto x : vecmap[e])
				g.set(x);

		auto slp = generate_slp(g);

		if (vecmap.size() > (N + 1) && is_odd_parity(g.weight()))
		{
			g.set(power_of_2(M) - 1);

			auto slp2 = generate_slp(g);

			if (slp2.topology.gateCount() < slp.topology.gateCount())
			{
//#pragma omp critical
//				{
//					std::cout << "MC decreased from " << slp.topology.gateCount() << " to " << slp2.topology.gateCount() << std::endl;
//				}

				slp = slp2;
			}
		}

//		if (slp.topology.gateCount() == M)
//		{
//#pragma omp critical
//			std::cout << "MC 6 observed" << std::endl;
//		}

		auto mc = reductions + slp.topology.gateCount();

		TopologySLP<N> slpex(slp.topology);

		bf_vector<N> inputs;
		
		inputs.reserve(slp.inputs.size());

		bf_tt<N> ffinal(init::zero);

		for (auto f : slp.inputs)
		{
			bf_tt<N> fi(init::zero);

			auto fa = f.anf();

			for (uint64_t i = 0; i < M; i++)
				if (fa.get(power_of_2(i)))
					fi ^= hwbits_serialized[i];

			if (fa.get(0))
				fi.complement();

			inputs.push_back(fi);
		}

		{
			auto fa = slp.ffinal.anf();

			for (uint64_t i = 0; i < M; i++)
				if (fa.get(power_of_2(i)))
					ffinal ^= hwbits_serialized[i];

			if (fa.get(0))
				ffinal.complement();
		}

		auto fgen = slpex.evaluateMask(inputs, slp.outputMask);

		fgen ^= ffinal;

		auto fsym = fs.to_bf();

		if (fgen == fsym)
		{
#pragma omp critical
			{
				//std::cout << "#fsym  : " << fsym << std::endl;

				//std::cout << "#E     : ";
				//for (auto x : Edecomp)
				//	std::cout << x << ' ';
				//std::cout << '\n';

				//std::cout << "#Sigma : ";
				//for (auto x : Sigmadecomp)
				//	std::cout << x << ' ';
				//std::cout << '\n';

				//std::cout << "mc : " << mc << std::endl;

				//slp.print();

				//std::cout << slp.topology.str() << std::endl;
				//std::cout << "inputs:" << '\n';
				//for (auto &l : slp.inputs)
				//	std::cout << l << '\n';
				//std::cout << "adjusted inputs:" << '\n';
				//for (auto &l : inputs)
				//	std::cout << l << '\n';
				//std::cout << ffinal << '\n';

				//slpex.str(inputs, ffinal, slp.outputMask);

				if (mc == (degree - 1))
					mc_tight++;

				++mcmap[mc];
			}
		}
		else
		{
			LOG_ERROR << "***SLP does not match the function" << LOG_ENDL;

#pragma omp critical
			{
				std::cout << "#E     : ";
				for (auto x : Edecomp)
					std::cout << x << ' ';
				std::cout << '\n';

				std::cout << "#Sigma : ";
				for (auto x : Sigmadecomp)
					std::cout << x << ' ';
				std::cout << '\n';

				slp.print();
			}

			LOG_DEBUG << "adjusted inputs:" << LOG_ENDL;
			for (auto &l : inputs)
				LOG_DEBUG << l << LOG_ENDL;
			LOG_DEBUG << ffinal << LOG_ENDL;

			LOG_DEBUG << "fsym : " << fsym << LOG_ENDL;
			LOG_DEBUG << "fgen : " << fgen << LOG_ENDL;
		}
	}

	std::cout << std::endl;
	std::cout << "mc distribution:" << std::endl;
	for (auto m : mcmap)
		std::cout << m.first << " : " << m.second << std::endl;

	std::cout << "tight : " << mc_tight << std::endl;
}



template<uint64_t N>
void symmetric_generate_slp()
{
	std::cout << "Generating SLPs for N = " << N << std::endl;

	std::vector<bf_vector<N>> hwbits(weight_bits(N));

	for (uint64_t i = 1; i <= N; i++)
		hwbits[0].push_back(bf_tt<N>("x" + std::to_string(i)));

	uint64_t current_bits{ N };
	uint64_t reductions{ 0 };


	// do an initial half-adder if N=22
	if (N == 22)
	{
		auto f1 = hwbits[0].back(); hwbits[0].pop_back();
		auto f2 = hwbits[0].back(); hwbits[0].pop_back();

		bf_tt<N> sum(f1);
		sum ^= f2;

		bf_tt<N> prod(f1);
		prod &= f2;

		hwbits[0].push_back(sum);
		hwbits[0 + 1].push_back(prod);

		reductions++;
	}

	for (uint64_t hw_coordinate = 0; hw_coordinate < hwbits.size() - 1; hw_coordinate++)
	{
		LOG_DEBUG << "#bits in partial hamming weight computation : " << current_bits << LOG_ENDL;

		if (current_bits <= 6)
		{
			LOG_DEBUG << "proceeding to the second phase" << LOG_ENDL;

			bf_vector<N> hwbits_serialized;

			for (auto &v : hwbits)
			{
				if (v.empty())
				{
					LOG_DEBUG << "<empty>" << LOG_ENDL;
				}
				else
					for (auto &f : v)
					{
						hwbits_serialized.push_back(f);

						if (N <= 11)
							LOG_DEBUG << f << LOG_ENDL;
					}

				LOG_DEBUG << LOG_ENDL;
			}

			std::vector<uint64_t> bitweights;
			uint64_t bitcount{ 0 };

			uint64_t current_weight{ 1 };
			uint64_t total_weight{ 0 };
			for (auto &v : hwbits)
			{
				bitcount += v.size();

				for (uint64_t i = 0; i < v.size(); i++)
					bitweights.push_back(current_weight);

				total_weight += current_weight * v.size();

				current_weight *= 2;
			}

			std::vector<std::vector<uint64_t>> vecmap(total_weight + 1);

			for (uint64_t mask = 0; mask < power_of_2(bitcount); mask++)
			{
				uint64_t sum{ 0 };
				for (uint64_t m = 0; m < bitcount; m++)
					if (get_bit(mask, m))
						sum += bitweights[m];

				//if(sum <= N)
				vecmap[sum].push_back(mask);
			}

			/*
			std::cout << "bitcount : " << bitcount << std::endl;
			std::cout << "weights  : ";
			for (auto x : bitweights)
				std::cout << x << ' ';
			std::cout << std::endl;

			std::cout << "weight classes:" << std::endl;
			for (uint64_t w = 0; w < vecmap.size(); w++)
			{
				std::cout << w << " : ";
				for (auto x : vecmap[w])
					std::cout << x << ' ';
				std::cout << std::endl;
			}
			*/

			auto minDegree = power_of_2(hw_coordinate);

			LOG_DEBUG << "min degree : " << minDegree << LOG_ENDL;
			LOG_DEBUG << "reductions : " << reductions << LOG_ENDL;

			if (bitcount == 2)
				produce_slp<N, 2>(hwbits_serialized, vecmap, reductions, minDegree);
			if (bitcount == 3)
				produce_slp<N, 3>(hwbits_serialized, vecmap, reductions, minDegree);
			else if (bitcount == 4)
				produce_slp<N, 4>(hwbits_serialized, vecmap, reductions, minDegree);
			else if (bitcount == 5)
				produce_slp<N, 5>(hwbits_serialized, vecmap, reductions, minDegree);
			else if (bitcount == 6)
				produce_slp<N, 6>(hwbits_serialized, vecmap, reductions, minDegree);
			else
			{
				std::cout << "skipping the case bitcount = " << bitcount << std::endl;
			}

			break;
		}
		else
		{
			LOG_DEBUG << "processing coordinate #" << hw_coordinate << " of the hamming weight function" << LOG_ENDL;

			while (hwbits[hw_coordinate].size() >= 3)
			{
				if ((hw_coordinate == 2) && (hwbits[hw_coordinate].size() == 3))
					break;

				auto f1 = hwbits[hw_coordinate].back(); hwbits[hw_coordinate].pop_back();
				auto f2 = hwbits[hw_coordinate].back(); hwbits[hw_coordinate].pop_back();
				auto f3 = hwbits[hw_coordinate].back(); hwbits[hw_coordinate].pop_back();

				bf_tt<N> sum(f1);
				sum ^= f2;
				sum ^= f3;
				hwbits[hw_coordinate].push_back(sum);

				bf_tt<N> p1(f1); p1 ^= f2;
				bf_tt<N> p2(f1); p2 ^= f3;
				p1 &= p2;
				p1 ^= f1;
				hwbits[hw_coordinate + 1].push_back(p1);

				reductions++;
				current_bits--;
			}
		}
	}

	std::cout << std::endl;
}


template<uint64_t N>
void symmetric_generate_slp_rec();

// termination condition
template<>
void symmetric_generate_slp_rec<26>() {
}

template<uint64_t N>
void symmetric_generate_slp_rec() {

	symmetric_generate_slp<N>();

	symmetric_generate_slp_rec<N + 1>();
}

int main(int argc, char *argv[])
{
	benchmark bm("main");

	// compute mc upper bounds for a particular N
	//symmetric_generate_slp<10>();

	// compute mc upper bounds from N=2 up to N=25
	symmetric_generate_slp_rec<2>();

	return 0;
}
