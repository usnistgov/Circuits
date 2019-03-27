#pragma once

#include "topology.h"
#include "bf.h"
#include <sstream>

// TODO: remove the namespace declaration from the header file
using namespace bfl;

template<uint64_t N>
class TopologySLP
{
public:
	TopologySLP(const Topology &topology) : _topology(topology) 
	{
		_gateCount = _topology.gateCount();

		_inputCount = _gateCount * 2;

		_terminalGates = _topology.getTerminalGates();

		_terminalOutputCount = _terminalGates.size();

		_intermediateOutputCount = _gateCount - _terminalOutputCount;

		compile();
	}
	
	~TopologySLP() {}

	uint64_t inputCount() const {

		return _inputCount;
	}

	uint64_t outputCount() const {

		return _gateCount;
	}

	uint64_t terminalOutputCount() const {

		return _terminalOutputCount;
	}

	uint64_t intermediateOutputCount() const {

		return _intermediateOutputCount;
	}

	uint64_t gateCount() const {
		
		return _gateCount;
	}

	const std::vector<uint64_t>& terminalOutputs() const {
		return _terminalOutputs;
	}

	const std::vector<uint64_t>& intermediateOutputs() const {
		return _intermediateOutputs;
	}

	std::vector<bf_tt<N>> evaluateCombinations(const bf_vector<N> &inputs)
	{
		bf_vector<N> result;

		uint64_t outputCount = _terminalOutputCount + _intermediateOutputCount;

		result.reserve(1ULL << outputCount);

		auto fout = evaluate(inputs);

		for (uint64_t mask = 1; mask < 1ULL << outputCount; mask++)
		{
			bf_tt<N> f(bfl::init::zero);

			for (uint64_t m = 0; m < outputCount; m++)
			{
				if ((mask >> m) & 1)
					f ^= fout[m];
			}

			result.push_back(f);
		}

		return result;
	}

	std::vector<bf_tt<N>> evaluateNondegenerateMasks(const std::vector<bf_tt<N>> &inputs) const
	{
		bf_vector<N> result;

		result.reserve(1ULL << _intermediateOutputCount);

		auto fout = evaluate(inputs);

		bf_tt<N> fterminals(init::zero);

		for (auto i : terminalOutputs())
			fterminals ^= fout[i];

		auto intermediates = intermediateOutputs();

		for (uint64_t mask = 0; mask < 1ULL << _intermediateOutputCount; mask++)
		{
			bf_tt<N> f(fterminals);

			for (uint64_t m = 0; m < _intermediateOutputCount; m++)
			{
				if ((mask >> m) & 1)
					f ^= fout[intermediates[m]];
			}

			result.push_back(f);
		}

		return result;
	}


	bf_vector<N> evaluate(const bf_vector<N> &inputs) const
	{
		std::vector<bf_tt<N>> variables(inputCount() + gateCount());

		std::copy_n(std::begin(inputs), inputCount(), std::begin(variables));
		//std::copy(std::begin(inputs), std::end(inputs), std::begin(variables));

		for (auto &line : _slp)
		{
			bf_tt<N> f1(init::zero);
			bf_tt<N> f2(init::zero);

			for (auto e : line.leftInputs)
				f1 ^= variables[e];

			for (auto e : line.rightInputs)
				f2 ^= variables[e];

			f1 &= f2;

			variables[line.outputIndex] = f1;
		}

		//std::cout << "output: " << _outputCount << std::endl;
		//std::cout << "inputs: " << inputs.size() << std::endl;
		//std::cout << "variables: " << _variables.size() << std::endl;

		bf_vector<N> result;
		result.reserve(outputCount());

		for (uint64_t i = inputCount(); i < variables.size(); i++)
			result.push_back(variables[i]);

		return result;
	}

	bf_tt<N> str(const bf_vector<N> &inputs, const bf_tt<N> &ffinal, uint64_t mask, bool printIntermediateValues = false) const
	{
		const uint64_t K = gateCount();

		std::string line;
		std::string topstr;
		std::string repstr;
		bf_vector<N> ai(K);
		bf_vector<N> li(inputs);

		for (int i = 0; i < K; i++)
		{
			std::ostringstream ostr1, ostr2;

			bf_tt<N> L(li[2 * i + 0]);
			bf_tt<N> R(li[2 * i + 1]);

			ostr1 << "a" << i << " = (";

			bool firstTerm{ false };

			for (auto j : _topology.gate(i)._leftInputs)
			{
				L ^= ai[j];

				if (!firstTerm)
					firstTerm = true;
				else
					ostr1 << "+";

				ostr1 << "a" << j;
			}

			if (!li[2 * i + 0].is_zero() || (!firstTerm))
			{
				if (firstTerm)
					ostr1 << "+";

				ostr1 << li[2 * i + 0].anf().str();
			}

			ostr1 << ") * (";

			firstTerm = false;

			for (auto j : _topology.gate(i)._rightInputs)
			{
				R ^= ai[j];

				if (!firstTerm)
					firstTerm = true;
				else
					ostr1 << "+";

				ostr1 << "a" << j;
			}

			if (!li[2 * i + 1].is_zero() || (!firstTerm))
			{
				if (firstTerm)
					ostr1 << "+";

				ostr1 << li[2 * i + 1].anf().str();
			}

			ostr1 << ")";

			ai[i] = L;
			ai[i] &= R;

			if(printIntermediateValues)
				ostr2 << "# a" << i << " = " << ai[i].anf().str();

			std::cout << std::left << std::setw(50) << ostr1.str() << " " << ostr2.str() << std::endl;
		}

		std::cout << "f  = ";

		bf_tt<N> fout(init::zero);
		bool firstTerm{ false };

		for (uint64_t g = 0; g < _gateCount; g++)
		{
			if ((mask >> g) & 1)
			{
				fout ^= ai[g];

				if (!firstTerm)
					firstTerm = true;
				else
					std::cout << "+";

				std::cout << "a" << g;
			}
		}

		if (!ffinal.is_zero())
		{
			fout ^= ffinal;

			std::cout << "+" << ffinal.anf().str();
		}

		if(printIntermediateValues)
			std::cout << " # f = " << fout.anf().str();

		std::cout << "\n" << std::endl;

		return fout;
	}

	void dump_slp() const {

		for (int i = 0; i < _slp.size(); i++)
		{
			std::cout << i << " : ";
			std::cout << "out = " << _slp[i].outputIndex << " ";
			std::cout << "left: ";
			for (auto j : _slp[i].leftInputs)
				std::cout << j << " ";
			std::cout << "right: ";
			for (auto j : _slp[i].rightInputs)
				std::cout << j << " ";

			std::cout << std::endl;
		}
	}

	bf_tt<N> evaluateIntermediateMask(const bf_vector<N> &inputs, uint64_t mask) const
	{
		bf_tt<N> result(init::zero);

		auto fout = evaluate(inputs);

		for (auto i : _terminalOutputs)
			result ^= fout[i];

		for (uint64_t m = 0; m < _intermediateOutputCount; m++)
		{
			if ((mask >> m) & 1)
				result ^= fout[_intermediateOutputs[m]];
		}

		return result;
	}

	bf_tt<N> evaluateMask(const bf_vector<N> &inputs, uint64_t mask) const
	{
		bf_tt<N> result(init::zero);

		auto fout = evaluate(inputs);

		for (uint64_t m = 0; m < outputCount(); m++)
		{
			if ((mask >> m) & 1)
				result ^= fout[m];
		}

		return result;
	}

private:

	void compile()
	{
		auto terminals = _topology.getTerminalGates();

		//auto layers = _topology.getLayers();

		uint64_t inputIndex{ 0 };

		//for (uint64_t l = 0; l < layers.size(); l++)
		{
			for(int g = 0; g < _topology.gateCount(); g++)
			//for (auto g : layers[l])
			{
				Line line;

				const typename Topology::TGate &gate = _topology.gate(g);

				line.leftInputs.push_back(inputIndex++);
				for (auto e : gate._leftInputs)
					line.leftInputs.push_back(inputCount() + e);

				line.rightInputs.push_back(inputIndex++);
				for (auto e : gate._rightInputs)
					line.rightInputs.push_back(inputCount() + e);

				line.outputIndex = inputCount() + g;

				bool isTerminalGate = std::find(std::begin(terminals), std::end(terminals), g) != std::end(terminals);

				if (isTerminalGate)
					_terminalOutputs.push_back(g);
				else
					_intermediateOutputs.push_back(g);

				_slp.push_back(line);
			}
		}
	}

	struct Line
	{
		std::vector<uint64_t> leftInputs;
		std::vector<uint64_t> rightInputs;
		uint64_t outputIndex;
	};

	std::vector<uint64_t>	_terminalOutputs;
	std::vector<uint64_t>	_intermediateOutputs;

	std::vector<Line>		_slp;

	const Topology&_topology;

	uint64_t _gateCount;
	uint64_t _inputCount;
	uint64_t _intermediateOutputCount;
	uint64_t _terminalOutputCount;
	std::vector<uint64_t> _terminalGates;
};

