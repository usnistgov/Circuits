#pragma once

#include <cstdint>
#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <map>
#include <set>
#include "bf.h"


class Topology
{
public:

	static const uint64_t SIZE_HINT = 8;

	enum class add_gate{always, if_normalized, normalize};

	Topology() {}
	~Topology() {}

	Topology(const std::string &str){
		set(str);
	}

	bool operator==(const Topology &rhs) const {

		if (gateCount() != rhs.gateCount())
			return false;

		if (degreeDistribution() != rhs.degreeDistribution())
			return false;

		//if (!checkEquivalenceLayers(rhs))
		//	return false;

		return checkEquivalence(rhs);
	}

	uint64_t gateCount() const {
		return _gates.size();
	}

	bool checkEquivalenceLayers(const Topology& rhs) const {

		auto layers1 = getLayers();
		auto layers2 = rhs.getLayers();

		if (layers1.size() != layers2.size())
			return false;

		for (uint64_t i = 0; i < layers1.size(); i++)
			if (layers1[i].size() != layers2[i].size())
				return false;

		return true;
	}

	std::vector<std::vector<uint64_t>> getLayers() const {

		std::vector < std::vector<uint64_t> > layers;

		layers.reserve(gateCount());

		for (uint64_t i = 0; i < _gates.size(); i++)
		{
			if (!_gates[i].hasLeftInputs() && !_gates[i].hasRightInputs())
			{
				if (layers.empty())
					layers.push_back(std::vector<uint64_t>());

				layers[0].push_back(i);
			}
			else
			{
				int depth{ -1 };

				for (int l = static_cast<int>(layers.size()) - 1; l >= 0; l--)
				{
					for (auto g : _gates[i]._leftInputs)
					{
						if (std::find(std::begin(layers[l]), std::end(layers[l]), g) != std::end(layers[l]))
						{
							depth = l;
							break;
						}
					}

					if (depth == l)
						break;
					
					for (auto g : _gates[i]._rightInputs)
					{
						if (std::find(std::begin(layers[l]), std::end(layers[l]), g) != std::end(layers[l]))
						{
							depth = l;
							break;
						}
					}

					if (depth == l)
						break;
				}

				depth++;

				if (depth >= layers.size())
				{
					layers.push_back(std::vector<uint64_t>());
					layers.back().push_back(i);
				}
				else
					layers[depth].push_back(i);
			}
		}

		return layers;
	}

	uint64_t depth() const {

		std::vector<uint64_t> depths(_gates.size(), 0);

		for (uint64_t i = 0; i < _gates.size(); i++)
		{
			uint64_t leftDepth{ 0 }, rightDepth{ 0 };

			for (uint64_t j = 0; j < _gates[i].leftInputCount(); j++)
			{
				auto d = depths[_gates[i]._leftInputs[j]];

				if (d > leftDepth)
					leftDepth = d;
			}

			for (uint64_t j = 0; j < _gates[i].rightInputCount(); j++)
			{
				auto d = depths[_gates[i]._rightInputs[j]];

				if (d > rightDepth)
					rightDepth = d;
			}

			depths[i] = 1 + std::max(leftDepth, rightDepth);
		}

		return *(std::max_element(std::begin(depths), std::end(depths)));
	}

	void set(const std::string &str) {

		_gates.clear();

		size_t gateStart{ 0 };

		while ((gateStart = str.find("(A", gateStart)) != std::string::npos)
		{
			auto leftStart = str.find_first_of(' ', gateStart);

			leftStart++;

			auto leftEnd = str.find(':', leftStart);

			std::string leftStr = str.substr(leftStart, leftEnd - leftStart - 1);

			auto rightStart = str.find_first_not_of(' ', leftEnd + 1);

			auto rightEnd = str.find_first_of(')', rightStart);

			std::string rightStr = str.substr(rightStart, rightEnd - rightStart);

			gateStart = rightEnd;

			uint64_t left{ 0 };
			uint64_t right{ 0 };

			if (leftStr != "L")
			{
				size_t posStart = 0, posEnd;

				do
				{
					posEnd = leftStr.find_first_of(' ', posStart);

					uint64_t gateIndex = (1ULL << std::atoi(leftStr.substr(posStart, posEnd - posStart).c_str()));

					left ^= gateIndex;

					posStart = posEnd + 1;

				} while (posEnd != std::string::npos);
			}

			if (rightStr != "L")
			{
				size_t posStart = 0, posEnd;

				do
				{
					posEnd = rightStr.find_first_of(' ', posStart);

					uint64_t gateIndex = (1ULL << std::atoi(rightStr.substr(posStart, posEnd - posStart).c_str()));

					right ^= gateIndex;

					posStart = posEnd + 1;

				} while (posEnd != std::string::npos);
			}


			addGate(left, right);
		}
	}

	std::string str() const {

		std::string ret;

		if (_gates.empty())
			ret = "<empty topology>";

		for (uint64_t i = 0; i < _gates.size(); i++)
		{
			ret += "(A" + std::to_string(i);

			if (_gates[i].leftInputCount() == 0)
				ret += " L";

			for (uint64_t j = 0; j < _gates[i].leftInputCount(); j++)
				ret += " " + std::to_string(_gates[i]._leftInputs[j]);

			ret += " :";
			
			if (!_gates[i].hasRightInputs())
				ret += " L";

			for (uint64_t j = 0; j < _gates[i].rightInputCount(); j++)
				ret += " " + std::to_string(_gates[i]._rightInputs[j]);

			ret += ')';

			if (i != (_gates.size() - 1))
				ret += ' ';
		}

		return ret;
	}

	bool addGate(uint64_t leftMask, uint64_t rightMask, add_gate ag = add_gate::always) {

		//if (!checkInclusion(leftMask, rightMask))
		//	return false;

		//if (!checkSymmetry(leftMask, rightMask))
		//	return false;

		if (ag != add_gate::always)
		{
			uint64_t normLeft{ leftMask }, normRight{ rightMask };

			normalizeMask(normLeft, normRight);

			if (ag == add_gate::if_normalized &&
				((normLeft != leftMask) || (normRight != rightMask)))
				return false;

			leftMask = normLeft;
			rightMask = normRight;
		}

		_gates.push_back(TGate(leftMask, rightMask));

		return true;
	}

	void insertGateFront(uint64_t leftMask, uint64_t rightMask) {

		for (uint64_t i = 0; i < _gates.size(); i++)
		{
			for (auto &x : _gates[i]._leftInputs)
				x++;

			if ((leftMask >> i) & 1)
				_gates[i]._leftInputs.insert(_gates[i]._leftInputs.begin(), 0);

			for (auto &x : _gates[i]._rightInputs)
				x++;

			if ((rightMask >> i) & 1)
				_gates[i]._rightInputs.insert(_gates[i]._rightInputs.begin(), 0);
		}

		_gates.insert(_gates.begin(), TGate());
	}

	struct TGate {

		enum Input : bool { Left, Right };

		TGate() {
			_leftInputs.reserve(SIZE_HINT);
			_rightInputs.reserve(SIZE_HINT);
		}

		TGate(uint64_t leftMask, uint64_t rightMask) : TGate() {

			for (uint64_t i = 0; i < word_bits<uint64_t>(); i++)
			{
				if ((leftMask >> i) & 1)
					addInput(i, TGate::Left);

				if ((rightMask >> i) & 1)
					addInput(i, TGate::Right);
			}
		}

		uint64_t getMask(Input inp) const {

			const std::vector<uint64_t> &inputs = (inp == Input::Left) ? _leftInputs : _rightInputs;

			uint64_t mask{ 0 };

			for (auto x : inputs)
				mask |= 1ULL << x;

			return mask;
		}

		bool addInput(uint64_t index, Input inp) {

			std::vector<uint64_t> &inputs = (inp == Input::Left) ? _leftInputs : _rightInputs;

			if (std::find(std::begin(inputs), std::end(inputs), index) != std::end(inputs))
				return false;

			inputs.push_back(index);

			return true;
		}

		bool inputContains(uint64_t index, Input inp) const {

			const std::vector<uint64_t> &inputs = (inp == Input::Left) ? _leftInputs : _rightInputs;

			return std::find(std::begin(inputs), std::end(inputs), index) != std::end(inputs);
		}

		uint64_t leftInputCount() const {
			return _leftInputs.size();
		}

		uint64_t rightInputCount() const {
			return _rightInputs.size();
		}

		bool hasLeftInputs() const {
			return leftInputCount() != 0;
		}

		bool hasRightInputs() const {
			return rightInputCount() != 0;
		}

		std::vector<uint64_t> _leftInputs;
		std::vector<uint64_t> _rightInputs;
	};

	const TGate& gate(uint64_t index) const {
		return _gates[index];
	}

	std::vector<uint64_t> getTerminalGates() const {

		std::vector<uint64_t> result;

		std::vector<uint8_t> used(gateCount(), 0);

		for (auto &g : _gates)
		{
			for (auto e : g._leftInputs)
				used[e] = 1;

			for (auto e : g._rightInputs)
				used[e] = 1;
		}

		for (uint64_t i = 0; i < used.size(); i++)
			if (!used[i])
				result.push_back(i);

		return result;
	}

	std::vector<uint64_t> getIntermediateGates() const {

		auto terminals = getTerminalGates();

		std::vector<uint64_t> intermediateGates;

		for (uint64_t i = 0; i < _gates.size(); i++)
			if (std::find(std::begin(terminals), std::end(terminals), i) == std::end(terminals))
				intermediateGates.push_back(i);

		return intermediateGates;
	}

	uint64_t xorCount() const {

		uint64_t count{ 0 };

		for (auto &g : _gates)
		{
			if (g.leftInputCount() > 1)
				count += g.leftInputCount() - 1;

			if (g.rightInputCount() > 1)
				count += g.rightInputCount() - 1;
		}

		return count;
	}

	bool reflect(uint64_t inputIndex) {

		auto gateIndex = inputIndex / 2;

		auto &g = gate(gateIndex);

		auto &inputsOther = (inputIndex % 2) ? g._leftInputs : g._rightInputs;

		if (inputsOther.empty())
			return false;

		for (auto &gt : _gates)
		{
			if (std::find(std::begin(gt._leftInputs), std::end(gt._leftInputs), gateIndex) != std::end(gt._leftInputs))
			{
				for (auto x : inputsOther)
				{
					auto pos = std::find(std::begin(gt._leftInputs), std::end(gt._leftInputs), x);
					
					if (pos == std::end(gt._leftInputs))
						gt._leftInputs.push_back(x);
					else
						gt._leftInputs.erase(pos);
				}

				std::sort(std::begin(gt._leftInputs), std::end(gt._leftInputs));
			}

			if (std::find(std::begin(gt._rightInputs), std::end(gt._rightInputs), gateIndex) != std::end(gt._rightInputs))
			{
				for (auto x : inputsOther)
				{
					auto pos = std::find(std::begin(gt._rightInputs), std::end(gt._rightInputs), x);

					if (pos == std::end(gt._rightInputs))
						gt._rightInputs.push_back(x);
					else
						gt._rightInputs.erase(pos);
				}

				std::sort(std::begin(gt._rightInputs), std::end(gt._rightInputs));
			}

			uint64_t lm = gt.getMask(TGate::Input::Left);
			uint64_t rm = gt.getMask(TGate::Input::Right);

			normalizeMask(lm, rm);

			gt = TGate(lm, rm);
		}

		return true;
	}

	std::vector<uint64_t> getLinearInputIndices() const {

		std::vector<uint64_t> vec;

		for (uint64_t i = 0; i < gateCount(); i++)
		{
			if (!_gates[i].hasLeftInputs())
				vec.push_back(2 * i + 0);

			if (!_gates[i].hasRightInputs())
				vec.push_back(2 * i + 1);
		}

		return vec;
	}

	void write_graph(const std::string &filename) {

		std::ofstream ofile(filename);

		if (!ofile)
			return;

		ofile << R"(digraph  "Topology")" << std::endl;
		ofile << "{" << std::endl;

		for (uint64_t i = 0; i < _gates.size(); i++)
		{
			ofile << "A" << std::to_string(i) << ' ';
			ofile << R"([shape = circle, color = blue3, style = filled, label = "", height = .3, width = .3])" << std::endl;
		}

		for (uint64_t i = 0; i < _gates.size(); i++)
		{
			if (_gates[i].hasLeftInputs())
			{
				for (auto x : _gates[i]._leftInputs)
					ofile << "A" << x << " -> " << "A" << i << std::endl;
			}
		}

		ofile << "}" << std::endl;
	}

private:

	std::map<uint64_t, uint64_t> degreeDistribution() const {

		std::map<uint64_t, uint64_t> dist;

		for (auto &g : _gates)
		{
			uint64_t leftInputs = g.leftInputCount();
			uint64_t rightInputs = g.rightInputCount();

			uint64_t key{ 0 };

			if (leftInputs <= rightInputs)
				key = (leftInputs << 32) | rightInputs;
			else
				key = (rightInputs << 32) | leftInputs;

			dist[key]++;
		}

		return dist;
	}

	uint64_t concatMasks(uint64_t leftMask, uint64_t rightMask) {

		return (leftMask << 32) | rightMask;
	}

	bool checkInclusion(uint64_t leftMask, uint64_t rightMask)
	{
		// eliminate forms of (f+L1)*(f+g+L2)
		uint64_t intersection = leftMask & rightMask;
		if (intersection != 0 && (intersection == leftMask || intersection == rightMask))
			return false;

		return true;
	}

	void normalizeMask(uint64_t &leftMask, uint64_t &rightMask)
	{
		// choose the minimum from the set {a.b, a.(a+b), (a+b).b}

		uint64_t masks[] = {
			concatMasks(leftMask, rightMask),
			concatMasks(leftMask, leftMask ^ rightMask),
			concatMasks(leftMask ^ rightMask, rightMask),
			concatMasks(rightMask, leftMask),
			concatMasks(leftMask ^ rightMask, leftMask),
			concatMasks(rightMask, leftMask ^ rightMask)
		};

		auto minMask = *(std::min_element(std::begin(masks), std::end(masks)));

		leftMask = minMask >> 32;
		rightMask = minMask & 0xffffffff;
	}

	bool checkSymmetry(uint64_t leftMask, uint64_t rightMask)
	{
		// choose the minimum from the set {a.b, a.(a+b), (a+b).b}
		uint64_t form1 = concatMasks(leftMask, rightMask);
		uint64_t form2 = concatMasks(leftMask, leftMask ^ rightMask);
		uint64_t form3 = concatMasks(leftMask ^ rightMask, rightMask);

		uint64_t form4 = concatMasks(rightMask, leftMask);
		uint64_t form5 = concatMasks(leftMask ^ rightMask, leftMask);
		uint64_t form6 = concatMasks(rightMask, leftMask ^ rightMask);

		if (form1 > form2 ||
			form1 > form3 ||
			form1 > form4 ||
			form1 > form5 ||
			form1 > form6)
			return false;

		return true;
	}

	bool checkEquivalence(const Topology &rhs) const {

		std::vector<uint64_t> perm(gateCount());

		std::iota(std::begin(perm), std::end(perm), 0);

		do
		{
			if (verifyMapping(perm, rhs))
				return true;

		} while (std::next_permutation(std::begin(perm), std::end(perm)));

		return false;
	}

	bool verifyMapping(const std::vector<uint64_t> &perm, const Topology &rhs) const {

		for (uint64_t i = 0; i < perm.size(); i++)
		{
			const TGate& gate1 = gate(i);
			const TGate& gate2 = rhs.gate(perm[i]);

			if (gate1.leftInputCount() == gate2.leftInputCount() && gate1.rightInputCount() == gate2.rightInputCount())
			{
				bool gateMaps{ true };

				for (auto e : gate1._leftInputs)
					if (!gate2.inputContains(perm[e], TGate::Left))
					{
						gateMaps = false;
						break;
					}

				for (auto e : gate1._rightInputs)
					if (!gate2.inputContains(perm[e], TGate::Right))
					{
						gateMaps = false;
						break;
					}

				if (gateMaps)
					continue;
			}

			if (gate1.leftInputCount() == gate2.rightInputCount() && gate1.rightInputCount() == gate2.leftInputCount())
			{
				bool gateMaps{ true };

				for (auto e : gate1._leftInputs)
					if (!gate2.inputContains(perm[e], TGate::Right))
					{
						gateMaps = false;
						break;
					}

				for (auto e : gate1._rightInputs)
					if (!gate2.inputContains(perm[e], TGate::Left))
					{
						gateMaps = false;
						break;
					}

				if (gateMaps)
					continue;
			}

			return false;
		}

		return true;
	}

	std::vector<TGate> _gates;
};

using TopologyVector = std::vector<Topology>;

TopologyVector load_topologies(const std::string &fileName);

TopologyVector load_topologies(uint64_t K);
