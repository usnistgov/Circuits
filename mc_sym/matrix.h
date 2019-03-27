#pragma once

#include <vector>
#include <array>
#include <set>
#include <cstdint>
#include <cassert>
#include <iostream>
#include "logger.h"
#include "utils.h"
#include "constant_time.h"

template<uint64_t Rows, uint64_t Cols = Rows, typename T = uint64_t>
class matrix
{
public:

	static const uint64_t WordBits = 8 * sizeof(T);

	static_assert(Cols <= WordBits, "number of columns too large");

	static const uint64_t ColMask = partial_mask<T, Cols>();

	matrix(){
	}

	matrix(std::initializer_list<T> elem) {
		assert(elem.size() <= Rows);
		std::copy(elem.begin(), elem.end(), _rows.begin());
	}

	~matrix(){}

	matrix(const matrix<Rows, Cols, T> &rhs) {
		_rows = rhs._rows;
	}

	T& operator[](uint64_t index) {
		assert(index < Rows);
		return _rows[index];
	}

	const T& operator[](uint64_t index) const {
		assert(index < Rows);
		return _rows[index];
	}

	bool operator==(const matrix<Rows, Cols, T>& rhs) {
		return _rows == rhs._rows;
	}

	void clear() {
		_rows.fill(0);
	}

	void complement() {
		for (auto &r : _rows)
			r = (~r) & ColMask;
	}

	T getRow(uint64_t index) const {
		assert(index < Rows);
		return _rows[index];
	}

	void setRow(uint64_t index, T value) {
		assert(index < Rows);
		assert(value < power_of_2(Cols));
		_rows[index] = value;
	}

	T getColumn(uint64_t index) const {

		assert(index < Cols);
		assert(Rows <= WordBits);

		T value{ 0 };

		for (uint64_t i = 0; i < Rows; i++)
		{
			value <<= 1;
			value ^= get_bit(_rows[Cols - 1 - i], index);
		}

		return value;
	}

	T entry(uint64_t r, uint64_t c) const {

		return static_cast<T>(get_bit(_rows, c));
	}

	void setColumn(uint64_t index, T value) {

		for (uint64_t i = 0; i < Rows; i++)
		{
			clear_bit(_rows[i], index);
			_rows[i] ^= get_bit(value, i) << index;
		}
	}

	uint64_t row_weight(uint64_t index) const {

		return popcnt(getRow(index));
	}

	uint64_t col_weight(uint64_t index) const {

		return popcnt(getColumn(index));
	}

	uint64_t col_common_terms(uint64_t col1, uint64_t col2) const {

		return popcnt(getColumn(col1) & getColumn(col2));
	}

	void print() const {

		for (int i = 0; i < Rows; i++)
		{
			for (int j = Cols - 1; j >= 0; j--)
				std::cout << ((_rows[i] >> j) & 1) << " ";

			std::cout << std::endl;
		}

		std::cout << std::endl;
	}

	static matrix<Rows, Cols, T> identity() {

		matrix<Rows, Cols, T> I;

		for (uint64_t i = 0; i < Rows; i++)
			I.setRow(i, unit_vector<T>(Cols - 1 - i));

		return I;
	}

	static matrix<Rows, Cols, T> ae_identity() {

		matrix<Rows, Cols, T> I;

		for (uint64_t i = 0; i < Rows; i++)
			I.setRow(i, unit_vector<T>(i));

		return I;
	}

protected:
	std::array<T, Rows> _rows;
};

class dmatrix
{
public:

	dmatrix(uint64_t cols = 0) 
		: _cols(cols) {
	}

	dmatrix(uint64_t rows, uint64_t cols) 
		: _cols(cols), _data(rows) {
	}

	uint64_t& operator[](uint64_t index) {
		return _data[index];
	}

	const uint64_t& operator[](uint64_t index) const {
		return _data[index];
	}

	uint64_t rows() const {
		return _data.size();
	}

	uint64_t cols() const {
		return _cols;
	}

	void set_rows(uint64_t rows) {
		_data.resize(rows);
	}

	void add_row(uint64_t r) {
		_data.push_back(r);
	}

private:
	uint64_t _cols;
	std::vector<uint64_t> _data;
};


uint64_t make_row_reduced(dmatrix &mat);

uint64_t rank(const dmatrix &mat);

template<uint64_t N, uint64_t M, typename T = uint64_t>
uint64_t make_row_reduced(matrix<N, M, T> &m)
{
	//auto col_set_mask = full_mask<uint64_t>();
	//for (uint64_t r = 0; r < N; r++)
	//	col_set_mask |= m[r];

	uint64_t rank{ 0 };

	for (uint64_t c = 0; c < M; c++)
	{
		//if (!get_bit(col_set_mask, c))
		//	continue;

		for (uint64_t r = rank; r < N; r++)
		{
			if (get_bit(m[r], c))
			{
				//ct::swap_if<uint64_t>(m[r], m[rank], r != rank);
				//if (r != rank)
					std::swap(m[r], m[rank]);

				for (uint64_t j = r + 1; j < N; j++)
//					ct::xor_if(m[j], m[rank], get_bit(m[j], c));
				if (get_bit(m[j], c))
					m[j] ^= m[rank];

				rank++;

				break;
			}
		}
	}

	return rank;
}

template<uint64_t N, uint64_t M, typename T = uint64_t>
uint64_t rank(const matrix<N, M, T> &m)
{
	auto rr(m);

	auto rank = make_row_reduced(rr);

	return rank;
}

template<uint64_t N, uint64_t M, typename T = uint64_t>
uint64_t is_invertible(const matrix<N, M, T> &m)
{
	return rank(m) == N;
}


template<uint64_t N>
class invertible_matrix_generator
{
public:

	invertible_matrix_generator() {
		for (uint64_t i = 0; i < N; i++)
			_m.setRow(i, 1ULL << i);
	}

	bool next()
	{
		if (nextVector(N - 1))
			return true;
		else
		{
			int64_t i;
			for (i = N - 2; i >= 0; i--)
			{
				calculateSpan(i - 1);

				if (nextVector(i))
				{
					for (int64_t j = i; j < N - 1; j++)
					{
						//std::set<u64> newset;

						//for (auto elem : _span)
						//{
						//	newset.insert(elem);
						//	newset.insert(elem ^ _matrix[j]);
						//}

						//_span.swap(newset);

						calculateSpan(j);

						_m[j + 1] = 0;
						nextVector(j + 1);
					}

					break;
				}
			}

			if (i == -1 && _m[0] == (1ULL << N))
				return false;
		}

		return true;
	}

	void print()
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = N - 1; j >= 0; j--)
				std::cout << ((_m[i] >> j) & 1) << " ";

			std::cout << std::endl;
		}
	}

	const matrix<N>& get_matrix() const {
		return _m;
	}

	uint64_t size() const {

		uint64_t s{ 1 };

		for (uint64_t i = 0; i < N; i++)
			s *= power_of_2(N) - power_of_2(i);

		return s;
	}

private:

	void calculateSpan(uint64_t rows)
	{
		_span.clear();

		for (uint64_t j = 0; j < (1ULL << (rows + 1)); j++)
		{
			uint64_t vec = 0;

			for (uint64_t k = 0; k < (rows + 1); k++)
				if ((j >> k) & 1)
					vec ^= _m[k];

			_span.insert(vec);
		}
	}

	bool nextVector(uint64_t row)
	{
		do
		{
			_m[row]++;

			if (_m[row] == (1ULL << N))
				return false;

		} while (_span.count(_m[row]) > 0);

		return true;
	}

	matrix<N> _m;
	std::set<uint64_t> _span;
};
