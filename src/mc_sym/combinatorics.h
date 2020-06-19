#pragma once
#include <cstdint>
#include <vector>
#include <set>
#include <memory>
#include <iostream>
#include <numeric>
#include <cassert>

#include "logger.h"

namespace combinatorics {

	uint64_t choose(uint64_t n, uint64_t k);

	uint64_t factorial(uint64_t n);

	class subset
	{
	public:

		class Iterator : public std::iterator<std::forward_iterator_tag, const std::vector<uint64_t>>
		{
		public:

			Iterator(subset &s, uint64_t index) : s_(s), index_(index) {
			}

			bool operator==(const Iterator &rhs) const {
				return index_ == rhs.index_;
			}

			bool operator!=(const Iterator &rhs) const {
				return index_ != rhs.index_;
			}

			const std::vector<uint64_t>& operator*() const {
				return s_();
			}

			Iterator operator++() {
				s_.next();
				index_++;
				return *this;
			}

		private:
			subset & s_;
			uint64_t index_;
		};

		subset(uint64_t objectCount, uint64_t chooseCount);
		~subset();

		Iterator begin() {
			return Iterator(*this, _stateIndex);
		}

		Iterator end() {
			return Iterator(*this, _terminateIndex);
		}

		uint64_t items() const {
			return _chooseCount;
		}

		void set_terminate_index(uint64_t index) {
			_terminateIndex = index;
		}

		virtual bool next();

		uint64_t count() const;

		uint64_t state_index() const {
			return _stateIndex;
		}

		uint64_t operator[](uint64_t index) const {
			return _items[index];
		}

		const std::vector<uint64_t>& operator()() const {
			return _items;
		}

		virtual void init();

		uint64_t size() const {
			return choose(_objectCount, _chooseCount);
		}

		bool seek(uint64_t state_index) {

			if (state_index >= size())
				return false;

			uint64_t rank{ state_index };
			uint64_t placed{ 0 };
			for (uint64_t i = 0; (i < _objectCount) && (placed != _chooseCount); i++)
			{
				auto c = choose(_objectCount - i - 1, _chooseCount - placed - 1);

				// check if item i is contained in the selection
				if (rank < c)
					_items[placed++] = i;
				else
					rank -= c;
			}

			assert(placed == _chooseCount);

			_stateIndex = state_index;
			
			return true;
		}

	protected:

		const uint64_t _objectCount;
		const uint64_t _chooseCount;
		uint64_t _stateIndex;
		uint64_t _terminateIndex;

		std::vector<uint64_t> _items;
	};

	class enumeration {

	public :
		enumeration(uint64_t items) : items_(items), limits_(items), state_(items) {

			init_state();
		}

		enumeration(std::initializer_list<uint64_t> il) : enumeration(il.size()) {

			set_limits(il);
		}

		void set_limit(uint64_t index, uint64_t value) {

			limits_[index] = value;
		}

		void set_limits(std::initializer_list<uint64_t> il) {

			std::copy(std::begin(il), std::end(il), std::begin(limits_));
		}

		uint64_t operator[](uint64_t index) const {

			return state_[index];
		}

		uint64_t size() const {

			return items_;
		}

		uint64_t count() const {
			
			uint64_t count{ 1 };
			
			for (auto x : limits_)
				count *= x;

			return count;
		}

		const std::vector<uint64_t>& state() const {

			return state_;
		}

		bool next() {
			
			state_[0]++;

			for (uint64_t i = 0; i < items_ - 1; i++)
			{
				if (state_[i] >= limits_[i])
				{
					for (uint64_t j = 0; j <= i; j++)
						state_[i] = 0;

					state_[i + 1]++;
				}
			}

			if (state_.back() >= limits_.back())
			{
				init_state();
				return false;
			}

			return true;
		}

		struct enum_iterator
		{
			enum_iterator(enumeration &en, bool end = false) : en_(en), end_(end) {}

			bool operator!=(const enum_iterator& rhs) {
				return end_ != rhs.end_;
			}

			const std::vector<uint64_t>& operator*() {
				return en_.state();
			}

			void operator++() {
				if (!en_.next())
					end_ = true;
			}

			enumeration& en_;
			bool end_;
		};

		enum_iterator begin() {
			return enum_iterator(*this);
		}

		enum_iterator end() {
			return enum_iterator(*this, true);
		}

	private:

		void init_state() {
			std::fill(std::begin(state_), std::end(state_), 0);
		}

		uint64_t items_;
		std::vector<uint64_t> limits_;
		std::vector<uint64_t> state_;
	};

	class permutation
	{
	public:

		permutation(uint64_t uItems);
		
		~permutation() = default;

		bool next();

	protected:

		void init();

		std::vector<uint64_t> _state;
	};

	
	class partitioned_subset
	{
	public:
		partitioned_subset(uint64_t object_count, uint64_t choose_count, uint64_t partition_count) {

			_partitions.reserve(partition_count);

			subset sub(object_count, choose_count);

			auto total_combins = choose(object_count, choose_count);

			auto partition_size = (total_combins + partition_count - 1) / partition_count;

			LOG_INFO << "partitioning C(" << object_count << ',' << choose_count << ") = " << total_combins << " into " << partition_count << " parts of size " << partition_size << LOG_ENDL;

			_partitions.push_back(sub);

			for (uint64_t i = 1; i < partition_count; i++)
			{
				for (uint64_t j = 0; (j < partition_size); j++)
					sub.next();

				if (sub.state_index() < _partitions[i - 1].state_index())
					break;

				_partitions.push_back(sub);
				_partitions[i - 1].set_terminate_index(sub.state_index());
			}

			auto generated_combins = std::accumulate(_partitions.begin(), _partitions.end(),
				0ULL,
				[](uint64_t a, subset &s) {
				return a + s.count();
			});

			LOG_VERIFY(generated_combins == total_combins);

			LOG_DEBUG << "partition sizes : ";
			for (auto &part : _partitions)
				LOG_DEBUG << part.count() << ' ';
			LOG_DEBUG << LOG_ENDL;
		}

		uint64_t size() const {
			return _partitions.size();
		}

		subset& operator[](uint64_t index) {
			return _partitions[index];
		}

		const subset& operator[](uint64_t index) const {
			return _partitions[index];
		}

	private:
		std::vector<subset> _partitions;
	};

	
}