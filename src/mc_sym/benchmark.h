#pragma once

#include <string>
#include <iostream>
#include <iomanip>
#include <cstdint>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>

namespace {

	static uint64_t take_time_stamp() {
		return std::chrono::high_resolution_clock::now().time_since_epoch().count();
	}

	static double elapsed(uint64_t since) {
		return (take_time_stamp() - since) * 1e-9;
	}
}

//! benchmarks the execution time of a structured code block
struct benchmark
{
	enum class report {none, begin, end, all};

	benchmark(const std::string &caption = "", report r = report::all, uint64_t iterations = 1) 
		: caption_(caption), report_(r), num_iterations_(iterations)
	{		
		if ((report_ == report::begin) || (report_ == report::all))
			print("started");

		start_ = take_time_stamp();
	}

	~benchmark() 
	{
		if ((report_ == report::end) || (report_ == report::all))	{

			double t = elapsed(start_);
			print("finished " + format_elapsed(t));
		}
	}

	void checkpoint(const std::string &name) const {

		double t = elapsed(start_);
		print("checkpoint '" + name + "' " + format_elapsed(t));
	}

	bool is_done() const {
		return (current_iteration_ >= num_iterations_);
	}

	bool is_done(double& e) {
		e = elapsed(start_);
		return is_done();
	}

	void next() {
		++current_iteration_;
	}

private:


	static std::string format_elapsed(double elapsed) {

		std::stringstream ss;

		ss << "[ " << std::fixed << std::setprecision(Precision) << elapsed << " sec ]";

		return ss.str();
	}

	void print(const std::string &str) const {

		std::cout << "benchmark '" << caption_ << "' " << str << std::endl;
	}

	static const int Precision = 4;

	std::string caption_;
	report report_;
	uint64_t start_;
	uint64_t num_iterations_;
	uint64_t current_iteration_{ 0 };
};

#define BENCHMARK(name) for(benchmark _bm(name); !_bm.is_done(); _bm.next())  
#define BENCHMARK_N(name, count) for(benchmark _bm(name, benchmark::report::end, count); !_bm.is_done(); _bm.next())
#define BENCHMARK_SAVE(name, elapsed) for(benchmark _bm(name, benchmark::report::none); !_bm.is_done(elapsed); _bm.next())

//! prints the progression of a computation, showing the percentage completed
struct progress_bar
{
public:

	progress_bar(const std::string &title = "", uint64_t units = 100, bool silent = false) 
		: _title(title), _units(units) , _silent(silent) {
	}

	~progress_bar() {

		if (!_silent)
		{
			std::cout << '\r';
			for (uint64_t i = 0; i < _title.length() + BarLength + 10; i++)
				std::cout << ' ';
			std::cout << '\r';
		}
	}

	bool set(uint64_t completed) {

		_completed = completed;

		auto newpos = get_pos(_completed);

		if (newpos != _currentpos)
		{
			_currentpos = newpos;

			showbar();

			return true;
		}

		return false;
	}

	void update(uint64_t units = 1) {

		set(_completed + units);
	}

	uint64_t current_pos() const {
		return _currentpos;
	}

private:

	uint64_t get_pos(uint64_t completed) const {

		return std::min(100ULL, (100ULL * completed) / _units);
	}

	void showbar() const {

		if (_silent)
			return;

		std::cout << '\r' << _title << ' ';

		auto pos = (_currentpos * BarLength) / 100;

		for (uint64_t x = 0; x <= BarLength; x++)
		{
			if (x <= pos)
				std::cout << FullCell;
			else
				std::cout << EmptyCell;
		}

		std::cout << ' ' << _currentpos << '%';
	}


	static const uint64_t BarLength = 20;

	static const uint8_t EmptyCell = 177;

	static const uint8_t FullCell = 219;

	std::string _title;
	uint64_t _units;
	uint64_t _completed = 0;
	uint64_t _currentpos = 0xff;
	bool _silent;
};
