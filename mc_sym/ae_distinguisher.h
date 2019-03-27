#pragma once

#include <cstdint>
#include <vector>
#include <set>

#include "ae_group.h"
#include "ae_helper.h"
#include "ae_local_connections.h"
#include "ae_indicators.h"
#include "bf.h"
#include "benchmark.h"

namespace ae {

	enum class distinguisher_t : uint32_t { local_connection, indicator };


	//! Checks whether the signatures of the local connections of two functions are equal
	//! Note: The signatures being equal does not mean the functions are affine equivalent
	template<uint64_t N>
	bool is_affine_equivalent_lc(const bf_tt<N> &f1, const bf_tt<N> &f2, const uint32_t max_distance)
	{
		assert(max_distance > 0);

		for (uint32_t distance = 1; distance <= max_distance; distance++)
		{
			auto ns = compute_neighbor_signatures_lc(f1, distance);

			for (const auto f : local_connections(f2, distance))
				if (ns.signatures.count(calculate_affine_invariants(f)) == 0)
					return false;
		}

		return true;
	}

	//! Checks whether the signatures of the indicators of two functions are equal
	//! Note: The signatures being equal does not mean the functions are affine equivalent
	template<uint64_t N>
	bool is_affine_equivalent_ind(const bf_tt<N> &f1, const bf_tt<N> &f2, const uint32_t max_distance)
	{
		for (uint32_t distance = 1; distance <= max_distance; distance++)
		{
			auto ns = compute_neighbor_signatures_ind(f1, distance);

			for (auto f : get_indicators(f2, distance))
				if (ns.signatures.count(calculate_affine_invariants(f)) == 1)
					return false;
		}

		return true;
	}


	struct neighbor_groups
	{
		uint64_t distance;
		std::vector<std::vector<uint64_t>> groupIndices;
	};

	struct ae_distinguisher
	{
		distinguisher_t type;
		std::vector<neighbor_groups> data;
	};


	template<uint64_t N>
	std::vector<uint64_t> compute_neighbor_groups_lc(const group_vector<N> &groupdata, const bf_tt<N> &f, uint64_t distance)
	{
		auto connections = get_local_connections(f, distance);

		auto indices = get_group_indices(connections, groupdata);

		return { indices.begin(), indices.end() };
	}

	template<uint64_t N>
	std::vector<std::vector<uint64_t>> compute_neighbor_groups_lc(const group_vector<N> &groupdata, const bf_vector<N> &fvec, uint64_t distance)
	{
		std::vector<std::vector<uint64_t>> vec;

		vec.reserve(fvec.size());

		for (auto f : fvec)
			vec.push_back(compute_neighbor_groups_lc(groupdata, f, distance));

		return vec;
	}


	template<uint64_t N>
	std::vector<uint64_t> compute_neighbor_groups_ind(const group_vector<N> &groupdata, const bf_tt<N> &f, uint64_t distance)
	{
		bf_set<N> indicators;

		bf_vector<N> indicator_vec;

		if (distance == 0)
			indicator_vec = get_indicators(f);
		else
		{
			for (auto g : local_connections(f, distance))
			{
				auto ind = get_indicators(g);

				for (auto e : ind)
					indicators.insert(e);
			}

			std::copy(std::begin(indicators), std::end(indicators), std::back_inserter(indicator_vec));
		}

		auto indices = get_group_indices(indicator_vec, groupdata);

		return { indices.begin(), indices.end() };
	}

	template<uint64_t N>
	std::vector<std::vector<uint64_t>> compute_neighbor_groups_ind(const group_vector<N> &groupdata, const bf_vector<N> &fvec, uint64_t distance)
	{
		std::vector<std::vector<uint64_t>> vec;

		vec.reserve(fvec.size());

		for (auto f : fvec)
			vec.push_back(compute_neighbor_groups_ind(groupdata, f, distance));

		return vec;
	}


	template<uint64_t N>
	bfl::bf_tt<N> find_representative_ind(const bfl::bf_tt<N> &f, const ae::group<N> &group, const ae::group_vector<N> &groups, const ae_distinguisher &aed)
	{
		assert(aed.type == distinguisher_t::indicator);

		const uint64_t functionsInGroup = group.functions.size();

		elimination_list el(functionsInGroup);

		for (uint64_t distIndex = 0; distIndex < aed.data.size(); distIndex++)
		{
			auto distance = aed.data[distIndex].distance;

			bfl::bf_vector<N> indicators;

			if (distance == 0)
				indicators = get_indicators(f);
			else
			{
				for (auto g : local_connections(f, distance))
				{
					auto tempind = get_indicators(g);
					std::copy(std::begin(tempind), std::end(tempind), std::back_inserter(indicators));
				}
			}

			for (auto &fi : indicators)
			{
				auto index = get_group_index(fi, groups);

				for (uint64_t i = 0; i < el.count(); i++)
				{
					if (!std::binary_search(std::begin(aed.data[distIndex].groupIndices[el.get(i)]), std::end(aed.data[distIndex].groupIndices[el.get(i)]), index))
					{
						if (el.remove(i))
							return group.functions[el.head()];
						else
							i--;
					}
				}
			}
		}

		const uint64_t lastDistance = aed.data.size() - 1;

		uint64_t min = aed.data[lastDistance].groupIndices[el.get(0)].size();
		uint64_t minIndex = 0;

		for (uint64_t i = 1; i < el.count(); i++)
			if (aed.data[lastDistance].groupIndices[el.get(i)].size() < min)
			{
				min = aed.data[lastDistance].groupIndices[el.get(i)].size();
				minIndex = i;
			}

		return group.functions[el.get(minIndex)];
	}

	template<uint64_t N>
	bfl::bf_tt<N> find_representative_lc(const bfl::bf_tt<N> &f, const ae::group<N> &group, const ae::group_vector<N> &groups, const ae_distinguisher &aed)
	{
		assert(aed.type == distinguisher_t::local_connection);

		const uint64_t functionsInGroup = group.functions.size();

		elimination_list el(functionsInGroup);

		for (uint64_t distIndex = 0; distIndex < aed.data.size(); distIndex++)
		{
			auto distance = aed.data[distIndex].distance;

			for (auto g : local_connections(f, distance))
			{
				auto index = get_group_index(g, groups);

				for (uint64_t i = 0; i < el.count(); i++)
				{
					if (!std::binary_search(std::begin(aed.data[distIndex].groupIndices[el.get(i)]), std::end(aed.data[distIndex].groupIndices[el.get(i)]), index))
					{
						if (el.remove(i))
							return group.functions[el.head()];
						else
							i--;
					}
				}
			}
		}

		const uint64_t lastDistance = aed.data.size() - 1;

		uint64_t min = aed.data[lastDistance].groupIndices[el.get(0)].size();
		uint64_t minIndex = 0;

		for (uint64_t i = 1; i < el.count(); i++)
			if (aed.data[lastDistance].groupIndices[el.get(i)].size() < min)
			{
				min = aed.data[lastDistance].groupIndices[el.get(i)].size();
				minIndex = i;
			}

		return group.functions[el.get(minIndex)];
	}


	template<uint64_t N>
	class class_distinguisher_base
	{
	public:

		class_distinguisher_base() {
		}

		class_distinguisher_base(const bf_vector<N> &reps) : reps_(reps) {
		}

		virtual bf_tt<N> get_representative(const bf_tt<N> &f) const = 0;

		virtual uint64_t get_representative_index(const bf_tt<N> &f) const {
			auto rep = get_representative(f);

			auto pos = std::lower_bound(std::begin(this->reps_), std::end(this->reps_), rep);

			assert(pos != std::end(this->reps_));

			return std::distance(std::begin(this->reps_), pos);
		}

		const bf_vector<N>& representatives() const {
			return reps_;
		}

		const group_vector<N>& groups() const {
			return groups_;
		}

	protected:

		typename group_vector<N>::const_iterator group_of(const bf_tt<N> &f) const {

			auto ai = calculate_affine_invariants(f);

			auto p = std::lower_bound(std::begin(groups_), std::end(groups_), group<N>(ai));

			assert(p != std::end(groups_) && p->ai == ai);

			return p;
		}

		bf_vector<N> reps_;
		group_vector<N> groups_;
	};

	template<uint64_t N>
	class class_distinguisher_basic : public class_distinguisher_base<N>
	{
	public:
		class_distinguisher_basic(const bf_vector<N> &reps) : class_distinguisher_base<N>(reps) {
			std::sort(std::begin(this->reps_), std::end(this->reps_));

			this->groups_ = generate_groups_fast(this->reps_);

			assert(std::all_of(std::begin(this->groups_), std::end(this->groups_), [](const group<N> &gr) {
				return gr.functions.size() == 1; }));
		}

		bf_tt<N> get_representative(const bf_tt<N> &f) const override {
			return this->group_of(f)->functions[0];
		}
	};

	template<uint64_t N>
	class class_distinguisher_advanced : public class_distinguisher_base<N>
	{
	public:

		class_distinguisher_advanced(const std::string &groupFile) {
			load_file(groupFile);
		}

		bf_tt<N> get_representative(const bf_tt<N> &f) const override {

			auto gr = this->group_of(f);

			if (gr->functions.size() == 1)
				return gr->functions[0];

			auto index = std::distance(std::begin(this->groups_), gr);

			return (distinguisher_[index].type == distinguisher_t::local_connection) ? 
				find_representative_lc(f, this->groups_[index], this->groups_, distinguisher_[index]) :
				find_representative_ind(f, this->groups_[index], this->groups_, distinguisher_[index]);
		}

	private:

		int load_file(const std::string &fileName) {

			benchmark bm("loading group file " + fileName, benchmark::report::end);

			std::ifstream ifile(fileName);

			if (!ifile.is_open())
				return -1;

			while (!ifile.eof())
			{
				std::string line;

				std::getline(ifile, line);

				if (!line.empty())
				{
					group<N> gr;
					ae_distinguisher aed;

					std::istringstream stream(line);

					uint64_t fnSize, group_id;
					std::string groupstr, sizestr, diststr;

					stream >> groupstr >> group_id >> sizestr >> fnSize;

					LOG_TRACE << "loading group " << group_id << " functions : " << fnSize << LOG_ENDL;

					gr.functions.reserve(fnSize);

					if (fnSize > 1)
					{
						stream >> diststr;

						aed.type = (diststr == "lc") ? distinguisher_t::local_connection : distinguisher_t::indicator;

						uint64_t distance;
						while (!stream.eof())
							if (stream >> distance)
								aed.data.push_back(neighbor_groups{ distance });
					}

					// load functions
					for (uint64_t i = 0; i < fnSize; i++)
					{
						std::getline(ifile, line);

						bf_tt<N> f;
						f.table().set_hex(line);

						gr.functions.push_back(f);

						this->reps_.push_back(f);

						if (i == 0)
							gr.ai = calculate_affine_invariants(f);
					}

					this->groups_.push_back(gr);

					distinguisher_.push_back(aed);
				}
			}

			LOG_TRACE << "groups : " << this->groups_.size() << LOG_ENDL;

			LOG_VERIFY(std::is_sorted(std::begin(this->groups_), std::end(this->groups_)));
			//std::sort(std::begin(this->groups_), std::end(this->groups_));

			std::sort(std::begin(this->reps_), std::end(this->reps_));

			progress_bar pb("processing groups", this->groups_.size());

#pragma omp parallel for schedule(dynamic)
			for (int g = 0; g < this->groups_.size(); g++)
			{
				LOG_TRACE << "processing group #" << g << LOG_ENDL;

				if (this->groups_[g].functions.size() > 1)
				{
					if (distinguisher_[g].type == distinguisher_t::local_connection)
					{
						for (auto &dd : distinguisher_[g].data)
							dd.groupIndices = compute_neighbor_groups_lc(this->groups_, this->groups_[g].functions, dd.distance);
					}
					else
					{
						for (auto &dd : distinguisher_[g].data)
							dd.groupIndices = compute_neighbor_groups_ind(this->groups_, this->groups_[g].functions, dd.distance);
					}
				}

#pragma omp critical
				pb.update();
			}

			return 0;
		}

		std::vector<ae_distinguisher> distinguisher_;
	};

	template<uint64_t N>
	std::unique_ptr<class_distinguisher_base<N>> class_distinguisher();

	template<>
	std::unique_ptr<class_distinguisher_base<6>> class_distinguisher();

	template<uint64_t N>
	std::unique_ptr<class_distinguisher_base<N>> class_distinguisher() {
		static_assert(N >= 2 && N <= 5, "");
		return std::unique_ptr<class_distinguisher_base<N>>(new class_distinguisher_basic<N>(load_representatives<N>()));
	}
}
