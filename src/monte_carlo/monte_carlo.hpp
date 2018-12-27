#pragma once
#include "point.hpp"
#include "matrix3.hpp"

#include <es_la/dense.hpp>
#include <es_la/io/matfile_writer.hpp>
#include <es_util/container/fenwick_tree.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <random>
#include <vector>
#include <utility>

#include <iostream>

// Monte-Carlo algorithm based on variable step size method
class Monte_carlo
{
private:
	using Index = unsigned int;
	static constexpr auto index_unavailable = static_cast<Index>(-1);
	static constexpr auto index_empty = static_cast<Index>(-2);

	static constexpr Index n_neighbours = 6;
	static constexpr Point3 neighbours[n_neighbours] =
		{{1, 0, 0}, {0 - 1u, 0, 0}, {0, 1, 0}, {0, 0 - 1u, 0}, {0, 0, 1}, {0, 0, 0 - 1u}};

	// Returns the twin neighbour index twin(i): it a site A has B as its
	// i-th neighbour, then the site B has A as its twin(i)-th neighbour
	Index twin_neighbour_index(Index index)
	{
		return index ^ 1;
	}

public:
	template<class Shape_predicate>
	void define_shape(Point3 size, Shape_predicate is_point_available)
	{
		indices_.clear();
		indices_.resize(size, index_unavailable);

		indices_.for_each(
			[&is_point_available](Index& index, const Point3& pt)
			{
				if (is_point_available(pt))
					index = index_empty;
			});
	}

	// Randomly uniformly distributes points in the lattice with a given filling factor
	void init_uniform(double filling)
	{
		assert(filling <= 1);

		time_ = 0;
		for (const auto& point : occupied_)
			indices_[point] = index_empty;

		const auto n_available = indices_.count_if(is_empty);
		const auto n_to_fill = static_cast<Index>(filling * n_available);
		assert(n_to_fill > 0);

		std::vector<Point3> available;
		available.reserve(n_available);

		indices_.for_each(
			[&available](Index index, const Point3& point)
			{
				if (is_empty(index))
					available.push_back(point);
			});

		occupied_.resize(n_to_fill);
		std::sample(available.begin(), available.end(), occupied_.begin(), n_to_fill, rnd_generator_);

		Index index = 0;
		for (const auto& point : occupied_)
			indices_[point] = index++;
	}

	template<class Rate_fn>
	double estimate_step_duration(Rate_fn rate_fn)
	{
		init_rates(rate_fn);
		const auto total_rate = rates_.sum();
		assert(total_rate > 0);

		return 1 / total_rate;
	}

	// Runs the simulation for the given amount of time, returns the number of steps performed
	template<class Rate_fn>
	std::pair<unsigned int, double> run(unsigned int max_n_steps, double max_duration, Rate_fn rate_fn)
	{
		std::exponential_distribution rnd_time_step{};
		std::uniform_real_distribution rnd_num{};

		init_rates(rate_fn);

		const auto start_time = time_;
		const auto end_time = start_time + max_duration;

		unsigned int n_steps = 0;

		while (time_ < end_time && n_steps < max_n_steps)
		{
			const auto total_rate = rates_.sum();
			assert(total_rate > 0);

			time_ += rnd_time_step_(rnd_generator_) / total_rate;

			const auto r = rnd_num_(rnd_generator_);
			const auto event_index = rates_.lower_bound(r * total_rate);
			assert(event_index < rates_.size());

			update_rates(event_index, rate_fn);
			++n_steps;
		}

		return {n_steps, time_ - start_time};
	}

	// Runs the simulation for the given number of steps, returns the elapsed time
	template<class Rate_fn>
	double run_steps(unsigned int n_steps, Rate_fn rate_fn)
	{
		init_rates(rate_fn);

		const auto start_time = time_;

		for (unsigned int i = 0; i < n_steps; ++i)
			do_step(rate_fn);

		return time_ - start_time;
	}

	Point3 grid_extents() const
	{
		return indices_.extents();
	}

	bool is_occupied(const Point3& pt) const
	{
		return is_occupied(indices_[pt]);
	}

	// template<class Solution>
	void write(const std::string& file_name/* , const Solution& sol */) const
	{
		const auto n_occupied = occupied_.size();
		std::vector<unsigned int> vec(n_occupied);

		la::Matfile_writer mat_file(file_name);

		for (std::size_t i = 0; i < n_occupied; ++i)
			vec[i] = occupied_[i].x;
		mat_file.write("x", vec);

		for (std::size_t i = 0; i < n_occupied; ++i)
			vec[i] = occupied_[i].y;
		mat_file.write("y", vec);

		for (std::size_t i = 0; i < n_occupied; ++i)
			vec[i] = occupied_[i].z;
		mat_file.write("z", vec);

		// std::vector<double> c(n_occupied);
		// for (std::size_t i = 0; i < n_occupied; ++i)
		// 	c[i] = sol({occupied_[i]});
		// mat_file.write("c", c);
	}

private:
	static bool is_empty(Index index)
	{
		return index == index_empty;
	}

	static bool is_occupied(Index index)
	{
		return index != index_empty && index != index_unavailable;
	}

	template<class Rate_fn>
	void init_rates(Rate_fn rate_fn)
	{
		const auto n_occupied = occupied_.size();
		assert(n_occupied > 0);

		rates_.reset(n_neighbours * n_occupied);

		for (std::size_t i = 0; i < n_occupied; ++i)
			for (Index k = 0; k < n_neighbours; ++k)
			{
				const auto pt_from = occupied_[i];
				const auto pt_to = pt_from + neighbours[k];
				const bool can_jump = indices_.contains(pt_to) && is_empty(indices_[pt_to]);
				if (can_jump)
				{
					const auto rate = rate_fn(pt_from, pt_to);
					rates_.add(n_neighbours * i + k, rate);
				}
			}
	}

	// Updates the rates list given an index of the event that has been chosen
	template<class Rate_fn>
	void update_rates(Index event_index, Rate_fn rate_fn)
	{
		const auto index = event_index / n_neighbours;
		const auto neighbour = event_index % n_neighbours;
		const auto pt_from = occupied_[index];
		const auto pt_to = pt_from + neighbours[neighbour];

		assert(is_occupied(indices_[pt_from]));
		assert(is_empty(indices_[pt_to]));

		// State of the point (pt_from) changes from occupied to unoccupied

		// Jumps from (pt_from) are now impossible:
		//   1) subtract their rates from (total_rate_diff)
		//   2) remove them from (rates) table (this is done automatically)
		for (Index k = 0; k < n_neighbours; ++k)
			rates_.set(n_neighbours * index + k, 0);

		// Jumps to (pt_from) are now possible:
		//   1) add their rates to (total_rate_diff)
		//   2) add them to (rates) table
		for (Index k = 0; k < n_neighbours; ++k)
		{
			const auto pt_n = pt_from + neighbours[k];

			// Jump from (pt_to) is not accounted for here, (indices[pt_to])
			// has not been updated yet and holds old value (index_empty), see below
			if (indices_.contains(pt_n) && is_occupied(indices_[pt_n]))
			{
				const auto rate = rate_fn(pt_n, pt_from);
				rates_.set(n_neighbours * indices_[pt_n] + twin_neighbour_index(k), rate);
			}
		}

		std::swap(indices_[pt_from], indices_[pt_to]);
		occupied_[index] = pt_to;

		// State of the point (pt_from) changes from unoccupied to occupied

		// Jumps from (pt_to) are now possible:
		//   1) add their rates to (total_rate_diff)
		//   2) add them to (rates) table
		for (Index k = 0; k < n_neighbours; ++k)
		{
			const auto pt_n = pt_to + neighbours[k];

			// Jump to (pt_from) is accounted for here, (indices[pt_from])
			// has been updated yet and holds new value (index_empty)
			if (indices_.contains(pt_n) && is_empty(indices_[pt_n]))
			{
				const auto rate = rate_fn(pt_to, pt_n);
				rates_.set(n_neighbours * index + k, rate);
			}
		}

		// Jumps to (pt_to) are now impossible:
		//   1) subtract their rates from (total_rate_diff)
		//   2) remove them from (rates) table
		for (Index k = 0; k < n_neighbours; ++k)
		{
			const auto pt_n = pt_to + neighbours[k];
			if (indices_.contains(pt_n) && is_occupied(indices_[pt_n]))
				rates_.set(n_neighbours * indices_[pt_n] + twin_neighbour_index(k), 0);
		}
	}

private:
	Matrix3<Index> indices_;
	std::vector<Point3> occupied_;
	es_util::Fenwick_tree<double> rates_;

	std::random_device rnd_device_{};
    std::mt19937 rnd_generator_{rnd_device_()};

	std::exponential_distribution<double> rnd_time_step_{};
	std::uniform_real_distribution<double> rnd_num_{};

	double time_ = 0;
};
