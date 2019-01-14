#pragma once
#include "../params.hpp"
#include "grid.hpp"
#include "write_occupied.hpp"
#include "rates_table.hpp"
#include "rate.hpp"

#include <algorithm>
#include <cassert>
#include <random>
#include <utility>

// Monte-Carlo algorithm based on variable step size method
class Monte_carlo
{
private:
	using Index = unsigned int;

public:
	Monte_carlo(Rate rate_fn) : rates_{std::move(rate_fn), grid_}
	{
		std::random_device rnd_device;
		rnd_generator_.seed(rnd_device());
	}

	template<typename Bool, class Boundary_predicate>
	void init(const Tensor<Bool>& available_points, Boundary_predicate is_boundary_point)
	{
		grid_.init(available_points, is_boundary_point);
	}

	// Randomly uniformly distributes points in the lattice with a given filling factor
	void fill_uniform(double filling)
	{
		grid_.fill_uniform(filling, rnd_generator_);
	}

	double estimate_step_duration()
	{
		rates_.init();
		const auto rate = rates_.total_rate();

		assert(rate > 0);
		return 1 / rate;
	}

	// Runs the simulation for the given amount of time, returns the number of steps performed
	std::pair<unsigned int, double> run(unsigned int max_n_steps, double max_duration)
	{
		std::exponential_distribution<double> rnd_exp;

		rates_.init();

		const auto start_time = time_;
		const auto end_time = start_time + max_duration;
		unsigned int n_steps = 0;

		while (time_ < end_time && n_steps < max_n_steps)
		{
			const auto int_rate = rates_.total_int_rate();
			const auto bnd_rate = rates_.total_bnd_rate();

			std::uniform_real_distribution<double> rnd_01;
			if (const auto r = (int_rate + bnd_rate) * rnd_01(rnd_generator_); r <= int_rate)
			{
				const auto next_event = rates_.int_event(r);
				update_int_events(next_event);
			}
			else
			{
				const auto next_event = rates_.bnd_event(r - int_rate);
				update_bnd_events(next_event);
			}

			time_ += rnd_exp(rnd_generator_) / (int_rate + bnd_rate);
			++n_steps;
		}

		return {n_steps, time_ - start_time};
	}

	const Grid<Index>& grid() const
	{
		return grid_;
	}

	void write(const std::string& file_name) const
	{
		write_occupied(file_name, grid_);
	}

private:
	// Updates the grid and the rates for the given index of the chosen internal event
	void update_int_events(std::size_t ev_index)
	{
		const auto [src, dest] = rates_.int_event_points(ev_index);
		assert(grid_.is_occupied(src));
		assert(grid_.is_empty(dest));

		constexpr auto is_src_empty = true; // `src` becomes empty
		rates_.update_rates_for_neighbours_as_dest(src, is_src_empty);
		rates_.update_rates_for_neighbours_as_src(src, is_src_empty);

		grid_.move_occupied(src, dest);

		constexpr auto is_dest_empty = false; // `dest` becomes occupied
		rates_.update_rates_for_neighbours_as_dest(dest, is_dest_empty);
		rates_.update_rates_for_neighbours_as_src(dest, is_dest_empty);
	}

	void update_bnd_events(std::size_t ev_index)
	{
		const auto point = rates_.bnd_event_point(ev_index);

		// If the point was occupied, it now becomes empty, and vice versa
		bool is_point_empty = grid_.is_occupied(point);

		if (is_point_empty)
		{
			rates_.remove_events(point);
			grid_.mark_empty(point);
		}
		else
		{
			rates_.add_events();
			grid_.mark_occupied(point);
		}

		rates_.update_rates_for_neighbours_as_dest(point, is_point_empty);
		rates_.update_rates_for_neighbours_as_src(point, is_point_empty);
	}

private:
	Grid<Index> grid_;
	Rates_table rates_;

	std::mt19937 rnd_generator_;

	double time_ = 0;
};
