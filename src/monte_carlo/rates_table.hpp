#pragma once
#include "../params.hpp"
#include "grid.hpp"
#include "neighbour.hpp"

#include <es_util/container/fenwick_tree.hpp>

#include <vector>

class Rates_table
{
private:
	using Index = unsigned int;

public:
	Rates_table(const Grid<Index>& grid) : grid_(grid)
	{}

	///////////////////////////////////////////////////////////////////////
	//* Access */

	double total_rate() const
	{
		return total_int_rate() + total_bnd_rate();
	}

	double total_int_rate() const
	{
		return int_rates_.sum();
	}

	double total_bnd_rate() const
	{
		return !bnd_rates_.is_empty() ? bnd_rates_.sum() : 0;
	}

	// Returns the index of the internal event specified by its given accumulated
	// rate, i.e. returns the smallest event index such that the sum of rates of
	// all events up to that one (inclusive) is not smaller than the given rate
	Index int_event(double accumulated_rate)
	{
		return int_rates_.lower_bound(accumulated_rate);
	}

	// Returns the index of the boundary event specified by its given accumulated
	// rate, i.e. returns the smallest event index such that the sum of rates of
	// all events up to that one (inclusive) is not smaller than the given rate
	Index bnd_event(double accumulated_rate)
	{
		return bnd_rates_.lower_bound(accumulated_rate);
	}

	//////////////////////////////////////////////////////////////////////
	//* Modifiers */

	template<class Rate_fn>
	void init(Rate_fn rate_fn)
	{
		init_internal(rate_fn);
		init_boundary(rate_fn);
	}

	// For the given point iterates over its neighbours and recalculates
	// 1) internal rates of hops to this point from its neighbours,
	// 2) boundary rates of hops to/from its boundary neighbours (if any)
	template<class Rate_fn>
	void update_rates_for_neighbours_as_src(const Point& dest, bool is_dest_empty, Rate_fn rate_fn)
	{
		for_each_neighbour(dest, [&, this](const Point& src, Neighbour_index n_index)
		{
			if (!contains(grid_.extents(), src))
				return;

			const auto is_src_occupied = grid_.is_occupied(src);
			if (is_src_occupied)
			{
				const auto ev_index = int_event_index(grid_.occupied_index(src), twin(n_index));
				int_rates_.set(ev_index, is_dest_empty ? rate_fn(src, dest) : 0);
			}

			if (grid_.is_boundary(src))
			{
				const auto ev_index = bnd_event_index(grid_.boundary_index(src));
				if (is_src_occupied)
				{
					constexpr auto eta = params::initial_filling;
					const auto delta_in_rate = (1 - eta) / eta * rate_fn(dest, src);
					bnd_rates_.add(ev_index, is_dest_empty ? -delta_in_rate : delta_in_rate);

					// if (bnd_rates_[ev_index] < 0)
					// {
					// 	std::cout << is_dest_empty << '\n';
					// 	std::cout << bnd_rates_[ev_index] << ' ' << delta_in_rate << '\n';
					// }
				}
				else
				{
					constexpr auto eta = params::initial_filling;
					const auto delta_out_rate = eta / (1 - eta) * rate_fn(dest, src);
					bnd_rates_.add(ev_index, is_dest_empty ? delta_out_rate : -delta_out_rate);

					// if (bnd_rates_[ev_index] < 0)
					// {
					// 	std::cout << is_dest_empty << '\n';
					// 	std::cout << bnd_rates_[ev_index] << ' ' << delta_out_rate << '\n';
					// }
				}
				// assert(bnd_rates_[ev_index] >= 0);	// TODO : remove
			}
		});
	}

	// For the given point iterates over its neighbours and recalculates
	// 1) internal rates of hops from this point to its neighbours,
	// 2) boundary rates of hops to/from this point (if it is a boundary one)
	template<class Rate_fn>
	void update_rates_for_neighbours_as_dest(const Point& src, bool is_src_empty, Rate_fn rate_fn)
	{
		if (!is_src_empty)
		{
			const auto index = grid_.occupied_index(src);
			for_each_neighbour(src, [&, this](const Point& dest, Neighbour_index n_index)
			{
				const auto ev_index = int_event_index(index, n_index);
				const auto is_dest_empty = contains(grid_.extents(), dest) && grid_.is_empty(dest);
				int_rates_.set(ev_index, is_dest_empty ? rate_fn(src, dest) : 0);
			});
		}

		if (grid_.is_boundary(src))
		{
			const auto ev_index = bnd_event_index(grid_.boundary_index(src));
			bnd_rates_.set(ev_index, boundary_rate(src, is_src_empty, rate_fn));
		}
	}

	std::pair<Point, Point> int_event_points(Index ev_index) const
	{
		assert(ev_index < n_int_events());

		const auto block_index = ev_index / *n_neighbours;
		const auto n_index = static_cast<Neighbour_index>(ev_index % *n_neighbours);

		const auto src = grid_.occupied_point(block_index);
		const auto dest = neighbour(src, n_index);

		assert(grid_.is_occupied(src));
		assert(grid_.is_empty(dest));

		return {src, dest};
	}

	Point bnd_event_point(Index ev_index) const
	{
		assert(ev_index < n_bnd_events());

		return grid_.boundary_point(ev_index);
	}

	void add_events()
	{
		for (Neighbour_index n_index{0}; n_index < n_neighbours; ++n_index)
			int_rates_.push(0);
	}

	void remove_events(const Point& point)
	{
		assert(grid_.is_occupied(point));

		if (grid_.n_occupied() > 1)
		{
			const auto point_index = grid_.occupied_index(point);
			const auto last_index = grid_.n_occupied() - 1;
			for (Neighbour_index n_index{0}; n_index < n_neighbours; ++n_index)
			{
				const auto src_ev_index = int_event_index(last_index, n_index);
				const auto dest_ev_index = int_event_index(point_index, n_index);
				int_rates_.set(dest_ev_index, int_rates_[src_ev_index]);
			};
		}

		for (Neighbour_index n_index{0}; n_index < n_neighbours; ++n_index)
			int_rates_.pop();
	}

private:
	static Index int_event_index(Index occup_index, Neighbour_index n_index)
	{
		return *n_neighbours * occup_index + *n_index;
	}

	// Returns the total number of internal events (some of them can have zero rate)
	Index n_int_events() const
	{
		return int_event_index(grid_.n_occupied(), Neighbour_index{0});
	}

	static Index bnd_event_index(Index bnd_index)
	{
		return bnd_index;
	}

	// Returns the total number of boundary events
	Index n_bnd_events() const
	{
		return bnd_event_index(grid_.n_boundary());
	}

	template<class Rate_fn>
	double boundary_rate(const Point& src, bool is_src_empty, Rate_fn&& rate_fn)
	{
		if (is_src_empty)
		{
			double out_rate = 0;
			for_each_neighbour(src, [&, this](const Point& dest, Neighbour_index)
			{
				if (contains(grid_.extents(), dest) && grid_.is_empty(dest))
					out_rate += rate_fn(src, dest);
			});

			return params::initial_filling / (1 - params::initial_filling) * out_rate;
		}
		else
		{
			double in_rate = 0;
			for_each_neighbour(src, [&, this](const Point& dest, Neighbour_index)
			{
				if (contains(grid_.extents(), dest) && grid_.is_occupied(dest))
					in_rate += rate_fn(dest, src);
			});

			return (1 - params::initial_filling) / params::initial_filling * in_rate;
		}
	}

	template<class Rate_fn>
	void init_internal(Rate_fn rate_fn)
	{
		std::vector<double> rates(n_int_events(), 0);

		for (auto& src : grid_.occupied_points())
		{
			const auto occup_index = grid_.occupied_index(src);
			for_each_neighbour(src, [&, this](const Point& dest, Neighbour_index n_index)
			{
				if (contains(grid_.extents(), dest) && grid_.is_empty(dest))
				{
					const auto ev_index = int_event_index(occup_index, n_index);
					rates[ev_index] = rate_fn(src, dest);
				}
			});
		}

		int_rates_.reset(std::move(rates));
	}

	template<class Rate_fn>
	void init_boundary(Rate_fn&& rate_fn)
	{
		std::vector<double> rates(n_bnd_events(), 0);

		for (auto& src : grid_.boundary_points())
		{
			const auto ev_index = bnd_event_index(grid_.boundary_index(src));
			rates[ev_index] = boundary_rate(src, grid_.is_empty(src), rate_fn);
		}

		bnd_rates_.reset(std::move(rates));
	}

private:
	const Grid<Index>& grid_;

	// Rates of hopping events inside the system
	es_util::Fenwick_tree<double> int_rates_;

	// Rates of in- and out-hopping events due to permeable boundaries
	es_util::Fenwick_tree<double> bnd_rates_;
};
