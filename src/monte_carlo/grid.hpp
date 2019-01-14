#pragma once
#include "point.hpp"
#include "tensor.hpp"

#include <es_util/iterator.hpp>

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>

template<typename Index>
class Grid
{
	static_assert(std::is_unsigned_v<Index>);

public:
	struct Site
	{
		Index occup_index;
		Index bnd_index;
	};

public:
	//////////////////////////////////////////////////////////////////////
	//* Capacity */

	const Point& extents() const
	{
		return sites_.extents();
	}

	// Returns the number of occupied points
	Index n_occupied() const
	{
		return static_cast<Index>(occupied_.size());
	}

	Index n_boundary() const
	{
		return static_cast<Index>(boundary_.size());
	}

	//////////////////////////////////////////////////////////////////////
	//* Initialization */

	template<typename Bool, class Boundary_predicate>
	void init(const Tensor<Bool>& available_points, Boundary_predicate is_boundary_point)
	{
		sites_.assign(available_points.extents(), Site{index_invalid, index_invalid});
		boundary_.clear();

		sites_.for_each([&, this](Site& site, const Point& point) {
			if (available_points[point])
			{
				site.occup_index = index_empty;
				if (is_boundary_point(point))
				{
					site.bnd_index = n_boundary();
					boundary_.push_back(point);
				}
			}
		});
	}

	// Randomly uniformly distributes points in the grid with a given filling factor
	template<class Uniform_rnd_generator>
	void fill_uniform(double filling, Uniform_rnd_generator& rnd_generator)
	{
		assert(0 < filling && filling < 1);

		for (const auto& point : occupied_)
			sites_[point].occup_index = index_empty;

		const auto n_available = sites_.count_if(is_empty_site);
		const auto n_to_fill = static_cast<Index>(filling * n_available);
		assert(n_to_fill > 0);

		std::vector<Point> available;
		available.reserve(n_available);

		sites_.for_each([this, &available](const Site& site, const Point& point) {
			if (is_empty_site(site))
				available.push_back(point);
		});

		occupied_.resize(n_to_fill);
		std::sample(
			available.begin(), available.end(), occupied_.begin(), n_to_fill, rnd_generator);

		Index index = 0;
		for (const auto& point : occupied_)
			sites_[point].occup_index = index++;
	}

	//////////////////////////////////////////////////////////////////////
	//* Element access */

	// Returns the index of the given occupied point in the list of occupied points
	Index occupied_index(const Point& point) const
	{
		const auto site = sites_[point];
		assert(is_occupied_site(site));

		return site.occup_index;
	}

	// Returns the point by its index in the list of occupied points
	const Point& occupied_point(Index index) const
	{
		return occupied_[index];
	}

	// Returns the index of the given occupied point in the list of boundary points
	Index boundary_index(const Point& point) const
	{
		const auto site = sites_[point];
		assert(is_boundary_site(site));

		return site.bnd_index;
	}

	// Returns the point by its index in the list of occupied points
	const Point& boundary_point(Index index) const
	{
		return boundary_[index];
	}

	//////////////////////////////////////////////////////////////////////
	//* Flags */

	// Checks whether the given point belongs to the grid (then it is either empty or occupied)
	bool is_valid(const Point& point) const
	{
		return is_valid_site(sites_[point]);
	}

	// Checks whether the given point is empty
	bool is_empty(const Point& point) const
	{
		return is_empty_site(sites_[point]);
	}

	// Checks whether the given point is occupied
	bool is_occupied(const Point& point) const
	{
		return is_occupied_site(sites_[point]);
	}

	bool is_boundary(const Point& point) const
	{
		return is_boundary_site(sites_[point]);
	}

	//////////////////////////////////////////////////////////////////////
	//* Iterators */

	auto occupied_points() const
	{
		return es_util::Iterable{occupied_};
	}

	auto boundary_points() const
	{
		return es_util::Iterable{boundary_};
	}

	//////////////////////////////////////////////////////////////////////
	//* Modifiers */

	// Moves a particle from the first occupied given point to the second
	// empty given point, the newly occupied point has the same index
	// in the list of occupied points as the old one
	void move_occupied(const Point& src, const Point& dest)
	{
		assert(is_occupied(src));
		assert(is_empty(dest));

		std::swap(sites_[src].occup_index, sites_[dest].occup_index);
		occupied_[sites_[dest].occup_index] = dest;
	}

	// Marks the given empty point as occupied
	void mark_occupied(const Point& point)
	{
		assert(is_empty(point));

		sites_[point].occup_index = n_occupied();
		occupied_.push_back(point);
	}

	// Marks the given occupied point as empty by moving the last point
	// in the list of occupied points into the position of the given point
	void mark_empty(const Point& point)
	{
		assert(is_occupied(point));

		auto& index = sites_[point].occup_index;
		const auto last_point = occupied_.back();

		occupied_[index] = last_point;
		occupied_.pop_back();

		sites_[last_point].occup_index = index_empty;
		std::swap(index, sites_[last_point].occup_index);
	}

private:
	// Dummy index value of `occup_index` to denote points outside the system boundary
	// and of `bnd_index` to denote non-boundary points
	static constexpr auto index_invalid = static_cast<Index>(-1);

	// Dummy index value of `occup_index` to denote empty points
	static constexpr auto index_empty = index_invalid - 1;

	static constexpr auto max_occup_index = index_empty - 1;

private:
	static bool is_valid_site(Site site)
	{
		return site.occup_index != index_invalid;
	}

	static bool is_empty_site(Site site)
	{
		return site.occup_index == index_empty;
	}

	static bool is_occupied_site(Site site)
	{
		return site.occup_index <= max_occup_index;
	}

	bool is_boundary_site(Site site) const
	{
		return site.bnd_index != index_invalid;
	}

private:
	Tensor<Site> sites_;

	// Occupied points in the grid
	std::vector<Point> occupied_;

	// Boundary points in the grid
	std::vector<Point> boundary_;
};
