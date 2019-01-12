#pragma once
#include "point.hpp"
#include "tensor.hpp"

#include <es_la/io/matfile_writer.hpp>
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

	template<class Shape_predicate, class Boundary_predicate>
	void define_shape(
		Point extents, Shape_predicate is_point_available, Boundary_predicate is_boundary_point)
	{
		sites_.assign(extents, Site{index_invalid, index_invalid});
		boundary_.clear();

		sites_.for_each([&, this](Site& site, const Point& pt)
		{
			if (is_point_available(pt))
			{
				site.occup_index = index_empty;
				if (is_boundary_point(pt))
				{
					site.bnd_index = n_boundary();
					boundary_.push_back(pt);
				}
			}
		});
	}

	// Randomly uniformly distributes points in the grid with a given filling factor
	template<class Uniform_rnd_generator>
	void init(double filling, Uniform_rnd_generator& rnd_generator)
	{
		assert(0 < filling && filling < 1);

		for (const auto& point : occupied_)
			sites_[point].occup_index = index_empty;

		const auto n_available = sites_.count_if(
			[this](const Site& site) { return is_empty(site); });
		const auto n_to_fill = static_cast<Index>(filling * n_available);
		assert(n_to_fill > 0);

		std::vector<Point> available;
		available.reserve(n_available);

		sites_.for_each([this, &available](const Site& site, const Point& point)
		{
			if (is_empty(site))
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
		assert(is_occupied(site));

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
		assert(is_boundary(site));

		return site.bnd_index;
	}

	// Returns the point by its index in the list of occupied points
	const Point& boundary_point(Index index) const
	{
		return boundary_[index];
	}

	// Checks whether the given point belongs to the grid (then it is either empty or occupied)
	bool is_valid(const Point& point) const
	{
		return is_valid(sites_[point]);
	}

	// Checks whether the given point is empty
	bool is_empty(const Point& point) const
	{
		return is_empty(sites_[point]);
	}

	// Checks whether the given point is occupied
	bool is_occupied(const Point& point) const
	{
		return is_occupied(sites_[point]);
	}

	bool is_boundary(const Point& point) const
	{
		return is_boundary(sites_[point]);
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

	void write_occupied(const std::string& file_name) const
	{
		la::Matfile_writer mat_file(file_name);

		const auto n = n_occupied();
		std::vector<unsigned int> x(n), y(n), z(n);

		for (Index index = 0; index < n; ++index)
		{
			const auto& point = occupied_[index];
			x[index] = point.x;
			y[index] = point.y;
			z[index] = point.z;
		}

		mat_file.write("x", x);
		mat_file.write("y", y);
		mat_file.write("z", z);
	}

	Site _site(Point pt) const
	{
		return sites_[pt];
	}

private:
	// Dummy index value of `occup_index` to denote points outside the system boundary
	// and of `bnd_index` to denote non-boundary points
	static constexpr auto index_invalid = static_cast<Index>(-1);

	// Dummy index value of `occup_index` to denote empty points
	static constexpr auto index_empty = index_invalid - 1;

	static constexpr auto max_occup_index = index_empty - 1;

private:
	static bool is_valid(Site site)
	{
		return site.occup_index != index_invalid;
	}

	static bool is_empty(Site site)
	{
		return site.occup_index == index_empty;
	}

	static bool is_occupied(Site site)
	{
		return site.occup_index <= max_occup_index;
	}

	bool is_boundary(Site site) const
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
