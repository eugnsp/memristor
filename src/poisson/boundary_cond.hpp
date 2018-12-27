#pragma once
#include "../params.hpp"

#include <es_fe/types.hpp>
#include <es_fe/mesh/mesh2.hpp>
#include <es_fe/mesh/algorithm/vertices_in_linestrip.hpp>

#include <es_geom/compare.hpp>
#include <es_geom/linestring.hpp>
#include <es_geom/point.hpp>

#include <cassert>
#include <cstddef>
#include <vector>

namespace poisson
{
class Dirichlet
{
public:
	Dirichlet(const es_fe::Mesh2& mesh, const es_geom::Linestring& ls)
	{
		vertices_ = es_fe::vertices_in_linestrip(ls, mesh);
	}

	static constexpr bool is_essential()
	{
		return true;
	}

	auto begin_vertex() const
	{
		return vertices_.begin();
	}

	auto end_vertex() const
	{
		return vertices_.end();
	}

protected:
	std::vector<es_fe::Vertex_index> vertices_;
};

class Dirichlet_const : public Dirichlet
{
public:
	using Dirichlet::Dirichlet;

	double value(const es_geom::Point&) const
	{
		return value_;
	}

	void set_value(double value)
	{
		value_ = value;
	}

private:
	double value_;
};

class Dirichlet_core : public Dirichlet
{
public:
	Dirichlet_core(const es_fe::Mesh2& mesh, const es_geom::Linestring& ls,
				   const std::vector<double>& core_potential) :
		Dirichlet(mesh, ls), core_potential_(core_potential)
	{
		// Remove the first and the last points, they belong to the electrodes
		vertices_.erase(vertices_.begin());
		vertices_.pop_back();
	}

	// Returns the boundary value at a given point using linear interpolation
	// from z-grid mid-points (i + 1/2), at which the core potential is defined
	double value(const es_geom::Point& pt) const
	{
		using params::grid_spacing;

		assert(es_geom::is_equal(pt.x(), 0));

		const auto z = pt.y() - grid_spacing / 2;
		if (z <= 0)
			return core_potential_.front();

		const auto index = static_cast<std::size_t>(z / grid_spacing);
		if (index + 1 >= core_potential_.size())
			return core_potential_.back();

		const auto alpha = z / grid_spacing - index;
		return (1 - alpha) * core_potential_[index] + alpha * core_potential_[index + 1];
	}

private:
	const std::vector<double>& core_potential_;
};
} // namespace poisson
