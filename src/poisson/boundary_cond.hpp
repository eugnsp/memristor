#pragma once
#include "../params.hpp"

#include <es_fe/types.hpp>
#include <es_fe/mesh/mesh2.hpp>
#include <es_fe/mesh/algorithm/vertices_in_linestrip.hpp>

#include <es_fe/geom/compare.hpp>
#include <es_fe/geom/linestring.hpp>
#include <es_fe/geom/point.hpp>

#include <cassert>
#include <cstddef>
#include <vector>

class Poisson_dirichlet
{
public:
	Poisson_dirichlet(const es_fe::Mesh2& mesh, const es_fe::Linestring& ls)
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

class Poisson_dirichlet_const : public Poisson_dirichlet
{
public:
	using Poisson_dirichlet::Poisson_dirichlet;

	double value(const es_fe::Point&) const
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

class Poisson_dirichlet_core : public Poisson_dirichlet
{
public:
	Poisson_dirichlet_core(
		const es_fe::Mesh2& mesh,
		const es_fe::Linestring& ls,
		const std::vector<double>& potential) :
		Poisson_dirichlet(mesh, ls),
		potential_(potential)
	{
		// Remove the first and the last points, they belong to the electrodes
		vertices_.erase(vertices_.begin());
		vertices_.pop_back();
	}

	// Returns the boundary value at a given point using linear interpolation
	// from z-grid mid-points (i + 1/2), at which the core potential is defined
	double value(const es_fe::Point& pt) const
	{
		assert(es_fe::is_geom_equal(pt.x(), 0));

		auto z = pt.y() / params::grid_spacing + .5;
		const auto n = potential_.size() - 2;

		const auto index = static_cast<std::size_t>(z);
		double alpha = z - index;

		// First and last points are special: they are z-grid points, not mid-points
		if (z < 1)
			alpha = 2 * alpha - 1;
		else if (z > n)
			alpha = 2 * alpha;

		assert(-.00001 < alpha && alpha < 1.00001);
		assert(index <= n);

		return (1 - alpha) * potential_[index] + alpha * potential_[index + 1];
	}

private:
	const std::vector<double>& potential_;
};
