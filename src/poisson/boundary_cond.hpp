#pragma once
#include "../params.hpp"
#include "element.hpp"

#include <es_fe/types.hpp>
#include <es_fe/mesh/mesh2.hpp>
#include <es_fe/boundary_cond.hpp>

#include <es_fe/geom/compare.hpp>
#include <es_fe/geom/linestring.hpp>
#include <es_fe/geom/point.hpp>

#include <cassert>
#include <cstddef>
#include <vector>

using Poisson_dirichlet_contact = es_fe::Uniform_boundary_cond<Poisson_element>;

class Poisson_dirichlet_core final : public es_fe::Boundary_cond<Poisson_element>
{
public:
	Poisson_dirichlet_core(
		const es_fe::Mesh2& mesh,
		const es_fe::Linestring& boundary,
		const std::vector<double>& potential) :
		Boundary_cond(mesh, boundary),
		potential_(potential)
	{
		// Remove the first and the last points, they belong to the electrodes
		this->vertices_.erase(this->vertices_.begin());
		this->vertices_.pop_back();
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
