#include "simulator.hpp"
#include "params.hpp"

#include <es_util/numeric.hpp>

#include <cassert>
#include <numeric>
#include <vector>

void Simulator::compute_filament_shape(std::vector<unsigned int>& radii) const
{
	assert(radii.size() == grid_size_z_);

	std::vector<unsigned int> occup(grid_center_xy_ + 1);
	std::vector<unsigned int> total(grid_center_xy_ + 1);

	for (auto z = 0u; z < grid_size_z_; ++z)
	{
		radii[z] = 0;

		std::fill(occup.begin(), occup.end(), 0);
		std::fill(total.begin(), total.end(), 0);

		for (auto x = 0u; x < grid_size_xy_; ++x)
			for (auto y = 0u; y < grid_size_xy_; ++y)
			{
				const auto xd = static_cast<double>(x) - grid_center_xy_;
				const auto yd = static_cast<double>(y) - grid_center_xy_;
				const auto rd = es_util::hypot(xd, yd);
				const auto r = static_cast<unsigned int>(std::round(rd));

				if (r <= grid_center_xy_)
				{
					occup[r] += mc_.grid().is_occupied({x, y, z});
					total[r] += 1;
				}
			}

		std::partial_sum(occup.begin(), occup.end(), occup.begin());
		std::partial_sum(total.begin(), total.end(), total.begin());

		assert(params::min_filament_radius > 0);
		for (auto r = grid_center_xy_; r >= params::min_filament_radius; --r)
			if (const auto frac = static_cast<double>(occup[r]) / total[r];
				frac >= params::filament_filling_threshold)
			{
				radii[z] = r;
				break;
			}
	}
}

void Simulator::compute_potential_and_heat(const std::vector<unsigned int>& filament_shape)
{
	assert(filament_shape.size() == grid_size_z_);

	// Compute resistivity
	std::vector<double> resistivity(grid_size_z_, 0);
	for (auto z = 0u; z < grid_size_z_; ++z)
	{
		if (metallic_regions_[z])
			continue;

		if (const auto radius = filament_shape[z]; radius > 0)
		{
			const auto area = es_util::math::pi * es_util::sq(params::grid_spacing * radius);
			resistivity[z] = params::filament_resistivity / area;
		}
		else
			resistivity[z] = params::grain_bnd_resistivity;
	}

	// Compute current
	const auto total_resistivity = std::accumulate(
		resistivity.begin() + 1, resistivity.end() - 1,
		(resistivity.front() + resistivity.back()) / 2);

	resistance_ = params::grid_spacing * total_resistivity;
	const auto current = bias_ / resistance_;

	// Compute heat source linear density
	assert(core_heat_.size() == grid_size_z_);
	assert(core_potential_.size() == grid_size_z_ + 1);

	for (auto z = 0u; z < grid_size_z_; ++z)
		core_heat_[z] = es_util::sq(current) * resistivity[z];

	// Compute potential
	double acc_resistivity = resistivity.front() / 2;

	core_potential_.front() = 0;
	for (auto z = 1u; z < grid_size_z_; ++z)
	{
		core_potential_[z] = current * params::grid_spacing * acc_resistivity;
		acc_resistivity += resistivity[z];
	}
	core_potential_.back() = bias_;
}
