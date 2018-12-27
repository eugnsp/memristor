#include "simulator.hpp"

#include "monte_carlo/point.hpp"
#include "params.hpp"

#include <es_math/const.hpp>
#include <es_util/numeric.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

std::vector<double> Simulator::compute_resistivity(
	const std::vector<unsigned int>& filament_shape) const
{
	std::vector<double> resist(grid_size_z_, 0);

	for (auto z = 0u; z < grid_size_z_; ++z)
	{
		if (is_metallic_region_[z])
			continue;

		if (const auto radius = filament_shape[z]; radius > 0)
		{
			const auto area = math::pi * es_util::sq(params::grid_spacing * radius);
			resist[z] = params::filament_resistivity / area;
		}
		else
			resist[z] = params::grain_bnd_resistivity;
	}

	return resist;
}

void Simulator::compute_potential_and_heating(const std::vector<unsigned int>& filament_shape)
{
	using params::grid_spacing;

	const auto resist = compute_resistivity(filament_shape);

	const auto total_resistance =
		grid_spacing * std::accumulate(resist.begin(), resist.end(), 0.);
	current_ = voltage_bias_ / total_resistance;

	core_potential_.resize(grid_size_z_ + 1);
	core_potential_.front() = 0;
	core_potential_.back() = voltage_bias_;

	double resistance = grid_spacing * resist.front() / 2;
	for (auto z = 1u; z < grid_size_z_; ++z)
	{
		resistance += grid_spacing * resist[z];
		core_potential_[z] = current_ * resistance;
	}

	core_heat_.resize(grid_size_z_);
	for (auto z = 0u; z < grid_size_z_; ++z)
		core_heat_[z] = es_util::sq(current_) * resist[z];
}
