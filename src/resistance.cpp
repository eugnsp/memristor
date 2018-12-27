#include "resistance.hpp"
#include "params.hpp"
#include "monte_carlo/monte_carlo.hpp"

#include <es_util/numeric.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

std::vector<unsigned int> filament_shape(const Monte_carlo& mc)
{
	assert(mc.grid_extents().x == mc.grid_extents().y);

	const auto size_xy = mc.grid_extents().x;
	const auto size_z = mc.grid_extents().z;
	const auto center_xy = (size_xy - 1) / 2;

	std::vector<unsigned int> occup(center_xy + 1);
	std::vector<unsigned int> total(center_xy + 1);

	std::vector<unsigned int> radii(size_z, 0);
	for (auto z = 0u; z < size_z; ++z)
	{
		std::fill(occup.begin(), occup.end(), 0);
		std::fill(total.begin(), total.end(), 0);

		for (auto x = 0u; x < size_xy; ++x)
			for (auto y = 0u; y < size_xy; ++y)
			{
				const auto xd = static_cast<double>(x) - center_xy;
				const auto yd = static_cast<double>(y) - center_xy;
				const auto rd = es_util::hypot(xd, yd);
				const auto r = static_cast<unsigned int>(std::round(rd));

				if (r <= center_xy)
				{
					occup[r] += mc.is_occupied({x, y, z});
					total[r] += 1;
				}
			}

		std::partial_sum(occup.begin(), occup.end(), occup.begin());
		std::partial_sum(total.begin(), total.end(), total.begin());

		assert(params::min_filament_radius > 0);
		for (auto r = center_xy; r >= params::min_filament_radius; --r)
			if (const auto frac = static_cast<double>(occup[r]) / total[r];
				frac >= params::filament_filling_threshold)
			{
				radii[z] = r;
				break;
			}
	}

	return radii;
}
