#pragma once
#include "../params.hpp"
#include "matrix3.hpp"

#include <es_geom/algorithm.hpp>
#include <es_geom/compare.hpp>
#include <es_geom/point.hpp>
#include <es_util/numeric.hpp>

#include <algorithm>
#include <cmath>

template<class Fe_solution_view>
void interpolate(const Fe_solution_view& from, Matrix3<double>& to)
{
	using params::grid_spacing;

	assert(to.extents().x == to.extents().y);
	const auto size_xy = to.extents().x;
	const auto size_z = to.extents().z;
	const auto center_xy = (size_xy - 1) / 2;

	for (auto& face : from.mesh().faces())
	{
		const auto br = es_geom::bounding_rect(face);

		const auto left = (br.left() - es_geom::delta) / grid_spacing;
		const auto right = (br.right() + es_geom::delta) / grid_spacing;
		const auto bottom = (br.bottom() - es_geom::delta) / grid_spacing;
		const auto top = (br.top() + es_geom::delta) / grid_spacing;

		const auto z_min = static_cast<unsigned int>(std::max(0., std::floor(bottom)));
		const auto z_max = std::min(size_z - 1, static_cast<unsigned int>(std::ceil(top)));
		const auto y_max = static_cast<unsigned int>(right);

		for (auto z = z_min; z <= z_max; ++z)
			for (auto y = center_xy - y_max; y <= center_xy + y_max; ++y)
			{
				const auto yd = static_cast<double>(y) - center_xy;
				const auto x_min =
					(std::abs(yd) >= left)
						? 0u
						: static_cast<unsigned int>(std::ceil(es_util::cathetus(left, yd)));
				const auto x_max = static_cast<unsigned int>(es_util::cathetus(right, yd));

				for (auto x = center_xy - x_max; x <= center_xy - x_min; ++x)
				{
					const auto xd = static_cast<double>(x) - center_xy;
					const auto r = es_util::hypot(xd, yd);

					const es_geom::Point pt{r * grid_spacing, z * grid_spacing};
					if (es_geom::contains(face, pt))
						to[{x, y, z}] = from(face, pt);
				}

				for (auto x = center_xy + x_min; x <= center_xy + x_max; ++x)
				{
					const auto xd = static_cast<double>(x) - center_xy;
					const auto r = es_util::hypot(xd, yd);

					const es_geom::Point pt{r * grid_spacing, z * grid_spacing};
					if (es_geom::contains(face, pt))
						to[{x, y, z}] = from(face, pt);
				}
			}
	}
}
