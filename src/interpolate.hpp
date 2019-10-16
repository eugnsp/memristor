#pragma once
#include "params.hpp"
#include "monte_carlo/tensor.hpp"

#include <esf/geometry.hpp>
#include <esu/numeric.hpp>

#include <algorithm>
#include <cmath>

template<class Fe_solution_view, typename T>
void interpolate(
	const Fe_solution_view& fe_solution, Tensor<T>& tensor)
{
	using params::grid_spacing;

	assert(tensor.extents().x == tensor.extents().y);
	const auto size_xy = tensor.extents().x;
	const auto size_z = tensor.extents().z;
	const auto center_xy = (size_xy - 1) / 2;

	for (auto& face : fe_solution.mesh().faces())
	{
		const auto br = bounding_rect(face);

		const auto left = (br.left() - esf::delta) / grid_spacing;
		const auto right = (br.right() + esf::delta) / grid_spacing;
		const auto bottom = (br.bottom() - esf::delta) / grid_spacing;
		const auto top = (br.top() + esf::delta) / grid_spacing;

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
						: static_cast<unsigned int>(std::ceil(esu::cathetus(left, yd)));
				const auto x_max = static_cast<unsigned int>(esu::cathetus(right, yd));

				for (auto x = center_xy - x_max; x <= center_xy - x_min; ++x)
				{
					const auto xd = static_cast<double>(x) - center_xy;
					const auto r = esu::hypot(xd, yd);

					const esf::Point2 pt{r * grid_spacing, z * grid_spacing};
					if (contains(face, pt))
						tensor(x, y, z) = fe_solution(pt, face);
				}

				for (auto x = center_xy + x_min; x <= center_xy + x_max; ++x)
				{
					const auto xd = static_cast<double>(x) - center_xy;
					const auto r = esu::hypot(xd, yd);

					const esf::Point2 pt{r * grid_spacing, z * grid_spacing};
					if (contains(face, pt))
						tensor(x, y, z) = fe_solution(pt, face);
				}
			}
	}
}
