#pragma once
#include "grid.hpp"

#include <esl/io/matfile_writer.hpp>

#include <vector>

template<typename Index>
inline void write_occupied(const std::string& file_name, const Grid<Index>& grid)
{
	esl::Matfile_writer mat_file(file_name);

	const auto n = grid.n_occupied();
	std::vector<Point::Type> x(n), y(n), z(n);

	for (Index index = 0; index < n; ++index)
	{
		const auto& point = grid.occupied_point(index);
		x[index] = point.x;
		y[index] = point.y;
		z[index] = point.z;
	}

	mat_file.write("nx", grid.extents().x);
	mat_file.write("ny", grid.extents().y);
	mat_file.write("nz", grid.extents().z);

	mat_file.write("x", x);
	mat_file.write("y", y);
	mat_file.write("z", z);
}
