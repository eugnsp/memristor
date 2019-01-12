#include "simulator.hpp"
#include "monte_carlo/tensor.hpp"
#include "monte_carlo/is_point_available.hpp"
#include "interpolate.hpp"
#include "params.hpp"
#include "tags_as_solution.hpp"

#include <es_fe/geom/compare.hpp>
#include <es_fe/mesh/io.hpp>
#include <es_fe/mesh/mesh2.hpp>
#include <es_fe/mesh/tools/mesh_filter.hpp>
#include <es_util/phys.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

void Simulator::init_meshes(const char* mesh_file)
{
	using namespace es_util::au::literals;

	const unsigned int physical_tag_index = 1;
	heat_mesh_ = es_fe::read_gmsh_mesh(mesh_file, 1_nm);
	heat_tags_ = es_fe::read_gmsh_tags(mesh_file, physical_tag_index);

	// The mesh for the Poisson equation is obtained from that for the heat equation
	// by excluding all elements with tag = CONTACT
	const auto tag_filter = [](auto tag) { return tag != params::Tags::CONTACT; };
	const auto face_filter = [this, tag_filter](const auto& face) {
		return tag_filter(heat_tags_[**face]);
	};

	poisson_mesh_ = es_fe::mesh_filter_by_faces(heat_mesh_, face_filter);
	poisson_tags_.resize(*poisson_mesh_.n_faces());
	std::copy_if(heat_tags_.begin(), heat_tags_.end(), poisson_tags_.begin(), tag_filter);
}

void Simulator::init_monte_carlo()
{
	const auto mc_grid_size = Point{grid_size_xy_, grid_size_xy_, grid_size_z_};

	mc_potential_.resize(mc_grid_size);
	mc_potential_.assign(0);

	mc_temp_.resize(mc_grid_size);
	mc_temp_.assign(0);

	Tensor<Is_point_available> available_points(mc_grid_size);
	interpolate(Tags_as_solution{poisson_mesh_, poisson_tags_}, available_points);

	const auto is_boundary_point = [this](const Point& point)
	{
		const auto z = point.z * params::grid_spacing;
		return es_fe::is_geom_equal(std::abs(z - system_height_), 0);
	};

	mc_.init(available_points, is_boundary_point);
	mc_.fill_uniform(params::initial_filling);
}

void Simulator::init()
{
	init_meshes("../mesh/mesh3.msh");

	const auto mesh_br = poisson_mesh_.bounding_rect();
	system_radius_ = mesh_br.width();
	system_height_ = mesh_br.height();

	grid_center_xy_ = static_cast<unsigned int>(system_radius_ / params::grid_spacing);
	grid_size_xy_ = 1 + 2 * grid_center_xy_;
	grid_size_z_ = 1 + static_cast<unsigned int>(system_height_ / params::grid_spacing);

	core_heat_.resize(grid_size_z_);
	core_potential_.resize(grid_size_z_ + 1);

	init_monte_carlo();
	find_forbidden_grid_regions();

	std::cout << "System diameter: " << es_util::au::to_nm(2 * system_radius_) << " nm\n"
			  << "System height:   " << es_util::au::to_nm(system_height_)
			  << " nm\n\n"
			  //   << "Poisson/heat FEM mesh file: " << mesh_file << '\n'
			  << poisson_mesh_ << '\n'
			  << std::endl;
}
