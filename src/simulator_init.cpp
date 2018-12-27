#include "simulator.hpp"

#include "monte_carlo/monte_carlo.hpp"
#include "params.hpp"

#include <es_la/dense.hpp>
#include <es_fe/mesh/io.hpp>
#include <es_fe/mesh/mesh2.hpp>
#include <es_fe/mesh/tools/mesh_filter.hpp>
#include <es_geom/algorithm.hpp>
#include <es_geom/compare.hpp>
#include <es_phys/atomic_units.hpp>
#include <es_util/numeric.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>

void Simulator::init_meshes(const char* mesh_file)
{
	using namespace es_phys::au::literals;

	const unsigned int physical_tag_index = 1;
	heat_mesh_ = es_fe::read_gmsh_mesh(mesh_file, 1_nm);
	heat_tags_ = es_fe::read_gmsh_tags(mesh_file, physical_tag_index);

	// The mesh for the Poisson equation is obtained from that for the heat equation
	// by excluding all elements with tag = CONTACT
	const auto tag_filter = [](auto tag) { return tag != Mesh_tags::CONTACT; };
	const auto face_filter = [this, tag_filter](const auto& face) {
		return tag_filter(heat_tags_[**face]);
	};

	poisson_mesh_ = es_fe::mesh_filter_by_faces(heat_mesh_, face_filter);
	poisson_tags_.resize(*poisson_mesh_.n_faces());
	std::copy_if(heat_tags_.begin(), heat_tags_.end(), poisson_tags_.begin(), tag_filter);
}

void Simulator::init()
{
	using namespace es_phys::au::literals;

	init_meshes("../mesh/mesh3.msh");

	const auto mesh_br = poisson_mesh_.bounding_rect();
	system_radius_ = mesh_br.width();
	system_height_ = mesh_br.height();

	grid_center_xy_ = static_cast<unsigned int>(system_radius_ / params::grid_spacing);
	grid_size_xy_ = 1 + 2 * grid_center_xy_;
	grid_size_z_ = 1 + static_cast<unsigned int>(system_height_ / params::grid_spacing);
	const auto mc_grid_size = Point3{grid_size_xy_, grid_size_xy_, grid_size_z_};

	mc_potential_.resize(mc_grid_size);
	mc_potential_.set_all(0);

	mc_temp_.resize(mc_grid_size);
	mc_temp_.set_all(0);

	find_forbidden_grid_regions();

	mc_.define_shape(
		mc_grid_size,
		// TODO : use mesh tags
		[this](Point3 pt) {
			const auto x = pt.x * params::grid_spacing - params::grid_spacing * grid_center_xy_;
			const auto y = pt.y * params::grid_spacing - params::grid_spacing * grid_center_xy_;
			const auto r = std::hypot(x, y);

			const auto z = pt.z * params::grid_spacing;

			if (r > 15_nm) // system boundary
				return false;
			if (z < 5_nm && r < 1_nm) // tip
				return false;
			if (es_util::sq(r / 1.5_nm) + es_util::sq((z - 20_nm) / 3_nm) < 1) // granule
				return false;

			return true;
		});

	mc_.init_uniform(params::initial_filling);

	std::cout << "System diameter: " << es_phys::au::to_nm(2 * system_radius_) << " nm\n"
			  << "System height:   " << es_phys::au::to_nm(system_height_)
			  << " nm\n\n"
			  //   << "Poisson/heat FEM mesh file: " << mesh_file << '\n'
			  << poisson_mesh_ << '\n'
			  << "Poisson/heat Monte-Carlo grid:\n"
			  << mc_potential_ << std::endl;
}
