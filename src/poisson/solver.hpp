#pragma once
#include "system.hpp"
#include "../params.hpp"

#include <es_math/const.hpp>
#include <es_la/solver/pardiso_solver.hpp>
#include <es_la/io/matfile_writer.hpp>
#include <es_geom/algorithm.hpp>
#include <es_fe/math/jacobian.hpp>
#include <es_fe/mesh/mesh2.hpp>
#include <es_fe/matrix_based/solver.hpp>

#include <es_fe/matrix_based/seq_assembler.hpp>
#include <es_fe/mesh/tools/mesh_face_colour_map.hpp>

#include <es_la/function.hpp>
#include <es_fe/var_list.hpp>
#include <es_fe/dof/layered_dof_mapper.hpp>
#include <es_fe/dof/tools.hpp>
#include <es_fe/element/lagrange.hpp>
#include <es_fe/math/matrix.hpp>
#include <es_fe/math/quadr.hpp>
#include <es_fe/io/matlab_writer.hpp>
#include <es_fe/io/vtk_writer_3d.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include <iostream>

namespace poisson
{
using Sp_solver = la::Pardiso_solver<la::Sparse_matrix<double, la::Symmetric_upper>>;

class Solver final : public es_fe::Matrix_based_solver<System, Sp_solver>
{
private:
	using Base = es_fe::Matrix_based_solver<System, Sp_solver>;
	using System = typename Base::System;

public:
	using Base::mesh;
	using Base::system;

public:
	Solver(const es_fe::Mesh2& mesh, const std::vector<unsigned int>& tags,
		const std::vector<double>& core_potential) :
		Base(mesh), tags_(tags), core_potential_(core_potential)
	{
		system().variable().set_name("phi");
	}

	void init()
	{
		const auto mesh_br = mesh().bounding_rect();
		const auto width = mesh_br.width();
		const auto height = mesh_br.height();

		system().variable().set_bnd_cond<0>(mesh(), es_geom::Linestring{{0, 0}, {width, 0}});
		system().variable().set_bnd_cond<1>(mesh(), es_geom::Linestring{{0, height}, {width, height}});
		system().variable().set_bnd_cond<2>(mesh(), es_geom::Linestring{{0, 0}, {0, height}}, core_potential_);

		Base::init();
		es_fe::compute_and_set_sparsity_pattern(system(), matrix_);
	}

	void solve()
	{
		system().variable().bnd_cond<0>().set_value(core_potential_.front());
		system().variable().bnd_cond<1>().set_value(core_potential_.back());

		Base::solve();
	}

private:
	virtual void set_bnd_values() override
	{
		system().variable().for_each_ess_bnd_cond(
			[this](const auto& bc)
			{
				for (auto vertex = bc.begin_vertex(); vertex != bc.end_vertex(); ++vertex)
				{
					typename System::template Var_vertex_dofs<0> vertex_dofs;
					system().dof_mapper().template vertex_dofs<0>(*vertex, vertex_dofs);
					//vertex_dof_indices.resize(vertex_dofs.rows(), vertex_dofs.cols());

					for (std::size_t j = 0; j < vertex_dofs.cols(); ++j)
						for (std::size_t i = 0; i < vertex_dofs.rows(); ++i)
						{
							assert(vertex_dofs(i, j).is_free == false);
							solution_[vertex_dofs(i, j).index] = bc.value(mesh().vertex(*vertex));
						}
				}
			});
	}

public:
	virtual void assemble() override
	{
   		for (const auto& face : mesh().faces())
   			assemble_on_face(face);
	}

	void assemble_on_face(const es_fe::Mesh2::Face_view& face)
	{
		const auto tag = tags_[**face];

		const auto eps = (tag == Mesh_tags::GRANULE || tag == Mesh_tags::METALLIC_GERM) ? 100. : 1.;
		const auto factor = -es_geom::area(face) * es_geom::center(face).x() * eps;

		using Quadr = es_fe::Quadr<1, 2>;
		const auto grads = es_fe::gradients<Element, Quadr>(es_fe::inv_transp_jacobian(face));
		const auto stiffness_matrix = es_fe::stiffness_matrix<Element, Quadr>(grads, factor);

		add_matrix_to_global(system().dofs(face), stiffness_matrix);
	}

	template<class Dofs, class Expr>
	void add_matrix_to_global(const Dofs& dofs, const Expr& matrix)
	{
		for (es_fe::Local_index c = 0; c < dofs.size(); ++c)
			if (dofs[c].is_free)
			{
				for (es_fe::Local_index r = 0; r <= c; ++r)
					if (dofs[r].is_free)
					{
						auto i1 = dofs[r].index;
						auto i2 = dofs[c].index;
						es_util::sort2(i1, i2);

						matrix_(i1, i2) += matrix(r, c);
					}
			}
			else
			{
				for (es_fe::Local_index r = 0; r < dofs.size(); ++r)
					if (dofs[r].is_free)
					{
						auto i1 = dofs[r].index;
						auto i2 = dofs[c].index;

						const auto m = matrix(r, c);
						rhs_[i1] -= m * solution_[i2];
					}
			}
	}

	virtual void after_solve() override
	{
		using namespace es_phys::au::literals;
		return;

		la::Vector_xd phi(*mesh().n_vertices(), 0);

		for (es_fe::Vertex_index vertex{0}; vertex < mesh().n_vertices(); ++vertex)
		{
			typename Base::System::template Var_vertex_dofs<0> vertex_dofs;
			system().dof_mapper().template vertex_dofs<0>(vertex, vertex_dofs);
			phi[*vertex] = es_phys::au::to_volt(solution_[vertex_dofs[0].index]);
		}

		es_fe::Matlab_writer m("p.mat", mesh(), 1_nm);
		m.write_vertex_field("data", phi);
	}

	void write(const std::string& file_name)
	{
		using namespace es_phys::au::literals;
//		return;

		la::Vector_xd phi(*mesh().n_vertices(), 0);

		for (es_fe::Vertex_index vertex{0}; vertex < mesh().n_vertices(); ++vertex)
		{
			typename Base::System::template Var_vertex_dofs<0> vertex_dofs;
			system().dof_mapper().template vertex_dofs<0>(vertex, vertex_dofs);
			phi[*vertex] = es_phys::au::to_volt(solution_[vertex_dofs[0].index]);
		}

		es_fe::Matlab_writer m(file_name, mesh(), 1_nm);
		m.write_vertex_field("data", phi);
	}


private:
	using Base::solution_;
	using Base::rhs_;
	using Base::matrix_;

	const std::vector<unsigned int>& tags_;
	const std::vector<double>& core_potential_;
};
}
