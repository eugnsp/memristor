#pragma once
#include "../params.hpp"
#include "system.hpp"

#include <esf/dof/tools.hpp>
#include <esf/geometry.hpp>
#include <esf/io/matlab_writer2.hpp>
#include <esf/math.hpp>
#include <esf/matrix_based/solver.hpp>
#include <esf/mesh/mesh2.hpp>
#include <esf/var_list.hpp>
#include <esl/dense.hpp>
#include <esl/io/matfile_writer.hpp>
#include <esl/sparse/solver/pardiso_solver.hpp>
#include <esu/algorithm.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include <iostream>

using Poisson_sp_solver = esl::Pardiso_solver<esl::Csr_matrix<double, esl::Symmetric_upper>>;

class Poisson_solver final : public esf::Matrix_based_solver<Poisson_system, Poisson_sp_solver>
{
private:
	using Base = esf::Matrix_based_solver<Poisson_system, Poisson_sp_solver>;
	using System = typename Base::System;
	using Element = typename System::template Var_t<0>::Element;

public:
	using Base::mesh;
	using Base::system;

public:
	Poisson_solver(
		const esf::Mesh<2>& mesh, const std::vector<unsigned int>& tags, const std::vector<double>& core_potential) :
		Base(mesh),
		tags_(tags), core_potential_(core_potential)
	{
		system().variable().set_name("phi");
	}

	void init()
	{
		const auto mesh_br = mesh().bounding_rect();
		const auto width = mesh_br.width();
		const auto height = mesh_br.height();

		system().variable().set_bnd_cond<0>(mesh(), esf::Linestring{{0, 0}, {width, 0}});
		system().variable().set_bnd_cond<1>(mesh(), esf::Linestring{{0, height}, {width, height}});
		system().variable().set_bnd_cond<2>(mesh(), esf::Linestring{{0, 0}, {0, height}}, core_potential_);

		Base::init();
		esf::compute_and_set_sparsity_pattern(system(), matrix_);
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
		system().variable().for_each_ess_bnd_cond([this](const auto& bc) {
			for (auto& vertex : bc.vertices())
			{
				typename System::template Var_vertex_dofs<0> vertex_dofs;
				system().dof_mapper().template vertex_dofs<0>(vertex, vertex_dofs);

				for (esf::Local_index i = 0; i < vertex_dofs.size(); ++i)
				{
					assert(!vertex_dofs[i].is_free);
					solution_[vertex_dofs[i].index] = bc.value(mesh().vertex(vertex));
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

	void assemble_on_face(const esf::Mesh<2>::Face_view& face)
	{
		const auto tag = tags_[**face];

		const bool is_metal = tag == params::Tags::GRANULE || tag == params::Tags::TIP;

		const auto metal_eps = 100.;
		const auto eps = is_metal ? metal_eps : 1.;
		const auto factor = -area(face) * center(face).x() * eps;

		using Quadr = esf::Quadr<1, 2>;
		const auto grads = esf::gradients<Element, Quadr>(esf::inv_transp_jacobian(face));
		const auto stiffness_matrix = esf::stiffness_matrix<Element, Quadr>(grads, factor);

		add_matrix_to_global(system().dof_mapper().dofs(face), stiffness_matrix);
	}

	template<class Dofs, class Expr>
	void add_matrix_to_global(const Dofs& dofs, const Expr& matrix)
	{
		for (esf::Local_index c = 0; c < dofs.size(); ++c)
			if (dofs[c].is_free)
			{
				for (esf::Local_index r = 0; r <= c; ++r)
					if (dofs[r].is_free)
					{
						auto [i1, i2] = esu::sorted(dofs[r].index, dofs[c].index);
						matrix_(i1, i2) += matrix(r, c);
					}
			}
			else
			{
				for (esf::Local_index r = 0; r < dofs.size(); ++r)
					if (dofs[r].is_free)
					{
						auto i1 = dofs[r].index;
						auto i2 = dofs[c].index;

						const auto m = matrix(r, c);
						rhs_[i1] -= m * solution_[i2];
					}
			}
	}

	void write(const std::string& file_name)
	{
		using namespace esu::au::literals;
		return;

		esl::Vector_xd phi(*mesh().n_vertices(), 0);

		for (esf::Vertex_index vertex{0}; vertex < mesh().n_vertices(); ++vertex)
		{
			typename Base::System::template Var_vertex_dofs<0> vertex_dofs;
			system().dof_mapper().template vertex_dofs<0>(vertex, vertex_dofs);
			phi[*vertex] = esu::au::to_volt(solution_[vertex_dofs[0].index]);
		}

		esf::Matlab_writer2 m(file_name, mesh(), 1_nm);
		m.write_vertex_field("data", phi);
	}

private:
	using Base::matrix_;
	using Base::rhs_;
	using Base::solution_;

	const std::vector<unsigned int>& tags_;
	const std::vector<double>& core_potential_;
};
