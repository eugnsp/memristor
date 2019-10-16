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
#include <esu/numeric.hpp>
#include <esu/phys.hpp>
#include <esu/algorithm.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include <iostream>

using Heat_sp_solver = esl::Pardiso_solver<esl::Csr_matrix<double, esl::Symmetric_upper>>;

class Heat_solver final : public esf::Matrix_based_solver<Heat_system, Heat_sp_solver>
{
private:
	using Base = esf::Matrix_based_solver<Heat_system, Heat_sp_solver>;
	using System = typename Base::System;
	using Element = typename System::template Var_t<0>::Element;

	static constexpr auto n_dofs = Element::n_total_face_dofs;

public:
	using Base::mesh;
	using Base::system;

public:
	Heat_solver(const esf::Mesh<2>& mesh, const std::vector<unsigned int>& tags,
		const std::vector<double>& core_heat_source) :
		Base(mesh),
		tags_(tags), core_heat_(core_heat_source)
	{
		system().variable().set_name("temp");
	}

	void init()
	{
		const auto mesh_br = mesh().bounding_rect();
		const auto width = mesh_br.width();
		const auto height = mesh_br.height();

		system().variable().set_bnd_cond<0>(
			mesh(), esf::Linestring{{0., 0.}, {width, 0.}, {width, height}, {0., height}}, params::temperature);

		Base::init();
		esf::compute_and_set_sparsity_pattern(system(), matrix_);
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
					solution_[vertex_dofs[i].index] = bc.value();
				}
			}

			if constexpr (Element::has_edge_dofs)
			{
				for (auto& halfedge : bc.halfedges())
				{
					typename System::template Var_edge_dofs<0> edge_dofs;
					system().dof_mapper().template edge_dofs<0>(edge(halfedge), edge_dofs);
					for (esf::Local_index i = 0; i < edge_dofs.size(); ++i)
					{
						assert(!edge_dofs[i].is_free);
						solution_[edge_dofs[i].index] = bc.value();
					}
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
		using Quadr = esf::Quadr<2 * (Element::order - 1) + 1, 2>;

		const auto tag = tags_[**face];

		const auto grads = esf::gradients<Element, Quadr>(esf::inv_transp_jacobian(face));
		const auto factor =
			area(face) *
			(tag == params::Tags::GRANULE ? 2 * params::thermal_conductivity : params::thermal_conductivity) *
			center(face).x();

		const auto matrix = esf::stiffness_matrix<Element, Quadr>(grads, factor);

		// auto vertex_circ = face.vertex_circ();
		// const auto x1 = vertex_circ->vertex().x();
		// const auto x2 = (++vertex_circ)->vertex().x();
		// const auto x3 = (++vertex_circ)->vertex().x();

		// constexpr auto n_dofs = Element::n_total_face_dofs;
		// la::Matrix_d<n_dofs, n_dofs> matrix;

		// for (esf::Local_index i = 0; i < n_dofs; ++i)
		// 	for (esf::Local_index j = 0; j < n_dofs; ++j)
		// 		matrix(i, j) = factor * Quadr::sum(
		// 			[i, j, &grads, x1, x2, x3](auto q)
		// 			{
		// 				const auto pt = Quadr::point(q);
		// 				const auto x = x1 * (1 - pt[0] - pt[1]) + x2 * pt[0] + x3 * pt[1];
		// 				return x * dot(grads(q, i), grads(q, j));
		// 			});

		add_matrix_to_global(system().dof_mapper().dofs(face), matrix);

		const auto br = bounding_rect(face);
		if (esf::is_geom_greater_equal(br.left(), params::heat_source_radius))
			return;

		auto vertex_circ = face.vertex_circ();
		const auto p1 = vertex_circ->vertex();
		const auto p2 = (++vertex_circ)->vertex();
		const auto p3 = (++vertex_circ)->vertex();

		const auto f = area(face) * center(face).x() /
					   (esu::math::pi * params::heat_source_radius * params::heat_source_radius);

		using Q = esf::Quadr<5, 2>;
		const auto vector = esf::load_vector<Element, Q>(
			[this, p1, p2, p3](auto q) {
				const auto pt = Q::point(q);
				const auto x = p1.x() * (1 - pt[0] - pt[1]) + p2.x() * pt[0] + p3.x() * pt[1];
				if (x > params::heat_source_radius)
					return 0.;

				const auto y = p1.y() * (1 - pt[0] - pt[1]) + p2.y() * pt[0] + p3.y() * pt[1];
				const auto z = y;

				const auto index = static_cast<std::size_t>(std::round(z / params::grid_spacing));
				return index < core_heat_.size() ? core_heat_[index] : 0;
			},
			f);

		add_vector_to_global(system().dof_mapper().dofs(face), vector);
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

	template<class Dofs, class Expr>
	void add_vector_to_global(const Dofs& dofs, const Expr& vector)
	{
		for (esf::Local_index r = 0; r < n_dofs; ++r)
			if (dofs[r].is_free)
				rhs_[dofs[r].index] += vector[r];
	}

	void write(const std::string& file_name) const
	{
		using namespace esu::au::literals;
		//		return;

		esl::Vector_xd temp(*mesh().n_vertices(), 0);

		for (esf::Vertex_index vertex{0}; vertex < mesh().n_vertices(); ++vertex)
		{
			typename Base::System::template Var_vertex_dofs<0> vertex_dofs;
			system().dof_mapper().template vertex_dofs<0>(vertex, vertex_dofs);
			temp[*vertex] = esu::au::to_kelvin(solution_[vertex_dofs[0].index]);
		}

		esf::Matlab_writer2 m(file_name, mesh(), 1_nm);
		m.write_vertex_field("data", temp);
	}

private:
	using Base::matrix_;
	using Base::rhs_;
	using Base::solution_;

	const std::vector<unsigned int>& tags_;
	const std::vector<double>& core_heat_;
};
