#pragma once
#include <es_fe/mesh/algorithm/vertices_in_linestrip.hpp>
#include <es_fe/mesh/mesh2.hpp>
#include <es_util/iterator.hpp>

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>

class Heat_dirichlet
{
public:
	Heat_dirichlet(const es_fe::Mesh2& mesh, const es_fe::Linestring& ls, double value) :
		value_(value)
	{
		std::tie(vertices_, halfedges_) = es_fe::vertices_and_halfedges_in_linestrip(ls, mesh);

		std::transform(
			halfedges_.begin(), halfedges_.end(), std::back_inserter(edges_), [](auto halfedge) {
				return edge(halfedge);
			});
	}

	static constexpr bool is_essential()
	{
		return true;
	}

	auto begin_vertex() const
	{
		return vertices_.begin();
	}

	auto end_vertex() const
	{
		return vertices_.end();
	}

	auto begin_halfedge() const
	{
		return halfedges_.begin();
	}

	auto end_halfedge() const
	{
		return halfedges_.end();
	}

	auto begin_edge() const
	{
		return edges_.begin();
	}

	auto end_edge() const
	{
		return edges_.end();
	}

	auto halfedges() const
	{
		return es_util::Iterable{begin_halfedge(), end_halfedge()};
	}

	double value() const
	{
		return value_;
	}

private:
	std::vector<es_fe::Vertex_index> vertices_;
	std::vector<es_fe::Edge_index> edges_;
	std::vector<es_fe::Halfedge_index> halfedges_;

	const double value_;
};
