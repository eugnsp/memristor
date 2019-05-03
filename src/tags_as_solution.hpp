#include <es_fe/mesh/mesh2.hpp>
#include <es_fe/geometry.hpp>

#include <cassert>
#include <vector>

class Tags_as_solution
{
public:
	using Value = unsigned int;

public:
	Tags_as_solution(const es_fe::Mesh<2>& mesh, const std::vector<Value>& tags) :
		mesh_(mesh), tags_(tags)
	{
		assert(tags_.size() == *mesh.n_faces());
	}

	Value operator()(const es_fe::Mesh<2>::Face_view& face, const es_fe::Point2&) const
	{
		return tags_[**face];
	}

	const es_fe::Mesh<2>& mesh() const
	{
		return mesh_;
	}

private:
	const es_fe::Mesh<2>& mesh_;
	const std::vector<Value>& tags_;
};
