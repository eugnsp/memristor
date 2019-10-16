#include <esf/mesh/mesh2.hpp>
#include <esf/geometry.hpp>

#include <cassert>
#include <vector>

class Tags_as_solution
{
public:
	using Value = unsigned int;

public:
	Tags_as_solution(const esf::Mesh<2>& mesh, const std::vector<Value>& tags) :
		mesh_(mesh), tags_(tags)
	{
		assert(tags_.size() == *mesh.n_faces());
	}

	Value operator()(const esf::Point2&, const esf::Mesh<2>::Face_view& face) const
	{
		return tags_[**face];
	}

	const esf::Mesh<2>& mesh() const
	{
		return mesh_;
	}

private:
	const esf::Mesh<2>& mesh_;
	const std::vector<Value>& tags_;
};
