#pragma once
#include "point3.hpp"

#include <es_la/io/matfile_writer.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

// 3D grid with first-dimension-major storage
template<class T>
class Grid3
{
private:
	using Container = std::vector<T>;

public:
	Grid3() = default;

	explicit Grid3(const Point3& extents, const T& init_value = T()) : extents_(extents)
	{
		assert(extents.x > 0 && extents.y > 0 && extents.z > 0);
		data_.resize(size(), init_value);
	}

	const Point3& extents() const
	{
		return extents_;
	}

	void resize(const Point3& extents, const T& value = T())
	{
		assert(extents.x > 0 && extents.y > 0 && extents.z > 0);

		extents_ = extents;
		data_.resize(size(extents_), value);
	}

	void set_all(const T& value)
	{
		data_.assign(data_.size(), value);
	}

	void clear()
	{
		data_.clear();
	}

	typename Container::reference operator[](const Point3& pt)
	{
		return data_[linear_index(pt)];
	}

	typename Container::const_reference operator[](const Point3& pt) const
	{
		return data_[linear_index(pt)];
	}

	const Container& data() const
	{
		return data_;
	}

	std::size_t size() const
	{
		return size(extents_);
	}

	// Checks whether a point belongs to the grid
	bool contains(const Point3& pt) const
	{
		return pt.x < extents_.x && pt.y < extents_.y && pt.z < extents_.z;
	}

	// Returns the number of points that satisfy the given predicate
	template<class Unary_predicate>
	std::size_t count_if(Unary_predicate pred) const
	{
		return std::count_if(data_.begin(), data_.end(), pred);
	}

	// Applies a given function object to all points of the grid traversing
	// the grid in the contiguous x-major order.
	//
	// Parameters:
	//  fn - the function object to be applied, the signature should be
	//       equivalent to the following: void fn(T&, Point);
	template<class Fn>
	void for_each(Fn fn)
	{
		auto it = data_.begin();

		for (auto z = 0u; z < extents_.z; ++z)
			for (auto y = 0u; y < extents_.y; ++y)
				for (auto x = 0u; x < extents_.x; ++x)
					fn(*it++, Point3{x, y, z});

		assert(it == data_.end());
	}

	// Returns approximate total size of memory in bytes occupied by the data structure
	std::size_t memory_size() const
	{
		return data_.capacity() * sizeof(T);
	}

private:
	std::size_t linear_index(const Point3& pt) const
	{
		assert(contains(pt));
		return pt.x + (pt.y + pt.z * static_cast<std::size_t>(extents_.y)) * extents_.x;
	}

	static std::size_t size(const Point3& extents)
	{
		return static_cast<std::size_t>(extents.x) * extents.y * extents.z;
	}

private:
	Point3 extents_{};
	Container data_;
};

template<class T>
void write(const std::string& file_name, const Grid3<T>& matrix)
{
	la::Matfile_writer mat_file(file_name);

	mat_file.write("nx", matrix.extents().x);
	mat_file.write("ny", matrix.extents().y);
	mat_file.write("nz", matrix.extents().z);
	mat_file.write("data", matrix.data());
}
