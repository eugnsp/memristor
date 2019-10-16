#pragma once
#include "point.hpp"

#include <esl/io/matfile_writer.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

// A tensor of rank 3 with first-dimension-major storage
template<class T>
class Tensor
{
	template<typename S>
	friend class Tensor;

private:
	using Container = std::vector<T>;
	using Coord = Point::Type;

public:
	Tensor() = default;

	explicit Tensor(const Point& extents, const T& init_value = T{}) : extents_(extents)
	{
		data_.resize(linear_size(extents_), init_value);
	}

	const Point& extents() const
	{
		return extents_;
	}

	// Returns the total number of elements in the tensor (the product of extents)
	std::size_t size() const
	{
		return linear_size(extents_);
	}

	// Returns the approximate total size of memory in bytes occupied by the tensor
	std::size_t memory_size() const
	{
		return data_.capacity() * sizeof(T);
	}

	//////////////////////////////////////////////////////////////////////
	/** Element access */

	typename Container::reference operator[](const Point& point)
	{
		return data_[linear_index(point)];
	}

	typename Container::const_reference operator[](const Point& point) const
	{
		return data_[linear_index(point)];
	}

	typename Container::reference operator()(Coord x, Coord y, Coord z)
	{
		return (*this)[{x, y, z}];
	}

	typename Container::const_reference operator()(Coord x, Coord y, Coord z) const
	{
		return (*this)[{x, y, z}];
	}

	//////////////////////////////////////////////////////////////////////
	/** Modifiers */

	template<typename S>
	Tensor& operator=(const Tensor<S>& other)
	{
		extents_ = other.extents();
		data_.assign(other.data_.begin(), other.data_.end());

		return *this;
	}

	void resize(const Point& extents, const T& value = T())
	{
		extents_ = extents;
		data_.resize(linear_size(extents_), value);
	}

	void assign(const T& value)
	{
		data_.assign(data_.size(), value);
	}

	void assign(const Point& extents, const T& value)
	{
		extents_ = extents;
		data_.assign(linear_size(extents), value);
	}

	void clear()
	{
		extents_.x = extents_.y = extents_.z = 0;
		data_.clear();
	}

	//////////////////////////////////////////////////////////////////////

	// Returns the number of elements that satisfy the given predicate
	template<class Unary_predicate>
	std::size_t count_if(Unary_predicate pred) const
	{
		return std::count_if(data_.begin(), data_.end(), pred);
	}

	// Applies the given function object to all elements traversing
	// the tensor in the contiguous first-dimension-major order;
	// the signature should be equivalent to the following:
	// `void fn(T&, const Point&);`
	template<class Fn>
	void for_each(Fn fn)
	{
		auto it = data_.begin();

		for (Coord z = 0; z < extents_.z; ++z)
			for (Coord y = 0; y < extents_.y; ++y)
				for (Coord x = 0; x < extents_.x; ++x)
					fn(*it++, Point{x, y, z});

		assert(it == data_.end());
	}

	void write(const std::string& file_name)
	{
		esl::Matfile_writer mat_file(file_name);

		mat_file.write("nx", extents_.x);
		mat_file.write("ny", extents_.y);
		mat_file.write("nz", extents_.z);
		mat_file.write("data", data_);
	}

private:
	static std::size_t linear_size(const Point& extents)
	{
		return static_cast<std::size_t>(extents.x) * extents.y * extents.z;
	}

	std::size_t linear_index(const Point& pt) const
	{
		assert(contains(extents_, pt));

		const auto x_dim = static_cast<std::size_t>(extents_.x);
		const auto y_dim = static_cast<std::size_t>(extents_.y);
		return pt.x + x_dim * (pt.y + y_dim * pt.z);
	}

private:
	Point extents_{};
	Container data_;
};
