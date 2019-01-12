#pragma once
#include "point.hpp"

#include <es_util/enum_class_index.hpp>

#include <cassert>

ES_UTIL_ENUM_CLASS_INDEX_TYPE(unsigned char, Neighbour_index)

inline static constexpr auto n_neighbours = Neighbour_index{6};

inline Neighbour_index twin(Neighbour_index neighbour)
{
	assert(neighbour < n_neighbours);

	// Should be compatible with `neighbours[]` definition in the `neighbour` function
	return static_cast<Neighbour_index>(*neighbour ^ 1);
}

// For the given point returns its neighbour specified the the given index
inline Point neighbour(const Point& point, Neighbour_index neighbour)
{
	constexpr auto m1 = static_cast<unsigned int>(-1);
	static constexpr Point neighbours[] = {{1, 0, 0},  {m1, 0, 0}, {0, 1, 0},
										   {0, m1, 0}, {0, 0, 1},  {0, 0, m1}};

	assert(neighbour < n_neighbours);
	return point + neighbours[*neighbour];
}

// Calls the given function for each neihgbour point of the given point;
// the function should have the signature equivalent to the following:
// `void fn(const Point&, Neighbour_index);`
template<class Fn>
void for_each_neighbour(const Point& point, Fn&& fn)
{
	for (Neighbour_index n_index{0}; n_index < n_neighbours; ++n_index)
	{
		const auto n_point = neighbour(point, n_index);
		fn(n_point, n_index);
	}
}
