#pragma once
#include <ostream>

// 3D point with integral coordinates
struct Point
{
	using Type = unsigned int;

	Type x;
	Type y;
	Type z;
};

inline Point operator+(Point pt1, const Point& pt2)
{
	pt1.x += pt2.x;
	pt1.y += pt2.y;
	pt1.z += pt2.z;

	return pt1;
}

inline bool operator==(const Point& pt1, const Point& pt2)
{
	return pt1.x == pt2.x && pt1.y == pt2.y && pt1.z == pt2.z;
}

inline bool operator!=(const Point& pt1, const Point& pt2)
{
	return !(pt1 == pt2);
}

// Checks whether a box with the given extents contains the given point
inline bool contains(const Point& extents, const Point& pt)
{
	return pt.x < extents.x && pt.y < extents.y && pt.z < extents.z;
}

inline std::ostream& operator<<(std::ostream& os, const Point& pt)
{
	os << '(' << pt.x << ", " << pt.y << ", " << pt.z << ')' << std::flush;
	return os;
}
