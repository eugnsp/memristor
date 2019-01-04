#pragma once

// 3D point with integral coordinates
struct Point3
{
	unsigned int x;
	unsigned int y;
	unsigned int z;
};

inline Point3 operator+(Point3 pt1, const Point3& pt2)
{
	pt1.x += pt2.x;
	pt1.y += pt2.y;
	pt1.z += pt2.z;

	return pt1;
}

inline bool operator==(const Point3& pt1, const Point3& pt2)
{
	return pt1.x == pt2.x && pt1.y == pt2.y && pt1.z == pt2.z;
}

inline bool operator!=(const Point3& pt1, const Point3& pt2)
{
	return !(pt1 == pt2);
}
