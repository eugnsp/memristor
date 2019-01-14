#include "../params.hpp"
#include "tensor.hpp"

// A proxy class to represent a flag indicating that a point
// should be included into a Monte-Carlo grid
class Is_point_available
{
public:
	void operator=(unsigned int tag)
	{
		is_inside_system = true;
		if (tag == params::Tags::TIP || tag == params::Tags::GRANULE)
			is_available = false;
	}

	explicit operator bool() const
	{
		return is_inside_system && is_available;
	}

private:
	bool is_inside_system = false;
	bool is_available = true;
};
