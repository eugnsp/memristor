#pragma once
#include "../params.hpp"
#include "point.hpp"
#include "tensor.hpp"

#include <cmath>

class Rate
{
public:
	Rate(const Tensor<double>& temp, const Tensor<double>& potential) :
		temp_(temp), potential_(potential)
	{}

	double operator()(const Point& src, const Point& dest) const
	{
		const auto temp = (temp_[src] + temp_[dest]) / 2;
		const auto delta = potential_[dest] - potential_[src];

		const auto z = -(params::activation_energy + delta) / temp;
		return params::debye_frequency * std::exp(z);
	}

private:
	const Tensor<double>& temp_;
	const Tensor<double>& potential_;
};
