#pragma once
#include "monte_carlo/monte_carlo.hpp"

#include <vector>


std::vector<unsigned int> filament_shape(const Monte_carlo&);

std::vector<double> core_resistivity(const std::vector<unsigned int>& filament_shape);
