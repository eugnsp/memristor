#pragma once
#include <esu/phys.hpp>

namespace params
{
using namespace esu::au::literals;

// Physical tags of mesh elements; the values should agree
// with those defined in the .geo/.msh mesh file
enum Tags : unsigned int
{
	DIELECTRIC = 1,
	TIP = 2,
	GRANULE = 3,
	CONTACT = 4
};

//////////////////////////////////////////////////////////////////////
/** Monte-Carlo simulation */

// Initial filling factor of O-vacancies grid
inline constexpr auto initial_filling = .06;

// O-vacancies grid spacing
inline constexpr auto grid_spacing = .25_nm;

// O-vacancy Debye frequency
inline constexpr auto debye_frequency = 1e13 / 1_sec;

// O-vacancy activation energy
inline constexpr auto activation_energy = 1.1_evolt;

// Maximum number of steps with fixed temperature and potential
inline constexpr auto steps_per_round = 25'000u;

// Maximum step duration
inline constexpr auto max_step_duration = .25e-3_sec;

//////////////////////////////////////////////////////////////////////
/** Filament shape calculation */

// Minimal radius of the filament (positive; in grid_delta units)
inline constexpr auto min_filament_radius = 4u;

// Filament formation filling threshold
inline constexpr auto filament_filling_threshold = .2;

//////////////////////////////////////////////////////////////////////
/** Resistance calculation */

// Filament volume resistivity (Ohm * m^2 / m = Ohm * m)
inline constexpr auto filament_resistivity = 2.4e-6_ohm * 1_m;

// Grain boundary linear resistivity (Ohm / m)
inline constexpr auto grain_bnd_resistivity = 5e3_ohm / 1_nm;

//////////////////////////////////////////////////////////////////////
/** Heat equation */

// Ambient temperature
inline constexpr auto temperature = 300_kelvin;

// Thermal conductivity
inline constexpr auto thermal_conductivity = 3_watt / (1_m * 1_kelvin);

// Heat source radius (grain boundary + filament)
inline constexpr auto heat_source_radius = 2_nm;

//////////////////////////////////////////////////////////////////////
/** Main simulation */

// Bias voltage sweep speed
inline constexpr auto bias_sweep_rate = 200_volt / 1_sec;

// Maximum external bias voltage
inline constexpr auto max_bias = 3.4_volt;

// Minimum external bias voltage
//inline constexpr auto min_bias = 3.4_volt;

// Maximum current
inline constexpr auto max_current = 25e-6_amp;
} // namespace params
