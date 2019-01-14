#pragma once
#include "params.hpp"
#include "poisson/solver.hpp"
#include "heat/solver.hpp"
#include "interpolate.hpp"
#include "monte_carlo/monte_carlo.hpp"
#include "monte_carlo/rate.hpp"

#include <es_la/dense.hpp>
#include <es_la/io/matfile_writer.hpp>
#include <es_fe/mesh/mesh2.hpp>
#include <es_util/phys.hpp>
#include <es_util/numeric.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>

inline std::string fname(std::string prefix, unsigned int index)
{
	char name[200];
	sprintf(name, "mat/%s%.3d.mat", prefix.c_str(), index);
	return std::string{name};
}

class Simulator
{
public:
	Simulator() : mc_(Rate{mc_temp_, mc_potential_})
	{}

	void run(int /* argc */, const char** /* argv */)
	{
		using namespace es_util::au::literals;

		// Step 1
		init();

		Heat_solver heat_solver{heat_mesh_, heat_tags_, core_heat_};
		heat_solver.init();

		Poisson_solver poisson_solver{poisson_mesh_, poisson_tags_, core_potential_};
		poisson_solver.init();

		// Step 2
		double external_bias = 0;
		double limiting_resistance = 0;
		double time = 0;
		int sweep = 0;

		std::vector<unsigned int> filament_shape(grid_size_z_);
		std::vector<double> biases, currents, time_steps;

		std::cout << std::setprecision(4);

		bias_ = external_bias;
		for (unsigned int i = 0; ; ++i)
		{
			if (i % 20 == 0)
				std::cout
					<< "-----------------------------------------------------------------------------------------------------------\n"
					<< "#    Time      Fil.vol   Ext.bias  Syst.bias Syst.res. Lim.res.  Current   MC est.      MC act.   Num.vacs.\n"
					<< "     [msec]    [arb.u.]  [V]       [V]       [kOhm]    [kOhm]    [uA]      [usec/step]  [steps]            \n"
					<< "-----------------------------------------------------------------------------------------------------------"
					<< std::endl;

			std::cout << std::left << std::setw(4) << i + 1 << ' ' << std::setw(9)
					  << es_util::au::to_sec(time) / 1e-3 << ' ' << std::flush;

			// Step 3
			compute_filament_shape(filament_shape);
			unsigned int filament_volume = 0;
			for (auto r : filament_shape)
				filament_volume += r * r;

			// Step 4
			compute_potential_and_heat(filament_shape);
			std::cout << std::setw(9) << filament_volume << ' ' << std::setw(9)
					  << es_util::au::to_volt(external_bias) << ' ' << std::setw(9)
					  << es_util::au::to_volt(bias_) << ' ' << std::setw(9)
					  << es_util::au::to_ohm(bias_ / current_) / 1e3 << ' ' << std::setw(9)
					  << es_util::au::to_ohm(limiting_resistance) / 1e3 << ' ' << std::setw(9)
					  << es_util::au::to_amp(current_) / 1e-6 << ' ' << std::flush;

			// Step 5
			heat_solver.solve();
			interpolate(heat_solver.solution_view(), mc_temp_);

			poisson_solver.solve();
			interpolate(poisson_solver.solution_view(), mc_potential_);

			// Step 6
			double time_step = params::min_step_duration;
			double est_step_duration = mc_.estimate_step_duration();
			std::cout << std::setw(9) << es_util::au::to_sec(est_step_duration) * 1e6 << "    "
					  << std::flush;

			if (est_step_duration < time_step)
			{
				// Step 7
				const auto mc_res = mc_.run(params::steps_per_round, time_step);
				time_step = mc_res.second;

				std::cout << std::setw(9) << mc_res.first << std::flush;
			}
			else
				std::cout << std::setw(9) << 0 << std::flush;

			std::cout << ' ' << std::setw(9) << mc_.grid().n_occupied() << std::endl;

			{
				biases.push_back(es_util::au::to_volt(bias_));
				currents.push_back(es_util::au::to_amp(current_) / 1e-6);
				time_steps.push_back(time_step);

				la::Matfile_writer mw("mat/iv.mat");
				mw.write("v", biases);
				mw.write("i", currents);
				mw.write("t", time_steps);

				mc_.write(fname("m", i));
				poisson_solver.write(fname("p", i));
				heat_solver.write(fname("h", i));
			}

			// Step 8
			external_bias += (sweep == 1 ? -1 : 1) * time_step * params::bias_sweep_rate;

			auto limiting_drop = std::abs(external_bias) / params::max_current * current_ - bias_;
			if (external_bias >= 0)
				limiting_drop = std::max(0., limiting_drop);
			else
				limiting_drop = std::min(0., limiting_drop);

			bias_ = external_bias - limiting_drop;
			limiting_resistance = limiting_drop / current_;

			time += time_step;

			if (external_bias >= params::max_bias)
				sweep = 1;
			if (sweep == 1 && external_bias < -params::max_bias)
				sweep = 2;
			if (sweep == 2 && external_bias > 0)
				break;
		}
	}

private:
	void init_meshes(const char* mesh_file);
	void init_metallic_regions();
	void init_monte_carlo();

	void init();

	void compute_filament_shape(std::vector<unsigned int>&) const;
	void compute_potential_and_heat(const std::vector<unsigned int>& filament_shape);

private:
	es_fe::Mesh2 poisson_mesh_;
	std::vector<unsigned int> poisson_tags_;

	es_fe::Mesh2 heat_mesh_;
	std::vector<unsigned int> heat_tags_;

	double bias_;
	double current_;

	std::vector<double> core_potential_;
	std::vector<double> core_heat_;

	double system_radius_;
	double system_height_;

	unsigned int grid_center_xy_;
	unsigned int grid_size_xy_;
	unsigned int grid_size_z_;

	Tensor<double> mc_temp_;
	Tensor<double> mc_potential_;

	Monte_carlo mc_;
	std::vector<bool> metallic_regions_;
};
