#pragma once
#include "heat/solver.hpp"
#include "interpolate.hpp"
#include "monte_carlo/monte_carlo.hpp"
#include "monte_carlo/rate.hpp"
#include "params.hpp"
#include "poisson/solver.hpp"
#include "printer.hpp"

#include <esl/dense.hpp>
#include <esl/io/matfile_writer.hpp>
#include <esf/mesh/mesh2.hpp>
#include <esu/phys.hpp>
#include <esu/numeric.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>

class Simulator
{
public:
	Simulator() : mc_(Rate{mc_temp_, mc_potential_})
	{}

	void run(int /* argc */, const char** /* argv */)
	{
		using namespace esu::au::literals;

		// Step 1
		init();

		Heat_solver heat_solver{heat_mesh_, heat_tags_, core_heat_};
		heat_solver.init();

		Poisson_solver poisson_solver{poisson_mesh_, poisson_tags_, core_potential_};
		poisson_solver.init();

		// Step 2
		double external_bias = 0;
		double lim_resistance = 0;
		double time = 0;
		unsigned int sweep = 0;

		std::vector<unsigned int> filament_shape(grid_size_z_);
		std::vector<double> biases, currents, time_steps;

		Printer printer{std::cout};

		bias_ = external_bias;
		for (unsigned int i = 0; ; ++i)
		{
			if (i % 20 == 0)
				printer.header();

			printer(i + 1, 1, 5);
			printer(time, 1e-3_sec);

			// Step 3
			compute_filament_shape(filament_shape);
			unsigned int filament_volume = 0;
			for (auto r : filament_shape)
				filament_volume += r * r;

			// Step 4
			compute_potential_and_heat(filament_shape);
			const auto current = bias_ / resistance_;

			printer(filament_volume);
			printer(external_bias, 1_volt);
			printer(bias_, 1_volt);
			printer(current, 1e-6_amp);
			printer(resistance_, 1e3_ohm);
			printer(lim_resistance, 1e3_ohm);

			// Step 5
			heat_solver.solve();
			interpolate(heat_solver.solution_view(), mc_temp_);

			poisson_solver.solve();
			interpolate(poisson_solver.solution_view(), mc_potential_);

			// Step 6
			double time_step = params::max_step_duration;
			double est_step_duration = mc_.estimate_step_duration();
			printer(est_step_duration, 1e-6_sec, 13);

			if (est_step_duration < time_step)
			{
				// Step 7
				const auto mc_res = mc_.run(params::steps_per_round, time_step);
				time_step = mc_res.second;
				printer(mc_res.first);
			}
			else
				printer(0);

			printer(mc_.grid().n_occupied());
			printer.endl();

			{
				biases.push_back(esu::au::to_volt(bias_));
				currents.push_back(esu::au::to_amp(current) / 1e-6);
				time_steps.push_back(time_step);

				esl::Matfile_writer mw("mat/iv.mat");
				mw.write("v", biases);
				mw.write("i", currents);
				mw.write("t", time_steps);

				auto name = std::to_string(i);
				if (name.length() < 3)
					name.insert(0, 3 - name.length(), '0');

				mc_.write("mat/m" + name + ".mat");
				poisson_solver.write("mat/p" + name + ".mat");
				heat_solver.write("mat/h" + name + ".mat");
			}

			// Step 8
			external_bias += (sweep == 1 ? -1 : 1) * time_step * params::bias_sweep_rate;

			lim_resistance = limiting_resistance(external_bias);
			bias_ = external_bias - lim_resistance * current;

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

	double limiting_resistance(double external_bias) const
	{
		const auto resistance = std::abs(external_bias) / params::max_current - resistance_;
		return std::max(0., resistance);
	}

private:
	esf::Mesh<2> poisson_mesh_;
	std::vector<unsigned int> poisson_tags_;

	esf::Mesh<2> heat_mesh_;
	std::vector<unsigned int> heat_tags_;

	double bias_;
	double resistance_;

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
