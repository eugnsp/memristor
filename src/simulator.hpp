#pragma once
#include "params.hpp"
#include "poisson_solver.hpp"
#include "heat_solver.hpp"
#include "interpolate.hpp"
#include "monte_carlo.hpp"

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

class Simulator
{
public:
	void run(int /* argc */, const char** /* argv */)
	{
		using namespace es_util::au::literals;

		// Step 1
		init();

		Heat_solver heat_solver{heat_mesh_, heat_tags_, core_heat_};
		heat_solver.init();

		Poisson_solver poisson_solver{poisson_mesh_, poisson_tags_, core_potential_};
		poisson_solver.init();

		mc_.init_uniform(params::initial_filling);

		// Step 2
		voltage_bias_ = 0;
		double time = 0;
		unsigned int i = 0;

		const auto rate_fn = [this](const auto& from, const auto& to) {
			return mc_rate_fn(from, to);
		};

		std::vector<unsigned int> fil_shape(grid_size_z_);

		std::vector<double> biases;
		std::vector<double> currents;

		bool sweep_direction = true;

		std::cout << "#    Time      Fil.size  Bias      Current   MC est.      MC act.\n"
				  << "     [msec]    [layers]  [V]       [uA]      [usec/step]  [steps]\n"
				  << "-----------------------------------------------------------------"
				  << std::endl
				  << std::setprecision(3);

		while (true)
		{
			std::cout << std::left << std::setw(4) << i + 1 << ' ' << std::setw(9)
					  << es_util::au::to_sec(time) * 1e3 << ' ' << std::flush;

			// Step 3
			compute_filament_shape(fil_shape);
			const auto fil_size =
				std::count_if(fil_shape.begin(), fil_shape.end(), [](auto r) { return r > 0; });

			// Step 4
			compute_potential_and_heat(fil_shape);
			std::cout << std::setw(9) << fil_size << ' ' << std::setw(9)
					  << es_util::au::to_volt(voltage_bias_) << ' ' << std::setw(9)
					  << es_util::au::to_amp(current_) * 1e6 << ' ' << std::flush;

			// Step 5
			heat_solver.solve();
			interpolate(heat_solver.solution_view(), mc_temp_);

			poisson_solver.solve();
			interpolate(poisson_solver.solution_view(), mc_potential_);

			// Step 6
			double time_step = params::min_step_duration;
			double est_step_duration = mc_.estimate_step_duration(rate_fn);
			std::cout << std::setw(9) << es_util::au::to_sec(est_step_duration) * 1e6 << "    "
					  << std::flush;

			if (est_step_duration < time_step)
			{
				// Step 7
				const auto mc_res = mc_.run(sweep_direction ? params::steps_per_round
															: 10 * params::steps_per_round,
											time_step, rate_fn);
				time_step = mc_res.second;

				std::cout << std::setw(9) << mc_res.first << std::endl;
			}
			else
				std::cout << 0 << std::endl;

			{
				biases.push_back(es_util::au::to_volt(voltage_bias_));
				currents.push_back(es_util::au::to_amp(current_) * 1e6);

				la::Matfile_writer mw("iv.mat");
				mw.write("v", biases);
				mw.write("i", currents);

				if (i % 5 == 0)
				{
					mc_.write("out/mc" + std::to_string(i / 5) + ".mat");
					poisson_solver.write("out/p" + std::to_string(i / 5) + ".mat");
					heat_solver.write("out/h" + std::to_string(i / 5) + ".mat");
				}
			}

			// Step 8
			if (sweep_direction)
				voltage_bias_ += time_step * params::bias_sweep_rate;
			else
				voltage_bias_ -= time_step * params::bias_sweep_rate;

			time += time_step;
			++i;

			if (fil_size > 40)
				sweep_direction = false;

			if (!sweep_direction && voltage_bias_ < 0)
				break;
			//		 mc_.write("mc" + std::to_string(i) + ".mat");
		}
	}

private:
	void init_meshes(const char* mesh_file);
	void init();

	void find_forbidden_grid_regions()
	{
		metallic_regions_.resize(grid_size_z_, false);

		for (const auto& face : poisson_mesh_.faces())
		{
			const auto tag = poisson_tags_[**face];
			if (tag != params::Tags::TIP && tag != params::Tags::GRANULE)
				continue;

			const auto br = bounding_rect(face);

			const auto z_min = static_cast<unsigned int>(
				std::floor(std::max(0., (br.bottom() - es_fe::delta) / params::grid_spacing)));
			const auto z_max =
				std::min(grid_size_z_ - 1, static_cast<unsigned int>(std::ceil(
											   (br.top() + es_fe::delta) / params::grid_spacing)));

			for (auto z = z_min; z <= z_max; ++z)
				metallic_regions_[z] = true;
		}

		// using T = monte_carlo::Point::Type;
		// min_radius_.resize(grid_size_z_, 0);

		// for (const auto& face : mesh_->faces())
		// {
		// 	const auto tag = mesh_tags_[**face];
		// 	if (tag != Mesh_tags::METALLIC_GERM && tag != Mesh_tags::GRANULE)
		// 		continue;

		// 	auto [r_min, z_min, r_max, z_max] =
		// 		bounding_rect_in_grid(es_geom::bounding_rect(face));

		// 	r_max = std::min(r_max, grid_center_xy_);
		// 	z_max = std::min(z_max, grid_size_z_ - 1);

		// 	for (T z = z_min; z <= z_max; ++z)
		// 		for (T r = r_min; r <= r_max; ++r)
		// 		{
		// 			const es_geom::Point pt{params::grid_delta * r, params::grid_delta * z};
		// 			if (es_geom::contains(face, pt))
		// 				min_radius_[z] = std::max(min_radius_[z], r + 1);
		// 		}
		// }
	}

	double mc_rate_fn(const Point3& from, const Point3& to)
	{
		const auto temp = (mc_temp_[from] + mc_temp_[to]) / 2;
		const auto delta_phi = mc_potential_[to] - mc_potential_[from];

		const auto rate =
			params::debye_frequency * std::exp(-(params::activation_energy + delta_phi) / temp);
		return rate;
	};

	void compute_filament_shape(std::vector<unsigned int>&) const;
	void compute_potential_and_heat(const std::vector<unsigned int>& filament_shape);

private:
	es_fe::Mesh2 poisson_mesh_;
	std::vector<unsigned int> poisson_tags_;

	es_fe::Mesh2 heat_mesh_;
	std::vector<unsigned int> heat_tags_;

	double voltage_bias_;
	double current_;

	std::vector<double> core_potential_;
	std::vector<double> core_heat_;

	double system_radius_;
	double system_height_;

	unsigned int grid_center_xy_;
	unsigned int grid_size_xy_;
	unsigned int grid_size_z_;

	Grid3<double> mc_temp_;
	Grid3<double> mc_potential_;

	Monte_carlo mc_;
	std::vector<bool> metallic_regions_;
};
