#pragma once
#include "params.hpp"
#include "poisson/solver.hpp"
#include "heat/solver.hpp"
#include "interpolate.hpp"
#include "monte_carlo/monte_carlo.hpp"

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
	sprintf(name, "out/%s%.3d.mat", prefix.c_str(), index);
	return std::string{name};
}

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

		// TODO : in init()?
		mc_.init(params::initial_filling);

		// Step 2
		double external_bias = 0;
		limiting_resistance_ = 0;
		bias_ = external_bias;
		double time = 0;
		unsigned int i = 0;

		const auto rate_fn = [this](const auto& from, const auto& to) {
			return mc_rate_fn(from, to);
		};

		std::vector<unsigned int> fil_shape(grid_size_z_);

		std::vector<double> biases;
		std::vector<double> currents;
		std::vector<double> time_steps;

		std::vector<double> lim_res(5, 0);

		int sweep = 0;

		const auto print_table_header = []
		{
			std::cout
				<< "-----------------------------------------------------------------------------------------------------------\n"
				<< "#    Time      Fil.vol   Ext.bias  Syst.bias Syst.res. Lim.res.  Current   MC est.      MC act.   Num.vacs.\n"
				<< "     [msec]    [arb.u.]  [V]       [V]       [kOhm]    [kOhm]    [uA]      [usec/step]  [steps]            \n"
				<< "-----------------------------------------------------------------------------------------------------------"
				<< std::endl;
		};

		std::cout << std::setprecision(4);
		print_table_header();

		while (true)
		{
			std::cout << std::left << std::setw(4) << i + 1 << ' ' << std::setw(9)
					  << es_util::au::to_sec(time) / 1e-3 << ' ' << std::flush;

			// Step 3
			compute_filament_shape(fil_shape);
			auto fil_size = 0u;
			for (auto r : fil_shape)
				fil_size += r * r;

			// Step 4
			compute_potential_and_heat(fil_shape);
			std::cout << std::setw(9) << fil_size << ' ' << std::setw(9)
					  << es_util::au::to_volt(external_bias) << ' ' << std::setw(9)
					  << es_util::au::to_volt(bias_) << ' ' << std::setw(9)
					  << es_util::au::to_ohm(bias_ / current_) / 1e3 << ' ' << std::setw(9)
					  << es_util::au::to_ohm(limiting_resistance_) / 1e3 << ' ' << std::setw(9)
					  << es_util::au::to_amp(current_) / 1e-6 << ' ' << std::flush;

			// Step 5
			heat_solver.solve();
			interpolate(heat_solver.solution_view(), mc_temp_);

			poisson_solver.solve();
			interpolate(poisson_solver.solution_view(), mc_potential_);

			// mc_temp_.write("t.mat");
			// mc_potential_.write("p.mat");
			// return;

			// Step 6
			double time_step = params::min_step_duration;
			double est_step_duration = mc_.estimate_step_duration(rate_fn);
			std::cout << std::setw(9) << es_util::au::to_sec(est_step_duration) * 1e6 << "    "
					  << std::flush;

			if (est_step_duration < time_step)
			{
				// Step 7
				const auto mc_res = mc_.run(params::steps_per_round, time_step, rate_fn);
				time_step = mc_res.second;

				std::cout << std::setw(9) << mc_res.first << std::flush;
			}
			else
				std::cout << std::setw(9) << 0 << std::flush;

			std::cout << ' ' << std::setw(9) << mc_.n_occupied() << std::endl;

			{
				biases.push_back(es_util::au::to_volt(bias_));
				currents.push_back(es_util::au::to_amp(current_) * 1e6);
				time_steps.push_back(time_step);

				la::Matfile_writer mw("iv.mat");
				mw.write("v", biases);
				mw.write("i", currents);
				mw.write("t", time_steps);

				const auto n_steps_write = 1;
				if (i % n_steps_write == 0/*  && est_step_duration < time_step */)
				{
					mc_.write(fname("m", i / n_steps_write));
					poisson_solver.write(fname("p", i / n_steps_write));
					heat_solver.write(fname("h", i / n_steps_write));
				}
			}

			// Step 8
			external_bias += (sweep == 1 ? -1 : 1) * time_step * params::bias_sweep_rate;

			// if (current_ >= params::max_current)
			// 	limiting_resistance_ = external_bias / params::max_current - bias_ / current_;
			// else
			// 	limiting_resistance_ = 0;

			// const auto avg = std::accumulate(lim_res.begin(), lim_res.end(), 0.) / 5;
			// lim_res.erase(lim_res.begin());
			// lim_res.push_back(limiting_resistance_);

			// bias_ = external_bias - params::max_current * limiting_resistance_;
			// bias_ = external_bias - params::max_current * limiting_resistance_;

			if (external_bias >= 0)
			{
				auto z = external_bias / params::max_current * current_ - bias_;
				if (z < 0)
					z = 0;
				bias_ = external_bias - z;
				limiting_resistance_ = z / current_;
			}
			else
			{
				auto z = external_bias / params::max_current * current_ + bias_;
				if (z < 0)
					z = 0;
				bias_ = external_bias + z;
				limiting_resistance_ = -z / current_;
			}

			time += time_step;
			++i;

			if (external_bias >= params::max_bias)
				sweep = 1;
			if (sweep == 1 && external_bias < -params::max_bias)
				sweep = 2;
			if (sweep == 2 && external_bias > 0)
				break;

			if (i % 20 == 0)
				print_table_header();
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
			const auto z_max = std::min(
				grid_size_z_ - 1, static_cast<unsigned int>(
									  std::ceil((br.top() + es_fe::delta) / params::grid_spacing)));

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

	double mc_rate_fn(const Point& from, const Point& to)
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

	double bias_;
	double limiting_resistance_;
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
