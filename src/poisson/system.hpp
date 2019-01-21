#pragma once
#include "boundary_cond.hpp"

#include <es_fe/var.hpp>
#include <es_fe/var_list.hpp>
#include <es_fe/system.hpp>
#include <es_fe/element/lagrange.hpp>
#include <es_fe/dof/dof_mapper.hpp>

#include <string>

using Poisson_element = es_fe::Lagrange<1, 2>;
using Poisson_var = es_fe::Var<
	Poisson_element,
	1,
	Poisson_dirichlet_const,
	Poisson_dirichlet_const,
	Poisson_dirichlet_core>;

class Poisson_system final : public es_fe::System<es_fe::Var_list<Poisson_var>, es_fe::Dof_mapper>
{
private:
	using Base = es_fe::System<es_fe::Var_list<Poisson_var>, es_fe::Dof_mapper>;

public:
	using Base::Base;

	virtual std::string name() const override
	{
		return "2D Poisson solver";
	}
};
