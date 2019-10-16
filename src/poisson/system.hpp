#pragma once
#include "boundary_cond.hpp"
#include "element.hpp"

#include <esf/var.hpp>
#include <esf/var_list.hpp>
#include <esf/system.hpp>
#include <esf/dof/dof_mapper.hpp>

#include <string>

using Poisson_var = esf::Var<
	Poisson_element,
	1,
	Poisson_dirichlet_contact,
	Poisson_dirichlet_contact,
	Poisson_dirichlet_core>;

class Poisson_system final : public esf::System<esf::Var_list<Poisson_var>, esf::Dof_mapper>
{
private:
	using Base = esf::System<esf::Var_list<Poisson_var>, esf::Dof_mapper>;

public:
	using Base::Base;

	virtual std::string name() const override
	{
		return "2D Poisson solver";
	}
};
