#pragma once
#include "boundary_cond.hpp"

#include <es_fe/var_list.hpp>
#include <es_fe/system.hpp>
#include <es_fe/element/lagrange.hpp>
#include <es_fe/dof/dof_mapper.hpp>

#include <string>

using Heat_element = es_fe::Lagrange<1, 2>;
using Heat_var = es_fe::Var<Heat_element, 1, Heat_dirichlet>;

class Heat_system final : public es_fe::System<es_fe::Var_list<Heat_var>, es_fe::Dof_mapper>
{
private:
	using Base = es_fe::System<es_fe::Var_list<Heat_var>, es_fe::Dof_mapper>;

public:
	using Base::Base;

	virtual std::string name() const override
	{
		return "2D heat solver";
	}
};
