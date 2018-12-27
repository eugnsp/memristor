#pragma once
#include "boundary_cond.hpp"

#include <es_fe/var_list.hpp>
#include <es_fe/system.hpp>
#include <es_fe/element/lagrange.hpp>
#include <es_fe/dof/dof_mapper.hpp>

#include <string>

namespace poisson
{
using Element = es_fe::Lagrange<1, 2>;
using Var = es_fe::Var<Element, 1, Dirichlet_const, Dirichlet_const, Dirichlet_core>;

class System final : public es_fe::System<es_fe::Var_list<Var>, es_fe::Dof_mapper>
{
private:
	using Base = es_fe::System<es_fe::Var_list<Var>, es_fe::Dof_mapper>;

public:
	using Base::Base;

	virtual std::string name() const override
	{
		return "2D Poisson solver";
	}
};
} // namespace poisson
