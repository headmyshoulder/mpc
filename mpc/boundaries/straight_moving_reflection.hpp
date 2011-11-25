#ifndef STRAIGHT_MOVING_REFLECTION_HPP_INCLUDED
#define STRAIGHT_MOVING_REFLECTION_HPP_INCLUDED

#include <cstdio>
#include <iostream>

#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/point_traits.hpp>
#include <mpc/mass_point.hpp>

MPC_NAMESPACE_BEGIN



template< class Point , size_t which , bool min_or_max >
class straight_moving_reflection
{
public:

	typedef point_traits< Point > point_traits;
	const static size_t dim = point_traits::dim;
	typedef typename point_traits::value_type value_type;
	typedef typename point_traits::value_type point_type;
	typedef mass_point< point_type > mass_point_type;


	straight_moving_reflection( value_type bound = 0.0 ) : m_bound( bound ) { }

	FUNC_DECL void operator()( mass_point_type &m , const point_type &bound_vel ) const
	{
		point_type &p = m.coor;
		if( larger( m_bound , p[which] ) )
		{
			p[which] = 2.0 * m_bound - p[which];
			m.vel[which] = -m.vel[which];
			m.vel += bound_vel;
		}
	}

	FUNC_DECL void operator()( mass_point_type &m , value_type t ) const
	{
		this->operator()( m );
	}


	//
	// getters
	//
	value_type get_bound( void ) const { return m_bound; }

	//
	// modifiers
	//
	void set_bound( value_type bound ) { m_bound = bound; }

protected:

	value_type m_bound ;

	bool larger( value_type one , value_type two ) const
	{
		return  min_or_max ? ( one > two ) : ( one < two );
	};
};




template< class Point , size_t which , bool min_or_max >
class straight_moving_reflection_periodic : straight_moving_reflection< Point , which , min_or_max >
{
	typedef straight_moving_reflection< Point , which , min_or_max > base_type;

public:

	typedef typename base_type::value_type value_type;
	typedef typename base_type::point_type point_type;
	typedef typename base_type::mass_point_type mass_point_type;

	straight_moving_reflection_periodic( value_type amplitude , value_type omega , value_type bound = 0.0 )
	: base_type( bound ) , m_amplitude( amplitude ) , m_omega( omega )
	{ }

	void operator()( mass_point_type &m , value_type t ) const
	{
		base_type::operator()( m , point_type( m_amplitude * cos( m_omega * t ) , 0.0 ) );
	}

private:

	value_type m_amplitude;
	value_type m_omega;
};



template< class Point , size_t which , bool min_or_max >
class straight_moving_reflection_constant : straight_moving_reflection< Point , which , min_or_max >
{
	typedef straight_moving_reflection< Point , which , min_or_max > base_type;

public:

	typedef typename base_type::value_type value_type;
	typedef typename base_type::point_type point_type;
	typedef typename base_type::mass_point_type mass_point_type;

	straight_moving_reflection_constant( value_type velocity , value_type bound = 0.0 )
	: base_type( bound ) , m_velocity( velocity )
	{ }

	void operator()( mass_point_type &m , value_type t ) const
	{
		base_type::operator()( m , point_type( m_velocity , 0.0 ) );
	}

private:

	value_type m_velocity;
};


MPC_NAMESPACE_END



#endif // STRAIGHT_MOVING_REFLECTION_HPP_INCLUDED
