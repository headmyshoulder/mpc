#ifndef STRAIGHT_REFLECTION_HPP_INCLUDED
#define STRAIGHT_REFLECTION_HPP_INCLUDED

#include <cstdio>
#include <iostream>

#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/point_traits.hpp>
#include <mpc/mass_point.hpp>

MPC_NAMESPACE_BEGIN



template< class Point , size_t which , bool min_or_max >
class straight_reflection
{
public:

	typedef point_traits< Point > point_traits;
	const static size_t dim = point_traits::dim;
	typedef typename point_traits::value_type value_type;
	typedef typename point_traits::point_type point_type;
	typedef mass_point< point_type > mass_point_type;


	straight_reflection( value_type bound = 0.0 ) : m_bound( bound ) { }

	FUNC_DECL void operator()( mass_point_type &m ) const
	{
		point_type &p = m.coor;
		if( larger( m_bound , p[which] ) )
		{
			p[which] = 2.0 * m_bound - p[which];
			m.vel[which] = -m.vel[which];
		}
	}

	FUNC_DECL void operator()( mass_point_type &m , value_type t ) const
	{
		this->operator()( m );
	}


	//
	// modifiers
	//
	value_type get_bound( void ) const { return m_bound; }

	//
	// getters
	//
	void set_bound( value_type bound ) { m_bound = bound; }

protected:

	value_type m_bound ;

	FUNC_DECL bool larger( value_type one , value_type two ) const
	{
		return  min_or_max ? ( one > two ) : ( one < two );
	};
};








MPC_NAMESPACE_END

#endif // STRAIGHT_REFLECTION_HPP_INCLUDED
