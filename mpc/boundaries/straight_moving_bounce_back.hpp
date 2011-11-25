#ifndef SRD_SERVICES_STRAIGHT_MOVING_BOUNCE_BACK_HPP_INCLUDED
#define SRD_SERVICES_STRAIGHT_MOVING_BOUNCE_BACK_HPP_INCLUDED

#include <cstdio>
#include <iostream>

#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/point_traits.hpp>
#include <mpc/mass_point.hpp>


MPC_NAMESPACE_BEGIN


/*
 * ToDo : check
 */
template< class Point , size_t which , bool min_or_max >
class straight_moving_bounce_back
{
public:

	typedef point_traits< Point > point_traits;
	const static size_t dim = point_traits::dim;
	typedef typename point_traits::value_type value_type;
	typedef typename point_traits::point_type point_type;
	typedef mass_point< point_type > mass_point_type;


	straight_moving_bounce_back( value_type bound = 0.0 , value_type dt = 0.025 ) : m_bound( bound ) , m_dt( dt ) { }

	FUNC_DECL void operator()( mass_point_type &m , const point_type &bound_vel ) const
	{
		point_type &p = m.coor;
		if( larger( m_bound , p[which] ) )
		{
			value_type frac = diff( m_bound , p[which] ) / m_dt / m.vel[which];
			check_frac_range( frac );
			frac *= 2.0 * m_dt;
			for( size_t i=0 ; i<dim ; ++i )
			{
				if( i == which ) p[i] = 2.0 * m_bound - p[i];
				else p[i] -= frac * m.vel[i];
			}
			m.vel = -m.vel;
			m.vel += bound_vel;
		}
	}

	//
	// geters
	//
	value_type get_bound( void ) const { return m_bound; }
	value_type get_dt( void ) const { return m_dt; }

	//
	// modifiers
	//
	void set_bound( value_type bound ) { m_bound = bound; }
	void set_dt( value_type dt ) { m_dt = dt; }

protected:

	value_type m_bound , m_dt;

	FUNC_DECL bool larger( value_type one , value_type two ) const
	{
		return  min_or_max ? ( one > two ) : ( one < two );
	};

	FUNC_DECL value_type diff( value_type one , value_type two ) const
	{
		return min_or_max ? two - one : one - two;
	}

	void check_frac_range( value_type &frac ) const
	{
		if( frac <= 0.0 )
		{
//			char str[512] = "";
//			sprintf( str , "check_frac_range() : frac < 0 : frac = %f" , frac );
			frac = 0.0;
			//            std::clog << str << std::endl;
		}
		if( frac >= 1.0 )
		{
//			char str[512] = "";
//			sprintf( str , "check_frac_range() : frac > 1 : frac = %f" , frac );
			frac = 1.0;
			//            std::clog << str << std::endl;
		}
	}
};




/*
 * ToDo : check
 */
template< class Point , size_t which , bool min_or_max >
class straight_moving_bounce_back_periodic : public straight_moving_bounce_back< Point , which , min_or_max >
{
	typedef straight_moving_bounce_back< Point , which , min_or_max > base_type;

public:

	typedef typename base_type::value_type value_type;
	typedef typename base_type::point_type point_type;
	typedef typename base_type::mass_point_type mass_point_type;

	straight_moving_bounce_back_periodic( value_type amplitude , value_type omega , value_type bound = 0.0 , value_type dt = 0.025 )
	: base_type( bound , dt ) , m_amplitude( amplitude ) , m_omega( omega )
	{ }

	FUNC_DECL void operator()( mass_point_type &m , value_type t ) const
	{
		base_type::operator()( m , point_type( m_amplitude * cos( m_omega * t ) , 0.0 ) );
	}

private:

	value_type m_amplitude;
	value_type m_omega;
};









template< class Point , size_t which , bool min_or_max >
class straight_moving_bounce_back_constant : public straight_moving_bounce_back< Point , which , min_or_max >
{
	typedef straight_moving_bounce_back< Point , which , min_or_max > base_type;

public:

	typedef typename base_type::value_type value_type;
	typedef typename base_type::point_type point_type;
	typedef typename base_type::mass_point_type mass_point_type;

	straight_moving_bounce_back_constant( value_type velocity , value_type bound = 0.0 , value_type dt = 0.025 )
	: base_type( bound , dt ) , m_velocity( velocity )
	{ }

	void operator()( mass_point_type &m , value_type t ) const
	{
		base_type::operator()( m , point_type( m_velocity , 0.0 ) );
	}

private:

	value_type m_velocity;
};


MPC_NAMESPACE_END

#endif // SRD_SERVICES_STRAIGHT_MOVING_BOUNCE_BACK_HPP_INCLUDED
