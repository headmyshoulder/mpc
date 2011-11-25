#ifndef SRD_SERVICES_PERIODIC_BOUNDARIES_HPP_INCLUDED
#define SRD_SERVICES_PERIODIC_BOUNDARIES_HPP_INCLUDED

#include <iostream>

#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/point_traits.hpp>
#include <mpc/mass_point.hpp>

MPC_NAMESPACE_BEGIN

template< class T , size_t which >
class periodic_boundary
{
public:

	FUNC_DECL FUNC_DECL periodic_boundary( T min = 0.0 , T max = 30.0 )
	: m_min( min ) , m_max( max ) { calc_shift(); }

	template< class MassPoint >
    FUNC_DECL void operator()( MassPoint &m ) const
    {
		typedef typename MassPoint::point_type point_type;
		typedef typename point_traits< point_type >::value_type value_type;

		point_type &p = m.coor;
		value_type diff1 = p[which] - m_min;
		value_type diff2 = p[which] - m_max;
		if( diff1 < 0.0 ) p[which] += m_shift;
		if( diff2 > 0.0 ) p[which] -= m_shift;
    }

	template< class MassPoint , class ValueType >
	FUNC_DECL void operator()( MassPoint &m , ValueType t ) const
	{
		this->operator()( m );
	}


	//
	// modifiers
	//
	FUNC_DECL void set_min( T min ) { m_min = min ; calc_shift() ; }
	FUNC_DECL void set_max( T max ) { m_max = max ; calc_shift() ; }
	FUNC_DECL void set_min_max( T min , T max ) { m_min = min ; m_max = max ; calc_shift(); }

    // getters
    FUNC_DECL T get_min( void ) const { return m_min; }
    FUNC_DECL T get_max( void ) const { return m_max; }

private:

    FUNC_DECL void calc_shift( void )
    {
    	m_shift = m_max - m_min ;
//    	if( m_shift < 0.0 )
//    		std::clog << "calc_shift() ; diff < 0.0" << std::endl;
	}


	T m_min , m_max , m_shift;
};


MPC_NAMESPACE_END


#endif //SRD_SERVICES_PERIODIC_BOUNDARIES_HPP_INCLUDED
