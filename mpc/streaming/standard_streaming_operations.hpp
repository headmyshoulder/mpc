/*
 * standard_streaming_operations.hpp
 *
 *  Created on: Dec 8, 2010
 *      Author: karsten
 */

#ifndef STANDARD_STREAMING_OPERATIONS_HPP_
#define STANDARD_STREAMING_OPERATIONS_HPP_

#include <boost/fusion/algorithm.hpp>
#include <boost/bind.hpp>
#include <boost/bind/apply.hpp>
#include <boost/ref.hpp>

#include <mpc/defines.hpp>

MPC_NAMESPACE_BEGIN



/*
 * instead of the struct cuda_streaming_operations can used, they work as well with std::vector ...
 */
struct standard_streaming_operations
{
	template< class Value , class BC >
	struct streaming_a
	{
		Value m_t , m_dt , m_mass;
		Value m_dtsq2_m , m_dt2_m;
		BC &m_bc;

		streaming_a( Value t , Value dt , Value mass , BC &bc )
		: m_t( t ) , m_dt( dt ) , m_mass( mass ) ,
		  m_dtsq2_m( 0.5 * dt *dt / mass ) , m_dt2_m( 0.5 * dt / mass ) ,
		  m_bc( bc ) { }

		template< class MassPoint , class Point >
		void operator()( MassPoint &mp , Point &f )
		{
			mp.coor += m_dt * mp.vel + m_dtsq2_m * f;
			mp.vel += m_dt2_m * f;
			boost::fusion::for_each( m_bc , boost::bind( boost::apply< void >() , _1 , boost::ref( mp ) , m_t ) );
		}
	};

	template< class Value , class BC >
	static streaming_a< Value , BC > make_streaming_a( Value t , Value dt , Value mass , BC &bc )
	{
		return streaming_a< Value , BC >( t , dt , mass , bc );
	}




	template< class Value >
	struct streaming_b
	{
		Value m_dt , m_mass;
		Value m_dt2_m;

		streaming_b( Value dt , Value mass )
		: m_dt( dt ) , m_mass( mass ) , m_dt2_m( 0.5 * dt / mass ) { }

		template< class MassPoint , class Point >
		void operator()( MassPoint &mp , const Point &f ) const
		{
			mp.vel += m_dt2_m * f;
		}
	};

	template< class Value >
	static streaming_b< Value > make_streaming_b( Value dt , Value mass )
	{
		return streaming_b< Value >( dt , mass );
	}
};

MPC_NAMESPACE_END

#endif /* STANDARD_STREAMING_OPERATIONS_HPP_ */
