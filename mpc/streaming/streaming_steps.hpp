/*
 * streaming.hpp
 *
 *  Created on: Dec 8, 2010
 *      Author: karsten
 */

#ifndef STREAMING_STEPS_HPP_
#define STREAMING_STEPS_HPP_

#include <mpc/defines.hpp>

MPC_NAMESPACE_BEGIN

/*
 * The const version exist for working with boost.range
 */
template< class Algebra , class StreamingOperations >
struct streaming
{

	template< class MassPoint , class Force , class Value , class BC >
	static void step_a( MassPoint &mp , Force &f , Value t , Value dt , Value mass , BC bc )
	{
		Algebra::for_each2( mp , f , StreamingOperations::make_streaming_a( t , dt , mass , bc ) );
	}

	template< class MassPoint , class Force , class Value , class BC >
	static void step_a( const MassPoint &mp , const Force &f , Value t , Value dt , Value mass , BC bc )
	{
		Algebra::for_each2( mp , f , StreamingOperations::make_streaming_a( t , dt , mass , bc ) );
	}






	template< class MassPoint , class Force , class Value >
	static void step_b( MassPoint &mp , Force &f , const Value dt , const Value mass )
	{
		Algebra::for_each2( mp , f , StreamingOperations::make_streaming_b( dt , mass ) );
	}

	template< class MassPoint , class Force , class Value >
	static void step_b( const MassPoint &mp , const Force &f , const Value dt , const Value mass )
	{
		Algebra::for_each2( mp , f , StreamingOperations::make_streaming_b( dt , mass ) );
	}

};

MPC_NAMESPACE_END

#endif /* STREAMING_STEPS_HPP_ */
