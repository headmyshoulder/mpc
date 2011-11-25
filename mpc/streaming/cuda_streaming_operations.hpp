/*
 * cuda_streaming_operations.hpp
 *
 *  Created on: Dec 8, 2010
 *      Author: karsten
 */

#ifndef CUDA_STREAMING_OPERATIONS_HPP_
#define CUDA_STREAMING_OPERATIONS_HPP_

#include <thrust/convert_tuple.hpp>
#include <thrust/tuple_for_each.hpp>

#include <mpc/defines.hpp>

MPC_NAMESPACE_BEGIN


template< class MassPoint >
struct invoke_func
{
	MassPoint &m_mp;

        __device__ __host__
	invoke_func( MassPoint &mp ) : m_mp( mp ) { }

	template< class Func >
	__device__ __host__
	void operator()( Func f ) const
	{
		f( m_mp );
	}
};

template< class MassPoint >
__device__ __host__
invoke_func< MassPoint > make_invoke_func( MassPoint &mp )
{
	return invoke_func< MassPoint >( mp );
}

struct cuda_streaming_operations
{
	template< class Value , class BC >
	struct streaming_a
	{
		Value m_t , m_dt , m_mass;
		Value m_dtsq2_m , m_dt2_m;
		BC m_bc;

		streaming_a( Value t , Value dt , Value mass , const BC &bc )
		: m_t( t ) , m_dt( dt ) , m_mass( mass ) ,
		  m_dtsq2_m( 0.5 * dt *dt / mass ) , m_dt2_m( 0.5 * dt / mass ) ,
		  m_bc( bc ) { }

		template< class Tuple >
		__device__ __host__
		void operator()( Tuple t )
		{
			thrust::get< 0 >( t ).coor += m_dt * thrust::get< 0 >( t ).vel + m_dtsq2_m * thrust::get< 1 >( t) ;
			thrust::get< 0 >( t ).vel += m_dt2_m * thrust::get< 1 >( t );
			thrust::tuple_for_each( m_bc , make_invoke_func( thrust::get< 0 >( t ) ) );
		}
	};

	template< class Value , class BC >
	static streaming_a< Value , typename thrust::result_of::as_tuple< BC >::type >
	make_streaming_a( Value t , Value dt , Value mass , BC &bc )
	{
		typedef typename thrust::result_of::as_tuple< BC >::type bc_type;
		return streaming_a< Value , bc_type >( t , dt , mass , thrust::as_tuple( bc ) );
	}



	template< class Value >
	struct streaming_b
	{
		Value m_dt , m_mass;
		Value m_dt2_m;

		streaming_b( Value dt , Value mass )
		: m_dt( dt ) , m_mass( mass ) , m_dt2_m( 0.5 * dt / mass ) { }

		template< class Tuple >
		__device__ __host__
		void operator()( Tuple t ) const
		{
			thrust::get<0>( t ).vel += m_dt2_m * thrust::get< 1 >( t );
		}
	};

	template< class Value >
	static streaming_b< Value > make_streaming_b( Value dt , Value mass )
	{
		return streaming_b< Value >( dt , mass );
	}
};


MPC_NAMESPACE_END


#endif /* CUDA_STREAMING_OPERATIONS_HPP_ */
