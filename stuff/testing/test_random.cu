/*
 * test_anderson.cu
 *
 *  Created on: Dec 15, 2010
 *      Author: karsten
 */


/*
 * streaming.cpp
 *
 *  Created on: Dec 6, 2010
 *      Author: karsten
 */

#include <iostream>
#include <vector>

// cuda include
#include <curand.h>
#include <curand_kernel.h>

// thrust includes
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/random.h>

#include <mpc/defines.hpp>

using namespace std;


struct gaussian
{
	template< class Tuple >
	FUNC_DECL
	void operator()( Tuple t )
	{
		curandState state = thrust::get< 0 >( t );
		thrust::get< 1 >( t ) = curand_normal_double( &state );
	}
};

int main( int argc , char **argv )
{
	typedef thrust::random::default_random_engine::result_type result_type;
	thrust::random::default_random_engine rng1( 10 ) , rng2( 10 );


//	result_type r10 = rng1() , r20 = rng2();
//	for( size_t i=0 ; i<10 ; ++i )
//	{
//		result_type r1 = rng1() , r2 = rng2();
//		cout << i << "\t" << r1 << "\t" << r2 << "\n";
//		if( i == 2 ) rng2.seed( r10 );
//		if( i == 5 ) rng2.seed( r1 );
//	}
//
//	curandState state;
//	curand_init( 123 , 1 , 1 , &state );
//
//	for( size_t i=0 ; i<10 ; ++i )
//		cout << curand( &state ) << "\n";

	const size_t n = 32;
	typedef std::vector< curandState > seed_vector_type1;
	typedef thrust::host_vector< curandState > seed_vector_type2;
	typedef thrust::device_vector< curandState > seed_vector_type3;

	typedef std::vector< double > value_vector_type1;
	typedef thrust::host_vector< double > value_vector_type2;
	typedef thrust::device_vector< double > value_vector_type3;

	seed_vector_type1 seeds1( n );
	for( size_t i=0 ; i<n ; ++i ) curand_init( 1856234 , 3 * i , 1 , &seeds1[i] );
	seed_vector_type2 seeds2 = seeds1;
	seed_vector_type3 seeds3 = seeds1;

	value_vector_type1 values1( n );
	value_vector_type2 values2( n );
	value_vector_type3 values3( n );


	thrust::for_each(
			thrust::make_zip_iterator( thrust::make_tuple( seeds1.begin() , values1.begin() ) ) ,
			thrust::make_zip_iterator( thrust::make_tuple( seeds1.end() , values1.end() ) ) ,
			gaussian() );

	thrust::for_each(
			thrust::make_zip_iterator( thrust::make_tuple( seeds2.begin() , values2.begin() ) ) ,
			thrust::make_zip_iterator( thrust::make_tuple( seeds2.end() , values2.end() ) ) ,
			gaussian() );


	thrust::for_each(
			thrust::make_zip_iterator( thrust::make_tuple( seeds3.begin() , values3.begin() ) ) ,
			thrust::make_zip_iterator( thrust::make_tuple( seeds3.end() , values3.end() ) ) ,
			gaussian() );

	value_vector_type2 values3_host = values3;
	for( size_t i=0 ; i<n ; ++i )
		cout << values1[i] << "\t" << values2[i] << "\t" << values3_host[i] << "\n";



	return 0;
}

