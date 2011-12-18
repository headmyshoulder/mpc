/*
 * streaming.cpp
 *
 *  Created on: Dec 6, 2010
 *      Author: karsten
 */

#include <iostream>
#include <vector>

// thrust includes
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

// fusion includes
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/container.hpp>

// algebra includes
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>

#include <mpc/mass_point.hpp>
#include <mpc/maxwell_velocity.hpp>
#include <mpc/boundaries.hpp>
#include <mpc/streaming.hpp>
#include <mpc/drand48_generator.hpp>


using namespace std;
using namespace mpc2;
namespace fusion = boost::fusion;

using boost::numeric::odeint::range_algebra;
using boost::numeric::odeint::thrust_algebra;

template< class Algebra , class StreamingOperations >
struct test_streaming
{
	typedef streaming< Algebra , StreamingOperations > streaming_type;

	template< class MassPoints , class Force , class Value , class BC >
	void operator()( MassPoints &mp , Force &f , Value &t , const Value dt , const Value mass , BC &bc , size_t num_of_streaming_steps )
	{
		for( size_t i=0 ; i<num_of_streaming_steps ; ++i , t+=dt )
		{
			streaming_type::step_a( mp , f , t , dt , mass , bc );
			streaming_type::step_b( mp , f , dt , mass );
		}
	}
};

const static size_t dim = 2;
typedef double value_type;
typedef point< value_type , dim > point_type;
typedef mass_point< point_type > mass_point_type;

typedef std::vector< mass_point_type > mass_point_vector_type1;
typedef std::vector< point_type > point_vector_type1;

typedef thrust::host_vector< mass_point_type > mass_point_vector_type2;
typedef thrust::host_vector< point_type > point_vector_type2;

typedef thrust::device_vector< mass_point_type > mass_point_vector_type3;
typedef thrust::device_vector< point_type > point_vector_type3;


const size_t n_solvent = 32;
const point_type range( 16.0 , 16.0 );
const point_type null_point( 0.0 , 0.0 );
const value_type sqrt_kbT_mass = 1.0;
const value_type mass = 1.0;


int main( int argc , char **argv )
{
	srand48( 325345777 );
	drand48_generator rng;


	// initialize the solvent particles
	mass_point_vector_type1 solvent1;
	point_vector_type1 force1;
	for( size_t i=0 ; i<n_solvent ; ++i )
	{
		mass_point_type mp;
		mp.coor[0] = drand48() * range[0];
		mp.coor[1] = drand48() * range[1];
		maxwell_velocity( mp.vel , sqrt_kbT_mass , rng );
		solvent1.push_back( mp );
		force1.push_back( point_type( 0.0 , -0.1 ) );
	}

	mass_point_vector_type2 solvent2 = solvent1;
	point_vector_type2 force2 = force1;

	mass_point_vector_type3 solvent3 = solvent1;
	point_vector_type3 force3 = force1;

	value_type t1 = 0.0 , t2 = 0.0 , t3 = 0.0;
	const value_type dt = 0.0025;

	fusion::vector
	<
		mpc2::periodic_boundary< value_type , 0 > ,
		mpc2::straight_bounce_back< point_type , 1 , true > ,
		mpc2::straight_bounce_back< point_type , 1 , false > ,
		mpc2::fixed_boundary_corrector< point_type >
	>
	bc
	(
		mpc2::periodic_boundary< value_type , 0 >( 0.0 , range[0] ) ,
		mpc2::straight_bounce_back< point_type , 1 , true >( 0.0 , dt ) ,
		mpc2::straight_bounce_back< point_type , 1 , false >( range[1] , dt ) ,
		mpc2::fixed_boundary_corrector< point_type >( null_point , range )
	);

	cout << "unset key" << endl;
	for( size_t i=0 ; i<10000 ; ++i )
	{
		test_streaming< range_algebra , standard_streaming_operations >()( solvent1 , force1 , t1 , dt , mass , bc , 4 );
		test_streaming< range_algebra , standard_streaming_operations >()( solvent2 , force2 , t2 , dt , mass , bc , 4 );
		test_streaming< thrust_algebra , cuda_streaming_operations >()( solvent3 , force3 , t3 , dt , mass , bc , 4 );
		mass_point_vector_type2 vec = solvent3;

		cout << "p '-' u 1:2 pt 7,'-' u 3:4 pt 7,'-' u 5:6 pt 7" << endl;
		for( size_t i=0 ; i<solvent1.size() ; ++i )
			cout << solvent1[i].coor << "\t" << solvent2[i].coor << "\t" << vec[i].coor << "\n";
		cout << "e" << endl;
	}



//	for( size_t i=0 ; i<1000 ; ++i )
//	{
//		test_streaming< standard_algebra , standard_streaming_operations >()( solvent1 , force1 , t1 , dt , mass , bc , 4 );
//		cout << "p '-' u 1:2 pt 7" << endl;
//		for( size_t i=0 ; i<solvent1.size() ; ++i )
//			cout << solvent1[i].coor << "\n";
//		cout << "e" << endl;
//	}
//	getchar();
//
//	for( size_t i=0 ; i<1000 ; ++i )
//	{
//		test_streaming< standard_algebra , standard_streaming_operations >()( solvent2 , force2 , t2 , dt , mass , bc , 4 );
//		cout << "p '-' u 1:2 pt 7" << endl;
//		for( size_t i=0 ; i<solvent2.size() ; ++i )
//			cout << solvent2[i].coor << "\n";
//		cout << "e" << endl;
//	}
//	getchar();
//
//	for( size_t i=0 ; i<1000 ; ++i )
//	{
//		test_streaming< thrust_algebra , cuda_streaming_operations >()( solvent3 , force3 , t3 , dt , mass , bc , 4 );
//		mass_point_vector_type2 vec = solvent3;
//		cout << "p '-' u 1:2 pt 7" << endl;
//		for( size_t i=0 ; i<vec.size() ; ++i )
//			cout << vec[i].coor << "\n";
//		cout << "e" << endl;
//	}
//	getchar();

	return 0;
}

