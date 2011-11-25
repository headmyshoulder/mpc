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
#include <utility>
#include <stdexcept>
#include <fstream>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/random.hpp>

// thrust includes
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>

// fusion includes
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/container.hpp>

// cuda include
#include <curand.h>
#include <curand_kernel.h>

// algebra includes
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>

#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/mass_point.hpp>
#include <mpc/boundaries.hpp>
#include <mpc/streaming.hpp>
#include <mpc/indexer.hpp>
#include <mpc/collision.hpp>
#include <mpc/analysis.hpp>

using namespace std;
using namespace mpc2;
namespace fusion = boost::fusion;

using boost::numeric::odeint::range_algebra;
using boost::numeric::odeint::thrust_algebra;


const static size_t dim = 2;
typedef double value_type;
typedef point< value_type , dim > point_type;
typedef mass_point< point_type > mass_point_type;

typedef std::vector< mass_point_type > mass_point_vector_type1;
typedef std::vector< point_type > point_vector_type1;
typedef std::vector< size_t > cell_type1;
typedef std::vector< point< size_t , dim > > index_type1;
typedef std::vector< curandState > rng_state_type1;


typedef thrust::host_vector< mass_point_type > mass_point_vector_type2;
typedef thrust::host_vector< point_type > point_vector_type2;
typedef thrust::host_vector< size_t > cell_type2;
typedef thrust::host_vector< point< size_t , dim > > index_type2;
typedef thrust::host_vector< curandState > rng_state_type2;


typedef thrust::device_vector< mass_point_type > mass_point_vector_type3;
typedef thrust::device_vector< point_type > point_vector_type3;
typedef thrust::device_vector< size_t > cell_type3;
typedef thrust::device_vector< point< size_t , dim > > index_type3;
typedef thrust::device_vector< curandState > rng_state_type3;


const point_type null_point( 0.0 , 0.0 );
const point_type spacing( 1.0 , 1.0 );

const point_type range( 30.0 , 15.0 );
// const point_type range( 15.0 , 8.0 );
const value_type density = 5;
const size_t n_solvent = size_t( range[0] * range[1] / spacing[0] / spacing[1] * density );
const size_t n_virtual_particles = range[0] / spacing[0] * density;
const size_t n = n_solvent + n_virtual_particles;

const point_type gravitation( 0.001 , 0.0 );


const value_type sqrt_kbT_mass = 1.0;
const value_type mass = 1.0;
const value_type dt = 0.025;




template
<
	class MassPointVector ,
	class PointVector ,
	class BC ,
	class Indexer ,
	class IndexVector ,
	class Rng ,
	class Value
>
void test_poiseuille
(
	MassPointVector &mp ,
	PointVector &force ,
	size_t num_of_particles , size_t num_of_virtual_particles ,
	BC &bc ,
	Indexer &indexer ,
	IndexVector &indices , IndexVector &permutated_indices ,
	IndexVector &begins ,
	IndexVector &ends ,
	Rng &rng ,
	Value &t , Value dt , Value mass ,
	size_t num_of_streaming_steps
)
{
	typedef streaming< thrust_algebra , cuda_streaming_operations > streaming_type;

	size_t num_of_all_particles = mp.size() , num_of_cells = indexer.number_of_cells();

	if( ( num_of_all_particles != ( num_of_particles + num_of_virtual_particles ) ) ||
		( force.size() != num_of_all_particles ) ||
		( indices.size() != num_of_all_particles ) ||
		( begins.size() != num_of_cells ) ||
		( ends.size() != num_of_cells ) ||
		( rng.size() != num_of_cells ) )
	{
		throw std::invalid_argument( "Sizes does not match!" );
	}

	for( size_t i=0 ; i<num_of_streaming_steps ; ++i , t+=dt )
	{
		streaming_type::step_a( std::make_pair( mp.begin() , mp.begin() + num_of_particles ) , force , t , dt , mass , bc );
		streaming_type::step_b( std::make_pair( mp.begin() , mp.begin() + num_of_particles ) , force , dt , mass );
//		streaming_type::step_a( mp , force , t , dt , mass , bc );
//		streaming_type::step_b( mp , force , dt , mass );
	}


	thrust::counting_iterator< size_t > search_begin( 0 );
	thrust::copy( search_begin , search_begin + num_of_all_particles , permutated_indices.begin() );

	thrust::transform( mp.begin() , mp.end() , indices.begin() , make_indexer_caller( indexer ) );
	thrust::sort_by_key( indices.begin() , indices.end() , permutated_indices.begin() );
	thrust::lower_bound( indices.begin() , indices.end() , search_begin , search_begin + num_of_cells , begins.begin() );
	thrust::upper_bound( indices.begin() , indices.end() , search_begin , search_begin + num_of_cells , ends.begin() );

	thrust::for_each(
			thrust::make_zip_iterator( thrust::make_tuple( begins.begin() , ends.begin() , rng.begin() ) ) ,
			thrust::make_zip_iterator( thrust::make_tuple( begins.end() , ends.end() , rng.end() ) ) ,
			anderson_plus_a_permutating< MassPointVector , IndexVector >( mp , permutated_indices )
			);
}







int main( int argc , char **argv )
{
	srand48( 325345777 );
	typedef boost::mt19937 rng_type;
	rng_type rng;


	typedef mpc2::regular_indexer< point_type ,
			with_start_point_tag ,
			check_and_throw_index_policy > regular_indexer_type;
	regular_indexer_type indexer( range , spacing );

	indexer.automatic_grid_shift( drand48 );
	clog << "Indexer start point : " << indexer.start() << endl;



	// initialize the solvent particles
	mass_point_vector_type1 solvent1;
	point_vector_type1 force1;
	for( size_t i=0 ; i<n ; ++i )
	{
		mass_point_type mp;
		mp.coor[0] = drand48() * range[0];
		mp.coor[1] = drand48() * range[1];
		maxwell_velocity( mp.vel , sqrt_kbT_mass , rng );
		solvent1.push_back( mp );
		force1.push_back( gravitation );
	}


	mass_point_vector_type2 solvent2 = solvent1;
	mass_point_vector_type3 solvent3 = solvent1;
	point_vector_type2 force2 = force1;
	point_vector_type3 force3 = force1;

	cell_type1 indices1( n ) , permutated_indices1( n );
	cell_type2 indices2( n ) , permutated_indices2( n );
	cell_type3 indices3( n ) , permutated_indices3( n );

	cell_type1 begins1( indexer.number_of_cells() ) , ends1( indexer.number_of_cells() );
	cell_type2 begins2( indexer.number_of_cells() ) , ends2( indexer.number_of_cells() );
	cell_type3 begins3( indexer.number_of_cells() ) , ends3( indexer.number_of_cells() );

	rng_state_type1 rng1( indexer.number_of_cells() );
	for( size_t i=0 ; i<indexer.number_of_cells() ; ++i ) curand_init( 1856234 , 0 , 0 , &rng1[i] );
	rng_state_type2 rng2 = rng1;
	rng_state_type3 rng3 = rng1;

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


	value_type t1 = 0.0 , t2 = 0.0 , t3 = 0.0;
//	test_poiseuille( solvent1 , force1 , n_solvent , n_virtual_particles , bc , indexer , indices1 , permutated_indices1 , begins1 , ends1 , rng1 , t1 , dt , mass , 4 );
//	test_poiseuille( solvent2 , force2 , n_solvent , n_virtual_particles , bc , indexer , indices2 , permutated_indices2 , begins2 , ends2 , rng2 , t2 , dt , mass , 4 );
//	test_poiseuille( solvent3 , force3 , n_solvent , n_virtual_particles , bc , indexer , indices3 , permutated_indices3 , begins3 , ends3 , rng3 , t3 , dt , mass , 4 );

    mpc2::velocity_dist_y< mass_point_type >	vel_dist( indexer.dimensions()[1] , range[1] , 0.0 , range[0] );

//	cout << "set term x11" << endl;
//	cout << "unset key" << endl;
//	cout << "set size ratio " << range[1] / range[0] << endl;
	for( size_t i=0 ; i<100000000 ; ++i )
	{
		// cout << "p [-1:" << range[0] + 1 << "][-1:" << range[1] + 1 << "] '-' w vec \n";
//		cout << "p [-1:" << range[0] + 1 << "][-1:" << range[1] + 1 << "] '-' pt 7 \n";
//		for( size_t i=0 ; i<n_solvent ; ++i )
//			cout << solvent1[i].coor << "\t" << solvent1[i].vel << "\n";
//		cout << "e" << endl;

//		double a = 1.001 * double(i);
//		for( size_t i=0 ; i<10000000 ; ++i ) a *= 1.0001;
//		clog << a << endl;

		vel_dist.add( solvent3 );

		test_poiseuille( solvent3 , force3 , n_solvent , n_virtual_particles , bc , indexer , indices3 , permutated_indices3 , begins3 , ends3 , rng3 , t3 , dt , mass , 4 );

		if( i  && ( !(i%100)))
		{

			clog << i << endl;
			ofstream fout( "dat/vel_dist.dat" );
			vel_dist.write( fout );
		}

	}


//	mass_point_vector_type2 solvent3_host = solvent3;
//
//	for( size_t i=0 ; i<indexer.number_of_cells() ; ++i )
//	{
//		cout << i << "\t";
//		cout << begins1[i] << "\t" << ends1[i] << "\t" << solvent1[i].coor << "\t" << solvent1[i].vel << "\t";
//		cout << begins2[i] << "\t" << ends2[i] << "\t" << solvent2[i].coor << "\t" << solvent2[i].vel  << "\t";
//		cout << begins3[i] << "\t" << ends3[i] << "\t" << solvent3_host[i].coor << "\t" << solvent3_host[i].vel  << "\n";
//	}

	return 0;
}

