#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <fstream>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

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

#include <omp.h>

// algebra includes
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>

#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/mass_point.hpp>
#include <mpc/boundaries.hpp>
#include <mpc/streaming.hpp>
#include <mpc/indexer.hpp>
#include <mpc/collision.hpp>
#include <mpc/analysis.hpp>
#include <mpc/drand48_generator.hpp>

using namespace std;
using namespace mpc2;
namespace fusion = boost::fusion;

using boost::numeric::odeint::thrust_algebra;


const static size_t dim = 2;
typedef double value_type;
typedef point< value_type , dim > point_type;
typedef mass_point< point_type > mass_point_type;

typedef VECTOR< mass_point_type > mass_point_vector_type;
typedef VECTOR< point_type > point_vector_type;
typedef VECTOR< size_t > cell_type;
typedef VECTOR< point< size_t , dim > > index_type;
typedef VECTOR< curandState > rng_state_type;


const point_type null_point( 0.0 , 0.0 );
const point_type spacing( 1.0 , 1.0 );

RANGE
DENSITY


// const point_type range( 15.0 , 8.0 );
// const value_type density = 5;
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
	OMP_NUM_THREADS

	srand48( 325345777 );
	typedef drand48_generator rng_type;
	rng_type rng;


	typedef mpc2::regular_indexer< point_type ,
			with_start_point_tag ,
			check_and_throw_index_policy > regular_indexer_type;
	regular_indexer_type indexer( range , spacing );

	indexer.automatic_grid_shift( drand48 );
//	clog << "Indexer start point : " << indexer.start() << endl;



	// initialize the solvent particles
	std::vector< mass_point_type > tmp_solvent;
	std::vector< point_type > tmp_force;
	for( size_t i=0 ; i<n ; ++i )
	{
		mass_point_type mp;
		mp.coor[0] = drand48() * range[0];
		mp.coor[1] = drand48() * range[1];
		maxwell_velocity( mp.vel , sqrt_kbT_mass , rng );
		tmp_solvent.push_back( mp );
		tmp_force.push_back( gravitation );
	}
	mass_point_vector_type solvent = tmp_solvent;
	point_vector_type force = tmp_force;


	cell_type indices( n ) , permutated_indices( n );
	cell_type begins( indexer.number_of_cells() ) , ends( indexer.number_of_cells() );

	std::vector< curandState > tmp_rng_state( indexer.number_of_cells() );
	for( size_t i=0 ; i<indexer.number_of_cells() ; ++i ) curand_init( 1856234 , 0 , 0 , &tmp_rng_state[i] );
	rng_state_type rng_state = tmp_rng_state;

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


	value_type t = 0.0 ;

	for( size_t i=0 ; i<200000 ; ++i )
	{
		test_poiseuille( solvent , force , n_solvent , n_virtual_particles , bc , indexer , indices , permutated_indices , begins , ends , rng_state , t , dt , mass , 4 );
		indexer.automatic_grid_shift( drand48 );
	}

	ofstream fout( OUTFILE );
	thrust::host_vector< mass_point_type > tmp = solvent;
	for( size_t i=0 ; i<tmp.size() ; ++i )
		fout << tmp[i].coor << "\t" << tmp[i].vel << "\n";

	return 0;
}
