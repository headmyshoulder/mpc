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

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

// thrust includes
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>

// cuda include
#include <curand.h>
#include <curand_kernel.h>

// algebra includes
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>


#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/mass_point.hpp>
#include <mpc/indexer.hpp>
#include <mpc/maxwell_velocity.hpp>
#include <mpc/collision.hpp>
#include <mpc/drand48_generator.hpp>

using namespace std;
using namespace mpc2;

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


const size_t n_solvent = 32;
const point_type range( 2.0 , 2.0 );
const point_type null_point( 0.0 , 0.0 );
const point_type spacing( 1.0 , 1.0 );
const value_type sqrt_kbT_mass = 1.0;



int main( int argc , char **argv )
{
	srand48( 325345777 );
	typedef drand48_generator rng_type;
	rng_type rng;


	typedef mpc2::regular_indexer< point_type ,
			with_start_point_tag ,
			check_and_throw_index_policy > regular_indexer_type;
	regular_indexer_type indexer( range , spacing );

	indexer.automatic_grid_shift( drand48 );
	clog << "Indexer start point : " << indexer.start() << endl;



	// initialize the solvent particles
	mass_point_vector_type1 solvent1;
	for( size_t i=0 ; i<n_solvent ; ++i )
	{
		mass_point_type mp;
		mp.coor[0] = drand48() * range[0];
		mp.coor[1] = drand48() * range[1];
		maxwell_velocity( mp.vel , sqrt_kbT_mass , rng );
		solvent1.push_back( mp );
	}

	mass_point_vector_type2 solvent2 = solvent1;
//	mass_point_vector_type3 solvent3 = solvent1;
	cell_type1 cell1( n_solvent );
	cell_type2 cell2( n_solvent );
//	cell_type3 cell3( n_solvent );

	cell_type1 begins1( indexer.number_of_cells() ) , ends1( indexer.number_of_cells() );
	cell_type2 begins2( indexer.number_of_cells() ) , ends2( indexer.number_of_cells() );
//	cell_type3 begins3( indexer.number_of_cells() ) , ends3( indexer.number_of_cells() );

	cout << 1 << endl;

	rng_state_type1 rng1( indexer.number_of_cells() );
	for( size_t i=0 ; i<indexer.number_of_cells() ; ++i ) curand_init( 1856234 , 0 , 0 , &rng1[i] );
	rng_state_type2 rng2 = rng1;
//	rng_state_type3 rng3 = rng1;


//	anderson( solvent1 , begins1 , ends1 , rng1 , indexer );
//	anderson( solvent2 , begins2 , ends2 , rng2 , indexer );
//	anderson( solvent3 , begins3 , ends3 , rng3 , indexer );
//
//	mass_point_vector_type2 solvent3_host = solvent3;
//
//	for( size_t i=0 ; i<indexer.number_of_cells() ; ++i )
//	{
//		cout << i << "\t";
//		cout << begins1[i] << "\t" << ends1[i] << "\t" << solvent1[i].coor << "\t" << solvent1[i].vel << "\t";
//		cout << begins2[i] << "\t" << ends2[i] << "\t" << solvent2[i].coor << "\t" << solvent2[i].vel  << "\t";
//		cout << begins3[i] << "\t" << ends3[i] << "\t" << solvent3_host[i].coor << "\t" << solvent3_host[i].vel  << "\n";
//	}

	cout << "set term x11" << endl;
	cout << "unset key" << endl;
	cout << "set size ratio " << range[1] / range[0] << endl;
	for( size_t i=0 ; i<100000 ; ++i )
	{
		cout << "p [-1:" << range[0] + 1 << "][-1:" << range[1] + 1 << "] '-' w vec \n";
		for( size_t i=0 ; i<n_solvent ; ++i )
			cout << solvent1[i].coor << "\t" << solvent1[i].vel << "\n";
		cout << "e" << endl;

		for( size_t i=0 ; i<rng1.size() ; ++i )
		{
			clog << "rng " << i << " " << rng1[i].d << " ";
			clog << rng1[i].v[0] << " " << rng1[i].v[1] << " " << rng1[i].v[2] << " " << rng1[i].v[3] << " " << rng1[i].v[4] << endl;
		}


		anderson( solvent1 , begins1 , ends1 , rng1 , indexer );
		getchar();
	}

	return 0;
}

