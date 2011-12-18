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

// thrust includes
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

// algebra includes
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>

#include <mpc/point.hpp>
#include <mpc/mass_point.hpp>
#include <mpc/indexer.hpp>
#include <mpc/maxwell_velocity.hpp>
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
typedef std::vector< size_t > cell_type1;
typedef std::vector< point< size_t , dim > > index_type1;

typedef thrust::host_vector< mass_point_type > mass_point_vector_type2;
typedef thrust::host_vector< size_t > cell_type2;
typedef thrust::host_vector< point< size_t , dim > > index_type2;

typedef thrust::device_vector< mass_point_type > mass_point_vector_type3;
typedef thrust::device_vector< size_t > cell_type3;
typedef thrust::device_vector< point< size_t , dim > > index_type3;


const size_t n_solvent = 32;
const point_type range( 16.0 , 16.0 );
const point_type null_point( 0.0 , 0.0 );
const point_type spacing( 1.0 , 1.0 );
const value_type sqrt_kbT_mass = 1.0;



int main( int argc , char **argv )
{
    srand48( 325345777 );
    drand48_generator rng;


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
	mass_point_vector_type3 solvent3 = solvent1;
	cell_type1 cell1( n_solvent );
	cell_type2 cell2( n_solvent );
	cell_type3 cell3( n_solvent );

	index_type1 indices1( n_solvent );
	index_type2 indices2( n_solvent );
	index_type3 indices3( n_solvent );



	typedef mpc2::regular_indexer< point_type ,
			with_start_point_tag ,
			check_and_throw_index_policy > regular_indexer_type;
	regular_indexer_type indexer( range , spacing );

	indexer.automatic_grid_shift( drand48 );

	std::transform( solvent1.begin() , solvent1.end() , cell1.begin() , make_indexer_caller( indexer ) );
	std::transform( solvent2.begin() , solvent2.end() , cell2.begin() , make_indexer_caller( indexer ) );
	thrust::transform( solvent3.begin() , solvent3.end() , cell3.begin() , make_indexer_caller( indexer ) );

	std::transform( solvent1.begin() , solvent1.end() , indices1.begin() ,  make_indexer_caller2( indexer ) );
	std::transform( solvent2.begin() , solvent2.end() , indices2.begin() , make_indexer_caller2( indexer ) );
	thrust::transform( solvent3.begin() , solvent3.end() , indices3.begin() , make_indexer_caller2( indexer ) );

	clog << indexer.start() << endl;


	cell_type2 cell3_host = cell3;
	index_type2 indices3_host = indices3;

	for( size_t i=0 ; i<cell1.size() ; ++i )
	{
		cout << solvent1[i].coor << "\t";
		cout << cell1[i] << "\t" << cell2[i] << "\t" << cell3_host[i] << "\t";
		cout << indices1[i] << "\t" << indices2[i] << "\t" << indices3_host[i] << "\n";
	}




	return 0;
}

