#include <iostream>
#include <fstream>
#include <vector>

// thrust includes
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

// includes from boost
#include <boost/random.hpp>

#include <mpc/point.hpp>
#include <mpc/mass_point.hpp>
#include <mpc/indexer.hpp>
#include <mpc/maxwell_velocity.hpp>
#include <mpc/analysis.hpp>

using namespace std;
using namespace mpc2;

const static size_t dim = 2;
typedef double value_type;
typedef point< value_type , dim > point_type;
typedef mass_point< point_type > mass_point_type;

typedef thrust::device_vector< point_type > point_vector_type;
typedef thrust::device_vector< mass_point_type > mass_point_vector_type;
typedef thrust::device_vector< size_t > cell_type;
typedef thrust::device_vector< point< size_t , dim > > index_type;

const size_t n_solvent = 256 * 256 *16;
const point_type range( 10.0 , 10.0 );
const value_type sqrt_kbT_mass = 1.0;
const size_t num_of_iterations = 1024;
const size_t num_of_bins = 100;


int main( int argc , char **argv )
{
	srand48( 325345777 );
	typedef boost::mt19937 rng_type;
	rng_type rng;

	typedef velocity_dist_y_device< point_vector_type , cell_type , mass_point_vector_type > vel_dist_type;

	vel_dist_type vel_dist( num_of_bins , range[1] , 0.0 , range[0] - 0.0 );

	// initialize the solvent particles
	mass_point_vector_type solvent( n_solvent );
	for( size_t i=0 ; i<n_solvent ; ++i )
	{
		mass_point_type mp;
		value_type y = drand48() * range[1];
		mp.coor[0] = drand48() * range[0];
		mp.coor[1] = y;
		maxwell_velocity( mp.vel , sqrt_kbT_mass , rng );
		mp.vel[0] += 4.0 * y * ( range[1] - y ) / range[1] / range[1] ;
		solvent[i] = mp;
	}

	for( size_t i=0 ; i<num_of_iterations ; ++i )
	{
		vel_dist.add( solvent );
	}


	std::vector< point_type > dist;
	vel_dist.get_result( dist );

	ofstream fout( "dat/result1.dat");
	for( size_t i=0 ; i<dist.size() - 1 ; ++i )
		fout << vel_dist.bin_mean( i ) << " " << dist[i] << endl;

	return 0;
}

