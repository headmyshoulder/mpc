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
#include <thrust/sort.h>
#include <thrust/binary_search.h>

// algebra includes
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>

#include <mpc/defines.hpp>
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
typedef std::vector< point_type > point_vector_type1;
typedef std::vector< size_t > cell_type1;
typedef std::vector< point< size_t , dim > > index_type1;

typedef thrust::host_vector< mass_point_type > mass_point_vector_type2;
typedef thrust::host_vector< point_type > point_vector_type2;
typedef thrust::host_vector< size_t > cell_type2;
typedef thrust::host_vector< point< size_t , dim > > index_type2;

typedef thrust::device_vector< mass_point_type > mass_point_vector_type3;
typedef thrust::device_vector< point_type > point_vector_type3;
typedef thrust::device_vector< size_t > cell_type3;
typedef thrust::device_vector< point< size_t , dim > > index_type3;


const size_t n_solvent = 32;
const point_type range( 16.0 , 16.0 );
const point_type null_point( 0.0 , 0.0 );
const point_type spacing( 8.0 , 8.0 );
const value_type sqrt_kbT_mass = 1.0;



template< class MassPointVector >
struct calculate_com
{
	typedef typename MassPointVector::value_type mass_point_type;
	typedef typename mass_point_type::point_type point_type;
	typedef typename point_traits< point_type >::value_type value_type;

	typedef typename MassPointVector::const_pointer const_pointer;

//	const MassPointVector &m_mp;
	const_pointer m_mp;

	calculate_com( const MassPointVector &mp ) : m_mp( mp.data() ) { }

	template< class Index >
	FUNC_DECL point_type operator()( Index begin , Index end ) const
	{
		point_type com;
		Index size = end - begin;
		for( Index i=begin ; i<end ; ++i )
		{
			mass_point_type cur = m_mp[i];
			com += cur.coor;
		}
		if( size ) com /= value_type( size );
		return com;
	}
};

template< class MassPointVector , class PointVector , class IndexVector , class Indexer >
void com( MassPointVector &mp , PointVector &com , IndexVector &begins , IndexVector &ends , const Indexer &indexer )
{
	size_t num_of_particles = mp.size() , num_of_cells = indexer.number_of_cells();
	IndexVector indices( num_of_particles );
	thrust::counting_iterator< size_t > search_begin( 0 );

	thrust::transform( mp.begin() , mp.end() , indices.begin() , make_indexer_caller( indexer ) );
	thrust::sort_by_key( indices.begin() , indices.end() , mp.begin() );
	thrust::lower_bound( indices.begin() , indices.end() , search_begin , search_begin + num_of_cells , begins.begin() );
	thrust::upper_bound( indices.begin() , indices.end() , search_begin , search_begin + num_of_cells , ends.begin() );
	thrust::transform( begins.begin() , begins.end() , ends.begin() , com.begin() , calculate_com< MassPointVector >( mp ) );
}





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
	mass_point_vector_type3 solvent3 = solvent1;
	cell_type1 cell1( n_solvent );
	cell_type2 cell2( n_solvent );
	cell_type3 cell3( n_solvent );

	cell_type1 begins1( indexer.number_of_cells() ) , ends1( indexer.number_of_cells() );
	cell_type2 begins2( indexer.number_of_cells() ) , ends2( indexer.number_of_cells() );
	cell_type3 begins3( indexer.number_of_cells() ) , ends3( indexer.number_of_cells() );
	point_vector_type1 com1( indexer.number_of_cells() );
	point_vector_type2 com2( indexer.number_of_cells() );
	point_vector_type3 com3( indexer.number_of_cells() );

	com( solvent1 , com1 , begins1 , ends1 , indexer );
	com( solvent2 , com2 , begins2 , ends2 , indexer );
	com( solvent3 , com3 , begins3 , ends3 , indexer );

	point_vector_type2 com3_host = com3;

	for( size_t i=0 ; i<indexer.number_of_cells() ; ++i )
	{
		cout << i << "\t";
		cout << begins1[i] << "\t" << ends1[i] << "\t" << com1[i] << "\t";
		cout << begins2[i] << "\t" << ends2[i] << "\t" << com2[i] << "\t";
		cout << begins3[i] << "\t" << ends3[i] << "\t" << com3_host[i] << "\n";
	}


//	for( size_t i=0 ; i<n_solvent ; ++i )
//	{
//		cout << i << "\t" << cell1[i] << "\t" << solvent1[i].coor << "\n";
//	}
//	cout << "\n\n";
//	for( size_t i=0 ; i<indexer.number_of_cells() ; ++i )
//	{
//		cout << i << "\t" << bucket_begin[i] << "\t" << bucket_end[i] << "\n";
//	}




//	cell_type2 cell3_host = cell3;
//	index_type2 indices3_host = indices3;
//
//	for( size_t i=0 ; i<cell1.size() ; ++i )
//	{
//		cout << solvent1[i].coor << "\t";
//		cout << cell1[i] << "\t" << cell2[i] << "\t" << cell3_host[i] << "\n";
//	}




	return 0;
}

