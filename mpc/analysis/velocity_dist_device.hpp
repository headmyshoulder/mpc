#ifndef SRD_SERVICES_ANALYSIS_VELOCITY_DIST_DEVICE_HPP_INCLUDED
#define SRD_SERVICES_ANALYSIS_VELOCITY_DIST_DEVICE_HPP_INCLUDED

#include <vector>
#include <ostream>
#include <numeric>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <thrust/host_vector.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/counting_iterator.h>

#include <mpc/defines.hpp>
#include <mpc/point_traits.hpp>


MPC_NAMESPACE_BEGIN


/*
 * ToDo: Generalize for arbitrary directions
 */
template< class PointVector , class CountVector , class MassPointVector >
class velocity_dist_y_device
{
public:

	typedef PointVector point_vector_type;
	typedef MassPointVector mass_point_vector_type;
	typedef CountVector count_vector_type;

	typedef typename point_vector_type::value_type point_type;
	typedef typename mass_point_vector_type::value_type mass_point_type;
	typedef typename count_vector_type::value_type size_type;
	typedef typename point_traits< point_type >::value_type value_type;
	const static size_t dim = point_traits< point_type >::dim;

	typedef typename boost::is_same< point_type , typename mass_point_type::point_type >::type xyz_type;
	BOOST_STATIC_ASSERT(( boost::is_same< point_type , typename mass_point_type::point_type >::type::value  ));

	velocity_dist_y_device( size_t num_of_bins , value_type height , value_type x_min , value_type x_max )
	: m_num_of_bins( num_of_bins ) , m_height( height ) , m_x_min( x_min ) , m_x_max( x_max ) ,
	  m_bins( num_of_bins + 1 ) , m_counts( num_of_bins + 1 ) ,
	  m_bin_min( num_of_bins ) , m_bin_mean( num_of_bins ) , m_bin_max( num_of_bins ) ,
	  m_indices() , m_permutated_indices() ,
	  m_keys( num_of_bins + 1 ) , m_av( num_of_bins + 1 ) ,
	  m_bucket_begin( num_of_bins + 1 ) , m_bucket_end( num_of_bins + 1 ) , m_bucket_sizes( num_of_bins + 1 )
	{
		thrust::fill( m_bins.begin() , m_bins.end() , point_type( 0.0 ) );
		thrust::fill( m_counts.begin() , m_counts.end() , size_type( 0.0 ) );
		value_type factor = m_height / value_type( m_num_of_bins );
		for( size_t i=0 ; i<m_num_of_bins ; ++i )
		{
			m_bin_min[i] = value_type(i) * factor;
			m_bin_mean[i] = ( value_type(i) + 0.5 ) * factor;
			m_bin_max[i] = value_type(i+1) * factor;
		}

	}

	void add( const mass_point_vector_type &mp )
	{
		if( m_indices.size() != mp.size() ) m_indices.resize( mp.size() );
		if( m_permutated_indices.size() != mp.size() ) m_permutated_indices.resize( mp.size() );

		thrust::counting_iterator<unsigned int> search_begin(0);

		// determine the index of each particle
		thrust::transform( mp.begin() , mp.end() , m_indices.begin() ,
				mp_to_indices( m_num_of_bins , m_height , m_x_min , m_x_max ) );

		// sort according to the indices
		thrust::copy( search_begin , search_begin + mp.size() , m_permutated_indices.begin() );
		thrust::sort_by_key( m_indices.begin() , m_indices.end() , m_permutated_indices.begin() );

		// calculate the center of mass velocity of each bins
		thrust::reduce_by_key( m_indices.begin() , m_indices.end() ,
				thrust::make_transform_iterator(
						thrust::make_permutation_iterator( mp.begin() , m_permutated_indices.begin() ) ,
						mp_to_vel() ) ,
				m_keys.begin() , m_av.begin() );

		// add the center of mass velocity to the bins
		thrust::transform( m_bins.begin() , m_bins.end() , m_av.begin() , m_bins.begin() ,
				thrust::plus< point_type >() );

		// find the beginning of each bucket's list of points
		thrust::lower_bound( m_indices.begin() , m_indices.end() ,
				search_begin , search_begin + m_num_of_bins + 1 ,
				m_bucket_begin.begin() );

		// find the end of each bucket's list of points
		thrust::upper_bound( m_indices.begin(), m_indices.end() ,
				search_begin , search_begin + m_num_of_bins + 1 ,
		        m_bucket_end.begin() );

		// take the difference between bounds to find each bucket’s size
		thrust::transform( m_bucket_end.begin() , m_bucket_end.end() ,
				m_bucket_begin.begin() , m_bucket_sizes.begin() ,
				thrust::minus<unsigned int>() );

		// add the bucket size to the counts
		thrust::transform( m_counts.begin() , m_counts.end() , m_bucket_sizes.begin() , m_counts.begin() ,
				thrust::plus< size_type >() );
	}


	void get_result( std::vector< point_type > &res )
	{
		thrust::host_vector< point_type > bins( m_bins );
		thrust::host_vector< size_type > counts( m_counts );

		assert( bins.size() == counts.size() );

		res.resize( bins.size() );
		for( size_t i=0 ; i<bins.size() ; ++i )
		{
			res[i] = ( ( counts[i] == 0 ) ? point_type( 0.0 , 0.0 ) : bins[i] / value_type( counts[i] ) );
		}
	}

	value_type bin_min( size_t i ) const { return m_bin_min[i]; }
	value_type bin_mean( size_t i ) const { return m_bin_mean[i]; }
	value_type bin_max( size_t i ) const { return m_bin_max[i]; }

	struct mp_to_vel : thrust::unary_function< mass_point_type , point_type >
	{
		FUNC_DECL
		point_type operator()( mass_point_type mp ) const { return mp.vel; }
	};

	struct mp_to_indices
	{
		size_t m_num_of_bins;
		value_type m_height , m_x_min , m_x_max;

		mp_to_indices( size_t num_of_bins , value_type height , value_type x_min , value_type x_max )
		: m_num_of_bins( num_of_bins ) , m_height( height ) , m_x_min( x_min ) , m_x_max( x_max ) { }

		FUNC_DECL
		size_type operator()( const mass_point_type &mp ) const
		{
			if( ( mp.coor[0] > m_x_min ) && ( mp.coor[0] < m_x_max ) )
			{
				return size_type( mp.coor[1] / m_height * value_type( m_num_of_bins ) );
			}
			else
			{
				return m_num_of_bins;
			}
		}
	};



private:

	size_t m_num_of_bins;
	value_type m_height;
	value_type m_x_min , m_x_max;
	point_vector_type m_bins;
	count_vector_type m_counts;
	std::vector< value_type > m_bin_min , m_bin_mean , m_bin_max;

	count_vector_type m_indices , m_permutated_indices;
	count_vector_type m_keys;
	point_vector_type m_av;
	count_vector_type m_bucket_begin , m_bucket_end , m_bucket_sizes;
};




MPC_NAMESPACE_END

#endif //SRD_SERVICES_ANALYSIS_VELOCITY_DIST_DEVICE_HPP_INCLUDED
