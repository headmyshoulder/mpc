#ifndef SRD_SERVICES_ANALYSIS_VELOCITY_DIST_HPP_INCLUDED
#define SRD_SERVICES_ANALYSIS_VELOCITY_DIST_HPP_INCLUDED

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
template< class MassPoint >
class velocity_dist_y
{

public:

	typedef MassPoint mass_point_type ;
	typedef typename mass_point_type::point_type point_type ;
	typedef typename mass_point_type::value_type value_type ;


private:

	size_t m_num_of_bins;
	value_type m_height;    
	value_type m_x_min , m_x_max;
	std::vector< point_type > m_bins;
	std::vector< size_t > m_counts;
	std::vector< value_type > m_bin_min , m_bin_mean , m_bin_max;

	size_t num_of_counts( void ) const
	{
		return std::accumulate( m_counts.begin() ,m_counts.end() , 0 );
	}

public:

	velocity_dist_y
	(
			size_t num_of_bins ,
			value_type height ,
			value_type x_min , value_type x_max
	)
	: m_num_of_bins( num_of_bins ) , m_height( height ) ,
	  m_x_min( x_min ) , m_x_max( x_max ) ,
	  m_bins( num_of_bins , point_type( 0.0 ) ) ,
	  m_counts( num_of_bins , 0 ) ,
	  m_bin_min( num_of_bins , 0.0 ) , m_bin_mean( num_of_bins , 0.0 ) ,
	  m_bin_max( num_of_bins , 0.0 )
	{
		value_type factor = m_height / value_type( m_num_of_bins );
		for( size_t i=0 ; i<m_num_of_bins ; ++i )
		{
			m_bin_min[i] = value_type(i) * factor;
			m_bin_mean[i] = ( value_type(i) + 0.5 ) * factor;
			m_bin_max[i] = value_type(i+1) * factor;
		}
	}

	template< class MassPointContainer >
	void add( const MassPointContainer &mass_points )
	{
		thrust::host_vector< mass_point_type > tmp_points = mass_points;

		for( size_t i=0 ; i<tmp_points.size() ; ++i )
		{
			const mass_point_type &m = tmp_points[i];
			const point_type &v = m.vel , &p = m.coor;

			if( ( m_x_min < p[0] ) && ( p[0] < m_x_max ) )
			{
				size_t bin_index = size_t( p[1] / m_height * value_type( m_num_of_bins ) );
				m_bins[bin_index] += v;
				m_counts[bin_index] ++;
			}
			if( p[1] >= m_height ) std::cerr << "Scheisse1!" << std::endl;
			if( p[1] <= 0.0 ) std::cerr << "Scheisse2!" << std::endl;
		}
	}

	value_type bin_min( size_t i ) const { return m_bin_min[i]; }
	value_type bin_mean( size_t i ) const { return m_bin_mean[i]; }
	value_type bin_max( size_t i ) const { return m_bin_max[i]; }


	point_type bin_value( size_t i ) const
	{
		size_t count = m_counts[i];
		return ( count ) ? m_bins[i] / value_type( count ) : 0.0;
	}

	value_type bin_density( size_t i ) const
	{
		size_t num = num_of_counts();
		if( num ) return value_type( m_counts[i] ) / value_type( num );
		else return 0.0;
	}

	void write( std::ostream &out ) const
	{
		for( size_t i=0 ; i<m_num_of_bins ; ++i )
		{
			out << bin_min(i) << "\t" << bin_mean(i) << "\t";
			out << bin_max(i) << "\t";
			out << bin_value(i) << "\t" << bin_density(i) << "\t";
			out << m_bins[i] << "\t" << m_counts[i] << std::endl;
		}
	}

	size_t size( void ) const { return m_num_of_bins; }
	point_type operator[]( size_t i ) const { return bin_value( i ); }
};

MPC_NAMESPACE_END

#endif //SRD_SERVICES_ANALYSIS_VELOCITY_DIST_HPP_INCLUDED
