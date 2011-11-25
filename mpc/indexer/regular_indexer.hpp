/*
 * regular_index_calculator.hpp
 *
 *  Created on: Oct 1, 2010
 *      Author: karsten
 */

#ifndef REGULAR_INDEX_CALCULATOR_NEW_HPP_
#define REGULAR_INDEX_CALCULATOR_NEW_HPP_


#include <boost/mpl/or.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>

#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/point_traits.hpp>
#include <mpc/mass_point.hpp>


MPC_NAMESPACE_BEGIN


/*
 * how to calculate the index
 */
struct start_point_tag {};
struct with_start_point_tag : start_point_tag {};
struct without_start_point_tag : start_point_tag {};


struct no_check_index_policy
{
	template< class Index >
	FUNC_DECL void check_index( const Index &index , const Index &dimensions ) const
	{
	}
};

struct check_and_throw_index_policy
{
	// ToDo : use point_traits
	template< class Index >
	FUNC_DECL void check_index( const Index &index , const Index &dimensions ) const
	{
		for( size_t i=0 ; i<Index::dim ; ++i )
		{
			#ifndef __CUDACC__
			bool bad = false;
			for( size_t i=0 ; i<dim ; ++i )
			{
				if( tmp[i] >= m_dimensions[i] ) bad = true;
			}
			if( bad )
			{
				throw std::domain_error( "check_and_throw_index_policy : index >= m_dimensions" );
			}
			#endif
		}
	}
};

struct no_periodic_boundaries_policy
{
	template< class Index >
	FUNC_DECL void treat_periodic_boundaries( const Index &range , Index &index )
	{
	}
};




/*
 * regular index base class
 * calculates the index according tho the CalcPolicy:
 *
 * with_start_point     : i = p.x / dim.x
 * without_start_point  : i = ( p.x - start.x ) / dim.x
 *
 * ToDo
 */
template<
	class Point ,
	class IndexCalcPolicy = with_start_point_tag ,
	class CheckIndexPolicy = no_check_index_policy ,
	class PeriodicBoundariesPolicy = no_periodic_boundaries_policy ,
	class SizeType = size_t
>
class regular_indexer
{

public:

	/*
	 * typedefs
	 */
	typedef Point point_type;
	const static size_t dim = point_traits< point_type >::dim;
	typedef typename point_traits< point_type >::value_type value_type;

	typedef SizeType size_type;
	typedef point< size_type , dim > index_type;

	typedef IndexCalcPolicy index_calc_policy;
	typedef CheckIndexPolicy index_check_policy_type;
	typedef PeriodicBoundariesPolicy periodic_boundaries_policy_type;



	/*
	 * checks if index calc policy has the right type
	 */
	BOOST_STATIC_ASSERT((
			boost::mpl::or_<
			typename boost::is_same< index_calc_policy , with_start_point_tag >::type ,
			typename boost::is_same< index_calc_policy , without_start_point_tag >::type
			>::type::value
	));



	/*
	 * constructors
	 */
	FUNC_DECL
	regular_indexer(
		const point_type &range ,
		const point_type &spacing ,
		const index_check_policy_type &index_check_policy = index_check_policy_type() ,
		const periodic_boundaries_policy_type &periodic_boundaries_policy = periodic_boundaries_policy_type() )
	: m_range() , m_spacing() , m_start( 0.0 ) ,
	  m_dimensions( 0 ) , m_number_of_cells( 0 ) ,
	  m_offsets() ,
	  m_index_check_policy( index_check_policy ) ,
	  m_periodic_boundaries_policy( periodic_boundaries_policy )
	{
		initialize( range , spacing );
	}



	/*
	 * modifiers
	 */
	FUNC_DECL void initialize( const point_type &range , const point_type &spacing )
	{
		m_range = range;
		m_spacing = spacing;
		for( size_t i=0 ; i<dim ; ++i ) m_dimensions[i] = size_type( m_range[i] / m_spacing[i] ) + 2;
		m_number_of_cells = 1;
		for( size_t i=0 ; i<dim ; ++i ) m_number_of_cells *= m_dimensions[i];

		initialize_offsets( m_dimensions , m_offsets , boost::mpl::size_t< dim >() );
	}

	FUNC_DECL void set_start( const point_type &p )
	{
		m_start = p;
	}

	template< class UniformRand >
	FUNC_DECL void automatic_grid_shift( UniformRand &rand )
	{
		for( size_t i=0 ; i<dim ; ++i ) m_start[i] = - rand() * m_spacing[i] ;
	}



	/*
	 * getters for the member variables
	 */
	point_type FUNC_DECL start( ) const { return m_start; }

	point_type FUNC_DECL spacing( void ) const { return m_spacing; }

	point_type FUNC_DECL range( void ) const { return m_range; }

	index_type FUNC_DECL dimensions( void ) const { return m_dimensions; }

	size_type FUNC_DECL number_of_cells( void ) const { return m_number_of_cells; }


	/*
	 * index calculation
	 */
	size_type FUNC_DECL calc_index_from_point( const point_type &p ) const
	{
	  index_type ind = calc_index_tuple_from_point( p );
	  size_type s = calc_index_from_index_tuple( ind );
	  return s;
	}

	index_type FUNC_DECL calc_index_tuple_from_point( const point_type &p ) const
	{
		index_type tmp = calc_index_tuple_from_point_impl( p , index_calc_policy() );
		m_index_check_policy.check_index( tmp , m_dimensions );
		return tmp;
	}

	size_type FUNC_DECL calc_index_from_index_tuple( const index_type &ind ) const
	{
		size_type index = 0;
		for( size_t i=0 ; i<dim-1 ; ++i ) index += m_offsets[i] * ind[i];
		return index + ind[dim-1];
	}

	index_type FUNC_DECL calc_index_tuple_from_index( size_type ind ) const
	{
		return index_type( ind / m_dimensions[1] , ind % m_dimensions[1] );
	}



	/*
	 * minimal and maximal indices
	 */
	index_type FUNC_DECL min_index( void ) const
	{
		index_type tmp;
		for( size_t i=0 ; i<dim ; ++i ) tmp[i] = min_index_impl( i , index_calc_policy() );
		return tmp;
	}

	index_type FUNC_DECL max_index( void ) const
	{
		index_type tmp;
		for( size_t i=0 ; i<dim ; ++i ) tmp[i] = max_index_impl( i , index_calc_policy() );
		return tmp;
	}



	point_type FUNC_DECL low_left( const index_type &ind ) const
	{
		return calc_low_left_from_index_impl( ind , index_calc_policy() );
	}

	point_type FUNC_DECL up_right( const index_type &ind ) const
	{
		return low_left( ind ) + m_spacing;
	}


protected:

	FUNC_DECL
	void initialize_offsets( const index_type &dimensions , size_type offsets[dim-1] , boost::mpl::size_t< 1 > )
	{
	}

	FUNC_DECL
	void initialize_offsets( const index_type &dimensions , size_type offsets[dim-1] , boost::mpl::size_t< 2 > )
	{
		offsets[0] = dimensions[1];
	}


	FUNC_DECL
	void initialize_offsets( const index_type &dimensions , size_type offsets[dim-1] , boost::mpl::size_t< 3 > )
	{
		offsets[0] = dimensions[1] * dimensions[2];
		offsets[1] = dimensions[2];
	}



	/*
	 * implementation of the indices according to the index_calc_policy
	 */
	index_type FUNC_DECL calc_index_tuple_from_point_impl( const point_type &p , without_start_point_tag ) const
	{
		index_type tmp;
		for( size_t i=0 ; i<dim ; ++i ) tmp[i] = size_type ( p[i] / m_spacing[i] );
		return tmp;
	}

	index_type FUNC_DECL calc_index_tuple_from_point_impl( const point_type &p , with_start_point_tag ) const
	{
		index_type tmp;
		for( size_t i=0 ; i<dim ; ++i ) tmp[i] = size_type ( ( p[i] - m_start[i] ) / m_spacing[i] );
		return tmp;
	}

	point_type FUNC_DECL calc_low_left_from_index_impl( const index_type &ind , without_start_point_tag ) const
	{
		point_type p;
		for( size_t i=0 ; i<dim ; ++i ) p[i] = value_type( ind[i] ) * m_spacing[i];
		return p;
	}

	point_type FUNC_DECL calc_low_left_from_index_impl( const index_type &ind , with_start_point_tag ) const
	{
		point_type p;
		for( size_t i=0 ; i<dim ; ++i ) p[i] = m_start[i] + value_type( ind[i] ) * m_spacing[i];
		return p;
	}


	size_type FUNC_DECL min_index_impl( const size_t i , start_point_tag ) const
	{
		return 0;
	}

	size_type FUNC_DECL max_index_impl( const size_t i , without_start_point_tag ) const
	{
		return size_type( m_range[i] / m_spacing[i] );
	}

	size_type FUNC_DECL max_index_impl( const size_t i , with_start_point_tag ) const
	{
		return size_type( ( m_range[i] - m_start[i] ) / m_spacing[i] );
	}


	/*
	 * member variables
	 */
	point_type m_range , m_spacing , m_start;
	index_type m_dimensions;
	size_type m_number_of_cells;
	size_type m_offsets[dim-1];

	index_check_policy_type m_index_check_policy;
	periodic_boundaries_policy_type m_periodic_boundaries_policy;
};

MPC_NAMESPACE_END


#endif /* REGULAR_INDEX_CALCULATOR_NEW_HPP_ */
