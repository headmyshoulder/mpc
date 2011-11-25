/*
 * fixed_boundary_corrector.hpp
 *
 *  Created on: Oct 6, 2010
 *      Author: karsten
 */

#ifndef FIXED_BOUNDARY_CORRECTOR_HPP_
#define FIXED_BOUNDARY_CORRECTOR_HPP_



#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/point_traits.hpp>
#include <mpc/mass_point.hpp>

MPC_NAMESPACE_BEGIN

template< class Point >
class fixed_boundary_corrector
{
public:

	typedef point_traits< Point > point_traits;
	const static size_t dim = point_traits::dim;
	typedef typename point_traits::value_type value_type;
	typedef typename point_traits::point_type point_type;
	typedef mass_point< point_type > mass_point_type;

	fixed_boundary_corrector( const point_type &min_bound , const point_type &max_bound )
	: m_min_bound( min_bound ) , m_max_bound( max_bound ) { }

	FUNC_DECL void operator()( mass_point_type &m ) const
	{
		const value_type eps = 0.001;
		for( size_t i=0 ; i<dim ; ++i )
		{
			if( m.coor[i] > m_max_bound[i] ) m.coor[i] = m_max_bound[i] - eps;
			if( m.coor[i] < m_min_bound[i] ) m.coor[i] = m_min_bound[i] + eps;
		}
	}

	FUNC_DECL void operator()( mass_point_type &m , value_type t ) const
	{
		this->operator()( m );
	}


	const point_type& min_bound( void ) const { return min_bound; }
	const point_type& max_bound( void ) const { return max_bound; }
	point_type& min_bound( void ) { return min_bound; }
	point_type& max_bound( void ) { return max_bound; }
	void set_min_bound( const point_type &min_bound ) { m_min_bound = min_bound; }
	void set_max_bound( const point_type &max_bound ) { m_max_bound = max_bound; }


private:

	point_type m_min_bound , m_max_bound;
};


MPC_NAMESPACE_END

#endif /* FIXED_BOUNDARY_CORRECTOR_HPP_ */
