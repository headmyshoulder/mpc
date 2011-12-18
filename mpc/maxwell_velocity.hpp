/*
 * maxwell_velocity.hpp
 *
 *  Created on: Jan 13, 2011
 *      Author: karsten
 */

#ifndef MAXWELL_VELOCITY_HPP_
#define MAXWELL_VELOCITY_HPP_

#include <curand_kernel.h>

#include <mpc/defines.hpp>
#include <mpc/point.hpp>
#include <mpc/point_traits.hpp>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

MPC_NAMESPACE_BEGIN


/*
 * assign maxwell velocity to the components of the point
 */
template< class Point , class Random >
void maxwell_velocity( Point &p , typename point_traits< Point >::value_type sqrt_kbT_mass , Random &random )
{
	typedef point_traits< Point > point_traits;
	typedef typename point_traits::value_type value_type;
	const static size_t dim = point_traits::dim;
	typedef boost::normal_distribution< value_type > gauss_type;
	typedef boost::variate_generator< Random& , gauss_type > generator_type;

	gauss_type gauss( 0.0 , sqrt_kbT_mass );
	generator_type gen( random , gauss );
	for( size_t i=0 ; i<dim ; ++i ) p[i] = gen();
}


template< class Point >
FUNC_DECL
void maxwell_velocity_curand( Point &p , typename point_traits< Point >::value_type sqrt_kbT_mass , curandState &random )
{
	const size_t dim = point_traits< Point >::dim;
	for( size_t i=0 ; i<dim ; ++i )
		p[i] = curand_normal( &random ) * sqrt_kbT_mass;
}



MPC_NAMESPACE_END

#endif /* MAXWELL_VELOCITY_HPP_ */
