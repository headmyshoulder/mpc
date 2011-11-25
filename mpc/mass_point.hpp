/*
 * mass_point.hpp
 *
 *  Created on: May 16, 2010
 *      Author: karsten
 */

#ifndef MPC_COMMON_MASS_POINT_HPP_
#define MPC_COMMON_MASS_POINT_HPP_

#include <mpc/point_traits.hpp>

#include <ostream>

MPC_NAMESPACE_BEGIN


//
// mass point type
//
template< class PointType >
class mass_point
{
public:

	typedef mass_point mass_point_type;
	typedef typename point_traits< PointType >::point_type point_type;
	typedef typename point_traits< PointType >::value_type value_type;
	static const size_t dim = point_traits< PointType >::dim;

	FUNC_DECL mass_point( void ) : coor() , vel() {}

	point_type coor;
	point_type vel;
};

template< class PointType >
std::ostream& operator<<( std::ostream& out , const mass_point< PointType > &mp )
{
	out << mp.coor << " " << mp.vel;
	return out;
}


MPC_NAMESPACE_END






#endif /* MPC_COMMON_MASS_POINT_HPP_ */
