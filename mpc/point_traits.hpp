/*
 * point_traits.hpp
 *
 *  Created on: Dec 17, 2010
 *      Author: karsten
 */

#ifndef POINT_TRAITS_HPP_
#define POINT_TRAITS_HPP_

#include <cstddef>

#include <mpc/point.hpp>

MPC_NAMESPACE_BEGIN


//
// point traits
//
template< class Point >
struct point_traits;

//
// point traits for mpc::point< T , Dim >
//
template< class T , size_t Dim >
struct point_traits< point< T , Dim > >
{
	typedef point< T , Dim > point_type;
	static const size_t dim = point_type::dim;
	typedef T value_type;
};


MPC_NAMESPACE_END

#endif /* POINT_TRAITS_HPP_ */
