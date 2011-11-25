/*
 * point.hpp
 *
 *  Created on: May 16, 2010
 *      Author: karsten
 */


#ifndef MPC_COMMON_POINT_HPP_INCLUDED
#define MPC_COMMON_POINT_HPP_INCLUDED

#include <ostream>

#include <boost/static_assert.hpp>

#include <mpc/defines.hpp>


MPC_NAMESPACE_BEGIN

namespace detail
{
	template< class T , size_t Dim > struct initialize;

	template< class T > struct initialize< T , 1 >
	{
		FUNC_DECL void operator()( T val[1] , T x , T y , T z )
		{
			val[0] = x;
		}
	};

	template< class T > struct initialize< T , 2 >
	{
		FUNC_DECL void operator()( T val[2] , T x , T y , T z )
		{
			val[0] = x; val[1] = y;
		}
	};

	template< class T > struct initialize< T , 3 >
	{
		FUNC_DECL void operator()( T val[3] , T x , T y , T z )
		{
			val[0] = x; val[1] = y; val[2] = z;
		}
	};
}



//
// the point type
//
template< class T , size_t Dim >
class point
{
public:

	const static size_t dim = Dim;
	typedef T value_type;
	typedef point< value_type , dim > point_type;

	FUNC_DECL inline point( void )
	{
		for( size_t i=0 ; i<dim ; ++i ) m_val[i] = 0.0;
	}

	FUNC_DECL inline point( value_type val )
	{
		for( size_t i=0 ; i<dim ; ++i ) m_val[i] = val;
	}

	FUNC_DECL inline point( value_type x , value_type y , value_type z = 0.0 )
	{
		detail::initialize< value_type , dim >()( m_val , x , y , z );
	}

	template< size_t i > FUNC_DECL inline T get( void ) const
	{
		BOOST_STATIC_ASSERT( i < dim );
		return m_val[i];
	}
	template< size_t i > FUNC_DECL inline T& get( void )
	{
		BOOST_STATIC_ASSERT( i < dim );
		return m_val[i];
	}

	FUNC_DECL inline T operator[]( size_t i ) const { return m_val[i]; }
	FUNC_DECL inline T& operator[]( size_t i ) { return m_val[i]; }

	FUNC_DECL inline void assign( value_type val )
	{
		for( size_t i=0 ; i<dim ; ++i ) m_val[i] = val;
	}



	FUNC_DECL inline point_type& operator+=( const point_type& p )
	{
		for( size_t i=0 ; i<dim ; ++i )
			m_val[i] += p[i];
		return *this;
	}

	FUNC_DECL inline point_type& operator-=( const point_type& p )
	{
		for( size_t i=0 ; i<dim ; ++i )
			m_val[i] -= p[i];
		return *this;
	}

	FUNC_DECL inline point_type& operator+=( const value_type& val )
	{
		for( size_t i=0 ; i<dim ; ++i )
			m_val[i] += val;
		return *this;
	}

	FUNC_DECL inline point_type& operator-=( const value_type& val )
	{
		for( size_t i=0 ; i<dim ; ++i )
			m_val[i] -= val;
		return *this;
	}

	FUNC_DECL inline point_type operator*=( const point_type &p )
	{
		for( size_t i=0 ; i<dim ; ++i )
			m_val[i] *= p[i];
		return *this;
	}

	FUNC_DECL inline point_type operator/=( const point_type &p )
	{
		for( size_t i=0 ; i<dim ; ++i )
			m_val[i] /= p[i];
		return *this;
	}

	FUNC_DECL inline point_type& operator*=( const value_type &val )
	{
		for( size_t i=0 ; i<dim ; ++i )
			m_val[i] *= val;
		return *this;
	}

	FUNC_DECL inline point_type& operator/=( const value_type &val )
	{
		for( size_t i=0 ; i<dim ; ++i )
			m_val[i] /= val;
		return *this;
	}

private:

	T m_val[dim];
};

//
// the  p1 + p2 operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator+( const point< T , Dim > &p1 , const point< T , Dim > &p2 )
{
	point< T , Dim > tmp( p1 );
	tmp += p2;
	return tmp;
}

//
// the p1 + val operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator+( const point< T , Dim > &p1 , const T &val )
{
	point< T , Dim > tmp( p1 );
	tmp += val;
	return tmp;
}

//
// the val + p1 operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator+( const T &val , const point< T , Dim > &p1 )
{
	point< T , Dim > tmp( p1 );
	tmp += val;
	return tmp;
}

//
// the  p1 - p2 operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator-( const point< T , Dim > &p1 , const point< T , Dim > &p2 )
{
	point< T , Dim > tmp( p1 );
	tmp -= p2;
	return tmp;
}

//
// the p1 - val operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator-( const point< T , Dim > &p1 , const T &val )
{
	point< T , Dim > tmp( p1 );
	tmp -= val;
	return tmp;
}

//
// the val - p1 operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator-( const T &val , const point< T , Dim > &p1 )
{
	point< T , Dim > tmp( val );
	tmp -= p1;
	return tmp;
}



//
// the p1 * val operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator*( const point< T , Dim > &p1 , const T &val )
{
	point< T , Dim > tmp( p1 );
	tmp *= val;
	return tmp;
}

//
// the val * p1 operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator*( const T &val , const point< T , Dim > &p1 )
{
	point< T , Dim > tmp( p1 );
	tmp *= val;
	return tmp;
}



//
// the p1 / val operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator/( const point< T , Dim > &p1 , const T &val )
{
	point< T , Dim > tmp( p1 );
	tmp /= val;
	return tmp;
}







//
// the + operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator+( const point< T , Dim > &p )
{
	point< T , Dim > tmp;
	for( size_t i=0 ; i<Dim ; ++i ) tmp[i] = p[i];
	return tmp;
}




//
// the - operator
//
template< class T , size_t Dim >
FUNC_DECL inline point< T , Dim > operator-( const point< T , Dim > &p )
{
	point< T , Dim > tmp;
	for( size_t i=0 ; i<Dim ; ++i ) tmp[i] = -p[i];
	return tmp;
}



//
// scalar product
//
template< class T , size_t Dim >
FUNC_DECL inline T scalar_prod( const point< T , Dim > &p1 , const point< T , Dim > &p2 )
{
	T tmp = 0.0;
	for( size_t i=0 ; i<Dim ; ++i ) tmp += p1[i] * p2[i];
	return tmp;
}



//
// norm
//
template< class T , size_t Dim >
FUNC_DECL inline T norm( const point< T , Dim > &p1 )
{
	return scalar_prod( p1 , p1 );
}




//
// absolute value
//
template< class T , size_t Dim >
FUNC_DECL inline T abs( const point< T , Dim > &p1 )
{
	return sqrt( norm( p1 ) );
}


//
// cross product for Dim = 2
template< class T >
FUNC_DECL inline T cross_product( const point< T , 2 > &p1 , const point< T , 2 > &p2 )
{
	return p1[0] * p2[1] - p1[1] * p2[0];
}





//
// output operator
//
template< class T , size_t Dim >
std::ostream& operator<<( std::ostream &out , const point< T , Dim > &p )
{
	if( Dim > 0 ) out << p[0];
	for( size_t i=1 ; i<Dim ; ++i ) out << " " << p[i];
	return out;
}




MPC_NAMESPACE_END


#endif // MPC_COMMON_POINT_HPP_INCLUDED
