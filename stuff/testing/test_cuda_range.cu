
#include <iostream>
#include <vector>
#include <utility>

// thrust includes
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/random.h>

#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>

using namespace std;
using boost::numeric::odeint::thrust_algebra;

struct initialize
{
	double m_a , m_b;
	initialize( double a , double b ) : m_a( a ) , m_b( b ) { }

	template< class Tuple >
	__device__ __host__
	void operator()( Tuple t ) const
	{
		thrust::get< 0 >( t ) = m_a;
		thrust::get< 1 >( t ) = m_b;
	}
};


int main( int argc , char **argv )
{
	typedef thrust::device_vector< double > vector_type;
	typedef vector_type::iterator iterator_type;
	typedef std::pair< iterator_type , iterator_type > range_type1;
	const size_t n = 32;
	vector_type vec1( n ) , vec2( n );

	thrust_algebra::for_each2( vec1 , vec2 , initialize( 1.0 , 2.0 ) );

	range_type1 r1( vec1.begin() , vec1.begin() + n / 2 );
	range_type1 r2( vec2.begin() , vec2.begin() + n / 2 );
	thrust_algebra::for_each2( r1 , r2  , initialize( 4.0 , 8.0 ) );
	thrust_algebra::for_each2(
		std::make_pair( vec1.begin() , vec1.begin() + 4 ) ,
		std::make_pair( vec2.begin() , vec2.begin() + 4 ) ,
		initialize( 16.0 , 32.0 ) );

	for( size_t i=0 ; i<n ; ++i )
		cout << vec1[i] << "\t" << vec2[i] << endl;


	return 0;
}

