/*
 * anderson.hpp
 *
 *  Created on: Dec 23, 2010
 *      Author: karsten
 */

#ifndef ANDERSON_PLUS_A_HPP_
#define ANDERSON_PLUS_A_HPP_

#include <mpc/defines.hpp>
#include <mpc/point_traits.hpp>
#include <mpc/maxwell_velocity.hpp>

MPC_NAMESPACE_BEGIN

template< class MassPointVector >
struct anderson_plus_a
{
	typedef MassPointVector mass_point_vector;
	typedef typename mass_point_vector::value_type mass_point_type;
	typedef typename mass_point_type::point_type point_type;
	typedef typename point_traits< point_type >::value_type value_type;
	const static size_t dim = point_traits< point_type >::dim;

	typedef mass_point_type* mass_point_pointer;

	mass_point_pointer m_mp;
	value_type m_sqrt_kbT_mass;


	anderson_plus_a( mass_point_vector &mp , const value_type &sqrt_kbT_mass = value_type( 1.0 ) )
	: m_mp( thrust::raw_pointer_cast( mp.data() ) ) , m_sqrt_kbT_mass( sqrt_kbT_mass )
	{ }

	template< class Tuple >
	FUNC_DECL
	void operator()( Tuple t ) const
	{
		typedef typename thrust::tuple_element< 0 , Tuple >::type size_type_ref;
		typedef typename thrust::tuple_element< 1 , Tuple >::type size_type_ref2;
		typedef typename thrust::tuple_element< 2 , Tuple >::type rng_state_type;
		BOOST_STATIC_ASSERT(( boost::is_same< size_type_ref , size_type_ref2 >::value ));
		BOOST_STATIC_ASSERT(( boost::is_reference< size_type_ref >::value ));
		BOOST_STATIC_ASSERT(( boost::is_reference< rng_state_type >::value ));

		typedef typename boost::remove_reference< size_type_ref >::type size_type;


		size_type begin = thrust::get< 0 >( t ) , end = thrust::get< 1 >( t );
		size_type size = end - begin;
		rng_state_type &rng = thrust::get< 2 >( t );

		if( size < 2 ) return;

		value_type inv_n = 1.0 / value_type( size );

		point_type com_vel( 0.0 ) , com_coor( 0.0 );
		value_type angular_momentum = 0.0;
		for( size_type i=begin ; i<end ; ++i )
		{
			mass_point_type m = m_mp[i];
			com_vel += m.vel;
			com_coor += m.coor;
			angular_momentum += cross_product( m.coor , m.vel );
		}
		com_vel *= inv_n;
		com_coor *= inv_n;
		angular_momentum -= cross_product( com_coor , com_vel ) / inv_n;


		value_type inertia_tensor = 0.0;
		point_type new_com_vel( 0.0 );
		for( size_type i=begin ; i<end ; ++i )
		{
			mass_point_type m = m_mp[i];
			maxwell_velocity_curand( m.vel , m_sqrt_kbT_mass , rng );
			new_com_vel += m.vel;
			inertia_tensor += norm( m.coor - com_coor );
			m_mp[i] = m;
		}
		new_com_vel *= inv_n;


		value_type angular_momentum_new = 0.0;
		for( size_type i=begin ; i<end ; ++i )
		{
			mass_point_type m = m_mp[i];
			m.vel -= new_com_vel;
			angular_momentum_new += cross_product( m.coor , m.vel );
			m_mp[i] = m;
		}


		value_type angular_momentum_diff = angular_momentum - angular_momentum_new;
		value_type omega = angular_momentum_diff / inertia_tensor ;
		for( size_type i=begin ; i<end ; ++i )
		{
			mass_point_type m = m_mp[i];
			m.vel[0] -= omega * ( m.coor[1] - com_coor[1] );
			m.vel[1] += omega * ( m.coor[0] - com_coor[0] );
			m.vel += com_vel;
			m_mp[i] = m;
		}
	}
};

template< class MassPointVector , class IndexVector >
struct anderson_plus_a_permutating
{
	typedef MassPointVector mass_point_vector;
	typedef IndexVector index_vector;
	typedef typename mass_point_vector::value_type mass_point_type;
	typedef typename index_vector::value_type index_type;
	typedef typename mass_point_type::point_type point_type;
	typedef typename point_traits< point_type >::value_type value_type;
	const static size_t dim = point_traits< point_type >::dim;

	typedef mass_point_type* mass_point_pointer;
	typedef index_type* index_pointer;


	mass_point_pointer m_mp;
	index_pointer m_indices;
	value_type m_sqrt_kbT_mass;


	anderson_plus_a_permutating( mass_point_vector &mp , index_vector &indices , const value_type &sqrt_kbT_mass = value_type( 1.0 ) )
	: m_mp( thrust::raw_pointer_cast( mp.data() ) ) , m_indices( thrust::raw_pointer_cast( indices.data() ) ) , m_sqrt_kbT_mass( sqrt_kbT_mass )
	{ }

	template< class Tuple >
	FUNC_DECL
	void operator()( Tuple t ) const
	{
		typedef typename thrust::tuple_element< 0 , Tuple >::type size_type_ref;
		typedef typename thrust::tuple_element< 1 , Tuple >::type size_type_ref2;
		typedef typename thrust::tuple_element< 2 , Tuple >::type rng_state_type;
		BOOST_STATIC_ASSERT(( boost::is_same< size_type_ref , size_type_ref2 >::value ));
		BOOST_STATIC_ASSERT(( boost::is_reference< size_type_ref >::value ));
		BOOST_STATIC_ASSERT(( boost::is_reference< rng_state_type >::value ));

		typedef typename boost::remove_reference< size_type_ref >::type size_type;


		size_type begin = thrust::get< 0 >( t ) , end = thrust::get< 1 >( t );
		size_type size = end - begin;
		rng_state_type &rng = thrust::get< 2 >( t );

		if( size < 2 ) return;

		value_type inv_n = 1.0 / value_type( size );

		point_type com_vel( 0.0 ) , com_coor( 0.0 );
		value_type angular_momentum = 0.0;
		for( size_type i=begin ; i<end ; ++i )
		{
//			mass_point_type m = m_mp[ m_indices[i] ];
			mass_point_type m( m_mp[ m_indices[i] ] );
			com_vel += m.vel;
			com_coor += m.coor;
			angular_momentum += cross_product( m.coor , m.vel );
		}
		com_vel *= inv_n;
		com_coor *= inv_n;
		angular_momentum -= cross_product( com_coor , com_vel ) / inv_n;

		value_type inertia_tensor = 0.0;
		point_type new_com_vel( 0.0 );
		for( size_type i=begin ; i<end ; ++i )
		{
//			mass_point_type m = m_mp[ m_indices[i] ];
			mass_point_type m( m_mp[ m_indices[i] ] );
			maxwell_velocity_curand( m.vel , m_sqrt_kbT_mass , rng );
			new_com_vel += m.vel;
			inertia_tensor += norm( m.coor - com_coor );
			m_mp[ m_indices[i] ] = m;
		}
		new_com_vel *= inv_n;


		value_type angular_momentum_new = 0.0;
		for( size_type i=begin ; i<end ; ++i )
		{
//			mass_point_type m = m_mp[ m_indices[i] ];
			mass_point_type m( m_mp[ m_indices[i] ] );
			m.vel -= new_com_vel;
			angular_momentum_new += cross_product( m.coor , m.vel );
			m_mp[ m_indices[i] ] = m;
		}


		value_type angular_momentum_diff = angular_momentum - angular_momentum_new;
		value_type omega = angular_momentum_diff / inertia_tensor ;
		for( size_type i=begin ; i<end ; ++i )
		{
//			mass_point_type m = m_mp[ m_indices[i] ];
			mass_point_type m( m_mp[ m_indices[i] ] );
			m.vel[0] -= omega * ( m.coor[1] - com_coor[1] );
			m.vel[1] += omega * ( m.coor[0] - com_coor[0] );
			m.vel += com_vel;
			m_mp[ m_indices[i] ] = m;
		}
	}
};



template< class MassPointVector , class IndexVector , class RngStateVector , class Indexer >
void anderson( MassPointVector &mp , IndexVector &begins , IndexVector &ends , RngStateVector &rng , const Indexer &indexer )
{
	size_t num_of_particles = mp.size() , num_of_cells = indexer.number_of_cells();
	IndexVector indices( num_of_particles );
	thrust::counting_iterator< size_t > search_begin( 0 );

	thrust::transform( mp.begin() , mp.end() , indices.begin() , make_indexer_caller( indexer ) );
	thrust::sort_by_key( indices.begin() , indices.end() , mp.begin() );
	thrust::lower_bound( indices.begin() , indices.end() , search_begin , search_begin + num_of_cells , begins.begin() );
	thrust::upper_bound( indices.begin() , indices.end() , search_begin , search_begin + num_of_cells , ends.begin() );
	thrust::for_each(
			thrust::make_zip_iterator( thrust::make_tuple( begins.begin() , ends.begin() , rng.begin() ) ) ,
			thrust::make_zip_iterator( thrust::make_tuple( begins.end() , ends.end() , rng.end() ) ) ,
			anderson_plus_a< MassPointVector >( mp )
			);
}

MPC_NAMESPACE_END


#endif /* ANDERSON_PLUS_A_HPP_ */
