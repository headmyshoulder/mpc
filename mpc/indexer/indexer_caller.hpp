/*
 * indexer_caller.hpp
 *
 *  Created on: Jan 13, 2011
 *      Author: karsten
 */

#ifndef INDEXER_CALLER_HPP_
#define INDEXER_CALLER_HPP_

#include <mpc/defines.hpp>

MPC_NAMESPACE_BEGIN

template< class Indexer >
struct indexer_caller
{
	typedef typename Indexer::size_type size_type;
	Indexer m_indexer;
	indexer_caller( const Indexer &indexer = Indexer() ) : m_indexer( indexer ) { }


	template< class MassPoint >
	FUNC_DECL size_type operator()( const MassPoint &m ) const
	{
		return m_indexer.calc_index_from_point( m.coor );
	}
};

template< class Indexer >
indexer_caller< Indexer > make_indexer_caller( const Indexer &indexer )
{
	return indexer_caller< Indexer >( indexer );
}

template< class Indexer >
struct indexer_caller2
{
	typedef typename Indexer::index_type index_type;
	Indexer m_indexer;
	indexer_caller2( const Indexer &indexer = Indexer() ) : m_indexer( indexer ) { }

	template< class MassPoint >
	FUNC_DECL index_type operator()( const MassPoint &m ) const
	{
		return m_indexer.calc_index_tuple_from_point( m.coor );
	}
};

template< class Indexer >
indexer_caller2< Indexer > make_indexer_caller2( const Indexer &indexer )
{
	return indexer_caller2< Indexer >( indexer );
}

MPC_NAMESPACE_END

#endif /* INDEXER_CALLER_HPP_ */
