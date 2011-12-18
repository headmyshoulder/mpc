/*
 * drand48_generator.hpp
 *
 *  Created on: Dec 17, 2011
 *      Author: karsten
 */

#ifndef DRAND48_GENERATOR_HPP_
#define DRAND48_GENERATOR_HPP_


#include <mpc/defines.hpp>

MPC_NAMESPACE_BEGIN


/*
 * Generator adaptor for drand48. drand48 can now be used with the random distributions
 */
struct drand48_generator
{
    typedef double result_type;
    result_type operator()( void ) const { return drand48(); }
    result_type min( void ) const { return 0.0; }
    result_type max( void ) const { return 1.0; }
};


MPC_NAMESPACE_END

#endif /* DRAND48_GENERATOR_HPP_ */
