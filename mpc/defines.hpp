/*
 * mpc_defines.hpp
 *
 *  Created on: Dec 8, 2010
 *      Author: karsten
 */

#ifndef MPC_DEFINES_HPP_
#define MPC_DEFINES_HPP_

#ifdef __CUDACC__
#define FUNC_DECL __host__ __device__
#else
#define FUNC_DECl
#endif

#define MPC_NAMESPACE_BEGIN namespace mpc2 {
#define MPC_NAMESPACE_END }

#endif /* MPC_DEFINES_HPP_ */
