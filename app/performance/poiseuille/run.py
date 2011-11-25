#! /usr/bin/python

import os
import time

from build_env import *
from helpers import *


dim_param_list = list()
dim_param_list += [ [   30 , 15 , 5 ] ]
dim_param_list += [ [   60 , 15 , 5 ] ]
dim_param_list += [ [  120 , 15 , 5 ] ]
dim_param_list += [ [  240 , 15 , 5 ] ]
dim_param_list += [ [  480 , 15 , 5 ] ]
dim_param_list += [ [  960 , 15 , 5 ] ]
dim_param_list += [ [  1920 , 15 , 5 ] ]
#dim_param_list += [ [   30 , 15 , 15 ] ]
#dim_param_list += [ [   60 , 15 , 15 ] ]
#dim_param_list += [ [  120 , 15 , 15 ] ]
#dim_param_list += [ [  240 , 15 , 15 ] ]
#dim_param_list += [ [  480 , 15 , 15 ] ]
#dim_param_list += [ [  960 , 480 , 15 ] ]

run_param_list = list()
# run_param_list  += [ [ "std::vector" , False , 0 ] ]
run_param_list += [ [ "thrust::device_vector" , False , 0 ]]
# run_param_list += [ [ "thrust::device_vector" , True , 1  ] ]
# run_param_list += [ [ "thrust::device_vector" , True , 2 ] ]
# run_param_list += [ [ "thrust::device_vector" , True , 4 ] ]

param_list = get_param_list( dim_param_list , run_param_list )

template_name = "poiseuille_template.cu"
template_str = get_template_string( template_name )
dir = "res_landau_gpu"

create_dir_structure( dir )


result = list()


for p in param_list :
    range_params = p[0]
    run_params = p[1]
    add_cxx_flags = run_params[1]
    
    print "Range       : " + str( range_params[0] ) + " " + str( range_params[1] )
    print "Density     : " + str( range_params[2] )
    print "Vector type : " + run_params[0]
    if run_params[1] != False :
        print "Compiling with OpenMP"
    if run_params[2] != 0 :
        print "Threads : " + str( run_params[2] )
    
    
    cfile_str = template_str;
    
    outfile_name = get_outfile_name( dir , range_params , run_params )
    cfile_str = replace_range( cfile_str , range_params )
    cfile_str = replace_vector( cfile_str , run_params )
    cfile_str = replace_outfile( cfile_str , outfile_name )
    cfile_str = replace_omp_num_threads( cfile_str , run_params )
    
  
    cfile_name = get_cfile_name( dir , range_params , run_params )
    write_cfile( cfile_str , cfile_name )
    
    executable_file = get_executable_name( dir , range_params , run_params )  
    
    print "src file : " + cfile_name
    print "exe file : " + executable_file
    print "out file : " + outfile_name
    
    call_str = get_system_string( cfile_name , executable_file , add_cxx_flags  )
    print "Compiling now!"
    print call_str
    compile_res = os.system( call_str + " 2> tmp.dat" )
    print "Ready!"
    
    print "Running now " + executable_file + "!" 
    t_start = time.time()
    os.system( executable_file )
    t_end = time.time()
    print "Ready!"
    
    dt = t_end - t_start
    print "Execution time : " + str( dt )
    print "\n\n\n\n"

    result += [ [ num_of_particles( range_params ) , range_params[0] , range_params[1] , range_params[2] , run_params[0] , run_params[1] , run_params[2] , dt ] ]

    f = open( dir + "/result.dat" , "w" )
    for res in result :
        line = ""
        for r in res :
            line += str( r ) + " " + "\n"
        f.write( line )
