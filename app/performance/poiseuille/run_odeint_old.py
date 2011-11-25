#! /usr/bin/python

import os
import time

from helpers import *

dim_param_list  = [ [   30 , 15 , 5 ] ]
dim_param_list += [ [   60 , 15 , 5 ] ]
dim_param_list += [ [  120 , 15 , 5 ] ]
dim_param_list += [ [  240 , 15 , 5 ] ]
dim_param_list += [ [  480 , 15 , 5 ] ]
dim_param_list += [ [  960 , 15 , 5 ] ]

executable_name = "/home/karsten/src/mpc/mpc/app/poiseuille_flow_2d/gcc-4.4/release/poiseuille_flow_2d_anderson_mpc_at_plus_a_performance"
dir = "res_old_home"
create_dir_structure_old( dir )

result = list()


for p in dim_param_list :
    
    print "Range       : " + str( p[0] ) + " " + str( p[1] )
    print "Density     : " + str( p[2] )
    
    outfile_name = get_outfile_name_old( dir , p )
    
    print "out file : " + outfile_name
    
    executable_file = executable_name + " " + str( p[0] ) + " " + str( p[1] ) + " " + str( p[2] ) + " " + outfile_name
    
    print "Running now " + executable_file + "!" 
    t_start = time.time()
    os.system( executable_file )
    t_end = time.time()
    print "Ready!"
    
    dt = t_end - t_start
    print "Execution time : " + str( dt )
    print "\n\n\n\n"

    result += [ [ num_of_particles( p ) , p[0] , p[1] , p[2] , dt ] ]

    f = open( dir + "/result.dat" , "w" )
    for res in result :
        line = ""
        for r in res :
            line += str( r ) + " "
        line += "\n"
        f.write( line )
