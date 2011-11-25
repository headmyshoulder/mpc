from build_env import *

import os

def num_of_particles( dim_param ):
    return dim_param[0] * dim_param[1] * dim_param[2] 


def estimate_memory( dim_param ):
    dim = 2
    bytes_per_float = 8
    var = 3
    return num_of_particles( dim_param ) * dim * bytes_per_float * var


def create_dir_structure( dir ):
    if dir[ len(dir) -1 ] == '/' :
        dir = dir[0:len(dir)-1]
    if not os.path.exists( dir ) :
        os.mkdir( dir )
    if not os.path.exists( dir + "/src" ) :
        os.mkdir( dir + "/src" )
    if not os.path.exists( dir + "/bin" ) :
        os.mkdir( dir + "/bin" )
    if not os.path.exists( dir + "/dat" ) :
        os.mkdir( dir + "/dat" )
        
def create_dir_structure_old( dir ) :
    if dir[ len(dir) -1 ] == '/' :
        dir = dir[0:len(dir)-1]
    if not os.path.exists( dir ) :
        os.mkdir( dir )
    if not os.path.exists( dir + "/dat" ) :
        os.mkdir( dir + "/dat" )
        
def get_file_name_base( range_params , run_params ) :
     
    val1 = '%(val)04d' % { "val" : range_params[0] }
    val2 = '%(val)04d' % { "val" : range_params[1] }
    val3 = '%(val)04d' % { "val" : range_params[2] }
    val4 = run_params[0].replace( "::" , "_" )
    fn = "poiseuille_" + val1 + "_" + val2 + "_" + val3 + "_" + val4
    if run_params[1] == True :
        fn += "_open_mp"
    if run_params[2] != 0 :
        fn += "_" + str( run_params[2] )
    return fn
        
def get_cfile_name( dir , range_params , run_params ) :
    if dir[ len(dir) -1 ] == '/' :
        dir = dir[0:len(dir)-1]
    return dir + "/src/" + get_file_name_base( range_params , run_params ) + ".cu"


def get_outfile_name( dir , range_params , run_params ) :
    if dir[ len(dir) -1 ] == '/' :
        dir = dir[0:len(dir)-1]
    return dir + "/dat/" + get_file_name_base( range_params , run_params ) + ".dat"


def get_outfile_name_old( dir , range_params ) :
    if dir[ len(dir) -1 ] == '/' :
        dir = dir[0:len(dir)-1]
    val1 = '%(val)04d' % { "val" : range_params[0] }
    val2 = '%(val)04d' % { "val" : range_params[1] }
    val3 = '%(val)04d' % { "val" : range_params[2] }
    base = "poiseuille_" + val1 + "_" + val2 + "_" + val3
    return dir + "/dat/" + base + ".dat"




def get_executable_name( dir , range_params , run_params ) :
    if dir[ len(dir) -1 ] == '/' :
        dir = dir[0:len(dir)-1]
    return dir + "/bin/" + get_file_name_base( range_params , run_params )


def get_param_list( dim_param_list , run_param_list ):
    param_list = list()
    for p1 in dim_param_list :
        for p2 in run_param_list :
            param_list += [[ p1 , p2 ]]
    return param_list
    

def get_template_string( template_name ):
    f = open( template_name , "r" )
    template_str = ""
    for line in f:
        template_str += line
    return template_str 


def replace_range( cfile_str , range_params ) :
    val1 = '%(val).1f' % { "val" : float( range_params[0] ) }
    val2 = '%(val).1f' % { "val" : float( range_params[1] ) }
    val3 = '%(val).1f' % { "val" : float( range_params[2] ) }
    range_str = "const point_type range( " + val1 + " , " + val2  + " );"
    density_str = "const value_type density = " + val3 + ";"   
    cfile_str = cfile_str.replace( "RANGE" , range_str )
    cfile_str = cfile_str.replace( "DENSITY" , density_str )
    return cfile_str


def replace_vector( cfile_str , run_params ) :
    cfile_str = cfile_str.replace( "VECTOR" , run_params[0] )
    return cfile_str


def replace_outfile( cfile_str , outfile_name ) :
    outfile_name = "\"" + outfile_name + "\""
    cfile_str = cfile_str.replace( "OUTFILE" , outfile_name )
    return cfile_str

def replace_omp_num_threads( cfile_str , run_params ) :
    omp_str = ""
    if run_params[1] == True :
        num_threads = 1
        if run_params[2] > 1 :
            num_threads = run_params[2]
        omp_str = "omp_set_num_threads( " + str( num_threads ) + " );"
    cfile_str = cfile_str.replace( "OMP_NUM_THREADS" , omp_str )
    return cfile_str

def write_cfile( cfile_str , cfile_name ) :
    f = open( cfile_name , "w" )
    f.write( cfile_str )


def get_system_string( cfile_name , executable_file , run_params ) :
    addCXXFLAGS = ""
    if run_params == True :
        addCXXFLAGS = OPENMP_OPT
    call = CXX + " " + INCLUDEDIR + " " + CXXFLAGS + " " + addCXXFLAGS + " " + LDLIBS + " " + cfile_name + " -o " + executable_file
    return call


    '''replace_range( cfile_str , range_params )
    replace_vector( cfile_str , run_params )
    replace_outfile( cfile_str , run_params )
    
    cfile_name = get_cfile_name( range_params , run_params )
    write_cfile( cfile_str , cfile_name )
    
    executable_file = get_executable_name( range_params , run_params )  
    call_str = get_system_string( cfile_name , executable_file , )'''