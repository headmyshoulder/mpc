OPENMP_OPT = "-Xcompiler -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP"

MPCDIR = "-I../../../"
INCLUDEDIR = MPCDIR + " " + "-I$BOOST_ROOT -I$CUDA_ROOT -I$THRUST_ROOT -I$ODEINT_NEW_ROOT -I$THRUST_EXT_ROOT"

CXX = "/usr/local/cuda/bin/nvcc"
CXXFLAGS = "-arch sm_13 -O3"
LDLIBS = "-lgomp -lcudart"

    
# /usr/local/cuda/bin/nvcc -I../../ -I$(BOOST_ROOT) -I$(CUDA_ROOT) -I$(THRUST_ROOT) -I$(ODEINT_NEW_ROOT) -I$(THRUST_EXT_ROOT) -arch sm_13 -O3 -lgomp -lcudart $(>) -o $(<)  
