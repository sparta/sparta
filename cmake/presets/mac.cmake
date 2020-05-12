# TODOs:
# set(CMAKE_OSX_SYSROOT "" CACHE STRING "")
#
# md5 of executable doesn't match. Not clear if this is due to .a and exe file
# name differences between the old build and this cmake build. Need to verify.

#################### BEGIN SPARTA OPTIONS ####################
#get_property(SPARTA_CXX_COMPILE_FLAGS_DOCSTRING CACHE SPARTA_CXX_COMPILE_FLAGS PROPERTY DOCSTRING)
set(SPARTA_CXX_COMPILE_FLAGS 
    -fPIC -DSPARTA_UNORDERED_MAP
    CACHE 
    STRING 
    "Compiler flags use when building .o files for spa_*"
    )
set(SPARTA_CMAKE_CONFIG_STRING
    mac
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    )
set(BUILD_MPI_TPL
    OFF
    CACHE
    BOOL
    "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
    )
set(PKG_MPI_STUBS
    ON
    CACHE
    BOOL
    "Enable or disable sparta mpi stubs package. Default: OFF.")
#################### END SPARTA OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_C_COMPILER "cc" CACHE STRING "")
set(CMAKE_CXX_COMPILER "c++" CACHE STRING "")
set(CMAKE_CXX_FLAGS "-O" CACHE STRING "")
set(CMAKE_AR "ar" CACHE STRING "")
set(CMAKE_SHARED_LINKER_FLAGS "-fPIC -shared" CACHE STRING "")
#################### END CMAKE OPTIONS ####################

#################### BEGIN KOKKOS OPTIONS ####################
set(Kokkos_ENABLE_OpenMP OFF CACHE STRING "")
#################### END   KOKKOS OPTIONS ####################
