#################### BEGIN SPARTA OPTIONS ####################
#get_property(SPARTA_CXX_COMPILE_FLAGS_DOCSTRING CACHE SPARTA_CXX_COMPILE_FLAGS PROPERTY DOCSTRING)
set(SPARTA_CMAKE_CONFIG_STRING
    chama
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    )
set(SPARTA_CXX_COMPILE_FLAGS 
    -fPIC
    CACHE 
    STRING 
    "Compiler flags use when building .o files for spa_*"
    )
#################### END   SPARTA OPTIONS ####################

#################### BEGIN TPL OPTIONS ####################
set(BUILD_MPI_TPL
    ON
    CACHE
    BOOL
    "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
    )
#################### END   TPL OPTIONS ####################

#################### BEGIN PACKAGE OPTIONS ####################
set(PKG_MPI_STUBS
    OFF
    CACHE
    BOOL
    "Enable or disable sparta mpi stubs package. Default: OFF."
    )
set(PKG_KOKKOS
    OFF
    CACHE
    BOOL
    "Enable or disable sparta kokkos package. Default: OFF."
    )
#################### END   PACKAGE OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
set(CMAKE_CXX_COMPILER "mpic++" CACHE STRING "")
set(CMAKE_CXX_FLAGS "-O2 -xsse4.2 -funroll-loops -fstrict-aliasing" CACHE STRING "")
set(CMAKE_AR "ar" CACHE STRING "")
#ARFLAGS = -rcsv
set(CMAKE_EXE_LINKER_FLAGS "-O -xsse4.2" CACHE STRING "")
set(CMAKE_SHARED_LINKER_FLAGS "-fPIC -shared" CACHE STRING "")
#################### END CMAKE OPTIONS ####################

#################### BEGIN MPI OPTIONS ####################
#################### END   MPI OPTIONS ####################

#################### BEGIN KOKKOS OPTIONS ####################
set(Kokkos_ENABLE_OpenMP OFF CACHE STRING "")
#################### END   KOKKOS OPTIONS ####################
